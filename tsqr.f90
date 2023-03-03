module reduction_tree_m
  use mpi_f08

  implicit none

  type reduction_tree_node
    integer :: l = 0
    integer :: r = 0
    integer :: parent = 0
    logical :: is_root = .false.
  end type

  type reduction_tree
    integer :: parent_of_leaf
    integer :: nnodes = 0
    type(reduction_tree_node), allocatable :: nodes(:)
    type(MPI_COMM) :: comm
  end type reduction_tree

  contains
    subroutine binaryReductionTree(comm, tree)
      type(MPI_COMM), intent(in) :: comm
      type(reduction_tree), intent(out) :: tree

      integer :: r, s
      integer :: c
      integer :: ier

      logical :: is_root

      call MPI_Comm_rank(comm, r, ier);
      call MPI_Comm_size(comm, s, ier);

      tree%comm = comm
      tree%parent_of_leaf = r - (mod(r,2))
      allocate(tree%nodes(s))

      if(r .eq. s-1 .and. mod(r,2) .eq. 0) tree%parent_of_leaf = r-2

      c = 0
      tree%nnodes = 0
      do while((mod(ishft(r, -c),2) .eq. 0) .and. (ishft(1, c) .lt. s))
        if(r+(ishft(1, c)) .lt. s) then
          is_root = .false.
          if((r .eq. 0) .and. (ishft(1, (c+1)) .ge. s)) is_root = .true.

          tree%nnodes = tree%nnodes + 1
          tree%nodes(tree%nnodes) = reduction_tree_node(&
            r,&
            r + ishft(1,c),&
            ishft(ishft(r, -(c+2)), (c+2)),&
            is_root&
          )
        end if
        c = c + 1
      end do
    end subroutine binaryReductionTree

end module reduction_tree_m


module tsqr
  use reduction_tree_m

  implicit none

  type qr_buffer
    real(8), allocatable :: l(:,:)
    real(8), allocatable :: r(:,:)
  end type qr_buffer


  contains
    subroutine tsqr_reduction_tree_reduce(Q, R, tree)
      real(8), intent(out)   :: Q(:,:)
      real(8), intent(inout) :: R(:,:)
      type(reduction_tree), intent(in) :: tree

      integer,            parameter :: TAG = 4733
      type(MPI_DATATYPE), parameter :: dt  = MPI_REAL8

      real(8), allocatable :: Qin(:,:)

      type(qr_buffer), allocatable :: buff_q(:)
      type(qr_buffer), allocatable :: buff_r(:)

      type(mpi_request), allocatable :: recv_req_l(:), recv_req_r(:)
      type(mpi_request) :: leaf_req

      integer :: i
      integer :: ier
      integer :: m, n

      m = size(Q, 1)
      n = size(Q, 2)

      allocate(buff_R(tree%nnodes))
      allocate(buff_Q(tree%nnodes))

      do i = 1, tree%nnodes
        allocate(buff_R(i)%l(n,n))
        allocate(buff_R(i)%r(n,n))

        allocate(buff_Q(i)%l(n,n))
        allocate(buff_Q(i)%r(n,n))
        buff_R(i)%l = 0.0d0
        buff_R(i)%r = 0.0d0

        buff_Q(i)%r = 0.0d0
        buff_Q(i)%l = 0.0d0
      end do

      allocate(recv_req_l(tree%nnodes))
      allocate(recv_req_r(tree%nnodes))

      allocate(Qin(n,n))

      Qin = 0.0d0

      ! recurse tree
      ! generate recv request for all nodes
      do i = 1, tree%nnodes
        ! receive left
        call mpi_Irecv(buff_R(i)%l, n*n, dt, tree%nodes(i)%l, &
          TAG, tree%comm, recv_req_l(i), ier)

        ! receive right
        call mpi_Irecv(buff_R(i)%r, n*n, dt, tree%nodes(i)%r, &
          TAG, tree%comm, recv_req_r(i), ier)
      end do

      ! send data to parent. Might be self
      call mpi_send(R, n*n, dt, tree%parent_of_leaf, &
        TAG, tree%comm, ier)

      ! actual reduce op
      do i = 1, tree%nnodes
        call mpi_wait(recv_req_l(i), MPI_STATUS_IGNORE, ier)
        call mpi_wait(recv_req_r(i), MPI_STATUS_IGNORE, ier)

        call reduce(buff_R(i)%l, buff_R(i)%r, buff_Q(i)%l, buff_Q(i)%r)

        if(.not. tree%nodes(i)%is_root) then
          call mpi_send(buff_R(i)%l, n*n, dt, tree%nodes(i)%parent, &
            TAG, tree%comm, ier)
        end if
      end do

      deallocate(recv_req_l)
      deallocate(recv_req_r)


      ! backpropagate tree
      ! generate recv request for all nodes
      allocate(recv_req_l(tree%nnodes))
      do i = tree%nnodes, 1, -1
        if(tree%nodes(i)%is_root) cycle

        call mpi_Irecv(Qin, n*n, dt, tree%nodes(i)%parent, &
          TAG, tree%comm, recv_req_l(i), ier)
      end do

      ! recv for leaf
      call mpi_Irecv(Qin, n*n, dt, tree%parent_of_leaf, &
        TAG, tree%comm, leaf_req, ier)

      if(tree%nnodes .gt. 0) then
        if(tree%nodes(tree%nnodes)%is_root) then
          ! root process needs to send its Q to two children
          R = buff_R(tree%nnodes)%l
          ! call backpropagate(QintermedL(tree%nnodes,:,:), &
          !   QintermedR(tree%nnodes,:,:))

          call mpi_send(buff_q(tree%nnodes)%l, n*n, dt, tree%nodes(tree%nnodes)%l, &
            TAG, tree%comm, ier)

          call mpi_send(buff_q(tree%nnodes)%r, n*n, dt, tree%nodes(tree%nnodes)%r, &
            TAG, tree%comm, ier)
        end if
      end if

      ! backpropagation
      do i = tree%nnodes, 1, -1
        if(tree%nodes(i)%is_root) cycle

        call mpi_wait(recv_req_l(i), MPI_STATUS_IGNORE, ier)

        call backpropagate(Qin, buff_q(i)%l, buff_q(i)%r)

        call mpi_send(buff_q(i)%l, n*n, dt, tree%nodes(i)%l, &
          TAG, tree%comm, ier)

        call mpi_send(buff_q(i)%r, n*n, dt, tree%nodes(i)%r, &
          TAG, tree%comm, ier)
      end do

      call mpi_wait(leaf_req, MPI_STATUS_IGNORE, ier)

      Q(1:n, 1:n) = Qin(:,:)

      deallocate(buff_R)
      deallocate(buff_Q)
      deallocate(recv_req_l)
    end subroutine tsqr_reduction_tree_reduce

    subroutine reduce(Rr, Rl, Q1, Q2)
      real(8), intent(inout) :: Rr(:,:), Rl(:,:)
      real(8), intent(out) :: Q1(:,:), Q2(:,:)

      real(8), allocatable :: R(:,:), tau(:), work(:)
      real(8) :: query(1)
      integer :: m, n, lda, ier, lwork

      n = size(Rr, 1)
      m = 2*n

      lda = m

      allocate(R(m,n))
      !stack Rr, Rl
      R(1:n,:) = Rr(:,:)
      R(n+1:,:) = Rl(:,:)

      allocate(tau(n))
      call dgeqrf(m, n, R, lda, tau, query, -1, ier)
      lwork = int(query(1))

      allocate(work(lwork))
      call dgeqrf(m, n, R, lda, tau, work, lwork, ier)

      Rr(:,:) = R(1:n,:)

      deallocate(work)

      allocate(work(n))
      call dorg2r(m, n, n, R, m, tau, work, ier)

      Q1(:,:) = R(1:n,:)
      Q2(:,:) = R(n+1:,:)

      deallocate(R)
      deallocate(tau)
      deallocate(work)

    end subroutine reduce

    subroutine backpropagate(Qi, Qa, Qb)
      real(8), intent(in) :: Qi(:,:)
      real(8), intent(inout) :: Qa(:,:), Qb(:,:)

      integer :: n
      real(8), allocatable :: work(:,:)

      n = size(Qi, 2)
      allocate(work(n,n))

      ! Qa = matmul(Qa, Qi)
      call dgemm('N', 'N', n, n, n, 1.0d0, Qa, n, Qi, n, 0.0d0, work, n)
      Qa = work

      ! Qb = matmul(Qb, Qi)
      call dgemm('N', 'N', n, n, n, 1.0d0, Qb, n, Qi, n, 0.0d0, work, n)
      Qb = work

      deallocate(work)
    end subroutine backpropagate

    subroutine pp_arr(Grid)
      real(8), intent(in) :: Grid(:,:)

      integer :: i,j, n,m
      n = size(Grid, 1)
      m = size(Grid, 2)

      write( * , "(*(g0.4))" ) ( (Grid(i,j)," ",j=1,m), new_line("A"), i=1,n)
    end subroutine

end module tsqr
