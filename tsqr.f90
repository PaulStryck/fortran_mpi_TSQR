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

    !=========
    !  Construct a binary reduction tree
    !  Will only work if size of comm is a power of 2
    !=========
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

  ! only ever needs to hold triangular matrices.
  !  thus use packed storage
  type packed_buffer
    real(8), allocatable :: l(:)
    real(8), allocatable :: r(:)
  end type packed_buffer

  type dense_buffer
    real(8), allocatable :: l(:,:)
    real(8), allocatable :: r(:,:)
  end type dense_buffer

  integer :: rnk


  contains
    subroutine tsqr_reduction_tree_reduce(Q, R, tree)
      real(8), intent(out)   :: Q(:,:)
      real(8), intent(inout) :: R(:,:)
      type(reduction_tree), intent(in) :: tree

      integer,            parameter :: TAG = 4733
      type(MPI_DATATYPE), parameter :: dt  = MPI_REAL8

      real(8), allocatable :: Qin(:) ! packed upper triag
      real(8), allocatable :: R_packed(:)

      type(packed_buffer), allocatable :: buff_Q(:)
      type(packed_buffer), allocatable :: buff_R(:)

      type(mpi_request), allocatable :: recv_req_l(:), recv_req_r(:)
      type(mpi_request) :: leaf_req

      integer :: i, j, k, l
      integer :: ier
      integer :: m, n
      integer :: s

      call mpi_comm_rank(MPI_COMM_WORLD, rnk)
      m = size(Q, 1)
      n = size(Q, 2)
      s = (n * n + n) / 2 ! use integer division, always even

      allocate(buff_R(tree%nnodes))
      allocate(buff_Q(tree%nnodes))

      do i = 1, tree%nnodes
        allocate(buff_R(i)%l(s))
        allocate(buff_R(i)%r(s))

        allocate(buff_Q(i)%l(s))
        allocate(buff_Q(i)%r(s))

        buff_R(i)%l = 0.0d0
        buff_R(i)%r = 0.0d0

        buff_Q(i)%l = 0.0d0
        buff_Q(i)%r = 0.0d0
      end do

      allocate(recv_req_l(tree%nnodes))
      allocate(recv_req_r(tree%nnodes))

      allocate(Qin(s))

      Q    = 0.0d0
      Qin  = 0.0d0

      ! forall (i=1:min(m,n)) Qin(i,i) = 1.0d0

      allocate(R_packed(s))
      l = 1
      do j = 1, n
        do i = 1, j
          R_packed(l) = R(i,j)
          l = l + 1
        end do
      end do

      ! recurse tree
      ! generate recv request for all nodes
      do i = 1, tree%nnodes
        ! receive left
        call mpi_Irecv(buff_R(i)%l, s, dt, tree%nodes(i)%l, &
          TAG, tree%comm, recv_req_l(i), ier)

        ! receive right
        call mpi_Irecv(buff_R(i)%r, s, dt, tree%nodes(i)%r, &
          TAG, tree%comm, recv_req_r(i), ier)
      end do

      ! send data to parent. Might be self
      call mpi_send(R_packed, s, dt, tree%parent_of_leaf, &
        TAG, tree%comm, ier)

      ! actual reduce op
      do i = 1, tree%nnodes
        call mpi_wait(recv_req_l(i), MPI_STATUS_IGNORE, ier)
        call mpi_wait(recv_req_r(i), MPI_STATUS_IGNORE, ier)

        call reduce(buff_R(i)%l, buff_R(i)%r, buff_Q(i)%l, buff_Q(i)%r, n)

        if(.not. tree%nodes(i)%is_root) then
          call mpi_send(buff_R(i)%l, s, dt, tree%nodes(i)%parent, &
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

        call mpi_Irecv(Qin, s, dt, tree%nodes(i)%parent, &
          TAG, tree%comm, recv_req_l(i), ier)
      end do

      ! recv for leaf
      call mpi_Irecv(Qin, s, dt, tree%parent_of_leaf, &
        TAG, tree%comm, leaf_req, ier)

      if(tree%nnodes .gt. 0) then
        if(tree%nodes(tree%nnodes)%is_root) then
          ! root contains global R
          R_packed = buff_R(tree%nnodes)%l
          ! call mpi_bcast(R_packed, s, dt, 0, tree%comm, ier)

          call mpi_send(buff_Q(tree%nnodes)%l, s, dt, tree%nodes(tree%nnodes)%l, &
            TAG, tree%comm, ier)

          call mpi_send(buff_Q(tree%nnodes)%r, s, dt, tree%nodes(tree%nnodes)%r, &
            TAG, tree%comm, ier)
        end if
      end if

      ! backpropagation
      do i = tree%nnodes, 1, -1
        if(tree%nodes(i)%is_root) cycle

        call mpi_wait(recv_req_l(i), MPI_STATUS_IGNORE, ier)

        call backpropagate(Qin, buff_Q(i)%l, buff_Q(i)%r, n)

        call mpi_send(buff_Q(i)%l, s, dt, tree%nodes(i)%l, &
          TAG, tree%comm, ier)

        call mpi_send(buff_Q(i)%r, s, dt, tree%nodes(i)%r, &
          TAG, tree%comm, ier)
      end do

      call mpi_wait(leaf_req, MPI_STATUS_IGNORE, ier)

      Q = 0.0d0
      call dtpttr('U', n, Qin, Q(1:n, 1:n), n, ier)
      call dtpttr('U', n, R_packed, R, n, ier)

      deallocate(buff_R)
      deallocate(buff_Q)
      deallocate(recv_req_l)
    end subroutine tsqr_reduction_tree_reduce

    subroutine reduce(Rl, Rr, Q1, Q2, n)
      real(8), intent(inout) :: Rl(:) ! packed upper triag
      real(8), intent(in)    :: Rr(:) ! packed upper triag
      real(8), intent(out)   :: Q1(:), Q2(:) ! packed upper triag

      integer, intent(in)    :: n

      real(8), allocatable :: R(:,:), tau(:), work(:)
      real(8) :: query(1)
      integer :: m, lda, ier, lwork
      integer :: i, j, l

      m = 2*n

      lda = m

      allocate(R(m,n))

      ! stack packed Rl, Rr
      R = 0.0d0
      l = 1
      do j = 1, n
        do i = 1, j
          R(i  ,j)  = Rl(l)
          R(i+n, j) = Rr(l)
          l = l + 1
        end do
      end do

      allocate(tau(n))
      call dgeqrf(m, n, R, lda, tau, query, -1, ier)
      lwork = int(query(1))

      allocate(work(lwork))
      call dgeqrf(m, n, R, lda, tau, work, lwork, ier)

      ! recover new upper triag as packed
      l = 1
      do j = 1, n
        do i = 1, j
          Rl(l) = R(i,j)
          l = l+1
        end do
      end do

      deallocate(work)

      allocate(work(n))
      call dorg2r(m, n, n, R, m, tau, work, ier)

      l = 1
      do j = 1, n
        do i = 1, j
          Q1(l) = R(i  ,j)
          Q2(l) = R(i+n,j)
          l = l+1
        end do
      end do

      deallocate(R)
      deallocate(tau)
      deallocate(work)

    end subroutine reduce

    subroutine backpropagate(Qi, Qa, Qb, n)
      real(8), intent(in)    :: Qi(:)        ! packed upper triag
      real(8), intent(inout) :: Qa(:), Qb(:) ! packed upper triag
      integer, intent(in) :: n

      integer :: ier

      real(8), allocatable :: A(:,:), B(:,:)

      allocate(A(n,n))
      allocate(B(n,n))
      A = 0.0d0
      B = 0.0d0

      ! expand Qi from packed to square
      call dtpttr('U', n, Qi, A, n, ier)

      ! Qa = matmul(Qa, Qi)
      B = 0.0d0
      call dtpttr('U', n, Qa, B, n, ier)
      call dtrmm('R', 'U', 'N', 'N', n, n, 1.0d0, A, n, B, n)
      call dtrttp('U', n, B, n, Qa, ier)
      ! call dgemm('N', 'N', n, n, n, 1.0d0, Qa, n, Qi, n, 0.0d0, work, n)
      ! Qa = work

      ! Qb = matmul(Qb, Qi)
      B = 0.0d0
      call dtpttr('U', n, Qb, B, n, ier)
      call dtrmm('R', 'U', 'N', 'N', n, n, 1.0d0, A, n, B, n)
      call dtrttp('U', n, B, n, Qb, ier)
      ! call dgemm('N', 'N', n, n, n, 1.0d0, Qb, n, Qi, n, 0.0d0, work, n)
      ! Qb = work

      deallocate(A)
      deallocate(B)
    end subroutine backpropagate

    subroutine pp_arr_p(P, n)
      real(8), intent(in) :: P(:)
      integer, intent(in) :: n
      real(8), allocatable :: Grid(:,:)
      integer :: ier

      allocate(Grid(n,n))
      call dtpttr('U', n, P, Grid, n, ier)
      call pp_arr(Grid)
    end subroutine
    subroutine pp_arr(Grid)
      real(8), intent(in) :: Grid(:,:)

      integer :: i,j, n,m
      n = size(Grid, 1)
      m = size(Grid, 2)

      write( * , "(*(g0.4))" ) ( (Grid(i,j)," ",j=1,m), new_line("A"), i=1,n)
    end subroutine

end module tsqr
