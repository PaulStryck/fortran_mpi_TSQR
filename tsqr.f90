module reduction_tree_m
  use mpi_f08
  implicit none

  type reduction_tree_node
    integer :: s1 = 0
    integer :: s2 = 0
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


  contains
    subroutine tsqr_reduction_tree_reduce(Q, R, tree)
      real(8), intent(inout) :: Q(:,:)
      real(8), intent(inout) :: R(:,:)
      type(reduction_tree), intent(in) :: tree

      real(8), allocatable :: R1(:,:), R2(:,:), Qin(:,:)

      real(8), allocatable :: QintermedA(:,:,:), QintermedB(:,:,:)
      TYPE(MPI_Request), allocatable :: recv_req1(:), recv_req2(:)
      type(mpi_request) :: leaf_req
      integer, parameter :: TAG = 4733
      type(MPI_DATATYPE) :: dt = MPI_REAL8

      integer :: i
      integer :: ier
      integer :: m, n
      integer :: rnk

      m = size(Q, 1)
      n = size(Q, 2)

      allocate(QintermedA(tree%nnodes, n, n))
      allocate(QintermedB(tree%nnodes, n, n))
      QintermedA = 0.0d0
      QintermedB = 0.0d0

      allocate(recv_req1(tree%nnodes))
      allocate(recv_req2(tree%nnodes))

      call mpi_comm_rank(MPI_COMM_WORLD, rnk, ier)

      allocate(R1(n,n), R2(n,n), Qin(n,n))

      R1  = 0.0d0
      R2  = 0.0d0
      Qin = 0.0d0

      ! recurse tree
      ! generate recv request for all nodes
      do i = 1, tree%nnodes
        ! receive left
        call mpi_Irecv(R1, n*n, dt, tree%nodes(i)%s1, TAG, &
          tree%comm, recv_req1(i), ier)

        ! receive right
        call mpi_Irecv(R2, n*n, dt, tree%nodes(i)%s2, TAG, &
          tree%comm, recv_req2(i), ier)
      end do

      ! send data to parent. Might be self
      call mpi_send(R, n*n, dt, tree%parent_of_leaf, TAG, &
        tree%comm, ier)

      ! actual reduce op
      do i = 1, tree%nnodes
        call mpi_wait(recv_req1(i), MPI_STATUS_IGNORE, ier)
        call mpi_wait(recv_req2(i), MPI_STATUS_IGNORE, ier)

        call reduce(R1, R2, QintermedA(i,:,:), QintermedB(i,:,:))

        if(.not. tree%nodes(i)%is_root) then
          call mpi_send(R1, n*n, dt, tree%nodes(i)%parent, &
            TAG, tree%comm, ier)
        end if
      end do

      deallocate(recv_req1)
      deallocate(recv_req2)

      allocate(recv_req1(tree%nnodes))
      ! backpropagate tree
      ! generate recv request for all nodes
      do i = tree%nnodes, 1, -1
        if(tree%nodes(i)%is_root) cycle

        call mpi_Irecv(Qin, n*n, dt, tree%nodes(i)%parent, &
          TAG, tree%comm, recv_req1(i), ier)
      end do

      ! recv for leaf
      call mpi_Irecv(Qin, n*n, dt, tree%parent_of_leaf, &
        TAG, tree%comm, leaf_req, ier)

      if(tree%nnodes .gt. 0) then
        if(tree%nodes(tree%nnodes)%is_root) then
          ! root process needs to send its Q to two children
          R = R1
          ! call backpropagate(QintermedA(tree%nnodes,:,:), &
          !   QintermedB(tree%nnodes,:,:))

          call mpi_send(QintermedA(tree%nnodes,:,:), n*n, dt, tree%nodes(tree%nnodes)%s1, &
            TAG, tree%comm, ier)

          call mpi_send(QintermedB(tree%nnodes,:,:), n*n, dt, tree%nodes(tree%nnodes)%s2, &
            TAG, tree%comm, ier)
        end if
      end if

      ! backpropagation
      do i = tree%nnodes, 1, -1
        if(tree%nodes(i)%is_root) cycle

        call mpi_wait(recv_req1(i), MPI_STATUS_IGNORE, ier)

        call backpropagate(Qin, QintermedA(i,:,:), QintermedB(i,:,:))
        call mpi_send(QintermedA(i,:,:), n*n, dt, tree%nodes(i)%s1, &
          TAG, tree%comm, ier)

        call mpi_send(QintermedB(i,:,:), n*n, dt, tree%nodes(i)%s2, &
          TAG, tree%comm, ier)
      end do

      call mpi_wait(leaf_req, MPI_STATUS_IGNORE, ier)
      Q = matmul(Q, Qin)
    end subroutine tsqr_reduction_tree_reduce

    subroutine reduce(R1, R2, Q1, Q2)
      real(8), intent(inout) :: R1(:,:), R2(:,:)
      real(8), intent(out) :: Q1(:,:), Q2(:,:)

      real(8), allocatable :: R(:,:), tau(:), work(:)
      real(8) :: query(1)
      integer :: m, n, lda, ier, lwork

      n = size(R1, 1)
      m = 2*n

      lda = m

      allocate(R(m,n))
      !stack R1, R2
      R(1:n,:) = R1(:,:)
      R(n+1:,:) = R2(:,:)

      allocate(tau(n))
      call dgeqrf(m, n, R, lda, tau, query, -1, ier)
      lwork = int(query(1))

      allocate(work(lwork))
      call dgeqrf(m, n, R, lda, tau, work, lwork, ier)

      R1(:,:) = R(1:n,:)

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

      Qa = matmul(Qa, Qi)
      Qb = matmul(Qb, Qi)
    end subroutine backpropagate

    subroutine pp_arr(Grid)
      real(8), intent(in) :: Grid(:,:)

      integer :: i,j, n,m
      n = size(Grid, 1)
      m = size(Grid, 2)

      write( * , "(*(g0.4))" ) ( (Grid(i,j)," ",j=1,m), new_line("A"), i=1,n)
    end subroutine

end module tsqr
