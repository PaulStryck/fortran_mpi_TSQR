program main
  use  mpi_f08
  use tsqr
  use reduction_tree_m

  implicit none

  integer :: ier
  integer :: nprocesses, thisprocess
  type(reduction_tree) :: tree
  integer, parameter :: n = 3
  integer, parameter :: m = 5

  real(8) :: Q(m,n)
  real(8) :: B(m,n)
  real(8) :: R(n,n)
  call random_number(Q)
  B(:,:) = Q(:,:)

  call qr_local(Q, R)

  call mpi_init(ier)
  call mpi_comm_size(MPI_COMM_WORLD, nprocesses, ier)
  call mpi_comm_rank(MPI_COMM_WORLD, thisprocess, ier)


  call binaryReductionTree(MPI_COMM_WORLD, tree)

  call tsqr_reduction_tree_reduce(Q, R, tree)
  call mpi_bcast(R, 3*3, MPI_REAL8, 0, MPI_COMM_WORLD, ier)

  call pp_arr(matmul(Q,R) - B)

  call mpi_finalize(ier)

  contains
    subroutine qr_local(A, B)
      real(8), intent(inout) :: A(:,:), B(:,:)
      integer :: mA, nA, lda, ier, lwork
      real(8), allocatable :: work(:), tau(:)
      real(8) :: query(1)
      integer :: i,j

      mA = size(A, 1)
      nA = size(A, 2)
      lda = mA

      allocate(tau(min(mA, nA)))

      call dgeqrf(mA, nA, A, lda, tau, query, -1, ier)
      lwork = int(query(1))

      allocate(work(lwork))
      call dgeqrf(mA, nA, A, lda, tau, work, lwork, ier)
      B = 0.0d0
      do j = 1, nA
        forall (i=j:nA) B(j,i) = A(j,i)
      end do

      deallocate(work)

      allocate(work(nA))
      call dorg2r(mA, nA, min(nA, mA), A, mA, tau, work, ier)

      deallocate(work)
      deallocate(tau)

    end subroutine qr_local

end program main
