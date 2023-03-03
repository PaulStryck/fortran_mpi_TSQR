program main
  use  mpi_f08
  use tsqr
  use reduction_tree_m

  implicit none

  integer :: nprocesses, thisprocess


  integer :: ier
  type(reduction_tree) :: tree
  integer, parameter :: n = 5
  integer, parameter :: m = 10

  real(8) :: B(m,n)
  real(8) :: Qred(m,n)
  real(8) :: Q(m,n)
  real(8) :: R(n,n)
  real(8) :: tau(n)
  real(8) :: x
  real(8) :: tic, toc, overall = 0.0d0

  real(8) :: query(1)
  real(8), allocatable :: work(:)
  integer :: lwork

  logical :: onmain

  call mpi_init(ier)
  call mpi_comm_size(MPI_COMM_WORLD, nprocesses, ier)
  call mpi_comm_rank(MPI_COMM_WORLD, thisprocess, ier)

  onmain = thisprocess.eq.0

  ! generate reduction tree
  call binaryReductionTree(MPI_COMM_WORLD, tree)

  Q = 0.0d0
  R = 0.0d0
  Qred = 0.0d0
  tau = 0.0d0

  call random_number(Q)

  ! backup Q for later verification
  B(:,:) = Q(:,:)

  ! local QR on each pocessor
  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(tic)

  call qr_local(Q, R, tau)

  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(toc)
  if(onmain) write(*,*) "initial qr", toc-tic
  overall = overall + toc-tic


  ! global R reduction
  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(tic)

  call tsqr_reduction_tree_reduce(Qred, R, tree)
  call mpi_bcast(R, n*n, MPI_REAL8, 0, MPI_COMM_WORLD, ier)

  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(toc)
  if(onmain) write(*,*) "tsqr", toc - tic
  overall = overall + toc - tic

  ! apply initial householder transformations to Q
  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(tic)
  call dormqr('L', 'N', m, n, n, Q, m, tau, Qred, m, query, -1, ier)
  lwork = int(query(1))
  allocate(work(lwork))
  call dormqr('L', 'N', m, n, n, Q, m, tau, Qred, m, work, lwork, ier)
  deallocate(work)
  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(toc)
  if(onmain) write(*,*) "Q", toc - tic
  overall = overall + toc-tic

  !Qred * R = B


  x = sum(abs(matmul(Qred,R) - B)) / (m*n)
  if (x .gt. 1d-13) write(*,*) "err on proc ", thisprocess, x


  if(onmain) write(*,*) "overall", overall



  call mpi_finalize(ier)

  contains
    subroutine qr_local(A, R, tau)
      real(8), intent(inout) :: A(:,:), R(:,:)
      real(8), intent(out) :: tau(:)
      integer :: mA, nA, lda, ier, lwork
      real(8), allocatable :: work(:)
      real(8) :: query(1)
      integer :: i,j

      mA = size(A, 1)
      nA = size(A, 2)
      lda = mA

      call dgeqrf(mA, nA, A, lda, tau, query, -1, ier)
      lwork = int(query(1))

      allocate(work(lwork))
      call dgeqrf(mA, nA, A, lda, tau, work, lwork, ier)
      R = 0.0d0
      do j = 1, nA
        forall (i=j:nA) R(j,i) = A(j,i)
      end do

      deallocate(work)

    end subroutine qr_local

end program main
