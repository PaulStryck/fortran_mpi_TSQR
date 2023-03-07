program main
  use  mpi_f08
  use gsqr
  use tsqr
  use reduction_tree_m

  implicit none

  integer, parameter :: n = 500
  integer, parameter :: m = 100000

  type(reduction_tree) :: tree

  real(8) :: A(m,n)
  real(8) :: Q(m,n)
  real(8) :: R(n,n)
  real(8) :: x

  logical :: onmain
  integer :: nprocesses, thisprocess
  integer :: ier

  real(8) :: tic, toc


  call mpi_init(ier)
  call mpi_comm_size(MPI_COMM_WORLD, nprocesses, ier)
  call mpi_comm_rank(MPI_COMM_WORLD, thisprocess, ier)

  onmain = thisprocess.eq.0

  ! generate reduction tree
  call mpiReductionTree(MPI_COMM_WORLD, tree, 0)

  Q = 1.0d0
  R  = 0.0d0

  call random_number(Q)

  ! backup Q for later verification
  A(:,:) = Q(:,:)

  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(tic)
  call tsqr_qr(Q, R, tree)
  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(toc)
  if(onmain) write(*,*) "tsqr", toc - tic

  x = sum(abs(matmul(Q,R) - A)) / (m*n)
  if (x .gt. 1d-13) write(*,*) "err on proc ", thisprocess, x

  Q = 0.0d0
  R = 0.0d0
  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(tic)
  call gsqr_mat(A, Q, R)
  call mpi_barrier(MPI_COMM_WORLD); call cpu_time(toc)
  if(onmain) write(*,*) "gs", toc - tic


  call mpi_finalize(ier)

end program main
