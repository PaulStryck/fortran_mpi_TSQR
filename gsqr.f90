module gsqr
  use mpi_f08

contains
  subroutine parallel_sum(x)
    real(8), intent(inout) :: x
    real(8) :: tmp

    integer :: ier

    tmp = 0.0d0

    call mpi_allreduce(x, tmp, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ier)
    x = tmp

  end subroutine parallel_sum

  subroutine gsqr_mat(A, Q, R)

    REAL(8), INTENT(IN) :: A(:,:)
    REAL(8), INTENT(OUT) :: Q(:,:), R(:,:)

    INTEGER :: nA, mA, i, j, k
    REAL(8), ALLOCATABLE :: v(:)

    mA = size(A,1)
    nA = size(A,2)

    allocate(v(mA))

    Q(:,:) = 0.0d0
    R(:,:) = 0.0d0
    v(:) = 0.0d0

    DO j=1, nA
      v(:) = A(:,j)
      DO k=1, 1
        DO i=1,j-1
          R(i,j) = DOT_PRODUCT(Q(:,i),v(:))
          CALL parallel_sum(R(i,j))
          v(:) = v(:) - R(i,j)*Q(:,i)
        END DO
      END DO

      R(j,j) = DOT_PRODUCT(v(:),v(:))
      CALL parallel_sum(R(j,j))
      R(j,j) = SQRT(R(j,j))
      Q(:,j) = v(:)/R(j,j)
    END DO

    deallocate(v)
  end subroutine gsqr_mat
end module gsqr
