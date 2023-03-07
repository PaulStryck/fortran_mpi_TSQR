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

    !=========
    !  Reverse engineer a reduction tree used by MPI_Reduce
    !=========
    subroutine mpiReductionTree(comm, tree, root)
      type(MPI_COMM), intent(in) :: comm
      type(reduction_tree), intent(out) :: tree
      integer, intent(in) :: root

      integer, parameter :: TAG = 1107
      integer :: rnk, sze ! rank and size of comm

      integer :: cnt      ! TAG + cnt will be used as mpi_tag to differentiate messages

      integer :: d1(2), d2(2) ! data buffers during reduction

      type(MPI_OP) :: op

      type(MPI_REQUEST), allocatable :: recv_reqs(:)


      cnt = 0

      tree%comm = comm

      call mpi_comm_rank(comm, rnk)
      call mpi_comm_size(comm, sze)

      allocate(recv_reqs(1))    ! this will frequently be extended

      ! references to members of tree%nodes(i) will be used as buffers
      !   thus they MUST NEVER change their memory address.
      ! allocate upper limit so it will never move in memory
      allocate(tree%nodes(sze))
      tree%nnodes = 0 ! keep track of actually used nodes


      call mpi_irecv(tree%parent_of_leaf, 1, MPI_INT, MPI_ANY_SOURCE,&
        TAG, comm, recv_reqs(Ubound(recv_reqs,1)))

      d1 = [rnk, cnt]
      d2 = [rnk, cnt]

      call mpi_op_create(insert_node, .true., op)
      call mpi_reduce(d1, d2, 2, MPI_INT, op, root, comm)

      if(rnk .eq. root) then
        call mpi_send(rnk, 1, MPI_INT, rnk, TAG+cnt, comm)
        if(sze .gt. 1) tree%nodes(tree%nnodes)%is_root = .true.
      end if

      call mpi_op_free(op)
      call mpi_waitall(size(recv_reqs), recv_reqs, MPI_STATUSES_IGNORE)

      contains

      ! will need to reference local tree, local recv_reqs, local cnt.
      !   MUST thus be bound to the outer subroutine to keep
      !   tree, recv_reqs, cnt in its scope
      subroutine insert_node(inv, outv, l, t)
        use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

        type(c_ptr), value :: inv, outv
        integer(c_int) :: l
        type(MPI_Datatype) :: t

        integer, pointer :: inv_i(:), outv_i(:)

        call c_f_pointer(inv,  inv_i,  [l])
        call c_f_pointer(outv, outv_i, [l])

        tree%nnodes = 1 + tree%nnodes
        tree%nodes(tree%nnodes)%l = inv_i(1)
        tree%nodes(tree%nnodes)%r = outv_i(1)
        tree%nodes(tree%nnodes)%parent = -1
        tree%nodes(tree%nnodes)%is_root = .false.

        call mpi_send(rnk, 1, MPI_INT, inv_i(1),  TAG+inv_i(2),  comm)
        call mpi_send(rnk, 1, MPI_INT, outv_i(1), TAG+outv_i(2), comm)

        recv_reqs = [recv_reqs, MPI_REQUEST(-1)]
        cnt = 1 + cnt

        call mpi_irecv(tree%nodes(tree%nnodes)%parent, 1, MPI_INT,&
          MPI_ANY_SOURCE, TAG+cnt, comm, recv_reqs(Ubound(recv_reqs, 1)))

        outv_i(1) = rnk
        outv_i(2) = cnt

      end subroutine

    end subroutine mpiReductionTree


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

  contains

    subroutine tsqr_qr(Q, R, tree)
      real(8), intent(inout) :: Q(:,:)
      real(8), intent(out)   :: R(:,:)
      type(reduction_tree), intent(in) :: tree

      real(8), allocatable :: tau(:)
      real(8), allocatable :: R_p(:), Qred_p(:), Qred(:,:)
      integer :: m, n, s

      integer :: lwork1, lwork2, ier
      real(8) :: query1(1), query2(1,1)
      real(8), allocatable :: work1(:), work2(:,:)


      m = size(Q, 1)
      n = size(Q, 2)
      s = (n*n+n) / 2

      allocate(tau(n))
      allocate(R_p(s))
      allocate(Qred_p(s))

      Qred_p = 0.0d0
      R_p = 0.0d0
      tau = 0.0d0

      ! compute local Q, R
      ! Q is overwritten by householder reflectors
      ! R_p contains packed upper triag R
      call qr_local(Q, R_p, tau)

      ! compute local Qs for global R
      call tsqr_reduction_tree_reduce(Qred_p, R_p, n, tree)
      call dtpttr('U', n, R_p, R, n, ier)

      ! apply householder reflectors
      if(.true.) then
        ! if enough memory for Qred use super fast dormqr
        allocate(Qred(m,n))
        Qred = 0.0d0
        call dtpttr('U', n, Qred_p, Qred(1:n, 1:n), n, ier)

        call dormqr('L', 'N', m, n, n, Q, m, tau, Qred, m, query1, -1, ier)
        lwork1 = int(query1(1))
        allocate(work1(lwork1))
        call dormqr('L', 'N', m, n, n, Q, m, tau, Qred, m, work1, lwork1, ier)
        deallocate(work1)

        Q(:,:) = Qred(:,:)

        deallocate(Qred)
      else
        ! if not enough memory, use slower inplace
        !   speed depends on compiler optimization
        !   but usally 4-5x slower than above method
        call dtpttr('U', n, Qred_p, Q, m, ier)

        call dormqr_inplace('L', 'N', m, n, n, Q, m, tau, query1, -1, query2, -1, ier)
        lwork1 = int(query1(1))
        lwork2 = int(query2(1,1))

        allocate(work1(lwork1))
        allocate(work2(m,lwork2))
        call dormqr_inplace('L', 'N', m, n, n, Q, m, tau, work1, lwork1, work2, lwork2, ier)
        deallocate(work1)
        deallocate(work2)
      end if

      deallocate(Qred_p)
      deallocate(R_p)
      deallocate(tau)

    end subroutine tsqr_qr


    subroutine dormqr_inplace(SIDE, TRANS, M, N, K, A, LDA, TAU, &
        WORK, LWORK, WORK2, LWORK2, INFO)
      ! .. Scalar Arguments ..
      character :: side, trans
      integer   :: info, k, lda, lwork, lwork2, m, n

      ! .. Array Arguments ..
      real(8) :: A(lda, *), tau(*), work(*),  WORK2(lda, lwork2)

      ! ====================================
      ! .. Parameters ..
      integer, parameter :: nbmax = 2
      integer, parameter :: ldt = nbmax + 1
      integer, parameter :: tsize = ldt * nbmax

      ! .. Local Scalars ..
      logical :: left, lquery, notran
      integer :: i, i1, i2, i3, ib, ic, iwt, jc, ldwork, &
        lwkopt, mi, nb, nbmin, ni, nq, nw, ii

      ! .. External Functions ..
      logical LSAME
      integer ILAENV
      external lsame, ilaenv

      ! .. External Subroutines ..
      external dlarfb, dlarft, dorm2r, xerbla

      ! .. intrinsic functions ..
      intrinsic max, min

      ! .. Executable Statements ..
      ! Test the input arguments
      info = 0
      left = lsame(side, 'L')
      notran = lsame(trans, 'N')
      lquery = lwork .eq. -1

      IF( left ) THEN
        nq = m
        nw = max( 1, n )
      END IF

      IF( .NOT.left ) THEN
        info = -1
      ELSE IF( .NOT.notran ) THEN
        info = -2
      ELSE IF( m.LT.0 ) THEN
        info = -3
      ELSE IF( n.LT.0 ) THEN
        info = -4
      ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
        info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
        info = -7
      ELSE IF( lwork.LT.nw .AND. .NOT.lquery ) THEN
        info = -12
      END IF

      IF( info.EQ.0 ) THEN
        ! Compute the workspace requirements
        nb = min( nbmax, ilaenv( 1, 'DORMQR', side // trans, m, n, k, -1 ) )
        lwkopt = nw*nb + tsize
        work( 1 )    = lwkopt
        work2( 1,1 ) = nb
      END IF

      IF( info.NE.0 ) THEN
        CALL xerbla( 'DORMQR', -info )
        RETURN
      ELSE IF( lquery ) THEN
        RETURN
      END IF

      ! quick return if possible
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 ) THEN
        work( 1 ) = 1
        work2( 1,1 ) = 1
        RETURN
      END IF

      nbmin = 2
      ldwork = nw
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
        IF( lwork.LT.lwkopt .OR. lwork2.lt.nb) THEN
          nb = min(lwork, (lwork-tsize) / ldwork)
          nbmin = max( 2, ilaenv( 2, 'DORMQR', side // trans, m, n, k, -1 ) )
        END IF
      END IF

      IF( nb.LT.nbmin .OR. nb.GE.k ) THEN
        ! Use unblocked code
        write(*,*) "TODO nb: ", nb, nbmin, k
        ! CALL dorm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work, iinfo )
      ELSE
        ! Use blocked code
        iwt = 1 + nw*nb

        i1 = ( ( k-1 ) / nb )*nb + 1
        i2 = 1
        i3 = -nb

        ni = n
        jc = 1

        DO i = i1, i2, i3
          ib = min( nb, k-i+1)

          call dlarft( 'Forward', 'Columnwise', nq-i+1, ib, a( i, i ),&
            lda, tau(i), work(iwt), ldt)
          ! copy block reflector to be used by dlarfb
          ! allocate(cc(lda, ib))
          work2 = 0.0d0
          work2 = a(1:lda, i:i+k)
          ! zero out used reflector elements
          do ii = i,i+ib-1
              a( ii+1:m,ii ) = 0.0d0
          end do

          ! H is applied to C(i:m,1:n)
          mi = m - i + 1
          ic = i

          ! Apply H
          CALL dlarfb( side, trans, 'Forward', 'Columnwise', mi, n-ic+1, &
            ib, work2(i,1), lda, work( iwt ), ldt, &
            a( ic, ic ), lda, work, ldwork )
        END DO
      END IF
      work(1) = lwkopt
      work2( 1,1) = nb
      RETURN
    end subroutine

    subroutine qr_local(A, R_p, tau)
      real(8), intent(inout) :: A(:,:), R_p(:)
      real(8), intent(out) :: tau(:)

      integer :: mA, nA, lda, ier, lwork
      real(8), allocatable :: work(:)
      real(8) :: query(1)

      mA = size(A, 1)
      nA = size(A, 2)
      lda = mA

      call dgeqrf(mA, nA, A, lda, tau, query, -1, ier)
      lwork = int(query(1))

      allocate(work(lwork))
      call dgeqrf(mA, nA, A, lda, tau, work, lwork, ier)
      deallocate(work)

      R_p = 0.0d0
      call dtrttp('U', nA, A(1:nA, 1:nA), nA, R_p, ier)


    end subroutine qr_local


    subroutine tsqr_reduction_tree_reduce(Q_p, R_p, n, tree)
      real(8), intent(out)   :: Q_p(:)
      real(8), intent(inout) :: R_p(:)
      integer, intent(in)    :: n
      type(reduction_tree), intent(in) :: tree

      integer,            parameter :: TAG = 4733
      type(MPI_DATATYPE), parameter :: dt  = MPI_REAL8


      type(packed_buffer), allocatable :: buff_Q(:)
      type(packed_buffer), allocatable :: buff_R(:)

      real(8), allocatable :: Q_leaf_buffer(:)

      type(mpi_request), allocatable :: recv_req_l(:), recv_req_r(:)
      type(mpi_request) :: leaf_req

      integer :: i
      integer :: ier
      integer :: s

      s = (n * n + n) / 2 ! use integer division, always even

      allocate(Q_leaf_buffer(s))
      Q_leaf_buffer = 0.0d0

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

      Q_p  = 0.0d0

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
      call mpi_send(R_p, s, dt, tree%parent_of_leaf, &
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

        call mpi_Irecv(buff_R(i)%l, s, dt, tree%nodes(i)%parent, &
          TAG, tree%comm, recv_req_l(i), ier)
      end do

      ! recv for leaf
      call mpi_Irecv(Q_leaf_buffer, s, dt, tree%parent_of_leaf, &
        TAG, tree%comm, leaf_req, ier)

      if(tree%nnodes .gt. 0) then
        if(tree%nodes(tree%nnodes)%is_root) then
          ! root contains global R
          R_p = buff_R(tree%nnodes)%l
          ! call mpi_bcast(R_p, s, dt, 0, tree%comm, ier)

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

        call backpropagate(buff_R(i)%l, buff_Q(i)%l, buff_Q(i)%r, n)

        call mpi_send(buff_Q(i)%l, s, dt, tree%nodes(i)%l, &
          TAG, tree%comm, ier)

        call mpi_send(buff_Q(i)%r, s, dt, tree%nodes(i)%r, &
          TAG, tree%comm, ier)
      end do

      call mpi_wait(leaf_req, MPI_STATUS_IGNORE, ier)
      call mpi_bcast(R_p, s, MPI_REAL8, 0, tree%comm, ier)

      ! Q_leaf_buffer has asynchronous attribute,
      !  but has been awaitet above.
      ! To satisfy the compiler, values must be copied to matrix
      !  without asynchronous attribute
      Q_p = Q_leaf_buffer


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

      m = 2*n

      lda = m

      allocate(R(m,n))

      ! stack packed Rl, Rr
      R = 0.0d0
      call dtpttr('U', n, Rl, R, lda, ier)
      call dtpttr('U', n, Rr, R(n+1:, 1:n), n, ier)

      allocate(tau(n))
      call dgeqrf(m, n, R, lda, tau, query, -1, ier)
      lwork = int(query(1))

      allocate(work(lwork))
      call dgeqrf(m, n, R, lda, tau, work, lwork, ier)

      ! recover new upper triag as packed
      call dtrttp('U', n, R, m, Rl, ier)

      deallocate(work)

      ! call dorg2r(m, n, n, R, m, tau, work, ier) ! unblocked
      call dorgqr(m, n, n, R, m, tau, query, -1, ier)
      lwork = int(query(1))
      allocate(work(lwork))
      call dorgqr(m, n, n, R, m, tau, work, lwork, ier) ! blocked

      call dtrttp('U', n, R, m, Q1, ier)
      call dtrttp('U', n, R(n+1:,1:n), n, Q2, ier)

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
      Grid = 0.0d0
      call dtpttr('U', n, P, Grid, n, ier)
      call pp_arr(Grid)
      deallocate(Grid)
    end subroutine
    subroutine pp_arr(Grid)
      real(8), intent(in) :: Grid(:,:)

      integer :: i,j, n,m
      n = size(Grid, 1)
      m = size(Grid, 2)

      write( * , "(*(g0.4))" ) ( (Grid(i,j)," ",j=1,m), new_line("A"), i=1,n)
    end subroutine

end module tsqr
