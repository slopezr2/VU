!###############################################################################
!
! numerical tools : linear algebra
!
! Wrapper to LAPack and/or LAPack95/MKL interfaces.
!
! For LAPack interface, see:
!   http://www.netlib.org/lapack/
!     select version
!       open "Explore LAPACK code"
!  
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "num.inc"
!
!###############################################################################



module Num_LinAlg

  use GO, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LinAlg_Sym_Solve
  public  ::  LinAlg_Sym_Factorize
  public  ::  LinAlg_Sym_FactorSolve
  public  ::  LinAlg_Pod_Eigen


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'Num_LinAlg'
  
  
contains
  
  
  !
  ! solve A X = B   with A symmetric
  !

  subroutine LinAlg_Sym_Solve( A, B, X, status )
  
#ifdef with_mkl
    use MKL95_LAPack, only : SySV
#else
#ifdef with_mkl_17
    use LAPack95, only : SySV
#else
#ifdef with_lapack95
    use F95_LAPack, only : LA_SySV
#endif
#endif
#endif

    ! --- in/out ------------------------------
    
    real, intent(in)        ::  A(:,:)   ! (n,n)
    real, intent(in)        ::  B(:,:)   ! (n,m)
    real, intent(out)       ::  X(:,:)   ! (n,m)
    integer, intent(out)    ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'LinAlg_Sym_Solve'
    
    ! --- local -------------------------------

    integer                 ::  n, m
    real, allocatable       ::  LU(:,:)    ! dims (np,np)
    integer                 ::  i
#ifdef with_lapack
    integer, allocatable    ::  ipiv(:)
    real, allocatable       ::  work(:)
    integer                 ::  lwork
#endif

    ! --- begin -------------------------------
    
    ! output shape:
    n = size(X,1)
    m = size(X,2)
    
    ! check ..
    if ( any(shape(A) /= (/n,n/)) ) then
      write (gol,'("input shape A (",i0,",",i0,") does not match output shape (",i0,",",i0,")")') &
                       shape(A), shape(X); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( any(shape(B) /= (/n,m/)) ) then
      write (gol,'("input shape B (",i0,",",i0,") does not match output shape (",i0,",",i0,")")') &
                       shape(B), shape(X); call goErr
      TRACEBACK; status=1; return
    end if

    ! solve X from  A X = B
    if ( n == 1 ) then

      ! A is scalar; devide by only element :
      X = B / A(1,1)

    else

      ! Solve  X from   A X = B ;
      ! input A is replaced by LU on ouput, thus use copy ...
      ! create storage for copy:
      allocate( LU(n,n), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! copy:
      LU  = A
      ! fill output X on input with right-hand-side B :
      X = B

#ifdef with_lapack
      ! work arrays:
      allocate( ipiv(n), stat=status )
      IF_NOTOK_RETURN(status=1)
      lwork = n*3
      allocate( work(lwork), stat=status )
      IF_NOTOK_RETURN(status=1)
#endif

      ! solve:
#if defined (with_mkl) || defined (with_mkl_17)
      call SySV( LU, X, info=status )
#else
#ifdef with_lapack95
      call LA_SySV( LU, X, info=status )
#else
#ifdef with_lapack
      ! solve:
      if ( kind(LU) == 4 ) then
        call sSySV( 'U', n, m, LU, n, ipiv, X, n, work, lwork, status )
      else if ( kind(LU) == 8 ) then
        call dSySV( 'U', n, m, LU, n, ipiv, X, n, work, lwork, status )
      else
        write (gol,'("unsupported real kind ",i0)') kind(LU); call goErr
        TRACEBACK; status=1; return
      end if
#else
      write (gol,'("code was not compiled with mkl, lapack95, or lapack")'); call goErr
      TRACEBACK; status=1; return
#endif
#endif
#endif
      if (status/=0) then
        write (gol,'("from solving A X = B")'); call goErr
        !write (gol,'("  A : ")'); call goErr
        !do i = 1, size(A,1)
        !  write (*,*) A(i,:)
        !end do
        !write (gol,'("  B : ")'); call goErr
        !do i = 1, size(B,1)
        !  write (*,*) B(i,:)
        !end do
        write (gol,'("  info  : ",i6)') status; call goErr
        TRACEBACK; status=1; return
      end if

#ifdef with_lapack
      ! clear:
      deallocate( ipiv, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( work, stat=status )
      IF_NOTOK_RETURN(status=1)
#endif

      ! clear:
      deallocate( LU, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! scalar or matrix

    ! ok
    status = 0
    
  end subroutine LinAlg_Sym_Solve
  
  ! *
  
  !
  ! solve A X = B   with A symmetric using:
  !  - factorization into triangular factors:
  !      A = U^T U
  !  - solve X from:
  !      U^T U X = B
  !

  subroutine LinAlg_Sym_Factorize( A, U, status )
  
#ifdef with_mkl
    use MKL95_LAPack, only : PoTrF
#else
#ifdef with_mkl_17
    use LAPack95, only : PoTrF
#else
#ifdef with_lapack95
    use F95_LAPack, only : LA_PoTrF
#endif
#endif
#endif

    ! --- in/out ------------------------------
    
    real, intent(in)        ::  A(:,:)   ! (n,n)
    real, intent(out)       ::  U(:,:)   ! (n,n)
    integer, intent(out)    ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'LinAlg_Sym_Factorize'
    
    ! --- local -------------------------------

    integer             ::  n

    ! --- begin -------------------------------
    
    ! output shape:
    n = size(A,1)
    
    ! check ..
    if ( any(shape(A) /= (/n,n/)) ) then
      write (gol,'("input matrix A (",i0,",",i0,") should be square")') shape(A); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( any(shape(U) /= (/n,n/)) ) then
      write (gol,'("output shape U (",i0,",",i0,") does not match input size ",i0)') shape(U), n; call goErr
      TRACEBACK; status=1; return
    end if

    ! factorize:
    if ( n == 1 ) then

      ! A i a scalar, factorize:
      U = sqrt(A)

    else

      ! copy input into output:
      U = A
      ! positive-definite triangular factorization:
#if defined (with_mkl) || defined (with_mkl_17)
      call    PoTrF( U, uplo='U', info=status )
#else
#ifdef with_lapack95
      call LA_PoTrF( U, uplo='U', info=status )
#else
#ifdef with_lapack
      if ( kind(U) == 4 ) then
        call sPoTrF( 'U', n, U, n, status )
      else if ( kind(U) == 8 ) then
        call dPoTrF( 'U', n, U, n, status )
      else
        write (gol,'("unsupported real kind ",i0)') kind(U); call goErr
        TRACEBACK; status=1; return
      end if
#else
      write (gol,'("code was not compiled with mkl, lapack95, or lapack")'); call goErr
      TRACEBACK; status=1; return
#endif
#endif
#endif
      if (status/=0) then
        write (gol,'("from factorization A = U^T U")'); call goErr
        write (gol,'("  info  : ",i6)') status; call goErr
        TRACEBACK; status=1; return
      end if

    end if ! scalar or matrix

    ! ok
    status = 0
    
  end subroutine LinAlg_Sym_Factorize
  
  ! *

  subroutine LinAlg_Sym_FactorSolve( U, B, X, status )
  
#ifdef with_mkl
    use MKL95_LAPack, only : PoTrS
#else
#ifdef with_mkl_17
    use LAPack95, only : PoTrS
#else
#ifdef with_lapack95
    use F95_LAPack, only : LA_PoTrS
#endif
#endif
#endif

    ! --- in/out ------------------------------
    
    real, intent(in)        ::  U(:,:)   ! (n,n)
    real, intent(in)        ::  B(:,:)   ! (n,m)
    real, intent(out)       ::  X(:,:)   ! (n,m)
    integer, intent(out)    ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'LinAlg_Sym_FactorSolve'
    
    ! --- local -------------------------------

    integer             ::  n, m

    ! --- begin -------------------------------
    
    ! output shape:
    n = size(X,1)
    m = size(X,2)
    
    ! check ..
    if ( any(shape(U) /= (/n,n/)) ) then
      write (gol,'("input shape U (",i0,",",i0,") does not match output shape (",i0,",",i0,")")') &
                       shape(U), shape(X); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( any(shape(B) /= (/n,m/)) ) then
      write (gol,'("input shape B (",i0,",",i0,") does not match output shape (",i0,",",i0,")")') &
                       shape(B), shape(X); call goErr
      TRACEBACK; status=1; return
    end if

    ! solve X from  U^T U X = B
    if ( n == 1 ) then

      ! U is scalar; devide by only element :
      X = B / (U(1,1)**2)

    else

      ! Solve  X from   U^T U X = B :
      ! fill output X on input with right-hand-side B :
      X = B
      ! positive-definite triangular solve:
#if defined (with_mkl) || defined (with_mkl_17)
      call    PoTrS( U, X, uplo='U', info=status )
#else
#ifdef with_lapack95
      call LA_PoTrS( U, X, uplo='U', info=status )
#else
#ifdef with_lapack
      if ( kind(U) == 4 ) then
        call sPoTrS( 'U', n, m, U, n, X, n, status )
      else if ( kind(U) == 8 ) then
        call dPoTrS( 'U', n, m, U, n, X, n, status )
      else
        write (gol,'("unsupported real kind ",i0)') kind(U); call goErr
        TRACEBACK; status=1; return
      end if
#else
      write (gol,'("code was not compiled with mkl, lapack95, or lapack")'); call goErr
      TRACEBACK; status=1; return
#endif
#endif
#endif
      if (status/=0) then
        write (gol,'("from solving U^T U X = B")'); call goErr
        write (gol,'("  info  : ",i6)') status; call goErr
        TRACEBACK; status=1; return
      end if

    end if ! scalar or matrix

    ! ok
    status = 0
    
  end subroutine LinAlg_Sym_FactorSolve
  
  ! *
  
  !
  ! Eigenvalue decomposition of symetric positive definite matrix:
  !   A Q = Q Lambda
  ! Since Q is orthogonal:
  !  A = Q Lambda Q^T
  !

  subroutine LinAlg_Pod_Eigen( A, lambda, Q, status )
  
#ifdef with_mkl
    use MKL95_LAPack, only : SyEV
#else
#ifdef with_mkl_17
    use LAPack95    , only : PoTrF
#else
#ifdef with_lapack95
    use F95_LAPack  , only : LA_PoTrF
#endif
#endif
#endif

    ! --- in/out ------------------------------
    
    real, intent(in)        ::  A(:,:)      ! (n,n) symetric positive definite
    real, intent(out)       ::  lambda(:)    ! (n)
    real, intent(out)       ::  Q(:,:)      ! (n,n)
    integer, intent(out)    ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'LinAlg_Sym_Figen'
    
    ! --- local -------------------------------

    integer             ::  n
#ifdef with_lapack
    integer             ::  lwork
    real, allocatable   ::  work(:)  ! (lwork)
#endif
    integer             ::  i
    
      !! testing ...
      !real, allocatable    ::  B(:,:)  ! (n,n)
      !real, allocatable    ::  D(:,:)  ! (n,n)
      !integer              ::  j

    ! --- begin -------------------------------
    
    ! output shape:
    n = size(A,1)
    
    ! check ..
    if ( any(shape(A) /= (/n,n/)) ) then
      write (gol,'("input matrix A (",i0,",",i0,") should be square")') shape(A); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( size(lambda) /= n ) then
      write (gol,'("output size lambda (",i0,") does not match input size ",i0)') size(lambda), n; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( any(shape(Q) /= (/n,n/)) ) then
      write (gol,'("output shape Q (",i0,",",i0,") does not match input size ",i0)') shape(Q), n; call goErr
      TRACEBACK; status=1; return
    end if

    ! factorize:
    if ( n == 1 ) then

      ! A i a scalar, factorize:
      lambda = A(1,1)
      Q = 1.0      

    else

      ! copy input into output:
      Q = A
      ! symetric eigenvalue decomposition:
#if defined (with_mkl) || defined (with_mkl_17)
      call    SyEV( Q, lambda, jobz='V', uplo='U', info=status )
#else
#ifdef with_lapack95
      call LA_PoTrF( Q, lambda, jobz='V', uplo='U', info=status )
#else
#ifdef with_lapack
      ! length of work array:
      lwork = 5*n
      ! storage:
      allocate( work(lwork), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! switch:
      if ( kind(Q) == 4 ) then
        call sSyEV( 'V', 'U', n, Q, n, lambda, work, lwork, status )
      else if ( kind(Q) == 8 ) then
        call dSyEV( 'V', 'U', n, Q, n, lambda, work, lwork, status )
      else
        write (gol,'("unsupported real kind ",i0)') kind(Q); call goErr
        TRACEBACK; status=1; return
      end if
      ! clear:
      deallocate( work, stat=status )
      IF_NOTOK_RETURN(status=1)
#else
      write (gol,'("code was not compiled with mkl, lapack95, or lapack")'); call goErr
      TRACEBACK; status=1; return
#endif
#endif
#endif
      if (status/=0) then
        write (gol,'("from eigenvalue system A Q = Q Lambda")'); call goErr
        write (gol,'("  info  : ",i6)') status; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! check ...
      if ( any(lambda <= 0.0) ) then
        write (gol,'("found non-positive eigenvalues:")'); call goErr
        do i = 1, n
          write (gol,*) i, lambda(i); call goErr
        end do
      end if       

    end if ! scalar or matrix
    
    !!
    !! testing ...
    !!
    !! storage:
    !allocate( B(n,n), stat=status )
    !IF_NOTOK_RETURN(status=1)
    !allocate( D(n,n), stat=status )
    !IF_NOTOK_RETURN(status=1)
    !!
    !write (gol,*) ''; call goPr
    !write (gol,*) 'A = '; call goPr
    !do i = 1, n
    !  write (gol,*) A(i,:); call goPr
    !end do
    !!
    !write (gol,*) ''; call goPr
    !write (gol,*) 'lambda = ', lambda; call goPr
    !!
    !write (gol,*) ''; call goPr
    !write (gol,*) 'Q = '; call goPr
    !do i = 1, n
    !  write (gol,*) Q(i,:); call goPr
    !end do
    !!
    !D = 0.0
    !do i = 1, n
    !  D(i,i) = lambda(i)
    !end do
    !!
    !write (gol,*) ''; call goPr
    !write (gol,*) 'abs( Q Lambda Q^T - A ) ~ O :'; call goPr
    !B = matmul( matmul( Q, D ), transpose(Q) )
    !do i = 1, n
    !  write (gol,*) abs( B(i,:) - A(i,:) ); call goPr
    !end do
    !!
    !write (gol,*) ''; call goPr
    !write (gol,*) 'Q^T Q ~ I :'; call goPr
    !B = matmul( transpose(Q), Q )
    !do i = 1, n
    !  write (gol,*) B(i,:); call goPr
    !end do
    !write (gol,*) ''; call goPr
    !!
    !! clear:
    !deallocate( B, stat=status )
    !IF_NOTOK_RETURN(status=1)
    !deallocate( D, stat=status )
    !IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LinAlg_Pod_Eigen


end module Num_LinAlg

