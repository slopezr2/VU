!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!######################################################################

module C3PO_Grid_Spg

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Grid_Spg
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Grid_Spg'
  

  ! --- types ----------------------------------------
  
  type T_Grid_Spg
    integer             ::  T     ! triangular truncation
    integer             ::  np    ! number of data points
  contains
    procedure   ::  Init_T    => Grid_Spg_Init_T
    procedure   ::  Init_np   => Grid_Spg_Init_np
    procedure   ::  Final     => Grid_Spg_Final
    procedure   ::  Pnm       => Grid_Spg_Pnm
    procedure   ::  FourrierCoeff => Grid_Spg_FourrierCoeff_pnm
    procedure   ::  Nabla     => Grid_Spg_Nabla
    procedure   ::  Eval_Lons_Ulat => Grid_Spg_Eval_Lons_Ulat
    procedure   ::  Eval_Lons_Xm_1 => Grid_Spg_Eval_Lons_Xm_1
!    !final       ::  Grid_Spg_Final
  end type T_Grid_Spg

! adhoc ...
#ifdef without_f2003
#define XTYPE type
#else
#define XTYPE class
#endif
    
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** sh file
  ! ***
  ! ********************************************************************


  ! Spectral number :  SN = (T+1)*(T+2)/2
  ! (half triangle  m=0,..,T , n=m,..,T  ;
  !  we assume real restult thus only m>=0 required)
  
  integer function SpN( T )
  
    ! --- in/out -----------------
    
    integer, intent(in)   ::  T   ! triangular truncation
    
    ! --- begin ------------------
    
    SpN = (T+1)*(T+2)/2
    
  end function SpN
  
  ! reverse:
  !   (T+1)(T+2)/2 = n
  !   T**2 + 3 T + (2-2n) = 0
  !   T = ( -3 + sqrt(8*n+1) )/2
  
  integer function SpT( n )
  
    ! --- in/out -----------------
    
    integer, intent(in)   ::  n   ! numper of coeff
    
    ! --- begin ------------------
    
    SpT = nint( -3.0 + sqrt(8*n+1.0) )/2
    
  end function SpT
  
  ! *

  subroutine Grid_Spg_Init_T( self, T, status )
  
    ! --- in/out ------------------------------------
    
    XTYPE(T_Grid_Spg), intent(out)      ::  self
    integer, intent(in)                 ::  T
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Spg_Init_T'
    
    ! --- begin ----------------------------------

    ! set Triangular Truncation
    if ( T < 0 ) then
      write (gol,'("tried to set triangular truncation ",i6)') T; call goErr
      TRACEBACK; status=1; return
    end if

    ! store triangular trunction:
    self%T = T
    
    ! number of complex data points:
    self%np = SpN( self%T )
    
    ! ok
    status = 0

  end subroutine Grid_Spg_Init_T
  
  ! *

  subroutine Grid_Spg_Init_np( self, np, status )
  
    ! --- in/out ------------------------------------
    
    XTYPE(T_Grid_Spg), intent(out)      ::  self
    integer, intent(in)                 ::  np
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Spg_Init_np'
    
    ! --- begin ----------------------------------

    ! set Triangular Truncation
    if ( np < 0 ) then
      write (gol,'("tried to set triangular coefficients ",i6)') np; call goErr
      TRACEBACK; status=1; return
    end if

    ! compute triangular trunction:
    self%T = SpT( np )
    
    ! store number of complex data points:
    self%np = np
    
    ! ok
    status = 0

  end subroutine Grid_Spg_Init_np


  ! *
  

  subroutine Grid_Spg_Final( self )
  
    ! --- in/out ------------------------------------
    
    XTYPE(T_Grid_Spg), intent(inout)       ::  self

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Spg_Final'
    
    ! --- begin ----------------------------------

    ! reset:
    self%T = -999
    self%np = -999

  end subroutine Grid_Spg_Final
  
  ! *
  
  
  ! compute associated legendre coeff for
  !
  !      d X           2  d X                         dX   dX
  !  ( -------- , (1-mu ) ---- )   =   r cos(theta) ( -- , -- )
  !    d lambda           d mu                        dx   dy
  !
  ! given associated legendre coeff of X(lambda,mu)
  !
  
  subroutine Grid_Spg_Nabla( self, X, NablaX, status )

    ! --- in/out -----------------------------
    
    XTYPE(T_Grid_Spg), intent(in)          ::  self
    real, intent(in)                       ::  X(:,:)        ! (2,np)
    real, intent(out)                      ::  NablaX(:,:,:) ! (2,np,ndir)
    integer, intent(out)                   ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Spg_Nabla'
    
    ! complex:
    integer, parameter  :: ireal = 1
    integer, parameter  :: iimag = 2
    
    ! directions:
    integer, parameter  :: ilambda = 1
    integer, parameter  :: imu     = 2
    
    ! --- local ------------------------
    
    integer   ::  m, n, T
    integer   ::  k
    
    ! --- begin ----------------------------

    ! extract triangular truncation
    T = self%T
    
    ! loop over all coeff    
    k = 0
    do m = 0, T
      do n = m, T
        k = k + 1

        ! d(ln ps)                                i m lon
        ! -------- = sum  {i m X(m,n)} P(mu;m,n) e
        !  d lon     m,n

        !NablaX(1)%c(k) = (0.0,1.0) * m * Xnm%c(k)
        NablaX(ireal,k,ilambda) = - m * X(iimag,k)
        NablaX(iimag,k,ilambda) =   m * X(ireal,k)

        !      2  d(ln ps)       ~                 i m lon
        ! (1-mu ) -------- = sum X(m,n) P(mu;m,n) e
        !           d mu     m,n
        !
        !  ~      = (n+2) eps(m,n+1) X(m,n+1)                           , m<T, m=n
        !  X(m,n) = (n+2) eps(m,n+1) X(m,n+1) - (n-1) eps(m,n) X(m,n-1) , m<T, m<n<T
        !         =                           - (n-1) eps(m,n) X(m,n-1) , m<T,   n=T
        !         =                          0.0                        , m=T, n=T
        
        if ( m < T ) then
          if ( n == m ) then
            !NablaX(2)%c(k) = (n+2) * epsi(m,n+1) * Xnm%c(k+1)
            NablaX(:,k,imu) = (n+2) * epsi(m,n+1) * X(:,k+1)
          else if ( n > m .and. n < T ) then
            !NablaX(2)%c(k) = (n+2) * epsi(m,n+1) * Xnm%c(k+1) - &
            !                 (n-1) * epsi(m,n  ) * Xnm%c(k-1)
            NablaX(:,k,imu) = (n+2) * epsi(m,n+1) * X(:,k+1) - &
                              (n-1) * epsi(m,n  ) * X(:,k-1)
          else ! n == T
            !NablaX(2)%c(k) = (n-1) * epsi(m,n)   * Xnm%c(k-1)
            NablaX(:,k,imu) = (n-1) * epsi(m,n)   * X(:,k-1)
          end if
        else
          !NablaX(2)%c(k) = 0.0
          NablaX(:,k,imu) = 0.0
        end if

      end do  ! n
    end do  ! m
    
    ! ok
    status = 0
    
  contains

    !          n*n - m*m  1/2
    ! epsi = ( --------- )
    !          4*n*n - 1

    real function epsi( m, n )

      ! --- in/out --------------------------------

      integer, intent(in) ::  m, n

      ! --- begin --------------------------------

      epsi = sqrt((n*n-m*m)*1.0/(4*n*n-1.0))

    end function epsi
    
  end subroutine Grid_Spg_Nabla
  

  ! ====
  
  !
  ! Evaluate associate Legendre functions at given latitude.
  ! Corner (T,T) is set to 0.0   .
  !
  ! Based on recurent formula:
  !
  !
  !   mu P(mu;m,n-1) = eps(m,n) P(mu;m,n) + eps(m,n-1) P(mu;m,n-2)
  !
  !                                          n^2 - m^2
  !      mu = sin(lat)  ,   eps(m,n) = sqrt( --------- )
  !                                          4n^2 - 1
  !
  !   P(m,n) = ( mu P(m,n-1) - eps(m,n-1) P(m,n-2) ) / eps(m,n)
  !
  !
  !   P(0,0) = 1
  !
  !
  !   P(0,1) = mu / eps(0,1) = mu sqrt(3)
  !
  !                                   
  !              2m+1   1/2  sqrt(1-mu^2)^m 
  !   P(m,m) = ( ----- )     -------------- (2m)!
  !              (2m)!            2^m m!
  !
  !                          sqrt(1-mu^2)^m 
  !          = sqrt((2m+1)!) --------------        
  !                             2^m m!
  !
  !                                    sqrt(1-mu^2)sqrt(1-mu^2)^(m-1)
  !          = sqrt((2m+1)2m(2(m-1))!) ------------------------------ 
  !                                         2 2^(m-1) m (m-1)!
  !
  !                                      cos(lat)
  !          = P(m-1,m-1) sqrt((2m+1)2m) --------
  !                                        2m
  !
  !          = P(m-1,m-1) sqrt((2m+1)) cos(lat) / sqrt(2m)
  !
  !   P(m,m+1) = mu P(m,m) / eps(m,m+1)
  !
  !

  subroutine Grid_Spg_Pnm( self, Pnm, lat, status )

    ! --- in/out ----------------------------------

    XTYPE(T_Grid_Spg), intent(in)   ::  self
    real, intent(out)               ::  Pnm(:)   ! (np)
    real, intent(in)                ::  lat      ! rad
    integer, intent(out)            ::  status

    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Spg_Pnm'
    
    ! --- local ------------------------------------

    integer   ::  m, n
    integer   ::  k
    real      ::  mu, rmu
    real      ::  fmm, fmmp
    real      ::  eps, eps1
    
    ! --- begin ---------------------------------------

    ! Use EMOS library:    
    !Pnm = 0.0
    !call jspleg1( Pnm, lat, T-1 )
    !return

    if ( size(Pnm) /= self%np ) then
      write (gol,'("wrong size of output array:")'); call goErr
      write (gol,'("  size(Pnm)   : ",i6)') size(Pnm); call goErr
      write (gol,'("  expected    : ",i6)') self%np; call goErr
      write (gol,'("  truncation  : ",i6)') self%T; call goErr
      TRACEBACK; status=1; return
    end if
    
    mu  = sin(lat)
    rmu = sqrt(1.0-mu*mu)    ! cos(lat)

    ! loop over columns:
    do m = 0, self%T-1

      if ( m == 0 ) then
      
        fmmp    = sqrt(3.0)
        k = 1; Pnm(k) = 1.0          ! (0,0)
        k = 2; Pnm(k) = fmmp*mu      ! (0,1)

      else

        fmm  = fmmp * rmu / sqrt( 2.0*m )
        fmmp = fmm * sqrt( 2*m + 3.0 )

        ! n = m and n = m+1

        k = k+1; Pnm(k) = fmm
        k = k+1; Pnm(k) = fmmp * mu

      endif

      ! leave if truncation is reached:
      if ( m == self%T-1 ) exit

      eps1  = 1.0 / sqrt( 2*m + 3.0 )

      do n = m+2, self%T

        eps = sqrt((n*n*1.0-m*m*1.0)/(4.0*n*n-1.0))
        k = k+1; Pnm(k) = ( mu*Pnm(k-1) - eps1*Pnm(k-2) ) / eps
        eps1 = eps

      end do

    end do  ! m
    
    ! set corner:
    k = k + 1; Pnm(k) = 0.0
    
    ! ok
    status = 0
    
  end subroutine Grid_Spg_Pnm
  
  
  ! ***
  
  
  pure subroutine Grid_Spg_FourrierCoeff_pnm( self, Xnm, Pnm, X, status )
  
    ! --- in/out -----------------------------------
    
    XTYPE(T_Grid_Spg), intent(in)       ::  self
    complex, intent(in)                 ::  Xnm(:)  ! (np)
    real, intent(in)                    ::  Pnm(:)  ! (np)
    complex, intent(out)                ::  X(0:self%T)
    integer, intent(out)                ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Spg_FourrierCoeff_pnm'
    
    ! --- local ---------------------------------
    
    integer        ::  m, n, kx, kp
    
    ! --- begin ----------------------------------

    !! check input ...
    !if ( size(Xnm) /= self%np ) then
    !  write (gol,'("mismatch between size of Xnm and specified truncation:")'); call goErr
    !  write (gol,'("  shi_in%T    : ",i6)') self%T; call goErr
    !  write (gol,'("  shi_in%np   : ",i6)') self%np; call goErr
    !  write (gol,'("  Xnm         : ",i6)') size(Xnm); call goErr
    !  TRACEBACK; status=1; return
    !end if
    !if ( size(Pnm) /= self%np ) then
    !  write (gol,'("mismatch between size of Pnm and specified truncation:")'); call goErr
    !  write (gol,'("  shi%T       : ",i6)') self%T; call goErr
    !  write (gol,'("  shi%np      : ",i6)') self%np; call goErr
    !  write (gol,'("  Pnm         : ",i6)') size(Pnm); call goErr
    !  TRACEBACK; status=1; return
    !end if
  
    ! loop over all coeff:
    kx = 0
    kp = 0
    do m = 0, self%T
      X(m) = (0.0,0.0)
      do n = m, self%T 
        kx = kx + 1
        kp = kp + 1
        X(m) = X(m) + Xnm(kx) * Pnm(kp)
      end do
    end do
    
    ! adhoc: ensure that X(0) is real:
    X(0) = real(X(0)) * (1.0,0.0)
    
    ! ok
    status = 0
    
  end subroutine Grid_Spg_FourrierCoeff_pnm
  
  !
  ! Evaluate x(lat) and x(-lat) given:
  !   Xnm      : complex Associated Legendre Func. coefficients
  !   Pnm      : Associated Legendre Functions evaluated at mu=sin(lat)
  !   K        : number of lon values on a complete latitude circle;
  !              distance between to lon values is thus 2pi/K;
  !   lon0     : western boundary of lon grid (rad)
  !   Kout     : actual number of lon points from lon0 to east .
  !
  ! Evaluation at oposite latitdues at the same time saves some
  ! computation time.
  !
  !  x(t0+k2pi/K)                             ,  k=0,1,..,K0-1
  !
  !         T         i m [t0 + k 2pi/K]
  !     =  sum  X(m) e                        ,  k=0,1,..,K0-1
  !        m=-T                  
  !
  !              L such that KL >= (2T+1)
  !              k=(j+1)/L  ,  j=kL-1   ,   J=KL
  !            
  !         T         i m t0  i m (j+1)2pi/KL
  !     =  sum  X(m) e       e                    ,  j=-1,L-1,2L-1,...,(K0-1)L-1
  !        m=-T                        
  !
  !         T         im[t0+2pi/J]   i m j 2pi/J
  !     =  sum  [X(m)e            ] e             ,  j=-1,L-1,2L-1,...,(K0-1)L-1
  !        m=-T                         
  !
  !   _            
  !   X(m) = X(m) exp(i phi)    ,   phi=m[t0+2pi/J]
  !        =  (a+ib)(cos(phi)+isin(phi))
  !        = [ a cos(phi) - b sin(phi) ]  +  i [ a sin(phi) + b cos(phi) ]
  !
  ! FFT99 uses that X(-m) = conjg(X(m)), and thus requires X(m) for m>=0 only;
  ! higher fourrier coefficients are zero:
  !    X(0)  X(1)  ...  X(T)  0  0  0  ...      
  !   |---    1+floor(J/2) elements     ---|
  !
  !         o   o   x   x   x   0   0 ... 0    0   x   x
  !  m =            0   1   2                     J-2 J-1
  !
  ! Number of K longitudes on complete latitude circle,
  ! while Kout evaluations are required.
  ! Note that Kout might be both smaller and larger than K.
  ! Routine FFT99 returns J+2 evaluations:
  !    x(j=-1) x(j=0)  ...  x(j=J-1=-1) x(j=J=0)
  ! thus if Kout exceeds K+2 then the result should be taken 
  ! cyclic from the output of FFT99.
  !
  
  ! Evaluate given spherical harmonic coefficients and latitude

  subroutine Grid_Spg_Eval_Lons_Ulat( self, rXnm, lat_rad, KK, lon_start, &
                                        Kout, llgrid, status )

    ! --- in/out -------------------------
    
    XTYPE(T_Grid_Spg), intent(in)   ::  self
    real, intent(in)                ::  rXnm(:,:)  ! (2,np)
    real, intent(in)                ::  lat_rad    ! [rad]
    integer, intent(in)             ::  KK
    real, intent(in)                ::  lon_start
    integer, intent(in)             ::  Kout
    real,intent(out)                ::  llgrid(Kout)
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Spg_Eval_Lons_Ulat'
    
    ! --- local --------------------------

    complex, allocatable   ::  Xnm(:)  ! (np)
    complex , allocatable  ::  X(:)    ! (0:T)
    real, allocatable      ::  Pnm(:)  ! (np)
    integer                ::  i
    
    ! --- begin --------------------------
    
    ! storage:
    allocate( X(0:self%T) )
    allocate( Xnm(1:self%np) )
    allocate( Pnm(1:self%np) )
    
    ! convert real coeffs to complex:
    do i = 1, self%np
      Xnm(i) = cmplx( rXnm(1,i), rXnm(2,i) )
    end do
    
    ! evaluate Legendre functions at given latitude:
    call self%Pnm( Pnm, lat_rad, status )
    IF_NOT_OK_RETURN(status=1)

    ! compute Fourrier coefficients
    ! (latitude is implicit defined in Pnm)
    call self%FourrierCoeff( Xnm, Pnm, X, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! evaluate at longitudes given Fourrier coefficients:
    call self%Eval_Lons_Xm_1( llgrid, X, KK, lon_start, Kout, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( X )
    deallocate( Xnm )
    deallocate( Pnm )
    
    ! ok
    status = 0
    
  end subroutine Grid_Spg_Eval_Lons_Ulat
  
  ! *
  
  ! Shift to first lon:
  !
  !  x(t0+k2pi/K)                             ,  k=0,1,..,KK-1
  !
  !         T         i m [t0 + k 2pi/K]
  !     =  sum  X(m) e                        ,  k=0,1,..,K0-1
  !        m=-T                  
  !
  !         T          i m t0   k 2pi/K
  !     =  sum  [X(m) e      ] e              ,  k=0,1,..,K0-1
  !        m=-T                  
  !
  ! PERFORMANCE
  !   This routine is about 7% slower than with emos fft ...
  !
  !  real(4), external ::  etime
  !  real ::  tarr(2), t0
  !  t0=etime(tarr)
  !  print *, 'time:',etime(tarr)-t0
    
  
  subroutine Grid_Spg_Eval_Lons_Xm_1( self, llgrid, Xm, KK, lon_start, Kout, status )

    use Binas, only : pi
    use C3PO_Grid_Singleton, only : fft
  
    ! --- in/out -------------------------
    
    XTYPE(T_Grid_Spg), intent(in)   ::  self
    complex, intent(in)             ::  Xm(0:self%T)
    integer, intent(in)             ::  KK
    real, intent(in)                ::  lon_start
    integer, intent(in)             ::  Kout
    real,intent(out)                ::  llgrid(Kout)
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Spg_Eval_Lons_Xm_1'
    
    ! --- local --------------------------

    integer                ::  j, k
    integer                ::  m
    
    complex                ::  X(0:self%T)

    integer                ::  LL, JJ
    complex, allocatable   ::  C(:), F(:)
    real                   ::  fac

    ! --- begin --------------------------
    
    ! shift over west most longitude:
    !   X(m) :=  X(m) exp( i m lon0 )
    do m = 0, self%T
      X(m) = Xm(m) * exp( (0.0,1.0)*m*lon_start )
    end do

    ! choose L such that JJ=KK*LL >= 2T+1
    LL = 1
    do
      if ( KK*LL >= 2*self%T+1 ) exit
      LL = LL * 2
    end do
    
    ! number of longitudes to be evaluated:
    JJ = KK*LL

    ! fft arrays:
    allocate( C(JJ) )
    allocate( F(JJ) )
    
    ! fill coeff array:
    C = (0.0,0.0)
    C(1) = X(0)
    do m = 1, self%T
      C(1+m)    = X(m)
      C(JJ+1-m) = conjg(X(m))
    end do

    ! Apply fast fourrier transform.
    !
    !            1    J-1       i m lambda(j)
    !  F(j) = ------- sum C(m) e
    !         sqrt(J) m=0
    !
    !  where
    !
    !    lambda(j) = j 2pi/J = 0, dlon, 2*dlon, ...
    !                                           ____      ____
    !    C = [ X(0), X(1), .., X(T), 0, ..., 0, X(T), .., X(1) ]
    !
    ! Coeff in C are complex conj, thus imag part of F is zero.returned array is real.
    !
    F = fft( C, inv=.true. )
    fac = sqrt(JJ*1.0)

    ! Extract result.
    ! Since 'only' KK lons are evaluated (covering complete circle),
    ! select values cyclic for k>KK.
    ! Select first each group of LL elements
    do k = 0, Kout-1
      j = mod(k,KK)*LL + 1    ! 1, LL+1, 2*LL+1, ...
      llgrid(k+1) = real(F(j))*fac
    end do
    
    ! done
    deallocate( C )
    deallocate( F )
    
    ! ok
    status = 0
    
  end subroutine Grid_Spg_Eval_Lons_Xm_1
  

end module C3PO_Grid_Spg
