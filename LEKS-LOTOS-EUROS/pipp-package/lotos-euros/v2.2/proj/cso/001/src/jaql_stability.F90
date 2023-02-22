#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if

!#######################################################################
!
! JAQL_Stability
!   Joint Air Quality Library
!     Stability tools
!
!
!#######################################################################


module JAQL_Stability

  use GO, only : gol, GoErr, GoPr
 
  implicit none
  

  ! --- in/out --------------------------------
  
  private
  
  public  ::  HalfLevelKzMSP
  public  ::  HalfLevelKzIFS
  public  ::  f_m_stability, f_h_stability
  public  ::  calc_ustar
  public  ::  Atmospheric_Resistance
  
  public  ::  exposure_new
  public  ::  CalcMonin_stabclass

  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'JAQL_Stability'


contains
  
  ! Half level Kz following original scheme with correction for large difference between cell height
  
  subroutine HalfLevelKzMSP( halt, monin_inv, blh, ustar, Kz, status )
    
    ! --- in/out ---------------------------------

    real, intent(in)        ::  halt(:) ! (1:nlev+1) half level altitude [m]
    real, intent(in)        ::  monin_inv   ! Inverse MoninObukhov length [1/m]
    real, intent(in)        ::  blh     ! boundary layer height [m]
    real, intent(in)        ::  ustar   ! friction velocity [m/s]
    real, intent(out)       ::  Kz(:)   ! (1:nlev+1) Half-Level Kz following MSP scheme
    integer, intent(out)    ::  status 
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/HalfLevelKzMSP'
    
    ! --- local ----------------------------------

    real    ::  z1, z2
    real    ::  hm
    real    ::  ufac
    real    ::  z1fac, z2fac
    real    ::  p
    real    ::  r12
    real    ::  alpha
    real    ::  qfac, q1fac, q2fac
    
    ! --- begin ----------------------------

    ! at the mixing layer
    kz(2) = 1.0

    ! above mixing layer
    kz(3:) = 0.1

    ! both 'h' and 'blh' are in m:
    z1    = 0.5 * halt(1)
    z2    = 0.5 * ( halt(2) + halt(1) )
    hm    = max( halt(2), blh )

    ufac  = ustar*0.4/0.74
    if (monin_inv > 0.0) then
      z1fac = (1.0-z1/hm)**0.5
      z2fac = (1.0-z2/hm)**0.5
      P     = (1.0-z2fac)*(1.0+z1fac)/( (1.0-z1fac)*(1.0+z2fac) )
      r12   = (log(P) + 2*(1.0+6.3*hm*monin_inv)*(1.0/z2fac-1.0/z1fac))/ufac
      if (r12 < 0)  print *, monin_inv, r12, log(P), 2*(1.0+6.3*hm*monin_inv)*(1.0/z2fac-1.0/z1fac), ufac
    else
      alpha=9.0*abs(monin_inv)
      qfac =(1.0+alpha*hm)**0.5
      q1fac=(1.0+alpha*z1)**0.5
      q2fac=(1.0+alpha*z2)**0.5
      r12  = (log((qfac+q2fac)*(qfac-q1fac)/((qfac+q1fac)*(qfac-q2fac)))/qfac - &
             log((1.0+q2fac)*(1.0-q1fac)/((1.0+q1fac)*(1.0-q2fac))) )/ufac
      if (r12<0) print *, r12, monin_inv, log((qfac+q2fac)*(qfac-q1fac)/((qfac+q1fac)*(qfac-q2fac)))/qfac, &
             log((1.0+q2fac)*(1.0-q1fac)/((1.0+q1fac)*(1.0-q2fac))), ufac
    endif
    ! kz changed by correction factor for dc/dz by MvL and MSP
    !  kz(i,j,1) = kz(i,j,1) * 0.5*h(i,j,2)/h(i,j,1)
    kz(1) = (z2-z1)/R12
    ! new expression for Kz to be checked, what is the reference here?
    ! TvN

    ! truncate Kz between surface and mixing layer to [1,100] :
    kz(1) = min( max( 1.0, kz(1) ), 100.0 )

    ! ok
    status = 0
  
  end subroutine HalfLevelKzMSP
  
  
  ! Half level Kz following IFS scheme
  ! Returns Kz profile at meteo half levels
  ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes
  
  subroutine HalfLevelKzIFS( temper, duudz, ri, monin, halt, blh, ustar, KzIFS, status ) 
    
    ! --- in/out ---------------------------------

    real, intent(in)        ::  temper  ! temperature [K]
    !real, intent(in)        ::  dse     ! dry static energy [J/kg]
    !real, intent(in)        ::  uu_ifs  ! wind profile [m/s]
    real, intent(in)        ::  duudz   ! Half level gradient [1/s]
    real, intent(in)        ::  ri      ! Richardson number [-]
    real, intent(in)        ::  monin   ! MoninObukhov length [m]
    real, intent(in)        ::  halt    ! Half level altitude [m]
    real, intent(in)        ::  blh     ! boundary layer height [m]
    real, intent(in)        ::  ustar   ! Friction velocity [m/s]
    real, intent(out)       ::  KzIFS   ! Half-Level Kz following IFS scheme [m2/s]
    integer, intent(out)    ::  status 
    
    ! --- const ---------------------------------- 
    
    character(len=*), parameter ::  rname = mname//'/HalfLevelKzIFS'

    ! --- local ----------------------------------
    
    character(len=32)       ::  regime
    real                    ::  lmix, lmix_lambda

    ! --- begin ----------------------------------
    
    ! which stability regime    
    regime = GetRegime( monin, ri, halt, blh )
    
    select case ( trim(regime) )
      
      case ( 'M-O' )
        
        ! asymptotic mixing length following eq. (3.54):
        lmix_lambda = 150 ! m

        ! mixing length:
        lmix = mixing_length(halt,lmix_lambda )

        ! following Monin-Obukhov similarities for heat transfer:
        KzIFS = MoninObukhov_Kz_heat_IFS(lmix,  duudz,  ri)
      
      case ( 'Louis' )
      
        ! asymptotic mixing length following eq. (3.54):
        if ( halt <= blh ) then
           lmix_lambda = 0.1 * blh ! [m]
        else 
           lmix_lambda = 30. ! [m]

        end if    

        ! mixing length:
        lmix = mixing_length( halt,lmix_lambda )

        ! following Louis scheme for heat transfer:
        KzIFS = Louis_Kz_heat_IFS(lmix,  duudz,  ri)

      case ( 'EDMF' )
        
        ! Following EDMF scheme ;
        ! use only surface driven diffusion, no in-cloud diffusion yet:
        KzIFS = EDMF_Kz_sfc_heat(ustar, monin,  blh,  halt )
      
      case default
        write( gol, '("Unknown stability regime found: ",a )' ) trim(regime) ; call GoErr
        TRACEBACK;status=1;return

    end select
    
    ! ok
    status = 0
    
  end subroutine HalfLevelKzIFS
  

  !
  ! * Calculation of stability regime
  !
  
  elemental function GetRegime( monin, ri, halt, blh ) result(regime)
    
    ! ---in/out---
    real, intent(in)   ::  monin
    real, intent(in)   ::  ri
    real, intent(in)   ::  halt
    real, intent(in)   ::  blh
    character(len=32)  ::  regime
        
    ! --- begin ---
    
    ! stable or unstable:
    if ( monin >= 0 ) then
      ! stable surface layer: check Richardson number:
      if (ri > 0 ) then
        ! use Louis Scheme:
        regime = 'Louis'
      else 
        ! use monin-obukhov similarity:
        regime = 'M-O'
      end if
    
    else    
      ! unstable surface layer; check height:
      if ( halt < blh ) then
        ! within mixing layer, use EDMF scheme
        regime = 'EDMF'
      else
        ! above mixing layer; check Richardson number:
        if ( ri > 0 ) then
          ! use Louis scheme:
          regime = 'Louis'
        else
          ! use monin-obukhov similarity:
          regime = 'M-O'
        end if
        
      end if ! within mixing layer?
    
    end if ! stable or unstable?
        
  end function GetRegime
  
  ! ***
  
  elemental function mixing_length(z,lamb) result( mix_L )
  
    ! Dimensionless shear function for momentum (3.19)
    ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes
    
    use binas, only : vkarman

    ! local 
    real,intent(in) :: z      ! height above surface [m]
    real,intent(in) :: lamb ! assymptotic length scale [m]
    real            :: mix_L  ! mixing length [m]

    ! mixing length scale :
    mix_L = 1.0 / ( 1.0/(vkarman * z) + 1.0/lamb )
        
  end function mixing_length  
  
  ! ***

  elemental function MoninObukhov_Kz_heat_IFS( lmix,duudz,ri ) result(Kz)
  
    ! Diffussion coefficient for heat using Monin-Obukov similarity (3.54)
    ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes

    ! --- in/out --------------------------- 

    real, intent(in)   :: lmix
    real, intent(in)   :: duudz
    real, intent(in)   :: ri
    real               :: Kz
    
    ! --- local ----------------------------------
    real ::  Phi_M  ! shear function for momentum [-] 
    real ::  Phi_H  ! shear function for heat [-]   

    ! --- begin ----------------------------   
    
    ! evaluate gradient functions using zeta=Ri:
    Phi_M = shear_function_Phi_M( ri )
    Phi_H = shear_function_Phi_H( ri )

    ! Diffusion coefficient
    kz = MoninObukhov_Kz_heat(Phi_M,Phi_H,lmix,duudz )

  end function MoninObukhov_Kz_heat_IFS
 
  ! ***
  
  elemental function  shear_function_Phi_M ( zeta ) result ( Phi_m )
  
    ! Dimensionless shear function for momentum (3.19)
    ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes
    
    ! local
    real,intent(in) :: zeta  ! either z/monin of Ri (-)
    real            :: Phi_M ! shear function

    ! shear function:
    Phi_m = ( 1.0 - 16 * zeta )**(-1.0/4.0)

  end function shear_function_Phi_M

  ! ***  
  
  elemental function shear_function_Phi_H( zeta ) result (Phi_H)
                     
    ! Dimensionless shear function for heat and moisture (3.19)
    ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes
    
    ! local
    real,intent(in) :: zeta ! either z/monin of Ri (-)
    real            :: Phi_H ! shear function
    
    ! shear function
    Phi_H = ( 1.0 - 16 * zeta )**(-1.0/2.0)

  end function shear_function_Phi_H
  
  ! ***
  
  elemental function shear_function_Phi_M0( zeta ) result(Phi_M0)
  
    ! Dimensionless shear function for momentum (3.19)
    ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes

    ! local
    real,intent(in) :: zeta ! either z/monin of Ri (-)
    real            :: Phi_M0 ! shear function

    ! shear function:
    Phi_m0 = ( 1.0 - 15 * zeta )**(-1.0/3.0)

  end function shear_function_Phi_M0

  ! ***  
  
  elemental function shear_function_Phi_H0 (zeta) result( Phi_H0 )
  
    ! Dimensionless shear function for heat and moisture (3.42)
    ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes

    ! local
    real,intent(in) :: zeta ! either z/monin of Ri (-)
    real            :: Phi_H0 ! shear function
        
    ! shear function
    Phi_h0 = ( 1.0 - 39 * zeta )**(-1.0/3.0)

  end function shear_function_Phi_H0
  
  ! ***
  
  elemental function MoninObukhov_Kz_heat(Phi_M,Phi_H,lmix,dU_dz) result( kz_heat )
  
    ! local 
    real,intent(in) :: Phi_M   ! momentum shear function [-]
    real,intent(in) :: Phi_H   ! momentum shear function [-]
    real,intent(in) :: lmix    ! mixing length [m]
    real,intent(in) :: dU_dz   ! half level windpseed gradient [m]
    real            :: kz_heat ! difusion coefficient
    
    kz_heat = lmix**2 / (Phi_M * Phi_H) * dU_dz   ! m2/s 

        
  end function MoninObukhov_Kz_heat  
  
  ! ***
  
  elemental function Louis_Kz_heat_IFS( lmix,duudz,ri_hlv ) result(Kz)
  
    ! Diffussion coefficient for heat using Louis scheme (3.55)
    ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes

    ! --- in/out --------------------------- 
    real,intent(in)   :: lmix
    real,intent(in)   :: duudz    
    real,intent(in)   :: ri_hlv    
    real              :: Kz   


    ! Diffusion coefficient
    Kz = Louis_Kz_heat(ri_hlv,lmix,duudz)


  end function Louis_Kz_heat_IFS
  
  ! *
  
  elemental function Louis_Kz_heat (Ri,lmix,dU_dz) result( kz_heat )
  
    ! local 
    real,intent(in) :: Ri       ! Richardson number [-]
    real,intent(in) :: dU_dz    ! vertical gradient of horizontal windspeed U [(m/s)/m]
    real,intent(in) :: lmix     ! mixing length [m]    
    real            :: f_LTG_H  ! function by Louis/Tiedtke/Geleyn following (3.56)
    real            :: kz_heat  ! difusion coefficient    
    
    real, parameter :: b = 5.
    real, parameter :: d = 1.    
    
    ! evaluate
    f_LTG_H = 1.0 / ( 1.0 + 2*b*Ri * (1.0 + d*Ri)**(1.0/2.0) )

    ! K_H from eq. (3.55):
    !           m2          1     (m/s)/m
    kz_heat = lmix**2 * f_LTG_H * dU_dz   ! (1/s)
        
  end function Louis_Kz_heat  
    
  ! ***
  
  elemental function EDMF_Kz_sfc_heat(ustar,monin,blh,z_hlv) result(Kz)
  
    ! Diffussion coefficient for surface driven diffustion of heat (3.41)
    ! described in IFS DOCUMENTATION Cy41r1 Part IV, physical processes

    ! --- in/out ---------------------------
    real, intent(in)       :: ustar
    real, intent(in)       :: monin
    real, intent(in)       :: blh
    real, intent(in)       :: z_hlv          
    real                   :: Kz
    
    !-----local--------------------------------
    
    real                   :: Phi_M0

    
    ! Shear function
    Phi_M0 =  shear_function_Phi_M0(z_hlv / monin )
    
    ! Diffusion coefficient
    Kz = EDMF_Kz_heat(ustar,Phi_M0,z_hlv,blh)


  end function EDMF_Kz_sfc_heat
  
  ! *
  
  elemental function EDMF_Kz_heat (u_star,Phi_M0,z,blh) result( kz_heat )
  
    use binas, only : vkarman
  
    ! local 
    real,intent(in) :: u_star   ! friction velocity [m/s]
    real,intent(in) :: Phi_M0   ! mixing length [m]
    real,intent(in) :: z        ! height above surface [m]
    real,intent(in) :: blh      ! boundary layer height [m] 
    real            :: kz_heat  ! difusion coefficient  [m2/s] 
    
    ! K_sfc_H from eq. (3.41), with extra factor z :
    !          1        m/s   /      1             m    1
    kz_heat = vkarman * u_star / Phi_M0 * z * ( 1.0 - z/blh )**2  ! [m2/s]
 
        
  end function EDMF_Kz_heat  
  
  !
  ! stability function for momentum (from Jacobsen-2005)
  !
  
  elemental function f_m_stability( z, znul, monin_inv ) result (f_m)
  
    use Binas, only : beta_m, gamma_m
    
    ! --- in/out ---------------------------------

    real, intent(in) :: z
    real, intent(in) :: znul
    real, intent(in) :: monin_inv
    real             :: f_m  ! output variable
     
    ! --- const ----------------------------------
    
    ! arbitrary small number to trap numerical errors
    real, parameter  ::  eps = 1.0e-6
    ! threshold for interval [delta,0.0]
    real, parameter  ::  delta = 1.0 - (1.0+eps)**4

    ! --- local ----------------------------------
    
    real  ::  phi_mz
    real  ::  phi_mznul
    
    ! --- begin ----------------------------------
    
    if ( monin_inv >= 0.0 ) then     ! stable case
      
      ! take stable case from formula (8.31, Jacobsen-2005)
      f_m = log(z/znul) + beta_m * (z-znul) * monin_inv

    else if ( monin_inv < 0.0 ) then   ! unstable case
      
      ! trap numerical errors for monin_inv too close to zero;
      ! use that z >= z0 :
      !
      !          phi_mz                    - 1 <        eps
      !     (1-gamma_m*z*monin_inv)**0.25  - 1 <        eps
      !      1-gamma_m*z*monin_inv             <     (1+eps)**4
      !        gamma_m*z*monin_inv             > 1 - (1+eps)**4
      !
      if ( gamma_m*znul*monin_inv > delta ) then
        
        ! take neutral case from formula (8.31, Jacobsen-2005)
        f_m = log(z/znul)
        
      else
        
        ! take unstable case from formula (8.31, Jacobsen-2005)
        phi_mz    = ( 1.0 - gamma_m*z   *monin_inv )**0.25
        phi_mznul = ( 1.0 - gamma_m*znul*monin_inv )**0.25
        f_m = ( log( ( phi_mz   -1 ) / (phi_mz    + 1 ) ) + 2*atan( phi_mz   ) ) - &
              ( log( ( phi_mznul-1 ) / (phi_mznul + 1 ) ) + 2*atan( phi_mznul) )

      end if
    
    end if
  
  end function f_m_stability  
  
  
  ! ***
  
  
  ! Stability function for heat 
  
  elemental function f_h_stability( z, znul, monin_inv ) result (f_h)
  
    use Binas, only : beta_h, gamma_h, pr_t
    
    ! --- in/out ---------------------------------
    
    real, intent(in) :: znul
    real, intent(in) :: z
    real, intent(in) :: monin_inv
    real             :: f_h
     
    ! --- const ----------------------------------
    
    ! arbitrary small number to trap numerical errors
    real, parameter  ::  eps = 1.0e-6
    ! threshold:
    real, parameter  ::  delta = 1.0 - (1.0+eps)**2
  
    ! --- local ----------------------------------
    
    real  ::  phi_hz
    real  ::  phi_hznul
    
    ! --- begin ----------------------------------
    
    if ( monin_inv >= 0.0 ) then     ! stable case
      
      ! take stable case from formula (8.38, Jacobsen-2005)
      f_h = pr_t * log(z/znul) + beta_h * (z-znul) * monin_inv
      
    else if ( monin_inv < 0.0 ) then   ! unstable case
      
      ! trap numerical errors for monin_inv too close to zero;
      ! use that z >= z0 :
      !
      !        phi_hz0                  -  1  <      eps
      !  (1-gamma_h*z0*monin_inv)**0.5  -  1  <      eps
      !   1-gamma_h*z0*monin_inv              <   (1+eps)**2
      !     gamma_h*z0*monin_inv              > 1-(1+eps)**2
      !
      if ( gamma_h*znul*monin_inv > delta ) then
        
        ! take neutral case from formula (8.38, Jacobsen-2005)
        f_h = pr_t * log(z/znul)
        
      else

        ! take unstable case from formula (8.38, Jacobsen-2005)
        phi_hz    = ( 1.0 - gamma_h*z   *monin_inv )**0.5
        phi_hznul = ( 1.0 - gamma_h*znul*monin_inv )**0.5
        f_h = pr_t * ( log( ( phi_hz   -1 ) / (phi_hz    + 1 ) ) - &
                       log( ( phi_hznul-1 ) / (phi_hznul + 1 ) )   )

      end if

    end if
  
  end function f_h_stability  
  
  
  ! function to calculate u* from stability for mass
  elemental function calc_ustar( kappa_stab, windmg, fm ) result ( ust )
           
    ! --- arguments -----------

    real, intent(in)    ::  kappa_stab
    real, intent(in)    ::  windmg
    real, intent(in)    ::  fm
    real                ::  ust
    
    !------------------------------
    
    ust = kappa_stab * windmg / fm
    
  end function calc_ustar  
  

  ! ***

  
  elemental function Atmospheric_Resistance( kappa_stab, ustar, fh ) result ( Ra )
  
    ! --- arguments -----------

    real, intent(in)  ::  kappa_stab
    real, intent(in)  ::  ustar
    real, intent(in)  ::  fh
    real              ::  Ra

    ! --- begin ---------------------------------
           
    Ra = fh / kappa_stab / ustar
    
  end function Atmospheric_Resistance
  
  
  ! ***
  
  subroutine CalcMonin_stabclass( stab_class_letter, z0, monin, status )
  
    ! in/out
    character(len=1), intent(in)  ::  stab_class_letter
    real, intent(in)              ::  z0
    real, intent(out)             ::  monin
    integer, intent(out)          ::  status
    
    ! --- const ----
    character(len=*), parameter ::  rname = mname//'/CalcMonin_stabclass'
    
    ! --- local ---
    real    ::  a, b
    real    ::  elinv
    
    ! inverse of monin-obuhkov length
    select case ( stab_class_letter )

      ! extremely unstable case
      case ( 'A', 'a' )
        a = -0.096
        b =  0.029
      ! moderately unstable case
      case ( 'B', 'b' )
        a = -0.037
        b =  0.029
      ! slightly unstable case
      case ( 'C', 'c' )
        a = -0.002
        b =  0.018
      ! neutral case
      case ( 'D', 'd' )
        !a =  0.0
        a = 1.e-6
        b =  0.0
      ! slightly stable case
      case ( 'E', 'e' )
        a =  0.004
        b = -0.018
      ! moderately stable case
      case ( 'F', 'f' )
        a =  0.035
        b = -0.036

      case default
        write( gol, '("Unknown exposure class: ", a )' ) trim(stab_class_letter)
    end select

    ! Golder calculation (1972)
    elinv = a + b * log(z0) 
    monin = 1.0/elinv

    ! ok
    status = 0

  end subroutine calcMonin_stabclass
    
  
  ! ***
      

  subroutine exposure_new( rad, wsurf, tcc, coszen, sd, h_local, stab_class_letter, stab_class, status )
    
    
    ! --- in/out ---
    real, intent(in)              ::  rad    ! radiation (W/m2)
    real, intent(in)              ::  wsurf  ! wind speed (m/s)
    real, intent(in)              ::  tcc    ! total cloud cover (0-1)
    real, intent(in)              ::  coszen ! cos of solar zenith angle (-1-1)
    real, intent(in)              ::  sd     ! snowdepth
    integer, intent(in)           ::  h_local ! local hour (need to check for low radiation (morning or evening) )
    character(len=1), intent(out) ::  stab_class_letter  ! letter of stability class
    integer, intent(out)          ::  stab_class ! number of stability class A=1, B=2, C=3, D=4, E=5, F=6 --> use for plotting
    integer, intent(out)          ::  status
    
    ! --- const ----
    
    character(len=*), parameter ::  rname = mname//'/exposure_new'
    
    character(len=1), parameter ::  day_classes_letter(5,5) = &
                                          reshape( (/'A','A','B','E','C', &
                                                     'A','B','C','D','D', &
                                                     'B','B','C','D','D', &
                                                     'C','C','D','D','D', &
                                                     'C','D','D','D','D'/), &
                                                     (/5,5/) )
    integer, parameter          ::  day_classes(5,5) = &
                                          reshape( (/1,1,2,5,4, &
                                                     1,2,3,4,4, &
                                                     2,2,3,4,4, &
                                                     3,3,4,4,4, &
                                                     3,4,4,4,4/), &
                                                     (/5,5/) )
    character(len=1), parameter ::  night_classes_letter(2,5) = &
                                          reshape( (/'F','F', &
                                                     'E','F', &
                                                     'D','E', &
                                                     'D','D', &
                                                     'D','D'/), &
                                                     (/2,5/) )
    integer, parameter          ::  night_classes(2,5) = &
                                          reshape( (/6,6, &
                                                     5,6, &
                                                     4,5, &
                                                     4,4, &
                                                     4,4/), &
                                                     (/2,5/) )
    ! --- local ---
    integer       ::  wsurf_index
    integer       ::  cloud_index
    integer       ::  rad_index
    
    ! --- begin
    
    ! wind index
    if ( wsurf < 2.0 ) then
      wsurf_index = 1
    else if ( wsurf >= 2.0 .and. wsurf < 3.0 ) then
      wsurf_index = 2
    else if ( wsurf >= 3.0 .and. wsurf < 5.0 ) then
      wsurf_index = 3
    else if ( wsurf >= 5.0 .and. wsurf < 6.0 ) then
      wsurf_index = 4
    else if ( wsurf >= 6.0 ) then
      wsurf_index = 5
    end if
    
    ! cloud cover index
    if ( tcc > 0.5 ) then
      cloud_index = 1
    else
      cloud_index = 2
    end if
    
    ! radiation index
    if ( rad >= 700 ) then
      rad_index = 1
    else if ( rad >= 350 .and. rad < 700 ) then
      rad_index = 2
    else if ( rad < 350 .and. rad >= 125 ) then
      rad_index = 3
    else
      if ( h_local < 12 ) then
        ! Low radiation in the morning (just after sunrise), so neutral or still stable with very low wind
        rad_index = 4
      else
        ! Low radiation in the evening (just before sunset), so neutral or still unstable with low wind
        rad_index = 5
      end if
    end if
    
    ! daytime
    if ( coszen > 0 ) then
      stab_class_letter = day_classes_letter(rad_index,wsurf_index)
      stab_class        = day_classes(rad_index,wsurf_index)
    ! nigthtime
    else
      stab_class_letter = night_classes_letter(cloud_index,wsurf_index)
      stab_class        = night_classes(cloud_index,wsurf_index)
    end if    
    
    
    ! ok
    status = 0
    
  end subroutine exposure_new
  
  
        

end module JAQL_Stability

!#############################################################################
!
! Python module used to test IFS diffusion scheme.
! Comment in this module was lost during translation to Fortran ..
! Therefore included here to save it for future.
!
! Some of the comment here should be ported to "le_data_varibles.F90",
! or better, included in 1D routines in this module.
!
!  From:                    To:
!
!  MoninObukovLength        Variables_Calc_MoninObukhovLengthIFS
!
!  ScaleParameters          Variables_Calc_ScaleParameterSpecificHumidity
!                           Variables_Calc_ScaleParameterDryStaticEnergy
!
!  BottomWindspeed2         Variables_Calc_BottomWindSpeed
!
!  VirtualTemperatureFlux   Variables_Calc_VirtualTemperatureFlux
!
!  Richardson_Number        Variables_Calc_RichardsonNumber
!
!
!######################################################################
!
!
!  #! /usr/bin/env python
!
!  """
!  Test IFS diffussion scheme.
!
!
!  Tracer diffusion
!  ----------------
!
!  From section 3.7:
!    "Tracers are difused in the same way as heat and moisture, but no mass 
!    flux term is used.
!    - The surface boundary condition consists of an externally specified flux. 
!    - The implicitness factor is set to 1, because a higher value is not necessary 
!      for stability. As for momentum, heat and moisture the implicit solver uses 
!      the dynamics term as source terms to obtain balance and small time step 
!      dependence for long time steps. It can be demonstrated that implicitness 
!      factors larger than 1 can lead to negative tracer concentrations
!      due to the combination with the dynamics source term."
!
!
!  Monin-Obukov length
!  -------------------    
!
!    Iterative search for L?
!    The Obukov length L depends on surface fluxes,
!    but the surface fluxes also depend on L .
!
!    From eq. (3.8) :
!       L = func( u*, Q0v(u*,s*,q*) )
!
!    If u* is known, then s* and q* are related to sshf and slhf:
!      s* = func( u*, sshf )
!      q* = func( u*, slhf )
!    thus effectively:
!      L = func( u*, sshf, slhf )
!
!    These 3 variables are currently available from the ECMWF 
!    operational meteo output:
!
!      Name                        Short Name      Units   Parameter ID
!      --------------------------  --------------  ------  -----------------
!      Friction velocity           zust            m s-1   228003  (  3.228)
!      Surface sensible heat flux  sshf            J m-2      146  (146.128)
!      Surface latent heat flux    slhf            J m-2      147  (147.128)
!
!
!  Atmospheric resistance
!  ----------------------
!
!  For deposition, derived quantities are needed.
!    - Ra is the 'atmospheric resistance' between z and ztop :
!                         ztop
!                           (   1
!          Ra([z,ztop])  =  | ----- ds        [s/m]
!                           ) Kz(s)
!                          s=z
!      Requires integrals over gradient functions Phi (see (3.19) and (3.21))
!      into stable profile functions Psi (see (3.20) and (3.22)).
!      See also header of "le_stability.F90".
!
!
!  References
!  ----------
!
!    IFS DOCUMENTATION / Cy41r1 / Operational implementation 12 May 2015 / PART IV: PHYSICAL PROCESSES
!    http://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf
!
!  """
!
!  #
!  # Units for better readability:
!  #   energy =        force          * distance 
!  #          = (mass * acceleration) * distance
!  #     J    =   kg        m/s2           m    
!  #     J    =   kg m2/s2
!  #   J/kg   =  m2/s2
!
!  class Diffusion( object ) :
!
!      """
!      Definition of IFS diffusion method.
!      """
!
!      def __init__( self ) :
!
!          """
!          Initialize diffusion process with surface and meteo properties.
!          """
!
!          # modules:
!          import binas
!
!          # Von Karman constant (0.4)
!          self.kappa = binas.vkarman  # 1 ?
!
!          # scale height of the boundary layer depth as used in
!          # and explained under eq. (3.18):
!          self.z_i   = 1000.0  # m
!
!          # acceleration of gravity:
!          self.g = binas.grav  # m/s2
!
!          # specific heat at constant presure of moist air:
!          self.c_p = binas.cp_air   # J/kg/K ?
!
!          # scale factor for virtual temperatures etc:
!          #  Rdry  : gas constant for dry air (287 J/K/kg)
!          #  Rvap  : gas constant for water vapour (461 J/K/kg)
!          self.epsilon = binas.Rvap / binas.Rdry - 1.0    # 1
!
!          # limitation of zeta=z/L as described at end of 
!          # introduction paragraphs of section 3.2:
!          self.zeta_max = 5.0
!
!          # no values yet:
!          regime = None
!
!      #enddef __init__
!
!      # *
!
!      def HalfLevel_Kz( self, z_hlv, z, u, v, T, q, 
!                         blh, u_star, sshf, slhf, rho ) :
!
!          """
!          Return Kz profile at half levels.
!
!          Number of layers is denoted by 'n' ;
!          the n full levels have logical numbering 1:n (stored 0:n-1),
!          the n+1 half levels have logical numbering 0:n or 0.5:n+0.5 (stored 0:n).
!
!          Arguments:
!             z_hlv[n+1] : height of half levels [m]
!             z[n]       : height            of full levels [m]
!             u[n]       : eastward  wind    at full levels [m/s]
!             v[n]       : northward wind    at full levels [m/s]
!             T[n]       : temperature       at full levels [K]
!             q[n]       : specific humidity at full levels [kg/kg]
!             blh        : boundary layer height [m]
!             u_star     : reference u* [m/s], valid for grass surface?
!             sshf       : surface sensible heat flux [J/m2/s]
!             slhf       : surface latent   heat flux [J/m2/s]
!             rho        : air density [kg/m3]
!
!          Return values:
!             Kz_hlv[n+1] : diffusion coefficients at half levels
!
!          """
!
!          # modules:
!          import numpy
!
!          # dims:
!          n = len(z)
!
!          # compute profile of dry static energy:
!          s = self.DryStaticEnergy( z, T )
!
!          # absolute windspeed:
!          UU = numpy.sqrt( u**2 + v**2 )
!          # replace bottom value:
!          UU[0] = numpy.sqrt( self.BottomWindspeed2( u[0], v[0], T[0], Q0v ) )
!
!          # derivative to z at half levels:
!          dU_dz_hlv = self.HalfLevelGradient( z_hlv, z, UU )
!
!          # storage for Kz at half levels:
!          Kz_hlv = numpy.zeros((n+1),float)
!
!          # obtain scale parameters for vertical profiles given fluxes:
!          s_star,q_star = self.ScaleParameters( u_star, sshf, slhf, T[0], rho )
!
!          # surface flux:
!          Q0v = self.VirtualTemperatureFlux( u_star, s_star, q_star, T[0] )
!
!          # Monin-Obukov length:
!          L = self.MoninObukovLength( u_star, T[0], Q0v )
!
!          # asymptotic mixing length following description below eq. (3.54):
!          lmix_lambda = 150.0  # m
!
!          # loop over half levels:
!          for ihlv in range(n+1) :
!
!              # bounds:
!              if (ihlv == 0) or (ihlv == n) :
!
!                  # no flux through sfc and top:
!                  Kz_hlv[ihlv] = 0.0
!
!              else :
!
!                  # mixing length:
!                  lmix = self.mixing_length( lmix_lambda, z_hlv[ihlv] )
!
!                  # get Richardson number at half level given value of lower and upper full levels:
!                  Ri_hlv = self.Richardson_Number( z[ihlv-1:ihlv], 
!                                                   u[ihlv-1:ihlv], v[ihlv-1:ihlv], 
!                                                   T[ihlv-1:ihlv], s[ihlv-1:ihlv] )
!
!                  # which regime?
!                  regime = self.GetRegime( L, Ri_hlv, blh, z_hlv[ihlv] )
!
!                  # swich:
!                  if regime == 'M-O' :
!
!                      # following Monin-Obukov similarities for heat transfer:
!                      Kz_hlv[ihlv] = self.MoninObukov_Kz_heat( Ri_hlv, lmix, dU_dz_hlv[ihlv] )
!
!                  elif regime == 'Louis' :
!
!                      # following Louis scheme for heat transfer:
!                      Kz_hlv[ihlv] = self.Louis_Kz_heat( Ri_hlv, lmix, dU_dz_hlv[ihlv] )
!
!                  elif regime == 'EDMF' :
!
!                      # following EDMF scheme ;
!                      # use only surface driven diffusion, no in-cloud diffusion yet:
!                      Kz_hlv[ihlv] = self.EDMF_Kz_sfc_heat( L, u_star, blh, z_hlv[ihlv] )
!
!                  else :
!
!                      print( 'unsupported regime "%s"' % regime )
!                      raise Exception
!
!                  #endif # regime
!
!              #endfor # sfc/top or upper half-level
!
!          #endfor # half levels           
!
!      #enddef Halflevel_Kz
!
!      # *
!
!      def GetRegime( self, L_MO, Ri, z_entr, z ) :
!
!          """
!          Set the boundary layer regime following Figure 3.1 .
!
!          Arguments:
!            L_MO     : Monin-Obukov lenght (only sign is used)
!            Ri       : bulk-Richardson number (only sign is used)
!            z_entr   : entrainment level [m]
!            z        : height abovef surface [m]
!          """
!
!          # stable or unstable?
!          if L_MO >= 0 :
!              # stable surface layer; check Richardson number:
!              if Ri > 0 :
!                  # use Louis scheme:
!                  regime = 'Louis'
!              else :
!                  # use Monin-Obukov similarity:
!                  regime = 'M-O'
!              #endif
!          else :
!              # unstable surface layer; check height:
!              if z < z_entr :
!                  # within mixed layer, use EDMF scheme:
!                  regime = 'EDMF'
!              else :
!                  # above mixed layer; check Richardson number:
!                  if Ri > 0 :
!                      # use Louis scheme:
!                      regime = 'Louis'
!                  else :
!                      # use Monin-Obukov similarity:
!                      regime = 'M-O'
!                  #endif
!              #endif # within mixed layer ?
!          #endif  # stable or unstable
!
!          # ok
!          return regime
!
!      #enddef SetRegime
!
!      # *
!
!      def ScaleParameters( self, u_star, sshf, slhf, Tn, rho ) :
!
!          """
!          Compute scale parameters s* and q* ;
!          these have the same function in vertical profiles 
!          for heat and moisture as u* has for velocity.
!
!          From eq. (3.7):
!            rho u* s* = Js   (=sshf [J/m2/s] )
!            rho u* q* = Jq   (~slhf)
!
!          Arguments:
!            u_star   : scale parameter for wind profile (friction velocity) [m/s]
!            sshf     : surface sensible heat flux [W/m2 = J/s/m2]
!            slhf     : surface latent   heat flux [W/m2 = J/s/m2]
!            Tn       : reference temperature [K] taken as near-surface temperature
!                       (here temperature of the lowest model level n)
!            rho      : air density [(kg air)/m3]
!
!          Return values:
!            s_star   : scale parameter for dry static energy [J/(kg air)]
!            q_star   : scale parameter for specific humidity [(kg water)/(kg air)]
!          """
!
!          # scale parameter for dry static energy:
!          #        J/m2/s  / (kg air)/m3 /  m/s
!          s_star =  sshf   /     rho     / u_star   # J/(kg air)
!
!          # scale parameter for specific humidity:
!          #        J/m2/s / ( J/(kg water)/K  K  ) / (kg air)/m3 /  m/s
!          q_star =  slhf  / (    self.c_p   * Tn ) /      rho    / u_star  # (kg water)/(kg air)
!
!          # ok
!          return s_star,q_star
!
!      #enddef ScaleParameters
!
!      # *
!
!      def VirtualTemperatureFlux( self, u_star, s_star, q_star, Tn ) :
!
!          """
!          Virtual temperature flux following eq. (3.8):
!
!                  u* s*
!            Q0v = ----- + epsilon Tn u* q*
!                   c_p
!
!          Parameters:
!             c_p     : specific heat at constant presure of moist air [J/kg/K]
!             epsilon : Rvap/Rdry - 1 [1]
!
!          Arguments:
!            u*    : friction velocity [m/s]
!            s*    : dry static energy profile scale parameter [J/(kg air)]
!            q*    : specific humidity profile scale parameter [(kg water)/(kg air)]
!            Tn    :  reference temperature [K] taken as near-surface temperature
!                     (here temperature of the lowest model level n)
!
!          Return value:
!            Q0v      : virtual temperature flux [K m/s]
!          """
!
!          # compute following eq. (3.8) :
!          #       m/s      J/kg  / J/kg/K          1         K     m/s      kg/kg
!          return u_star * s_star / self.c_p + self.epsilon * Tn * u_star * q_star  # K m/s
!
!      #enddef VirtualTemperatureFlux
!
!      # *
!
!      def MoninObukovLength( self, u_star, Tn, Q0v ) :
!
!          """
!          Monin-Obukov length following eq. (3.8):
!
!                     (u*)^3
!            L =  - -----------
!                   kappa g
!                   ------- Q0v
!                     Tn
!
!          Parameters:
!            kappa    : Von Karman constant (=0.4)
!            g        : gravity acceleration [m/s2]
!
!          Arguments:
!            u*    :  friction velocity [m/s]
!            Tn    :  reference temperature [K] taken as near-surface temperature
!                     (here temperature of the lowest model level n)
!            Q0v   :  virtual temperature flux in the surface layer [K m/s]
!
!          Return values:
!            L     :  Monin-Obukov length [m]
!          """
!
!          # compute L:
!          #          m3/s3   / (       1       m/s2  / K   K m/s )
!          return - u_star**3 / ( self.kappa * self.g / Tn * Q0v  )  # m
!
!      #enddef MoninObukovLength
!
!      # *
!
!      def BottomWindspeed2( self, u_n, v_n, T_n, Q0v ) :
!
!          """
!          Absolute windspeed (squared) at the bottom
!          model layer following eq. (3.17):
!
!             |U_n|^2 = u_n^2 + v_n^2 + w_*^2
!
!          with w_* the free convection velocity scale defined by:
!
!                         g       1/3
!            w_* = ( z_i --- Q0v )
!                        T_n
!
!          Note that Q0v could be negative (which gives L>0, stable conditions).
!          Therefore it is better to formulate the convection velocity scale as:
!
!                 3        g
!            (w_*)  = z_i --- Q0v
!                         T_n
!          and
!                    |      g      | 2/3
!            w_*^2 = | z_i --- Q0v |
!                    |     T_n     |
!
!          Parameters:
!             z_i   : scale height of the boundary layer depth (1000 m)
!             g     : acceleration of gravity
!
!          Arguments:
!             u_n, v_n  : horizontal wind components at lowest model layer
!             T_n       : temperature                at lowest model layer
!             Q0v       : virtual temperature flux in the surface layer [K m/s ?]
!          """
!
!          # free convection velocity scale:
!          #               m         m/s2 /  K    Km/s
!          w_star = abs( self.z_i * self.g / T_n * Q0v )**(1.0/3.0) # m/s
!
!          # squared absolute windspeed at lowest model layer:
!          return u_n**2 + v_n**2 + w_star**2   # J/kg
!
!      #enddef BottomWindspeed2    
!
!      # *
!
!      def DryStaticEnergy( self, z, T ) :
!
!          """
!          Compute dry static energy following eq. (3.5):
!             s = g z + c_p T
!
!          Parameters:
!             c_p     : specific heat at constant presure of moist air [J/kg/K]
!             g       : accelartion of gravity [m/s2]
!
!          Arguments:
!             z   : height above surface [m]
!             T   : temperature [K]
!
!          Return value:
!             s   : dry static energy [J/kg]
!          """
!
!          # compute:
!          #       m/s2    m    J/kg/K    K
!          return self.g * z + self.c_p * T
!
!      #enddef DryStaticEnergy
!
!      # *
!
!      def Richardson_Number( self, z, u, v, T, s ) :
!
!          """
!          Richardson number at half level following eq. (3.52);
!          note that in IFS level order 'k' is upper and 'k+1' is lower:
!
!                                            {(Delta sv)/(c_p T)}_{k+1/2}
!             Ri_{k+1/2} = g (z_k - z_{k+1}) ----------------------------
!                                                |Delta_U|_{k+1/2}^2
!
!          with the length of wind vector change between the layers:
!
!             |Delta_U|_{k+1/2}^2 = (u_k - u_{k+1})^2 + (v_k - v_{k+1})^2
!
!          and:
!
!              Delta sv                   2( s_k - s_{k+1} )
!             (--------)_{k+1/2} = ---------------------------------- + epsilon (q_k - q_{k+1})
!               c_p T              (s_k - g z_k + s_{k+1} - g z_{k+1}
!
!          The dry static energy s is defined in eq. (3.5) as:
!
!            s = g z + c_p T
!
!          which simplifies the above formula to a form in better agreement with
!          the notation:
!
!              Delta sv                s_k - s_{k+1}
!             (--------)_{k+1/2} = --------------------- + epsilon (q_k - q_{k+1})
!               c_p T              c_p (T_k + T_{k+1})/2
!
!          Rename 'k+1' to 'low', 'k+1/2' to 'hlv', and 'k' to 'upp' to avoid confusion:
!
!                                        {(Delta sv)/(c_p T)}_hlv
!             Ri_hlv = g (z_upp - z_low) ------------------------
!                                             |Delta_U|_hlv^2
!
!             |Delta_U|_hlv^2 = (u_upp - u_low)^2 + (v_upp - v_low)^2
!
!              Delta sv            s_upp - s_low
!             (--------)_hlv = --------------------- + epsilon (q_upp - q_low)
!               c_p T          c_p (T_upp + T_low)/2
!
!          Parameters:
!             c_p     : specific heat at constant presure of moist air [J/kg/K]
!             epsilon : R_vap / R_dry - 1 [1]
!             g       : accelartion of gravity  [m/s2]
!
!          Arguments:
!             z[2]    : height            of lower and upper layer [m]
!             u[2]    : eastward  wind    in lower and upper layer [m/s]
!             v[2]    : northward wind    in lower and upper layer [m/s]
!             T[2]    : temperature       in lower and upper layer [K]
!             q[2]    : specific humidity in lower and upper layer [kg/kg]
!             s[2]    : dry static energy in lower and upper layer [K]
!
!          """
!
!          # length of wind vector change:
!          Delta_U_2 = (u[1] - u[0])**2 + (v[1] - v[0])**2  # m2/s2
!
!          # energy term:
!          #             J/kg         J/kg/K        K                  1           kg/kg
!          Ds_cpT = (s[1] - s[0]) / (self.c_p * (T[1]+T[0])/2) + self.epsilon * (q[1]-q[0])  # 1
!
!          # combine:
!          #       m/s2         m            1    /   m2/s2
!          return self.g * (z[1] - z[0]) * Ds_cpT / Delta_U_2  # 1
!
!      #enddef Richardson_Number
!
!      # *
!
!      def Phi_M( self, zeta ) :
!
!          """
!          Dimensionless gradient function for momentum.
!          Following eq. (3.19) :
!            Phi_M(zeta) = ( 1 - 16 zeta)^{-1/4}
!
!          Arguments:
!            zeta = z / L < 0 with:
!              z   : height above surface [m]
!              L   : Monin-Obukov length [m] (<0 for unstable conditions)
!          """
!
!          # check ..
!          if zeta >= 0 :
!              print( 'expecting zeta<0, but argument is:', zeta )
!              raise Exception
!          #endif
!
!          # shear function:
!          return ( 1.0 - 16 * zeta )**(-1.0/4.0)
!
!      #enddef Phi_M
!
!      # *
!
!      def Phi_H( self, zeta ) :
!
!          """
!          Dimensionless gradient function for heat (and moisture).
!          Following eq. (3.19) :
!            Phi_M(zeta) = ( 1 - 16 zeta)^{-1/2}
!
!          Arguments:
!            zeta = z / L < 0 with:
!              z   : height above surface [m]
!              L   : Monin-Obukov length [m] (<0 for unstable conditions)
!          """
!
!          # check ..
!          if zeta >= 0 :
!              print( 'expecting zeta<0, but argument is:', zeta )
!              raise Exception
!          #endif
!
!          # shear function:
!          return ( 1.0 - 16 * zeta )**(-1.0/2.0)
!
!      #enddef Phi_H
!
!      # *
!
!      def Phi_M0( self, zeta ) :
!
!          """
!          Dimensionless gradient function for momentum.
!          Following eq. (3.42) :
!            Phi_M0(zeta) = ( 1 - 15 zeta)^{-1/3}
!
!          Arguments:
!            zeta = z / L < 0 with:
!              z   : height above surface [m]
!              L   : Monin-Obukov length [m] (<0 for unstable conditions)
!          """
!
!          # check ..
!          if zeta >= 0 :
!              print( 'expecting zeta<0, but argument is:', zeta )
!              raise Exception
!          #endif
!
!          # shear function:
!          return ( 1.0 - 15 * zeta )**(-1.0/3.0)
!
!      #enddef Phi_M0
!
!      # *
!
!      def Phi_H0( self, zeta ) :
!
!          """
!          Dimensionless gradient function for heat (and moisture).
!          Following eq. (3.42) :
!            Phi_H0(zeta) = ( 1 - 39 zeta)^{-1/3}
!
!          Arguments:
!            zeta = z / L < 0 with:
!              z   : height above surface [m]
!              L   : Monin-Obukov length [m] (<0 for unstable conditions)
!          """
!
!          # check ..
!          if zeta >= 0 :
!              print( 'expecting zeta<0, but argument is:', zeta )
!              raise Exception
!          #endif
!
!          # shear function:
!          return ( 1.0 - 39 * zeta )**(-1.0/3.0)
!
!      #enddef Phi_H0
!
!      # *
!
!      def mixing_length( self, lamb, z ) :
!
!          """
!          Mixing length scale used in the surface layer,
!          bounded by an asymptotic length scale lambda.
!          Following Eq. (3.53):
!            1      1         1
!            - = ------- + ------
!            l   kappa z   lambda
!          thus:
!            l = 1 / ( 1/(kappa z) + 1/lambda )
!
!          Parameters:
!            kappa    : Von Karman constant (0.4)
!
!          Arguments:
!            lambda   : assymptotic lenght scale [m]
!            z        : height above surface [m]
!
!          Return value:
!            lmix     : mixing length [m]
!          """
!
!          # mixing length scale:
!          return 1.0 / ( 1.0/(self.kappa * z) + 1.0/lamb )
!
!      #enddef mixing_length
!
!      # *
!
!      def MoninObukov_Kz_heat( self, Ri, lmix, dU_dz ) :
!
!          """
!          Diffussion coefficient for heat using Monin-Obukov similarity;
!          for unstable layers with Ri<0. Following eq. (3.54):
!
!                     l^2      | dU |
!            K_H = ----------- | -- |
!                  Phi_M Phi_H | dz |
!
!          with elements:
!            l       : mixing length [m] computed from z and an asymptotic mixing length;
!            Phi_M   : gradient function for momentum following eq. (3.19) with zeta=Ri
!            Phi_H   : gradient function for heat     following eq. (3.19) with zeta=Ri
!
!          Arguments:
!            Ri      : Richardson number (<0) [1]
!            lmix    : mixing length scale [m], using asymptotic length scale lambda=150m
!            dU_dz   : vertical gradient of horizontal windspeed U [(m/s)/m]
!
!          Return value:
!            Kz      : diffusion coefficient [m2/s]
!          """
!
!          # check ..
!          if Ri >= 0 :
!              print( 'expecting Ri<0, but argument is:', Ri )
!              raise Exception
!          #endif
!
!          # evaluate gradient functions using zeta=Ri:
!          Phi_M = self.Phi_M( Ri )
!          Phi_H = self.Phi_H( Ri )
!
!          # K_H from eq. (3.54):
!          #        m2    / (  1       1  )  (m/s)/m
!          return lmix**2 / (Phi_M * Phi_H) * dU_dz   # m2/s
!
!      #enddef MoninObukov_Kz_heat
!
!      # *
!
!      def Louis_Kz_heat( self, Ri, lmix, z, dU_dz ) :
!
!          """
!          Diffussion coefficient for heat using Louis scheme;
!          for stable layers with Ri>0. Following eq. (3.55):
!
!                                  | dU |
!            K_H = l^2 f_LTG_H(Ri) | -- |
!                                  | dz |
!
!          with elements:
!            l       : mixing length [m] computed from z and an asymptotic mixing length;
!            f_LTG_H : function by Louis/Tiedtke/Geleyn following (3.56)
!
!          Arguments:
!            Ri      : Richardson number (<0) [1]
!            lmix    : mixing length scale [m], using asymptotic length scale lambda=150m
!            dU_dz   : vertical gradient of horizontal windspeed U [(m/s)/m]
!
!          Return value:
!            Kz      : diffusion coefficient [m2/s]
!          """
!
!          # check ..
!          if Ri <= 0 :
!              print( 'expecting Ri>0, but argument is:', Ri )
!              raise Exception
!          #endif
!
!          # function parameters:
!          b = 5.0
!          d = 1.0
!          # evaluate:
!          f_LTG_H = 1.0 / ( 1.0 + 2*b*Ri * (1.0 + d*Ri)**(1.0/2.0) )  # 1
!
!          # K_H from eq. (3.55):
!          #        m2          1    (m/s)/m
!          return lmix**2 * f_LTG_H * dU_dz   # m2/s
!
!      #enddef Louis_Kz_heat
!
!      # *
!
!      def EDMF_Kz_sfc_heat( self, L, u_star, blh, z ) :
!
!          """
!          Diffussion coefficient for surface driven diffustion of heat
!          in unstable boundary layer in the EDMF scheme.
!          Following eq. (3.41):
!
!                        kappa u*
!            K_sfc_H = ----------- z (1-z/blh)^2
!                      Phi_H0(z/L)
!
!          Note: In the original version of the IFS manual, the factor 'z' was omitted.
!          For the original equations, see Eq. (1) in (Troen and Mahtr, 1986) via:
!            http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.461.9396&rep=rep1&type=pdf
!
!          Arguments:
!            L       : Monin-Obukov length [m]
!            u_star  : scale parameter of velocity profile (friction velocity) [m/s]
!            blh     : boundary layer height [m]
!            z       : height above surface [m]
!
!          Return values:
!            Kz      : diffusion coefficient [m2/s]
!          """
!
!          # check ..
!          if z >= blh :
!              print( 'expecting z <= blh, but arguments are ', z, blh )
!              raise Exception
!          #endif
!
!          # K_sfc_H from eq. (3.41), with extra factor z :
!          #           1        m/s   /      1             m         1
!          return self.kappa * u_star / self.Phi_M0(z/L) * z * ( 1.0 - z/blh )**2  # m2/s
!
!      #enddef EDMF_Kz_sfc_heat
!
!
!  #endclass Diffusion
!
