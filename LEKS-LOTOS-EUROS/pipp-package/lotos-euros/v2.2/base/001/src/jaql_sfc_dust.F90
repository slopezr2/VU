!#######################################################################
!
! JAQL_SFC_Dust
!   Joint Air Quality Library
!   Surface modelling
!   Dust parameterizations
!
! Reports:
!
!   For a full description of dust emisisons as used in LOTOS-EUROS,
!   see the report of the BOP project:
!
!   [Schaap 2009]
!     Schaap M; Manders AMM; Hendriks ECJ; Cnossen JM; Segers AJS; Denier van der Gon HAC;
!     Jozwicka M; Sauter F; Velders G; Matthijsen J; Builtjes PJH
!     "Regional modelling of particulate matter for the Netherlands"
!     PBL Report no.    500099008
!     http://www.pbl.nl/bibliotheek/rapporten/500099008.pdf
!
!   or the summary report:
!
!   [Schaap et al., 2010]
!    "Mineral Dust as a component of particulate matter"
!    http://www.pbl.nl/en/publications/2010/Mineral-Dust-component-particulate-matter
!
! References:
!
!   [Alfaro&Gomes 2001]
!     St\'ephae C. Alfaro, Laurent Gomes
!     Modeling mineral aerosol production by wind erosion:
!       Emission intensities and aerosol size distribution in source areas
!     JGR 106 D16, p 18.075-18.084
!
!   [Marticorena and Bergametti, 1995]
!     Modeling the atmoshperic dust cycle: 1. Design of a soil-derived dust emission scheme
!     JGR, vol. 100, no. D8, p. 16,415-16,430 .
!
!#######################################################################


module JAQL_SFC_Dust

  implicit none
  

  ! --- in/out --------------------------------
  
  private
  
  public    ::  nmode
  
  public    ::  mode_gmean, mode_gsigma, mode_massfrac
  public    ::  Friction_Velocity_Threshold
  public    ::  Friction_Velocity_Threshold_Factor_Soilwater
  public    ::  Friction_Velocity_Threshold_Factor_Soilwater_Demo
  public    ::  Friction_Velocity_Threshold_Scale_Roughness
  public    ::  JAQL_SFC_Dust_Emis


  ! --- const --------------------------------
  
  ! * aerosol size distribution
  
  ! number of aerosol modes:
  integer, parameter  ::  nmode = 3

  ! size distribution used by Zender et al. (2003)
  !real, parameter    :: mode_gmean (nmode) =   (/ 0.832e-6,    4.82e-6,    19.38e-6    /)  ! geometric mean mass diameter [m]
  !real, parameter    :: mode_gsigma(nmode) =   (/ 2.1,         1.9,        1.6         /)  ! geometric standard deviation [-]
  !real, parameter   :: mode_massfrac(nmode) =  (/ 0.036,       0.957,      0.007       /)  ! geometric mass fraction      [-]

  ! size distribution used by Claquin 1999
  !real, parameter    :: mode_gmean (nmode) =   (/ 0.011e-6,    2.54e-6,    42.1e-6     /)  ! geometric mean mass diameter [m]
  !real, parameter    :: mode_gsigma(nmode) =   (/ 1.89,        2.0,        2.13        /)  ! geometric standard deviation [-]
  !real, parameter   :: mode_massfrac(nmode) =  (/ 2.6e-6,      0.78,       0.22        /)  ! geometric mass fraction      [-]
  
  ! size distribution used by DEAD model (seems to be best if angstrom exponent is compared between LE and MODIS)
  ! best of these three options, but far from perfect (S.L. Janson, July 2014)
  real, parameter    :: mode_gmean (nmode) =    (/ 0.2e-6,      1.67e-6,    11.6e-6     /)  ! geometric mean mass diameter [m]
  real, parameter    :: mode_gsigma(nmode) =    (/ 1.75,        1.76,       1.7         /)  ! geometric standard deviation [-]
  real, parameter    :: mode_massfrac(nmode) =  (/ 0.0008,      0.0092,     0.99        /)  ! geometric mass fraction      [-]

contains


  ! ===============================================================

  
  !
  ! Lowest u* needed to relase a particle from the surface.
  !
  
  elemental function Friction_Velocity_Threshold( Dp, rho_p, rho_a ) result( ustar_thr_u )
  
    use Binas, only : grav  ! gravity acceleration (m/s2)
  
    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  Dp            ! particle diameter (m)
    real, intent(in)      ::  rho_p         ! particle density (kg/m3), i.e. 2.65 g/cm3
    real, intent(in)      ::  rho_a         ! air density (kg/m3), i.e. 0.00123 g/cm3
    real                  ::  ustar_thr_u   ! friction velocity threshold (m/s)
    
    ! --- const ----------------------------------
    
    real, parameter   ::  a = 1331.0   ! cm**(-x)
    real, parameter   ::  b = 0.38     ! -
    real, parameter   ::  x = 1.56     ! -
    
   ! --- local ----------------------------------
    
    real      ::  Dp_cm
    real      ::  grav_cms2
    real      ::  rho_p_gcm3
    real      ::  rho_a_gcm3
    real      ::  Re
    real      ::  K
    
    ! --- begin ----------------------------------
    
    ! In the parameterization all distances are in cm,
    ! including the gravity acceleration etc.
    ! Therefore convert to cm first:
    Dp_cm      = Dp   * 1e2    ! cm
    grav_cms2  = grav * 1e2    ! cm/s2
    rho_p_gcm3 = rho_p * 1e-3  ! g/cm3
    rho_a_gcm3 = rho_a * 1e-3  ! g/cm3
    
    ! estimate of Reynolds number; minium for Dp==0 is b :
    Re = a * Dp_cm**x + b   ! -
    
    ! factor:
    !            g/cm3        cm/s2      cm    g/cm3          
    K = sqrt( rho_p_gcm3 * grav_cms2 * Dp_cm / rho_a_gcm3 ) &   ! cm/s
        !       g/cm**0.5/s2  g/cm3        cm/s2      cm**2.5
        * sqrt( 1 + 0.006/( rho_p_gcm3 * grav_cms2 * Dp_cm**2.5 ) )  ! 1
    
    ! set threshold given Reynodls number;
    ! note that original formulation a test 0.03<Re<10 is used,
    ! but since Re >= b = 0.38 the situaion Re <= 0.03 never occures:
    if ( Re < 10.0 ) then
      !               1?   cm/s             1
      ustar_thr_u = 0.129 * K / sqrt( 1.928 * Re**0.092 - 1.0 )   ! cm/s
    else
      !               1?   cm/s             1
      ustar_thr_u = 0.129 * K * ( 1.0 - 0.0858*exp( - 0.0617*(Re-10.0) ) )  ! cm/s
    end if
    
    ! convert to m:
    ustar_thr_u = ustar_thr_u * 1.0e-2  ! m/s
    
    
    ustar_thr_u = ustar_thr_u * 0.66  ! Tegen correction
    
  
  end function Friction_Velocity_Threshold
  
  
  ! Soil moisture inhibition after Fecan et al. 1999 .
  !
  ! The minimum friction velocity u* that is required to elevate
  ! dust from the surface is described as the "small-scale friction velocity u*s"
  ! times a factor:
  !
  !  u*_threshold  =  u*_s  f_sw
  !
  ! The factor f_sw depends on the soil water content and the clay fraction.
  ! If the soil contains water, it will be more difficult to elevate
  ! particles from it and a higher u* is required.
  ! However, a soil with a lot of clay will absorp more water before
  ! this effect is reached; therefore, the amount of soil water needs
  ! to exceed a threshold depending on the clay fraction.
  !
  ! Sketch of the result:
  !
  !         |                               
  !    f_sw |                                   *
  !         |                          *
  !         |                    *
  !         |                *
  !         |              *
  !       1 +**************
  !         |
  !         +-------------+------------------------------+
  !         0       sw_threshold(clay)           sw ->   1
  !                        
  
  elemental function Friction_Velocity_Threshold_Factor_Soilwater( sw, clay_fraction ) result( f_sw )
  
    ! --- in/out -----------------------------
    
    real, intent(in)      ::  sw              ! soil water content (kg water)/(kg soil)
    real, intent(in)      ::  clay_fraction   ! 0-1
    real                  ::  f_sw
    
    ! --- in/out -----------------------------
    
    real        ::  sw_threshold  ! (kg water)/(kg soil)
    
    ! --- begin ------------------------------
    
    !
    ! NOTE: the original formulation is in % for both 
    !       clay_fraction and the result;
    !       this converts from/to [0,1] instead of %;
    !       actual values in [0.00,0.31] :
    !   sw_theshold = [ 0.0014 (100 c)**2 + 0.17 (100 c) ]/100
    !               = [     14 c**2       +    17  c     ]/100
    !               =     0.14 c**2       + 0.17 c
    !
    !! threshold soil moisture dependence upon clay fraction
    !sw_threshold = 0.14*(clay_fraction)**2 + 0.17*clay_fraction  ! [0,1]
    
    !! according Mokhtari et al. (2012) the factor 3 represents the best emission rates
    !sw_threshold = 3 * ( 0.14*(clay_fraction)**2 + 0.17*clay_fraction )
    
    ! FIXED: following eq. (7) in Mokhtari et al. (2012) :
    !   - factor should be 0.0014
    !   - truncation to interval
    !sw_threshold = min( max( 0.053, 3 * ( 0.17*clay_fraction + 0.0014*(clay_fraction)**2 ) ), 0.15 )
    
    ! Fixed again: original formulation is with clay_fraction in %, and result is in % too ...
    ! soil water threshold in [0,1] (actual range [0.00,0.31]):
    sw_threshold = ( 0.0014 * (clay_fraction*100.0)**2 + 0.17*(clay_fraction*100.0) )/100.0
    ! apply Mokhtari et al. scaling and min/max values:
    sw_threshold = min( max( 0.053, 3*sw_threshold ) , 0.15 )
    
    ! only increased value if water content above threshold:
    if ( sw > sw_threshold) then
        ! BUG? what is the factor 100 doing? parameterization in % ?
        !f_sw = ( 1.0 + 1.21 *( 100*(sw-sw_threshold) )**0.68 )**0.5
        ! Try corrected:
        !f_sw = sqrt( 1.0 + 1.21 * (sw-sw_threshold)**0.68 )
        ! Original was ok, parameterization is indeed in %,
        ! so need to convert the sw values first ...
        f_sw = ( 1.0 + 1.21 *( 100*(sw-sw_threshold) )**0.68 )**0.5
    else
        f_sw = 1.0
    endif
    
  end function Friction_Velocity_Threshold_Factor_Soilwater


  ! *
  
  subroutine Friction_Velocity_Threshold_Factor_Soilwater_Demo( fu, filename, status )
  
    ! --- in/out -----------------------------
    
    integer, intent(in)             ::  fu        ! file unit
    character(len=*), intent(in)    ::  filename  ! target file
    integer, intent(out)            ::  status
    
    ! --- const ------------------------------
    
    ! clay fractions similar to (Fecan et al., 1999, Fig. 5.):
    integer, parameter    ::  nclayfrac = 4
    real, parameter       ::  clayfracs(nclayfrac) = (/ 0.00, 0.18, 0.42, 1.00 /) ! [0-1]
    
    ! number of soil water values:
    integer, parameter    ::  nsw = 100
    
    ! --- local ------------------------------
    
    integer       ::  iclayfrac
    integer       ::  isw
    real          ::  sw
    real          ::  factor
    
    ! --- begin ------------------------------
    
    ! create new file:
    open( unit=fu, file=trim(filename), form='formatted', iostat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - opening unit ",i0," for file : ",a)') fu, trim(filename)
      status=1; return
    end if
    ! header:
    write (fu,'("clay fraction,gravimetric soil moisture,friction velocity factor")')
    write (fu,'("0-1,(kg water)/(kg soil),1")')
    
    ! loop over clay fractions:
    do iclayfrac = 1, nclayfrac
      ! loop over soil water fractions:
      do isw = 0, nsw
        ! current:
        sw = float(isw)/nsw
        ! evaluate:
        factor = Friction_Velocity_Threshold_Factor_Soilwater( sw, clayfracs(iclayfrac) )
        ! write record:
        write (fu,*) clayfracs(iclayfrac), sw, factor
      end do
    end do
  
    ! close:
    close( unit=fu, iostat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - closing unit ",i0," for file : ",a)') fu, trim(filename)
      status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine Friction_Velocity_Threshold_Factor_Soilwater_Demo


  !
  ! Required friction velocity for dust elevation following:
  !
  !   Marticorena and Bergametti (1995)
  !   Modeling the atmoshperic dust cycle: 1. Design of a soil-derived dust emission scheme
  !   JGR, vol. 100, no. D8, p. 16,415-16,430 .
  !
  ! A parameterization is defined for the ratio between the 
  ! small-scale friction velocity u*_s and the 
  ! threshold of the surface friction velocity u* :
  !
  !            u*_s           ln(z0/z0s)
  !   f_eff = ------ = 1 - ( ------------- )      (eq.17)
  !           u*_thr         ln(delta/z0s)
  !
  ! where:
  !
  !   z0             roughness length
  !   z0s            small-scale roughness length 
  !                     Van Loon et al. 2003         : 1.0e-5
  !                     Marticorena et al, JGR, 1997 : 4.0e-6\
  !                     LOTOS-EUROS                  : 30e-6
  !   delta/z0s  = a (X/z0s)**p         (Eq.18)
  !   a ~ 0.35
  !   p ~ 0.8
  !   X ~ 0.1 m      (in article everything in cm!)
  !
  ! Sketch of the result, see Marticorena 1995 Figure 2:
  !
  !       |
  !     1 *
  !       |    *
  !       |        *
  !       |             *
  ! f_eff |                 *
  !       |                     *
  !       |                          *
  !       +--------+--------+--------+----*---+
  !      1e-6    1e-5    1e-4     1e-3      1e-2
  !                                    z0 (m) -->
  !
  
  elemental function Friction_Velocity_Threshold_Scale_Roughness &
              ( z0, z0s ) result ( f_eff )
  
    ! --- in/out -----------------------------
    
    real, intent(in)      ::  z0
    real, intent(in)      ::  z0s
    real                  ::  f_eff
    
    ! --- const ------------------------------
    
    real, parameter  ::  a = 0.35
    real, parameter  ::  p = 0.8
    real, parameter  ::  X = 0.10  ! m , 10 cm in the article
    
    ! --- begin ------------------------------
    
    ! first guess value:
    f_eff = 1.0 - ( log(z0/z0S) / log(a*(X/z0s)**p) )
    
    ! at high z0 values f_eff becomes negative,
    ! yielding high emissions for the wrong reason;
    ! therefore set a minimum value:
    if ( f_eff < 1e-6 ) f_eff = 1e-6
    if ( f_eff > 1.0 ) f_eff = 1.0
    
  end function Friction_Velocity_Threshold_Scale_Roughness
  
  
  ! ***
  
  
  !
  ! Total dust flux
  !
  ! AJS: This should not be a function of issd !
  !   Used to compute average soil particle diameter for Alpha.
  !   Strange anayway, soil-sizes are defined as mass distribution,
  !   that is not the same as a diameter distribution ...
  !

  subroutine JAQL_SFC_Dust_Emis(issd, Dp, rho_s, rho_a, ustar, ustar_thr_fac, em )
    
    use JAQL_SFC_SOIL, only      :  ssd_gmean, ssd_mfrac
    ! --- in/out ---------------------------------

    real, intent(in)          ::  Dp              ! soil particle diameter (m)
    integer, intent(in)       ::  issd
    real, intent(in)          ::  rho_s           ! soil density (kg/m3)
    real, intent(in)          ::  rho_a           ! air density (kg/m3)
    real, intent(in)          ::  ustar           ! friction velocity (m/s)
    real, intent(in)          ::  ustar_thr_fac   ! friction velocity threshold factor 
    !real, intent(out)         ::  em(nmode)       ! emission per aerosol mode (kg/m2/s)
    real, intent(out)         ::  em              ! total emission

    ! --- const ----------------------------------

    real, parameter ::  Dd = 6.7e-6 ! [m]

    ! --- local ----------------------------------
    
    real            ::  ustar_thr_u
    real            ::  ustar_thr
    real            ::  Fh
    !real           ::  alpha(nmode)
    real            ::  alpha
    real            ::  ustar_t_Dd
    real            ::  Dp_s
    
    ! --- begin ----------------------------------
  
    ! uncorrect u* treshold:
    ustar_thr_u = Friction_Velocity_Threshold( Dp, rho_s, rho_a )  ! m/s

    ! total threshold:
    ustar_thr = ustar_thr_u * ustar_thr_fac  ! m/s

    ! horizontal sand flux:
    Fh = JAQL_SFC_Dust_Fh( rho_a, ustar, ustar_thr )  ! (kg air)/m/s

    ! sandblasting efficiency:
    ! ~ u* required to compute alpha based on 
    !   average aerosol particle diameter (Dd ~ 6.7 micrometer)
    ustar_t_Dd = Friction_Velocity_Threshold( Dd, rho_s, rho_a )
    ! average soil diameter:  
    !   CORRECT? gmean/gsigma is for mass distribution!
     Dp_s =   ssd_gmean(issd, 1) * ssd_mfrac(issd, 1) +  &
              ssd_gmean(issd, 2) * ssd_mfrac(issd, 2) +  &
              ssd_gmean(issd, 3) * ssd_mfrac(issd, 3)   ! [m]
    ! ~ compute efficiency:
    alpha = JAQL_SFC_Dust_Alpha( Dp_s, rho_s, rho_a, ustar, ustar_t_Dd ) ! (kg soil)/(kg air)/m

    ! combine:
    !    (kg air)/m/s    (kg soil)/(kg air)/m
    em =    Fh        *          alpha          ! (kg soil)/m2/s
    
  end subroutine JAQL_SFC_Dust_Emis
  

  ! *
  
  
  !
  ! Horizontal sand flux.
  !  
  
  elemental function JAQL_SFC_Dust_Fh( rho_a, ustar, ustar_thr ) result( Fh )
    
    use Binas, only : grav  ! m/s2
    
    ! --- in/out ---------------------------------

    real, intent(in)          ::  rho_a           ! air density [(kg air)/m3]
    real, intent(in)          ::  ustar           ! friction velocity (m/s)
    real, intent(in)          ::  ustar_thr       ! friction velocity threshold (m/s)
    real                      ::  Fh              ! (kg air)/m/s

    ! --- local ----------------------------------
    
    real      ::  R
    
    ! --- begin ----------------------------------
    
    ! only flux if u* above threshold:
    if ( ustar <= ustar_thr ) then
      ! no flux:
      Fh = 0.0
    else
      ! ratio between threshold and friction velocity:
      R = ustar_thr / ustar
      ! horizontal flux:
      !    (kg air)/m3 / m/s2   (m/s)**3
      Fh =     rho_a   / grav * ustar**3 * (1-R) * (1+R)**2  ! (kg air)/m/s
    end if
    
  end function JAQL_SFC_Dust_Fh
  
  
  ! *
  
  
  !
  ! Aerosol emission efficiency due to horizontal sand flux.
  ! Method acquired by using Shao et al., 1996; as shown in Mokhtari, 2012
  !  
  ! AJS: What is the factor 0.001 in beta ?
  !
  ! AJS: Why minimum value for alpha ?
  !
  
  pure function JAQL_SFC_Dust_Alpha( Dp_s, rho_s, rho_a, ustar, ustar_t ) result( alpha )
    
    use Binas,              only :  grav  ! m/s2
    
    ! --- in/out ---------------------------------------
    
    real, intent(in)            ::  Dp_s       ! soil particle diameter (m)
    real, intent(in)            ::  rho_s      ! soil density (kg/m3)
    real, intent(in)            ::  rho_a      ! air  density (kg/m3)
    real, intent(in)            ::  ustar      ! (m/s)
    real, intent(in)            ::  ustar_t    ! (m/s)
    real                        ::  alpha      ! (kg soil)/(kg air)/m
    
    ! --- const ---------------------------------------

    ! constant:
    real, parameter             ::  gamma = 2.5

    ! fixed average aerosol size diameter:
    real, parameter             ::  Dd = 6.7e-6 ! um

    ! --- local ---------------------------------------
    
    real                        ::  Dd_mm
    real                        ::  Dp_s_mm
    real                        ::  Beta
    
    ! --- begin ---------------------------------------
    
    ! convert:
    Dd_mm   = Dd   * 1e3  ! mm
    Dp_s_mm = Dp_s * 1e3  ! mm
    
    ! scale factor:
    Beta = max( 0.0, 0.001 * (log(Dp_s_mm) * 0.125e-4 + 0.328e-4) * exp(-140.7*Dd_mm + 0.37) )

    !    (kg soil)/(kg air)/m  = (kg soil)/m3 / (kg air)/m3    1       1     m/s2 /   (m/s)**2
    alpha = max( 1e-4, 2.0/3.0 *    rho_s     /     rho_a   * Beta * Gamma * grav / (ustar_t)**2)

  end function JAQL_SFC_Dust_Alpha


end module JAQL_SFC_Dust
