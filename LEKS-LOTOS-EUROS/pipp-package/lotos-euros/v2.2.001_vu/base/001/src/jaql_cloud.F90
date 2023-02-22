!#######################################################################
!
! JAQL_Cloud
!   Joint Air Quality Library
!   Cloud chemistry and scavening function
!
!
!#######################################################################


module JAQL_Cloud

  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  get_hplus
  public  ::  dissociation_henry
  public  ::  henry_func
  public  ::  Rh_Growth
  
  public  ::  CloudCoverOverhead
  

  ! --- const --------------------------------
  
  
  

contains
  

  ! ====================================================================
  
  !
  ! Calculation of hplus, from CAMx
  !
  ! Function initializes result [H+] with value -999.0 ("NO_DATA"),
  ! but this is immediately reset to a value of at least 1.0e-10  (ph=10).
  !
  ! Seems that the routine also returns negative [H+],
  ! as well as [H+] values corresponding to a pH << 0  (around -3, -4) .
  !
  
  pure function get_hplus( ch_so2, ch_so4, ch_nh3, ch_nh4, ch_hno3, ch_no3, co2, &
                               dkco2, dkso2, dknh3, dkhso3, dkh2o, &
                               hkco2, hkso2, hknh3, &
                               xl ) result ( hplus )
  
    !--- in/out ------------------------------------------------------------------

    real, intent(in)  :: ch_so2, ch_so4, ch_nh3, ch_nh4, ch_hno3, ch_no3, co2  ! concentration in ppb (co2 in ppm?)
    real, intent(in)  :: dkco2, dkso2, dknh3, dkhso3, dkh2o
    real, intent(in)  :: hkco2, hkso2, hknh3
    real, intent(in)  :: xl                    ! volume of water in cloud
    real              :: hplus                 ! concentration H+ (mol/l)

    !--- local ------------------------------------------------------------------

    integer, parameter :: maxiter = 10    ! maximal number of iterations
    real, parameter    :: NO_DATA = -999.0

    !--- local ------------------------------------------------------------------

    real           :: zph(1:maxiter+1)    ! pH of cloud water (for check only) xxx can be removed later 

    ! variables for iteration:
    integer :: iter     ! iteration count
    real    :: x1       ! H+ concentration due to strong acids
    real    :: x2       ! H+ concentration due to total NHx
    real    :: a2       ! H+ concentration due to SO2
    real    :: x3       ! combined dissolution and solubility constant for CO2 * co2 xxx ??
    real    :: a1       ! integration constant  ?
    real    :: z        ! ??
    real    :: a        ! first coefficient in quadratic equation for H+
    real    :: b        ! second coefficient in quadratic equation for H+
    real    :: d        ! third coefficient in quadratic equation for H+
    real    :: discr    ! discriminant of quadratic equation

    !--- begin ------------------------------------------------------------------

    ! Initialisation (NO DATA):
    hplus = NO_DATA
    zph   = NO_DATA

    ! Compute initial guess for H+ concentration (mol/l), only taking into account 
    ! the strong acids SO4a, NH4a, HNO3, NO3a and assuming that all mass of these 
    ! acids enter the cloud droplets:
    ! xxx if (2.0*ch_so4 - ch_nh4 + ch_hno3 + ch_no3)/xl) < 1 -> pH < 0
    ! xxx test 1e tijdstap -> pH = -0.12 (regel 19539 screen output)
    hplus  = max(1.e-10,(2.0*ch_so4 - ch_nh4 + ch_hno3 + ch_no3)/xl)
    zph(1) = -log10(hplus)

    ! Compute H+ concentration, iteratively for pH > 4.5.
    ! Now we take also into account the contributions of HSO3, SO3, NH4(particle), HCO3;
    ! iteration is needed, since the amount of mass that enters a cloud droplet depends
    ! on the pH that we are computing.
    do iter = 1,maxiter

       ! Only continue iteration, if pH > 4.52 (-log10(3e-5) = 4.52)
       if  (hplus .lt. 3e-5) then

          ! Compute initial H+ concentrations 
          ! x1: strong acids
          ! x2: total NHx
          ! a2: SO2
          x1 = (2.0*ch_so4 - ch_nh4 + ch_hno3 + ch_no3)/xl
          x2 = (ch_nh3 + ch_nh4)/xl
          a2 = ch_so2/xl

          ! x3 is combined dissolution and solubility constant for CO2 * co2 xxx ??
          x3 = dkco2*hkco2*co2

          ! xxx underlying formulas needed 
          ! a1: integration constant  ?
          !
          !      dkh2o        1
          ! a1 = ----- (1 + ----- ) 
          !      dknh3      hknh3
          !
          ! 
          !                  a2
          ! z  = ---------------------------------
          !      [H+]          1      dkhso3
          !      ----- {(1 + -----) + ------ + 1}
          !      dkso2       hkso2     [H+]
          !
          !
          a1 = dkh2o/dknh3*(1.+1./hknh3) 
          z  = a2/(hplus/dkso2*(1.+1./hkso2)+dkhso3/hplus+1.)

          !                               2
          ! Solve quadratic equation a(H+) + b(H+) + d = 0;
          !                   a1
          ! xxx with a = 1 + ----, b = -x1-z, d = -x3 - 2.*dkhso3*z
          !                   H+

          a      = 1.+x2/(a1+hplus)
          b      = -x1-z
          d      = -x3-2.*dkhso3*z
          discr  = max(0.,(b*b-4.*a*d))
          hplus = max(1.e-10,(-b+sqrt(discr))/(2.*a))

          ! pH (for check only)
          zph(iter+1) = -log10(hplus)
       endif
    end do  ! iter

    ! TEST VOOR HET ZETTEN VAN pH op 7
    ! staat nu weer uit: met deze test goede resultaten!!!
    ! nu gaan checken voor de pH-waarden

    !hplus=10**(-5.5)

    ! Write pH (check) xxx remove later:
    !if (hplus .eq. NO_DATA) then
    !   write(*,'(12(1x,f8.2))') zph,NO_DATA
    !else
    !   write(*,'(12(1x,f8.2))') zph,-log10(hplus)
    !endif

  end function get_hplus

  !-----------------------------------------------------------------------------
  !
  ! This routine calculates the Henry and dissociation equilibria for the cloud chemistry:
  ! usage is in calculating the pH-value of the cloud water and the SO2 -> SO4 transformation.
  !
  ! coefficients are declared in the module header
  ! dkh2o          ! dissociation constant water
  ! hkco2          ! dimensionless Henry's constant CO2 
  ! dkco2          ! Dissociation constant CO2
  ! hkso2          ! dimensionless Henry's constant for sulfur dioxide
  ! dknh3          ! dissociation constant ammonia
  ! hknh3          ! dimensionless Henry's constant ammonia
  ! hkh2o2         ! dimensionless Henry's constant for hydroperoxide
  ! hko3           ! dimensionless Henry's constant for ozone
  ! dkso2          ! Dissociation constant for SO2
  ! dkhso3         ! Dissociation constant for HSO3-
  ! ztr            ! Temperature-related parameter (1/Kelvin)
                   ! ztr = (1/T - 1/Tref)

  pure subroutine dissociation_henry( temp_loc, ztr, &
                                   dkco2, dkso2, dknh3, dkhso3, dkh2o, &
                                   hkco2, hkso2, hknh3, hko3, hkh2o2 )

    ! -- in/out ---------------------------------------------------------------

    real, intent(in)  :: temp_loc       ! Temperature in grid cell in Kelvin
    real, intent(out) :: ztr            ! Temperature-related parameter (1/Kelvin)
                                        ! ztr = (1/T - 1/Tref)
    real, intent(out)  :: dkco2, dkso2, dknh3, dkhso3, dkh2o
    real, intent(out)  :: hkco2, hkso2, hknh3, hko3, hkh2o2

    ! -- local ---------------------------------------------------------------

    real,parameter :: rg=0.08314        ! R/100, with R: gas constant xxx klopt dit? zo ja gebruik dan ook R/100
    real           :: rt                ! temp*rg
    real,parameter :: Tref=298.         ! reference temperature (K) ; 298 K = 25 C

    ! -- begin ---------------------------------------------------------------

    ztr = (1./temp_loc - 1./Tref)
    rt  = temp_loc*rg

    ! xxx references ??
    dkh2o  = 1.01e-14*exp(-6706.0 *ztr) ! H2O <=> Hplus + SO3--
    hkco2  = 3.4e-2*exp(2420.*ztr)      ! dimensionless Henry coefficient CO2
    dkco2  = 4.5E-7*exp(-1000.*ztr)     ! CO2aq <=> HCO3- + Hplus
    hkso2  = 1.2*exp(3120.*ztr)*rt      ! dimensionless Henry coefficient SO2
    dknh3  = 1.8e-5*exp(-450.*ztr)      ! NH3 <=> NH4+ + OH-
    hknh3  = 76.0*exp(3400.*ztr)*rt     ! dimensionless Henry coefficient NH3
    hkh2o2 = 8.3e4*exp(7400.*ztr)*rt    ! dimensionless Henry coefficient H2O2
    hko3   = 1.1e-2*exp(2300.*ztr)*rt   ! dimensionless Henry coefficient O3
    dkso2  = 1.7e-2*exp(2090.*ztr)      ! SO2 <=> HSO3- + Hplus
    dkhso3 = 6.6e-8*exp(1510.*ztr)      ! HSO3- <=> SO3-- + Hplus

  end subroutine dissociation_henry
  
  
  pure function henry_func( ispec, hlaw0, Tfact, T, pH ) result( hlaw )
    
    use indices, only : ispec_nh3, ispec_hno3, ispec_so2
    
    ! ---in/out-----------------------------
    integer, intent(in)   ::  ispec  ! Tracer
    real, intent(in)      ::  hlaw0  ! Baseline Henry's Law constant @298K (M/atm)
    real, intent(in)      ::  Tfact  ! temperature factor
    real, intent(in)      ::  T      ! ambient temperature (K)
    real, intent(in)      ::  pH     ! pH of liquid 
    
    real                  ::  hlaw   ! Adjusted Henry's Law constant (M/atm)
    
    ! ---local------------------------------
    
    real                  ::  diss1, diss2
    
    ! --- begin ----------------------------
    
    hlaw = hlaw0*exp(Tfact*(1./298. - 1./T))
    if ( ispec == ispec_nh3 ) then
       diss1 = 10.**(-189.1/T - 4.117)
       diss2 = 10.**(-5839.5/T - 9.7618*alog(T) + 61.206)
       hlaw = hlaw*(1. + (diss1/diss2)*10.**(-pH))
    elseif (ispec == ispec_HNO3 ) then
       diss1 = 15.4
       hlaw = hlaw*(1. + diss1/(10.**(-pH)))
    elseif (ispec == ispec_so2 ) then
       diss1 = 10.**(853./T)/54950.
       diss2 = 10.**(621.9/T)/1.897e+9
       hlaw = hlaw*(1. + diss1/(10.**(-pH)) + diss1*diss2/(10.**(-2.*pH)))
    endif
         
  end function
  
  !
  ! Increase of extinction efficency (or cross-section)
  ! with increasing relative humidity.
  ! Based on (Veefkind et al, 1996)
  !
  !    Veefkind J.P., J.C.H. van der Hage, and H. M. ten Brink, 1996. 
  !    Nephelometer derived and directly measured aerosol optical thickness 
  !    of the atmospheric boundary layer,
  !    Atmos. Res., 41, 217-228.
  !

  elemental real function Rh_Growth( rh )
  
    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  rh   ! relative humidity [%]
    
    ! --- const ----------------------------------
    
    ! fit parameters growth funtion
    real, parameter :: ap(0:6) = (/0.196581266 ,  0.100650483  , -0.00389645252 , 5.63596671e-5, &
                                   6.4988432e-8, -7.95507834e-9,  4.95445298e-11 /)
  
    ! --- begin ----------------------------------
    
    Rh_Growth = ap(0) + ap(1)*Rh    + ap(2)*Rh**2 + ap(3)*Rh**3 +  &
                        ap(4)*Rh**4 + ap(5)*Rh**5 + ap(6)*Rh**6

  end function Rh_Growth
  
  
  !
  ! Given profile of cloud cover (0-1 per cell), compute "cloud cover overhead", 
  ! i.e. fraction of cell which has any cloud cover in a layer above.
  ! Uses ECMWF's "maximum random overlap" approach, as described in:
  !
  !   IFS Documentation Cy38r1
  !     Part IV: Physical Processes
  !       Chapter 7: Clouds and large scale precipitation
  !         7.2.3 Defintion of the source and sink terms
  !
  ! Levels in profile are ordered surface to top.
  !
  ! Code based on TM5's "phys_cloudcover" module, reversed layer order.
  !
  
  subroutine CloudCoverOverhead( n, cc, cco, status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(in)   ::  n      ! number of levels
    real, intent(in)      ::  cc(n)  ! cloud cover profile (sfc->top) [0-1]
    real, intent(out)     ::  cco(n) ! cloud cover overhead (sfc->top) [0-1]
    integer, intent(out)  ::  status

    ! --- const ----------------------------------

    ! minimum cloud fraction to avoid division by zero:
    real       ::  cltres = 1.0e-4
    
    ! --- local -----------------------------------

    real      ::   zclear, zcloud
    integer   ::   k
    
    ! --- begin -----------------------------------

    ! init flags:
    zclear = 1.0
    zcloud = 0.0
    ! loop over levels from top to bottom:
    do k = n, 1, -1
      zclear = zclear * (1.0-max(cc(k),zcloud)) / (1.0-min(zcloud,1.0-cltres))
      zcloud = cc(k)
      cco(k) = 1.0 - zclear
    end do
    
    ! ok
    status = 0
    
  end subroutine CloudCoverOverhead

    

end module JAQL_Cloud
