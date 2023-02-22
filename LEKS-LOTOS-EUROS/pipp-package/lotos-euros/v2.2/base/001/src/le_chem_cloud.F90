!#################################################################
!
! NAME
!
!   LE_Chem_Cloud  -  in-cloud chemistry
!
!
! DESCRIPTION
!
!   Module cloudchem for cloud chemistry.
!   Aplies SO2 -> SO4 conversion in clouds, depending on pH.
!
!
! HISTORY
!
!   2011 sep, Eric van der Swaluw, RIVM
!     Original code.
!   2011 oct, Arjo Segers, TNO
!     Adapted for LOTOS-EUROS v1.7.9 .
!     Moved arrays from Dims module to here.
!     Added module init and done routines.
!   2013-05, Arjo Segers, TNO
!     Changed threshold for volume of water in cloud.
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') mname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Chem_Cloud

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------

  private
  
  public  ::  LE_Chem_Cloud_Init, LE_Chem_Cloud_Done
  public  ::  wetS
  
  
  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Chem_Cloud'

  logical, parameter :: fast_slow = .true.     ! false: fast oxidation (all SO2 -> SO4)
                                               ! true : explicit oxidation scheme

  real, parameter :: NO_DATA         = -999.0  ! NO DATA value

  real, parameter :: cloud_threshold = 0.05    ! threshold for cloud cover; 
                                               ! if cloud cover is below threshold,
                                               ! there is no SO2-> SO4 conversion 

  ! threshold for xl (volume of water in cloud)
  !real, parameter :: xl_threshold    = 1.0e-10 ! threshold for xl; if xl is below threshold,
  !                                             ! no SO2-> SO4 conversion xxx value oke ????
  real, parameter :: xl_threshold    = 1.0e-5  ! new value after problems with RACMO meteo;
                                               ! alternative would be to change the cloud conversion
  
  real,parameter  :: xm_air=28.964e-3          ! mass of air (kg mol-1) xxx use from Binas ?



contains


  ! ====================================================================


  subroutine LE_Chem_Cloud_Init( status )
  
    use Dims         , only : nx, ny, nz
    use Indices      , only : nspec
    use Indices      , only : i_so2, i_so4a_f, i_o3, i_h2o2, i_hno3, i_no3a_f, i_nh3, i_nh4a_f
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ---------------------------------

    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Chem_Cloud_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'icc', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'iclwc', status )
    IF_NOTOK_RETURN(status=1)

    ! check if all involved tracers are in use:
    if ( any( (/ i_so2, i_so4a_f, i_o3, i_h2o2, i_hno3, i_no3a_f, i_nh3, i_nh4a_f /) < 0 ) ) then
      write (gol,'("cloud chemistry requires that following tracers are enabled:")'); call goErr
      write (gol,'("  SO2, SO4a_f, O3, H2O2, HNO3, NO3a_f, NH3, NH4a_f")'); call goErr
      write (gol,'("  Define both cbm4 and sia in selected tracer groups")'); call goErr
      TRACEBACK; status=1; return
    end if


    ! ok
    status = 0

  end subroutine LE_Chem_Cloud_Init


  ! ***


  subroutine LE_Chem_Cloud_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Chem_Cloud_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine LE_Chem_Cloud_Done


  ! ***


  !
  ! cloud chemistry SO2 -> SO4;
  ! ref: CAMx xxx aanvullen
  !

  subroutine wetS( dt, c, status )

    use Dims   , only : nx, ny, nz
    use Indices, only : nspec
    use Dims   , only: pH
    use JAQL   , only: get_hplus, dissociation_henry
    
    use LE_Data      , only : LE_Data_GetPointer

    !--- in/out ------------------------------------------------------------------

    real, intent(in)     :: dt          ! time step (s) 
    real, intent(inout)  :: c(:,:,:,:)  ! (nx,ny,nz,nspec) concentrations 
                                        ! (gas in ppb, aerosol in ug/m3)
    integer, intent(out) :: status
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/wetS'
    
    real,parameter :: co2 = 320.0e-6       ! assumed CO2 concentration (320 ppm)

    !--- local ------------------------------------------------------------------

    ! washing below-cloud and in-cloud :
    real, allocatable         ::  washbc(:)
    real, allocatable         ::  washic(:)

    ! Concentrations in a grid cell
    real           :: ch_so2              ! concentration SO2 in ppb
    real           :: ch_hno3             ! concentration HNO3 in ppb
    real           :: ch_nh3              ! concentration NH3 in ppb
    real           :: ch_so4              ! concentration SO4 in ppb
    real           :: ch_nh4              ! concentration NH4 in ppb
    real           :: ch_no3              ! concentration NO3 in ppb
    real           :: ch_sulf             ! concentration sulfate in ppb ( = ch_so4 xxx is dat inderdaad hetzelfde?)
    real           :: ch_o3               ! concentration O3 in ppb
    real           :: ch_h2o2             ! concentration H2O2 in ppb

    ! meteo parameters in a grid cell
    real           :: zclwc               ! cloud liquid water density in m3 water/m3 air
    real           :: xl                  ! volume of water in cloud
    real           :: xliq                ! zclwc/(cloud cover)

    ! Dissociation and Henry coefficients in a grid cell
    real           :: ztr                 ! Temperature-related parameter (1/Kelvin)
                                          ! ztr = (1/T - 1/Tref)
    real           :: dkh2o               ! dissociation constant water
    real           :: hkco2               ! dimensionless Henry's constant CO2 
    real           :: dkco2               ! Dissociation constant CO2
    real           :: hkso2               ! dimensionless Henry's constant for SO2
    real           :: dkso2               ! Dissociation constant for SO2
    real           :: dknh3               ! dissociation constant ammonia
    real           :: hknh3               ! dimensionless Henry's constant ammonia
    real           :: hkh2o2              ! dimensionless Henry's constant for hydroperoxide
    real           :: hko3                ! dimensionless Henry's constant for ozone
    real           :: dkhso3              ! Dissociation constant for HSO3-

    real           :: hplus               ! concentration H+ (mol/l)
    integer        :: ix, iy, iz          ! grid cell indices
    !real           :: hsave(nx,ny,nz)     ! hsave om weg te schrijven
    
    ! meteo data:
    real, pointer        ::   temp(:,:,:)   ! (lon,lat,1)
    real, pointer        ::   dens(:,:,:)   ! (lon,lat,alt)
    real, pointer        :: icc(:,:,:)   ! (lon,lat,alt)        
    real, pointer        ::  iclwc(:,:,:)   ! (lon,lat,alt)                

    ! --- begin ---------------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 't', temp, status, check_units='K' )    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')     
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'icc', icc, status, check_units ='1')     
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'iclwc', iclwc, status, check_units ='kg/kg')     
    IF_NOTOK_RETURN(status=1)        

    ! washing:
    allocate( washbc(nspec) )
    allocate( washic(nspec) )

    ! Reset pH
    pH = NO_DATA

    ! Loop over grid cells in horizontal:
    
    do ix = 1,nx
       do iy = 1,ny
         ! cloud (tcc) is only use to compare to threshold, but not used any further..
         ! not available from NC meteo, so :
         ! temporary trick for now..
         ! Loop over vertical layers:
         do iz = 1,nz
  
          ! Check whether we have clouds:
          ! now test on 3D overhead cloud coverage:
          if ( icc(ix,iy,iz) > cloud_threshold ) then 

                ! Get meteo parameters in current grid cell;
                ! zclwc, xl, xliq (defined in module header):
                ! now use 3D overhead cloud coverage:
                call get_meteo( temp(ix,iy,iz), dens(ix,iy,iz), &
                                  iclwc(ix,iy,iz), icc(ix,iy,iz), &
                                  zclwc, xl, xliq )

                ! Check threshold for xl (xxx see also original code Sabine -> cwmin)
                if (xl > xl_threshold) then

                   ! Get concentrations in current grid cell (all in ppb);
                   ! ch_so2, ch_hno3, ch_nh3, ch_so4, ch_nh4, ch_no3, ch_sulf, ch_o3, ch_h2o2
                   ! (defined in module header):
                   call get_conc( c(ix,iy,iz,:), dens(ix,iy,iz), &
                                  ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2, &
                                  ch_so4, ch_sulf, ch_nh4, ch_no3 )

                   ! Compute dissociation and Henry coefficients;
                   ! dkh2o, hkco2, dkco2, hkso2, dknh3, hknh3, hkh2o2, hko3, dkso2, dkhso3
                   ! (defined in module header):
                   call dissociation_henry( temp(ix,iy,iz), ztr, &
                                             dkco2, dkso2, dknh3, dkhso3, dkh2o, &
                                             hkco2, hkso2, hknh3, hko3, hkh2o2 )

!                   ! Get H+ concentration (mol/l):
!                   call get_hplus( ch_so2, ch_so4, ch_nh3, ch_nh4, ch_hno3, ch_no3, co2, &
!                                   dkco2, dkso2, dknh3, dkhso3, dkh2o, &
!                                   hkco2, hkso2, hknh3, &
!                                   xl, hplus)
                   hplus = get_hplus(ch_so2, ch_so4, ch_nh3, ch_nh4, ch_hno3, ch_no3, co2, &
                                     dkco2, dkso2, dknh3, dkhso3, dkh2o, &
                                     hkco2, hkso2, hknh3, &
                                     xl)
                                     
                   pH(ix,iy,iz)=-log10(hplus)

                   ! Update concentrations due to wet chemistry sulfate production:
                   call sulfwetchem( dt, temp(ix,iy,iz), dens(ix,iy,iz), hplus, &
                                       c(ix,iy,iz,:), icc(ix,iy,iz), &
                                       ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2, &
                                       ch_so4, ch_sulf, ch_nh4, ch_no3, &
                                       dkco2, dkso2, dknh3, dkhso3, dkh2o, &
                                       hkco2, hkso2, hknh3, hko3, hkh2o2, &
                                       xl, xliq, ztr )

                endif ! xl > xl_threshold

            endif ! cloud > cloud_threshold
          enddo ! iz
       enddo ! iy
    enddo ! ix

    ! EXTRA WEGSCHRIJVEN VOOR pH
    !open(12,file='hplus.dat',status='unknown')
    !do ix = 1, nx
    !  do iy = 1, ny
    !    write(12,'(100(1pe18.10))')hsave(ix,iy,2),hsave(ix,iy,3)
    !  end do
    !end do
    !close(12)

    ! clear:
    deallocate( washbc )
    deallocate( washic )
    
    ! ok
    status = 0

  end subroutine wetS


  !-----------------------------------------------------------------------------
  !
  ! Test routine for cloud chemistry SO2 -> SO4;
  ! in this test all SO2 is instantaneously converted to SO4.
  !
  ! factor 0.0409*96.0 is used to convert from ppb to ug/m3 (96 is molecular weight SO4)          
  ! icc: cloud coverage in grid cell (fraction)
  !

  subroutine wetS_all( c, status )
  
    use Indices, only : i_so2, i_so4a_f
    use LE_Data, only : LE_Data_GetPointer

    real, intent(inout)  :: c(:,:,:,:)   ! (nx,ny,nz,nspec)  concentration
    integer, intent(out) :: status
    
    real, pointer        :: icc(:,:,:)   ! (lon,lat,alt)    

    call LE_Data_GetPointer( 'icc', icc, status, check_units ='1')     
    IF_NOTOK_RETURN(status=1)

    c(:, :, :, i_so2   ) = c(:, :, :, i_so2   ) * (1.0-icc(:, :, :))
    c(:, :, :, i_so4a_f) = c(:, :, :, i_so4a_f) + (icc(:, :, :)*c(:, :, :, i_so2))*0.0409*96.0
    
    ! ok
    status = 0

  end subroutine wetS_all


!  !-----------------------------------------------------------------------------
!  !
!  ! Get PH in cloud (for one grid cell)
!  ! ref: CAMx xxx aanvullen
!  !
!  ! from data:
!  ! icc: cloud coverage in grid cell (fraction)
!  !
!
!  pure subroutine get_hplus( ch_so2, ch_so4, ch_nh3, ch_nh4, ch_hno3, ch_no3, co2, &
!                               dkco2, dkso2, dknh3, dkhso3, dkh2o, &
!                               hkco2, hkso2, hknh3, &
!                               xl, hplus )
!  
!    use Indices, only : nspec
!
!    !--- in/out ------------------------------------------------------------------
!
!    real, intent(in)  :: ch_so2, ch_so4, ch_nh3, ch_nh4, ch_hno3, ch_no3, co2  ! concentration in ppb (co2 in ppm?)
!    real, intent(in)  :: dkco2, dkso2, dknh3, dkhso3, dkh2o
!    real, intent(in)  :: hkco2, hkso2, hknh3
!    real, intent(in)  :: xl                    ! volume of water in cloud
!    real, intent(out) :: hplus                 ! concentration H+ (mol/l)
!
!    !--- local ------------------------------------------------------------------
!
!    integer, parameter :: maxiter = 10    ! maximal number of iterations
!
!    !--- local ------------------------------------------------------------------
!
!    real           :: zph(1:maxiter+1)    ! pH of cloud water (for check only) xxx can be removed later 
!
!    ! variables for iteration:
!    integer :: iter     ! iteration count
!    real    :: x1       ! H+ concentration due to strong acids
!    real    :: x2       ! H+ concentration due to total NHx
!    real    :: a2       ! H+ concentration due to SO2
!    real    :: x3       ! combined dissolution and solubility constant for CO2 * co2 xxx ??
!    real    :: a1       ! integration constant  ?
!    real    :: z        ! ??
!    real    :: a        ! first coefficient in quadratic equation for H+
!    real    :: b        ! second coefficient in quadratic equation for H+
!    real    :: d        ! third coefficient in quadratic equation for H+
!    real    :: discr    ! discriminant of quadratic equation
!
!    !--- begin ------------------------------------------------------------------
!
!    ! Initialisation (NO DATA):
!    hplus = NO_DATA
!    zph   = NO_DATA
!
!    ! Compute initial guess for H+ concentration (mol/l), only taking into account 
!    ! the strong acids SO4a, NH4a, HNO3, NO3a and assuming that all mass of these 
!    ! acids enter the cloud droplets:
!    ! xxx if (2.0*ch_so4 - ch_nh4 + ch_hno3 + ch_no3)/xl) < 1 -> pH < 0
!    ! xxx test 1e tijdstap -> pH = -0.12 (regel 19539 screen output)
!    hplus  = max(1.e-10,(2.0*ch_so4 - ch_nh4 + ch_hno3 + ch_no3)/xl)
!    zph(1) = -log10(hplus)
!
!    ! Compute H+ concentration, iteratively for pH > 4.5.
!    ! Now we take also into account the contributions of HSO3, SO3, NH4(particle), HCO3;
!    ! iteration is needed, since the amount of mass that enters a cloud droplet depends
!    ! on the pH that we are computing.
!    do iter = 1,maxiter
!
!       ! Only continue iteration, if pH > 4.52 (-log10(3e-5) = 4.52)
!       if  (hplus .lt. 3e-5) then
!
!          ! Compute initial H+ concentrations 
!          ! x1: strong acids
!          ! x2: total NHx
!          ! a2: SO2
!          x1 = (2.0*ch_so4 - ch_nh4 + ch_hno3 + ch_no3)/xl
!          x2 = (ch_nh3 + ch_nh4)/xl
!          a2 = ch_so2/xl
!
!          ! x3 is combined dissolution and solubility constant for CO2 * co2 xxx ??
!          x3 = dkco2*hkco2*co2
!
!          ! xxx underlying formulas needed 
!          ! a1: integration constant  ?
!          !
!          !      dkh2o        1
!          ! a1 = ----- (1 + ----- ) 
!          !      dknh3      hknh3
!          !
!          ! 
!          !                  a2
!          ! z  = ---------------------------------
!          !      [H+]          1      dkhso3
!          !      ----- {(1 + -----) + ------ + 1}
!          !      dkso2       hkso2     [H+]
!          !
!          !
!          a1 = dkh2o/dknh3*(1.+1./hknh3) 
!          z  = a2/(hplus/dkso2*(1.+1./hkso2)+dkhso3/hplus+1.)
!
!          !                               2
!          ! Solve quadratic equation a(H+) + b(H+) + d = 0;
!          !                   a1
!          ! xxx with a = 1 + ----, b = -x1-z, d = -x3 - 2.*dkhso3*z
!          !                   H+
!
!          a      = 1.+x2/(a1+hplus)
!          b      = -x1-z
!          d      = -x3-2.*dkhso3*z
!          discr  = max(0.,(b*b-4.*a*d))
!          hplus = max(1.e-10,(-b+sqrt(discr))/(2.*a))
!
!          ! pH (for check only)
!          zph(iter+1) = -log10(hplus)
!       endif
!    end do  ! iter
!
!    ! TEST VOOR HET ZETTEN VAN pH op 7
!    ! staat nu weer uit: met deze test goede resultaten!!!
!    ! nu gaan checken voor de pH-waarden
!
!    !hplus=10**(-5.5)
!
!    ! Write pH (check) xxx remove later:
!    !if (hplus .eq. NO_DATA) then
!    !   write(*,'(12(1x,f8.2))') zph,NO_DATA
!    !else
!    !   write(*,'(12(1x,f8.2))') zph,-log10(hplus)
!    !endif
!
!  end subroutine get_hplus


  !**********************************************************************
  !
  !wetS - aqueous phase chemistry of sulfur  (and other)
  !programmed by Ad Jeuken (KNMI), Frank Dentener (IMAU) and Martijn Schaap (IMAU)
  !
  !purpose
  !-------
  !oxidation of SO2 and uptake of other gases in the aqueous phase
  !
  !interface
  !---------
  !call  wetS(c,dt,zclwc,y,zph,c4)
  !c     concentration array
  !dt    chemistry timestep
  !zclwc dimensionless liquid water content
  !zph   pH of cloudwater
  !c4    budget accumulator
  !
  !method
  !------
  !implicit solution of oxidation of SO2
  !
  !external
  !--------
  !none
  !
  !reference
  !---------
  !-
  !**********************************************************************

  pure subroutine sulfwetchem( dt, temp, dens, hplus, c, icc, &
                           ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2, &
                           ch_so4, ch_sulf, ch_nh4, ch_no3, &
                           dkco2, dkso2, dknh3, dkhso3, dkh2o, &
                           hkco2, hkso2, hknh3, hko3, hkh2o2, &
                           xl, xliq, ztr )
  
    use Indices, only : specmolm
    use Indices, only : nspec
    use Indices, only : i_o3, i_h2o2, i_so2, i_so4a_f, i_hno3

    !--- in/out ------------------------------------------------------------------

    real, intent(in)    :: dt                    ! time step (s)
    real, intent(in)    :: temp                  ! temperature (K)
    real, intent(in)    :: dens                  ! air density (kg/m3)
    real, intent(in)    :: hplus                 ! H+ concentration (mol/l)
    real, intent(in)    :: icc                   ! Fraction cloud in cell
    real, intent(inout) :: c(nspec)              ! concentration (gas in ppb, aerosol in ug/m3)
    real, intent(inout) :: ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2  ! gas in ppb
    real, intent(inout) :: ch_so4, ch_sulf, ch_nh4, ch_no3  ! aerosol in ppb
    real, intent(in)    :: dkco2, dkso2, dknh3, dkhso3, dkh2o
    real, intent(in)    :: hkco2, hkso2, hknh3, hko3, hkh2o2
    real, intent(in)    :: xl                    ! volume of water in cloud
    real, intent(in)    :: xliq                  ! zclwc/(cloud cover)
    real, intent(in)    :: ztr                   ! Temperature-related parameter (1/Kelvin)
                                                 ! ztr = (1/T - 1/Tref)

    !--- local  ------------------------------------------------------------------

    real,parameter     :: cwmin= 0.0015
    real,parameter     :: ccmin = 0.01
    integer, parameter :: nreacw = 3              ! number of wet reactions (xxx 2 reactions are used)
    integer, parameter :: KSO2HP = 1              ! index into rw for oxidation by H2O2
    integer, parameter :: KSO2O3 = 2              ! index into rw for oxidation by O3
    real,parameter     :: avo=6.03e20             ! Avogadro's number xxx is not used ?!?!

    real               :: x1,x2,x3,b1,b2,disc     ! help variables
    real               :: so2x                    ! SO2 concentration
    real               :: dso2                    ! changes in concentrations
    real               :: xso2o3a, xso2o3b        ! xxx ?

    real               :: phs4                    ! effective dissolvation of S(IV)
    real               :: phso2                   ! effective dissolvation of SO2
    real               :: phh2o2                  ! effective dissolvation of H2O2
    real               :: phozone                 ! effective dissolvation of O3

    real               :: xcov                    ! cloud cover xxx or icc ??

    real,dimension(nreacw) :: rw                  ! reaction rates

    !---------------------------------------------------------------------------------
    ! wet phase reactions
    ! S&P Seinfeld & Pandis: Atmospheric chemistry and physics, Wiley & sons, 1997
    !---------------------------------------------------------------------------------

    ! xxx FS check of er wolken zijn gebeurt in routine wetS;
    ! xxx die check heb ik hier weggelaten.
    ! xxx de checks zijn wel anders:
    ! xxx in wetS       : cloud > cloud_threshold = 0.1  en xl > xl_threshold = 1e-10
    ! xxx in sulfwetchem: cwc >= cwmin = 0.0015 
    ! xxx nog beslissen hoe we dat doen

    ! test write all concentrations (xxx can be removed later):
    ! xxx check: sulf2 en sulf4 moeten hetzelfde getal geven
    ! xxx (som van de concentraties in ppb voor en na de wetchem)
    ! ES write(*,'(a,8(1x,e12.4))') 'sulf1 ',c(i_so2),c(i_hno3),c(i_nh3),c(i_o3),c(i_h2o2),c(i_so4a),c(i_nh4a),c(i_no3a)
    ! ES write(*,'(a,8(1x,e12.4))') 'sulf2 ',c(i_so2)+c(i_o3)+c(i_h2o2)+ch_sulf ! mass balance (ppb)

    ! S(IV): total dissolved sulphur in oxidiation state 4:
    ! [S(IV)] = [SO2.H2O] + [HSO3-] + [SO3--] (S&P 6.39)
    !
    ! Compute phase factor ratio of aqueous phase to gas phase concentration;
    ! (S&P 6.40, with Ks1 = dkso2, Ks2 = dkhso3; see S&P 6.34, 6.35) 
    ! xxx afleiding it in de buurt, maar preciese afleiding met xliq ??
    ! xxx dimensies oke ??
    phs4    = hkso2 *(1. + dkso2/hplus + dkhso3*dkso2/hplus/hplus)*xliq
    phso2   = hkso2 *xliq
    phh2o2  = hkh2o2*xliq
    phozone = hko3  *xliq

    ! rate constants xxx afleiding?
    rw(KSO2HP) = 8e4*exp(-3560.*ztr)/(0.1 + hplus)
    XSO2O3A    = 4.39e11*exp(-4131/temp) + 2.56e3*exp(-966/temp)  !S(IV)
    XSO2O3B    = 2.56e3*exp(-966/temp)/hplus                    !divide by [H+]!S(IV)

    !  Make rate constants dimensionless by multiplying by (1./xliq/avo=6e20)
    !  multiply with the fractions of concentrations residing in the aqueous phase
    ! xxx afleiding ?
    rw(KSO2HP) = rw(KSO2HP)/xl*phso2/(1.+phs4)*phh2o2/(1.+phh2o2)
    rw(KSO2O3) = (XSO2O3A + XSO2O3B)/xl*phs4/(1.+phs4)*phozone/(1.+phozone)

    so2x=ch_so2
    if (fast_slow) then
       !
       ! oxidation of S(IV) by O3
       !
       so2x = ch_so2
       xcov = icc ! 3-dimensial cloud-coverage

       ! Compute dso2: change in [SO2], from differential equation:
       ! d[SO2]/dt = -Rw*[O3]*[SO2], with solution [SO2](t+dt) = [SO2](t)*exp(-Rw*[O3]*dt).
       ! From this follows d[SO2] = [SO2](t+dt) - [SO2](t) = [SO2](t)*[exp(-Rw*[O3]*dt)-1].
       ! Note that dso2 is always < 0 
       x1   = min(100.,rw(kso2o3)*ch_o3*dt)
       dso2 = ch_so2*xcov*(exp(-x1)-1.) ! only applied to xcov part of cloud
       dso2 = max(-ch_o3*xcov,dso2)     ! limit to O3 availability

       !! debug ...
       !print*,(-dso2/ch_so2)

       ! Update SO2, SO4 and O3 concentrations:
       ch_so2  = ch_so2  + dso2
       ch_sulf = ch_sulf - dso2
       c(i_o3) = ch_o3   + dso2

       !
       ! oxidation of S(IV) by H2O2
       !
       ! Here we explicitly solve the differential equation:
       ! y' = P-Q*y-R*y*y (P and Q are 0 => b3=0.)
       ! 
       ! b1*y^2 + b2*y + b3 = 0, with b1 = Rw, b2 = Rw*([H2O2]-[SO2]), b3 = 0 xxx afleiding ??

       so2x  = ch_so2                               ! note that ch_so2 has already been updated
                                                    ! due to oxidation by O3
       b1    = rw(kso2hp)
       b2    = b1*(ch_h2o2-so2x)
       disc  = min(100.,sqrt(b2*b2))                ! for b3=0, disc = b2 (disc = sqrt(discriminant))
       x1    = (b2-disc)/(-2.*b1)                   ! x1 = 0
       x2    = (b2+disc)/(-2.*b1)                   ! x2 = 2*b2/(-2Rw) = -b2/Rw
       x3    = (so2x-x1)/(so2x-x2)*exp(-disc*dt)    ! x3 = ([SO2]/([SO2]+(b2/Rw)))*exp(-b2*dt)
       if ( x3 == 1.0 ) then
         dso2 = 0.0
       else
         so2x  = (x1-x2*x3)/(1.-x3)                   ! xxx ?? new SO2 concentration
         dso2  = (so2x - ch_so2)*xcov
       end if
       dso2  = max(dso2,-ch_h2o2*xcov)              ! limit to H2O2 availability

       ! Update H2O2 concentration:
       ! hier zat factor /1.e+3 nog in weggehaald voor cchem011
       c(i_h2o2) = (ch_h2o2 + dso2)
    else
       ! TEST TEST TEST TEST TEST TEST
       ! maximal conversion of SO2 -> SO4 in clouds
       ! xxx compare with wetS_all
       xcov = icc ! 3-dimensional cloud-coverage
       dso2 = -xcov*ch_so2
    endif ! fast_slow

    ! Update SO2 (ppb) and SO4 concentrations again and convert aerosol concentration to ug/m3:
    c(i_so2)  =  ch_so2  + dso2
    ch_sulf   = (ch_sulf - dso2)
    c(i_so4a_f) =  ch_sulf * dens / xm_air * (specmolm(i_so4a_f))

    ! test write all concentrations (xxx can be removed later):
    ! ES write(*,'(a,8(1x,e12.4))') 'sulf3 ',c(i_so2),c(i_hno3),c(i_nh3),c(i_o3),c(i_h2o2),c(i_so4a),c(i_nh4a),c(i_no3a)
    ! ES write(*,'(a,8(1x,e12.4))') 'sulf4 ',c(i_so2)+c(i_o3)+c(i_h2o2)+ch_sulf ! mass balance (ppb)


    !--------------------------------------------------------
    ! xxx rest of routine not needed anymore? of toch meenemen?
    ! done somewhere else in RCG:
    !  c4(i,j,k,2) = c4(i,j,k,2) + dso2
    !
    ! NH3 uptake in cloud droplets is limited by H2SO4 availability
    ! no HNO3 is considered at this point
    ! assume instantaneous uptake of NH3 incloud  only in cloudy part
    !
    !sb  dnh3 = max((2.*ch_sulf-ch_nh4)*xcov,0.)
    !sb  dnh3 = max(-ch_nh3*xcov,-dnh3)
    !sb  c(i,j,k,i_nh3) = ch_nh3 + dnh3                 ! dnh3 is loss of NH3
    !sb  c(i,j,k,i_nh4) = (ch_nh4-dnh3)*0.0409*18
    !sb  c4(i,j,k,3) = c4(i,j,k,3) + dnh3
    !--------------------------------------------------------

  end subroutine sulfwetchem


  !-----------------------------------------------------------------------------
  ! Get concentrations in the current grid cell and
  ! convert aerosol concentrations from ug/m3 to ppb

  pure subroutine get_conc( c, dens, &
                       ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2, &
                       ch_so4, ch_sulf, ch_nh4, ch_no3 )

    use Indices, only : nspec
    use Indices, only : i_so2, i_o3, i_h2o2, i_hno3, i_so4a_f, i_nh3, i_no3a_f, i_nh4a_f
    use Indices, only : specmolm          ! mol mass of elements (kg/mol)

    ! -- in/out ---------------------------------------------------------------

    real, intent(in)  :: c(nspec)     ! concentration (gas in ppb, aerosol in ug/m3)
    real, intent(in)  :: dens         ! air density (kg/m3)
    real, intent(out) :: ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2  ! gas in ppb
    real, intent(out) :: ch_so4, ch_sulf, ch_nh4, ch_no3  ! aerosol in ppb

    ! -- local  ---------------------------------------------------------------

    ! gasesous components:
    ch_so2  = c(i_so2)
    ch_hno3 = c(i_hno3)
    ch_nh3  = c(i_nh3)
    ch_o3   = c(i_o3)
    ch_h2o2 = c(i_h2o2)

    ! aerosols:
    ch_so4  = c(i_so4a_f)/dens * xm_air / specmolm(i_so4a_f)
    ch_sulf = ch_so4
    ch_nh4  = c(i_nh4a_f)/dens * xm_air / specmolm(i_nh4a_f)
    ch_no3  = c(i_no3a_f)/dens * xm_air / specmolm(i_no3a_f)

  end subroutine get_conc

  !-----------------------------------------------------------------------------
  ! Get meteo parameters in the current grid cell.
  !
  ! Meteo parameters are declared in module header
  ! zclwc               ! cloud liquid water density in m3 water/m3 air
  ! xl                  ! volume of water in cloud
  ! xliq                ! zclwc/(cloud cover)
  !

  pure subroutine get_meteo( temp, dens, clwc, icc, &
                        zclwc, xl, xliq )

    ! -- in/out ---------------------------------------------------------------

    real, intent(in)  :: temp         ! temperature (K)
    real, intent(in)  :: dens         ! air density (kg/m3)
    real, intent(in)  :: clwc         ! cloud liquid water content (kg/kg)
    real, intent(in)  :: icc          ! cloud thickness (vertical fraction cell)
    real, intent(out) :: zclwc        ! cloud liquid water density in m3 water/m3 air
    real, intent(out) :: xl           ! volume of water in cloud
    real, intent(out) :: xliq         ! zclwc/(cloud cover)

    ! -- local  ---------------------------------------------------------------

    real, parameter :: R=8.31451           ! Gas constant xxx use from Binas ??

    ! -- begin  ---------------------------------------------------------------

    ! Compute zclwc = cloud liquid water density in m3 water/m3 air.
    ! The unit in the code for zlwc to be used is m3 water/m3 air;
    ! initially clwc is in kg/kg, so we have to multiply with 10^-3 assuming a 
    ! water density of 1000 kg/m3 and we have to multiply with dens to 
    ! account for the density of air:

    !ES clwc nemen we vast
    !zclwc = 3.e-4*1.e-3*dens
    zclwc = clwc*1.e-3*dens

    ! Compute xl = volume of water in cloud 
    ! xl is used to convert ppb (nmol of species/mol air) concentrations to 
    ! concentrations in water (mol H+/l) 
    ! hplus = c/xl -> xl = c/hplus -> [xl] = (nmol of species/mol air) / (mol H+/l)
    ! ideal gas law: p*V = n*R*T <=> V = n*R*T/p
    ! p = 10^5 hPa
    ! f1 = conversion factor ppb   = 10^-9 mol/mol air
    ! f2 = conversion factor liter = 10^-3 m3
    ! 1/(f1*f2*p) = 10^7
    ! xxx in sulfwetchem: xl = zclwc/cloud*R*temp*10e7 welke is correct??
    ! xxx in getph      : xl = zclwc*R*temp*10e7       welke is correct??
    xl = zclwc*R*1.e+7*temp 

    ! Compute xliq = zclwcoud
    xliq = zclwc

  end subroutine get_meteo

!  !-----------------------------------------------------------------------------
!  !
!  ! This routine calculates the Henry and dissociation equilibria for the cloud chemistry:
!  ! usage is in calculating the pH-value of the cloud water and the SO2 -> SO4 transformation.
!  !
!  ! coefficients are declared in the module header
!  ! dkh2o          ! dissociation constant water
!  ! hkco2          ! dimensionless Henry's constant CO2 
!  ! dkco2          ! Dissociation constant CO2
!  ! hkso2          ! dimensionless Henry's constant for sulfur dioxide
!  ! dknh3          ! dissociation constant ammonia
!  ! hknh3          ! dimensionless Henry's constant ammonia
!  ! hkh2o2         ! dimensionless Henry's constant for hydroperoxide
!  ! hko3           ! dimensionless Henry's constant for ozone
!  ! dkso2          ! Dissociation constant for SO2
!  ! dkhso3         ! Dissociation constant for HSO3-
!  ! ztr            ! Temperature-related parameter (1/Kelvin)
!                   ! ztr = (1/T - 1/Tref)
!
!  pure subroutine dissociation_henry( temp_loc, ztr, &
!                                   dkco2, dkso2, dknh3, dkhso3, dkh2o, &
!                                   hkco2, hkso2, hknh3, hko3, hkh2o2 )
!
!    ! -- in/out ---------------------------------------------------------------
!
!    real, intent(in)  :: temp_loc       ! Temperature in grid cell in Kelvin
!    real, intent(out) :: ztr            ! Temperature-related parameter (1/Kelvin)
!                                        ! ztr = (1/T - 1/Tref)
!    real, intent(out)  :: dkco2, dkso2, dknh3, dkhso3, dkh2o
!    real, intent(out)  :: hkco2, hkso2, hknh3, hko3, hkh2o2
!
!    ! -- local ---------------------------------------------------------------
!
!    real,parameter :: rg=0.08314        ! R/100, with R: gas constant xxx klopt dit? zo ja gebruik dan ook R/100
!    real           :: rt                ! temp*rg
!    real,parameter :: Tref=298.         ! reference temperature (K) ; 298 K = 25 C
!
!    ! -- begin ---------------------------------------------------------------
!
!    ztr = (1./temp_loc - 1./Tref)
!    rt  = temp_loc*rg
!
!    ! xxx references ??
!    dkh2o  = 1.01e-14*exp(-6706.0 *ztr) ! H2O <=> Hplus + SO3--
!    hkco2  = 3.4e-2*exp(2420.*ztr)      ! dimensionless Henry coefficient CO2
!    dkco2  = 4.5E-7*exp(-1000.*ztr)     ! CO2aq <=> HCO3- + Hplus
!    hkso2  = 1.2*exp(3120.*ztr)*rt      ! dimensionless Henry coefficient SO2
!    dknh3  = 1.8e-5*exp(-450.*ztr)      ! NH3 <=> NH4+ + OH-
!    hknh3  = 76.0*exp(3400.*ztr)*rt     ! dimensionless Henry coefficient NH3
!    hkh2o2 = 8.3e4*exp(7400.*ztr)*rt    ! dimensionless Henry coefficient H2O2
!    hko3   = 1.1e-2*exp(2300.*ztr)*rt   ! dimensionless Henry coefficient O3
!    dkso2  = 1.7e-2*exp(2090.*ztr)      ! SO2 <=> HSO3- + Hplus
!    dkhso3 = 6.6e-8*exp(1510.*ztr)      ! HSO3- <=> SO3-- + Hplus
!
!  end subroutine dissociation_henry
!

end module LE_Chem_Cloud
