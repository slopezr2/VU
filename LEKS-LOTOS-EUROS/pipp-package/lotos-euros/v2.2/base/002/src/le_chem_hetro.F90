!###############################################################################
!
! Introduction
! ------------
!
! Heterogeneous reaction rate of tracer GAS with the water around
! an aerosol AER.
! Tracer GAS is one of :
!   - N2O5
!   - HNO3
! Aerosol AER is one of : 
!   - seasalt (actually NaCl)
!   - (NH4)HSO4 
!
! Rate depends on the probability that a gas molecule meets
! an aerosol of a certain size. The number density distribution 
! as a function of the aerosol size is important therefore.
! Here a the shape of the number density distribution is kept
! constant (median and width), the actual number per size bin
! depends on the amount of aerosol 'moles' per m3 air.
! Result is a reaction rate with unit [1/(ppb AER)/s] .
!
! Since computation is expensive, reaction rates are stored
! in a look-up-table with indices for pressure and increment
! factor due to water uptake.
!
!
! Reaction summary
! ----------------
!
! In terms of atoms:
!
!    (NH4)SO4a (s) + H2O (l) + N2O5 (g)  ->  (NH4)SO4a (s) + 2 HNO3 (g)
!
!     2 NaCl (s)   + H2O (l) + N2O5 (g)  ->   2 NaCl (s)   + 2 NO3a- (aq) + 2 H+ (aq)
!
!       NaCl (s)   + H2O (l) + HNO3 (g)  ->     NaCl (s)   +   NO3a- (aq) +   H+ (aq) + H2O(l)
!
! In terms of reactants, products, and katalysators;
! reaction rates is expressed for this:
!
!    (NH4)SO4a (s) + N2O5 (g)  ->  (NH4)SO4a (s) + 2 HNO3 (g)
!
!     seasalt (s)  + N2O5 (g)  ->   seasalt (s)  + 2 NO3a- (aq)
!
!     seasalt (s)  + HNO3 (g)  ->   seasalt (s)  +   NO3a- (aq)
!
!
! Method
! ------
!
! Reactions:
!
! * {NaCl (s) + H2O} + HNO3 (g) + air (g) <--> NaNO3 (s) + HCl (g) + air (g)
!       particle        a          b
!
! * {NaCl (s) + H2O} + N2O5 (g) + air (g) [<--> NaCl (s) + 2HNO3 (g) + air (g) ]
!      particle         a          b    [ intermediate product is directly converted into]
!
!                                          <--> 2NaNO3 (s) + 2HCl (g) + air (g)
!
! * {(NH4)HSO4 (s) + H2O} + N2O5 (g) + air (g) <--> (NH4)HSO4 (s) + 2HNO3 (g) + air (g)
!        particle            a          b
!
! Kn = labda_ab / aerosol_radius          = Knudsen number [-]
!                                           where:  labda_ab        = mean free path [m]
!                                                   aerosol radius  = radius of aerosol particle AER [m]
!
! labda_ab = 3 * D_ab / c_a               = mean free path (average distance traveled by molecule 'a' before it encounters another molecule ('a' or) 'b' [m]
!                                           where:  D_ab            = diffusivity of gas 'a' in 'b' [m2 s-1]
!                                                   c_a             = mean speed of gas molecules 'a' [m s-1]
!
! c_a = [(8*R*T)/(pi*M_a)]**0.5           = mean speed of gas molecules 'a' [m s-1]
!                                           where:  R               = molar gas constant [J mol-1 K-1 == kg m2 s-2 mol-1 K-1]
!                                                   T               = temperature [K]
!                                                   pi              = 3.14159
!                                                   M_a             = molecular weight of 'a' [kg mol-1]
!
!        1e-3 * (temp**1.75) * (1/(M_a*1e3)+1/(M_b*1e3))**0.5
! D_ab = ----------------------------------------------------*1e-4  = diffusivity of gas 'a' in 'b' [m2 s-1]
!            pres_atm  *  ( (v_a)**(1/3) + (v_b)**(1/3) )**2
!                                           where:  temp            = temperature [K]
!                                                   M_a             = molecular weight of 'a' [kg mol-1]
!                                                   M_b             = molecular weight of 'b' [kg mol-1]
!                                                   pres_atm        = pressure [atm]
!                                                   v_a             = atomic (molecular) diffusion volume of 'a' [?] Fuller, Schettler, Giddings relation (1966) from Perry's Chemical Engineers' Handbook (1984)
!                                                   v_b             = atomic (molecular) diffusion volume of 'b' [?] Fuller, Schettler, Giddings relation (1966) from Perry's Chemical Engineers' Handbook (1984)
!
!
! History
!   1996-04, Fank Dentener, IMAU
!     Original code for TM3.
!     Use Whitby sulfate distribution, and Fuchs' rate expression
!     to integrate rate coefficient on aerosol distribution.
!     Dentener (1993) Ph.D. thesis
!   2001-05, Martijn Schaap, IMAU/TNO
!     Included in LOTOS.
!   2011-06, Roy Wichink Kruit, TNO
!     Extension to other hetregeneous reactions.
!   2011-09, Arjo Segers, TNO
!     Generalized implementation for LOTOS-EUROS v1.7.7 .
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Chem_Hetro

  use GO, only : gol, goPr, goErr
  
  use Binas, only : xm_H, xm_N, xm_O, xm_S, xm_Na, xm_seasalt

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  nreac_het
  public  ::  ireac_N2O5_NH4HSO4a_f
  public  ::  ireac_N2O5_ss_f
  public  ::  ireac_N2O5_ss_c
  public  ::  ireac_HNO3_ss_f
  public  ::  ireac_HNO3_ss_c
  public  ::  gamma_het
  public  ::  LUT_Rk_Het
  
  public  ::  LE_Chem_Hetro_Init, LE_Chem_Hetro_Done
  public  ::  LUT_Rk_Het_Get
  public  ::  increment_factor_aerh2o
  public  ::  increment_factor_rh


  ! --- const ------------------------------------

  ! module name:
  character(len=*), parameter  ::  mname = 'LE_Chem_Hetro'


  ! ***

  ! number of hetrogenious reactions:
  integer, parameter     :: nreac_het = 5
  ! indices:
  integer, parameter     :: ireac_N2O5_NH4HSO4a_f = 1
  integer, parameter     :: ireac_N2O5_ss_f   = 2
  integer, parameter     :: ireac_N2O5_ss_c   = 3
  integer, parameter     :: ireac_HNO3_ss_f   = 4
  integer, parameter     :: ireac_HNO3_ss_c   = 5

  ! ***
  
  ! set parameters for the heterogeneous routine
  real, parameter     :: gamma_het(nreac_het) = &
           (/ 0.05, &       !               HNO3    formation       on        sulphate aerosol
              0.05, &       ! direct fine   nitrate formation: N2O5 on fine   sea salt aerosol
              0.10, &       ! direct coarse nitrate formation: N2O5 on coarse sea salt aerosol
              0.05, &       !        fine   nitrate formation: HNO3 on fine   sea salt aerosol
              0.10 /)       !        coarse nitrate formation: HNO3 on coarse sea salt aerosol

  ! atomic diffusion volumes
  ! (Fuller, Schettler, Giddings relation (1966) 
  ! from Perry's Chemical Engineers' Handbook (1984))
  real, parameter  :: v_H    = 1.98   ! H atom
  real, parameter  :: v_N    = 5.69   ! N atom
  real, parameter  :: v_O    = 5.48   ! O atom
  ! molecular diffusion volumes 
  ! (Fuller, Schettler, Giddings relation (1966) 
  ! from Perry's Chemical Engineers' Handbook (1984))
  real, parameter  :: v_air  = 20.1   ! air

  ! *  properties of N2O5
  
  real, parameter  :: v_N2O5   =  v_N * 2 +  v_O * 5  ! molecular diffusion volume

  ! * properties of HNO3
  
  real, parameter  :: v_HNO3   =  v_H +  v_N +  v_O * 3 ! molecular diffusion volume

  ! * properties of sea salt

  ! constants for aerosol size distribution of sea salt:
  ! ~ fine mode  : mean radius = 0.1 um; log(sigma)~log(1.9)=0.279
  !   (Jaenicke R., 1993, Ed. Hobbs)
  real,     parameter    :: apar_ss_f(3)  = (/1., 0.1e-6, 0.279/)
  ! ~ coarse mode: mean radius = 1.0 um; log(sigma)~log(2.0)=0.301 
  !  (Haldis, 2009)
  real,     parameter    :: apar_ss_c(3)  = (/1., 1.0e-6, 0.301/)

  !   cons_ss = max( c(:,:,:,i_Na_f) * Na_to_SeaSalt, 1e-6 )
  !   cons_ss = max( c(:,:,:,i_Na_c) * Na_to_SeaSalt, 1e-6 )

  !! relative increase of aerosol radius due to water uptake
  !finc = increment_factor_rh( conc_AER(i,j,k), rho_AER, rh(i,j,k) )

  ! * proprties of (NH4)HSO4

  ! aerosol size distribution of nh4hso4: 
  !   mean radius = 0.034 um; log(sigma)~log(2)=0.301
  real, parameter  :: apar_NH4HSO4a_f(3) = (/1., 0.034e-6, 0.301/)

  !! nh4hso4 concentrations:
  !conc_nh4hso4 = max( c(i,j,k,i_so4a)*xm_nh4hso4/xm_so4, 1.e-6 )

  !! relative increase of aerosol radius due to water uptake
  !finc = increment_factor_aerh2o( conc_AER(i,j,k), rho_AER, aerh2o(i,j,k), h2odens)

  ! * look-up-table axis

  ! number of intervals for integration over pressure range (11 in range 1-0 atm)
  integer, parameter    ::  pres_interval = 11

  !! pressure [Pa]
  !real, parameter       ::  pres          = 1.0e5

  ! number of intervals for integration over relative increase of aerosol radius due to water uptake
  integer, parameter   :: finc_interval = 13
  ! relative increase of aerosol radius due to water uptake intervals
  real, parameter      :: finc_table(finc_interval) = (/1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3., 4., 5., 10., 15./)


  ! --- type --------------------------------------
  
  type T_LUT_Rk_Het
    ! lookup table values:
    real, allocatable       ::  value(:,:)
  end type T_LUT_Rk_Het
  
  
  ! --- var --------------------------------------
  
  ! lookup tables:
  type(T_LUT_Rk_Het)        ::  LUT_Rk_Het(nreac_het)



contains


  ! ====================================================================


  subroutine LE_Chem_Hetro_Init( status )
  
    use Binas, only : xm_N2O5
    use Binas, only : xm_HNO3
    use Binas, only : xm_NH4HSO4, rho_NH4HSO4a
    use Binas, only : xm_seasalt, rho_seasalt

    ! --- in/out ---------------------------------
    
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Chem_Hetro_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! fill table for N2O5 + fine (NH4)HSO4 aerosol :
    call LUT_Rk_Het_Init( LUT_Rk_Het(ireac_N2O5_NH4HSO4a_f), &
                               gamma_het(ireac_N2O5_NH4HSO4a_f), &
                               xm_N2O5, v_N2O5, &
                               xm_NH4HSO4, rho_NH4HSO4a, apar_NH4HSO4a_f, &
                               status )
    IF_NOTOK_RETURN(status=1)
  
    ! fill table for N2O5 + fine seasalt aerosol :
    call LUT_Rk_Het_Init( LUT_Rk_Het(ireac_N2O5_ss_f), &
                               gamma_het(ireac_N2O5_ss_f), &
                               xm_N2O5, v_N2O5, &
                               xm_seasalt, rho_seasalt, apar_ss_f, &
                               status )
    IF_NOTOK_RETURN(status=1)
  
    ! fill table for N2O5 + coarse seasalt aerosol :
    call LUT_Rk_Het_Init( LUT_Rk_Het(ireac_N2O5_ss_c), &
                               gamma_het(ireac_N2O5_ss_c), &
                               xm_N2O5, v_N2O5, &
                               xm_seasalt, rho_seasalt, apar_ss_c, &
                               status )
    IF_NOTOK_RETURN(status=1)
  
    ! fill table for HNO3 + fine seasalt aerosol :
    call LUT_Rk_Het_Init( LUT_Rk_Het(ireac_HNO3_ss_f), &
                               gamma_het(ireac_HNO3_ss_f), &
                               xm_HNO3, v_HNO3, &
                               xm_seasalt, rho_seasalt, apar_ss_f, &
                               status )
    IF_NOTOK_RETURN(status=1)
  
    ! fill table for HNO3 + coarse seasalt aerosol :
    call LUT_Rk_Het_Init( LUT_Rk_Het(ireac_HNO3_ss_c), &
                               gamma_het(ireac_HNO3_ss_c), &
                               xm_HNO3, v_HNO3, &
                               xm_seasalt, rho_seasalt, apar_ss_c, &
                               status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Chem_Hetro_Init


  ! ***


  subroutine LE_Chem_Hetro_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Chem_Hetro_Done'

    ! --- local ----------------------------------
    
    integer       ::  ireac

    ! --- begin ----------------------------------
    
    ! loop over hetrogeneous reactions:
    do ireac = 1, nreac_het
      call LUT_Rk_Het_Done( LUT_Rk_Het(ireac), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! ok
    status = 0

  end subroutine LE_Chem_Hetro_Done


  ! ====================================================================


  !*************************************************************************
  ! 1 particle (unity)/cm3, radius and log(sigma) measurements from Whitby
  ! the radius is assumed to be 'dry radius
  ! We take aerosol size from Whitby accumulation mode (1978)
  !  Numbermean radius: 0.034um, sigma=2, 1 (unity) particles cm-3
  ! Molecular weight NH4HSO4 115 g/mol
  ! aerosol density of dry NH4HSO4 1.8 E3 kg/m3= 1.8 gcm-3
  ! temperature is not a determining factor is implicitly accounted for
  ! as a function of pressure.
  ! temperature is assumed to follow an adiabatic lapse rate:
  ! (T2/T1)=(P2/P1)^{(g-1)/g} function of pressure with g=Cp/Cv=ca. 1.40
  !*************************************************************************

  subroutine LUT_Rk_Het_Init( lut, gamma, xm_GAS, v_GAS, &
                                     xm_AER, rho_AER, apar_AER, &
                                     status )

    use Binas, only : xm_air, Rgas, pi, twopi

    ! --- in/out ---------------------------------
    
    type(T_LUT_Rk_Het), intent(out)     ::  lut
    real, intent(in)                    ::  gamma
    real, intent(in)                    ::  xm_GAS  ! molecular weight [kg mol-1]
    real, intent(in)                    ::   v_GAS  ! molecular diffusion volume
    real, intent(in)                    ::  xm_AER    ! aerosol molecular weight [kg mol-1]
    real, intent(in)                    ::  rho_AER    ! aerosol density [kg m-3]
    real, intent(in)                    ::  apar_AER(3) ! aerosol size distribution parameters
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LUT_Rk_Het_Init'

    real, parameter      :: vent          = 1.0        ! what is this? seems to be useless
    real, parameter      :: fln10         = 2.302585   ! Constant used for 10-based logarithms in reaction rates => ln(x)=2.3025log(x)
    real, parameter      :: pres_atm_std  = 1.         ! standard pressure == 1 [atm]
    real, parameter      :: temp_std      = 288.       ! standard temperature == 288 [K]
    real, parameter      :: hcr           = 1.40       ! heat capacity ratio == cp/cv == 1.4
    real, parameter      :: conc          = 1e-9       ! volume mixing ratio of 1 ppbv (mol mol-1)
    integer, parameter   :: nt            = 4          ! quantities of integration (e.g.number, surface, volume and rate coefficient)
    integer, parameter   :: nint1         = 2000       ! number of integration intervals over aerosol size distribution (1000 in range 0-1 um; 1000 in range 1-101 um)

    ! --- local ----------------------------------

    integer   ::  ip
    real      ::  pres_atm
    real      ::  temp
    real      ::  d_GAS
    real      ::  c_a
    real      ::  labda_ab
    real      ::  rho_b
    real      ::  aervol
    integer   ::  m
    integer   ::  iaero
    real      ::  rx1, rx2
    real      ::  fn1,fn2
    real      ::  fr1,fr2
    real      ::  fa1,fa2
    real      ::  fv1,fv2
    real      ::  zlogs
    real      ::  rmean
    real      ::  xfac
    real      ::  xkn
    real      ::  xlab
    real      ::  qi(nt)
    real      ::  qi30

    ! integration stepsizes [m],0-1 um,1-101 um
    real      :: rint(nint1)

    ! --- begin ----------------------------------
    
    ! init table:
    allocate( lut%value(finc_interval,pres_interval) )

    ! fill integegration stepsizes;
    !~ first half: 1000 steps of 0.001 um between 0-1 um
    do m = 1, nint1/2
      rint(m) = 0.001E-6
    end do
    !~ seconc half: 1000 steps of 0.1 um between 1-101 um
    do m = nint1/2+1, nint1
      rint(m) = 0.1E-6
    end do

    ! pressure from 1 to 0 atm
    do ip = 1,pres_interval

      ! atmosphere (minimum is 0.001 atm) GAS rename because of confusion with pres = 100000 and see next point
      pres_atm = max(0.001,1.1-ip*0.1)
      ! this estimate of temp is not very accurate [K]
      temp     = max(210.,temp_std*(pres_atm/pres_atm_std)**((hcr-1)/hcr))
      ! diffusivity of GAS in air [m2 s-1]
      ! (Fuller, Schettler, Giddings relation (1966) from Perry's Chemical Engineers' Handbook (1984))
      d_GAS    = 1.e-3*(temp**1.75)*((1./(xm_air*1.e3)+1./(xm_GAS*1.e3))**0.5)/(pres_atm*((v_air)**(1./3.)+v_GAS**(1./3.))**2.)*1.e-4
      ! mean molecular speed [m s-1] (Maxwell-Boltzmann speed distribution)
      c_a      = sqrt((8.*Rgas*temp)/(pi*xm_GAS))
      ! free molecular path length [m]
      labda_ab = 3.*d_GAS/c_a
      ! air density of gas 'b' [mol m-3]
      ! rho_b = (atm)/(RT) = (10^5 N m-2)/(8.314 N m mol-1 K-1 * K) = mol m-3
      rho_b    = pres_atm/(Rgas*temp)*1.e5
      ! dimensionless aerosol volume fraction:
      ! (mol aerosol mol-1 air) *(mol air m-3 air) * (kg aerosol mol-1 aerosol) / (kg aerosol m-3 aerosol) => [m3 aerosol m-3 air]
      ! aervol is the volume of 1 pbbv dry AER at temp and pres_atm
      aervol   = conc*rho_b*xm_AER/rho_AER
      
      ! dummy value to avoid warnings:
      qi30 = 1.0

      ! aerosol increased radius loop (1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3., 4., 5., 10., 15.)
      do iaero=1,finc_interval
        ! initialize rx1
        rx1  = 0.0
        ! initialize qi
        qi(:)= 0.0
        ! mean aerosol radius including the increase of aerosol radius due to water uptake [m]
        rmean= apar_AER(2) * finc_table(iaero)
        zlogs= apar_AER(3)  ! log(sigma)
        ! actual integration
        ! loop over 1 to nint1 (number of integration intervals)
        do m=1,nint1-1
          ! radius 1
          rx1=rx1+rint(m)
          ! radius 2
          rx2=rx1+rint(m+1)
          ! Knudsen number
          xkn= 3.*d_GAS/(rx1*c_a)
          xlab=(xkn*4./3.+0.71)/(xkn+1.)
          fn1=(log10(rx1/rmean))**2
          fn2=(log10(rx2/rmean))**2
          ! number integration
          fn1= exp(-fn1/(2.*zlogs**2))/rx1
          fn2= exp(-fn2/(2.*zlogs**2))/rx2
          ! surface integration
          fa1= rx1*rx1*fn1
          fa2= rx2*rx2*fn2
          ! volume integration
          fv1= rx1*rx1*rx1*fn1
          fv2= rx2*rx2*rx2*fn2
          ! reactivity integration
          fr1=1./(1.+xkn*(xlab+(4.*(1.-gamma)/3./gamma)))*rx1*fn1
          fr2=1./(1.+xkn*(xlab+(4.*(1.-gamma)/3./gamma)))*rx2*fn2
          ! EULER INTEGRATION
          qi(1)=qi(1)+rint(m)/2.*(fn1+fn2)
          qi(2)=qi(2)+rint(m)/2.*(fa1+fa2)
          qi(3)=qi(3)+rint(m)/2.*(fv1+fv2)
          qi(4)=qi(4)+rint(m)/2.*(fr1+fr2)
        enddo ! m=1,nint1-1

        ! constant integration factor
        xfac=apar_AER(1)*1.e6/fln10/sqrt(twopi)/zlogs
        ! conversion cm3=>m3 number
        qi(1)=qi(1)*xfac*1.e-6
        ! m=>cm surface
        qi(2)=qi(2)*4.*pi*xfac*1.e-2
        ! volume
        qi(3)=qi(3)*4./3.*pi*xfac
        ! reactivity at dry volume qi30
        qi(4)=qi(4)*4.*pi*d_GAS*xfac*vent
        ! store dry volume for smallest radius:
        if ( iaero == 1 ) qi30 = qi(3)
        ! lookup table of removal coefficients corrected to the real dry aerosol k (aervol/qi30)
        lut%value(iaero,ip) = aervol/qi30*qi(4)
      end do !iaero
    end do !ip

    ! ok
    status = 0

  end subroutine LUT_Rk_Het_Init


  ! ***


  subroutine LUT_Rk_Het_Done( lut, status )

    ! --- in/out ---------------------------------

    type(T_LUT_Rk_Het), intent(inout)   ::  lut
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LUT_Rk_Het_Done'

    ! --- begin ----------------------------------

    ! clear:
    deallocate( lut%value )

    ! ok
    status = 0

  end subroutine LUT_Rk_Het_Done


  ! ***


  pure subroutine LUT_Rk_Het_Get( lut, finc, pres, rk )

    ! --- in/out ---------------------------------

    type(T_LUT_Rk_Het), intent(in)      ::  lut
    real, intent(in)                    ::  finc
    real, intent(in)                    ::  pres
    real(8), intent(out)                ::  rk

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LUT_Rk_Het_Get'

    ! --- local ----------------------------------

    integer  :: l
    integer  :: np1,np2
    integer  :: nr1,nr2
    real     :: px
    real     :: hx
    real     :: hp1,hp2

    ! --- begin ----------------------------------

    ! interpolation from lookup table:
    np1=min(pres_interval,1+nint(10.-pres/10000.))          ! pressure
    np1=max(1,np1)
    np2=min(pres_interval,np1+1)

    ! lower bound of finc array
    nr1=1
    do l=1,finc_interval
      if(finc.ge.finc_table(l)) nr1=l
    enddo

    ! upper bound of finc
    nr2=min(nr1+1,finc_interval)

    ! pressure interpolation fraction:
    px=((pres_interval-np1)*10000.-pres)/10000.

    ! aerosol increment interpolation fraction:
    hx=1.
    if (nr1.ne.nr2) hx=(finc-finc_table(nr1))/(finc_table(nr2)-finc_table(nr1))

    ! bi-linear interpolation:
    hp1 = px * lut%value(nr1,np2) + (1.0-px) * lut%value(nr1,np1)
    hp2 = px * lut%value(nr2,np2) + (1.0-px) * lut%value(nr2,np1)
    
    ! fill reaction rate:
    rk = hx * hp2 + (1.0-hx) * hp1     ! [s ppbv-1]

  end subroutine LUT_Rk_Het_Get



  ! ***
  

  !--------------------------------------------------------------------------------------------------------
  ! Return relative increase of aerosol radius due to water uptake
  ! by Roy Wichink Kruit @ TNO June 2011
  !
  !               |--| delta_r_water (increment of aerosol radius due to aerosol water)
  !            |-----| total_radius (= aerosol_radius * increment_factor)
  !            |--| aerosol radius
  !        ________
  !       /        \
  !      /   /---\  \
  !     /   /-----\  \
  !     \   \-----/  /
  !      \   \---/  /
  !       \________/
  !
  !
  ! aerosol_mass = aerosol_density * aerosol_volume 
  !              = aerosol_density * 4/3 * pi * aerosol_radius**3
  !
  ! water_mass   = water_density * water_volume     
  !              = water_density * ( 4/3 * pi * total_radius**3 -  4/3 * pi * aerosol_radius**3 )
  !              = water_density * ( 4/3 * pi * (aerosol_radius * increment_factor)**3 -  4/3 * pi * aerosol_radius**3 )
  !              = water_density * ( 4/3 * pi * aerosol_radius**3 ) * ( increment_factor)**3 - 1 )
  !              = water_density *   aerosol_mass/aerosol_density   * ( increment_factor)**3 - 1 )
  !
  ! f_inc3 = increment_factor**3 = (water_mass/water_density)*(aerosol_density/aerosol_mass) + 1
  ! increment_factor = [(water_mass/water_density)*(aerosol_density/aerosol_mass) + 1]**(1/3)
  !
  ! If we consider 1 m3 of air, the water_mass and the aerosol_mass can be replaced 
  ! by the water_concentration and aerosol_concentration respectively
  !
  !--------------------------------------------------------------------------------------------------------

  function increment_factor_aerh2o(aer_conc, aer_dens, h2o_conc, h2o_dens) result( finc )

    real, intent(in)  :: aer_conc   ! (kg aerosol)/(m3 air)
    real, intent(in)  :: aer_dens   ! (kg aerosol)/(m3 aerosol)
    real, intent(in)  :: h2o_conc   ! (kg water)/(m3 air)
    real, intent(in)  :: h2o_dens   ! (kg water)/(m3 water)
    real              :: finc

    real              :: f_inc3     ! third power of increment_factor

    ! see the comment above ...
    f_inc3 = (h2o_conc/h2o_dens)*(aer_dens/aer_conc) + 1

    ! return value:
    finc = f_inc3 ** (1.0/3.0)

  end function increment_factor_aerh2o
  
  ! *
  
  ! idem using relative humidity:
  
  function increment_factor_rh( aer_conc, aer_dens, rh ) result( finc )
  
    use Binas, only : pi
  
    ! --- in/out ---------------------------------
    
    real, intent(in)  ::  aer_conc   ! (kg aerosol)/(m3 air)
    real, intent(in)  ::  aer_dens   ! (kg aerosol)/(m3 aerosol)
    real, intent(in)  ::  rh         ! %
    real              ::  finc
    
    ! --- local ----------------------------------
    
    real              ::  Rhf
    real              ::  fwet, rwet, rdry
    
    ! --- begin ----------------------------------
    
    ! Rh in %; Rh-0.01 because of problems with dividing by 0
    !Rhf = ( Rh - 0.01 )/100.0   ! [0,1]
    ! original Rh is sometimes larger than 100% ...
    Rhf = ( min(max(0.0,Rh),100.0) - 0.01 )/100.0   ! [0,1]

    ! wetness factor:
    fwet = (4.0/3.7)*(((2.0-rhf)/(1.0-rhf))**(1./3.))
    if (fwet<1) fwet=1
    if (fwet>4) fwet=4

    ! radius of wet and dry aerosol:
    !    rho 4/3 pi R**3 = mass(per volume)
    rwet = ( (3.0*aer_conc) / (4.0*pi*aer_dens) )**(1./3.) * fwet
    rdry = ( (3.0*aer_conc) / (4.0*pi*aer_dens) )**(1./3.) 
    
    ! relative increase of aerosol radius due to water uptake
    finc = rwet / rdry
    
    ! ... thus 'finc' is the same as 'fwet' ?

  end function increment_factor_rh


end module LE_Chem_Hetro
