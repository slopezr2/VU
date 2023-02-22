!###############################################################################
!
! Boundary-layer deposition resistance (Rb) following Zhang (2001).
!
! Not implemented yet:
!  o correction for vegetation height
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

module LE_DryDepos_Aero_Zhang

  use GO, only : gol, goPr, goErr
  
  implicit none


  ! --- in/out --------------------------------
  
  private
  
  public  ::  LE_DryDepos_Aero_Zhang_vd


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_DryDepos_Aero_Zhang'
  


contains


  ! ===============================================================
  

  !in: landuse class, particle size, surface temperature, temperature
  !    density, snowcover, rain, Ra (=ftop), ustar
  !out: vd, vd+lu, vs
  !input: particle size (diameter) in m.

  subroutine LE_DryDepos_Aero_Zhang_vd( vd, Rs, vs1, partsize, slipcor, &
                                       nwet, tsurf, dens1, viscos1, &
                                       luc, lai, ftop_lu, ustar_lu, status )

    use Binas        , only : grav, pi
    use dims         , only : nz
    use LE_Landuse_Data, only : ilu_water_sea, ilu_water_inland, ilu_other, ilu_desert, ilu_ice, ilu_arable
    use LE_Landuse_Data, only : alfa_lu, gamma_lu, A_lu
   
    ! --- in/out ---------------------------------
    
    real, intent(out)        ::  vd          ! deposition velocity (m/s)
    real, intent(out)        ::  Rs          ! Sedimentaion resistance (s/m)
    real, intent(in)         ::  vs1         ! sedimentation velocity in lowest layer
    real, intent(in)         ::  partsize    ! particle diameter (m)
    real, intent(in)         ::  slipcor     ! slip correction factor
    integer, intent(in)      ::  nwet        ! 1=rain, 9=snowcover
    real, intent(in)         ::  tsurf       ! surface temperature (K)
    real, intent(in)         ::  dens1       ! air density (kg/m3) in lowest layer
    real, intent(in)         ::  viscos1     ! air viscosity in lowest layer
    integer, intent(in)      ::  luc         ! depac landuse code
    real, intent(in)         ::  lai         ! leaf area index
    real, intent(in)         ::  ftop_lu     ! atmospheric resistnace Ra
    real, intent(in)         ::  ustar_lu    ! friction velocity u*    
    integer, intent(out)     ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_DryDepos_Aero_Zhang_vd'
    
    real, parameter   ::  beta     = 2.0
    real, parameter   ::  epsilon0 = 3.0
    real, parameter   ::  kb       = 1.38066e-23

    ! --- local ----------------------------------

    real              ::  kinvisc, diff_part
    real              ::  schmidt,stokes, Ebrown, Eimpac, Einterc, Reffic

    ! --- begin ----------------------------------

    ! kinetic viscosity & diffusivity
    kinvisc = viscos1 / dens1 ! only needed at surface
    diff_part = kb * tsurf * slipcor / (3*pi*viscos1*partsize)

    ! Schmidt number
    schmidt = kinvisc / diff_part

    ! calculate collection efficiencie E
    Ebrown = Schmidt**(-gamma_lu(luc)) ! Brownian diffusion
    ! determine Stokes number, interception efficiency 
    ! and sticking efficiency R (1 = no rebound)
    
    ! Icy surface  
    if ( luc == ilu_ice .or. nwet.eq.10 .or. (luc==ilu_water_inland .and. tsurf <= 269.15) ) then
      stokes=vs1*ustar_lu**2/(grav*kinvisc)
      Einterc=0.0
      Reffic=exp(-Stokes**0.5)

    ! Wet surface
    else if ( nwet.eq.9 .or. luc == ilu_water_sea .or. luc == ilu_water_inland ) then !snow cover, sneeuw uit meteo nog meenemen! or.snow....
      stokes=vs1*ustar_lu**2/(grav*kinvisc)
      Einterc=0.0
      Reffic=1.0
      
    ! Dry surface without any vegetation (tundra,desert or fallow arable land )
    else if ( luc == ilu_other .or. luc == ilu_desert .or. ( luc == ilu_arable .and. lai < 0.01 ) ) then 
      stokes=vs1*ustar_lu**2/(grav*kinvisc)
      Einterc=0.0
      Reffic=exp(-Stokes**0.5)
    else
      stokes=vs1*ustar_lu/(grav*A_lu(luc)*1.e-3)
      Einterc=0.5*(partsize/(A_lu(luc)*1e-3))**2
      Reffic=exp(-Stokes**0.5)
    end if
    ! when surface is wet all particles do not rebound:
    if(nwet.eq.1) Reffic = 1.0
    
    ! Snow surface or icy surface
    if ( nwet.eq.9 .or. nwet.eq.10 ) then
      Eimpac = ( stokes / (50.0 + stokes ) ) **beta
    else
      ! determine impaction efficiency:
      Eimpac = ( stokes / (alfa_lu(luc)+stokes) )**beta
    end if

    ! sedimentation resistance:
    Rs = 1.0 / ( epsilon0 * ustar_lu * (Ebrown+Eimpac+Einterc) * Reffic )

    ! deposition velocity according to Seinfeld and Pandis (= SP, 2006; eq 19.7):
    !
    !            1
    !    vd = ------------------ + vs
    !         Ra + Rs + Ra*Rs*vs
    !
    !    where: Rs = Rb (in SP, 2006)
    !
    !
    vd = 1.0 / ( ftop_lu + Rs + ftop_lu*Rs*vs1) + vs1

    ! ok
    status = 0

  end subroutine LE_DryDepos_Aero_Zhang_vd



end module LE_DryDepos_Aero_Zhang

