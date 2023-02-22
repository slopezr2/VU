module JAQL_drydeposition

  implicit none


  ! --- in/out --------------------------------

  private

  public  ::  Free_Path_Length
  public  ::  Slip_Correction_Factor
  public  ::  Sedimentation_Velocity
  public  ::  Viscosity


contains


  ! ====================================================================


  ! free path length
  ! should be improved, 3D field

  elemental function Free_Path_Length( tsurf, psurf ) result ( freepathlen )

    ! --- in/out ---------------------------------

    real, intent(in)        ::  tsurf      ! surface temperature [K]
    real, intent(in)        ::  psurf      ! surface pressure [Pa]
    real                    ::  freepathlen

    ! --- begin ----------------------------------

    ! free path length:
    freepathlen = 2.332e-5 * tsurf/psurf

  end function Free_Path_Length


  ! ***


  ! compute sedimentaiton velocity

  elemental function Viscosity( temp ) result ( visc )

    ! --- in/out ---------------------------------

    real, intent(in)        ::  temp       ! layer temperature
    real                    ::  visc

    ! --- begin ----------------------------------

    ! viscosity:
    visc = 1.496e-6 * temp**1.5 / (temp+120.0)

  end function Viscosity


  ! ***


  ! compute slip correction factor
  ! E.g. Seinfeld & Pandis 2006, Eq. 9.34 p.40

  elemental function Slip_Correction_Factor( freepathlen, partsize ) result ( slipcor )

    ! --- in/out ---------------------------------

    real, intent(in)        ::  freepathlen
    real, intent(in)        ::  partsize   ! diameter in m
    real                    ::  slipcor

    ! --- begin ----------------------------------

    ! slip correction factor
    slipcor = 1.0 + (2*freepathlen/partsize)*(1.257+0.4*exp(-0.55*partsize/freepathlen))

    !!>>> compatibility with v1.9.000 --------------+
    !!                                              v
    !slipcor = 1.0 + (2*freepathlen/partsize)*(1.257*0.4*exp(-0.55*partsize/freepathlen))
    !!<<<

  end function Slip_Correction_Factor


  ! *


  ! compute sedimentaiton velocity

  elemental function Sedimentation_Velocity( rhopart, partsize, slipcor, visc ) result ( vs )

    use Binas, only : grav

    ! --- in/out ---------------------------------

    real, intent(in)        ::  rhopart    ! kg/m3
    real, intent(in)        ::  partsize   ! m
    real, intent(in)        ::  slipcor    ! m
    real, intent(in)        ::  visc       ! viscosity
    real                    ::  vs

    ! --- begin ----------------------------------

    !viscosity & sedimentation velocity
    vs = rhopart * (partsize**2) * grav * slipcor / (18*visc)

  end function Sedimentation_Velocity


end module JAQL_drydeposition
