!###############################################################################
!
! Landuse_LSM  -  Make land sea mask
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!
#include "le.inc"
!
!###############################################################################


module LE_LandUse_LSM

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  Landuse_LSM_Init, Landuse_LSM_Done


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'Landuse_LSM'


  ! --- var --------------------------------------

contains


  ! ========================================================================


  subroutine Landuse_LSM_Init( rcF, status )

    use GO, only : TrcFile, ReadRc
    use GO, only : GoVarValue
    use Grid        , only : TllGridInfo, Init, Done, IndexFractions
    use LE_Grid     , only : lli
    use dims, only : nx, ny
    use LE_Landuse_Data, only : waterbodies, nwaterbody_type
    use LE_LandUse_File_nc, only : T_Landuse_File_nc, Init, Done, Get
    
    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_LandUse_LSM_Init'
  
    ! --- local ----------------------------------
    
    character(len=512)          ::  query
    character(len=256)          ::  fname
    ! read from world water body file:
    type(T_Landuse_File_nc)     ::  wwb
    type(TllGridInfo)       ::  lliX

    integer, allocatable    ::  waterbodiesX(:,:)
    real                    ::  west_bound, east_bound, dlon
    real                    ::  south_bound, north_bound, dlat
    integer                 ::  nlon, nlat
    
    ! --- begin ----------------------------------
    
    call ReadRc( RcF, 'my.landsea.waterbody.file', query, status )
    IF_NOTOK_RETURN(status=1)
    
    
    ! name of waterbody file:
    call goVarValue( trim(query), ';', 'file', '=', fname, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open database:
    call Init( wwb, trim(fname), status )
    IF_NOTOK_RETURN(status=1)

    ! extract grid parameters:
    call Get( wwb, status, &
                west_bound=west_bound, east_bound=east_bound, nlon=nlon, &
                south_bound=south_bound, north_bound=north_bound, nlat=nlat )
    IF_NOTOK_STOP

    ! spacing:
    dlon = (  east_bound -  west_bound ) / nlon
    dlat = ( north_bound - south_bound ) / nlat

    ! setup storage:
    allocate( waterbodiesX(nlon,nlat) )

    ! read:
    call Get( wwb, status, wwb=waterbodiesX )
    IF_NOTOK_STOP

    ! close:
    call Done( wwb, status )
    IF_NOTOK_STOP
  
    ! info ...
    write (gol,'("  convert to LE grid ...")'); call goPr
    
    ! setup input grid:
    call Init( lliX, west_bound+0.5*dlon, dlon, nlon, &
                    south_bound+0.5*dlat, dlat, nlat, status )
    IF_NOTOK_STOP

    ! convert from 2D field with indices to 3D field with fraction per cell for each index:
    call IndexFractions( lliX, waterbodiesX, lli, 0, nwaterbody_type-1, waterbodies, status )
    IF_NOTOK_RETURN(status=1)

    ! done:
    call Done( lliX, status )
    IF_NOTOK_STOP

    ! clear:
    if ( allocated( waterbodiesX   ) ) deallocate( waterbodiesX    )
        
    ! ok
    status = 0

  end subroutine Landuse_LSM_Init


  ! ***


  subroutine Landuse_LSM_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_LandUse_LSM_Done'

    ! --- begin ----------------------------------

    ! clear:
    !deallocate( seafraction )

    ! ok
    status = 0

  end subroutine Landuse_LSM_Done


end module LE_LandUse_LSM

