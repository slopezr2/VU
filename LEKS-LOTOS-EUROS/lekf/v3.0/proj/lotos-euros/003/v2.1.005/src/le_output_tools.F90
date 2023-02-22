!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Output_Tools

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_Create, NF90_Close
  use NetCDF, only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
  use NetCDF, only : NF90_Put_Var, NF90_Put_Att
  use NetCDF, only : NF90_NOCLOBBER, NF90_CLOBBER, NF90_GLOBAL, NF90_UNLIMITED
  use NetCDF, only : NF90_REAL, NF90_INT, NF90_DOUBLE
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  implicit none

  ! --- in/out -----------------------------

  private
  
  public :: LE_Output_Define_Dims_Lon_Lat
  public :: LE_Output_Define_Dims_Lev
  public :: LE_Output_Define_Dims_Time
 
  public :: LE_Output_Define_Vars_Lon_Lat
  public :: LE_Output_Define_Vars_Lev
  public :: LE_Output_Define_Vars_Time
 
  public :: LE_Output_Put_Var_Lon_Lat
  public :: LE_Output_Put_Var_Lev
  public :: LE_Output_Put_Var_Time
  
  public :: LE_Output_Put_Var_Domains
  public :: LE_Output_Put_Var_4D
  
  !public  ::  itim_init_output_conc, itim_init_output_conc_com, itim_init_output_conc_rc
  
  
  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Tools'

    
  ! --- var --------------------------------
  
  !! timers:
  !integer   ::  itim_init_output_conc, itim_init_output_conc_com, itim_init_output_conc_rc



contains 


  !============================================================================
    

  subroutine LE_Output_Define_Dims_Lon_Lat(ncid,dimid_lon,dimid_lat, &
                                           nlon, nlat, &
                                           convention, status)
      
    ! --- in/out ----------------------
    integer, intent(inout)        ::  ncid, dimid_lon, dimid_lat
    integer, intent(in)           ::  nlon, nlat
    character(len=*), intent(in)  ::  convention
    integer, intent(out)          ::  status
 
    ! --- const ------------------------------  
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Define_Dims_Lon_Lat'
    
    ! --- begin ------------------------------
    
    ! Default dimension names for comparison with paraview program 'longitude', 'latitude'
    ! for racmo 'X', 'Y' needed
    select case ( convention )
      case( 'racmo' )
        status = NF90_Def_Dim( ncid, 'X', nlon, dimid_lon )
        IF_NF90_NOTOK_RETURN(status=1)

        status = NF90_Def_Dim( ncid, 'Y', nlat, dimid_lat )
        IF_NF90_NOTOK_RETURN(status=1)
      
      case default
        status = NF90_Def_Dim( ncid, 'longitude', nlon, dimid_lon )
        IF_NF90_NOTOK_RETURN(status=1)

        status = NF90_Def_Dim( ncid, 'latitude', nlat, dimid_lat )
        IF_NF90_NOTOK_RETURN(status=1)
    end select
    ! ok
    status = 0

  end subroutine LE_Output_Define_Dims_Lon_Lat

! -----------------------------------------------------------------------------  
  
  subroutine LE_Output_Define_Dims_Lev(ncid,dimid_lev, nlev, convention, status, dimid_hlev)
  
    ! --- in/out ----------------------
    integer, intent(inout)            ::  ncid, dimid_lev
    integer, intent(in)               ::  nlev
    character(len=*), intent(in)      ::  convention
    integer, intent(out)              ::  status
    integer, intent(inout), optional  ::  dimid_hlev
 
    ! --- const ------------------------------  
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Define_Dims_Lev'
    
    ! --- begin ------------------------------
    
    ! Default dimension names for comparison with paraview program 'level'
    ! for racmo 'Z' needed
    select case ( convention )
      case ( 'racmo' )
        status = NF90_Def_Dim( ncid, 'Z', nlev, dimid_lev )
        IF_NF90_NOTOK_RETURN(status=1)
      case default
        status = NF90_Def_Dim( ncid, 'level', nlev, dimid_lev )
        IF_NF90_NOTOK_RETURN(status=1)
        ! half level dimension?
        if ( present(dimid_hlev) .and. nlev > 1 ) then
          status = NF90_Def_Dim( ncid, 'hlevel', nlev+1, dimid_hlev )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
    end select
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Define_Dims_Lev  

! -----------------------------------------------------------------------------  
  
  subroutine LE_Output_Define_Dims_Time(ncid, dimid_time, status ) 
  
    ! --- in/out ----------------------
    integer, intent(inout)        ::  ncid, dimid_time
    integer, intent(out)          ::  status
 
    ! --- const ------------------------------  
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Define_Dims_Time'

    ! no special default for racmo  
    status = NF90_Def_Dim( ncid, 'time', NF90_UNLIMITED, dimid_time )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Define_Dims_Time  
  
! -----------------------------------------------------------------------------  
  
  subroutine LE_Output_Define_Vars_Lon_Lat(ncid, varid_lon, varid_lat, &
                                   dimid_lon, dimid_lat, convention, &
                                   status)
  
    ! --- in/out -------------------------------
    integer, intent(inout)      ::  ncid, varid_lon, varid_lat
    integer, intent(in)         ::  dimid_lon, dimid_lat
    character(len=*), intent(in)::  convention
    integer, intent(out)        ::  status
 
    ! --- const ------------------------------  
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Define_Vars_Lon_Lat'
    
    ! --- local --------------------------------   
    integer                     :: varid

    ! --- begin --------------------------------
    
    ! for paraview default variables must have the same name as the dimensions for longitude, latitude    
    
    select case ( convention )
      case ( 'racmo' )
        ! longitude     
        status = NF90_Def_Var( ncid, 'lon', NF90_REAL, (/dimid_lon,dimid_lat/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'longitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'longitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Lon' )
        IF_NF90_NOTOK_RETURN(status=1)
        varid_lon = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'lat', NF90_REAL, (/dimid_lon,dimid_lat/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'latitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'latitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Lat' )
        IF_NF90_NOTOK_RETURN(status=1)
        varid_lat = varid
    case default
        ! longitude     
        status = NF90_Def_Var( ncid, 'longitude', NF90_REAL, dimid_lon, varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'longitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'longitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Lon' )
        IF_NF90_NOTOK_RETURN(status=1)
        varid_lon = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'latitude', NF90_REAL, dimid_lat, varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'latitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'latitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Lat' )
        IF_NF90_NOTOK_RETURN(status=1)
        varid_lat = varid
    end select
    
    ! ok
    status = 0    

  end subroutine LE_Output_Define_Vars_Lon_Lat
  
! -----------------------------------------------------------------------------    
  
  subroutine LE_Output_Define_Vars_Lev(ncid, varid_lev, dimid_lev, &
                                       levtype, convention, status)
  
    ! --- in/out -------------------------------
    integer, intent(inout)      ::  ncid, varid_lev
    integer, intent(in)         ::  dimid_lev
    character(len=*), intent(in)::  levtype
    character(len=*), intent(in)::  convention   
    integer, intent(out)        ::  status
 
    ! --- const ------------------------------  
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Define_Vars_Lev'
    
    ! --- local --------------------------------   
    integer                     :: varid

    ! --- begin --------------------------------
    
    ! for paraview default variables must have the same name as the dimensions for lon, lat and level
    
    select case ( convention )
      
      case( 'racmo' )
        ! level
        status = NF90_Def_Var( ncid, 'lev', NF90_REAL, dimid_lev, varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'layer' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Height' )
        IF_NF90_NOTOK_RETURN(status=1)
        varid_lev = varid
      case default
        ! which type of output levels ?
        select case ( trim(levtype) )
          case ( 'levels', 'half_levels', 'input_levels', 'input_halflevels' )
            status = NF90_Def_Var( ncid, 'level', NF90_INT, dimid_lev, varid )
            IF_NF90_NOTOK_RETURN(status=1)
            status = NF90_Put_Att( ncid, varid, 'long_name', 'layer' )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 'heights' )
            status = NF90_Def_Var( ncid, 'level', NF90_REAL, dimid_lev, varid )
            IF_NF90_NOTOK_RETURN(status=1)
            status = NF90_Put_Att( ncid, varid, 'long_name', 'layer height relative to orography' )
            IF_NF90_NOTOK_RETURN(status=1)
            status = nf90_put_att( ncid, varid, 'standard_name', 'height' )
            IF_NF90_NOTOK_RETURN(status=1)
            status = NF90_Put_Att( ncid, varid, 'units', 'm' )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 'elevations' )
            status = NF90_Def_Var( ncid, 'level', NF90_REAL, dimid_lev, varid )
            IF_NF90_NOTOK_RETURN(status=1)
            status = NF90_Put_Att( ncid, varid, 'long_name', 'layer height relative to sea level' )
            IF_NF90_NOTOK_RETURN(status=1)
            status = nf90_put_att( ncid, varid, 'standard_name', 'altitude' )
            IF_NF90_NOTOK_RETURN(status=1)
            status = NF90_Put_Att( ncid, varid, 'units', 'm' )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported level type : ",a)') trim(levtype); call goErr
            TRACEBACK; status=1; return
        end select
        ! other:
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Height' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! fill:
        varid_lev = varid
    end select
    
    ! ok 
    status = 0
    
  end subroutine LE_Output_Define_Vars_Lev
  
! -----------------------------------------------------------------------------    
  
  subroutine LE_Output_Define_Vars_Time(ncid, varid_time, varid_time_dtg, &
                                        dimid_time, &
                                        convention, t0, status)
  
    ! --- in/out -------------------------------
    integer, intent(inout)      ::  ncid, varid_time, varid_time_dtg
    integer, intent(in)         ::  dimid_time
    character(len=*), intent(in)::  convention   
    type(TDate), intent(in)     ::  t0     
    integer, intent(out)        ::  status
 
    ! --- const ------------------------------  
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Define_Vars_Time'
    
    ! --- local --------------------------------   
    integer                     :: varid
    character(len=128)          :: time_since_name
    
    ! --- begin --------------------------------
       
    select case ( convention )
      case( 'racmo' ) ! racmo also has dtg output
        ! dtg
        status = NF90_Def_Var( ncid, 'dtg', NF90_INT, dimid_time, varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'date-time' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'datum-time' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'yyyymmddhh' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'calendar', 'standard' )
        IF_NF90_NOTOK_RETURN(status=1)
        varid_time_dtg = varid
    end select
        
    ! time
    status = NF90_Def_Var( ncid, 'time', NF90_DOUBLE, dimid_time, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid, 'standard_name', 'time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid, 'long_name', 'time' )
    IF_NF90_NOTOK_RETURN(status=1)
    write (time_since_name,'("seconds since ",i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2," UTC")') &
                    t0%year, t0%month, t0%day, t0%hour,t0%min, t0%sec  
    status = NF90_Put_Att( ncid, varid, 'units', trim(time_since_name) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid, 'calendar', 'standard' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Time' )
    IF_NF90_NOTOK_RETURN(status=1)
    varid_time = varid
    
!    ! extra variable time_day
!    status = NF90_Def_Var( leo%ncid, 'time_day', NF90_REAL, leo%dimid_time, varid )
!    IF_NF90_NOTOK_RETURN(status=1)
!    status = NF90_Put_Att( leo%ncid, varid, 'standard_name', 'time' )
!    IF_NF90_NOTOK_RETURN(status=1)
!    status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'time' )
!    IF_NF90_NOTOK_RETURN(status=1)
!    status = NF90_Put_Att( leo%ncid, varid, 'units', 'yyyymmdd+dayfraction' )
!    IF_NF90_NOTOK_RETURN(status=1)
!    status = NF90_Put_Att( leo%ncid, varid, 'calendar', 'standard' )
!    IF_NF90_NOTOK_RETURN(status=1)
!    leo%varid_time_day = varid

    !ok
    status = 0
    
  end subroutine LE_Output_Define_Vars_Time
  
! -----------------------------------------------------------------------------    
  
  subroutine LE_Output_Put_Var_Lon_Lat( ncid, varid_lon, varid_lat,  &
                                         nlon, nlat, lons, lats, &
                                         convention, status )
                                       
    use Dims, only    : nx, ny
    
    ! --- in/out -------------------------------
    integer, intent(inout)      ::  ncid
    integer, intent(in)         ::  varid_lon, varid_lat
    integer, intent(in)         ::  nlon, nlat
    real, intent(in)            ::  lons(nlon), lats(nlat)
    character(len=*)            ::  convention
    integer, intent(out)        ::  status
    
    ! --- const ------------------------------  

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Put_Var_Lon_Lat'

    ! --- local --------------------------------

    integer                     ::  i, j

    ! --- begin --------------------------------
    
    select case ( convention )
      case ( 'racmo' )
        do j = 1, nlat
          status = NF90_Put_Var( ncid, varid_lon, lons, &
                                 start=(/1,j/), count=(/nlon,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end do
        do i = 1, nlon
          status = NF90_Put_Var( ncid, varid_lat, lats, &
                                 start=(/i,1/), count=(/1,nlat/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end do
      
      case default
        status = NF90_Put_Var( ncid, varid_lon, lons )
        IF_NF90_NOTOK_RETURN(status=1)
        
        status = NF90_Put_Var( ncid, varid_lat, lats )
        IF_NF90_NOTOK_RETURN(status=1)
    end select
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Put_Var_Lon_Lat
  
! -----------------------------------------------------------------------------    
  
  subroutine LE_Output_Put_Var_Lev(ncid, varid_lev, nlev, levtype, ilev, heights, status)
                                       
    
    ! --- in/out -------------------------------
    integer, intent(inout)      ::  ncid
    integer, intent(in)         ::  varid_lev
    integer, intent(in)         ::  nlev
    character(len=*)            ::  levtype
    integer, intent(in)         ::  ilev(1:nlev)
    real, intent(in)            ::  heights(1:nlev)
    integer, intent(out)        ::  status
    
    ! --- const ------------------------------  
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Put_Var_Lev'

    ! --- local --------------------------------
   
    ! --- begin --------------------------------
    
    select case ( levtype )
      case ( 'levels', 'half_levels', 'input_levels', 'input_halflevels' )
        status = NF90_Put_Var( ncid, varid_lev, real(ilev(1:nlev)) )
        IF_NF90_NOTOK_RETURN(status=1)
      case ( 'heights', 'elevations' )
        status = NF90_Put_Var( ncid, varid_lev, heights(1:nlev) )
        IF_NF90_NOTOK_RETURN(status=1)
      case default
        write (gol,'("unsupported level type : ",a)') trim(levtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Put_Var_Lev
  
! -----------------------------------------------------------------------------    
  
  subroutine LE_Output_Put_Var_Time(ncid, varid_time, varid_time_dtg, &
                                    time, time6, convention, itrec, status)
                                       
    
    ! --- in/out -------------------------------
    integer, intent(inout)      ::  ncid
    integer, intent(in)         ::  varid_time, varid_time_dtg
    integer, intent(in)         ::  time    
    integer, intent(in)         ::  time6(6)
    character(len=*), intent(in)::  convention
    integer, intent(in)         ::  itrec
    integer, intent(out)        ::  status
    
    ! --- const ------------------------------  
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Put_Var_Time'

    ! --- local --------------------------------
    integer                     :: dtg
    ! --- begin --------------------------------
    
    ! time variable
    status = NF90_Put_Var( ncid, varid_time, (/time/), start = (/itrec/) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    select case ( convention )
      case ( 'racmo' )
        !20-1-2010 time6 factor changed
        dtg = time6(1)*1000000 + time6(2)*10000 + time6(3)*100 + time6(4)
        ! write time record:
        status = NF90_Put_Var( ncid, varid_time_dtg, (/dtg/), start=(/itrec/) )
        IF_NF90_NOTOK_RETURN(status=1)
    end select
   
!    ! extra variable time_day   
!    status = NF90_Put_Var( leo%ncid, leo%varid_time_day, (/time/), start=(/leo%itrec/) )
!    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Put_Var_Time


  ! ***


  !
  ! Collect local 2D domains on root, put into netcdf file.
  ! Use ilev<0 for variable without level dimension.
  ! Use irec<0 for variable without time dimension.
  !

  subroutine LE_Output_Put_Var_Domains( ncid, varid, ilev, irec, values, status, bb_start, bb_count, ilab )
  
#ifdef with_netcdf
    use NetCDF , only : NF90_Put_Var
#endif
    use GO     , only : goc
    use LE_Grid, only : dom

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  ncid
    integer, intent(in)       ::  varid
    integer, intent(in)       ::  ilev
    integer, intent(in)       ::  irec
    real, intent(in)          ::  values(:,:)
    integer, intent(out)      ::  status
    integer, intent(in), optional :: bb_start(:)
    integer, intent(in), optional :: bb_count(:)
    integer, intent(in), optional :: ilab
    

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Put_Var_Domains'

    ! --- local ----------------------------------
    
    integer               ::  shp(2)
    real, allocatable     ::  glb_values(:,:)
    integer               ::  i1, i2, j1, j2, c1, c2
    ! --- begin ----------------------------------
    
    ! global shape:
    call dom%Get( status, glb_shp=shp )
    IF_NOTOK_RETURN(status=1)

    ! storage:
    if ( goc%root ) then
      allocate( glb_values(shp(1),shp(2)), stat=status )
      IF_NOTOK_RETURN(status=1)
    else
      allocate( glb_values(1,1), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! collect:
    call dom%GatherV( values, glb_values, status )
    IF_NOTOK_RETURN(status=1)
    
    ! Bounding box? (start and count indices)
    if ( present( bb_start ) ) then
      i1 = bb_start(1)
      i2 = i1 + bb_count(1)-1
      j1 = bb_start(2)
      j2 = j1 + bb_count(2)-1
      c1 = bb_count(1)
      c2 = bb_count(2)
    ! default full grid
    else 
      i1 = 1
      i2 = shp(1)
      j1 = 1
      j2 = shp(2)
      c1 = shp(1)
      c2 = shp(2)
    end if
    
#ifdef with_netcdf
    ! root?
    if ( goc%root ) then
    
      ! without time dimension?
      if ( irec < 0 ) then

        ! without level dimension?
        if ( ilev < 0 ) then
          ! write:
          status = NF90_Put_Var( ncid, varid, glb_values(i1:i2,j1:j2) )
          IF_NF90_NOTOK_RETURN(status=1)
        else
          ! write level:
          status = NF90_Put_Var( ncid, varid, glb_values(i1:i2,j1:j2),  &
                                 start=(/1,1,ilev/), count=(/c1,c2,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if ! level dim?

      else

        ! without level dimension?
        if ( ilev < 0 ) then
          ! With label dimension
          if ( present(ilab) ) then
            ! write:
            status = NF90_Put_Var( ncid, varid, glb_values(i1:i2,j1:j2), &
                                   start=(/1,1,ilab,irec/), count=(/c1,c2,1,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
          else 
            ! write:
            status = NF90_Put_Var( ncid, varid, glb_values(i1:i2,j1:j2), &
                                   start=(/1,1,irec/), count=(/c1,c2,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
        else
          ! With label dimension
          if ( present(ilab) ) then
            ! write level:
            status = NF90_Put_Var( ncid, varid, glb_values(i1:i2,j1:j2), &
                                   start=(/1,1,ilev,ilab,irec/), count=(/c1,c2,1,1,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
          else 
            ! write level:
            status = NF90_Put_Var( ncid, varid, glb_values(i1:i2,j1:j2), &
                                   start=(/1,1,ilev,irec/), count=(/c1,c2,1,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
        end if ! level dim?

      end if ! time dim?
      
    end if  ! root
#endif

    ! storage:
    deallocate( glb_values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Put_Var_Domains


  ! ***


  !
  ! Collect local (x,y,:,:) domains on root, put into netcdf file.
  ! Use irec<0 for variable without time dimension.
  !

  subroutine LE_Output_Put_Var_4D( ncid, varid, irec, values, status, bb_start, bb_count )
  
#ifdef with_netcdf
    use NetCDF , only : NF90_Put_Var
#endif
    use GO     , only : goc
    use LE_Grid, only : dom

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  ncid
    integer, intent(in)       ::  varid
    integer, intent(in)       ::  irec
    real, intent(in)          ::  values(:,:,:,:)
    integer, intent(out)      ::  status
    integer, intent(in), optional :: bb_start(:)
    integer, intent(in), optional :: bb_count(:)
    

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Put_Var_Domains'

    ! --- local ----------------------------------
    
    integer               ::  shp(2)
    integer               ::  nz, nv
    real, allocatable     ::  glb_values(:,:,:,:)  ! (nx,ny,nz,nv)
    integer               ::  i1, i2, j1, j2, c1, c2
    ! --- begin ----------------------------------
    
    ! global shape:
    call dom%Get( status, glb_shp=shp )
    IF_NOTOK_RETURN(status=1)
    ! other dims:
    nz = size(values,3)
    nv = size(values,4)

    ! storage:
    if ( goc%root ) then
      allocate( glb_values(shp(1),shp(2),nz,nv), stat=status )
      IF_NOTOK_RETURN(status=1)
    else
      allocate( glb_values(1,1,1,1), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! collect:
    call dom%GatherV( values, glb_values, status )
    IF_NOTOK_RETURN(status=1)
    
    ! Bounding box? (start and count indices)
    if ( present( bb_start ) ) then
      i1 = bb_start(1)
      i2 = i1 + bb_count(1)-1
      j1 = bb_start(2)
      j2 = j1 + bb_count(2)-1
      c1 = bb_count(1)
      c2 = bb_count(2)
    ! default full grid
    else 
      i1 = 1
      i2 = shp(1)
      j1 = 1
      j2 = shp(2)
      c1 = shp(1)
      c2 = shp(2)
    end if
    
#ifdef with_netcdf
    ! root?
    if ( goc%root ) then
    
      ! without time dimension?
      if ( irec < 0 ) then

        ! write:
        status = NF90_Put_Var( ncid, varid, glb_values(i1:i2,j1:j2,:,:) )
        IF_NF90_NOTOK_RETURN(status=1)

      else

        ! write slab:
        status = NF90_Put_Var( ncid, varid, glb_values(i1:i2,j1:j2,:,:), &
                               start=(/1,1,1,1,irec/), count=(/c1,c2,nz,nv,1/) )
        IF_NF90_NOTOK_RETURN(status=1)

      end if ! time dim?
      
    end if  ! root
#endif

    ! storage:
    deallocate( glb_values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Put_Var_4D

end module LE_Output_Tools
