!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
!   C3PO_File_Ugg  : access files with data defined on univerdsal-gaussian-grid
!
! EXAMPLE
!
!    use GO     , only : TDate, NewDate
!    use C3PO   , only : T_File_Ugg
!    use C3PO   , only : T_Grid_Ugg
!    use C3PO   , only : T_Levs_Hyb
!
!    ! variables:
!    type(T_File_Ugg)      ::  infile
!    type(T_Grid_Ugg)      ::  ugg
!    type(T_Levs_Hyb)      ::  levs
!    type(TDate)           ::  t
!    integer               ::  irec
!    real, allocatable     ::  values(:,:)  ! (nlon,nlat,nlev)
!    integer               ::  status
!
!    ! open file:
!    call infile%Open( 'temperature.nc', status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! init grid definition from file:
!    call infile%Get_Grid( ugg, status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! get level definition
!    call infile%Get_Levs( levs, status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! target time:
!    t = NewDate( 2000, 01, 02, 00, 00 )
!    
!    ! index of time record:
!    call infile%Inq_TimeRecord( 'time', t, irec, status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! storage:
!    allocate( values_in(ugg%nlon,ugg%nlat,levs%nlev) )
!
!    ! read data slab:
!    call infile%Get_Var( trim(nametype), trim(name), values_in, units_in, status, &
!                              start=(/1,1,irec/), count=(/ugg%nlon,ugg%nlat,levs%nlev,1/) )
!    IF_NOT_OK_RETURN(status=1)
!      
!    ! clear input field:
!    deallocate( values_in )
!    
!    ! clear input grid definition:
!    call ugg%Done( status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! clear input level definition:
!    call levs%Done( status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! close input file:
!    call infile%Close( status )
!    IF_NOT_OK_RETURN(status=1)
!
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
#include "c3po.inc"
!
!######################################################################

module C3PO_File_Ugg

  use GO    , only : gol, goPr, goErr
  use NetCDF, only : NF90_NOERR, NF90_StrError
  
  use C3PO_File_Nc, only : T_File_Nc, maxdim

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_File_Ugg
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_File_Ugg'

  ! ugg parameters:
  character(len=*), parameter  ::  dimname_nv   = 'nv'
  character(len=*), parameter  ::  varext_bnds  = '_bnds'
  integer, parameter           ::  nv = 2
  

  ! --- types ----------------------------------------
  
  type, extends(T_File_Nc) :: T_File_Ugg
    ! ugg dims:
    integer                 ::  dimid_lat
    integer                 ::  dimid_lon
    integer                 ::  dimid_nv
    ! ugg variables:
    integer                 ::  varid_longitude
    integer                 ::  varid_latitude
    integer                 ::  varid_longitude_bnds
    integer                 ::  varid_latitude_bnds
  contains
    procedure   ::  Get_Grid                => File_Ugg_Get_Grid
!    procedure   ::  Def_Grid                => File_Ugg_Def_Grid
!    procedure   ::  Put_Grid                => File_Ugg_Put_Grid
!    procedure   ::  Copy_Variable_Selection => File_Ugg_Copy_Variable_Selection
  end type T_File_Ugg

  
  

contains


  ! ********************************************************************
  ! ***
  ! *** ugg file
  ! ***
  ! ********************************************************************

  !
  ! Read grid definition from file.
  ! Provide variable id of sample variable,
  ! from this the dimensions and corresponding coordinates are read.
  !
  ! Optional arguments:
  !
  !   subset  : (in) integer array (/i1,i2,j1,j2/) to select only subdomain;
  !             use negative values for all
  !
  !   ugg_to  : (in) target grid after interpolation ; 
  !             this limits the defined 'ugg' to only the necessary domain ;
  !             should be combined with arguments:
  !   start_ind(2), count_ind(2)  :  (out) 2D start/count indices for reading 
  !                                  slab of data from netcdf file
  !

  subroutine File_Ugg_Get_Grid( self, vid, ugg, status, &
                                  ugg_to, start_ind, count_ind, &
                                  subset )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att, NF90_ENOTATT
        
    use Num          , only : Interp_Lin, EXTRAPOL_LINEAR
    use Num          , only : Interval, Interval_Modulo
    use C3PO_Grid_ugg, only : T_Grid_ugg
    use GO           , only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_File_Ugg), intent(in)           ::  self
    integer, intent(in)                     ::  vid
    class(T_Grid_ugg), intent(out)          ::  ugg
    integer, intent(out)                    ::  status
    class(T_Grid_ugg), intent(in), optional ::  ugg_to
    integer, intent(out), optional          ::  start_ind(2)
    integer, intent(out), optional          ::  count_ind(2)
    integer, intent(in), optional           ::  subset(4)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/File_Ugg_Get_Grid'
    
    ! --- local ----------------------------------

    integer                         ::  ndims
    integer, allocatable            ::  dimids(:)
    integer                         ::  idim
    integer                         ::  dimid
    character(len=32), allocatable  ::  dimnames(:)
    character(len=32)               ::  dimname
    character(len=32)               ::  varname
    integer                         ::  varid
    integer                         ::  varid_lon, varid_lat
    character(len=16)               ::  dim_case

    character(len=128)              ::  coordinates
    character(len=32)               ::  coorname
    character(len=64)               ::  standard_name
    real, allocatable               ::  xlons(:,:)
    real, allocatable               ::  xlats(:,:)

    integer                         ::  nlon, nlat, nv
    real, allocatable               ::  lats(:)
    real, allocatable               ::  lats_bnds(:,:)
    real, allocatable               ::  lons(:)
    real, allocatable               ::  lons_bnds(:,:)
    real, allocatable               ::  lons2(:)
    real, allocatable               ::  lons2_bnds(:,:)
    real, allocatable               ::  ii(:)
    integer                         ::  i, j
    
    real                            ::  bbox(4)
    integer                         ::  idir, jdir
    integer                         ::  i1, i2, j1, j2
    
    ! --- begin ----------------------------------
    
    ! empty dimensions
    dim_case = 'None'

    ! get  number of dimensions:
    status = NF90_Inquire_Variable( self%ncid, vid, ndims=ndims )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( dimids(ndims), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( dimnames(ndims), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! get dimension id's:
    status = NF90_Inquire_Variable( self%ncid, vid, dimids=dimids )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    !
    ! check for coordinates attribute ;
    ! this maps from rotated to unrotated locations:
    !
    !   ! variable along rotated coordinates (1D),
    !   ! attribute describes the actual coordinates:
    !   float HGT(Time, south_north, west_east) ;
    !     HGT:coordinates = "XLONG XLAT" ;
    !
    !   ! actual coordinates (2D):
    !   float XLAT(Time, south_north, west_east) ;
    !     XLAT:units = "degree_north" ;
    !   float XLONG(Time, south_north, west_east) ;
    !     XLONG:units = "degree_east" ;

    ! select first dimension:
    dimid = dimids(1)
    ! get dimension name and length:
    status = NF90_Inquire_Dimension( self%ncid, dimid, name=dimnames(1), len=nlon )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! select second dimension:
    dimid = dimids(2)
    ! get dimension name and length:
    status = NF90_Inquire_Dimension( self%ncid, dimid, name=dimnames(2), len=nlat )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! check if grid variables are defined in coordinate attribute
    !
    call self%Inq_VarCoordinates( vid, coordinates, status )
    !~ coordinates attribute defined:
    if ( status == 0 ) then
      
      ! loop over coordinates in attribute:
      do
        ! leave if empty:
        if ( len_trim(coordinates) == 0 ) exit
        ! extract next:
        call goReadFromLine( coordinates, coorname, status, sep=' ' )
        IF_NF90_NOT_OK_RETURN(status=1)
        
        ! get variable id for this coordinate:
        status = NF90_INQ_VarID( self%ncid, coorname, varid )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! try to get standard name ...
        status = NF90_Get_Att( self%ncid, varid, 'standard_name', standard_name )
        !~ attribute not found ?
        if ( status == NF90_ENOTATT ) then
          ! try to guess standard_name based on coordinate name:
          select case ( trim(coorname) )
            !~ longitudes:
            case ( 'lon', 'longitude', 'XLONG', 'XLONG_U', 'XLONG_V' )
              ! set standard name:
              standard_name = 'longitude'
            !~ latitudes:
            case ( 'lat', 'latitude', 'XLAT', 'XLAT_U', 'XLAT_V' )
              ! set standard name:
              standard_name = 'latitude'
            !~ other:
            case ( 'time', 'XTIME' )
              ! to be skipped ...
              standard_name = 'None'
            !~
            case default
              write (gol,'("unsupported coordinate name `",a,"`")') trim(coorname); call goErr
              TRACEBACK; status=1; return
          end select
        !~ other error?
        else if ( status /= NF90_NOERR ) then
          ! some error ...
          gol=NF90_StrError(status); call goErr
          TRACEBACK; status=1; return
        end if

        ! switch:
        select case ( trim(standard_name) )

          !~ longitudes:
          case ( 'longitude' )
            ! Check dimensions
            status = NF90_Inquire_Variable( self%ncid, varid, ndims=ndims )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! define which dimension case for grid definition
            if ( ndims == 1 ) then
              dim_case = '1D'
            else if ( ndims == 2 ) then
              dim_case = '2D'
            else if ( ndims == 3 ) then
              dim_case = '3D'
            else
              write (gol,'("unsupported number of dimensions in longitude coordinate `",a,"`")') trim(coorname); call goErr
              TRACEBACK; status=1; return
            end if
            ! save variable id for grid definitions
            varid_lon = varid

          !~ latitudes:
          case ( 'latitude' )
            ! save variable id for grid definitions
            varid_lat = varid

          !~ other
          case ( 'None', 'height' )
            ! no lon or lat variable (level, time?); do nothing ..
          !~
          case default
            write (gol,'("unsupported standard name `",a,"` for coordinate `",a,"`")') &
                    trim(standard_name), trim(coorname); call goErr
            TRACEBACK; status=1; return
        end select

      end do ! coordinates
      
    !~ no coordinates attribute defined:
    else if ( status /= 0 ) then

      ! check whether coordinate variable name is identical to dimension name
      varname = dimnames(1)
      status = NF90_INQ_VarID( self%ncid, trim(varname), varid )
      if ( status /= NF90_NOERR ) then
        write (gol,'("could not find coordinate variable for dimension `",a,"`")') trim(varname); call goErr
        write (gol,'("  file: ",a)') trim(self%filename); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! get number of dimensions:
      status = NF90_Inquire_Variable( self%ncid, varid, ndims=ndims )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! dimension case for grid definition (should be 1, because it is a coordinate axis )
      if ( ndims == 1 ) then
        dim_case = '1D'
      else 
        write( gol, '("Dimension of supposed coordinate variable is not 1, but: ", i0)') ndims; call GoErr
        TRACEBACK; status=1; return
      end if

      ! save variable id for grid definitions:
      varid_lon = varid

      ! get variable id for other dimension:
      varname = dimnames(2)
      status = NF90_INQ_VarID( self%ncid, trim(varname), varid )
      IF_NOT_OK_RETURN(status=1)
      ! save variable id:
      varid_lat = varid
    
    end if

    ! Select case for dimension of grid variables
    select case ( trim(dim_case) )
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 1-dimensional
      case ( '1D' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! number of boundary values:
        nv = 2

        ! ~ longitudes

        ! variable with longitudes:
        allocate( lons(nlon), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! read:
        status = NF90_Get_Var( self%ncid, varid_lon, lons )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! some files have lons (...,358,359,0,1,2,...)
        if ( lons(1) > lons(nlon) ) then
          ! shift:
          do i = 1, nlon
            if ( lons(i) > lons(nlon) ) then
              lons(i) = lons(i) - 360.0
            end if
          end do
        end if

        ! range:
        allocate( ii(nlon), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill:
        do i = 1, nlon
          ii(i) = i
        end do
        ! storage:
        allocate( lons_bnds(nv,nlon), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! interpolate boundaries:
        call Interp_Lin( ii, lons, ii-0.5, lons_bnds(1,:), status, extrapol=EXTRAPOL_LINEAR )
        IF_NOT_OK_RETURN(status=1)
        call Interp_Lin( ii, lons, ii+0.5, lons_bnds(2,:), status, extrapol=EXTRAPOL_LINEAR )
        IF_NOT_OK_RETURN(status=1)
        ! clear:
        deallocate( ii, stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ latitudes

        ! storage:
        allocate( lats(nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! read:
        status = NF90_Get_Var( self%ncid, varid_lat, lats )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! range:
        allocate( ii(nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill:
        do i = 1, nlat
          ii(i) = i
        end do
        ! storage:
        allocate( lats_bnds(nv,nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! interpolat boundaries:
        call Interp_Lin( ii, lats, ii-0.5, lats_bnds(1,:), status, extrapol=EXTRAPOL_LINEAR )
        IF_NOT_OK_RETURN(status=1)
        call Interp_Lin( ii, lats, ii+0.5, lats_bnds(2,:), status, extrapol=EXTRAPOL_LINEAR )
        IF_NOT_OK_RETURN(status=1)
        ! clear:
        deallocate( ii, stat=status )
        IF_NOT_OK_RETURN(status=1)
        
        ! ~ start and count indices

        ! by default define grid on full input file
        i1 = 1
        i2 = nlon
        j1 = 1
        j2 = nlat

        ! subset defined?
        if ( present(subset) ) then
          if ( subset(1) > 0 ) i1 = subset(1)
          if ( subset(2) > 0 ) i2 = subset(2)
          if ( subset(3) > 0 ) j1 = subset(3)
          if ( subset(4) > 0 ) j2 = subset(4)
        end if

        ! Find indices of boundary of target domain to save memory
        if ( present( ugg_to ) ) then
        
          ! this will reset the i1,i2,j1,j2 values,
          ! and thus ignore the 'subset' argument ...
          if ( present(subset) ) then
            write (gol,'("argument `subset` not supported in combination with `ugg_to`")'); call goErr
            TRACEBACK; status=1; return
          end if

          ! bounding box for target grid:
          call ugg_to%GetBoundingBox( bbox, status )
          IF_NOT_OK_RETURN(status=1)

          ! shift longitude from [0,360] to [-360,0] if target longitudes
          ! are defined for negative longitudes only;
          ! check if east of target domain is left of first longitude:
          if ( bbox(2) < lons(1) ) then
            ! shit 360 degrees to left:
            lons      = lons      - 360.0
            lons_bnds = lons_bnds - 360.0
          end if

          ! find west index:
          call Interval_Modulo( lons, bbox(1), 360.0, i1, status )
          if ( status /= 0 ) then
            write (gol,'("could not find west ",f8.2," of target bounding box in longitude axis:")') bbox(1); call goErr
            do i = 1, nlon
              write (gol,'("  ",i4,f8.2)') i, lons(i); call goErr
            end do
            TRACEBACK; status=1; return
          end if

          ! find east index:
          call Interval_Modulo( lons, bbox(2), 360.0, i2, status )
          if ( status /= 0 ) then
            write (gol,'("could not find east ",f8.2," of target bounding box in longitude axis:")') bbox(2); call goErr
            do i = 1, nlon
              write (gol,'("  ",i4,f8.2)') i, lons(i); call goErr
            end do
            TRACEBACK; status=1; return
          end if

          ! need to shift slab [180,360] to [-180,0] ?
          if ( i1 > i2 ) then
            ! shift the 'i1' index to a negative value,
            ! the read routine should handle this:
            i1 = i1 - nlon
            ! shift the lon arrays to new index space;
            ! storage for copy:
            allocate( lons2(nlon), stat=status )
            IF_NOT_OK_RETURN(status=1)
            allocate( lons2_bnds(nv,nlon), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! copy:
            lons2      = lons
            lons2_bnds = lons_bnds
            ! clear original:
            deallocate( lons, stat=status )
            IF_NOT_OK_RETURN(status=1)
            deallocate( lons_bnds, stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! allocate again with new index space:
            allocate( lons(1-nlon/2:nlon/2), stat=status )
            IF_NOT_OK_RETURN(status=1)
            allocate( lons_bnds(nv,1-nlon/2:nlon/2), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! reset original:
            do i = 1-nlon/2, 0
              lons     (  i) = lons2     (  i+nlon) - 360.0
              lons_bnds(:,i) = lons2_bnds(:,i+nlon) - 360.0
            end do
            do i = 1, nlon/2
              lons     (  i) = lons2     (  i)
              lons_bnds(:,i) = lons2_bnds(:,i)
            end do
            !! testing ...
            !write (gol,*) 'xxx shifted bounds:'; call goPr
            !do i = lbound(lons,1), ubound(lons,1)
            !  write (gol,*) 'xxx  ', i, lons(i), lons_bnds(:,i); call goPr
            !end do
            !write (gol,*) 'xxx break'; call goPr
            !TRACEBACK; status=1; return
          end if

          ! Find indices, south-north
          if ( lats(1) < lats(nlat) ) then
            j1=0
            call Interval( lats, bbox(3), j1, status )
            if ( status /= 0 ) then
              write (gol,'("could not find interval holding bbox south ",f10.4," in lats range [",f10.4,",",f10.4,"]")') &
                                bbox(3), lats(1), lats(nlat); call goErr
              TRACEBACK; status=1; return
            end if
            j2 = j1
            call Interval( lats, bbox(4), j2, status )
            if ( status /= 0 ) then
              write (gol,'("could not find interval holding bbox north ",f10.4," in lats range [",f10.4,",",f10.4,"]")') &
                                bbox(4), lats(1), lats(nlat); call goErr
              TRACEBACK; status=1; return
            end if
          else
            j1=0
            call Interval( -1.0*lats, -1.0*bbox(4), j1, status )
            if ( status /= 0 ) then
              write (gol,'("could not find interval holding bbox north ",f10.4," in lats range [",f10.4,",",f10.4,"]")') &
                                bbox(4), lats(1), lats(nlat); call goErr
              TRACEBACK; status=1; return
            end if
            j2 = j1
            call Interval( -1.0*lats, -1.0*bbox(3), j2, status )
            if ( status /= 0 ) then
              write (gol,'("could not find interval holding bbox south ",f10.4," in lats range [",f10.4,",",f10.4,"]")') &
                                bbox(3), lats(1), lats(nlat); call goErr
              TRACEBACK; status=1; return
            end if
          end if

          ! Add 2 extra cells at boundary to be sure
          i1 = max( lbound(lons,dim=1), i1-2 )
          i2 = min( ubound(lons,dim=1), i2+2 )
          j1 = max(1   ,j1-2)
          j2 = min(nlat,j2+2)

          ! Fill in start and count indices to read from netCDF file ;
          ! expected to be present ...
          if ( .not. ( present(start_ind) .and. present(count_ind) ) ) then
            write( gol, '("routine called with with optional `ugg_to` argument for target grid,")' ) ; call goErr
            write( gol, '("but without `start` and `count` arguments to return index range for GetVar")' ) ; call GoErr
            TRACEBACK;status=1;return
          end if
          ! fill start and count:
          start_ind = (/i1,j1/)
          count_ind = (/i2-i1+1,j2-j1+1/)

        else

          ! no target grid, thus no index range:
          if ( present(start_ind) .or. present(count_ind) ) then
            write( gol, '("routine called with with optional `start` or `count` arguments,")' ) ; call goErr
            write( gol, '("but without `ugg_to` arguments to define target grid")' ) ; call GoErr
            TRACEBACK;status=1;return
          end if

        end if  ! ugg_to present

        ! define grid:
        call ugg%Init( lons(i1:i2), lons_bnds(:,i1:i2), lats(j1:j2), lats_bnds(:,j1:j2), status )
        IF_NOT_OK_RETURN(status=1)   

        ! clear:
        deallocate( lons, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( lons_bnds, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( lats, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( lats_bnds, stat=status )
        IF_NOT_OK_RETURN(status=1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2-dimensional, 
      ! or 3-dimensional with last dim (time?) ignored:
      case ( '2D', '3D' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        ! grid arrays:
        allocate( xlons(nlon,nlat), stat=status)
        IF_NOT_OK_RETURN(status=1)
        allocate( xlats(nlon,nlat), stat=status)
        IF_NOT_OK_RETURN(status=1)
    
        ! read from 2D or 3D arrays:
        if ( trim(dim_case) == '2D' ) then
          ! read:
          status = NF90_Get_Var( self%ncid, varid_lon, xlons )
          IF_NF90_NOT_OK_RETURN(status=1)
          status = NF90_Get_Var( self%ncid, varid_lat, xlats )
          IF_NF90_NOT_OK_RETURN(status=1)
        else if ( trim(dim_case) == '3D' ) then
          ! read:
          status = NF90_Get_Var( self%ncid, varid_lon, xlons, &
                                   start=(/1,1,1/), count=(/nlon,nlat,1/) )
          IF_NF90_NOT_OK_RETURN(status=1)
          status = NF90_Get_Var( self%ncid, varid_lat, xlats, &
                                   start=(/1,1,1/), count=(/nlon,nlat,1/) )
          IF_NF90_NOT_OK_RETURN(status=1)
        else
          write (gol,'("could not read coordinate arrays for dim_case `",a,"`")') trim(dim_case); call goErr
          TRACEBACK; status=1; return
        end if

        ! by default all:
        i1 = 1
        i2 = nlon
        j1 = 1
        j2 = nlat

        ! subset defined?
        if ( present(subset) ) then
          if ( subset(1) > 0 ) i1 = subset(1)
          if ( subset(2) > 0 ) i2 = subset(2)
          if ( subset(3) > 0 ) j1 = subset(3)
          if ( subset(4) > 0 ) j2 = subset(4)
        end if

        ! Find indices of boundary of target domain to save memory
        if ( present( ugg_to ) ) then
        
          ! this will reset the i1,i2,j1,j2 values,
          ! and thus ignore the 'subset' argument ...
          if ( present(subset) ) then
            write (gol,'("argument `subset` not supported in combination with `ugg_to`")'); call goErr
            TRACEBACK; status=1; return
          end if

          ! bounding box for target grid:
          call ugg_to%GetBoundingBox( bbox, status )
          IF_NOT_OK_RETURN(status=1)
    
          ! quick scan to select a small i-range
          ! that contains the target box;
          !~ direction:
          if ( xlons(1,1) < xlons(nlon,1) ) then
            idir = 1
          else
            idir = -1
          end if
          !~ loop over columns:
          do i = 1, nlon
            ! all left of target ?
            if ( all( idir*xlons(i,:) <= idir*bbox(1) ) ) then
              ! reset:
              i1 = max( i, i1 )
            end if
            ! all right of traget?
            if ( all( idir*xlons(i,:) >= idir*bbox(2) ) ) then
              ! reset:
              i2 = min( i, i2 )
              ! found end, leave:
              exit
            end if
          end do ! i

          ! quick scan to select a small j-range
          ! that contains the target point ;
          !~ direction:
          if ( xlats(1,1) < xlats(1,nlat) ) then
            jdir = 1
          else
            jdir = -1
          end if
          !~ loop over columns:
          do j = 1, nlat
            ! all below target ?
            if ( all( jdir*xlats(i1:i2,j) <= jdir*bbox(3) ) ) then
              ! reset:
              j1 = max( j, j1 )
            end if
            ! all above traget?
            if ( all( jdir*xlats(i1:i2,j) >= jdir*bbox(4) ) ) then
              ! reset:
              j2 = min( j, j2 )
              ! found end, leave:
              exit
            end if
          end do ! j

          ! Fill in start and count indices to read from netCDF file ;
          ! expected to be present ...
          if ( .not. ( present(start_ind) .and. present(count_ind) ) ) then
            write( gol, '("routine called with with optional `ugg_to` argument for target grid,")' ) ; call goErr
            write( gol, '("but without `start` and `count` arguments to return index range for GetVar")' ) ; call GoErr
            TRACEBACK;status=1;return
          end if
          ! fill start and count:
          start_ind = (/i1,j1/)
          count_ind = (/i2-i1+1,j2-j1+1/)

        else

          ! no target grid, thus no index range:
          if ( present(start_ind) .or. present(count_ind) ) then
            write( gol, '("routine called with with optional `start` or `count` arguments,")' ) ; call goErr
            write( gol, '("but without `ugg_to` arguments to define target grid")' ) ; call GoErr
            TRACEBACK;status=1;return
          end if

        end if  ! ugg_to present

        ! define grid:
        call ugg%Init( xlons(i1:i2,j1:j2), xlats(i1:i2,j1:j2), status )
        IF_NOT_OK_RETURN(status=1)   

        ! clear:
        deallocate( xlons, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( xlats, stat=status )
        IF_NOT_OK_RETURN(status=1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! other dim_case
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~
        write( gol, '("Unknown dimension case for grid definition: ", a)' ) trim(dim_case); call GoErr
        TRACEBACK;status=1;return

    end select
    
    ! clear:
    deallocate( dimids, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( dimnames, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine File_Ugg_Get_Grid
  

  ! ***
  
  
!  subroutine File_Ugg_Def_Grid( self, ugg, status )
!
!    use NetCDF, only : NF90_SHORT, NF90_INT, NF90_FLOAT
!    use NetCDF, only : NF90_Def_Dim
!    use NetCDF, only : NF90_Def_Var
!    use NetCDF, only : NF90_Put_Att
!
!    use C3PO_Grid_RGG, only : T_Grid_ugg
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_File_Ugg), intent(inout)    ::  self
!    class(T_Grid_ugg), intent(in)       ::  ugg
!    integer, intent(out)                ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/File_Ugg_Def_Grid'
!    
!    ! --- local ----------------------------------
!    
!    integer               ::  varid
!    character(len=1024)   ::  line
!
!    ! --- begin ----------------------------------
!
!    ! define dimension:
!    status = NF90_Def_Dim( self%ncid, dimname_lat, ugg%nlat, self%dimid_ulat )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! define dimension:
!    status = NF90_Def_Dim( self%ncid, dimname_lon, ugg%nulon, self%dimid_ulon )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! define dimension:
!    status = NF90_Def_Dim( self%ncid, dimname_rgrid, ugg%npoint, self%dimid_rgrid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! define dimension:
!    status = NF90_Def_Dim( self%ncid, dimname_nv, ugg%nv, self%dimid_nv )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! ~ index
!    
!    ! add index variable:
!    status = NF90_Def_Var( self%ncid, dimname_rgrid, NF90_INT, &
!                              (/self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'reduced grid index' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    write (line,'(a," ",a)') dimname_lat, dimname_lon
!    status = NF90_Put_Att( self%ncid, varid, 'compress', trim(line) )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    write (line,'("zero based index that specifies location in regular (",a,",",a,") grid : ", &
!                & " j = rgrid(n) / ",a," ; i = rgrid(n) - ",a," * j")') &
!             dimname_lat, dimname_lon, dimname_lon, dimname_lon
!    IF_NF90_NOT_OK_RETURN(status=1)
!    status = NF90_Put_Att( self%ncid, varid, 'description', trim(line) )
!    !
!    self%varid_irgrid = varid
!        
!    ! ~ centers
!    
!    ! longitude values for each point:
!    status = NF90_Def_Var( self%ncid, varname_lon, NF90_FLOAT, &
!                              (/self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'longitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'longitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_east' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'bounds', varname_lon//varext_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_longitude = varid
!    
!    ! latitude values for each point:
!    status = NF90_Def_Var( self%ncid, varname_lat, NF90_FLOAT, &
!                              (/self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'latitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'bounds', varname_lat//varext_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_latitude = varid
!
!    ! ~ boundary values
!    
!    ! longitude values for each point:
!    status = NF90_Def_Var( self%ncid, varname_lon//varext_bnds, NF90_FLOAT, &
!                              (/self%dimid_nv,self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'longitude cell boundaries' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_east' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_longitude_bnds = varid
!    
!    ! latitude values for each point:
!    status = NF90_Def_Var( self%ncid, varname_lat//varext_bnds, NF90_FLOAT, &
!                              (/self%dimid_nv,self%dimid_rgrid/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude cell boundaries' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_latitude_bnds = varid
!
!    ! ~ latitude bands (extra info)
!    
!    ! longitude values for each point:
!    status = NF90_Def_Var( self%ncid, dimname_lat, NF90_FLOAT, &
!                              (/self%dimid_ulat/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'standard_name', 'latitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'bounds', dimname_lat//varext_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_band_lats = varid
!    
!    ! latitude values for each point:
!    status = NF90_Def_Var( self%ncid, dimname_lat//varext_bnds, NF90_FLOAT, &
!                              (/self%dimid_nv,self%dimid_ulat/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'latitude cell boundaries' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', 'degrees_north' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_band_lats_bnds = varid
!    
!    ! number of longitudes per row:
!    status = NF90_Def_Var( self%ncid, dimname_lat//'_nlon', NF90_SHORT, &
!                              (/self%dimid_ulat/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'number of longitudes on latitude band' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_band_nlon = varid
!
!    ! index in reduced grid:
!    status = NF90_Def_Var( self%ncid, dimname_lat//'_i0', NF90_INT, &
!                              (/self%dimid_ulat/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'zero-based index in reduced grid of first point on latitude band' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    !
!    self%varid_band_i0 = varid
!    
!    ! ok
!    status = 0
!    
!  end subroutine File_Ugg_Def_Grid
!
!
!  ! ***
!  
!  
!  subroutine File_Ugg_Put_Grid( self, ugg, status )
!
!    use NetCDF, only : NF90_Put_Var
!
!    use C3PO_Grid_RGG, only : T_Grid_ugg
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_File_Ugg), intent(inout)    ::  self
!    class(T_Grid_ugg), intent(in)       ::  ugg
!    integer, intent(out)                ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/File_Ugg_Put_Grid'
!    
!    ! --- local ----------------------------------
!    
!    ! --- begin ----------------------------------
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_irgrid, ugg%irgrid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_longitude, ugg%longitude )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_latitude, ugg%latitude )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_longitude_bnds, ugg%longitude_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_latitude_bnds, ugg%latitude_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_band_lats, ugg%band_lats )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_band_lats_bnds, ugg%band_lats_bnds )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_band_nlon, ugg%band_nlon )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_band_i0, ugg%band_i0 )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    
!    ! ok
!    status = 0
!    
!  end subroutine File_Ugg_Put_Grid
!
!
!  ! ***
!  
!  
!  subroutine File_Ugg_Copy_Variable_Selection( self, varname, outfile, selection, status )
!
!    use NetCDF, only : NF90_INT, NF90_FLOAT, NF90_DOUBLE
!    use NetCDF, only : NF90_Inquire_Dimension
!    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
!    use NetCDF, only : NF90_Get_Var, NF90_Put_Var
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_File_Ugg), intent(in)     ::  self
!    character(len=*), intent(in)      ::  varname
!    class(T_File_Ugg), intent(inout)  ::  outfile
!    integer, intent(in)               ::  selection(:)  ! 1-based indices in rgrid
!    integer, intent(out)              ::  status
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/File_Ugg_Copy_Variable_Selection'
!    
!    ! --- local ----------------------------------
!    
!    integer             ::  varid, varid_out
!    integer             ::  xtype
!    integer             ::  ndim, idim
!    integer             ::  dimids(maxdim)
!    integer             ::  shp(maxdim)
!    logical             ::  rgridim(maxdim)
!    character(len=64)   ::  dimname
!    integer             ::  nselected
!    integer             ::  i
!    
!    integer(4), allocatable   ::  data_i4_1d(:), data_i4_1d_out(:)
!    real(4), allocatable      ::  data_r4_1d(:), data_r4_1d_out(:)
!    real(8), allocatable      ::  data_r8_1d(:), data_r8_1d_out(:)
!    
!    integer(4), allocatable   ::  data_i4_2d(:,:), data_i4_2d_out(:,:)
!    real(4), allocatable      ::  data_r4_2d(:,:), data_r4_2d_out(:,:)
!    real(8), allocatable      ::  data_r8_2d(:,:), data_r8_2d_out(:,:)
!    
!    integer(4), allocatable   ::  data_i4_3d(:,:,:), data_i4_3d_out(:,:,:)
!    real(4), allocatable      ::  data_r4_3d(:,:,:), data_r4_3d_out(:,:,:)
!    real(8), allocatable      ::  data_r8_3d(:,:,:), data_r8_3d_out(:,:,:)
!    
!    ! --- begin ----------------------------------
!    
!    ! id of input variable:
!    status = NF90_INQ_VarID( self%ncid, varname, varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    
!    ! obtain info:
!    status = NF90_Inquire_Variable( self%ncid, varid, xtype=xtype, &
!                                       ndims=ndim, dimids=dimids )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! obtain shape:
!    do idim = 1, ndim
!      ! get length:
!      status = NF90_Inquire_Dimension( self%ncid, dimids(idim), len=shp(idim) )
!      IF_NF90_NOT_OK_RETURN(status=1)
!    end do
!    
!    ! init flags for rgrid dim:
!    rgridim = .false.
!    do idim = 1, ndim
!      ! get name:
!      status = NF90_Inquire_Dimension( self%ncid, dimids(idim), name=dimname )
!      IF_NF90_NOT_OK_RETURN(status=1)
!      ! compare:
!      rgridim(idim) = trim(dimname) == trim(dimname_rgrid)
!    end do
!    ! check ...
!    if ( count(rgridim) > 1 ) then
!      write (gol,'("could not select points for ",i2," rgrid dimensions")') count(rgridim); call goErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! id of output variable:
!    status = NF90_INQ_VarID( outfile%ncid, varname, varid_out )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    
!    ! number of selected grid points:
!    nselected = size( selection )
!    
!    ! switch on rank:
!    select case ( ndim )
!
!      !~ 1D
!      case ( 1 )
!        ! create storage and read:
!        select case ( xtype )
!        
!          !~ i4
!          case ( NF90_INT )
!            ! storage:
!            allocate( data_i4_1d(shp(1)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_i4_1d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! storage:
!              allocate( data_i4_1d_out(nselected) )
!              ! copy:
!              do i = 1, nselected
!                data_i4_1d_out(i) = data_i4_1d(selection(i))
!              end do
!            else
!              ! storage:
!              allocate( data_i4_1d_out(shp(1)) )
!              ! copy:
!              data_i4_1d_out = data_i4_1d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_1d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_i4_1d )
!            deallocate( data_i4_1d_out )
!        
!          !~ r4
!          case ( NF90_FLOAT )
!            ! storage:
!            allocate( data_r4_1d(shp(1)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r4_1d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! storage:
!              allocate( data_r4_1d_out(nselected) )
!              ! copy:
!              do i = 1, nselected
!                data_r4_1d_out(i) = data_r4_1d(selection(i))
!              end do
!            else
!              ! storage:
!              allocate( data_r4_1d_out(shp(1)) )
!              ! copy:
!              data_r4_1d_out = data_r4_1d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_1d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r4_1d )
!            deallocate( data_r4_1d_out )
!        
!          !~ r8
!          case ( NF90_DOUBLE )
!            ! storage:
!            allocate( data_r8_1d(shp(1)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r8_1d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! storage:
!              allocate( data_r8_1d_out(nselected) )
!              ! copy:
!              do i = 1, nselected
!                data_r8_1d_out(i) = data_r8_1d(selection(i))
!              end do
!            else
!              ! storage:
!              allocate( data_r8_1d_out(shp(1)) )
!              ! copy:
!              data_r8_1d_out = data_r8_1d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_1d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r8_1d )
!            deallocate( data_r8_1d_out )
!              
!          !~ unkown
!          case default
!            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
!            TRACEBACK; status=1; return
!        end select
!
!      !~ 2D
!      case ( 2 )
!        ! create storage and read:
!        select case ( xtype )
!        
!          !~ i4
!          case ( NF90_INT )
!            ! storage:
!            allocate( data_i4_2d(shp(1),shp(2)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_i4_2d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_i4_2d_out(nselected,shp(2)) )
!                ! copy:
!                do i = 1, nselected
!                  data_i4_2d_out(i,:) = data_i4_2d(selection(i),:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_i4_2d_out(shp(1),shp(2)) )
!              ! copy:
!              data_i4_2d_out = data_i4_2d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_2d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_i4_2d )
!            deallocate( data_i4_2d_out )
!        
!          !~ r4
!          case ( NF90_FLOAT )
!            ! storage:
!            allocate( data_r4_2d(shp(1),shp(2)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r4_2d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_r4_2d_out(nselected,shp(2)) )
!                ! copy:
!                do i = 1, nselected
!                  data_r4_2d_out(i,:) = data_r4_2d(selection(i),:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_r4_2d_out(shp(1),shp(2)) )
!              ! copy:
!              data_r4_2d_out = data_r4_2d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_2d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r4_2d )
!            deallocate( data_r4_2d_out )
!        
!          !~ r8
!          case ( NF90_DOUBLE )
!            ! storage:
!            allocate( data_r8_2d(shp(1),shp(2)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r8_2d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_r8_2d_out(nselected,shp(2)) )
!                ! copy:
!                do i = 1, nselected
!                  data_r8_2d_out(i,:) = data_r8_2d(selection(i),:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_r8_2d_out(shp(1),shp(2)) )
!              ! copy:
!              data_r8_2d_out = data_r8_2d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_2d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r8_2d )
!            deallocate( data_r8_2d_out )
!              
!          !~ unkown
!          case default
!            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
!            TRACEBACK; status=1; return
!        end select
!
!      !~ 3D
!      case ( 3 )
!        ! create storage and read:
!        select case ( xtype )
!        
!          !~ i4
!          case ( NF90_INT )
!            ! storage:
!            allocate( data_i4_3d(shp(1),shp(2),shp(3)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_i4_3d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_i4_3d_out(nselected,shp(2),shp(3)) )
!                ! copy:
!                do i = 1, nselected
!                  data_i4_3d_out(i,:,:) = data_i4_3d(selection(i),:,:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_i4_3d_out(shp(1),shp(2),shp(3)) )
!              ! copy:
!              data_i4_3d_out = data_i4_3d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_i4_3d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_i4_3d )
!            deallocate( data_i4_3d_out )
!        
!          !~ r4
!          case ( NF90_FLOAT )
!            ! storage:
!            allocate( data_r4_3d(shp(1),shp(2),shp(3)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r4_3d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_r4_3d_out(nselected,shp(2),shp(3)) )
!                ! copy:
!                do i = 1, nselected
!                  data_r4_3d_out(i,:,:) = data_r4_3d(selection(i),:,:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_r4_3d_out(shp(1),shp(2),shp(3)) )
!              ! copy:
!              data_r4_3d_out = data_r4_3d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r4_3d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r4_3d )
!            deallocate( data_r4_3d_out )
!        
!          !~ r8
!          case ( NF90_DOUBLE )
!            ! storage:
!            allocate( data_r8_3d(shp(1),shp(2),shp(3)) )
!            ! read:
!            status = NF90_Get_Var( self%ncid, varid, data_r8_3d )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! need to select ?
!            if ( any(rgridim) ) then
!              ! which ?
!              if ( rgridim(1) ) then
!                ! storage:
!                allocate( data_r8_3d_out(nselected,shp(2),shp(3)) )
!                ! copy:
!                do i = 1, nselected
!                  data_r8_3d_out(i,:,:) = data_r8_3d(selection(i),:,:)
!                end do
!              else
!                ! not yet ...
!                write (gol,'("unsupported rgrid dimensions : ",2l2)') rgridim(1:ndim); call goErr
!                TRACEBACK; status=1; return
!              end if
!            else
!              ! storage:
!              allocate( data_r8_3d_out(shp(1),shp(2),shp(3)) )
!              ! copy:
!              data_r8_3d_out = data_r8_3d
!            end if
!            ! write:
!            status = NF90_Put_Var( outfile%ncid, varid_out, data_r8_3d_out )
!            IF_NF90_NOT_OK_RETURN(status=1)
!            ! clear:
!            deallocate( data_r8_3d )
!            deallocate( data_r8_3d_out )
!              
!          !~ unkown
!          case default
!            write (gol,'("unsupported xtype ",i6)') xtype; call goErr
!            TRACEBACK; status=1; return
!        end select
!        
!      !~
!      case default
!        write (gol,'("unsupported ndim ",i6)') ndim; call goErr
!        TRACEBACK; status=1; return
!    end select
!    
!    ! ok
!    status = 0
!    
!  end subroutine File_Ugg_Copy_Variable_Selection


end module C3PO_File_Ugg
