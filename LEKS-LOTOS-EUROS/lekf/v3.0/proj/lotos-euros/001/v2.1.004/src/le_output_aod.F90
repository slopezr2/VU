!###############################################################################
!
! NAME
!
!   LE_Output_aod  -  LOTOS-EUROS output of aod fields
!
! HISTORY
!
!   2008 oct, Arjo Segers, TNO
!
!
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

module LE_Output_aod

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate
  
#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif
  
  use LE_Output_Common, only : T_LE_Output_Common
  
  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  T_LE_Output_aod
  
  public  ::  LE_Output_AOD_Init, LE_Output_AOD_Done
  public  ::  LE_Output_AOD_PutOut


  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LE_Output_aod'
  
    
  ! --- types ------------------------------
  
  type T_LE_Output_aod
    ! name for this file:
    character(len=16)           ::  typ
    character(len=16)           ::  name
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
    ! adhoc convention, might change per file ...
    character(len=32)           ::  convention
    ! replace existing files ?
    logical                     ::  replace
    ! file opened ?
    logical                     ::  opened
    ! pack output ?
    logical                     ::  pack_output
    ! current time range:
    type(TDate)                 ::  tr(2)
    ! time resolution:
    real                        ::  dhour
    ! time record counter:
    integer                     ::  itrec
    ! state name:
    character(len=16)           ::  state
    ! file name:
    character(len=256)          ::  fname
    ! file handle:
    integer                     ::  ncid
    ! dimension handles:
    integer                     ::  dimid_lon
    integer                     ::  dimid_lat
    integer                     ::  dimid_lev
    integer                     ::  dimid_time
    integer                     ::  varid_time_dtg
    ! dimension variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_time
    ! output variables:
    character(len=32)           ::  units
    real                        ::  unitconv
    integer                     ::  varid
    character(len=32)           ::  varname
    ! other
    logical                     ::  write_biascorr
    ! grads ctl file ?
    logical                     ::  grads_ctl
    character(len=256)          ::  grads_ctl_file
    character(len=256)          ::  grads_ctl_base
    integer                     ::  grads_ctl_nt
    type(TDate)                 ::  grads_ctl_t1
    type(TIncrDate)             ::  grads_ctl_dt
    ! bounding box
    integer                     ::  i1, i2, ni
    integer                     ::  j1, j2, nj
    real                        ::  westb, southb
  end type T_LE_Output_aod

  
contains


  ! ====================================================
  
  
  subroutine LE_Output_AOD_Init( leo, rcfile, rckey, typ, name, state, status )

    use GO     , only : TrcFile, Init, Done, ReadRc
    use GO     , only : goMatchValues
    use GO     , only : AnyDate
    use Dims   , only : nx, ny, nz
    use Indices, only : nspec
    !use Dims   , only : nspec_maxaccum
    use Indices, only : specname, specunit
    use LE_Grid, only : ugg
    use LE_Output_Common, only : Init
  
    ! --- in/out --------------------------------
    
    type(T_LE_Output_aod), intent(out)    ::  leo
    character(len=*), intent(in)          ::  rcfile
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  state
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Output_AOD_Init'
    
    ! --- local ---------------------------------
    
    type(TrcFile)         ::  rcF
    character(len=64)     ::  basekey

    character(len=512)    ::  key
    real                  ::  west, south, east, north
    real, pointer         ::  ff(:,:)

    ! --- begin ---------------------------------
    
    ! store name:
    leo%typ  = typ
    leo%name = name
    
    ! init common stuff:
    call Init( leo%com, rcfile, rckey, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open rcfile:
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)

    ! replace existing files ?
    call ReadRc( rcF, trim(rckey)//'.replace', leo%replace, status )
    IF_NOTOK_RETURN(status=1)

    ! write GrADS ctl file ?
    call ReadRc( rcF, trim(rckey)//'.ctl', leo%grads_ctl, status )
    IF_NOTOK_RETURN(status=1)

    ! base key:
    write (basekey,'(a,".",a,".",a)') trim(rckey), trim(typ), trim(name)

    ! output convention:
    call ReadRc( rcF, trim(basekey)//'.convention', leo%convention, status )
    IF_NOTOK_RETURN(status=1)

    ! output time resolution:
    call ReadRc( rcF, trim(basekey)//'.dhour', leo%dhour, status )
    IF_NOTOK_RETURN(status=1)
    
    ! pack output?
    call ReadRC( rcF, trim(basekey)//'.pack_output', leo%pack_output, status )
    IF_NOTOK_RETURN(status=1)

    ! write bias corrected aod ?
    call ReadRc( rcF, trim(basekey)//'.biascorr', leo%write_biascorr, status )
    IF_NOTOK_RETURN(status=1)

    ! state name:
    leo%state = trim(state)
    
    ! files not open yet:
    leo%opened = .false.
    
    ! no time range set yet:
    leo%tr(1) = AnyDate()
    leo%tr(2) = AnyDate()

    ! init GrADS stuff:
    if ( leo%grads_ctl ) then
      ! no times written yet:
      leo%grads_ctl_nt = 0
    end if

    ! bounding box
    call ReadRc( rcF, trim(basekey)//'.bounding_box', key, status )
    IF_NOTOK_RETURN(status=1)
    ! empty?
    if (len_trim(key) == 0) then
      ! full domain
      leo%i1 = 1
      leo%i2 = nx
      leo%ni = nx
      leo%j1 = 1
      leo%j2 = ny
      leo%nj = ny
      leo%westb = ugg%longitude_bnds(1,1,1)
      leo%southb = ugg%latitude_bnds(1,1,1)
    else
      
      select case ( trim(ugg%type) ) 
        
        case ( 'cartesian-regular') 
          ! read domain from key
          read(key,*,iostat=status) west, east, south, north
          if(status/=0) then
            write (gol,'("could not read domain from key: ",a)') trim(key); call goErr
            TRACEBACK; status=1; return
          endif

          ! Check if bounding box is in run domain
          if ( west < ugg%longitude_bnds_1d(1,1) .or. east > ugg%longitude_bnds_1d(2,ugg%nlon) .or. &
               south < ugg%latitude_bnds_1d(1,1) .or. north > ugg%latitude_bnds_1d(2,ugg%nlat) ) then
            write( gol, '("Bounding box domain is (partly) outside run domain")' ) ; call goErr
            write( gol, '("Run domain: ", 4f8.2)' ) ugg%longitude_bnds_1d(1,1),ugg%longitude_bnds_1d(2,ugg%nlon),ugg%latitude_bnds_1d(1,1),ugg%latitude_bnds_1d(2,ugg%nlat); call goErr
            write( gol, '("Bounding Box domain: ", 4f8.2)' ) west, east, south, north ; call goErr
            TRACEBACK;status=1;return
          endif

          ! for safety
          nullify(ff)
          ! get cell range covered by box
          call ugg%GetDistribution(west,east,south,north,leo%i1,leo%i2,leo%j1,leo%j2,ff,status)
          IF_NOTOK_RETURN(status=1)
          !clear, fractions not used
          if ( associated(ff) ) deallocate(ff)
          ! set shape
          leo%ni = leo%i2-leo%i1+1
          leo%nj = leo%j2-leo%j1+1
          ! set west/south bounds
          leo%westb  = ugg%longitude_bnds_1d(1,leo%i1)
          leo%southb = ugg%latitude_bnds_1d(1,leo%j1)
        case default 
          write( gol, '("Definition of bounding box not clear for grid-type: ", a)' ) trim(ugg%type) ; call goErr
          TRACEBACK;status=1;return
      end select
    end if
    
    ! close
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Output_AOD_Init
  
  
  ! ***
  

  subroutine LE_Output_AOD_Done( leo, status )
  
#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif
    use LE_Output_Common, only : Done
  
    ! --- in/out --------------------------------
    
    type(T_LE_Output_aod), intent(inout)   ::  leo
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Output_AOD_Done'

    character(len=256) :: commandline
    
    ! --- begin ---------------------------------
    
    ! file opened ?
    if ( leo%opened ) then
      ! close:
#ifdef with_netcdf
      status = NF90_Close( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)

      ! pack output?
      if ( leo%pack_output ) then
        write( commandline, '(a," " a)' ) trim(leo%com%packscript), trim(leo%fname)
        call SYSTEM( trim(commandline) )
      end if      
#else
      stop 'not compiled with netcdf support'
#endif
      ! reset flag:
      leo%opened = .false.
      ! write GrADS ctl file if necessary:
      call Write_GrADS_Ctl( leo, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! done with common stuff ...
    call Done( leo%com, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Output_AOD_Done
  
  
  ! ***
  
  
  subroutine LE_Output_AOD_PutOut( leo, t, c, status )

    use GO     , only : goc
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get
    use GO     , only : operator(+), operator(-), operator(<), iTotal
    use GO     , only : wrtgol, Precisely
    use LE_Grid, only : ugg, glb_ugg
    use C3PO   , only : T_Grid_NcDef

#ifdef with_netcdf
    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT
#endif

    use LE_CF_Conventions, only : LE_CF_names

    use Dims   , only : nx, ny, nz, nspec
    use Indices, only : specname, specunit
    use LE_AOD , only : LE_AOD_Calc3d

    use LE_Output_Common, only : PutOut_GlobalAttributes
    
    use LE_Output_Tools , only : LE_Output_Define_Dims_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Define_Dims_Lev
    use LE_Output_Tools , only : LE_Output_Define_Dims_Time

    use LE_Output_Tools , only : LE_Output_Define_Vars_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Define_Vars_Lev
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time

    use LE_Output_Tools , only : LE_Output_Put_Var_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Put_Var_Lev
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains

    use LE_BiasCorr, only : biascorr__aod__factor
    
    ! --- in/out --------------------------------
    
    type(T_LE_Output_aod), intent(inout)    ::  leo
    type(TDate), intent(in)                 ::  t
    real, intent(in)                        ::  c (nx,ny,nz,nspec)
    integer, intent(out)                    ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Output_AOD_PutOut'
    
    ! --- local ---------------------------------
    
    integer               ::  time6(6)
    integer               ::  time
#ifdef with_netcdf    
    integer               ::  cmode
#endif
!    integer               ::  itr, ilev
    type(T_Grid_NcDef)    ::  gncd
    integer               ::  varid
    type(TDate)           ::  t0
    integer               ::  i, j
    real, allocatable     ::  pat(:,:)
    real, allocatable     ::  field(:,:,:), field2(:,:,:)
    
    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units
    character(len=512)    ::  comment
    character(len=256)    :: commandline
    
    ! --- begin ---------------------------------
    
#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif

    ! for multiples of dhour only ...
    if ( .not. Precisely(t,leo%dhour,'hour')  ) then
      status=0; return
    end if

    ! current time not in time range ?
    if ( (t < leo%tr(1)) .or. (leo%tr(2) < t) ) then

      ! file opened ?
      if ( leo%opened ) then

        ! root?
        if ( goc%root ) then
          ! close:
#ifdef with_netcdf    
          status = NF90_Close( leo%ncid )
          IF_NF90_NOTOK_RETURN(status=1)
#endif
          ! pack output?
          if ( leo%pack_output ) then
            write( commandline, '(a," " a)' ) trim(leo%com%packscript), trim(leo%fname)
            call SYSTEM( trim(commandline) )
          end if

          ! write GrADS ctl file if necessary:
          call Write_GrADS_Ctl( leo, status )
          IF_NOTOK_RETURN(status=1)

        end if ! root

        ! reset flag:
        leo%opened = .false.

      end if ! opened

      ! extract time fields:
      call Get( t, time6=time6 )

      ! set time range [00,24) for this day:
      leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
      leo%tr(2) = leo%tr(1) + IncrDate( day=1 ) - IncrDate(mili=1)

      ! new file name:
      write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
                trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                trim(leo%name), time6(1:3)
      if ( len_trim(leo%state) > 0 ) write (leo%fname,'(a,"_",a)') trim(leo%fname), trim(leo%state)
      write (leo%fname,'(a,".nc")') trim(leo%fname)

      ! root?
      if ( goc%root ) then

#ifdef with_netcdf    
        ! set creation mode flag:
        if ( leo%replace ) then
          cmode = NF90_CLOBBER       ! overwrite existing files
        else
          cmode = NF90_NOCLOBBER     ! do not overwrite existing files
        end if

        ! create file:
        status = NF90_Create( trim(leo%fname), cmode, leo%ncid )
        if ( status /= NF90_NOERR ) then
          write (gol,'("creating file : ")'); call goErr
          write (gol,'("  file name  : ",a)') trim(leo%fname); call goErr
          write (gol,'("  nf90 error : ",a)') trim(nf90_strerror(status)); call goErr
          TRACEBACK; status=1; return
        end if
#else
        ! dummy ..
        leo%ncid = -1
#endif
      
        ! reset flag:
        leo%opened = .true.

        ! write global attributes:
        call PutOut_GlobalAttributes( leo%com, leo%ncid, status )
        IF_NOTOK_RETURN(status=1)
      
        ! grid dimensions/variables
        call glb_ugg%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                    dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )
                                    !subset=(/leo%i1,leo%i2,leo%j1,leo%j2/) )
        IF_NOTOK_RETURN(status=1)

        ! time dimensions /variable
        call LE_Output_Define_Dims_Time(leo%ncid, leo%dimid_time, status)
        IF_NOTOK_RETURN(status=1)

        ! time since t0
        t0 = leo%com%t0      
        ! time variables
        call LE_Output_Define_Vars_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                        leo%dimid_time, trim(leo%com%CF_convention), t0, status)
        IF_NOTOK_RETURN(status=1)

      end if ! root
      
      !
      ! define variables:
      !
      
      ! variable name:
      leo%varname = 'AOD'
      if ( leo%write_biascorr ) leo%varname = trim(leo%varname)//'_biascorr'

      ! no comment yet
      comment = ''
      ! get names following CF conventions;
      ! store conversion factor for later usage:
      call LE_CF_names( &
                   trim(leo%varname), '1', &
                   cf_standard_name, cf_long_name, cf_units, &
                   leo%unitconv, comment, &
                   status )
      IF_NOTOK_RETURN(status=1)

      ! store units for later usage (GrADS ctl file, conversions):
      leo%units = trim(cf_units)

      ! root only:
      if ( goc%root ) then

#ifdef with_netcdf

        ! define variable:
        status = NF90_Def_Var( leo%ncid, trim(leo%varname), NF90_REAL, &
                                 (/leo%dimid_lon,leo%dimid_lat,leo%dimid_time/), varid )
        IF_NF90_NOTOK_RETURN(status=1)

        ! write attributes:
        if ( len_trim(cf_standard_name) > 0 ) then
          status = nf90_put_att( leo%ncid, varid, 'standard_name', trim(cf_standard_name) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
        status = nf90_put_att( leo%ncid, varid, 'long_name', trim(cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( leo%ncid, varid, 'units', trim(cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time level latitude longitude')
        IF_NF90_NOTOK_RETURN(status=1)
        call ugg%DefCoor_NetCDF( gncd, varid, status )
        IF_NF90_NOTOK_RETURN(status=1)
        ! add comment?
        if ( len_trim(comment) > 0 ) then
          status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
#endif
        
        ! store variable id:
        leo%varid = varid

        ! end defintion mode:

#ifdef with_netcdf    
        status = NF90_EndDef( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif

      end if ! root

      ! no records written yet:
      leo%itrec = 0
    
    end if  ! new file
    
    ! next time record:
    leo%itrec = leo%itrec + 1

    ! root only:
    if ( goc%root ) then
    
      ! GrADS time counter:
      if ( leo%grads_ctl ) then
        ! increase counter:
        leo%grads_ctl_nt = leo%grads_ctl_nt + 1
        ! set times if necessary:
        if ( leo%grads_ctl_nt == 1 ) then
          leo%grads_ctl_t1 = t
          leo%grads_ctl_dt = IncrDate(day=1)   ! dummy ...
        end if
        if ( leo%grads_ctl_nt == 2 ) then
          leo%grads_ctl_dt = t - leo%grads_ctl_t1
        end if
      end if

      ! write dimension data only once ...
      if ( leo%itrec == 1 ) then

        ! write grid to netCDF file
        call glb_ugg%PutGrid_NetCDF( gncd, status )
        IF_NOTOK_RETURN(status=1)

      end if  ! first record

      ! date up to seconds:
      call Get( t, time6=time6 )

      ! time since t0
      t0 = leo%com%t0
      time = iTotal( t - t0, 'sec' )

      ! write time record:
#ifdef with_netcdf
      call LE_Output_Put_Var_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                   time, time6, trim(leo%com%CF_convention), leo%itrec, status )
      IF_NOTOK_RETURN(status=1)
#endif

    end if  ! root
    
    ! 2d field:
    allocate( pat(1:nx,1:ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! 3d field:
    allocate( field(1:nx,1:ny,1:nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( field2(1:nx,1:ny,1:nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! compute AOD field:
    call LE_AOD_Calc3d( c, field, field2, status )
    IF_NOTOK_RETURN(status=1)

    ! column:
    pat = sum(field,3)
    ! unit conversion:
    pat = pat * leo%unitconv

    ! apply bias correction ?
    if ( leo%write_biascorr ) then
      ! apply bias correction factor, might be 1.0
      pat = pat * biascorr__aod__factor
    end if

    ! write 2D field, without level, with time index:
    call LE_Output_Put_Var_Domains( leo%ncid, leo%varid, -1, leo%itrec, &
                                      pat, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( pat, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( field, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( field2, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Output_AOD_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use LE_Grid, only : glb_ugg
    
    ! --- in/out ---------------------------------

    type(T_LE_Output_aod), intent(inout)    ::  leo
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Write_GrADS_Ctl'

    ! --- local ----------------------------------

    type(T_GrADS_Ctl)     ::  ctl
    integer               ::  l
    character(len=512)    ::  line

    ! --- begin ----------------------------------

    ! write ctl file ?
    if ( leo%grads_ctl ) then

      ! ctl file name:
      write (leo%grads_ctl_file,'(a,"_",a,"_",a)') &
                trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
      if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_file,'(a,"_",a)') trim(leo%grads_ctl_file), trim(leo%state)
      write (leo%grads_ctl_file,'(a,".ctl")') trim(leo%grads_ctl_file)

      ! template for data files:
      write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2")') &
                trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
      if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_base,'(a,"_",a)') trim(leo%grads_ctl_base), trim(leo%state)
      write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)

      ! open ctl file:
      call GrADS_Ctl_Init( ctl, trim(leo%com%outdir), trim(leo%grads_ctl_file), status )
      IF_NOTOK_RETURN(status=1)
      ! comment ...
      call GrADS_Ctl_Comment( ctl, '', status )
      call GrADS_Ctl_Comment( ctl, 'GrADS Data Descriptor File', status )
      call GrADS_Ctl_Comment( ctl, '', status )
      ! data set:
      call GrADS_Ctl_DSet( ctl, trim(leo%grads_ctl_base), status )
      IF_NOTOK_RETURN(status=1)
      ! title:
      write (line,'("model: ",a,"; expid: ",a)') trim(leo%com%model), trim(leo%com%expid)
      call GrADS_Ctl_Title( ctl, trim(line), status )
      IF_NOTOK_RETURN(status=1)
      ! write xdef/ydef from grid/projection definition
      call glb_ugg%WriteCtlProjection( ctl, status)
      IF_NOTOK_RETURN(status=1)
      ! zdef:
      call GrADS_Ctl_ZDef( ctl, (/0/), status )
      IF_NOTOK_RETURN(status=1)
      ! tdef:
      call GrADS_Ctl_TDef( ctl, leo%grads_ctl_nt, &
                               (/leo%grads_ctl_t1%year,leo%grads_ctl_t1%month,leo%grads_ctl_t1%day,leo%grads_ctl_t1%hour,leo%grads_ctl_t1%min/), &
                               (/                    0,                     0,leo%grads_ctl_dt%day,leo%grads_ctl_dt%hour,leo%grads_ctl_dt%min/), &
                               status )
      IF_NOTOK_RETURN(status=1)
      ! number of variables:
      call GrADS_Ctl_Vars( ctl, 1, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over tracers to be written:
      do l = 1, 1
        ! set variable description:
        write (line,'(a," [",a,"]")') 'AOD', trim(leo%units)
        ! add variable description:
        call GrADS_Ctl_Var( ctl, trim(leo%varname), 1, 't,y,x', trim(line), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! end of variables section:
      call GrADS_Ctl_EndVars( ctl, status )
      IF_NOTOK_RETURN(status=1)
      ! close ctl file:
      call GrADS_Ctl_Done( ctl, status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! ok
    status = 0

  end subroutine Write_GrADS_Ctl


end module LE_Output_aod
