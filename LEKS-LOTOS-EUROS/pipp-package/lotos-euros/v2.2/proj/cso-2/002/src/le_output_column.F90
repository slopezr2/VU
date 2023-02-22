!###############################################################################
!
! NAME
!
!   LE_Output_Column  -  LOTOS-EUROS output of concentration fields
!
! HISTORY
!
!   2007 may, Arjo Segers, TNO
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

module LE_Output_Column

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_Column

  public  ::  LE_Output_Column_Init, LE_Output_Column_Done
  public  ::  LE_Output_Column_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Column'


  ! --- types ------------------------------

  type T_LE_Output_Column
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
    ! current time range:
    type(TDate)                 ::  tr(2)
    ! time resolution:
    real                        ::  dhour
    ! collect: daily, instant
    character(len=32)           ::  collect
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
    integer                     ::  dimid_time
    ! dimension variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_time
    integer                     ::  varid_time_day
    integer                     ::  varid_time_dtg
    integer                     ::  varid_dh
    ! tracer variables:
    integer                     ::  ntr
    integer, pointer            ::  itr(:)
    character(len=32), pointer  ::  name_tr(:)
    character(len=32), pointer  ::  unit_tr(:)
    real, pointer               ::  unitconv(:)
    integer, pointer            ::  varid_tr(:)
    ! level selection:
    character(len=16)           ::  levtype
    integer                     ::  nlev
    ! top presssure:
    real                        ::  ptop
    ! grads ctl file ?
    logical                     ::  grads_ctl
    character(len=256)          ::  grads_ctl_file
    character(len=256)          ::  grads_ctl_base
    integer                     ::  grads_ctl_nt
    type(TDate)                 ::  grads_ctl_t1
    type(TIncrDate)             ::  grads_ctl_dt
    ! local box
    integer                     ::  i1, i2, ni
    integer                     ::  j1, j2, nj
    ! bounding box
    integer                     ::  bb_i1, bb_i2, bb_ni
    integer                     ::  bb_j1, bb_j2, bb_nj
    real                        ::  westb, southb
  end type T_LE_Output_Column


contains


  ! ====================================================


  subroutine LE_Output_Column_Init( leo, rcF, rckey, typ, name, state, status )

    use GO     , only : TrcFile
    use GO     , only : goMatchValues
    use GO     , only : AnyDate
    use Indices, only : specname, specunit
    use Dims   , only : nspec_all
    use Dims   , only : nx,ny
    use LE_Grid, only : glb_ugg

    use LE_Data         , only : LE_Data_Enable
    use LE_Data_Common  , only : nlev, nlev_top

    ! --- in/out --------------------------------

    type(T_LE_Output_Column), intent(out)   ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  state
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Column_Init'

    ! --- local ---------------------------------

    character(len=32)     ::  basekey
    character(len=1024)   ::  tracer_names
    integer               ::  itr
    integer               ::  i
  
    character(len=512)    ::  key
    real                  ::  west, south, east, north
    real, pointer         ::  ff(:,:)

    ! --- begin ---------------------------------

    ! store name:
    leo%typ  = typ
    leo%name = name

    ! init common stuff:
    call leo%com%Init( rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)

    ! replace existing files ?
    call rcF%Get( trim(rckey)//'.replace', leo%replace, status )
    IF_NOTOK_RETURN(status=1)

    ! write GrADS ctl file ?
    call rcF%Get( trim(rckey)//'.ctl', leo%grads_ctl, status )
    IF_NOTOK_RETURN(status=1)

    ! base key:
    write (basekey,'(a,".",a,".",a)') trim(rckey), trim(typ), trim(name)

    ! output convention:
    call rcF%Get( trim(basekey)//'.convention', leo%convention, status )
    IF_NOTOK_RETURN(status=1)

    ! level type:
    call rcF%Get( trim(basekey)//'.levtype', leo%levtype, status )
    IF_NOTOK_RETURN(status=1)
    ! layers or heights ?
    select case ( trim(leo%levtype) )
      !~ model levels:
      case ( 'all' )
        leo%nlev = nlev
      !~ model levels and top levels:
      case ( 'top' )
        leo%nlev = nlev_top
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! top pressure:
    call rcF%Get( trim(basekey)//'.ptop', leo%ptop, status )
    IF_NOTOK_RETURN(status=1)

    ! collect daily or instant
    call rcF%Get( trim(basekey)//'.collect', leo%collect, status )
    IF_NOTOK_RETURN(status=1)

    ! output time resolution:
    call rcF%Get( trim(basekey)//'.dhour', leo%dhour, status )
    IF_NOTOK_RETURN(status=1)

    ! tracer names:
    call rcF%Get( trim(basekey)//'.fields', tracer_names, status )
    IF_NOTOK_RETURN(status=1)

    ! state name:
    leo%state = trim(state)

    ! setup storage for tracer fields:
    allocate( leo%itr     (nspec_all) )
    allocate( leo%name_tr (nspec_all) )
    allocate( leo%unit_tr (nspec_all) )
    allocate( leo%unitconv(nspec_all) )
    allocate( leo%varid_tr(nspec_all) )

    ! match tracer names:
    call goMatchValues( tracer_names, specname, &
                          leo%ntr, leo%name_tr, leo%itr, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected tracers for drydep output:")'); call goPr
    do i = 1, leo%ntr
      itr = leo%itr(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3, a10," ",a10,")")') &
                  i, leo%name_tr(i), &
                  itr, trim(specname(itr)), trim(specunit(itr)); call goPr
    end do

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

    ! full domain
    leo%i1 = 1
    leo%i2 = nx
    leo%ni = nx
    leo%j1 = 1
    leo%j2 = ny
    leo%nj = ny

    ! bounding box
    call rcF%Get( trim(basekey)//'.bounding_box', key, status )
    IF_NOTOK_RETURN(status=1)

    ! empty?
    if (len_trim(key) == 0) then
      ! full domain
      leo%bb_i1 = 1
      leo%bb_i2 = glb_ugg%nlon
      leo%bb_ni = glb_ugg%nlon
      leo%bb_j1 = 1
      leo%bb_j2 = glb_ugg%nlat
      leo%bb_nj = glb_ugg%nlat      
      leo%westb = glb_ugg%longitude_bnds(1,1,1)
      leo%southb = glb_ugg%latitude_bnds(1,1,1)
    else
      
      select case ( trim(glb_ugg%type) ) 
        
        case ( 'cartesian-regular') 
          ! read domain from key
          read(key,*,iostat=status) west, east, south, north
          if(status/=0) then
            write (gol,'("could not read domain from key: ",a)') trim(key); call goErr
            TRACEBACK; status=1; return
          endif

          ! Check if bounding box is in run domain
          if ( west < glb_ugg%longitude_bnds_1d(1,1) .or. east > glb_ugg%longitude_bnds_1d(2,glb_ugg%nlon) .or. &
               south < glb_ugg%latitude_bnds_1d(1,1) .or. north > glb_ugg%latitude_bnds_1d(2,glb_ugg%nlat) ) then
            write( gol, '("Bounding box domain is (partly) outside run domain")' ) ; call goErr
            write( gol, '("Run domain: ", 4f8.2)' ) glb_ugg%longitude_bnds_1d(1,1),glb_ugg%longitude_bnds_1d(2,glb_ugg%nlon),glb_ugg%latitude_bnds_1d(1,1),glb_ugg%latitude_bnds_1d(2,glb_ugg%nlat); call goErr
            write( gol, '("Bounding Box domain: ", 4f8.2)' ) west, east, south, north ; call goErr
            TRACEBACK;status=1;return
          endif

          ! for safety
          nullify(ff)
          ! get cell range covered by box
          call glb_ugg%GetDistribution(west,east,south,north,leo%bb_i1,leo%bb_i2,leo%bb_j1,leo%bb_j2,ff,status)
          IF_NOTOK_RETURN(status=1)
          !clear, fractions not used
          if ( associated(ff) ) deallocate(ff)
          ! set shape
          leo%bb_ni = leo%bb_i2-leo%bb_i1+1
          leo%bb_nj = leo%bb_j2-leo%bb_j1+1
          ! set west/south bounds
          leo%westb  = glb_ugg%longitude_bnds_1d(1,leo%bb_i1)
          leo%southb = glb_ugg%latitude_bnds_1d(1,leo%bb_j1)
        case default 
          write( gol, '("Definition of bounding box not clear for grid-type: ", a)' ) trim(glb_ugg%type) ; call goErr
          TRACEBACK;status=1;return
      end select
    end if

    ! enable data:
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'hp', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Column_Init


  ! ***


  subroutine LE_Output_Column_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif

    ! --- in/out --------------------------------

    type(T_LE_Output_Column), intent(inout)   ::  leo
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Column_Done'

    character(len=256) :: commandline
    ! --- begin ---------------------------------

#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif

    ! file opened ?
    if ( leo%opened ) then
      ! close:
#ifdef with_netcdf
      status = NF90_Close( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)      
#endif
      ! reset flag:
      leo%opened = .false.
      ! write GrADS ctl file if necessary:
      call Write_GrADS_Ctl( leo, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! clear storage for tracer fields:
    deallocate( leo%itr      )
    deallocate( leo%name_tr  )
    deallocate( leo%unit_tr  )
    deallocate( leo%unitconv )
    deallocate( leo%varid_tr )

    ! done with common stuff ...
    call leo%com%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Column_Done


  ! ***


  subroutine LE_Output_Column_PutOut( leo, t, c, status )

    use Binas  , only : xm_air, grav
    use LE_Comm, only : goc_dom
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>), rTotal, dTotal, iTotal
    use GO     , only : wrtgol, Precisely, MidNight
    use LE_Grid, only : glb_ugg
    use C3PO   , only : T_Grid_NcDef

#ifdef with_netcdf
    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF , only : NF90_NETCDF4
    use NetCDF , only : NF90_GLOBAL
    use NetCDF , only : NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT
#endif

    use LE_CF_Conventions, only : LE_CF_names

    use Dims            , only : nx, ny, nz, nspec
    ! point to meteo data:
    use LE_Data         , only : LE_Data_GetPointer
    use Indices         , only : specname, specunit, specmolm
    use Indices         , only : accum_n, accum_ii, accum_ww, accum_ppb_to_ugm3
    use LE_Bound        , only : caloft
    use LE_Output_Tools , only : LE_Output_Define_Dims_Time
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains

    ! --- in/out --------------------------------

    type(T_LE_Output_Column), intent(inout)   ::  leo
    type(TDate), intent(in)                   ::  t
    real, intent(in)                          ::  c(nx,ny,nz,nspec)
    integer, intent(out)                      ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Column_PutOut'

    ! --- local ---------------------------------

    integer               ::  time6(6)
    integer               ::  time
#ifdef with_netcdf
    integer               ::  cmode
#endif
    integer               ::  itr, iz
    integer               ::  varid
    type(T_Grid_NcDef)    ::  gncd
    type(TDate)           ::  t0
    integer               ::  i, j
    integer               ::  l, la
    real, allocatable     ::  pat(:,:)
    real, allocatable     ::  col(:,:)
    real, allocatable     ::  dh(:,:)
    real, allocatable     ::  dm(:,:)
    real, allocatable     ::  dm_sum(:,:)
    real, allocatable     ::  convfact(:,:)

    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units, cf_tracer_name
    character(len=512)    ::  comment
    character(len=32)     ::  varname
    character(len=32)     ::  afac
    
    character(len=256) :: commandline

    ! meteo data:
    real, pointer          ::  h_m(:,:,:)   ! (lon,lat,1)
    real, pointer          ::  hp(:,:,:)    ! (lon,lat,0:nlev)
    real, pointer          ::  dens(:,:,:)   ! (lon,lat,alt)    

    ! --- begin ----------------------------
    
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m')    
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOTOK_RETURN(status=1)   
    call LE_Data_GetPointer( 'hp', hp, status, check_units ='Pa')    
    IF_NOTOK_RETURN(status=1)

#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif

    ! for multiples of dhour only ...
    if ( .not. Precisely(t,leo%dhour,'hour')  ) then
      status=0; return
    end if
    
    ! dummy ...
    cf_units = 'no-cf-units'

    ! current time not in time range ?
    if ( (t < leo%tr(1)) .or. (leo%tr(2) < t) ) then

      ! extract time fields:
      call Get( t, time6=time6 )

      ! daily or less ?
      select case ( trim(leo%collect) )
        ! collect daily for [00,24)
        case ( 'daily' )
          ! set time range [00,24) for this day:
          leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
          leo%tr(2) = leo%tr(1) + IncrDate( day=1 ) - IncrDate(mili=1)
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), time6(1:3)
          if ( len_trim(leo%state) > 0 ) write (leo%fname,'(a,"_",a)') trim(leo%fname), trim(leo%state)
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! collect daily for (00,24]
        case ( 'daily24' )
          ! set time range (00,24] for this day:
          if ( Midnight(t) ) then
            leo%tr(1) = t - IncrDate(day=1)
            leo%tr(2) = t
          else
            leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
            leo%tr(2) = leo%tr(1) + IncrDate( day=1 )
            leo%tr(1) = leo%tr(1) + IncrDate(mili=1)
          end if
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), leo%tr(1)%year, leo%tr(1)%month, leo%tr(1)%day
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! files with instant fields:
        case ( 'instant' )
          ! set time range for current instant time:
          leo%tr(1) = t
          leo%tr(2) = t
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), time6(1:5)
          if ( len_trim(leo%state) > 0 ) write (leo%fname,'(a,"_",a)') trim(leo%fname), trim(leo%state)
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! unknonw ...
        case default
          write (gol,'("unsupported collect key : ",a)') trim(leo%collect); call goErr
          TRACEBACK; status=1; return
      end select

      ! root only:
      if ( goc_dom%root ) then

#ifdef with_netcdf
        ! set creation mode flag:
        if ( leo%replace ) then
          cmode = NF90_CLOBBER       ! overwrite existing files
        else
          cmode = NF90_NOCLOBBER     ! do not overwrite existing files
        end if

        ! enable large file support:
        cmode = or( cmode, NF90_NETCDF4 )

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
        call leo%com%PutOut_GlobalAttributes( leo%ncid, status )
        IF_NOTOK_RETURN(status=1)

#ifdef with_netcdf

        ! grid dimensions/variables
        call glb_ugg%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                    dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat,  &
                                    subset=(/leo%bb_i1,leo%bb_i2,leo%bb_j1,leo%bb_j2/) )
        IF_NOTOK_RETURN(status=1)

        ! time dimensions
        call LE_Output_Define_Dims_Time(leo%ncid, leo%dimid_time, status)
        IF_NOTOK_RETURN(status=1)

        ! time since t0
        t0 = leo%com%t0      
        ! time variables
        call LE_Output_Define_Vars_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                        leo%dimid_time, trim(leo%com%CF_convention), t0, status)
        IF_NOTOK_RETURN(status=1)

#endif
  
      end if ! root

      ! loop over tracers to be written:
      do l = 1, leo%ntr

        ! global tracer index
        itr = leo%itr(l)

        ! CF standard name for concentration/mixing ratio/column:

        ! no comment yet
        comment = ''

        ! set names etc following some convention:
        select case ( trim(leo%convention) )
        
          case ( 'cf' )

            ! standard variable name:
            varname = trim(leo%name_tr(l))//'_column'

            ! get names following CF conventions;
            ! store conversion factor for later usage:
            call LE_CF_names( &
                         specname(itr), specunit(itr), &
                         cf_standard_name, cf_long_name, cf_units, &
                         leo%unitconv(l), comment, &
                         status, cf_tracer_name=cf_tracer_name )
            IF_NOTOK_RETURN(status=1)
            
            ! reset to new name and unit ...
            cf_standard_name = 'column_integrated_'//trim(cf_tracer_name)
            cf_long_name     = trim(cf_standard_name)
            cf_units = 'ug m-2'
            ! reset to dummy value, not used:
            leo%unitconv(l) = -999.9
        
          case ( 'x' )

            ! standard variable name:
            varname = 'x'//trim(leo%name_tr(l))

            ! get names following CF conventions;
            ! store conversion factor for later usage:
            call LE_CF_names( &
                         specname(itr), specunit(itr), &
                         cf_standard_name, cf_long_name, cf_units, &
                         leo%unitconv(l), comment, &
                         status, cf_tracer_name=cf_tracer_name )
            IF_NOTOK_RETURN(status=1)
            
            ! reset gasses to new name and unit ...
            if ( trim(specunit(itr)) == 'ppb' ) then
              ! column mixing ratio's:
              cf_standard_name = ''
              cf_long_name     = 'column mixing ratio of '//trim(cf_tracer_name)
              cf_units         = 'ppb'
            else
              ! reset to new name and unit ...
              cf_standard_name = ''
              cf_long_name     = 'column_integrated_'//trim(cf_tracer_name)
              cf_units         = 'ug m-2'
            end if
            ! reset to dummy value, not used:
            leo%unitconv(l)  = -999.9
            
          case ( 'megapoli' )

            ! standard variable name:
            varname = trim(leo%name_tr(l))//'_column'

            ! get names following CF conventions;
            ! store conversion factor for later usage:
            call LE_CF_names( &
                         specname(itr), specunit(itr), &
                         cf_standard_name, cf_long_name, cf_units, &
                         leo%unitconv(l), comment, &
                         status, cf_tracer_name=cf_tracer_name )
            IF_NOTOK_RETURN(status=1)
            
            ! reset to new name and unit ...
            cf_standard_name = 'column_integrated_'//trim(cf_tracer_name)
            cf_long_name     = trim(cf_standard_name)
            cf_units = 'ug m-2'
            ! reset to dummy value, not used:
            leo%unitconv(l) = -999.9
            
          case default
          
            write (gol,'("unsupported output convention : ",a)') trim(leo%convention); call goErr
            TRACEBACK; status=1; return
            
        end select
        
        ! accumulated specs ?
        if ( accum_n(itr) > 1 ) then
          if ( len_trim(comment) > 0 ) comment = trim(comment)//' ; '
          do la = 1, accum_n(itr)
            if ( la > 1 ) comment = trim(comment)//' + '
            if ( accum_ppb_to_ugm3(itr,la) ) then
              write (afac,'(f5.2)') accum_ww(itr,la)/xm_air
              comment = trim(comment)//trim(specname(accum_ii(itr,la)))//'*'//trim(afac)//'*dens'
            else
              write (afac,'(f5.2)') accum_ww(itr,la)
              comment = trim(comment)//trim(specname(accum_ii(itr,la)))//'*'//trim(afac)
            end if
          end do
        end if

        ! store units for later usage (GrADS ctl file, conversions):
        leo%unit_tr(l) = trim(cf_units)

        ! root only:
        if ( goc_dom%root ) then

#ifdef with_netcdf
          ! define variable:
          status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                   (/leo%dimid_lon,leo%dimid_lat,leo%dimid_time/), varid, &
                                   deflate_level=leo%com%deflate_level )
          if (status/=0) then
            write (gol,'("defining variable : ",a)') trim(leo%name_tr(l)); call goErr
            TRACEBACK; status=1; return
          end if

          ! write attributes:
          if ( len_trim(cf_standard_name) > 0 ) then
            status = nf90_put_att( leo%ncid, varid, 'standard_name', trim(cf_standard_name) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
          status = nf90_put_att( leo%ncid, varid, 'long_name', trim(cf_long_name) )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, varid, 'units', trim(cf_units) )
          IF_NF90_NOTOK_RETURN(status=1)
          call glb_ugg%DefCoor_NetCDF( gncd, varid, status )
          IF_NF90_NOTOK_RETURN(status=1)
          ! molemass:
          status = nf90_put_att( leo%ncid, varid, 'molemass', specmolm(itr) )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, varid, 'molemass_unit', 'kg mole-1' )
          IF_NF90_NOTOK_RETURN(status=1)
          ! optional comment:
          if ( len_trim(comment) > 0 ) then
            status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
#endif

          ! store variable id:
          leo%varid_tr(l) = varid
          
        end if  ! root

      end do  ! written tracers

      ! root only:
      if ( goc_dom%root ) then

        ! end defintion mode:
#ifdef with_netcdf
        status = NF90_EndDef( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif

      end if  ! root

      ! no records written yet:
      leo%itrec = 0

    end if

    ! next time record:
    leo%itrec = leo%itrec + 1

    ! root only:
    if ( goc_dom%root ) then
    
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
    allocate( col(1:nx,1:ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( dh(1:nx,1:ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( convfact(1:nx,1:ny), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! air masses:
    allocate( dm(1:nx,1:ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( dm_sum(1:nx,1:ny), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over all tracer to be written:
    do l = 1, leo%ntr
      ! global tracer index:
      itr = leo%itr(l)
      
      ! init column sum:
      col = 0.0
      ! init air mass sum:
      dm_sum = 0.0
      
      ! loop over layers:
      do iz = 1, leo%nlev

        ! init 2d field to zero:
        pat = 0.0

        ! add contributions:
        do la = 1, accum_n(itr)
          ! conversion ?
          convfact = 1.0
          if ( accum_ppb_to_ugm3(itr,la) ) convfact = dens(1:nx,1:ny,iz)/xm_air
          ! add:
          if ( iz <= nz ) then
            pat(1:nx,1:ny) = pat(1:nx,1:ny) + c     (1:nx,1:ny,iz,accum_ii(itr,la)) * convfact * accum_ww(itr,la)
          else
            pat(1:nx,1:ny) = pat(1:nx,1:ny) + caloft(1:nx,1:ny,iz,accum_ii(itr,la)) * convfact * accum_ww(itr,la)
          end if
        end do

        ! layer heigh in m:
        if ( iz == 1 ) then
          dh = h_m(:,:,iz)  ! m
        else
          dh = h_m(:,:,iz) - h_m(:,:,iz-1)   ! m
        end if
        
        ! air mass in layer:
        !        kg/m3             m
        dm = dens(1:nx,1:ny,iz) * dh  ! (kg air)/layer
        ! add contribution:
        dm_sum = dm_sum + dm

        ! specials ...
        select case ( trim(specunit(itr))//' -> '//trim(leo%unit_tr(l)) )
          !~ other notations ...
          case ( 'ug/m3 -> ug m-2' )
            ! integerate over height:
            pat = pat * dh
          !~ volume mixing ratio to mass concentration:
          case ( 'ppb -> ug m-2' )
            ! ug/m3        =      ppb  m    * ((mole tr/mole air)/ppb) * (kg tr/mole tr) *   (kg air/m3 air)  / (kg air/mole air) * ug/kg
            pat(1:nx,1:ny) = pat(1:nx,1:ny) *       1e-9               *  specmolm(itr)  * dens(1:nx,1:ny,iz) /     xm_air        * 1e9
            ! integrate over height:
            pat = pat * dh
          !~ column mixing ratio:
          case ( 'ppb -> ppb' )
            ! mass weighted:
            pat = pat * dm  ! ppb*kg
          !~ unkown ...
          case default
            write (gol,'("do not know how to convert `",a,"` from `",a,"` to `",a,"`")') &
                    trim(specname(itr)), trim(specunit(itr)), trim(leo%unit_tr(l)); call goErr
            TRACEBACK; status=1; return
        end select
        
        ! take fractions if below ptop:
        do j = 1, ny
          do i = 1, nx
            ! switch:
            if ( hp(i,j,iz-1) < leo%ptop ) then
              ! bottom above ptop:
              pat(i,j) = 0.0
            else if ( hp(i,j,iz) < leo%ptop ) then
              ! take fraction:
              pat(i,j) = pat(i,j) * (hp(i,j,iz-1) - leo%ptop)/(hp(i,j,iz-1) - hp(i,j,iz))
            end if
          end do ! i
        end do ! j

        ! add contribution:
        col = col + pat
        
      end do  ! layers
      
      ! postprocess?
      select case ( trim(specunit(itr))//' -> '//trim(leo%unit_tr(l)) )
        !~ no postprocessing needed:
        case ( 'ug/m3 -> ug m-2', 'ppb -> ug m-2' )
          ! skip
        !~ column mixing ratio:
        case ( 'ppb -> ppb' )
          ! mass weighted:
          !     ppb*kg /   kg
          col =   col  / dm_sum
        !~ unkown ...
        case default
          write (gol,'("do not know how to convert `",a,"` from `",a,"` to `",a,"`")') &
                  trim(specname(itr)), trim(specunit(itr)), trim(leo%unit_tr(l)); call goErr
          TRACEBACK; status=1; return
      end select

      ! write 2D field, without level, with time index:
      call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_tr(l), -1, leo%itrec, &
                                          col(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                          bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
      IF_NOTOK_RETURN(status=1)
    
    end do  ! tracers

    ! clear:
    deallocate( pat, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( col, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( dh, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( dm, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( dm_sum, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( convfact, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! root?
    if ( goc_dom%root ) then

      ! next time exceeds interval ?
      if ( t+IncrDate(hour=int(leo%dhour)) > leo%tr(2) ) then
        ! close:
#ifdef with_netcdf
        status = NF90_Close( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif
        ! reset flag:
        leo%opened = .false.
        ! write GrADS ctl file if necessary:
        call Write_GrADS_Ctl( leo, status )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
      
    end if  ! root

    ! ok
    status = 0

  end subroutine LE_Output_Column_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use LE_Grid, only : glb_ugg

    ! --- in/out ---------------------------------

    type(T_LE_Output_Column), intent(inout)   ::  leo
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

      ! daily or less ?
      select case ( trim(leo%collect) )
        ! collect daily for [00,24)
        case ( 'daily', 'daily24' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_base,'(a,"_",a)') trim(leo%grads_ctl_base), trim(leo%state)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! files with instant fields:
        case ( 'instant' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2_%h2%n2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_base,'(a,"_",a)') trim(leo%grads_ctl_base), trim(leo%state)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! unknonw ...
        case default
          write (gol,'("unsupported collect key : ",a)') trim(leo%collect); call goErr
          TRACEBACK; status=1; return
      end select

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
      call GrADS_Ctl_Vars( ctl, leo%ntr, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over tracers to be written:
      do l = 1, leo%ntr
        ! set variable description:
        write (line,'(a,"_column [",a,"]")') trim(leo%name_tr(l)), trim(leo%unit_tr(l))
        ! add variable description:
        call GrADS_Ctl_Var( ctl, trim(leo%name_tr(l))//'_column', 0, 't,y,x', trim(line), status )
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


end module LE_Output_Column
