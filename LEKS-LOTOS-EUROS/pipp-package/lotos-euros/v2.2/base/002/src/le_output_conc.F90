!###############################################################################
!
! NAME
!
!   LE_Output_Conc  -  LOTOS-EUROS output of concentration fields
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

module LE_Output_Conc

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_Conc

  public  ::  LE_Output_Conc_Init, LE_Output_Conc_Done
  public  ::  LE_Output_Conc_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Conc'

  ! maximum number of levels in output:
  integer, parameter  ::  maxlev = 20


  ! --- types ------------------------------

  type T_LE_Output_Conc
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
    integer                     ::  dimid_lev
    integer                     ::  dimid_hlev
    integer                     ::  dimid_time
    ! dimension variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_lev
    integer                     ::  varid_time
    !integer                     ::  varid_time_day
    integer                     ::  varid_time_dtg
    integer                     ::  varid_dh
    integer                     ::  varid_height
    ! extra meteo variables for 3d-output (temperature, humidity, half-level-pres, half-level-alt)
    integer                     ::  varid_T
    integer                     ::  varid_Q
    integer                     ::  varid_hp
    integer                     ::  varid_halt
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
    integer, pointer            ::  ilev(:)
    real                        ::  heights(maxlev)
    ! heights distribution (interp or sample)
    character(len=16)           ::  height_dist
    ! put out halo cells ?
    logical                     ::  halo
    integer                     ::  nhx, nhy, nhz
    ! put out height/layer thickness field too ?
    logical                     ::  addheight
    logical                     ::  addthickness
    ! grads ctl file ?
    logical                     ::  grads_ctl
    character(len=256)          ::  grads_ctl_file
    character(len=256)          ::  grads_ctl_base
    integer                     ::  grads_ctl_nt
    type(TDate)                 ::  grads_ctl_t1
    type(TIncrDate)             ::  grads_ctl_dt
    ! indices    
    integer                     ::  i1, i2, ni
    integer                     ::  j1, j2, nj
    !bounding box
    integer                     ::  bb_i1, bb_i2, bb_ni
    integer                     ::  bb_j1, bb_j2, bb_nj
    real                        ::  westb, southb

  end type T_LE_Output_Conc


contains


  ! ====================================================


  subroutine LE_Output_Conc_Init( leo, rcF, rckey, typ, name, state, status )

    use GO     , only : TrcFile
    use GO     , only : goMatchValues, goSplitString, goReadFromLine
    use GO     , only : AnyDate
    !use GO     , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use Dims   , only : nx, ny, nz
    
    use Indices         , only : nspec, nspec_all
    use Indices         , only : specname, specunit
    use LE_Grid         , only : glb_ugg
    use LE_Data         , only : LE_Data_Enable
    use LE_Data         , only : LE_Data_GetPointer
    use LE_Data_Common  , only : nlev, mixlayer_top
    use LE_Data_Common  , only : nlev_top

    ! --- in/out --------------------------------

    type(T_LE_Output_Conc), intent(out)   ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  state
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Conc_Init'

    ! --- local ---------------------------------

    character(len=64)     ::  basekey
    character(len=1024)   ::  tracer_names
    integer               ::  itr
    character(len=256)    ::  level_names
    integer               ::  ilev
    integer               ::  i, i1, i2

    character(len=512)    ::  key
    real                  ::  west, south, east, north
    real, pointer         ::  ff(:,:)

    ! --- begin ---------------------------------

    !! start timing:
    !call GO_Timer_Start( itim_init_output_conc, status )
    !IF_NOTOK_RETURN(status=1)

    ! store name:
    leo%typ  = typ
    leo%name = name

    !! start timing:
    !call GO_Timer_Start( itim_init_output_conc_com, status )
    !IF_NOTOK_RETURN(status=1)

    ! init common stuff:
    call leo%com%Init( rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)

    !! end timing:
    !call GO_Timer_End( itim_init_output_conc_com, status )
    !IF_NOTOK_RETURN(status=1)

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

    ! collect daily or instant
    call rcF%Get( trim(basekey)//'.collect', leo%collect, status )
    IF_NOTOK_RETURN(status=1)

    ! output time resolution:
    call rcF%Get( trim(basekey)//'.dhour', leo%dhour, status )
    IF_NOTOK_RETURN(status=1)

    ! tracer names:
    call rcF%Get( trim(basekey)//'.fields', tracer_names, status )
    IF_NOTOK_RETURN(status=1)

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
    write (gol,'("selected tracers for conc output:")'); call goPr
    do i = 1, leo%ntr
      itr = leo%itr(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3, a10," ",a10,")")') &
                  i, leo%name_tr(i), &
                  itr, trim(specname(itr)), trim(specunit(itr)); call goPr
    end do
    write (gol,'("(transported tracers up to nr. ",i3,")")') nspec; call goPr

    ! enable data:
    call LE_Data_Enable( 'oro', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    
    ! extra meteo output
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'q', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'hp', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'halt', status )
    IF_NOTOK_RETURN(status=1)

    ! level type:
    call rcF%Get( trim(basekey)//'.levtype', leo%levtype, status )
    IF_NOTOK_RETURN(status=1)
    ! level descriptions:
    call rcF%Get( trim(basekey)//'.levels', level_names, status )
    IF_NOTOK_RETURN(status=1)
    ! layers or heights ?
    select case ( trim(leo%levtype) )
      !~ model levels:
      case ( 'levels' )
        ! setup storage for level indices (surface + levels + upper boundary)
        !allocate( leo%ilev(1+nz+1) )
        allocate( leo%ilev(1:100) )
        ! switch:
        if ( trim(level_names) == 'all' ) then
          ! all:
          leo%nlev = nz
          do i = 1, leo%nlev
            leo%ilev(i) = i
          end do
        else if ( trim(level_names) == 'top' ) then
          ! all:
          leo%nlev = nlev_top
          do i = 1, leo%nlev
            leo%ilev(i) = i
          end do
        else if ( index(trim(level_names),':') > 0 ) then
          ! extract range:
          call goReadFromLine( level_names, i1, status, sep=':' )
          IF_NOTOK_RETURN(status=1)
          call goReadFromLine( level_names, i2, status )
          IF_NOTOK_RETURN(status=1)
          ! count:
          leo%nlev = i2 - i1 + 1
          ! store indices:
          do i = 1, leo%nlev
            leo%ilev(i) = i1 - 1 + i
          end do
        else
          ! set selected level indices:
          call goMatchValues( level_names, 0, nlev_top, leo%nlev, leo%ilev, status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! info ...
        write (gol,'("selected levels for conc output:")'); call goPr
        do i = 1, leo%nlev
          ilev = leo%ilev(i)
          write (gol,'("  ",i3,"  ",i3)') i, ilev; call goPr
        end do
      !~ height levels:
      case ( 'heights', 'elevations' )
        ! extract heights (m):
        call goSplitString( trim(level_names), leo%nlev, leo%heights, status )
        IF_NOTOK_RETURN(status=1)
        ! need to store 3D boundary fields ?
        ! compare with level definition:
        if ( any( leo%heights(1:leo%nlev) > mixlayer_top(nlev) ) ) then
          ! requested height above model
          write( gol, '("Requested output height above model top")' ) ; call goErr
          write( gol, '("Model top: ", f7.4)' ) mixlayer_top(nlev) ; call goErr
          write( gol, '("Maximum requested height: ", f7.4)' ) maxval(leo%heights(1:leo%nlev)) ; call goErr
          TRACEBACK;status=1;return
        end if
        ! dummy ...
        allocate( leo%ilev(leo%nlev) )
        leo%ilev = -999
        ! read height distibution (interp or sample)
        call rcF%Get( trim(basekey)//'.height_dist', leo%height_dist, status)
        IF_NOTOK_RETURN(status=1)
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! write halo cells ?
    call rcF%Get( trim(basekey)//'.halo', leo%halo, status )
    IF_ERROR_RETURN(status=1)
    if ( leo%halo ) then
      ! not yet ...
      write (gol,'("output of halo not supported yet in domain decomposition ...")'); call goErr
      TRACEBACK; status=1; return
      ! for future:
      leo%nhx = 1
      leo%nhy = 1
    else
      leo%nhx = 0
      leo%nhy = 0
    end if
    write (gol,'("halo cells in x and y direction : ",2i3)') leo%nhx, leo%nhy; call goPr

    ! write height field to conc file ?
    call rcF%Get( trim(basekey)//'.addheight', leo%addheight, status )
    IF_ERROR_RETURN(status=1)

    ! write layer thickness field to conc file ?
    call rcF%Get( trim(basekey)//'.addthickness', leo%addthickness, status )
    IF_ERROR_RETURN(status=1)

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
      ! check HALO not used
      if (leo%halo ) then
        write (gol,'("used halo and bounding box together, this is not possible")'); call goErr
        TRACEBACK; status=1; return
      endif
      
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

    !! end timing:
    !call GO_Timer_End( itim_init_output_conc, status )
    !IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Conc_Init


  ! ***


  subroutine LE_Output_Conc_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif

    ! --- in/out --------------------------------

    type(T_LE_Output_Conc), intent(inout)   ::  leo
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Conc_Done'

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

    ! layers or heights ?
    select case ( trim(leo%levtype) )
      !~ model levels:
      case ( 'levels' )
        ! clear storage for level indices:
        deallocate( leo%ilev )
      !~ height levels:
      case ( 'heights', 'elevations' )
        ! clear storage for level indices:
        deallocate( leo%ilev )
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! done with common stuff ...
    call leo%com%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Conc_Done


  ! ***


  subroutine LE_Output_Conc_PutOut( leo, t, c, cg, status )

    use Binas  , only : xm_air
    use GO     , only : goc
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>), rTotal, dTotal, iTotal
    use GO     , only : wrtgol, Precisely
    use Num    , only : Interp_Lin
    use Num    , only : Interval
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

    use LE_CF_Conventions   , only : LE_CF_names

    use Dims            , only : nx, ny, nz, nspec
    use LE_Data_Common  , only : nlev_top
    use Indices         , only : specname, specunit, specmolm
    use Indices         , only : accum_n, accum_ii, accum_ww, accum_ppb_to_ugm3
    use LE_Data         , only : LE_Data_GetPointer
    use LE_Bound        , only : bc_west, bc_east, bc_south, bc_north, caloft
    use LE_Output_Tools , only : LE_Output_Define_Dims_Lev
    use LE_Output_Tools , only : LE_Output_Define_Dims_Time

    use LE_Output_Tools , only : LE_Output_Define_Vars_Lev
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time

    use LE_Output_Tools , only : LE_Output_Put_Var_Lev
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
    
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains

    ! --- in/out --------------------------------

    type(T_LE_Output_Conc), intent(inout)   ::  leo
    type(TDate), intent(in)                 ::  t
    real, intent(in)                        ::  c (nx,ny,nz,nspec)
    real, intent(in)                        ::  cg(nx,ny,nspec)
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Conc_PutOut'

    ! --- local ---------------------------------

    integer               ::  time6(6)
    integer               ::  time
#ifdef with_netcdf
    integer               ::  cmode
#endif
    integer               ::  itr, ilev, iz
    integer               ::  varid
    type(T_Grid_NcDef)    ::  gncd
    type(TDate)           ::  t0
    integer               ::  i, j, k, l, la
    integer               ::  ilayer
    real, allocatable     ::  field(:,:,:)  ! nx,ny,0:nlev_top
    real, allocatable     ::  pat(:,:)
    real, allocatable     ::  convfact(:,:)  ! (nx,ny)

    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units, cf_tracer_name
    character(len=512)    ::  comment
    character(len=32)     ::  varname
    character(len=32)     ::  afac

    real                  ::  cg_o3(nx,ny)

    real, allocatable     ::  hh(:,:,:)
    real, allocatable     ::  hsamp(:,:)
    integer               ::  ilast
    character(len=256)    ::  bt_units

    character(len=256) :: commandline
    
    ! meteo data:
    real, pointer        ::  oro(:,:,:)    ! (lon,lat,1)
    real, pointer        ::  dens(:,:,:)   ! (lon,lat,lev)    
    real, pointer        ::  h_m(:,:,:)    ! (lon,lat,lev)
  
    real, pointer        ::  temp(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::  Q(:,:,:)      ! (lon,lat,lev)
    real, pointer        ::  hp(:,:,:)     ! (lon,lat,hlev)
    real, pointer        ::  halt(:,:,:)   ! (lon,lat,hlev)

    ! --- begin ---------------------------------

#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif
    
    ! point to meteo data:
    call LE_Data_GetPointer( 'oro', oro, status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)

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
      if ( goc%root ) then

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

        ! level dimension
        call LE_Output_Define_Dims_Lev(leo%ncid, leo%dimid_lev, leo%nlev, trim(leo%com%CF_convention), status, leo%dimid_hlev)
        IF_NOTOK_RETURN(status=1)

        ! time dimensions
        call LE_Output_Define_Dims_Time(leo%ncid, leo%dimid_time, status)
        IF_NOTOK_RETURN(status=1)

        ! level variables
        call LE_Output_Define_Vars_Lev(leo%ncid, leo%varid_lev, leo%dimid_lev, &
                                         trim(leo%levtype), trim(leo%com%CF_convention), status)
        IF_NOTOK_RETURN(status=1)

        ! time since t0
        t0 = leo%com%t0      
        ! time variables
        call LE_Output_Define_Vars_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                        leo%dimid_time, trim(leo%com%CF_convention), t0, status)
        IF_NOTOK_RETURN(status=1)


        ! add layer heights as 4D field ?
        if ( leo%addheight ) then
           !! info ...
           !write (gol,'("LE_Output_Conc: add height field definition ")'); call goPr
           ! open:
           status = NF90_Def_Var( leo%ncid, 'altitude', NF90_REAL, &
                                  (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), varid, &
                                  deflate_level=leo%com%deflate_level )
           IF_NF90_NOTOK_RETURN(status=1)
           status = nf90_put_att( leo%ncid, varid, 'standard_name', 'altitude' )
           IF_NF90_NOTOK_RETURN(status=1)
           status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'layer height relative to sea level' )
           IF_NF90_NOTOK_RETURN(status=1)
           status = NF90_Put_Att( leo%ncid, varid, 'units', 'm' )
           IF_NF90_NOTOK_RETURN(status=1)
           leo%varid_height = varid
        end if

        ! add layer thickness as 4D field ?
        if ( leo%addthickness ) then
           ! info ...
           write (gol,'("LE_Output_Conc: add layer thickness field definition ")'); call goPr
           ! open:
           status = NF90_Def_Var( leo%ncid, 'cell_thickness', NF90_REAL, &
                                  (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), varid, &
                                  deflate_level=leo%com%deflate_level )
           IF_NF90_NOTOK_RETURN(status=1)
           status = nf90_put_att( leo%ncid, varid, 'standard_name', 'cell_thickness' )
           IF_NF90_NOTOK_RETURN(status=1)
           status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'vertical extent of the model grid cell' )
           IF_NF90_NOTOK_RETURN(status=1)
           status = NF90_Put_Att( leo%ncid, varid, 'units', 'm' )
           IF_NF90_NOTOK_RETURN(status=1)
           leo%varid_dh = varid
        end if

        !status = NF90_Def_Var( leo%ncid, 'time_day', NF90_REAL, leo%dimid_time, varid )
        !IF_NF90_NOTOK_RETURN(status=1)
        !status = NF90_Put_Att( leo%ncid, varid, 'standard_name', 'time' )
        !IF_NF90_NOTOK_RETURN(status=1)
        !status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'time' )
        !IF_NF90_NOTOK_RETURN(status=1)
        !status = NF90_Put_Att( leo%ncid, varid, 'description', 'yyyymmdd+dayfraction' )
        !IF_NF90_NOTOK_RETURN(status=1)
        !status = NF90_Put_Att( leo%ncid, varid, 'calendar', 'standard' )
        !IF_NF90_NOTOK_RETURN(status=1)
        !status = NF90_Put_Att( leo%ncid, varid, '_CoordinateAxisType', 'Time' )
        !IF_NF90_NOTOK_RETURN(status=1)
        !leo%varid_time_day = varid

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
            varname = trim(leo%name_tr(l))

            ! get names following CF conventions;
            ! store conversion factor for later usage:
            call LE_CF_names( &
                         specname(itr), specunit(itr), &
                         cf_standard_name, cf_long_name, cf_units, &
                         leo%unitconv(l), comment, &
                         status )
            IF_NOTOK_RETURN(status=1)

          case ( 'megapoli' )

            ! standard variable name:
            varname = trim(leo%name_tr(l))

            ! get names following CF conventions;
            ! store conversion factor for later usage:
            call LE_CF_names( &
                         specname(itr), specunit(itr), &
                         cf_standard_name, cf_long_name, cf_units, &
                         leo%unitconv(l), comment, &
                         status, cf_tracer_name=cf_tracer_name )
            IF_NOTOK_RETURN(status=1)

            ! reset to new unit ...
            cf_standard_name = 'mass_concentration_of_'//trim(cf_tracer_name)//'_in_air'
            cf_long_name     = trim(cf_standard_name)
            cf_units = 'ug m-3'
            leo%unitconv(l) = -999.9    ! this enables conversion from ppb to ug/m3 below ...

          case default

            write (gol,'("unsupported output convention : ",a)') trim(leo%convention); call goErr
            TRACEBACK; status=1; return

        end select

        ! store units for later usage (GrADS ctl file, conversions):
        leo%unit_tr(l) = trim(cf_units)

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
        
        ! root?
        if ( goc%root ) then

#ifdef with_netcdf
          ! define variable:
          status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                 (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), varid, &
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
          IF_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time level latitude longitude')
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
      
      ! define extra variables for 3d -output
      if ( goc%root ) then
        if ( leo%nlev > 1 ) then
          ! temperature:
          status = nf90_Def_Var( leo%ncid, 'T', NF90_REAL, &
                                 (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), leo%varid_T, &
                                 deflate_level=leo%com%deflate_level )
          if ( status/= 0) then
            write( gol,'("defining variable Temperature: ")' ) ; call goErr
            TRACEBACK; status=1; return
          end if 
          status = nf90_put_att( leo%ncid, leo%varid_T, 'long_name', 'Temperature' )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, leo%varid_T, 'units', 'K' )
          IF_NF90_NOTOK_RETURN(status=1)
          call glb_ugg%DefCoor_NetCDF( gncd, leo%varid_T, status )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, leo%varid_T, '_CoordinateAxes', 'time level latitude longitude')
          IF_NF90_NOTOK_RETURN(status=1)
          ! specific humidity:
          status = nf90_Def_Var( leo%ncid, 'Q', NF90_REAL, &
                                 (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), leo%varid_Q, &
                                 deflate_level=leo%com%deflate_level )
          if ( status/= 0) then
            write( gol,'("defining variable Specific humidity: ")' ) ; call goErr
            TRACEBACK; status=1; return
          end if 
          status = nf90_put_att( leo%ncid, leo%varid_Q, 'long_name', 'Specific humidity' )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, leo%varid_Q, 'units', 'kg/kg' )
          IF_NF90_NOTOK_RETURN(status=1)
          call glb_ugg%DefCoor_NetCDF( gncd, leo%varid_Q, status )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, leo%varid_Q, '_CoordinateAxes', 'time level latitude longitude')
          IF_NF90_NOTOK_RETURN(status=1)
          ! Half level pressures:
          status = nf90_Def_Var( leo%ncid, 'hp', NF90_REAL, &
                                 (/leo%dimid_lon,leo%dimid_lat,leo%dimid_hlev,leo%dimid_time/), leo%varid_hp, &
                                 deflate_level=leo%com%deflate_level )
          if ( status/= 0) then
            write( gol,'("defining variable half level pressure: ")' ) ; call goErr
            TRACEBACK; status=1; return
          end if 
          status = nf90_put_att( leo%ncid, leo%varid_hp, 'long_name', 'Half-level-pressures' )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, leo%varid_hp, 'units', 'Pa' )
          IF_NF90_NOTOK_RETURN(status=1)
          call glb_ugg%DefCoor_NetCDF( gncd, leo%varid_hp, status )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, leo%varid_hp, '_CoordinateAxes', 'time hlevel latitude longitude')
          IF_NF90_NOTOK_RETURN(status=1)
          ! Half level altitudes:
          status = nf90_Def_Var( leo%ncid, 'halt', NF90_REAL, &
                                 (/leo%dimid_lon,leo%dimid_lat,leo%dimid_hlev,leo%dimid_time/), leo%varid_halt, &
                                 deflate_level=leo%com%deflate_level )
          if ( status/= 0) then
            write( gol,'("defining variable half level altitude: ")' ) ; call goErr
            TRACEBACK; status=1; return
          end if 
          status = nf90_put_att( leo%ncid, leo%varid_halt, 'long_name', 'Half-level-altitudes' )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, leo%varid_halt, 'units', 'm' )
          IF_NF90_NOTOK_RETURN(status=1)
          call glb_ugg%DefCoor_NetCDF( gncd, leo%varid_halt, status )
          IF_NF90_NOTOK_RETURN(status=1)
          status = nf90_put_att( leo%ncid, leo%varid_halt, '_CoordinateAxes', 'time hlevel latitude longitude')
          IF_NF90_NOTOK_RETURN(status=1)
        end if  ! nlev > 1
      end if  ! root
      ! end defintion mode:

      ! root?
      if ( goc%root ) then
        ! end defintion mode:
#ifdef with_netcdf
        status = NF90_EndDef( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif
      end if ! root

      ! no records written yet:
      leo%itrec = 0

    end if

    ! next time record:
    leo%itrec = leo%itrec + 1

    ! root ?
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

#ifdef with_netcdf
        ! write level indices
        call LE_Output_Put_Var_Lev(leo%ncid, leo%varid_lev, leo%nlev, &
                                     trim(leo%levtype), leo%ilev, leo%heights, status)
        IF_NOTOK_RETURN(status=1)                                   
#endif

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

    ! 3d field from surface to aloft:
    allocate( field(leo%i1-leo%nhx:leo%i2+leo%nhx,leo%j1-leo%nhy:leo%j2+leo%nhy,0:nlev_top) )
    ! 2d field:
    allocate( pat(leo%i1-leo%nhx:leo%i2+leo%nhx,leo%j1-leo%nhy:leo%j2+leo%nhy) )
    ! 3d conversion field:
    allocate( convfact(nx,ny) )

    ! add layer heights as 4D field ?
    if ( leo%addheight ) then
      ! info ...
      write (gol,'("LE_Output_Conc: add height field values ")'); call goPr
      ! loop over levels:
      do k = 1, leo%nlev
        ! global level index:
        ilev = leo%ilev(k)
        if ( ilev <= 0 ) then
           pat(leo%i1:leo%i2,leo%j1:leo%j2) = 2.0  ! surface level at 2 m
        else if ( ilev <= nz ) then
           pat(leo%i1:leo%i2,leo%j1:leo%j2) = h_m(leo%i1:leo%i2,leo%j1:leo%j2,ilev) ! m
        else
           pat = 2.0e5  ! 200 km
        end if

        ! add orography to model height
        pat(leo%i1:leo%i2,leo%j1:leo%j2) = pat(leo%i1:leo%i2,leo%j1:leo%j2) + oro(leo%i1:leo%i2,leo%j1:leo%j2,1) 
        ! fill halo cells if necessary:
        if ( leo%halo ) then
           do i = 1, leo%nhx
              pat( 1-i,1:ny) = pat( 1,1:ny)
              pat(nx+i,1:ny) = pat(nx,1:ny)
           end do
           do j = 1, leo%nhy
              pat(1:nx, 1-j) = pat(1:nx, 1)
              pat(1:nx,ny+j) = pat(1:nx,ny)
           end do
        end if

        !status = NF90_Put_Var( leo%ncid, leo%varid_height, pat, &
        !     start=(/1,1,k,leo%itrec/), count=(/leo%ni+2*leo%nhx,leo%nj+2*leo%nhy,1,1/) )
        !IF_NF90_NOTOK_RETURN(status=1)

        ! write 2D field with level and time index:
        call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_height, k, leo%itrec, &
                                          pat(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                          bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
        IF_NOTOK_RETURN(status=1)

      end do
    end if

    ! add layer thickness as 4D field ?
    if ( leo%addthickness ) then
      ! info ...
      write (gol,'("LE_Output_Conc: add layer thickness field values ")'); call goPr
      do k = 1, leo%nlev
        ! layer thickness:
        if ( k == 1 ) then
          pat = h_m(:,:,k) ! m
        else
          pat = h_m(:,:,k) - h_m(:,:,k-1) ! m
        end if
        ! write 2d field:
!#ifdef with_netcdf
!        status = NF90_Put_Var( leo%ncid, leo%varid_dh, pat, &
!                       start=(/1,1,k,leo%itrec/), count=(/nx+2*leo%nhx,ny+2*leo%nhy,1,1/) )
!        IF_NF90_NOTOK_RETURN(status=1)
!#endif
        ! write 2D field with level and time index:
        call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dh, k, leo%itrec, &
                                          pat(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                          bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! dummy storage to avoid warnings:
    allocate( hh(1,1,1), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( hsamp(1,1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! height samples ?
    if ( (trim(leo%levtype) == 'heights') .or. &
         (trim(leo%levtype) == 'elevations') ) then
      ! clear dummy storage:
      deallocate( hh, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! half level heights:
      allocate( hh(1-leo%nhx:nx+leo%nhx,1-leo%nhy:ny+leo%nhy,0:nlev_top), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! lowest is orography:
      hh(1:nx,1:ny,0) = oro(1:nx,1:ny,1)  ! m
      ! interp or samples
      if (trim(leo%height_dist) == 'interp') then
        ! mid of layers:
        hh(1:nx,1:ny,1) = oro(1:nx,1:ny,1) + h_m(1:nx,1:ny,1) * 0.5
        do k = 2, nlev_top
          hh(1:nx,1:ny,k) = oro(1:nx,1:ny,1) + ( h_m(1:nx,1:ny,k-1) + h_m(1:nx,1:ny,k) ) * 0.5
        end do
      else if (trim(leo%height_dist) == 'sample') then
        ! top of layers
        do k = 1, nlev_top
          hh(1:nx,1:ny,k) = oro(1:nx,1:ny,1) + h_m(1:nx,1:ny,k)  ! m
        end do
      else
        write (gol,'(" unsupported height distribution type : ",a)') trim(leo%height_dist); call goErr
        write (gol,'(" type must be either interp or sample")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! halo cells:
      if ( leo%halo ) then
        do i = 1, leo%nhx
          hh( 1-i,:,:) = hh( 1,:,:)
          hh(nx+i,:,:) = hh(nx,:,:)
        end do
        do j = 1, leo%nhy
          hh(:, 1-j,:) = hh(:, 1,:)
          hh(:,ny+j,:) = hh(:,ny,:)
        end do
      end if
      ! clear dummy storage:
      deallocate( hsamp, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! sample height:
      allocate( hsamp(1-leo%nhx:nx+leo%nhx,1-leo%nhy:ny+leo%nhy), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! loop over all tracer to be written:
    do l = 1, leo%ntr

      ! global tracer index:
      itr = leo%itr(l)

      ! loop over all model levels, including surface and aloft:
      do ilev = 0, nlev_top

        ! init 2d field to zero:
        pat = 0.0

        ! extract 2d field:
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if ( ilev == 0 ) then  ! -- extract surface field:
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add contributions:
          do la = 1, accum_n(itr)
            ! conversion needed ?
            convfact = 1.0
            if ( accum_ppb_to_ugm3(itr,la) ) convfact = dens(:,:,1)/xm_air  ! (mole air)/m3
            ! core:
            pat(leo%i1:leo%i2,leo%j1:leo%j2) = pat(leo%i1:leo%i2,leo%j1:leo%j2) &
                         + cg(leo%i1:leo%i2,leo%j1:leo%j2,accum_ii(itr,la)) &
                                * convfact(leo%i1:leo%i2,leo%j1:leo%j2) * accum_ww(itr,la)
            ! halo:
            if ( leo%halo ) then
              do i = 1, leo%nhx
                pat( 1-i,1:ny) = pat( 1-i,1:ny) + bc_west (:,1,accum_ii(itr,la)) * convfact( 1,:) * accum_ww(itr,la)
                pat(nx+i,1:ny) = pat(nx+i,1:ny) + bc_east (:,1,accum_ii(itr,la)) * convfact(nx,:) * accum_ww(itr,la)
              end do
              do j = 1, leo%nhy
                pat(1:nx, 1-j) = pat(1:nx, 1-j) + bc_south(:,1,accum_ii(itr,la)) * convfact(:, 1) * accum_ww(itr,la)
                pat(1:nx,ny+j) = pat(1:nx,ny+j) + bc_north(:,1,accum_ii(itr,la)) * convfact(:,ny) * accum_ww(itr,la)
              end do
            end if
          end do

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else if ( ilev > nz ) then ! -- extract from upper boundary:
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add contributions:
          do la = 1, accum_n(itr)
            ! conversion needed ?
            convfact = 1.0
            if ( accum_ppb_to_ugm3(itr,la) ) convfact = dens(:,:,ilev)/xm_air
            ! core:
            pat(leo%i1:leo%i2,leo%j1:leo%j2) = pat(leo%i1:leo%i2,leo%j1:leo%j2) &
                    + caloft(leo%i1:leo%i2,leo%j1:leo%j2,ilev,accum_ii(itr,la)) &
                       * convfact(leo%i1:leo%i2,leo%j1:leo%j2) * accum_ww(itr,la)
          end do

          ! fill halo cells if necessary;
          ! no boundary conditions available for this layer,
          ! so copy from core:
          if ( leo%halo ) then
            do i = 1, leo%nhx
              pat( 1-i,1:ny) = pat( 1,1:ny)
              pat(nx+i,1:ny) = pat(nx,1:ny)
            end do
            do j = 1, leo%nhy
              pat(1:nx, 1-j) = pat(1:nx, 1)
              pat(1:nx,ny+j) = pat(1:nx,ny)
            end do
          end if

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else ! -- extract level:
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add contributions:
          do la = 1, accum_n(itr)
            ! conversion needed ?
            convfact = 1.0
            if ( accum_ppb_to_ugm3(itr,la) ) convfact = dens(:,:,ilev)/xm_air
            ! add 2D field:
            pat(leo%i1:leo%i2,leo%j1:leo%j2) = pat(leo%i1:leo%i2,leo%j1:leo%j2) &
                         + c(leo%i1:leo%i2,leo%j1:leo%j2,ilev,accum_ii(itr,la)) &
                                * convfact(leo%i1:leo%i2,leo%j1:leo%j2) * accum_ww(itr,la)
            ! update halo?
            if ( leo%halo ) then
              do i = 1, leo%nhx
                pat( 1-i,1:ny) = pat( 1-i,1:ny) + bc_west (:,ilev,accum_ii(itr,la)) * convfact( 1,:) * accum_ww(itr,la)
                pat(nx+i,1:ny) = pat(nx+i,1:ny) + bc_east (:,ilev,accum_ii(itr,la)) * convfact(nx,:) * accum_ww(itr,la)
              end do
              do j = 1, leo%nhy
                pat(1:nx, 1-j) = pat(1:nx, 1-j) + bc_south(:,ilev,accum_ii(itr,la)) * convfact(:, 1) * accum_ww(itr,la)
                pat(1:nx,ny+j) = pat(1:nx,ny+j) + bc_north(:,ilev,accum_ii(itr,la)) * convfact(:,ny) * accum_ww(itr,la)
              end do
            end if
          end do

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end if ! level
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! unit conversion:
        if ( leo%unitconv(l) > 0.0 ) then
          ! single factor ...
          pat = pat * leo%unitconv(l)
        else
          ! specials ...
          select case ( trim(specunit(itr))//' -> '//trim(leo%unit_tr(l)) )
            !~ other notations ...
            case ( 'ug/m3 -> ug m-3' )
              ! no further conversion
            !~ volume mixing ratio to mass concentration:
            case ( 'ppb -> ug m-3' )
              ! 1-nz, index in density array
              iz = min( max( 1, ilev ), nz )
              ! ug/m3 =  ppb * ((mole tr/mole air)/ppb) *  (kg tr/mole tr) *             (kg air/m3 air)         /  (kg air/mole air) * ug/kg
              pat     = pat  *       1e-9               *   specmolm(itr)  * dens(leo%i1:leo%i2,leo%j1:leo%j2,iz) /     xm_air        * 1e9
              if ( leo%halo ) then
                stop 'please implement unit convesion for halo cells'
              end if
            !~ unkown ...
            case default
              write (gol,'("do not know how to convert `",a,"` from `",a,"` to `",a,"`")') &
                      trim(specname(itr)), trim(specunit(itr)), trim(leo%unit_tr(l)); call goErr
              TRACEBACK; status=1; return
          end select
        end if

        ! fill halo corners if necessry:
        if ( leo%halo ) then
          do i = 1, leo%nhx
            do j = 1, leo%nhy
              pat( 1-i, 1-j) = ( pat( 1, 1-j) + pat( 1-i, 1) )/2.0
              pat(nx+i, 1-j) = ( pat(nx, 1-j) + pat(nx+i, 1) )/2.0
              pat( 1-i,ny+j) = ( pat( 1-i,ny) + pat( 1,ny+j) )/2.0
              pat(nx+i,ny+j) = ( pat(nx+i,ny) + pat(nx,ny+j) )/2.0
            end do
          end do
        end if

        ! store:
        field(:,:,ilev) = pat

      end do   ! model layers

      ! * vertical

      ! which output levels ?
      select case ( trim(leo%levtype) )

        !~ selected model levels:
        case ( 'levels' )

          ! loop over selected layer:
          do k = 1, leo%nlev
            ! global level index:
            ilev = leo%ilev(k)
            ! extract 2D field:
            pat = field(:,:,ilev)
!            ! write concentrations:
!#ifdef with_netcdf
!            status = NF90_Put_Var( leo%ncid, leo%varid_tr(l), pat, &
!                           start=(/1,1,k,leo%itrec/), count=(/leo%ni+2*leo%nhx,leo%nj+2*leo%nhy,1,1/) )
!            IF_NF90_NOTOK_RETURN(status=1)
!#endif
            ! write 2D field with level and time index:
            call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_tr(l), k, leo%itrec, &
                                            pat(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                            bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
            IF_NOTOK_RETURN(status=1)
          end do

        !~ height levels:
        case ( 'heights', 'elevations' )

          ! dummy:
          ilast = 1
          ! loop over target levels:
          do k = 1, leo%nlev

            ! sample height:
            if ( trim(leo%levtype) == 'heights' ) then
              ! relative to orography:
              hsamp = hh(:,:,0) + leo%heights(k)
            else if ( trim(leo%levtype) == 'elevations' ) then
              ! absolute, interpolation will take surface value if below orography:
              hsamp = leo%heights(k)
            else
              write (gol,'("no sample height defined for level type : ",a)') trim(leo%levtype); call goErr
              TRACEBACK; status=1; return
            end if

            ! height distribution (interp or samples)
            if (trim(leo%height_dist) == 'interp') then

              ! loop over horizontal cells:            
              do j = leo%j1-leo%nhy, leo%j2+leo%nhy
                do i = leo%i1-leo%nhx, leo%i2+leo%nhx
                  ! below or above model top ?
                  if ( hsamp(i,j) <= hh(i,j,nz) ) then
                    ! vertical interpolation of model layers;
                    ! below surface is taken from lowest value:
                    call Interp_Lin( hh(i,j,:), field(i,j,:), max(hh(i,j,0),hsamp(i,j)), &
                                         pat(i,j), ilast, status )
                    IF_NOTOK_RETURN(status=1)
                  else
                    write( gol, '("interpolation height above model top!")' ) ; call GoErr
                    TRACEBACK;status=1;return
                  end if
                end do  ! i
              end do  ! j

            else if (trim(leo%height_dist) == 'sample') then

              ! loop over horizontal cells:            
              do j = leo%j1-leo%nhy, leo%j2+leo%nhy
                do i = leo%i1-leo%nhx, leo%i2+leo%nhx
                  ! below or above model top ?
                  if ( hsamp(i,j) <= hh(i,j,nz) ) then
                    ! get layer index;
                    ! note that return value is in 1,..,nz+2
                    ! return status -1 if sample below lowest value (ok, sample cg)
                    call Interval( hh(i,j,:), hsamp(i,j), ilayer, status )
                    IF_ERROR_RETURN(status=1)
                    ! reset to 0,..,nlev_top
                    ilayer = ilayer -1
                    ! extract:
                    pat(i,j) = field(i,j,ilayer)
                  else
                    write( gol, '("sampling height above model top!")' ) ; call GoErr
                    TRACEBACK;status=1;return                    
                  end if
                end do
              end do

            else
              write (gol,'(" unsupported height distribution type : ",a)') trim(leo%height_dist); call goErr
              write (gol,'(" type must be either interp or sample")'); call goErr
              TRACEBACK; status=1; return
            end if

            ! write 2D field with level and time index:
            call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_tr(l), k, leo%itrec, &
                                            pat(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                            bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
            IF_NOTOK_RETURN(status=1)

          end do  ! heights

        !~ unknown ...
        case default
          write (gol,'("unsupported level type : ",a)') trim(leo%levtype); call goErr
          TRACEBACK; status=1; return

      end select

    end do   ! tracers
    
    ! Add extra meteo variables for 3d-output
    if ( leo%nlev > 1 ) then
      ! temperature
      call LE_Data_GetPointer( 't', temp, status, check_units='K')
      IF_NOTOK_RETURN(status=1)
      ! specific humidity
      call LE_Data_GetPointer( 'q', Q, status, check_units='kg/kg')
      IF_NOTOK_RETURN(status=1)
      ! half-level-pressure
      call LE_Data_GetPointer( 'hp', hp, status, check_units='Pa')
      IF_NOTOK_RETURN(status=1)
      ! half-level-altitude
      call LE_Data_GetPointer( 'halt', halt, status, check_units='m')
      IF_NOTOK_RETURN(status=1)

      ! loop over selected layer (full levels):
      do k = 1, leo%nlev
        ! global level index:
        ilev = leo%ilev(k)
        ! extract 2D field:
        pat = temp(:,:,ilev)
        ! write concentrations:
        call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_T, k, leo%itrec, &
                                        pat(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                        bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
        IF_NOTOK_RETURN(status=1)
        
        ! extract 2D field:
        pat = Q(:,:,ilev)
        ! write concentrations:
        call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_Q, k, leo%itrec, &
                                        pat(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                        bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
        IF_NOTOK_RETURN(status=1)
      end do
      
      ! loop over selected layer (half levels):
      do k = 0, leo%nlev
        ! global level index:
        if ( k == 0 ) then
          ilev = 0
        else 
          ilev = leo%ilev(k)
        end if
        ! extract 2D field:
        pat = hp(:,:,ilev)
        ! write half-level pressures:
        call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_hp, k+1, leo%itrec, &
                                        pat(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                        bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
        IF_NOTOK_RETURN(status=1)
        ! extract 2D field:
        pat = halt(:,:,ilev)
        ! write  half-level altitudes:
        call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_halt, k+1, leo%itrec, &
                                        pat(leo%i1:leo%i2,leo%j1:leo%j2), status, &
                                        bb_start=(/leo%bb_i1,leo%bb_j1/), bb_count=(/leo%bb_ni,leo%bb_nj/) )
        IF_NOTOK_RETURN(status=1)
      end do
      
    end if  ! (levo%nlev > 0)  
      
    
    ! clear:
    deallocate( pat )
    deallocate( field )
    deallocate( convfact )
    if ( trim(leo%levtype) == 'heights' ) then
      deallocate( hh )
      deallocate( hsamp )
    end if
    
    ! root?
    if ( goc%root ) then
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
        IF_NOTOK_RETURN(status=1)
      end if
    end if ! root

    ! ok
    status = 0

  end subroutine LE_Output_Conc_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use LE_Grid, only : glb_ugg
    
    ! --- in/out ---------------------------------

    type(T_LE_Output_Conc), intent(inout)   ::  leo
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
        case ( 'daily' )
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
      call glb_ugg%WriteCtlProjection( ctl, status )
      IF_NOTOK_RETURN(status=1)
      ! zdef:
      select case ( trim(leo%levtype) )
        case ( 'levels' )
          call GrADS_Ctl_ZDef( ctl, leo%ilev(1:leo%nlev), status )
          IF_NOTOK_RETURN(status=1)
        case ( 'heights', 'elevations' )
          call GrADS_Ctl_ZDef( ctl, nint(leo%heights(1:leo%nlev)), status )
          IF_NOTOK_RETURN(status=1)
        case default
          write (gol,'("do not know how to write zdef for levtyp : ",a)') trim(leo%levtype); call goPr
          TRACEBACK; status=1; return
      end select
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
        write (line,'(a," [",a,"]")') trim(leo%name_tr(l)), trim(leo%unit_tr(l))
        ! add variable description:
        call GrADS_Ctl_Var( ctl, trim(leo%name_tr(l)), leo%nlev, 't,z,y,x', trim(line), status )
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


end module LE_Output_Conc
