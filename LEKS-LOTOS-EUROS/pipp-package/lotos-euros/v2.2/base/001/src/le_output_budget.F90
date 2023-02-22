!###############################################################################
!
! NAME
!
!   LE_Output_Budget  -  LOTOS-EUROS output of deposition/fluxes
!
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

module LE_Output_Budget

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_Budget

  public  ::  LE_Output_Budget_Init, LE_Output_Budget_Done
  public  ::  LE_Output_Budget_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Budget'


  ! maximum number of supported data sets:
  integer, parameter  ::  ndep = 40


  ! --- types ------------------------------

  type T_LE_Output_Budget
    ! name for this file:
    character(len=16)           ::  typ
    character(len=16)           ::  name
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
    ! replace existing files ?
    logical                     ::  replace
    ! file opened ?
    logical                     ::  opened
    ! current time range:
    type(TDate)                 ::  tr(2)
    ! time resolution:
    ! time record counter:
    integer                     ::  itrec
    ! file name:
    character(len=256)          ::  fname
    ! file handle:
    integer                     ::  ncid
    ! dimension handles:
    integer                     ::  dimid_lon
    integer                     ::  dimid_lat
    integer                     ::  dimid_lev
    integer                     ::  dimid_time
    ! dimension variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_lev
    integer                     ::  varid_time
    integer                     ::  varid_time_dtg
    integer                     ::  varid_time_day
    ! tracer variables:
    integer                     ::  ndep
    integer                     ::  ndep_total
    integer, pointer            ::  idat(:)
    character(len=10), pointer   ::  name_dat(:)
    character(len=32), pointer  ::  unit_dat(:)
    real, pointer               ::  unitconv(:)
    integer, pointer            ::  varid_dat(:,:)
    ! landuse variables:
    integer                     ::  nlu
    integer, pointer            ::  ilu(:)
    character(len=32), pointer  ::  name_lu(:)

    ! level selection:
    integer                     ::  nlev

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

  end type T_LE_Output_Budget

  type T_LE_Dat
    character(len=32)      ::  name
    character(len=32)      ::  unit
    integer                ::  rank
    logical                ::  const
    logical                ::  landuse_dependent
  end type T_LE_Dat


  ! --- var --------------------------------------

  type(T_LE_Dat)     ::  le_dat(ndep)


contains


  ! ====================================================


  subroutine SetDat( d, name, unit, rank, const, landuse_dependent )

    ! --- in/out ----------------------------------

    type(T_LE_Dat), intent(out)     ::  d
    character(len=*), intent(in)    ::  name
    character(len=*), intent(in)    ::  unit
    integer, intent(in)             ::  rank
    logical, intent(in)             ::  const
    logical, intent(in)             ::  landuse_dependent

    ! --- begin ----------------------------------

    d%name = name
    d%unit = unit
    d%rank = rank
    d%const = const
    d%landuse_dependent = landuse_dependent

  end subroutine SetDat


  ! ====================================================


  subroutine LE_Output_Budget_Init( leo, rcF, rckey, typ, name, status )

    use GO     , only : TrcFile
    use GO     , only : goMatchValues
    use GO     , only : AnyDate
    use Dims   , only : nx, ny
    use Indices, only : ispec_noya, ispec_soxa, ispec_nhxa
    use Indices, only : ispec_o3, ispec_nh3, ispec_so2
    use Indices, only : specname
    use LE_Grid, only : ugg
    use LE_Landuse_Data , only : nlu, lu_name, lu_name_abbr
    use LE_Landuse_Data , only : ludep_output_whole_grid
    use LE_Data, only : LE_Data_Enable

    use LE_Output_Common, only : Init

    ! --- in/out --------------------------------

    type(T_LE_Output_Budget), intent(out) ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Budget_Init'

    ! rckey extensions:

    ! --- local ---------------------------------

    character(len=32)     ::  basekey
    character(len=1024)   ::  selected_names
    character(len=1024)   ::  field_names
    integer               ::  idat,ilu
    integer               ::  i
    integer               ::  ispec
    character(len=8)      ::  sname

    character(len=512)    ::  key
    real                  ::  west, south, east, north
    real, pointer         ::  ff(:,:)
    ! --- begin ---------------------------------

    ! store name:
    leo%typ  = typ
    leo%name = name

    ! init common stuff:
    call Init( leo%com, rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)

    ! replace existing files ?
    call rcF%Get( trim(rckey)//'.replace', leo%replace, status )
    IF_NOTOK_RETURN(status=1)

    ! write GrADS ctl file ?
    call rcF%Get( trim(rckey)//'.ctl', leo%grads_ctl, status )
    IF_NOTOK_RETURN(status=1)

    ! base key:
    write (basekey,'(a,".",a,".",a)') trim(rckey), trim(typ), trim(name)

    ! tracer names:
    call rcF%Get( trim(basekey)//'.fields', field_names, status )
    IF_NOTOK_RETURN(status=1)


    ! define supported data:
    !   name   :  used in rcfile to identify this field;
    !   unit   :  in which the corresponding field in LE is expressed;
    !             the model units are converted to 'cf'-conventions
    !             (SI-units) on output;
    !   rank   :  2 for 2D field, and 3 for 3D
    !   const  :  logical, T for constant fields
    !                        name        unit       rank  const  , landuse 
    call SetDat( LE_Dat( 1), 'drysox', 'ug m-2'   ,   2, .false. , .false. )
    call SetDat( LE_Dat( 2), 'drynoy', 'ug m-2'   ,   2, .false. , .false. )
    call SetDat( LE_Dat( 3), 'drynhx', 'ug m-2'   ,   2, .false. , .false. )
    call SetDat( LE_Dat( 4), 'wetsox', 'ug m-2'   ,   2, .false. , .false. )
    call SetDat( LE_Dat( 5), 'wetnoy', 'ug m-2'   ,   2, .false. , .false. )
    call SetDat( LE_Dat( 6), 'wetnhx', 'ug m-2'   ,   2, .false. , .false. )
    call SetDat( LE_Dat( 7), 'o3max' , 'ppb'      ,   2, .false. , .false. )
    call SetDat( LE_Dat( 8), 'nh3ave', 'ppb'      ,   2, .false. , .false. )
    call SetDat( LE_Dat( 9), 'so2ave', 'ppb'      ,   2, .false. , .false. )
    call SetDat( LE_Dat(10), 'o3f'   , 'umol m-2' ,   2, .false. , .true. )
    call SetDat( LE_Dat(11), 'o3e'   , 'umol m-2' ,   2, .false. , .true. )
    call SetDat( LE_Dat(12), 'drysox_lu', 'ug m-2',   2, .false. , .true. )
    call SetDat( LE_Dat(13), 'drynoy_lu', 'ug m-2',   2, .false. , .true. )
    call SetDat( LE_Dat(14), 'drynhx_lu', 'ug m-2',   2, .false. , .true. )
    
    ! setup storage for tracer fields:
    allocate( leo%idat     (ndep) )
    allocate( leo%name_dat (ndep) )
    allocate( leo%unit_dat (ndep) )
    allocate( leo%unitconv (ndep) )

    ! match tracer names:
    call goMatchValues( field_names, LE_Dat(:)%name, &
                          leo%ndep, leo%name_dat, leo%idat, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected fields for data output:")'); call goPr
    do i = 1, leo%ndep
      idat = leo%idat(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3, a10," ",a10,")")') &
                  i, leo%name_dat(i), &
                  idat, trim(LE_Dat(idat)%name), trim(LE_Dat(idat)%unit); call goPr
    end do

    ! info ...
    write (gol,'("check if required tracers are enabled ...")'); call goPr
    ! safety check, arrays always exist but might have remained zero
    ! if the accumulated species were not enabled ;
    ! loop over all variables to be written:
    do i = 1, leo%ndep
      ! set (accumulated) spec index used for this output:
      select case ( trim(leo%name_dat(i)) )
        case ( 'drysox_lu', 'drysox', 'wetsox' ) ; ispec = ispec_soxa ; sname = 'SOx'
        case ( 'drynoy_lu', 'drynoy', 'wetnoy' ) ; ispec = ispec_noya ; sname = 'NOy'
        case ( 'drynhx_lu', 'drynhx', 'wetnhx' ) ; ispec = ispec_nhxa ; sname = 'NHx'
        case ( 'o3max', 'o3f', 'o3e'           ) ; ispec = ispec_o3   ; sname = 'O3'
        case ( 'nh3ave'                        ) ; ispec = ispec_nh3  ; sname = 'NH3'
        case ( 'so2ave'                        ) ; ispec = ispec_so2  ; sname = 'SO2'
        case default
          write (gol,'("no check on species implemented for budget output `",a,"`")') &
                   trim(leo%name_dat(i)); call goErr
          TRACEBACK; status=1; return
      end select
      ! check ...
      if ( ispec < 0 ) then
        write (gol,'("budget output `",a,"` needs spec `",a,"`")') &
            trim(leo%name_dat(i)), trim(sname); call goErr
        TRACEBACK; status=1; return
      else
        write (gol,'("  budget `",a,"` requires `",a,"` : ok")') &
            trim(leo%name_dat(i)), trim(sname); call goPr
      end if
    end do  ! budget outputs
    ! landuse names:
    call rcF%Get( trim(basekey)//'.landuses', selected_names, status )
    IF_NOTOK_RETURN(status=1)
    ! setup storage for tracer fields:
    allocate( leo%ilu     (nlu) )
    allocate( leo%name_lu (nlu) )
    ! match tracer names:
    call goMatchValues( selected_names, lu_name_abbr, &
                          leo%nlu, leo%name_lu, leo%ilu, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected landuses for budget output:")'); call goPr
    do i = 1, leo%nlu
      ilu = leo%ilu(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3," ",a,")")') &
                  i, leo%name_lu(i), ilu, trim(lu_name(ilu)); call goPr
      ludep_output_whole_grid(ilu) = .true.                  
    end do

    ! storage for written variables:
    allocate( leo%varid_dat(leo%ndep,leo%nlu) )

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
    call rcF%Get( trim(basekey)//'.bounding_box', key, status )
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
      ! not yet ...
      write (gol,'("no output subset supported for domain decomposition yet")'); call goErr
      TRACEBACK; status=1; return
      
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

    ! enable data:
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Budget_Init


  ! ***


  subroutine LE_Output_Budget_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif
    use LE_Output_Common, only : Done

    ! --- in/out --------------------------------

    type(T_LE_Output_Budget), intent(inout)   ::  leo
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Budget_Done'

    character(len=256) :: commandline
    ! --- begin ---------------------------------

    ! file opened ?
    if ( leo%opened ) then
      ! close:
#ifdef with_netcdf
      status = NF90_Close( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#else
    stop 'not compiled with netcdf support'
#endif
      ! reset flag:
      leo%opened = .false.
      ! write GrADS ctl file if necessary:
      call Write_GrADS_Ctl( leo, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! clear storage for tracer fields:
    deallocate( leo%idat      )
    deallocate( leo%name_dat  )
    deallocate( leo%unit_dat  )
    deallocate( leo%unitconv  )
    deallocate( leo%varid_dat )



    ! done with common stuff ...
    call Done( leo%com, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Budget_Done


  ! ***


  subroutine LE_Output_Budget_PutOut( leo, t, bud, status )

    use Binas  , only : xm_air
    use GO     , only : goc
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>), rTotal, dTotal, iTotal
    use GO     , only : wrtgol, Precisely, MidNight
    use LE_Grid, only : ugg, glb_ugg
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

    use Dims   , only : nx, ny, nz
    use Indices, only : accum_ii, accum_ww, accum_n, accum_ppb_to_ugm3
    use Indices, only : ispec_noya, ispec_soxa, ispec_nhxa
    !use indices, only: !i_o3 ! i_so4a, i_so2, i_no2, i_no, i_no3a, i_nh3, i_nh4a, i_hno3
    !use LE_GrOutput, only : bud_depos, cdry, o3max, bud_stomata, bud_stomata_exd !,cwet
    use LE_Budget         , only : T_Budget
    use LE_Budget_Param   , only : ibud_sox, ibud_noy, ibud_nhx
    use LE_Budget_DryDepos, only : ex_drydepo
    use LE_Budget_WetDepos, only : ex_wetdepo

    use LE_Output_Common, only : PutOut_GlobalAttributes

    use LE_Output_Tools , only : LE_Output_Define_Dims_Time
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
        
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains
    
    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out --------------------------------

    type(T_LE_Output_Budget), intent(inout)   ::  leo
    type(TDate), intent(in)                   ::  t
    type(T_Budget), intent(in)                ::  bud
    integer, intent(out)                      ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Budget_PutOut'

    ! --- local ---------------------------------

    integer               ::  time6(6)
    integer               ::  time
#ifdef with_netcdf
    integer               ::  cmode
#endif
    integer               ::  idat, ilu
    type(T_Grid_NcDef)    ::  gncd
    integer               ::  varid
    type(TDate)           ::  t0, tdum
    integer               ::  i, j, l, n
    integer               ::  ibud_spec
    real                  ::  pat(nx,ny)
    real                  ::  convfact(nx,ny)

    character(len=64)     ::  varname
    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units
    character(len=512)    ::  comment

    character(len=256) :: commandline

    ! meteo data:
    real, pointer          ::  dens(:,:,:)   ! (lon,lat,alt)    

    ! --- begin ----------------------------
    
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOTOK_RETURN(status=1)

#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif

    ! extract time fields:
    tdum=t-IncrDate(day=1)
    call Get( tdum, time6=time6 )

    ! set time range [00,24) for this day:
    leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3),hour=00 )

    !leo%tr(2) = leo%tr(1) + IncrDate( day=1 ) - IncrDate(mili=1)
    ! new file name:
    write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
              trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
              trim(leo%name), time6(1:3)
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
      ! dummmy ...
      leo%ncid = -1
#endif

      ! reset flag:
      leo%opened = .true.

      ! write global attributes:
      call PutOut_GlobalAttributes( leo%com, leo%ncid, status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_netcdf
      
      ! grid dimensions/variables
      call glb_ugg%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                  dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )
                                  !subset=(/leo%i1,leo%i2,leo%j1,leo%j2/) )
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
      leo%ndep_total = 0
    end if  ! root

    ! loop over data fields to be written:
    do l = 1, leo%ndep

      ! global tracer index
      idat = leo%idat(l)

      ! loop over landuses to be written
      do n = 1, leo%nlu
        ! not landuse dependent? skip after first
        if ( (.not. le_dat(idat)%landuse_dependent) .and. (n > 1) ) cycle
          ! landuse index
          ilu = leo%ilu(n)

          ! CF standard name for concentration/mixing ratio/column:

          ! initial comment:
          comment = ''

          ! variable name:
          varname = trim(leo%name_dat(l) )
          ! add landuse name in necessary:
          if ( le_dat(idat)%landuse_dependent ) then
            varname = trim(varname)//'_'//trim(leo%name_lu(n))
          end if

          ! get names following CF conventions;
          ! store conversion factor for later usage:
          call LE_CF_names( &
                       LE_Dat(idat)%name, LE_Dat(idat)%unit, &
                       cf_standard_name, cf_long_name, cf_units, &
                       leo%unitconv(l), comment, &
                       status )
          IF_NOTOK_RETURN(status=1)

          ! store units for later usage (GrADS ctl file):
          leo%unit_dat(l) = trim(cf_units)
          
          ! root?
          if ( goc%root ) then
            ! define variable:
#ifdef with_netcdf
            status = NF90_Def_Var( leo%ncid, varname, NF90_REAL, &
                               (/leo%dimid_lon,leo%dimid_lat,leo%dimid_time/), varid )
            IF_NF90_NOTOK_RETURN(status=1)
#endif

           ! write attributes:
#ifdef with_netcdf
            status = nf90_put_att( leo%ncid, varid, 'standard_name', trim(cf_standard_name) )
            IF_NF90_NOTOK_RETURN(status=1)
            status = nf90_put_att( leo%ncid, varid, 'long_name', trim(cf_long_name) )
            IF_NF90_NOTOK_RETURN(status=1)
            ! write units:
            status = nf90_put_att( leo%ncid, varid, 'units', trim(cf_units) )
            IF_NF90_NOTOK_RETURN(status=1)
            call glb_ugg%DefCoor_NetCDF( gncd, varid, status )
            IF_NOTOK_RETURN(status=1)
            ! add comment:
            if ( len_trim(comment) > 0 ) then
              status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
              IF_NF90_NOTOK_RETURN(status=1)
            end if
#endif
    
            ! store variable id:
            leo%varid_dat(l,n) = varid
          
            ! total number of data variables
            leo%ndep_total = leo%ndep_total + 1
          end if  ! root
          
        end do ! landuses
      end do  ! written tracers

      if ( goc%root ) then
#ifdef with_netcdf
        ! end defintion mode:
        status = NF90_EndDef( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
#endif

      ! no records written yet:
      leo%itrec = 0

      ! next time record:
      leo%itrec = leo%itrec + 1
    
    if ( goc%root ) then
      ! GrADS time counter:
      if ( leo%grads_ctl ) then
        ! increase counter:
        leo%grads_ctl_nt = leo%grads_ctl_nt + 1
        ! set times if necessary:
        if ( leo%grads_ctl_nt == 1 ) then
          if ( MidNight(t) ) then
            leo%grads_ctl_t1 = t - IncrDate(day=1)
          else
            leo%grads_ctl_t1 = t
          end if
          leo%grads_ctl_dt = IncrDate(day=1)   ! dummy ...
        end if
        if ( leo%grads_ctl_nt == 2 ) then
          if (MidNight(t) ) then
            leo%grads_ctl_dt = t - (leo%grads_ctl_t1 + IncrDate(day=1))
          else 
            leo%grads_ctl_dt = t - leo%grads_ctl_t1
          end if
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
    end if  ! root ?

    ! loop over all tracer to be written:
    do l = 1, leo%ndep
      ! global tracer index:
      idat = leo%idat(l)

      ! constant fields written only once:
      ! loop over landuses to be written:
      do n = 1, leo%nlu
        ! not landuse dependend ? then skip after first:
        if ( (.not. le_dat(idat)%landuse_dependent) .and. (n > 1) ) cycle
        ! landuse index:
        ilu = leo%ilu(n)

        ! extract 2d field:
        select case ( trim(LE_Dat(idat)%name) )

          case('drysox'); pat = bud%drydepos%ex_day(:,:,ibud_sox,ex_drydepo)
          case('drynoy'); pat = bud%drydepos%ex_day(:,:,ibud_noy,ex_drydepo)
          case('drynhx'); pat = bud%drydepos%ex_day(:,:,ibud_nhx,ex_drydepo)
          case('wetsox'); pat = bud%wetdepos%ex_day(:,:,ibud_sox,ex_wetdepo)
          case('wetnoy'); pat = bud%wetdepos%ex_day(:,:,ibud_noy,ex_wetdepo)
          case('wetnhx'); pat = bud%wetdepos%ex_day(:,:,ibud_nhx,ex_wetdepo)
          case('o3max' ); pat = bud%o3max
          case('nh3ave'); pat = bud%drydepos%cnh3_ave_prev
          case('so2ave'); pat = bud%drydepos%cso2_ave_prev
          case('o3f'   ); pat = bud%drydepos%o3flx_stomata_day_lu(:,:,ilu)
          case('o3e'   ); pat = bud%drydepos%o3exc_stomata_day_lu(:,:,ilu)
          case('drysox_lu')
            pat = 0.0
            ! loop over different components in SOx, with different weights
            do ibud_spec = 1, accum_n(ispec_soxa)
              convfact = 1.0
              if (accum_ppb_to_ugm3(ispec_soxa,ibud_spec) ) convfact = dens(:,:,1)/xm_air
              pat = pat + bud%drydepos%ex_day_lu(:,:,accum_ii(ispec_soxa,ibud_spec),ilu,ex_drydepo) * convfact * accum_ww(ispec_soxa,ibud_spec)
            end do            
          case('drynoy_lu')
            pat = 0.0
            ! loop over different components in NOy, with different weights
            do ibud_spec = 1, accum_n(ispec_noya)
              convfact = 1.0
              if (accum_ppb_to_ugm3(ispec_noya,ibud_spec) ) convfact = dens(:,:,1)/xm_air
              pat = pat + bud%drydepos%ex_day_lu(:,:,accum_ii(ispec_noya,ibud_spec),ilu,ex_drydepo) * convfact * accum_ww(ispec_noya,ibud_spec)
            end do              
          case('drynhx_lu')
            pat = 0.0
            ! loop over different components in NHx, with different weights
            do ibud_spec = 1, accum_n(ispec_nhxa) 
              convfact = 1.0
              if (accum_ppb_to_ugm3(ispec_nhxa,ibud_spec) ) convfact = dens(:,:,1)/xm_air
              pat = pat + bud%drydepos%ex_day_lu(:,:,accum_ii(ispec_nhxa,ibud_spec),ilu,ex_drydepo) * convfact * accum_ww(ispec_nhxa,ibud_spec)
            end do

        case default
            write (gol,'("unsupported 2d data : ",a)') trim(LE_Dat(idat)%name); call goErr
            TRACEBACK; status=1; return
        end select
        ! unit conversion:
        pat = pat * leo%unitconv(l)

        ! write data:
#ifdef with_netcdf
        call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l,n), -999, leo%itrec, &
                                        pat, status )
        IF_NOTOK_RETURN(status=1)                                          
#endif
      end do ! landuses
    end do ! data fields

    ! root?
    if ( goc%root ) then
#ifdef with_netcdf
      ! close
      status = NF90_Close( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#endif
      ! reset flag:
      leo%opened = .false.
      ! write GrADS ctl file if necessary:
      call Write_GrADS_Ctl( leo, status )
      IF_NOTOK_RETURN(status=1)

    end if  ! root

    ! ok
    status = 0

  end subroutine LE_Output_Budget_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use LE_Grid, only : glb_ugg

    ! --- in/out ---------------------------------

    type(T_LE_Output_Budget), intent(inout)    ::  leo
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Write_GrADS_Ctl'

    ! --- local ----------------------------------

    type(T_GrADS_Ctl)     ::  ctl
    integer               ::  l, idat, ilu, n
    character(len=512)    ::  line
    character(len=64)     ::  varname

    ! --- begin ----------------------------------

    ! write ctl file ?
    if ( leo%grads_ctl ) then

      ! ctl file name:
      write (leo%grads_ctl_file,'(a,"_",a,"_",a)') &
                trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
      write (leo%grads_ctl_file,'(a,".ctl")') trim(leo%grads_ctl_file)



          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! files with instant fields:


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
      call GrADS_Ctl_ZDef( ctl, (/1/), status )
      IF_NOTOK_RETURN(status=1)
      ! tdef:
      call GrADS_Ctl_TDef( ctl, leo%grads_ctl_nt, &
                               (/leo%grads_ctl_t1%year,leo%grads_ctl_t1%month,leo%grads_ctl_t1%day,leo%grads_ctl_t1%hour,leo%grads_ctl_t1%min/), &
                               (/                    0,                     0,leo%grads_ctl_dt%day,leo%grads_ctl_dt%hour,leo%grads_ctl_dt%min/), &
                               status )
      IF_NOTOK_RETURN(status=1)
      ! number of variables:
      call GrADS_Ctl_Vars( ctl, leo%ndep_total, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over data fields to be written:
      do l = 1, leo%ndep
        ! global tracer index
        idat = leo%idat(l)
        ! loop over landuses to be written:
        do n = 1, leo%nlu
          ! not landuse dependend ? then skip after first:
          if ( (.not. le_dat(idat)%landuse_dependent) .and. (n > 1) ) cycle
          ! landuse index:
          ilu = leo%ilu(n)
          varname = trim(leo%name_dat(l))
          ! add landuse name if necessary:
          if ( le_dat(idat)%landuse_dependent ) then
            varname = trim(varname)//'_'//trim(leo%name_lu(n))
          end if
          
          ! set variable lineiption:
          write (line,'(a," [",a,"]")') trim(varname), trim(leo%unit_dat(l))
          ! define variable:

              if ( LE_Dat(idat)%const ) then
                call GrADS_Ctl_Var( ctl, trim(varname), 1, 'y,x', trim(line), status )
                IF_NOTOK_RETURN(status=1)
              else
                call GrADS_Ctl_Var( ctl, trim(varname), 1, 't,y,x', trim(line), status )
                IF_NOTOK_RETURN(status=1)
              end if
        end do ! landuses      

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



end module LE_Output_Budget
