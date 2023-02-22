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

module LE_Output_Budget_Label

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_Budget_Label

  public  ::  LE_Output_Budget_Label_Init, LE_Output_Budget_Label_Done
  public  ::  LE_Output_Budget_Label_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Budget_Label'


  ! maximum number of supported data sets:
  integer, parameter  ::  ndep = 40


  ! --- types ------------------------------

  type T_LE_Output_Budget_Label
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
    integer                     ::  dimid_label
    integer                     ::  dimid_labelname_len
    ! dimension variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_lev
    integer                     ::  varid_time
    integer                     ::  varid_time_dtg
    integer                     ::  varid_time_day
    integer                     ::  varid_labelname
    ! tracer variables:
    integer                     ::  ndep
    integer                     ::  ndep_total
    integer, pointer            ::  idat(:)
    integer                     ::  ntr
    integer, pointer            ::  itr(:)
    integer, pointer            ::  itr_glob(:)
    character(len=10), pointer  ::  name_dat(:)
    character(len=32), pointer  ::  unit_dat(:)
    character(len=32), pointer  ::  name_tr(:)
    real, pointer               ::  unitconv(:)
    integer, pointer            ::  varid_dat(:,:,:)
    ! landuse variables:
    integer                     ::  nlu
    integer, pointer            ::  ilu(:)
    integer, pointer            ::  ilu_glob(:)
    character(len=32), pointer  ::  name_lu(:)
    character(len=32), pointer  ::  name_lu_glob(:)

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

  end type T_LE_Output_Budget_Label

  type T_LE_Dat
    character(len=32)      ::  name
    character(len=32)      ::  unit
    integer                ::  rank
    logical                ::  const
    logical                ::  landuse_dependent
    logical                ::  tracer_dependent
  end type T_LE_Dat


  ! --- var --------------------------------------

  type(T_LE_Dat)     ::  le_dat(ndep)


contains


  ! ====================================================


  subroutine SetDat( d, name, unit, rank, const, landuse_dependent, tracer_dependent )

    ! --- in/out ----------------------------------

    type(T_LE_Dat), intent(out)     ::  d
    character(len=*), intent(in)    ::  name
    character(len=*), intent(in)    ::  unit
    integer, intent(in)             ::  rank
    logical, intent(in)             ::  const
    logical, intent(in)             ::  landuse_dependent
    logical, intent(in)             ::  tracer_dependent

    ! --- begin ----------------------------------

    d%name = name
    d%unit = unit
    d%rank = rank
    d%const = const
    d%landuse_dependent = landuse_dependent
    d%tracer_dependent  = tracer_dependent

  end subroutine SetDat


  ! ====================================================


  subroutine LE_Output_Budget_Label_Init( leo, rcF, rckey, typ, name, status )

    use GO     , only : TrcFile
    use GO     , only : goMatchValues
    use GO     , only : AnyDate
    use Dims   , only : nx, ny
    
    use Indices, only : nspec, specname
    use LE_Grid, only : ugg
    use LE_Landuse_Data, only : nlu, lu_name, lu_name_abbr
    use LE_Landuse_Data, only : ludep_output_whole_grid
    use LE_Data, only : LE_Data_Enable
      
    use SA_Labeling, only : labelled_specs_budget_output_names, landuse_budget_output_names
    use SA_Labeling, only : nbudget_label_output_tracer, nbudget_label_output_lu
        
    ! --- in/out --------------------------------

    type(T_LE_Output_Budget_Label), intent(out) ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Budget_Label_Init'

    ! rckey extensions:

    ! --- local ---------------------------------

    character(len=64)     ::  basekey
    character(len=256)    ::  selected_names
    character(len=256)    ::  field_names
    integer               ::  idat,ilu
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
    !                        name        unit       rank  const  , landuse , tracer
    call SetDat( LE_Dat(1), 'drydep'   , 'ug m-2',   2, .false. , .false., .true. )
    call SetDat( LE_Dat(2), 'drydep_lu', 'ug m-2',   2, .false. , .true. , .true. )
    call SetDat( LE_Dat(3), 'wetdep'   , 'ug m-2',   2, .false. , .false., .true. )
    
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

    ! tracer names:
    call rcF%Get( trim(basekey)//'.tracers', selected_names, status )
    IF_NOTOK_RETURN(status=1)
    ! setup storage for tracer fields:
    allocate( leo%itr     (nspec) )
    allocate( leo%itr_glob(nspec) )
    allocate( leo%name_tr (nspec) )
    
    ! match tracer names:
    call goMatchValues( selected_names, labelled_specs_budget_output_names(1:nbudget_label_output_tracer), &
                          leo%ntr, leo%name_tr, leo%itr, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! match tracer names:
    call goMatchValues( selected_names, specname, &
                          leo%ntr, leo%name_tr, leo%itr_glob, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! landuse names:
    call rcF%Get( trim(basekey)//'.landuses', selected_names, status )
    IF_NOTOK_RETURN(status=1)
    ! setup storage for tracer fields:
    allocate( leo%ilu     (nlu) )
    allocate( leo%ilu_glob(nlu) )
    allocate( leo%name_lu (nlu) )
    allocate( leo%name_lu_glob(nlu) )
    ! match tracer names:
    call goMatchValues( selected_names, landuse_budget_output_names(1:nbudget_label_output_lu), &
                          leo%nlu, leo%name_lu, leo%ilu, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! match tracer names:
    call goMatchValues( selected_names, lu_name_abbr, &
                          leo%nlu, leo%name_lu_glob, leo%ilu_glob, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected landuses for budget output:")'); call goPr
    do i = 1, leo%nlu
      ilu = leo%ilu_glob(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3," ",a,")")') &
                  i, leo%name_lu(i), ilu, trim(lu_name(ilu)); call goPr
      ludep_output_whole_grid(ilu) = .true.                  
    end do

    ! storage for written variables:
    allocate( leo%varid_dat(leo%ndep,leo%ntr,leo%nlu) )

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

  end subroutine LE_Output_Budget_Label_Init


  ! ***


  subroutine LE_Output_Budget_Label_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif

    ! --- in/out --------------------------------

    type(T_LE_Output_Budget_Label), intent(inout)   ::  leo
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Budget_Label_Done'

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
    call leo%com%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Budget_Label_Done


  ! ***


  subroutine LE_Output_Budget_Label_PutOut( leo, t, status )

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
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR
#endif

    use LE_CF_Conventions   , only : LE_CF_names

    use Dims   , only : nx, ny, nz
    
    use indices, only : specname, specunit, specmolm
    !use LE_GrOutput, only : bud_depos, cdry, o3max, bud_stomata, bud_stomata_exd !,cwet
    use LE_LandUse_Data, only : ilu_grass, ilu_arable, ilu_permanent_crops, &
                                ilu_coniferous_forest, ilu_deciduous_forest, &
                                ilu_water_sea, ilu_water_inland, ilu_urban, ilu_other, ilu_desert, &
                                ilu_ice, ilu_wheat, ilu_beech, ilu_spruce, &
                                ilu_semi_natural_veg

    use LE_Output_Tools , only : LE_Output_Define_Dims_Time
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_Time    
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains
    
    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer
    
    use SA_Labeling     , only : SA_nlabel, SA_Label_names
    use SA_Labeling     , only : labelname_len

    use SA_Labeling     , only : SA_budget_drydepos_day
    use SA_Labeling     , only : SA_budget_wetdepos_day
    use SA_Labeling     , only : SA_budget_drydepos_lu_day
    
    ! --- in/out --------------------------------

    type(T_LE_Output_Budget_Label), intent(inout)   ::  leo
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Budget_Label_PutOut'

    ! --- local ---------------------------------

    integer               ::  time6(6)
    integer               ::  time
#ifdef with_netcdf
    integer               ::  cmode
#endif
    integer               ::  idat, ilu, itr, itr_glob
    type(T_Grid_NcDef)    ::  gncd
    integer               ::  varid
    type(TDate)           ::  t0, tdum
    integer               ::  i, j, l, m, n
    integer               ::  labellen, SA_ilabel, iz
    integer               ::  ibud_spec
    real                  ::  pat(nx,ny,SA_nlabel)

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
      call leo%com%PutOut_GlobalAttributes( leo%ncid, status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_netcdf
      
      ! grid dimensions/variables
      call glb_ugg%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                  dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )
      IF_NOTOK_RETURN(status=1)
      
      ! time dimensions
      call LE_Output_Define_Dims_Time(leo%ncid, leo%dimid_time, status)
      IF_NOTOK_RETURN(status=1)

      ! label dimension
      status = NF90_Def_Dim( leo%ncid, 'label', SA_nlabel, leo%dimid_label )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! len_labelname dimension
      status = NF90_Def_Dim( leo%ncid, 'labelname_len', labelname_len, leo%dimid_labelname_len )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! time since t0
      t0 = leo%com%t0      
      ! time variables
      call LE_Output_Define_Vars_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                      leo%dimid_time, trim(leo%com%CF_convention), t0, status)
      IF_NOTOK_RETURN(status=1)

      ! labelname variable
      status = NF90_Def_Var( leo%ncid, 'labelnames', NF90_CHAR, (/leo%dimid_labelname_len, leo%dimid_label/), varid )
      IF_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( leo%ncid, varid, 'standard_name', 'Names of the label defnitions')
      IF_NOTOK_RETURN(status=1)
      leo%varid_labelname = varid      
#endif
      leo%ndep_total = 0
    end if  ! root

    ! loop over data fields to be written:
    do l = 1, leo%ndep

      ! global tracer index
      idat = leo%idat(l)
        
      ! loop over tracers to be written
      do m = 1, leo%ntr
        ! not tracer dependent? then skip after first:
        if ( (.not. le_dat(idat)%tracer_dependent) .and. (m > 1 ) ) cycle
        ! tracer index
        itr = leo%itr(m)
        ! global tracer index
        itr_glob = leo%itr_glob(m)
          
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
          ! add tracer name if necessary:
          if ( le_dat(idat)%tracer_dependent ) then
            varname = trim(varname)//'_'//trim(leo%name_tr(m))
          end if
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
                                   (/leo%dimid_lon,leo%dimid_lat,leo%dimid_label,leo%dimid_time/), varid, &
                                   deflate_level=leo%com%deflate_level )
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
            leo%varid_dat(l,m,n) = varid
          
            ! total number of data variables
            leo%ndep_total = leo%ndep_total + 1
          end if  ! root
          
        end do ! landuses
      end do  ! written tracers
    end do  ! depositions
#ifdef with_netcdf
    if ( goc%root ) then
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

        ! write labelnames
        do SA_ilabel = 1, SA_nlabel
          labellen = len(trim( SA_Label_Names(SA_ilabel)) )
          status = NF90_Put_Var(leo%ncid, leo%varid_labelname, trim(SA_Label_Names(SA_ilabel)), start=(/1,SA_ilabel/), count=(/labellen,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end do

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

    ! loop over all deposition types to be written:
    do l = 1, leo%ndep
      ! global tracer index:
      idat = leo%idat(l)
      
      ! constant fields, written only once:
      if ( LE_Dat(idat)%const .and. (l > 1) ) cycle
      
      ! loop over tracers to be written
      do m = 1, leo%ntr
        ! not tracer dependent? then skip after first
        if ( (.not. le_dat(idat)%tracer_dependent) .and. (m > 1) ) cycle
        ! tracer index
        itr = leo%itr(m)
        ! labelled tracer indec
        itr_glob = leo%itr_glob(m)
      
        ! constant fields written only once:
        ! loop over landuses to be written:
        do n = 1, leo%nlu
          ! not landuse dependend ? then skip after first:
          if ( (.not. le_dat(idat)%landuse_dependent) .and. (n > 1) ) cycle
          ! landuse index:
          ilu = leo%ilu(n)

          ! extract 2d field:
          select case ( trim(LE_Dat(idat)%name) )

            case('drydep_lu') ; pat = SA_budget_drydepos_lu_day(:,:,itr,ilu,:)
            case('drydep') ; pat = SA_budget_drydepos_day(:,:,itr,:)
            case('wetdep') ; pat = SA_budget_wetdepos_day(:,:,itr,:)
                
            case default
              write( gol, '("Unknown data type, a" )' ) trim(LE_Dat(idat)%name) ; call goErr
              TRACEBACK;status=1;return
              
          end select
              
#ifdef with_netcdf
          do SA_ilabel = 1, SA_nlabel
            call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l,m,n), -999, leo%itrec, &
                                            pat(:,:,SA_ilabel), status, ilab=SA_ilabel )
            IF_NOTOK_RETURN(status=1)                                          
          end do
#endif
        end do ! landuses
      end do ! tracers
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

  end subroutine LE_Output_Budget_Label_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use SA_Labeling, only : SA_nlabel
    use SA_Labeling, only : SA_Label_Names, SA_Short_Label_Names
    use LE_Grid, only : glb_ugg

    ! --- in/out ---------------------------------

    type(T_LE_Output_Budget_Label), intent(inout)    ::  leo
    integer, intent(out)                             ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Write_GrADS_Ctl'

    ! --- local ----------------------------------

    type(T_GrADS_Ctl)     ::  ctl
    integer               ::  l, idat, itr, ilu, n, m
    integer               ::  SA_ilabel, Nvars
    character(len=512)    ::  dim_descr
    character(len=512)    ::  line
    character(len=64)     ::  varname
    character(len=128)    ::  comment_line
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
      call GrADS_Ctl_Comment( ctl, 'Label definitions', status )
      call GrADS_Ctl_Comment( ctl, '', status )
      do SA_ilabel = 1, SA_nlabel
        write( comment_line , '(a," : ", a)' ) trim(SA_Short_Label_names(SA_ilabel)), trim(SA_Label_Names(SA_ilabel))
        call GrADS_Ctl_Comment( ctl, trim(comment_line) , status )
      end do
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
      Nvars = leo%ndep_total * SA_nlabel      
      call GrADS_Ctl_Vars( ctl, Nvars, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over data fields to be written:
      do l = 1, leo%ndep
        ! global tracer index
        idat = leo%idat(l)
        
        ! loop over tracers to be written:
        do m = 1, leo%ntr
          ! not tracer dependend ? then skip after first:
          if ( (.not. le_dat(idat)%tracer_dependent) .and. (m > 1) ) cycle
          ! tracer index:
          itr = leo%itr(m)          
                   
          ! loop over landuses to be written:
          do n = 1, leo%nlu
            ! not landuse dependend ? then skip after first:
            if ( (.not. le_dat(idat)%landuse_dependent) .and. (n > 1) ) cycle
            ! landuse index:
            ilu = leo%ilu(n)
            varname = trim(leo%name_dat(l))
            ! add tracer name if necessary:
            if ( le_dat(idat)%tracer_dependent ) then
              varname = trim(varname)//'_'//trim(leo%name_tr(m))
            end if
            ! add landuse name if necessary:
            if ( le_dat(idat)%landuse_dependent ) then
              varname = trim(varname)//'_'//trim(leo%name_lu(n))
            end if
            
            ! loop over labels:
            do SA_ilabel = 1, SA_nlabel
              
              ! set variable lineiption:
              write (line,'(a," [",a,"]")') trim(varname), trim(leo%unit_dat(l))
              
              if ( LE_Dat(idat)%const ) then
                ! set dimension description labelnr, y,x
                write (dim_descr, '(i2,",y,x")') SA_ilabel-1  
              else
                ! set dimension description t, labelnr, y,x
                write (dim_descr, '("t, ",i2,",y,x")') SA_ilabel-1  
              end if
              
              ! add variable description:              
              call GrADS_Ctl_Var( ctl, trim(varname), 1, trim(dim_descr), trim(line), status, label_name=SA_Short_Label_Names(SA_ilabel) )
              IF_NOTOK_RETURN(status=1)
            end do ! labels
            
          end do ! landuses      
        end do ! tracers
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



end module LE_Output_Budget_Label
