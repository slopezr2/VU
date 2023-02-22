!###############################################################################
!
! NAME
!
!   LE_Output_vd  -  LOTOS-EUROS output of deposition velocities, diagnostic calculations on different heights
!
!   Relevant variables:
!
!     Ra(dh,ilu)         : Atmoshperic resistance over interval [z0,z0+dh], 
!                          where dh represents an observation height
!                          and htop the height where the observation profile starts
!                          Depends on landuse since z0 depends on it too.
!
!     Rb(ispec,ilu)      : Laminar boundary layer resistance.
!                          Always landuse dependend 
!                          For O3 fluxes following (McNaughton and Van der Hurk, 1995)
!
!     Rc_eff(ispec,ilu)  : Effective total surface resistance, depends on vegetation (landuse) and tracer.
!
! HISTORY
!
!   2012 feb, Arjo Segers, TNO
!   2013 Aug, Roy Wichink Kruit, TNO
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

module LE_Output_vd_diag

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf  
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_vd_diag

  public  ::  LE_Output_vd_diag_Init, LE_Output_vd_diag_Done
  public  ::  LE_Output_vd_diag_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_vd_diag'


  ! maximum number of supported data sets:
  integer, parameter  ::  ndat = 25
  
  ! maximum number of levels in output:
  integer, parameter  ::  maxlev = 10

  ! --- types ------------------------------

  type T_LE_Output_vd_diag
    ! name for this file:
    character(len=32)           ::  typ
    character(len=32)           ::  name
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
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
    integer                     ::  dimid_time
    ! dimension variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_lev
    integer                     ::  varid_time
    !integer                     ::  varid_time_day
    integer                     ::  varid_time_dtg
    ! data variables:
    integer                     ::  ndat
    integer, pointer            ::  idat(:)
    character(len=32), pointer  ::  name_dat(:)
    character(len=32), pointer  ::  unit_dat(:)
    real, pointer               ::  unitconv(:)
    integer, pointer            ::  varid_dat(:,:,:)  ! (ndat,ntr,nlu)
    integer                     ::  ndat_total  ! total data variables (needed for ctl file)
    ! tracer variables:
    integer                     ::  ntr
    integer, pointer            ::  itr(:)
    character(len=32), pointer  ::  name_tr(:)
    ! landuse variables:
    integer                     ::  nlu
    integer, pointer            ::  ilu(:)
    character(len=32), pointer  ::  name_lu(:)
    ! level selection:
    character(len=16)           ::  levtype
    integer                     ::  nlev
    integer, pointer            ::  ilev(:)
    real                        ::  heights(maxlev)
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
  end type T_LE_Output_vd_diag

  type T_LE_Dat
    character(len=32)      ::  name
    character(len=32)      ::  unit
    integer                ::  rank
    logical                ::  const
    logical                ::  tracer_dependend
    logical                ::  landuse_dependend
  end type T_LE_Dat


  ! --- var --------------------------------------

  type(T_LE_Dat)     ::  le_dat(ndat)


contains


  ! ====================================================


  subroutine SetDat( d, name, unit, rank, const, tracer_dependend, landuse_dependend )

    ! --- in/out ----------------------------------

    type(T_LE_Dat), intent(out)     ::  d
    character(len=*), intent(in)    ::  name
    character(len=*), intent(in)    ::  unit
    integer, intent(in)             ::  rank
    logical, intent(in)             ::  const
    logical, intent(in)             ::  tracer_dependend
    logical, intent(in)             ::  landuse_dependend

    ! --- begin ----------------------------------

    d%name = name
    d%unit = unit
    d%rank = rank
    d%const = const
    d%tracer_dependend  = tracer_dependend
    d%landuse_dependend = landuse_dependend

  end subroutine SetDat


  ! ====================================================


  subroutine LE_Output_vd_diag_Init( leo, rcF, rckey, typ, name, state, status )

    use GO              , only : TrcFile
    use GO              , only : goMatchValues, goSplitString
    use GO              , only : AnyDate
    use Dims            , only : nx, ny, nz
    use Indices         , only : nspec, specname
    use LE_Landuse_Data , only : nlu, lu_name, lu_name_abbr
    use LE_Landuse_Data , only : ludep_output_whole_grid
    use LE_Grid         , only : ugg
    use LE_Output_Common, only : Init
    use LE_Data, only : LE_Data_Enable

    ! --- in/out --------------------------------

    type(T_LE_Output_vd_diag), intent(out)::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  state
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_vd_diag_Init'

    ! rckey extensions:

    ! --- local ---------------------------------

    character(len=64)     ::  basekey
    character(len=1024)   ::  selected_names
    integer               ::  idat, itr, ilu
    character(len=32)     ::  level_names
    integer               ::  ilev
    integer               ::  i

    character(len=512)    ::  key
    real                  ::  west, south, east, north
    real, pointer         ::  ff(:,:)

    ! --- begin ---------------------------------

    ! store name:
    leo%typ  = typ
    leo%name = name

    ! state name:
    leo%state = trim(state)

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
    
    ! collect daily or instant
    call rcF%Get( trim(basekey)//'.collect', leo%collect, status )
    IF_NOTOK_RETURN(status=1)

    ! output time resolution:
    call rcF%Get( trim(basekey)//'.dhour', leo%dhour, status )
    IF_NOTOK_RETURN(status=1)

    ! define supported data:
    !   name   :  used in rcfile to identify this field;
    !   unit   :  in which the corresponding field in LE is expressed;
    !             the model units are converted to 'cf'-conventions
    !             (SI-units) on output;
    !   rank   :  2 for 2D field, and 3 for 3D
    !   const  :  logical, T for constant fields in time
    !   tracer :  logical, T for multiple components
    !   landuse:  logical, T for landuse dependent
    !                        name              unit                 rank  const , tracer, landuse
    call SetDat( LE_Dat( 1), 'Ra'            , 's/m'                , 3, .false., .false., .true. )
    call SetDat( LE_Dat( 2), 'vd'            , 'm/s'                , 3, .false., .true. , .true. )
    call SetDat( LE_Dat( 3), 'conc_sfc_lu'   , 'ug m-3'             , 2, .false., .true. , .true. )

    ! variable names:
    call rcF%Get( trim(basekey)//'.fields', selected_names, status )
    IF_NOTOK_RETURN(status=1)    
    ! setup storage for tracer fields:
    allocate( leo%idat     (ndat) )
    allocate( leo%name_dat (ndat) )
    allocate( leo%unit_dat (ndat) )
    allocate( leo%unitconv (ndat) )
    ! match variable names:
    call goMatchValues( selected_names, LE_Dat(:)%name, &
                          leo%ndat, leo%name_dat, leo%idat, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected fields for vd output:")'); call goPr
    do i = 1, leo%ndat
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
    allocate( leo%name_tr (nspec) )
    ! match tracer names:
    call goMatchValues( selected_names, specname, &
                          leo%ntr, leo%name_tr, leo%itr, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected tracers for vd output:")'); call goPr
    do i = 1, leo%ntr
      itr = leo%itr(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3, a10,")")') &
                  i, leo%name_tr(i), itr, trim(specname(itr)); call goPr
    end do

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
    write (gol,'("selected landues for vd output:")'); call goPr
    do i = 1, leo%nlu
      ilu = leo%ilu(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3," ",a,")")') &
                  i, leo%name_lu(i), ilu, trim(lu_name(ilu)); call goPr
      ludep_output_whole_grid(ilu) = .true.                  
    end do

    ! storage for written variables:
    allocate( leo%varid_dat(leo%ndat,leo%ntr,leo%nlu) )
    
    ! level type::
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
        allocate( leo%ilev(1+nz+1) )
        ! set selected level indices:
        call goMatchValues( level_names, 0, nz+1, leo%nlev, leo%ilev, status )
        IF_NOTOK_RETURN(status=1)
        write (gol,'("selected levels for conc output:")'); call goPr
        do i = 1, leo%nlev
          ilev = leo%ilev(i)
          write (gol,'("  ",i3,"  ",i3)') i, ilev; call goPr
        end do
      !~ height levels:
      case ( 'heights' )
        ! extract heights:
        call goSplitString( trim(level_names), leo%nlev, leo%heights, status )
        IF_NOTOK_RETURN(status=1)
        ! dummy indices:
        allocate( leo%ilev(leo%nlev) )
        leo%ilev = -999
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select

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

    ! enable data:
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_vd_diag_Init


  ! ***


  subroutine LE_Output_vd_diag_Done( leo, status )

#ifdef with_netcdf  
    use NetCDF          , only : NF90_Close
#endif
    use LE_Output_Common, only : Done

    ! --- in/out --------------------------------

    type(T_LE_Output_vd_diag), intent(inout)   ::  leo
    integer, intent(out)                       ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_vd_diag_Done'

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
      IF_NF90_NOTOK_RETURN(status=1)
    end if

    ! clear storage for tracer fields:
    deallocate( leo%idat      )
    deallocate( leo%name_dat  )
    deallocate( leo%unit_dat  )
    deallocate( leo%unitconv  )
    deallocate( leo%varid_dat )

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
    call Done( leo%com, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_vd_diag_Done


  ! ***


  subroutine LE_Output_vd_diag_PutOut( leo, t, c, bud, status )

    use Binas  , only : xm_air, xm_o3
    use Binas  , only : kappa_stab
    use GO     , only : goc
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>), rTotal, iTotal
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

    use Dims            , only : nx, ny, nz, nspec
    use LE_DryDepos     , only : Rb
    use LE_DryDepos     , only : Rc_eff
    use LE_DryDepos     , only : Rs
    use LE_DryDepos     , only : vs
    use LE_DryDepos     , only : mix2ground_lu
    use JAQL_Stability  , only : Atmospheric_Resistance, f_h_stability 
    use LE_DryDepos     , only : spectype_depos
    use LE_DryDepos     , only : SPECTYPE_DEPOS_AERO_ZHANG, SPECTYPE_DEPOS_GAS_DEPAC
    use LE_Landuse_Data , only : z0m_depac, z0h_depac
    use LE_Landuse_Data , only : ustar_lu
    use Indices         , only : specname, specunit, specmolm
    use LE_Landuse_Data , only : lu_name
    use LE_Landuse_Data , only : nlu
    use LE_Budget_DryDepos, only : T_DryDepos_Budget, ex_drydepo

    use LE_Output_Common, only : PutOut_GlobalAttributes

    use LE_Output_Tools , only : LE_Output_Define_Dims_Lev
    use LE_Output_Tools , only : LE_Output_Define_Dims_Time

    use LE_Output_Tools , only : LE_Output_Define_Vars_Lev
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time

    use LE_Output_Tools , only : LE_Output_Put_Var_Lev
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains

    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out --------------------------------

    type(T_LE_Output_vd_diag), intent(inout)  ::  leo
    type(TDate), intent(in)                   ::  t
    real, intent(in)                          ::  c(nx,ny,nz,nspec)
    type(T_DryDepos_Budget), intent(in)       ::  bud   
    integer, intent(out)                      ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_vd_diag_PutOut'

    ! --- local ---------------------------------

    integer               ::  time6(6)
    integer               ::  time
#ifdef with_netcdf  
    integer               ::  cmode
#endif
    integer               ::  idat, itr, ilu, ilev, ilev1
    integer               ::  varid
    type(T_Grid_NcDef)    ::  gncd
    type(TDate)           ::  t0
    integer               ::  i, j, k, l, m, n, iz
    real                  ::  pat(nx,ny)
    real                  ::  Ra(nx,ny)
    real                  ::  fh(nx,ny)
    real, allocatable     ::  cground(:,:,:,:)

    character(len=64)     ::  varname
    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units
    character(len=512)    ::  comment

    character(len=256) :: commandline

    ! meteo data:
    real, pointer         ::  dens(:,:,:)   ! (lon,lat,alt)    
    real, pointer         ::  monin_inv(:,:,:)   ! (lon,lat,alt)        

    ! --- begin ----------------------------
    
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'monin_inv', monin_inv, status, check_units ='1/m' )
    IF_NOTOK_RETURN(status=1)

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
        ! dummmy ...
        leo%ncid = -1
#endif

        ! reset flag:
        leo%opened = .true.

        ! write global attributes:
        call PutOut_GlobalAttributes( leo%com, leo%ncid, status )
        IF_NOTOK_RETURN(status=1)

#ifdef with_netcdf
        !
        ! define dimensions:
        !

        ! grid dimensions/variables
        call glb_ugg%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                    dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )
                                    !subset=(/leo%i1,leo%i2,leo%j1,leo%j2/) )
        IF_NOTOK_RETURN(status=1)

        ! level dimension
        call LE_Output_Define_Dims_Lev(leo%ncid, leo%dimid_lev, leo%nlev, trim(leo%com%CF_convention), status)
        IF_NOTOK_RETURN(status=1)

        ! time dimensions
        call LE_Output_Define_Dims_Time(leo%ncid, leo%dimid_time, status)
        IF_NOTOK_RETURN(status=1)

        !
        ! define variables:
        !

        ! level variables
        call LE_Output_Define_Vars_Lev( leo%ncid, leo%varid_lev, leo%dimid_lev, &
                                          trim(leo%levtype), trim(leo%com%CF_convention), status)
        IF_NOTOK_RETURN(status=1)

        ! time since t0
        t0 = leo%com%t0      
        ! time variables
        call LE_Output_Define_Vars_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                        leo%dimid_time, trim(leo%com%CF_convention), t0, status)
        IF_NOTOK_RETURN(status=1)
#endif

      end if  ! root
      
      ! total number of data variables
      leo%ndat_total = 0
      ! loop over data fields to be written:
      do l = 1, leo%ndat
        ! variable index:
        idat = leo%idat(l)
        
        ! loop over tracers to be written:
        do m = 1, leo%ntr
          ! not tracer dependend ? then skip after first:
          if ( (.not. le_dat(idat)%tracer_dependend) .and. (m > 1) ) cycle
          ! tracer index:
          itr = leo%itr(m)          
          
          ! loop over landuses to be written:
          do n = 1, leo%nlu
            ! not landuse dependend ? then skip after first:
            if ( (.not. le_dat(idat)%landuse_dependend) .and. (n > 1) ) cycle
            ! landuse index:
            ilu = leo%ilu(n)
            
            ! variable name:
            varname = trim(leo%name_dat(l))
            ! add tracer name if necessary:
            if ( le_dat(idat)%tracer_dependend ) then
              varname = trim(varname)//'_'//trim(leo%name_tr(m))
            end if
            ! add landuse name if necessary:
            if ( le_dat(idat)%landuse_dependend ) then
              varname = trim(varname)//'_'//trim(leo%name_lu(n))
            end if

            ! initial comment:
            select case ( trim(LE_Dat(idat)%name) )
              case default
                comment = ''
            end select
            
            ! long name:
            select case ( trim(leo%name_dat(l)) )
              case ( 'area' )
                cf_long_name = 'grid cell area'
              case ( 'Ra' )
                write (cf_long_name,'("atmospheric resistance [z0,z0+h] above ",a)') trim(lu_name(ilu))
              case ( 'vd' )
                write (cf_long_name,'("deposition velocity above ",a," for ",a)') trim(lu_name(ilu)), trim(specname(itr))
              case ( 'conc_sfc_lu' )
                write (cf_long_name,'("Surface concentration of ",a," over ",a)') trim(specname(itr)), trim(lu_name(ilu))
              case default
                write (gol,'("could not set longname for variable `",a,"`")') trim(leo%name_dat(l)); call goErr
                TRACEBACK; status=1; return
            end select
            
            ! same units as in model:
            cf_units = trim(LE_Dat(idat)%unit)
            leo%unitconv(l) = -999.9

            ! store units for later usage (GrADS ctl file):
            leo%unit_dat(l) = trim(cf_units)

            ! root only:
            if ( goc%root ) then

#ifdef with_netcdf
              ! define variable:
              select case ( LE_Dat(idat)%rank )
                case ( 2 )
                  if ( LE_Dat(idat)%const ) then
                    status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                             (/leo%dimid_lon,leo%dimid_lat/), varid )
                    IF_NF90_NOTOK_RETURN(status=1)
                    status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'latitude longitude')
                    IF_NF90_NOTOK_RETURN(status=1)
                  else
                    status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                             (/leo%dimid_lon,leo%dimid_lat,leo%dimid_time/), varid )
                    IF_NF90_NOTOK_RETURN(status=1)
                    status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time latitude longitude')
                    IF_NF90_NOTOK_RETURN(status=1)
                  end if
                case ( 3 )
                  if ( LE_Dat(idat)%const ) then
                    status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                             (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev/), varid )
                    IF_NF90_NOTOK_RETURN(status=1)
                    status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'level latitude longitude')
                    IF_NF90_NOTOK_RETURN(status=1)
                  else
                    status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                             (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), varid )
                    IF_NF90_NOTOK_RETURN(status=1)
                    status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time level latitude longitude')
                    IF_NF90_NOTOK_RETURN(status=1)
                  end if
                case default
                  write (gol,'("unsupported data rank : ",i4)') LE_Dat(idat)%rank; call goErr
                  TRACEBACK; status=1; return
              end select

              ! write attributes:
              !status = nf90_put_att( leo%ncid, varid, 'standard_name', trim(cf_standard_name) )
              !IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( leo%ncid, varid, 'long_name', trim(cf_long_name) )
              IF_NF90_NOTOK_RETURN(status=1)
              ! write units:
              status = nf90_put_att( leo%ncid, varid, 'units', trim(cf_units) )
              IF_NF90_NOTOK_RETURN(status=1)
              call ugg%DefCoor_NetCDF( gncd, varid, status )
              IF_NF90_NOTOK_RETURN(status=1)
              ! add comment:
              if ( len_trim(comment) > 0 ) then
                status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
                IF_NF90_NOTOK_RETURN(status=1)
              end if
#endif

              ! store variable id:
              leo%varid_dat(l,m,n) = varid

              ! total number of data variables
              leo%ndat_total = leo%ndat_total + 1
              
            end if  ! root
            
          end do   ! landuses

        end do   ! tracers

      end do  ! written variables

      ! root only:
      if ( goc%root ) then

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
        call LE_Output_Put_Var_Lev( leo%ncid, leo%varid_lev, leo%nlev, &
                                   trim(leo%levtype), leo%ilev, leo%heights, status)
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

!    ! sample ?    
!    if ( (trim(leo%levtype) == 'heights'   ) .or. &
!         (trim(leo%levtype) == 'elevations') ) then
!      ! lowest is orography:
!      hh(:,:,0) = oro(:,:)  ! m
!      ! mid of layers:
!      hh(:,:,1) = oro(:,:) + h(:,:,1)*1e3 * 0.5
!      do k = 2, nz
!        hh(:,:,k) = oro(:,:) + ( h(:,:,k-1) + h(:,:,k) )*1e3 * 0.5
!      end do
!      ! dummy for aloft ...
!      hh(:,:,nz+1) = h(:,:,nz)*1e3 + 1000.0   ! 1000 m above top
!    end if

    ! loop over all tracer to be written:
    do l = 1, leo%ndat
      ! global tracer index:
      idat = leo%idat(l)
  
      if ( trim(LE_Dat(idat)%name) == 'conc_sfc_lu' ) then
        !
        allocate ( cground( nx,ny,nspec,nlu) )
        call mix2ground_lu(c, cground, status )
        IF_NOTOK_RETURN(status=1)
      end if

        
      ! constant fields written only once:
      if ( LE_Dat(idat)%const .and. (leo%itrec > 1) ) cycle

      ! loop over tracers to be written:
      do m = 1, leo%ntr
        ! not tracer dependend ? then skip after first:
        if ( (.not. le_dat(idat)%tracer_dependend) .and. (m > 1) ) cycle
        ! tracer index:
        itr = leo%itr(m)
        
        ! loop over landuses to be written:
        do n = 1, leo%nlu
          ! not landuse dependend ? then skip after first:
          if ( (.not. le_dat(idat)%landuse_dependend) .and. (n > 1) ) cycle
          ! landuse index:
          ilu = leo%ilu(n)

          ! 2d or 3d ?
          select case ( LE_Dat(idat)%rank )

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case ( 2 )
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              pat = 0.0                    
              ! extract 2d field:
              select case ( trim(LE_Dat(idat)%name) )
                
                case ( 'conc_sfc_lu' )
                  
                  ! pick concentration at ground level for each lu-class                  
                  pat = cground(:,:,itr,ilu)
                  
                  ! unit conversion; 'cdry' is 'concentation * height', 
                  ! thus 'ppb m' or 'ug/m3 m'
                  select case ( trim(specunit(itr))//' -> '//trim(LE_Dat(idat)%unit) )
                  
                    case ( 'ug/m3 -> ug m-3' )                                     
                      ! do nothing
                      
                    case ( 'ppb -> ug m-3' )
                      ! surface layer
                      iz = 1
                      ! ug/m3     =  ppb    * ((mole tr/mole air)/ppb) * (kg tr/mole tr) *   (kg air/m3 air)  / (kg air/mole air) * ug/kg
                      pat         =  pat    *       1e-9               *  specmolm(itr)  * dens(1:nx,1:ny,iz) /     xm_air        * 1e9
                      
                    case default
                      write (gol,'("do not know how to convert `",a,"` from `",a," m` to `",a,"`")') &
                              trim(specname(itr)), trim(specunit(itr)), trim(LE_Dat(idat)%unit); call goErr
                      TRACEBACK; status=1; return
                  end select
                      
                case default
                  write (gol,'("unsupported 2d data : ",a)') trim(LE_Dat(idat)%name); call goErr
                  TRACEBACK; status=1; return

              end select
              ! write data:
              if ( LE_Dat(idat)%const ) then
                ! write 2D field, no level, no time index:
                call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l,m,n), -1, -1, &
                                                  pat, status )
                IF_NOTOK_RETURN(status=1)
              else
                ! write 2D field, no level, wit time index:
                call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l,m,n), -1, leo%itrec, &
                                                  pat, status )
                IF_NOTOK_RETURN(status=1)
              end if

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case ( 3 )
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
              ! which output levels ?
              select case ( trim(leo%levtype) )
    
                !~ height levels:
                case ( 'heights' )
    
                  ! loop over target levels:
                  do k = 1, leo%nlev
                    ! extract 2D field:
                    select case ( trim(LE_Dat(idat)%name) )
          
                      ! atmospheric resistance (for gaseous tracers)
                      case ( 'Ra' )
                        ! first stability function for heat
                        fh = f_h_stability( z0h_depac(ilu)+leo%heights(k), z0h_depac(ilu), monin_inv(:,:,1) )
                        ! calc atmospheric resistance
                        Ra =  Atmospheric_Resistance( kappa_stab, ustar_lu(:,:,ilu), fh )                                                        
                        ! set in output variable
                        pat = Ra
                      ! deposition velocity (for gaseous tracers):
                      case ( 'vd' )
                        ! first stability function for heat
                        fh = f_h_stability( z0h_depac(ilu)+leo%heights(k), z0h_depac(ilu), monin_inv(:,:,1) )
                        ! calc atmospheric resistance
                        Ra =  Atmospheric_Resistance( kappa_stab, ustar_lu(:,:,ilu), fh )                                                                                                       

                        ! combine:
                        ! deposition type of tracer?
                        select case ( spectype_depos (itr) )
                          
                          case ( SPECTYPE_DEPOS_GAS_DEPAC )
                           ! gaseous tracer (vd = 1/(Ra+Rb+Rc) )
                            where( Rc_eff (:,:,itr,ilu) < 0.0 )
                              ! negative effective resistance 
                              ! (surface concentration larger than atmospheric concentration --> netto emission --> negative vd )
                              pat = 1.0 / (  Rc_eff(:,:,itr,ilu) - Rb(:,:,itr,ilu) - Ra )
                            elsewhere( Rc_eff( :,:,itr,ilu) >= 0.0 )
                              pat = 1.0 / (  Rc_eff(:,:,itr,ilu) + Rb(:,:,itr,ilu) + Ra )
                            endwhere
                            
                          case ( SPECTYPE_DEPOS_AERO_ZHANG )
                            ! aerosol tracer (vd = 1/(Ra + Rs) + vs
                            pat = 1.0 / ( Ra + Rs(:,:,itr,ilu) ) + vs(:,:,itr)
                          
                          ! unknown ..
                          case default
                            write (gol,'("unsupported deposition type : ",i4)') spectype_depos (itr); call goErr
                            TRACEBACK; status=1; return
                          end select
                          
                      ! unknown ..
                      case default
                        write (gol,'("unsupported 3d data : ",a)') trim(LE_Dat(idat)%name); call goErr
                        TRACEBACK; status=1; return
                    end select
                    ! write data:
                    if ( LE_Dat(idat)%const ) then
                      ! write 2D field for current level, no time index:
                      call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l,m,n), k, -1, &
                                                        pat, status )
                      IF_NOTOK_RETURN(status=1)
                    else
                      ! write 2D field for current level, with time index:
                      call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l,m,n), k, leo%itrec, &
                                                        pat, status )
                      IF_NOTOK_RETURN(status=1)
                    end if
                  end do   ! heights
    
                !~ unknown ...
                case default
                  write (gol,'("unsupported level type : ",a)') trim(leo%levtype); call goErr
                  TRACEBACK; status=1; return
    
              end select

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            case default
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              write (gol,'("unsupported data rank : ",i4)') LE_Dat(idat)%rank; call goErr
              TRACEBACK; status=1; return

          end select

        end do   ! landuses

      end do   ! tracers
      
      if ( trim(LE_Dat(idat)%name) == 'conc_sfc_lu' ) then
        ! done
        deallocate ( cground )
      end if

    end do   ! variables

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

  end subroutine LE_Output_vd_diag_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use LE_Grid, only : glb_ugg

    ! --- in/out ---------------------------------

    type(T_LE_Output_vd_diag), intent(inout)::  leo
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Write_GrADS_Ctl'

    ! --- local ----------------------------------

    type(T_GrADS_Ctl)     ::  ctl
    integer               ::  l, idat
    character(len=512)    ::  line
    character(len=64)     ::  varname
    integer               ::  itr, ilu, m, n

    ! --- begin ----------------------------------

    ! write ctl file ?
    if ( leo%grads_ctl ) then

      ! ctl file name:
      write (leo%grads_ctl_file,'(a,"_",a,"_",a)') &
                trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
      write (leo%grads_ctl_file,'(a,".ctl")') trim(leo%grads_ctl_file)

      ! daily or less ?
      select case ( trim(leo%collect) )
        ! collect daily for [00,24)
        case ( 'daily', 'daily24' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! files with instant fields:
        case ( 'instant' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2_%h2%n2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
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
      call GrADS_Ctl_Vars( ctl, leo%ndat_total, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over data fields to be written:
      do l = 1, leo%ndat
        ! global tracer index
        idat = leo%idat(l)
        
        ! loop over tracers to be written:
        do m = 1, leo%ntr
          ! not tracer dependend ? then skip after first:
          if ( (.not. le_dat(idat)%tracer_dependend) .and. (m > 1) ) cycle
          ! tracer index:
          itr = leo%itr(m)          
          
          ! loop over landuses to be written:
          do n = 1, leo%nlu
            ! not landuse dependend ? then skip after first:
            if ( (.not. le_dat(idat)%landuse_dependend) .and. (n > 1) ) cycle
            ! landuse index:
            ilu = leo%ilu(n)
            
            ! variable name:
            varname = trim(leo%name_dat(l))
            ! add tracer name if necessary:
            if ( le_dat(idat)%tracer_dependend ) then
              varname = trim(varname)//'_'//trim(leo%name_tr(m))
            end if
            ! add landuse name if necessary:
            if ( le_dat(idat)%landuse_dependend ) then
              varname = trim(varname)//'_'//trim(leo%name_lu(n))
            end if

            ! set variable lineiption:
            write (line,'(a," [",a,"]")') trim(varname), trim(leo%unit_dat(l))
            ! define variable:
            select case ( LE_Dat(idat)%rank )
              case ( 2 )
                if ( LE_Dat(idat)%const ) then
                  call GrADS_Ctl_Var( ctl, trim(varname), 1, 'y,x', trim(line), status )
                  IF_NOTOK_RETURN(status=1)
                else
                  call GrADS_Ctl_Var( ctl, trim(varname), 1, 't,y,x', trim(line), status )
                  IF_NOTOK_RETURN(status=1)
                end if
              case ( 3 )
                if ( LE_Dat(idat)%const ) then
                  call GrADS_Ctl_Var( ctl, trim(varname), leo%nlev, 'z,y,x', trim(line), status )
                  IF_NOTOK_RETURN(status=1)
                else
                  call GrADS_Ctl_Var( ctl, trim(varname), leo%nlev, 't,z,y,x', trim(line), status )
                  IF_NOTOK_RETURN(status=1)
                end if
              case default
                write (gol,'("unsupported data rank : ",i4)') LE_Dat(idat)%rank; call goErr
                TRACEBACK; status=1; return
            end select
          end do
        end do
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



end module LE_Output_vd_diag
