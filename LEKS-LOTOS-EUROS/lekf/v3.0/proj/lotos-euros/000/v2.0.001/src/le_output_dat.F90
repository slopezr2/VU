!###############################################################################
!
! NAME
!
!   LE_Output_dat  -  LOTOS-EUROS output of 3D fields
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

module LE_Output_Dat

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_dat

  public  ::  LE_Output_Dat_Init, LE_Output_Dat_Done
  public  ::  LE_Output_Dat_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_dat'


  ! maximum number of supported data sets:
  integer, parameter  ::  ndat = 200

  ! maximum number of levels in output:
  integer, parameter  ::  maxlev = 10


  ! --- types ------------------------------

  type T_LE_Dat
    character(len=32)      ::  name
    character(len=32)      ::  unit
    integer                ::  rank
    logical                ::  const
  end type T_LE_Dat

  type T_LE_Output_dat
    ! name for this file:
    character(len=16)           ::  typ
    character(len=16)           ::  name
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
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
    ! collect: daily, instant
    character(len=32)           ::  collect
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
    !integer                     ::  varid_time_day
    integer                     ::  varid_time_dtg
    ! database with supported variables:
    type(T_LE_Dat)              ::  LE_Dat(ndat)
    ! tracer variables:
    integer                     ::  ndat
    integer, pointer            ::  idat(:)
    character(len=32), pointer  ::  name_dat(:)
    character(len=32), pointer  ::  unit_dat(:)
    real, pointer               ::  unitconv(:)
    integer, pointer            ::  varid_dat(:)
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
  end type T_LE_Output_dat


contains


  ! ====================================================


  subroutine SetDat( d, name, unit, rank, const )

    ! --- in/out ----------------------------------

    type(T_LE_Dat), intent(out)     ::  d
    character(len=*), intent(in)    ::  name
    character(len=*), intent(in)    ::  unit
    integer, intent(in)             ::  rank
    logical, intent(in)             ::  const

    ! --- begin ----------------------------------

    d%name = name
    d%unit = unit
    d%rank = rank
    d%const = const

  end subroutine SetDat


  ! ====================================================


  subroutine LE_Output_Dat_Init( leo, rcfile, rckey, typ, name, status )

    use GO     , only : TrcFile, Init, Done, ReadRc
    use GO     , only : goMatchValues, goSplitString, goReadFromLine
    use GO     , only : AnyDate
    use Grid   , only : GetDistribution
    use Dims   , only : nx, ny, nz
    use Dims   , only : runF
    use LE_Grid, only : lli
    use LE_Output_Common, only : Init
    use LE_Data, only : LE_Data_Enable

    ! --- in/out --------------------------------

    type(T_LE_Output_dat), intent(out)    ::  leo
    character(len=*), intent(in)          ::  rcfile
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Dat_Init'

    ! rckey extensions:

    ! --- local ---------------------------------

    type(TrcFile)         ::  rcF
    character(len=64)     ::  basekey
    character(len=1024)   ::  field_names
    integer               ::  idat
    character(len=32)     ::  level_names
    integer               ::  ilev
    integer               ::  i, i1, i2

    character(len=512)    ::  key
    real                  ::  west, south, east, north
    real, pointer         ::  ff(:,:)

    integer               ::  imf

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

    ! collect daily or instant
    call ReadRc( rcF, trim(basekey)//'.collect', leo%collect, status )
    IF_NOTOK_RETURN(status=1)

    ! output time resolution:
    call ReadRc( rcF, trim(basekey)//'.dhour', leo%dhour, status )
    IF_NOTOK_RETURN(status=1)

    ! pack output?
    call ReadRC( rcF, trim(basekey)//'.pack_output', leo%pack_output, status )
    IF_NOTOK_RETURN(status=1)

    ! tracer names:
    call ReadRc( rcF, trim(basekey)//'.fields', field_names, status )
    IF_NOTOK_RETURN(status=1)

    ! level type::
    call ReadRc( rcF, trim(basekey)//'.levtype', leo%levtype, status )
    IF_NOTOK_RETURN(status=1)

    ! level descriptions:
    call ReadRc( rcF, trim(basekey)//'.levels', level_names, status )
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
          call goMatchValues( level_names, 0, nz+1, leo%nlev, leo%ilev, status )
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
        ! extract heights:
        call goSplitString( trim(level_names), leo%nlev, leo%heights, status )
        IF_NOTOK_RETURN(status=1)
        ! dummy ...
        allocate( leo%ilev(leo%nlev) )
        leo%ilev = -999
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! init empty:
    do i = 1, size(leo%LE_Dat)
      call SetDat( leo%LE_Dat(i), '', '', 0, .false. )
    end do

    ! define supported data:
    !   name   :  used in rcfile to identify this field;
    !   unit   :  in which the corresponding field in LE is expressed;
    !             the model units are converted to 'cf'-conventions
    !             (SI-units) on output;
    !   rank   :  2 for 2D field, and 3 for 3D
    !   const  :  logical, T for constant fields
    !                        name        unit            rank  const
    i = 1
    call SetDat( leo%LE_Dat(i), 'area'    , 'km2'           , 2, .true.  ); i=i+1
    call SetDat( leo%LE_Dat(i), 'oro'     , 'm'             , 2, .true.  ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lsm'     , '1'             , 2, .true.  ); i=i+1
    call SetDat( leo%LE_Dat(i), 'h'       , 'km'            , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'dh'      , 'm'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rhumid'  , '%'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'temper'  , 'K'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'tcc'     , '1'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'wsurf'   , 'm s-1'         , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'usurf'   , 'm s-1'         , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'vsurf'   , 'm s-1'         , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'tsurf'   , 'K'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'dsurf'   , 'K'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'sst'     , 'K'             , 2, .false. ); i=i+1    
    call SetDat( leo%LE_Dat(i), 'srh'     , '%'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'blh'     , 'km'            , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rad'     , 'W m-2'         , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rain'    , 'mm 3hr-1'      , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'snow'    , 'mm 3hr-1'      , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'sd'      , 'm'             , 2, .false. ); i=i+1
    !call SetDat( leo%LE_Dat(i), 'ctop'    , '-'             , 2, .false. ); i=i+1
    !call SetDat( leo%LE_Dat(i), 'cbot'    , '-'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'u'       , 'km min-1'      , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'v'       , 'km min-1'      , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'kz'      , 'm2 s-1'        , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'dens'    , 'kg m-3'        , 3, .false. ); i=i+1

    !call SetDat( leo%LE_Dat(i), 'vd_o3'   , 'm/s'           , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'z0m'     , 'm'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'monin'   , 'm'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'ustar'   , 'm s-1'         , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'zenith'  , '??'            , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'expcls'  , '  '            , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'snowcover', '1'            , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'swg'     , 'kg/kg'         , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'us_emdust_ara', 'm s-1'    ,2 , .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'us_emdust_des', 'm s-1'    ,2 , .false. ); i=i+1
    ! for deposition
    !call SetDat( leo%LE_Dat(i), 'ra '     , 's/m'           , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'vsrais'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'vsracs'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'vsrcos'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'vdrais'  , '-'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'vdracs'  , '-'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'vdrcoa'  , '-'             , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rdnucs'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rdaits'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rdaccs'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rdcoas'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rdaiti'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rdacci'  , '-'             , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'rdcoai'  , '-'             , 3, .false. ); i=i+1                                               

    ! for new additions
   ! pressure not explicitly defined in the model
    call SetDat( leo%LE_Dat(i), 'psurf'   , 'Pa'            , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'slpress'   , 'Pa'          , 2, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), 'pressure', 'Pa'            , 3, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), '', ''      , 3, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), '', ''      , 3, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), '', ''      , 3, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), '', ''      , 3, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), '', ''      , 2, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), '', ''      , 2, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), '', ''      , 2, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), 'vs_na_f', 'm/s'     , 3, .false. ); i=i+1
!    call SetDat( leo%LE_Dat(i), 'vs_na_c', 'm/s'     , 3, .false. ); i=i+1

    ! in-cloud chemistry
    call SetDat( leo%LE_Dat(i), 'pH'        , '-'           , 3,  .false.); i=i+1

    ! for co2 bio emis
    call SetDat( leo%LE_Dat(i), 'radd'      , 'W m-2'       , 2, .false. ); i=i+1

    ! country fraction for NL
    call SetDat( leo%LE_Dat(i), 'cnt_NLD'   , '1'         , 2, .true. ); i=i+1

    ! sibcasa gpp flux for co2
    call SetDat( leo%LE_Dat(i), 'gpp_co2'   , '(kg C)/m2/s' , 2, .false. ); i=i+1

    ! landuse fractions:
    call SetDat( leo%LE_Dat(i), 'lu_grs'   , '1'         , 2, .true. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lu_ara'   , '1'         , 2, .true. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lu_crp'   , '1'         , 2, .true. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lu_cnf'   , '1'         , 2, .true. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lu_dcf'   , '1'         , 2, .true. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lu_wat'   , '1'         , 2, .true. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lu_urb'   , '1'         , 2, .true. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lu_tun'   , '1'         , 2, .true. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'lu_des'   , '1'         , 2, .true. ); i=i+1 
    call SetDat( leo%LE_Dat(i), 'lu_ice'   , '1'         , 2, .true. ); i=i+1 
    call SetDat( leo%LE_Dat(i), 'lu_sav'   , '1'         , 2, .true. ); i=i+1 
    call SetDat( leo%LE_Dat(i), 'lu_trf'   , '1'         , 2, .true. ); i=i+1 
    call SetDat( leo%LE_Dat(i), 'lu_wai'   , '1'         , 2, .true. ); i=i+1 
    call SetDat( leo%LE_Dat(i), 'lu_med'   , '1'         , 2, .true. ); i=i+1 

    ! needed to create soil water average
    call SetDat( leo%LE_Dat(i), 'slt'       , '1'        , 2, .true.  ); i=i+1
    call SetDat( leo%LE_Dat(i), 'swvl1'     , 'm3/m3'    , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'swvl2'     , 'm3/m3'    , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'swvl3'     , 'm3/m3'    , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'swvl4'     , 'm3/m3'    , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'smi2'      , '1'        , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'smi3'      , '1'        , 2, .false. ); i=i+1

    ! pressure output
    call SetDat( leo%LE_Dat(i), 'mpressure' , 'Pa'       , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'hpressure' , 'Pa'       , 3, .false. ); i=i+1

#ifdef with_pollen
    ! accumulated heatsum
    call SetDat( leo%LE_Dat(i), 'accum_heatsum_polb', 'degree days', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'accum_heatsum_polo', 'degree days', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'rel_amt_polb_left' , '1', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'rel_amt_polo_left' , '1', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'rel_amt_polg_left' , '1', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'start_day_polg'    , '-', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'end_day_polg'    , '-', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'len_day_polg'    , '-', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'max_avail_polg'  , '-', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'rest_avail_polg'  , '-', 2, .false.); i=i+1
    call SetDat( leo%LE_Dat(i), 'gamma_field'  , '-', 2, .false.); i=i+1
#endif

    ! cloud coverage
    call SetDat( leo%LE_Dat(i), 'icc'   , '1', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'bcc'   , '1', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'clwc'  , '1', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'iclwc' , '1', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'praini', '1', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'raini' , 'mm/h', 3, .false. ); i=i+1
    
    ! volume(fluxes)
    call SetDat( leo%LE_Dat(i), 'ovol'  , 'km3'    , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'vol'   , 'km3'    , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'afluxx', 'km3/min', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'afluxy', 'km3/min', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'afluxz', 'km3/min', 3, .false. ); i=i+1
    
    ! cfl numbers
    call SetDat( leo%LE_Dat(i), 'cflx'  , '1'      , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'cfly'  , '1'      , 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'cflz'  , '1'      , 3, .false. ); i=i+1
    
#ifdef with_radiation
    ! AOD
    call SetDat( leo%LE_Dat(i), 'tau'        , '-', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'extinction' , '-', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'ssa'        , '-', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'asy'        , '-', 3, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'AOD2'       , '-', 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'angstrom'   , '-', 2, .false. ); i=i+1
#endif

    ! setup storage for tracer fields:
    allocate( leo%idat     (ndat) )
    allocate( leo%name_dat (ndat) )
    allocate( leo%unit_dat (ndat) )
    allocate( leo%unitconv (ndat) )
    allocate( leo%varid_dat(ndat) )

    ! match tracer names:
    call goMatchValues( field_names, leo%LE_Dat(:)%name, &
                          leo%ndat, leo%name_dat, leo%idat, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected fields for data output:")'); call goPr
    do i = 1, leo%ndat
      idat = leo%idat(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3, a10," ",a10,")")') &
                  i, leo%name_dat(i), &
                  idat, trim(leo%LE_Dat(idat)%name), trim(leo%LE_Dat(idat)%unit); call goPr
    end do

    ! * enable meteo
    
    ! needed for output at heights or altitudes:
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over all tracer to be written:
    do i = 1, leo%ndat
      ! global tracer index:
      idat = leo%idat(i)
      ! switch
      select case ( trim(leo%LE_Dat(idat)%name) )
        case ( 'oro', 'lsm', 'dh', 'tsurf', 'dsurf', 'srh', 'blh', 'dens', &
               'icc', 'bcc', 'clwc', 'iclwc', 'raini', 'praini', &
               'slt', 'swvl1', 'swvl2', 'swvl3', 'smi2', 'smi3' )
          call LE_Data_Enable( trim(leo%LE_Dat(idat)%name), status )
          IF_NOTOK_RETURN(status=1)
        case ( 'sst' )
          call LE_Data_Enable( 'sstk', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'usurf' )
          call LE_Data_Enable( 'u10', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'vsurf' )
          call LE_Data_Enable( 'v10', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'snow' )
          call LE_Data_Enable( 'sf', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'sd' )
          call LE_Data_Enable( 'sd', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'swg' )
          call LE_Data_Enable( 'swg', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'tcc' )
          call LE_Data_Enable( 'occ', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'rad' )
          call LE_Data_Enable( 'ssr', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'radd' )
          call LE_Data_Enable( 'ssrd', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'hpressure' )
          call LE_Data_Enable( 'hp', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'mpressure' )
          call LE_Data_Enable( 'hp', status )
          IF_NOTOK_RETURN(status=1)
          call LE_Data_Enable( 'p', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'temper' )
          call LE_Data_Enable( 't', status )
          IF_NOTOK_RETURN(status=1)
        case ( 'rhumid' )
          call LE_Data_Enable( 'rh', status )
          IF_NOTOK_RETURN(status=1)
      end select
    end do  ! output variables

    ! *

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
      leo%westb = runF%westb
      leo%southb = runF%southb
    else
      ! read domain from key
      read(key,*,iostat=status) west, east, south, north
      if(status/=0) then
        write (gol,'("could not read domain from key: ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      endif

      ! Check if bounding box is in run domain
      if ( west < runF%westb .or. east > runF%eastb .or. &
           south < runF%southb .or. north > runF%northb ) then
        write( gol, '("Bounding box domain is (partly) outside run domain")' ) ; call goErr
        write( gol, '("Run domain: ", 4f8.2)' ) runF%westb,runF%eastb,runF%southb,runF%northb ; call goErr
        write( gol, '("Bounding Box domain: ", 4f8.2)' ) west, east, south, north ; call goErr
        TRACEBACK;status=1;return
      endif

      ! for safety
      nullify(ff)
      ! get cell range covered by box
      call GetDistribution(lli,west,east,south,north,leo%i1,leo%i2,leo%j1,leo%j2,ff,status)
      IF_NOTOK_RETURN(status=1)
      !clear, fractions not used
      if ( associated(ff) ) deallocate(ff)
      ! set shape
      leo%ni = leo%i2-leo%i1+1
      leo%nj = leo%j2-leo%j1+1
      ! set west/south bounds
      leo%westb  = lli%blon_deg(leo%i1-1)
      leo%southb = lli%blat_deg(leo%j1-1)
    end if

    ! close:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Dat_Init


  ! ***


  subroutine LE_Output_Dat_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif
    use LE_Output_Common, only : Done

    ! --- in/out --------------------------------

    type(T_LE_Output_dat), intent(inout)   ::  leo
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Dat_Done'

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
        ! nothing to be done ..
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

  end subroutine LE_Output_Dat_Done


  ! ***


  subroutine LE_Output_Dat_PutOut( leo, t, status )

    use Binas  , only : xm_air
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>), rTotal, dTotal, iTotal
    use GO     , only : wrtgol, Precisely, MidNight
    use GO     , only : goMatchValue
    use Num    , only : Interp_Lin

#ifdef with_netcdf
    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT
#endif

    use LE_CF_Conventions   , only : LE_CF_names

    use Dims   , only : nx, ny, nz
    use Dims   , only : runF  
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Meteo_Data, only : ovolume, volume, afluxx, afluxy, afluxz
    use LE_Meteo_Data, only : cflx, cfly, cflz

    use Dims   , only : Kz
    !use Dims   , only : vd
    !use Dims   , only : ra, rb
    use Dims   , only : monin, ustar, zenith, expcls, coszen
#ifdef with_m7    
    use LE_M7  , only : rdrym7modes ,rwetm7modes !,waterm7modes, densm7modes
    !use LE_M7  , only : vsmass,vdmass
#endif

    !! testing co2 ...
    !use LE_DryDepos, only : gpp_co2

    !use deposition      , only: vd_lu, ftop_lu
    !use LE_Sedim, only : vs
    !use LE_Chem_Cloud   , only : pH


    use indices ,  only: i_no, i_no2,i_o3, i_ppm_f, i_ppm_c, i_so2, i_nh3
    !use indices ,  only: i_Na_f, i_Na_c

    use LE_Emis, only : emis_set, max_emis

    use LE_Country, only : country_map_code, country_map_frac

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

    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_Data, only : ilu_grass, ilu_arable, ilu_permanent_crops, &
                                ilu_coniferous_forest, ilu_deciduous_forest, &
                                ilu_water_sea, ilu_urban, ilu_other, ilu_desert, &
                                ilu_ice, ilu_savanna, ilu_tropical_forest, &
                                ilu_water_inland, ilu_mediterrean_scrub
    use LE_Landuse_Data, only : ustar_lu_dust_emis         
    use LE_Landuse_Data, only : z0m                        
#ifdef with_radiation
    use LE_Radiation,     only : tau, extinction, asy, ssa, angstrom                            
#endif

    ! --- in/out --------------------------------

    type(T_LE_Output_dat), intent(inout)  ::  leo
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Dat_PutOut'

    ! --- local ---------------------------------

    integer               ::  time6(6)
    integer               ::  time
#ifdef with_netcdf
    integer               ::  cmode
#endif
    integer               ::  idat, ilev, ilev1, ilevn
    integer               ::  varid
    type(TDate)           ::  t0
    real, allocatable     ::  lons(:), lats(:)
    integer               ::  i, j, k, l
    real                  ::  pat(nx,ny)
    real                  ::  field(nx,ny,0:nz+1)

    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units
    character(len=512)    ::  comment

    real                  ::  hh(nx,ny,0:nz+1)
    real                  ::  hsamp(nx,ny)
    integer               ::  ilast

    integer               ::  icountry
    integer               ::  iemis
    integer               ::  itr

    character(len=256) :: commandline

    integer               ::  imf
    logical               ::  match

    ! meteo data:
    real, pointer        ::   oro(:,:,:)      ! (lon,lat,1)
    real, pointer        ::   hpres(:,:,:)    ! (lon,lat,lev)
    real, pointer        ::   h_m(:,:,:)      ! (lon,lat,lev)
    real, pointer        ::   pdata(:,:,:)    ! (lon,lat,1)

    ! --- begin ---------------------------------

#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif

    ! for multiples of dhour only ...
    if ( .not. Precisely(t,leo%dhour,'hour')  ) then
      status=0; return
    end if

    ! point to meteo data:
    call LE_Data_GetPointer( 'oro', oro, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)

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
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! unknonw ...
        case default
          write (gol,'("unsupported collect key : ",a)') trim(leo%collect); call goErr
          TRACEBACK; status=1; return
      end select

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
      call PutOut_GlobalAttributes( leo%com, leo%ncid, status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_netcdf
      !
      ! define dimensions:
      !

      ! lon/lat dimensions
      call LE_Output_Define_Dims_Lon_Lat(leo%ncid, leo%dimid_lon, leo%dimid_lat, &
                                         leo%ni, leo%nj, trim(leo%com%CF_convention), status)
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

      ! lon/lat variables
      call LE_Output_Define_Vars_Lon_Lat(leo%ncid, leo%varid_lon, leo%varid_lat, &
                                         leo%dimid_lon, leo%dimid_lat, trim(leo%com%CF_convention), status)
      IF_NOTOK_RETURN(status=1)

      ! level variables
      call LE_Output_Define_Vars_Lev(leo%ncid, leo%varid_lev,leo%dimid_lev, &
                                       trim(leo%levtype), trim(leo%com%CF_convention), status)
      IF_NOTOK_RETURN(status=1)

      ! time since t0
      t0 = leo%com%t0
      ! time variables
      call LE_Output_Define_Vars_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                      leo%dimid_time, trim(leo%com%CF_convention), t0, status)
      IF_NOTOK_RETURN(status=1)
#endif

      ! loop over data fields to be written:
      do l = 1, leo%ndat

        ! global tracer index
        idat = leo%idat(l)

        ! CF standard name for concentration/mixing ratio/column:

        ! initial comment:
        select case ( trim(leo%LE_Dat(idat)%name) )
          case ( 'mpressure' )
            comment = 'pressure at mid of layer (full level)'
          case ( 'hpressure' )
            comment = 'pressure at top of layer (upper half level)'
          case default
            comment = ''
        end select

        ! get names following CF conventions;
        ! store conversion factor for later usage:
        call LE_CF_names( &
                     leo%LE_Dat(idat)%name, leo%LE_Dat(idat)%unit, &
                     cf_standard_name, cf_long_name, cf_units, &
                     leo%unitconv(l), comment, &
                     status )
        IF_NOTOK_RETURN(status=1)

        ! define variable:
        select case ( leo%LE_Dat(idat)%rank )
          case ( 2 )
            if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf
              status = NF90_Def_Var( leo%ncid, trim(leo%name_dat(l)), NF90_REAL, &
                                       (/leo%dimid_lon,leo%dimid_lat/), varid )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'latitude longitude')
              IF_NF90_NOTOK_RETURN(status=1)
#endif
            else
#ifdef with_netcdf
              status = NF90_Def_Var( leo%ncid, trim(leo%name_dat(l)), NF90_REAL, &
                                       (/leo%dimid_lon,leo%dimid_lat,leo%dimid_time/), varid )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time latitude longitude')
              IF_NF90_NOTOK_RETURN(status=1)
#endif
            end if
          case ( 3 )
            if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf
              status = NF90_Def_Var( leo%ncid, trim(leo%name_dat(l)), NF90_REAL, &
                                       (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev/), varid )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'level latitude longitude')
              IF_NF90_NOTOK_RETURN(status=1)
#endif
            else
#ifdef with_netcdf
              status = NF90_Def_Var( leo%ncid, trim(leo%name_dat(l)), NF90_REAL, &
                                       (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), varid )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time level latitude longitude')
              IF_NF90_NOTOK_RETURN(status=1)
#endif
            end if
          case default
            write (gol,'("unsupported data rank : ",i4)') leo%LE_Dat(idat)%rank; call goErr
            TRACEBACK; status=1; return
        end select

        ! write attributes:
#ifdef with_netcdf
        if ( len_trim(cf_standard_name) > 0 ) then
          status = nf90_put_att( leo%ncid, varid, 'standard_name', trim(cf_standard_name) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
        status = nf90_put_att( leo%ncid, varid, 'long_name', trim(cf_long_name) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! write units:
        status = nf90_put_att( leo%ncid, varid, 'units', trim(cf_units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! specials ...
        select case ( trim(leo%LE_Dat(idat)%name) )
          case ( 'dens' )
            ! moleweight of air:
            status = nf90_put_att( leo%ncid, varid, 'moleweight_air', xm_air )
            IF_NF90_NOTOK_RETURN(status=1)
            status = nf90_put_att( leo%ncid, varid, 'moleweight_air_unit', 'kg mole-1' )
            IF_NF90_NOTOK_RETURN(status=1)
        end select
        ! store units for later usage (GrADS ctl file):
        leo%unit_dat(l) = trim(cf_units)
        ! add comment:
        if ( len_trim(comment) > 0 ) then
          status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
#endif

        ! store variable id:
        leo%varid_dat(l) = varid

      end do  ! written tracers

      ! end defintion mode:

#ifdef with_netcdf
      status = NF90_EndDef( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#endif

      ! no records written yet:
      leo%itrec = 0

    end if

    ! next time record:
    leo%itrec = leo%itrec + 1

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

      ! write longitudes:
      allocate( lons(leo%i1:leo%i2) )
      do i = leo%i1, leo%i2
        lons(i) = runF%westb + (i-0.5)*runF%dlon
      end do

      allocate( lats(leo%j1:leo%j2) )
      do j = leo%j1, leo%j2
        lats(j) = runF%southb + (j-0.5)*runF%dlat
      end do
#ifdef with_netcdf

      call LE_Output_Put_Var_Lon_Lat( leo%ncid, leo%varid_lon, leo%varid_lat, &
                                      leo%ni, leo%nj, lons, lats, &
                                      trim(leo%com%CF_convention), status)
      IF_NOTOK_RETURN(status=1)

#endif
      deallocate( lons )
      deallocate( lats )

#ifdef with_netcdf
      ! write level indices
      call LE_Output_Put_Var_Lev(leo%ncid, leo%varid_lev, leo%nlev, &
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

    ! sample ?
    if ( (trim(leo%levtype) == 'heights'   ) .or. &
         (trim(leo%levtype) == 'elevations') ) then
      ! lowest is orography:
      hh(:,:,0) = oro(:,:,1)  ! m
      ! mid of layers:
      hh(:,:,1) = oro(:,:,1) + h_m(:,:,1) * 0.5  ! m
      do k = 2, nz
        hh(:,:,k) = oro(:,:,1) + ( h_m(:,:,k-1) + h_m(:,:,k) ) * 0.5  ! m
      end do
      ! dummy for aloft ...
      hh(:,:,nz+1) = h_m(:,:,nz) + 1000.0   ! 1000 m above top
    end if

    ! loop over all tracer to be written:
    do l = 1, leo%ndat

      ! global tracer index:
      idat = leo%idat(l)

      ! constant fields written only once:
      if ( leo%LE_Dat(idat)%const .and. (leo%itrec > 1) ) cycle

      ! 2d or 3d ?
      select case ( leo%LE_Dat(idat)%rank )

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 2 )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! init flag:
          match = .false.

          ! extract 2d field:
          select case ( trim(leo%LE_Dat(idat)%name) )
            case ( 'oro' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'oro', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            case ( 'lsm' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'lsm', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            case ( 'slt' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'slt', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            case ( 'area' )
              do j = 1, ny
                pat(:,j) = runF%dx(j) * runF%dy
              end do
            case ( 'psurf' )
              ! point to half-level pressures:
              call LE_Data_GetPointer( 'hp', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,0)  ! surface pressure
            case ( 'wsurf' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'wsurf', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'tsurf' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'tsurf', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'dsurf' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'dsurf', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'sst' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'sstk', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'srh' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'srh', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'blh' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'blh', pdata, status, check_units='m' )    
              IF_NOTOK_RETURN(status=1)
              ! select:
              if ( trim(leo%LE_Dat(idat)%unit) == 'km' ) then
                pat = pdata(1:nx,1:ny,1) * 1.0e-3  ! m to km
              else
                write (gol,'("unsupported units `",a,"` for blh")') trim(leo%LE_Dat(idat)%unit); call goErr
              end if
            !
            case ( 'usurf' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'u10', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'vsurf' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'v10', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'snow' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'sf', pdata, status, check_units='m/s' )
              IF_NOTOK_RETURN(status=1)
              ! select:
              if ( trim(leo%LE_Dat(idat)%unit) == 'mm 3hr-1' ) then
                pat = pdata(1:nx,1:ny,1) * 1.0e3 * 3*3600.0  ! m/s -> mm/3hr
              else
                write (gol,'("unsupported units `",a,"` for `raini`")') trim(leo%LE_Dat(idat)%unit); call goErr
              end if
            !
            case ( 'sd' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'sd', pdata, status, check_units='m' )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'swg' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'swg', pdata, status, check_units='kg/kg' )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'tcc' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'occ', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'rad' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'ssr', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'radd' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'ssrd', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'rain' )

              ! point to meteo data:
              call LE_Data_GetPointer( 'rain', pdata, status, check_units='m/s' )
              IF_NOTOK_RETURN(status=1)
              ! select:
              if ( trim(leo%LE_Dat(idat)%unit) == 'mm 3hr-1' ) then
                pat = pdata(1:nx,1:ny,1) * 1.0e3 * 3*3600.0  ! m/s -> mm/3hr
              else
                write (gol,'("unsupported units `",a,"` for `raini`")') trim(leo%LE_Dat(idat)%unit); call goErr
              end if
            !
            case ( 'snowcover' )
              ! point to meteo data:
              call LE_Data_GetPointer( 'snowcover', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'swvl1', 'swvl2', 'swvl3', 'swvl4' )
              ! point to meteo data:
              call LE_Data_GetPointer( trim(leo%LE_Dat(idat)%name), pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)
            !
            case ( 'smi2', 'smi3' )
              ! point to meteo data:
              call LE_Data_GetPointer( trim(leo%LE_Dat(idat)%name), pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )
              IF_NOTOK_RETURN(status=1)
              ! select:
              pat = pdata(1:nx,1:ny,1)

            !case ( 'ctop'  ) ;  pat = ctop
            !case ( 'cbot'  ) ;  pat = cbot
            case ( 'monin'  ) ;  pat = monin
            case ( 'ustar'  ) ;  pat = ustar
            case ( 'z0m'    ) ;  pat = z0m
            case ( 'us_emdust_ara' ); pat=ustar_lu_dust_emis(:,:,ilu_arable)
            case ( 'us_emdust_des' ); pat=ustar_lu_dust_emis(:,:,ilu_desert)
            case ( 'zenith'  ) ;  pat = zenith
            case ( 'expcls'  ) ;  pat = expcls
            case ( 'coszen'  ) ; pat = coszen
#ifdef with_radiation
            case ( 'AOD2'    )
               pat=0.
                do ilev = 1, nz
                    pat =pat+ tau(:,:,ilev,10)
                end do
            case ( 'angstrom' ) ; pat = angstrom
#endif
            !case ( 'ra'    ) ;  pat = ra

#ifdef with_pollen
            !~ accumulated heatsum for ripening of Birch pollen
            case ( 'accum_heatsum_polb')
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                ! check if emission set is pollen
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                    ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_b' ) then
                      ! extract heatsum
                      pat = emis_set(iemis)%pollen%emp(itr)%Birch_Pollen%heatsum
                      exit
                    end if
                  end do
                end if
              end do
            !~ accumulated heatsum for ripening of Birch pollen          
            case ( 'accum_heatsum_polo') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                ! check if emission set is pollen
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                    ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_o' ) then
                      ! extract heatsum 
                      pat = emis_set(iemis)%pollen%emp(itr)%Olive_Pollen%heatsum
                      exit
                    end if
                  end do
                end if
              end do
            !~ Relative amount of Birch pollen left in the catkins
            case ( 'rel_amt_polb_left')
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_b' ) then
                      where (emis_set(iemis)%pollen%emp(itr)%Birch_Pollen%max_avail_grains > 0.0 )
                        pat = emis_set(iemis)%pollen%emp(itr)%Birch_Pollen%rest_avail_grains / &
                               emis_set(iemis)%pollen%emp(itr)%Birch_Pollen%max_avail_grains
                      elsewhere
                        pat = 0.0
                      end where
                    end if
                  end do
                end if
              end do

            !~ Relative amount of Olive pollen left 
            case ( 'rel_amt_polo_left') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_o' ) then
                      where (emis_set(iemis)%pollen%emp(itr)%Olive_Pollen%max_avail_grains > 0.0 )
                        pat = emis_set(iemis)%pollen%emp(itr)%Olive_Pollen%rest_avail_grains / &
                               emis_set(iemis)%pollen%emp(itr)%Olive_Pollen%max_avail_grains
                      elsewhere
                        pat = 0.0
                      end where
                    end if
                  end do
                end if
              end do

            !~ Relative amount of Grass pollen left 
            case ( 'rel_amt_polg_left') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_g' ) then
                      where (emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%max_avail_grains > 0.0 )
                        pat = emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%rest_avail_grains / &
                               emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%max_avail_grains
                      elsewhere
                        pat = 0.0
                      end where
                    end if
                  end do
                end if
              end do
            !~ start day of grass pollen
            case ( 'start_day_polg') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_g' ) then
                      pat = emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%startday_em
                      exit
                    end if
                  end do
                end if
              end do
            !~ end day of grass pollen
            case ( 'end_day_polg') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_g' ) then
                      pat = emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%endday_em
                      exit
                    end if
                  end do
                end if
              end do
            !~ length of grass pollen season
            case ( 'len_day_polg') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_g' ) then
                      pat = emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%endday_em-&
                            emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%startday_em
                      exit
                    end if
                  end do
                end if
              end do
            !~ max avialable grass pollen
            case ( 'max_avail_polg') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_g' ) then
                      pat = emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%max_avail_grains
                      exit
                    end if
                  end do
                end if
              end do
            !~ rest avialable grass pollen
            case ( 'rest_avail_polg') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_g' ) then
                      pat = emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%rest_avail_grains
                      exit
                    end if
                  end do
                end if
              end do
            !~ rest avialable grass pollen
            case ( 'gamma_field') 
              ! loop over emission input sets and find pollen set
              do iemis = 1, max_emis
                if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
                  ! loop over different pollen tracer emissions
                  do itr = 1, emis_set(iemis)%pollen%ntr
                  ! check if emitted pollen tracer is Birch
                    if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer) == 'pol_g' ) then
                      pat = emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%gamma_field
                      exit
                    end if
                  end do
                end if
              end do
#endif
            case ( 'cnt_NLD' )
              ! get index of country, extract 3-letter code from variable:
              call goMatchValue( leo%LE_Dat(idat)%name(5:7), country_map_code, icountry, status )
              IF_NOTOK_RETURN(status=1)
              ! copy:
              pat = country_map_frac(:,:,icountry)

            ! landuse fractions:
            case ( 'lu_grs' ) ; pat = lu_fracs(:,:,ilu_grass            )
            case ( 'lu_ara' ) ; pat = lu_fracs(:,:,ilu_arable           )
            case ( 'lu_crp' ) ; pat = lu_fracs(:,:,ilu_permanent_crops  )
            case ( 'lu_cnf' ) ; pat = lu_fracs(:,:,ilu_coniferous_forest)
            case ( 'lu_dcf' ) ; pat = lu_fracs(:,:,ilu_deciduous_forest )
            case ( 'lu_wat' ) ; pat = lu_fracs(:,:,ilu_water_sea        )
            case ( 'lu_urb' ) ; pat = lu_fracs(:,:,ilu_urban            )
            case ( 'lu_tun' ) ; pat = lu_fracs(:,:,ilu_other            )
            case ( 'lu_des' ) ; pat = lu_fracs(:,:,ilu_desert           )  
            case ( 'lu_ice' ) ; pat = lu_fracs(:,:,ilu_ice              )  
            case ( 'lu_sav' ) ; pat = lu_fracs(:,:,ilu_savanna          )  
            case ( 'lu_trf' ) ; pat = lu_fracs(:,:,ilu_tropical_forest  )  
            case ( 'lu_wai' ) ; pat = lu_fracs(:,:,ilu_water_inland     )  
            case ( 'lu_med' ) ; pat = lu_fracs(:,:,ilu_mediterrean_scrub)  


            !case ( 'vdrais' ) ; pat = vdmass(:,:,2)
            !case ( 'vdracs' ) ; pat = vdmass(:,:,3)
            !case ( 'vdrcoa' ) ; pat = vdmass(:,:,6)

            case default
              write (gol,'("unsupported 2d data : ",a)') trim(leo%LE_Dat(idat)%name); call goErr
              TRACEBACK; status=1; return
          end select

          ! unit conversion:
          pat = pat * leo%unitconv(l)

          ! write data:
          if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf
            status = NF90_Put_Var( leo%ncid, leo%varid_dat(l), pat(leo%i1:leo%i2,leo%j1:leo%j2) )
            IF_NF90_NOTOK_RETURN(status=1)
#endif
          else
#ifdef with_netcdf
            status = NF90_Put_Var( leo%ncid, leo%varid_dat(l), pat, &
                           start=(/1,1,leo%itrec/), count=(/leo%ni,leo%nj,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
#endif
          end if

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 3 )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! loop over all model levels, including surface and aloft:
          do ilev = 0, nz+1
            ! corresponding model layer, extract surface (0) from layer 1:
            ilev1 = max(ilev,1)
            ilevn = min(ilev1,nz)
            ! set to zero by default, might be used for ilev==nz+1
            pat = 0.0
            ! extract 2d field:
            select case ( trim(leo%LE_Dat(idat)%name) )
              case ( 'h' )
                if ( ilev == 0 ) then
                  pat = 2.0e-3  ! surface layer at 2 m, convert to km
                else if ( ilev <= nz ) then
                  pat = h_m(:,:,ilev)*1e-3  ! km
                else
                  pat = 200.0  ! km
                end if
              ! upper half level pressure:
              case ( 'hpressure' )
                ! point to meteo data:
                call LE_Data_GetPointer( 'hp', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) then
                  pat = pdata(1:nx,1:ny,ilev)  ! Pa, valid for ilev=0,..,nz
                else
                  pat = 0.0  ! Pa
                end if
              ! mid level pressure:
              case ( 'mpressure' )
                ! point to meteo data:
                call LE_Data_GetPointer( 'hp', hpres, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! point to meteo data:
                call LE_Data_GetPointer( 'p', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev == 0 ) then
                  pat = hpres(:,:,0)  ! use surface pressure for this
                else if ( ilev <= nz ) then
                  pat = pdata(1:nx,1:ny,ilev)  ! Pa
                else
                  pat = 0.0  ! Pa
                end if
              ! layer thickness:
              case ( 'dh' )
                ! point to meteo data:
                call LE_Data_GetPointer( 'dh', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev == 0 ) then
                  pat = 0.0  ! m
                else if ( ilev <= nz ) then
                  pat = pdata(1:nx,1:ny,ilev)  ! m
                else
                  pat = 2000.0  ! m
                end if
              ! temperature:
              case ( 'temper' )
                ! point to meteo data:
                call LE_Data_GetPointer( 't', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) then
                  pat = pdata(1:nx,1:ny,ilev1)
                else
                  pat = 293.0
                end if
              case ( 'dens' )
                ! point to meteo data:
                call LE_Data_GetPointer( 'dens', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) then
                  pat = pdata(1:nx,1:ny,ilev1)
                else
                  pat = pdata(1:nx,1:ny,nz)     ! <-- density at top layer
                end if
              case ( 'ovol' )
                !! point to meteo data:
                !call LE_Data_GetPointer( 'ovol', pdata, status, check_units='m3' )    
                !IF_NOTOK_RETURN(status=1)
                !! select:
                !pat = pdata(1:nx,1:ny,ilevn) * 1e-9  ! km3
                pat = ovolume(1:nx,1:ny,ilevn)
              case ( 'vol' )
                !! point to meteo data:
                !call LE_Data_GetPointer( 'vol', pdata, status, check_units='m3' )    
                !IF_NOTOK_RETURN(status=1)
                !! select:
                !pat = pdata(1:nx,1:ny,ilevn) * 1e-9  ! km3
                pat = volume(1:nx,1:ny,ilevn)
              case ( 'afluxx' )
                !! point to meteo data:
                !call LE_Data_GetPointer( 'uflux', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                !IF_NOTOK_RETURN(status=1)
                !! average over u-edges:
                !pat = 0.5 * ( pdata(0:nx-1,:,ilevn) + pdata(1:nx,:,ilevn) )
                ! average over x-edges:
                pat = 0.5 * ( afluxx(0:nx-1,1:ny,ilevn) + afluxx(1:nx,1:ny,ilevn) )
              case ( 'afluxy' )
                !! point to meteo data:
                !call LE_Data_GetPointer( 'vflux', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                !IF_NOTOK_RETURN(status=1)
                !! average over v-edges:
                !pat = 0.5 * ( pdata(:,0:ny-1,ilevn) + pdata(:,1:ny,ilevn) )
                ! average over y-edges:
                pat = 0.5 * ( afluxy(1:nx,0:ny-1,ilevn) + afluxy(1:nx,1:ny,ilevn) )
              case ( 'afluxz' )
                !! point to meteo data:
                !call LE_Data_GetPointer( 'wflux', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                !IF_NOTOK_RETURN(status=1)
                !! average over half levels:
                !pat = 0.5 * ( pdata(:,:,ilevn-1) + pdata(:,:,ilevn) )
                ! average over z-edges:
                pat = 0.5 * ( afluxz(1:nx,1:ny,ilevn-1) + afluxz(1:nx,1:ny,ilevn) )
              case ( 'cflx' )
                pat = cflx(1:nx,1:ny,ilevn)
              case ( 'cfly' )
                pat = cfly(1:nx,1:ny,ilevn)
              case ( 'cflz' )
                pat = cflz(1:nx,1:ny,ilevn)
              case ( 'rhumid' )
                ! point to meteo data:
                call LE_Data_GetPointer( 'rh', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) pat = pdata(1:nx,1:ny,ilev1)
              case ( 'kz' )
                ! flux     = - Kz dc/dz
                ! mol/m2/s = m2/s (mol/m3)/m
                if ( ilev <= nz ) pat = Kz(:,:,ilev1)  ! m2/s ?

#ifdef with_m7
              ! M7 varriables
              !case ( 'vsrais' )
              !  if ( ilev <= nz ) pat = vsmass(:,:,ilev1,2)  !ppb ?
              !case ( 'vsracs' )
              !  if ( ilev <= nz ) pat = vsmass(:,:,ilev1,3)  !ppb ?
              !case ( 'vsrcos' )
              !  if ( ilev <= nz ) pat = vsmass(:,:,ilev1,4)  !ppb ?
              case ( 'rdnucs' )
                if ( ilev <= nz ) pat = rdrym7modes(:,:,ilev1,1) !ppb ?
              case ( 'rdaits' )
                if ( ilev <= nz ) pat = rdrym7modes(:,:,ilev1,2) !ppb ?
              case ( 'rdaccs' )
                if ( ilev <= nz ) pat = rdrym7modes(:,:,ilev1,3) !ppb ?
              case ( 'rdcoas' )
                if ( ilev <= nz ) pat = rdrym7modes(:,:,ilev1,4) !ppb ?
              case ( 'rdaiti' )
                if ( ilev <= nz ) pat = rwetm7modes(:,:,ilev1,5)  !ppb ?
                case ( 'rdacci' )
                if ( ilev <= nz ) pat = rwetm7modes(:,:,ilev1,6)  !ppb ?
              case ( 'rdcoai' )
                if ( ilev <= nz ) pat = rwetm7modes(:,:,ilev1,7)  !ppb ?
#endif                
              case ( 'icc'  )
                ! point to meteo data:
                call LE_Data_GetPointer( 'icc', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) pat = pdata(1:nx,1:ny,ilev1)
              case ( 'bcc'  )
                ! point to meteo data:
                call LE_Data_GetPointer( 'bcc', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) pat = pdata(1:nx,1:ny,ilev1)
              case ( 'clwc'   )
                ! point to meteo data:
                call LE_Data_GetPointer( 'clwc', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) pat = pdata(1:nx,1:ny,ilev1)
              case ( 'iclwc'  )
                ! point to meteo data:
                call LE_Data_GetPointer( 'iclwc', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) pat = pdata(1:nx,1:ny,ilev1)
              case ( 'praini'  )
                ! point to meteo data:
                call LE_Data_GetPointer( 'praini', pdata, status, check_units=trim(leo%LE_Dat(idat)%unit) )    
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( ilev <= nz ) pat = pdata(1:nx,1:ny,ilev1)
              case ( 'raini'  )
                ! point to meteo data:
                call LE_Data_GetPointer( 'raini', pdata, status, check_units='m/s' )
                IF_NOTOK_RETURN(status=1)
                ! select:
                if ( trim(leo%LE_Dat(idat)%unit) == 'mm/h' ) then
                  if ( ilev <= nz ) pat = pdata(1:nx,1:ny,ilev1) * 1.0e3 * 3600.0
                else
                  write (gol,'("unsupported units `",a,"` for `raini`")') trim(leo%LE_Dat(idat)%unit); call goErr
                end if
#ifdef with_radiation
              case ( 'tau' )
                if ( ilev <= nz ) then
                  pat = tau(:,:,ilev1,10) ! ~~visible band, 550 nm
                else
                  pat = tau(:,:,nz, 10)     ! <-- aod top layer
                end if
				                
              case ( 'extinction' )
                if ( ilev <= nz ) then
                  pat = extinction(:,:,ilev1,10)
                else
                  pat = extinction(:,:,nz,10)     ! <-- extinction top layer
                end if
              case ( 'ssa' )
                if ( ilev <= nz ) then
                  pat = ssa(:,:,ilev1,10)
                else
                  pat = ssa(:,:,nz,10)     ! <-- ssa top layer
                end if
              case ( 'asy' )
                !print *,'write asy'
                if ( ilev <= nz ) then
                  pat = asy(:,:,ilev1,10)
                else
                  pat = asy(:,:,nz,10)     ! <-- asy top layer
                end if
#endif

              case default
                write (gol,'("unsupported 3d data : ",a)') trim(leo%LE_Dat(idat)%name); call goErr
                TRACEBACK; status=1; return
            end select
            ! unit conversion:
            pat = pat * leo%unitconv(l)
            ! store:
            field(:,:,ilev) = pat

          end do  ! model layers

          ! * vertical selection or interpolation:

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
                ! write data:
                if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf
                  status = NF90_Put_Var( leo%ncid, leo%varid_dat(l), pat(leo%i1:leo%i2,leo%j1:leo%j2), &
                                 start=(/1,1,k/), count=(/nx,ny,1/) )
                  IF_NF90_NOTOK_RETURN(status=1)
#endif
                else
#ifdef with_netcdf
                  status = NF90_Put_Var( leo%ncid, leo%varid_dat(l), pat(leo%i1:leo%i2,leo%j1:leo%j2), &
                                 start=(/1,1,k,leo%itrec/), count=(/leo%ni,leo%nj,1,1/) )
                 IF_NF90_NOTOK_RETURN(status=1)
#endif
                end if
              end do   ! selected model layers

            !~ height levels:
            case ( 'heights', 'elevations' )

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
                ! loop over horizontal cells:
                do j = 1, ny
                  do i = 1, nx
                    ! vertical interpolation:
                    call Interp_Lin( hh(i,j,:), field(i,j,:), hsamp(i,j), pat(i,j), ilast, status )
                    IF_NOTOK_RETURN(status=1)
                  end do  ! i
                end do  ! j
                ! write data:
                if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf
                  status = NF90_Put_Var( leo%ncid, leo%varid_dat(l), pat(leo%i1:leo%i2,leo%j1:leo%j2), &
                                 start=(/1,1,k/), count=(/leo%ni,leo%nj,1/) )
                  IF_NF90_NOTOK_RETURN(status=1)
#endif
                else
#ifdef with_netcdf
                  status = NF90_Put_Var( leo%ncid, leo%varid_dat(l), pat(leo%i1:leo%i2,leo%j1:leo%j2), &
                                 start=(/1,1,k,leo%itrec/), count=(/leo%ni,leo%nj,1,1/) )
                 IF_NF90_NOTOK_RETURN(status=1)
#endif
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

          write (gol,'("unsupported data rank : ",i4)') leo%LE_Dat(idat)%rank; call goErr
          TRACEBACK; status=1; return

      end select

    end do   ! variables
    ! next time exceeds interval ?
    if ( t+IncrDate(hour=int(leo%dhour)) > leo%tr(2) ) then
      ! close:
#ifdef with_netcdf
      status = NF90_Close( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)

      ! pack output?
      if ( leo%pack_output ) then
        write( commandline, '(a," " a)' ) trim(leo%com%packscript), trim(leo%fname)
        call SYSTEM( trim(commandline) )
      end if
#endif
      ! reset flag:
      leo%opened = .false.
      ! write GrADS ctl file if necessary:
      call Write_GrADS_Ctl( leo, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine LE_Output_Dat_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : runF
    use dims, only : nx, ny
    use GrADS_Ctl

    ! --- in/out ---------------------------------

    type(T_LE_Output_dat), intent(inout)    ::  leo
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Write_GrADS_Ctl'

    ! --- local ----------------------------------

    type(T_GrADS_Ctl)     ::  ctl
    integer               ::  l, idat
    character(len=512)    ::  line

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
      ! xdef:
      call GrADS_Ctl_XYDef( ctl, 'XDEF', nx, runF%westb +0.5*runF%dlon, runF%dlon, status )
      IF_NOTOK_RETURN(status=1)
      ! ydef:
      call GrADS_Ctl_XYDef( ctl, 'YDEF', ny, runF%southb+0.5*runF%dlat, runF%dlat, status )
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
      call GrADS_Ctl_Vars( ctl, leo%ndat, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over data fields to be written:
      do l = 1, leo%ndat
        ! global tracer index
        idat = leo%idat(l)
        ! set variable lineiption:
        write (line,'(a," [",a,"]")') trim(leo%name_dat(l)), trim(leo%unit_dat(l))
        ! define variable:
        select case ( leo%LE_Dat(idat)%rank )
          case ( 2 )
            if ( leo%LE_Dat(idat)%const ) then
              call GrADS_Ctl_Var( ctl, trim(leo%name_dat(l)), 1, 'y,x', trim(line), status )
              IF_NOTOK_RETURN(status=1)
            else
              call GrADS_Ctl_Var( ctl, trim(leo%name_dat(l)), 1, 't,y,x', trim(line), status )
              IF_NOTOK_RETURN(status=1)
            end if
          case ( 3 )
            if ( leo%LE_Dat(idat)%const ) then
              call GrADS_Ctl_Var( ctl, trim(leo%name_dat(l)), leo%nlev, 'z,y,x', trim(line), status )
              IF_NOTOK_RETURN(status=1)
            else
              call GrADS_Ctl_Var( ctl, trim(leo%name_dat(l)), leo%nlev, 't,z,y,x', trim(line), status )
              IF_NOTOK_RETURN(status=1)
            end if
          case default
            write (gol,'("unsupported data rank : ",i4)') leo%LE_Dat(idat)%rank; call goErr
            TRACEBACK; status=1; return
        end select
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



end module LE_Output_dat
