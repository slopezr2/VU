!###############################################################################
!
! NAME
!   LE_Emis_TNO_Base - access to TNO base emission file
!
! FILE CONTENT
!
!   A base emisison file contains the following data:
!    - country codes and names
!    - emisison categories
!    - emisison types (areao or point)
!    - emisison grid defintion
!    - base emissions; each is assigned:
!        type (point or area)
!        category
!        country
!        grid cell
!   A cocktail of tracers is emitted: NOx, SOx, NH3, CO, CH4, VOC, PM fine, PM coarse (or total)
!
!   Format: NetCDF
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_TNO_Base

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -------------------------------
  
  private
  
!  public  ::  N_SOURCETYPE, SOURCETYPE_LABEL, SOURCETYPE_AREA, SOURCETYPE_POINT

  public  ::  T_Emis_TNO_Base, T_Country_Info, T_SourceCat_Info, T_SourceType_Info, T_Emis_Info, T_Emis_Gathered
  
  public  ::  LE_Emis_TNO_Base_Init, LE_Emis_TNO_Base_Done
  public  ::  LE_Emis_TNO_Base_Summary
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_TNO_Base'
  
  ! number of category codes:
  integer, parameter  ::  N_CATCODE = 4
  integer, parameter  ::  I_SECTOR_CODE    = 1
  integer, parameter  ::  I_SUBSECTOR_CODE = 2
  integer, parameter  ::  I_FUEL_CODE      = 3
  integer, parameter  ::  I_SUBFUEL_CODE   = 4

  ! maximum lengths:
  integer, parameter  ::  LEN_NAME  = 32
  integer, parameter  ::  LEN_UNITS = 32
  
  
  ! --- local --------------------------------
  
  ! * raw data, gathered (lon,lat,country,category,type)
  
  ! raw data:
  type T_Emis_Gathered_All
    ! grid definition:
    integer                                 ::  nlon, nlat
    real, allocatable                       ::  lonbs(:,:)   ! (nlon,2)
    real, allocatable                       ::  latbs(:,:)   ! (nlat,2)
    ! country info:
    integer                                 ::  ncountry
    integer                                 ::  len_country_code ! actual length
    character(len=3), allocatable           ::  country_code(:)  ! country code (e.g. 'ALB', 'BEL')
    integer                                 ::  len_country_name ! actual length
    character(len=100), allocatable         ::  country_name(:)  ! country name
    ! emission categories
    integer                                 ::  ncat         ! number of emission categories
    character(len=8), allocatable           ::  cat_code(:)  ! code of each emission category
    integer                                 ::  len_cat_code ! actual maximum length of cat names
    character(len=100), allocatable         ::  cat_name(:)  ! name of each emission category 
    integer                                 ::  len_cat_name ! actual maximum length of cat names
    ! type info:
    integer                                 ::  ntyp
    character(len=2)                        ::  typs
    ! sources info:
    integer                                 ::  nsource
    integer, allocatable                    ::  icats(:)   ! (nsource)
    integer, allocatable                    ::  ityps(:)   ! (nsource)
    integer, allocatable                    ::  ilons(:)   ! (nsource)
    integer, allocatable                    ::  ilats(:)   ! (nsource)
    integer, allocatable                    ::  icnts(:)   ! (nsource)
    ! for point source use the exact location:
    real, allocatable                       ::  plons(:)   ! (nsource)
    real, allocatable                       ::  plats(:)   ! (nsource)
    ! number of emisison variables:
    integer                                 ::  nemis
    ! emitted components, units:
    character(len=LEN_NAME), allocatable    ::  name(:)   ! (nemis)
    character(len=LEN_UNITS), allocatable   ::  units(:)  ! (nemis)
    ! values:
    real, allocatable                       ::  emis(:,:)   ! (nsource,nemis)
  end type T_Emis_Gathered_All
  
  ! * country info
  
  type T_Country_Info
    ! ISO code:
    character(len=3)          ::  code
    ! long name:
    character(len=100)        ::  name
  end type T_Country_Info
  
  ! * source category info
  
  type T_SourceCat_Info
    ! name:
    character(len=100)        ::  name
    ! short code:
    character(len=8)          ::  code
  end type T_SourceCat_Info
  
  ! * sourcetype info
  
  type T_SourceType_Info
    ! single character code:
    character(len=1)          ::  code
  end type T_SourceType_Info
  
  ! * emission info
  
  type T_Emis_Info
    ! emitted component:
    character(len=LEN_NAME)   ::  name
    ! units:
    character(len=LEN_UNITS)  ::  units
  end type T_Emis_Info
  
  ! * base emisisons, gather (lon,lat,country)
  
  ! base emissions for single category/sourcetype combination:
  type T_Emis_Gathered
    ! current sector/fuel codes:
    integer                             ::  icat
    ! source type: area or point
    integer                             ::  ityp
    ! number of countries:
    integer                             ::  ncountry
    ! number of sources:
    integer                             ::  nsource
    ! x, y, and country indices in model grid:
    integer, allocatable                ::  ix(:)  ! (nsource)
    integer, allocatable                ::  iy(:)  ! (nsource)
    integer, allocatable                ::  ic(:)  ! (nsource)
    ! number of emisison variables:
    integer                             ::  nemis
    ! emissions:
    real, allocatable                   ::  emis(:,:)  ! (nsource,nemis)
  end type T_Emis_Gathered
  
  ! all data:
  type T_Emis_TNO_Base
    ! countries:
    integer                                 ::  ncountry
    type(T_Country_Info), allocatable       ::  cnti(:)  ! (ncountry)
    ! categories:
    integer                                 ::  ncat
    type(T_SourceCat_Info), allocatable     ::  cati(:)  ! (ntyp)
    ! source types:
    integer                                 ::  ntyp
    type(T_SourceType_Info), allocatable    ::  typi(:)  ! (ntyp)
    ! emitted tracers:
    integer                                 ::  nemis
    type(T_Emis_Info), allocatable          ::  emi(:)  ! (nemis)
    ! base emissions for single category/sourcetype combination,
    ! regridded to model grid
    type(T_Emis_Gathered), allocatable      ::  emg(:,:)   ! (ncat,ntyp)
  end type T_Emis_TNO_Base


contains


  ! ====================================================================
  ! ===
  ! === country info
  ! ===
  ! ====================================================================


  subroutine Country_Info_Init( cnti, code, name, status )

    ! --- in/out ------------------------------

    type(T_Country_Info), intent(inout)       ::  cnti
    character(len=*), intent(in)              ::  code
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Country_Info_Init'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! store:
    cnti%code = trim(code)
    cnti%name = trim(name)

    ! ok
    status = 0

  end subroutine Country_Info_Init
  
  
  ! ***
  
  
  subroutine Country_Info_Done( cnti, status )

    ! --- in/out ------------------------------

    type(T_Country_Info), intent(inout)    ::  cnti
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Country_Info_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------

    ! ok
    status = 0

  end subroutine Country_Info_Done


  ! ====================================================================
  ! ===
  ! === source category info
  ! ===
  ! ====================================================================


  subroutine SourceCat_Info_Init( cati, code, name, status )

    ! --- in/out ------------------------------

    type(T_SourceCat_Info), intent(inout)     ::  cati
    character(len=*), intent(in)              ::  code
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/SourceCat_Info_Init'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! store name:
    cati%name = trim(name)
    ! store code:
    cati%code = trim(code)

    ! ok
    status = 0

  end subroutine SourceCat_Info_Init
  
  
  ! ***
  
  
  subroutine SourceCat_Info_Done( cati, status )

    ! --- in/out ------------------------------

    type(T_SourceCat_Info), intent(inout)    ::  cati
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/SourceCat_Info_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------

    ! ok
    status = 0

  end subroutine SourceCat_Info_Done


  ! ====================================================================
  ! ===
  ! === source type info
  ! ===
  ! ====================================================================


  subroutine SourceType_Info_Init( typi, code, status )

    ! --- in/out ------------------------------

    type(T_SourceType_Info), intent(inout)    ::  typi
    character(len=1), intent(in)              ::  code
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/SourceType_Info_Init'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! store:
    typi%code  = code

    ! ok
    status = 0

  end subroutine SourceType_Info_Init
  
  
  ! ***
  
  
  subroutine SourceType_Info_Done( typi, status )

    ! --- in/out ------------------------------

    type(T_SourceType_Info), intent(inout)    ::  typi
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/SourceType_Info_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------

    ! ok
    status = 0

  end subroutine SourceType_Info_Done


  ! ====================================================================
  ! ===
  ! === emission info
  ! ===
  ! ====================================================================


  subroutine Emis_Info_Init( emi, name, units, status )

    ! --- in/out ------------------------------

    type(T_Emis_Info), intent(inout)    ::  emi
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  units
    integer, intent(out)                ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Emis_Info_Init'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! store:
    emi%name  = trim(name)
    emi%units = trim(units)

    ! ok
    status = 0

  end subroutine Emis_Info_Init
  
  
  ! ***
  
  
  subroutine Emis_Info_Done( emi, status )

    ! --- in/out ------------------------------

    type(T_Emis_Info), intent(inout)    ::  emi
    integer, intent(out)                ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Emis_Info_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------

    ! ok
    status = 0

  end subroutine Emis_Info_Done
  

  ! ====================================================================
  ! ===
  ! === gathered emissions
  ! ===
  ! ====================================================================

    
  subroutine Emis_Gathered_All_Init( ema, filename, status )

    use MDF, only : MDF_Open, MDF_Close
    use MDF, only : MDF_Inquire
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Inquire_Variable, MDF_Get_Var
    use MDF, only : MDF_Get_Att
    use MDF, only : MDF_NETCDF, MDF_READ, MDF_GLOBAL

    use LE_Emis_Tools, only : ShortSNAP_to_Code
    use LE_Emis_Tools, only : MDF_Get_StrArr

    ! --- in/out ------------------------------

    type(T_Emis_Gathered_All), intent(out)    ::  ema
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Emis_Gathered_All_Init'
    
    ! --- local -------------------------------
    
    integer                 ::  hid
    integer                 ::  dimid
    integer                 ::  varid
    integer, allocatable    ::  varids(:)
    integer                 ::  icat
    integer                 ::  iemis
    integer                 ::  iemis_pm25, iemis_pm10
    integer                 ::  nvar
    character(len=32)       ::  varname
    character(len=512)      ::  long_name
    
    character(len=8)        ::  format_id
    character(len=32)       ::  varname_longitude_bounds
    character(len=32)       ::  varname_latitude_bounds
    character(len=32)       ::  varname_country_name_len
    character(len=32)       ::  varname_emis_cat_name_len
    character(len=32)       ::  varname_emis_cat_name
    character(len=32)       ::  varname_emis_cat_shortsnap
    character(len=32)       ::  varname_emis_cat_code_len
    character(len=32)       ::  varname_emis_cat_code
    character(len=32)       ::  varname_longitude_source
    character(len=32)       ::  varname_latitude_source
    character(len=32)       ::  attname_units
    
    integer, allocatable    ::  shortsnap(:)

    ! --- begin -------------------------------

    ! open file:
    call MDF_Open( trim(filename), MDF_NETCDF, MDF_READ, hid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! file format: a format_id might be present ("1.3"),
    ! if not, assume that it is "1.2" :
    call MDF_Get_Att( hid, MDF_GLOBAL, 'format_id', format_id, status )
    if ( status /= 0 ) then
      format_id = '1.2'
      write (gol,'("IGNORE - attribute `format_id` not found, assume value `",a,"`")') trim(format_id); call goPr
    end if
    
    ! set some format specific values:
    select case ( trim(format_id) )
      !~ old:
      case ( '1.2' )
        varname_longitude_bounds   = 'grid_longitude_borders'
        varname_latitude_bounds    = 'grid_latitude_borders'
        varname_country_name_len   = 'country_name'
        varname_emis_cat_name_len  = 'emis_cat_name'
        varname_emis_cat_name      = 'emis_cat_name'
        varname_emis_cat_code_len  = ''
        varname_emis_cat_code      = ''
        varname_emis_cat_shortsnap = 'emis_cat_lesnap'
        varname_longitude_source   = 'longitude'
        varname_latitude_source    = 'latitude'
        attname_units              = 'unit'
      !~ latest
      case ( '1.3', '1.4' )
        varname_longitude_bounds   = 'longitude_bounds'
        varname_latitude_bounds    = 'latitude_bounds'
        varname_country_name_len   = 'country_name_len'
        varname_emis_cat_name_len  = 'emis_cat_name_len'
        varname_emis_cat_name      = 'emis_cat_name'
        varname_emis_cat_code_len  = ''
        varname_emis_cat_code      = ''
        varname_emis_cat_shortsnap = 'emis_cat_shortsnap'
        varname_longitude_source   = 'longitude_source'
        varname_latitude_source    = 'latitude_source'
        attname_units              = 'units'
      !~ emission model output
      case ( '1.5' )
        varname_longitude_bounds   = 'longitude_bounds'
        varname_latitude_bounds    = 'latitude_bounds'
        varname_country_name_len   = 'country_name_len'
        varname_emis_cat_name_len  = 'emis_cat_name_len'
        varname_emis_cat_name      = 'emis_sector_name'
        varname_emis_cat_code_len  = ''
        varname_emis_cat_code      = ''
        varname_emis_cat_shortsnap = 'sectorID'
        varname_longitude_source   = 'longitude_source'
        varname_latitude_source    = 'latitude_source'
        attname_units              = 'units'
      !~ latest
      case ( '1.6' )
        varname_longitude_bounds   = 'longitude_bounds'
        varname_latitude_bounds    = 'latitude_bounds'
        varname_country_name_len   = 'country_name_len'
        varname_emis_cat_name_len  = 'emis_cat_name_len'
        varname_emis_cat_name      = 'emis_cat_name'
        varname_emis_cat_code_len  = 'emis_cat_code_len'
        varname_emis_cat_code      = 'emis_cat_code'
        varname_emis_cat_shortsnap = ''
        varname_longitude_source   = 'longitude_source'
        varname_latitude_source    = 'latitude_source'
        attname_units              = 'units'
      !~ unknown ...
      case default
        write (gol,'("unsupported format_id `",a,"`")') trim(format_id); call goErr
        TRACEBACK; status=1; return
    end select
      

    !
    ! * read grid arrays
    !

    ! number of longitudes:
    call MDF_Inq_DimID( hid, 'longitude', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ema%nlon )
    IF_NOTOK_RETURN(status=1)
    ! number of latitudes:
    call MDF_Inq_DimID( hid, 'latitude', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ema%nlat )
    IF_NOTOK_RETURN(status=1)

    ! storage for grid cell boundaries:
    allocate( ema%lonbs(2,ema%nlon) )
    allocate( ema%latbs(2,ema%nlat) )
    ! fill grid:
    call MDF_Inq_VarID( hid, trim(varname_longitude_bounds), varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%lonbs, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inq_VarID( hid, trim(varname_latitude_bounds), varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%latbs, status )
    IF_NOTOK_RETURN(status=1)
    
    !
    ! * country arrays
    !

    ! number of countries:
    call MDF_Inq_DimID( hid, 'country', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ema%ncountry )
    IF_NOTOK_RETURN(status=1)
    ! length of country codes:
    ema%len_country_code = 3
    ! length of country names:
    call MDF_Inq_DimID( hid, trim(varname_country_name_len), dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ema%len_country_name )
    IF_NOTOK_RETURN(status=1)
    ! storage:
    allocate( ema%country_code(ema%ncountry) )
    allocate( ema%country_name(ema%ncountry) )
    ! read country codes:
    call MDF_Get_StrArr( hid, 'country_id'  , ema%len_country_code, ema%ncountry, ema%country_code, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_StrArr( hid, 'country_name', ema%len_country_name, ema%ncountry, ema%country_name, status )
    IF_NOTOK_RETURN(status=1)
    
    
    !
    ! * source category arrays
    !
    
    ! number of categories:
    call MDF_Inq_DimID( hid, 'emis_cat', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ema%ncat )
    IF_NOTOK_RETURN(status=1)
    ! maximum length of name:
    call MDF_Inq_DimID( hid, trim(varname_emis_cat_name_len), dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ema%len_cat_name )
    IF_NOTOK_RETURN(status=1)

    ! storage:
    allocate( ema%cat_code(ema%ncat) )
    allocate( ema%cat_name(ema%ncat) )

    ! read emission category codes
    if ( len_trim(varname_emis_cat_shortsnap) > 0 ) then
      !~ LE SNAP
      allocate( shortsnap(ema%ncat), stat=status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Inq_VarID( hid, trim(varname_emis_cat_shortsnap), varid, status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Get_Var( hid, varid, shortsnap, status )
      IF_NOTOK_RETURN(status=1)
      do icat = 1, ema%ncat
        call ShortSNAP_to_Code( shortsnap(icat), ema%cat_code(icat), status )
        IF_NOTOK_RETURN(status=1)
      end do
      deallocate( shortsnap, stat=status )
      IF_NOTOK_RETURN(status=1)
    !  
    else if ( len_trim(varname_emis_cat_code) > 0 ) then
      ! maximum length of name:
      call MDF_Inq_DimID( hid, trim(varname_emis_cat_code_len), dimid, status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Inquire_Dimension( hid, dimid, status, length=ema%len_cat_code )
      IF_NOTOK_RETURN(status=1)
      ! character codes:
      call MDF_Get_StrArr( hid, trim(varname_emis_cat_code), ema%len_cat_code, ema%ncat, ema%cat_code, status )
      IF_NOTOK_RETURN(status=1)
    !
    else
      write (gol,'("no varname for cat code defined")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! read emission category names:
    call MDF_Get_StrArr( hid, trim(varname_emis_cat_name), ema%len_cat_name, ema%ncat, ema%cat_name, status )
    IF_NOTOK_RETURN(status=1)
    ! adhoc fixes ...
    do icat = 1, ema%ncat
      select case ( ema%cat_name(icat) )
        case ( '"Residential,_commercial_and_other_combustion"' )
          ema%cat_name(icat) = 'Residential,_commercial_and_other_combustion'
      end select
    end do

    !
    ! * source type arrays
    !

    ! number of source types:
    call MDF_Inq_DimID( hid, 'source_type', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ema%ntyp )
    IF_NOTOK_RETURN(status=1)

    ! source types character codes ('p' for point, 'a' for area):
    call MDF_Inq_VarID( hid, 'source_type_code', varid, status )
    IF_NOTOK_RETURN(status=1)
    ! read into character array:
    call MDF_Get_Var( hid, varid, ema%typs(1:ema%ntyp), status )
    IF_NOTOK_RETURN(status=1)

    !
    ! * read source arrays
    !

    ! number of sources:
    call MDF_Inq_DimID( hid, 'source', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ema%nsource )
    IF_NOTOK_RETURN(status=1)
    ! storage:
    allocate( ema%icats(ema%nsource) )
    allocate( ema%ityps(ema%nsource) )
    allocate( ema%ilons(ema%nsource) )
    allocate( ema%ilats(ema%nsource) )
    allocate( ema%icnts(ema%nsource) )
    ! indices in categories array:
    call MDF_Inq_VarID( hid, 'emission_category_index', varid , status)
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%icats, status )
    IF_NOTOK_RETURN(status=1)
    ! indices in sourcetypes array:
    call MDF_Inq_VarID( hid, 'source_type_index', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%ityps, status )
    IF_NOTOK_RETURN(status=1)
    ! indices in lons array:
    call MDF_Inq_VarID( hid, 'longitude_index', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%ilons, status )
    IF_NOTOK_RETURN(status=1)
    ! indices in lats array:
    call MDF_Inq_VarID( hid, 'latitude_index', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%ilats, status )
    IF_NOTOK_RETURN(status=1)
    ! indices in countries array:
    call MDF_Inq_VarID( hid, 'country_index', varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%icnts, status )
    IF_NOTOK_RETURN(status=1)
    
    ! check ...
    if ( any(ema%ityps < 1) .or. any(ema%ityps > ema%ntyp) ) then
      write (gol,'("found source type indices out of range 1 .. ",i4)') ema%ntyp; call goErr
      write (gol,*) '  min, max : ', minval(ema%ityps), maxval(ema%ityps); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( ema%plons(ema%nsource) )
    allocate( ema%plats(ema%nsource) )
    ! longitudes (used for point sources):
    call MDF_Inq_VarID( hid, trim(varname_longitude_source), varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%plons, status )
    IF_NOTOK_RETURN(status=1)
    ! latitudes (used for point sources):
    call MDF_Inq_VarID( hid, (varname_latitude_source), varid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Get_Var( hid, varid, ema%plats, status )
    IF_NOTOK_RETURN(status=1)
    
    !
    ! * read emissions
    !

    ! info ...
    write (gol,'("  Search emission variables: ")'); call goPr
    ! emitted components; init counter:
    ema%nemis = 0
    ! get total number of data sets:
    call MDF_Inquire( hid, status, nVariables=nvar )
    IF_NOTOK_RETURN(status=1)
    ! storage for variable id's:
    allocate( varids(nvar) )
    ! loop over all variables:
    do varid = 1, nvar
      ! extract name:
      call MDF_Inquire_Variable( hid, varid, status, name=varname )
      IF_NOTOK_RETURN(status=1)
      ! try to get a the long_name attribute, might fail:
      call MDF_Get_Att( hid, varid, 'long_name', long_name, status )
      if (status/=0) cycle
      ! check if it is long enough to hold 'emission of ...'
      if ( len_trim(long_name) < 12 ) cycle
      ! compare:
      if ( long_name(1:12) == 'emission of ' ) then
        ! info ...
        write (gol,'("    found `",a,"` (",a,")")') trim(varname), trim(long_name); call goPr
        ! increase counter:
        ema%nemis = ema%nemis + 1
        ! store:
        varids(ema%nemis) = varid
      end if  ! long name says it is an emission ...
    end do
    ! check ...
    if ( ema%nemis == 0 ) then
      write (gol,'("could not find emission variables in: ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage for emissions:
    allocate( ema%name (ema%nemis) )
    allocate( ema%units(ema%nemis) )
    allocate( ema%emis(ema%nsource,ema%nemis) )
    ! loop over emissions:
    do iemis = 1, ema%nemis
      ! current:
      varid = varids(iemis)
      ! extract name:
      call MDF_Inquire_Variable( hid, varid, status, name=ema%name(iemis) )
      IF_NOTOK_RETURN(status=1)
      ! units are stored as attribute:
      call MDF_Get_Att( hid, varid, trim(attname_units), ema%units(iemis), status )
      IF_NOTOK_RETURN(status=1)
      ! read record:
      call MDF_Get_Var( hid, varid, ema%emis(:,iemis), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! clear:
    deallocate( varids )

    !
    ! * convert total pm10 (0-10um) to coarse fraction (2.5-10um)
    !
    
    ! search for index of coarse mode emission:
    iemis_pm10 = -1
    do iemis = 1, ema%nemis
      select case ( trim(ema%name(iemis)) )
        case ( 'PPM10', 'pm10', 'pmco' )
          iemis_pm10 = iemis
          exit
        case ( 'ch4', 'co', 'co2', 'nh3', 'nmvoc', 'nox', 'pm25', 'pm2_5', 'so2', &
               'pm1', 'bc_1', 'ec_1-25', 'ec_25-10', 'oc_25', 'oc_25-10', &
               'so4ks', 'so4as', &
               'bcks', 'bcas', 'bcki', 'ocks', 'ocas', 'ocki', &
               'duas', 'duai', &
               'nais', 'nacs', 'naii', 'naci', &
               'pb', 'cd', 'cu_scaled_ex_bw' )
          ! not pm10
        case default
          write (gol,'("do not know whether tracer `",a,"` represents coarse mode aerosol or something else")') trim(ema%name(iemis)); call goErr
          TRACEBACK; status=1; return
      end select
    end do
    ! found ? then continue ...
    if ( iemis_pm10 > 0 ) then
      ! search for index of fine mode emission:
      iemis_pm25 = -1
      do iemis = 1, ema%nemis
        select case ( trim(ema%name(iemis)) )
          case ( 'PPM25', 'pm2_5', 'pm25' )
            iemis_pm25 = iemis
            exit
        end select
      end do
      ! check ...
      if ( iemis_pm25 < 0 ) then
        write (gol,'("could not find fine mode pm in emitted tracers :")'); call goErr
        do iemis = 1, ema%nemis
          write (gol,'("  ",a)') trim(ema%name(iemis)); call goErr
        end do
        TRACEBACK; status=1; return
      end if
      ! info ...
      write (gol,'("  Convert `",a,"` to `pm25_pm10` using `",a,"`")') &
              trim(ema%name(iemis_pm10)), trim(ema%name(iemis_pm25)); call goPr
      ! substract:  [2.5,10]          [0,10]          -         [0,2.5]
      ema%emis(:,iemis_pm10) = ema%emis(:,iemis_pm10) - ema%emis(:,iemis_pm25)
      ! change name:
      ema%name(iemis_pm10) = 'pm25_pm10'
    end if

    !
    ! * close
    !

    ! close file:
    call MDF_close( hid, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! * end
    !

    ! ok
    status = 0

  end subroutine Emis_Gathered_All_Init
  
  
  ! ***
  
  
  subroutine Emis_Gathered_All_Done( ema, status )

    ! --- in/out ------------------------------

    type(T_Emis_Gathered_All), intent(inout)  ::  ema
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Emis_Gathered_All_Done'
    
    ! --- local -------------------------------
       
    ! --- begin -------------------------------

    ! clear grid:
    deallocate( ema%lonbs )
    deallocate( ema%latbs )
    ! clear countries:
    deallocate( ema%country_code )
    deallocate( ema%country_name )
    ! clear categories:
    deallocate( ema%cat_code )
    deallocate( ema%cat_name )
    ! clear source info:
    deallocate( ema%icats )
    deallocate( ema%ityps )
    deallocate( ema%ilons )
    deallocate( ema%ilats )
    deallocate( ema%icnts )
    ! clear point locations:
    deallocate( ema%plons )
    deallocate( ema%plats )
    ! clear emisisons:
    deallocate( ema%name  )
    deallocate( ema%units )
    deallocate( ema%emis  )

    ! ok
    status = 0

  end subroutine Emis_Gathered_All_Done
  
  
  ! ====================================================================
  ! ===
  ! === collection of regridded emissions
  ! ===
  ! ====================================================================


  subroutine Emis_Gathered_Init( emg, ema, icat, ityp, status )
  
    use Dims   , only : nx, ny, nz
    use LE_Grid, only : ugg

    ! --- in/out ------------------------------

    type(T_Emis_Gathered), intent(out)      ::  emg
    type(T_Emis_Gathered_All), intent(in)   ::  ema
    integer, intent(in)                     ::  icat
    integer, intent(in)                     ::  ityp
    integer, intent(out)                    ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Emis_Gathered_Init'
    
    ! --- local -------------------------------

    integer, allocatable          ::  indx(:,:,:)   ! (nx,ny,ncountry)
    integer                       ::  iloop
    integer                       ::  isrc
    integer                       ::  ilon, ilat, icnt
    real                          ::  west, east, south, north
    real                          ::  lon, lat
    integer                       ::  i1, i2, j1, j2
    real, pointer                 ::  ff(:,:)
    integer                       ::  ix, iy
    integer                       ::  isource
    integer                       ::  iemis
    integer                       ::  n

    ! --- begin -------------------------------
    
    ! fill category codes index:
    emg%icat = icat
    ! fill source type code:
    emg%ityp = ityp
    
    ! number of countries:
    emg%ncountry = ema%ncountry

    ! number of emitted tracers:
    emg%nemis = ema%nemis
    
    ! array to map x/y/country tupples to 1D index:
    allocate( indx(nx,ny,emg%ncountry) )
    
    ! init weights:
    nullify( ff )

    ! no sources filled yet:
    emg%nsource = 0

    ! no 1D indices yet:
    indx = 0

    ! fill in 2 loops:
    !   loop 1 : count the number of sources, allocate
    !   loop 2 : fill
    do iloop = 1, 2

      ! loop over all original sources:
      do isrc = 1, ema%nsource

        ! not this category or type ? then skip:
        if ( ema%icats(isrc) /= icat ) cycle
        if ( ema%ityps(isrc) /= ityp ) cycle

        ! original indices:
        ilon = ema%ilons(isrc)
        ilat = ema%ilats(isrc)
        icnt = ema%icnts(isrc)

        ! distribute over model grid cells;
        ! depends on source type:
        select case ( ema%typs(ityp:ityp) )
          !~~~~~~~~~~~~~~~
          case ( 'A', 'a' )
          !~~~~~~~~~~~~~~~
            ! extract bounding box:
            west  = ema%lonbs(1,ilon)
            east  = ema%lonbs(2,ilon)
            south = ema%latbs(1,ilat)
            north = ema%latbs(2,ilat)
            ! distribution of area over target grid:
            call ugg%GetDistribution( west, east, south, north, &
                                    i1, i2, j1, j2, ff, status )
            IF_ERROR_RETURN(status=1)
            ! not in domain ? then skip:
            if ( status < 0 ) cycle
          !~~~~~~~~~~~~~~~
          case ( 'P', 'p' )
          !~~~~~~~~~~~~~~~
            ! extact point:
            lon = ema%plons(isrc)
            lat = ema%plats(isrc)
            ! grid cell:
            call ugg%GetLocation( lon, lat, i1, j1, status )
            IF_ERROR_RETURN(status=1)
            ! not in domain ? then skip:
            if ( status < 0 ) cycle
            ! fill up to full distribution: assign all to single cell:
            i2 = i1
            j2 = j1
            if ( associated(ff) ) deallocate( ff )
            allocate( ff(i1:i2,j1:j2) )
            ff = 1.0
          !~~~~~~~~~~~~~~~
          case default
          !~~~~~~~~~~~~~~~
            write (gol,'("unsupported source type : ",a)') ema%typs(ityp:ityp); call goErr
            TRACEBACK; status=1; return
        end select
        ! loop over target cells:
        do iy = j1, j2
          do ix = i1, i2
            !~~ counting loop ?
            if ( iloop == 1 ) then
              ! reset maping index, used as a flag later on:
              indx(ix,iy,icnt) = 1
            !~~ fill loop
            else
              ! source index:
              isource = indx(ix,iy,icnt)
              ! store:
              emg%ix(isource) = ix
              emg%iy(isource) = iy
              emg%ic(isource) = icnt
              ! loop over emitted tracers:
              do iemis = 1, ema%nemis
                ! add emission contributions to this cell:
                emg%emis(isource,iemis) = emg%emis(isource,iemis) &
                                  + ema%emis(isrc,iemis) * ff(ix,iy)
              end do  ! tracers
            end if
            !~~
          end do  ! i
        end do  ! j

      end do   ! original sources

      !~~ counting loop ?
      if ( iloop == 1 ) then
        ! fill source indices:
        n = 0
        do icnt = 1, emg%ncountry
          do iy = 1, ny
            do ix = 1, nx
              ! used ?
              if ( indx(ix,iy,icnt) > 0 ) then
                ! incease index:
                n = n + 1
                ! store:
                indx(ix,iy,icnt) = n
              end if  ! coordinates in use
            end do   ! ix
          end do  ! iy
        end do  ! country
        ! store maximum:
        emg%nsource = n
        ! any sources found ?
        if ( n > 0 ) then
          ! info ...
          write (gol,'("    gathered sources for cat ",i2," typ ",i1," : ",i8)') icat, ityp, n; call goPr
          ! storage for coordinate indices:
          allocate( emg%ix(n) )
          allocate( emg%iy(n) )
          allocate( emg%ic(n) )
          ! storage for emissions:
          allocate( emg%emis(n,ema%nemis) )
          ! set to zero, contributions will be added:
          emg%emis = 0.0
        end if
      end if
      !~~

    end do   ! loop
    
    ! clear:
    if ( associated(ff) ) deallocate( ff )

    ! ok
    status = 0
    
  end subroutine Emis_Gathered_Init


  ! ***
  

  subroutine Emis_Gathered_Done( emg, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Gathered), intent(inout)    ::  emg
    integer, intent(out)                    ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Emis_Gathered_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! any sources ?
    if ( emg%nsource > 0 ) then
      ! clear:
      deallocate( emg%ix )
      deallocate( emg%iy )
      deallocate( emg%ic )
      deallocate( emg%emis )
    end if

    ! ok
    status = 0

  end subroutine Emis_Gathered_Done
  

  ! ====================================================================
  ! ===
  ! === collection of regridded emissions
  ! ===
  ! ====================================================================


  subroutine LE_Emis_TNO_Base_Init( emb, filename, status )
  
    use Grid   , only : GetDistribution, GetLocation
    use Dims   , only : nx, ny, nz

    ! --- in/out ------------------------------

    type(T_Emis_TNO_Base), intent(out)  ::  emb    
    character(len=*), intent(in)        ::  filename
    integer, intent(out)                ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Base_Init'
    
    ! --- local -------------------------------

    type(T_Emis_Gathered_All)     ::  ema 
    integer                       ::  icnt
    integer                       ::  icat
    integer                       ::  ityp
    integer                       ::  iemis
    
    ! --- begin -------------------------------
    
    ! read all data:
    call Emis_Gathered_All_Init( ema, filename, status )
    IF_NOTOK_RETURN(status=1)
    
    ! number of countries:
    emb%ncountry = ema%ncountry
    ! storage for source country info:
    allocate( emb%cnti(ema%ncountry) )
    ! loop over emitted tracers:
    do icnt = 1, emb%ncountry
      ! init info:
      call Country_Info_Init( emb%cnti(icnt), ema%country_code(icnt), &
                                    ema%country_name(icnt), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! number of categories:
    emb%ncat = ema%ncat
    ! storage for source category info:
    allocate( emb%cati(ema%ncat) )
    ! loop over categories:
    do icat = 1, emb%ncat
      ! init info:
      call SourceCat_Info_Init( emb%cati(icat), ema%cat_code(icat), ema%cat_name(icat), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! number of source types:
    emb%ntyp = ema%ntyp
    ! storage for source type info:
    allocate( emb%typi(ema%ntyp) )
    ! loop over source types:
    do ityp = 1, emb%ntyp
      ! init info:
      call SourceType_Info_Init( emb%typi(ityp), ema%typs(ityp:ityp), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! storage for gathered emisisons:
    allocate( emb%emg(emb%ncat,emb%ntyp) )

    ! number of emitted tracers:
    emb%nemis = ema%nemis
    ! storage for emission info:
    allocate( emb%emi(ema%nemis) )
    ! loop over emitted tracers:
    do iemis = 1, emb%nemis
      ! init info:
      call Emis_Info_Init( emb%emi(iemis), trim(ema%name(iemis)), trim(ema%units(iemis)), status )
      IF_NOTOK_RETURN(status=1)
    end do   

    ! init per base emissions;
    ! loop over source types:
    do ityp = 1, emb%ntyp
      ! loop over categories:
      do icat = 1, emb%ncat

        ! gather and regrid for this category/type:
        call Emis_Gathered_Init( emb%emg(icat,ityp), ema, icat, ityp, status )
        IF_NOTOK_RETURN(status=1)

      end do  ! categories
    end do  ! source types
    
    ! clear input data:
    call Emis_Gathered_All_Done( ema, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_TNO_Base_Init


  ! ***
  

  subroutine LE_Emis_TNO_Base_Done( emb, status )
      
    ! --- in/out ------------------------------
    
    type(T_Emis_TNO_Base), intent(inout)  ::  emb
    integer, intent(out)                  ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Base_Done'
    
    ! --- local -------------------------------
    
    integer     ::  icnt
    integer     ::  icat
    integer     ::  ityp
    integer     ::  iemis
    
    ! --- begin -------------------------------
    
    ! loop over countries:
    do icnt = 1, emb%ncountry
      ! clear:
      call Country_Info_Done( emb%cnti(icnt), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! loop over source types:
    do ityp = 1, emb%ntyp
      ! clear:
      call SourceType_Info_Done( emb%typi(ityp), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! loop over emitted tracers:
    do iemis = 1, emb%nemis
      ! clear:
      call Emis_Info_Done( emb%emi(iemis), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! loop over source types:
    do ityp = 1, emb%ntyp
      ! loop over categories:
      do icat = 1, emb%ncat
        ! clear gathered emissions:
        call Emis_Gathered_Done( emb%emg(icat,ityp), status )
        IF_NOTOK_RETURN(status=1)
      end do  ! categories
    end do  ! types
    ! clear:
    deallocate( emb%emg )

    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Base_Done
  

  ! ***
  
  ! write summary tables with emisisons totals
  
  subroutine LE_Emis_TNO_Base_Summary( emb, basename, status )
  
    use GO, only : goGetFU

    ! --- in/out ---------------------------------
    
    type(T_Emis_TNO_Base), intent(inout)    ::  emb
    character(len=*), intent(in)            ::  basename
    integer, intent(out)                    ::  status
    
    !integer, intent(in)    :: lines
    !real   , intent(in)    :: fac_unit ! factor for unit conversion
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Base_Summary'
    
    ! --- local ----------------------------------

    real, allocatable      ::  emission_tab(:,:,:)
    integer                ::  isrc
    integer                ::  icountry, icat, ityp, iemis
    integer                ::  fu
    character(len=512)     ::  fname
    character(len=1024)    ::  line
    character(len=32)      ::  field
    character(len=1)       ::  sep
    
    ! --- begin ----------------------------------

    ! storage:
    allocate( emission_tab(emb%ncountry,emb%ncat,emb%nemis) )

    ! add emission for each record to totals for each country 
    ! and each enission category
    emission_tab(:,:,:) = 0.0
    ! loop over categories:
    do icat = 1, emb%ncat
      ! loop over source types:
      do ityp = 1, emb%ntyp
        ! loop over all sources:
        do isrc = 1, emb%emg(icat,ityp)%nsource
          ! get the country number
          icountry = emb%emg(icat,ityp)%ic(isrc)
          ! loop over emitted tracers:
          do iemis = 1, emb%nemis
            ! add contribution:
            emission_tab(icountry,icat,iemis) = &
               emission_tab(icountry,icat,iemis) + &
               emb%emg(icat,ityp)%emis(isrc,iemis)
          end do
        end do  ! sources
      end do  ! source types
    end do  ! categories
    
    ! seperation character:
    sep = ';'
    
    ! * total per countriy and category
    
    ! output file:
    write (fname,'(a,"_country_category.csv")') trim(basename)
    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open:
    open( unit=fu, file=trim(fname), form='formatted', iostat=status )
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(5a)') 'country', &
                            sep, 'sector', sep, 'subsector'
    do iemis = 1, emb%nemis
      field = trim(emb%emi(iemis)%name)
      line = trim(line)//sep//trim(field)
    end do
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(5a)') 'iso3', sep, 'index', sep, 'index'
    do iemis = 1, emb%nemis
      field = trim(emb%emi(iemis)%units)
      line = trim(line)//sep//trim(field)
    end do
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! loop over countries:
    do icountry = 1, emb%ncountry
      ! loop over categories:
      do icat = 1, emb%ncat
        ! write record:
        write (line,'(3a)') trim(emb%cnti(icountry)%code), sep, emb%cati(icat)%code
        do iemis = 1, emb%nemis
          write (field,'(e12.4)') emission_tab(icountry,icat,iemis)
          line = trim(line)//sep//trim(field)
        end do
        write (fu,'(a)',iostat=status) trim(line)
        IF_NOTOK_RETURN(status=1)
      end do  ! categories
    end do  ! countries
    ! close:
    close( unit=fu, iostat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! * total per country
    
    ! output file:
    write (fname,'(a,"_country.csv")') trim(basename)
    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open:
    open( unit=fu, file=trim(fname), form='formatted', iostat=status )
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(a)') 'country'
    do iemis = 1, emb%nemis
      field = trim(emb%emi(iemis)%name)
      line = trim(line)//sep//trim(field)
    end do
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(a)') 'iso3'
    do iemis = 1, emb%nemis
      field = trim(emb%emi(iemis)%units)
      line = trim(line)//sep//trim(field)
    end do
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! loop over countries:
    do icountry = 1, emb%ncountry
      ! write record:
      write (line,'(a,4(a,i2.2))') trim(emb%cnti(icountry)%code)
      do iemis = 1, emb%nemis
        write (field,'(e12.4)') sum(emission_tab(icountry,:,iemis))
        line = trim(line)//sep//trim(field)
      end do
      write (fu,'(a)',iostat=status) trim(line)
      IF_NOTOK_RETURN(status=1)
    end do  ! countries
    ! close:
    close( unit=fu, iostat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! * country info
    
    ! output file:
    write (fname,'(a,"_countries.csv")') trim(basename)
    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open:
    open( unit=fu, file=trim(fname), form='formatted', iostat=status )
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(3a)') 'code', sep, 'name'
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! loop over countries:
    do icountry = 1, emb%ncountry
      ! write record:
      write (line,'(3a)') trim(emb%cnti(icountry)%code), sep, trim(emb%cnti(icountry)%name)
      write (fu,'(a)',iostat=status) trim(line)
      IF_NOTOK_RETURN(status=1)
    end do  ! countries
    ! close:
    close( unit=fu, iostat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! * category info
    
    ! output file:
    write (fname,'(a,"_categories.csv")') trim(basename)
    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open:
    open( unit=fu, file=trim(fname), form='formatted', iostat=status )
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(5a)') 'sector', sep, 'subsector', &
                          sep, 'category'
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! loop over categories:
    do icat = 1, emb%ncat
      ! write record:
      write (line,'(3a)') emb%cati(icat)%code, sep, trim(emb%cati(icat)%name)
      write (fu,'(a)',iostat=status) trim(line)
      IF_NOTOK_RETURN(status=1)
    end do  ! categories
    ! close:
    close( unit=fu, iostat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! *

    ! clear:
    deallocate( emission_tab )

    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Base_Summary


end module LE_Emis_TNO_Base
