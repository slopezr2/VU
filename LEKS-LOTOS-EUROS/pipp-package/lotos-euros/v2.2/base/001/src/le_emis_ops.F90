!###############################################################################
!
! NAME
!   LE_Emis_OPS - routines to read and work with emission files from OPS-input
!
! EMISSION FILES
!
! NOx_OPS_bronbestand.brn 
!
! x,y locations of the sources (either rdm or lon/lat)
!
! snr    x(m)    y(m)    q(g/s) hc(MW)  h(m)   r(m)  s(m)  dv cat area ps  component
!2011  22.080  44.130 2.835E+00 20.000 100.0      0  50.0   1  04 891   0  NOX         
!2011  20.312  44.819 7.824E-02  0.500   4.0      0   2.0   3  08 891   0  NOX         
!2011  19.000  42.300 3.371E-03 20.000 100.0      0  50.0   1  09 891   0  NOX         
!2011  19.650  45.817 1.054E-02 20.000 100.0      0  50.0   1  09 891   0  NOX         
!2011  21.800  44.500 1.251E-03 20.000 100.0      0  50.0   1  09 891   0  NOX         
!2011  22.200  44.400 3.694E-03 20.000 100.0      0  50.0   1  07 891   0  NOX         
!2011  18.750  40.125 1.670E-01  0.000  10.0  91398   5.0   2  02   8   0  NOX         
!2011  18.750  40.875 1.045E+00  0.000  10.0  90888   5.0   2  02   8   0  NOX         
!2011  20.250  39.375 6.699E-02  0.000  10.0  91897   5.0   2  02   8   0  NOX         
!2011  20.250  40.125 3.126E+00  0.000  10.0  91398   5.0   2  02   8   0  NOX        
!2011   8.938  52.094 1.027E-01 10.000  50.0   6828  25.0   1  03 276   0  NOX         
!2011   8.938  52.156 1.171E-01 10.000  50.0   6823  25.0   1  03 276   0  NOX         
!2011   8.938  52.531 3.047E-02 10.000  50.0   6794  25.0   1  03 276   0  NOX         
!2011   8.938  52.594 3.079E-02 10.000  50.0   6789  25.0   1  03 276   0  NOX
!2011 -597364  176204 3.799E-02  2.100  30.0   5000  15.0   03840 9993840           NOx  

! The text files should have comment explaining the format and content.
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

module LE_Emis_OPS

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile

  use LE_Emis_HeightDistr, only : T_Emis_HeightDistr
  use LE_Emis_Time_Prof  , only : T_Emis_Time_Prof
  use LE_Emis_Composition, only : T_Emis_Composition
  
!#ifdef with_netcdf
!  use NetCDF
!#endif

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Emis_OPS

  public  ::  LE_Emis_OPS_Init, LE_Emis_OPS_Done
  public  ::  LE_Emis_OPS_Setup


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_OPS'

  ! maximum number of emitted tracers:
  integer, parameter            ::  max_pollutant = 10
  ! maximum number of sources per tracer
  integer, parameter            ::  max_source=9000000
  ! maximum number of daily patterns
  integer, parameter            ::  max_pattern = 40
  ! maximum number of model tracers composited from pollutant
  integer, parameter            ::  max_spec = 10
  ! maximum number of countries defined
  integer, parameter            ::  max_country = 70
  ! maximum number of SNAP sectors defined
  integer, parameter            ::  max_cat = 20
  ! maximum number of defined time types
  integer, parameter            ::  max_time_type = 5
  
  
!  ! maximum number of defined time types
!  integer, parameter                ::  max_time_type = 5
!
!  ! for temporary files that are closed directly after use
!  integer, parameter :: u_tmp  = 61

  include "OPS_PS_to_SNAP_tables.inc"
  
  include "OPS_Change_PS_Code.inc"
  
  ! --- local --------------------------------
  
  ! emission properties per pollutant
  type T_Emis_OPS_Pol
    ! name of the pollutant
    character(len=32)                 ::  name
    ! number of sources for this pollutant
    integer                           ::  nsource
    ! characteristics of the emissions
    real, allocatable                 ::  lon(:)  ! degree_east
    real, allocatable                 ::  lat(:)  ! degree_north
    integer, allocatable              ::  ilon(:) ! lon_index in LE-grid
    integer, allocatable              ::  ilat(:) ! lat_index in LE-grid
    real, allocatable                 ::  Q(:)    ! Source strength
    real, allocatable                 ::  hc(:)   ! Heat in the source
    real, allocatable                 ::  h(:)    ! source height
    real, allocatable                 ::  s(:)    ! vertical spread of area sources
    integer, allocatable              ::  r(:)    ! source size (assumed square) 
    integer, allocatable              ::  dv(:)   ! temporal variation code
    integer, allocatable              ::  icat(:) ! index of source sector
    integer, allocatable              ::  area(:) ! area code of each source 
    character(len=3), allocatable     ::  country(:) ! country of source
    integer, allocatable              ::  icountry(:) ! index of country code
    integer, allocatable              ::  ps(:)   ! code for particle distribution
    
    ! how to distirbute pollutants over tracers
    integer                           ::  nspec
    integer                           ::  ispecs(max_spec)
    real                              ::  compfrac(max_spec)
    
    ! scenario factor of the pollutant
    real                              ::  scenfac
    
    ! emission composition
    type(T_Emis_Composition)          ::  emcomp
    
    ! unit of emission
    character(len=10)                 ::  unit
    
    ! molecular mass per pollutant
    real                              ::  xm
  end type T_Emis_OPS_Pol

  ! emission data base
  type T_Emis_OPS

    ! label assigned to this emission:
    character(len=32)                 :: label

    ! year for which this data set was made:
    ! (in future, implement this with a time range [t1,t2])
    integer                           ::  year
    ! settings:
    character(len=512)                ::  rcfile
    character(len=64)                 ::  rckey
   
    ! arrays for country code:
    integer                           ::  ncountry         ! number of countries
    integer                           ::  country_code(max_country)
    character(len=3)                  ::  country_ISO3(max_country)  ! country code (e.g. 'ALB', 'BEL')
    character(len=100)                ::  country_name(max_country)  ! country name

    ! height distribution
    logical                           ::  with_height_distribution
    type(T_Emis_HeightDistr)          ::  hdistr
    
    ! time factors country/snap dependent
    type(T_Emis_Time_Prof)            ::  time_prof


    ! description of time factors to be used:
    !   'unity' | 'mdh' | 'time_prof' | 'hourly-gridded'
    integer                           ::  ntimetype
    character(len=16)                 ::  timetype(max_time_type)
    ! emission distribution:
    character(len=32)                 ::  distr_nam   ! 'SNAP15'

    ! emission categories
    integer                           ::  ncat         ! number of emission categories
    character(len=8)                  ::  cat_code(max_cat)  ! code-number of each emission category
    character(len=100)                ::  cat_nam(max_cat)   ! name of each emission category

    ! number of emitted components:
    integer                           ::  npol
    character(len=10)                 ::  pol_names(max_pollutant)
    ! names:
    type(T_Emis_OPS_Pol), allocatable ::  pol(:)
    
    ! Composition of PM in 6 classes
    real                              ::  pm_comp(max_country,max_cat, 6)
    ! fraction of ec in pm25
    real, allocatable                 ::  ec_frac(:,:)
    
    ! temporal variations
    real                              ::  daily_pattern(0:max_pattern,0:23)
    
    ! country/SNAP totals
    real, allocatable                 ::  country_snap_total(:,:,:)  ! pollutant,country,snap
    
    real, allocatable                 ::  annual_total_gridded(:,:,:) ! pollutant, lon,lat

  end type T_Emis_OPS
  


contains


  ! ===============================================================


  subroutine LE_Emis_OPS_Init( emo, label, rcF, rckey, t, status )

    use GO, only : TrcFile
    use GO, only : GO_Timer_Def
    use GO, only : TDate
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    type(T_Emis_OPS), intent(out)       ::  emo
    character(len=*), intent(in)        ::  label
    type(TrcFile), intent(in)           ::  rcF
    character(len=*), intent(in)        ::  rckey
    type(TDate), intent(in)             ::  t       ! start time
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_Init'

    ! --- local ------------------------------

    ! --- begin -------------------------------

    ! no data Init yet ...
    emo%year = -1

    ! store label:
    emo%label = trim(label)

    ! store info on settings:
    emo%rcfile = trim(rcF%fname)
    emo%rckey = trim(rckey)

    ! init
    emo%country_code = 0
    emo%cat_code     = ''
    
    ! read for current year:
    call LE_Emis_OPS_Read( emo, t%year, status )
    IF_NOTOK_RETURN(status=1)

    ! enable data:
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lon', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lat', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_OPS_Init


  ! ***


  subroutine LE_Emis_OPS_Done( emo, status )

    ! --- in/out ------------------------------

    type(T_Emis_OPS), intent(inout)    ::  emo
    integer, intent(out)               ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! any data loaded ?
    if ( emo%year > 0 ) then
      ! cleanup:
      call LE_Emis_OPS_Clear( emo, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine LE_Emis_OPS_Done


  ! ===============================================================


  subroutine LE_Emis_OPS_Read( emo, year, status )

    use Binas          , only : ae, deg2rad
    use Binas          , only : xm_N, xm_O, xm_S, xm_H, xm_dummy
    use GO             , only : ReadRc, Init, Done, ReadRc
    use GO             , only : goLoCase, goMatchValue, goReadFromLine
    use GO             , only : goReplace, goVarValue, goSplitString
    use GO             , only : goGetFU
    use LE_grid        , only : ugg
    
!    use LE_Logging     , only : u_log, u_err
    use indices        , only : nspec, specname
    use indices        , only : i_no, i_no2, i_so2, i_nh3 

    use LE_Emis_Time_Prof  , only : LE_Emis_Time_Prof_Init
    use LE_Emis_OPS_PM_Composition, only : LE_Emis_OPS_PM_Composition_Init
    use LE_Emis_OPS_EC_Fractions, only : LE_Emis_OPS_EC_Fractions_Init
    use LE_Emis_OPS_Output, only : LE_Emis_OPS_Write_Annual_Total
    
    use LE_Data  , only : LE_Data_GetPointer
    use LE_Config, only : outputdir
    
    ! --- in/out ------------------------------

    type(T_Emis_OPS), intent(inout)     ::  emo
    integer, intent(in)                 ::  year
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_Read'
    integer, parameter          ::  maxcol=20
    integer, parameter        ::  max_query = 5 ! maximum amount of defined profiles

    ! --- local ------------------------------

    type(TrcFile)                   ::  rcF
    character(len=512)              ::  list
    character(len=512)              ::  pm_composition_file
    character(len=512)              ::  ec_fraction_file
    character(len=512)              ::  base_file
    character(len=1024)             ::  query_list
    character(len=512)              ::  query
    character(len=512)              ::  queries(max_query)
    integer                         ::  iquery
    integer                         ::  ipol
    integer                         ::  iline
    
    character(len=1)                ::  sep, comment, charc
    
    integer                         ::  snr, x,y,r, dv, area,ps
    character(len=8)                ::  cat
    integer                         ::  cat_nr
    real                            ::  lon, lat, Q, hc, h, s
    integer                         ::  ilon, ilat
    integer                         ::  ilon1, ilon2, ilat1, ilat2
    real                            ::  west, east, south, north
    real, pointer                   ::  fracs(:,:)
    character(len=12)               ::  comp
    integer                         ::  isource
    logical                         ::  indomain
    integer                         ::  icat, icountry
    integer                         ::  ichange
    integer                         ::  iter

    logical                         ::  exist
    integer                         ::  fu
    
    character(len=1024)             ::  line
    integer                         ::  nheader
    character(len=512)              ::  headers(maxcol)
    character(len=64)               ::  header
    character(len=64)               ::  field
    
    real                            ::  cat_frac(11)
    integer                         ::  ips, icat_frac
    
    character(len=64)               ::  model, expid
    character(len=1024)             ::  basename
    character(len=1024)             ::  fname
    
    ! meteo data:
    real, pointer        ::  lons(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  lats(:,:,:)   ! (lon,lat,1)

    ! --- begin -------------------------------

    call LE_Data_GetPointer( 'lon', lons, status, check_units ='degrees_east' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lat', lats, status, check_units ='degrees_north' )
    IF_NOTOK_RETURN(status=1)

    ! store year:
    emo%year = year

    !
    ! settings
    !
    
    ! init file:
    call Init( rcF, trim(emo%rcfile), status )
    IF_NOTOK_RETURN(status=1)

    ! read pm-distribution file (this is also used for country/SNAP information)
    call ReadRc( rcF, trim(emo%rckey)//'.pm_composition', pm_composition_file, status )
    IF_NOTOK_RETURN(status=1)    
    call LE_Emis_OPS_PM_composition_Init( pm_composition_file, &
                                          max_country, max_cat, &
                                          emo%country_ISO3, emo%country_code, emo%ncountry, &
                                          emo%cat_code, emo%ncat, &
                                          emo%pm_comp, &
                                          status )
    IF_NOTOK_RETURN(status=1)   
        
    allocate ( emo%ec_frac( emo%ncountry, emo%ncat ) )
    emo%ec_frac = -1
    ! read ec fraction file
    call ReadRc( rcF, trim(emo%rckey)//'.ec_fraction', ec_fraction_file, status )
    IF_NOTOK_RETURN(status=1)
    call LE_Emis_OPS_EC_fractions_Init ( ec_fraction_file, &
                                         emo%country_ISO3(1:emo%ncountry), emo%ncountry, &
                                         emo%cat_code(1:emo%ncat), emo%ncat, &
                                         emo%ec_frac, &
                                         status )
    IF_NOTOK_RETURN(status=1)

    ! read in emission input for neceassary pollutants    
    call ReadRc( rcF, trim(emo%rckey)//'.pollutants', list, status )
    IF_NOTOK_RETURN(status=1)
    
    call goSplitString( list, emo%npol, emo%pol_names, status )
    IF_NOTOK_RETURN(status=1)
    
    !--------------------------------------------------------
    write(*,*) '  Read height distribution'
    !--------------------------------------------------------

    ! key to set or read distribution:
    call ReadRc( rcF, trim(emo%rckey)//'.height_distribution', query, status )
    IF_NOTOK_RETURN(status=1)
    ! if not empty than this is enabled ...
    emo%with_height_distribution = len_trim(query) > 0

    !! check ...
    !if ( emt%with_stack_height .and. emt%with_height_distribution ) then
    !  write (gol,'("could not use both stack-height table and height-distribution table")'); call goErr
    !  TRACEBACK; status=1; return
    !end if

    ! defined ?
    if ( emo%with_height_distribution ) then
      ! setup:
      call emo%hdistr%Init( query, emo%cat_code(1:emo%ncat), status )
      IF_NOTOK_RETURN(status=1)
    end if

    !--------------------------------------------------------
    write(*,*) '  Setup time profiles ...'
    !--------------------------------------------------------

    ! set dummy flag for safety:
    emo%timetype = 'None'

    ! read query:
    call ReadRc( rcF, trim(emo%rckey)//'.time_profiles', query_list, status )
    IF_NOTOK_RETURN(status=1)
    
    call GoSplitString( query_list, emo%ntimetype, queries, status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over different timetypes
    do iquery = 1, emo%ntimetype
      
      ! 
      query = trim(queries(iquery) )
      
      call goVarValue( trim(query), ';', 'type', '=', emo%timetype(iquery) , status )
      IF_ERROR_RETURN(status=1)
    
      ! switch:
      select case ( trim(emo%timetype(iquery)) )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'hourly' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! info ...
          write (*,*) '  read hourly time factors: Country/Snap dependent'

          ! replace some values:
          call goReplace( query, '%{year}', '(i4.4)', emo%year, status )
          IF_NOTOK_RETURN(status=1)

          ! read:
          call LE_Emis_Time_Prof_Init( emo%time_prof, query, emo%year, &
                                         emo%cat_code(1:emo%ncat), emo%country_ISO3(1:emo%ncountry), status )
          IF_NOTOK_RETURN(status=1)
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          write (gol,'("unsupported time profile type : ",a)') trim(emo%timetype(iquery)); call goErr
          TRACEBACK; status=1; return

      end select  ! timeprof type
    
    end do ! queries of defined time profiles


    ! --------------------------------------
    write( *,*) ' Read emissions'
    ! --------------------------------------
    
    ! total emissions per country/SNAP/Pollutant
    allocate(emo%country_snap_total(emo%npol, emo%ncountry,emo%ncat) )
    emo%country_snap_total = 0.0
    
    ! total annual gridded emissions
    allocate(emo%annual_total_gridded(emo%npol, ugg%nlon, ugg%nlat) )
    emo%annual_total_gridded = 0.0
    
    ! allocate emission properties
    allocate(emo%pol(emo%npol) )
    ! loop over pollutants    
    do ipol = 1, emo%npol
      
      ! Fill in emissions for this pollutant   
      emo%pol(ipol)%name = trim(emo%pol_names(ipol))
      
      ! Scenario factor per pollutant
      call ReadRc( rcF, trim(emo%rckey)//'.'//trim(emo%pol(ipol)%name)//'.scenfac', emo%pol(ipol)%scenfac, status )
      IF_NOTOK_RETURN(status=1)
            
      ! inout file per pollutant
      call ReadRc( rcF, trim(emo%rckey)//'.'//trim(emo%pol(ipol)%name), query, status )
      IF_NOTOK_RETURN(status=1)
      
      ! file
      base_file = ''
      call goVarValue( trim(query), ';', 'file', '=', base_file, status )
      IF_ERROR_RETURN(status=1)
      
      ! replace year string for actual year
      call goReplace( base_file, '%{year}', '(i4.4)', emo%year, status )
      IF_NOTOK_RETURN(status=1)

      ! comment lines
      comment = '#' 
      call goVarValue( trim(query), ';', 'comment', '=', comment, status )
      IF_ERROR_RETURN(status=1)
            
      ! seperation character
      sep = ';'
      call goVarValue( trim(query), ';', 'sep', '=', sep, status )
      IF_ERROR_RETURN(status=1)
      if ( sep == 't' ) then
        sep = achar(9) ! tab character
      end if
      
      ! file should be present
      inquire( file=trim(base_file), exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found : ",a)') trim(base_file); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! new file unit
      call goGetFU( fu, status )
      IF_NOTOK_RETURN(status=1)
      
      ! read twice, first to define number of sources      
      do iter = 1,2
        ! open file:
        open( fu, file=trim(base_file), status='old', form='formatted', iostat=status )
        if (status/=0) then
          write (gol,'("opening file : ",a)') trim(base_file); call goErr
          TRACEBACK; status=1; return
        end if

        ! no sources yet 
        isource = 0

        ! line counter
        iline = 0

        do 
          ! read line
          read ( fu, '(a)', iostat=status) line
          if ( status/= 0 ) then
            write (gol,'("reading header line from file : ",a)') trim(base_file); call goErr
            TRACEBACK; status=1; return
          end if
          ! empty ? then skip:
          if ( len_trim(line) == 0 ) cycle
          ! comment ? then skip:
          if ( line(1:1) == comment ) cycle
          ! found non-comment line, leave loop:
          exit
        end do

        ! read rest of the file
        do
          iline = iline + 1
          ! read line
          read( fu, '(a)', iostat=status ) line
          if ( status < 0 ) then
            ! end of file reached
            exit
          else if (status > 0 ) then
            write( gol, '("Error reading line ", i8, "in: ", a)') iline, trim(base_file) ; call goErr
            TRACEBACK;status=1;return
          end if
          if ( line(1:4) == ' snr' .or. line (1:4) == ' ssn' ) cycle
          ! check if it is lon/lat or x/y format
          if ( line(9:9) == '.' ) then
            ! lon lat format
            read(line, '(I4, 2F8.3, E10.3, F7.3, F6.1, I7, F6.1, 4I4, 2x, A12)', iostat=status) snr, lon,lat,Q,hc,h,r,s,dv,cat_nr,area,ps,comp
            if ( status /= 0 ) then
              write( gol, '("Error reading in file: ", a)') trim(base_file) ; call goErr
              write( gol, '("Error reading: ", a)') trim(line) ; call goErr
              TRACEBACK;status=1;return
            endif

          else 
            ! x/y rdm format
            read(line, '(I4, 2I8, E10.3, F7.3, F6.1, I7, F6.1, 4I4, 2x, A12)', iostat=status) snr,x,y,Q,hc,h,r,s,dv,cat_nr,area,ps,comp
            if ( status /= 0 ) then
              write( gol, '("Error reading in file: ", a)') trim(base_file) ; call goErr
              write( gol, '("Error reading: ", a)') trim(line) ; call goErr
              TRACEBACK;status=1;return
            endif

            ! convert rdm coordinates to lon/lat 
            call rdm2geo( x,y,lon,lat, status )
            IF_NOTOK_RETURN(status=1)

          end if

          ! check if this source is in model domain
          call ugg%InDomain( lon, lat, indomain, status )
          IF_NOTOK_RETURN(status=1)        
          if ( .not. indomain ) then
            !print *, 'OUTSIDE: ', lon,lat
            ! next source
            cycle
          end if

          ! match category code
          if ( cat_nr > 10 ) then
            ! ps code
            ! PM10
            ! check whether combination of country and PS-code should be changed (sea shipping NLD--> sea shipping NOS )
            do ichange = 1, n_change_PS_Code
              if ( change_PS_Code(1,1,ichange) == area .and. change_PS_Code(2,1,ichange) == cat_nr ) then
                area = change_PS_Code(1,2,ichange)
                cat_nr = change_PS_Code(2,2,ichange)
              end if            
            end do

            cat_frac = 0.0
            if ( trim(emo%pol(ipol)%name) == 'PM10' ) then
              call goMatchValue( cat_nr, ps_codes_PM10, ips, status) 
              IF_NOTOK_RETURN(status=1)
              cat_frac = Translate_PM10(:,ips)
            else if ( trim(emo%pol(ipol)%name) == 'NOx' ) then
              call goMatchValue( cat_nr, ps_codes_NOx, ips, status) 
              IF_NOTOK_RETURN(status=1)
              cat_frac = Translate_NOx(:,ips)
            else if ( trim(emo%pol(ipol)%name) == 'NH3' ) then
              call goMatchValue( cat_nr, ps_codes_NH3, ips, status) 
              IF_NOTOK_RETURN(status=1)
              cat_frac = Translate_NH3(:,ips)
            else if ( trim(emo%pol(ipol)%name) == 'SO2' ) then
              call goMatchValue( cat_nr, ps_codes_SO2, ips, status) 
              IF_NOTOK_RETURN(status=1)
              cat_frac = Translate_SO2(:,ips)
            else
              write( gol, '("Unknown tracer for snap distribution of PS-code")' ) ; call goErr
              write( gol, '("search for tracer: ", a)') trim(emo%pol(ipol)%name) ; call goErr
              TRACEBACK;status=1;return
            end if
          else
            ! init
            cat_frac = 0.0
            ! convert:
            write (cat,'(i0)') cat_nr
            ! category code  
            call goMatchValue( cat, emo%cat_code, icat, status )
            IF_NOTOK_RETURN(status=1)
            cat_frac(icat) = 1.0
          endif


          ! macth country code
          call goMatchValue( area, emo%country_code, icountry, status )
          IF_NOTOK_RETURN(status=1)

          ! For area source, distribute emission over several grid cells.
          ! for this execrcise area sources are taken on lon/lat square around the center
          if ( r > 0.0 ) then

            select case ( trim(emo%country_ISO3(icountry)) )

              case ( 'NLD', 'NOS' )
                ! RDM-coordinate grid
                ! longitude of west and east boundary defined as shift from center lon
                west = (lon*deg2rad - 0.5*r / (ae*cos(lat*deg2rad)) )/deg2rad
                east = (lon*deg2rad + 0.5*r / (ae*cos(lat*deg2rad)) )/deg2rad

                ! latitude of south and north boundary defined as shift from center lat
                south = (lat*deg2rad - 0.5*r/ae) / deg2rad
                north = (lat*deg2rad + 0.5*r/ae) / deg2rad

              case ( 'BEL', 'GBR', 'DEU', 'LUX', 'DNK', 'FRA' )
                ! original area sources on 1/8x1/16 lon/lat
                west = lon - 0.0625
                east = lon + 0.0625
                south = lat - 0.03125
                north = lat + 0.03125     

              case ( 'ALB', 'ARM', 'AUT', 'AZE', 'BGR', 'BIH', 'BLR', 'CHE', 'CYP', &
                     'CZE', 'ESP', 'EST', 'FIN', 'GEO', 'GRC', 'HRV', 'HUN', 'IRL', &
                     'ITA', 'LTU', 'LVA', 'MDA', 'MKD', 'MLT', 'NOR', 'POL', 'PRT', &
                     'ROU', 'RUS', 'SVK', 'SVN', 'SWE', 'TUR', 'UKR', 'YUG' )
                ! original area sources on 3/2x3/4 lon/lat
                west = lon - 0.75
                east = lon + 0.75
                south = lat - 0.375
                north = lat + 0.375     
              case default
                !
                write ( gol, '("Unknown emission resolution for country: ",a)') trim(emo%country_ISO3(icountry))
                TRACEBACK;status=1;return

            end select

            ! calculate contributions to each LE-gridcell
            call ugg%GetDistribution( west, east, south, north, &
                                  ilon1, ilon2, ilat1, ilat2, fracs, &
                                  status )

            IF_NOTOK_RETURN(status=1)

          else ! point source, only one target cell

            ! indices in LE-grid
            call ugg%GetLocation( lon, lat, ilon, ilat, status )
            IF_NOTOK_RETURN(status=1)

            ! only one target cell, with all emission in this cell          
            ilon1 = ilon
            ilon2 = ilon
            ilat1 = ilat
            ilat2 = ilat
            if ( associated(fracs)) nullify(fracs)
            allocate(fracs(ilon:ilon,ilat:ilat) )
            fracs(ilon,ilat) = 1.0

          end if

          ! loop over target cells to distribute emissions
          do ilon = ilon1, ilon2
          do ilat = ilat1, ilat2
            do icat_frac = 1, 11
              if ( cat_frac(icat_frac) /= 0.0 ) then
                icat = icat_frac
  !              if ( icat_frac == 11 ) then 
  !                ! use this further on as category 8
  !                icat = 8
  !              end if

                if (iter == 1 ) then
                  ! count:
                  isource = isource + 1
                else if (iter == 2 ) then
                  ! fill:

                  isource = isource + 1
                  emo%pol(ipol)%lon(isource)      = lons(ilon,ilat,1)  
                  emo%pol(ipol)%lat(isource)      = lats(ilon,ilat,1)  
                  emo%pol(ipol)%ilon(isource)     = ilon               
                  emo%pol(ipol)%ilat(isource)     = ilat               
                  emo%pol(ipol)%Q(isource)        = Q*fracs(ilon,ilat)*cat_frac(icat)*emo%pol(ipol)%scenfac
                  emo%pol(ipol)%hc(isource)       = hc                 
                  emo%pol(ipol)%h(isource)        = h                  
                  emo%pol(ipol)%r(isource)        = r                  
                  emo%pol(ipol)%s(isource)        = s                  
                  emo%pol(ipol)%dv(isource)       = dv                       
                  emo%pol(ipol)%icat(isource)     = icat               
                  emo%pol(ipol)%area(isource)     = area               
                  emo%pol(ipol)%icountry(isource) = icountry    
                  emo%pol(ipol)%country(isource)  = trim(emo%country_ISO3(icountry))       
                  emo%pol(ipol)%ps(isource)       = ps       

                  emo%country_snap_total(ipol,icountry,icat) = &
                      emo%country_snap_total(ipol,icountry,icat) + emo%pol(ipol)%Q(isource)
                  emo%annual_total_gridded(ipol,ilon,ilat) = &
                      emo%annual_total_gridded(ipol,ilon,ilat) + emo%pol(ipol)%Q(isource)

                endif

              end if ! positive fraction for this category
            end do ! fraction categories
          end do
          end do                           

        end do

        ! clean
        if ( associated(fracs)) nullify(fracs)
        
        if ( iter == 1 ) then
          ! number of sources for this pollutant
          emo%pol(ipol)%nsource = isource
          
          allocate( emo%pol(ipol)%lon( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%lat( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%ilon( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%ilat( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%Q( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%hc( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%h( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%s( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%r( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%dv( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%icat( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%area( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%icountry( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%country( emo%pol(ipol)%nsource ) )
          allocate( emo%pol(ipol)%ps( emo%pol(ipol)%nsource ) )
        end if
        
        ! close:
        close( unit=fu, iostat=status )
        IF_NOTOK_RETURN(status=1)
  
      end do ! two times read of file
         
      ! pollutant distribution over different model tracers
      select case ( trim(emo%pol(ipol)%name) )
        
        case ( 'NOx', 'NOX', 'nox' )
          emo%pol(ipol)%nspec = 2
          emo%pol(ipol)%ispecs(1:2)   = (/i_no, i_no2/)
          emo%pol(ipol)%compfrac(1:2) = (/0.97,0.03/)
          emo%pol(ipol)%xm            = xm_N + 2*xm_O
          
        case ( 'SOx', 'sox', 'SO2', 'so2' )
          emo%pol(ipol)%nspec = 1
          emo%pol(ipol)%ispecs(1)     = i_so2
          emo%pol(ipol)%compfrac(1)   = 1.0
          emo%pol(ipol)%xm            = xm_S + 2*xm_O

        case ( 'PM10' )

          ! read Emission compositions
          call ReadRc( rcF, trim(emo%rckey)//'.'//trim(emo%pol(ipol)%name)//'.composition', query, status )
          IF_NOTOK_RETURN(status=1)
      
          ! initialize given query:
          ! Get composition for fine part of PM10
          call emo%pol(ipol)%emcomp%Init( trim(query), 'pm25', &
                                         emo%cat_code, emo%country_ISO3, &
                                         specname(1:nspec), year, status )
          IF_NOTOK_RETURN(status=1)

          ! Distribution of PM10 is done by composition files
          emo%pol(ipol)%nspec = 0

        case ( 'NH3' )
          ! Distribution of PM10 is done by composition files
          emo%pol(ipol)%nspec = 1
          emo%pol(ipol)%ispecs(1)     = i_nh3
          emo%pol(ipol)%compfrac(1)   = 1.0
          emo%pol(ipol)%xm            = xm_N + 3*xm_H

        case default
        
          write( gol, '("Do not know how to distribute pollutant: ",a , " over model tracers")' ) trim(emo%pol(ipol)%name); call goErr
          TRACEBACK; status=1; return
      end select                  

      ! state units of input file (should actually be read from header)
      emo%pol(ipol)%unit = 'g/s'
    end do ! pollutants
    
    !-------------------------------------------
    write (*,*) ' summary tables'
    !-------------------------------------------

    ! output directory etc:
    call ReadRc( rcF, 'le.output.model', model, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'le.output.expid', expid, status )
    IF_NOTOK_RETURN(status=1)

    ! base name of files:  /output/LE_RUNID_emission_macc_2007
    write (basename,'(a,"/",a,"_",a,"_emis_",a,"_",i4.4)') &
               trim(outputdir), trim(model), trim(expid), &
               trim(emo%label), emo%year

    ! seperation character:
    sep = ';'
    
    ! output file:
    write (fname,'(a,"_country_SNAP.csv")') trim(basename)
    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open:
    open( unit=fu, file=trim(fname), form='formatted', iostat=status )
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(3a)') 'country', &
                            sep, 'sector'
    do ipol = 1, emo%npol
      field = trim(emo%pol(ipol)%name)
      line = trim(line)//sep//trim(field)
    end do
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(3a)') 'iso3', sep, 'index'
    do ipol = 1, emo%npol
      field = trim(emo%pol(ipol)%unit)
      line = trim(line)//sep//trim(field)
    end do
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! loop over countries
    do icountry = 1, emo%ncountry
      ! loop over categories
      do icat = 1, emo%ncat
        ! write record:
        write (line, '(a,a,i2.2)') &
                       trim(emo%country_ISO3(icountry)), &
                       sep, emo%cat_code(icat)
        !
        do ipol = 1, emo%npol
          write( field, '(e12.4)') emo%country_snap_total(ipol,icountry,icat)
          line = trim(line)//sep//trim(field)
        end do
        write (fu,'(a)',iostat=status) trim(line)
        IF_NOTOK_RETURN(status=1)
      end do  ! categories
    end do  ! countries
    ! close:
    close( unit=fu, iostat=status )
    IF_NOTOK_RETURN(status=1)

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
    do ipol = 1, emo%npol
      field = trim(emo%pol(ipol)%name)
      line = trim(line)//sep//trim(field)
    end do
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! write header:
    write (line,'(a)') 'iso3'
    do ipol = 1, emo%npol
      field = trim(emo%pol(ipol)%unit)
      line = trim(line)//sep//trim(field)
    end do
    write (fu,'(a)',iostat=status) trim(line)
    IF_NOTOK_RETURN(status=1)
    ! loop over countries
    do icountry = 1, emo%ncountry
      ! write record:
      write (line, '(a)') &
                     trim(emo%country_ISO3(icountry))
      !
      do ipol = 1, emo%npol
        write( field, '(e12.4)') sum(emo%country_snap_total(ipol,icountry,:))
        line = trim(line)//sep//trim(field)
      end do
      write (fu,'(a)',iostat=status) trim(line)
      IF_NOTOK_RETURN(status=1)
    end do  ! countries
    ! close:
    close( unit=fu, iostat=status )
    IF_NOTOK_RETURN(status=1)

    ! write out ncfile with annual total emissions
    call LE_Emis_OPS_Write_Annual_Total( emo%npol, emo%pol_names(1:emo%npol), emo%annual_total_gridded, emo%pol(1)%unit, basename, status )
    IF_NOTOK_RETURN(status=1)
    
!    ! fill in time profile patterns
!    emo%daily_pattern(0,0:23)  = (/1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00, &
!                                   1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/)
!    emo%daily_pattern(1,0:23)  = (/0.73,0.73,0.69,0.69,0.68,0.68,1.00,1.00,1.29,1.29,1.31,1.31, &
!                                   1.24,1.24,1.21,1.21,1.09,1.09,0.97,0.97,0.93,0.93,0.86,0.86/)
!    emo%daily_pattern(2,0:23)  = (/0.33,0.33,0.33,0.33,0.35,0.35,0.80,0.80,1.50,1.50,1.55,1.55, &
!                                   1.20,1.20,1.16,1.16,1.22,1.22,1.35,1.35,1.45,1.45,0.77,0.77/)
!    emo%daily_pattern(3,0:23)  = (/0.24,0.24,0.16,0.16,0.23,0.23,1.50,1.50,1.75,1.75,1.21,1.21, &
!                                   1.27,1.27,1.54,1.54,1.90,1.90,1.12,1.12,0.60,0.60,0.48,0.48/)
!
!    emo%daily_pattern(31,0:23) = (/0.24,0.24,0.16,0.16,0.23,0.23,1.50,1.50,1.75,1.75,1.21,1.21, &
!                                   1.27,1.27,1.54,1.54,1.90,1.90,1.12,1.12,0.60,0.60,0.48,0.48/)
!    emo%daily_pattern(32,0:23) = (/0.24,0.24,0.16,0.16,0.23,0.23,1.50,1.50,1.75,1.75,1.21,1.21, &
!                                   1.27,1.27,1.54,1.54,1.90,1.90,1.12,1.12,0.60,0.60,0.48,0.48/)
!    emo%daily_pattern(33,0:23) = (/0.24,0.24,0.16,0.16,0.23,0.23,1.50,1.50,1.75,1.75,1.21,1.21, &
!                                   1.27,1.27,1.54,1.54,1.90,1.90,1.12,1.12,0.60,0.60,0.48,0.48/)
    ! ok
    status = 0

  end subroutine LE_Emis_OPS_Read


  ! ***


  subroutine LE_Emis_OPS_Clear( emo, status )


    ! --- in/out ------------------------------

    type(T_Emis_OPS), intent(inout)    ::  emo
    integer, intent(out)               ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_Clear'

    ! --- local -------------------------------

    integer     ::  ipol

    ! --- begin -------------------------------
    
    emo%country_code = 0
    emo%cat_code     = ''

    ! clear emission info
    do ipol = 1, emo%npol    
      deallocate( emo%pol(ipol)%lon )
      deallocate( emo%pol(ipol)%lat )
      deallocate( emo%pol(ipol)%ilon )
      deallocate( emo%pol(ipol)%ilat )
      deallocate( emo%pol(ipol)%Q )
      deallocate( emo%pol(ipol)%hc )
      deallocate( emo%pol(ipol)%h )
      deallocate( emo%pol(ipol)%s )
      deallocate( emo%pol(ipol)%r )
      deallocate( emo%pol(ipol)%dv )
      deallocate( emo%pol(ipol)%icat )
      deallocate( emo%pol(ipol)%area )
      deallocate( emo%pol(ipol)%country )
      deallocate( emo%pol(ipol)%icountry )
      deallocate( emo%pol(ipol)%ps )
    end do
    deallocate( emo%pol )
    deallocate( emo%country_snap_total)
    deallocate( emo%ec_frac)
    deallocate( emo%annual_total_gridded)
    ! ok
    status = 0

  end subroutine LE_Emis_OPS_Clear


  ! ***


  subroutine LE_Emis_OPS_Setup( emo, emis_a, t1, t2, status )

    use GO       , only : TDate, Get, operator(+), operator(-), operator(/)
!    use GO       , only : GO_Timer_Start, GO_Timer_End
    use GO       , only : calc_DayNumber
    use Num      , only : IntervalSum
    use dims     , only : nx, ny, nz, nspec
    use LE_Data, only : LE_Data_GetPointer
    use LE_Time  , only : local_time2
    use indices  , only : specunit
    use indices  , only : i_ppm_1, i_ppm_1_25, i_ppm_25_4, i_ppm_4_10
    use indices  , only : i_ec_1, i_ec_1_25

    ! --- in/out ---------------------------

    type(T_Emis_OPS), intent(inout)       ::  emo
    real, intent(inout)                   ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_Setup'

    ! conversion parameters 
    real, parameter   ::  g_per_sec_to_g_per_min = 60.
    
    real, parameter   ::  g_per_sec_to_ug_per_min = 60.0*1e6
    
    ! --- local ----------------------------
    
    ! time definitions
    type(TDate)           ::  tmid
    integer               ::  yy, mm, dd, hh
    integer               ::  yyh, mmh, ddh, hhh, iday, ihour
    integer               ::  julian_day
    
    ! source properties
    integer               ::  ipol
    integer               ::  isource
    integer               ::  ilon, ilat, iz
    integer               ::  ispec
    integer               ::  icat, icountry

!    integer              ::  temp_index
!    real                 ::  temp_fac_CO, temp_fac_VOC
    
    ! height profile of emissions
    real, allocatable     ::  hd_profiles(:,:,:,:)  ! (nx,ny,nz,ncat)
    real, allocatable     ::  profile(:)  ! height profile of emission in LE-grid cells    
    real, allocatable     ::  hhb(:)
    integer               ::  ilast
    real                  ::  emtop
    
    ! timepattern
    integer               ::  itimetyp
    real                  ::  timefact
    
    ! compositionof each pollutant
    integer               ::  icompspec
    real                  ::  compfact
    integer               ::  icomp
    
    ! unit conversion
    character(len=128)    ::  conversion
    real                  ::  unit_conv
    
    ! emission
    real                  ::  delta_emis(nz)
    ! meteo data:
    real, pointer        ::  h_m(:,:,:)   ! (lon,lat,nz)

    ! --- begin ----------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 'h', h_m, status,check_units='m')
    IF_NOTOK_RETURN(status=1)

    ! extract time parameters:
    tmid = t1 + (t2-t1)/2
    call Get( tmid, year=yy, month=mm, day=dd, hour=hh )
    
    ! other year then loaded ?
    if ( tmid%year /= emo%year ) then
      ! cleanup previous data:
      call LE_Emis_OPS_Clear( emo, status )
      IF_NOTOK_RETURN(status=1)
      ! read for current year:
      call LE_Emis_OPS_Read( emo, tmid%year, status )
      IF_NOTOK_RETURN(status=1)
    end if
 
    ! injection profile:
    allocate( profile(1:nz) )

    ! height distribution ?
    if ( emo%with_height_distribution ) then
      ! lookup table:
      allocate( hd_profiles(nx,ny,nz,emo%ncat) )
      ! storage for model layer heighs:
      allocate( hhb(0:nz) )
      ! emission top:
      emtop = emo%hdistr%heightb(emo%hdistr%nlay)
      ! init :
      hd_profiles = 0.0
      ! loop over grid cells:
      do ilon = 1, nx
        do ilat = 1, ny
          ! extra model layer heights:
          hhb(0)    = 0.0  ! m
          hhb(1:nz) = h_m(ilon,ilat,1:nz)   ! m
          ! loop over categories:
          do icat = 1, emo%ncat
            ! loop over model layers:
            ilast = -1
            do iz = 1, nz
              ! fraction of profile in this interval:
              call IntervalSum( emo%hdistr%heightb, emo%hdistr%fraction(icat,:), &
                                 hhb(iz-1), min(hhb(iz),emtop), &
                                 hd_profiles(ilon,ilat,iz,icat), ilast, status )
              IF_NOTOK_RETURN(status=1)
              ! end ?
              if ( hhb(iz) > emtop ) exit
            end do  ! model layers
          end do  ! categories
        end do ! ix
      end do ! iy
      ! clear:
      deallocate( hhb )

    end if  ! with height distribution
    
    ! loop over pollutants 
    do ipol = 1, emo%npol
      
      ! loop over sources for this pollutant
      do isource = 1, emo%pol(ipol)%nsource 
        
        ! init
        delta_emis = 0.0
        
        ! extract
        icountry = emo%pol(ipol)%icountry(isource)
        icat     = emo%pol(ipol)%icat(isource)
        
        ! horizontal location
        ilon = emo%pol(ipol)%ilon(isource)
        ilat = emo%pol(ipol)%ilat(isource)
                
!        ! calculate local time of the source
!        call local_time2( trim(emo%pol(ipol)%country(isource)), &
!                          yy , mm , dd , hh , &   ! uct
!                          yyh, mmh, ddh, hhh, &   ! local time  [yyh-mmh-ddh,hhh:00:00,yyh-mmh-ddh,hhh+1:00:00)
!                          iday, &                  ! day-of-the-week number
!                          status )
!        IF_NOTOK_RETURN(status=1)

        ! loop over time types
        do itimetyp = 1, emo%ntimetype
          ! switch:
          select case ( trim(emo%timetype(itimetyp) ) )
              !~ profile for each hour, category, and country:
              case ( 'hourly' )
                ! set the time factor (equal for all components)
                julian_day = calc_DayNumber('gregorian',yy,mm,dd)
                ihour = (julian_day-1)*24+hh + 1
                ! Snap/country dependent time factor
                timefact = emo%time_prof%profile(ihour,icat,icountry)
                if ( timefact < 0 ) then
                  write( gol, '("invalid timefactor for category/country: ",i0," ",i0)' ) icat,icountry ; call goErr
                  TRACEBACK; status=1; return
                end if 
              case default
                write (gol,'("unsupported timetype `",a,"`")') trim(emo%timetype(itimetyp)); call goErr
                TRACEBACK; status=1; return
            end select
        end do        
        
        ! fill profile from height distribution table ?
        if ( emo%with_height_distribution ) then
          ! extract:
          profile = hd_profiles(ilon,ilat,:,icat)
        else 
          ! everythin in lowest layer
          profile = 0.0
          profile(1) = 1.0
        end if

        ! composition of PM10 tracers
        if ( emo%pol(ipol)%name == 'PM10') then

          ! unit conversion (to units of ppm_1, assume that all ppm and ec units are identical )
          conversion = trim(emo%pol(ipol)%unit)//' ; '//trim(specunit(i_ppm_1))
          select case ( trim(conversion) )       
            
            case ( 'g/s ; ppb' )
              unit_conv = g_per_sec_to_g_per_min/ ( emo%pol(ipol)%xm*1e3 )  ! mol/min
              
            case ( 'g/s ; ug/m3', 'g/s ; ug m-3' )
              unit_conv = g_per_sec_to_ug_per_min      ! ug/min
              
            case default
              write( gol, '(" do not know how to convert emissions: ", a )' ) trim(conversion) ; call goErr
              TRACEBACK;status=1;return

          end select

          ! Fine fraction divided over five subtracers
          do icomp = 1, emo%pol(ipol)%emcomp%ncomp
            ! target tracer:
            ispec = emo%pol(ipol)%emcomp%itracer(icomp)

            ! not a tracer in this run ? then skip:
            if ( ispec < 0 ) cycle
            
            ! extract:         # fraction fine subtracers             * Fraction fine in total PM
            compfact = emo%pol(ipol)%emcomp%frac(icat,icountry,icomp) * (1.0-emo%pm_comp(icountry,icat,3)-emo%pm_comp(icountry,icat,4) )
          
            !
            delta_emis = emo%pol(ipol)%Q(isource) * timefact * compfact * &
                            profile * unit_conv
          end do                            
!         
          ! Coarse Fractions of PM (pm25_4 and pm4_10) 
          ! ppm_2.5_4.0 (2.5-4.0 um)
          delta_emis = emo%pol(ipol)%Q(isource) * timefact * &
                          emo%pm_comp(icountry,icat,3) * & 
                            profile * unit_conv 
          emis_a(ilon,ilat,:,i_ppm_25_4) = emis_a(ilon,ilat,:,i_ppm_25_4) + delta_emis
        
          ! ppm_4.0-10.0 (4.0-10.0 um)
          delta_emis = emo%pol(ipol)%Q(isource) * timefact * &
                          emo%pm_comp(icountry,icat,4) * & 
                            profile * unit_conv 
          emis_a(ilon,ilat,:,i_ppm_4_10) = emis_a(ilon,ilat,:,i_ppm_4_10) + delta_emis
        ! Non PM-tracers
        else   
                                
          ! loop over model tracers composited per pollutant
          do icompspec = 1, emo%pol(ipol)%nspec

            ! find factors  
            ispec = emo%pol(ipol)%ispecs(icompspec)
            compfact = emo%pol(ipol)%compfrac(icompspec)

            ! unit conversion
            conversion = trim(emo%pol(ipol)%unit)//' ; '//trim(specunit(ispec))
            select case ( trim(conversion) )       

              case ( 'g/s ; ppb' )
                unit_conv = g_per_sec_to_g_per_min/ (emo%pol(ipol)%xm*1e3)   ! mol/min

              case ( 'g/s ; ug/m3', 'g/s ; ug m-3' )
                unit_conv = g_per_sec_to_ug_per_min      ! ug/min

              case default
                write( gol, '(" do not know how to convert emissions: ", a )' ) trim(conversion) ; call goErr
                TRACEBACK;status=1;return

              end select

            ! calculate emission
            delta_emis = emo%pol(ipol)%Q(isource) * timefact * compfact * profile * unit_conv
            ! add:
            emis_a(ilon,ilat,:,ispec) = emis_a(ilon,ilat,:,ispec) + delta_emis       
          end do ! composited tracers
        
        end if  
      end do ! source per pollutant
    end do ! pollutants         

    ! done
    deallocate( profile)
    
    ! ok
    status = 0

  end subroutine LE_Emis_OPS_Setup

! ---------------------------------

  subroutine rdm2geo( x,y ,lon,lat, status )
  
    integer, intent(in)     ::  x, y
    real, intent(out)       ::  lon, lat
    integer, intent(out)    ::  status
    
    ! constant
    character(len=*), parameter ::  rname = mname//'/rdm2geo'
    
    ! local
    real     ::  xx, yy
    
    ! begin
    xx = (x - 155000.)/100000.
    yy = (y - 463000.)/100000.

    lat = 3236.033*yy - 32.592*xx**2 - 0.247*yy**2 - 0.850*xx**2*yy - 0.065*yy**3 + 0.005*xx**4 - 0.017*xx**2*yy**2
    lon = 5261.305*xx + 105.979*xx*yy + 2.458*xx*yy**2 - 0.819*xx**3 + 0.056*xx*yy**3 - 0.056*xx**3*yy

    lon = lon/3600. + 5. + 23./60. + 15./3600. + 0.5/3600.
    lat = lat/3600. + 52. + 9./60. + 22./3600. + 0.178/3600.
    
    ! ok
    status = 0
  
  end subroutine rdm2geo



end module LE_Emis_OPS

