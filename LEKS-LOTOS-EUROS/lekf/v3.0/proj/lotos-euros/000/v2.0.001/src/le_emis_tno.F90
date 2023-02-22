!###############################################################################
!
! NAME
!   LE_Emis_TNO - routines to read and work with emission file in TNO format
!
! EMISSION FILES
!
! * Base emissions:
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
!   Format: NetCDF
!
! * Time profiles (distribution of emission in time, for each category)
!   Assigned per month in the year, day in the week, and hour in the day.
!
!   Supported rcfile settings:
!   !~ no profile, use unity factors:
!   le.emis.<set>.time_profiles  :  type=unity
!   !~ old month/day/hour profiles in text file:
!   le.emis.<set>.time_profiles  :  file=/path/time_var_emis.txt;type=mdh
!   !~ netcdf files with hourly factors for each year;
!   !  a key '%{year}' is replaced by the actual year:
!   le.emis.<set>.time_profiles  :  file=/path/time_profiles_%{year}.nc;type=hourly
!   !~ netcdf files with hourly factors for each year specified per component;
!   !  a key '%{year}' is replaced by the actual year:
!   le.emis.<set>.time_profiles  :  file=/path/time_profiles_%{year}.nc;type=hourly-comp
!   !  a key '%{year}' is replaced by the actual year:
!   le.emis.<set>.time_profiles  :  file=/path/time_profiles_gridded_%{year}.nc;type=hourly-gridded
!
! * VOC-profiles (distribution of VOC-emission over different species, for each category)
!   Assigned per category.
!   Format: text file
!
! * Scenario-profiles (used to reduce or increase emissions, for each country, for each category)
!   Assigned per category and country.
!   Format: text file
!
! * Emisison composition, e.g. distribution of NOx emisisons of NO and NO2, etc.
!   Assigned per category and country.
!   Format: text file
!
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

module LE_Emis_TNO

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile

#ifdef with_netcdf
  use NetCDF
#endif
  use LE_Emis_TNO_Base   , only : T_Emis_TNO_Base
  use LE_Emis_Composition, only : T_Emis_Composition
  use LE_Emis_Scenario   , only : T_Emis_Scenario
  use LE_Emis_HeightDistr, only : T_Emis_HeightDistr
  use LE_Emis_Time_Prof  , only : T_Emis_Time_Prof
  use LE_Emis_Time_Prof_Comp, only : T_Emis_Time_Prof_Comp
  use LE_Emis_Time_Prof_Grid, only : T_Emis_Time_Prof_Grid

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Emis_TNO

  public  ::  LE_Emis_TNO_Init, LE_Emis_TNO_Done
  public  ::  LE_Emis_TNO_Setup


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_TNO'

  ! maximum number of emitted tracers:
  integer, parameter            ::  max_emis = 20
  ! maximum number of defined time types
  integer, parameter                ::  max_time_type = 5

  ! for temporary files that are closed directly after use
  integer, parameter :: u_tmp  = 61


  ! --- local --------------------------------

  ! emission data base
  type T_Emis_TNO

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
    character(len=3), allocatable     ::  country_code(:)  ! country code (e.g. 'ALB', 'BEL')
    character(len=100), allocatable   ::  country_name(:)  ! country name

    ! emission distribution:
    character(len=32)                 ::  distr_nam   ! 'SNAP15'

    ! emission categories
    integer                           ::  ncat         ! number of emission categories
    integer, allocatable              ::  cat_code(:)  ! code-number of each emission category
    character(len=100), allocatable   ::  cat_nam(:)   ! name of each emission category

    ! number of emitted components:
    integer                           ::  nemis
    ! names:
    character(len=32)                 ::  emis_name(max_emis)
    ! molecule masses:
    real                              ::  xm(max_emis)  ! kg/mole

    ! base emissions:
    type(T_Emis_TNO_Base)             ::  emb

    ! scenario factors for each category, country, emitted species
    type(T_Emis_Scenario)             ::  scenario

    ! emisison compositions:
    type(T_Emis_Composition)          ::  emcomp(max_emis)

    ! skip flags per model tracer:
    logical, allocatable              ::  skip_tracer(:)   ! (nspec)

    ! VOC profiles
    integer                           ::  nvoc               ! number of VOC-species that are emitted
    real, pointer                     ::  vocprof(:,:)       ! vocprof(1:ncat,1:nvoc)

    !! stack height
    !logical                           ::  with_stack_height
    !real, allocatable                 ::  hstack(:) ! hstack(1:ncat)
    ! height distribution:
    logical                           ::  with_height_distribution
    type(T_Emis_HeightDistr)          ::  hdistr

    ! timing:
    integer                           ::  itim_height_distr

    ! description of time factors to be used:
    !   'unity' | 'mdh' | 'time_prof' | 'hourly-gridded'
    integer                           ::  ntimetype
    character(len=16)                 ::  timetype(max_time_type)
    ! time dependency factors
    real, allocatable                 ::  imonthdp(:,:)
    real, allocatable                 ::  idaydp(:,:)
    real, allocatable                 ::  ihourdp(:,:)
    ! time factors country/snap dependent
    type(T_Emis_Time_Prof)            ::  time_prof

    ! time factors country/snap/component dependent
    type(T_Emis_Time_Prof_Comp)       ::  time_prof_comp

    ! time factors grid dependent
    type(T_Emis_Time_Prof_Grid)       ::  time_prof_grid
    
    ! temperature variations for VOC and CO
    real, pointer                     ::  temp_var_VOC(:,:)
    real, pointer                     ::  temp_var_CO(:,:)

  end type T_Emis_TNO



contains


  ! ===============================================================


  subroutine LE_Emis_TNO_Init( emt, label, rcF, rckey, t, status )

    use GO, only : TrcFile
    use GO, only : GO_Timer_Def
    use GO, only : TDate
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    type(T_Emis_TNO), intent(out)       ::  emt
    character(len=*), intent(in)        ::  label
    type(TrcFile), intent(in)           ::  rcF
    character(len=*), intent(in)        ::  rckey
    type(TDate), intent(in)             ::  t       ! start time
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Init'

    ! --- local ------------------------------

    ! --- begin -------------------------------

    ! no data Init yet ...
    emt%year = -1

    ! store label:
    emt%label = trim(label)

    ! store info on settings:
    emt%rcfile = trim(rcF%fname)
    emt%rckey = trim(rckey)

    ! define timer:
    call GO_Timer_Def( emt%itim_height_distr, 'emis_height_distr', status )
    IF_NOTOK_RETURN(status=1)

    ! read for current year:
    call LE_Emis_TNO_Read( emt, t%year, status )
    IF_NOTOK_RETURN(status=1)

    ! enable data:
    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Init


  ! ***


  subroutine LE_Emis_TNO_Done( emt, status )

    ! --- in/out ------------------------------

    type(T_Emis_TNO), intent(inout)    ::  emt
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! any data loaded ?
    if ( emt%year > 0 ) then
      ! cleanup:
      call LE_Emis_TNO_Clear( emt, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Done


  ! ===============================================================


  subroutine LE_Emis_TNO_Read( emt, year, status )

    use Binas          , only : xm_H, xm_C, xm_N, xm_O, xm_S
    use GO             , only : ReadRc, Init, Done, ReadRc
    use GO             , only : goLoCase, goMatchValue, goReadFromLine
    use GO             , only : goReplace, goVarValue, goSplitString
    use GO             , only : goc
    use LE_Logging     , only : u_log, u_err
    use indices        , only : nspec, specname
    use dims           , only : nx, ny
    use dims           , only : runF
    use LE_IO_Tools    , only : debopt
    use LE_IO_Tools    , only : io_read_table_country_cat_data
    use LE_IO_Tools    , only : io_read_table1
    use LE_Emis_TNO_Base, only : LE_Emis_TNO_Base_Init
    use LE_Emis_TNO_Base, only : LE_Emis_TNO_Base_Summary
    use LE_Emis_Composition, only : LE_Emis_Composition_Init
    use LE_Emis_Scenario   , only : LE_Emis_Scenario_Init
    use LE_EMis_HeightDistr, only : LE_Emis_HeightDistr_Init
    use LE_Emis_Time_Prof  , only : LE_Emis_Time_Prof_Init
    use LE_Emis_Time_Prof_Comp, only : LE_Emis_Time_Prof_Comp_Init
    use LE_Emis_Time_Prof_Grid, only : LE_Emis_Time_Prof_Grid_Init
#ifdef with_labeling
    use SA_Labeling, only         : SA_Emis_Sectors
#endif

    ! --- in/out ------------------------------

    type(T_Emis_TNO), intent(inout)     ::  emt
    integer, intent(in)                 ::  year
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Read'
    integer, parameter        ::  max_query = 5 ! maximum amount of defined profiles

    ! --- local ------------------------------

    type(TrcFile)                   ::  rcF
    character(len=512)              ::  base_file
    character(len=512)              ::  aux_dir

    integer                         ::  icountry
    integer                         ::  icat
    integer                         ::  iemis
    integer                         ::  ispec

    character(len=250)              ::  fnam
    logical                         ::  exist
    character(len=500)              ::  colheader
    character(len=50), allocatable  ::  emis_tempprof(:)
    real, allocatable               ::  temp_var(:,:)
    character(len=50)               ::  fmt
    integer                         ::  ii
    integer                         ::  i0
    integer                         ::  irec
    logical                         ::  firstread
    integer                         ::  ncol
    integer                         ::  icol
    integer                         ::  t0, t1, dt, temp
    real                            ::  aver

    ! emission composition:
    character(len=512)              ::  query
    character(len=1024)             ::  query_list
    integer                         ::  iquery, nquery
    character(len=512)              ::  queries(max_query)
    ! skip list:
    character(len=32)               ::  emspec

    character(len=512)              ::  outdir
    character(len=512)              ::  basename
    character(len=32)               ::  model, expid

    logical                         ::  trac_dep  ! timeprofiles component dependent?

    ! --- begin -------------------------------

    ! store year:
    emt%year = year

    ! check implementation:
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read `TNO` emissions;")'); call goErr
    write (gol,'("define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    !
    ! settings
    !

    ! init file:
    call Init( rcF, trim(emt%rcfile), status )
    IF_NOTOK_RETURN(status=1)

    ! read emission distribution:
    call ReadRc( rcF, trim(emt%rckey)//'.distribution', emt%distr_nam , status )
    IF_NOTOK_RETURN(status=1)

    ! read name of data file:
    call ReadRc( rcF, trim(emt%rckey)//'.base.file', base_file, status )
    IF_NOTOK_RETURN(status=1)

    ! replace some values:
    call goReplace( base_file, '%{year}', '(i4.4)', emt%year, status )
    IF_NOTOK_RETURN(status=1)

    ! directory with auxilary files:
    call ReadRc( rcF, trim(emt%rckey)//'.aux.dir', aux_dir, status )
    IF_NOTOK_RETURN(status=1)

    !-------------------------------------------
    write (gol,'("  Base emissions")'); call goPr
    !-------------------------------------------

    ! read:
    call LE_EMIS_TNO_Base_Init( emt%emb, trim(base_file), status )
    IF_NOTOK_RETURN(status=1)

    !-------------------------------------------
    write (gol,'("  summary tables")'); call goPr
    !-------------------------------------------
    
    ! root only ...
    if ( goc%root ) then

      ! output directory etc:
      call ReadRc( rcF, 'le.output.outdir', outdir, status )
      IF_NOTOK_RETURN(status=1)
      call ReadRc( rcF, 'le.output.model', model, status )
      IF_NOTOK_RETURN(status=1)
      call ReadRc( rcF, 'le.output.expid', expid, status )
      IF_NOTOK_RETURN(status=1)

      ! base name of files:  /output/LE_RUNID_emission_macc_2007
      write (basename,'(a,"/",a,"_",a,"_emis_",a,"_",i4.4)') &
                 trim(outdir), trim(model), trim(expid), &
                 trim(emt%label), emt%year

      ! write:
      call LE_Emis_TNO_Base_Summary( emt%emb, basename, status )
      IF_NOTOK_RETURN(status=1)
      
    end if

    !------------------------------------------------------
    write (gol,'("  Read country codes and make indices for each country")'); call goPr
    !------------------------------------------------------

    ! Country file consists of iso3 country_code ('ALB')
    ! followed by long country name ('ALBANIA')
    ! note: the long name is not used here:

    ! number of countries:
    emt%ncountry = emt%emb%ncountry
    ! storage:
    allocate( emt%country_code(emt%ncountry) )
    allocate( emt%country_name(emt%ncountry) )
    ! loop over contries:
    do icountry = 1, emt%ncountry
      ! copy:
      emt%country_code(icountry) = trim(emt%emb%cnti(icountry)%code)
      emt%country_name(icountry) = trim(emt%emb%cnti(icountry)%name)
    end do

    !-------------------------------------------
    write (gol,'("  Read definition of emission categories")'); call goPr
    !-------------------------------------------

    ! Emission categories file consists of integer codes and name

    ! number of categories:
    emt%ncat = emt%emb%ncat
    ! storage:
    allocate( emt%cat_code(emt%ncat) )
    allocate( emt%cat_nam (emt%ncat) )
    ! loop over categories:
    do icat = 1, emt%ncat
      ! copy:
      emt%cat_code(icat) = emt%emb%cati(icat)%code1
      emt%cat_nam (icat) = trim(emt%emb%cati(icat)%name)
    end do

#ifdef with_labeling
    call SA_Emis_Sectors(emt%cat_nam, emt%ncat, status, tag=emt%label)
    IF_NOTOK_RETURN(status=1)
#endif

    !--------------------------------------------------
    write (gol,'("  Read base emissions")'); call goPr
    !--------------------------------------------------

    ! number of emitted tracers:
    emt%nemis = emt%emb%nemis
    ! check ...
    if ( emt%nemis > max_emis ) then
      write (gol,'("number of emission records (",i2,") exceeds maximum (",i2,");")') emt%nemis, max_emis; call goErr
      write (gol,'("increase parameter `max_emis` in module `",a,"`")') trim(mname); call goErr
      TRACEBACK; status=1; return
    end if
    ! loop over emitted tracers:
    do iemis = 1, emt%nemis
      ! copy:
      emt%emis_name(iemis) = trim(emt%emb%emi(iemis)%name)
    end do

    !-------------------------------------------
    write (gol,'("  Set molecule masses")'); call goPr
    !-------------------------------------------

    ! loop over emisisons:
    do iemis = 1, emt%nemis
      ! switch:
      select case ( trim(emt%emis_name(iemis)) )
        case ( 'so2', 'sox' ) ; emt%xm(iemis) = xm_S + xm_O * 2  ! kg/mole
        case ( 'nox' ) ; emt%xm(iemis) = xm_N + xm_O * 2  ! kg/mole  NO2 equivalent
        case ( 'nh3' ) ; emt%xm(iemis) = xm_N + xm_H * 3  ! kg/mole
        case ( 'co'  ) ; emt%xm(iemis) = xm_C + xm_O      ! kg/mole
        case ( 'co2' ) ; emt%xm(iemis) = xm_C + xm_O * 2  ! kg/mole
        case ( 'ch4' ) ; emt%xm(iemis) = xm_C + xm_H * 4  ! kg/mole
        case ( 'so4ns', 'so4ks', 'so4as', 'so4cs' )
                         emt%xm(iemis) = xm_S + xm_O * 4  ! kg/mole
        case default   ; emt%xm(iemis) = -999.9   ! to be trapped ...
      end select
    end do

    !-------------------------------------------
    write (gol,'("  Read emission compositions")'); call goPr
    !-------------------------------------------

    ! loop over emisisons:
    do iemis = 1, emt%nemis
      ! component fractions query:
      call ReadRc( rcF, trim(emt%rckey)//'.composition.'//trim(emt%emis_name(iemis)), query, status )
      IF_NOTOK_RETURN(status=1)
      ! replace some values:
      call goReplace( query, '%{year}', '(i4.4)', emt%year, status )
      IF_NOTOK_RETURN(status=1)
      ! initialize given query:
      call LE_Emis_Composition_Init( emt%emcomp(iemis), trim(query), &
                                     trim(emt%emis_name(iemis)), &
                                     emt%cat_code, emt%country_code, &
                                     specname(1:nspec), status )
      IF_NOTOK_RETURN(status=1)
    end do

    !--------------------------------------------------------
    write (gol,'("   Read skip list")'); call goPr
    !--------------------------------------------------------

    ! storage per model tracer:
    allocate( emt%skip_tracer(nspec) )
    ! by default skip nothing:
    emt%skip_tracer = .false.

    ! list with model tracer names to be skipped:
    call ReadRc( rcF, trim(emt%rckey)//'.skip.species', query, status )
    IF_NOTOK_RETURN(status=1)
    ! loop over list items:
    do
      ! leave if empty:
      if ( len_trim(query) == 0 ) exit
      ! extract:
      call goReadFromLine( query, emspec, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      ! index:
      call goMatchValue( trim(emspec), specname(1:nspec), ispec, status, quiet=.true. )
      IF_ERROR_RETURN(status=1)
      ! found ?
      if ( ispec > 0 ) then
        ! reset flag:
        emt%skip_tracer(ispec) = .true.
        ! info ...
        write (gol,'("    skip tracer ",a," ...")') trim(specname(ispec)); call goPr
      end if
    end do


    !--------------------------------------------------
    write (gol,'("  Read scenario factors")'); call goPr
    !--------------------------------------------------

    ! scenario factors file:
    call ReadRc( rcF, trim(emt%rckey)//'.scenario', query, status )
    IF_NOTOK_RETURN(status=1)
    ! replace some values:
    call goReplace( query, '%{year}', '(i4.4)', emt%year, status )
    IF_NOTOK_RETURN(status=1)

    ! read table:
    call LE_Emis_Scenario_Init( emt%scenario, query, &
                                  emt%cat_code, emt%country_code, &
                                  emt%emis_name(1:emt%nemis), &
                                  status )
    IF_NOTOK_RETURN(status=1)


    !--------------------------------------------------------
    write (gol,'("  Read VOC profiles")'); call goPr
    !--------------------------------------------------------

    ! for safety ...
    nullify( emt%vocprof )

    ! fixed file name ...
    write (fnam,'(a,"/voc_profiles_CBM4_",a,".txt")') trim(aux_dir), trim(emt%distr_nam)
    ! check ...
    inquire( file=trim(fnam), exist=exist )
    if ( .not. exist ) then
      write (gol,'("    WARNING - VOC profile file not found :")'); call goPr
      write (gol,'("    WARNING -  ",a)') trim(fnam); call goPr
      write (gol,'("    WARNING - continue, assume this is an emission file without VOC ...")'); call goPr
    else
      ! number of components:
      emt%nvoc = 10
      ! allocate array:
      allocate( emt%vocprof(emt%ncat,emt%nvoc) )

      ! CBM4: 10 emitted species  OLE PAR TOL XYL FORM ALD KET ACET ETH UNR
      ! read table with VOC splits;
      ! row    -> category
      ! column -> 10 VOC species
      firstread = .true.
      emt%nvoc      = 10
      colheader = 'code category_name OLE PAR TOL XYL FORM ALD KET ACET ETH UNR'
      call io_read_table1( 'VOC profiles CBM4', '[mol/kg VOC]', emt%ncat, emt%nvoc, &
                          colheader, emt%cat_code, emt%cat_nam, &
                          firstread, fnam, irec, u_tmp, u_log, u_err, debopt, &
                          emt%vocprof )
      ! close file
      close(u_tmp)
    end if

    !!--------------------------------------------------------
    !write (gol,'("   Read stack heights")'); call goPr
    !!--------------------------------------------------------

    ! In future: use effective heights form emission database if possible ..

    !! storage:
    !allocate( emt%hstack(emt%ncat) )

    !! default values: inject at surface:
    !emt%hstack = 0.0

    !! name of stack height file:
    !call ReadRc( rcF, trim(emt%rckey)//'.stack_height.file', query, status )
    !IF_NOTOK_RETURN(status=1)
    !! if not empty than this is enabled ...
    !emt%with_stack_height = len_trim(query) > 0

    !! read if necessary:
    !if ( emt%with_stack_height ) then
    !  ! read table with stack heights;
    !  ! row    -> category
    !  ! column -> stack height (1 column)
    !  firstread = .true.
    !  call io_read_table1('stack heights for point sources','[m]',emt%ncat,1, &
    !                      'code category_name stack_height', emt%cat_code, emt%cat_nam, &
    !                      firstread, trim(query), irec, &
    !                      u_tmp,u_log,u_err,debopt, &
    !                      emt%hstack )
    !  ! Close file
    !  close(u_tmp)
    !end if

    !--------------------------------------------------------
    write (gol,'("  Read height distribution")'); call goPr
    !--------------------------------------------------------

    ! key to set or read distribution:
    call ReadRc( rcF, trim(emt%rckey)//'.height_distribution', query, status )
    IF_NOTOK_RETURN(status=1)
    ! if not empty than this is enabled ...
    emt%with_height_distribution = len_trim(query) > 0

    !! check ...
    !if ( emt%with_stack_height .and. emt%with_height_distribution ) then
    !  write (gol,'("could not use both stack-height table and height-distribution table")'); call goErr
    !  TRACEBACK; status=1; return
    !end if

    ! defined ?
    if ( emt%with_height_distribution ) then
      ! setup:
      call LE_Emis_HeightDistr_Init( emt%hdistr, query, emt%cat_code, status )
      IF_NOTOK_RETURN(status=1)
    end if

    !--------------------------------------------------------
    write (gol,'("   Setup time profiles ...")'); call goPr
    !--------------------------------------------------------

    ! set dummy flag for safety:
    emt%timetype = 'None'

    ! read query:
    call ReadRc( rcF, trim(emt%rckey)//'.time_profiles', query_list, status )
    IF_NOTOK_RETURN(status=1)
    
    call GoSplitString( query_list, emt%ntimetype, queries, status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over different timetypes
    do iquery = 1, emt%ntimetype
      ! extract
      query = queries(iquery)
      
      ! read time factors file type;
      call goVarValue( trim(query), ';', 'type', '=', emt%timetype(iquery), status )
      IF_NOTOK_RETURN(status=1)
      ! switch:
      select case ( trim(emt%timetype(iquery)) )

      ! special value or file specification ?
        case ('unity' )

          ! info ...
          write (gol,'("     no time profiles, use unity factor ...")'); call goPr

        case ( 'mdh' )

          ! info ...
          write (gol,'("     read month/day/hour profiles ...")'); call goPr

          ! Allocate arrays:
          allocate( emt%imonthdp(emt%ncat,  12) )
          allocate( emt%idaydp  (emt%ncat,   7) )
          allocate( emt%ihourdp (emt%ncat,1:24) )  ! [00:00 - 01:00) is first record

          !  name of table file, raise error message if not defined (status<0):
          call goVarValue( trim(query), ';', 'file', '=', fnam, status )
          if ( status < 0 ) then
            write (gol,'("no keyword `file` in time profile query : ",a)') trim(query); call goErr
            TRACEBACK; status=1; return
          end if
          IF_NOTOK_RETURN(status=1)

          ! read table with 12 monthly factors
          ! row    -> category
          ! column -> month
          firstread = .true.
          ncol      = 12
          colheader = 'code category_name jan feb mar apr may jun jul aug sep oct nov dec'
          call io_read_table1('monthly time profiles emissions','[-]', emt%ncat, ncol, &
                              colheader, emt%cat_code, emt%cat_nam, &
                              firstread, fnam, irec, u_tmp, u_log, u_err, debopt, &
                              emt%imonthdp )
          ! check ...
          do icat = 1, emt%ncat
            aver = sum(emt%imonthdp(icat,:)) / real(ncol)
            if ( abs(aver-1.0) > 1.0e-4 ) then
              write (gol,'("category ",i2.2," month profile has not average 1.0 but ",f12.7)') icat, aver; call goErr
              TRACEBACK; status=1; return
            end if
          end do
          ! reset flag:
          firstread = .false.

          ! read table with 7 daily factors
          ! row    -> category
          ! column -> day in week
          ncol      = 7
          colheader = 'code category_name mon tue wed thu fri sat sun'
          call io_read_table1('daily time profiles emissions','[-]',emt%ncat,ncol, &
                              colheader,emt%cat_code,emt%cat_nam, &
                              firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
                              emt%idaydp )
          ! check ...
          do icat = 1, emt%ncat
            aver = sum(emt%idaydp(icat,:)) / real(ncol)
            if ( abs(aver-1.0) > 1.0e-4 ) then
              write (gol,'("category ",i2.2," day profile has not average 1.0 but ",f12.7)') icat, aver; call goErr
              TRACEBACK; status=1; return
            end if
          end do


          ! read table with 24 hourly factors
          ! row    -> category
          ! column -> hour
          ncol      = 24
          colheader = 'code category_name 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24'
          call io_read_table1('hourly time profiles emissions','[-]',emt%ncat,ncol, &
                              colheader, emt%cat_code, emt%cat_nam, &
                              firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
                              emt%ihourdp(:,1:ncol) )
          ! skip this, first record in file is valid for [00:00-01:00)
          ! set hour 0 equal to hour 24
          !emt%ihourdp(:,0) = emt%ihourdp(:,24)
          ! check ...
          do icat = 1, emt%ncat
            aver = sum(emt%ihourdp(icat,1:ncol)) / real(ncol)
            if ( abs(aver-1.0) > 1.0e-4 ) then
              write (gol,'("category ",i2.2," (",a,") hour profile has not average 1.0 but ",f12.7)') &
                              icat, trim(emt%cat_nam(icat)), aver; call goErr
              write (gol,'("values:")'); call goErr
              write (gol,*) emt%ihourdp(icat,1:ncol); call goErr
              TRACEBACK; status=1; return
            end if
          end do

          if (debopt > 1) then
             write(*,*) 'em/monthdp 7 ',emt%imonthdp(7,:)
             write(*,*) 'em/idaydp 7 ',emt%idaydp(7,:)
             write(*,*) 'em/ihourdp 7 ',emt%ihourdp(7,:)
          endif

          ! Close file
          close(u_tmp)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'hourly' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! info ...
          write (gol,'("  read hourly time factors: Country/Snap dependent")'); call goPr

          ! replace some values:
          call goReplace( query, '%{year}', '(i4.4)', emt%year, status )
          IF_NOTOK_RETURN(status=1)

          ! read:
          call LE_Emis_Time_Prof_Init( emt%time_prof, query, emt%year, &
                                         emt%cat_code, emt%country_code, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'hourly-comp' )
        !--------------------------------------------------------

          ! info ...
          write (gol,'("   Read time factors: Country/Snap/Component dependent")'); call goPr

          ! replace some values:
          call goReplace( query, '%{year}', '(i4.4)', emt%year, status )
          IF_NOTOK_RETURN(status=1)

          ! read:
          call LE_Emis_Time_Prof_Comp_Init( emt%time_prof_comp, query, emt%year, &
                                            emt%cat_code, emt%country_code, emt%emis_name(1:emt%nemis), status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'hourly-gridded' )
        !--------------------------------------------------------

          ! info ...
          write (gol,'("   Read time factors: Grid dependent")'); call goPr

          ! replace some values:
          call goReplace( query, '%{year}', '(i4.4)', emt%year, status )
          IF_NOTOK_RETURN(status=1)

          ! read:
          call LE_Emis_Time_Prof_Grid_Init( emt%time_prof_grid, query, &
                                            emt%cat_code, emt%year, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          write (gol,'("unsupported time profile type : ",a)') trim(emt%timetype(iquery)); call goErr
          TRACEBACK; status=1; return

      end select  ! timeprof type
    
    end do ! queries of defined time profiles
    
    !--------------------------------------------------------
    write (gol,'("   Read temperature dependent emission factors")'); call goPr
    !--------------------------------------------------------

    ! temperature dependent factors to be applied to VOC and CO emissions
    !
    ! extension to more species not trivial ; check code !!
    ! maybe extension to array temp_var(ncat,t0:t1,nemis,ndistr)

    ! by default clean:
    nullify( emt%temp_var_VOC )
    nullify( emt%temp_var_CO  )

    ! input file:
    write (fnam,'(a,"/temperature_var_emis_",a,".txt")') trim(aux_dir), trim(emt%distr_nam)
    ! check ...
    inquire( file=trim(fnam), exist=exist )
    if ( .not. exist ) then
      write (gol,'("    WARNING - VOC profile file not found :")'); call goPr
      write (gol,'("    WARNING -  ",a)') trim(fnam); call goPr
      write (gol,'("    WARNING - continue, assume this is an emission file without VOC ...")'); call goPr
    else

      ! Set species to read temp. dep. emission factors for:
      allocate(emis_tempprof(2))
      emis_tempprof = (/'VOC','CO '/)

      firstread = .true.

      ! Set temperature range and step size (integers !)
      t0   = -30
      t1   = 40
      dt   = 2
      ncol = (t1 - t0)/dt + 1

      ! Allocate arrays:
      allocate( emt%temp_var_VOC(emt%ncat,t0:t1) )
      allocate( emt%temp_var_CO (emt%ncat,t0:t1) )

      ! temporary storage:
      allocate( temp_var(emt%ncat,ncol) )

      ! Loop over species:
      do ii = 1,size(emis_tempprof)

         ! read table with temperature factors for temperatures ranging from t0=-30 to t1=40 degrees C
         ! with step dt = 2 degrees C
         ! row    -> category
         ! column -> temperature

         ! make column header with temperatures
         colheader = 'code category_name '
         i0 = len_trim(colheader) + 2
         do temp = t0,t1,dt
            write(colheader(i0:i0+4),'(i4)') temp
            i0 = i0 + 4
         enddo

         select case ( trim(emis_tempprof(ii)) )

          case ( 'VOC')

            ! Initialise:
            emt%temp_var_VOC = -999.

            ! Read data:
            call io_read_table1('temperature dependent emission factors'//trim(emis_tempprof(ii)), &
                                '[-]',emt%ncat,ncol, &
                                colheader,emt%cat_code,emt%cat_nam, &
                                firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
                                temp_var)
            firstread = .false.

            ! load data into array
            icol = 1
            do temp = t0,t1,dt
               emt%temp_var_VOC(:,temp) = temp_var(:,icol)
               icol = icol + 1
            enddo

         case ( 'CO' )

            ! Initialise:
            emt%temp_var_CO = -999.

            ! Read data:
            call io_read_table1('temperature dependent emission factors'//trim(emis_tempprof(ii)), &
                                '[-]',emt%ncat,ncol, &
                                colheader,emt%cat_code,emt%cat_nam, &
                                firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
                                temp_var)

            ! load data into array
            icol = 1
            do temp = t0,t1,dt
               emt%temp_var_CO(:,temp) = temp_var(:,icol)
               icol = icol + 1
            enddo

         case default

             write(u_err,*) ' unknown species for temperature dependent emission factors'
             stop

         end select

      enddo

      ! close:
      close(u_tmp)

      ! clear:
      deallocate( temp_var )

      ! interpolate for the odd values: only even values are read in from file
      !do i=-29,39,2
      do icat = 1, emt%ncat
         do temp = t0+1,t1-1,dt
           if (emt%temp_var_VOC(icat,temp-1) > 0.0 .AND. emt%temp_var_VOC(icat,temp+1) > 0) &
               emt%temp_var_VOC(icat,temp) = 0.5*(emt%temp_var_VOC(icat,temp-1) + emt%temp_var_VOC(icat,temp+1))
           if (emt%temp_var_CO(icat,temp-1) > 0.0 .AND. emt%temp_var_CO(icat,temp+1) > 0) &
               emt%temp_var_CO(icat,temp) = 0.5*(emt%temp_var_CO(icat,temp-1) + emt%temp_var_CO(icat,temp+1))
         enddo
      enddo

      ! write to log file
      if (debopt .gt. 0) then
         write(u_log,*) '-----------------------------------------------------------------------'
         write(u_log,*) 'temperature dependent emission factors VOC for ',trim(emt%distr_nam)
         fmt = '(i6,1x,a50,   (e12.5))'
         write(fmt(12:14),'(i3)') ncol
         do icat = 1,emt%ncat
            write(u_log,fmt) emt%cat_code(icat),trim(emt%cat_nam(icat)),emt%temp_var_VOC(icat,:)
         enddo
         write(u_log,*) 'temperature dependent emission factors CO for ',trim(emt%distr_nam)
         do icat = 1,emt%ncat
            write(u_log,fmt) emt%cat_code(icat),trim(emt%cat_nam(icat)),emt%temp_var_CO(icat,:)
         enddo
         write(u_log,*) '-----------------------------------------------------------------------'
      endif

    end if  ! file present


    !-------------------------------------------
    ! done
    !-------------------------------------------

    ! close:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Read


  ! ***


  subroutine LE_Emis_TNO_Clear( emt, status )

    use LE_Emis_TNO_Base   , only : LE_Emis_TNO_Base_Done
    use LE_Emis_Composition, only : LE_Emis_Composition_Done
    use LE_Emis_Scenario   , only : LE_Emis_Scenario_Done
    use LE_EMis_HeightDistr, only : LE_Emis_HeightDistr_Done
    use LE_Emis_Time_Prof  , only : LE_Emis_Time_Prof_Done
    use LE_Emis_Time_Prof_Comp, only : LE_Emis_Time_Prof_Comp_Done
    use LE_Emis_Time_Prof_Grid, only : LE_Emis_Time_Prof_Grid_Done

#ifdef with_labeling
    use SA_Labeling        , only : SA_Emis_Sectors_Done
#endif

    ! --- in/out ------------------------------

    type(T_Emis_TNO), intent(inout)    ::  emt
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Clear'

    ! --- local -------------------------------

    integer     ::  iemis

    ! --- begin -------------------------------

    ! clear country arrays:
    deallocate( emt%country_code )
    deallocate( emt%country_name )

    ! cllear category code arrays:
    deallocate( emt%cat_code )
    deallocate( emt%cat_nam  )

    ! clear profiles:
    if ( associated(emt%vocprof     ) ) deallocate( emt%vocprof      )
    if ( associated(emt%temp_var_VOC) ) deallocate( emt%temp_var_VOC )
    if ( associated(emt%temp_var_CO ) ) deallocate( emt%temp_var_CO  )

    if (allocated(emt%imonthdp) ) deallocate( emt%imonthdp, emt%idaydp, emt%ihourdp )

    ! done with base emissions:
    call LE_Emis_TNO_Base_Done( emt%emb, status )
    IF_NOTOK_RETURN(status=1)

    ! done with scenaro factors:
    call LE_Emis_Scenario_Done( emt%scenario, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over emisisons:
    do iemis = 1, emt%nemis
      ! done with compositions:
      call LE_Emis_Composition_Done( emt%emcomp(iemis), status )
      IF_NOTOK_RETURN(status=1)
    end do

    !! done with stack heights:
    !deallocate( emt%hstack )

    ! done with height distributions:
    if ( emt%with_height_distribution ) then
      call LE_Emis_HeightDistr_Done( emt%hdistr, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! done with country/snap/tracer dependent time profiles
    call LE_Emis_Time_Prof_Done(emt%time_prof, status)
    IF_NOTOK_RETURN(status=1)

    ! done with country/snap/tracer dependent time profiles
    call LE_Emis_Time_Prof_Comp_Done(emt%time_prof_comp, status)
    IF_NOTOK_RETURN(status=1)
    
    ! done with grid dependent time profiles
    call LE_Emis_Time_Prof_Grid_Done(emt%time_prof_grid, status)
    IF_NOTOK_RETURN(status=1)
    
#ifdef with_labeling
    ! done with sector labeling:
    call SA_Emis_Sectors_Done(status, tag=emt%label )
    IF_NOTOK_RETURN(status=1)
#endif

    ! clear skip list:
    deallocate( emt%skip_tracer )

    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Clear


  ! ***


  subroutine LE_Emis_TNO_Setup( emt, emis_a, t1, t2, status )

    use GO       , only : TDate, Get, operator(+), operator(-), operator(/)
    use GO       , only : GO_Timer_Start, GO_Timer_End
    use GO       , only : calc_DayNumber
    use Binas    , only : Avog
    use Num      , only : IntervalSum
    use dims     , only : runF!, outF
    use dims     , only : nx, ny, nz, nspec
    use LE_Time  , only : local_time2
    use indices  , only : specname, specunit
    use indices  , only : ispec_nh3, ispec_co, ispec_ch4
    use indices  , only : ispec_eth, ispec_ole, ispec_par, ispec_ald, ispec_form, ispec_xyl, ispec_tol
    use le_emis_tables, only : nh3fac, nh3_month
#ifdef with_labeling
    use SA_Labeling, only : SA_Emis_Setup_TNO
#endif

    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out ---------------------------

    type(T_Emis_TNO), intent(inout)    ::  emt
    real, intent(inout)                   ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Setup'

    ! number of hours in a year of 365 days:
    real, parameter  ::  hour_per_year = 8760.0

    ! this parameter used to be hard-coded, tried to explain:
    !  (g/kg)/(min/hour) = 1000.0/60.0 ~ 16.67
    real,parameter   ::  kg_per_hour_to_g_per_min = 16.67  !  (g/min)/(kg/hour)

    ! conversion from kg/year to ug/min :
    !                                                ug/kg /   hour/year   / min/hour
    real, parameter  ::  kg_per_year_to_ug_per_min = 1.0e9 / hour_per_year /  60.0  ! (ug/min)/(kg/year)

    ! conversion from #/year to #/min :
    !                                               #  /   hour/year   / min/hour
    real, parameter  ::  n_per_year_to_n_per_min = 1.0 / hour_per_year /  60.0  ! (#/min)/(#/year)

    ! --- local ----------------------------

    type(TDate)          ::  tmid
    integer              ::  yy, mm, dd, hh
    integer              ::  iemis
    integer              ::  icat
    integer              ::  icountry
    integer              ::  icomp
    integer              ::  ispec
    integer              ::  yyh, mmh, ddh, hhh, iday
    integer              ::  ix, iy
    real, allocatable    ::  emis(:)
    integer              ::  temp_index
    real                 ::  temp_fac_CO, temp_fac_VOC
    real                 ::  fact_nh3
    real,allocatable     ::  time_fact(:)
    real                 ::  compfrac
    character(len=128)   ::  conversion

    integer              ::  ihour, julian_day
    integer              ::  iz
    !real                 ::  hstack1
    real, allocatable    ::  profile(:)
    real, allocatable    ::  hhb(:)
    integer              ::  ilast
    real                 ::  emtop
    ! pre-computed height distribution profiles:
    real, allocatable    ::  hd_profiles(:,:,:,:)  ! (nx,ny,nz,ncat)
    !real, allocatable    ::  emis_a_cols(:,:)   ! (nz,nspec)

    integer               ::  ityp
    integer               ::  isource
    integer               ::  iicat
    integer               ::  itimetyp    

    real, allocatable     ::  delta_emis(:)
    ! status argument for parallel computing
    integer               ::  status_par

    ! meteo data:
    real, pointer               ::  tsurf(:,:,:)   ! (lon,lat,1)
    real, pointer               ::    h_m(:,:,:)   ! (lon,lat,lev)

    ! --- begin ----------------------------

    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)

    ! extract time parameters:
    tmid = t1 + (t2-t1)/2
    call Get( tmid, year=yy, month=mm, day=dd, hour=hh )

    ! other year then loaded ?
    if ( tmid%year /= emt%year ) then
      ! cleanup previous data:
      call LE_Emis_TNO_Clear( emt, status )
      IF_NOTOK_RETURN(status=1)
      ! read for current year:
      call LE_Emis_TNO_Read( emt, tmid%year, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! array for all emitted components:
    allocate( emis(emt%nemis) )
    
    ! array for time factors per emitted component
    allocate( time_fact(emt%nemis) )

    ! injection profile:
    allocate( profile(1:nz) )


    ! height distribution ?
    if ( emt%with_height_distribution ) then
      ! start timer:
      call GO_Timer_Start( emt%itim_height_distr, status )
      IF_NOTOK_RETURN(status=1)
      ! lookup table:
      allocate( hd_profiles(nx,ny,nz,emt%ncat) )
      ! storage for model layer heighs:
      allocate( hhb(0:nz) )
      ! emission top:
      emtop = emt%hdistr%heightb(emt%hdistr%nlay)
      ! init :
      hd_profiles = 0.0
      ! loop over grid cells:
      do iy = 1, ny
        do ix = 1, nx
          ! extra model layer heights:
          hhb(0)    = 0.0  ! m
          hhb(1:nz) = h_m(ix,iy,1:nz)   ! m
          ! loop over categories:
          do icat = 1, emt%ncat
            ! loop over model layers:
            ilast = -1
            do iz = 1, nz
              ! fraction of profile in this interval:
              call IntervalSum( emt%hdistr%heightb, emt%hdistr%fraction(icat,:), &
                                 hhb(iz-1), min(hhb(iz),emtop), &
                                 hd_profiles(ix,iy,iz,icat), ilast, status )
              IF_NOTOK_RETURN(status=1)
              ! end ?
              if ( hhb(iz) > emtop ) exit
            end do  ! model layers
          end do  ! categories
        end do ! ix
      end do ! iy
      ! clear:
      deallocate( hhb )
      ! end timer:
      call GO_Timer_End( emt%itim_height_distr, status )
      IF_NOTOK_RETURN(status=1)
    end if  ! with height distribution

    !! emissions per gridcell:
    !allocate( emis_a_cols(nz,nspec) )

    ! OpenMP notes:
    !  * Around 'isource' loop requires reduction of a
    !    full 4D array 'emis_a_loc(nx,ny,nz,npsec)',
    !    but this fails because the array seems to be too large
    !    (with a smaller array it seems to work properly).
    !  * Around the 'iemis' loop is problematic too.
    !    Some tracers could be the result of more than one emission,
    !    for example SO4a originates from both 'SOx' and 'PM'-with-composition.
    !    This would require a local copy of 'emis_a' on each thread which
    !    is then reduced by sum afterwards; these arrays contain many zero's,
    !    and therefore the computation might actually be slowed down.

    ! loop over categories:
    do icat = 1, emt%emb%ncat
      ! loop over source types:
      do ityp = 1, emt%emb%ntyp

        !! type code:
        !stype = emt%emb%typi(ityp)%code
        !! get source height and emission layer
        !select case ( stype )
        !  case ( 'A', 'a' )
        !    iz      = 1
        !    hstack1 = 0.0
        !  case ( 'P', 'p' )
        !    iz      = 1    ! use lowest layer if stack_height is zero
        !    hstack1 = emt%hstack(icat)
        !  case default
        !    write (gol,'("unsupported source type : ",a)') stype; call goErr
        !    TRACEBACK; status=1; return
        !end select
        ! use lowest layer if stack_height is zero:
        iz = 1

        ! default profile: inject everything in layer iz:
        profile = 0.0
        profile(iz) = 1.0

        ! loop over sources:
        do isource = 1, emt%emb%emg(icat,ityp)%nsource

          ! expand coordinate:
          ix       = emt%emb%emg(icat,ityp)%ix(isource)
          iy       = emt%emb%emg(icat,ityp)%iy(isource)
          icountry = emt%emb%emg(icat,ityp)%ic(isource)

          ! compute actual time for emission based on the country code
          ! also get day-of-the-week number
          call local_time2( emt%country_code(icountry), &
                              yy , mm , dd , hh , &   ! uct
                              yyh, mmh, ddh, hhh, &   ! local time  [yyh-mmh-ddh,hhh:00:00,yyh-mmh-ddh,hhh+1:00:00)
                              iday, &                  ! day-of-the-week number
                              status )
          IF_NOTOK_RETURN(status=1)
          
          ! loop over time types
          do itimetyp = 1, emt%ntimetype
            ! switch:
            select case ( trim(emt%timetype(itimetyp) ) )
              !~ unity factors:
              case ( 'unity' )
                time_fact(:) = 1.0
              !~ month/day/hour profiles: 
              case ( 'mdh' )
                ! set the time factor (equal for all components)
                time_fact(:) = emt%imonthdp(icat,mmh  ) * &
                             emt%idaydp  (icat,iday ) * &
                             emt%ihourdp (icat,hhh+1)  ! at local hour 00 (hhh=0), the time profile for the first hour is valid
              !~ profile for each hour, category, and country:
              case ( 'hourly' )
                ! set the time factor (equal for all components)
                julian_day = calc_DayNumber('gregorian',yy,mm,dd)
                ihour = (julian_day-1)*24+hh + 1
                ! Snap/country dependent time factor
                time_fact(:) = emt%time_prof%profile(ihour,icat,icountry)
              case ( 'hourly-comp' )
                ! set the time factor (component dependent)
                julian_day = calc_DayNumber('gregorian',yy,mm,dd)
                ihour = (julian_day-1)*24+hh + 1
                ! Snap/country dependent time factor
                time_fact(:) = emt%time_prof_comp%profile(ihour,icat,icountry,:)
              case ( 'hourly-gridded' )
                ! set the timefactor for grid dependency
                do iicat = 1, emt%time_prof_grid%ncat 
                  ! match icat with categories in grid dependent file
                  if ( emt%time_prof_grid%icats(iicat) == icat ) then
                    ! set the time factor (component dependent)
                    julian_day = calc_DayNumber('gregorian',yy,mm,dd)
                    ihour = (julian_day-1)*24+hh + 1
                    ! fill in timefactor from grid dependent fileS
                    time_fact(:) = emt%time_prof_grid%profile(ihour,iicat,ix,iy)
                    exit
                  end if
                end do
              !~ strange ...
              case default
                write (gol,'("unsupported timetype `",a,"`")') trim(emt%timetype(itimetyp)); call goErr
                TRACEBACK; status=1; return
            end select
            
          end do ! timetypes

          ! check ...
          if ( any(time_fact < 0.0) ) then
            write (gol,'("BUG - found negative time factor : ",e16.8)') time_fact; call goErr
            TRACEBACK; status=1; return
          end if

          !! stack heights defined ?
          !if ( emt%with_stack_height ) then
          !  ! determine the layer in which the emission takes place
          !  if ( hstack1 > 0.0 ) then
          !    do k = 1, nz
          !      if ( hstack1 < h(ix,iy,k)*1000.0 ) then
          !        iz = k
          !        exit
          !      end if
          !    end do
          !  end if
          !  ! reset distribution profile now that iz might be changed:
          !  profile = 0.0
          !  profile(iz) = 1.0
          !endif

          ! fill profile from height distribution table ?
          if ( emt%with_height_distribution ) then
            ! extract:
            profile = hd_profiles(ix,iy,:,icat)
          end if

          ! check ...
          if ( abs(sum(profile)-1.0) > 1.0e-4 ) then
            write (gol,'("sum over profile not one:")'); call goErr
            write (gol,*) profile; call goErr
            TRACEBACK; status=1; return
          end if

          ! get temperature dependent VOC and CO factors
          if ( associated(emt%temp_var_CO) .or. associated(emt%temp_var_VOC) ) then
            ! check ...
            if ( .not. (associated(emt%temp_var_CO) .and. associated(emt%temp_var_VOC)) ) then
              write (gol,'("both temp_var_CO and temp_var_VOC should be allocated")'); call goErr
              TRACEBACK; status=1; return
            end if
            ! lookup value:
            temp_index = int(tsurf(ix,iy,1)-273.15)
            temp_index = max(-30,min(40,temp_index))
            temp_fac_CO  = emt%temp_var_CO(icat,temp_index)
            temp_fac_VOC = emt%temp_var_VOC(icat,temp_index)
            if (temp_fac_VOC < 0.0) then
               print *,'warning: trying to use an uninitialized temperature factor for VOC!'
               stop
            endif
            if (temp_fac_CO < 0.0) then
               print *,'warning: trying to use an uninitialized temperature factor for CO!'
               stop
            endif
          else
            ! dummy ...
            temp_fac_CO  = 1.0
            temp_fac_VOC = 1.0
          end if

          ! Units:
          !   emt%base_emis%emis  : kg/year
          !   emis_a for gases    : mol/min
          !   emis_a for aerosols : ug/min

          ! extract emissions:
          emis = emt%emb%emg(icat,ityp)%emis(isource,:)

          ! apply scenario factor:
          emis = emis * emt%scenario%factor(icat,icountry,:)

          !! emission profile for all specs;
          !! each tread will make a private copy,
          !! the copies will be reduced by summing:
          !!emis_a_cols = 0.0

          ! status for parallel computing
          status_par = 0

          ! loop over all emitted tracers:
          !>>> don't use OpenMP in this way, see above !
          !xOMP parallel &
#ifndef __GFORTRAN__
          !xOMP   default( none ) &
          !xOMP   shared ( nh3fac, nh3_month ) &
          !xOMP   shared ( specname, specunit ) &
#endif
          !xOMP   shared ( emt ) &
          !xOMP   shared ( icat, icountry ) &
          !xOMP   shared ( mmh, hhh, nz, yy ) &
          !xOMP   shared ( ix, iy ) &
          !xOMP   shared ( emis, emis_a ) &  ! <-- use array delta_emis_a with recuction:+
          !xOMP   shared ( profile ) &
          !xOMP   shared ( time_fact ) &
          !xOMP   shared ( temp_fac_CO, temp_fac_VOC ) &
          !xOMP   private( gol ) &
          !xOMP   private( iemis, ispec ) &
          !xOMP   private( fact_nh3 ) &
          !xOMP   private( compfrac ) &
          !xOMP   private( conversion ) &
          !xOMP   private( delta_emis ) &
          !xOMP   private( status ) &
          !xOMP   reduction( + : status_par )
          !xOMP   do
          !<<<
          do iemis = 1, emt%nemis

            ! check number of components in composition;
            ! might be zero, for example if composition description included 'skip=T':
            if ( emt%emcomp(iemis)%ncomp == 0 ) cycle

            ! only positive emisisons:
            if ( emis(iemis) <= 0.0 ) cycle

            ! help array delta_emis
            if ( .not. allocated(delta_emis) ) allocate( delta_emis(1:nz) )

            ! which ?
            select case ( trim(emt%emis_name(iemis)) )

              !=================================================
              case ( 'NH3', 'nh3' )
              !=================================================

                ! NH3 emissions
                ispec = ispec_nh3
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    !thelp=tsurf(ix,iy)-273.15
                    !fact_nh3 = nh3fac(hhh)/100.0 * max(1.0,0.25*(thelp-20.0))
                    fact_nh3 = (nh3fac(hhh)/100.0) * (nh3_month(mmh)/100.0)

                    ! always inject in lowest layer, first set all layers to zero:
                    delta_emis = 0.0
                    ! use Snap/country/tracer dependent time factor
                    !delta_emis(1) =  emis(iemis) * time_fact(iemis) * kg_per_hour_to_g_per_min / (emt%xm(iemis)*1e3*hour_per_year)
                    ! For NH3 use special factor (fact_nh3)
                    delta_emis(1) =  emis(iemis) * fact_nh3  * kg_per_hour_to_g_per_min / (emt%xm(iemis)*1e3*hour_per_year)
                    ! add:
                    emis_a(ix,iy,1,ispec) = emis_a(ix,iy,1,ispec) + delta_emis(1)
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if


              !=================================================
              case ( 'CO', 'co' )
              !=================================================

                ! CO emission
                ispec = ispec_co
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    ! fill contribution:
                    delta_emis = profile * &
                             emis(iemis) * time_fact(iemis) * kg_per_hour_to_g_per_min &
                             / (emt%xm(iemis)*1e3*hour_per_year) * temp_fac_CO
                    ! add:
                    emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if


              !=================================================
              case ( 'VOC', 'nmvoc' )
              !=================================================

                ! VOCs need to be divided according to a voc-split profile
                !  emis      : kg/year
                !  vocprof   : mol/(kg VOC)
                !
                !    emis(iemis) * emt%vocprof(icat,5) * kg_per_hour_to_g_per_min/(1000.*hour_per_year)
                !    kg/year        mol/kg                 g/min hour/kg      /(g/kg  hour/year)
                !    kg             mol/kg                  /min
                !    mol/min

                !
                ! To be done: implement the voc split as a component file.
                ! Include a composit column for:
                !    PAR emission = PAR + 3*ACET + 4*KET
                ! Use country code '***' to set these for all countries,
                ! not implemented yet in LE_Emis_Composition.
                !

                ! check ..
                if ( .not. associated(emt%vocprof) ) then
                  write (gol,'("No VOC-profile allocated (file was not found ?)")'); call goErr
                  TRACEBACK; status_par=status_par+1; cycle
                end if

                ! ETH emission
                ispec = ispec_eth
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    ! fill contribution:
                    delta_emis = profile * &
                                 emis(iemis)*emt%vocprof(icat,9)*time_fact(iemis)*kg_per_hour_to_g_per_min/(1000.*hour_per_year)*temp_fac_VOC
                    ! add:
                    emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if

                ! OLE emission
                ispec = ispec_ole
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    ! fill contribution:
                    delta_emis = profile * &
                                 emis(iemis)*emt%vocprof(icat,1)*time_fact(iemis)*kg_per_hour_to_g_per_min/(1000.*hour_per_year)*temp_fac_VOC
                    ! add:
                    emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if

                ! PAR emission = PAR + 3*ACET + 4*KET
                ispec = ispec_par
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    ! fill contribution:
                    delta_emis = profile * &
                                  emis(iemis)* &
                                  ( emt%vocprof(icat,2) + 3*emt%vocprof(icat,8) + 4*emt%vocprof(icat,7) )*  &
                                   time_fact(iemis)*kg_per_hour_to_g_per_min/(1000.*hour_per_year)*temp_fac_VOC
                    ! add:
                    emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if

                ! ALD emission
                ispec = ispec_ald
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    ! fill contribution:
                    delta_emis = profile * &
                                 emis(iemis)*emt%vocprof(icat,6)*time_fact(iemis)*kg_per_hour_to_g_per_min/(1000.*hour_per_year)*temp_fac_VOC
                    ! add:
                    emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if

                ! FORM emission
                ispec = ispec_form
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    ! fill contribution:
                    delta_emis = profile * &
                         emis(iemis)*emt%vocprof(icat,5)*time_fact(iemis)*kg_per_hour_to_g_per_min/(1000.*hour_per_year)*temp_fac_VOC
                    ! add:
                    emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if

                ! XYL emission
                ispec = ispec_xyl
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    ! fill contribution:
                    delta_emis = profile * &
                                 emis(iemis)*emt%vocprof(icat,4)*time_fact(iemis)*kg_per_hour_to_g_per_min/(1000.*hour_per_year)*temp_fac_VOC
                    ! add:
                    emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if

                ! TOL emission
                ispec = ispec_tol
                ! enabled ?
                if ( ispec > 0 ) then
                  ! not explicitly skipped ?
                  if ( .not. emt%skip_tracer(ispec) ) then
                    ! fill contribution:
                    delta_emis = profile * &
                                 emis(iemis)*emt%vocprof(icat,3)*time_fact(iemis)*kg_per_hour_to_g_per_min/(1000.*hour_per_year)*temp_fac_VOC
                    ! add:
                    emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis
#ifdef with_labeling
                    call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                    if (status /= 0 ) then
                      write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                    end if
#endif
                  end if
                end if

              !=================================================
              case default
              !=================================================

                ! loop over components in composition;
                ! might be zero, for example if composition
                ! description included 'skip=T':
                do icomp = 1, emt%emcomp(iemis)%ncomp
                  ! target tracer:
                  ispec = emt%emcomp(iemis)%itracer(icomp)
                  ! not a tracer in this run ? then skip:
                  if ( ispec < 0 ) cycle
                  ! explicitly skip ?
                  if ( emt%skip_tracer(ispec) ) cycle
                  ! extract:
                  compfrac = emt%emcomp(iemis)%frac(icat,icountry,icomp)
                  ! check for undefined values ...
                  if ( compfrac < 0.0 ) then
                    write (gol,'("found undefined component fraction for emission:")'); call goErr
                    write (gol,'("  emission   : ",a)') trim(emt%emis_name(iemis)); call goErr
                    write (gol,'("  category   : ",i6)') emt%cat_code(icat); call goErr
                    write (gol,'("  country    : ",a)') emt%country_code(icountry); call goErr
                    write (gol,'("  component  : ",a)') specname(ispec); call goErr
                    TRACEBACK; status_par=status_par+1; cycle
                  end if
                  ! units of emissions and target tracer:
                  conversion = trim(emt%emb%emi(iemis)%units)//' ; '//trim(specunit(ispec))
                  ! distribute emission (in kg) over components:
                  select case ( trim(conversion) )

                    ! convert to ppb/min if tracer unit is volume mixing ratio:
                    case ( 'kg year**-1 ; ppb', 'kg year-1 ; ppb' )
                      ! check ...
                      if ( emt%xm(iemis) < 0.0 ) then
                        write (gol,'("no molecule mass defined for emitted tracer `",a,"`")') &
                               trim(emt%emis_name(iemis)); call goErr
                        TRACEBACK; status_par=status_par+1; cycle
                      end if
                      ! use Snap/country dependent time factor
                      delta_emis = profile * &
                                 compfrac * emis(iemis) * time_fact(iemis) * kg_per_hour_to_g_per_min &
                                 / (emt%xm(iemis)*1e3*hour_per_year)        ! mol/min

                    ! convert to ug/min if targer tracer unit is mass concentration:
                    case ( 'kg year**-1 ; ug/m3', 'kg year-1 ; ug/m3' )
                      ! emission increment:
                      delta_emis = profile * &
                                   compfrac * emis(iemis) * time_fact(iemis) * kg_per_year_to_ug_per_min  ! ug/min

                    ! convert to mlc/min if targer tracer unit is number concentration:
                    case ( 'kg year**-1 ; mlc/cm3', 'kg year-1 ; mlc/cm3' )
                      ! check ...
                      if ( emt%xm(iemis) < 0.0 ) then
                        write (gol,'("no molecule mass defined for emitted tracer `",a,"`")') &
                               trim(emt%emis_name(iemis)); call goErr
                        TRACEBACK; status_par=status_par+1; cycle
                      end if
                      ! emission increment:
                      !      emis(iemis) * n_per_year_to_n_per_min  / emt%xm(iemis) * Avog 
                      !        kg/year   * (mol/min)/(mol/year)     /  (kg/mol)     * (mlc/mol)
                      !        kg/year     mol/min year/mol            mol/kg         mlc/mol
                      !        mlc/min
                      !
                      delta_emis = profile * &
                                   compfrac * emis(iemis) * time_fact(iemis) * n_per_year_to_n_per_min & ! kg/min
                                   / emt%xm(iemis) * Avog   ! mlc/min

                    ! convert to #/min if targer tracer unit is number concentration:
                    case ( '1 year**-1 ; 1/cm3', '1 year-1 ; 1/cm3')
                      ! emission increment:
                      delta_emis = profile * &
                                   compfrac * emis(iemis) * time_fact(iemis) * n_per_year_to_n_per_min  ! #/min

                    ! unknown ...
                    case default
                      write (gol,'("unsupported emis;tracer units pair `",a,"` for tracer ",i3," (",a,")")') &
                             trim(trim(conversion)), ispec, trim(specname(ispec)); call goErr
                      TRACEBACK; status_par=status_par+1; cycle
                  end select

                  ! add:
                  emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis   ! mol/min, ug/min, mlc/min, #/min
#ifdef with_labeling
                  ! update labels if necessary:
                  call SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis,status, tag=emt%label)
                  if (status /= 0 ) then
                    write (gol, '(" SA_Emis_Setup_TNO")' ); call goErr
                    TRACEBACK; status_par=status_par+1; cycle
                  end if
#endif

                end do  ! components

            end select   ! switch over emitted tracers

            ! clear:
            deallocate( delta_emis )

          end do   ! emitted components
          !>>> don't use OpenMP in this way, see above
          !xOMP   end do
          !xOMP end parallel
          !<<<

          ! check on errors:
          if ( status_par /= 0 ) then
            write (gol,'("error status returned from loop over sources")'); call goErr
            TRACEBACK; status=1; return
          end if

          !! add contribution to emission array:
          !emis_a(ix,iy,:,:) = emis_a(ix,iy,:,:) + emis_a_cols

        end do   ! sources

      end do  ! source types
    end do  ! categories

    ! clear:
    deallocate( emis )
    deallocate( profile )
    deallocate( time_fact )
    ! height distribution ?
    if ( emt%with_height_distribution ) then
      deallocate( hd_profiles )
    end if
    !deallocate( emis_a_cols )

    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Setup


end module LE_Emis_TNO

