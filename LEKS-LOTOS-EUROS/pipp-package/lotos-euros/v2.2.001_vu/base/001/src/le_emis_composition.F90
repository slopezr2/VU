!###############################################################################
!
! LE_Emis_Composition - LOTOS-EUROS emission routines
!
! SUMMARY
!
!  Provides a split of an emitted bulk into different components.
!  Input:
!   - list with supported country codes: 
!        (/'ALB,'ARM',.../)
!   - list with supported cat codes: 
!        (/01,02,..,75/)
!   - list with tracer names, '-' for unsupported:
!        (/'no2','no','o3','-',../)
!  Output:
!   - number of target components over which bulk is distributed:
!       ncomp
!   - original tracer index for a component:
!       itracer(1:ncomp)
!   - lookup table with fraction as a function of category,
!     country, and target component index:
!       frac(1:ncat,1:ncountry,1:ncomp)
!
!  USAGE
!
!    type(T_Emis_Composition)   :  emcomp
!
!    ! init table using query, see below for definitions:
!    call Emis_Composition_Init( emcomp, &
!                    'file=PM_compositionV2.csv;format=csv', &
!                    'ppm_f', &
!                    country_codes, cat_codes, specnames, status )
!    if (status/=0) stop
!
!    ! loop over emissions:
!    do ..
!
!      ! obtain country index:
!      icountry = ..
!      ! obtain cat index:
!      icat = ..
!      ! obtain bulk emission:
!      embulk = ..
!
!      ! loop over target components:
!      do icomp = 1, emcomp%ncomp
!
!        ! model index of target specie:
!        itr = emcomp%itracer(icomp)
!
!        ! check ...
!        if ( emcomp%frac(icat,icountry,icomp) < 0.0 ) then
!          print *, 'found undefined value'
!          stop
!        end if
!
!        ! compute component emission:
!        emis(ix,iy,iz,itr) = embulk * emcomp%frac(icat,icountry,icomp)
!
!      end do   ! components
!
!    end do  ! emission record
!
!    ! done:
!    call Emis_Composition_Done( emcomp, status )
!    if (status/=0) stop
!      
! COMPOSITION QUERIES
!
!  A query is a ';'-seperated list with 'key=value' settings.
!  Examples:
!     ''
!     'skip=T'
!     'target_tracer=ec_f'
!     'file=PM_composition.csv;format=csv;comment=#;sep=,;default_tracer=ppm_f;allow_missing_tracers=F'
!
!  If the 'skip' key is present, then the emission is not
!  assigned to any tracer at all.
!    'skip=F'
!
!  If the 'file' key is not specified, then the emission will be simply assign
!  to a single model tracer. By default it is supposed to have the same
!  name as the emitted tracer; a different name could be specified through
!  the value assigned to the 'target_tracer' key.
!
!  To read a table from a file, use the 'file' key for the filename.
!  The format is defined below.
!  Set the 'allow_missing_tracers' flag to 'T' if it is no problem that
!  an emitted component could not be assigned to a model tracer,
!  for example because of running a model with limitted chemistry.
!  Specify a default_tracer to assign the remaining fraction to a certain
!  model tracer; if not provided, sum of fractions should be 1.0 .
!     
!
! INPUT TABLE
!
!    Should be a .csv file, probably saved from an Excel sheet:
!
!      ISO3;cat;%EC;%OC (Full Molecular Mass);%SO4;%Na;%Other mineral
!      ALB;1;0.084395256;0.058440165;0.095846356;0.015;0.746318224
!      ARM;1;0.092480775;0.126437726;0.128324079;0.015;0.63775742
!      AUT;1;0.07384835;0.13867002;0.043577376;0.015;0.728904254
!      :
!
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

module LE_Emis_Composition

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Composition
  public  ::  Emis_Composition_Init, Emis_Composition_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Composition'
  
  
  ! --- types --------------------------------
  
  ! storage for emission composition:
  type T_Emis_Composition
    ! dimensions:
    integer                   ::  ncountry
    integer                   ::  ncat
    integer                   ::  ncomp
    ! model index of target component:
    integer, allocatable      ::  itracer(:)   ! (ncomp)
    ! fraction assigned to component:
    real, allocatable         ::  frac(:,:,:)  ! (ncat,ncountry,ncomp)
  contains
    procedure   ::  Init                 => Emis_Composition_Init
    procedure   ::  Done                 => Emis_Composition_Done
  end type
  
  ! *
  
  type T_EmisComp_MassSplit
    ! dimensions:
    integer                           ::  ncountry
    integer                           ::  ncat
    ! number of substancances:
    integer                           ::  nsubstance
    ! substance names:
    character(len=32), allocatable    ::  substance_code(:)  ! (nsubstance)
    ! fraction assigned to substances:
    real, allocatable                 ::  fraction(:,:,:)  ! (ncat,ncountry,nsubstance)  [(kg substance)/(kg nmvoc)]
  contains
    procedure   ::  Init          => EmisComp_MassSplit_Init
    procedure   ::  Done          => EmisComp_MassSplit_Done
  end type T_EmisComp_MassSplit
  
  ! *
  
  type T_EmisComp_MolarMass
    ! number of substancances:
    integer                           ::  nsubstance
    ! substance code:
    character(len=8), allocatable     ::  code(:)  ! (nsubstance)
    ! substance names:
    character(len=32), allocatable    ::  name(:)  ! (nsubstance)
    ! molar weight:
    real, allocatable                 ::  xm(:)  ! (nsubstance)  [(kg substance)/(mole substance)]
    character(len=32)                 ::  xm_units
  contains
    procedure   ::  Init          => EmisComp_MolarMass_Init
    procedure   ::  Done          => EmisComp_MolarMass_Done
  end type T_EmisComp_MolarMass
  
  ! *
  
  type T_EmisComp_CarbonBonds
    ! number of substancances:
    integer                           ::  nsubstance
    ! substance code:
    character(len=8), allocatable     ::  substance_code(:)  ! (nsubstance)
    ! substance names:
    character(len=32), allocatable    ::  substance_name(:)  ! (nsubstance)
    ! number of carbon bonds:
    integer                           ::  nbond
    ! carbon bond code:
    character(len=8), allocatable     ::  bond_code(:)  ! (nbond)
    ! assignments:
    real, allocatable                 ::  bonds(:,:)  ! (nsubstance,nbond)  [(mole carbon bond)/(mole substance)]
  contains
    procedure   ::  Init          => EmisComp_CarbonBonds_Init
    procedure   ::  Done          => EmisComp_CarbonBonds_Done
  end type T_EmisComp_CarbonBonds



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine Emis_Composition_Init( self, query, emitted_tracer, &
                                       cat_codes, country_codes, tracer_names, year, &
                                       status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue
    use LE_Emis_Tools, only : Get_ShortSNAP, ShortSNAP_to_Code

    ! --- in/out ------------------------------
    
    class(T_Emis_Composition), intent(out)    ::  self
    character(len=*), intent(in)              ::  query
    character(len=*), intent(in)              ::  emitted_tracer
    character(len=*), intent(in)              ::  cat_codes(:)      ! (/'01.00','02.00',../)
    character(len=*), intent(in)              ::  country_codes(:)  ! (/'ALB',.../)
    character(len=*), intent(in)              ::  tracer_names(:)     ! (/'no2','no',../)
    integer, intent(in)                       ::  year
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Emis_Composition_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 10
    
    ! --- local -------------------------------
    
    logical                 ::  skip
    character(len=32)       ::  target_tracer
    logical                 ::  allow_missing_tracers
    
    character(len=1024)     ::  filename
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    character(len=32)       ::  fileformat
    character(len=32)       ::  default_tracer
    logical                 ::  exist

    integer                 ::  ifield_country
    integer                 ::  ifield_cat
    integer                 ::  ifield_sectorid
    logical                 ::  with_snap
    integer                 ::  itracer_emitted
    integer                 ::  icomp_emitted
    integer                 ::  itracer_default
    integer                 ::  icomp_default
    integer, allocatable    ::  icomps(:)   ! (nfield)

    integer                 ::  fu
    character(len=1024)     ::  line
    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    integer                 ::  nfield
    integer                 ::  ifield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    integer                 ::  iline
    integer                 ::  icountry
    integer                 ::  shortsnap
    character(len=8)        ::  cat_code
    integer                 ::  icat
    character(len=64)       ::  component
    integer                 ::  itracer
    real                    ::  frac
    real                    ::  fsum
    
    character(len=1024)           ::  mass_split_file
    character(len=1024)           ::  molar_mass_file
    character(len=1024)           ::  carbon_bonds_file
    type(T_EmisComp_MassSplit)    ::  mass_split
    type(T_EmisComp_MolarMass)    ::  molar_mass
    type(T_EmisComp_CarbonBonds)  ::  carbon_bonds
    integer                       ::  icomp
    integer                       ::  isubstance

    ! --- begin -------------------------------
    
    ! store currently known dimensions:
    self%ncountry = size(country_codes)
    self%ncat     = size(cat_codes)
    
    ! skip this emission ?
    skip = .false.
    call goVarValue( trim(query), ';', 'skip', '=', skip, status )
    IF_ERROR_RETURN(status=1)
    
    ! extract filename, by default empty:
    filename = ''
    call goVarValue( trim(query), ';', 'file', '=', filename, status )
    IF_ERROR_RETURN(status=1)
    
    ! 3 level tables:
    mass_split_file = ''
    call goVarValue( trim(query), ';', 'mass_split_file', '=', mass_split_file, status )
    IF_ERROR_RETURN(status=1)
    molar_mass_file = ''
    call goVarValue( trim(query), ';', 'molar_mass_file', '=', molar_mass_file, status )
    IF_ERROR_RETURN(status=1)
    carbon_bonds_file = ''
    call goVarValue( trim(query), ';', 'carbon_bonds_file', '=', carbon_bonds_file, status )
    IF_ERROR_RETURN(status=1)

    ! flag Santiago:
    !allow_missing_tracers = .false.
    allow_missing_tracers = .true.
    call goVarValue( trim(query), ';', 'allow_missing_tracers', '=', allow_missing_tracers, status )
    IF_ERROR_RETURN(status=1)
    
    ! seperation character:
    sep = ';'
    call goVarValue( trim(query), ';', 'sep', '=', sep, status )
    IF_ERROR_RETURN(status=1)
    ! comment character:
    comment = '#'
    call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    IF_ERROR_RETURN(status=1)
    ! format description:
    fileformat = 'csv'
    call goVarValue( trim(query), ';', 'format', '=', fileformat, status )
    IF_ERROR_RETURN(status=1)

    ! default tracer ?
    default_tracer = ''
    call goVarValue( trim(query), ';', 'default_tracer', '=', default_tracer, status )
    IF_ERROR_RETURN(status=1)

    ! how to fill ?
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( skip ) then   ! do not use ...
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      ! no emitted components:
      self%ncomp = 0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if ( trim(query) == 'special' ) then
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
      ! dummy, but indicate that emission should be distributed over
      ! more than zero components:
      self%ncomp = 999

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if ( len_trim(mass_split_file) > 0 ) then
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
      ! Fill fractions using 3 tables:
      !  (mol carbon bond)/(kg NMVOC) =
      !       * mass_split         (kg substance)/(kg NMVOC)
      !       / molar_mass         /((kg substance)/(mole substance))
      !       * carbon_bond        (mole carbon bond)/(mole substance)
      
      ! read mass split table:
      call mass_split%Init( mass_split_file, cat_codes, country_codes, year, status )
      IF_ERROR_RETURN(status=1)
      
      ! read molar mass table:
      call molar_mass%Init( molar_mass_file, mass_split%substance_code, status )
      IF_ERROR_RETURN(status=1)
      
      ! read molar mass table:
      call carbon_bonds%Init( carbon_bonds_file, mass_split%substance_code, status )
      IF_ERROR_RETURN(status=1)
      
      ! copy dimensions:
      self%ncountry = mass_split%ncountry
      self%ncat     = mass_split%ncat
      self%ncomp    = carbon_bonds%nbond
      
      ! storage for mapping to model tracers:
      allocate( self%itracer(self%ncomp), stat=status )
      IF_ERROR_RETURN(status=1)
      ! loop:
      do icomp = 1, self%ncomp
        ! match carbon bond with tracer names:
        call goMatchValue( trim(carbon_bonds%bond_code(icomp)), tracer_names, self%itracer(icomp), status )
        IF_NOTOK_RETURN(status=1)
      end do
      
      ! storage for composition:
      allocate( self%frac(self%ncat,self%ncountry,self%ncomp), stat=status )
      IF_ERROR_RETURN(status=1)
      ! loop:
      do icat = 1, self%ncat
        do icountry = 1, self%ncountry
          do icomp = 1, self%ncomp
            ! init sum:
            self%frac(icat,icountry,icomp) = 0.0   ! (mole carbon bond)/(kg NMVOC)
            ! loop over intermediate substances:
            do isubstance = 1, mass_split%nsubstance
              ! add contribution:
              self%frac(icat,icountry,icomp) = self%frac(icat,icountry,icomp) + &    ! (mole carbon bond)/(kg NMVOC)
                           mass_split%fraction(icat,icountry,isubstance) & !  (kg substance)/(kg NMVOC)
                           / molar_mass%xm(isubstance) &                   !  / ((kg substance)/(mole substance))
                           * carbon_bonds%bonds(isubstance,icomp)          !  * (mole carbon bond)/(mole substance)
            end do ! substances
          end do ! icomp
        end do ! icountry
      end do ! icat     
      
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if ( len_trim(filename) == 0 ) then   ! no file, fill with identity
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      ! target tracer, by default the same name as the emitted tracer:
      target_tracer = trim(emitted_tracer)
      ! extract from the key if defined in there:
      call goVarValue( trim(query), ';', 'target_tracer', '=', target_tracer, status )
      IF_ERROR_RETURN(status=1)

      ! assign everthing to target tracer; get the index:
      call goMatchValue( trim(target_tracer), tracer_names, itracer_emitted, status, quiet=.true. )
      IF_ERROR_RETURN(status=1)
      if ( itracer_emitted < 0 ) then
        if ( allow_missing_tracers ) then
          write (gol,'("WARNING - target tracer `",a,"` not found; allowed ...")') trim(target_tracer); call goPr
        else
          write (gol,'("emitted tracer `",a,"` not found and `allow_missing_tracers` flag not set ...")') trim(target_tracer); call goErr
          TRACEBACK; status=1; return
        end if
      end if

      ! single component:
      self%ncomp = 1
      ! index of emitted component:
      icomp_emitted = 1    

      ! storage for tracer indices:
      allocate( self%itracer(self%ncomp) )
      ! assign to emitted component:
      self%itracer = -1
      self%itracer(icomp_emitted) = itracer_emitted

      ! storage for fractions:
      allocate( self%frac(self%ncat,self%ncountry,self%ncomp) )
      ! by default assign everything to emitted component:
      self%frac = 0.0
      self%frac(:,:,icomp_emitted) = 1.0
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else  ! read table from file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! file should be present:
      inquire( file=trim(filename), exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! switch per format:
      select case ( trim(fileformat) )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'csv' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! new file unit:
          call goGetFU( fu, status )
          IF_NOTOK_RETURN(status=1)

          ! open file:
          open( fu, file=trim(filename), status='old', form='formatted', iostat=status )
          if (status/=0) then
            write (gol,'("opening file : ",a)') trim(filename); call goErr
            TRACEBACK; status=1; return
          end if

          ! line counter:          
          iline = 0

          ! read header line after first comment:
          do
            ! read line:
            read (fu,'(a)',iostat=status) line
            if (status/=0) then
              write (gol,'("reading header line from file : ",a)') trim(filename); call goErr
              TRACEBACK; status=1; return
            end if
            ! empty ? then skip:
            if ( len_trim(line) == 0 ) cycle
            ! comment ? then skip:
            if ( line(1:1) == comment ) cycle
            ! found non-comment line, leave loop:
            exit
          end do
          ! split:
          call goSplitString( line, nheader, headers, status, sep=sep )
          IF_NOTOK_RETURN(status=1)

          ! storage for tracer indices, just ensure that it is large enough:
          allocate( self%itracer(nheader) )
          ! not assigned yet:
          self%itracer = -1

          ! set variables from header; first dummy values:
          ifield_country = -1
          ifield_cat = -1
          with_snap  = .false.
          ifield_sectorid = -1
          self%ncomp = 0
          ! storage to map from field index to component index:
          allocate( icomps(nheader) ) ; icomps = -1
          ! loop over fields:
          do ifield = 1, nheader
            ! current:
            header = headers(ifield) 
            ! which column ?
            select case ( trim(header) )
              !~ country code:
              case ( 'ISO3', 'country' )
                ! store field index:
                ifield_country = ifield
              !~ category code:
              case ( 'SNAP', 'SNAP1' )
                ! store field index:
                ifield_cat = ifield
                ! flag:
                with_snap = .true.
              !~ category code:
              case ( 'GNFR_Category', 'GNFR_Sector' )
                ! store field index:
                ifield_cat = ifield
              !~ sector_id, used in emission from emis-model :
              case ( 'sector_id' )
                ! store field index:
                ifield_sectorid = ifield
              !~ not used ...
              case ( 'Year' )
                ! ignore ...
              !~ tracer ...
              case default
                ! increase counter:
                self%ncomp = self%ncomp + 1
                ! assign to field index:
                icomps(ifield) = self%ncomp
                ! set name of target component; depends on emitted tracer:
                select case ( trim(emitted_tracer) )
                  !~ NOx emissions
                  case ( 'nox' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_NO'  ) ; component = 'no'
                      case ( 'Fraction_NO2' ) ; component = 'no2'
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ SOx emissions
                  case ( 'sox', 'so2' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_SO2'  ) ; component = 'so2'
                      case ( 'Fraction_SO4a' ) ; component = 'so4a_f'
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ fine mode aerosol
                  case ( 'ppm_f', 'pm2_5', 'pm25' )
#ifdef with_m7
                    !Santiago
                    select case ( trim(header) )
                      case ( 'Fraction_EC'                      , 'EC_fine'     ) ; component = 'bcks'
                      case ( 'Fraction_OC (Full Molecular Mass)', 'OC_fine'     ) ; component = 'ocks'
                      case ( 'Fraction_SO4'                     , 'SO4_fine'    ) ; component = 'so4ks'
                      case ( 'Fraction_Na'                      , 'Na_fine'     ) ; component = 'naks'
                      case ( 'Fraction_Other mineral'           , 'OthMin_fine' ) ; component = 'ssks'
                     
                
#else
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_EC'                      , 'EC_fine'     ) ; component = 'ec_f'
                      case ( 'Fraction_OC (Full Molecular Mass)', 'OC_fine'     ) ; component = 'pom_f'
                      case ( 'Fraction_SO4'                     , 'SO4_fine'    ) ; component = 'so4a_f'
                      case ( 'Fraction_Na'                      , 'Na_fine'     ) ; component = 'na_f'
                      case ( 'Fraction_Other mineral'           , 'OthMin_fine' ) ; component = 'ppm_f'
#endif
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ coarse mode aerosol
                  case ( 'ppm_c', 'pm25_pm10' )
                    ! switch per emitted component:
                    select case ( trim(header) )
#ifdef with_m7
                      case ( 'Fraction_EC'                      , 'EC_coarse'     ) ; component = 'bccs'
                      case ( 'Fraction_OC (Full Molecular Mass)', 'OC_coarse'     ) ; component = 'occs'
                      case ( 'Fraction_SO4'                     , 'SO4_coarse'    ) ; component = 'so4cs'
                      case ( 'Fraction_Na'                      , 'Na_coarse'     ) ; component = 'nacs'
                      case ( 'Fraction_Other mineral'           , 'OthMin_coarse' ) ; component = 'sscs'

#else
                      case ( 'Fraction_EC'                      , 'EC_coarse'     ) ; component = 'ec_c'
                      case ( 'Fraction_OC (Full Molecular Mass)', 'OC_coarse'     ) ; component = 'pom_c'
                      case ( 'Fraction_SO4'                     , 'SO4_coarse'    ) ; component = 'so4a_c'
                      case ( 'Fraction_Na'                      , 'Na_coarse'     ) ; component = 'na_c'
                      case ( 'Fraction_Other mineral'           , 'OthMin_coarse' ) ; component = 'ppm_c'
#endif
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ lead
                  case ( 'pb' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_Pb_f'                     ) ; component = 'pb_f'
                      case ( 'Fraction_Pb_c'                     ) ; component = 'pb_c'
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ cadmium
                  case ( 'cd' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_Cd_f'                     ) ; component = 'cd_f'
                      case ( 'Fraction_Cd_c'                     ) ; component = 'cd_c'
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ unknown ...
                  case default
                    write (gol,'("unsupported emitted tracer `",a,"`")') trim(emitted_tracer); call goErr
                    TRACEBACK; status=1; return
                end select  ! tracer
                ! search:
                call goMatchValue( component, tracer_names, itracer, status )
                IF_ERROR_RETURN(status=1)
                if ( itracer < 0 ) then
                  if ( allow_missing_tracers ) then
                    write (gol,'("WARNING - emitted tracer `",a,"` not found; allowed ...")') &
                                   trim(component); call goPr
                  else
                    write (gol,'("emitted tracer `",a,"` not found and `allow_missing_tracers` flag not set ...")') &
                                   trim(component); call goErr
                    TRACEBACK; status=1; return
                  end if
                end if
                ! store model tracer index:
                self%itracer(self%ncomp) = itracer
              !~
            end select  ! header
          end do  ! fields

          ! check index:
          if ( ifield_cat < 0 .and. ifield_sectorid < 0) then
            write (gol,'("category index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return
          end if

          ! check index:
          if ( ifield_country < 0 ) then
            write (gol,'("country index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return
          end if

          ! index for default tracer:
          itracer_default = -1
          ! default tracer specified ?
          if ( len_trim(default_tracer) > 0 ) then
            ! search:
            call goMatchValue( trim(default_tracer), tracer_names, itracer_default, status )
            IF_ERROR_RETURN(status=1)
            if ( itracer_default < 0 ) then
              if ( allow_missing_tracers ) then
                write (gol,'("WARNING - emitted tracer `",a,"` not found; allowed ...")') &
                               trim(default_tracer); call goPr
              else
                write (gol,'("emitted tracer `",a,"` not found and `allow_missing_tracers` flag not set ...")') &
                               trim(default_tracer); call goErr
                TRACEBACK; status=1; return
              end if
            end if
            ! increase counter:
            self%ncomp = self%ncomp + 1
            ! set index for default tracer:
            icomp_default = self%ncomp
            ! store model tracer index:
            self%itracer(icomp_default) = itracer_default
          end if

          ! storage for fractions:
          allocate( self%frac(self%ncat,self%ncountry,self%ncomp) )
          ! by default assign a warning value ...
          self%frac = -999.9

          ! loop over records:
          do

            ! increase record counter:
            iline = iline + 1
            ! try to read line:
            read (fu,'(a)',iostat=status) line
            ! eof ?
            if (status<0) exit
            ! error ?
            if (status>0) then
              write (gol,'("reading line ",i6," from file : ",a)') iline, trim(filename); call goErr
              TRACEBACK; status=1; return
            end if

            ! empty ? then skip:
            if ( len_trim(line) == 0 ) cycle
            ! comment ? then skip:
            if ( line(1:1) == comment ) cycle

            ! split into records:
            call goSplitString( line, nfield, fields, status, sep=sep )
            IF_NOTOK_RETURN(status=1)
            ! check ...
            if ( nfield /= nheader ) then
              write (gol,'("number of fields (",i2,") in line ",i2," :")') nfield, iline; call goPr
              write (gol,'("  ",a)') trim(line); call goErr
              write (gol,'("fields:")'); call goErr
              do ifield = 1, nfield
                write (gol,'(i6," : ",a)') ifield, trim(fields(ifield)); call goErr
              end do
              write (gol,'("differs from number of headers ",i2," in file ",a)') nheader, trim(filename); call goErr
              TRACEBACK; status=1; return
            end if

            ! get country code:
            call goMatchValue( fields(ifield_country), country_codes, icountry, status )
            IF_NOTOK_RETURN(status=1)
            
            if ( ifield_sectorid > 0 ) then
              ! convert category code to integer (1,2,71,..)
              call Get_ShortSNAP( fields(ifield_sectorid), shortsnap, status )
              IF_NOTOK_RETURN(status=1)
              ! convert to '01.00' etc:
              call ShortSNAP_to_Code( shortsnap, cat_code, status )
              IF_NOTOK_RETURN(status=1)
              ! search ...
              call goMatchValue( cat_code, cat_codes, icat, status )
              IF_NOTOK_RETURN(status=1)
            else
              ! switch:
              if ( with_snap ) then
                ! convert category code to integer (1,2,71,..)
                call Get_ShortSNAP( fields(ifield_cat), shortsnap, status )
                IF_NOTOK_RETURN(status=1)
                ! convert to '01.00' etc:
                call ShortSNAP_to_Code( shortsnap, cat_code, status )
                IF_NOTOK_RETURN(status=1)
              else
                ! copy:
                cat_code = trim(fields(ifield_cat))
              end if
              ! search ...
              call goMatchValue( cat_code, cat_codes, icat, status )
              IF_NOTOK_RETURN(status=1)
            end if
            
            ! check ...
            if ( any(self%frac(icat,icountry,:) >= 0.0) ) then
              write (gol,'("found second record for country `",a,"` category ",i4)') &
                              country_codes(icountry), cat_codes(icat); call goErr
              write (gol,'("line    : ",i6)') iline; call goErr
              write (gol,'("file    : ",a)') trim(filename); call goErr
              TRACEBACK; status=1; return
            end if              

            ! default tracer specified ?
            if ( itracer_default > 0 ) then
              ! initialize values, assign everything to emitted if possible:
              self%frac(icat,icountry,:            ) = 0.0
              self%frac(icat,icountry,icomp_default) = 1.0
            end if

            ! loop over fields:
            do ifield = 1, nfield

              ! no component ? then skip:
              if ( icomps(ifield) < 0 ) cycle

              ! current field characters:
              field = fields(ifield)            
              ! read fraction:
              read (field,*,iostat=status) frac
              if (status/=0) then
                write (gol,'("reading fraction from `",a,"`")') field; call goErr
                TRACEBACK; status=1; return
              end if

              ! check for strange values ...
              if ( (frac < 0.0) .or. (frac > 1.0) ) then
                write (gol,*) 'found strange fraction: ', frac; call goErr
                write (gol,'("fields  : ",a)') trim(line); call goErr
                write (gol,'("line    : ",i6)') iline; call goErr
                write (gol,'("file    : ",a)') trim(filename); call goErr
                TRACEBACK; status=1; return
              end if

              ! store:
              self%frac(icat,icountry,icomps(ifield)) = frac

              ! default tracer specified ?
              if ( itracer_default > 0 ) then
                ! remove contribution of current component:
                self%frac(icat,icountry,icomp_default) = &
                      self%frac(icat,icountry,icomp_default) - frac
              end if

            end do  ! fields

            ! current sum:
            fsum = sum(self%frac(icat,icountry,:))
            ! check ...
            if ( abs(fsum-1.0) > 1.0e-4 ) then
              write (gol,'("difference in fractions is too large : ")'); call goErr
              write (gol,'("  country  : ",i6," ",a)') icountry, country_codes(icountry); call goErr
              write (gol,'("  category : ",i6," ",a)') icat, cat_codes(icat); call goErr
              write (gol,'("  records  : ")'); call goErr
              do ifield = 1, nfield
                write (gol,*) ifield, ' : ', fields(ifield); call goErr
              end do
              write (gol,'("  extracted : ")'); call goErr
              do ifield = 1, self%ncomp
                write (gol,*) ifield, ' : ', self%frac(icat,icountry,ifield); call goErr
              end do
              write (gol,'("  sum      : ")'); call goErr
              write (gol,*) fsum; call goErr
              TRACEBACK; status=1; return
            end if
            ! scale all emisisons:
            self%frac(icat,icountry,:) = self%frac(icat,icountry,:) / fsum

          end do  ! lines
          
          !! check ...
          !if ( any(self%frac < 0.0) ) then
          !  write (gol,'("some values were not set in fraction table ...")'); call goErr
          !  write (gol,'("first undefined record : ")'); call goErr
          !  do icat = 1, self%ncat
          !    do icountry = 1, self%ncountry
          !      if ( any(self%frac(icat,icountry,:) < 0.0) ) then
          !        write (gol,'("  category code : ",i6)') cat_codes(icat); call goErr
          !        write (gol,'("  country code  : ",a)') trim(country_codes(icountry)); call goErr
          !        TRACEBACK; status=1; return
          !      end if
          !    end do
          !  end do
          !end if

          ! clear:
          deallocate( icomps )

          ! close
          close( fu, iostat=status )
          if (status/=0) then
            write (gol,'("closing file : ",a)') trim(filename); call goErr
            TRACEBACK; status=1; return
          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          write (gol,'("unsupported format : ",a)') trim(fileformat); call goErr
          TRACEBACK; status=1; return

      end select

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! ok
    status = 0

  end subroutine Emis_Composition_Init
  
  
  ! ***
  

  subroutine Emis_Composition_Done( self, status )
  
    ! --- in/out ------------------------------
    
    class(T_Emis_Composition), intent(inout)  ::  self
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Emis_Composition_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    if ( allocated(self%itracer) ) deallocate( self%itracer )
    if ( allocated(self%frac   ) ) deallocate( self%frac    )

    ! ok
    status = 0

  end subroutine Emis_Composition_Done



  ! ===============================================================
  ! ===
  ! === mass split table
  ! ===
  ! ===============================================================
  

  subroutine EmisComp_MassSplit_Init( self, filename, cat_codes, country_codes, year, status )

    use GO            , only :  goGetFU
    use GO            , only :  goVarValue, goSplitString, goMatchValue

    ! --- in/out ------------------------------
    
    class(T_EmisComp_MassSplit), intent(out)  ::  self
    character(len=*), intent(in)              ::  filename
    character(len=*), intent(in)              ::  cat_codes(:)      ! (/'A','B',../)
    character(len=*), intent(in)              ::  country_codes(:)  ! (/'ALB',.../)
    integer, intent(in)                       ::  year
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/EmisComp_MassSplit_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 100
    
    ! --- local -------------------------------
    
    character(len=1)        ::  comment
    character(len=1)        ::  sep

    integer                 ::  ifield_year
    integer                 ::  ifield_country
    integer                 ::  ifield_cat
    integer                 ::  offset_substance

    integer                 ::  fu
    character(len=1024)     ::  line
    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    integer                 ::  nfield
    integer                 ::  ifield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    integer                 ::  iline
    integer                 ::  rec_year
    integer                 ::  icountry
    integer                 ::  cat_code
    integer                 ::  icat
    integer                 ::  isubstance
    real                    ::  frac
    real                    ::  fsum
    
    ! --- begin -------------------------------
    
    ! store currently known dimensions:
    self%ncountry = size(country_codes)
    self%ncat     = size(cat_codes)
    
    ! seperation character:
    sep = ';'
    !call goVarValue( trim(query), ';', 'sep', '=', sep, status )
    !IF_ERROR_RETURN(status=1)

    ! comment character:
    comment = '#'
    !call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    !IF_ERROR_RETURN(status=1)

    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( fu, file=trim(filename), status='old', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! line counter:          
    iline = 0

    ! read header line after first comment:
    do
      ! read line:
      read (fu,'(a)',iostat=status) line
      if (status/=0) then
        write (gol,'("reading header line from file : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle
      ! found non-comment line, leave loop:
      exit
    end do
    ! split:
    call goSplitString( line, nheader, headers, status, sep=sep )
    IF_NOTOK_RETURN(status=1)

    ! set variables from header; first dummy values:
    ifield_year      = -1
    ifield_country   = -1
    ifield_cat       = -1
    offset_substance = -1
    ! loop over fields:
    do ifield = 1, nheader
      ! current:
      header = headers(ifield) 
      ! which column ?
      select case ( trim(header) )
        !~ year:
        case ( 'Year' )
          ! store field index:
          ifield_year = ifield
        !~ country code:
        case ( 'ISO3', 'country' )
          ! store field index:
          ifield_country = ifield
        !~ category code:
        case ( 'GNFR_Category' )
          ! store field index:
          ifield_cat = ifield
        !~ substance ...
        case default
          ! not initialized yet?
          if ( offset_substance < 0 ) then
            ! set offset:
            offset_substance = ifield - 1
            ! assumed length:
            self%nsubstance = nheader - offset_substance
            ! storage:
            allocate( self%substance_code(self%nsubstance), stat=status )
            IF_NOTOK_RETURN(status=1)
          end if
          ! target substance:
          isubstance = ifield - offset_substance
          ! store:
          self%substance_code(isubstance) = trim(header)
      end select  ! header
    end do  ! fields

    ! check index:
    if ( ifield_cat < 0 ) then
      write (gol,'("category index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
      write (gol,'("  ",a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if

    ! check index:
    if ( ifield_country < 0 ) then
      write (gol,'("country index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
      write (gol,'("  ",a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage for fractions:
    allocate( self%fraction(self%ncat,self%ncountry,self%nsubstance), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! by default assign a warning value ...
    self%fraction = -999.9

    ! loop over records:
    do

      ! increase record counter:
      iline = iline + 1
      ! try to read line:
      read (fu,'(a)',iostat=status) line
      ! eof ?
      if (status<0) exit
      ! error ?
      if (status>0) then
        write (gol,'("reading line ",i6," from file : ",a)') iline, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle

      ! split into records:
      call goSplitString( line, nfield, fields, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( nfield /= nheader ) then
        write (gol,'("number of fields (",i2,") in line ",i2," :")') nfield, iline; call goPr
        write (gol,'("  ",a)') trim(line); call goErr
        write (gol,'("fields:")'); call goErr
        do ifield = 1, nfield
          write (gol,'(i6," : ",a)') ifield, trim(fields(ifield)); call goErr
        end do
        write (gol,'("differs from number of headers ",i2," in file ",a)') nheader, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! extract year:
      field = trim(fields(ifield_year))
      read (field,*,iostat=status) rec_year
      if (status/=0) then
        write (gol,'("reading year from `",a,"`")') field; call goErr
        TRACEBACK; status=1; return
      end if
      ! filter ...
      if ( rec_year /= year ) cycle

      ! get country code:
      call goMatchValue( fields(ifield_country), country_codes, icountry, status )
      IF_NOTOK_RETURN(status=1)

      ! index in category code list:
      call goMatchValue( fields(ifield_cat), cat_codes, icat, status )
      IF_NOTOK_RETURN(status=1)

      ! check ...
      if ( any(self%fraction(icat,icountry,:) >= 0.0) ) then
        write (gol,'("found second record for country `",a,"` category `",a,"`")') &
                        trim(country_codes(icountry)), trim(cat_codes(icat)); call goErr
        write (gol,'("line    : ",i6)') iline; call goErr
        write (gol,'("file    : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! loop over substance fields:
      do ifield = offset_substance+1, nfield

        ! current field characters:
        field = fields(ifield)            
        ! read fraction:
        read (field,*,iostat=status) frac
        if (status/=0) then
          write (gol,'("reading fraction from `",a,"`")') field; call goErr
          TRACEBACK; status=1; return
        end if

        ! check for strange values ...
        if ( (frac < 0.0) .or. (frac > 1.0) ) then
          write (gol,*) 'found strange fraction: ', frac; call goErr
          write (gol,'("fields  : ",a)') trim(line); call goErr
          write (gol,'("line    : ",i6)') iline; call goErr
          write (gol,'("file    : ",a)') trim(filename); call goErr
          TRACEBACK; status=1; return
        end if

        ! store:
        self%fraction(icat,icountry,ifield-offset_substance) = frac

      end do  ! fields

      ! current sum:
      fsum = sum(self%fraction(icat,icountry,:))
      ! any value at al?
      if ( fsum > 0.0 ) then
        ! sum should be close to 1.0 ...
        if ( abs(fsum-1.0) > 1.0e-4 ) then
          write (gol,'("difference in fractions is too large : ")'); call goErr
          write (gol,'("  country  : ",i6," ",a)') icountry, country_codes(icountry); call goErr
          write (gol,'("  category : ",i6," ",a)') icat, cat_codes(icat); call goErr
          write (gol,'("  records  : ")'); call goErr
          do ifield = 1, nfield
            write (gol,*) '    ', trim(headers(ifield)), ' : ', fields(ifield); call goErr
          end do
          write (gol,'("  extracted : ")'); call goErr
          do isubstance = 1, self%nsubstance
            write (gol,*) isubstance, ' : ', self%fraction(icat,icountry,isubstance); call goErr
          end do
          write (gol,'("  sum      : ")'); call goErr
          write (gol,*) fsum; call goErr
          TRACEBACK; status=1; return
        end if
        ! scale to eliminate rounding errors:
        self%fraction(icat,icountry,:) = self%fraction(icat,icountry,:) / fsum
      end if

    end do  ! lines

    ! close
    close( fu, iostat=status )
    if (status/=0) then
      write (gol,'("closing file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( any(self%fraction < 0.0) ) then
      !! while testing, just reset ...
      !self%fraction = max( 0.0, self%fraction )
      ! info ..
      write (gol,'("mass split fraction not defined for some categories/countries:")'); call goErr
      ! loop:
      do icat = 1, self%ncat
        do icountry = 1, self%ncountry
          if ( any( self%fraction(icat,icountry,:) < 0.0 ) ) then
            write (gol,*) '  ', trim(cat_codes(icat)), ' ', trim(country_codes(icountry)); call goErr
          end if
        end do
      end do
      write (gol,'("file: ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine EmisComp_MassSplit_Init
  
  
  ! ***
  

  subroutine EmisComp_MassSplit_Done( self, status )
  
    ! --- in/out ------------------------------
    
    class(T_EmisComp_MassSplit), intent(inout)    ::  self
    integer, intent(out)                          ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/EmisComp_MassSplit_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate( self%substance_code, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%fraction, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine EmisComp_MassSplit_Done



  ! ===============================================================
  ! ===
  ! === molar mass table
  ! ===
  ! ===============================================================
  

  subroutine EmisComp_MolarMass_Init( self, filename, substance_codes, status )

    use GO            , only :  goGetFU
    use GO            , only :  goVarValue, goSplitString, goMatchValue

    ! --- in/out ------------------------------
    
    class(T_EmisComp_MolarMass), intent(out)  ::  self
    character(len=*), intent(in)              ::  filename
    character(len=*), intent(in)              ::  substance_codes(:)   ! (/'101','102',../)
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/EmisComp_MolarMass_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 3
    
    ! --- local -------------------------------
    
    character(len=1)        ::  comment
    character(len=1)        ::  sep

    integer                 ::  ifield_code
    integer                 ::  ifield_name
    integer                 ::  ifield_xm

    integer                 ::  fu
    character(len=1024)     ::  line
    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    integer                 ::  nfield
    integer                 ::  ifield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    integer                 ::  iline
    integer                 ::  isubstance
    real                    ::  value
    
    ! --- begin -------------------------------
    
    ! store currently known dimensions:
    self%nsubstance = size(substance_codes)

    ! storage:
    allocate( self%code(self%nsubstance), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%name(self%nsubstance), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%xm(self%nsubstance), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! not assigned yet:
    self%code = ''
    self%name = ''
    self%xm = -999.9

    ! seperation character:
    sep = ';'
    !call goVarValue( trim(query), ';', 'sep', '=', sep, status )
    !IF_ERROR_RETURN(status=1)

    ! comment character:
    comment = '#'
    !call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    !IF_ERROR_RETURN(status=1)

    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( fu, file=trim(filename), status='old', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! line counter:          
    iline = 0

    ! read header line after first comment:
    do
      ! read line:
      read (fu,'(a)',iostat=status) line
      if (status/=0) then
        write (gol,'("reading header line from file : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle
      ! found non-comment line, leave loop:
      exit
    end do
    ! split:
    call goSplitString( line, nheader, headers, status, sep=sep )
    IF_NOTOK_RETURN(status=1)

    ! set variables from header; first dummy values:
    ifield_code = -1
    ifield_name = -1
    ifield_xm   = -1
    ! loop over fields:
    do ifield = 1, nheader
      ! current:
      header = headers(ifield) 
      ! which column ?
      select case ( trim(header) )
        !~ description:
        case ( 'substance code' )
          ! store field index:
          ifield_code = ifield
        !~ description:
        case ( 'substance name' )
          ! store field index:
          ifield_name = ifield
        !~ value:
        case ( 'molecular weight' )
          ! store field index:
          ifield_xm = ifield
        !~ unknown ...
        case default
          write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
          TRACEBACK; status=1; return
      end select  ! header
    end do  ! fields

    ! check index:
    if ( ifield_code < 0 ) then
      write (gol,'("substance code index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
      write (gol,'("  ",a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if

    ! check index:
    if ( ifield_name < 0 ) then
      write (gol,'("substance name index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
      write (gol,'("  ",a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if

    ! check index:
    if ( ifield_xm < 0 ) then
      write (gol,'("molar mass index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
      write (gol,'("  ",a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if

    ! loop over records:
    do

      ! increase record counter:
      iline = iline + 1
      ! try to read line:
      read (fu,'(a)',iostat=status) line
      ! eof ?
      if (status<0) exit
      ! error ?
      if (status>0) then
        write (gol,'("reading line ",i6," from file : ",a)') iline, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle

      ! split into records:
      call goSplitString( line, nfield, fields, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( nfield /= nheader ) then
        write (gol,'("number of fields (",i2,") in line ",i2," :")') nfield, iline; call goPr
        write (gol,'("  ",a)') trim(line); call goErr
        write (gol,'("fields:")'); call goErr
        do ifield = 1, nfield
          write (gol,'(i6," : ",a)') ifield, trim(fields(ifield)); call goErr
        end do
        write (gol,'("differs from number of headers ",i2," in file ",a)') nheader, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! get code:
      field = fields(ifield_code)
      ! get index:
      call goMatchValue( field, substance_codes, isubstance, status )
      IF_NOTOK_RETURN(status=1)
      ! store:
      self%code(isubstance) = trim(field)
      
      ! get name:
      field = fields(ifield_name)
      ! store:
      self%name(isubstance) = trim(field)
      
      ! value field:
      field = fields(ifield_xm)
      ! read value:
      read (field,*,iostat=status) value
      if (status/=0) then
        write (gol,'("reading fraction from `",a,"`")') field; call goErr
        TRACEBACK; status=1; return
      end if

      ! check for strange values ...
      if ( value < 0.0 ) then
        write (gol,*) 'found strange molar mass: ', value; call goErr
        write (gol,'("fields  : ",a)') trim(line); call goErr
        write (gol,'("line    : ",i6)') iline; call goErr
        write (gol,'("file    : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! store, convert from g/mole to kg/mole:
      self%xm(isubstance) = value/1.0e3
      self%xm_units = 'kg/mole'

    end do  ! lines

    ! close
    close( fu, iostat=status )
    if (status/=0) then
      write (gol,'("closing file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( any(self%xm < 0.0) ) then
      ! info ..
      write (gol,'("molar masses not defined for some substance codes:")'); call goErr
      ! loop:
      do isubstance = 1, self%nsubstance
        write (gol,'("  ",i3," ",a4,f12.4," ",a)') isubstance, &
                trim(substance_codes(isubstance)), self%xm(isubstance), trim(self%name(isubstance)); call goErr
      end do  ! substances
      write (gol,'("file: ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine EmisComp_MolarMass_Init
  
  
  ! ***
  

  subroutine EmisComp_MolarMass_Done( self, status )
  
    ! --- in/out ------------------------------
    
    class(T_EmisComp_MolarMass), intent(inout)    ::  self
    integer, intent(out)                          ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/EmisComp_MolarMass_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate( self%code, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%name, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%xm, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine EmisComp_MolarMass_Done



  ! ===============================================================
  ! ===
  ! === carbon bonds table
  ! ===
  ! ===============================================================
  

  subroutine EmisComp_CarbonBonds_Init( self, filename, substance_codes, status )

    use GO            , only :  goGetFU
    use GO            , only :  goVarValue, goSplitString, goMatchValue

    ! --- in/out ------------------------------
    
    class(T_EmisComp_CarbonBonds), intent(out)  ::  self
    character(len=*), intent(in)                ::  filename
    character(len=*), intent(in)                ::  substance_codes(:)   ! (/'101','102',../)
    integer, intent(out)                        ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/EmisComp_CarbonBonds_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 20
    
    ! --- local -------------------------------
    
    character(len=1)        ::  comment
    character(len=1)        ::  sep

    integer                 ::  ifield_code
    integer                 ::  ifield_name
    integer                 ::  offset_bond

    integer                 ::  fu
    character(len=1024)     ::  line
    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    integer                 ::  nfield
    integer                 ::  ifield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    integer                 ::  iline
    integer                 ::  isubstance
    integer                 ::  ibond
    real                    ::  value
    
    ! --- begin -------------------------------
    
    ! store currently known dimensions:
    self%nsubstance = size(substance_codes)

    ! storage:
    allocate( self%substance_code(self%nsubstance), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%substance_name(self%nsubstance), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! not assigned yet:
    self%substance_code = ''
    self%substance_name = ''

    ! seperation character:
    sep = ';'
    !call goVarValue( trim(query), ';', 'sep', '=', sep, status )
    !IF_ERROR_RETURN(status=1)

    ! comment character:
    comment = '#'
    !call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    !IF_ERROR_RETURN(status=1)

    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( fu, file=trim(filename), status='old', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! line counter:          
    iline = 0

    ! read header line after first comment:
    do
      ! read line:
      read (fu,'(a)',iostat=status) line
      if (status/=0) then
        write (gol,'("reading header line from file : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle
      ! found non-comment line, leave loop:
      exit
    end do
    ! split:
    call goSplitString( line, nheader, headers, status, sep=sep )
    IF_NOTOK_RETURN(status=1)

    ! set variables from header; first dummy values:
    ifield_code = -1
    ifield_name = -1
    offset_bond = -1
    ! loop over fields:
    do ifield = 1, nheader
      ! current:
      header = headers(ifield) 
      ! which column ?
      select case ( trim(header) )
        !~ description:
        case ( 'substance code' )
          ! store field index:
          ifield_code = ifield
        !~ description:
        case ( 'substance name' )
          ! store field index:
          ifield_name = ifield
        !~ carbon bonds ...
        case default
          ! not initialized yet?
          if ( offset_bond < 0 ) then
            ! set offset:
            offset_bond = ifield - 1
            ! assumed length:
            self%nbond = nheader - offset_bond
            ! storage:
            allocate( self%bond_code(self%nbond), stat=status )
            IF_NOTOK_RETURN(status=1)
          end if
          ! target substance:
          ibond = ifield - offset_bond
          ! store:
          self%bond_code(ibond) = trim(header)
      end select  ! header
    end do  ! fields

    ! check index:
    if ( ifield_code < 0 ) then
      write (gol,'("substance code index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
      write (gol,'("  ",a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if

    ! check index:
    if ( ifield_name < 0 ) then
      write (gol,'("substance name index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
      write (gol,'("  ",a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage for bonds:
    allocate( self%bonds(self%nsubstance,self%nbond), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! by default assign a warning value ...
    self%bonds = -999.9

    ! loop over records:
    do

      ! increase record counter:
      iline = iline + 1
      ! try to read line:
      read (fu,'(a)',iostat=status) line
      ! eof ?
      if (status<0) exit
      ! error ?
      if (status>0) then
        write (gol,'("reading line ",i6," from file : ",a)') iline, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle

      ! split into records:
      call goSplitString( line, nfield, fields, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( nfield /= nheader ) then
        write (gol,'("number of fields (",i2,") in line ",i2," :")') nfield, iline; call goPr
        write (gol,'("  ",a)') trim(line); call goErr
        write (gol,'("fields:")'); call goErr
        do ifield = 1, nfield
          write (gol,'(i6," : ",a)') ifield, trim(fields(ifield)); call goErr
        end do
        write (gol,'("differs from number of headers ",i2," in file ",a)') nheader, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! get code:
      field = fields(ifield_code)
      ! get index, table might have more entries than used in mass split:
      call goMatchValue( field, substance_codes, isubstance, status, quiet=.true. )
      IF_ERROR_RETURN(status=1)
      ! not needed? then skip:
      if ( isubstance < 0 ) cycle
      ! store:
      self%substance_code(isubstance) = trim(field)
      
      ! get name:
      field = fields(ifield_name)
      ! store:
      self%substance_name(isubstance) = trim(field)

      ! check ...
      if ( any(self%bonds(isubstance,:) >= 0.0) ) then
        write (gol,'("found second record for substance `",a,"`")') &
                        trim(substance_codes(isubstance)); call goErr
        write (gol,'("line    : ",i6)') iline; call goErr
        write (gol,'("file    : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! loop over carbon bond fields:
      do ifield = offset_bond+1, nfield
      
        ! value field:
        field = fields(ifield)
        ! read value:
        read (field,*,iostat=status) value
        if (status/=0) then
          write (gol,'("reading carbon bond ratio from `",a,"`")') field; call goErr
          TRACEBACK; status=1; return
        end if

        ! check for strange values ...
        if ( value < 0.0 ) then
          write (gol,*) 'found strange carbon bond ratio: ', value; call goErr
          write (gol,'("fields  : ",a)') trim(line); call goErr
          write (gol,'("line    : ",i6)') iline; call goErr
          write (gol,'("file    : ",a)') trim(filename); call goErr
          TRACEBACK; status=1; return
        end if

        ! index in bond array:
        ibond = ifield - offset_bond
        ! store
        self%bonds(isubstance,ibond) = value
        
      end do  ! bond fields

    end do  ! lines

    ! close
    close( fu, iostat=status )
    if (status/=0) then
      write (gol,'("closing file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( any(self%bonds < 0.0) ) then
      ! info ..
      write (gol,'("carbon bond ratios not defined for some substance codes:")'); call goErr
      ! loop:
      do isubstance = 1, self%nsubstance
        write (gol,'("  ",i3," ",a4" (",a,")")') isubstance, &
                trim(substance_codes(isubstance)), trim(self%substance_name(isubstance)); call goErr
        do ibond = 1, self%nbond
          write (gol,'("    ",a8," ",f8.2)') trim(self%bond_code(ibond)), self%bonds(isubstance,ibond); call goErr
        end do  ! bonds
      end do  ! substances
      write (gol,'("file: ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine EmisComp_CarbonBonds_Init
  
  
  ! ***
  

  subroutine EmisComp_CarbonBonds_Done( self, status )
  
    ! --- in/out ------------------------------
    
    class(T_EmisComp_CarbonBonds), intent(inout)    ::  self
    integer, intent(out)                            ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/EmisComp_CarbonBonds_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate( self%substance_code, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%substance_name, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%bonds, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine EmisComp_CarbonBonds_Done




end module LE_Emis_Composition
