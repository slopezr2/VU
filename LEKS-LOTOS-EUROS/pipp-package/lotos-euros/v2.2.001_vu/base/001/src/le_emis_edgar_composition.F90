!###############################################################################
!
! LE_Emis_EDGAR_Composition - LOTOS-EUROS emission routines
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
!    call LE_Emis_EDGAR_Composition_Init( emcomp, &
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
!    call LE_Emis_EDGAR_Composition_Done( emcomp, status )
!    if (status/=0) stop
!      
! COMPOSITION QUERIES
!
!  A query is a ';'-seperated list with 'key=value' settings.
!  Examples:
!     ''
!     'skip=T'
!     'target_tracer=ec_f'
!     'file=PM_composition.csv;format=csv;comment=#;sep=,;default_tracer=ppm_f;allow_missing_tracers=F;nocheck=F'
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
!  If "nocheck" is true, the common check if the sum of all fractions is 1 and if all fractions are between   
!     0 and 1 is not performed (for NMVOC, for which the composition table gives mol tracer / kg NMVOC instead of fraction).
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

module LE_Emis_EDGAR_Composition

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Composition
  public  ::  LE_Emis_Composition_Init, LE_Emis_Composition_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_EDGAR_Composition'
  
  
  ! --- types --------------------------------
  
  ! storage for emission composition:
  type T_Emis_Composition
    ! dimensions:
    integer                   ::  ncountry
    integer                   ::  ncat
    integer                   ::  ncomp
    ! model index of target component:
    integer, allocatable      ::  itracer(:)   ! (icomp)
    ! fraction assigned to component:
    real, allocatable         ::  frac(:,:,:)  ! (icat,icountry,icomp)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_Composition_Init( emcomp, query, emitted_tracer, &
                                       cat_codes, country_codes, tracer_names, &
                                       status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue

    ! --- in/out ------------------------------
    
    type(T_Emis_Composition), intent(out)     ::  emcomp
    character(len=*), intent(in)              ::  query
    character(len=*), intent(in)              ::  emitted_tracer
    character(len=*), intent(in)              ::  cat_codes(:)      ! (/1,2,..,75/)
    character(len=*), intent(in)              ::  country_codes(:)  ! (/'ALB',.../)
    character(len=*), intent(in)              ::  tracer_names(:)     ! (/'no2','no',../)
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_EDGAR_Composition_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 11
    
    ! --- local -------------------------------
    
    logical                 ::  skip
    character(len=32)       ::  target_tracer
    logical                 ::  allow_missing_tracers
    
    character(len=512)      ::  filename
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    character(len=32)       ::  fileformat
    character(len=32)       ::  default_tracer
    logical                 ::  exist
    logical                 ::  nocheck

    integer                 ::  ifield_country
    integer                 ::  ifield_cat
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
 !   character(len=64)       ::  cat_code
    integer                 ::  icat
    character(len=64)       ::  component
    integer                 ::  itracer
    real                    ::  frac
    real                    ::  fsum
    
    ! --- begin -------------------------------
    
    ! store currently known dimensions:
    emcomp%ncountry = size(country_codes)
    emcomp%ncat     = size(cat_codes)
    
    ! skip this emission ?
    skip = .false.
    call goVarValue( trim(query), ';', 'skip', '=', skip, status )
    IF_ERROR_RETURN(status=1)
    
    ! target tracer, by default the same name as the emitted tracer:
    target_tracer = trim(emitted_tracer)
    ! extract from the key if defined in there:
    call goVarValue( trim(query), ';', 'target_tracer', '=', target_tracer, status )
    IF_ERROR_RETURN(status=1)
    
    ! extract filename, by default empty:
    filename = ''
    call goVarValue( trim(query), ';', 'file', '=', filename, status )
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

    ! do not check sums?
    nocheck = .false.
    call goVarValue( trim(query), ';', 'nocheck', '=', nocheck, status )
    IF_ERROR_RETURN(status=1)

    ! how to fill ?
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( skip ) then   ! do not use ...
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      ! no emitted components:
      emcomp%ncomp = 0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if ( len_trim(filename) == 0 ) then   ! no file, fill with identity
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
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
      emcomp%ncomp = 1
      ! index of emitted component:
      icomp_emitted = 1    

      ! storage for tracer indices:
      allocate( emcomp%itracer(emcomp%ncomp) )
      ! assign to emitted component:
      emcomp%itracer = -1
      emcomp%itracer(icomp_emitted) = itracer_emitted

      ! storage for fractions:
      allocate( emcomp%frac(emcomp%ncat,emcomp%ncountry,emcomp%ncomp) )
      ! by default assign everything to emitted component:
      emcomp%frac = 0.0
      emcomp%frac(:,:,icomp_emitted) = 1.0
    
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
          allocate( emcomp%itracer(nheader) )
          ! not assigned yet:
          emcomp%itracer = -1

          ! set variables from header; first dummy values:
          ifield_country = -1
          ifield_cat = -1
          emcomp%ncomp = 0
          ! storage to map from field index to component index:
          allocate( icomps(nheader) ) ; icomps = -1
          ! loop over fields:
          do ifield = 1, nheader
            ! current:
            header = headers(ifield) 
            ! which column ?
            select case ( trim(header) )
              !~ country code:
              case ( 'Zone' )
                ! store field index:
                ifield_country = ifield
              !~ category code:
              case ( 'IPCC_emis_cat', 'MEIC_emis_cat' )
                ! store field index:
                ifield_cat = ifield
              case ( 'SNAP' )
                ! just skip
              !~ tracer ...
              case default
                ! increase counter:
                emcomp%ncomp = emcomp%ncomp + 1
                ! assign to field index:
                icomps(ifield) = emcomp%ncomp
                ! set name of target component; depends on emitted tracer:
                select case ( trim(emitted_tracer) )
                  !~ NOx emissions
                  case ( 'nox', 'NOx' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_NO'  ) ; component = 'no'
                      case ( 'Fraction_NO2' ) ; component = 'no2'
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ SOx emissions
                  case ( 'sox', 'so2', 'SO2' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_SO2'  ) ; component = 'so2'
#ifdef with_m7
                      case ( 'Fraction_so4ks' ) ; component = 'so4ks'
                      case ( 'Fraction_so4as' ) ; component = 'so4as'
#else
                      case ( 'Fraction_SO4a' ) ; component = 'so4a_f'
#endif

                      case default
                        write (gol,'("SOX unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  
                  case ( 'oc', 'OC' )
                    ! switch per emitted component:
                    ! write(gol, '("Jianbing test: print the emission header ", a)' ) trim(header); call GoPr
                    select case ( trim(header) )
#ifdef with_m7
                      case ( 'Fraction_ocks' ) ; component = 'ocks'
                      case ( 'Fraction_ocki' ) ; component = 'ocki'
#endif
                      case default
                        write (gol,'("OC unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ fine mode aerosol
                  case ( 'ppm_f', 'pm2_5', 'PM2.5_bio', 'PM2.5_fossil' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_EC'                       ) ; component = 'ec_f'
                      case ( 'Fraction_OC (Full Molecular Mass)' ) ; component = 'pom_f'
                      case ( 'Fraction_SO4'                      ) ; component = 'so4a_f'
                      case ( 'Fraction_Na'                       ) ; component = 'na_f'
                      case ( 'Fraction_Other mineral'            ) ; component = 'ppm_f'
                      case ( 'Fraction_PPM'                      ) ; component = 'ppm_f'
#ifdef with_m7
                      case ( 'Fraction_duas' ) ; component = 'duas'
                      case ( 'Fraction_ducs' ) ; component = 'ducs'
                      case ( 'Fraction_duai' ) ; component = 'duai'
                      case ( 'Fraction_duci' ) ; component = 'duci'
                      case ( 'Fraction_ocks' ) ; component = 'ocks'
                      case ( 'Fraction_ocas' ) ; component = 'ocas'
                      case ( 'Fraction_occs' ) ; component = 'occs'
                      case ( 'Fraction_ocki' ) ; component = 'ocki'
                      case ( 'Fraction_bcas' ) ; component = 'bcas'
                      case ( 'Fraction_bccs' ) ; component = 'bccs'
#endif                    
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ coarse mode aerosol
                  case ( 'ppm_c', 'pm25_pm10', 'PM25_PM10')
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_EC'                       ) ; component = 'ec_c'
                      case ( 'Fraction_OC (Full Molecular Mass)' ) ; component = 'pom_c'
                      case ( 'Fraction_SO4'                      ) ; component = 'so4a_c'
                      case ( 'Fraction_Na'                       ) ; component = 'na_c'
                      case ( 'Fraction_Other mineral'            ) ; component = 'ppm_c'
                      case ( 'Fraction_PPM'                      ) ; component = 'ppm_c'
#ifdef with_m7
                      case ( 'Fraction_ducs'                     ) ; component = 'ducs'
                      case ( 'Fraction_duci'                     ) ; component = 'duci'
                      case ( 'Fraction_occs'                     ) ; component = 'occs'
                      case ( 'Fraction_ocki'                     ) ; component = 'ocki' ! this is wrong, temporal assumpation 

#endif
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ PM10
                  case ( 'PM10' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'Fraction_EC'                       ) ; component = 'ec_f'
                      case ( 'Fraction_PPMf'                     ) ; component = 'ppm_f'
                      case ( 'Fraction_PPMc'                     ) ; component = 'ppm_c'
#ifdef with_m7
                      case ( 'Fraction_ducs'                     ) ; component = 'ducs'
                      case ( 'Fraction_duci'                     ) ; component = 'duci'
#endif                        
                      case default
                        write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
                        TRACEBACK; status=1; return
                    end select
                  !~ unknown ...
                  !~ NMVOC
                  case ( 'NMVOC' )
                    ! switch per emitted component:
                    select case ( trim(header) )
                      case ( 'ETH'                               ) ; component = 'eth'
                      case ( 'OLE'                               ) ; component = 'ole'
                      case ( 'PAR'                               ) ; component = 'par'
                      case ( 'ALD'                               ) ; component = 'ald'
                      case ( 'FORM'                              ) ; component = 'form'
                      case ( 'XYL'                               ) ; component = 'xyl'
                      case ( 'TOL'                               ) ; component = 'tol'
                      case ( 'UNR'                               )  
                        ! just skip ...
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
                emcomp%itracer(emcomp%ncomp) = itracer
              !~
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
            emcomp%ncomp = emcomp%ncomp + 1
            ! set index for default tracer:
            icomp_default = emcomp%ncomp
            ! store model tracer index:
            emcomp%itracer(icomp_default) = itracer_default
          end if

          ! storage for fractions:
          allocate( emcomp%frac(emcomp%ncat,emcomp%ncountry,emcomp%ncomp) )
          ! by default assign a warning value ...
          emcomp%frac = -999.9

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

            ! search category code ...
            call goMatchValue( fields(ifield_cat), cat_codes, icat, status )
            IF_NOTOK_RETURN(status=1)
            
            ! check ...
            if ( any(emcomp%frac(icat,icountry,:) >= 0.0) ) then
              write (gol,'("found second record for country `",a,"` category ",i4)') &
                              country_codes(icountry), cat_codes(icat); call goErr
              write (gol,'("line    : ",i6)') iline; call goErr
              write (gol,'("file    : ",a)') trim(filename); call goErr
              TRACEBACK; status=1; return
            end if              

            ! default tracer specified ?
            if ( itracer_default > 0 ) then
              ! initialize values, assign everything to emitted if possible:
              emcomp%frac(icat,icountry,:            ) = 0.0
              emcomp%frac(icat,icountry,icomp_default) = 1.0
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

              ! check for strange values ?
              if (.not. nocheck ) then
                if ( (frac < 0.0) .or. (frac > 1.0) ) then
                  write (gol,*) 'found strange fraction: ', frac; call goErr
                  write (gol,'("fields  : ",a)') trim(line); call goErr
                  write (gol,'("line    : ",i6)') iline; call goErr
                  write (gol,'("file    : ",a)') trim(filename); call goErr
                  TRACEBACK; status=1; return
                end if
              end if

              ! store:
              emcomp%frac(icat,icountry,icomps(ifield)) = frac

              ! default tracer specified ?
              if ( itracer_default > 0 ) then
                ! remove contribution of current component:
                emcomp%frac(icat,icountry,icomp_default) = &
                      emcomp%frac(icat,icountry,icomp_default) - frac
              end if

            end do  ! fields


            ! check for sum equal 1 ?
            if (.not. nocheck ) then
              ! current sum:
              fsum = sum(emcomp%frac(icat,icountry,:))
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
                do ifield = 1, emcomp%ncomp
                  write (gol,*) ifield, ' : ', emcomp%frac(icat,icountry,ifield); call goErr
                end do
                write (gol,'("  sum      : ")'); call goErr
                write (gol,*) fsum; call goErr
                TRACEBACK; status=1; return
              end if
              ! scale all emisisons:
              emcomp%frac(icat,icountry,:) = emcomp%frac(icat,icountry,:) / fsum
            end if

          end do  ! lines
          
          !! check ...
          !if ( any(emcomp%frac < 0.0) ) then
          !  write (gol,'("some values were not set in fraction table ...")'); call goErr
          !  write (gol,'("first undefined record : ")'); call goErr
          !  do icat = 1, emcomp%ncat
          !    do icountry = 1, emcomp%ncountry
          !      if ( any(emcomp%frac(icat,icountry,:) < 0.0) ) then
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

  end subroutine LE_Emis_Composition_Init
  
  
  ! ***
  

  subroutine LE_Emis_Composition_Done( emcomp, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Composition), intent(inout)   ::  emcomp
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_EDGAR_Composition_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    if ( allocated(emcomp%itracer) ) deallocate( emcomp%itracer )
    if ( allocated(emcomp%frac   ) ) deallocate( emcomp%frac    )

    ! ok
    status = 0

  end subroutine LE_Emis_Composition_Done




end module LE_Emis_EDGAR_Composition
