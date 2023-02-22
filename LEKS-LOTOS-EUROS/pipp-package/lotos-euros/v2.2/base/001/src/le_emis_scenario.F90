!###############################################################################
!
! LE_Emis_Scenario - LOTOS-EUROS emission scenario factor
!
!
!  Provides structure to stare emission scenario factors.
!  Input:
!   - list with supported country codes: 
!        (/'ALB,'ARM',.../)
!   - list with supported cat codes: 
!        (/01,02,..,75/)
!   - list with emitted tracer names:
!        (/'nox','co','pm10',../)
!  Output:
!   - lookup table with factor as a function of category,
!     country, and target component index:
!       factor(1:ncat,1:ncountry,1:ncomp)
!
!  Input file is a text file with ';'-seperated fields.
!  Comment lines start with '#' .
!  The first record contains the column headers.
!  First column has the 3-letter ISO country code, 
!  should match with one of the values in the provided list.
!  Second column has the SNAP code for the category, 
!  should match with one of the values in the provided list.
!  Other columns contain the factors for a specific emitted component,
!  the names should be the same as the provided names.
!
!    #
!    # Scenario factors
!    #
!    #
!    iso3; snap;   ch4;    co;   nh3; nmvoc;   nox;   so2;   pm10; pm2_5
!    ALB ;    1;   1.0;   1.0;   1.0;   1.0;   1.0;   1.0;    1.0;   1.0
!    ALB ;    2;   1.0;   1.0;   1.0;   1.0;   1.0;   1.0;    1.0;   1.0
!    ALB ;    3;   1.0;   1.0;   1.0;   1.0;   1.0;   1.0;    1.0;   1.0
!    :
!
!  Usage:
!
!    type(T_Emis_Scenario)   :  scen
!
!    !
!    ! Read table given provided file query.
!    ! With an empty query, all factor are set to 1.0 .
!    ! The query has the following fields and defaults:
!    !  - filename        # full path
!    !  - format=csv      # file format (optional)
!    !  - sep=;           # field seperation character (optional)
!    !  - comment=#       # line comment character (optional)
!    ! Provide country codes, cat codes, and tracer names.
!    !
!    call LE_Emis_Scenario_Init( scen, 'file=emisfac.csv', &
!                                 country_codes, cat_codes, specnames, status )
!    if (status/=0) stop
!
!    ...
!    ! apply factor:
!    em = em * scen%factor(icat,icountry,icomp)
!    ...
!
!    ! done:
!    call LE_Emis_Scenario_Done( scen, status )
!    if (status/=0) stop
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

module LE_Emis_Scenario

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Scenario
  public  ::  LE_Emis_Scenario_Init, LE_Emis_Scenario_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Scenario'
  
  
  ! --- types --------------------------------
  
  ! storage for emission composition:
  type T_Emis_Scenario
    ! dimensions:
    integer                   ::  ncountry
    integer                   ::  ncat
    integer                   ::  ncomp
    ! factor to be applied:
    real, allocatable         ::  factor(:,:,:)  ! (icat,icountry,icomp)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_Scenario_Init( scen, query, &
                                       cat_codes, country_codes, comp_names, &
                                       status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue
    use LE_Emis_Tools, only : Get_ShortSNAP, ShortSNAP_to_Code

    ! --- in/out ------------------------------
    
    type(T_Emis_Scenario), intent(out)        ::  scen
    character(len=*), intent(in)              ::  query
    character(len=*), intent(in)              ::  cat_codes(:)      ! (/'01.00','02.00',../)
    character(len=*), intent(in)              ::  country_codes(:)  ! (/'ALB',.../)
    character(len=*), intent(in)              ::  comp_names(:)     ! (/'no2','no',../)
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Scenario_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 12
    
    ! --- local -------------------------------
    
    character(len=512)      ::  filename
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    character(len=32)       ::  fileformat
    logical                 ::  exist

    integer                 ::  ifield_country
    integer                 ::  ifield_cat
    integer                 ::  ifield_sectorid
    integer, allocatable    ::  icomps(:)   ! (nfield)

    integer                 ::  fu
    integer                 ::  iline
    character(len=1024)     ::  line
    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    integer                 ::  nfield
    integer                 ::  ifield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    integer                 ::  shortsnap
    character(len=8)        ::  cat_code
    integer                 ::  icountry
    integer                 ::  icat
    integer                 ::  icomp
    real                    ::  factor
    
    ! --- begin -------------------------------
    
    ! extract dimensions:
    scen%ncat     = size(cat_codes)
    scen%ncountry = size(country_codes)
    scen%ncomp    = size(comp_names)
    
    ! storage:
    allocate( scen%factor(scen%ncat,scen%ncountry,scen%ncomp) )
    ! fill with dummy value to be able to check for missing data:
    scen%factor = -999.9
    
    ! empty ?
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( len_trim(query) == 0 ) then   ! empty
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
      ! all unity ...
      scen%factor = 1.0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else  ! read table from file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! extract filename:
      call goVarValue( trim(query), ';', 'file', '=', filename, status )
      if ( status /= 0 ) then
        write (gol,'("could not read key `file=xxx` from emission scenario query : ",a)') trim(query); call goErr
        TRACEBACK; status=1; return
      end if
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

          ! set variables from header; first dummy values:
          ifield_country = -1
          ifield_cat = -1
          ifield_sectorid = -1
          ! storage to map from field index to component index:
          allocate( icomps(nheader) ) ; icomps = -1
          ! loop over fields:
          do ifield = 1, nheader
            ! current:
            header = headers(ifield) 
            ! which column ?
            select case ( trim(header) )
              !~ country code:
              case ( 'ISO3', 'iso3' )
                ! store field index:
                ifield_country = ifield
              !~ category code:
              case ( 'SNAP', 'snap', 'snap1' )
                ! store field index:
                ifield_cat = ifield
              case ( 'SectorID', 'sector_id' )
                ! store field index:
                ifield_sectorid = ifield
              !~ emitted component ...
              case default
                ! search:
                call goMatchValue( trim(header), comp_names, icomp, status )
                IF_NOTOK_RETURN(status=1)
                ! assign to field index:
                icomps(ifield) = icomp
              !~
            end select  ! header
          end do  ! fields

          ! check index:
          if ( ifield_cat < 0 .and. ifield_sectorid < 0 ) then
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
              ! search ...
              call goMatchValue(  fields(ifield_cat), cat_codes, icat, status )
              IF_NOTOK_RETURN(status=1)
            end if
            
            ! check ...
            if ( any(scen%factor(icat,icountry,:) >= 0.0) ) then
              write (gol,'("found second record for country `",a,"` category `",a,"`")') &
                              country_codes(icountry), trim(cat_codes(icat)); call goErr
              write (gol,'("line    : ",i6)') iline; call goErr
              write (gol,'("file    : ",a)') trim(filename); call goErr
              TRACEBACK; status=1; return
            end if              

            ! loop over fields:
            do ifield = 1, nfield

              ! current component:
              icomp = icomps(ifield)
              ! no component ? then skip:
              if ( icomp < 0 ) cycle

              ! current field characters:
              field = fields(ifield)            
              ! read factor:
              read (field,*,iostat=status) factor
              if (status/=0) then
                write (gol,'("reading fraction from `",a,"`")') field; call goErr
                TRACEBACK; status=1; return
              end if

              ! check for strange values ...
              if ( factor < 0.0 ) then
                write (gol,*) 'found strange factor: ', factor; call goErr
                write (gol,'("fields  : ",a)') trim(line); call goErr
                write (gol,'("line    : ",i6)') iline; call goErr
                write (gol,'("file    : ",a)') trim(filename); call goErr
                TRACEBACK; status=1; return
              end if

              ! store:
              scen%factor(icat,icountry,icomp) = factor

            end do  ! fields

          end do  ! lines
          
          ! check ...
          if ( any(scen%factor < 0.0) ) then
            write (gol,'("some scenario factors were not set ...")'); call goErr
            write (gol,'("first undefined record : ")'); call goErr
            do icat = 1, scen%ncat
              do icountry = 1, scen%ncountry
                do icomp = 1, scen%ncomp
                  if ( scen%factor(icat,icountry,icomp) < 0.0 ) then
                    write (gol,'("  category code : `",a,"`")') trim(cat_codes(icat)); call goErr
                    write (gol,'("  country code  : ",a)') trim(country_codes(icountry)); call goErr
                    write (gol,'("  component     : ",a)') trim(comp_names(icomp)); call goErr
                    TRACEBACK; status=1; return
                  end if
                end do
              end do
            end do
          end if

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

  end subroutine LE_Emis_Scenario_Init
  
  
  ! ***
  

  subroutine LE_Emis_Scenario_Done( scen, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Scenario), intent(inout)      ::  scen
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Scenario_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate( scen%factor )

    ! ok
    status = 0

  end subroutine LE_Emis_Scenario_Done




end module LE_Emis_Scenario
