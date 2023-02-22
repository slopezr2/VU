!###############################################################################
!
! LE_Emis_TempVar - LOTOS-EUROS emission temperature variability for VOC and CO
!
!
!  Provides structure to store temperature-dependent emission factors.
!  Input:
!   - list with supported cat codes: 
!        (/01,02,..,75/)
!  Output:
!   - number of temp categories:
!       ntemp
!   - lookup table with factor as a function of category and layer:
!       fraction(1:ncat,1:ntemp)
!
!  Input file is a text file with ';'-seperated fields.
!  Comment lines start with '#' .
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_EDGAR_TempVar

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_TempVar
  public  ::  LE_Emis_TempVar_Init, LE_Emis_TempVar_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_TempVar'
  
  
  ! --- types --------------------------------
  
  ! storage for emission composition:
  type T_Emis_TempVar
    ! dimensions:
    integer                   ::  ncat
    integer                   ::  ntemp
    ! temperature values:
    real, allocatable         ::  temper(:)
    ! factor to be applied:
    real, allocatable         ::  fraction(:,:)  ! (icat,itemp)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_TempVar_Init( temp_var, query, cat_codes, status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue

    ! --- in/out ------------------------------
    
    type(T_Emis_TempVar), intent(out)         ::  temp_var
    character(len=*), intent(in)              ::  query
    character(len=*), intent(in)              ::  cat_codes(:)      
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_TempVar_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 50
    
    ! --- local -------------------------------
    
    character(len=512)      ::  filename
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    character(len=32)       ::  fileformat
    logical                 ::  exist

    integer                 ::  ifield_cat

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
    integer                 ::  icat
    integer                 ::  itemps(maxcol)
    integer                 ::  itemp
    real                    ::  fraction
    real                    ::  fsum
    real                    ::  hhb(0:maxcol)
    
    ! --- begin -------------------------------
    
    ! extract dimensions:
    temp_var%ncat     = size(cat_codes)
    ! extract filename:
    filename = 'no-file-defined'
    call goVarValue( trim(query), ';', 'file', '=', filename, status )
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
          ! increase line counter:
          iline = iline + 1
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

        ! no temps found yet:
        temp_var%ntemp = 0

        ! set variables from header; first dummy values:
        ifield_cat = -1
        ! init mapping from  field index to layer index:
        itemps(:) = -1
        ! loop over fields:
        do ifield = 1, nheader
          ! testing ...
          ! current:
          header = headers(ifield) 
          ! which column ?
          select case ( trim(header) )
            !~ category code:
            case ( 'IPCC_emis_cat', 'MEIC_emis_cat' )
            ! store field index:
              ifield_cat = ifield
            case ( 'SNAP', 'snap' )
               ! just skip ...
            !~ category name:
            case ( 'category' )
              ! just skip ...
            !~ layer top
            case default
              ! increase counter:
              temp_var%ntemp = temp_var%ntemp + 1
              ! check ...
              if ( temp_var%ntemp > maxcol ) then
                write (gol,'("number of temperatures ",i3," exceeds maximum;")') temp_var%ntemp; call goErr
                write (gol,'("please increase parameter maxcol in module ",a)') trim(mname); call goErr
                TRACEBACK; status=1; return
              end if
              ! read value:
              read (header,*,iostat=status) hhb(temp_var%ntemp)
              ! assign layer index to column:
              itemps(ifield) = temp_var%ntemp
            !~
          end select  ! header
        end do  ! fields

        ! storage:
        allocate( temp_var%fraction(temp_var%ncat,temp_var%ntemp), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init with dummy:
        temp_var%fraction=-999.0
        ! check index:
        if ( ifield_cat < 0 ) then
          write (gol,'("category index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
          write (gol,'("  ",a)') trim(line); call goErr
          TRACEBACK; status=1; return
        end if

        ! storage for layer boundaries
        allocate( temp_var%temper(temp_var%ntemp), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! copy:
        temp_var%temper = hhb(1:temp_var%ntemp) 

        ! loop over records:
        do

          ! increase line counter:
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

          ! search ...
          call goMatchValue( fields(ifield_cat), cat_codes, icat, status )
          IF_NOTOK_RETURN(status=1)

          ! check ...
          if ( any(temp_var%fraction(icat,:) >= 0.0) ) then
            write (gol,'("found second record for category ",a)') &
                            cat_codes(icat); call goErr
            write (gol,'("line    : ",i6)') iline; call goErr
            write (gol,'("file    : ",a)') trim(filename); call goErr
            TRACEBACK; status=1; return
          end if              

          ! loop over fields:
          do ifield = 1, nfield

            ! current layer:
            itemp = itemps(ifield)
            ! no component ? then skip:
            if ( itemp < 0 ) cycle

            ! current field characters:
            field = fields(ifield)            
            ! read factor:
            read (field,*,iostat=status) fraction
            if (status/=0) then
              write (gol,'("reading fraction from `",a,"`")') field; call goErr
              TRACEBACK; status=1; return
            end if

            ! store:
            temp_var%fraction(icat,itemp) = fraction

          end do  ! fields

        end do  ! lines

        ! close
        close( fu, iostat=status )
        if (status/=0) then
          write (gol,'("closing file : ",a)') trim(filename); call goErr
          TRACEBACK; status=1; return
        end if

        ! check ...
        if ( any(temp_var%fraction < 0.0) ) then
          write (gol,'("some height distribution factors were not set ...")'); call goErr
          write (gol,'("first undefined record : ")'); call goErr
          do icat = 1, temp_var%ncat
            do itemp = 1, temp_var%ntemp
              if ( temp_var%fraction(icat,itemp) < 0.0 ) then
                write (gol,'("  category code : ",a)') cat_codes(icat); call goErr
                write (gol,'("  layer         : ",i6)') itemp; call goErr
                TRACEBACK; status=1; return
              end if
            end do
          end do
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        write (gol,'("unsupported format : ",a)') trim(fileformat); call goErr
        TRACEBACK; status=1; return

    end select

    ! ok
    status = 0

  end subroutine LE_Emis_TempVar_Init
  
  
  ! ***
  

  subroutine LE_Emis_TempVar_Done( temp_var, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_TempVar), intent(inout)   ::  temp_var
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_TempVar_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate( temp_var%fraction )

    ! ok
    status = 0

  end subroutine LE_Emis_TempVar_Done




end module LE_Emis_EDGAR_TempVar
