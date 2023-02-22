!###############################################################################
!
! LE_Emis_HeightDistr - LOTOS-EUROS emission height distribution
!
!
!  Provides structure to store emission height distribution.
!  Input:
!   - list with supported cat codes: 
!        (/01,02,..,75/)
!  Output:
!   - number of height layers:
!       nlay
!   - layer boundaries in m :
!       heightb(0:nlay)
!   - lookup table with factor as a function of category and layer:
!       fraction(1:ncat,1:nlay)
!
!  Input file is a text file with ';'-seperated fields.
!  Comment lines start with '#' .
!  Columns:
!    snap
!      Category code.
!    category
!      Category name, only used for interpretation of this file.
!    0.0
!      Fraction of emission that should be injected into lowest layer.
!    90.0 170.0 ...
!      Fractions of emission that should be injected between two height levels. 
!      The column header indicates the top level (in meters), 
!      the header of the previous column the bottom level.
!      Thus in this example, column '90' has the fractions for [0m,90m],
!      column '170' the fractions for [90m,170m] .
!
!    # emission height distribution
!    snap; category                  ;    0;   90;  170;  310;  470;  710;  990
!       1; public power stations     ; 0.00; 0.00; 0.00; 0.08; 0.46; 0.29; 0.17
!       2; comm. inst. combustion    ; 0.00; 0.50; 0.50; 0.00; 0.00; 0.00; 0.00
!       3; industrial combustion     ; 0.00; 0.00; 0.04; 0.19; 0.41; 0.30; 0.06
!       4; production processes      ; 0.00; 0.90; 0.10; 0.00; 0.00; 0.00; 0.00
!       5; extraction fossil fuel    ; 0.00; 0.90; 0.10; 0.00; 0.00; 0.00; 0.00
!       6; solvents                  ; 1.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00
!       7; road transport            ; 1.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00
!       8; other mobile              ; 1.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00
!       9; waste                     ; 0.00; 0.10; 0.15; 0.40; 0.35; 0.00; 0.00
!      10; agriculture               ; 1.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00
!
!  Usage:
!
!    type(T_Emis_HeightDistr)   :  hd
!
!    !
!    ! Read table given provided file query.
!    ! With an empty query, all emissions are assigned to the emitted tracer.
!    ! The query has the following fields and defaults:
!    !  - filename        # full path
!    !  - format=csv      # file format (optional)
!    !  - sep=;           # field seperation character (optional)
!    !  - comment=#       # line comment character (optional)
!    ! Provide name of emitted tracer to be able to assign
!    ! fractions to the correct tracers .
!    ! Provide country codes, cat codes, and tracer names.
!    ! Set 'allow_missing_tracers' flag if it is ok that some emissions
!    ! could not be assigned to a tracer (for example if only a limitted
!    ! model is used); if not set, all tracers should exist.
!    !
!    call LE_Emis_HeightDistr_Init( hd, 'file=emisfac.csv', cat_codes, status )
!    if (status/=0) stop
!
!    ! distribute emisisons:
!    em(1:nz) = ... , hd%heightb, hd%fraction(icat,:), ...
!
!    ! done:
!    call LE_Emis_HeightDistr_Done( hd, status )
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

module LE_Emis_EDGAR_HeightDistr

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_HeightDistr
  public  ::  LE_Emis_HeightDistr_Init, LE_Emis_HeightDistr_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_HeightDistr'
  
  
  ! --- types --------------------------------
  
  ! storage for emission composition:
  type T_Emis_HeightDistr
    ! dimensions:
    integer                   ::  ncat
    integer                   ::  nlay
    ! height layer boundaries:
    real, allocatable         ::  heightb(:)
    ! factor to be applied:
    real, allocatable         ::  fraction(:,:)  ! (icat,ilay)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_HeightDistr_Init( hd, query, cat_codes, status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue

    ! --- in/out ------------------------------
    
    type(T_Emis_HeightDistr), intent(out)     ::  hd
    character(len=*), intent(in)              ::  query
    character(len=*), intent(in)              ::  cat_codes(:)      
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_HeightDistr_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 10
    
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
    real                    ::  hhb(0:maxcol)
    integer                 ::  ilays(maxcol)
    integer                 ::  ilay
    real                    ::  fraction
    real                    ::  fsum
    
    ! --- begin -------------------------------
    
    ! extract dimensions:
    hd%ncat     = size(cat_codes)
    
    ! empty ?
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( len_trim(query) == 0 ) then   ! empty
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
      ! single layer:
      hd%nlay = 1
      ! storage:
      allocate( hd%heightb(0:hd%nlay) )
      allocate( hd%fraction(hd%ncat,1:hd%nlay) )
      ! dummy layer of 1 cm ...
      hd%heightb(0) = 0.00 ! m
      hd%heightb(1) = 0.01 ! m
      ! assign everything to this layer:
      hd%fraction(:,:) = 1.0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else  ! read table from file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
          
          ! no layers found yet:
          hd%nlay = 0
          ! set bottom:
          hhb(0) = 0.0   ! m

          ! set variables from header; first dummy values:
          ifield_cat = -1
          ! init mapping from  field index to layer index:
          ilays(:) = -1
          ! loop over fields:
          do ifield = 1, nheader
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
                hd%nlay = hd%nlay + 1
                ! check ...
                if ( hd%nlay > maxcol ) then
                  write (gol,'("number of layers ",i3," exceeds maximum;")') hd%nlay; call goErr
                  write (gol,'("please increase parameter maxcol in module ",a)') trim(mname); call goErr
                  TRACEBACK; status=1; return
                end if
                ! read value:
                read (header,*,iostat=status) hhb(hd%nlay)
                ! assign layer index to column:
                ilays(ifield) = hd%nlay
              !~
            end select  ! header
          end do  ! fields

          ! check index:
          if ( ifield_cat < 0 ) then
            write (gol,'("category index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return
          end if

          ! storage for layer boundaries
          allocate( hd%heightb(0:hd%nlay) )
          ! copy:
          hd%heightb(0:hd%nlay) = hhb(0:hd%nlay) 
          
          ! storage for layer boundaries
          allocate( hd%fraction(1:hd%ncat,1:hd%nlay) )
          ! fill with dummy value to trap undefined elements:
          hd%fraction = -999.9

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
            if ( any(hd%fraction(icat,:) >= 0.0) ) then
              write (gol,'("found second record for category `",a,"`")') &
                              trim(cat_codes(icat)); call goErr
              write (gol,'("line    : ",i6)') iline; call goErr
              write (gol,'("file    : ",a)') trim(filename); call goErr
              TRACEBACK; status=1; return
            end if              

            ! loop over fields:
            do ifield = 1, nfield

              ! current layer:
              ilay = ilays(ifield)
              ! no component ? then skip:
              if ( ilay < 0 ) cycle

              ! current field characters:
              field = fields(ifield)            
              ! read factor:
              read (field,*,iostat=status) fraction
              if (status/=0) then
                write (gol,'("reading fraction from `",a,"`")') field; call goErr
                TRACEBACK; status=1; return
              end if

              ! check for strange values ...
              if ( (fraction < 0.0) .or. (fraction > 1.0) ) then
                write (gol,*) 'found strange layer fraction: ', fraction; call goErr
                write (gol,'("fields  : ",a)') trim(line); call goErr
                write (gol,'("line    : ",i6)') iline; call goErr
                write (gol,'("file    : ",a)') trim(filename); call goErr
                TRACEBACK; status=1; return
              end if

              ! store:
              hd%fraction(icat,ilay) = fraction

            end do  ! fields
            
            ! row sum should be one:
            fsum = sum(hd%fraction(icat,:))
            if ( abs(fsum-1.0) > 1.0e-4 ) then
              write (gol,*) 'sum of fraction is not 1.0 but: ', fsum; call goErr
              write (gol,'("values  : ")'); call goErr
              write (gol,*) hd%fraction(icat,:); call goErr
              write (gol,'("line    : ",i6)') iline; call goErr
              write (gol,'("file    : ",a)') trim(filename); call goErr
              TRACEBACK; status=1; return
            end if

          end do  ! lines
          
          ! close
          close( fu, iostat=status )
          if (status/=0) then
            write (gol,'("closing file : ",a)') trim(filename); call goErr
            TRACEBACK; status=1; return
          end if

          ! check ...
          if ( any(hd%fraction < 0.0) ) then
            write (gol,'("some height distribution factors were not set ...")'); call goErr
            write (gol,'("first undefined record : ")'); call goErr
            do icat = 1, hd%ncat
              do ilay = 1, hd%nlay
                if ( hd%fraction(icat,ilay) < 0.0 ) then
                  write (gol,'("  category index : ",i0)') icat; call goErr
                  write (gol,'("  category code  : `",a,"`")') trim(cat_codes(icat)); call goErr
                  write (gol,'("  layer          : ",i6)') ilay; call goErr
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

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! ok
    status = 0

  end subroutine LE_Emis_HeightDistr_Init
  
  
  ! ***
  

  subroutine LE_Emis_HeightDistr_Done( hd, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_HeightDistr), intent(inout)   ::  hd
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_HeightDistr_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate( hd%fraction )
    deallocate( hd%heightb )

    ! ok
    status = 0

  end subroutine LE_Emis_HeightDistr_Done




end module LE_Emis_EDGAR_HeightDistr
