!###############################################################################
!
! NAME
!   LE_Emis_OPS_EC_Fraction - Fraction of EC in PM25 for OPS emissions
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

module LE_Emis_OPS_EC_Fractions

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile


  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  LE_Emis_OPS_EC_Fractions_Init
  public  ::  LE_Emis_OPS_EC_Fractions_Done

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_OPS_EC_Fractions'
  
contains  
  
  ! ========================
  
  subroutine LE_Emis_OPS_EC_Fractions_Init( query, country_ISO3, ncountry, cat_codes, ncat, ec_frac, status )
  
    use GO, only : GoMatchValue, goVarValue, goGetFu
    use GO, only : goSplitString
    ! --- in/out ---
    
    character(len=*), intent(in)              ::  query
    character(len=3), intent(in)              ::  country_ISO3(1:ncountry)
    integer, intent(in)                       ::  ncountry
    character(len=*), intent(in)              ::  cat_codes(1:ncat)
    integer, intent(in)                       ::  ncat
    real, intent(inout)                       ::  ec_frac(ncountry,ncat)
    integer, intent(out)                      ::  status
    
    ! --- const ---
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_EC_Fractions_Init'
    integer, parameter  ::  max_col = 5
        
    ! --- local ---
    
    logical             ::  exist
    character(len=512)  ::  filename
    integer             ::  fu
    character(len=256)  ::  line
    integer             ::  iline
    character(len=1)    ::  comment, sep
    
    character(len=100)  ::  headers(max_col)
    character(len=100)  ::  header
    integer             ::  nheader
    character(len=100)  ::  fields(max_col)
    integer             ::  ifield, nfield
    integer             ::  ifield_cnt, ifield_cat, ifield_frac
    character(len=3)    ::  country_code
    character(len=8)    ::  cat_code
    integer             ::  icountry, icat
        
    
    ! --- begin ---
    
    ! extract filename, by default empty:
    filename = ''
    call goVarValue( trim(query), ';', 'file', '=', filename, status )
    IF_ERROR_RETURN(status=1)
    ! seperation character:
    sep = ';'
    call goVarValue( trim(query), ';', 'comment', '=', sep, status )
    IF_ERROR_RETURN(status=1)
    ! comment character:
    comment = '#'
    call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    IF_ERROR_RETURN(status=1)

    ! file should be present:
    inquire( file=trim(filename), exist=exist  )
    if ( .not. exist ) then
      write (gol, '("file not found : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( fu, file=trim(filename), status='old', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! line counter
    iline = 0
    
    ! read headerline
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
    
    ! init
    ifield_cnt  = -1
    ifield_cat  = -1
    ifield_frac = -1
    do ifield = 1, nheader
      header = headers(ifield)
      
      select case ( trim(header) )
        case ( 'ISO3' )
          ifield_cnt = ifield
        case ( 'SNAP' )
          ifield_cat = ifield 
        case ( 'EC/PM2.5' )
          ifield_frac = ifield
        case default
          write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
          TRACEBACK; status=1; return
        end select
    end do


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
      call goMatchValue( fields(ifield_cnt), country_ISO3, icountry, status )
      if ( status /= 0 ) then
        ! country not found in PM_composition
        ! dont fill EC-fractions
        cycle
      endif
      
      ! get categroy code
      read ( fields(ifield_cat), * ) cat_code
      call goMatchValue( cat_code, cat_codes, icat, status )
      IF_NOTOK_RETURN(status=1)
         
      ! fill in ec fraction
      read ( fields(ifield_frac), * ) ec_frac(icountry,icat)            
      
    end do          
      
    ! ok
    status = 0
    
  end subroutine LE_Emis_OPS_EC_Fractions_Init
  
  subroutine LE_Emis_OPS_EC_Fractions_Done( status )

    ! --- in/out ---
    integer, intent(out)                      ::  status
    
    ! --- const ---
    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_EC_Fractions_Done'
    
    ! --- local ---
    
    ! --- begin ---
    
    ! ok
    status = 0
  
  end subroutine LE_Emis_OPS_EC_Fractions_Done
  
  
end module LE_Emis_OPS_EC_Fractions
