!#################################################################
!
! NAME
!
!   MAORI_DataFile_CSV  -  read observation data from CSV files.
!
! DATA FILE
!
!   Example of data file:
!
!    ---[201504150000_dust_c.csv]--------------------------
!    compound, concentration, location
!    dust_c,   092,             China001
!    dust_c,   553,             China002
!    dust_c,   397,             China003
!    :
!    :
!    ---------------------------------------------------------------------
!
! USAGE
!
!   use MAORI_DataFile_CSV
!
!   type(T_MAORI_DataFile_CSV)        ::  file
!   character(len=1024)                  ::  keys
!   character(len=32)                    ::  code
!   character(len=32)                    ::  comp
!   real                                 ::  value
!   character(len=32)                    ::  units
!
!   ! Specify header names to be used for needed variables:
!   !    code, comp, value, units
!   ! If no units are defined in the file, enclose value by double quotes:
!   !    units="ug/m3"
!   ! If units are included in the header as "compound(ug/m3)" use:
!   !    units=(UNITS)
!   ! If no component name is in the file, the calling program should handle:
!   !    comp=None
!   ! To filter on selected component only, specify the name to filter on:
!   !    filter_comp=dust_c
!   ! If the file has multiple time records, filter on time and specify headers for columns with time values:
!   !    filter_time=T;year=YEAR;month=MONTH;day=DAY;hour=HOUR
!   ! Optional:
!   !     comment=#
!   !     sep=,          : default is ';'
!   ! Specification for above example:
!   keys = 'code=id;name=longname;longitude=lon;latitude=lat;altitude=alt;sep=,'
!
!   ! open file:
!   call file%Init( '201504150000_dust_c.csv', status )
!   if (status/=0) stop
!
!   ! loop over records:
!   do
!     ! read next record:
!     call self%ReadRecord( code, comp, value, units, status )
!     if (status<0) exit ! eof
!     if (status/=0) stop
!   end do
!  
!   ! close:
!   call self%Done( status )
!   if (status/=0) stop
!  
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!#################################################################

module MAORI_DataFile_CSV

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_DataFile_CSV


  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_DataFile_CSV'
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_DataFile_CSV
    ! file name:
    character(len=512)          ::  fname
    ! seperation character:
    character(len=1)            ::  sep
    character(len=1)            ::  comment
    ! file unit:
    integer                     ::  fu
    ! line number:
    integer                     ::  iline
    ! column indices:
    integer                     ::  icol_code
    integer                     ::  icol_comp
    integer                     ::  icol_value
    integer                     ::  icol_units
    ! filter on time?
    logical                     ::  filter_time
    integer                     ::  icol_year, icol_month, icol_day, icol_hour
    ! filter on component?
    character(len=32)           ::  filter_comp
    ! fixed values:
    character(len=32)           ::  units
    character(len=32)           ::  comp
  contains
    procedure :: Init         => CSV_Init
    procedure :: Done         => CSV_Done
    procedure :: ReadRecord   => CSV_ReadRecord
  end type T_MAORI_DataFile_CSV


contains


  ! ======================================================================


  subroutine CSV_Init( self, fname, keys, status )
  
    use GO, only : goGetFU
    use GO, only : goVarValue, goSplitString, goMatchValue
    
    ! --- in/out --------------------------------
    
    class(T_MAORI_DataFile_CSV), intent(out)      ::  self
    character(len=*), intent(in)                  ::  fname
    character(len=*), intent(in)                  ::  keys
    integer, intent(out)                          ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/CSV_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=4000)   ::  line
    character(len=128)    ::  key
    character(len=128)    ::  headers(20)
    character(len=128)    ::  units(20)
    integer               ::  nheader
    integer               ::  i, k, l
    
    ! --- begin ----------------------------------
    
    ! set special characers:
    self%comment = '#'
      call goVarValue( keys, ';', 'comment', '=', self%comment, status )
      IF_ERROR_RETURN(status=1)
    self%sep     = ';'
      call goVarValue( keys, ';', 'sep', '=', self%sep, status )
      IF_ERROR_RETURN(status=1)
    
    ! store file name:
    self%fname = trim(fname)
    
    ! check ...
    inquire( file=trim(self%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("data file not found:")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( self%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    !! testing ...
    !write (gol,'(a,": open file: ",a)') rname, trim(self%fname); call goPr
    
    ! open file:
    open( self%fu, file=trim(self%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening data file:")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! init line counter:
    self%iline = 0

    ! loop over header lines:
    do
      ! next line ...
      self%iline = self%iline + 1
      ! read line:
      read (self%fu,'(a)',iostat=status) line
      if ( status/=0 ) then
        write (gol,'("reading line from EMEP data file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(self%fname); call goErr
        write (gol,'("  line   : ",i6)') self%iline; call goErr
        TRACEBACK; status=1; return
      end if
      ! skip blank lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == self%comment ) cycle
      ! first non-blank line, quit:
      exit
    end do   ! lines
    
    ! split:
    call goSplitString( line, nheader, headers, status, sep=self%sep )
    IF_NOTOK_RETURN(status=1)
    ! fill units if in headers:  "concentration(ug/m3)"
    do i = 1, nheader
      ! length:
      l = len_trim(headers(i))
      ! ends with ')' ?
      if ( headers(i)(l:l) == ')' ) then
        ! search start:
        k = index( headers(i)(1:l), '(' )
        if ( k < 1 ) then
          write (gol,'("header `",a,"` ends with `)` but start `(` not found ...")') trim(headers(i)); call goErr
          write (gol,'("  file : ",a)') trim(self%fname); call goErr
          TRACEBACK; status=1; return
        end if
        ! extract units:
        units(i) = headers(i)(k+1:l-1)
        ! limit header:
        headers(i) = headers(i)(1:k-1)
      else
        ! no units in header ..
        units(i) = ''
      end if
    end do
    
    ! search column with station code:
    call goVarValue( keys, ';', 'code', '=', key, status, verbose=.true. )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( key, headers(1:nheader), self%icol_code, status )
    IF_NOTOK_RETURN(status=1)

    ! by default no component column:
    key = 'None'
    ! search column with component name:
    call goVarValue( keys, ';', 'comp', '=', key, status, verbose=.false. )
    IF_ERROR_RETURN(status=1)
    if ( trim(key) == 'None' ) then
      ! dumy index:
      self%icol_comp = -999
      ! store value, this will be returned:
      self%comp = trim(key)
    else
      ! name of comp column:
      call goMatchValue( key, headers(1:nheader), self%icol_comp, status )
      IF_NOTOK_RETURN(status=1)
      ! dummy:
      self%comp = '*'
    end if

    ! search column with data values:
    call goVarValue( keys, ';', 'value', '=', key, status, verbose=.true. )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( key, headers(1:nheader), self%icol_value, status )
    IF_NOTOK_RETURN(status=1)

    ! search column with units:
    call goVarValue( keys, ';', 'units', '=', key, status, verbose=.true. )
    IF_NOTOK_RETURN(status=1)
    !~ units not in file, specified in settings:
    if ( key(1:1) == '"' ) then
      ! dummy:
      self%icol_units = -999
      ! fill actual value, strip quotes:
      l = len_trim(key)
      self%units = key(2:l-1)
    !~ units in header as "concentration(units)", extracted above:
    else if ( key(1:1) == '(' ) then
      ! dummy:
      self%icol_units = -999
      ! copy:
      self%units = trim(units(self%icol_value))
      ! check ..
      if ( len_trim(self%units) == 0 ) then
        write (gol,'("units key `",a,"` implies units in header of `",a,"` column, but extracted value is empty ...")') &
                 trim(key), trim(headers(self%icol_value)); call goErr
        TRACEBACK; status=1; return
      end if
    !~ in file column:
    else
      ! search column:
      call goMatchValue( key, headers(1:nheader), self%icol_units, status )
      IF_NOTOK_RETURN(status=1)
      ! dummy ...
      self%units = ''
    end if

    ! by default no filter on time record:
    self%filter_time = .false.
    ! enable filter on time records?
    call goVarValue( keys, ';', 'filter_time', '=', self%filter_time, status, verbose=.false. )
    IF_ERROR_RETURN(status=1)
    ! enabled?
    if ( self%filter_time ) then
      ! search columns with date value:
      call goVarValue( keys, ';', 'year', '=', key, status, verbose=.true. )
      IF_NOTOK_RETURN(status=1)
      call goMatchValue( key, headers(1:nheader), self%icol_year, status )
      IF_NOTOK_RETURN(status=1)
      ! search columns with date value:
      call goVarValue( keys, ';', 'month', '=', key, status, verbose=.true. )
      IF_NOTOK_RETURN(status=1)
      call goMatchValue( key, headers(1:nheader), self%icol_month, status )
      IF_NOTOK_RETURN(status=1)
      ! search columns with date value:
      call goVarValue( keys, ';', 'day', '=', key, status, verbose=.true. )
      IF_NOTOK_RETURN(status=1)
      call goMatchValue( key, headers(1:nheader), self%icol_day, status )
      IF_NOTOK_RETURN(status=1)
      ! search columns with date value:
      call goVarValue( keys, ';', 'hour', '=', key, status, verbose=.true. )
      IF_NOTOK_RETURN(status=1)
      call goMatchValue( key, headers(1:nheader), self%icol_hour, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! by default no filter on components:
    self%filter_comp = ''
    ! enable filter on components:
    call goVarValue( keys, ';', 'filter_comp', '=', self%filter_comp, status, verbose=.false. )
    IF_ERROR_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine CSV_Init
  
  
  ! ***
  
  
  subroutine CSV_Done( self, status )
  
    ! --- in/out --------------------------------
    
    class(T_MAORI_DataFile_CSV), intent(inout)   ::  self
    integer, intent(out)                            ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/CSV_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( self%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing data file:")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine CSV_Done
  

  ! ***
  
  
  subroutine CSV_ReadRecord( self, code, comp, value, units, status, t )

    use GO, only : goSplitString
    use GO, only : TDate, NewDate, operator(/=), wrtgol

    ! --- in/out --------------------------------
    
    class(T_MAORI_DataFile_CSV), intent(inout)      ::  self
    character(len=*), intent(out)                   ::  code
    character(len=*), intent(out)                   ::  comp
    real, intent(out)                               ::  value
    character(len=*), intent(out)                   ::  units
    integer, intent(out)                            ::  status
    type(TDate), intent(in), optional               ::  t

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/CSV_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=4000)   ::  line
    integer                ::  nvalue
    character(len=128)     ::  values(20)
    integer                ::  year, month, day, hour
    type(TDate)            ::  trec
    
    ! --- begin ----------------------------------
    
    ! loop until first non-empty line or end-of-file:
    do
      ! next line:
      self%iline = self%iline + 1
      ! read line:
      read (self%fu,'(a)',iostat=status) line
      if (status<0) then
        status=-1; return   ! eof
      end if
      if ( status/=0 ) then
        write (gol,'("reading line from CSV loclist file:")'); call goErr
        write (gol,'("  file   : ",a)') trim(self%fname); call goErr
        write (gol,'("  line   : ",i6)') self%iline; call goErr
        TRACEBACK; status=1; return
      end if

      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == self%comment ) cycle
      ! skip if empty:
      if ( len_trim(line) == 0 ) cycle

      ! split:
      call goSplitString( line, nvalue, values, status, sep=self%sep )
      IF_NOTOK_RETURN(status=1)

      ! filter on time records?
      if ( self%filter_time ) then
        ! extract time values:
        read (values(self%icol_year ),'(i6)') year
        read (values(self%icol_month),'(i6)') month
        read (values(self%icol_day  ),'(i6)') day
        read (values(self%icol_hour ),'(i6)') hour
        ! fill:
        trec = NewDate( year, month, day, hour, 0, 0 )
        ! check ..
        if ( .not. present(t) ) then
          write (gol,'("time filter on csv records enabled, but no time arguement past")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! skip if not the same:
        if ( trec /= t ) then
          !! testing ..
          !call wrtgol( '  x skip record, time ', trec, ' not equal to ', t ); call goPr
          ! next record:
          cycle
        end if
      end if  ! time filter

      ! extract station code:
      code = trim(values(self%icol_code))

      ! extract component:
      if ( self%icol_comp < 0 ) then
        ! copy from structure (read from settings):
        comp = trim(self%comp)
      else
        ! copy from column:
        comp = trim(values(self%icol_comp))
      end if
      ! filter on component?
      if ( len_trim(self%filter_comp) > 0 ) then
        ! skip if not the same:
        if ( trim(comp) /= trim(self%filter_comp) ) then
          !! testing ..
          !write (gol,'("  x skip record, component `",a,"` not equal to `",a,"`")') trim(comp), trim(self%filter_comp); call goPr
          ! next record:
          cycle
        end if
      end if

      ! read concentration value:
      read (values(self%icol_value),*,iostat=status) value
      if ( status /= 0 ) then
        write (gol,'("could not read value from column ",i0," : ",a)') &
                        self%icol_value, trim(values(self%icol_value)); call goErr
        write (gol,'("  file    : ",a)') trim(self%fname); call goErr
        write (gol,'("  line    : ",i6)') self%iline; call goErr
        write (gol,'("  content : ",a)') trim(line); call goErr
        TRACEBACK; status=1; return
      end if

      ! units:
      if ( self%icol_units < 0 ) then
        ! copy fixed value:
        units = trim(self%units)
      else
        ! read:
        units = trim(values(self%icol_units))
      end if
      
      !! testing ...
      !write (gol,*) '  x value ', value, ' ', trim(units); call goPr
      
      ! record extracted, leave:
      exit
      
    end do ! lines

    ! ok
    status = 0
    
  end subroutine CSV_ReadRecord
  

end module MAORI_DataFile_CSV

