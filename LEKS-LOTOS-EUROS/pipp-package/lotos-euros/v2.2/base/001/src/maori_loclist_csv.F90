!#######################################################################
!
! NAME
!
!   MAORI_LocList_CSV  -  read list of station locations and meta data
!
! DATA FILE
!
!   Example of csv file:
!
!    id,      lon, lat, alt, longname
!    STN001, 5.23,52.4, 2.0, Station upon Countryside
!     :
!
! USAGE
!
!   use MAORI_LocList_CSV_file
!
!   integer                         ::  status
!   character(len=128)              ::  keys
!   type(T_MAORI_LocList_CSV)       ::  file
!   integer                         ::  nrec
!   integer                         ::  irec
!   character(len=16)               ::  code
!   character(len=128)              ::  name
!   character(len=128)              ::  meta
!   real                            ::  lon, lat, alt
!
!   ! Specify header names to be used for needed variables:
!   !    code, name, longitude, latitude, altitude
!   ! If no station height provided, use:
!   !    altitude=None
!   ! Optional:
!   !     comment=#
!   !     sep=,          : default is ';'
!   ! Specification for above example:
!   keys = 'code=id;name=longname;longitude=lon;latitude=lat;altitude=alt;sep=,'
!
!   ! open file, count records:
!   call file%Init( 'loclist.csv', keys, nrec, status )
!   if (status/=0) stop
!
!   ! loop over records:
!   do irec = 1, nrec
!     ! read location, fill meta data with station name etc:
!     call file%ReadRecord( name, code, lon, lat, alt, meta, status )
!     if (status<0) exit   ! eof
!     if (status/=0) stop
!   end do
!
!   call file%Done( status )
!   if (status/=0) stop
!  
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_LocList_CSV

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_MAORI_LocList_CSV
 

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'MAORI_LocList_CSV'
  
  ! value for no data ...
  real, parameter   ::  CSV_nodata = -999.9
  
  
  ! --- types ----------------------------------
   
  type T_MAORI_LocList_CSV
    ! file name:
    character(len=512)          ::  fname
    ! comment character:
    character(len=1)            ::  comment
    ! field seperation character:
    character(len=1)            ::  sep
    ! file unit:
    integer                     ::  fu
    ! line number:
    integer                     ::  iline
    ! column indices:
    integer                     ::  icol_name
    integer                     ::  icol_code
    integer                     ::  icol_lon
    integer                     ::  icol_lat
    integer                     ::  icol_alt
  contains
    procedure :: Init         => MAORI_LocList_CSV_Init
    procedure :: Done         => MAORI_LocList_CSV_Done
    procedure :: ReadRecord   => MAORI_LocList_CSV_ReadRecord
  end type T_MAORI_LocList_CSV
  


contains


  ! ======================================================================


  subroutine MAORI_LocList_CSV_Init( self, fname, keys, nrec, status )
  
    use GO, only : goGetFU
    use GO, only : goVarValue, goSplitString, goMatchValue
    
    ! --- in/out --------------------------------
    
    class(T_MAORI_LocList_CSV), intent(out)   ::  self
    character(len=*), intent(in)              ::  fname
    character(len=*), intent(in)              ::  keys
    integer, intent(out)                      ::  nrec
    integer, intent(out)                      ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/CSV_Init'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    character(len=1024)   ::  line
    character(len=128)    ::  key
    character(len=128)    ::  headers(10)
    integer               ::  nheader
    
    ! --- begin ----------------------------------
    
    ! store file name:
    self%fname = trim(fname)
    
    ! set special characers:
    self%comment = '#'
      call goVarValue( keys, ';', 'comment', '=', self%comment, status )
      IF_ERROR_RETURN(status=1)
    self%sep     = ';'
      call goVarValue( keys, ';', 'sep', '=', self%sep, status )
      IF_ERROR_RETURN(status=1)
    
    ! check ...
    inquire( file=trim(self%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("CSV loclist file not found:")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! free file unit:
    call goGetFU( self%fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open file:
    open( self%fu, file=trim(self%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening CSV loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    self%iline = 0
    
    ! record counter:
    nrec = 0
    do
      ! next line:
      self%iline = self%iline + 1
      ! read line:
      read (self%fu,'(a)',iostat=status) line
      if (status<0) exit   ! eof
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
      ! increase counter:
      nrec = nrec + 1
    end do
    ! do not count header line:
    nrec = nrec - 1
    
    ! close file:
    close( self%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing CSV loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! re-open file:
    open( self%fu, file=trim(self%fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening CSV loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if
    ! reset counter:
    self%iline = 0
    ! loop until header line is read:
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
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do
    
    ! split:
    call goSplitString( line, nheader, headers, status, sep=self%sep )
    IF_NOTOK_RETURN(status=1)
    ! search column:
    call goVarValue( keys, ';', 'code', '=', key, status, verbose=.true. )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( key, headers(1:nheader), self%icol_code, status )
    IF_NOTOK_RETURN(status=1)
    ! search column:
    call goVarValue( keys, ';', 'name', '=', key, status, verbose=.true. )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( key, headers(1:nheader), self%icol_name, status )
    IF_NOTOK_RETURN(status=1)
    ! search column:
    call goVarValue( keys, ';', 'longitude', '=', key, status, verbose=.true. )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( key, headers(1:nheader), self%icol_lon, status )
    IF_NOTOK_RETURN(status=1)
    ! search column:
    call goVarValue( keys, ';', 'latitude', '=', key, status, verbose=.true. )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( key, headers(1:nheader), self%icol_lat, status )
    IF_NOTOK_RETURN(status=1)
    ! search column:
    call goVarValue( keys, ';', 'altitude', '=', key, status, verbose=.true. )
    IF_NOTOK_RETURN(status=1)
    if ( trim(key) == 'None' ) then
      self%icol_alt = -999
    else
      call goMatchValue( key, headers(1:nheader), self%icol_alt, status )
      IF_NOTOK_RETURN(status=1)
    end if


    ! ok
    status = 0
    
  end subroutine MAORI_LocList_CSV_Init
  
  
  ! ***
  
  
  subroutine MAORI_LocList_CSV_Done( self, status )
  
    ! --- in/out --------------------------------
    
    class(T_MAORI_LocList_CSV), intent(inout)   ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/CSV_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    close( self%fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing CSV loclist file:")'); call goErr
      write (gol,'("  ",a)') trim(self%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine MAORI_LocList_CSV_Done
  

  ! ***
  
  
  subroutine MAORI_LocList_CSV_ReadRecord( self, name, code, lon, lat, alt, meta, status )
  
    use GO, only : goSplitString
    
    ! --- in/out --------------------------------
    
    class(T_MAORI_LocList_CSV), intent(inout)   ::  self
    character(len=*), intent(out)               ::  name
    character(len=*), intent(out)               ::  code
    real, intent(out)                           ::  lon, lat   ! [degree]
    real, intent(out)                           ::  alt        ! altitude above sea-level [m]
    character(len=*), intent(out)               ::  meta
    integer, intent(out)                        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/CSV_ReadRecord'
    
    ! --- local ----------------------------------
    
    character(len=1024)    ::  line
    integer                ::  nvalue
    character(len=128)     ::  values(10)
    
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
      ! not empty ? then leave:
      if ( len_trim(line) > 0 ) exit
    end do
    
    ! split:
    call goSplitString( line, nvalue, values, status, sep=self%sep )
    IF_NOTOK_RETURN(status=1)
    ! copy:
    name = trim(values(self%icol_name))
    code = trim(values(self%icol_code))
    ! read:
    read (values(self%icol_lon),*,iostat=status) lon
    if ( status /= 0 ) then
      write (gol,'("could not read lon value from column ",i0," : ",a)') &
                      self%icol_lon, trim(values(self%icol_lon)); call goErr
      write (gol,'("  file   : ",a)') trim(self%fname); call goErr
      write (gol,'("  line   : ",i6)') self%iline; call goErr
    end if
    ! read:
    read (values(self%icol_lat),*) lat
    if ( status /= 0 ) then
      write (gol,'("could not read lat value from column ",i0," : ",a)') &
                      self%icol_lat, trim(values(self%icol_lat)); call goErr
      write (gol,'("  file   : ",a)') trim(self%fname); call goErr
      write (gol,'("  line   : ",i6)') self%iline; call goErr
    end if
    ! altitude:
    if ( self%icol_alt < 0 ) then
      alt = 0.0
    else
      read (values(self%icol_alt),*) alt
    end if
    
    ! no meta data:
    meta = ''
    
    ! ok
    status = 0
    
  end subroutine MAORI_LocList_CSV_ReadRecord
  

end module MAORI_LocList_CSV

