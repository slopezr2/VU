!#################################################################
!
! call ReadRc( rcfile, 'test.flag', l, status [,default=.false.] )
!
! return status :
!   <0  : key not found, value set to default
!    0  : key found and value read without errors
!   >0  : some errors
!
! Search for extended keys:
!
!   call ReadRc( rcfile, 'test', (/'*  ','all','b  '/), flag, status, default=.true. )
!
! will search for (dots are inserted automatically):
!
!     test.*       :  F 
!     test.all     :  F 
!     test.b       :  T 
!
! The last found key overwrites all previous values.
!  
!#################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#################################################################

module GO_Rc

  use GO_Print, only : gol, goErr

  implicit none

  ! --- in/out ---------------------

  private
  
  public  ::  TrcFile

  public  ::  Init, Done
  public  ::  ReadRc


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'GO_Rc'

  ! maximum line length in rc file:
  integer, parameter     ::  buflen = 2000

  ! --- types ---------------------------------
  
  ! key/value pair:
  type T_KeyAndValue
    ! key and value:
    character(len=:), allocatable   ::  key
    character(len=:), allocatable   ::  value
  contains
    procedure ::  Init        =>  KeyAndValue_Init
    procedure ::  Done        =>  KeyAndValue_Done
  end type T_KeyAndValue

  ! file:
  type TrcFile
    ! input file:
    character(len=:), allocatable       ::  fname
    ! key/value pairs:
    integer                             ::  n
    type(T_KeyAndValue), allocatable    ::  element(:)
    !
  contains
    procedure ::  Init        =>  RcFile_Init
    procedure ::  Done        =>  RcFile_Done
    procedure ::  FindIndex   =>  RcFile_FindIndex
    !
    procedure ::                  RcFile_Get_i
    procedure ::                  RcFile_Get_r
    procedure ::                  RcFile_Get_l
    procedure ::                  RcFile_Get_s
    generic   ::  Get         =>  RcFile_Get_i, &
                                  RcFile_Get_r, &
                                  RcFile_Get_l, &
                                  RcFile_Get_s
    !
  end type TrcFile


  ! --- interfaces -------------------------------------

  interface Init
    module procedure RcFile_Init
  end interface

  interface Done
    module procedure RcFile_Done
  end interface

  interface ReadRc
    module procedure RcFile_Get_i
    module procedure RcFile_Get_i1
    module procedure RcFile_Get_r
    module procedure RcFile_Get_l
    module procedure RcFile_Get_s
  end interface


contains


  ! ================================================================
  ! ===
  ! === key and value pairs
  ! ===
  ! ================================================================


  subroutine KeyAndValue_Init( self, key, value, status )

    use GO_String, only : goSplitLine, goTab2Space
    use GO_File  , only : TTextFile, Init, Done, ReadLine

   ! --- in/out ---------------------------

    class(T_KeyAndValue), intent(out)   ::  self
    character(len=*), intent(in)        ::  key
    character(len=*), intent(in)        ::  value
    integer, intent(out)                ::  status
    
    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/KeyAndValue_Init'
    
    ! --- local --------------------------

    ! --- begin ---------------------------
    
    ! storage:
    allocate( character(len=len_trim(key)) :: self%key, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    self%key = trim(key)

    ! storage:
    allocate( character(len=len_trim(value)) :: self%value, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    self%value = trim(value)

    ! ok
    status = 0

  end subroutine KeyAndValue_Init


  ! ***


  subroutine KeyAndValue_Done( self,  status )

    ! --- in/out ---------------------------

    class(T_KeyAndValue), intent(inout)     ::  self
    integer, intent(out)                    ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/KeyAndValue_Done'
    
    ! --- local ---------------------------

    ! --- begin ---------------------------
    
    ! clear:
    deallocate( self%key, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%value, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine KeyAndValue_Done




  ! ================================================================
  ! ===
  ! === RcFile
  ! ===
  ! ================================================================


  subroutine RcFile_Init( self, filename, status )

    use GO_String, only : goSplitLine, goTab2Space
    use GO_File  , only : TTextFile, Init, Done, ReadLine

   ! --- in/out ---------------------------

    class(TrcFile), intent(out)      ::  self
    character(len=*), intent(in)      ::  filename
    integer, intent(out)              ::  status
    
    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Init'
    
    ! --- local --------------------------
    
    logical                 ::  exist
    type(TTextFile)         ::  file
    character(len=buflen)   ::  line
    integer                 ::  i
    integer                 ::  j
    character(len=buflen)   ::  key
    character(len=buflen)   ::  value

    ! --- begin ---------------------------

    ! file not present ?
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("rcfile not found :")'); call goErr
      write (gol,'("  ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( character(len=len_trim(filename)) :: self%fname, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    self%fname = trim(filename)

    ! first pass: count number of key/value lines:
    self%n = 0
    ! open commented text file:
    call Init( file, trim(self%fname), status, status='old', comment='!' )
    IF_NOTOK_RETURN(status=1)
    ! scan all lines 
    do
      ! read next non empty, non comment line:
      call ReadLine( file, line, status )
      if (status<0) exit  ! end of file
      IF_NOTOK_RETURN(status=1)
      ! increase counter:
      self%n = self%n + 1
    end do
    ! close:
    call Done( file, status )
    IF_NOTOK_RETURN(status=1)
    
    ! any values ?
    if ( self%n > 0 ) then
    
      ! storage:
      allocate( self%element(self%n), stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! open commented text file:
      call Init( file, trim(self%fname), status, status='old', comment='!' )
      IF_NOTOK_RETURN(status=1)
      ! loop over lines: 
      do i = 1, self%n
      
        ! read next non empty, non comment line:
        call ReadLine( file, line, status )
        IF_NOTOK_RETURN(status=1)
      
        ! replace tabs:
        call goTab2Space( line )

        ! split at colon:
        call goSplitLine( line, key, ':', value, status )
        IF_NOTOK_RETURN(status=1)
        
        ! check on doubles?
        if ( i > 1 ) then
          ! loop over current elements:
          do j = 1, i-1
            ! match?
            if ( trim(self%element(j)%key) == trim(key) ) then
              write (gol,'("found key `",a,"` at least twice in file:")') trim(key); call goErr
              write (gol,'("  ",a)') trim(self%fname); call goErr
              TRACEBACK; status=1; return
            end if  ! same?
          end do  ! previous elements
        end if ! check on doubles
        
        ! init key/value pair:
        call self%element(i)%Init( key, value, status )
        IF_NOTOK_RETURN(status=1)
        
      end do
      ! close:
      call Done( file, status )
      IF_NOTOK_RETURN(status=1)
      
    end if  ! n > 0
    
    ! ok
    status = 0

  end subroutine RcFile_Init


  ! ***


  subroutine RcFile_Done( self,  status )

    ! --- in/out ---------------------------

    class(TrcFile), intent(inout)    ::  self
    integer, intent(out)              ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Done'
    
    ! --- local ---------------------------
    
    integer     ::  i

    ! --- begin ---------------------------
    
    ! any elements?
    if ( self%n > 0 ) then
      ! loop over elements:
      do i = 1, self%n
        ! done:
        call self%element(i)%Done( status )
        IF_NOTOK_RETURN(status=1)
      end do ! elements
      ! clear:
      deallocate( self%element, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if ! n > 0

    ! clear:
    deallocate( self%fname, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine RcFile_Done


  ! ***


  ! 
  ! Return index of element with specified key,
  ! or negative number if not found.
  !
  
  subroutine RcFile_FindIndex( self, key, ind, status )

    ! --- in/out ---------------------------

    class(TrcFile), intent(in)       ::  self
    character(len=*)                  ::  key
    integer, intent(out)              ::  ind
    integer, intent(out)              ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_FindIndex'
    
    ! --- local ---------------------------
    
    integer     ::  i

    ! --- begin ---------------------------
    
    ! init return value as not found:
    ind = -999
    
    ! any elements?
    if ( self%n > 0 ) then
      ! loop over elements:
      do i = 1, self%n
        ! compare:
        if ( trim(key) == self%element(i)%key ) then
          ! fill return value:
          ind = i
          ! leave:
          exit
        end if
      end do ! elements
    end if ! any elements
    
    ! ok
    status = 0

  end subroutine RcFile_FindIndex
  
  
  ! ***


  subroutine RcFile_Get_i( self, key, i, status, default )

    use GO_Print, only : gol, goErr, goPr

    ! --- in/out ----------------------------

    class(TrcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    integer, intent(out)                        ::  i
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_i'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOTOK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        i = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      read (self%element(ind)%value,*,iostat=status) i
      if ( status /= 0 ) then
        write (gol,'("while reading integer:")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(self%element(ind)%value); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_i

  ! *

  subroutine RcFile_Get_i1( self, key, i, status, default )

    use GO_Print, only : gol, goErr, goPr

    ! --- in/out ----------------------------

    class(TrcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    integer, intent(out)                        ::  i(:)
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_i1'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOTOK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        i = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      read (self%element(ind)%value,*,iostat=status) i
      if ( status /= 0 ) then
        write (gol,'("while reading integer:")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(self%element(ind)%value); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_i1


  ! ***


  subroutine RcFile_Get_r( self, key, r, status, default )

    use GO_Print, only : gol, goErr, goPr

    ! --- in/out ----------------------------

    class(TrcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    real, intent(out)                           ::  r
    integer, intent(out)                        ::  status
    
    real, intent(in), optional                  ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_r'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOTOK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        r = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      read (self%element(ind)%value,*,iostat=status) r
      if ( status /= 0 ) then
        write (gol,'("while reading real :")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(self%element(ind)%value); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_r
  

  ! ***
  
  
  subroutine RcFile_Get_l( self, key, l, status, default )

    use GO_Print, only : gol, goErr, goPr

    ! --- in/out ----------------------------

    class(TrcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    logical, intent(out)                        ::  l
    integer, intent(out)                        ::  status
    
    logical, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_l'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOTOK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with warning:
      if ( present(default) ) then
        l = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      read (self%element(ind)%value,*,iostat=status) l
      if ( status /= 0 ) then
        write (gol,'("while reading logical :")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(self%element(ind)%value); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_l


  ! ***


  subroutine RcFile_Get_s( self, key, s, status, default )
  
    use GO_Print, only : gol, goErr, goPr

    ! --- in/out ----------------------------

    class(TrcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    character(len=*), intent(out)               ::  s
    integer, intent(out)                        ::  status
    
    character(len=*), intent(in), optional      ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_s'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOTOK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        s = trim(default)
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(self%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      s = trim(self%element(ind)%value)
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_s


end module GO_Rc
