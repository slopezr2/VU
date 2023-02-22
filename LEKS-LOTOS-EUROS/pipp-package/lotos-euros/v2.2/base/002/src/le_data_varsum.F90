!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Data_VarSum

  use GO  , only : gol, goPr, goErr

  ! --- in/out -----------------------------------
  
  private

  public  ::  T_VarSum
  
  public  ::  VarSum_Init
  public  ::  VarSum_Done
  

  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Data_VarSum'


  ! --- types --------------------------------
  
  ! linear combination of variables:
  type T_VarSum
    ! number variables in sum:
    integer                           ::  n
    ! variable names:
    character(len=32), allocatable    ::  name(:)  ! (n)
    ! weights:
    real, allocatable                 ::  w(:)     ! (n)
  contains
    procedure :: Init          => VarSum_Init
    procedure :: Done          => VarSum_Done
    procedure :: Show          => VarSum_Show
    procedure :: Get           => VarSum_Get
    procedure :: GetElement    => VarSum_GetElement
    procedure :: Enable        => VarSum_Enable
  end type T_VarSum
  

contains


  ! ====================================================================
  ! ===
  ! === VarSum
  ! ===
  ! ====================================================================


  ! Line has form of weighted sum:
  !   a + 0.1 * b - c / 3
  ! Analyse this up to the form:
  !   1.0 * a  +  0.1 * b  +  (-1.0/3) * c

  subroutine VarSum_Init( self, input, status, verbose )

    use Binas, only : xm_Na, xm_seasalt
    use GO   , only : goIsNum
    
    ! --- in/out ---------------------------------
    
    class(T_VarSum), intent(out)                  ::  self
    character(len=*), intent(in)                  ::  input
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  verbose
      
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VarSum_Init'
    
    ! assumed maximum number of variables involved:
    integer, parameter    ::  nmax = 10
      
    ! --- local ----------------------------------
    
    logical                ::  verb
    character(len=64)      ::  name
    real                   ::  w
    integer                ::  i, j
    character(len=1)       ::  c, ce
    character(len=64)      ::  element
    character(len=64)      ::  part
    logical                ::  isnum
    real                   ::  factor
    logical                ::  inv
    logical                ::  neg
      
    ! --- begin ----------------------------------

    ! info?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! init storage:
    allocate( self%name(nmax), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%w(nmax), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! no data yet:
    self%n = 0
    
    ! Init factor and name: 
    element = ''
    neg = .false.
    ! loop over characters:
    do i = 1, len_trim(input)
      ! current:
      c = input(i:i)
      ! skip whitespace:
      if ( c == ' ' ) cycle
      ! end of sum element?
      if ( (i == len_trim(input)) .or. &
           ( (len_trim(element) > 0) .and. ((c == '+') .or. (c == '-')) ) ) then
        ! add final character:
        if ( i == len_trim(input) ) element = trim(element)//c
        ! possible forms of element:
        !   vname
        !   0.1 * vname
        ! init target values:
        name = ''
        w = 1.0
        ! init part:
        part = ''
        inv = .false.
        ! loop over characters in element:
        do j = 1, len_trim(element)
          ! current:
          ce = element(j:j)
          ! end of part?
          if ( (j == len_trim(element)) .or. (ce == '*') .or. (ce == '/') ) then
            ! add final character if this is end of element:
            if ( j == len_trim(element) ) part = trim(part)//ce
            ! check for special names:
            select case ( trim(part) )
              !~ labeled constants:
              case ( 'xm_Na' )
                factor = xm_Na
                isnum = .true.
              case ( 'xm_seasalt' )
                factor = xm_seasalt
                isnum = .true.
              !~ might be a numeric value ...
              case default
                ! check characters:
                isnum = goIsNum( part )
                ! numberic value?
                if ( isnum ) then
                  ! read factor:
                  read (part,*,iostat=status) factor
                  if ( status /= 0 ) then
                    write (gol,'("part `",a,"` of element `",a,"` interpreted as numeric, but could not read value")') &
                                   trim(part), trim(element); call goErr
                    TRACEBACK; status=1; return
                  end if
                end if  ! read factor
            end select
            ! numberic value?
            if ( isnum ) then
              ! update weight:
              if ( inv ) then
                w = w / factor
              else
                w = w * factor
              end if
            else
              ! this is the variable name, check if still empty ...
              if ( len_trim(name) > 0 ) then
                write (gol,'("from element `",a,"` extracted name `",a,"` but found name `",a,"` before")') &
                                trim(element), trim(part), trim(name); call goErr
                TRACEBACK; status=1; return
              end if
              ! store:
              name = trim(part)
              ! could not invert variables yet ...
              if ( inv ) then
                write (gol,'("element `",a,"` interpreted as division by name `",a,"` which is not supported yet ..")') &
                               trim(element), trim(name); call goErr
                TRACEBACK; status=1; return
              end if
            end if
            ! reset:
            part = ''
            ! should next factor be inverted?
            inv = ce == '/'
          else
            ! add to part:
            part = trim(part)//ce
          end if
        end do ! j (element chars)
        ! swap sign?
        if ( neg ) w = -1.0 * w
        ! store:
        self%n = self%n + 1
        self%name(self%n) = trim(name)
        self%w   (self%n) = w
        ! reset:
        name = ''
        w    = 1.0
        element = ''
        neg = c == '-'
      else
        ! add to element:
        element = trim(element)//c
      end if
    end do ! i (input chars)
    
    ! info ...
    if ( verb ) then
      write (gol,'("interpreted variable sum description: ",a)') trim(input); call goPr
      do i = 1, self%n
        if ( len_trim(self%name(i)) == 0 ) then
          write (gol,'("  element: ",f16.4)') self%w(i); call goPr
        else
          write (gol,'("  element: ",f16.4," * ",a)') self%w(i), trim(self%name(i)); call goPr
        end if
      end do
    end if

    ! ok
    status = 0
    
  end subroutine VarSum_Init


  ! ***


  subroutine VarSum_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_VarSum), intent(inout)                ::  self
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VarSum_Done'
      
    ! --- begin ----------------------------------
    
    ! any storage?
    if ( self%n > 0 ) then
      ! clear:
      deallocate( self%name, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%w, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine VarSum_Done


  ! ***


  subroutine VarSum_Show( self, prompt, status )
    
    ! --- in/out ---------------------------------
    
    class(T_VarSum), intent(in)                   ::  self
    character(len=*), intent(in)                  ::  prompt
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VarSum_Done'
      
    ! --- local ----------------------------------
      
    integer     ::  i
      
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a,"number of elements: ",i0)') prompt, self%n; call goPr
    ! any?
    if ( self%n > 0 ) then
      ! loop:
      do i = 1, self%n
        write (gol,'(a,"..",i2," ",f16.4," ",a)') prompt, i, self%w(i), trim(self%name(i)); call goPr
      end do ! i
    end if
    
    ! ok
    status = 0
    
  end subroutine VarSum_Show


  ! ***


  subroutine VarSum_Get( self, status, n, molemass )
  
    use Binas, only : xm_H, xm_C, xm_N, xm_O, xm_S
    
    ! --- in/out ---------------------------------
    
    class(T_VarSum), intent(in)                   ::  self
    integer, intent(out)                          ::  status
    integer, intent(out), optional                ::  n
    real, intent(out), optional                   ::  molemass  ! kg/mole

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VarSum_Get'
      
    ! --- local ----------------------------------
      
    integer     ::  i
    real        ::  xm
      
    ! --- begin ----------------------------------
    
    ! return number of elements?
    if ( present(n) ) n = self%n
    
    ! interpreted as molemass ..
    if ( present(molemass) ) then
      ! Formula:
      !     C + H * 4
      !     64.8
      ! init sum:
      molemass = 0.0
      ! loop over contributions:
      do i = 1, self%n
        ! empty atom if molemass is inserted directly:
        if ( len_trim(self%name(i)) == 0 ) then
          xm = 1.0e-3  ! kg/mole
        else
          ! element:
          select case ( trim(self%name(i)) )
            case ( 'H' ) ; xm = xm_H  ! kg/mole
            case ( 'C' ) ; xm = xm_C  ! kg/mole
            case ( 'N' ) ; xm = xm_N  ! kg/mole
            case ( 'O' ) ; xm = xm_O  ! kg/mole
            case ( 'S' ) ; xm = xm_S  ! kg/mole
            case default
              write (gol,'("unsupported element `",a,"` for molemass")') trim(self%name(i)); call goErr
              TRACEBACK; status=1; return
          end select
        end if
        ! add contribution:
        molemass = molemass + self%w(i) * xm  ! kg/mole
      end do ! i
    end if ! molemass        
    
    ! ok
    status = 0
    
  end subroutine VarSum_Get


  ! ***


  subroutine VarSum_GetElement( self, i, status, &
                                  pdata, w, name, units, &
                                  check_units, check_lbo, check_ubo )

    use LE_Data_VData, only : LE_Data_VData_GetPointer
    use LE_Data_VData, only : LE_Data_VData_InqVar
    
    ! --- in/out ---------------------------------
    
    class(T_VarSum), intent(in)                   ::  self
    integer, intent(in)                           ::  i
    integer, intent(out)                          ::  status
    real, pointer, optional                       ::  pdata(:,:,:)
    real, intent(out), optional                   ::  w
    character(len=*), intent(out), optional       ::  name
    character(len=*), intent(out), optional       ::  units
    character(len=*), intent(in), optional        ::  check_units
    integer, intent(in), optional                 ::  check_lbo(3), check_ubo(3)

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VarSum_GetElement'
      
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (i < 1) .or. (i > self%n) ) then
      write (gol,'("index ",i0," out of range 1:",i0)') i, self%n; call goErr
      TRACEBACK; status=1; return
    end if

    ! return pointer to data?
    if ( present(pdata) ) then
      ! get pointer to variable data corresponding to name, perform checks:
      call LE_Data_VData_GetPointer( trim(self%name(i)), pdata, status, &
                                      units=units, check_units=check_units, &
                                     check_lbo=check_lbo, check_ubo=check_ubo )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! return weight?
    if ( present(w) ) w = self%w(i)
    
    ! return name?
    if ( present(name) ) name = trim(self%name(i))
    
    ! ok
    status = 0
    
  end subroutine VarSum_GetElement


  ! ***


  subroutine VarSum_Enable( self, status )

    use LE_Data_VData, only : LE_Data_VData_Enable
    
    ! --- in/out ---------------------------------
    
    class(T_VarSum), intent(in)                   ::  self
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/VarSum_Enable'
      
    ! --- local ----------------------------------
    
    integer      ::  i
      
    ! --- begin ----------------------------------
    
    ! loop over contributions:
    do i = 1, self%n
      ! defined?
      if ( len_trim(self%name(i)) > 0 ) then
        ! enable variable:
        call LE_Data_VData_Enable( trim(self%name(i)), status )
        IF_NOTOK_RETURN(status=1)
      end if
    end do ! contributions
    
    ! ok
    status = 0
    
  end subroutine VarSum_Enable

  

end module LE_Data_VarSum
