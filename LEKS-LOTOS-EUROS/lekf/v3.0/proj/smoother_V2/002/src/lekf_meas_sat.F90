!#################################################################
!
! NAME
!   LEKF_Meas_Sat  -  interface to satellite data
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
!###############################################################################

module LEKF_Meas_Sat

  use GO,  only : gol, goPr, goErr
  use Num, only : T_Random

  implicit none


  ! --- in/out ----------------------------

  private

  public  ::  T_LEKF_Meas_Sat
  public  ::  meas_sat

  public  ::  T_LEKF_Meas_Sat_Set


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LEKF_Meas_Sat'


  ! --- types --------------------------------

  !
  ! analysis info per sat set:
  !
  type T_LEKF_Meas_Sat_Set
    !
    ! name for this set, should match one of 'leo%file_sat(:)' names:
    character(len=32)     ::  name
    !
    ! analyse this set ?
    logical               ::  analyse
    !
    ! the radius of influence for analysis;
    ! "wild guess" to be sophisticated in future
    real                  ::  rho
    !
    ! screening factor: measurements are rejected
    ! if square of departure exceeds factor times variance:
    real                 ::  screening_factor
    !
    ! random number generator:
    type(T_Random)       ::  rnd
    !
  contains
    procedure :: Init            => LEKF_Meas_Sat_Set_Init
    procedure :: Done            => LEKF_Meas_Sat_Set_Done
  end type T_LEKF_Meas_Sat_Set

  !
  ! main object with multiple sets
  !
  type T_LEKF_Meas_Sat
    ! number of sets:
    integer                                   ::  nset
    ! corresponding sat output number:
    integer, allocatable                      ::  iout(:)  ! (nset)
    ! info per set:
    type(T_LEKF_Meas_Sat_Set), pointer        ::  set(:)  ! (nset)
  contains
    procedure :: Init            => LEKF_Meas_Sat_Init
    procedure :: Done            => LEKF_Meas_Sat_Done
  end type T_LEKF_Meas_Sat


  ! --- var ----------------------------------

  ! data:
  type(T_LEKF_Meas_Sat)     ::  meas_sat


contains


  ! ====================================================================
  ! ===
  ! === set
  ! ===
  ! ====================================================================


  subroutine LEKF_Meas_Sat_Set_Init( self, rcF, rcbase, name, seed0, status )
  
    use GO, only : TrcFile

    ! --- in/out ---------------------------------
    
    class(T_LEKF_Meas_Sat_Set), intent(out)     ::  self
    type(TrcFile), intent(in)                   ::  rcF
    character(len=*), intent(in)                ::  rcbase
    character(len=*), intent(in)                ::  name
    integer, intent(in)                         ::  seed0
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LEKF_Meas_Sat_Set_Init'
    
    ! --- local ----------------------------------
    
    character(len=64)   ::  rckey
    
    ! --- begin ----------------------------------
    
    ! store:
    self%name = trim(name)
    
    ! prefix:
    rckey = trim(rcbase)//'.'//trim(name)
    
    ! analyse?
    call rcF%Get( trim(rckey)//'.analyse', self%analyse, status )
    IF_NOTOK_RETURN(status=1)
    
    ! localization length scale:
    call rcF%Get( trim(rckey)//'.rho', self%rho, status )
    IF_NOTOK_RETURN(status=1)

    ! screening factor:
    call rcF%Get( trim(rckey)//'.alfa', self%screening_factor, status )
    IF_NOTOK_RETURN(status=1)

    ! init random number generator:
    call self%rnd%Init(  seed0, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LEKF_Meas_Sat_Set_Init


  ! ***


  subroutine LEKF_Meas_Sat_Set_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_LEKF_Meas_Sat_Set), intent(inout)     ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LEKF_Meas_Sat_Set_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! done with random number generator:
    call self%rnd%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LEKF_Meas_Sat_Set_Done


  ! ====================================================================
  ! ===
  ! === sets
  ! ===
  ! ====================================================================


  subroutine LEKF_Meas_Sat_Init( self, rcF, rcbase, status )
  
    use GO        , only : TrcFile
    use GO        , only : goSplitString, goMatchValue
    use LEKF_Data , only : leo

    ! --- in/out ---------------------------------
    
    class(T_LEKF_Meas_Sat), intent(out)         ::  self
    type(TrcFile), intent(in)                   ::  rcF
    character(len=*), intent(in)                ::  rcbase
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LEKF_Meas_Sat_Init'
    
    ! --- local ----------------------------------

    character(len=1024)     ::  line
    character(len=32)       ::  names(10)
    integer                 ::  iset
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (gol,'(a,": init sat measurements ...")') rname; call goPr
    
    ! list with names:
    call rcF%Get( trim(rcbase)//'.names', line, status )
    IF_NOTOK_RETURN(status=1)
    ! split:
    call goSplitString( line, self%nset, names, status )
    IF_NOTOK_RETURN(status=1)
    ! info ..
    write (gol,'(a,":   number of sets: ",i0)') rname, self%nset; call goPr

    ! storage:
    allocate( self%iout(self%nset), stat=status )   
    IF_NOTOK_RETURN(status=1)
    allocate( self%set(self%nset), stat=status )   
    IF_NOTOK_RETURN(status=1)
    ! loop:
    do iset = 1, self%nset
      ! info ..
      write (gol,'(a,":   set ",i0," : ",a)') rname, iset, trim(names(iset)); call goPr
      write (gol, *) 'trim(names(iset)',trim(names(iset)); call goPr
      write (gol, *) 'leo%file_cso',leo%file_cso; call goPr
      ! write (gol, *) 'self%iout(iset)', self%iout(iset); call goPr
      ! search in output:
      call goMatchValue( trim(names(iset)), leo%file_cso, self%iout(iset), status )
      IF_NOTOK_RETURN(status=1)
      ! init, provide set number as random seed:
      call self%set(iset)%Init( rcF, rcbase, names(iset), iset, status )
      IF_NOTOK_RETURN(status=1)
	  write (gol, *) 'self%iout(iset)', self%iout(iset); call goPr
    end do  ! sets
    
    !
    ! * setup state vectors
    !
    
    !
!    call self%sat_state%Init( status, key='ensXXX', &
!                                       description='simulated '//trim(names(i))//' retrievals' )
!    IF_NOTOK_RETURN(status=1)
    
    !
    ! *
    !

    ! ok
    status = 0
    
  end subroutine LEKF_Meas_Sat_Init


  ! ***


  subroutine LEKF_Meas_Sat_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_LEKF_Meas_Sat), intent(inout)         ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LEKF_Meas_Sat_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (gol,'(a,": done with sat measurements ...")'); call goPr

    ! clear:
    deallocate( self%iout, stat=status )   
    IF_NOTOK_RETURN(status=1) 
    deallocate( self%set, stat=status )   
    IF_NOTOK_RETURN(status=1) 

    ! ok
    status = 0
    
  end subroutine LEKF_Meas_Sat_Done


end module LEKF_Meas_SAT
