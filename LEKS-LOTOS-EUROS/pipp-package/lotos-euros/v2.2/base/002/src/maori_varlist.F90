!#######################################################################
!
! NAME
!
!   MAORI - Model And Output Routine Interface
!
! DESCRIPTION
!
!   Mode independent output routines.
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_VarList

  use GO, only : gol, goPr, goErr

  use MAORI_Param, only : MAORI_LEN_NAME, MAORI_LEN_LONGNAME, MAORI_LEN_UNIT, MAORI_LEN_LINE
  use MAORI_Param, only : MAORI_MAX_VAR

  implicit none


  ! --- in/out ---------------------------

  private

  public  ::  T_MAORI_Var, T_MAORI_VarList
  public  ::  MAORI_VarList_Init, MAORI_VarList_Done
  public  ::  MAORI_VarList_Inq
  public  ::  MAORI_Var_Inq


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI_VarList'


  ! --- types ------------------------------

  ! single variable:
  type T_MAORI_Var
    ! variable name and unit:
    character(len=MAORI_LEN_NAME)       ::  name
    character(len=MAORI_LEN_LONGNAME)   ::  cf_standard_name
    character(len=MAORI_LEN_LONGNAME)   ::  cf_long_name
    character(len=MAORI_LEN_UNIT)       ::  cf_units
    real                                ::  cf_units_conversion_factor
    character(len=MAORI_LEN_LINE)       ::  comments
    logical                             ::  halflevel
  end type T_MAORI_Var

  ! list
  type T_MAORI_VarList
    ! number of elements:
    integer                         ::  n
    ! list:
    type(T_MAORI_Var)               ::  var(MAORI_MAX_VAR)
  end type T_MAORI_VarList


contains


  ! ==================================================================


  subroutine MAORI_Var_Init( var, name, status )

    use LE_CF_Conventions, only : LE_CF_names
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Var), intent(out)        ::  var
    character(len=*), intent(in)          ::  name
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Var_Init'

    ! --- local -----------------------------------

    ! --- begin -----------------------------------

    ! store name:
    var%name = trim(name)
    
    ! fill other fields following cf conventions;
    ! use '<cf>' for current unit to avoid error messages ..
    call LE_CF_names( name, '<cf>', &
                                 var%cf_standard_name, &
                                 var%cf_long_name, &
                                 var%cf_units, var%cf_units_conversion_factor, &
                                 var%comments, &
                                 status, &
                                 halflevel=var%halflevel )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Var_Init


  ! ***


  subroutine MAORI_Var_Done( var, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Var), intent(inout)    ::  var
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Var_Done'

    ! --- begin -----------------------------------

    ! clear:
    var%name = ''

    ! ok
    status = 0

  end subroutine MAORI_Var_Done



  ! ==================================================================


  subroutine MAORI_VarList_Init( varlist, key, status )

    use GO, only : goSplitString

    ! --- in/out ---------------------------------

    type(T_MAORI_VarList), intent(out)    ::  varlist
    character(len=*), intent(in)          ::  key
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_VarList_Init'

    ! --- local -----------------------------------

    character(len=MAORI_LEN_NAME)   ::  varnames(MAORI_MAX_VAR)
    integer                         ::  ivar

    ! --- begin -----------------------------------

    ! store list:
    call goSplitString( key, varlist%n, varnames, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over extracted fields:
    do ivar = 1, varlist%n
      ! init variable:
      call MAORI_Var_Init( varlist%var(ivar), varnames(ivar), status )
      IF_NOTOK_RETURN(status=1)
    end do   ! variables

    ! ok
    status = 0

  end subroutine MAORI_VarList_Init


  ! ***


  subroutine MAORI_VarList_Done( varlist, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_VarList), intent(inout)    ::  varlist
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_VarList_Done'

    ! --- local -----------------------------------

    integer                         ::  ivar

    ! --- begin -----------------------------------

    ! loop over variables:
    do ivar = 1, varlist%n
      ! init variable:
      call MAORI_Var_Done( varlist%var(ivar), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! clear:
    varlist%n = 0

    ! ok
    status = 0

  end subroutine MAORI_VarList_Done



  ! ***


  subroutine MAORI_VarList_Inq( varlist, status, nvar )

    ! --- in/out ---------------------------------

    type(T_MAORI_VarList), intent(in)       ::  varlist
    integer, intent(out)                    ::  status

    integer, intent(out), optional          ::  nvar

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_VarList_Inq'

    ! --- begin ----------------------------------

    ! return values:
    if ( present(nvar) ) nvar = varlist%n

    ! ok
    status = 0

  end subroutine MAORI_VarList_Inq


  ! ***


  subroutine MAORI_Var_Inq( varlist, ivar, status, &
                                 name, unit )

    ! --- in/out ---------------------------------

    type(T_MAORI_VarList), intent(in)       ::  varlist
    integer, intent(in)                     ::  ivar
    integer, intent(out)                    ::  status

    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  unit

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Var_Inq'

    ! --- begin -----------------------------------

    ! check:
    if ( (ivar < 1) .or. (ivar > varlist%n) ) then
      write (gol,'("var index ",i6," outside valid range 1, ..",i6)') ivar, varlist%n; call goErr
      TRACEBACK; status=1; return
    end if

    ! return values:
    if ( present(name) ) name = trim(varlist%var(ivar)%name)
    if ( present(unit) ) unit = trim(varlist%var(ivar)%cf_units)

    ! ok
    status = 0

  end subroutine MAORI_Var_Inq


end module MAORI_VarList
