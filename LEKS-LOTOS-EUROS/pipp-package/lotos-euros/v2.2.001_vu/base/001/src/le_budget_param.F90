!###############################################################################
!
! LE_Budget_Param  -  parameters shared by dry- and wet-Budgetition
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

module LE_Budget_Param

  use Indices, only : i_SOxa, i_NOya, i_NHxa

  implicit none


  ! --- in/out --------------------------------
  
  private

  public  ::  nbud
  public  ::  ibud_sox, ibud_noy, ibud_nhx
  public  ::  ispec_bud
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Budget_Param'


  ! --- const ------------------------------------
  
  ! Set groups of species for which to compute mass budgets:
  integer, parameter :: nbud = 3
  ! indices:
  integer, parameter :: ibud_SOx = 1  ! mass budget of SOx
  integer, parameter :: ibud_NOy = 2  ! mass budget of NOy
  integer, parameter :: ibud_NHx = 3  ! mass budget of NHx

  ! corresponding accumulated tracers; 
  ! will be filled in LE_Budget_Init :
  integer  ::  ispec_bud(nbud)   ! (/ i_SOxa, i_NOya, i_NHxa /)


end module LE_Budget_Param

