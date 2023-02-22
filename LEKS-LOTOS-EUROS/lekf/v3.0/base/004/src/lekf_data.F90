!###############################################################################
!
! NAME
!
!   LEKF_Data  -  data module for LEKF
!
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

module LEKF_Data

  use LE_Output, only : T_LE_Output
#ifdef with_kf_meas_maori
  use LE_MAORI , only : T_MAORI_Data
#endif
  
  implicit none


  ! --- in/out -----------------------------
  
  public
  
  
  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LEKF_Data'
  
  
  ! --- var --------------------------------
  
  ! storage of output variables:
  type(T_LE_Output)           ::  leo
#ifdef with_kf_meas_maori
  type(T_MAORI_Data)          ::  mad
#endif

  ! flags:
  logical                     ::  lekfo_replace
  logical                     ::  lekfo_with_xi

end module LEKF_Data
    
