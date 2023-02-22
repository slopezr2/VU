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

  use LE_Output, only : T_LE_Output, T_LE_OutputState
#ifdef with_kf_meas_maori
  use LE_MAORI , only : T_MAORI_Data, T_MAORI_Output
#endif
  
  implicit none


  ! --- in/out -----------------------------
  
  public
  
  
  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LEKF_Data'
  
  
  ! --- var --------------------------------
  
  ! storage of output variables:
  type(T_LE_Output)        ::  leo
  type(T_LE_OutputState)   ::  leos_xb
  type(T_LE_OutputState)   ::  leos_xf, leos_sf  ! forecast
  type(T_LE_OutputState)   ::  leos_xa, leos_sa  ! analysis
  
#ifdef with_kf_meas_maori
  ! storage of output variables:
  type(T_MAORI_Data)       ::  mad
  type(T_MAORI_Output)     ::  mao_xb           ! model
  type(T_MAORI_Output)     ::  mao_xf, mao_sf   ! filter mean and std.dev., forecast
  type(T_MAORI_Output)     ::  mao_xa, mao_sa   ! filter mean and std.dev., analysis
#endif

  ! flags:
  logical                  ::  lekfo_replace


end module LEKF_Data
    
