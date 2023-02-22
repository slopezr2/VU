!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#include "c3po.inc"
!
!######################################################################

module C3PO_Grid_Tools

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  ll_area_deg_m2
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Grid_Tools'    
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** cell area
  ! ***
  ! ********************************************************************
    
  
  ! given rectangle [west,east]x[south,north] in radians,
  ! compute area in rad^2
  
  real function ll_area( west, east, south, north )
  
    ! --- in/out -------------------------------------
    
    real, intent(in)  ::  west, east, south, north    ! rad
    
    ! --- begin ------------------------------------

    ll_area = ( sin(max(north,south)) - sin(min(north,south)) ) * abs(east-west)
    
  end function ll_area
  
  
  ! *
    
  
  ! given rectangle [west,east]x[south,north] in degrees,
  ! compute area in m^2
  
  real function ll_area_deg_m2( west, east, south, north )

    use Binas, only : deg2rad, ae
  
    ! --- in/out -------------------------------------
    
    real, intent(in)  ::  west, east, south, north    ! deg
    
    ! --- begin ------------------------------------
    
    ll_area_deg_m2 = ll_area( west*deg2rad, east*deg2rad, south*deg2rad, north*deg2rad ) * ae**2
    
  end function ll_area_deg_m2


end module C3PO_Grid_Tools
