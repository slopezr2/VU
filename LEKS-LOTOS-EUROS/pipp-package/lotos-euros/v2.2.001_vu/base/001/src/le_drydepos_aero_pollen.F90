!###############################################################################
!
! Boundary-layer deposition resistance (Rb) for Pollen following Sofiev et al, 2006.
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

module LE_DryDepos_Aero_Pollen

  use GO, only : gol, goPr, goErr
  
  implicit none


  ! --- in/out --------------------------------
  
  private
  
  public  ::  LE_DryDepos_Aero_Pollen_vd


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_DryDepos_Aero_Pollen'
  


contains


  ! ===============================================================
  

  ! in : vs, Ra (=ftop)
  ! out: vd, vd
  !input: particle size (diameter) in m.

  subroutine LE_DryDepos_Aero_Pollen_vd( vd, vs, Ra, status )

    use Binas        , only : grav
   
    ! --- in/out ---------------------------------
    
    real, intent(out)        ::  vd       ! deposition velocity (m/s)
    real, intent(in)         ::  vs       ! sedimentation velocity (u in formula 6, form Sofiev, 2006)
    real, intent(in)         ::  Ra       ! atmospheric resistance
    integer, intent(out)     ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_DryDepos_Aero_Pollen_vd'
    

    ! --- local ----------------------------------

    real              ::  Rgrav

    ! --- begin ----------------------------------
        
    ! Gravitational resistance  [Eq (8), Sofiev 2006]
    Rgrav = 1 / vs
        
    ! deposition velocity :
    !             1          1
    !    vd = ---------- + -----
    !         Ra + Rgrav   Rgrav
    !
    vd = 1.0 / ( Ra + Rgrav ) + vs

    ! ok
    status = 0

  end subroutine LE_DryDepos_Aero_Pollen_vd



end module LE_DryDepos_Aero_Pollen

