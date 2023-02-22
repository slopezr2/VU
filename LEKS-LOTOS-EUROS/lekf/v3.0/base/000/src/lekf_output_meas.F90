!###############################################################################
!
! NAME
!
!   KF_Output_Meas  -  write measurements
!
! HISTORY
!
!   2007 may, Arjo Segers, TNO
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

module LEKF_Output_Meas

  use GO, only : gol, goPr, goErr
  
#ifdef with_kf_meas_pm
  use LEKF_Meas_PM    , only : meas_output_pm
#endif
#ifdef with_kf_meas_pm_m24
  use LEKF_Meas_PM_m24, only : meas_output_pm_m24
#endif
#ifdef with_kf_meas_aod
  use LEKF_Meas_AOD, only : meas_output_aod
#endif
#ifdef with_kf_meas_ground
  use LEKF_Meas_Ground, only : meas_output_ground
#endif
#ifdef with_kf_meas_omi_trc
  use LEKF_Meas_OMI_TRC, only : meas_output_omi_trc
#endif

  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  T_LEKF_Output_Meas


  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'KF_Output_Meas'
  
    
  ! --- types ------------------------------
  
  type T_LEKF_Output_Meas
    ! to avoid errors about empty structures ...
    integer                      ::  dummy
    ! output for specific measurements:
#ifdef with_kf_meas_pm
    type(meas_output_pm)         ::  mo_pm
#endif
#ifdef with_kf_meas_pm_m24
    type(meas_output_pm_m24)     ::  mo_pm_24
#endif
#ifdef with_kf_meas_aod
    type(meas_output_aod)        ::  mo_aod
#endif
#ifdef with_kf_meas_ground
    type(meas_output_ground)     ::  mo_ground
#endif
#ifdef with_kf_meas_omi_trc
    type(meas_output_omi_trc)    ::  mo_omi_trc
#endif
  contains
    procedure   ::  Init        =>  kfo_meas_Init
    procedure   ::  Done        =>  kfo_meas_Done
    procedure   ::  PutOut      =>  kfo_meas_PutOut
  end type T_LEKF_Output_Meas

  
contains


  ! ====================================================
  
  
  subroutine kfo_meas_Init( self, rcfile, rckey, status )

#ifdef with_kf_meas_pm
    use LEKF_Meas_PM    , only : Init
#endif
#ifdef with_kf_meas_pm_m24
    use LEKF_Meas_PM_m24, only : Init
#endif
#ifdef with_kf_meas_aod
    use LEKF_Meas_AOD, only : Init
#endif
#ifdef with_kf_meas_ground
    use LEKF_Meas_Ground, only : Init
#endif
  
    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_Meas), intent(out)  ::  self
    character(len=*), intent(in)            ::  rcfile
    character(len=*), intent(in)            ::  rckey
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/kfo_meas_Init'
    
    ! --- local ---------------------------------
    
    ! --- begin ---------------------------------
    
    ! dummy to have at least some output ...
    self%dummy = 0
    
#ifdef with_kf_meas_pm
    call Init( self%mo_pm   , rcfile, rckey, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_pm_m24
    call Init( self%mo_pm_24, rcfile, rckey, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_aod
    call Init( self%mo_aod, rcfile, rckey, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_ground
    call Init( self%mo_ground, rcfile, rckey, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_omi_trc
    call self%mo_omi_trc%Init( rcfile, rckey, status )
    IF_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine kfo_meas_Init
  
  
  ! ***
  

  subroutine kfo_meas_Done( self, status )
  
#ifdef with_kf_meas_pm
    use LEKF_Meas_PM    , only : Done
#endif
#ifdef with_kf_meas_pm_m24
    use LEKF_Meas_PM_m24, only : Done
#endif
#ifdef with_kf_meas_aod
    use LEKF_Meas_AOD, only : Done
#endif
#ifdef with_kf_meas_ground
    use LEKF_Meas_Ground, only : Done
#endif

    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_Meas), intent(inout)  ::  self
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/kfo_meas_Done'
    
    ! --- begin ---------------------------------
    
#ifdef with_kf_meas_pm
    call Done( self%mo_pm   , status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_pm_m24
    call Done( self%mo_pm_24, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_aod
    call Done( self%mo_aod, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_ground
    call Done( self%mo_ground, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_omi_trc
    call self%mo_omi_trc%Done( status )
    IF_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine kfo_meas_Done
  
  
  ! ***
  
  
  subroutine kfo_meas_PutOut( self, key, t, status, last )
  
    use GO     , only : TDate
#ifdef with_kf_meas_pm
    use LEKF_Meas_PM    , only : PutOut
#endif
#ifdef with_kf_meas_pm_m24
    use LEKF_Meas_PM_m24, only : PutOut
#endif
#ifdef with_kf_meas_aod
    use LEKF_Meas_AOD, only : PutOut
#endif
#ifdef with_kf_meas_ground
    use LEKF_Meas_Ground, only : PutOut
#endif

    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_Meas), intent(inout)  ::  self
    character(len=*), intent(in)              ::  key
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status
    logical, intent(in), optional             ::  last
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/kfo_meas_PutOut'
    
    ! --- local ---------------------------------
    
    ! --- begin ---------------------------------
    
#ifdef with_kf_meas_pm
    call PutOut( self%mo_pm   , key, t, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_pm_m24
    call PutOut( self%mo_pm_24, key, t, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_aod
    call PutOut( self%mo_aod, key, t, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_ground
    call PutOut( self%mo_Ground, key, t, status )
    IF_NOTOK_RETURN(status=1)
#endif
#ifdef with_kf_meas_omi_trc
    call self%mo_omi_trc%PutOut( key, t, status, last=last )
    IF_NOTOK_RETURN(status=1)
#endif

    ! ok
    status = 0
    
  end subroutine kfo_meas_PutOut


end module LEKF_Output_Meas
