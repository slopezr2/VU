!###############################################################################
!
! LOTOS-EUROS model data
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

module LE_Data

  use GO            , only : gol, goPr, goErr
  use LE_Data_VData , only : LE_Data_Inquire    => LE_Data_VData_Inquire
  use LE_Data_VData , only : LE_Data_InqVar     => LE_Data_VData_InqVar
  use LE_Data_VData , only : LE_Data_Enable     => LE_Data_VData_Enable
  use LE_Data_VData , only : LE_Data_GetPointer => LE_Data_VData_GetPointer
  use LE_Data_VData , only : LE_Data_Setup      => LE_Data_VData_Setup
  use LE_Data_VarSum, only : T_VarSum
  use LE_Data_Common, only : LE_Data_CompareUnits
  use LE_Data_Common, only : levels_type
  
  implicit none


  ! --- in/out -----------------------------------
  
  private

  public  ::  LE_Data_Init, LE_Data_Done
  public  ::  LE_Data_Inquire
  public  ::  LE_Data_InqVar
  public  ::  LE_Data_Enable
  public  ::  LE_Data_GetPointer
  public  ::  LE_Data_Setup
  public  ::  LE_Data_CompareUnits
  
  public  ::  T_VarSum
  
  public  ::  levels_type


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Data'
  

  ! --- var -------------------------------------- 
  

contains


  ! ====================================================================
  ! ===
  ! === module
  ! ===
  ! ====================================================================
  

  subroutine LE_Data_Init( rcF, status )
  
    use GO            , only : TrcFile
    use LE_Data_Common, only : LE_Data_Common_Init
    use LE_Data_Calc  , only : LE_Data_Calc_Init
    use LE_Data_VData , only : LE_Data_VData_Init
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), pointer                  ::  rcF
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Data_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------
    
    ! init common variables:
    call LE_Data_Common_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! init calculations:
    call LE_Data_Calc_Init( status )
    IF_NOTOK_RETURN(status=1)
    
    ! init variables:
    call LE_Data_VData_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Data_Init
  
  
  ! ***
  
  
  subroutine LE_Data_Done( status )
  
    use LE_Data_Common, only : LE_Data_Common_Done
    use LE_Data_Calc  , only : LE_Data_Calc_Done
    use LE_Data_VData , only : LE_Data_VData_Done
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! done with calculations:
    call LE_Data_Calc_Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with common variables:
    call LE_Data_Common_Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with variables:
    call LE_Data_VData_Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Data_Done
  
  
end module LE_Data

