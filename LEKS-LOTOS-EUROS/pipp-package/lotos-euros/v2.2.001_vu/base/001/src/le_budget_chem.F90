!
! Unused budget codes, kept here for future inspiration.
!
!!###############################################################################
!!
!! Chemistry budgets
!!
!!###############################################################################
!!
!#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
!#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!!
#include "le.inc"
!!
!!###############################################################################
!
!module LE_Budget_Chem
!
!  use GO, only : gol, goPr, goErr
!
!  implicit none
!
!
!  ! --- in/out --------------------------------
!  
!  private
!
!  public  ::  Chem_Budget_Init, Chem_Budget_Done
!  
!
!  ! --- const --------------------------------
!
!  character(len=*), parameter ::  mname = 'LE_Budget_Chem'
!
!  integer, parameter :: nproces = 5
!  integer, parameter :: nreac = 8
!
!
!  ! --- var --------------------------------------
!  
!  real, allocatable, dimension(:,:,:,:) :: budget
!  real, allocatable, dimension(:,:,:) :: budget_chem
!
!
!contains
!
!
!
!  ! ========================================================================
!  
!
!  subroutine Chem_Budget_Init( status )
!  
!    use Dims, only : nx, ny, nz
!    use Indices, only : nspec
!  
!    ! --- in/out ---------------------------------
!    
!    integer, intent(out)                    ::  status
!
!    ! --- const --------------------------------
!    
!    character(len=*), parameter   :: rname = mname//'/Chem_Budget_Init'
!    
!    ! --- local ------------------------------------------
!    
!    ! --- begin ------------------------------------------
!
!    ! storage:
!    allocate( budget(nx,ny,nspec,nproces) ) 
!    allocate( budget_chem(nx,ny,nreac) )
!
!    ! ok
!    status = 0
!    
!  end subroutine Chem_Budget_Init
!  
!  
!  ! ***
!  
!  
!  subroutine Chem_Budget_Done( status )
!  
!    ! --- in/out ---------------------------------
!    
!    integer, intent(out)                    ::  status
!    
!    ! --- const -------------------------------
!    
!    character(len=*), parameter ::  rname = mname//'/Chem_Budget_Done'
!    
!    ! --- begin ----------------------------------
!
!    ! clear:
!    deallocate( budget ) 
!    deallocate( budget_chem )
!
!    ! ok
!    status = 0
!    
!  end subroutine Chem_Budget_Done
!
!
!  
!end module LE_Budget_Chem
