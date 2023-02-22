!###############################################################################
!
! LE_Bound - LOTOS-EUROS boundary condition routines
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

module LE_Bound_Common

  use GO     , only : gol, goErr, goPr

  implicit none
  

  ! --- in/out --------------------------------

  private

  public  ::  LE_Bound_Common_Init, LE_Bound_Common_Done

  public  ::  bc_west, bc_east, bc_south, bc_north
  public  ::  caloft
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_Common'


  ! --- var --------------------------------

  ! storage for boundary fields:
  real, allocatable     ::  bc_west (  :,:,:)  ! (   ny,nlev           ,nspec)
  real, allocatable     ::  bc_east (  :,:,:)  ! (   ny,nlev           ,nspec)
  real, allocatable     ::  bc_north(:,  :,:)  ! (nx,   nlev           ,nspec)
  real, allocatable     ::  bc_south(:,  :,:)  ! (nx,   nlev           ,nspec)
  real, allocatable     ::  caloft  (:,:,:,:)  ! (nx,ny,nlev+1:nlev_top,nspec)



contains


  ! ====================================================================
  ! ===
  ! === module init/done
  ! ===
  ! ====================================================================


  subroutine LE_Bound_Common_Init( status )
  
    use Dims          , only : nx, ny
    use LE_Data_Common, only : nlev, nlev_top
    use Indices       , only : nspec
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Bound_Common_Init'
      
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! storage:
    allocate( bc_west(ny,nlev,nspec), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_east(ny,nlev,nspec), stat=status ) 
    IF_NOTOK_RETURN(status=1)
    allocate( bc_north(nx,nlev,nspec), stat=status ) 
    IF_NOTOK_RETURN(status=1)
    allocate( bc_south(nx,nlev,nspec), stat=status ) 
    IF_NOTOK_RETURN(status=1)
    allocate( caloft(nx,ny,nlev+1:nlev_top,nspec), stat=status ) 
    IF_NOTOK_RETURN(status=1)

    ! force all boundary concentrations to zero:
    bc_west  = 0.0
    bc_east  = 0.0
    bc_north = 0.0
    bc_south = 0.0
    caloft   = 0.0
    
    ! ok
    status = 0
    
  end subroutine LE_Bound_Common_Init


  ! ***


  subroutine LE_Bound_Common_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Bound_Common_Done'
      
    ! --- local ----------------------------------
    
    integer       ::  itracer
      
    ! --- begin ----------------------------------

    ! clear:
    deallocate( bc_west, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( bc_east, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( bc_north, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( bc_south, stat=status ) 
    IF_NOTOK_RETURN(status=1)
    deallocate( caloft, stat=status ) 
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Bound_Common_Done


end module LE_Bound_Common


