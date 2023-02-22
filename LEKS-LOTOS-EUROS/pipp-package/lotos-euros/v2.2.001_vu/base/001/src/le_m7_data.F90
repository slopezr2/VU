!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_M7_Data

  use GO, only : gol, goPr, goErr

  ! M7 parameters:
  use mo_aero_m7, only : n_m7numb => nmod     ! m7 number modes
  use mo_aero_m7, only : n_m7mass => naermod  ! m7 mass modes

  implicit none


  ! --- in/out -----------------------------------

  private

  public    ::  LE_M7_Data_Init, LE_M7_Data_Done
  
  public    ::  n_m7numb, n_m7mass
  public    ::  densm7modes
  public    ::  waterm7modes
  public    ::  rdrym7modes
  public    ::  rwetm7modes
  public    ::  vsmass, vsnumber
!  public    ::  vdmass, vdnumber
!  public    ::  vdtable


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_M7_Data'


  ! --- var --------------------------------------
  
  !M7 arrays for density, water content and deposition velocieties
  real, allocatable   ::  densm7modes (:,:,:,:)
  real, allocatable   ::  waterm7modes(:,:,:,:)
  real, allocatable   ::  rdrym7modes (:,:,:,:)
  real, allocatable   ::  rwetm7modes (:,:,:,:)
  real, allocatable   ::  vsmass(:,:,:,:), vsnumber(:,:,:,:)
!  real, allocatable   ::  vdmass(:,:,:), vdnumber(:,:,:)


contains


  ! ====================================================================
  ! ===
  ! === module init/done
  ! ===
  ! ====================================================================


  subroutine LE_M7_Data_Init( status )

    use dims   , only : nx, ny, nz
    use indices, only : n_m7

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_Init'

    ! --- begin ----------------------------------

    ! enabled ?
    if ( n_m7 > 0 ) then

      !m7 extras,dry and wet radius, density and water content of number modes
      !Note that only for soluble modes wet and dry radius are different.
      allocate( rdrym7modes (nx,ny,nz,4       ) )
      allocate( rwetm7modes (nx,ny,nz,n_m7numb) ) 
      allocate( densm7modes (nx,ny,nz,n_m7numb) ) 
      allocate( waterm7modes(nx,ny,nz,n_m7numb) ) 
      allocate( vsmass      (nx,ny,nz,n_m7numb) ) 
      allocate( vsnumber    (nx,ny,nz,n_m7numb) )
!      allocate( vdmass      (nx,ny   ,n_m7numb) ) 
!      allocate( vdnumber    (nx,ny   ,n_m7numb) ) 

    end if

    ! ok
    status = 0

  end subroutine LE_M7_Data_Init


  ! ***


  subroutine LE_M7_Data_Done( status )

    use indices, only : n_m7

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_Data_Done'

    ! --- begin ----------------------------------

    ! enabled ?
    if ( n_m7 > 0 ) then

      ! clear:
      deallocate( rdrym7modes ) 
      deallocate( rwetm7modes )
      deallocate( densm7modes )
      deallocate( waterm7modes )
      deallocate( vsmass )
      deallocate( vsnumber )
!      deallocate( vdmass )
!      deallocate( vdnumber ) 

    end if

    ! ok
    status = 0

  end subroutine LE_M7_Data_Done


end module LE_M7_Data

