!###############################################################################
!
! NAME
!   LE_Bound_Clim_Steady
!
! DESCRIPTION
!
!   Apply steady-state assumption to O3/NO/NO2 in boundary conditions.
!
! HISTORY
!   Taken from form 'logan' module in 'bound_logan.F90' .
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

module LE_Bound_Clim_Steady

  use GO, only : gol, goErr, goPr
  
  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  bc_steady
  

  ! --- const ------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_Clim_Steady'


contains

  
  ! ==========================================================
  
  
  !
  !   spec_filled(1:nspec)  : set to .true. if bc is filled, unchanged otherwise
  !

  subroutine bc_steady( rk1rk3, spec_filled )

    use LE_Logging, only : ident1
    use dims, only : nx, ny, nz, nspec
    use indices
    use LE_Data_Common , only : nlev
    use LE_Bound_Common, only : bc_west, bc_east, bc_south, bc_north, caloft

    ! --- in/out ----------------------------
    
    real, intent(in)            ::  rk1rk3(nx,ny)
    logical, intent(inout)      ::  spec_filled(nspec)
    
    ! --- local -----------------------------
    
    integer :: i,j,k
    
    ! --- begin -----------------------------
    
    if ( all((/i_o3,i_no,i_no2/) > 0) ) then

      !print *, ident1,'putting BC values for O3/NO2/NO in steady state'

      ! west and east
      do j=1,ny
      do k=1,nz
        ! apply photo steady-state to obtain O3/NO2/NO
        call getsteady(rk1rk3(1,j),  bc_west(j,k,i_o3), bc_west(j,k,i_no2),bc_west(j,k,i_no) )
        call getsteady(rk1rk3(nx,j), bc_east(j,k,i_o3), bc_east(j,k,i_no2),bc_east(j,k,i_no) )
      enddo
      enddo

      ! north and south
      do i=1,nx
      do k=1,nz
        ! apply photo steady-state to obtain O3/NO2/NO
        call getsteady(rk1rk3(i,1),  bc_south(i,k,i_o3), bc_south(i,k,i_no2),bc_south(i,k,i_no) )
        call getsteady(rk1rk3(i,ny), bc_north(i,k,i_o3), bc_north(i,k,i_no2),bc_north(i,k,i_no) )
      enddo
      enddo

      ! caloft
      do i=1,nx
      do j=1,ny
        call getsteady(rk1rk3(i,j),  caloft(i,j,nlev+1,i_o3), caloft(i,j,nlev+1,i_no2), caloft(i,j,nlev+1,i_no) )
      enddo
      enddo

      ! set flags:
      spec_filled(i_o3 ) = .true.
      spec_filled(i_no2) = .true.
      spec_filled(i_no ) = .true.
      
    end if

  end subroutine
  
  
  ! ***
  

  subroutine getsteady(rk1rk3,o3,no2,no)

    real :: nox, ox, discr, o3, no2, no, rk1rk3

    ox  = no2 + o3
    nox = no+no2
    discr = max( 0.0, (nox+ox+rk1rk3)**2-4.0*nox*ox)
    no2 = max(0.0, 0.5*(nox+ox+rk1rk3 - sqrt(discr)))
    no2 = min(no2,nox,ox)
    o3  = max(0.0, ox-no2)
    no  = max(0.0,nox-no2)

  end subroutine


end module LE_Bound_Clim_Steady

