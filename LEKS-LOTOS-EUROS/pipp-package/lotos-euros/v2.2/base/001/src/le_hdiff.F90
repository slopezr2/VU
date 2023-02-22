!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module hdiff

  implicit none

contains


  ! add horizontal diffusion to the concentration field

  subroutine hordif( c, dt, bc_w, bc_e, bc_s, bc_n )

    use dims , only : nx, ny, nz, khx, khy, volume
    use input, only : runF

    ! --- in/out --------------------------------------------
    
    real, intent(inout)   ::  c(nx,ny,nz)
    real, intent(in)      ::  bc_w(ny,nz), bc_e(ny,nz), bc_s(nx,nz), bc_n(nx,nz)
    real, intent(in)      ::  dt

    ! --- local ----------------------------------------------
    
    real      ::  fluxx(0:nx,ny,nz), fluxy(nx,0:ny,nz)
    integer   ::  i,j,k
    
    ! --- begin ----------------------------------------------

    ! compute advective fluxes
    do j = 1, ny
      do i = 1, nx-1
        fluxx(i,j,:)= dt*khx(i,j,:)*(c(i+1,j,:)-c(i,j,:))/runF%dx(j)
      end do
      ! at the boundaries...
      fluxx(0,j,:)  =  dt*khx(0 ,j,:)*(c(1, j,:)-bc_w(j,:))/runF%dx(j)
      fluxx(nx,j,:) =  dt*khx(nx,j,:)*(bc_e(j,:)-c(nx,j,:))/runF%dx(j)
    end do

    ! add the diffusive fluxes in y-direction to the already
    ! compute advective fluxes
    do i = 1, nx
      do j = 1, ny-1
         fluxy(i,j,:)= dt*khy(i,j,:)*(c(i,j+1,:)-c(i,j,:))/runF%dy
      end do
      fluxy(i,0,:)  =  dt*khy(i,0 ,:)*(c(i, 1,:)-bc_s(i,:))/runF%dy
      fluxy(i,ny,:) =  dt*khy(i,ny,:)*(bc_n(i,:)-c(i,ny,:))/runF%dy
    end do

    ! here update concentrations
    c(1:nx,1:ny,:) = c(1:nx,1:ny,:) + (fluxx(1:nx,1:ny,:)-fluxx(0:nx-1,1:ny,:) + &
                                       fluxy(1:nx,1:ny,:)-fluxy(1:nx,0:ny-1,:) )/volume(1:nx,1:ny,:)

  end subroutine hordif
  

end module hdiff
