!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_HDiff

  use GO, only : gol, goPr, goErr
  implicit none

  character(len=*), parameter   ::  mname = 'LE_HDiff'

contains

  !  *** 
  !  Mass-diffusion coefficients in x/y direction
  !  ***
  
  subroutine difcoef(kh,status)
 
    use dims   , only : nx, ny, nz
    use dims   , only : khx, khy
    use LE_Grid, only : ugg
    use LE_Meteo_Data, only : afluxx, afluxy, afluxz
    use LE_Data      , only : LE_Data_GetPointer

    ! --- in/out ---
    
    real, intent(in)     ::  kh(0:nx+1,0:ny+1,1:nz)
    integer, intent(out) :: status 
    
    ! --- const ---
    
    character(len=*), parameter   ::  rname = mname//'/difcoef'

    ! the hordif constant is set for a grif of 100kmx100km.To set the correct scale dependent
    ! difcof we need to use the following formula: K = Kref * sqrt(dx/dxref). The same formula is applied
    ! to difmin and difmax.    
    real, parameter     ::  dx_ref = 1.e6
    real, parameter     ::  dy_ref = 1.e6 
    real, parameter     ::  difcof = 9000.0
    real, parameter     ::  difmin = 1e1
    real, parameter     ::  difmax = 1e5
    
    ! --- local ---
    integer             ::  i, j, k    
    real                ::  xrate, yrate
    
    ! 'meteo'-input
    real, pointer       ::  dx(:,:,:)
    real, pointer       ::  dy(:,:,:)
    real, pointer       ::  dh(:,:,:)

    ! --- begin ----------------------------------------------

    call LE_Data_GetPointer( 'dxu', dy, status, check_units = 'm')    
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_GetPointer( 'dxv', dx, status, check_units = 'm')    
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_GetPointer( 'dh' , dh, status, check_units = 'm' )
    IF_NOTOK_RETURN(status=1)
    
    ! scale rates
    xrate = sqrt( maxval(dx)/dx_ref )
    yrate = sqrt( maxval(dy)/dy_ref )

    ! for each layer
    do k=1,nz
    
      ! now make the diffusive mass fluxes at cell edges
      ! x-direction
      do j=1,ny
        do i=0,nx
          khx(i,j,k)=0.5*(kh(i,j,k)+kh(i+1,j,k))*dh(i,j,k)*dy(i,j,1)*difcof*xrate

          ! limit:
          khx(i,j,k)= min(khx(i,j,k), xrate*difmax)
          khx(i,j,k)= max(khx(i,j,k), xrate*difmin)
        enddo

      enddo

      ! y-direction
      do i=1,nx
        do j=0,ny
          khy(i,j,k)=0.5*(kh(i,j,k)+kh(i,j+1,k))*dh(i,j,k)*dx(i,j,1)*difcof*yrate
          ! limit:
          khy(i,j,k)= min(khy(i,j,k), yrate*difmax)
          khy(i,j,k)= max(khy(i,j,k), yrate*difmin)
        enddo
      enddo

    enddo  ! levels

    ! ok
    status = 0

  end subroutine  ! difcoef

  
  ! ***
  ! apply horizontal diffusion to the concentration field
  ! ***
  
  subroutine hordif( c, dt, bc_w, bc_e, bc_s, bc_n,ispec, status )

    use dims         , only : nx, ny, nz, khx, khy
    use LE_Meteo_Data, only : volume
    use LE_Data      , only : LE_Data_GetPointer
#ifdef with_labeling
    use SA_Labeling, only : SA_Hdiff_Apply
#endif

    character(len=*), parameter   ::  rname = mname//'/hordif'

    ! --- in/out --------------------------------------------
    
    real, intent(inout)   ::  c(nx,ny,nz)
    real, intent(in)      ::  bc_w(ny,nz), bc_e(ny,nz), bc_s(nx,nz), bc_n(nx,nz)
    real, intent(in)      ::  dt
    integer, intent(in)   ::  ispec
    integer, intent(out)  ::  status 
    
    ! --- local ----------------------------------------------
    
    real      ::  fluxx(0:nx,ny,nz), fluxy(nx,0:ny,nz)
    integer   ::  i,j,iz
    
    real, pointer        ::  dx(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::  dy(:,:,:)   ! (lon,lat,lev)
    
    ! --- begin ----------------------------------------------

    call LE_Data_GetPointer( 'dxu', dy, status, check_units = 'm')    
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_GetPointer( 'dxv', dx, status, check_units = 'm')    
    IF_NOTOK_RETURN(status=1)        
    
    ! compute advective fluxes
    do j = 1, ny
      do i = 1, nx-1
        fluxx(i,j,:)= dt * khx(i,j,:)*(c(i+1,j,:)-c(i,j,:))/dx(i,j,1)
      end do
      ! at the boundaries...
      fluxx(0,j,:)  = dt * khx(0 ,j,:)*(c(1, j,:)-bc_w(j,:))/dx(1,j,1)
      fluxx(nx,j,:) = dt * khx(nx,j,:)*(bc_e(j,:)-c(nx,j,:))/dx(nx,j,1)
    end do

    ! add the diffusive fluxes in y-direction to the already
    ! compute advective fluxes
    do i = 1, nx
      do j = 1, ny-1
         fluxy(i,j,:)= dt * khy(i,j,:)*(c(i,j+1,:)-c(i,j,:))/dy(i,j,1)
      end do
      fluxy(i,0,:)  =  dt * khy(i,0 ,:)*(c(i, 1,:)-bc_s(i,:))/dy(i,1,1)
      fluxy(i,ny,:) =  dt * khy(i,ny,:)*(bc_n(i,:)-c(i,ny,:))/dy(i,ny,1)
    end do
#ifdef with_labeling
    ! Apply labeling on Hdiff
    do iz = 1, nz
      call SA_HDiff_Apply( c(1:nx,1:ny,iz), fluxx(0:nx,1:ny,iz), dx(1:nx,1:ny,1), fluxy(1:nx,0:ny,iz), dy(1:nx,1:ny,1), iz, ispec, status )
      IF_NOTOK_RETURN(status=1)
    end do
#endif       

    ! here update concentrations
    do iz = 1, nz
      c(1:nx,1:ny,iz) = c(1:nx,1:ny,iz) + &
                          (fluxx(1:nx,1:ny,iz)-fluxx(0:nx-1,1:ny,iz)) / dx(1:nx,1:ny,1) + &
                          (fluxy(1:nx,1:ny,iz)-fluxy(1:nx,0:ny-1,iz)) / dy(1:nx,1:ny,1)
    end do                        

    ! ok
    status = 0
    
  end subroutine hordif
    
end module LE_Hdiff
