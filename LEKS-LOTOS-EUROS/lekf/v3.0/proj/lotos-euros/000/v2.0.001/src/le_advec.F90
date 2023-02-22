!###############################################################################
!
! NAME
!   LE_Advec   -   advection
!
! METHOD
!
!   Following Walcek (2000?)
!
! DISCUSSION
!
!  Volume-mixing-ratio's or mass-mixing ratio's ?
!
!  In LOTOS-EUROS all units are in terms of volume:
!    volume       : cell volume [(m3 air)]
!    aflux[xyz]   : air volume fluxes throug cell interfaces [(m3 air)/s]
!    c            : volume mixing ratio's  [(m3 tracer)/(m3 air)]
!
!  Advection scheme is logically defined in terms of masses:
!    mass         : cell air mass [(kg air)]
!    mflux[xyz]   : air mass fluxes through cell interfaces [(kg air)/s]
!    c            : mass mixing ratio's
!
!  The later is much easer to understand, and does not need the
!  fuzy temporary densities.
!
! HISTORY
!
!   20??, ??, ??
!     Original.
!   2011-10, Arjo Segers, TNO
!     Use loop over list with advected tracers instead of a loop
!     from 1 to 'nadvect' ; this allows tracers to be in any order.
!   2011-10, Arjo Segers, TNO
!     Finally fixed the bug in the testing if flow in both edges is negative.
!     This bug has been discovered in the past by many users already
!     (Elja Huijbrechts, Richard Kranenburg, Arjo Segers), but for some
!     reason the fix never made it into the versions.
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

module LE_Advec

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public   ::  LE_Advec_Init, LE_Advec_Done
  public   ::  LE_Advec_Get_NStep
  public   ::  LE_Advec_Apply
  
  
  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Advec'

  ! maximum allowed courant number for the numerical scheme
  real, parameter :: courant_max = 1.0

  !the initial density
  real, parameter :: d0=1.0
  
  
  ! --- var --------------------------------------

  real, allocatable  :: inv_volume(:,:,:)

  !! timers:
  !integer   ::  itim_adv_x
  !integer   ::  itim_adv_y
  !integer   ::  itim_adv_z
  !integer   ::  itim_adv_dif


contains


  ! ====================================================================
  
  
  subroutine LE_Advec_Init( status )

    use dims, only : nx, ny, nz
    !use LE_Timers, only : Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ LE_Advec_Init'

    ! number of halo cells:
    integer, parameter  ::  nh = 2
    
    ! --- begin ----------------------------------
    
    !! define timers:
    !call Timer_Def( itim_adv_x  , 'advection x' )
    !call Timer_Def( itim_adv_y  , 'advection y' )
    !call Timer_Def( itim_adv_z  , 'advection z' )
    !call Timer_Def( itim_adv_dif, 'advection dif' )

    ! ok
    status = 0

  end subroutine LE_Advec_Init


  ! ***
  
  
  subroutine LE_Advec_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ LE_Advec_Done'

    ! --- begin ----------------------------------
    
    ! ok
    status = 0
    
  end subroutine LE_Advec_Done


  ! ***
  
  
  ! Return minimum number of advections steps
  ! that are necessary within a time step of 'dt_min'
  
  subroutine LE_Advec_Get_NStep( dt_min, nstep, status )

    use Dims, only : nx, ny, nz
    use LE_Meteo_Data, only : volume
    
    ! --- in/out ---------------------------------

    real, intent(in)                ::  dt_min  ! minutes
    integer, intent(out)            ::  nstep
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Advec_Get_NStep'

    ! --- local ----------------------------------
    
    real       ::  dta

    ! --- begin ----------------------------------

    ! storage:
    allocate( inv_volume(nx, ny, nz) )

    ! pre-compute 1.0/volume since devision is much more expensive than multiplication ...
    inv_volume = 1.0 / volume

    ! check courant condition, store cfls numbers:
    call courant( nx,ny,nz, dt_min, dta, nstep, store=.true. )

    ! clear:
    deallocate( inv_volume )

    ! ok
    status = 0
    
  end subroutine LE_Advec_Get_NStep


  ! ***


  subroutine LE_Advec_Apply( c, dt, status )

    use Indices, only : n_advected, ispecs_advected
    use dims, only : nx, ny, nz, nspec, &
        bc_west, bc_east, bc_north, bc_south, caloft, runF, outF
    use LE_Logging, only : ident2
    use LE_Meteo_Data, only : volume
    use LE_Meteo_Data, only : afluxx, afluxy, afluxz
#ifdef with_hdiff
    use khx, khy
    use hdiff, only : hordif
#endif
#ifdef with_labeling  
    use SA_Labeling, only : SA_Advec_Setup
#endif    
    !use LE_Timers, only : Timer_Start, Timer_End
    
    ! --- begin ------------------------------
    
    real, intent(in)       ::  dt
    real, intent(inout)    ::  c(nx,ny,nz,nspec)
    integer, intent(out)   ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Advec_Apply'
    
    ! --- local ------------------------------
    
    integer   ::  istep, ix, iy, iz
    integer   ::  itr, k
    integer, parameter  ::  nh = 2

    ! actual time step and number of steps for advection:
    real      ::  dta
    integer   ::  nstep_a
    
#ifdef with_hdiff
    ! actual time step and number of steps for diffusion:
    real      ::  dtd
    integer   ::  nstep_d
#endif    
    
    real, allocatable  :: ch(:,:,:)
    real, allocatable  :: dens0(:,:,:)
    real, allocatable  :: dens1(:,:,:)
    real, allocatable  :: inv_d2(:,:,:)
    
#ifdef _OPENMP    
    ! status argument for parallel computing
    integer            :: status_par
#endif    

    ! --- begin ----------------------------------

    ! info:
    if ( .not. outF%suppress ) then
      write (gol,'(a,"<advection>")') ident2; call goPr
    end if

    ! status for parallel computing
#ifdef _OPENMP    
    status_par = 0
#endif    
    !AJS: routine is called for each mode with first flag, thus moved to kf_lotdriver
    !if (runF%first) call init_advec

    ! storage:
    allocate( inv_volume  (nx, ny, nz) )

    ! pre-compute 1.0/volume since devision is much more expensive than multiplication ...
    inv_volume = 1.0 / volume

    ! check courant condition
#ifdef with_hdiff
    call courant( nx,ny,nz, dt, dta, nstep_a, dtd, nstep_d )
#else
    call courant( nx,ny,nz, dt, dta, nstep_a )
#endif

    ! safety check:
    if ( nstep_a /= 1 ) then
      write (gol,'("something wrong with operator splitting step,")'); call goErr
      write (gol,'("more than 1 sub-step required for advection : ",i6)') nstep_a; call goErr
      TRACEBACK; status=1; return
    end if

#ifdef with_labeling
    call SA_Advec_Setup(status) 
    IF_NOTOK_RETURN(status=1)    
#endif    
    !$OMP parallel &
#ifndef __GFORTRAN__
    !$OMP   default ( none ) &
    !$OMP   shared ( ispecs_advected ) &
#endif
    !$OMP   shared ( nx, ny, nz ) &
    !$OMP   shared ( nstep_a, dta ) &
    !$OMP   shared ( inv_volume ) &
    !$OMP   shared ( bc_west, bc_east, bc_south, bc_north, caloft ) &
    !$OMP   shared ( afluxx, afluxy, afluxz ) &
    !$OMP   shared ( c ) &
    !$OMP   shared ( status_par ) &
    !$OMP   shared ( gol ) &
    !$OMP   private( k, istep ) &
    !$OMP   private( itr ) &
    !$OMP   private( ch ) &
    !$OMP   private( dens0, dens1 ) &
    !$OMP   private( inv_d2 ) &
    !$OMP   private( ix,iy,iz ) &
    !$OMP   private( status ) 

    allocate( ch(1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh) )
    allocate( dens0(nx, ny, nz) )
    allocate( dens1(nx, ny, nz) )
    allocate( inv_d2(nx, ny, nz) )
    
    ! loop over advected tracers:
    !$OMP do
    do itr = 1, n_advected
      ! current index:
      k = ispecs_advected(itr)

      ! copy the concentration vector into the help array ch
      !ch = 0.0
      ch(1:nx,1:ny,1:nz) = c(1:nx,1:ny,1:nz,k)

      ! put bc's in dummy cells of stage vector and help conc. array
      ! west
      ch( 0,1:ny,1:nz)       = bc_west(1:ny,1:nz,k)
      ch(-1,1:ny,1:nz)       = bc_west(1:ny,1:nz,k)
      ! east
      ch(nx+1,1:ny,1:nz)     = bc_east(1:ny,1:nz,k)
      ch(nx+2,1:ny,1:nz)     = bc_east(1:ny,1:nz,k)
      ! south
      ch(1:nx,-1,1:nz)       = bc_south(1:nx,1:nz,k)
      ch(1:nx, 0,1:nz)       = bc_south(1:nx,1:nz,k)
      ! north
      ch(1:nx,ny+1,1:nz)     = bc_north(1:nx,1:nz,k)
      ch(1:nx,ny+2,1:nz)     = bc_north(1:nx,1:nz,k)
      ! upper
      ch(1:nx,1:ny,nz+1)  = caloft(:,:,k)
      ch(1:nx,1:ny,nz+2)  = caloft(:,:,k)
      ! lower boundary conditions are in the loop with the number of time steps!

      ! do the necessary number of time steps
      do istep = 1, nstep_a

        ! do the x-direction

        do iz = 1, nz
           do iy = 1, ny
              do ix = 1, nx

                 ! set initial 'density' (fraction of air density):
                 dens0(ix,iy,iz) = d0

                 ! update 'density' ; aflux* in (m3 air)/s through interface:
                 dens1(ix,iy,iz) = dens0(ix,iy,iz) + dta*(afluxx(ix-1,iy,iz) &
                                   - afluxx(ix,iy,iz)) * inv_volume(ix,iy,iz)
                 ! pre-compute inverse:
                 inv_d2(ix,iy,iz) = 1.0 / dens1(ix,iy,iz)

              enddo
           enddo
        enddo

        ! perform advection
        !call Timer_Start( itim_adv_x )
        call advecX(nx,ny,nz, dens0, dens1, ch, inv_d2,dta,k, status )
#ifdef _OPENMP        
        if (status /= 0 ) then  
          write (gol, '(" Error in Advection in x-direction")' ); call goErr
          status_par = status_par +1
        end if
#else
        IF_NOTOK_RETURN(status=0)
#endif                
        !call Timer_End( itim_adv_x )

        ! do the y-direction

        do iz = 1, nz
           do iy = 1, ny
              do ix = 1, nx

                 ! update density:
                 dens0(ix,iy,iz) = dens1(ix,iy,iz) + dta*(afluxy(ix,iy-1,iz) &
                                   - afluxy(ix,iy,iz)) * inv_volume(ix,iy,iz)
                 ! pre-compute inverse:
                 inv_d2(ix,iy,iz) = 1.0 / dens0(ix,iy,iz)
              enddo
           enddo
        enddo

        ! perform advection
        !call Timer_Start( itim_adv_y )
        call advecY( nx,ny,nz,dens1, dens0, ch, inv_d2, dta, k, status )
#ifdef _OPENMP        
        if (status /= 0 ) then  
          write (gol, '(" Error in Advection in y-direction")' ); call goErr
          status_par = status_par +1
        end if       
#else
        IF_NOTOK_RETURN(status=0)
#endif                
        !call Timer_End( itim_adv_y )

        ! do the z-direction

        ! set the lower boundary conditions here. They are kind of unrealistic, 
        ! but the scheme uses boundary cells to determine minimum and maximum 
        ! concentrations at the cell faces.
        ! Therefore we give the actual values of the lowest layers
        ch(:,:,0)           = ch(:,:,1)
        ch(:,:,-1)          = ch(:,:,1)

        do iz = 1, nz
           do iy = 1, ny
              do ix = 1, nx
                !
                !update density
                !
                 dens1(ix,iy,iz) = dens0(ix,iy,iz) + dta*(afluxz(ix,iy,iz-1) &
                                   - afluxz(ix,iy,iz)) * inv_volume(ix,iy,iz)
                !
                !pre-compute inverse:
                !
                 inv_d2(ix,iy,iz) = 1.0 / dens1(ix,iy,iz)
              enddo
           enddo
        enddo

        ! perform advection
        !call Timer_Start( itim_adv_z )
        call advecZ( nx,ny,nz,dens0, dens1, ch, inv_d2, dta, k, status )
#ifdef _OPENMP        
        if (status /= 0 ) then  
          write (gol, '(" Error in Advection in z-direction")' ); call goErr
          status_par = status_par +1
        end if       
#else
        IF_NOTOK_RETURN(status=0)
#endif                
        !call Timer_End( itim_adv_z )

        ! end loop over #advection steps
      end do

      ! put ch into c, the global concentration array
      c(1:nx,1:ny,1:nz,k) = ch(1:nx,1:ny,1:nz)

    end do    ! end loop over species
    !$OMP end do
    
    
    ! clear:
    deallocate( ch )
    deallocate( dens0 )
    deallocate( dens1 )
    deallocate( inv_d2 )

    !$OMP end parallel
#ifdef _OPENMP    
    if (status_par /= 0 ) then
      write(gol, '("Error in Advection routine")'); call goErr
      TRACEBACK; status=1;return
    endif
#endif

#ifdef with_hdiff
    ! now do the horizontal diffusion
    !call Timer_Start( itim_adv_dif )
    !$OMP parallel &
    !$OMP   default( none ) &
    !$OMP   private( ispec, istep )
    !$OMP   do
    ! loop over advected tracers:
    do itr = 1, n_advected
      ! current index:
      ispec = ispecs_advected(itr)
      do istep = 1, nstep_d
        ! do horizontal diffusion for the compound
        call hordif( c(1:nx,1:ny,1:nz,ispec), dtd, &
                     bc_west (1:ny,1:nz,ispec), bc_east(1:ny,1:nz,ispec), &
                     bc_south(1:nx,1:nz,ispec),bc_north(1:nx,1:nz,ispec) )
      end do  ! diffusion steps
    end do  ! specs
    !$OMP   end do
    !$OMP end parallel
    !call Timer_End( itim_adv_dif )
#endif
    
    ! clear:
    deallocate( inv_volume )

    ! ok
    status = 0

  end subroutine LE_Advec_Apply


  ! ***
  
  
  subroutine advecX( nx,ny,nz,d1, d2, ch, inv_d2,dta, k, status )                           

    use LE_Meteo_Data, only : afluxx, volume
#ifdef with_labeling
    use SA_Labeling, only : SA_Advec_x
#endif        
    
    ! --- in/out ---------------------------------
        
    integer, parameter  ::  nh = 2
    integer, intent(in)   ::  nx, ny, nz    
    real, intent(in)    ::   d1(nx,ny,nz), d2(nx,ny,nz)
    real, intent(in)    ::   dta
    real, intent(inout) ::   ch(1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in)    ::   inv_d2(nx,ny,nz)        
    integer, intent(out)    ::  status        
    integer, intent(in)     ::  k
    
    ! --- const ----------------------------------
    character(len=*), parameter :: rname = mname//'/advecX'
    
    ! --- local ----------------------------------
    
    logical, dimension(:,:,:), allocatable :: is_extreme
    real, dimension(:,:,:), allocatable    :: cgues
    real, dimension(:,:,:), allocatable    :: fluxx

    integer :: ix, iy, iz
    real :: c1, c2, u, cour, fcour, alpha, cg
    
    ! --- begin ----------------------------------

    ! storage:
    allocate ( is_extreme(-1:nx+1,ny,nz) )
    allocate ( cgues(nx,ny,nz) )
    allocate ( fluxx(0:nx,ny,nz) )

    ! initialize flux
    fluxx = 0.0

    ! determine the extremes
    is_extreme = .false.
    ! NOTE: attemps to openmp these loops slowed down the whole routine ...
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      c1 = ch(ix-1,iy,iz)
      c2 = ch(ix+1,iy,iz)
      if (ch(ix,iy,iz) >= max(c1, c2) .OR. ch(ix,iy,iz) <= min(c1,c2)) is_extreme(ix,iy,iz) = .true.
    enddo
    enddo
    enddo

    !
    ! Loop over cell interfaces from left to right:
    ! fill all Qf (tracer vmr in flux through interface)
    ! assuming linear distribution  (here: fluxx)
    !
    ! Loop over cells from left to right:
    ! for all inflow-only or outflow-only cells,
    ! replace Qf at both sides given vmr and air-mass-flux, new vmr.  
    !
    ! Loop over cells from left to right:
    ! for rightflow-only cells, compute new vmr given previously computed fluxes;
    ! limit mixing ratio, if necessary re-compute right flux given left flux
    ! and mass change.
    !
    ! Loop over cells from right to left:
    ! for left-only cells (bug: test for left outflow) compute new vmr;
    ! limit mixing ratio, if necessary re-compute left flux given right flux
    ! and mass change.
    !
    !                <--     <--     -->     -->     -->     <--     <--     -->  
    !                 | left  |outflow| right | right |inflow | left  |outflow|
    !              ---+-------+-------+-------+-------+-------+-------+-------+----
    !
    ! 0,nx           F,a     F,a     a,F     a,F     a,F     F,a     F,a     a,F
    !
    ! 1,nx                    F       F               F       F       F       F
    !                             Q                       Q               Q
    !
    ! 1,nx                                Q > F
    !                                             Q > F                             a)
    !
    ! nx,1,-1                                                        (F < Q    )    b)
    !                                                         F < Q
    !                        (F < Q     )                                           b)
    !                 F < Q
    !
    ! a) BUG: limiting of Q forces computation of a new flux F,
    !    but this doesn't seem to lead to an update of Q in the right cell
    ! 
    ! b) Test wrongly implemented ?
    !          
    !----------
    !
    !  up loop:                   .   F
    !                                     Q > F
    !                                             Q > F
    !                                                                     .   F
    !
    !  down loop: 
    !                                                                 F < Q
    !                                                         F < Q
    !                                                     Q                      
    !                         F < Q
    !                 F < Q
    !


    ! compute the fluxes: here still as a concentration at the cell face
    do iz=1,nz
    do iy=1,ny
    ! loop over cell interfaces:
    do ix=0,nx
      ! use volume flux through this interface for 'wind':
      u = afluxx(ix,iy,iz)
      ! courant number using volume left from the interface:
      !AJS: bug? should be 'upwind' courrant number with different volume given sign:
      !AJS:   cour = dta *     u(ix)  / volume(ix  )   , u >= 0.0
      !AJS:   cour = dta * abs(u(ix)) / volume(ix+1)   , u <  0.0
      cour = dta*abs(u) * inv_volume(max(1,ix),iy,iz)

      !
      !                     -->  
      !                   o     o
      !              o
      !
      !                     <--       *
      !                   *     *
      !    |     |     |     |     |     |
      !    +-----+-----+-----+-----+-----+
      !      i-2   i-1    i    i+1   i+2
      !               i-1    i
      !               imh
      !
      ! Average flux through interface during timestep.
      ! Assume linear tracer distribution within downwind cell
      ! with gradient computed using central difference
      ! and extra slope factor alfa:
      !
      !                Q(i+1) - Q(i-1)
      !    gradient =  --------------- alfa     , u_i >= 0
      !                     2 dx
      !
      !                Q(i+2) - Q(i  )
      !    gradient =  --------------- alfa     , u_i <  0
      !                     2 dx
      !
      ! which defines the mixing ratio within the cell:
      !
      !    Q(x) = Q(x_i) + gradient * (x-x_i)
      !
      ! The 'upwind' courant number is the faction of the cell flowing out.
      ! Given its linear distribution within the cell, the average
      ! concentration in the interval flowing out is the concentration
      ! halfway the interval:
      !
      !    xf  =  x_{i}  +dx/2 - c*dx/2  =  x_{i}   + (1-c)*dx/2   ,  u_i >= 0
      !        =  x_{i+1}-dx/2 + c*dx/2  =  x_{i+1} - (1-c)*dx/2   ,  u_i <  0
      !
      ! with concentration:
      !
      !    Q(xf) = Q(x_i) + gradient * (1-c)*dx/2
      !          = Q(x_i) - gradient * (1-c)*dx/2
      !     
      
      ! average concentration ((m3 tracer)/(m3 air) or (ug tracer)/(m3 air))
      fcour= (1.0-cour)/4.0
      alpha = 1.0
      if (u >=0.0) then
        if (is_extreme(ix-1,iy,iz)) alpha = max(1.5,1.2+0.6*cour)
        if (is_extreme(ix+1,iy,iz)) alpha = (1.75-0.45*cour)
        fluxx(ix,iy,iz) = ch(ix,iy,iz) + (ch(ix+1,iy,iz) - ch(ix-1,iy,iz))*fcour*alpha
      else
        if (is_extreme(ix+1,iy,iz)) alpha = max(1.5,1.2+0.6*cour)
        if (is_extreme(ix-1,iy,iz)) alpha = (1.75-0.45*cour)
        fluxx(ix,iy,iz) = ch(ix+1,iy,iz) + (ch(ix,iy,iz) - ch(ix+2,iy,iz))*fcour*alpha
      endif

      ! bound cell face concentration
      !AJS: BUG? Flux through interface should be within minimum/maximum
      !AJS: BUG? of three cells if cell has inflow only ..
      !AJS: if ( (u(imh) > 0.0) .and. (u(imh+1) < 0.0) ) then
      !AJS:   c1 = min(        min(ch(i-1),ch(i)),ch(i))
      !AJS:   c2 = max(        max(ch(i-1),ch(i)),ch(i))
      !AJS: else if ( (u(imh-1) > 0.0) .and. (u(imh) < 0.0) ) then
      !AJS:   c1 = min(ch(i-2),min(ch(i-1)),ch(i)      )
      !AJS:   c2 = max(ch(i-2),max(ch(i-1)),ch(i)      )
      !AJS: else
      !AJS:   c1 =             max(ch(i-1),ch(i)       )
      !AJS:   c2 =             max(ch(i-1),ch(i)       )
      !AJS: end if
      c1 = min( ch(ix,iy,iz), ch(ix+1,iy,iz) )
      c2 = max( ch(ix,iy,iz), ch(ix+1,iy,iz) )
      fluxx(ix,iy,iz) = max(c1, min( fluxx(ix,iy,iz), c2) )
      ! make a flux out of it ....
      fluxx(ix,iy,iz) = u*dta*fluxx(ix,iy,iz)*d0

    enddo
    enddo

    enddo

    ! now the exceptions: inflow or outflow only cells
    do iz=1,nz

    do iy=1,ny
    do ix=1,nx
      if (afluxx(ix,iy,iz) <= 0.0 .AND. afluxx(ix-1,iy,iz) >= 0.0) then 
        ! inflow only
        fluxx(ix-1,iy,iz) = dta*afluxx(ix-1,iy,iz)*ch(ix-1,iy,iz)*d0
        fluxx(ix  ,iy,iz) = dta*afluxx(ix  ,iy,iz)*ch(ix+1,iy,iz)*d0
        cgues(ix,iy,iz) = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
               (fluxx(ix-1,iy,iz) - fluxx(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
      else if (afluxx(ix,iy,iz) >= 0.0 .AND. afluxx(ix-1,iy,iz) <= 0.0) then 
        ! outflow only
        fluxx(ix-1,iy,iz) = dta*afluxx(ix-1,iy,iz)*ch(ix,iy,iz)*d0
        fluxx(ix  ,iy,iz) = dta*afluxx(ix  ,iy,iz)*ch(ix,iy,iz)*d0
        cgues(ix,iy,iz) = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
               (fluxx(ix-1,iy,iz) - fluxx(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
      endif
    enddo
    enddo

    enddo


    ! intermediate solution
    !cgues(1:nx,1:ny,1:nz) = (ch(1:nx,1:ny,1:nz)*d1(1:nx,1:ny,1:nz) + &
    !                        (fluxx(0:nx-1,1:ny,1:nz) - fluxx(1:nx,1:ny,1:nz))/volume(1:nx,1:ny,1:nz) )/ &
    !                        d2(1:nx,1:ny,1:nz)

    ! and now the limiting of the fluxes, first going up
    do iz=1,nz

    do iy=1,ny
    do ix=1,nx
       ! look for positive winds
       if (afluxx(ix,iy,iz) > 0.0 .AND. afluxx(ix-1,iy,iz) >=0.0) then
         cg = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
               (fluxx(ix-1,iy,iz) - fluxx(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
             c1 = min( ch(ix-1,iy,iz), ch(ix,iy,iz) )
             c2 = max( ch(ix-1,iy,iz), ch(ix,iy,iz) )
         if (cg < c1 .OR. cg > c2) then
             cg = min(c2, max(cg, c1) )
             fluxx(ix,iy,iz) = fluxx(ix-1,iy,iz) + &
                               (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz))*volume(ix,iy,iz)
         endif
         cgues(ix,iy,iz) = cg
       endif
    enddo
    enddo

    enddo

    ! now the limiting going down
    do iz=1,nz

    do iy=1,ny
    do ix=nx,1,-1
       ! look for negative winds
       !BUG? if (afluxx(ix-1,iy,iz) < 0.0) then
       if ( (afluxx(ix,iy,iz) < 0.0) .and. (afluxx(ix-1,iy,iz) < 0.0) ) then
         cg = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
               (fluxx(ix-1,iy,iz) - fluxx(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
         c1 = min( ch(ix,iy,iz), ch(ix+1,iy,iz) )
         c2 = max( ch(ix,iy,iz), ch(ix+1,iy,iz) )
         if (cg < c1 .OR. cg > c2) then
             cg = min(c2, max(cg, c1) )
             fluxx(ix-1,iy,iz) = fluxx(ix,iy,iz) - &
                               (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz))*volume(ix,iy,iz)
         endif
         cgues(ix,iy,iz) = cg
       endif
    enddo
    enddo

    enddo

#ifdef with_labeling    
    call SA_Advec_x( ch(1:nx,1:ny,1:nz), fluxx, k, status)
    IF_NOTOK_RETURN(status=1)    
#endif        

    ! copy cgues to c
    ch(1:nx,1:ny,1:nz) = cgues(1:nx,1:ny,1:nz) 

    ! clear:
    deallocate ( is_extreme )
    deallocate ( cgues )
    deallocate ( fluxx )
    
    ! ok
    status = 0
    
  end subroutine advecX


  ! ***----------------------------------------------
  
  
  subroutine advecY( nx,ny,nz,d1, d2, ch, inv_d2,dta,k, status)

    use LE_Meteo_Data, only : afluxy, volume
#ifdef with_labeling
    use SA_Labeling, only : SA_Advec_y
#endif        

    ! --- in/out ---------------------------------

    integer, parameter  ::  nh = 2
    
    integer, intent(in)   ::  nx, ny, nz    
    real, intent(in)     ::  d1(nx,ny,nz), d2(nx,ny,nz)
    real, intent(in)     ::  dta
    real, intent(inout) ::   ch(1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in)    ::   inv_d2(nx,ny,nz)
    integer, intent(in)   ::  k
    integer, intent(out)  ::  status  
    
    ! --- const ----------------------------------
    character(len=*), parameter :: rname = mname//'/advecY'

    ! --- local ----------------------------------

    logical, dimension(:,:,:), allocatable :: is_extreme
    real, dimension(:,:,:), allocatable    :: cgues
    real, dimension(:,:,:), allocatable    :: fluxy

    integer :: ix, iy, iz
    real :: c1, c2, u, cour, fcour, alpha, cg
    
    ! --- begin ----------------------------------

    ! storage:
    allocate ( is_extreme(nx,-1:ny+1,nz) )
    allocate ( cgues(nx,ny,nz) )
    allocate ( fluxy(nx,0:ny,nz) )

    ! initialize flux
    fluxy = 0.0
 
    ! determine the extremes
    is_extreme = .false.
    do iz=1,nz

    do iy=1,ny
    do ix=1,nx
      c1 = ch(ix,iy-1,iz)
      c2 = ch(ix,iy+1,iz)
      if (ch(ix,iy,iz) >= max(c1, c2) .OR. ch(ix,iy,iz) <= min(c1,c2)) is_extreme(ix,iy,iz) = .true.
    enddo
    enddo

    enddo

    ! compute the fluxes: here still as a concentration at the cell face
    do iz=1,nz

    do iy=0,ny
    do ix=1,nx
      u = afluxy(ix,iy,iz)
      cour = dta*abs(u) * inv_volume(ix,max(1,iy),iz)
      fcour= (1.0-cour)/4.0
      alpha = 1.0
      if (u >=0.0) then
        if (is_extreme(ix,iy-1,iz)) alpha = max(1.5,1.2+0.6*cour)
        if (is_extreme(ix,iy+1,iz)) alpha = (1.75-0.45*cour)
        fluxy(ix,iy,iz) = ch(ix,iy,iz) + (ch(ix,iy+1,iz) - ch(ix,iy-1,iz))*fcour*alpha
      else
        if (is_extreme(ix,iy+1,iz)) alpha = max(1.5,1.2+0.6*cour)
        if (is_extreme(ix,iy-1,iz)) alpha = (1.75-0.45*cour)
        fluxy(ix,iy,iz) = ch(ix,iy+1,iz) + (ch(ix,iy,iz) - ch(ix,iy+2,iz))*fcour*alpha
      endif
      ! bound cell face concentration
      c1 = min( ch(ix,iy,iz), ch(ix,iy+1,iz) ) 
      c2 = max( ch(ix,iy,iz), ch(ix,iy+1,iz) ) 
      fluxy(ix,iy,iz) = max(c1, min( fluxy(ix,iy,iz), c2) )
      ! make a flux out of it ....
      fluxy(ix,iy,iz) = u*dta*fluxy(ix,iy,iz)*d0

    enddo
    enddo

    enddo

    ! now the exceptions: inflow or outflow only cells
    do iz=1,nz

    do iy=1,ny
    do ix=1,nx
      if (afluxy(ix,iy,iz) <= 0.0 .AND. afluxy(ix,iy-1,iz) >= 0.0) then 
        ! inflow only
        fluxy(ix,iy-1,iz) = dta*afluxy(ix,iy-1,iz)*ch(ix,iy-1,iz)*d0
        fluxy(ix,iy  ,iz) = dta*afluxy(ix,iy  ,iz)*ch(ix,iy+1,iz)*d0
        cgues(ix,iy,iz)   = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
                     (fluxy(ix,iy-1,iz) - fluxy(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
      else if (afluxy(ix,iy,iz) >= 0.0 .AND. afluxy(ix,iy-1,iz) <= 0.0) then 
        ! outflow only
        fluxy(ix,iy-1,iz) = dta*afluxy(ix,iy-1,iz)*ch(ix,iy,iz)*d0
        fluxy(ix,iy  ,iz) = dta*afluxy(ix,iy,iz)*ch(ix,iy,iz)*d0
        cgues(ix,iy,iz)   = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
                     (fluxy(ix,iy-1,iz) - fluxy(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
      endif
    enddo
    enddo

    enddo

    ! and now the limiting of the fluxes, first going up
    do iz=1,nz

    do iy=1,ny
    do ix=1,nx
       ! look for positive winds
       if (afluxy(ix,iy,iz) > 0.0 .AND. afluxy(ix,iy-1,iz) >=0.0) then
         cg = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
              (fluxy(ix,iy-1,iz) - fluxy(ix,iy,iz))/volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
         c1 = min( ch(ix,iy-1,iz), ch(ix,iy,iz) )
         c2 = max( ch(ix,iy-1,iz), ch(ix,iy,iz) )
         if (cg < c1 .OR. cg > c2) then
             cg = min(c2, max(cg, c1) )
             fluxy(ix,iy,iz) = fluxy(ix,iy-1,iz) + &
                               (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz))*volume(ix,iy,iz)
         endif
         cgues(ix,iy,iz) = cg
       endif
    enddo
    enddo

    enddo

    ! now the limiting going down
    do iz=1,nz

    do iy=ny,1,-1
    do ix=1,nx
       ! look for negative winds
       !BUG? if (afluxy(ix,iy-1,iz) < 0.0) then
       if ( (afluxy(ix,iy,iz) < 0.0) .and. (afluxy(ix,iy-1,iz) < 0.0) ) then
         cg = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
              (fluxy(ix,iy-1,iz) - fluxy(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
         c1 = min( ch(ix,iy,iz), ch(ix,iy+1,iz) )
         c2 = max( ch(ix,iy,iz), ch(ix,iy+1,iz) )
         if (cg < c1 .OR. cg > c2) then
             cg = min(c2, max(cg, c1) )
             fluxy(ix,iy-1,iz) = fluxy(ix,iy,iz) - &
                               (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz))*volume(ix,iy,iz)
         endif
         cgues(ix,iy,iz) = cg
       endif
    enddo
    enddo

    enddo

#ifdef with_labeling
    call SA_Advec_y( ch(1:nx,1:ny,1:nz), fluxy, k, status )
    IF_NOTOK_RETURN(status=1)
#endif
    
    ! copy cgues to c
    ch(1:nx,1:ny,1:nz) = cgues(1:nx,1:ny,1:nz) 

    ! clear:
    deallocate ( is_extreme )
    deallocate ( cgues )
    deallocate ( fluxy )

    ! ok
    status = 0
    
  end subroutine advecY


  ! ***------------------------------------------------------
  
  
  subroutine advecZ( nx, ny, nz, d1, d2, ch, inv_d2, dta, k, status)

    use LE_Meteo_Data, only : afluxz, volume
#ifdef with_labeling
    use SA_Labeling, only : SA_Advec_z
#endif     

    ! --- in/out ---------------------------------
    
    integer, parameter  ::  nh = 2

    integer, intent(in)   ::  nx, ny, nz    
    real, intent(in)      :: d1(nx,ny,nz), d2(nx,ny,nz)
    real, intent(in)      ::  dta
    real, intent(inout)   ::   ch(1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in)      ::   inv_d2(nx,ny,nz)
    integer, intent(in)     ::  k
    integer, intent(out)    ::  status  

    ! --- const ----------------------------------
    character(len=*), parameter :: rname = mname//'/advecZ'
    
    ! --- local ----------------------------------
    
    real, dimension(:,:,:), allocatable    :: cgues
    real, dimension(:,:,:), allocatable    :: fluxz

    integer :: ix, iy, iz
    real :: c1, c2, u, cour, fcour, alpha, cg

    ! --- begin ----------------------------------
    
    ! storage:
    allocate ( cgues(nx,ny,nz) )
    allocate ( fluxz(nx,ny,0:nz) )

    ! initialize flux
    fluxz = 0.0

    ! personal communication: put alpha=5 for vertical to reduce diffusion
    alpha = 5.0


    ! compute the fluxes: here still as a concentration at the cell face
    do iz=0,nz

    do iy=1,ny
    do ix=1,nx
      u = afluxz(ix,iy,iz)
      cour = dta*abs(u) * inv_volume(ix,iy,max(1,iz))
      fcour= (1.0-cour)/4.0
      if (u >=0.0) then
        fluxz(ix,iy,iz) = ch(ix,iy,iz) + (ch(ix,iy,iz+1) - ch(ix,iy,iz-1))*fcour*alpha
      else
        fluxz(ix,iy,iz) = ch(ix,iy,iz+1) + (ch(ix,iy,iz) - ch(ix,iy,iz+2))*fcour*alpha
      endif
      ! bound cell face concentration
      c1 = min( ch(ix,iy,iz), ch(ix,iy,iz+1) ) 
      c2 = max( ch(ix,iy,iz), ch(ix,iy,iz+1) ) 
      fluxz(ix,iy,iz) = max(c1, min( fluxz(ix,iy,iz), c2) )
      ! make a flux out of it ....
      fluxz(ix,iy,iz) = u*dta*fluxz(ix,iy,iz)*d0

    enddo
    enddo

    enddo

    ! now the exceptions: inflow or outflow only cells
    do iz=1,nz

    do iy=1,ny
    do ix=1,nx
        if (afluxz(ix,iy,iz) <= 0.0 .AND. afluxz(ix,iy,iz-1) >= 0.0) then 
          ! inflow only
          fluxz(ix,iy,iz-1) = dta*afluxz(ix,iy,iz-1)*ch(ix,iy,iz-1)*d0
          fluxz(ix,iy,iz  ) = dta*afluxz(ix,iy,iz  )*ch(ix,iy,iz+1)*d0
          cgues(ix,iy,iz)   = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
                            (fluxz(ix,iy,iz-1) - fluxz(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
        else if (afluxz(ix,iy,iz) >= 0.0 .AND. afluxz(ix,iy,iz-1) <= 0.0) then 
          ! outflow only
          fluxz(ix,iy,iz-1) = dta*afluxz(ix,iy,iz-1)*ch(ix,iy,iz)*d0
          fluxz(ix,iy,iz  ) = dta*afluxz(ix,iy,iz)*ch(ix,iy,iz)*d0
          cgues(ix,iy,iz)   = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
                            (fluxz(ix,iy,iz-1) - fluxz(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
        endif
      enddo  ! x
      enddo  ! y

      ! upper boundary
      if ( iz == nz ) then
        ! put upstream at the boundaries (because of u(0)=0.0 no need for the flux at the bottom here)

        do iy=1,ny
        do ix=1,nx
          if (afluxz(ix,iy,nz) <= 0.0) then
             fluxz(ix,iy,nz) = dta*afluxz(ix,iy,nz)*ch(ix,iy,nz+1)*d0
          endif
        enddo  ! x
        enddo  ! y

      end if
    enddo  ! z


    ! and now the limiting of the fluxes, first going up
    do iz=1,nz

    do iy=1,ny
    do ix=1,nx
       ! look for positive winds
       if (afluxz(ix,iy,iz) > 0.0 .AND. afluxz(ix,iy,iz-1) >=0.0) then
         cg = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
              (fluxz(ix,iy,iz-1) - fluxz(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
         c1 = min( ch(ix,iy,iz-1), ch(ix,iy,iz) )
         c2 = max( ch(ix,iy,iz-1), ch(ix,iy,iz) )
         if (cg < c1 .OR. cg > c2) then
             cg = min(c2, max(cg, c1) )
             fluxz(ix,iy,iz) = fluxz(ix,iy,iz-1) + &
                               (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz))*volume(ix,iy,iz)
         endif
         cgues(ix,iy,iz) = cg
       endif
    enddo
    enddo

    enddo

    ! now the limiting going down
    do iz=nz,1,-1

    do iy=1,ny
    do ix=1,nx
       ! look for negative winds
       !BUG? if (afluxz(ix,iy,iz-1) < 0.0) then
       if ( (afluxz(ix,iy,iz) < 0.0) .and. (afluxz(ix,iy,iz-1) < 0.0) ) then
         cg = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
              (fluxz(ix,iy,iz-1) - fluxz(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)
         c1 = min( ch(ix,iy,iz), ch(ix,iy,iz+1) )
         c2 = max( ch(ix,iy,iz), ch(ix,iy,iz+1) )
         if (cg < c1 .OR. cg > c2) then
             cg = min(c2, max(cg, c1) )
             fluxz(ix,iy,iz-1) = fluxz(ix,iy,iz) - &
                               (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz))*volume(ix,iy,iz)
             if (iz==1) stop 'unacceptable flux adjustment'
         endif
         cgues(ix,iy,iz) = cg
       endif
    enddo
    enddo

    enddo

#ifdef with_labeling 
    call SA_Advec_z( ch(1:nx,1:ny,1:nz), fluxz, k, status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! copy cgues to c
    ch(1:nx,1:ny,1:nz) = cgues(1:nx,1:ny,1:nz) 

    ! clear:
    deallocate ( cgues )
    deallocate ( fluxz )

    ! ok
    status = 0
    
  end subroutine advecZ


  ! ***----------------------------------------------------------


  ! check courant condition for dt. If necessary,
  ! compute how many substeps are needed to satisfy
  ! the courant condition

#ifdef with_hdiff
  subroutine courant( nx, ny,nz, dt, dta, nstep_a, dtd, nstep_d, store )
#else
  subroutine courant( nx, ny,nz, dt, dta, nstep_a, store )
#endif

    use LE_Meteo_Data, only : afluxx, afluxy, afluxz
    use LE_Meteo_Data, only : cflx, cfly, cflz  ! (nx,ny,nz)
#ifdef with_hdiff
    use dims, only : khx, khy
#endif
    use dims, only : runF, outF
    use LE_Logging, only : ident2

    ! --- in/out ---------------------------------
    
    integer, intent(in)   ::  nx, ny, nz
    real, intent(in)      ::  dt
    
    ! actual time steps and number of steps for advection:
    real, intent(out)     ::  dta
    integer, intent(out)  ::  nstep_a

#ifdef with_hdiff
    ! actual time steps and number of steps for diffusion:
    integer, intent(out)  ::  nstep_d
    real, intent(out)     ::  dtd
#endif

    logical, intent(in), optional   ::  store

    ! --- local ----------------------------------
    
    real    :: cour_max, courx_max, coury_max, courz_max, courx, coury,courz
    integer :: i,j,k
    logical :: do_store

    ! --- begin ----------------------------------
    
    ! store cfl numbers?
    do_store = .false.
    if ( present(store) ) do_store = store
    
    cour_max = 0.0
    courx_max = 0.0
    coury_max = 0.0
    courz_max = 0.0

    do i=1,nx
    do j=1,ny
    do k=1,nz
      ! for the Walcek scheme, each direction is separate !
      courx = dt*max( abs( afluxx(i-1,j,k) ), abs( afluxx(i,j,k) ) ) * inv_volume(i,j,k)
      coury = dt*max( abs( afluxy(i,j-1,k) ), abs( afluxy(i,j,k) ) ) * inv_volume(i,j,k)
      courz = dt*max( abs( afluxz(i,j,k-1) ), abs( afluxz(i,j,k) ) ) * inv_volume(i,j,k)
      ! store?
      if ( do_store ) then
        cflx(i,j,k) = courx
        cfly(i,j,k) = coury
        cflz(i,j,k) = courz
      end if
      ! update maximum:
      cour_max = max(cour_max, courx, coury, courz)
      ! the 1D courant numbers
      courx_max = max(courx_max, courx)
      coury_max = max(coury_max, coury)
      courz_max = max(courz_max, courz)
    enddo
    enddo
    enddo

    if (cour_max > courant_max) then
      nstep_a = int( cour_max / courant_max) + 1
      dta = dt/nstep_a
    else
      dta = dt
      nstep_a = 1
    endif

    !if (.NOT.outF%suppress) then
    !if ( do_store ) then
    !   write (gol,*) ident2, '  >>> #advection steps:   ', nstep_a; call goPr
    !   write (gol,*) ident2, '  >>> courant number was: ', cour_max; call goPr
    !   write (gol,*) ident2, '  >>> courant number xyz: ', courx_max, coury_max, courz_max; call goPr
    !end if
  
#ifdef with_hdiff
    ! check time step limit for diffusion
    cour_max = 0.0

    do i=0,nx
    do j=1,ny
    do k=1,nz
       cour = abs( khx(i,j,k)*dt)/runF%dx(j) * inv_volume(max(1,i),j,k)
       cour_max = max(cour_max, cour)
    enddo
    enddo
    enddo

    do i=1,nx
    do j=0,ny
    do k=1,nz
       cour = abs( khy(i,j,k)*dt)/runF%dy * inv_volume(i,max(1,j),k)
       cour_max = max(cour_max, cour)
    enddo
    enddo
    enddo

    nstep_d = 1
    dtd     = dt
    if (cour_max > 0.25) then
      ! stability for hor diffusion violated
      nstep_d = int (cour_max/0.25) + 1
      dtd = dt/nstep_d
    endif

    if (nstep_d > nstep_a .AND. (.NOT.outF%suppress) ) then
      write (*,*) ident2, '  >>> diffusion dominates advection'
      write (*,*) ident2, '  >>> number of steps required for hor diffusion: ', nstep_d
    endif
#endif

  end subroutine courant


!  ! ===================================================================
!  ! ===
!  ! === Walcek scheme following original paper.
!  ! ===
!  ! === This implementation was never tested yet,
!  ! === stored here for future inspiration.
!  ! ===
!  ! ===================================================================
!  
!  
!  ! 
!  ! This subroutine calculates change in mass-mixing-ratio Q during time
!  ! step dt due to advection along a grid with n cells.
!  ! Air mass fluxes phi are specified at cell faces 0:n .
!  ! Due to the mass fluxes, the air mass m in the cells are changed.
!  !
!  !                     -->  
!  !                   o     o
!  !              o
!  !
!  !                     <--       *
!  !                   *     *
!  !    |     |     |     |     |     |
!  !    +-----+-----+-----+-----+-----+
!  !      i-2   i-1    i    i+1   i+2
!  !               i-1    i
!  !
!  ! Define 1D density:
!  ! 
!  !   rho_i = m_i / dx_i     [(kg air)/m]
!  !
!  ! The 'upwind' courant number is the faction of the cell flowing out:
!  ! 
!  !   cour = dt *  phi_{i}    / m_i   [0-1]   , phi_{i}   >= 0
!  !        = dt * |phi_{i-1}| / m_i   [0-1]   , phi_{i-1} <  0
!  !
!  ! Average flux through interface during timestep.
!  ! Assume linear tracer distribution within downwind cell
!  ! with gradient computed using central difference
!  ! and extra slope factor alfa:
!  !
!  !                Q(x_{i+1}) - Q(x_{i-1})
!  !    gradient =  ----------------------- alfa     , phi_i >= 0
!  !                          2 dx
!  !
!  !                Q(x_{i+2}) - Q(x_{i  })
!  !    gradient =  ----------------------- alfa     , phi_i <  0
!  !                          2 dx
!  !
!  ! which defines the mixing ratio within the cell:
!  !
!  !    Q(x) = Q(x_i) + gradient * (x-x_i)
!  !
!  ! Given its linear distribution within the cell, the average
!  ! mixing-ratio in the interval flowing out (with length cour*dx)
!  ! is the mixing-ratio halfway the interval at location:
!  !
!  !    xf  =  x_{i}  +dx/2 - cour*dx/2  =  x_{i}   + (1-cour)*dx/2   ,  phi_i >= 0
!  !        =  x_{i+1}-dx/2 + cour*dx/2  =  x_{i+1} - (1-cour)*dx/2   ,  phi_i <  0
!  !
!  ! The mixing-ratio in xf is:
!  !
!  !    Qf = Q(x_i) + gradient * (1-cour)*dx/2
!  !       = Q(x_i) - gradient * (1-cour)*dx/2
!  ! 
!  ! The distance dx in the gradient canceles the dx at the right,
!  ! such that Qf can be computed without dx:
!  !
!  !    Qf = Q(x_i) + [Q(i+1)-Q(i-1)]*alfa*(1-cour)/4   ,  phi_i >= 0
!  !       = Q(x_i) - [Q(i+2)-Q(i  )]*alfa*(1-cour)/4   ,  phi_i <  0
!  ! 
!  ! NOTES: 
!  !  o In LOTOS-EUROS implementation:
!  !    personel communication: put alpha=5 for vertical to reduce diffusion
!  !
! 
!  subroutine Walcek_1D_paper( n, phi, m_t0, Q_t0, dt, m_t1, Q_t1 )
!  
!    ! --- in/out ------------------------------------------
!    
!    integer, parameter    ::  nh = 1  ! number of halo cells
!    
!    integer, intent(in)   ::  n                ! number of grid cells
!    real, intent(in)      ::  phi(0:n)         ! air mass flux through cell faces  [(kg air)/s]
!    real, intent(in)      ::  m_t0(1-nh:n+nh)  ! initial air mass in (halo) cells  [(kg air)]
!    real, intent(in)      ::  Q_t0(1-nh:n+nh)  ! initial mass-mixing-ratios in (halo) cells [(kg tracer)/(kg air)]
!    real, intent(in)      ::  dt               ! time step  [s]
!    real, intent(out)     ::  m_t1(1-nh:n+nh)  ! final air mass in (halo) cells  [(kg air)]
!    real, intent(out)     ::  Q_t1(1-nh:n+nh)  ! final mass-mixing-ratios in (halo) cells [(kg tracer)/(kg air)]
!    
!    ! --- local -------------------------------------------
!    
!    integer          ::  i
!    logical          ::  extreme(1-nh:n+nh)
!    real             ::  Q_t1_min(1:n)
!    real             ::  Q_t1_max(1:n)
!    real             ::  Flux(0:n)    ! tracer flux [(kg tracer)/m2]
!    real             ::  Qf
!    real             ::  cour
!    real             ::  alfa
!    
!    ! --- begin -------------------------------------------
!    
!    ! local extreme ?
!    extreme(1-nh:1-1) = .false.  ! no extemes in halo cells
!    extreme(n+1:n+nh) = .false.  ! no extemes in halo cells
!    do i = 1, n
!      extreme(i) = ( Q_t0(i) >= max(Q_t0(i-1),Q_t0(i+1)) ) .or. &
!                   ( Q_t0(i) <= min(Q_t0(i-1),Q_t0(i+1)) )
!    end do
!    
!    ! tracer mixing ratio limits at new time:
!    do i = 1, n
!      if ( (u(i-1) >= 0.0) .and. (u(i) < 0.0) ) then
!        ! inflow only:
!        Q_t1_min(i) = min( min( Q_t0(i-1), Q_t0(i) ), Q_t0(i+1) )  ! line after Eq. 7
!        Q_t1_max(i) = max( max( Q_t0(i-1), Q_t0(i) ), Q_t0(i+1) )  ! line after Eq. 7
!      else if ( u(i-1) >= 0.0 ) then
!        ! inflow from left:
!        Q_t1_min(i) = min( Q_t0(i-1), Q_t0(i) )  ! Eq. 7a
!        Q_t1_max(i) = max( Q_t0(i-1), Q_t0(i) )  ! Eq. 7a
!      else if ( u(i) < 0.0 )
!        ! inflow from right:
!        Q_t1_min(i) = min( Q_t0(i), Q_t0(i+1) )  ! Eq. 7b
!        Q_t1_max(i) = max( Q_t0(i), Q_t0(i+1) )  ! Eq. 7b
!      else
!        ! outflow only does not change tracer mixing ratio ...
!        Q_t1_min(i) = Q_t0(i)
!        Q_t1_max(i) = Q_t0(i)
!      end if
!    end do
!    
!    ! air mass at end of time step:
!    do i = 1, n
!      m_t1(i) = m_t0(i) + ( phi(i-1) - phi(i) )*dt  ! (kg air)
!    end do
!    
!    ! Update tracer mixing ratios, taking limits into account.
!    ! First loop from left to right computing:
!    !   o right tracer flux (F) for outflow only cell;
!    !   o new mixing ratio (Q) for cells with flow to the right;
!    !     left inflow and mixing ratio change defines (>) right tracer flux (F)
!    ! Then loop from right to left computing:
!    !   o new mixing ratio (Q) for cells with inflow only;
!    !   o new mixing ratio (Q) for cells with outflow to left;
!    !     right in- or outflow and mixing ratio change defines (>) left tracer flux (F)
!    !
!    !                <--     <--     -->     -->     -->     <--     <--     -->  
!    !                 | left  |outflow| right | right |inflow | left  |outflow|
!    !              ---+-------+-------+-------+-------+-------+-------+-------+----
!    !
!    !  up loop:                   .   F
!    !                                     Q > F
!    !                                             Q > F
!    !                                                                     .   F
!    !
!    !  down loop: 
!    !                                                                 F < Q
!    !                                                         F < Q
!    !                                                     Q                      
!    !                         F < Q
!    !                 F < Q
!    !
!    
!    ! first going up for cells with outflow to right;
!    ! for each cell i, a right side tracer outflow Flux(i) becomes defined if present;
!    ! thus, left side tracer inflow Flux(i-1) is always previously defined for i>1,
!    ! but for first cell i=1 left side inflow Flux(0) should be defined first if necessary:
!    if ( phi(0) >= 0.0 ) Flux(0) = Q_t0(0) * phi(0) * dt  ! [(kg tracer)]
!    ! loop over cells from left to right:
!    do i = 1, n
!
!      ! inflow from right ? later ...
!      if ( phi(i) < 0.0 ) cycle
!
!      ! cell has outflow to right side;
!      ! also outflow at left side ?
!      if ( phi(i-1) < 0.0 ) then
!
!        ! outflow only cell, so mixing ratio probably remains the same,
!        ! but might be subject to limiting however;
!        ! new mixing ratio will be set in the downward loop,
!        ! compute right tracer outflow now:
!        Flux(i) = Q_t0(i) * phi(i) * dt  ! [(kg tracer)]
!
!      else
!
!        ! inflow from left and outflow to right;
!        ! tracer inflow flux Flux(i-1) is defined in this loop for the previous cell
!
!        ! 'upwind' Courant number; defines fraction of cell blown out:
!        cour = phi(i) * dt / m_t0(i)
!
!        ! extra scale factor for subgrid gradient to reduce numerical diffusion;
!        ! parameterisation is a result of experiments with advection of a
!        ! cone-shaped concentration:
!        if ( extreme(i-1) ) then
!          alfa = max( 1.5, 1.2+0.6*cour )   ! Eq. 10b
!        else if ( extreme(i+1) ) then
!          alfa = 1.75 - 0.45*cour           ! Eq. 10a
!        else
!          alfa = 1.0
!        end if
!        ! average mixing ratio in outflow flux:
!        Qf = Q_t0(i) + (Q_t0(i+1)-Q_t0(i-1)) * alfa * (1.0-cour)*0.25   ! Eq. 4a
!        ! limit to mixing ratios at each side of the face:
!        Qf = min( max( min(Q_t0(i),Q_t0(i+1)), Qf ), max(Q_t0(i),Q_t0(i+1)) )  ! Eq. 6a
!
!        ! first guess of new mixing ratio:
!        !          (kg tracer old)  from left     to right         (kg air new)
!        Qguess = ( Q_t0(i)*m_t0(i) + Flux(i-1) - cour*Qf*m_t0(i) ) / m_t1(i) ! [(kg tracer)/(kg air)]
!        ! limit to within original range of mixing ratios:
!        Q_t1(i) = min( max( Q_t1_min(i), Qguess ), Q_t1_max(i) )
!        
!        ! adjust right tracer flux since limitation changed total tracer mass:
!        !   right outflow = left inflow + mass decrease
!        Flux(i) = Flux(i-1) + Q_t0(i)*m_t0(i) -  Q_t1(i)*m_t1(i)  ! [(kg tracer)]
!        
!      end if
!    
!    end do
!
!    ! next going down for cells with outflow to left;
!    ! for each cell i, a left side tracer outflow Flux(i-1) becomes defined if present;
!    ! thus, right side tracer inflow Flux(i) is always previously defined for i<n ;
!    ! for 'first' cell i=n, define right side inflow Flux(n) if present:
!    if ( phi(n) < 0.0 ) Flux(n) = Q_t0(n+1) * phi(n) * dt  ! [(kg tracer)]
!    ! loop over cells from right to left:
!    do i = n, 1, -1
!
!      ! inflow from left ?
!      if ( u(i-1) >= 0.0 ) then
!
!        ! also inflow from right ?
!        if ( u(i) < 0.0 ) then
!          ! inflow only cell; compute first guess using mass balance:
!          Qguess = ( Q_t0(i)*m_t0(i) + Flux(i-1) - Flux(i) ) / m_t1(i)
!          ! left inflow tracer flux was calculated in upwards loop,
!          ! right inflow tracer flux was calculated in previous step of this loop;
!          ! since these can not be changed anymore, the new mixing ratio should
!          ! not be constrained since this might violate the mass balance ...
!          !Q_t1(i) = min( max( Q_t1_min(i), Qguess ), Q_t1_max(i) )
!          Q_t1(i) = Qguess
!          ! test ...
!          if ( (Q_t1(i) < Q_t1_min(i)) .or. (Q_t1(i) > Q_t1_max(i)) ) then
!            print *, 'CHECK - New mass-mixing-ratio in inflow-only cell violates limits:'
!            print *, 'CHECK -   cells           : ', i-1, i, i+1
!            print *, 'CHECK -   old air mass    : ', m_t0(i-1:i+1)
!            print *, 'CHECK -   new air mass    : ', m_t1(i-1:i+1)
!            print *, 'CHECK -   mixing ratios   : ', Q_t0(i-1:i+1)
!            print *, 'CHECK -   fluxes          : ', Flux(i-1:i)
!            print *, 'CHECK -   new mixing rat  : ', Q_t1(i)
!            print *, 'CHECK -   limits          : ', Q_t1_min(i), Q_t1_max(i)
!            stop 'CHECK - break in Walcek_1D'
!          end if
!        end if
!
!      else
!        
!        ! outflow to left;
!        ! right tracer outflow Flux(i) defined in upwards loop,
!        ! right tracer inflow Flux(i) is defined in this loop for previous cell
!
!        ! 'upwind' Courant number; defines fraction of cell blown out:
!        cour = abs(phi(i-1)) * dt / m_t0(i)
!
!        ! average mixing ratio in outflow flux:
!        if ( phi(i) >= 0.0 ) then
!          ! outflow only does not change mixing ratio;
!          ! therefore, no need to assume subgrid distribution
!          ! thus mixing ratio in outflow flux is same as cell average:
!          Qf = Q_t0(i)
!        else
!          ! outflow to left, inflow from right
!          ! assume linear subgrid mixing ratios;
!          ! extra scale factor for subgrid gradient to reduce numerical diffusion;
!          ! parameterisation is a result of experiments with advection of a
!          ! cone-shaped concentration:
!          if ( extreme(i-1) ) then
!            alfa = 1.75 - 0.45*cour           ! Eq. 10a
!          else if ( extreme(i+1) ) then
!            alfa = max( 1.5, 1.2+0.6*cour )   ! Eq. 10b
!          else
!            alfa = 1.0
!          end if
!          ! average mixing ratio in outflow flux:
!          Qf = Q_t0(i) - (Q_t0(i+1)-Q_t0(i-1)) * alfa * (1.0-cour)*0.25   ! Eq. 4a
!          ! limit to mixing ratios at each side of the face:
!          Qf = min( max( min(Q_t0(i-1),Q_t0(i)), Qf ), max(Q_t0(i-1),Q_t0(i)) )  ! Eq. 6a
!        end if
!
!        ! first guess of new mixing ratio:
!        !          (kg tracer old)       to left         to right      (kg air new)
!        Qguess = ( Q_t0(i)*D_t0(i) - cour*Qf*mf_t1(i-1) - Flux(i) ) / m_t1(i) ! [(kg tracer)/(kg air)]
!        ! limit to within original range of mixing ratios:
!        Q_t1(i) = min( max( Q_t1_min(i), Qguess ), Q_t1_max(i) )
!        
!        ! adjust left tracer flux since limitation changed total tracer mass:
!        !   left inflow = mass increase + right outflow
!        Flux(i-1) = Q_t1(i)*m_t1(i) - Q_t0(i)*m_t0(i) + Flux(i)  ! [(kg tracer)]
!        
!      end if
!    
!    end do
!  
!  end subroutine Walcek_1D_paper
!
!
!  ! ***
!
!
!  subroutine Walcek_1D( n, phi, m_t0, Q_t0, dt, m_t1, Q_t1 )
!  
!    ! --- in/out ------------------------------------------
!    
!    integer, parameter    ::  nh = 2  ! number of halo cells
!    
!    integer, intent(in)   ::  n                ! number of grid cells
!    real, intent(in)      ::  phi(0:n)         ! air mass flux through cell faces  [(kg air)/s]
!    real, intent(in)      ::  m_t0(1-nh:n+nh)  ! initial air mass in (halo) cells  [(kg air)]
!    real, intent(in)      ::  Q_t0(1-nh:n+nh)  ! initial mass-mixing-ratios in (halo) cells [(kg tracer)/(kg air)]
!    real, intent(in)      ::  dt               ! time step  [s]
!    real, intent(out)     ::  m_t1(1-nh:n+nh)  ! final air mass in (halo) cells  [(kg air)]
!    real, intent(out)     ::  Q_t1(1-nh:n+nh)  ! final mass-mixing-ratios in (halo) cells [(kg tracer)/(kg air)]
!    
!    ! --- local -------------------------------------------
!    
!    integer          ::  i
!    logical          ::  extreme(1-nh:n+nh)
!!    real             ::  Q_t1_min(1:n)
!!    real             ::  Q_t1_max(1:n)
!    real             ::  Flux(0:n)    ! tracer flux [(kg tracer)/m2]
!    real             ::  Qf
!    real             ::  cour
!    real             ::  alfa
!    
!    ! --- begin -------------------------------------------
!    
!    ! local extreme ?
!    extreme(1-nh:1-1) = .false.  ! no extemes in halo cells
!    extreme(n+1:n+nh) = .false.  ! no extemes in halo cells
!    do i = 1, n
!      extreme(i) = ( Q_t0(i) >= max(Q_t0(i-1),Q_t0(i+1)) ) .or. &
!                   ( Q_t0(i) <= min(Q_t0(i-1),Q_t0(i+1)) )
!    end do
!    
!!    ! tracer mixing ratio limits at new time:
!!    do i = 1, n
!!      if ( (u(i-1) >= 0.0) .and. (u(i) < 0.0) ) then
!!        ! inflow only:
!!        Q_t1_min(i) = min( min( Q_t0(i-1), Q_t0(i) ), Q_t0(i+1) )  ! line after Eq. 7
!!        Q_t1_max(i) = max( max( Q_t0(i-1), Q_t0(i) ), Q_t0(i+1) )  ! line after Eq. 7
!!      else if ( u(i-1) >= 0.0 ) then
!!        ! inflow from left:
!!        Q_t1_min(i) = min( Q_t0(i-1), Q_t0(i) )  ! Eq. 7a
!!        Q_t1_max(i) = max( Q_t0(i-1), Q_t0(i) )  ! Eq. 7a
!!      else if ( u(i) < 0.0 )
!!        ! inflow from right:
!!        Q_t1_min(i) = min( Q_t0(i), Q_t0(i+1) )  ! Eq. 7b
!!        Q_t1_max(i) = max( Q_t0(i), Q_t0(i+1) )  ! Eq. 7b
!!      else
!!        ! outflow only does not change tracer mixing ratio ...
!!        Q_t1_min(i) = Q_t0(i)
!!        Q_t1_max(i) = Q_t0(i)
!!      end if
!!    end do
!    
!    ! air mass at end of time step:
!    do i = 1, n
!      m_t1(i) = m_t0(i) + ( phi(i-1) - phi(i) )*dt  ! (kg air)
!    end do
!
!    !
!    ! Loop over cell interfaces from left to right:
!    ! fill right tracer flux (tracer vmr in flux through interface)
!    ! assuming linear distribution  (here: fluxx)
!    !
!    ! Loop over cells from left to right:
!    ! for all outflow-only cells,
!    ! replace Qf at both sides given vmr and air-mass-flux, new vmr.  
!    !
!    ! Loop over cells from left to right:
!    ! for rightflow-only cells, compute new vmr given previously computed fluxes;
!    ! limit mixing ratio, if necessary re-compute right flux given left flux
!    ! and mass change.
!    !
!    ! Loop over cells from right to left:
!    ! for left-only cells (bug: test for left outflow) compute new vmr;
!    ! limit mixing ratio, if necessary re-compute left flux given right flux
!    ! and mass change.
!    !
!    !                <--     <--     -->     -->     -->     <--     <--     -->  
!    !                 | left  |outflow| right | right |inflow | left  |outflow|
!    !              ---+-------+-------+-------+-------+-------+-------+-------+----
!    !
!    ! 0,nx           F,a     F,a     a,F     a,F     a,F     F,a     F,a     a,F
!    !
!    ! 1,nx                    F       F                               F       F
!    !                             Q                                       Q
!    !
!    ! 1,nx                                Q > F
!    !                                             Q > F
!    !
!    ! nx,1,-1                                                                    
!    !                                                         F < Q
!    !                                                     Q
!    !                 F < Q
!    !
!    
!    ! Loop over cell interfaces from left to right;
!    ! fill all Qf (tracer vmr in flux through interface)
!    ! assuming linear distribution 
!    do i = 0, n
!    
!      ! left or right flow through this interface ?
!      if ( phi(i) >= 0.0 ) then
!      
!        ! flow from left to right ...
!        
!        ! 'upwind' Courant number; defines fraction of air leaving left cell:
!        cour = phi(i) * dt / m_t0(i-1)  ! [0-1]
!
!        ! extra scale factor for subgrid gradient to reduce numerical diffusion;
!        ! parameterisation is a result of experiments with advection of a
!        ! cone-shaped concentration:
!        if ( extreme(i-1) ) then
!          alfa = max( 1.5, 1.2+0.6*cour )   ! Eq. 10b
!        else if ( extreme(i+1) ) then
!          alfa = 1.75 - 0.45*cour           ! Eq. 10a
!        else
!          alfa = 1.0
!        end if
!        
!        ! average mixing ratio in outflow flux:
!        Qf = Q_t0(i) + (Q_t0(i+1)-Q_t0(i-1)) * alfa * (1.0-cour)*0.25   ! Eq. 4a
!        
!      else
!      
!        ! flow from right to left ...
!        
!        ! 'upwind' Courant number; defines fraction of air leaving right cell:
!        cour = phi(i) * dt / m_t0(i)
!
!        ! extra scale factor for subgrid gradient to reduce numerical diffusion;
!        ! parameterisation is a result of experiments with advection of a
!        ! cone-shaped concentration:
!        if ( extreme(i-1) ) then
!          alfa = max( 1.5, 1.2+0.6*cour )   ! Eq. 10b
!        else if ( extreme(i+1) ) then
!          alfa = 1.75 - 0.45*cour           ! Eq. 10a
!        else
!          alfa = 1.0
!        end if
!        
!        ! average mixing ratio in outflow flux:
!        Qf = Q_t0(i+1) - (Q_t0(i+2)-Q_t0(i)) * alfa * (1.0-cour)*0.25   ! Eq. 4a
!        
!      end if
!    
!      ! limit to mixing ratios at each side of the face:
!      Qf = min( max( min(Q_t0(i),Q_t0(i+1)), Qf ), max(Q_t0(i),Q_t0(i+1)) )  ! Eq. 6a
!
!      ! tracer flux through surface:
!      !    (kg tracer)/(kg air) * (kg air)/s * s
!      Flux(i) = Qf * phi(i) * dt  ! [(kg tracer)]
!      
!    end do
!
!    ! Loop over cells from left to right:
!    ! for all inflow-only or outflow-only cells,
!    ! replace Qf at both sides given vmr and air-mass-flux, new vmr.  
!    do i = 1, n
!    
!      ! inflow- or outflow only ?
!      if ( (phi(i-1) >= 0.0) .and. (phi(i) <= 0.0) ) then
!        
!        ! inflow only ...
!        
!        ! tracer fluxes through surfaces:
!        Flux(i-1) = Q_t0(i-1) * phi(i-1) * dt  ! [(kg tracer)]
!        Flux(i  ) = Q_t0(i+1) * phi(i  ) * dt  ! [(kg tracer)]
!        
!        ! first guess of new tracer mixing ratio:
!        Qguess(i) = ( Q_t0(i)*m_t0(i) + Flux(i-1) - Flux(i) ) / m_t1(i)   ! [(kg tracer)/(kg air)]
!        
!      else
!      
!        ! outflow only ...
!        
!        ! tracer fluxes through surfaces:
!        Flux(i-1) = Q_t0(i) * phi(i-1) * dt  ! [(kg tracer)]
!        Flux(i  ) = Q_t0(i) * phi(i  ) * dt  ! [(kg tracer)]
!    
!        ! new tracer mixing ratio:
!        Q_t1(i) = ( Q_t0(i)*m_t0(i) + Flux(i-1) - Flux(i) ) / m_t1(i)   ! [(kg tracer)/(kg air)]
!      
!      end if
!      
!    end do
!
!    ! Loop over cells from left to right:
!    ! for rightflow-only cells, compute new vmr given previously computed fluxes;
!    ! limit mixing ratio, if necessary re-compute right flux given left flux and mass change.
!    do i = 1, n
!    
!      ! right-only flow:
!      if ( (phi(i) >= 0.0) .and. (phi(i) > 0.0) ) then
!      
!        ! first guess of new tracer mixing ratio:
!        Qgues = ( Q_t0(i)*m_t0(i) + Flux(i-1) - Flux(i) ) / m_t1(i)
!        
!        ! limit to range of original mixing ratios in left and this cell:
!        Qf = min( max( min(Q_t0(i-1),Q_t0(i)), Qguess ), max(Q_t0(i-1),Q_t0(i)) )  ! Eq. 6a
!        
!    
!    end do
!
!
!! --------------
!
!
!    ! next going down for cells with outflow to left;
!    ! for each cell i, a left side tracer outflow Flux(i-1) becomes defined if present;
!    ! thus, right side tracer inflow Flux(i) is always previously defined for i<n ;
!    ! for 'first' cell i=n, define right side inflow Flux(n) if present:
!    if ( phi(n) < 0.0 ) Flux(n) = Q_t0(n+1) * phi(n) * dt  ! [(kg tracer)]
!    ! loop over cells from right to left:
!    do i = n, 1, -1
!
!      ! inflow from left ?
!      if ( u(i-1) >= 0.0 ) then
!
!        ! also inflow from right ?
!        if ( u(i) < 0.0 ) then
!          ! inflow only cell; compute first guess using mass balance:
!          Qguess = ( Q_t0(i)*m_t0(i) + Flux(i-1) - Flux(i) ) / m_t1(i)
!          ! left inflow tracer flux was calculated in upwards loop,
!          ! right inflow tracer flux was calculated in previous step of this loop;
!          ! since these can not be changed anymore, the new mixing ratio should
!          ! not be constrained since this might violate the mass balance ...
!          !Q_t1(i) = min( max( Q_t1_min(i), Qguess ), Q_t1_max(i) )
!          Q_t1(i) = Qguess
!          ! test ...
!          if ( (Q_t1(i) < Q_t1_min(i)) .or. (Q_t1(i) > Q_t1_max(i)) ) then
!            print *, 'CHECK - New mass-mixing-ratio in inflow-only cell violates limits:'
!            print *, 'CHECK -   cells           : ', i-1, i, i+1
!            print *, 'CHECK -   old air mass    : ', m_t0(i-1:i+1)
!            print *, 'CHECK -   new air mass    : ', m_t1(i-1:i+1)
!            print *, 'CHECK -   mixing ratios   : ', Q_t0(i-1:i+1)
!            print *, 'CHECK -   fluxes          : ', Flux(i-1:i)
!            print *, 'CHECK -   new mixing rat  : ', Q_t1(i)
!            print *, 'CHECK -   limits          : ', Q_t1_min(i), Q_t1_max(i)
!            stop 'CHECK - break in Walcek_1D'
!          end if
!        end if
!
!      else
!        
!        ! outflow to left;
!        ! right tracer outflow Flux(i) defined in upwards loop,
!        ! right tracer inflow Flux(i) is defined in this loop for previous cell
!
!        ! 'upwind' Courant number; defines fraction of cell blown out:
!        cour = abs(phi(i-1)) * dt / m_t0(i)
!
!        ! average mixing ratio in outflow flux:
!        if ( phi(i) >= 0.0 ) then
!          ! outflow only does not change mixing ratio;
!          ! therefore, no need to assume subgrid distribution
!          ! thus mixing ratio in outflow flux is same as cell average:
!          Qf = Q_t0(i)
!        else
!          ! outflow to left, inflow from right
!          ! assume linear subgrid mixing ratios;
!          ! extra scale factor for subgrid gradient to reduce numerical diffusion;
!          ! parameterisation is a result of experiments with advection of a
!          ! cone-shaped concentration:
!          if ( extreme(i-1) ) then
!            alfa = max( 1.5, 1.2+0.6*cour )   ! Eq. 10b
!          else if ( extreme(i+1) ) then
!            alfa = 1.75 - 0.45*cour           ! Eq. 10a
!          else
!            alfa = 1.0
!          end if
!          ! average mixing ratio in outflow flux:
!          Qf = Q_t0(i) - (Q_t0(i+1)-Q_t0(i-1)) * alfa * (1.0-cour)/4.0   ! Eq. 4a
!          ! limit to mixing ratios at each side of the face:
!          Qf = min( max( min(Q_t0(i-1),Q_t0(i)), Qf ), max(Q_t0(i-1),Q_t0(i)) )  ! Eq. 6a
!        end if
!
!        ! first guess of new mixing ratio:
!        !          (kg tracer old)       to left         to right      (kg air new)
!        Qguess = ( Q_t0(i)*D_t0(i) - cour*Qf*mf_t1(i-1) - Flux(i) ) / m_t1(i) ! [(kg tracer)/(kg air)]
!        ! limit to within original range of mixing ratios:
!        Q_t1(i) = min( max( Q_t1_min(i), Qguess ), Q_t1_max(i) )
!        
!        ! adjust left tracer flux since limitation changed total tracer mass:
!        !   left inflow = mass increase + right outflow
!        Flux(i-1) = Q_t1(i)*m_t1(i) - Q_t0(i)*m_t0(i) + Flux(i)  ! [(kg tracer)]
!        
!      end if
!    
!    end do
!  
!  end subroutine Walcek_1D


end module LE_Advec

