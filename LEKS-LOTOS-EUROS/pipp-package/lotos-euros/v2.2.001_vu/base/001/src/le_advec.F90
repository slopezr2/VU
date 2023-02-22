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

  !! timers:
  !integer   ::  itim_adv_x
  !integer   ::  itim_adv_y
  !integer   ::  itim_adv_z
  !integer   ::  itim_adv_dif


contains


  ! ====================================================================


  subroutine LE_Advec_Init( status )

    use GO, only : GO_Timer_Def

    ! --- in/out ---------------------------------

    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ LE_Advec_Init'

    ! number of halo cells:
    integer, parameter  ::  nh = 2

    ! --- begin ----------------------------------

    !! define timers:
    !call GO_Timer_Def( itim_adv_x  , 'advection x', status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_adv_y  , 'advection y', status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_adv_z  , 'advection z', status )
    !IF_NOTOK_RETURN(status=1)
    !call GO_Timer_Def( itim_adv_dif, 'advection dif', status )
    !IF_NOTOK_RETURN(status=1)

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

    use GO           , only : goc
    use Dims, only : nx, ny, nz
    use LE_Meteo_Data, only : volume

    ! --- in/out ---------------------------------

    real, intent(in)                ::  dt_min  ! minutes
    integer, intent(out)            ::  nstep
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Advec_Get_NStep'

    ! --- local ----------------------------------

    real, allocatable  ::  inv_volume(:,:,:)
    real               ::  dta

    ! --- begin ----------------------------------

    ! storage for inverted volume, incl halo:
    allocate( inv_volume(0:nx+1,0:ny+1,1:nz+1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! pre-compute 1.0/volume since devision is much more expensive than multiplication ...
    inv_volume = 1.0 / volume

    ! check courant condition, store cfl numbers:
    call courant_a( nx,ny,nz, inv_volume, dt_min, dta, nstep, store=.true. )

    ! round to even number:
    nstep = int( ceiling(0.5*nstep) * 2 )

    ! maximum over all domains:
    call goc%AllReduce( 'max', nstep, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( inv_volume, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Advec_Get_NStep


  ! ***


  subroutine LE_Advec_Apply( c, dt, status )

    use GO           , only : goc
    use GO           , only : GO_Timer_Start, GO_Timer_End
    use Indices      , only : n_advected, ispecs_advected
    use Dims         , only : nx, ny, nz, nspec
    use Dims         , only : runF, outF
    use LE_Logging   , only : ident2
    use LE_Grid      , only : dom
    use LE_Bound     , only : bc_west, bc_east, bc_north, bc_south
    use LE_Bound     , only : caloft  ! (nx,ny,nz+1:nz_top,nspec)
    use LE_Meteo_Data, only : volume
    use LE_Meteo_Data, only : afluxx, afluxy, afluxz
#ifdef with_hdiff
    use Dims         , only : khx, khy
    use hdiff        , only : hordif
#endif

!          use LE_Grid, only : debug_i, debug_j

    ! --- begin ------------------------------

    real, intent(in)       ::  dt
    real, intent(inout)    ::  c(nx,ny,nz,nspec)
    integer, intent(out)   ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Advec_Apply'

    ! number of halo cells needed:
    integer, parameter  ::  nh = 2

    ! --- local ------------------------------

    integer   ::  istep, ix, iy, iz
    integer   ::  itr, k

    ! actual time step and number of steps for advection:
    real      ::  dta
    integer   ::  nstep_a

#ifdef with_hdiff
    ! actual time step and number of steps for diffusion:
    real      ::  dtd
    integer   ::  nstep_d
#endif

    ! pointer array since ExchangeHalo needs this argument
    ! in order to keep index space:
    real, pointer      ::  ch(:,:,:)
    real, pointer      ::  volume_loc(:,:,:)

    ! local arrays:
    real, allocatable  ::  dens0(:,:,:)
    real, allocatable  ::  dens1(:,:,:)
    real, allocatable  ::  inv_d2(:,:,:)
    real, allocatable  ::  inv_volume(:,:,:)
    real, allocatable  ::  volume0(:,:,:)

    logical            ::  has_west, has_east, has_south, has_north
    integer            ::  ih
    integer            ::  i, j
    
!        integer       ::  qi, qj
!        integer       ::  off(2)

    ! status argument within OpenMP loop:
    integer            ::  status_par
    

    ! --- begin ----------------------------------

    !! info:
    !if ( .not. outF%suppress ) then
    !  write (gol,'(a,"<advection>")') ident2; call goPr
    !end if

    !AJS: routine is called for each mode with first flag, thus moved to kf_lotdriver
    !if (runF%first) call init_advec

    ! 
    ! * cfl check
    !
    
    ! storage:
    allocate( inv_volume(0:nx+1,0:ny+1,1:nz+1), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! pre-compute 1.0/volume since devision is much more expensive than multiplication ...
    inv_volume = 1.0 / volume

    ! check courant condition
    call courant_a( nx,ny,nz, inv_volume, dt, dta, nstep_a )
#ifdef with_hdiff
    call courant_d( nx,ny,nz, inv_volume, dt, dtd, nstep_d )
    if (nstep_d > nstep_a .AND. (.NOT.outF%suppress) ) then
      write (*,*) ident2, '  >>> diffusion dominates advection'
      write (*,*) ident2, '  >>> number of steps required for hor diffusion: ', nstep_d
    endif
#endif
    ! clear:
    deallocate( inv_volume, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! minimum/maximum over all domains:
    call goc%AllReduce( 'min', dta, status )
    IF_NOTOK_RETURN(status=1)
    call goc%AllReduce( 'max', nstep_a, status )
    IF_NOTOK_RETURN(status=1)

    ! safety check:
    if ( nstep_a /= 1 ) then
      !write (gol,'("something wrong with operator splitting step,")'); call goErr
      !write (gol,'("more than 1 sub-step required for advection : ",i6)') nstep_a; call goErr
      !TRACEBACK; status=1; return
      write (gol,'("WARNING - operator splitting step not ok yet (temporal interpol?), ",i0," adv. steps needed")') nstep_a; call goPr
    end if

    !
    ! * advection steps
    !

    ! volume copy:
    allocate( volume0(0:nx+1,0:ny+1,1:nz+1), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! save current volume, later will be updated
    ! after advection of first tracer:
    volume0 = volume

    ! local domain at border?
    call dom%Get( status, has_west=has_west, has_east=has_east, has_south=has_south, has_north=has_north )
    IF_NOTOK_RETURN(status=1)

    ! init loop status:
    status_par = 0

    !xOMP parallel &
#ifndef __GFORTRAN__
    !xOMP   default ( none ) &
    !xOMP   shared ( ispecs_advected ) &
#endif
    !xOMP   shared ( dom ) &
    !xOMP   shared ( nx, ny, nz ) &
    !xOMP   shared ( nstep_a, dta ) &
    !xOMP   shared ( volume ) &
    !xOMP   shared ( volume0 ) &
    !xOMP   shared ( bc_west, bc_east, bc_south, bc_north, caloft ) &
    !xOMP   shared ( has_west, has_east, has_south, has_north ) & 
    !xOMP   shared ( afluxx, afluxy, afluxz ) &
    !xOMP   shared ( c ) &
    !xOMP   shared ( status_par ) &
    !xOMP   shared ( gol ) &
    !xOMP   private( k, istep ) &
    !xOMP   private( itr ) &
    !xOMP   private( ch ) &
    !xOMP   private( dens0, dens1 ) &
    !xOMP   private( inv_d2 ) &
    !xOMP   private( inv_volume ) &
    !xOMP   private( volume_loc ) &
    !xOMP   private( ix,iy,iz ) &
    !xOMP   private( status )

    ! private arrays:
    allocate( ch(1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh), stat=status )
    if (status/=0) status_par = status_par + 1
    allocate( dens0(nx, ny, nz), stat=status )
    if (status/=0) status_par = status_par + 1
    allocate( dens1(nx, ny, nz), stat=status )
    if (status/=0) status_par = status_par + 1
    allocate( inv_d2(nx, ny, nz), stat=status )
    if (status/=0) status_par = status_par + 1
    ! volumes incl 1 halo:
    allocate( volume_loc(0:nx+1,0:ny+1,1:nz+1), stat=status )
    if (status/=0) status_par = status_par + 1
    allocate( inv_volume(0:nx+1,0:ny+1,1:nz+1), stat=status )
    if (status/=0) status_par = status_par + 1

    ! loop over advected tracers:
    !xOMP do
    do itr = 1, n_advected
      ! current index:
      k = ispecs_advected(itr)


      ! copy initial volume:
      volume_loc = volume0

      ! copy the concentration vector into the help array ch
      !ch = 0.0
      ch(1:nx,1:ny,1:nz) = c(1:nx,1:ny,1:nz,k)

      ! put bc's in dummy cells of stage vector and help conc. array
      do ih = 1, nh
        ! horizontal:
        if ( has_west  ) ch( 1-ih, 1:ny,1:nz) = bc_west (1:ny,1:nz,k)
        if ( has_east  ) ch(nx+ih, 1:ny,1:nz) = bc_east (1:ny,1:nz,k)
        if ( has_south ) ch( 1:nx, 1-ih,1:nz) = bc_south(1:nx,1:nz,k)
        if ( has_north ) ch( 1:nx,ny+ih,1:nz) = bc_north(1:nx,1:nz,k)
        ! upper:
        ch(1:nx,1:ny,nz+ih) = caloft(:,:,nz+1,k)
      end do ! halo cells

      ! do the necessary number of time steps
      do istep = 1, nstep_a

        ! 
        ! ~ x-direction
        !

        ! exchange concentrations in halo cells in x direction:
        call dom%ExchangeHaloX( nh, ch, status )
        if (status /= 0 ) then
          write (gol, '(" Error from exchange of concentration halo cells")' ); call goErr
          status_par = status_par + 1
        end if
        ! advection steps have changed the volumes;
        ! exchange single cell of halo cells in x direction:
        call dom%ExchangeHaloX( 1, volume_loc, status )
        if (status /= 0 ) then
          write (gol, '(" Error from exchange of volume halo cells")' ); call goErr
          status_par = status_par + 1
        end if
        ! re-compute inverse:
        inv_volume = 1.0 / volume_loc

        ! fill densities:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx
              ! set initial 'density' (fraction of air density):
              dens0(ix,iy,iz) = d0

              ! update 'density' ; aflux* in (m3 air)/s through interface:
              !                        1            s                   m3/s                           1/m3
              dens1(ix,iy,iz) = dens0(ix,iy,iz) + dta*( afluxx(ix-1,iy,iz) - afluxx(ix,iy,iz) ) * inv_volume(ix,iy,iz)
              ! pre-compute inverse:
              inv_d2(ix,iy,iz) = 1.0 / dens1(ix,iy,iz)
            end do  ! ix
          end do  ! iy
        end do  ! iz
        
        !! start timing:
        !call GO_Timer_Start( itim_adv_x, status )
        !if (status/=0) status_par = status_par + 1

        ! perform advection:
        call advecX( nx,ny,nz, volume_loc, inv_volume, &
                     dens0, dens1, k, ch, inv_d2, dta, status )
        if ( status /= 0 ) then
          write (gol, '(" Error in Advection in x-direction")' ); call goErr
          status_par = status_par + 1
        end if

        ! update volume:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx
              ! add volume change due to x-flux (might be zero):
              !          km3                    km3         min                    km3/s
              volume_loc(ix,iy,iz) = volume_loc(ix,iy,iz) + dta*( afluxx(ix-1,iy,iz) - afluxx(ix,iy,iz) )
            end do  ! ix
          end do  ! iy
        end do  ! iz
        ! re-compute inverse:
        inv_volume = 1.0 / volume_loc

        !! end timing:
        !call GO_Timer_End( itim_adv_x, status )
        !if (status/=0) status_par = status_par + 1

        !
        ! ~ y-direction
        !

        ! exchange concentrations in halo cells in y direction:
        call dom%ExchangeHaloY( nh, ch, status )
        if (status /= 0 ) then
          write (gol, '(" Error from exchange of concentration halo cells")' ); call goErr
          status_par = status_par + 1
        end if
        

        ! advection steps have changed the volumes;
        ! exchange single cell of halo cells in y direction:
        call dom%ExchangeHaloY( 1, volume_loc, status )
        if (status /= 0 ) then
          write (gol, '(" Error from exchange of volume halo cells")' ); call goErr
          status_par = status_par + 1
        end if
        ! re-compute inverse:
        inv_volume = 1.0 / volume_loc

        ! fill densities:
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
        
        !! start timing:
        !call GO_Timer_Start( itim_adv_y, status )
        !if (status/=0) status_par = status_par + 1

        ! perform advection:
        call advecY( nx,ny,nz, volume_loc, inv_volume, &
                       dens1, dens0, k, ch, inv_d2, dta, status )
        if (status /= 0 ) then
          write (gol, '(" Error in Advection in y-direction")' ); call goErr
          status_par = status_par +1
        end if

        ! update volume:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx
              ! add volume change due to y-flux (might be zero):
              !          km3                    km3         min                    km3/s
              volume_loc(ix,iy,iz) = volume_loc(ix,iy,iz) + dta*( afluxy(ix,iy-1,iz) - afluxy(ix,iy,iz) )
            end do  ! ix
          end do  ! iy
        end do  ! iz
        ! re-compute inverse:
        inv_volume = 1.0 / volume_loc

        !! end timing:
        !call GO_Timer_End( itim_adv_y, status )
        !if (status/=0) status_par = status_par + 1

        !
        ! ~ z-direction
        !

        ! set the lower boundary conditions here. They are kind of unrealistic,
        ! but the scheme uses boundary cells to determine minimum and maximum
        ! concentrations at the cell faces.
        ! Therefore we give the actual values of the lowest layers
        do ih = 1, nh
          ch(:,:,1-ih) = ch(:,:,1)
        end do

        ! fill densities:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx
              ! update density
              dens1(ix,iy,iz) = dens0(ix,iy,iz) + dta*(afluxz(ix,iy,iz-1) &
                                   - afluxz(ix,iy,iz)) * inv_volume(ix,iy,iz)
              ! pre-compute inverse:
              inv_d2(ix,iy,iz) = 1.0 / dens1(ix,iy,iz)
            end do ! ix
          end do ! iy
        end do ! iz

        !! start timing:
        !call GO_Timer_Start( itim_adv_z, status )
        !if (status/=0) status_par = status_par + 1

        ! perform advection
        call advecZ( nx,ny,nz, volume_loc, inv_volume, &
                      dens0, dens1, k, ch, inv_d2, dta, status )
        if (status /= 0 ) then
          write (gol, '(" Error in Advection in z-direction")' ); call goErr
          status_par = status_par + 1
        end if

        ! update volume:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx
              ! add volume change due to z-flux (might be zero):
              volume_loc(ix,iy,iz) = volume_loc(ix,iy,iz) + &
                        dta*( afluxz(ix,iy,iz-1) - afluxz(ix,iy,iz) )
            end do  ! ix
          end do  ! iy
        end do  ! iz

        !! end timing:
        !call GO_Timer_End( itim_adv_z, status )
        !if (status/=0) status_par = status_par + 1

      end do  ! end loop over #advection steps

      ! put ch into c, the global concentration array
      c(1:nx,1:ny,1:nz,k) = ch(1:nx,1:ny,1:nz)

      ! reset global volume, changed in the same way for every tracer,
      ! copy if first is done to avoid two threads changing the same array:
      if ( itr == 1 ) then
        volume = volume_loc
      end if

    end do    ! end loop over species
    !xOMP end do

    ! clear:
    deallocate( ch, stat=status )
    if (status/=0) status_par = status_par + 1
    deallocate( dens0, stat=status )
    if (status/=0) status_par = status_par + 1
    deallocate( dens1, stat=status )
    if (status/=0) status_par = status_par + 1
    deallocate( inv_d2, stat=status )
    if (status/=0) status_par = status_par + 1
    deallocate( inv_volume, stat=status )
    if (status/=0) status_par = status_par + 1
    deallocate( volume_loc, stat=status )
    if (status/=0) status_par = status_par + 1
    !xOMP end parallel
    if ( status_par /= 0 ) then
      write(gol, '("non-zero return status from OpenMP loop: ",i0)') status_par; call goErr
      TRACEBACK; status=1;return
    endif

#ifdef with_hdiff
    ! now do the horizontal diffusion
    !call GO_Timer_Start( itim_adv_dif, status )
    !IF_NOTOK_RETURN(status=0)
    !xOMP parallel &
    !xOMP   default( none ) &
    !xOMP   private( ispec, istep )
    !xOMP   do
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
    !xOMP   end do
    !xOMP end parallel
    !call GO_Timer_End( itim_adv_dif, status )
    !IF_NOTOK_RETURN(status=0)
#endif

    ! clear:
    deallocate( volume0, stat=status )
    IF_NOTOK_RETURN(status=0)

    ! ok
    status = 0

  end subroutine LE_Advec_Apply


  ! ***


  subroutine advecX( nx,ny,nz, volume, inv_volume, d1, d2, k, &
                          ch, inv_d2, dta, status )

    use LE_Grid      , only : dom
    use LE_Meteo_Data, only : afluxx

    ! --- const ----------------------------------

    character(len=*), parameter :: rname = mname//'/advecX'

    ! number halo cells in concentration array:
    integer, parameter  ::  nh = 2

    ! --- in/out ---------------------------------

    integer, intent(in)           ::  nx, ny, nz
    real, intent(in), pointer     ::  volume(:,:,:)    ! (0:nx+1,0:ny+1,1:nz+1)
    real, intent(in)              ::  inv_volume(0:nx+1,0:ny+1,1:nz+1)
    real, intent(in)              ::  d1(nx,ny,nz)
    real, intent(in)              ::  d2(nx,ny,nz)
    real, intent(in)              ::  dta
    integer, intent(in)           ::  k
    real, intent(inout), pointer  ::  ch(:,:,:)        ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in)              ::  inv_d2(nx,ny,nz)
    integer, intent(out)          ::  status

    ! --- local ----------------------------------

    integer               ::  cart_size, cart_coord
    integer               ::  coord

    logical, pointer      ::  is_extreme(:,:,:)
    real, allocatable     ::  cgues(:,:,:)
    real, pointer         ::  fluxx(:,:,:)

    integer               ::  ix, iy, iz
    real                  ::  c1, c2
    real                  ::  u
    real                  ::  cour, fcour
    real                  ::  alpha
    real                  ::  cg
    
      !integer   ::  qi, qj, qk

    ! --- begin ----------------------------------

    ! Info on domain layout in first dimension:
    call dom%GetDim( 1, status, cart_size=cart_size, cart_coord=cart_coord )
    IF_NOTOK_RETURN(status=1)

    ! storage for new concentrations:
    allocate ( cgues(nx,ny,nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for tracer flux:
    allocate ( fluxx(0:nx,ny,nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! initialize flux
    fluxx = 0.0

    ! storage for extremes including 2 halo cells in x-direction:
    allocate ( is_extreme(-1:nx+2,ny,nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! init, only outside the domain this remains:
    is_extreme = .false.
    ! loop over x cells incl x halo:
    do iz = 1, nz
      do iy = 1, ny
        do ix = 1, nx
          ! concentrations in neighbour cells:
          c1 = ch(ix-1,iy,iz)
          c2 = ch(ix+1,iy,iz)
          ! concentrations here above/below max/min of neighbours?
          is_extreme(ix,iy,iz) = (ch(ix,iy,iz) >= max(c1,c2)) .or. (ch(ix,iy,iz) <= min(c1,c2))
        end do
      end do
    end do
    ! exchange halo cells:
    call dom%ExchangeHaloX( 2, is_extreme, status )
    IF_NOTOK_RETURN(status=1)

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


    ! compute the fluxes;
    ! loop over cell interfaces in x direction:
    do iz = 1, nz
      do iy = 1, ny
        do ix = 0, nx
          ! use volume flux through this interface for 'wind':
          u = afluxx(ix,iy,iz)
          ! courant number using volume left from the interface:
          !AJS: bug? should be 'upwind' courrant number with different volume given sign:
          !AJS:   cour = dta *     u(ix)  / volume(ix  )   , u >= 0.0
          !AJS:   cour = dta * abs(u(ix)) / volume(ix+1)   , u <  0.0
          ! FIX for domain decomposition:
          ! only skip first interface if at boundary:
          if ( cart_coord == 0 ) then
            cour = dta * abs(u) * inv_volume(max(1,ix),iy,iz)
          else
            cour = dta * abs(u) * inv_volume(      ix ,iy,iz)
          end if

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
          !
          ! cell -1     0     1     2      3
          ! ix             0     1     2     ...
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
          ! NOTE: loop over edges 0,1,..,nx  ,  using extremes {-1,1},{0,2},...,{nx-1,nx+1} ;
          !   extremes are based on concentrations in surrouding cells,
          !   thus concentrations needed for -2,-1,0,..,nx+1,nx+2  ; strange, assymetric!
          fcour = (1.0-cour)/4.0
          alpha = 1.0
          if ( u >= 0.0 ) then
            !AJS: strange .. assymetric!
            if (is_extreme(ix-1,iy,iz)) alpha = max(1.5,1.2+0.6*cour)
            if (is_extreme(ix+1,iy,iz)) alpha = (1.75-0.45*cour)
            !!AJS: probably this was needed:
            !if (is_extreme(ix  ,iy,iz)) alpha = max(1.5,1.2+0.6*cour)
            !if (is_extreme(ix+1,iy,iz)) alpha = (1.75-0.45*cour)
            ! 'flux' through interface, but actually concentration at interface:
            fluxx(ix,iy,iz) = ch(ix,iy,iz) + (ch(ix+1,iy,iz) - ch(ix-1,iy,iz))*fcour*alpha
          else
            !AJS: strange .. assymetric!
            if (is_extreme(ix+1,iy,iz)) alpha = max(1.5,1.2+0.6*cour)
            if (is_extreme(ix-1,iy,iz)) alpha = (1.75-0.45*cour)
            !!AJS: probably this was needed:
            !if (is_extreme(ix+1,iy,iz)) alpha = max(1.5,1.2+0.6*cour)
            !if (is_extreme(ix  ,iy,iz)) alpha = (1.75-0.45*cour)
            ! 'flux' through interface, but actually concentration at interface:
            fluxx(ix,iy,iz) = ch(ix+1,iy,iz) + (ch(ix,iy,iz) - ch(ix+2,iy,iz))*fcour*alpha
          end if

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
          ! concentration range:
          c1 = min( ch(ix,iy,iz), ch(ix+1,iy,iz) )
          c2 = max( ch(ix,iy,iz), ch(ix+1,iy,iz) )
          ! bound cell face concentration:
          fluxx(ix,iy,iz) = max( c1, min( fluxx(ix,iy,iz), c2 ) )
          ! make a flux out of it ....
          fluxx(ix,iy,iz) = u * dta * fluxx(ix,iy,iz) * d0

        end do ! ix (interfaces)
      end do ! iy
    end do ! iz
    
    ! *

    ! First handle exceptions: inflow or outflow only cells.
    ! In decomposed domain, the new concentrations 'cguess' that are filled
    ! will be the same independend of the decomposition, since it only depends
    ! on the left/right tracer fluxes 'fluxx' values that are recomputed from the
    ! left/right airmass fluxes in 'afluxx', and the halo concentrations.
    ! However, the tracers fluxes 'fluxx' might have been updated for the last
    ! cell in the domain on the west, and not being update for the first cell
    ! in this domain, and a difference over the domains could arise.
    ! Therefore perform this step from left to right, and copy the 'fluxx' 
    ! from the domain on the west.
    ! Loop over domains from west to east:
    do coord = 0, cart_size-1
    
      ! receive values from the domain on the west:
      call dom%TransferSlab( coord, 'from-west', 1, fluxx, status )
      IF_NOTOK_RETURN(status=1)
      
      ! is local domain in this column?
      if ( coord == cart_coord ) then
      
        ! loop over cells:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx

              !~ inflow only?
              if ( (afluxx(ix,iy,iz) <= 0.0) .and. (afluxx(ix-1,iy,iz) >= 0.0) ) then
                ! (re)set tracer fluxes:
                fluxx(ix-1,iy,iz) = dta * afluxx(ix-1,iy,iz) * ch(ix-1,iy,iz) * d0
                fluxx(ix  ,iy,iz) = dta * afluxx(ix  ,iy,iz) * ch(ix+1,iy,iz) * d0
                ! update concentartion:
                cgues(ix,iy,iz) = ( ch(ix,iy,iz) * d1(ix,iy,iz) + &
                                    (fluxx(ix-1,iy,iz) - fluxx(ix,iy,iz)) * inv_volume(ix,iy,iz) ) &
                                  * inv_d2(ix,iy,iz)

              !~ outflow only?
              else if (afluxx(ix,iy,iz) >= 0.0 .AND. afluxx(ix-1,iy,iz) <= 0.0) then
                ! (re)set tracer fluxes:
                fluxx(ix-1,iy,iz) = dta*afluxx(ix-1,iy,iz)*ch(ix,iy,iz)*d0
                fluxx(ix  ,iy,iz) = dta*afluxx(ix  ,iy,iz)*ch(ix,iy,iz)*d0
                ! update concentartion:
                cgues(ix,iy,iz) = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
                       (fluxx(ix-1,iy,iz) - fluxx(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)

              end if ! inflow or outflow only

            end do  ! ix
          end do ! iy
        enddo ! iz
        
      end if  ! local domain in this cart dim
    
      ! fluxes on left side might have changed, send to previous domain:
      call dom%TransferSlab( coord, 'to-west', 1, fluxx, status )
      IF_NOTOK_RETURN(status=1)
      
    end do  ! domains
    
    ! *
    
    ! Limit fluxes, first going up.
    ! This might change the 'fluxx' on the right side of the cell,
    ! therefore obtain latest value from domain on the west.
    ! Loop over domains from west to east:
    do coord = 0, cart_size-1
    
      ! receive values from the domain on the west:
      call dom%TransferSlab( coord, 'from-west', 1, fluxx, status )
      IF_NOTOK_RETURN(status=1)
      
      ! is local domain in this column?
      if ( coord == cart_coord ) then

        ! loop over cells:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx

              ! look for positive winds
              if ( (afluxx(ix,iy,iz) > 0.0) .and. (afluxx(ix-1,iy,iz) >= 0.0) ) then

                ! first guess:
                cg = ( ch(ix,iy,iz) * d1(ix,iy,iz) &
                       + (fluxx(ix-1,iy,iz) - fluxx(ix,iy,iz)) * inv_volume(ix,iy,iz) ) &
                     * inv_d2(ix,iy,iz)

                ! inital concentration range:
                c1 = min( ch(ix-1,iy,iz), ch(ix,iy,iz) )
                c2 = max( ch(ix-1,iy,iz), ch(ix,iy,iz) )

                ! cg should remain within [c1,c2]:
                if ( cg < c1 .or. cg > c2 ) then
                  ! limit:
                  cg = min( c2, max( cg, c1 ) )
                  ! change flux at right side:
                  fluxx(ix,iy,iz) = fluxx(ix-1,iy,iz) + &
                                   (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz)) * volume(ix,iy,iz)
                end if
                
                ! store:
                cgues(ix,iy,iz) = cg
                
              end if ! positive winds

            end do ! ix
          end do  ! iy
        end do  ! iz
        
      end if  ! local domain in this cart dim
      
    end do  ! domains

    ! *
    
    ! Loop over domains from east to west:
    do coord = cart_size-1, 0, -1
    
      ! receive values from the domain on the east:
      call dom%TransferSlab( coord, 'from-east', 1, fluxx, status )
      IF_NOTOK_RETURN(status=1)
      
      ! is local domain in this column?
      if ( coord == cart_coord ) then

        ! now the limiting going down
        do iz = 1, nz
          do iy = 1, ny
            do ix = nx, 1, -1
            
              ! look for negative winds:
              if ( (afluxx(ix,iy,iz) < 0.0) .and. (afluxx(ix-1,iy,iz) < 0.0) ) then

                ! first guess:
                cg = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
                     (fluxx(ix-1,iy,iz) - fluxx(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)

                ! inital concentration range:
                c1 = min( ch(ix,iy,iz), ch(ix+1,iy,iz) )
                c2 = max( ch(ix,iy,iz), ch(ix+1,iy,iz) )

                ! cg should remain within [c1,c2]:
                if (cg < c1 .OR. cg > c2) then
                  ! limit:
                  cg = min(c2, max(cg, c1) )
                  ! change flux at left side:
                  fluxx(ix-1,iy,iz) = fluxx(ix,iy,iz) - &
                                   (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz))*volume(ix,iy,iz)
                end if

                ! store:
                cgues(ix,iy,iz) = cg

              end if  ! negative winds

            end do ! ix
          end do  ! iy
        end do  ! iz
        
      end if  ! local domain in this cart dim
      
    end do  ! domains

    ! *


    ! copy cgues to c
    ch(1:nx,1:ny,1:nz) = cgues(1:nx,1:ny,1:nz)

    ! clear:
    deallocate ( is_extreme, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate ( cgues, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate ( fluxx, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advecX


  ! ***----------------------------------------------


  subroutine advecY( nx,ny,nz, volume, inv_volume, d1, d2, k, &
                          ch, inv_d2, dta, status)

    use LE_Grid      , only : dom
    use LE_Meteo_Data, only : afluxy

    ! --- const ----------------------------------

    character(len=*), parameter :: rname = mname//'/advecY'

    ! number halo cells in concentration array:
    integer, parameter  ::  nh = 2

    ! --- in/out ---------------------------------

    integer, intent(in)           ::  nx, ny, nz
    real, intent(in), pointer     ::  volume(:,:,:)   ! (0:nx+1,0:ny+1,1:nz+1)
    real, intent(in)              ::  inv_volume(0:nx+1,0:ny+1,1:nz+1)
    real, intent(in)              ::  d1(nx,ny,nz)
    real, intent(in)              ::  d2(nx,ny,nz)
    real, intent(in)              ::  dta
    integer, intent(in)           ::  k
    real, intent(inout), pointer  ::  ch(:,:,:)       ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in)              ::  inv_d2(nx,ny,nz)
    integer, intent(out)          ::  status

    ! --- local ----------------------------------

    integer               ::  cart_size, cart_coord
    integer               ::  coord

    logical, pointer      ::  is_extreme(:,:,:)
    real, allocatable     ::  cgues(:,:,:)
    real, pointer         ::  fluxy(:,:,:)

    integer               ::  ix, iy, iz
    real                  ::  c1, c2
    real                  ::  u
    real                  ::  cour, fcour
    real                  ::  alpha
    real                  ::  cg
    
    ! --- begin ----------------------------------

    ! Info on domain layout in second dimension:
    call dom%GetDim( 2, status, cart_size=cart_size, cart_coord=cart_coord )
    IF_NOTOK_RETURN(status=1)

    ! storage for new concentrations:
    allocate ( cgues(nx,ny,nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for tracer flux:
    allocate ( fluxy(nx,0:ny,nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! initialize flux
    fluxy = 0.0

    ! storage for extremes including y halo:
    allocate ( is_extreme(nx,-1:ny+2,nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! init, this value remains in halo:
    is_extreme = .false.
    ! loop over x cells incl y halo:
    do iz = 1, nz
      do iy = 1, ny
        do ix = 1,nx
          ! concentrations in neighbour cells:
          c1 = ch(ix,iy-1,iz)
          c2 = ch(ix,iy+1,iz)
          ! concentrations here above/below max/min of neighbours?
          is_extreme(ix,iy,iz) = (ch(ix,iy,iz) >= max(c1,c2) .or. ch(ix,iy,iz) <= min(c1,c2))
        end do
      end do
    end do
    ! exchange halo cells:
    call dom%ExchangeHaloY( 2, is_extreme, status )
    IF_NOTOK_RETURN(status=1)

    ! compute the fluxes;
    ! loop over cell interfaces in y direction:
    do iz = 1, nz
      do iy = 0, ny
        do ix = 1, nx
          ! use volume flux through this interface for 'wind':
          u = afluxy(ix,iy,iz)
          ! courant number using volume left from the interface:
          !AJS: bug? should be 'upwind' courrant number with different volume given sign:
          !AJS:   cour = dta *     u(ix)  / volume(iy  )   , u >= 0.0
          !AJS:   cour = dta * abs(u(ix)) / volume(iy+1)   , u <  0.0
          ! FIX for domain decomposition:
          ! only skip first interface if at boundary:
          if ( cart_coord == 0 ) then
            cour = dta * abs(u) * inv_volume(ix,max(1,iy),iz)
          else
            cour = dta * abs(u) * inv_volume(ix,      iy ,iz)
          end if

          ! average concentration ((m3 tracer)/(m3 air) or (ug tracer)/(m3 air))
          ! NOTE: loop over edges 0,1,..,ny  ,  using extremes {-1,1},{0,2},...,{ny-1,ny+1} ;
          !   extremes are based on concentrations in surrouding cells,
          !   thus concentrations needed for -2,-1,0,..,ny+1,ny+2  ; strange, assymetric!
          fcour = (1.0-cour)/4.0
          alpha = 1.0
          if (u >=0.0) then
            !AJS: strange .. assymetric!
            if (is_extreme(ix,iy-1,iz)) alpha = max(1.5,1.2+0.6*cour)
            if (is_extreme(ix,iy+1,iz)) alpha = (1.75-0.45*cour)
            !!AJS: probably this was needed:
            !if (is_extreme(ix,iy  ,iz)) alpha = max(1.5,1.2+0.6*cour)
            !if (is_extreme(ix,iy+1,iz)) alpha = (1.75-0.45*cour)
            ! 'flux' through interface, but actually concentration at interface:
            fluxy(ix,iy,iz) = ch(ix,iy,iz) + (ch(ix,iy+1,iz) - ch(ix,iy-1,iz))*fcour*alpha
          else
            !AJS: strange .. assymetric!
            if (is_extreme(ix,iy+1,iz)) alpha = max(1.5,1.2+0.6*cour)
            if (is_extreme(ix,iy-1,iz)) alpha = (1.75-0.45*cour)
            !!AJS: probably this was needed:
            !if (is_extreme(ix,iy+1,iz)) alpha = max(1.5,1.2+0.6*cour)
            !if (is_extreme(ix,iy  ,iz)) alpha = (1.75-0.45*cour)
            ! 'flux' through interface, but actually concentration at interface:
            fluxy(ix,iy,iz) = ch(ix,iy+1,iz) + (ch(ix,iy,iz) - ch(ix,iy+2,iz))*fcour*alpha
          end if

          ! concentration range:
          c1 = min( ch(ix,iy,iz), ch(ix,iy+1,iz) )
          c2 = max( ch(ix,iy,iz), ch(ix,iy+1,iz) )
          ! bound cell face concentration:
          fluxy(ix,iy,iz) = max( c1, min( fluxy(ix,iy,iz), c2 ) )
          ! make a flux out of it ....
          fluxy(ix,iy,iz) = u * dta * fluxy(ix,iy,iz) * d0

        end do ! ix
      end do ! iy (interfaces)
    end do ! iz
    
    ! *

    ! First handle exceptions: inflow or outflow only cells.
    ! Loop over domains from south to north:
    do coord = 0, cart_size-1
    
      ! receive values from the domain on the west:
      call dom%TransferSlab( coord, 'from-south', 1, fluxy, status )
      IF_NOTOK_RETURN(status=1)
      
      ! is local domain in this column?
      if ( coord == cart_coord ) then
      
        ! loop over cells:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx

              !~ inflow only?
              if ( (afluxy(ix,iy,iz) <= 0.0) .and. (afluxy(ix,iy-1,iz) >= 0.0) ) then
                ! (re)set tracer fluxes:
                fluxy(ix,iy-1,iz) = dta * afluxy(ix,iy-1,iz) * ch(ix,iy-1,iz) * d0
                fluxy(ix,iy  ,iz) = dta * afluxy(ix,iy  ,iz) * ch(ix,iy+1,iz) * d0
                ! update concentartion:
                cgues(ix,iy,iz)   = ( ch(ix,iy,iz) * d1(ix,iy,iz) + &
                                      (fluxy(ix,iy-1,iz) - fluxy(ix,iy,iz)) * inv_volume(ix,iy,iz) ) &
                                    * inv_d2(ix,iy,iz)

              !~ outflow only?
              else if ( (afluxy(ix,iy,iz) >= 0.0) .and. (afluxy(ix,iy-1,iz) <= 0.0) ) then
                ! (re)set tracer fluxes:
                fluxy(ix,iy-1,iz) = dta * afluxy(ix,iy-1,iz) * ch(ix,iy,iz) * d0
                fluxy(ix,iy  ,iz) = dta * afluxy(ix,iy  ,iz) * ch(ix,iy,iz) * d0
                ! update concentartion:
                cgues(ix,iy,iz)   = ( ch(ix,iy,iz)*d1(ix,iy,iz) + &
                                      (fluxy(ix,iy-1,iz) - fluxy(ix,iy,iz)) * inv_volume(ix,iy,iz) ) &
                                    * inv_d2(ix,iy,iz)
              end if ! inflow or outflow only

            end do  ! ix
          end do ! iy
        enddo ! iz
        
      end if  ! local domain in this cart dim
    
      ! fluxes on left side might have changed, send to previous domain:
      call dom%TransferSlab( coord, 'to-south', 1, fluxy, status )
      IF_NOTOK_RETURN(status=1)
      
    end do  ! domains
    
    ! *
    
    ! Limit fluxes, first going up.
    ! Loop over domains from south to north:
    do coord = 0, cart_size-1
    
      ! receive values from the domain on the west:
      call dom%TransferSlab( coord, 'from-south', 1, fluxy, status )
      IF_NOTOK_RETURN(status=1)
      
      ! is local domain in this column?
      if ( coord == cart_coord ) then

        ! loop over cells:
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx

              ! look for positive winds:
              if (afluxy(ix,iy,iz) > 0.0 .AND. afluxy(ix,iy-1,iz) >=0.0) then

                ! first guess:
                cg = ( ch(ix,iy,iz) * d1(ix,iy,iz) &
                       + (fluxy(ix,iy-1,iz) - fluxy(ix,iy,iz)) / volume(ix,iy,iz) ) &
                     * inv_d2(ix,iy,iz)

                ! inital concentration range:
                c1 = min( ch(ix,iy-1,iz), ch(ix,iy,iz) )
                c2 = max( ch(ix,iy-1,iz), ch(ix,iy,iz) )

                ! cg should remain within [c1,c2]:
                if ( cg < c1 .or. cg > c2 ) then
                  ! limit:
                  cg = min( c2, max( cg, c1 ) )
                  ! change flux at right side:
                  fluxy(ix,iy,iz) = fluxy(ix,iy-1,iz) + &
                                    (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz)) * volume(ix,iy,iz)
                end if
                
                ! store:
                cgues(ix,iy,iz) = cg
                
              end if ! positive winds

            end do ! ix
          end do  ! iy
        end do  ! iz
        
      end if  ! local domain in this row
      
    end do  ! domain rows
    
    ! *
    
    ! Loop over domains from east to west:
    do coord = cart_size-1, 0, -1
    
      ! receive values from the domain on the north:
      call dom%TransferSlab( coord, 'from-north', 1, fluxy, status )
      IF_NOTOK_RETURN(status=1)
      
      ! is local domain in this column?
      if ( coord == cart_coord ) then

        ! now the limiting going down
        do iz = 1, nz
          do iy = ny, 1, -1
            do ix = 1, nx

              ! look for negative winds:
              if ( (afluxy(ix,iy,iz) < 0.0) .and. (afluxy(ix,iy-1,iz) < 0.0) ) then

                ! first guess:
                cg = (ch(ix,iy,iz)*d1(ix,iy,iz) + &
                     (fluxy(ix,iy-1,iz) - fluxy(ix,iy,iz)) * inv_volume(ix,iy,iz) ) * inv_d2(ix,iy,iz)

                ! inital concentration range:
                c1 = min( ch(ix,iy,iz), ch(ix,iy+1,iz) )
                c2 = max( ch(ix,iy,iz), ch(ix,iy+1,iz) )

                ! cg should remain within [c1,c2]:
                if (cg < c1 .OR. cg > c2) then
                  ! limit:
                  cg = min(c2, max(cg, c1) )
                  ! change flux at left side:
                  fluxy(ix,iy-1,iz) = fluxy(ix,iy,iz) - &
                                    (ch(ix,iy,iz)*d1(ix,iy,iz)-cg*d2(ix,iy,iz))*volume(ix,iy,iz)
                end if

                ! store:
                cgues(ix,iy,iz) = cg

              end if  ! negative winds

            end do ! ix
          end do  ! iy
        end do  ! iz
        
      end if  ! local domain in this cart dim
      
    end do  ! domains
    
    ! *


    ! copy cgues to c
    ch(1:nx,1:ny,1:nz) = cgues(1:nx,1:ny,1:nz)

    ! clear:
    deallocate ( is_extreme, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate ( cgues, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate ( fluxy, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine advecY


  ! ***------------------------------------------------------


  subroutine advecZ( nx, ny, nz, volume, inv_volume, d1, d2, k, ch, inv_d2, dta, status)

    use LE_Meteo_Data, only : afluxz

    ! --- const ----------------------------------

    character(len=*), parameter :: rname = mname//'/advecZ'

    integer, parameter  ::  nh = 2

    ! --- in/out ---------------------------------

    integer, intent(in)           ::  nx, ny, nz
    real, intent(in), pointer     ::  volume(:,:,:)     ! (0:nx+1,0:ny+1,1:nz+1)
    real, intent(in)              ::  inv_volume(0:nx+1,0:ny+1,1:nz+1)
    real, intent(in)              ::  d1(nx,ny,nz)
    real, intent(in)              ::  d2(nx,ny,nz)
    real, intent(in)              ::  dta
    integer, intent(in)           ::  k
    real, intent(inout), pointer  ::  ch(:,:,:)         ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in)              ::  inv_d2(nx,ny,nz)
    integer, intent(out)          ::  status

    ! --- local ----------------------------------

    real, allocatable     ::  cgues(:,:,:)
    real, pointer         ::  fluxz(:,:,:)

    integer               ::  ix, iy, iz
    real                  ::  c1, c2
    real                  ::  u
    real                  ::  cour, fcour
    real                  ::  alpha
    real                  ::  cg

    ! --- begin ----------------------------------

    ! storage for new concentrations:
    allocate ( cgues(nx,ny,nz), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for tracer flux:
    allocate ( fluxz(nx,ny,0:nz), stat=status )
    IF_NOTOK_RETURN(status=1)
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


    ! copy cgues to c
    ch(1:nx,1:ny,1:nz) = cgues(1:nx,1:ny,1:nz)

    ! clear:
    deallocate ( cgues )
    deallocate ( fluxz )

    ! ok
    status = 0

  end subroutine advecZ


  ! ----------------------------------------------------------


  ! check courant condition for dt. If necessary,
  ! compute how many substeps are needed to satisfy
  ! the courant condition

  subroutine courant_a( nx, ny,nz, inv_volume, dt, dta, nstep_a, store )

    use LE_Meteo_Data, only : afluxx, afluxy, afluxz
    use LE_Meteo_Data, only : cflx, cfly, cflz  ! (nx,ny,nz)
    use dims, only : runF, outF
    use LE_Logging, only : ident2

    ! --- in/out ---------------------------------

    integer, intent(in)   ::  nx, ny, nz
    real, intent(in)      ::  inv_volume(0:nx+1,0:ny+1,1:nz+1)
    real, intent(in)      ::  dt

    ! actual time steps and number of steps for advection:
    real, intent(out)     ::  dta
    integer, intent(out)  ::  nstep_a

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
  end subroutine courant_a
  
  ! *

#ifdef with_hdiff
  subroutine courant_d( nx, ny,nz, inv_volume, dt, dtd, nstep_d, store )

    use dims, only : khx, khy
    use dims, only : runF, outF
    use LE_Logging, only : ident2

    ! --- in/out ---------------------------------

    integer, intent(in)   ::  nx, ny, nz
    real, intent(in)      ::  inv_volume(0:nx+1,0:ny+1,1:nz+1)
    real, intent(in)      ::  dt

    ! actual time steps and number of steps for diffusion:
    integer, intent(out)  ::  nstep_d
    real, intent(out)     ::  dtd

    ! --- local ----------------------------------

    real    :: cour
    real    :: cour_max
    integer :: i,j,k

    ! --- begin ----------------------------------
    ! check time step limit for diffusion
    cour_max = 0.0

    do i=0,nx
    do j=1,ny
    do k=1,nz
       cour = abs( khx(i,j,k)*dt)/runF%dx(j) * inv_volume(i,j,k)
       cour_max = max(cour_max, cour)
    enddo
    enddo
    enddo

    do i=1,nx
    do j=0,ny
    do k=1,nz
       cour = abs( khy(i,j,k)*dt)/runF%dy * inv_volume(i,j,k)
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

  end subroutine courant_d
#endif


end module LE_Advec

