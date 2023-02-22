!###############################################################################
!
! Sedimentation of aerosols.
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

module LE_Sedim

  use GO, only : gol, goPr, goErr
  
  implicit none


  ! --- in/out --------------------------------
  
  private
  
  public  ::  slinnfac
  
  public  ::  LE_Sedim_Init, LE_Sedim_Done
  public  ::  LE_Sedim_Apply


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Sedim'
  

  ! --- var --------------------------------------

  ! factor for sedimentation velocity, 1.0 for non-M7 tracers:
  real, allocatable   ::  slinnfac(:)     ! (nspec)
  
  

contains



  ! ========================================================================
  

  subroutine LE_Sedim_Init( status )
  
    use Dims            , only : nx, ny, nz
    use Indices         , only : nspec
    use Indices         , only : n_aerosol, tracer_is_numberconc
    use LE_Data, only : LE_Data_Enable
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Sedim_Init'
    
    ! --- local ------------------------------------------

    integer   ::  ispec

    ! --- begin ------------------------------------------

    ! factor used for computation of sedimentation velocity
    allocate( slinnfac(nspec) )
    ! loop over all tracers:
    do ispec = 1, nspec
         slinnfac(ispec) = 1.0
    end do

    ! enable data:
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'p', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Sedim_Init
  
  
  ! ***
  
  
  subroutine LE_Sedim_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Sedim_Done'
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( slinnfac )

    ! ok
    status = 0
    
  end subroutine LE_Sedim_Done
  
  
  ! ***

  
  subroutine LE_Sedim_Apply( c, bud, update_bud, dt_sec, status )

    use LE_Logging        , only : ident2
    use dims              , only : outF
    use Dims              , only : nx, ny, nz, nspec
    use LE_Meteo_Data     , only : freepathlen
    use LE_Meteo_Data     , only : viscos
    use Indices           , only : n_aerosol, ispecs_aerosol
    use Indices           , only : specmode, AEROSOL_FINE_MODES, AEROSOL_COARSE_MODE
    use Indices           , only : specmode, AEROSOL_FF_MODES, AEROSOL_CC_MODE,AEROSOL_CCC_MODE
    use Indices           , only : AEROSOL_ULTRA_FINE_MODE, AEROSOL_ULTRA_FINE_FINE_MODE, AEROSOL_FINE_MEDIUM_MODE, AEROSOL_MEDIUM_COARSE_MODE
    use Indices           , only : specname
    use Indices           , only : rhopart
    use LE_Particle_Data  , only : partsize
    use LE_Particle_Data  , only : slipcor
    use LE_Particle_Data  , only : LE_Particle_Data_Update
    use JAQL_drydeposition, only : Sedimentation_Velocity
    use LE_Budget_DryDepos, only : T_DryDepos_Budget
    ! point to meteo data:
    use LE_Data      , only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------------
    
    real, intent(inout)                      ::  c(nx,ny,nz,nspec)
    type(T_DryDepos_Budget), intent(inout)   ::  bud                ! budgets, a.o. per hour
    logical, intent(in)                      ::  update_bud
    real, intent(in)                         ::  dt_sec
    integer, intent(out)                     ::  status
    
    
    ! --- const ----------------------------------
    character(len=*), parameter     :: rname = mname//'/LE_Sedim_Apply'
    
    ! --- local ----------------------------------

    integer             ::  i_aerosol
    integer             ::  ispec
    real, allocatable   ::  vs(:,:,:)   ! (nx,ny,nz) sedimentation velocity
    ! meteo data:
    real, pointer       ::  temp(:,:,:)   ! (lon,lat,1)    
    real, pointer       ::  press(:,:,:)   ! (lon,lat,1)        
    
    ! --- begin ----------------------------------

    call LE_Data_GetPointer( 't', temp, status, check_units ='K')          
    IF_NOTOK_RETURN(status=1)  
    call LE_Data_GetPointer( 'p', press, status, check_units ='Pa')        
    IF_NOTOK_RETURN(status=1)

    !! info ...
    !if (.NOT.outF%suppress) then
    !  write (gol,'(a,"<sedim>")') ident2; call goPr
    !end if
    
    ! storage:
    allocate( vs(nx,ny,nz) )
    
    !>>> now called in LE_Driver before call to this routine:
    ! update particle sizes, slip correction factors etc (if necessary):
    !call LE_Particle_Data_Update( c, status )
    !IF_NOTOK_RETURN(status=1)
    !<<<
    
    ! loop over aerosol tracers:
    do i_aerosol = 1, n_aerosol
      ! global index:
      ispec = ispecs_aerosol(i_aerosol)

      ! switch:
      select case ( specmode(ispec) )

        ! fine modes:
        case ( AEROSOL_FINE_MODES, AEROSOL_FF_MODES, &
               AEROSOL_ULTRA_FINE_MODE, AEROSOL_ULTRA_FINE_FINE_MODE)

          ! no sedimentation of fine mode particles, only coarse ..

        ! coarse mode:
        case ( AEROSOL_COARSE_MODE, AEROSOL_CC_MODE, AEROSOL_CCC_MODE, &
               AEROSOL_FINE_MEDIUM_MODE, AEROSOL_MEDIUM_COARSE_MODE )

          ! sedimentation velocity:
          vs = slinnfac(ispec) * Sedimentation_Velocity( rhopart(ispec), &
                                                         partsize(:,:,:,i_aerosol), &
                                                         slipcor (:,:,:,i_aerosol), &
                                                         viscos )
          ! apply sedimentation to tracer:
          call LE_Sedim_Apply_comp( c(:,:,:,ispec), bud, update_bud, &
                                       vs, dt_sec, ispec, status )              
          IF_NOTOK_RETURN(status=1)
          
        ! unknown ...
        case default
          write (gol,'("unexpected specmode ",i6," for aaerosol tracer with ispec ",i6)') specmode(ispec), ispec; call goErr
          TRACEBACK;status=1; return
          
      end select

    end do  ! aerosols
                                                               
    ! clear:
    deallocate( vs )

    ! ok
    status = 0

  end subroutine LE_Sedim_Apply
  
  
  ! ***


  ! apply sedimentation
  
  subroutine LE_Sedim_Apply_comp( c, bud, update_bud, vsed, dt_sec, ispec, status )

    use Dims, only : nx, ny, nz, nspec
    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer
    
    use LE_Budget_DryDepos, only : T_DryDepos_Budget, ex_drydepo
    use LE_Landuse_Data   , only : nlu
    
    real, intent(inout)   ::  c(nx,ny,nz)
    type(T_DryDepos_Budget), intent(inout)  ::  bud                ! budgets, a.o. per hour
    logical, intent(in)                     ::  update_bud
    real, intent(in)                        ::  vsed(nx,ny,nz)
    real, intent(in)                        ::  dt_sec   ! s
    integer, intent(in)                     ::  ispec
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------
    character(len=*), parameter     :: rname = mname//'/LE_Sedim_Apply_comp'
    
    ! --- local --------------------

    integer   ::  i,j,k, nsteps
    real      ::  courmax, flux(nz+1)
    integer   ::  status_par    
    integer   ::  ilu
    
    integer   ::  iz
    real      ::  dh(nx,ny,nz)
    
    ! meteo data:
    real, pointer               ::  h_m(:,:,:)   ! (lon,lat,1)

    ! --- begin ----------------------------
    
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m')    
    IF_NOTOK_RETURN(status=1)

    ! Calculate cell heights
    dh(:,:,1) = h_m(:,:,1)
    do iz = 2, nz
      dh(:,:,iz) = h_m(:,:,iz) - h_m(:,:,iz-1)
    end do
    
    ! AJS: why a global courrant number ?
    !      Should be sufficient to do this per column.
    ! compute the courant condition for the mass in the mode
    !                 m/s    min s/min/  km    / (m/km)

    !                    m/s          s    /   m 
    courmax = maxval( vsed(:,:,:) * dt_sec / dh(:,:,:) )  ! 1
    nsteps = int(courmax) + 1
    
    ! init status used for parrallel loop:
    status_par = 0
    
    !$OMP parallel &
    !$OMP default( none ) &
    !$OMP shared( nsteps ) &    
    !$OMP shared( nx, ny, nz , nlu) &    
    !$OMP shared( bud, update_bud ) &
    !$OMP shared( c ) &    
    !$OMP shared( vsed, dt_sec ) &
    !$OMP shared( dh ) &
    !$OMP shared( ispec ) &
    !$OMP shared( gol ) &
    !$OMP private( i, j, ilu ) &
    !$OMP private( flux ) &
    !$OMP private( status ) &
    !$OMP reduction( + : status_par )
    !$OMP do 
    ! time integration
    do i=1,nx
    do j=1,ny


      ! loop over intermediate time steps:
      do k=1,nsteps
        flux(1)=0.0
        ! is the flux out of this layer z to z-1
        flux(2:nz) = dt_sec/real(nsteps) * vsed(i,j,2:nz) * c(i,j,2:nz)
        flux(nz+1) = 0.0
        
        ! update concentration:
        c(i,j,1:nz) = c(i,j,1:nz) + ( flux(2:nz+1) - flux(1:nz) )/dh(i,j,1:nz)

        
        ! update budgets ?
        if ( update_bud ) then
          ! loop over landuse classes:
          do ilu = 1, nlu
            ! landuse independents
            ! Fluxes negative to budgets because depositions are described to be negative
			! In this case results are not influenced, because sedimentation is not calculated in surface layer
			! Base settings: Sedimentation in surface layer is calculated in dry-deposition routine (according to Seinfeld Pandis)
            bud%ex_hour_lu(i,j,ispec,ilu, ex_drydepo ) = bud%ex_hour_lu(i,j,ispec,ilu, ex_drydepo ) - flux(1)
          end do
          ! actual sedimentation
          bud%ex_hour(i,j,ispec,ex_drydepo ) = bud%ex_hour(i,j,ispec,ex_drydepo ) - flux(1)
        end if
        
      enddo
    enddo
    enddo
    !$OMP end do
    !$OMP end parallel

    ! check ...   
    if ( status_par /= 0 ) then
      write (gol,'("Non-zero status returned from applicaiton of sedimentation")'); call goErr
      TRACEBACK; status=1;return
    endif
    
    ! ok
    status = 0

  end subroutine  LE_Sedim_Apply_comp


end module LE_Sedim

