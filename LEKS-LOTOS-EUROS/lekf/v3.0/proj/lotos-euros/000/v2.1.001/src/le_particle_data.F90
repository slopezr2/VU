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

module LE_Particle_Data

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out --------------------------------

  private

#ifdef with_m7
  public  ::  im7mode
#endif
  public  ::  partsize
  public  ::  slipcor

  public  ::  LE_Particle_Data_Init, LE_Particle_Data_Done, LE_Particle_Data_Update
  public  ::  partsize_modes


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Particle_Data'


  ! --- var --------------------------------------

#ifdef with_m7
  ! Original M7 mode: 1..7, or -999 for other tracers:
  integer, allocatable  ::  im7mode(:)  ! (nspec)
#endif
  
  ! average particle size per cell and component:
  real, allocatable     ::  partsize(:,:,:,:)  ! (nx,ny,nz,n_aerosol)
  
  real   :: partsize_modes(5)
  

  ! slip correction factor
  real, allocatable     ::  slipcor(:,:,:,:)     ! (nx,ny,nz,n_aerosol)



contains



  ! ========================================================================
  ! ===
  ! ===  particle size stuff
  ! ===
  ! ========================================================================


  subroutine LE_Particle_Data_Init( status )
  
#ifdef with_m7    
    use Indices, only : nspec, specmode, tracer_is_soluble
    use Indices, only : n_m7, ispecs_m7
    use Indices, only : AEROSOL_NUCL_MODE, AEROSOL_AITKEN_MODE, AEROSOL_ACCUM_MODE
    use Indices, only : specmode, NO_AEROSOL_MODE
    use Indices, only : AEROSOL_COARSE_MODE
#endif

    ! --- in/out ---------------------------------

    integer, intent(out)        ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_Particle_Data_Init'

    ! --- local ------------------------------------------
    
    integer   ::  k
    integer   ::  ispec

    ! --- begin ------------------------------------------

#ifdef with_m7
    ! original M7 mode numbers for each spec:
    allocate( im7mode(nspec) )
    ! dummy by default:
    im7mode(:) = -999
    ! any M7 tracers ?
    if ( n_m7 > 0 ) then
      ! loop over all M7 tracers:
      do k = 1, n_m7
        ! global index:
        ispec = ispecs_m7(k)
        ! switch:
        select case ( specmode(ispec) )
          case ( NO_AEROSOL_MODE )
            ! for example H2SO4, has no M7 mode
          case ( AEROSOL_NUCL_MODE )
            im7mode(ispec) = 1
          case ( AEROSOL_AITKEN_MODE )
            if ( tracer_is_soluble(ispec) )then
              im7mode(ispec) = 2
            else
              im7mode(ispec) = 5
            end if
          case ( AEROSOL_ACCUM_MODE )
            if ( tracer_is_soluble(ispec) )then
              im7mode(ispec) = 3
            else
              im7mode(ispec) = 6
            end if
          case ( AEROSOL_COARSE_MODE )
            if ( tracer_is_soluble(ispec) )then
              im7mode(ispec) = 4
            else
              im7mode(ispec) = 7
            end if
          case default
            write (gol,'("unsupported specmode ",i6," for M7 tracer")') specmode(ispec); call goErr
            TRACEBACK;status=1;return
        end select
      end do  ! M7 tracers
    end if  ! any M7 tracers
#endif

    ! init size stuff:
    call LE_Particle_Size_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! init slipcor stuff:
    call SlipCor_Init( status )
    IF_NOTOK_RETURN(status=1)

     ! ok
    status = 0

  end subroutine LE_Particle_Data_Init


  ! ***


  subroutine LE_Particle_Data_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Particle_Data_Done'

    ! --- begin ----------------------------------

    ! done with size stuff:
    call LE_Particle_Size_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with slipcor stuff:
    call SlipCor_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Particle_Data_Done

  ! *

  subroutine LE_Particle_Data_Update( c, status )

#ifdef with_m7
    use LE_M7, only : LE_M7_Updatemodes
#endif
    ! --- in/out ---------------------------------

    real, intent(inout)       ::  c(:,:,:,:)   ! (nx,ny,nz,nspec)
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Particle_Data_Update'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef with_m7
    ! compute new aerosol wet-radius etc:
    call LE_M7_Updatemodes( c, status, warning=.true. )
    IF_NOTOK_RETURN(status=1)
#endif

    ! update size stuff, for M7 tracers this depends on wet-radius:
    call LE_Particle_Size_Update( c, status )
    IF_NOTOK_RETURN(status=1)

    ! update slip correction, this uses the updated particle size!
    call SlipCor_Update( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Particle_Data_Update



  ! ========================================================================
  ! ===
  ! ===  particle size stuff
  ! ===
  ! ========================================================================


  subroutine LE_Particle_Size_Init( status )

    use Dims       , only : nx, ny, nz, nspec
    use Indices    , only : specname
    use Indices    , only : n_aerosol, ispecs_aerosol, ispec2aerosol
    use Indices    , only : rhopart, rhopart_default, rhopart_dust
    use Indices    , only : ispec_dust_ff, ispec_dust_f
    use Indices    , only : ispec_dust_ccc, ispec_dust_cc, ispec_dust_c
#ifdef with_m7
    use Indices    , only : tracer_is_m7
    use Indices    , only : AEROSOL_NUCL_MODE, AEROSOL_AITKEN_MODE, AEROSOL_ACCUM_MODE
#endif    
#ifdef with_pollen
    use Indices    , only : AEROSOL_POLLEN_MODE
    use Indices    , only : rhopart_pollen
    use Indices    , only : ispec_pol_b, ispec_pol_g, ispec_pol_o
#endif
    use Indices    , only : specmode, NO_AEROSOL_MODE, AEROSOL_FINE_MODES, AEROSOL_COARSE_MODE, AEROSOL_All_MODES
    use Indices    , only : AEROSOL_FF_MODES, AEROSOL_CC_MODE, AEROSOL_CCC_MODE
    use Indices    , only : AEROSOL_ULTRA_FINE_MODE, AEROSOL_ULTRA_FINE_FINE_MODE, AEROSOL_FINE_MEDIUM_MODE, AEROSOL_MEDIUM_COARSE_MODE
    use Indices    , only : i_Na_c

    ! --- in/out ---------------------------------

    integer, intent(out)        ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_Particle_Size_Init'

    ! Effective particle size, depends on location and component.
    ! This is necessary since the aerosol mass is defined for a
    ! size bin [0,2.5] etc, while we do not know the average size.
    ! A representation with size and number as used by M7 would
    ! not have this problem.
    ! For most tracers the same particle size is assumed everywhere,
    ! only for sea-salt a smaller size is assumed over land.
    
       
   
	! Default: fine and coarse, 0.7e-6 and 8e-6 m diameter (wet/dry), used for deposition
    ! and radiation
    ! for dust and sea salt: additional size classes

  	real, parameter   ::  partsize_ff   = 0.33e-6
    real, parameter   ::  partsize_fine   = 0.7e-6
    real, parameter   ::  partsize_ccc = 3.0e-6
    real, parameter   ::  partsize_cc = 5.0e-6
    real, parameter   ::  partsize_coarse = 8.0e-6

#ifdef with_pollen        
    real, parameter   ::  partsize_pollen_birch = 22.0e-6
    real, parameter   ::  partsize_pollen_grass = 32.0e-6
    real, parameter   ::  partsize_pollen_olive = 28.0e-6
#endif
    
    ! Values for tracers used in PPM with OPS
    real, parameter   ::  partsize_ultra_fine      = 0.7e-6
    real, parameter   ::  partsize_ultra_fine_fine = 0.7e-6
    real, parameter   ::  partsize_fine_medium     = 8.0e-6
    real, parameter   ::  partsize_medium_coarse   = 8.0e-6

    ! --- local ------------------------------------------

    integer   ::  ispec, i_aerosol, k

    ! --- begin ------------------------------------------

    ! storage:
    allocate( partsize(nx,ny,nz,n_aerosol) )
    allocate( rhopart(nspec) )
    rhopart = -999.0

    partsize_modes(1)=partsize_ff
    partsize_modes(2)=partsize_fine
    partsize_modes(3)=partsize_ccc
    partsize_modes(4)=partsize_cc
    partsize_modes(5)=partsize_coarse

    ! default no size at all ...
    partsize = 0.0

    ! loop over aerosols:
    do i_aerosol = 1, n_aerosol
      ! global index:
      ispec = ispecs_aerosol(i_aerosol)
        ! fine/coarse modes ; switch on mode:
        select case ( specmode(ispec) )
          case ( AEROSOL_FINE_MODES )
            partsize(:,:,:,i_aerosol) = partsize_fine  ! m
          case ( AEROSOL_FF_MODES )
            partsize(:,:,:,i_aerosol) = partsize_ff  ! m
          case ( AEROSOL_CCC_MODE )
            partsize(:,:,:,i_aerosol) = partsize_ccc   ! m
          case ( AEROSOL_CC_MODE )
            partsize(:,:,:,i_aerosol) = partsize_cc   ! m
          case ( AEROSOL_COARSE_MODE )
            partsize(:,:,:,i_aerosol) = partsize_coarse   ! m
#ifdef with_m7
            ! Overwrite with size for smallest coarse mode
            partsize(:,:,:,i_aerosol) = partsize_ccc   ! m
#endif

#ifdef with_pollen
          case ( AEROSOL_POLLEN_MODE )
            select case ( trim(specname(ispec) ))
              case ('pol_b' )
                partsize(:,:,:,i_aerosol) = partsize_pollen_birch
              case ('pol_g' )
                partsize(:,:,:,i_aerosol) = partsize_pollen_grass
              case ('pol_o' )
                partsize(:,:,:,i_aerosol) = partsize_pollen_olive
              case default
                write( gol, '("No particle size defined for pollen tracer ",i3," (",a,")")') &
                            ispec, trim(specname(ispec)); call goErr
                TRACEBACK;status=1;return
            end select
#endif          
          case (AEROSOL_ULTRA_FINE_MODE )
            partsize(:,:,:,i_aerosol) = partsize_ultra_fine ! m
          case (AEROSOL_ULTRA_FINE_FINE_MODE )
            partsize(:,:,:,i_aerosol) = partsize_ultra_fine_fine ! m
          case (AEROSOL_FINE_MEDIUM_MODE )
            partsize(:,:,:,i_aerosol) = partsize_fine_medium ! m
          case (AEROSOL_MEDIUM_COARSE_MODE )
            partsize(:,:,:,i_aerosol) = partsize_medium_coarse ! m  
            
#ifdef with_m7      
          ! Particle size based on M7 is calculated according to actual particle mass 
          ! and number concentration in dry deposition/sedimentation section
          case ( NO_AEROSOL_MODE )
            ! for example H2SO4
          case ( AEROSOL_NUCL_MODE, AEROSOL_AITKEN_MODE, AEROSOL_ACCUM_MODE )
            partsize(:,:,:,i_aerosol) = partsize_fine  ! m
#endif
          ! undefined particle sizes used for accumulated tracers (tpm10):
          case ( AEROSOL_All_MODES )
            ! undefined

          case default
            write (gol,'("no particle mode for aerosol tracer ",i3," (",a,")")') &
                     ispec, trim(specname(ispec)); call goErr
            TRACEBACK; status=1; return
        end select
        
        ! fill in aerosol density
        select case (ispec )
          case ( ispec_dust_ff, ispec_dust_f, ispec_dust_ccc, ispec_dust_cc, ispec_dust_c )
            rhopart(ispec) = rhopart_dust
#ifdef with_pollen            
          case ( ispec_pol_b, ispec_pol_g, ispec_pol_o )
            rhopart(ispec) = rhopart_pollen
#endif          
          case default
            rhopart(ispec) = rhopart_default
        end select
      
    end do  ! specs
    
    ! No need for this when more size classes are takeninto account    
    ! Special: sodium particles.
    ! Actaually sea salt particles, sodium as tracer.
    ! Coarse mode over land is smaller.
    !ispec = i_Na_c
    !if ( ispec > 0.0 ) then
    !  ! mapping:
    !  i_aerosol = ispec2aerosol(ispec)
    !  ! fill lowest layer:
    !  where ( seafraction(:,:) > 0.5 )
    !    partsize(:,:,1,i_aerosol) = partsize_coarse  ! m  ; sea: large particles
    !  elsewhere
    !    partsize(:,:,1,i_aerosol) = 4.0e-6           ! m  ; land: smaller particles
    !  end where
    !  ! copy to higher layers:
    !  do k = 2, nz
    !    partsize(:,:,k,i_aerosol) = partsize(:,:,1,i_aerosol)
    !  end do
    !end if
    
    ! ok
    status = 0

  end subroutine LE_Particle_Size_Init


  ! ***


  subroutine LE_Particle_Size_Done( status )

    use Indices, only : rhopart
    
    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Particle_Size_Done'

    ! --- begin ----------------------------------

    ! clear:
    deallocate( partsize )
    deallocate( rhopart )

    ! ok
    status = 0

  end subroutine LE_Particle_Size_Done

  ! *

  subroutine LE_Particle_Size_Update( c, status )

#ifdef with_m7
    use Indices   , only : specmode, tracer_is_soluble, tracer_is_numberconc
    use Indices   , only : n_m7, ispecs_m7, ispec2aerosol, tracer_is_m7
    use LE_M7     , only : sigmaln
    use LE_M7_Data, only : rwetm7modes
#endif
    use dims      , only : nx, ny, nz, nspec

    implicit none

    ! --- in/out ---------------------------------

    real, intent(inout)         ::  c(nx,ny,nz,nspec)
    integer, intent(out)        ::  status

    ! --- const -------------------------------

    character(len=*), parameter  ::  rname = mname//'/LE_Particle_Size_Update'

    ! --- local ----------------------------------

#ifdef with_m7
    integer     ::  k, ispec
    integer     ::  im7
    integer     ::  i_aerosol
#endif

    ! --- begin ----------------------------------

#ifdef with_m7
    ! update M7 particle sizes ?
    ! applied to complete array, so check if needed:
    if ( n_m7 > 0 ) then

      ! loop over all M7 tracers:
      do k = 1, n_m7
        ! global index:
        ispec = ispecs_m7(k)
        ! M7 mode number:
        im7 = im7mode(ispec)
        ! skip non-aerosols (H2SO4):
        if ( im7 <= 0 ) cycle
        ! index in 'n_aerosol' arrays:
        i_aerosol = ispec2aerosol(ispec)
        ! init particle size:
        partsize(:,:,:,i_aerosol) = 2.0 * rwetm7modes(:,:,:,im7) * 1e-2 ! cm->m,  radius->diam
        ! number concentration or mass ?
        if ( tracer_is_numberconc(ispec) ) then
          partsize(:,:,:,i_aerosol) = partsize(:,:,:,i_aerosol) * exp(       sigmaln(im7)**2 )
        else
          partsize(:,:,:,i_aerosol) = partsize(:,:,:,i_aerosol) * exp( 3.0 * sigmaln(im7)**2 )
        end if
      end do

    end if  ! any M7 tracers
#endif

    ! ok
    status = 0

  end subroutine LE_Particle_Size_Update



  ! ========================================================================
  ! ===
  ! ===  particle size stuff
  ! ===
  ! ========================================================================


  subroutine SlipCor_Init( status )

    use Dims   , only : nx, ny, nz
    use Indices, only : n_aerosol

    ! --- in/out ---------------------------------

    integer, intent(out)        ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/SlipCor_Init'

    ! --- local ------------------------------------------

    ! --- begin ------------------------------------------

    ! storage:
    allocate( slipcor(nx,ny,nz,n_aerosol) )
    ! dummy values:
    slipcor = 0.0

     ! ok
    status = 0

  end subroutine SlipCor_Init


  ! ***


  subroutine SlipCor_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/SlipCor_Done'

    ! --- begin ----------------------------------

    ! clear:
    deallocate( slipcor )

    ! ok
    status = 0

  end subroutine SlipCor_Done

  ! *

  subroutine SlipCor_Update( status )

    use LE_Meteo_Data, only : freepathlen
    use indices      , only : n_aerosol, ispecs_aerosol
    use indices      , only : specmode
    !use Indices      , only : AEROSOL_FINE_MODES, AEROSOL_COARSE_MODE
    use JAQL_drydeposition, only : Slip_Correction_Factor

    use indices      , only : ispec_dust_c

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/SlipCor_Update'

    ! --- local ----------------------------------

    integer :: i_aerosol

    ! --- begin ----------------------------------

    ! loop over aerosol tracers:
    do i_aerosol = 1, n_aerosol
      ! trap undefined (zero) particle sizes;
      ! in this case the aerosol load is probably zero (M7):
      where ( partsize(:,:,:,i_aerosol) > 0.0 )
        slipcor(:,:,:,i_aerosol) = Slip_Correction_Factor( freepathlen, partsize(:,:,:,i_aerosol) )
      elsewhere
        slipcor(:,:,:,i_aerosol) = 0.0
      end where
     
    end do  ! specs

    ! ok
    status = 0

  end subroutine SlipCor_Update



end module LE_Particle_Data

