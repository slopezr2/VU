!###############################################################################
!
! NAME
!
!   LEKF_Driver  -  kf interface to lotos-euros model
!
! CODE TREE
!
!   driver
!
!     output
!       output_dc
!       output_meas
!
!         meas
!           meas_ground
!           meas_omi
!           meas_modis
!
!       output_data
!
!            state
!
!              dims
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
!###############################################################################


module LEKF_Driver

  use GO, only : gol, goPr, goErr

  use dims, only : runF, outF

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  LEKF_Arguments
  public  ::  LEKF_Driver_Init_State
  public  ::  LEKF_TimeStep_Run
  public  ::  LEKF_CalcMeas
  public  ::  TruncateStates

  !public  ::  biascorr_kf_surface_ozone

  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LEKF_Driver'


  ! --- var --------------------------------

  !character(len=64)  ::  biascorr_kf_surface_ozone = 'none'


contains


  ! ===================================================================
  ! ===
  ! === arguments
  ! ===
  ! ===================================================================


  subroutine LEKF_Arguments( rcfile, status )

    use GO, only : goArgCount, goGetArg

    ! --- in/out ----------------------------------

    character(len=*), intent(out)     ::  rcfile
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LEKF_Arguments'

    ! --- local -----------------------------------

    integer               ::  narg
    integer               ::  iarg
    character(len=1024)   ::  line

    ! --- begin -----------------------------------

    ! on root only, since some mpirun version do not parse
    ! all arguments to each executable:

    ! number of arguments:
    call goArgCount( narg, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    if ( narg == 0 ) then
      write (gol,'("no arguments found ...")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! defaults:
    rcfile = 'None'

    ! loop over arguments:
    iarg = 0
    do
      ! next:
      iarg = iarg + 1
      ! get argument:
      call goGetArg( iarg, line, status )
      IF_NOTOK_RETURN(status=1)
      ! not filled yet ?
      if ( trim(rcfile) == 'None' ) then
        rcfile = trim(line)
      else
        write (gol,'("unsupported argument : ",a)') trim(line); call goErr
        TRACEBACK; status=1; return
      end if
      ! last one is processed now ?
      if ( iarg == narg ) exit
    end do

    ! ok
    status = 0

  end subroutine LEKF_Arguments


  ! ===================================================================
  ! ===
  ! === driver routines
  ! ===
  ! ===================================================================


  subroutine LEKF_Driver_Init_State( x, t, status )

    use GO          , only : TDate
    use LE_Driver   , only : LE_State_Init
    use LE_DryDepos , only : LE_DryDepos_Setup_vd
    use LE_DryDepos , only : mix2ground
    use LE_BiasCorr , only : LE_BiasCorr_Fill
    use LEKF_state  , only : TState
    use LEKF_state  , only : bud0

    ! --- in/out -------------------------

    type(TState), intent(inout)     ::  x
    type(TDate), intent(in)         ::  t
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_Driver_Init_State'

    ! --- local --------------------------

    ! --- begin ---------------------------

    ! setup concentrations:
    call LE_State_Init( x%c, t, status )
    IF_NOTOK_RETURN(status=1)

    ! re-compute deposition velocities if these depend on concentrations;
    ! also the adhoc factor o3fac is applied if necessary ;
    ! use overal budget as input, used for NH3 compensation point:
    call LE_DryDepos_Setup_vd( x%c, bud0%drydepos, .true., t, status )
    IF_NOTOK_RETURN(status=1)

    ! surface concentrations:
    call mix2ground( x%c, x%cg, status )
    IF_NOTOK_RETURN(status=1)

    ! fill bias corrected tracers:
    call LE_BiasCorr_Fill( x%c, x%cg, status )
    IF_NOTOK_RETURN(status=1)
    ! ok
    status = 0

  end subroutine LEKF_Driver_Init_State


  ! ***


  subroutine LEKF_TimeStep_Run( t, dt, x, is_model, status, xb )

    use GO         , only : TDate, TIncrDate
    use Num        , only : gasdev
    use dims       , only : runF
    use dims       , only : nx,ny,nz,nspec
    use dims       , only : emis_a
    use indices
    use LE_Bound   , only : bc_west, bc_east, bc_south, bc_north
    use LE_Bound   , only : caloft
    use LE_DryDepos, only : vd_o3fac
    use LE_DryDepos, only : LE_DryDepos_Setup_vd
    use LE_DryDepos, only : mix2ground
    use LE_Driver  , only : LE_TimeStep_Run

    use LEKF_noise , only : nnoise, noise_name
    use LEKF_noise , only : dfc => disturb_factor
    use LEKF_State , only : TState
    use LEKF_State , only : bud0

    ! --- in/out ----------------------------

    type(TDate), intent(in)       ::  t
    type(TIncrDate), intent(in)   ::  dt
    type(TState), intent(inout)   ::  x
    logical, intent(in)           ::  is_model
    integer, intent(out)          ::  status

    type(TState), intent(in), optional   ::  xb

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_TimeStep_Run'

    ! --- local -----------------------------

    ! backup of model parameters to be changed:
    real, allocatable    ::  emis_save    (:,:,:,:)    ! dims (nx,ny,nz,nspec)
    real, allocatable    ::  bc_west_save (:,:,:)      ! dims (ny,nz,nspec)
    real, allocatable    ::  bc_east_save (:,:,:)      ! dims (ny,nz,nspec)
    real, allocatable    ::  bc_north_save(:,:,:)      ! dims (nx,nz,nspec)
    real, allocatable    ::  bc_south_save(:,:,:)      ! dims (nx,nz,nspec)
    real, allocatable    ::  caloft_save  (:,:,:,:)    ! dims (nx,ny,1,nspec)

    integer                 ::  nvoc
    integer                 ::  inmvoc

    integer                 ::  iz
    integer                 ::  ispec, icomp, k
    integer                 ::  inoise
    integer                 ::  ilu

    ! --- begin -------------------------------

    !! info ...
    !write (gol,'("KF:     <timestep>")'); call goPr

    !! only needed for additive noise ...
    !if ( (.not. is_model) .and. (.not. present(xb)) ) then
    !  write (gol,'("non-model step without argument xb passed")'); call goErr
    !  TRACEBACK; status=1; return
    !end if

    ! create temp storage
    allocate( emis_save    (nx,ny,nz,nspec) )
    allocate( bc_west_save (ny,nz,nspec) )
    allocate( bc_east_save (ny,nz,nspec) )
    allocate( bc_north_save(nx,nz,nspec) )
    allocate( bc_south_save(nx,nz,nspec) )
    allocate( caloft_save  (nx,ny,1,nspec) )

    ! ensure positive concentrations at start ...
    x%c  = max( 0.0, x%c  )

    !! DO NOT USE THIS! ADDITIVE NOISE MIGHT BE NEGATIVE ....
    !x%dc = max( 0.0, x%dc )

    ! backup model parameters:
    emis_save            = emis_a
    bc_west_save         = bc_west
    bc_east_save         = bc_east
    bc_north_save        = bc_north
    bc_south_save        = bc_south
    caloft_save(:,:,1,:) = caloft(:,:,nz+1,:)

    ! loop over noise types:
    do inoise = 1, nnoise

      ! try types:
      select case ( trim(noise_name(inoise)) )

        ! for testing ...
        case ( 'None' )
        
          ! no change at all, ensemble spread should be zero
          ! if this is the only noise parameter !
        ! surface emissions and adjacent boundary:
        case ( 'NOx_emis' )

          ! loop over levels:
          do iz = 1, nz
            ! scale emissions:
            emis_a(:,:,iz,i_no  ) = emis_a(:,:,iz,i_no  ) * x%dc(:,:,inoise,1) * dfc(inoise)
            emis_a(:,:,iz,i_no2 ) = emis_a(:,:,iz,i_no2 ) * x%dc(:,:,inoise,1) * dfc(inoise)
            ! scale lowest boundary layer:
            if ( iz == 1 ) then
               bc_west (:,iz,i_no  )  = bc_west (:,iz,i_no  ) * x%dc( 1, :,inoise,1) * dfc(inoise)
               bc_west (:,iz,i_no2 )  = bc_west (:,iz,i_no2 ) * x%dc( 1, :,inoise,1) * dfc(inoise)
               bc_east (:,iz,i_no  )  = bc_east (:,iz,i_no  ) * x%dc(nx, :,inoise,1) * dfc(inoise)
               bc_east (:,iz,i_no2 )  = bc_east (:,iz,i_no2 ) * x%dc(nx, :,inoise,1) * dfc(inoise)
               bc_south(:,iz,i_no  )  = bc_south(:,iz,i_no  ) * x%dc( :, 1,inoise,1) * dfc(inoise)
               bc_south(:,iz,i_no2 )  = bc_south(:,iz,i_no2 ) * x%dc( :, 1,inoise,1) * dfc(inoise)
               bc_north(:,iz,i_no  )  = bc_north(:,iz,i_no  ) * x%dc( :,ny,inoise,1) * dfc(inoise)
               bc_north(:,iz,i_no2 )  = bc_north(:,iz,i_no2 ) * x%dc( :,ny,inoise,1) * dfc(inoise)
            end if
          end do  ! levels

        ! surface emissions and adjacent boundary:
        case ( 'VOC_emis' )
        
          ! loop over nmvocs:
          do inmvoc = 1, n_nmvoc
             ! le spec index:
             ispec = ispecs_nmvoc(inmvoc)
               ! loop over levels:
               do iz = 1, nz
                 ! scale  emissions:
                 emis_a(:,:,iz,ispec) = emis_a(:,:,iz,ispec) * x%dc(:,:,inoise,1) * dfc(inoise)
                 ! scale lowest boundary layer:
                 if ( iz == 1 ) then
                    bc_west (:,iz,ispec)  = bc_west (:,iz,ispec) * x%dc( 1, :,inoise,1) * dfc(inoise)
                    bc_east (:,iz,ispec)  = bc_east (:,iz,ispec) * x%dc(nx, :,inoise,1) * dfc(inoise)
                    bc_south(:,iz,ispec)  = bc_south(:,iz,ispec) * x%dc( :, 1,inoise,1) * dfc(inoise)
                    bc_north(:,iz,ispec)  = bc_north(:,iz,ispec) * x%dc( :,ny,inoise,1) * dfc(inoise)
                 end if
              end do  ! levels
          end do  ! nmvocs
          
                ! surface emissions and adjacent boundary:
        ! case ( 'HCHO_emis' )
        
          ! ! loop over levels:
          ! do iz = 1, nz
            ! ! scale emissions:
            ! emis_a(:,:,iz,i_iso  ) = emis_a(:,:,iz,i_iso  ) * x%dc(:,:,inoise,1) * dfc(inoise)
            ! emis_a(:,:,iz,i_form ) = emis_a(:,:,iz,i_form ) * x%dc(:,:,inoise,1) * dfc(inoise)
            ! ! scale lowest boundary layer:
            ! if ( iz == 1 ) then
               ! bc_west (:,iz,i_iso  )  = bc_west (:,iz,i_iso  ) * x%dc( 1, :,inoise,1) * dfc(inoise)
               ! bc_west (:,iz,i_form )  = bc_west (:,iz,i_form ) * x%dc( 1, :,inoise,1) * dfc(inoise)
               ! bc_east (:,iz,i_iso  )  = bc_east (:,iz,i_iso  ) * x%dc(nx, :,inoise,1) * dfc(inoise)
               ! bc_east (:,iz,i_form )  = bc_east (:,iz,i_form ) * x%dc(nx, :,inoise,1) * dfc(inoise)
               ! bc_south(:,iz,i_iso  )  = bc_south(:,iz,i_iso  ) * x%dc( :, 1,inoise,1) * dfc(inoise)
               ! bc_south(:,iz,i_form )  = bc_south(:,iz,i_form ) * x%dc( :, 1,inoise,1) * dfc(inoise)
               ! bc_north(:,iz,i_iso  )  = bc_north(:,iz,i_iso  ) * x%dc( :,ny,inoise,1) * dfc(inoise)
               ! bc_north(:,iz,i_form )  = bc_north(:,iz,i_form ) * x%dc( :,ny,inoise,1) * dfc(inoise)
            ! end if
          ! end do  ! levels

        ! upper boundary condition
        case ( 'o3_top' )

          ! scale top boundary condition values:
          caloft(:,:,nz+1,i_o3) = caloft(:,:,nz+1,i_o3) * x%dc(:,:,inoise,1) * dfc(inoise)

        !! extra error:
        !case ( 'o3_xerr', 'upm_f_xerr', 'upm_c_xerr' )
        !  ! not yet, added after emis etc

        ! surface emissions and adjacent boundary:
        case ( 'SOx_emis' )

          ! loop over levels:
          do iz = 1, nz
            ! scale emissions:
            emis_a(:,:,iz,i_so2   ) = emis_a(:,:,iz,i_so2   ) * x%dc(:,:,inoise,1) * dfc(inoise)
            emis_a(:,:,iz,i_so4a_f) = emis_a(:,:,iz,i_so4a_f) * x%dc(:,:,inoise,1) * dfc(inoise)
            ! scale lowest boundary layer:
            if ( iz == 1 ) then
              bc_west (:,iz,i_so2   ) = bc_west (:,iz,i_so2   ) * x%dc( 1, :,inoise,1) * dfc(inoise)
              bc_west (:,iz,i_so4a_f) = bc_west (:,iz,i_so4a_f) * x%dc( 1, :,inoise,1) * dfc(inoise)
              bc_east (:,iz,i_so2   ) = bc_east (:,iz,i_so2   ) * x%dc(nx, :,inoise,1) * dfc(inoise)
              bc_east (:,iz,i_so4a_f) = bc_east (:,iz,i_so4a_f) * x%dc(nx, :,inoise,1) * dfc(inoise)
              bc_south(:,iz,i_so2   ) = bc_south(:,iz,i_so2   ) * x%dc( :, 1,inoise,1) * dfc(inoise)
              bc_south(:,iz,i_so4a_f) = bc_south(:,iz,i_so4a_f) * x%dc( :, 1,inoise,1) * dfc(inoise)
              bc_north(:,iz,i_so2   ) = bc_north(:,iz,i_so2   ) * x%dc( :,ny,inoise,1) * dfc(inoise)
              bc_north(:,iz,i_so4a_f) = bc_north(:,iz,i_so4a_f) * x%dc( :,ny,inoise,1) * dfc(inoise)
            end if
          end do  ! levels

        ! surface emissions and adjacent boundary:
        case ( 'NH3_emis' )

          ! loop over levels:
          do iz = 1, nz
            ! scale emissions:
            emis_a(:,:,iz,i_nh3 ) = emis_a(:,:,iz,i_nh3 ) * x%dc(:,:,inoise,1) * dfc(inoise)
            ! scale lowest boundary layer:
            if ( iz == 1 ) then
              bc_west (:,iz,i_nh3 )  = bc_west (:,iz,i_nh3 ) * x%dc( 1, :,inoise,1) * dfc(inoise)
              bc_east (:,iz,i_nh3 )  = bc_east (:,iz,i_nh3 ) * x%dc(nx, :,inoise,1) * dfc(inoise)
              bc_south(:,iz,i_nh3 )  = bc_south(:,iz,i_nh3 ) * x%dc( :, 1,inoise,1) * dfc(inoise)
              bc_north(:,iz,i_nh3 )  = bc_north(:,iz,i_nh3 ) * x%dc( :,ny,inoise,1) * dfc(inoise)
            end if
          end do  ! levels

        ! surface emissions and adjacent boundary:
        case ( 'PM25_BC_emis' )

          ! loop over levels:
          do iz = 1, nz
            ! scale emissions:
            emis_a(:,:,iz,i_ppm_f) = emis_a(:,:,iz,i_ppm_f) * x%dc(:,:,inoise,1) * dfc(inoise)
            emis_a(:,:,iz,i_ec_f ) = emis_a(:,:,iz,i_ec_f ) * x%dc(:,:,inoise,1) * dfc(inoise)
            ! scale lowest boundary layer:
            if ( iz == 1 ) then
              bc_west (:,iz,i_ppm_f) = bc_west (:,iz,i_ppm_f) * x%dc( 1, :,inoise,1) * dfc(inoise)
              bc_west (:,iz,i_ec_f ) = bc_west (:,iz,i_ec_f ) * x%dc( 1, :,inoise,1) * dfc(inoise)
              bc_east (:,iz,i_ppm_f) = bc_east (:,iz,i_ppm_f) * x%dc(nx, :,inoise,1) * dfc(inoise)
              bc_east (:,iz,i_ec_f ) = bc_east (:,iz,i_ec_f ) * x%dc(nx, :,inoise,1) * dfc(inoise)
              bc_south(:,iz,i_ppm_f) = bc_south(:,iz,i_ppm_f) * x%dc( :, 1,inoise,1) * dfc(inoise)
              bc_south(:,iz,i_ec_f ) = bc_south(:,iz,i_ec_f ) * x%dc( :, 1,inoise,1) * dfc(inoise)
              bc_north(:,iz,i_ppm_f) = bc_north(:,iz,i_ppm_f) * x%dc( :,ny,inoise,1) * dfc(inoise)
              bc_north(:,iz,i_ec_f ) = bc_north(:,iz,i_ec_f ) * x%dc( :,ny,inoise,1) * dfc(inoise)
            end if
          end do  ! levels

        ! deposition velocities
        case ( 'o3_vd' )
             vd_o3fac(:,:) = x%dc(:,:,inoise,1) * dfc(inoise)
 
        ! dust boundary conditions:
        case ( 'dust_bc' )

          ! loop over levels:
          do iz = 1, nz
            ! scale boundary layers:
            bc_west (:,iz,i_dust_f)  = bc_west (:,iz,i_dust_f) * x%dc( 1, :,inoise,1) * dfc(inoise)
            bc_west (:,iz,i_dust_c)  = bc_west (:,iz,i_dust_c) * x%dc( 1, :,inoise,1) * dfc(inoise)
            bc_east (:,iz,i_dust_f)  = bc_east (:,iz,i_dust_f) * x%dc(nx, :,inoise,1) * dfc(inoise)
            bc_east (:,iz,i_dust_c)  = bc_east (:,iz,i_dust_c) * x%dc(nx, :,inoise,1) * dfc(inoise)
            bc_south(:,iz,i_dust_f)  = bc_south(:,iz,i_dust_f) * x%dc( :, 1,inoise,1) * dfc(inoise)
            bc_south(:,iz,i_dust_c)  = bc_south(:,iz,i_dust_c) * x%dc( :, 1,inoise,1) * dfc(inoise)
            bc_north(:,iz,i_dust_f)  = bc_north(:,iz,i_dust_f) * x%dc( :,ny,inoise,1) * dfc(inoise)
            bc_north(:,iz,i_dust_c)  = bc_north(:,iz,i_dust_c) * x%dc( :,ny,inoise,1) * dfc(inoise)
          end do  ! levels

        ! other ...
        case default

          write (gol,'("unsupported noise name `",a,"`")') trim(noise_name(inoise)); call goErr
          TRACEBACK; status=1; return

      end select

    end do  ! noise

    ! TESTING: fill upm boundary conditions:
    ! REJECTED: this gives clouds over atlantic ;
    !   concentrations in there are not too large, but clearly visible;
    !   might be due to adding all kind of aerosols (incl soluble) into uppm
    !do k = 1, n_unspecified
    !  ! global index:
    !  ispec = ispecs_unspecified(k)
    !  ! enabled ?
    !  if ( ispec > 0 ) then
    !    ! clear:
    !    bc_west (:,:,ispec) = 0.0
    !    bc_east (:,:,ispec) = 0.0
    !    bc_south(:,:,ispec) = 0.0
    !    bc_north(:,:,ispec) = 0.0
    !    ! loop over all:
    !    do icomp = 1, nspec
    !      ! skip natural emissions:
    !      if ( tracer_is_dust   (icomp) ) cycle
    !      if ( tracer_is_seasalt(icomp) ) cycle
    !      ! no unspecified tracers, otherwise it includes itself:
    !      if ( tracer_is_unspecified(icomp) ) cycle
    !      ! select aerosols from same size mode:
    !      if ( ispec == ispec_upm_f ) then
    !        if ( specmode(icomp) /= AEROSOL_FINE_MODES ) cycle
    !      else if ( ispec == ispec_upm_c ) then
    !        if ( specmode(icomp) /= AEROSOL_COARSE_MODE ) cycle
    !      else
    !        write (gol,'("unsupported ispec ",i0," (",a,")")') ispec, trim(specname(ispec)); call goErr
    !        TRACEBACK; status=1; return
    !      end if
    !      ! check ...
    !      if ( trim(specunit(icomp)) /= trim(specunit(ispec)) ) then
    !        write (gol,'("emission source units `",a,"` for `",a,"`")') trim(specunit(icomp)), trim(specunit(icomp)); call goPr
    !        write (gol,'("  does not match with target units `",a,"` for `",a,"`")') trim(specunit(ispec)), trim(specunit(ispec)); call goPr
    !        TRACEBACK; status=1; return
    !      end if
    !      ! add contribution:
    !      bc_west (:,:,ispec) = bc_west (:,:,ispec) + bc_west (:,:,icomp)
    !      bc_east (:,:,ispec) = bc_east (:,:,ispec) + bc_east (:,:,icomp)
    !      bc_south(:,:,ispec) = bc_south(:,:,ispec) + bc_south(:,:,icomp)
    !      bc_north(:,:,ispec) = bc_north(:,:,ispec) + bc_north(:,:,icomp)
    !    end do  ! spec
    !  end if  ! enabled
    !end do  ! upm tracers

    ! * testing: additive noise (bias correction?)
    ! only for ensemble members, might need xb ...
    if ( .not. is_model ) then
      ! loop over noise types:
      do inoise = 1, nnoise

        ! try types:
        select case ( trim(noise_name(inoise)) )

          ! already done ...
          case ( 'None' )
          case ( 'NOx_emis', 'VOC_emis', 'SOx_emis', 'NH3_emis', 'PM25_BC_emis' ,'HCHO_emis' )
          case ( 'o3_top', 'o3_vd' )
          case ( 'dust_bc' )

          !! extra error:
          !case ( 'o3_xerr', 'upm_f_xerr', 'upm_c_xerr' )
          !
          !  ! switch:
          !  select case ( trim(noise_name(inoise)) )
          !    case ( 'o3_xerr'    ) ; ispec = i_o3
          !    case ( 'upm_f_xerr' ) ; ispec = i_upm_f
          !    case ( 'upm_c_xerr' ) ; ispec = i_upm_c
          !    case default
          !      write (gol,'("unsupported noise name `",a,"`")') trim(noise_name(inoise)); call goErr
          !      TRACEBACK; status=1; return
          !  end select
          !
          !  ! enabled ?
          !  if ( ispec > 0 ) then
          !
          !    ! check ...
          !    if ( (.not. is_model) .and. (.not. present(xb)) ) then
          !      write (gol,'("additive noise requires that xb is passed as optional argument ...")'); call goErr
          !      TRACEBACK; status=1; return
          !    end if

          !    ! new: values around 0.0 with sigma as amplitude,
          !    ! negative values are removed below;
          !    ! surface and boundary layer only:
          !    do iz = 1, 2
          !
          !      ! ~ additive, zero mean
          !      !! add to concentrations:
          !      !x%c(:,:,iz,ispec) = x%c(:,:,iz,ispec) + x%dc(:,:,inoise,1)
          !      !! add to boundary conditions:
          !      !bc_west (:,iz,ispec) = bc_west (:,iz,ispec) + x%dc(1 ,: ,inoise,1)
          !      !bc_east (:,iz,ispec) = bc_east (:,iz,ispec) + x%dc(nx,: ,inoise,1)
          !      !bc_south(:,iz,ispec) = bc_south(:,iz,ispec) + x%dc(: ,1 ,inoise,1)
          !      !bc_north(:,iz,ispec) = bc_north(:,iz,ispec) + x%dc(: ,ny,inoise,1)
          !
          !      ! ~ fraction of background
          !      ! add to concentrations:
          !      x%c(:,:,iz,ispec) = x%c(:,:,iz,ispec) + xb%c(:,:,iz,ispec) * x%dc(:,:,inoise,1)
          !      !! add to boundary conditions:
          !      !! REJECTED, see above
          !      !bc_west (:,iz,ispec) = bc_west (:,iz,ispec) * ( 1.0 + x%dc(1 ,: ,inoise,1) )
          !      !bc_east (:,iz,ispec) = bc_east (:,iz,ispec) * ( 1.0 + x%dc(nx,: ,inoise,1) )
          !      !bc_south(:,iz,ispec) = bc_south(:,iz,ispec) * ( 1.0 + x%dc(: ,1 ,inoise,1) )
          !      !bc_north(:,iz,ispec) = bc_north(:,iz,ispec) * ( 1.0 + x%dc(: ,ny,inoise,1) )
          !
          !    end do
          !
          !  end if  ! spec enabled

          ! other ...
          case default
            write (gol,'("unsupported noise name `",a,"`")') trim(noise_name(inoise)); call goErr
            TRACEBACK; status=1; return
        end select

      end do  ! noise

    end if ! not the model run

#ifdef skip_timestep

    ! dummy ...
    x%c = x%c * max( 0.01, 1.0+gasdev()*0.10 )

#else

    ! model timestep;
    ! deposition velocities are recomputed before depos, 
    ! incl. application of o3fac if necessary ;
    ! update budgets only for xb:
    call LE_TimeStep_Run( t, dt, x%c, x%cg, x%aerh2o, &
                            bud0, is_model, &
                            status )
    IF_NOTOK_RETURN(status=1)

#endif

    ! positive concentrations only ...
    where ( x%c <0.0 ) x%c = 0.0

    ! reset the emissions and deposition
    emis_a             = emis_save
    bc_west            = bc_west_save
    bc_east            = bc_east_save
    bc_north           = bc_north_save
    bc_south           = bc_south_save
    caloft(:,:,nz+1,:) = caloft_save(:,:,1,:)
    ! reset depositon field and factor:
    vd_o3fac = 1.0

    ! clear:
    deallocate( emis_save )
    deallocate( bc_west_save )
    deallocate( bc_east_save )
    deallocate( bc_north_save )
    deallocate( bc_south_save )
    deallocate( caloft_save )

    ! ok
    status = 0

  end subroutine LEKF_TimeStep_Run


  ! ***


  subroutine LEKF_CalcMeas( x, t, status )

    use GO         , only : TDate
    use dims       , only : nx,ny,nz,nspec
    use LEKF_State , only : TState

#ifdef with_kf_meas_sat
    use LEKF_Data        , only : leo
#endif
#ifdef with_kf_meas_modis
    use LEKF_Meas_MODIS, only : CalcMeas_MODIS
#endif
#ifdef with_kf_meas_maori
    use LE_MAORI  , only : LE_MAORI_State_Put
    use LEKF_Data , only : mad
#endif

    ! --- in/out ----------------------------

    type(TState), intent(inout)   ::  x
    type(TDate), intent(in)       ::  t
    integer, intent(out)          ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LEKF_CalcMeas'

    ! --- local -----------------------------
    
    integer       ::  k

    ! --- begin -------------------------------

#ifdef with_kf_meas_aod
    ! re-compute the AOD fields:
    call calc_AOD( x%c, x%AOD )
#endif

#ifdef with_kf_meas_sat
    ! any sat output ?
	write(gol,*) 'writing x%nsat',x%nsat; call GoPr
    if ( x%nsat > 0 ) then
      ! loop over files to be written:
	  
      do k = 1, x%nsat
        ! simulate retrievals (see "le_output.F90" for original):
        call x%sat_state(k)%Setup( leo%cso_data(k), x%c, status )
        IF_NOTOK_RETURN(status=1)
      end do  ! output files
    end if  ! any output?
#endif

#ifdef with_kf_meas_modis
    ! simulate AOD field:
    call CalcMeas_MODIS( x%c, x%modis, x%nmodis, status )
    IF_NOTOK_RETURN(status=1)
#endif

#ifdef with_kf_meas_maori
    ! put model simulations into MAORI arrays:
    call LE_MAORI_State_Put( x%mas, mad, t, x%c, x%cg, x%aerh2o, status )
    IF_NOTOK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine LEKF_CalcMeas


  ! ***


  subroutine TruncateStates( status )

    use LE_DryDepos , only : mix2ground
    use LEKF_dims   , only : nx, ny, nz, nspec
    use LEKF_state  , only : nmodes
    use LEKF_state  , only : Ens
    use LEKF_state  , only : Ens_Mean_and_Sigma
    use LEKF_Noise  , only : ClipDC
#ifdef with_kf_meas_modis
    use LEKF_Meas_MODIS, only : CalcMeas_MODIS
#endif

    ! --- in/out -----------------------------

    integer, intent(out)          ::  status

    ! --- const ------------------------------

    character(len=*), parameter   ::  rname = mname//'/TruncateStates'

    ! --- local ------------------------------

    integer     ::  j
    integer     ::  inoise

    ! --- begin ------------------------------

    ! truncate:
    do j = 1, nmodes
       Ens(j)%c       = max( 0.0, Ens(j)%c       )
       Ens(j)%aerh2o  = max( 0.0, Ens(j)%aerh2o  )
      !Ens(j)%aod     = max( 0.0, Ens(j)%aod     )
      ! clip dc values to [min,max] range:
      call ClipDC( Ens(j)%dc, status )
      IF_NOTOK_RETURN(status=1)
#ifdef with_kf_meas_modis
       Ens(j)%modis   = max( 0.0, Ens(j)%modis   )
#endif
    end do  ! modes

    ! update measurement interpolations:
    do j = 1, nmodes

#ifdef with_kf_meas_modis
       ! simulate AOD field:
       call CalcMeas_MODIS( Ens(j)%c, Ens(j)%modis, Ens(j)%nmodis, status )
       IF_NOTOK_RETURN(status=1)
#endif

    end do

    ! new sample mean/stdv:
    call Ens_Mean_and_Sigma( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine TruncateStates


end module LEKF_Driver

