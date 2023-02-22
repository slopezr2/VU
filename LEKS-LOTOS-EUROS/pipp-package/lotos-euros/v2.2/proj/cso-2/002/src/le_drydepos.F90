!#################################################################
!
! Deposition
!
! Version including aerosol sedimentation following Zhang.
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

module LE_DryDepos

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private
  
  public  ::  spectype_depos
  public  ::  SPECTYPE_DEPOS_GAS_DEPAC, SPECTYPE_DEPOS_AERO_ZHANG
  
  public  ::  LE_DryDepos_Init, LE_DryDepos_Done
  public  ::  LE_DryDepos_Apply
  public  ::  LE_DryDepos_Setup_vd
  public  ::  mix2ground
  public  ::  mix2ground_lu
  public  ::  mix2canopytop
  
  public  ::  Rb                ! quasi-laminar (leaf) boundary layer resistance
  public  ::  Rc_eff            ! effective canopy resistance (s/m) assuming that Ccomp_tot=0 (which might
                                ! not be true for ammonia) and might become negative in case of emission 
  public  ::  Rc_total          ! Rc als used for vd
  public  ::  Rs                ! sedimentation resistance (for aerosols)                              
!  public  ::  Boundary_Resistance_lu_tr
!  public  ::  Canopy_Surface_Resistance_lu_tr
!  public  ::  Stomatal_Resistance_lu_tr

  public  ::  vd_lu
  
  public  ::  vs    ! sedimentation velocity

  ! ozone flux specials:
  public  ::  vd_zcanopytop_lu
  public  ::  vd_htop_lu
  public  ::  Cc_lu
  public  ::  Fdry_total_lu
  public  ::  Fdry_stom_lu
  public  ::  Fdry_w_lu
  public  ::  Fdry_soil_lu
  public  ::  Cstom
  public  ::  Cw
  public  ::  Csoil
  public  ::  gstom
  public  ::  gw
  public  ::  gsoil_eff  
  public  ::  frac_sto_o3_lu
  public  ::  rbc_o3_lu
  public  ::  czcanopytop_o3_lu
  
  public  ::  vd_o3fac ! assim parameter field


  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_DryDepos'

  ! Dry deposition of species is treated in different ways:
  ! no depos at all
  integer, parameter   :: SPECTYPE_DEPOS_NOT            = 0
  ! gasses following DEPAC 311
  integer, parameter   :: SPECTYPE_DEPOS_GAS_DEPAC      = 1
  !! aerosols following Wesely (1985), Erisman (1994), and Ruijgrok (1994)
  integer, parameter   :: SPECTYPE_DEPOS_AERO_WER       = 2
  ! aerosols following Zhang
  integer, parameter   :: SPECTYPE_DEPOS_AERO_ZHANG     = 3
  ! aerosols for M7
  integer, parameter   :: SPECTYPE_DEPOS_AERO_M7_MASS   = 4
  integer, parameter   :: SPECTYPE_DEPOS_AERO_M7_NUMBER = 5
  ! pollen deposition following Sofiev et al, 2006 (towards numerical...)
  integer, parameter   :: SPECTYPE_DEPOS_AERO_POLLEN    = 6
  ! labels:
  character(len=*), parameter  ::  SPECTYPE_DEPOS_LABEL(0:6) = &
         (/ 'NOT           ', &
            'GAS_DEPAC     ', 'AERO_WER      ', 'AERO_ZHANG    ', &
            'AERO_M7_MASS  ', 'AERO_M7_NUMBER', 'POLLEN        ' /)


  ! --- var --------------------------------------

  ! Resistance parameters
  real   , allocatable  ::  Rb(:,:,:,:)        ! (nx,ny,nspec,nlu)
  real   , allocatable  ::  Rc_eff(:,:,:,:)    ! (nx,ny,nspec,nlu)
  real   , allocatable  ::  Rc_total(:,:,:,:)  ! (nx,ny,nspec,nlu)
  real   , allocatable  ::  Rs(:,:,:,:)

  ! deposition velocity (m/s) per land use class
  real   , allocatable  ::  vd_lu(:,:,:,:)
  
  ! sedimentation velocity:
  real   , allocatable  ::  vs(:,:,:)

  ! adhoc factor for assim:
  real, allocatable     ::  vd_o3fac(:,:)  ! (nx,ny)

  ! extra output:
  logical               ::  extra_output
  real, allocatable     ::  Cc_lu(:,:,:,:)
  real, allocatable     ::  Fdry_total_lu(:,:,:,:)
  real, allocatable     ::  Fdry_stom_lu(:,:,:,:)
  real, allocatable     ::  Fdry_w_lu(:,:,:,:)
  real, allocatable     ::  Fdry_soil_lu(:,:,:,:)
  real, allocatable     ::  Cstom(:,:,:,:)
  real, allocatable     ::  Cw(:,:,:,:)
  real, allocatable     ::  Csoil(:,:,:,:)
  real, allocatable     ::  gstom(:,:,:,:)
  real, allocatable     ::  gw(:,:,:,:)
  real, allocatable     ::  gsoil_eff(:,:,:,:)
  ! only for ozone:
  real, allocatable     ::  vd_zcanopytop_lu(:,:,:)
  real, allocatable     ::  vd_htop_lu(:,:,:)
  real, allocatable     ::  czcanopytop_o3_lu(:,:,:)
  real, allocatable     ::  fac_surface_area_2_PLA(:,:,:)
  real, allocatable     ::  frac_sto_o3_lu(:,:,:)
  real, allocatable     ::  rbc_o3_lu(:,:,:)

  ! type of deposition (not, gas_depac, gas_similar, aero_f, aero_c)
  integer, allocatable :: spectype_depos(:)   ! (nspec)
  ! depends on concentrations ?
  logical, allocatable ::  vd_depends_on_c(:)  ! (nspec)

  ! Compensation point (= 0, but for NH3):
  real, allocatable, dimension(:,:,:,:) :: ccomp_tot_lu ! ppb (nx,ny,nspec,nlu_depac)
 
  integer, allocatable  ::  ns_regime(:,:)

  ! solubility/reactivity parameters   
  real, allocatable :: henry0(:)
  real, allocatable :: Tfact(:)  
  real, allocatable :: reactivity(:)
  
  ! timers:
  integer   ::  itim_vdepos
  integer   ::  itim_vdepos_vs, itim_vdepos_rb, itim_vdepos_rc
#ifdef with_m7
  integer   ::  itim_vdepos_m7
#endif  



contains


  ! ========================================================================


  subroutine LE_DryDepos_Init( rcF, status )

    use GO, only : TrcFile, ReadRc
    use GO, only : GO_Timer_Def
    use GO     , only : GoSplitString, GoGetFu, GoMatchValue
    use dims         , only : nx, ny, nspec
    use Indices      , only : specmode, NO_AEROSOL_MODE, AEROSOL_FINE_MODE, AEROSOL_COARSE_MODE
    use Indices      , only : i_so2, i_hno3, i_nh3, i_no, i_no2, i_o3, i_n2o5, i_no3, i_h2o2, i_co, &
                              i_vbs_pog1, i_vbs_pog2, i_vbs_pog3, i_vbs_pog4, i_vbs_pog5, &
                              i_vbs_pog6, i_vbs_pog7, i_vbs_pog8, i_vbs_pog9, &
                              i_vbs_sisog1, i_vbs_sisog2, i_vbs_sisog3, i_vbs_sisog4, i_vbs_sisog5, &
                              i_vbs_sisog6, i_vbs_sisog7, i_vbs_sisog8, &
                              i_vbs_asog1, i_vbs_asog2, i_vbs_asog3, &
                              i_vbs_asog4, i_vbs_asog5, i_vbs_asog6, &
                              i_vbs_bsog1, i_vbs_bsog2, i_vbs_bsog3, &
                              i_vbs_bsog4, i_vbs_bsog5, i_vbs_bsog6
#ifdef with_m7    
    use Indices      , only : tracer_is_m7, tracer_is_numberconc
#endif
#ifdef with_pollen
    use Indices      , only : tracer_is_pollen
#endif
    use Indices      , only : specname
    use LE_Landuse_Data, only : nlu
    use LE_DryDepos_Gas_GammaWater, only : LE_DryDepos_Gas_GammaWater_Init
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_DryDepos_Init'

    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 10

    ! --- local ---------------------------------
    
    integer                    ::  k
    integer                    ::  ispec
    character(len=1024)        ::  fname_template
    logical                    ::  apply

    character(len=512)      ::  wetscav_datafile
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    integer                 ::  fu
    logical                 ::  exist
    
    real                    ::  value
    integer                 ::  iline, iheader!, ispec
    character(len=1024)     ::  line

    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    
    integer                 ::  nfield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    integer                 ::  ifield
    
    integer                 ::  ifield_tracer
    integer                 ::  ifield_lbnd
    integer                 ::  ifield_hlaw
    integer                 ::  ifield_Tfac
    integer                 ::  ifield_diffrat
    integer                 ::  ifield_rscale
    integer                 ::  ifield_reactivity

    ! --- begin ----------------------------------
    
    ! define timers:
    call GO_Timer_Def( itim_vdepos  , 'dry deposition velocities', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_vdepos_vs, 'vdepos vsed', status )
    IF_NOTOK_RETURN(status=1)
#ifdef with_m7    
    call GO_Timer_Def( itim_vdepos_m7, 'vdepos M7 table', status )
    IF_NOTOK_RETURN(status=1)
#endif
    call GO_Timer_Def( itim_vdepos_rb, 'vdepos rb', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_vdepos_rc, 'vdepos rc', status )
    IF_NOTOK_RETURN(status=1)

    ! * ozone flux specials

    ! read flag:
    call ReadRc( rcF, 'deposition.extra_output', extra_output, status )
    IF_NOTOK_RETURN(status=1)

    ! * deposition method
    
    ! storage for deposition method per tracer:
    allocate( spectype_depos(nspec) )
    allocate( vd_depends_on_c(nspec) )

    ! Default no deposition:
    spectype_depos = SPECTYPE_DEPOS_NOT
    ! default only depends on meteo, not on concentrations:
    vd_depends_on_c = .false.

    ! set depos type:
    if (i_hno3>0) spectype_depos(i_hno3) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_no  >0) spectype_depos(i_no  ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_no2 >0) spectype_depos(i_no2 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_o3  >0) spectype_depos(i_o3  ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_n2o5>0) spectype_depos(i_n2o5) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_no3 >0) spectype_depos(i_no3 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_h2o2>0) spectype_depos(i_h2o2) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_co  >0) spectype_depos(i_co  ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog1 >0) spectype_depos(i_vbs_pog1 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog2 >0) spectype_depos(i_vbs_pog2 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog3 >0) spectype_depos(i_vbs_pog3 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog4 >0) spectype_depos(i_vbs_pog4 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog5 >0) spectype_depos(i_vbs_pog5 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog6 >0) spectype_depos(i_vbs_pog6 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog7 >0) spectype_depos(i_vbs_pog7 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog8 >0) spectype_depos(i_vbs_pog8 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_pog9 >0) spectype_depos(i_vbs_pog9 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_sisog1 >0) spectype_depos(i_vbs_sisog1 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_sisog2 >0) spectype_depos(i_vbs_sisog2 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_sisog3 >0) spectype_depos(i_vbs_sisog3 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_sisog4 >0) spectype_depos(i_vbs_sisog4 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_sisog5 >0) spectype_depos(i_vbs_sisog5 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_sisog6 >0) spectype_depos(i_vbs_sisog6 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_sisog7 >0) spectype_depos(i_vbs_sisog7 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_sisog8 >0) spectype_depos(i_vbs_sisog8 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_asog1 >0) spectype_depos(i_vbs_asog1 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_asog2 >0) spectype_depos(i_vbs_asog2 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_asog3 >0) spectype_depos(i_vbs_asog3 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_asog4 >0) spectype_depos(i_vbs_asog4 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_asog5 >0) spectype_depos(i_vbs_asog5 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_asog6 >0) spectype_depos(i_vbs_asog6 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_bsog1 >0) spectype_depos(i_vbs_bsog1 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_bsog2 >0) spectype_depos(i_vbs_bsog2 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_bsog3 >0) spectype_depos(i_vbs_bsog3 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_bsog4 >0) spectype_depos(i_vbs_bsog4 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_bsog5 >0) spectype_depos(i_vbs_bsog5 ) = SPECTYPE_DEPOS_GAS_DEPAC
    if (i_vbs_bsog6 >0) spectype_depos(i_vbs_bsog6 ) = SPECTYPE_DEPOS_GAS_DEPAC

    ! depends on SO2/NH3 ratio:
    if (i_so2 >0) then
      spectype_depos (i_so2) = SPECTYPE_DEPOS_GAS_DEPAC
      vd_depends_on_c(i_so2) = .true.
    end if
    ! uses compensation point:
    if (i_nh3 >0) then
      spectype_depos (i_nh3) = SPECTYPE_DEPOS_GAS_DEPAC
      vd_depends_on_c(i_nh3) = .true.
    end if
    
    ! aerosol deposition:
    do ispec = 1, nspec
      
      if ( specmode(ispec) == NO_AEROSOL_MODE ) then
        ! skip gasses, already set above
        cycle
      else
        ! following Zhang scheme:
        spectype_depos(ispec) = SPECTYPE_DEPOS_AERO_ZHANG
        ! aerosol size is constant for fine/course definition:
        vd_depends_on_c(ispec) = .false.
        
#ifdef with_m7      
        ! M7 or not ?
        if ( tracer_is_m7(ispec) ) then
          ! different for mass and number:
          if ( tracer_is_numberconc(ispec) ) then
            spectype_depos(ispec) = SPECTYPE_DEPOS_AERO_M7_NUMBER
          else
            spectype_depos(ispec) = SPECTYPE_DEPOS_AERO_M7_MASS
          end if
          ! vd depends on particle diameter,
          ! for M7 this depends on mass and number ratio:
          vd_depends_on_c(ispec) = .true.
        end if ! m7 type
#endif

#ifdef with_pollen        
        ! Pollen or not?
        if ( tracer_is_pollen(ispec) ) then
          ! 
          spectype_depos(ispec)  = SPECTYPE_DEPOS_AERO_POLLEN
          ! aerosol size is constant by definition
          vd_depends_on_c(ispec) = .false.        
        end if  ! aerosol type
#endif
                
      end if  ! specmode      
    end do ! specs
        
    ! *
    
    ! Allocate and initialise arrays:
    allocate( ccomp_tot_lu(nx,ny,nspec,nlu) ) ; ccomp_tot_lu  = 0.0
     
    ! storage:
    allocate( Rb(nx,ny,nspec,nlu) )
    allocate( Rc_eff(nx,ny,nspec,nlu) )
    allocate( Rc_total(nx,ny,nspec,nlu) )
    allocate( Rs(nx,ny,nspec,nlu) )
    
    allocate( vd_lu(nx,ny,nspec,nlu) )

    allocate( vs(nx,ny,nspec) )

    ! storage:
    allocate( vd_o3fac(nx,ny) )
    ! init to unity:
    vd_o3fac = 1.0
    
    if ( extra_output ) then
      allocate( Cc_lu(nx,ny,nspec,nlu) )
      allocate( Fdry_total_lu(nx,ny,nspec,nlu) )
      allocate( Fdry_stom_lu(nx,ny,nspec,nlu) )
      allocate( Fdry_w_lu(nx,ny,nspec,nlu) )
      allocate( Fdry_soil_lu(nx,ny,nspec,nlu) )
      allocate( Cstom(nx,ny,nspec,nlu) )
      allocate( Cw(nx,ny,nspec,nlu) )
      allocate( Csoil(nx,ny,nspec,nlu) )
      allocate( gstom(nx,ny,nspec,nlu) )
      allocate( gw(nx,ny,nspec,nlu) )
      allocate( gsoil_eff(nx,ny,nspec,nlu) )
      ! only ozone:
      allocate( vd_zcanopytop_lu(nx,ny,nlu) )
      allocate( vd_htop_lu(nx,ny,nlu) )
      allocate( czcanopytop_o3_lu     (nx,ny,nlu) )
      allocate( fac_surface_area_2_PLA(nx,ny,nlu) )
      allocate( frac_sto_o3_lu        (nx,ny,nlu) )
      allocate( rbc_o3_lu             (nx,ny,nlu) )
    end if
    
    
    allocate( ns_regime(nx,ny) )

    allocate ( henry0(nspec) )
    allocate ( Tfact(nspec) )
    allocate ( reactivity(nspec) )
    
    ! init to zero:
    vd_lu = 0.0
    vs = 0.0
    
    if ( extra_output ) then
      Cc_lu         = 0.0
      Fdry_total_lu = 0.0
      Fdry_stom_lu  = 0.0
      Fdry_w_lu     = 0.0
      Fdry_soil_lu  = 0.0
      Cstom         = 0.0
      Cw            = 0.0
      Csoil         = 0.0
      gstom         = 0.0
      gw            = 0.0
      gsoil_eff     = 0.0
      ! by default no ozone deposition:
      vd_zcanopytop_lu       = 0.0
      vd_htop_lu             = 0.0
      frac_sto_o3_lu         = 0.0
      rbc_o3_lu              = 0.0
      czcanopytop_o3_lu      = 0.0
      fac_surface_area_2_PLA = 0.0
    end if

    ! for all gaseous species, independ of meteo:
    Rb = 0.0
    Rs = 0.0

    henry0 = 0.0
    Tfact  = 0.0
    reactivity  = 0.0

    ! * info

    !! Write non-deposition species:
    !write(*,*) '------------------------------------------------------'
    !do ispec = 1,nspec
    !   if (spectype_depos(ispec) .eq. spectype_depos_not) &
    !      write(*,*) 'no deposition for species ',specname(ispec)
    !enddo
    !write(*,*) '------------------------------------------------------'
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("dry deposition scheme per tracer:")'); call goPr
    write (gol,'("------ --------------- ----- --------------- ----")'); call goPr
    write (gol,'(" ispec        specname dtype           label depc")'); call goPr
    write (gol,'("------ --------------- ----- --------------- ----")'); call goPr
    do ispec = 1, nspec
      !if ( spectype_depos(ispec) == SPECTYPE_DEPOS_NOT ) cycle
      write (gol,'(i6,a16,i6,a16,l5)') &
              ispec, specname(ispec), &
              spectype_depos(ispec), SPECTYPE_DEPOS_LABEL(spectype_depos(ispec)), &
              vd_depends_on_c(ispec); call goPr
    end do
    write (gol,'("------ --------------- ----- --------------- ----")'); call goPr
    write (gol,'("")'); call goPr
    
    ! * gamma water
    
    ! fill auxilary array:
    call LE_DryDepos_Gas_GammaWater_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! * meteo

    ! enable data:
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
!    call LE_Data_Enable( 'ssr', status )
!    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'ssrd', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'rain', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'sd', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'hp', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'sstk', status )
    IF_NOTOK_RETURN(status=1)
    !call LE_Data_Enable( 'srh', status )
    !IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'rh', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lat', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lai_lu', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'sai_lu', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'ustar', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'monin_inv', status )
    IF_NOTOK_RETURN(status=1)



    ! depac?
    if ( any( spectype_depos == SPECTYPE_DEPOS_GAS_DEPAC ) ) then
      ! enable meteo:
      call LE_Data_Enable( 'smi2', status )
      IF_NOTOK_RETURN(status=1)
      call LE_Data_Enable( 'smi3', status )
      IF_NOTOK_RETURN(status=1)

      ! data file for wet scavenging parameters
      ! warning: same file read twice: here and in wet_depos_camx !rj
      call ReadRc( rcF, 'le.wet_depos.CAMx.datafile', wetscav_datafile, status )
      IF_NOTOK_RETURN(status=1)
      
      ! file should be present:
      inquire( file=trim(wetscav_datafile), exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found : ",a)') trim(wetscav_datafile); call goErr
        TRACEBACK; status=1; return
      end if

      ! new file unit:
      call goGetFU( fu, status )
      IF_NOTOK_RETURN(status=1)
    
      ! open file:
      open( fu, file=trim(wetscav_datafile), status='old', form='formatted', iostat=status )
      if (status/=0) then
        write (gol,'("opening file : ",a)') trim(wetscav_datafile); call goErr
        TRACEBACK; status=1; return
      end if
    
      ! line counter:          
      iline = 0
      
      ! comment character:
      comment = '#'
      
      ! seperation character:
      sep = ';'
    
      ! read header line after first comment:
      do
        ! read line:
        read (fu,'(a)',iostat=status) line
        if (status/=0) then
          write (gol,'("reading header line from file : ",a)') trim(wetscav_datafile); call goErr
          TRACEBACK; status=1; return
        end if
        ! empty ? then skip:
        if ( len_trim(line) == 0 ) cycle
        ! comment ? then skip:
        if ( line(1:1) == comment ) cycle
        ! found non-comment line, leave loop:
        exit
      end do
    
      ! split:
      call goSplitString( line, nheader, headers, status, sep=sep )
      IF_NOTOK_RETURN(status=1)

      ! dummy values to avoid warnings ...
      ifield_tracer     = -999
      ifield_lbnd       = -999
      ifield_hlaw       = -999
      ifield_Tfac       = -999
      ifield_diffrat    = -999
      ifield_reactivity = -999
      ifield_Rscale     = -999 

      ! loop over headers
      do iheader = 1, nheader
        ! current:
        header = headers(iheader)
        ! which column ?
        select case ( trim(header) )
          ! tracer
          case ( 'Tracer' )
            ifield_tracer = iheader
          ! lower bnd
          case ( 'lower_bnd' )
            ifield_lbnd = iheader
          ! Henry laws constant
          case ( 'H_Law' )
            ifield_hlaw = iheader
          ! Temperature Fact
          case ( 'T_Fact' )
            ifield_Tfac = iheader
          ! Diffusion ratio
          case ( 'Diffrat' )
            ifield_diffrat = iheader
          ! reactivity
          case ( 'Reactivity' )
            ifield_reactivity = iheader
          ! Rscale
          case ( 'Rscale' )
            ifield_Rscale = iheader 
          ! unknown 
          case default
            write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
            TRACEBACK; status=1; return
        end select        
      end do 
    
      ! loop over records:
      do
    
        ! increase record counter:
        iline = iline + 1
        ! try to read line:
        read (fu,'(a)',iostat=status) line
        ! eof ?
        if (status<0) exit
        ! error ?
        if (status>0) then
          write (gol,'("reading line ",i6," from file : ",a)') iline, trim(wetscav_datafile); call goErr
          TRACEBACK; status=1; return
        end if
    
        ! empty ? then skip:
        if ( len_trim(line) == 0 ) cycle
        ! comment ? then skip:
        if ( line(1:1) == comment ) cycle
    
        ! split into records:
        call goSplitString( line, nfield, fields, status, sep=sep )
        IF_NOTOK_RETURN(status=1)
        ! check ...
        if ( nfield /= nheader ) then
          write (gol,'("number of fields (",i2,") in line ",i2," :")') nfield, iline; call goPr
          write (gol,'("  ",a)') trim(line); call goErr
          write (gol,'("fields:")'); call goErr
          do ifield = 1, nfield
            write (gol,'(i6," : ",a)') ifield, trim(fields(ifield)); call goErr
          end do
          write (gol,'("differs from number of headers ",i2," in file ",a)') nheader, trim(wetscav_datafile); call goErr
          TRACEBACK; status=1; return
        end if
      
        ! get tracer code
        call goMatchValue( trim(fields(ifield_tracer)), specname, ispec, status, quiet=.true. )
        if ( status /= 0 ) then
          write ( gol, '("Tracer in data file:, ",a, ", not in model run")' ) fields(ifield_tracer) ; call goPr
          cycle
        end if
        
        do ifield = 1, nfield
          
          ! skip tracer field
          if (ifield == ifield_tracer ) cycle
    
          field = fields(ifield)
          ! read field:
          read (field,*,iostat=status) value
          if (status/=0) then
            write (gol,'("reading fraction from `",a,"`")') field; call goErr
            TRACEBACK; status=1; return
          end if
          
          !~ henrys law
          if ( ifield == ifield_hlaw ) then
            henry0(ispec) = value
          !~ Temperature factor
          else if ( ifield == ifield_Tfac ) then
            Tfact(ispec)  = value
          ! reactivity
          else if ( ifield == ifield_reactivity ) then
            reactivity(ispec)   = value
    
          end if ! which header
        end do ! fields
      end do ! lines
      
      ! close file
      close( fu, iostat=status)
      if (status/=0) then
        write (gol,'("closing file : ",a)') trim(wetscav_datafile); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! which tracers are not filled with henry constant?
      do ispec = 1, nspec
        if ( henry0(ispec) == 0.0 ) then
          write( gol, '(" No henry diffusion coeff for tracer: ",a)') specname(ispec) ; call goPr
        endif
      enddo  

    end if ! depac

    ! *

    ! ok
    status = 0

  end subroutine LE_DryDepos_Init


  ! ***


  subroutine LE_DryDepos_Done( status )

    use LE_DryDepos_Gas_GammaWater, only : LE_DryDepos_Gas_GammaWater_Done

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_DryDepos_Done'

    ! --- begin ----------------------------------

    ! clear:
    deallocate( Rb )
    deallocate( Rc_eff )
    deallocate( Rc_total )
    deallocate( Rs )

    deallocate( vd_lu )
    deallocate( vs )

    if ( extra_output ) then
      deallocate( Cc_lu )
      deallocate( Fdry_total_lu )
      deallocate( Fdry_stom_lu )
      deallocate( Fdry_w_lu )
      deallocate( Fdry_soil_lu )
      deallocate( Cstom )
      deallocate( Cw )
      deallocate( Csoil )
      deallocate( gstom )
      deallocate( gw )
      deallocate( gsoil_eff )
      ! clear:
      deallocate( vd_zcanopytop_lu )
      deallocate( vd_htop_lu )
      deallocate( czcanopytop_o3_lu )
      deallocate( frac_sto_o3_lu )
      deallocate( fac_surface_area_2_PLA )
      deallocate( rbc_o3_lu )
    end if

    ! clear:
    deallocate( spectype_depos )
    deallocate( vd_depends_on_c )
    
    deallocate( henry0 )
    deallocate( Tfact )
    deallocate( reactivity )
       
    ! done with gammawater:
    call LE_DryDepos_Gas_GammaWater_Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_DryDepos_Done


  ! ***


  ! Compute mass budgets and new concentration due to deposition process;
  ! a compensation point makes it possible to have emission instead of deposition.
  ! Currently, only NH3 uses a compensation point, but depos has been made general
  ! in such a way that compensation points for other species are possible.
  !
  ! There are three time intervals for which budgets are kept;
  ! dimension (nx,ny,nspec,nex)
  !  bud%ex_day : daily budget
  !  bud%ex_hour: hourly budget
  !  bud_ex_now : budget for current time step
  ! 
  ! For dry deposition, separate budgets are computed per land use class:
  ! dimension (nx,ny,nspec,nlu_depac,nex)
  !  bud%ex_day_lu : daily budget per land use class
  !  bud%ex_hour_lu: hourly budget per land use class
  !  bud_ex_now_lu : budget for current time step per land use class
  !
  ! There are budgets for nex types of exchange processes; currently nex = 3:
  !  ex_drydepo = 1
  !  ex_dryemis = 2
  !  ex_wetdepo = 3
  !  ex_dry = (/ ex_drydepo, ex_dryemis /)
  !  bud_ex...(...,ex_drydepo): dry deposition
  !  bud_ex...(...,ex_dryemis): dry emission, in case there is a compensation point
  !  bud_ex...(...,ex_wetdepo): wet deposition
  !
  ! The "now" budgets are only present in subroutine depos; 
  ! hourly budgets are updated in subroutine depos;
  ! daily budgets are updated in le_groutput_update.

  subroutine LE_DryDepos_Apply( c, dt_min, bud, update_bud, status )

    use Binas             , only : xm_air, xm_C, xm_CO2
    use LE_Logging        , only : ident2
    use dims              , only : nx, ny, nz, nspec
    use dims              , only : rk1rk3
    use dims              , only : outF, runF
    use indices
    use LE_LandUse_Data   , only : Ra_h_htop_zcanopytop_lu
    use LE_LandUse_Data   , only : nlu, lu_fracs
    use LE_LandUse_Data   , only : fst_th_Y
    use LE_Budget_DryDepos, only : T_DryDepos_Budget
    use LE_Budget_DryDepos, only : nex, ex_drydepo, ex_dryemis
#ifdef with_labeling
    use SA_Labeling       , only : SA_DryDepos_Fractions, SA_Comp_point_budget
    use SA_Labeling       , only : SA_DryDepos_Fractions_lu_day
#endif    
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out -----------------------------------

    real, intent(inout)                      ::  c(nx,ny,nz,nspec)  ! atmospheric concentration (ppb or ug/m3)
    real, intent(in)                         ::  dt_min             ! time step (min)
    type(T_DryDepos_Budget), intent(inout)   ::  bud                ! budgets, a.o. per hour
    logical, intent(in)                      ::  update_bud
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_DryDepos_Apply'

    ! --- local ------------------------------------

    integer               ::  i, j
    integer               ::  k
    integer               ::  lu

    ! xxx is bud_ex_now also needed in wetdep?
    real, allocatable     ::  bud_ex_now(:,:,:,:)
    real, allocatable     ::  bud_ex_now_lu(:,:,:,:)

    ! pre-computed values:
    real, allocatable     ::  dh(:,:)
    real, allocatable     ::  inv_dh(:,:)
    real, allocatable     ::  dt_dh(:,:)
    ! for o3 stomata flux
    real, allocatable     ::  fst(:,:)

    ! meteo data:
    real, pointer          ::  h_m(:,:,:)    ! (lon,lat,1)
    real, pointer          ::  dens(:,:,:)   ! (lon,lat,alt)    
!    real, pointer          ::  rad(:,:,:)    ! (lon,lat,alt)    
    real, pointer          ::  radd(:,:,:)   ! (lon,lat,alt)

    ! --- begin ----------------------------
    
    call LE_Data_GetPointer( 'h' , h_m , status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3' )
    IF_NOTOK_RETURN(status=1)
!    call LE_Data_GetPointer( 'ssr' , rad , status, check_units='J/m2/s' )
!    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'ssrd', radd, status, check_units='J/m2/s' )
    IF_NOTOK_RETURN(status=1)

    !! info:
    !if ( .not. outF%suppress ) then
    !  write (gol,'(a,"<deposition>")') ident2; call goPr
    !end if

    ! storage:
    allocate( bud_ex_now   (nx,ny,nspec,nex) )
    allocate( bud_ex_now_lu(nx,ny,nspec,nlu) )
    allocate( dh(nx,ny) )
    allocate( inv_dh(nx,ny) )
    allocate( dt_dh(nx,ny) )

    ! pre-compute expressions used multiple times:
    !  tmp = 60.0*deltat/(h(:,:,1)*1000.0)
    dh = h_m(:,:,1)    ! m
    inv_dh = 1.0 / dh
    dt_dh = 60.0 * dt_min / dh

    ! * stomatal flux for O3
    
    ! ozone in use ?
    if ( i_o3 > 0 .and. extra_output ) then

      ! calculate landuse dependent concentration at canopy top
      call mix2canopytop(c(:,:,:,i_o3),czcanopytop_o3_lu)
     
      ! initialize stomatal flux at every time step
      allocate( fst(nx,ny) )
      fst = 0.0

      ! o3 stomatal flux exceedence calculated in order to infer the AFstY (Accumulated stomatal flux over threshold Y)
      ! for O3 (partly according to Simpson et al., 2007; Mills et al., 2004).
      ! BUT with DIFFERENT PARAMETRISATION FOR RB) and applying the constant flux approach (meaning that here we assume
      ! that the flux in layer 1 is equal to the flux at canopy height). ecjh. dd 31/03/2009.

      ! budget update per land use:
      do lu = 1, nlu
        
        
        where ( rbc_o3_lu(:,:,lu) > 0. )
          ! stomatal flux at current time step, units in nmol.m-2:
          ! Fstom (ctop) =         conc(ctop)               *     1/(rb+rc)        *   gs/(gs+gext+gsoil)   * factor surface area to PLA     * second in timestep
          !                                                                             frac_sto_o3_lu
          !        mol air /m3 air * 1e-9 mol O3/mol air    *       m/s            *            -           *          -                     * s
          fst    = 0.0409e+3       * czcanopytop_o3_lu(:,:,lu) * 1/ rbc_o3_lu(:,:,lu) * frac_sto_o3_lu(:,:,lu) * fac_surface_area_2_PLA(:,:,lu) * 60.0*dt_min 
        end where
        
        ! update budgets ?
        if ( update_bud ) then

          ! accumulated stomatal flux for o3
          where (radd(:,:,1) > 50.0 )
            bud%o3flx_stomata_hour_lu(:,:,lu) = bud%o3flx_stomata_hour_lu(:,:,lu) + fst
          end where

          ! accumulated o3 stomatal flux exceedence above threshold Y, 
          ! units in nmol.m-2 (over time period in seconds)
          if ( fst_th_Y(lu) >= 0.0 ) then
            where ( (fst > 60.0*dt_min*fst_th_Y(lu)) .and. (radd(:,:,1) > 50.0) )
              bud%o3exc_stomata_hour_lu(:,:,lu) = &
                  bud%o3exc_stomata_hour_lu(:,:,lu) + fst - 60.0*dt_min*fst_th_Y(lu)   ! nmol/m2
            end where
          end if
          
        end if ! update budgets

      enddo  ! depac landuse classes

      !! debug ...
      !print*,'frac_sto_o3(50,50,:)= ',frac_sto_o3(50,50,:)
      !print*,'fst(50,50,:)= ',fst(50,50,:)
      !print*,'acc_fst(50,50,:)= ',acc_fst(50,50,:)
      !print*,'acc_fst_exd(50,50,:)= ',acc_fst_exd(50,50,:)

    end if   ! O3 in use
    
    ! *

    ! initialize on zero for present deposition budgets    
    bud_ex_now = 0.0
    ! loop over species:
    do k = 1, nspec
      ! only the deposed tracers ...
      if ( spectype_depos(k) == SPECTYPE_DEPOS_NOT ) cycle

      ! loop over grid cells:
      !$OMP parallel &
      !$OMP   default( none ) &
      !$OMP   shared( k ) &
      !$OMP   shared( nx, ny ) &
      !$OMP   shared( nlu ) &
      !$OMP   shared( lu_fracs ) &
      !$OMP   shared( update_bud ) &
      !$OMP   shared( bud_ex_now ) &
      !$OMP   shared( bud_ex_now_lu ) &
      !$OMP   shared( bud ) &
      !$OMP   shared( ccomp_tot_lu ) &
      !$OMP   shared( vd_lu ) &
      !$OMP   shared( c ) &
      !$OMP   shared( dh, inv_dh, dt_dh ) &
      !$OMP   private( i, j ) &
      !$OMP   private( lu )
      !$OMP   do SCHEDULE( DYNAMIC )
      do j = 1, ny
        do i = 1, nx

          ! Loop over land use classes:
          do lu = 1,nlu
      
            ! Initialise mass exchange budgets of current time step:
            bud_ex_now_lu(i,j,k,lu) = 0.0

            ! Compute mass flux integrated over this time step in surface layer, 
            ! i.e. exchanged mass per area in this time step (ppb.m or ug/m2), for each land use type:
            ! Note: bud_ex_now_lu < 0, if c > ccomp_tot (deposition)
            !       bud_ex_now_lu > 0, if c < ccomp_tot (emission)
            ! compute mass budget for each land use class:
            bud_ex_now_lu(i,j,k,lu) = - (c(i,j,1,k) - ccomp_tot_lu(i,j,k,lu)) * &
                                       (1.0 - exp(-vd_lu(i,j,k,lu)*dt_dh(i,j) ))*dh(i,j)

            ! 1. Add contribution of this time step and this land use class to 
            !    total grid cell budgets (deposition or emission) in current time step
            ! 2. Add contribution of this time step and this land use class to 
            !    hourly land use class budgets (deposition or emission) 
            ! Note: depos < 0, emis > 0 xxx afspreken
            if ( bud_ex_now_lu(i,j,k,lu) <= 0.0 ) then
              bud_ex_now    (i,j,k,ex_drydepo) = bud_ex_now(i,j,k,ex_drydepo) + &
                                                 bud_ex_now_lu(i,j,k,lu)*lu_fracs(i,j,lu)
            else
              bud_ex_now    (i,j,k,ex_dryemis) = bud_ex_now(i,j,k,ex_dryemis) + &
                                                    bud_ex_now_lu(i,j,k,lu)*lu_fracs(i,j,lu)
            end if

            ! update budgets ?
            if ( update_bud ) then
              if ( bud_ex_now_lu(i,j,k,lu) <= 0.0 ) then
                bud%ex_hour_lu(i,j,k,lu,ex_drydepo) = bud%ex_hour_lu(i,j,k,lu,ex_drydepo) + &
                                                      bud_ex_now_lu(i,j,k,lu)
              else
                bud%ex_hour_lu(i,j,k,lu,ex_dryemis) = bud%ex_hour_lu(i,j,k,lu,ex_dryemis) + &
                                                      bud_ex_now_lu(i,j,k,lu)
              end if
            end if
            
          end do  ! landuse classes
          
          ! update budgets ?
          if ( update_bud ) then
            ! Add contribution of this time step to hourly grid cell budgets (drydepo and dryemis):
            bud%ex_hour(i,j,k,:) = bud%ex_hour(i,j,k,:) + bud_ex_now(i,j,k,:)
          end if  ! update budgets

          ! Compute new concentration, add contribution of all exchange processes:
          c(i,j,1,k) = c(i,j,1,k) + sum(bud_ex_now(i,j,k,:)) * inv_dh(i,j)

        end do   ! i
      end do  ! j
      !$OMP   end do
      !$OMP end parallel

    end do   ! loop over species

#ifdef with_labeling
    ! calculate fractions in deposition budgets
    call SA_DryDepos_Fractions( bud_ex_now, status )
    IF_NOTOK_RETURN(status=1)
    ! calculate fractions of new landuse dependent deposition budgets
    call SA_DryDepos_Fractions_lu_day( bud_ex_now_lu, status )
    IF_NOTOK_RETURN(status=1)
#endif         

    !! sec. inorg. aerosol enabled ?
    !if ( n_sia > 0 ) then
    !  ! limit the aerosols
    !  ! AJS: This was probalby necessary to avoid problems with issoropia
    !  !  or some other sia routine. Problem: this is an artifical source!
    !  if (i_so4a_f>0) where (c(:,:,1,i_so4a_f) < 1.0e-4) c(:,:,1,i_so4a_f) = 1.0e-4
    !  if (i_nh4a_f>0) where (c(:,:,1,i_nh4a_f) < 1.0e-4) c(:,:,1,i_nh4a_f) = 1.0e-4
    !  if (i_no3a_f>0) where (c(:,:,1,i_no3a_f) < 1.0e-4) c(:,:,1,i_no3a_f) = 1.0e-4
    !  ! AJS: don't do this for course SO4a and NO3a yet;
    !  ! keep it zero, such that results with and without these new tracers
    !  ! are the same.
    !  !where (c(:,:,1,i_so4a_c) < 1.0e-4) c(:,:,1,i_so4a_c) = 1.0e-4
    !  !where (c(:,:,1,i_no3a_c) < 1.0e-4) c(:,:,1,i_no3a_c) = 1.0e-4
    !end if

    ! Update summed NH3 concentration (needed for average concentration for DEPAC):
    if ( i_nh3 > 0 .and. i_so2 > 0 ) then
      ! update budgets ?
      if ( update_bud ) then
        bud%cnh3_sum  = bud%cnh3_sum + c(:,:,1,i_nh3) ! layer 1 
        bud%cnh3_nsum = bud%cnh3_nsum + 1
        bud%cso2_sum  = bud%cso2_sum + c(:,:,1,i_so2) ! layer 1 
        bud%cso2_nsum = bud%cso2_nsum + 1
        !if (outF%debug_print) then
        !  ix = outF%debug_ix
        !  iy = outF%debug_iy
        !  if (ispec .eq. i_nh3) then
        !    write(*,'(a,1x,a6,4(1x,i4),3(1x,e12.5))') &
        !            'drydep sum; spec,ix,iy,lu,cnh3_nsum,c,cnh3_sum,cnh3_sum/cnh3_nsum: ', &
        !            specname(i_nh3),ix,iy,0,cnh3_nsum,c(ix,iy,1,i_nh3),cnh3_sum(ix,iy),cnh3_sum(ix,iy)/cnh3_nsum
        !  endif
        !endif
      end if  ! update budgets
#ifdef with_labeling
      ! update budgets for compensation point in labeling
      ! .false. to state that it is not end of the month
      call SA_Comp_point_budget(  .false., status, conc_sfc = c(:,:,1,:) )
      IF_NOTOK_RETURN(status=1)
#endif      
    endif  ! NH3 and SO2 present

    !! Debug output:
    !if (outF%debug_print) then
    !  ix = outF%debug_ix
    !  iy = outF%debug_iy
    !  do k = 1,nspec
    !    if (spectype_depos(k) .ne. spectype_depos_not) then
    !      do lu = 1,nlu_depac
    !        if (ispec .eq. i_nh3) then
    !          write(*,'(a,1x,a6,3(1x,i4),1x,e12.5)') &
    !                  'drydep end1; spec,ix,iy,lu,bud_ex_now_lu: ', &
    !                  specname(k),ix,iy,lu,bud_ex_now_lu(ix,iy,k,lu)
    !          write(*,'(a,1x,a6,3(1x,i4),2(1x,e12.5))') &
    !                  'drydep end2; spec,ix,iy,lu,drydepo_hour,dryemis_hour: ', &
    !                  specname(k),ix,iy,lu,bud_ex_hour_lu(ix,iy,k,lu,ex_dry)
    !        endif
    !      enddo
    !      if (ispec .eq. i_nh3) then
    !        write(*,'(a,1x,a6,3(1x,i4),2(1x,e12.5))') &
    !                'drydep end3; spec,ix,iy,lu,drydepo_now,dryemis_now: ', &
    !                 specname(k),ix,iy,0,bud_ex_now(ix,iy,k,ex_dry)
    !        write(*,'(a,1x,a6,3(1x,i4),2(1x,e12.5))') &
    !                'drydep end4; spec,ix,iy,0,drydepo_hour,dryemis_hour: ', &
    !                specname(k),ix,iy,0,bud_ex_hour(ix,iy,k,ex_dry)
    !        write(*,'(a,1x,a6,3(1x,i4),1x,e12.5)') &
    !                'drydep end5; spec,ix,iy,lu,c(1): ', &
    !                specname(k),ix,iy,0,c(ix,iy,1,k)
    !      endif
    !    endif
    !  enddo
    !endif
    !
    !! xxx debug printing
    !outF%debug_print = .false.

    ! clear:
    deallocate( bud_ex_now )
    deallocate( bud_ex_now_lu )
    deallocate( dh )
    deallocate( inv_dh )
    deallocate( dt_dh )
    
    ! ok
    status = 0

  end subroutine LE_DryDepos_Apply

  ! ***


  ! Get deposition velocity and compensation point
  !
  ! get Ra, Rb
  ! loop over grid cells
  ! |  get meteo variables
  ! |  loop over grid cells
  ! |  |   1. compute vd and ccomp_tot for DEPAC species (if needed)
  ! |  |   2. compute vd for fine aerosol (if needed)
  ! |  |   3. compute vd for coarse aerosol (if needed)
  ! |  |   4. loop over species
  ! |  |      |   retrieve vd (and compensation point) for either DEPAC species, similar species, 
  ! |  |      |    fine aerosol of coarse aerosol
  ! |  |      end loop over species
  ! |  end loop over land use classes
  ! end loop over grid cells

  subroutine LE_DryDepos_Setup_vd( c, bud, only_if_depends_on_c, t, status )

    use Binas                     , only : xm_air
    use JAQL_drydeposition        , only : Sedimentation_Velocity

    use LE_Logging                , only : ident2
    use GO                        , only : TDate, DayNumber, Get
    use GO                        , only : goUpCase
    use GO                        , only : GO_Timer_Start, GO_Timer_End
    use Dims                      , only : nx, ny, nz, nspec
    use Dims                      , only : lst,coszen
    use Dims                      , only : outF, runF
    use indices
    use LE_Landuse_Data           , only : nlu, lu_fracs, lu_name
    use LE_Landuse_Data           , only : Ld               ! leaf dimension
    use LE_Landuse_Data           , only : Ra_h_htop_z0_lu
    use LE_Landuse_Data           , only : Ra_h_zcanopytop_z0_lu
    use LE_Landuse_Data           , only : u_zcanopytop_lu
    use LE_Landuse_Data           , only : ustar_lu
    use LE_Landuse_Data           , only : ilu_grass, ilu_arable, ilu_permanent_crops
    use LE_Landuse_Data           , only : ilu_coniferous_forest, ilu_deciduous_forest
    use LE_Landuse_Data           , only : ilu_water_sea, ilu_urban, ilu_other, ilu_desert
    use LE_Landuse_Data           , only : ilu_ice, ilu_savanna,ilu_tropical_forest, ilu_water_inland
    use LE_Landuse_Data           , only : ilu_mediterrean_scrub
    use LE_Landuse_Data           , only : ustar_lu
    use LE_Landuse_Data           , only : ludep_output_whole_grid
    use LE_Landuse_Data           , only : with_ozone_specials

    use LE_Particle_Data          , only : partsize
    use LE_Particle_Data          , only : slipcor
    !use LE_Particle_Data          , only : LE_Particle_Data_Update
    use LE_Sedim                  , only : slinnfac
    use LE_Meteo_data             , only : viscos

#ifdef with_m7
    use LE_M7_drydepo             , only : LE_M7_DryDepo_Setup, LE_M7_DryDepo_vd,m7_vd_table
    use LE_M7                     , only : rwetm7modes, densm7modes
#endif
    use LE_DryDepos_Aero_WER      , only : vdpart_fine, vdpart_coarse
    use LE_DryDepos_Aero_Zhang    , only : LE_DryDepos_Aero_Zhang_vd
#ifdef with_pollen    
    use LE_DryDepos_Aero_Pollen   , only : LE_DryDepos_Aero_Pollen_vd
#endif
    use LE_DryDepos_Gas_DEPAC     , only : DryDepos_Gas_DEPAC
    use LE_DryDepos_Gas_DEPAC     , only : Get_N_S_Regime, iratns_default
    use LE_DryDepos_Gas_GammaWater, only : gamma_water  ! (nx,ny) gamma soil for water
    use LE_Budget_DryDepos        , only : T_DryDepos_Budget

    use LE_Data                   , only : LE_Data_GetPointer
    use dims                      , only : substantial_snowdepth
   
    use JAQL                      , only : henry_func
    use Binas                     , only : Rgas

    !use JAQL, only : JAQL_BL_Rb_McNaughton_VdHurk
   
    ! --- in/out ---------------------------------

    ! atmospheric concentration (ppb (gas) or ug/m3 (aerosol))
    ! changed by m7 remoding
    real, intent(in)                      ::  c(nx,ny,nz,nspec)
    ! budgets, used for NH3 compensation point:
    type(T_DryDepos_Budget), intent(in)   ::  bud
    ! only if depends on c ?
    logical, intent(in)                   ::  only_if_depends_on_c
    ! time:
    type(TDate), intent(in)               ::  t
    ! return status:
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_DryDepos_Setup_vd'
    
    ! --- local ---------------------------------

    integer :: i                    ! index of grid cell (x-direction)
    integer :: j                    ! index of grid cell (y-direction)
    integer :: k                    ! index of land use class
    integer :: ispec                ! species index
    integer :: iday                 ! day of year
    integer :: iratns               ! index for NH3/SO2 ratio;
                                    !   iratns = 1: low NH3/SO2
                                    !   iratns = 2: high NH3/SO2
                                    !   iratns = 3: very low NH3/SO2
    integer :: nwet                 ! wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow; nwet=10 -> ice
    real    :: glrad                ! global radiation (W/m2)
    real    :: ps                   ! surface pressure (Pa)
    real    :: ts                   ! surface temperature (C)
    real    :: tsea                 ! sea surface temperature (K)
    real    :: rh_surf              ! relative humidity at surface
    real    :: ol                   ! 1/L, L = Monin Obukhov length (m)
    real    :: ust                  ! friction velocity (m/s)
    real    :: lat                  ! latitude (degrees)
    real    :: coszen1              ! cos(zenith angle) = sin(solar elevation)
    real    :: ppb_to_ugm3          ! conversion factor

    real    :: Rc_tot               ! total canopy resistance (s/m)

    character(len=100)  :: specnam1 ! species name
    real, allocatable   ::  diffc(:)       ! (nspec)
    real    :: catm                 ! actual atmospheric concentration (ug/m3)
    real    :: c_ave_prev_nh3       ! atmospheric concentration average over previous period (ug/m3)
    real    :: c_ave_prev_so2       ! atmospheric concentration average over previous period (ug/m3)
    real    :: ccomp_tot_lu1        ! compensation point from DEPAC (ug/m3)
    
    real    ::  Rc_temp_o3_sto      ! canopy resistance of ozone due to stomata
    real    ::  Rc_temp             ! canopy resistance
    
    real    ::  gstot, gs, gscu, gsoeff
    real    ::  gc_tot
    real    ::  hlaw                ! Henry's law constant (M atm-1) 
    real    ::  react               ! reactivity (-)
    ! default pH of rain water, somewhat acid due to presence of CO2: !rj: to be improved for water on surface 
    real, parameter             ::  pH_rain = 5.0

    integer             ::  imode
    integer             ::  ispecN
    integer             ::  i_aerosol
#ifdef with_pollen    
    integer             ::  i_pollen
#endif    
    real    ::  smi
    
    integer   ::  status_par
    
    ! meteo data:
    real, pointer               ::  rain (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  tsurf(:,:,:)   ! (lon,lat,1)
    real, pointer               ::  rh   (:,:,:)   ! (lon,lat,lev)
    real, pointer               ::  temp (:,:,:)   ! (lon,lat,lev)
    real, pointer               ::  hpres(:,:,:)   ! (lon,lat,lev)
    real, pointer               ::  dens (:,:,:)   ! (lon,lat,lev)
    real, pointer               ::  sst  (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  radd (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  sd   (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  smi2 (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  smi3 (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  lats (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  ustar (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  monin_inv (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  z0h_lu    (:,:,:)   ! (lon,lat,nlu)
    real, pointer               ::  lai_lu    (:,:,:)   ! (lon,lat,nlu)
    real, pointer               ::  sai_lu    (:,:,:)   ! (lon,lat,nlu)

    ! --- begin ----------------------------
    
    ! start timing:
    call GO_Timer_Start( itim_vdepos, status )
    IF_NOTOK_RETURN(status=1)
    
    !! info ...
    !write (gol,'("compute deposition velocities ...")'); call goPr
    !if ( any(ludep_output_whole_grid) ) then
    !  write (gol,'("  depositions on whole grid for landuses:")'); call goPr
    !  do k = 1, nlu
    !    if ( ludep_output_whole_grid(k) ) then
    !      write (gol,'("    ",a)') trim(lu_name(k)); call goPr
    !    end if
    !  end do
    !else
    !  write (gol,'("  depositions only for present landuses")'); call goPr
    !end if
          
    ! get pointer to meteo:
    call LE_Data_GetPointer( 'rain', rain, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 't', temp, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'hp', hpres, status, check_units ='Pa' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'sstk', sst, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'ssrd', radd, status, check_units='J/m2/s' )
    IF_NOTOK_RETURN(status=1)
    !call LE_Data_GetPointer( 'srh', sRH, status, check_units='%' )
    !IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'rh' , RH , status, check_units='%' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'sd', sd, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lat', lats, status, check_units='degrees_north' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'ustar', ustar, status, check_units='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'monin_inv', monin_inv, status, check_units='1/m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'z0h_lu', z0h_lu, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lai_lu', lai_lu, status, check_units='m2/m2' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'sai_lu', sai_lu, status, check_units='m2/m2' )
    IF_NOTOK_RETURN(status=1)


    ! depac?
    if ( any( spectype_depos == SPECTYPE_DEPOS_GAS_DEPAC ) ) then
      ! get pointer to meteo:
      call LE_Data_GetPointer( 'smi2', smi2, status, check_units='1' )
      IF_NOTOK_RETURN(status=1)
      call LE_Data_GetPointer( 'smi3', smi3, status, check_units='1' )
      IF_NOTOK_RETURN(status=1)
    end if ! depac

    ! set the day number in the year
    iday = DayNumber( t )

    ! init sum:
    gstot = 0.0

    !>>> now called in LE_Driver/LE_TimeStep_Run ,
    !    thus all derived aerosol parameters (particle size, slip correction)
    !    are assumed to be up-to-date now:
    !!! update particle sizes, slip correction factors etc (if necessary):
    !call LE_Particle_Data_Update( c, status )
    !IF_NOTOK_RETURN(status=1)
    !<<<

#ifdef with_m7
    ! start timing:
    call GO_Timer_Start( itim_vdepos_m7, status )
    IF_NOTOK_RETURN(status=1)
    ! setup M7 look-up table if necessary:
    call LE_M7_DryDepo_Setup( status )
    IF_NOTOK_RETURN(status=1)
    ! end timing:
    call GO_Timer_End( itim_vdepos_m7, status )
    IF_NOTOK_RETURN(status=1)
#endif
    
    ! ~ sedimentation velocity

    ! start timing:
    call GO_Timer_Start( itim_vdepos_vs, status )
    IF_NOTOK_RETURN(status=1)

    ! any aerosol deposition following Zhang ?
    if ( any( spectype_depos == SPECTYPE_DEPOS_AERO_ZHANG ) &
#ifdef with_pollen        
        .or. any( spectype_depos == SPECTYPE_DEPOS_AERO_POLLEN ) & 
#endif
        ) then

      ! loop over aerosol tracers:
      do i_aerosol = 1, n_aerosol
        ! global index:
        ispec = ispecs_aerosol(i_aerosol)
        ! sedimentation velocity in lowest layer:
        vs(1:nx,1:ny,ispec) = slinnfac(ispec) * Sedimentation_Velocity( rhopart(ispec), &
                                                        partsize(1:nx,1:ny,1,i_aerosol), &
                                                        slipcor (1:nx,1:ny,1,i_aerosol), &
                                                        viscos  (1:nx,1:ny,1) )
      end do  ! aerosol tracers
    end if

    ! end timing:
    call GO_Timer_End( itim_vdepos_vs, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ~ Rc

    ! start timing:
    call GO_Timer_Start( itim_vdepos_rc, status )
    IF_NOTOK_RETURN(status=1)
   
    ! parameter indicating the NH3/SO2 ratio regime ;
    ! both SO2 and NH3 present ?
    if ( (i_so2 > 0) .and. (i_nh3 > 0) ) then
      ! set parameter:
      ns_regime = Get_N_S_Regime( c(:,:,1,i_nh3), c(:,:,1,i_so2) )
    else
      ! default low ratio:
      ns_regime = iratns_default
    end if

    ! storage:
    allocate( diffc(nspec), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! use default value for HNO3, N2O5, NO3, H2O2
    diffc            = 0.11e-4
    ! overwrite with known coefficients of diffusivity from Massman (1998)
    if ( ispec_no2 > 0 ) diffc(ispec_no2) = 0.136e-4 
    if ( ispec_no  > 0 ) diffc(ispec_no ) = 0.199e-4
    if ( ispec_o3  > 0 ) diffc(ispec_o3 ) = 0.144e-4
    if ( ispec_co  > 0 ) diffc(ispec_co ) = 0.176e-4
    if ( ispec_so2 > 0 ) diffc(ispec_so2) = 0.112e-4
    if ( ispec_ch4 > 0 ) diffc(ispec_ch4) = 0.191e-4
    if ( ispec_nh3 > 0 ) diffc(ispec_nh3) = 0.191e-4

    ! loop over all species:
    do ispec = 1, nspec

      ! no need to update ?
      if ( only_if_depends_on_c .and. (.not. vd_depends_on_c(ispec)) ) cycle

      ! switch ...
      select case ( spectype_depos(ispec) )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( SPECTYPE_DEPOS_NOT )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! just skip ..

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( SPECTYPE_DEPOS_GAS_DEPAC )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! Compute Rc (DEPAC):
          specnam1 = goUpCase(trim(specname(ispec))) ! DEPAC expects upper case name
          
          ! init parallel status:
          status_par = 0
    
          ! loop over cells:
          !$OMP parallel &
#ifndef __GFORTRAN__
          !$OMP   default( none ) &
          !$OMP   shared( Ld ) &
          !$OMP   shared( specmolm ) &
#endif
          !$OMP   shared( ispec, specnam1 ) &
          !$OMP   shared( nx,ny ) &
          !$OMP   shared( nlu ) &
          !$OMP   shared( i_so2, i_nh3, i_o3 ) &
          !$OMP   shared( c ) &
          !$OMP   shared( ns_regime ) &
          !$OMP   shared( tsurf, sd, rain ) &
          !$OMP   shared( rh, radd, hpres, temp, sst, monin_inv, ustar, coszen ) &
          !$OMP   shared( dens) &
          !$OMP   shared( bud ) &
          !$OMP   shared( lu_fracs ) &
          !$OMP   shared( diffc ) &
          !$OMP   shared( u_zcanopytop_lu ) &
          !$OMP   shared( ustar_lu, z0h_lu ) &
          !$OMP   shared( ludep_output_whole_grid ) &
          !$OMP   shared( smi2, smi3 ) &
          !$OMP   shared( lats ) &
          !$OMP   shared( extra_output, with_ozone_specials ) &
          !$OMP   shared( iday ) &
          !$OMP   shared( lai_lu, sai_lu ) &
          !$OMP   shared( Ra_h_htop_z0_lu, Rb, Rc_eff, Rc_total ) &
          !$OMP   shared( gamma_water ) &
          !$OMP   shared( gw, gstom, gsoil_eff, cw, cstom, csoil ) &
          !$OMP   shared( ccomp_tot_lu ) &
          !$OMP   shared( vd_lu ) &
          !$OMP   shared( Cc_lu, Fdry_total_lu, Fdry_stom_lu, Fdry_w_lu, Fdry_soil_lu ) &
          !$OMP   shared( frac_sto_o3_lu, fac_surface_area_2_PLA ) &
          !$OMP   shared( vd_htop_lu, vd_zcanopytop_lu ) &
          !$OMP   shared( rbc_o3_lu ) &
          !$OMP   shared( Ra_h_zcanopytop_z0_lu ) &
          !$OMP   shared( henry0, Tfact, reactivity ) &

          !$OMP   private( i, j ) &
          !$OMP   private( iratns ) &
          !$OMP   private( nwet ) &
          !$OMP   private( rh_surf, glrad, ps, ts, tsea, ol, ust, coszen1 ) &
          !$OMP   private( lat ) &
          !$OMP   private( ppb_to_ugm3 ) &
          !$OMP   private( catm ) &
          !$OMP   private( c_ave_prev_nh3, c_ave_prev_so2 ) &
          !$OMP   private( k ) &
          !$OMP   private( smi ) &
          !$OMP   private( Rc_tot, ccomp_tot_lu1 ) &
          !$OMP   private( status ) &
          !$OMP   private( hlaw, react ) &
          !$OMP   reduction( + : status_par )
          !$OMP   do
          do j = 1, ny
            do i = 1 ,nx

              ! parameter indicating the NH3/SO2 ratio regime:
              iratns = ns_regime(i,j)

              !dry surface by default
              nwet = 0
              ! if it's freezing we assume dry surface or snow cover
              if (tsurf(i,j,1) < 273.15) then
                nwet = 0
                ! icheck for snow cover (assume below -4 C, icy surface)
                if ( sd(i,j,1) > substantial_snowdepth .and. tsurf(i,j,1) >  269.15 ) nwet=9
                if ( sd(i,j,1) > substantial_snowdepth .and. tsurf(i,j,1) <= 269.15 ) nwet=10
                                
              else
                ! wet, in case of rain
                if ( rain(i,j,1) > 0.0 ) nwet=1
                ! icheck for snow cover
                if ( sd(i,j,1) > substantial_snowdepth ) nwet=9
                !add wetness as function of RH or fog in the future
              endif

              ! Other meteo and grid variables:
              ! ~ relative humidity in lowest layer:
              rh_surf =  rh(i,j,1)  ! lowest layer of 3D field
              !rh_surf = srh(i,j,1)  ! based on t2m and d2m; should be better?
              ! ~ radiation:
              glrad   = radd(i,j,1)        
              ! ~ surface pressure:
              ps      = hpres(i,j,0)
              ! ~ surface temperature in degrees Celcius:
              ts      = temp(i,j,1) - 273.15 ! degrees celcius

              ! ~ sea-surface temperature:
              tsea = sst(i,j,1)  ! K
              ! trap cases where tsea over land is undefined;
              ! value is probably 0.0 then, reset to 2m temperature:
              if ( tsea < 200.0 ) tsea = tsurf(i,j,1)

              ! ~ stability parameters:
              ol      = monin_inv(i,j,1)
              ust     = ustar(i,j,1)
              ! ~ latitude :
              lat = lats(i,j,1)
              ! ~ cosine of solar zenith angle:
              coszen1 = coszen(i,j)

              ! Factor used for conversion from ppb to ug/m3 :
              !   ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
              !    c           1e-9                 xm_tracer         1e9       /       xm_air            dens
              ! thus:
              !    c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
              ! Use density at lowest layer:
              ppb_to_ugm3 =  dens(i,j,1) / xm_air  ! (mole air)/m3

              ! atmospheric concentration in DEPAC is requested in ug/m3:
              !  ug/m3              ppb          (ug/m3)/ppb/(kg/mole)     kg/mole
              catm           = c(i,j,1,ispec)         * ppb_to_ugm3 *   specmolm(ispec)  ! in ug/m3

              ! for compensation point:
              if ( ispec == i_nh3 ) then
                c_ave_prev_nh3 = bud%cnh3_ave_prev(i,j) * ppb_to_ugm3 * specmolm(i_nh3)  ! in ug/m3                  
                c_ave_prev_so2 = bud%cso2_ave_prev(i,j) * ppb_to_ugm3 * specmolm(i_so2)  ! in ug/m3                  
              else
                c_ave_prev_nh3 = 0.0
                c_ave_prev_so2 = 0.0
              endif

              ! calculate effective henry law constant !rj
              hlaw = henry_func( ispec, henry0(ispec), Tfact(ispec), tsurf(i,j,1), pH_rain ) 
              react=reactivity(ispec)

              ! Loop over land use classes:
              do k = 1, nlu

                ! skip if land use class is not present in current grid cell:
                ! not skip, if landuse dependent output is wanted on the whole grid
                if ( lu_fracs(i,j,k) <= 0.0 .and. ( .not.(ludep_output_whole_grid(k)) )) cycle

                ! get boundary layer resistance:
                call Get_Rb_cell( (k == ilu_water_sea) .or. (k == ilu_water_inland), Ld(k), &
                                  z0h_lu(i,j,k), ustar_lu(i,j,k), u_zcanopytop_lu(i,j,k), &
                                  diffc(ispec), &
                                  Rb(i,j,ispec,k) )

                ! determine till which depth the soil moisture index (SMI) needs to be taken into account
                ! smi2 => 0-28 cm
                ! smi3 => 0-100 cm
                ! which type of landuse, to determine soil moisture index
                if ( k == ilu_grass .or. &
                     k == ilu_other) then
                     ! .or. &
                     !k == ilbg_semi_natural_grass ) then
                  smi = smi2(i,j,1)
                else if ( k == ilu_arable .or. &
                          k == ilu_permanent_crops .or. &
                          k == ilu_coniferous_forest .or. &
                          k == ilu_deciduous_forest .or. &
                          k == ilu_water_sea .or. &
                          k == ilu_urban .or. &
                          k == ilu_desert .or. &
                          k == ilu_ice .or. &
                          k == ilu_tropical_forest .or. &
                          k == ilu_water_inland ) then
                          !.or. &
                          !k == ilu_spruce ) then
                  smi = smi3(i,j,1)
                else
                  !write (gol, '("Unknown landuse class: ", i2)' ) k ; call GoPr
                  !TRACEBACK; status=1; return
                  status_par = status_par + 1
                end if
               
                if ( ispec == i_o3 ) then
                  ! For ozone extra calculation of stomatal ozone fluxes, to put in budget)

                  if ( extra_output .or. with_ozone_specials ) then
                    ! Extra output, conductances and compensation points: external leaf, effective soil and stomata
                    call DryDepos_Gas_DEPAC( specnam1, iday, lat, ts, ust, glrad, coszen1, &
                             rh_surf, lai_lu(i,j,k), sai_lu(i,j,k), nwet, k, iratns, &
                             Rc_tot, ccomp_tot_lu1, &
                             hlaw, react, status, &
                             p=ps, tsea=tsea, smi=smi, &
                             c_ave_prev_nh3=c_ave_prev_nh3, c_ave_prev_so2=c_ave_prev_so2, &
                             catm=catm, gamma_soil_water=gamma_water(i,j), &
                             Ra=Ra_h_htop_z0_lu(i,j,k), Rb=Rb(i,j,ispec,k), Rc_eff=Rc_eff(i,j,ispec,k), &
                             calc_stom_o3flux       = .true. , &
                             frac_sto_o3_lu         = frac_sto_o3_lu(i,j,k), &
                             fac_surface_area_2_PLA = fac_surface_area_2_PLA(i,j,k), &
                             gw_out        = gw(i,j,ispec,k), &
                             gstom_out     = gstom(i,j,ispec,k), &
                             gsoil_eff_out = gsoil_eff(i,j,ispec,k), &
                             cw_out        = cw(i,j,ispec,k), &
                             cstom_out     = cstom(i,j,ispec,k), &
                             csoil_out     = csoil(i,j,ispec,k) )
                    if ( status /= 0 ) then
                      write (gol,'("cell (",i0,",",i0,") landuse ",i0)') i, j, k; call goErr
                      TRACEBACK; status=1; return
                    end if
                            
                    if (Rc_tot.le.0.0) Rc_tot = 9999.
                    ! ~ velocity:
                    vd_htop_lu(i,j,k) = 1.0/(Ra_h_htop_z0_lu(i,j,k)      +Rb(i,j,ispec,k)+Rc_tot)
                    vd_zcanopytop_lu(i,j,k) = 1.0/(Ra_h_zcanopytop_z0_lu(i,j,k)+Rb(i,j,ispec,k)+Rc_tot)
                    !  total surface resistance (Rb+Rc) needed for C at canopy top :
                    rbc_o3_lu(i,j,k) = Rb(i,j,ispec,k) + Rc_tot

                    ! Compute Cc (see le_drydepos_gas_depac.F90) and component fluxes
                    Cc_lu(i,j,ispec,k) = (catm/(Ra_h_htop_z0_lu(i,j,k)+Rb(i,j,ispec,k))+ &
                                          Cstom(i,j,ispec,k)*gstom(i,j,ispec,k)+Cw(i,j,ispec,k)*gw(i,j,ispec,k)+Csoil(i,j,ispec,k)*gsoil_eff(i,j,ispec,k))/ &
                                          (1/(Ra_h_htop_z0_lu(i,j,k)+Rb(i,j,ispec,k))+gstom(i,j,ispec,k)+gw(i,j,ispec,k)+gsoil_eff(i,j,ispec,k)) 
                    Fdry_total_lu(i,j,ispec,k)  = -(catm-Cc_lu(i,j,ispec,k))/(Ra_h_htop_z0_lu(i,j,k)+Rb(i,j,ispec,k))
                    Fdry_stom_lu(i,j,ispec,k)   = -(Cc_lu(i,j,ispec,k)-Cstom(i,j,ispec,k))*gstom(i,j,ispec,k)
                    Fdry_w_lu(i,j,ispec,k)      = -(Cc_lu(i,j,ispec,k)-Cw(i,j,ispec,k))*gw(i,j,ispec,k)
                    Fdry_soil_lu(i,j,ispec,k)   = -(Cc_lu(i,j,ispec,k)-Csoil(i,j,ispec,k))*gsoil_eff(i,j,ispec,k)

                  else ! only calculation of deposition resistances
                    call DryDepos_Gas_DEPAC( specnam1, iday, lat, ts, ust, glrad, coszen1, &
                             rh_surf, lai_lu(i,j,k), sai_lu(i,j,k), nwet, k, iratns, &
                             Rc_tot, ccomp_tot_lu1, &
                             hlaw, react, status, &
                             p=ps, tsea=tsea, smi=smi, &
                             c_ave_prev_nh3=c_ave_prev_nh3, c_ave_prev_so2=c_ave_prev_so2, &
                             catm=catm, gamma_soil_water=gamma_water(i,j), &
                             Ra=Ra_h_htop_z0_lu(i,j,k), Rb=Rb(i,j,ispec,k), Rc_eff=Rc_eff(i,j,ispec,k) )
                    if ( status /= 0 ) then
                      write (gol,'("cell (",i0,",",i0,") landuse ",i0)') i, j, k; call goErr
                      TRACEBACK; status=1; return
                    end if
                  end if 
                 
                 else ! other tracers

                  if ( extra_output ) then
                    ! extra output, conductances and compensation points: external leaf, effective soil and stomata
                    call DryDepos_Gas_DEPAC( specnam1, iday, lat, ts, ust, glrad, coszen1, &
                                 rh_surf, lai_lu(i,j,k), sai_lu(i,j,k), nwet, k, iratns, &
                                 Rc_tot, ccomp_tot_lu1, &
                                 hlaw, react, status, &
                                 p=ps, tsea=tsea, smi=smi, &
                                 c_ave_prev_nh3=c_ave_prev_nh3, c_ave_prev_so2=c_ave_prev_so2, &
                                 catm=catm, gamma_soil_water=gamma_water(i,j), &
                                 Ra=Ra_h_htop_z0_lu(i,j,k), Rb=Rb(i,j,ispec,k), Rc_eff=Rc_eff(i,j,ispec,k), &
                                 gw_out        = gw(i,j,ispec,k), &
                                 gstom_out     = gstom(i,j,ispec,k), &
                                 gsoil_eff_out = gsoil_eff(i,j,ispec,k), &
                                 cw_out        = cw(i,j,ispec,k), &
                                 cstom_out     = cstom(i,j,ispec,k), &
                                 csoil_out     = csoil(i,j,ispec,k) )
                    if ( status /= 0 ) then
                      write (gol,'("cell (",i0,",",i0,") landuse ",i0)') i, j, k; call goErr
                      TRACEBACK; status=1; return
                    end if

                    ! Compute Cc (see le_drydepos_gas_depac.F90) and component fluxes
                    Cc_lu(i,j,ispec,k) = (catm/(Ra_h_htop_z0_lu(i,j,k)+Rb(i,j,ispec,k))+ &
                                          Cstom(i,j,ispec,k)*gstom(i,j,ispec,k)+Cw(i,j,ispec,k)*gw(i,j,ispec,k)+Csoil(i,j,ispec,k)*gsoil_eff(i,j,ispec,k))/ &
                                          (1/(Ra_h_htop_z0_lu(i,j,k)+Rb(i,j,ispec,k))+gstom(i,j,ispec,k)+gw(i,j,ispec,k)+gsoil_eff(i,j,ispec,k)) 
                    Fdry_total_lu(i,j,ispec,k)  = -(catm-Cc_lu(i,j,ispec,k))/(Ra_h_htop_z0_lu(i,j,k)+Rb(i,j,ispec,k))
                    Fdry_stom_lu(i,j,ispec,k)   = -(Cc_lu(i,j,ispec,k)-Cstom(i,j,ispec,k))*gstom(i,j,ispec,k)
                    Fdry_w_lu(i,j,ispec,k)      = -(Cc_lu(i,j,ispec,k)-Cw(i,j,ispec,k))*gw(i,j,ispec,k)
                    Fdry_soil_lu(i,j,ispec,k)   = -(Cc_lu(i,j,ispec,k)-Csoil(i,j,ispec,k))*gsoil_eff(i,j,ispec,k)
                  else
                     call DryDepos_Gas_DEPAC( specnam1, iday, lat, ts, ust, glrad, coszen1, &
                                 rh_surf, lai_lu(i,j,k), sai_lu(i,j,k), nwet, k, iratns, &
                                 Rc_tot, ccomp_tot_lu1, &
                                 hlaw, react, status, &
                                 p=ps, tsea=tsea, smi=smi, &
                                 c_ave_prev_nh3=c_ave_prev_nh3, c_ave_prev_so2=c_ave_prev_so2, &
                                 catm=catm, gamma_soil_water=gamma_water(i,j), &
                                 Ra=Ra_h_htop_z0_lu(i,j,k), Rb=Rb(i,j,ispec,k), Rc_eff=Rc_eff(i,j,ispec,k) )
                    if ( status /= 0 ) then
                      write (gol,'("cell (",i0,",",i0,") landuse ",i0)') i, j, k; call goErr
                      TRACEBACK; status=1; return
                    end if

                  end if                                   

                end if ! ozone or not
  
                ! Convert compensation point back to ppb:
                ccomp_tot_lu(i,j,ispec,k) = ccomp_tot_lu1/(ppb_to_ugm3*specmolm(ispec))

                ! Set Rc to large value (no pathway) if Rc negative:
                if (Rc_tot .le. 0.0) Rc_tot = 9999.

                ! store:
                Rc_total(i,j,ispec,k) = Rc_tot

                ! Compute exchange velocity for land use class k:
                vd_lu(i,j,ispec,k) = 1.0 / ( Ra_h_htop_z0_lu(i,j,k)+ Rb(i,j,ispec,k) + Rc_tot )

              end do   ! landuse classes

            end do  ! i
          end do  ! j
          !$OMP end do
          !$OMP end parallel

          ! check status ...
          if ( status_par /= 0 ) then
            write (gol,'("error status returned from parallel loop over tracers : ",i6)') status_par; call goErr
            TRACEBACK; status=1; return
          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( SPECTYPE_DEPOS_AERO_ZHANG )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! init parallel status:
          status_par = 0
    
          ! loop over cells:
          !$OMP parallel &
#ifndef __GFORTRAN__
          !$OMP   default( none ) &
          !$OMP   shared( ispec2aerosol ) &
#endif
          !$OMP   shared( ispec ) &
          !$OMP   shared( nx,ny ) &
          !$OMP   shared( tsurf, sd, rain ) &
          !$OMP   shared( partsize ) &
          !$OMP   shared( Rs, vs, slipcor, viscos ) &
          !$OMP   shared( dens) &
          !$OMP   shared( Ra_h_htop_z0_lu) &
          !$OMP   shared( ustar_lu) &
          !$OMP   shared( nlu ) &
          !$OMP   shared( lu_fracs ) &
          !$OMP   shared( ludep_output_whole_grid ) &
          !$OMP   shared( vd_lu ) &
          !$OMP   shared( lai_lu ) &
          !$OMP   private( i, j ) &
          !$OMP   private( i_aerosol ) &
          !$OMP   private( k ) &
          !$OMP   private( nwet ) &
          !$OMP   private( status ) &
          !$OMP   reduction( + : status_par )
          !$OMP   do
          do j = 1, ny
            do i = 1 ,nx

              !dry surface by default
              nwet = 0
              ! if it's freezing we assume dry surface or snow cover
              if (tsurf(i,j,1) < 273.15) then
                nwet = 0
                ! icheck for snow cover
                if ( sd(i,j,1) > substantial_snowdepth .and. tsurf(i,j,1) >  269.15 ) nwet=9
                if ( sd(i,j,1) > substantial_snowdepth .and. tsurf(i,j,1) <= 269.15 ) nwet=10
              else
                ! wet, in case of rain
                if ( rain(i,j,1) > 0.0 ) nwet=1
                !add wetness as function of RH or fog in the future
                ! icheck for snow cover (assume below -4 C, icy surface)
                if ( sd(i,j,1) > substantial_snowdepth ) nwet=9
              endif

              ! aerosol index:
              i_aerosol = ispec2aerosol(ispec)
              
              ! For M7 tracer the particle size could be zero,
              ! if it could not be computed due to lack of aerosol mass and number;
              ! this leads to problems with division by zero when computing Schmidt nubmer.
              ! Does not occur for the fine/coarse tracers, for these always
              ! a particle size is defined independent if aerosol is present or not.
              ! Only call Zhang vd routine for non-zero particle size:
              if ( partsize(i,j,1,i_aerosol) > 0.0  ) then

                ! Aerosol deposition velocities.
                ! Loop over land use classes:
                do k = 1, nlu
                  ! skip if land use class is not present in current grid cell:
                  if ( lu_fracs(i,j,k) <= 0.0 .and. ( .not.(ludep_output_whole_grid(k)) )) cycle
                  ! for this landuse class only ;
                  ! sedimentation velocity in lowest layer is input:

                  call LE_DryDepos_Aero_Zhang_vd( vd_lu(i,j,ispec,k), Rs(i,j,ispec,k), &
                                   vs(i,j,ispec), &
                                   partsize(i,j,1,i_aerosol), &
                                   slipcor (i,j,1,i_aerosol), &
                                   nwet, tsurf(i,j,1), dens(i,j,1), viscos(i,j,1), &
                                   k, lai_lu(i,j,k), Ra_h_htop_z0_lu(i,j,k), ustar_lu(i,j,k), &
                                   status )
                  status_par = status_par + status
                end do   ! landuse classes
                
              else
              
                ! M7 zero particle size, set vd to zero for all landuse classes:
                vd_lu(i,j,ispec,:) = 0.0
                
              end if
            
            end do  ! i
          end do  ! j
          !$OMP end do
          !$OMP end parallel

          ! check status ...
          if ( status_par /= 0 ) then
            write (gol,'("error status returned from parallel loop over tracers : ",i6)') status_par; call goErr
            TRACEBACK; status=1; return
          end if

#ifdef with_pollen
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( SPECTYPE_DEPOS_AERO_POLLEN )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
          ! init parallel status:
          status_par = 0
    
          ! loop over cells:
          !$OMP parallel &
#ifndef __GFORTRAN__
          !$OMP   default( none ) &
          !$OMP   shared( ispec2aerosol ) &
#endif
          !$OMP   shared( ispec ) &
          !$OMP   shared( nx,ny ) &
          !$OMP   shared( nlu ) &
          !$OMP   shared( lu_fracs ) &
          !$OMP   shared( ludep_output_whole_grid ) &
          !$OMP   shared( vd_lu ) &
          !$OMP   shared( partsize ) &
          !$OMP   shared( vs ) &
          !$OMP   shared( Ra_h_htop_z0_lu ) &
          !$OMP   private( i, j ) &
          !$OMP   private( k ) &
          !$OMP   private( i_aerosol ) &
          !$OMP   private( status ) &
          !$OMP   reduction( + : status_par )
          !$OMP   do
          do j = 1, ny
            do i = 1 ,nx

              ! aerosol index:
              i_aerosol = ispec2aerosol(ispec)
              
              ! avoid particle size 0
              if ( partsize(i,j,1,i_aerosol) > 0.0 ) then
                
                ! pollen deposition velocities
                ! Loop over land use classes:
                do k = 1, nlu
                  ! skip if land use class is not present in current grid cell:
                  if ( lu_fracs(i,j,k) <= 0.0 .and. ( .not.(ludep_output_whole_grid(k)) )) cycle
                  ! for this landuse class only ;
                  ! sedimentation velocity in lowest layer is input:
                  call LE_DryDepos_Aero_Pollen_vd( vd_lu(i,j,ispec,k), &
                                   vs(i,j,ispec), &
                                   Ra_h_htop_z0_lu(i,j,k), &
                                   status )
                  !IF_NOTOK_RETURN(status=1)
                  status_par = status_par + status
                end do   ! landuse classes            
              
              end if ! partsize check

            end do  ! i
          end do  ! j
          !$OMP end do
          !$OMP end parallel

          ! check status ...
          if ( status_par /= 0 ) then
            write (gol,'("error status returned from parallel loop over tracers : ",i6)') status_par; call goErr
            TRACEBACK; status=1; return
          end if
#endif

#ifdef with_m7              
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( SPECTYPE_DEPOS_AERO_M7_MASS )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
          ! init parallel status:
          status_par = 0
    
          ! loop over cells:
          !$OMP parallel &
          !$OMP   default( none ) &
          !$OMP   shared( ispec ) &
          !$OMP   shared( nx,ny ) &
          !$OMP   shared( nlu ) &
          !$OMP   shared( c) &
          !$OMP   shared( lu_fracs ) &
          !$OMP   shared( ludep_output_whole_grid ) &!
          !$OMP   shared( vd_lu ) &         
          !$OMP   shared( m7_mode, m7_numberconc_ispecs, rwetm7modes, densm7modes, m7_vd_table ) &   
          !$OMP   private( i_aerosol, imode, ispecN ) &
          !$OMP   private( i, j ) &
          !$OMP   private( k ) &
          !$OMP   private( status ) &
          !$OMP   reduction( + : status_par )
          !$OMP   do
          do j = 1, ny
            do i = 1 ,nx

              ! NOTE: some overhead accepted here ; the vd only depends
              ! on the mode number, but is here re-computed for each M7 mass tracer
            
              ! M7 mode number:
              imode = m7_mode(ispec)
              ! tracer index with corresponding number density:
              ispecN = m7_numberconc_ispecs(imode)
            
              ! Loop over land use classes:
              do k = 1, nlu
               ! skip if land use class is not present in current grid cell:
                if ( lu_fracs(i,j,k) <= 0.0 .and. ( .not.(ludep_output_whole_grid(k)) )) cycle
                ! fill both mass and number velocity:
                call LE_M7_DryDepo_vd( &
                            vd_lu(i,j,ispec ,k)     , &   ! out: vd for aerosol mass
                            vd_lu(i,j,ispecN,k)     , &   ! out: vd for related aerosol number
                            imode                   , &   ! mode number
                            rwetm7modes(i,j,1,imode), &   ! wet radius
                            densm7modes(i,j,1,imode), &   ! density
                            c(i,j,1,ispecN)         , &   ! number densitity of this mode
                            m7_vd_table(i,j,:,k)    , &   ! lookup table for this landuse
                            status )
                !IF_NOTOK_RETURN(status=1)
                status_par = status_par + status
              end do

            end do  ! i
          end do  ! j
          !$OMP end do
          !$OMP end parallel

          ! check status ...
          if ( status_par /= 0 ) then
            write (gol,'("error status returned from parallel loop over tracers : ",i6)') status_par; call goErr
            TRACEBACK; status=1; return
          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( SPECTYPE_DEPOS_AERO_M7_NUMBER )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! skip, already done together with the aerosol mass
#endif              

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! error ...
          write (gol,'("unsupported spectype_depos ",i6," for tracer ",i6," ",a)') &
                     spectype_depos(ispec), ispec, trim(specname(ispec)); call goErr
          TRACEBACK; status=1; return

      end select

      ! adhoc factor used for assimilation:
      if ( ispec == ispec_o3 ) then
        ! loop over land use classes:
        do k = 1, nlu
          ! apply factor:
          vd_lu(:,:,ispec_o3,k) = vd_lu(:,:,ispec_o3,k) * vd_o3fac
        end do
      end if

    end do  ! tracers
    
    ! clear:
    deallocate( diffc, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! end timing:
    call GO_Timer_End( itim_vdepos_rc, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ~

    ! end timing:
    call GO_Timer_End( itim_vdepos, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_DryDepos_Setup_vd


  ! ***


  !
  ! Compute quasi-laminar boundary layer resistance as a function of landuse and tracer.
  ! If the leaf dimension (Ld) is known, the formulation by McNaughton and van der Hurk (1995) is used:
  !                  
  !  Rb = thk/diffc(ispec) x 150 x sqrt[Ld / V(h)]
  !
  ! Here: - V(h) (at canopy top h) is the landuse dependent windspeed at the top of the canopy.
  !       - thk/diffc(ispec) is the ratio of the diffusivity of heat and the diffusivity of the component 
  !         under consideration (e.g. for O3 the ratio is 1.3)
  !       - Ld is the leaf dimension
  ! If the leaf dimension is not known or no leaves are present (e.g. sea), the original EMEP formulation 
  ! (Simpson et al, 2003) is used. 
  !
  ! Usage:
  !   call Get_Rb_cell( is_water, Ld(k), &
  !                     z0h_lu(i,j,k), ustar_lu(i,j,k), u_zcanopytop_lu(i,j,k), &
  !                     diffc(ispec), &
  !                     Rb(i,j,ispec,k) )
  !

  pure subroutine Get_Rb_cell( is_water, Ld, z0h, ustar, u_zcanopytop, diffc, Rb )

    use Binas,            only : kappa_stab       ! von Karman constant

    ! --- in/out ---------------------------------
    
    logical, intent(in)   ::  is_water
    real, intent(in)      ::  Ld            ! leaf dimension
    real, intent(in)      ::  z0h           ! roughness length for heat
    real, intent(in)      ::  ustar         ! friction velocity
    real, intent(in)      ::  u_zcanopytop  ! wind speed at development height htop
    real, intent(in)      ::  diffc         ! coefficient of diffusivity
    real, intent(out)     ::  Rb            ! boundary layer resistance

    ! --- const ----------------------------------
   
    ! thermal diffusivity of dry air 20 C
    real, parameter       ::  thk = 0.19e-4         ! http://en.wikipedia.org/wiki/Thermal_diffusivity  (@ T=300K)
    
    !  --- begin ---------------------------------
    
    if ( Ld > 0.0 ) then
      ! use McNaughton and van der Hurk (1995)
      Rb = (thk/diffc) * 150.0 * sqrt( Ld /max(0.01,u_zcanopytop) )
    else
      ! use Simpson et al. (2003) 
      if ( is_water ) then
        Rb = 1.0 / (kappa_stab*max(0.01,ustar)) * log(z0h/diffc*kappa_stab*max(0.01,ustar))
        Rb = max(Rb,0.0)
      else
        Rb = 5.0 / max(0.01,ustar) * (thk/diffc)**0.67
      end if
    end if

  end subroutine Get_Rb_cell

      
  ! ====================================================================

  
  ! mix2ground computes a ground concentration (i.e. at measuring height)
  ! assuming a constant flux from height htop to the ground.
  ! htop is the height at which the deposition velocity is computed;
  ! for nz > 3, layer height of layer 1 = 25 m (see module stability): 
  ! htop=max(5.0, min(htopmx, 1./elinv, h(i,j,1)*500.0)) ! thus min(50, L, 12.5)
  ! Note that the measuring height is currently set in the module stability
  ! at 2.5 m: xxx 3.5 m is better?
  !  real, parameter :: hrfub = 2.5

  ! FS, 2009-11-24 update mix2ground for compensation point.
  !     New documentation.
  !     Placed inside module deposition, since it relates closely to the deposition scheme used.

  subroutine mix2ground( c, cg, status )
  
    use Dims, only      : nx, ny, nz, nspec
    use dims, only      : outF
    use LE_Logging, only : ident2
    use LE_Landuse_Data, only   : nlu, lu_fracs
    use LE_Landuse_Data, only   : Ra_h_htop_zobs_lu
    
    ! --- in/out ---------------------------------

    real, intent(in)      ::  c(nx,ny,nz,nspec) ! atmospheric concentration
    real, intent(out)     ::  cg(nx,ny,nspec)   ! concentration at measuring height
    integer, intent(out)  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/mix2ground'

    ! --- local ----------------------------------
    
    real    :: c_href_lu(nx,ny)  ! concentration at reference level (measuring height) for 1 land use class in ppb
    real    :: c_href_ave(nx,ny) ! concentration at reference level averaged over land use classes in ppb
    integer :: ispec             ! species index
    integer :: lu                ! land use class index

    ! --- begin ----------------------------------

    !if (.NOT.outF%suppress) print *, ident2,'<computing ground concentrations>'

    ! flux F = -vd*(Catm - Ccomp), 
    ! vd   : deposition velocity (m/s), 
    ! Catm : atmospheric concentration (ppb or ug/m3)
    ! Ccomp: compensation point (ppb or ug/m3)
    !
    ! For the relevant part htop - href = measuring height, the flux is
    ! -(Catm - C(href))/Ra(href).
    !
    ! Constant flux from htop to ground implies:
    ! -Catm + C(href) = Ra(href)*F = -Ra(href)*vd*(Catm - Ccomp),
    ! so C(href) = Catm - Ra(href)*vd*(Catm - Ccomp).
    ! and C(href) = Ccomp + (Catm-Ccomp)(1-Ra(href)*vd)
    !
    ! Note1: if we have no compensation point, we have:
    ! C(href) = Catm -Ra(href)*vd*Catm = Catm*(1-Ra(href)*vd).
    ! Note2: Ra(href), i.e. the resistance between htop and href = measuring height,
    ! is computed in module stability, variable ikz.
    !
    ! xxx In subroutine stability a grid averaged z0 is used; it would be better to compute
    ! htop and ikz for each land use class separately.

    ! Loop over species:
    do ispec = 1,nspec 

      ! Check whether species is deposited or not:
      if ( spectype_depos(ispec) == SPECTYPE_DEPOS_NOT ) then
        cg(:,:,ispec) = c(:,:,1,ispec) ! same as atmospheric concentration
      else if (( spectype_depos(ispec) == SPECTYPE_DEPOS_AERO_ZHANG) &
           .or.( spectype_depos(ispec) == SPECTYPE_DEPOS_AERO_WER) &
#ifdef with_pollen
           .or.( spectype_depos(ispec) == SPECTYPE_DEPOS_AERO_POLLEN) &
#endif          
            ) then
        !old method, no compensation point, average deposition velocity
        !cg(:,:,ispec) = c(:,:,1,ispec)*(1.0-ikz*vd(:,:,ispec))

        ! Initialise:
        c_href_ave = 0.0

        ! loop over land use classes:
        do lu = 1,nlu
          ! Compute concentration at reference level (measuring height) for current land use type:
          c_href_lu = ccomp_tot_lu(:,:,ispec,lu) + &                    
                      (c(:,:,1,ispec) - ccomp_tot_lu(:,:,ispec,lu)) * (1+ vs(:,:,ispec)*Ra_h_htop_zobs_lu(:,:,lu) - Ra_h_htop_zobs_lu(:,:,lu)*vd_lu(:,:,ispec,lu))
          ! Add contribution to average concentration:
          c_href_ave = c_href_ave + lu_fracs(:,:,lu)*c_href_lu
        enddo

        ! Assign average concentration to output array:
        cg(:,:,ispec) = c_href_ave
      else
        !old method, no compensation point, average deposition velocity
        !cg(:,:,ispec) = c(:,:,1,ispec)*(1.0-ikz*vd(:,:,ispec))

        ! Initialise:
        c_href_ave = 0.0

        ! loop over land use classes:
        do lu = 1,nlu
          ! Compute concentration at reference level (measuring height) for current land use type:
          c_href_lu = ccomp_tot_lu(:,:,ispec,lu) + &                    
                      (c(:,:,1,ispec) - ccomp_tot_lu(:,:,ispec,lu)) * (1.0 - Ra_h_htop_zobs_lu(:,:,lu)*vd_lu(:,:,ispec,lu))
          ! Add contribution to average concentration:
          c_href_ave = c_href_ave + lu_fracs(:,:,lu)*c_href_lu
        enddo

        ! Assign average concentration to output array:
        cg(:,:,ispec) = c_href_ave
      end if

    end do
    
    ! ok
    status = 0

  end subroutine mix2ground
  
  ! ====================================================================
  
  subroutine mix2ground_lu( c, cg, status )
  
    use Dims, only      : nx, ny, nz, nspec
    use LE_Logging, only : ident2
    use LE_Landuse_Data, only   : nlu
    use LE_Landuse_Data, only   : Ra_h_htop_zobs_lu

    ! --- in/out ---------------------------------

    real, intent(in)      ::  c(nx,ny,nz,nspec) ! atmospheric concentration
    real, intent(out)     ::  cg(nx,ny,nspec,nlu)   ! concentration at measuring height
    integer, intent(out)  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/mix2ground_lu'

    ! --- local ----------------------------------
    
    real    :: c_href_lu(nx,ny)  ! concentration at reference level (measuring height) for 1 land use class in ppb
    integer :: ispec             ! species index
    integer :: lu                ! land use class index

    ! --- begin ----------------------------------

    !if (.NOT.outF%suppress) print *, ident2,'<computing ground concentrations>'

    ! flux F = -vd*(Catm - Ccomp), 
    ! vd   : deposition velocity (m/s), 
    ! Catm : atmospheric concentration (ppb or ug/m3)
    ! Ccomp: compensation point (ppb or ug/m3)
    !
    ! For the relevant part htop - href = measuring height, the flux is
    ! -(Catm - C(href))/Ra(href).
    !
    ! Constant flux from htop to ground implies:
    ! -Catm + C(href) = Ra(href)*F = -Ra(href)*vd*(Catm - Ccomp),
    ! so C(href) = Catm - Ra(href)*vd*(Catm - Ccomp).
    ! and C(href) = Ccomp + (Catm-Ccomp)(1-Ra(href)*vd)
    !
    ! Note1: if we have no compensation point, we have:
    ! C(href) = Catm -Ra(href)*vd*Catm = Catm*(1-Ra(href)*vd).
    ! Note2: Ra(href), i.e. the resistance between htop and href = measuring height,
    ! is computed in module stability, variable ikz.
    !
    ! xxx In subroutine stability a grid averaged z0 is used; it would be better to compute
    ! htop and ikz for each land use class separately.

    ! Loop over species:
    do ispec = 1,nspec 

      ! Check whether species is deposited or not:
      if ( spectype_depos(ispec) == SPECTYPE_DEPOS_NOT ) then
        do lu = 1, nlu
          cg(:,:,ispec,lu) = c(:,:,1,ispec) ! same as atmospheric concentration
        end do
      else if (( spectype_depos(ispec) == SPECTYPE_DEPOS_AERO_ZHANG) &
           .or.( spectype_depos(ispec) == SPECTYPE_DEPOS_AERO_WER)   &
#ifdef with_pollen
           .or.( spectype_depos(ispec) == SPECTYPE_DEPOS_AERO_POLLEN)&
#endif      
           ) then
        ! loop over land use classes:
        do lu = 1,nlu
          ! Compute concentration at reference level (measuring height) for current land use type:
          c_href_lu = ccomp_tot_lu(:,:,ispec,lu) + &
                       (c(:,:,1,ispec) - ccomp_tot_lu(:,:,ispec,lu)) * (1+ vs(:,:,ispec)*Ra_h_htop_zobs_lu(:,:,lu) - Ra_h_htop_zobs_lu(:,:,lu)*vd_lu(:,:,ispec,lu))          
          ! Assign average concentration to output array:
          cg(:,:,ispec,lu) = c_href_lu
        enddo
      else 
        ! loop over land use classes:
        do lu = 1,nlu
          ! Compute concentration at reference level (measuring height) for current land use type:
          c_href_lu = ccomp_tot_lu(:,:,ispec,lu) + &
                      (c(:,:,1,ispec) - ccomp_tot_lu(:,:,ispec,lu)) * (1.0 - Ra_h_htop_zobs_lu(:,:,lu)*vd_lu(:,:,ispec,lu))
          ! Assign average concentration to output array:
          cg(:,:,ispec,lu) = c_href_lu
        enddo

      end if

    end do
    
    ! ok
    status = 0

  end subroutine mix2ground_lu
  
  ! ====================================================================

  
  subroutine mix2canopytop( c, czcanopytop_o3_lu )
  
    use indices
    use Dims, only                : nx, ny, nz, nspec
    use LE_Logging, only          : ident2
    use LE_Landuse_Data, only     : nlu
    use LE_Landuse_Data, only     : Ra_h_htop_zcanopytop_lu

    ! --- in/out ---------------------------------

    real, intent(in)    :: c(nx,ny,nz)             ! atmospheric concentration of ozone
    real, intent(out)   :: czcanopytop_o3_lu(nx,ny,nlu)  ! concentration at canopy top

    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/mix2canopytop'

    ! --- local ----------------------------------
    
    ! Local variables:
    integer :: lu                               ! land use class index

    ! --- begin ----------------------------------

    !if (.NOT.outF%suppress) print *, ident2,'<computing canopy top concentration>'

    ! flux F = -vd_lu(htop)*(C(htop) - Ccomp), 
    ! vd_lu   : deposition velocity (m/s), 
    ! C(htop) : atmospheric concentration (ppb or ug/m3)
    !
    ! Constant flux from htop to canopytop implies:
    ! -C(htop) + C(zcanopytop) = Ra_h_htop_zcanopytop_lu*F = -Ra_h_htop_zcanopytop_lu*vd_lu*(C(htop) - Ccomp),
    ! so C(zcanopytop) = C(htop) - Ra_h_htop_zcanopytop_lu*vd_lu*(C(htop) - Ccomp).
    !
    ! For the relevant part htop - zcanopytop = canopy top height, the flux is
    ! -(C(htop) - C(zcanopytop))/Ra_h_htop_zcanopytop_lu.
    !
    ! Note1: if we have no compensation point, we have:
    ! C(zcanopytop) = C(htop) -Ra_h_htop_zcanopytop_lu*vd_lu*C(htop) = C(htop)*(1-Ra_htop_zcanopytop_lu)*vd_lu).
    ! Note2: Ra_h_htop_zcanopytop_lu, i.e. the resistance between htop and zcanopytop = canopy top height,
    ! is computed in module stability, variable Ra_h_htop_zcanopytop_lu.
    !
    ! RWK, 2013-09-03 created for special ozon output.
         
    ! loop over land use classes:
    do lu = 1,nlu
      czcanopytop_o3_lu(:,:,lu)= c(:,:,1)  * (1.0 - Ra_h_htop_zcanopytop_lu(:,:,lu)*vd_htop_lu(:,:,lu))
    enddo

  end subroutine mix2canopytop

end module LE_DryDepos

