!#################################################################
!
! Source-Apportion Labeling module
!
! LIBRARY
!
!   Requires the Sparse-BLAS library:
!     http://math.nist.gov/spblas/
!
!   Current implementation based on the version included in the
!   Math-Kernel-Library (MKL) which comes with the Intel compiler suite.
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#################################################################

module SA_Labeling

  use NetCDF,           only :  NF90_StrError, NF90_NOERR
  use GO,               only :  gol, goErr, goPr
  use GO,               only :  TDate

  use LE_Meteo_Data,    only :  ovolume
  use dims,             only :  kz
  use dims,             only :  nx,ny,nz,nspec
  use dims,             only :  runF

  use indices

  implicit none


  ! --- in/out ------------------------------

  private

  public  ::  SA_nlabel
  public  ::  SA_Label_Names, SA_Short_Label_Names, labelname_len
  public  ::  SA_frac
  public  ::  SA_frac_adv
  public  ::  SA_frac_hdiff
  public  ::  SA_Emis
  public  ::  SA_conc_comp_point, SA_frac_comp_point
  public  ::  SA_budget_wetdepos_hour, SA_budget_wetdepos_day
  public  ::  SA_budget_drydepos_hour, SA_budget_drydepos_day, SA_budget_drydepos_lu_day
  public  ::  Label_Aloft
  public  ::  Labelled_specs, labelled_specs_names, nspec_labelled
  public  ::  labelled_specs_budget_output_names, landuse_budget_output_names
  public  ::  nbudget_label_output_tracer, nbudget_label_output_lu
  public  ::  labelled_specs_drydep_output_names, labelled_specs_wetdep_output_names
  public  ::  ndrydep_label_output_tracer, nwetdep_label_output_tracer
  public  ::  SA_Synchronize_Conc

  public  ::  SA_Labeling_Init, SA_Labeling_Done
  public  ::  SA_Restart_Restore
  public  ::  SA_Frac_Init
  public  ::  SA_Adjust_Setup, SA_Adjust_Conc, SA_Adjust_Fractions
  public  ::  SA_Advec_Setup
  public  ::  SA_Advec_x, SA_Advec_y, SA_Advec_z
  public  ::  SA_HDiff_Setup, SA_HDiff_Apply
  public  ::  SA_Bound_Setup, SA_Bound_Done
  public  ::  SA_Chem_Gas_Setup, SA_Chem_Gas_Iter, SA_Chem_Gas_Step, SA_Chem_SIA, SA_Chem_Cloud
  public  ::  SA_Emis_Setup_Dust, SA_Emis_Setup_Natural, SA_Emis_Setup_Fire, SA_Emis_Setup_TNO, SA_Emis_Setup_EDGAR, SA_Emis_Setup_Soil_NOx
  public  ::  SA_Emis_Sectors, SA_Emis_Sectors_Done
  public  ::  SA_Emis_Reset, SA_Emis_Delta_c
  public  ::  SA_Sedim_Setup, SA_Sedim_Conc
  public  ::  SA_DryDepos_Fractions, SA_DryDepos_Reset
  public  ::  SA_WetDepos_Fractions, SA_WetDepos_Reset
  public  ::  SA_DryDepos_Fractions_lu_day, SA_Depos_Reset_day
  public  ::  SA_Comp_point_Init, SA_Comp_point_Done
  public  ::  SA_Comp_point_budget
  public  ::  SA_Vdiff_Setup, SA_Vdiff_Conc, SA_Vdiff_Fractions



  ! --- const ------------------------------

  character(len=*), parameter             ::  mname = 'SA_Labeling'

  ! For numerical reasons concentrations may not be too small
  ! Therefore a numerical limit concentration is set (numlim_conc) for which numerical problems won't occur

  ! assume at least one molecule per m3
  ! top of the model at least 700 hPa
  ! maximum temperature 330 K
  ! Gas constant 8.314472
  ! Avogadro's number 6.022 e23
  ! Then the volume of 1 mol gas at least is 8.314472 * 330 / 700 = 0.039 m3
  ! 1 mol gas has 6.022e23 molecules, thus at most 6.022e23/0.039 = 1.54e25 molecules per m3
  ! Thus a concentration of 1 molecule per m3 means 1 molecule per 1.54e25 molecules
  ! This is equal to 1.54e-16 ppb
  !
  ! for the tracers which are set in micrograms 1 molecule per m3
  ! compares with at least 1/6.022e23  = 1.66e-24 g / m3
  !                                    = 1.66e-21 ug / m3

  ! For clarity reasons, the numerical limit for all concentrations is set to 1e-21

  real, parameter                         ::  numlim_conc = 1e-21
  integer, parameter                      ::  max_emis_set = 5
  
  ! --- types ------------------------------
  
  ! emission sectors ( store for all emission sets)
  type TEmis_sectors
    character(len=100),allocatable,dimension(:) ::  emis_sector_names
    integer                                     ::  nsector
    integer                                     ::  isector_resuspension, isector_agriculture    
  end type

  ! emission countries ( store for all emission sets)
  type TEmis_countries
    character(len=32),allocatable,dimension(:)  ::  emis_country_names
    integer                                     ::  ncountry    
  end type
  

  ! --- local ------------------------------
  integer                                 ::  SA_nlabel
  integer                                 ::  def_label
  integer                                 ::  Label_Bound, Label_Init
  integer                                 ::  Label_Aloft, Label_Natural

  integer, pointer                        ::  Labelled_specs(:)
  character(len=16), pointer              ::  labelled_specs_names(:)
  integer                                 ::  nspec_labelled

  character(len=512), pointer             ::  SA_Label_Names(:)
  character(len=512), pointer             ::  SA_Short_Label_Names(:)
  integer                                 ::  labelname_len

  logical                                 ::  First_def_step

  ! Concentrations and fractions
  real, allocatable, dimension(:,:,:,:,:) ::  SA_frac
  real, allocatable, dimension(:,:,:,:)   ::  LE_conc

  ! Emissions
  real, allocatable, dimension(:,:,:,:,:) ::  SA_emis         ! Added emission in each time step
  real, allocatable, dimension(:,:,:,:,:) ::  SA_flux_emis
  
  ! Emission sets
  character(len=64)                       ::  emis_set_names(max_emis_set)
  integer                                 ::  nemis_set      ! real number of emission sets
  type(TEmis_Sectors)                     ::  Emis_Sectors( max_emis_set )
  type(TEmis_Countries)                   ::  Emis_Countries( max_emis_set )

  ! Advection
  real, allocatable                       ::  SA_frac_adv(:,:,:,:,:)
  
  ! Horizontal diffusion
  real, allocatable                       ::  SA_frac_hdiff(:,:,:,:,:)
  
  ! Adjust
  real, allocatable, dimension(:)         ::  SA_caloft_frac
  real, allocatable, dimension(:,:)       ::  int_m_SA
  real, allocatable, dimension(:)         ::  SA_mass
  real, allocatable, dimension(:,:)       ::  SA_conc_adj

  ! Vert diffusion
  real, allocatable, dimension(:,:,:)     ::  SA_conc_vdiff

  ! Sedimentation
  real, allocatable, dimension(:,:)       ::  SA_frac_sedim
  real, allocatable, dimension(:)         ::  SA_conc_sedim
  
  ! Dry Deposition fluxes on hourly basis
  real, allocatable, dimension(:,:,:,:)   ::  SA_budget_drydepos_hour
  ! Wet Deposition fluxes on hourly basis
  real, allocatable, dimension(:,:,:,:)   ::  SA_budget_wetdepos_hour 
  
  ! Budgets on daily basis (also landuse dependent available)
  ! Calculate only tracers/landuses which are necessary for the output
  integer, parameter                      ::  max_budget_fields = 5
  integer                                 ::  nbudget_fields              ! number of different output types for budget-label
  character(len=32)                       ::  budget_fields(max_budget_fields) ! list of output type for budget-label 
  
  integer                                 ::  ndrydep_label_output_tracer ! number of tracers in dry deposition output
  integer, pointer                        ::  labelled_specs_drydep_output(:)
  character(len=16), pointer              ::  labelled_specs_drydep_output_names(:)
  integer, pointer                        ::  labelled_specs_drydep_output_numbers(:)

  integer                                 ::  nwetdep_label_output_tracer ! number of tracers in wet deposition output
  integer, pointer                        ::  labelled_specs_wetdep_output(:)
  character(len=16), pointer              ::  labelled_specs_wetdep_output_names(:)
  integer, pointer                        ::  labelled_specs_wetdep_output_numbers(:)

  integer                                 ::  nbudget_label_output_tracer ! number of tracers in budget output (daily)
  integer, pointer                        ::  labelled_specs_budget_output(:)
  character(len=16), pointer              ::  labelled_specs_budget_output_names(:)
  integer, pointer                        ::  labelled_specs_budget_output_numbers(:)
  
  integer                                 ::  nbudget_label_output_lu     ! number of landuses in budget output (daily)
  integer, pointer                        ::  landuse_budget_output(:)
  character(len=16), pointer              ::  landuse_budget_output_names(:)
  integer, pointer                        ::  landuse_budget_output_numbers(:)
  
  real, allocatable, dimension(:,:,:,:)   ::  SA_budget_drydepos_day    ! (nx,ny,nbudget_label_output_tracer,nlabel)
  real, allocatable, dimension(:,:,:,:,:) ::  SA_budget_drydepos_lu_day ! (nx,ny,nbudget_label_output_tracer,nbudget_label_output_lu,nlabel)
  real, allocatable, dimension(:,:,:,:)   ::  SA_budget_wetdepos_day    ! (nx,ny,nbudget_label_output_tracer,nlabel)
    
  ! compensation point
  real, allocatable, dimension(:,:,:,:)   ::  SA_frac_comp_point 
  real, allocatable, dimension(:,:,:,:)   ::  SA_conc_comp_point
  integer                                 ::  iexdry_emis, iexdry_depo
  
  ! Chemistry
  integer, allocatable, dimension(:)      ::  array_cols, array_rows
  real, allocatable, dimension(:)         ::  array_prod_total
  real, allocatable, dimension(:)         ::  array_prod_total_backup
  real                                    ::  dt_total, dt_total_backup
  real, allocatable, dimension(:)         ::  Labelled_Chem_Specs
  integer                                 ::  SA_n_nonzeros
  logical                                 ::  First_chem_step

  ! Boundary conditions
  integer                                 ::  cache_nlon, cache_nlat, cache_nlev
  logical                                 ::  cache_setup
  integer, allocatable, dimension(:,:,:)  ::  cache_ij_west, cache_ij_south
  integer, allocatable, dimension(:,:,:)  ::  cache_ij_east, cache_ij_north

  ! properties of boundary condition file:
  character(len=512)                      ::  arch_path, arch_model, arch_expid, arch_name

contains

!
!----------------Definition of the labels----------------------------------------------
!

  subroutine SA_Label_Definition(icat,icountry, emis_set_name,def_label,  status, country_name)

    use GO, only : GoMatchValue
    
    ! --- in/out ------------------
    integer, intent(in)                     ::  icat, icountry
    character(len=*), intent(in)            ::  emis_set_name
    integer, intent(out)                    ::  def_label
    integer, intent(out)                    ::  status
    character(len=*), intent(in), optional  :: country_name
    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Label_definition'

    ! --- local ------------------------------
    integer               ::  SA_ilabel
    integer               ::  iset
    ! --- begin -------------------

    if (First_def_step) then
      ! define the names of the labels, for the variable label name in the nc output file
      ! called in Initialization part of the labeling
      SA_Label_Names(1) = 'NLD'
      SA_Label_Names(2) = 'Abroad'
      !SA_Label_Names(..)      = '..'
      !SA_Label_Names(..)      = '..'

     ! define short label names, at most 10 characters such that for the grads/ctl the names are distinct
      SA_short_label_names(1) = 'NLD'
      SA_short_label_names(2)  = 'Abr'
      !SA_short_label_names(..)  = '..'
      !SA_short_label_names(..)  = '..'
      
      ! Maximum label length (need for dimension in nc-file)
      labelname_len = 1
      do SA_ilabel = 1, SA_nlabel
        labelname_len = max( labelname_len, len(trim(SA_Label_Names(SA_ilabel)) ) )
      end do
      
      ! labels defined, dummy label number = -1
      First_def_step = .false.
      def_label = -1
      
    else 
      
      ! match name of emission set
      call goMatchValue( trim(emis_set_name), emis_set_names, iset, status )
      IF_NOTOK_RETURN(status=1)      
      
      ! --------------------------------------------------------
      ! Important NOTE, if countries are labelled and 
      !  agricultural or resuspension dust is in simulations
      ! --------------------------------------------------------
      !
      ! Some countries defined in general country map, are not in emission input 
      ! like Greenland, Greenland Sea, Gibraltar, etc..
      ! Resuspension or agricultural dust from those areas are attached to countries from general map
      ! and therefore do not have a valid country number in the emission set
      ! those situations are tracked by country number -999, and optional argument country_name is filled
      !
      ! ---> Dust emissions from those areas should be attached to one of the defined labels
      ! For example Gibraltar --> Spain
      ! 
      ! if ( icountry == -999 ) then
      !   if trim(country_name) == 'GIB' 
      !     def_label =  <label_number_spain>
      !   elif ....
      ! 
      ! end if
       
      
      ! Define here the labels       
                
      ! default it can be a function of countries and categories or make combinations
      if ( icountry == -999 ) then
        ! undefined countries for dust emissions
        def_label = 2
      else if ( Emis_countries(iset)%emis_country_names(icountry) == 'NLD' ) then
        def_label = 1
      else 
        def_label = 2
      end if
      
      !if ( Emis_Sectors(iset)%emis_sector_names(icat) == 'Industry' ) then
      !  def_label = 1
      !else 
      !  def_label = 2
      !end if

    end if

    ! ok
    status = 0

  end subroutine SA_Label_Definition
  
!
! --------------- Initialization of the Labeling routines
!

  subroutine SA_Labeling_Init( rcF, status)

    use GO                , only :  AnyDate
    use GO                , only :  goMatchValues
    use GO                , only :  goVarValue
    use GO                , only :  TrcFile, ReadRc
    use Indices           , only :  specname

    use LE_Data           , only :  LE_Data_Enable

    ! --- in/out ---------------------------

    type(TRcFile), intent(in)           ::  rcF
    integer,intent(out)                 ::  status

    ! --- const ----------------------------
    character(len=*), parameter :: rname = mname//'/SA_Labeling_Init'

    ! --- local ----------------------------

    integer, pointer                ::  labelled_specs_numbers(:)
    character(len=512)              ::  labelled_specs_list
    integer                         ::  no_label, no_extra_labels
    integer                         ::  ispec
    integer                         ::  dummy_label
    character(len=5)                ::  dummy_set_name
    character(len=512)              ::  key

    ! --- begin ----------------------------

    ! label definition
    call ReadRc( rcF, 'labels.nlabel', no_label, status )
    IF_NOTOK_RETURN(status=1)

    ! number of extra labels (usually 4: natural, bound, init, and aloft)
    no_extra_labels = 4

    SA_nlabel = no_label + no_extra_labels

    Label_Natural = no_label + 1
    Label_Bound = no_label + 2
    Label_Init = no_label + 3
    Label_Aloft = no_label + 4

    ! Which species should be labelled ?
    call ReadRc( rcF, 'labels.labelled.specs', labelled_specs_list, status )
    IF_NOTOK_RETURN(status=1)

    ! allocate
    allocate( labelled_specs_names  (nspec), stat = status)
    IF_NOTOK_RETURN(status=1)
    allocate( labelled_specs_numbers(nspec), stat = status)
    IF_NOTOK_RETURN(status=1)
    allocate (Labelled_specs        (nspec), stat = status)
    IF_NOTOK_RETURN(status=1)

    ! match labelled species within list of all species
    call goMatchValues(labelled_specs_list, specname,  &
                       nspec_labelled, labelled_specs_names, labelled_specs_numbers, &
                       status)
    IF_NOTOK_RETURN(status=1)

    Labelled_specs = -1

    do ispec = 1, nspec_labelled
      Labelled_specs(labelled_specs_numbers(ispec)) = ispec
    enddo

    ! Allocate and initialise LE_concentration in this module
    allocate( LE_conc(nx,ny,nz,nspec_labelled), stat = status)
    IF_NOTOK_RETURN(status=1)

    LE_conc = 0.0

    ! Allocate and initialize Fractions Array with Source apportionments
    allocate( SA_frac(nx,ny,nz,nspec_labelled,SA_nlabel), stat = status)
    IF_NOTOK_RETURN(status=1)

    SA_frac = 0.0
    
    allocate( SA_frac_comp_point(nx,ny,nspec_labelled,SA_nlabel), stat=status) 
    IF_NOTOK_RETURN(status=1)
    SA_frac_comp_point = 0.0
    SA_frac_comp_point(:,:,:,Label_Natural) = 1.0

    allocate( SA_conc_comp_point(nx,ny,nspec_labelled,SA_nlabel), stat=status) 
    IF_NOTOK_RETURN(status=1)
    SA_conc_comp_point = 0.0
    
    allocate( SA_Label_Names(SA_nlabel), stat = status)
    IF_NOTOK_RETURN(status=1)

    allocate( SA_Short_Label_Names(SA_nlabel), stat = status)
    IF_NOTOK_RETURN(status=1)

    SA_Label_Names(Label_Natural) = 'Natural'
    SA_Label_Names(Label_Bound)   = 'Boundary'
    SA_Label_Names(Label_Init)    = 'Init'
    SA_Label_Names(Label_Aloft)   = 'Aloft'

    SA_Short_Label_Names(Label_Natural) = 'Nat'
    SA_Short_Label_Names(Label_Bound)   = 'Bound'
    SA_Short_Label_Names(Label_Init)    = 'Init'
    SA_Short_Label_Names(Label_Aloft)   = 'Aloft'
  
    ! No emission sets yet, fill with empty space
    emis_set_names = ''
    nemis_set = 0
    
    ! Get names from Label definition
    First_def_Step = .true.
    dummy_set_name = 'dummy'
    call SA_Label_definition( 1,1,trim(dummy_set_name),dummy_label,status )
    IF_NOTOK_RETURN(status=1)
    
    ! enable data:
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)

    ! Initialize variable array for the different processes
    call SA_Emis_Init(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Advec_Init(status)
    IF_NOTOK_RETURN(status=1)
    
    call SA_Hdiff_Init(status)
    IF_NOTOK_RETURN(status=1)
    
    call SA_Adjust_Init(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Vdiff_Init(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Sedim_Init(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Depos_Init(rcF, status)
    IF_NOTOK_RETURN(status=1)
    
    !
    ! Read Boundary conditions from a boundary run?
    !
    ! Read location and properties of boundary files
    call ReadRc( rcF, 'le.bound.le.path', arch_path, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'le.bound.le.key.labeling', key, status )
    IF_NOTOK_RETURN(status=1)

    arch_model = 'LE'
    call goVarValue( key, ';', 'model', '=', arch_model, status )
    IF_ERROR_RETURN(status=1)

    arch_expid = 'base'
    call goVarValue( key, ';', 'expid', '=', arch_expid, status )
    IF_ERROR_RETURN(status=1)
   
    arch_name = 'label-3d'
    call goVarValue( key, ';', 'name', '=', arch_name, status )
    IF_ERROR_RETURN(status=1)

    ! boundary cells not yet defined
    cache_setup = .false.

    deallocate( labelled_specs_numbers, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Labeling_Init

! --------------------------------------------------------------------------------------

  subroutine SA_Labeling_Done( status)

    ! --- in/out -------------------

    integer, intent(out)                    ::  status

    ! --- const --------------------
    character(len=*), parameter   ::  rname = mname//'/SA_Labeling_Done'

    ! --- local --------------------

    ! --- begin --------------------

    ! close help arrays for different processes
    call SA_Emis_Done(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Advec_Done(status)
    IF_NOTOK_RETURN(status=1)
    
    call SA_HDiff_Done(status)
    IF_NOTOK_RETURN(status=1)
    
    call SA_Adjust_Done(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Vdiff_Done(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Sedim_Done(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Depos_Done(status)
    IF_NOTOK_RETURN(status=1)

    call SA_Chem_Done(status)
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Labeling_Done

! --------------------------------------------------------------------------------------

  subroutine SA_Frac_Init(c,status)
    ! State the fractions for the initial condition

    ! --- in/out -------------------
    real, intent(in)        ::  c(nx,ny,nz,nspec)
    integer, intent(out)    ::  status

    ! --- local --------------------
    integer                 ::  ispec, SA_ispec

    ! --- begin --------------------

    ! set the Lotos-Euros concentration in this module (only for the labelled species)
    do ispec = 1, nspec
      SA_ispec = Labelled_specs(ispec)
      if (SA_ispec > 0) then
        LE_conc(:,:,:,SA_ispec) = c(:,:,:,ispec)
      endif
    enddo

    ! positive initial concentrations get initial label 1
    where( LE_conc >= 0.0 )
      SA_frac(:,:,:,:,Label_Init) = 1.0
    endwhere

    ! ok
    status = 0

  end subroutine SA_Frac_Init

!
! ----------------Synchronization-------------------------------------------------------
!

  subroutine SA_Synchronize_Conc(c,status)
    
    ! --- in/out --------------------
    real, intent(in)      ::  c(nx,ny,nz,nspec)
    integer, intent(out)  ::  status

    ! --- const ---------------------
    character(len=*), parameter :: rname = mname//'/SA_Synchronize_Conc'

    ! --- local ---------------------
    integer               ::  ispec, ix,iy,iz
    integer               ::  SA_ispec

    ! --- begin ---------------------
    
    !$OMP parallel &
    !$OMP default( none ) &
    !$OMP shared( Labelled_specs ) &
    !$OMP shared( LE_conc, c) &
    !$OMP private( ispec, SA_ispec )
    !$OMP do SCHEDULE( DYNAMIC )
    do ispec = 1, nspec      
      SA_ispec = Labelled_specs(ispec)

      if (SA_ispec > 0 ) then
        do iz = 1, nz
        do iy = 1, ny
        do ix = 1, nx
          if ( c(ix,iy,iz,ispec) < -1e-4  ) then
            write( gol, '("Found negative concentration: ",e12.5)' ) c(ix,iy,iz,ispec) ; call goErr
            write( gol, '("ix,iy,iz,spec: ", i0," ",i0," ",i0, " ",a)' ) ix,iy,iz,trim(specname(ispec)); call goErr
            TRACEBACK;status=1;return
          else if ( c(ix,iy,iz,ispec) < 0.0 .and. c(ix,iy,iz,ispec) > -1.e-4 ) then
            LE_conc(ix,iy,iz,SA_ispec) = 0.0
          else
            LE_conc(ix,iy,iz,SA_ispec) = c(ix,iy,iz,ispec)
          end if
        end do
        end do
        end do
      endif

    enddo
    !$OMP end do
    !$OMP end parallel

    ! ok
    status = 0

  end subroutine SA_Synchronize_Conc

!
! ----------------EMISSIONS-------------------------------------------------------------
!
  subroutine SA_Emis_Init(status)

    ! --- in/out -----------
    integer, intent(out) ::  status

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Init'

    ! --- begin ------------------------------

    ! allocate and initialize
    allocate( SA_emis(nx,ny,nz,nspec_labelled,SA_nlabel), stat=status )
    IF_NOTOK_RETURN(status=1)

    SA_emis = 0.0

    ! ok
    status = 0

  end subroutine SA_Emis_Init

! --------------------------------------------------------------------------------------

  subroutine SA_Emis_Done(status)

    ! --- in/out ------------------------
    integer, intent(out) ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Done'

    ! --- begin -------------------------

    ! deallocate arrays at end of modelrun
    deallocate(SA_emis, stat=status)
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Emis_Done

! --------------------------------------------------------------------------------------

  subroutine SA_Emis_Sectors( emis_set_name, em_sectors, nsector, em_countries, ncountry, status )
    
    use GO, only : GoMatchValue
    
    ! --- in/out ------------------------
    
    character(len=*), intent(in)    ::  emis_set_name  ! name of emission input set
    character(len=*), intent(in)    ::  em_sectors(nsector)
    integer, intent(in)             ::  nsector
    character(len=*), intent(in)    ::  em_countries(ncountry)
    integer, intent(in)             ::  ncountry
    integer, intent(out)            ::  status

    ! --- local -------------------------
    integer                         ::  isector, icountry, iisector
    integer                         ::  iem_set

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Sectors'

    ! --- begin -------------------------
    
    ! Find matching emission set ( check for double naming )
    call GoMatchValue( trim( emis_set_name), emis_set_names, iem_set, status, quiet=.true. )
    ! not found:
    if ( status < 0 ) then
      ! increase numbers
      nemis_set = nemis_set + 1
      iem_set = nemis_set
      emis_set_names(iem_set) = trim(emis_set_name)
      ! Allocate sector/country array for this set
      allocate( Emis_Sectors(iem_set)%emis_sector_names(nsector), stat=status)
      IF_NOTOK_RETURN(status=1)
      allocate( Emis_Countries(iem_set)%emis_country_names(ncountry), stat=status)
      IF_NOTOK_RETURN(status=1)
    else if ( status > 0 ) then
      write (gol, '("Unknown Error found, when try to match emission set: ",a)' ) trim(emis_set_name) ; call goErr
      TRACEBACK;status=1;return
    end if        
    
    ! Fill sector names for this set
    do isector = 1, nsector
      Emis_Sectors(iem_set)%emis_sector_names(isector) = trim(em_sectors(isector))
    end do
    ! Fill country names for this set
    do icountry = 1, ncountry
      Emis_Countries(iem_set)%emis_country_names(icountry) = trim(em_countries(icountry))
    end do
    
    ! Define sector number for dust-resuspension and dust-agriculture
    do isector =  1, nsector
      ! sector number for road transport needed to state dust resuspension
      if (em_sectors(isector) == 'Road_transport' .or. &
          em_sectors(isector) == 'Road transport' .or. &
          em_sectors(isector) == 'Road transport, tyre, brake and road wear' .or. &
          em_sectors(isector) == 'RoadTransport non-exhaust' .or. &
          em_sectors(isector) == 'RoadTransport non-exhaust-' .or. &
          em_sectors(isector) == 'Road transport, wear, LDV' ) then
        Emis_Sectors(iem_set)%isector_resuspension = isector
        exit
      endif
      if (isector == nsector) then
        write (gol, '("source: traffic, not found, dust resuspension cannot be stated in a specific sector")' ) ; call goErr
        write (gol, '("Available sources : ")' ) ; call goErr
        do iisector = 1, nsector 
          write( gol, '("Sector ",i0,": ",a)') iisector, trim(em_sectors(iisector)) ; call goErr
        end do
        TRACEBACK; status=1; return
      endif
    end do

    do isector =  1, nsector
      ! sector number for agriculture needed to state dust agriculture
      if (em_sectors(isector) == 'Agriculture' .or. &
          em_sectors(isector) == 'Other_agriculture' .or. &
          em_sectors(isector) == 'Other agriculture' .or. &
          em_sectors(isector) == 'Agricultural Other-' .or. &
          em_sectors(isector) == 'Agricultural Other' ) then
        Emis_Sectors(iem_set)%isector_agriculture = isector
        exit
      endif
      if (isector == nsector) then
        write (gol, '("source: agriculture, not found, dust agriculture cannot be stated in a specific sector")' ) ; call goErr
        write (gol, '("Available sources : ")' ) ; call goErr
        do iisector = 1, nsector 
          write( gol, '("Sector ",i0,": ",a)') iisector, trim(em_sectors(iisector)) ; call goErr
        end do
        TRACEBACK; status=1; return
      endif
    end do

    ! ok
    status = 0

  end subroutine SA_Emis_Sectors

! --------------------------------------------------------------------------------------

  subroutine SA_Emis_Sectors_Done(emis_set_name, status)
    
    use GO, only : GoMatchValue

    ! --- in/out -----------------
    character(len=*), intent(in)  ::  emis_set_name
    integer, intent(out)          ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Sectors_Done'

    ! --- local ------------------
    integer                     ::  iem_set

    ! --- begin ------------------

    ! Find matching emission set name
    call GoMatchValue( trim(emis_set_name), emis_set_names, iem_set, status )
    IF_NOTOK_RETURN(status=1)

    ! clear naming:
    emis_set_names(iem_set) = ''
    
    ! Clear emission sectors
    deallocate(Emis_Sectors(iem_set)%emis_sector_names, stat=status)
    IF_NOTOK_RETURN(status=1)
    ! Clear emission countries
    deallocate(Emis_Countries(iem_set)%emis_country_names, stat=status)
    IF_NOTOK_RETURN(status=1)
   
    ! ok
    status = 0

  end subroutine SA_Emis_Sectors_Done

! --------------------------------------------------------------------------------------

  subroutine SA_Emis_Setup_TNO(ix,iy,ispec,icat,icountry,delta_emis_a, emis_set_name, status)

    ! Collect the delta_emissions in each time step, on each location for each country and categorie.
    ! This is called in le_emis_type_base.f90

    ! --- in/out -----------

    integer, intent(in)           ::  ix,iy
    integer, intent(in)           ::  ispec, icat, icountry
    real, intent(in)              ::  delta_emis_a(1:nz)
    character(len=*), intent(in)  ::  emis_set_name           ! name of emission set
    integer, intent(out)          ::  status

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Setup_TNO'

    ! --- local ------------------------------
    integer              ::  SA_ispec
    integer              ::  SA_ispec_emis

    ! --- begin ------------------------------


    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    ! Component labelled ?
    if (SA_ispec > 0) then
      
      ! Define the label for this country and category
      call SA_Label_Definition(icat,icountry, emis_set_name, def_label, status)
      IF_NOTOK_RETURN(status=1)

      ! Add the delta emissions to the right label
      SA_emis(ix,iy,:,SA_ispec,def_label) = SA_emis(ix,iy,:,SA_ispec,def_label) + delta_emis_a

    end if  ! spec labelled?

    ! ok
    status = 0

  end subroutine SA_Emis_Setup_TNO

! --------------------------------------------------------------------------------------

  subroutine SA_Emis_Setup_EDGAR(ix,iy,ispec,icat,izone,delta_emis_a, emis_set_name, status)

    ! Collect the delta_emissions in each time step, on each location for each country and categorie.
    ! This is called in le_emis_type_base.f90

    ! --- in/out -----------
    integer, intent(in)           ::  ix,iy
    integer, intent(in)           ::  ispec, icat, izone
    real, intent(in)              ::  delta_emis_a(1:nz)
    character(len=*), intent(in)  ::  emis_set_name
    integer, intent(out)          ::  status

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Setup_EDGAR'

    ! --- local ------------------------------
    integer              ::  SA_ispec

    ! --- begin ------------------------------


    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    ! Component labelled ?
    if (SA_ispec > 0) then
      ! Define the label for this timezone and category
      call SA_Label_Definition(icat,izone, trim(emis_set_name), def_label, status)
      IF_NOTOK_RETURN(status=1)

      ! Add the delta emissions to the right label
      SA_emis(ix,iy,:,SA_ispec,def_label) = SA_emis(ix,iy,:,SA_ispec,def_label) + delta_emis_a
    endif

    ! ok
    status = 0

  end subroutine SA_Emis_Setup_EDGAR

!--------------------------------------------------------------------------------------

  subroutine SA_Emis_Setup_Natural(ix,iy,ispec,delta_emis, status)

    ! --- in/out ---------
    integer, intent(in)   ::  ix,iy,ispec
    real,intent(in)       ::  delta_emis
    integer, intent(out)  ::  status

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Setup_natural'

    ! --- local ------------------------------
    integer               ::  SA_ispec

    ! --- begin ------------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    ! Natural label
    def_label = Label_Natural

    ! Component labelled ?
    if (SA_ispec > 0) then
      ! Add the delta emissions to the right label
      SA_emis(ix,iy,1,SA_ispec,def_label) = SA_emis(ix,iy,1,SA_ispec,def_label) + delta_emis
    endif

    ! ok
    status = 0

  end subroutine SA_Emis_Setup_natural

! ---------------------------------------------------------------------------------------

  subroutine SA_Emis_Setup_Dust(ix,iy,ispec,delta_emis, source, status)

    use LE_country, only : country_map_n, country_map_frac, country_map_code
    use GO        , only : GoMatchValue
    ! --- in/out ---------
    integer, intent(in)   ::  ix,iy,ispec
    real,intent(in)       ::  delta_emis
    character(len=*)      ::  source
    integer, intent(out)  ::  status

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Setup_Dust'

    ! --- local ------------------------------
    integer               ::  SA_ispec
    integer               ::  icountry, iicountry
    real                  ::  fraction
    integer               ::  icat
    character(len=32)     ::  country
    integer               ::  iem_set

    ! --- begin ------------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    ! Component labelled ?
    if (SA_ispec > 0) then
      
      ! Find first index filled with emissions (to determine sectornumbers for agricultural and resuspension dust)
      do iem_set = 1, nemis_set
        if ( .not. (emis_set_names(iem_set) == '' ) ) then
          exit
        endif
      end do
      ! Dust from wind, natural emission
      if ( source == 'wind' ) then
        def_label = Label_Natural
        SA_emis(ix,iy,1,SA_ispec,def_label) = SA_emis(ix,iy,1,SA_ispec,def_label) + delta_emis

      ! Dust from resuspension or agriculture
      else if ( trim(source) == 'resuspension' .or. trim(source) == 'agriculture' ) then        
        ! Pick corresponding category from first available set
        if ( trim(source) == 'resuspension'  ) then
          icat = Emis_Sectors(iem_set)%isector_resuspension
        else if ( trim(source) == 'agriculture' ) then
          icat = Emis_Sectors(iem_set)%isector_agriculture
        else
          write(gol, '("Incorrect dust source sector found: ",a)' ) trim( source ) ; call GoErr
          TRACEBACK;status=1;return
        end if

        ! Check which country from general country map
        do icountry = 1, country_map_n
          ! fraction of country in this cell
          fraction = country_map_frac(ix,iy,icountry)
          
          ! country in this grid cell?
          if (fraction > 0.0 ) then
            ! pick country name
            country = country_map_code(icountry)
            ! match with country names in first available emission set
            call GoMatchValue( trim(country), Emis_Countries(iem_set)%emis_country_names, iicountry, status, quiet=.true. )
            ! country not found, because country not available in emission set 1: <-- for example Gibraltar not in emission set
            if ( status < 0 ) then
              iicountry = -999
            else if ( status > 0 ) then
              write( gol, '("Unknown error found for matching country")' ) ; call goErr
              TRACEBACK;status=1;return
            endif
            ! Get label number ( icat, icountry do now match with 1st available emission set )
            call SA_Label_Definition(icat,iicountry,emis_set_names(iem_set), def_label,status,country_name=trim(country))
            IF_NOTOK_RETURN(status=1)

            ! Add the delta emissions to the correct label
            SA_emis(ix,iy,1,SA_ispec,def_label) = SA_emis(ix,iy,1,SA_ispec,def_label) + fraction * delta_emis
          end if

        end do  ! countries

      else
        write (gol, '("incorrect dust source, must be wind, resuspension of agriculture")' ) ; call goErr
        TRACEBACK; status=1; return
      endif

    endif

    ! ok
    status = 0

  end subroutine SA_Emis_Setup_Dust

! ---------------------------------------------------------------------------------------

  subroutine SA_Emis_Setup_Fire(ispec,delta_emis, status)

    ! --- in/out ---------
    integer, intent(in)   ::  ispec
    real,intent(in)       ::  delta_emis(nx,ny,nz)
    integer, intent(out)  ::  status

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Setup_Fire'

    ! --- local ------------------------------
    integer               ::  SA_ispec

    ! --- begin ------------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    ! Natural label for fire emissions
    def_label = Label_Natural

    ! Component labelled ?
    if (SA_ispec > 0) then
      ! Add the delta emissions to the right label
      SA_emis(:,:,:,SA_ispec,def_label) = SA_emis(:,:,:,SA_ispec,def_label) + delta_emis
    endif

    ! ok
    status = 0

  end subroutine SA_Emis_Setup_Fire

! ---------------------------------------------------------------------------------------

  subroutine SA_Emis_Setup_Soil_NOx(emis_level,ispec,delta_emis, status)

    ! --- in/out ---------
    integer, intent(in)   ::  emis_level, ispec
    real,intent(in)       ::  delta_emis(nx,ny)
    integer, intent(out)  ::  status

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Setup_Soil_NOx'

    ! --- local ------------------------------
    integer               ::  SA_ispec

    ! --- begin ------------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    ! Natural label for soil nox emissions
    def_label = Label_Natural

    ! Component labelled ?
    if (SA_ispec > 0) then
      ! Add the delta emissions to the right label
      SA_emis(:,:,emis_level,SA_ispec,def_label) = SA_emis(:,:,emis_level,SA_ispec,def_label) + delta_emis
    endif

    ! ok
    status = 0

  end subroutine SA_Emis_Setup_Soil_NOx

! ---------------------------------------------------------------------------------------

  subroutine SA_Emis_Delta_c(ispec,factor, status)

    ! For each specie, the emissions are converted to a concentration

    ! -- in/out ----
    real, intent(in)          ::  factor(nx,ny,nz)    ! factor factor is stated in emis routine ppb_fac or ugm3_fac
    integer, intent(in)       ::  ispec     ! species
    integer, intent(out)      ::  status

    ! --- local ----
    integer                   ::  ix, iy, iz, SA_ispec
    real                      ::  denum, inv_denum
    real                      ::  SA_flux_emis(SA_nlabel)

    ! --- const ------------------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_Delta_c'

    ! --- begin ------------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    ! Component labelled ?
    if (SA_ispec > 0)  then

      ! Calculate new fractions after addition of the emissions
      do iz = 1,nz
      do iy = 1,ny
      do ix = 1,nx
        ! flux
        SA_flux_emis = factor(ix,iy,iz) * SA_emis(ix,iy,iz,SA_ispec,:)

        ! calculate denumerator only once
        denum = sum(SA_flux_emis) + LE_conc(ix,iy,iz,SA_ispec)

        if ( denum > numlim_conc ) then
          inv_denum = 1.0 / denum

          SA_frac(ix,iy,iz,SA_ispec,:) = (SA_frac(ix,iy,iz,SA_ispec,:) * LE_conc(ix,iy,iz,SA_ispec) + SA_flux_emis(:) ) * inv_denum
        else
          SA_frac(ix,iy,iz,SA_ispec,:) = 0.0
        endif

      enddo
      enddo
      enddo

    endif

    ! ok
    status = 0

  end subroutine SA_Emis_Delta_c

! -------------------------------------------------------------------------------------

  subroutine SA_Emis_Reset(status)

    ! --- in/out -----------------
    integer, intent(out) ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Emis_reset'

    ! --- begin ------------------

    ! reset emissions to zero otherwise emissions of the next time step includes also emissions from this timestep
    SA_emis = 0.0

    ! ok
    status = 0

  end subroutine SA_Emis_Reset

!
! ------------------------ADVECTION-----------------------------------------------------
!

  subroutine SA_Advec_Init(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Advec_Init'

    ! --- begin -------------------------

    ! Allocate help variable for the fractions during the advection process
    allocate( SA_frac_adv(SA_nlabel,0:nx+1,0:ny+1,0:nz+1,nspec_labelled), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! initialize Fractions of the boundary conditions. Zero inside domain and one at the boundary
    ! These default labels can be override by the boundary setup if a run is used with a larger domain
    SA_frac_adv = 0.0

    ! x-boundary
    SA_frac_adv(Label_Bound,0   ,1:ny,1:nz,:) = 1
    SA_frac_adv(Label_Bound,nx+1,1:ny,1:nz,:) = 1

    ! y-boundary
    SA_frac_adv(Label_Bound,1:nx,0   ,1:nz,:) = 1
    SA_frac_adv(Label_Bound,1:nx,ny+1,1:nz,:) = 1

    ! z-boundary (aloft)
    SA_frac_adv(Label_Bound,1:nx,1:ny,nz+1,:) = 1

    ! ok
    status = 0

  end subroutine SA_Advec_Init

! ---------------------------------------------------------------------------------------

  subroutine SA_Advec_Done(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Advec_Done'

    ! --- begin -------------------------

    deallocate(SA_frac_adv, stat=status)
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Advec_Done

! ---------------------------------------------------------------------------------------

  subroutine SA_Advec_Setup(k,status)

    ! --- in/out ------------------------
    
    integer, intent(in)       ::  k
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Advec_Setup'
    
    ! --- local ------------------
    
    integer                   ::  SA_ispec, SA_ilabel
    
    ! --- begin -------------------------
    
    SA_ispec = Labelled_specs(k) 
    if ( SA_ispec > 0 ) then
      !set fractions inside domain in the help variable for advection process
      do SA_ilabel = 1, SA_nlabel
        SA_frac_adv(SA_ilabel,1:nx,1:ny,1:nz,SA_ispec) = SA_frac(1:nx,1:ny,1:nz,SA_ispec,SA_ilabel)
      end do
    end if

    ! ok
    status = 0

  end subroutine SA_Advec_setup

! ---------------------------------------------------------------------------------------

  subroutine SA_Advec_x(ch, fluxx, volume_loc, ispec, status)

    use dims, only :  nx, ny, nz

    ! --- in/out ------------------------

    real, intent(in), pointer ::  ch(:,:,:)         ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in), pointer ::  fluxx(:,:,:)      ! (0:nx,1:ny,2:nz)
    real, intent(in), pointer ::  volume_loc(:,:,:) ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    integer, intent(in)       ::  ispec
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Advec_x'

    ! --- local ------------------

    integer                   ::  ix,iy,iz
    integer                   ::  SA_ispec, SA_ilabel
    real                      ::  denum, inv_denum
    real, allocatable         ::  SA_flux_adv(:,:,:,:)
    real, allocatable         ::  LE_conc_abs(:,:,:)

    ! --- begin -------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    if (SA_ispec > 0) then

      ! storage:
      allocate( SA_flux_adv(SA_nlabel,nx,ny,nz), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( LE_conc_abs(1:nx,1:ny,1:nz), stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! init:
      SA_flux_adv = 0.0
      ! state mass as help in this routine
      LE_conc_abs = ch(1:nx,1:ny,1:nz)*volume_loc(1:nx,1:ny,1:nz)

      ! loop over cells to make fluxes in each cell
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx
        ! calculate mass fluxes per cell per label (fluxes through cell edges are taken as fraction from donor cell)
        ! negative winds
        if      (fluxx(ix-1,iy,iz) <= 0 .and. fluxx(ix,iy,iz) <= 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxx(ix-1,iy,iz) * SA_frac_adv(:,ix  ,iy,iz,SA_ispec) - &
                                    fluxx(ix  ,iy,iz) * SA_frac_adv(:,ix+1,iy,iz,SA_ispec)
        ! inflow cells
        else if (fluxx(ix-1,iy,iz) > 0 .and. fluxx(ix,iy,iz) <= 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxx(ix-1,iy,iz) * SA_frac_adv(:,ix-1,iy,iz,SA_ispec) - &
                                    fluxx(ix  ,iy,iz) * SA_frac_adv(:,ix+1,iy,iz,SA_ispec)
        ! outflow cells
        else if (fluxx(ix-1,iy,iz) <= 0 .and. fluxx(ix,iy,iz) > 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxx(ix-1,iy,iz) * SA_frac_adv(:,ix  ,iy,iz,SA_ispec) - &
                                    fluxx(ix  ,iy,iz) * SA_frac_adv(:,ix  ,iy,iz,SA_ispec)
        ! positive winds
        else if (fluxx(ix-1,iy,iz) > 0 .and. fluxx(ix,iy,iz) > 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxx(ix-1,iy,iz) * SA_frac_adv(:,ix-1,iy,iz,SA_ispec) - &
                                    fluxx(ix  ,iy,iz) * SA_frac_adv(:,ix  ,iy,iz,SA_ispec)
        endif

        ! calculate new fractions after x-advection
        denum = sum(SA_flux_adv(:,ix,iy,iz)) + LE_conc_abs(ix,iy,iz)
        if ( denum > numlim_conc ) then
          inv_denum = 1.0 / denum
          SA_frac(ix,iy,iz,SA_ispec,:) = &
              (SA_frac_adv(:,ix,iy,iz,SA_ispec) * LE_conc_abs(ix,iy,iz) + SA_flux_adv(:,ix,iy,iz) ) * inv_denum

        else
          SA_frac(ix,iy,iz,SA_ispec,:) = 0.0
        endif

      end do
      end do
      end do
      
      ! Set the calculated fractions in the help array for the advection routine
      ! check for negative fractions (numerical errors)
      ! Some compilers have trouble with where loops in OPEN-MP code
      ! changed to nx,ny,nz loop
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx
        do SA_ilabel = 1, SA_nlabel
          if (SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) < 0.0 ) then
            SA_frac_adv(SA_ilabel,ix,iy,iz,SA_ispec) = 0.0
            SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) = 0.0
          else 
            SA_frac_adv(SA_ilabel,ix,iy,iz,SA_ispec) = SA_frac(ix,iy,iz,SA_ispec,SA_ilabel)
          end if
        end do
      end do
      end do
      end do
      
      ! clear:
      deallocate( SA_flux_adv, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( LE_conc_abs, stat=status )
      IF_NOTOK_RETURN(status=1)

    endif

    !ok
    status = 0

  end subroutine SA_Advec_x

! ---------------------------------------------------------------------------------------

  subroutine SA_Advec_y(ch, fluxy, volume_loc, ispec , status)

    use dims, only : nx, ny, nz

    ! --- in/out ------------------------
    real, intent(in), pointer ::  ch(:,:,:)         ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in), pointer ::  fluxy(:,:,:)      ! (1:nx,0:ny,1:nz)
    real, intent(in), pointer ::  volume_loc(:,:,:) ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    integer, intent(in)       ::  ispec
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Advec_y'

    ! --- local ------------------
    integer                   ::  ix,iy,iz
    integer                   ::  SA_ispec, SA_ilabel
    real                      ::  denum, inv_denum
    real, allocatable         ::  SA_flux_adv(:,:,:,:)
    real, allocatable         ::  LE_conc_abs(:,:,:)

    ! --- begin -------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    if (SA_ispec > 0) then

      allocate( SA_flux_adv(SA_nlabel,nx,ny,nz), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( LE_conc_abs(1:nx,1:ny,1:nz), stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! init:
      SA_flux_adv = 0.0
      ! state mass as help in this routine
      LE_conc_abs = ch(1:nx,1:ny,1:nz)*volume_loc(1:nx,1:ny,1:nz)

      ! loop over cells to make fluxes in each cell
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx
        ! calculate mass fluxes per cell per label (fluxes through cell edges are taken as fraction from donor cell)
        ! negative winds
        if      (fluxy(ix,iy-1,iz) <= 0 .and. fluxy(ix,iy,iz) <= 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxy(ix,iy-1,iz) * SA_frac_adv(:,ix,iy  ,iz,SA_ispec) - &
                                    fluxy(ix,iy  ,iz) * SA_frac_adv(:,ix,iy+1,iz,SA_ispec)
        ! inflow cells
        else if (fluxy(ix,iy-1,iz) > 0 .and. fluxy(ix,iy,iz) <= 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxy(ix,iy-1,iz) * SA_frac_adv(:,ix,iy-1,iz,SA_ispec) - &
                                    fluxy(ix,iy  ,iz) * SA_frac_adv(:,ix,iy+1,iz,SA_ispec)
        ! outflow cells
        else if (fluxy(ix,iy-1,iz) <= 0 .and. fluxy(ix,iy,iz) > 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxy(ix,iy-1,iz) * SA_frac_adv(:,ix,iy  ,iz,SA_ispec) - &
                                    fluxy(ix,iy  ,iz) * SA_frac_adv(:,ix,iy  ,iz,SA_ispec)
        ! positive winds
        else if (fluxy(ix,iy-1,iz) > 0 .and. fluxy(ix,iy,iz) > 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxy(ix,iy-1,iz) * SA_frac_adv(:,ix,iy-1,iz,SA_ispec) - &
                                    fluxy(ix,iy  ,iz) * SA_frac_adv(:,ix,iy  ,iz,SA_ispec)
        endif

        ! calculate new fractions after x-advection
        denum = sum(SA_flux_adv(:,ix,iy,iz)) + LE_conc_abs(ix,iy,iz)
        if ( denum > numlim_conc ) then
          inv_denum = 1.0 / denum

          SA_frac(ix,iy,iz,SA_ispec,:) = &
                (SA_frac_adv(:,ix,iy,iz,SA_ispec) * LE_conc_abs(ix,iy,iz) + SA_flux_adv(:,ix,iy,iz) ) * inv_denum
        else
          SA_frac(ix,iy,iz,SA_ispec,:) = 0.0
        endif

      end do
      end do
      end do

      ! Set the calculated fractions in the help array for the advection routine
      ! check for negative fractions (numerical errors)
      ! Some compilers have trouble with where loops in OPEN-MP code
      ! changed to nx,ny,nz loop
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx
        do SA_ilabel = 1, SA_nlabel
          if (SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) < 0.0 ) then
            SA_frac_adv(SA_ilabel,ix,iy,iz,SA_ispec) = 0.0
            SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) = 0.0
          else 
            SA_frac_adv(SA_ilabel,ix,iy,iz,SA_ispec) = SA_frac(ix,iy,iz,SA_ispec,SA_ilabel)
          end if
        end do
      end do
      end do
      end do

      ! clear
      deallocate( SA_flux_adv, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( LE_conc_abs, stat=status )
      IF_NOTOK_RETURN(status=1)

    endif


    !ok
    status = 0

  end subroutine SA_Advec_y

! ---------------------------------------------------------------------------------------

  subroutine SA_Advec_z(ch, fluxz, volume_loc, ispec, status)

    use dims, only : nx, ny, nz

    ! --- in/out ------------------------
    real, intent(in), pointer ::  ch(:,:,:)         ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    real, intent(in), pointer ::  fluxz(:,:,:)      ! (1:nx,1:ny,0:nz)
    real, intent(in), pointer ::  volume_loc(:,:,:) ! (1-nh:nx+nh,1-nh:ny+nh,1-nh:nz+nh)
    integer, intent(in)       ::  ispec
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Advec_z'

    ! --- local ------------------
    integer                   ::  ix,iy,iz
    integer                   ::  SA_ispec, SA_ilabel
    real                      ::  denum, inv_denum
    real, allocatable         ::  SA_flux_adv(:,:,:,:)
    real, allocatable         ::  LE_conc_abs(:,:,:)

    ! --- begin -------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    if (SA_ispec > 0) then

      allocate ( SA_flux_adv(SA_nlabel,nx,ny,nz), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate ( LE_conc_abs(nx,ny,nz), stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! init:
      SA_flux_adv = 0.0
      ! state mass as help in this routine
      LE_conc_abs = ch(1:nx,1:ny,1:nz)*volume_loc(1:nx,1:ny,1:nz)

      ! loop over cells to make fluxes in each cell
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx

        ! calculate mass fluxes per cell per label (fluxes through cell edges are taken as fraction from donor cell)
        ! negative winds
        if      (fluxz(ix,iy,iz-1) <= 0 .and. fluxz(ix,iy,iz) <= 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxz(ix,iy,iz-1) * SA_frac_adv(:,ix,iy,iz  ,SA_ispec) - &
                                    fluxz(ix,iy,iz  ) * SA_frac_adv(:,ix,iy,iz+1,SA_ispec)
        ! inflow cells
        else if (fluxz(ix,iy,iz-1) > 0 .and. fluxz(ix,iy,iz) <= 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxz(ix,iy,iz-1) * SA_frac_adv(:,ix,iy,iz-1,SA_ispec) - &
                                    fluxz(ix,iy,iz  ) * SA_frac_adv(:,ix,iy,iz+1,SA_ispec)
        ! outflow cells
        else if (fluxz(ix,iy,iz-1) <= 0 .and. fluxz(ix,iy,iz) > 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxz(ix,iy,iz-1) * SA_frac_adv(:,ix,iy,iz  ,SA_ispec) - &
                                    fluxz(ix,iy,iz  ) * SA_frac_adv(:,ix,iy,iz  ,SA_ispec)
        ! positive winds
        else if (fluxz(ix,iy,iz-1) > 0 .and. fluxz(ix,iy,iz) > 0) then
          SA_flux_adv(:,ix,iy,iz) = fluxz(ix,iy,iz-1) * SA_frac_adv(:,ix,iy,iz-1,SA_ispec) - &
                                    fluxz(ix,iy,iz  ) * SA_frac_adv(:,ix,iy,iz  ,SA_ispec)
       endif

        ! calculate new fractions after x-advection
        denum = sum(SA_flux_adv(:,ix,iy,iz)) + LE_conc_abs(ix,iy,iz)
        if ( denum > numlim_conc ) then
          inv_denum = 1.0 / denum

          SA_frac(ix,iy,iz,SA_ispec,:) = &
                (SA_frac_adv(:,ix,iy,iz,SA_ispec) * LE_conc_abs(ix,iy,iz) + SA_flux_adv(:,ix,iy,iz) ) * inv_denum
        else
          SA_frac(ix,iy,iz,SA_ispec,:) = 0.0
        endif

      end do
      end do
      end do

      ! Set the calculated fractions in the help array for the advection routine
      ! check for negative fractions (numerical errors)
      ! Some compilers have trouble with where loops in OPEN-MP code
      ! changed to nx,ny,nz loop
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx
        do SA_ilabel = 1, SA_nlabel
          if (SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) < 0.0 ) then
            SA_frac_adv(SA_ilabel,ix,iy,iz,SA_ispec) = 0.0
            SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) = 0.0
          else 
            SA_frac_adv(SA_ilabel,ix,iy,iz,SA_ispec) = SA_frac(ix,iy,iz,SA_ispec,SA_ilabel)
          end if
        end do
      end do
      end do
      end do

      deallocate( LE_conc_abs, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( SA_flux_adv, stat=status )
      IF_NOTOK_RETURN(status=1)

    endif

    !ok
    status = 0

  end subroutine SA_Advec_z


!
! ---------------Horizontal diffusion ---------------------------------------------------
!

  subroutine SA_Hdiff_Init(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Hdiff_Init'

    ! --- begin -------------------------

    ! Allocate help variable for the fractions during the advection process
    allocate( SA_frac_hdiff(SA_nlabel,0:nx+1,0:ny+1,1:nz,nspec_labelled), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! initialize Fractions of the boundary conditions. Zero inside domain and one at the boundary
    ! These default labels can be override by the boundary setup if a run is used with a larger domain
    SA_frac_hdiff = 0.0

    ! x-boundary
    SA_frac_hdiff(Label_Bound,0   ,1:ny,1:nz,:) = 1
    SA_frac_hdiff(Label_Bound,nx+1,1:ny,1:nz,:) = 1

    ! y-boundary
    SA_frac_hdiff(Label_Bound,1:nx,0   ,1:nz,:) = 1
    SA_frac_hdiff(Label_Bound,1:nx,ny+1,1:nz,:) = 1

    ! ok
    status = 0

  end subroutine SA_Hdiff_Init

! ---------------------------------------------------------------------------------------

  subroutine SA_Hdiff_Done(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Hdiff_Done'

    ! --- begin -------------------------

    deallocate(SA_frac_hdiff, stat=status)
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Hdiff_Done

! ---------------------------------------------------------------------------------------

  subroutine SA_Hdiff_Setup(k,status)

    ! --- in/out ------------------------
    
    integer, intent(in)       ::  k
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Hdiff_Setup'
    
    ! --- local ------------------
    
    integer                   ::  SA_ispec, SA_ilabel
    
    ! --- begin -------------------------
    
    SA_ispec = Labelled_specs(k) 
    if ( SA_ispec > 0 ) then
      !set fractions inside domain in the help variable for advection process
      do SA_ilabel = 1, SA_nlabel
        SA_frac_hdiff(SA_ilabel,1:nx,1:ny,1:nz,SA_ispec) = SA_frac(1:nx,1:ny,1:nz,SA_ispec,SA_ilabel)
      end do
    end if

    ! ok
    status = 0

  end subroutine SA_Hdiff_setup

! ---------------------------------------------------------------------------------------

  subroutine SA_Hdiff_Apply(ch, fluxx, dx, fluxy, dy, iz, ispec, status)

    use dims, only :  nx, ny, nz

    ! --- in/out ------------------------

    real, intent(in)          ::  ch(1:nx,1:ny)    ! (1:nx,1:ny)
    real, intent(in)          ::  fluxx(0:nx,1:ny) ! (0:nx,1:ny)
    real, intent(in)          ::  dx(1:nx,1:ny)          ! (1:nx,1:ny)
    real, intent(in)          ::  fluxy(1:nx,0:ny) ! (1:nx,0:ny)
    real, intent(in)          ::  dy(1:nx,1:ny)          ! (1:nx,1:ny)
    integer, intent(in)       ::  iz
    integer, intent(in)       ::  ispec
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Hdiff_Apply'

    ! --- local ------------------

    integer                   ::  ix,iy
    integer                   ::  SA_ispec, SA_ilabel
    real                      ::  denum, inv_denum
    real, allocatable         ::  SA_flux_hdiff(:,:,:)
    real, allocatable         ::  LE_conc_abs(:,:)

    ! --- begin -------------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)
    !write(gol, '("dx: ", i0," ",i0)' ) shape(dx,1), shape(dx,2) ; call goPr
    !write(gol, '("dy: ", i0," ",i0)' ) shape(dy,1), shape(dy,2) ; call goPr
    if (SA_ispec > 0) then

      ! storage:
      allocate( SA_flux_hdiff(SA_nlabel,nx,ny), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( LE_conc_abs(1:nx,1:ny), stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! init:
      SA_flux_hdiff = 0.0
      ! state mass as help in this routine
      LE_conc_abs = ch(1:nx,1:ny)

      ! loop over cells to make fluxes in each cell
      do iy = 1, ny
      do ix = 1, nx
        ! calculate mass fluxes per cell per label (fluxes through cell edges are taken as fraction from donor cell)
        ! negative winds
        if      (fluxx(ix-1,iy) <= 0 .and. fluxx(ix,iy) <= 0) then
          SA_flux_hdiff(:,ix,iy) = ( fluxx(ix-1,iy) * SA_frac_hdiff(:,ix  ,iy,iz,SA_ispec) - &
                                     fluxx(ix  ,iy) * SA_frac_hdiff(:,ix+1,iy,iz,SA_ispec) ) / dx(ix,iy)
        ! inflow cells
        else if (fluxx(ix-1,iy) > 0 .and. fluxx(ix,iy) <= 0) then
          SA_flux_hdiff(:,ix,iy) = ( fluxx(ix-1,iy) * SA_frac_hdiff(:,ix-1,iy,iz,SA_ispec) - &
                                     fluxx(ix  ,iy) * SA_frac_hdiff(:,ix+1,iy,iz,SA_ispec) ) / dx(ix,iy)
        ! outflow cells
        else if (fluxx(ix-1,iy) <= 0 .and. fluxx(ix,iy) > 0) then
          SA_flux_hdiff(:,ix,iy) = ( fluxx(ix-1,iy) * SA_frac_hdiff(:,ix  ,iy,iz,SA_ispec) - &
                                     fluxx(ix  ,iy) * SA_frac_hdiff(:,ix  ,iy,iz,SA_ispec) ) / dx(ix,iy)
        ! positive winds
        else if (fluxx(ix-1,iy) > 0 .and. fluxx(ix,iy) > 0) then
          SA_flux_hdiff(:,ix,iy) = ( fluxx(ix-1,iy) * SA_frac_hdiff(:,ix-1,iy,iz,SA_ispec) - &
                                     fluxx(ix  ,iy) * SA_frac_hdiff(:,ix  ,iy,iz,SA_ispec) ) / dx(ix,iy)
        endif

        ! negative winds
        if      (fluxy(ix,iy-1) <= 0 .and. fluxy(ix,iy) <= 0) then
          SA_flux_hdiff(:,ix,iy) = SA_flux_hdiff(:,ix,iy) + &
                                   ( fluxy(ix,iy-1) * SA_frac_hdiff(:,ix,iy  ,iz,SA_ispec) - &
                                     fluxy(ix,iy  ) * SA_frac_hdiff(:,ix,iy+1,iz,SA_ispec) ) /dy(ix,iy)
        ! inflow cells
        else if (fluxy(ix,iy-1) > 0 .and. fluxy(ix,iy) <= 0) then
          SA_flux_hdiff(:,ix,iy) = SA_flux_hdiff(:,ix,iy) + &
                                   ( fluxy(ix,iy-1) * SA_frac_hdiff(:,ix,iy-1,iz,SA_ispec) - &
                                     fluxy(ix,iy  ) * SA_frac_hdiff(:,ix,iy+1,iz,SA_ispec) ) /dy(ix,iy)
        ! outflow cells
        else if (fluxy(ix,iy-1) <= 0 .and. fluxy(ix,iy) > 0) then
          SA_flux_hdiff(:,ix,iy) = SA_flux_hdiff(:,ix,iy) + &
                                   ( fluxy(ix,iy-1) * SA_frac_hdiff(:,ix,iy  ,iz,SA_ispec) - &
                                     fluxy(ix,iy  ) * SA_frac_hdiff(:,ix,iy  ,iz,SA_ispec) ) /dy(ix,iy)
        ! positive winds
        else if (fluxy(ix,iy-1) > 0 .and. fluxy(ix,iy) > 0) then
          SA_flux_hdiff(:,ix,iy) = SA_flux_hdiff(:,ix,iy) + &
                                   ( fluxy(ix,iy-1) * SA_frac_hdiff(:,ix,iy-1,iz,SA_ispec) - &
                                     fluxy(ix,iy  ) * SA_frac_hdiff(:,ix,iy  ,iz,SA_ispec) ) /dy(ix,iy)
        endif

        ! calculate new fractions after x-advection
        denum = sum(SA_flux_hdiff(:,ix,iy)) + LE_conc_abs(ix,iy)
        if ( denum > numlim_conc ) then
          inv_denum = 1.0 / denum
          SA_frac(ix,iy,iz,SA_ispec,:) = &
              (SA_frac_hdiff(:,ix,iy,iz,SA_ispec) * LE_conc_abs(ix,iy) + SA_flux_hdiff(:,ix,iy) ) * inv_denum

        else
          SA_frac(ix,iy,iz,SA_ispec,:) = 0.0
        endif

      end do
      end do
      
!      ! Set the calculated fractions in the help array for the advection routine
!      ! check for negative fractions (numerical errors)
!      ! Some compilers have trouble with where loops in OPEN-MP code
!      ! changed to nx,ny,nz loop
!      do iy = 1, ny
!      do ix = 1, nx
!        do SA_ilabel = 1, SA_nlabel
!          if (SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) < 0.0 ) then
!            SA_frac_hdiff(SA_ilabel,ix,iy,iz,SA_ispec) = 0.0
!            SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) = 0.0
!          else 
!            SA_frac_hdiff(SA_ilabel,ix,iy,iz,SA_ispec) = SA_frac(ix,iy,iz,SA_ispec,SA_ilabel)
!          end if
!        end do
!      end do
!      end do
      
      ! clear:
      deallocate( SA_flux_hdiff, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( LE_conc_abs, stat=status )
      IF_NOTOK_RETURN(status=1)

    endif

    !ok
    status = 0

  end subroutine SA_Hdiff_Apply
  
  
!
! ----------------ADJUST LAYER DEPTHS----------------------------------------------------
!

  subroutine SA_Adjust_Init(status)

    ! --- in/out -----------------------
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Adjust_Init'

    ! --- begin ------------------------

    ! Allocate and initialize help arrays for Adjust routine
    allocate(SA_caloft_frac (1:SA_nlabel)      , stat=status)
    IF_NOTOK_RETURN(status=1)
    allocate(int_m_SA       (0:nz+1, SA_nlabel), stat=status)
    IF_NOTOK_RETURN(status=1)
    allocate(SA_mass        (1:SA_nlabel)      , stat=status)
    IF_NOTOK_RETURN(status=1)
    allocate(SA_conc_adj    (nz,SA_nlabel)     , stat=status)
    IF_NOTOK_RETURN(status=1)

    SA_caloft_frac(:)              = 0.0
    SA_caloft_frac(Label_Aloft)    = 1.0

    int_m_SA                       = 0.0
    SA_mass                        = 0.0

    ! ok
    status = 0

  end subroutine SA_Adjust_Init

! ---------------------------------------------------------------------------------------

  subroutine SA_Adjust_Done(status)

    ! --- in/out -----------------------
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Adjust_Done'

    ! --- begin ------------------------

    ! clear help arrays
    deallocate(SA_caloft_frac, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(int_m_SA      , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(SA_mass       , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(SA_conc_adj   , stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Adjust_Done

! ---------------------------------------------------------------------------------------

  subroutine SA_Adjust_Setup(i,j,l,k,caloft,SA_old_volume, status)

    ! --- in/out -----------------------
    integer, intent(in)         ::  i, j, l, k
    real, intent(in), optional  ::  caloft      ! Concentration at top (for cell i,j and spec k)
    real, intent(in), optional  ::  SA_old_volume
    integer, intent(out)        ::  status

    ! --- const ------------------------
    character(len=*), parameter :: rname = mname//'/SA_Adjust_Setup'

    ! --- local ------------------------
    integer                     :: SA_ispec

    SA_ispec = Labelled_specs(k)
    if ( SA_ispec > 0 ) then
      if ( l > nz ) then ! aloft layer
        ! Integrated mass upto top layer
        int_m_SA(l,:) = int_m_SA(l-1,:)
        int_m_SA(l,Label_Aloft) = SA_old_volume*caloft * SA_caloft_frac (Label_Aloft) + int_m_SA(l-1,Label_aloft)
        
        ! initialize total mass in column
        SA_mass = 0.0
      else
        ! integrated mass upto current layer
        int_m_SA(l,:) = ovolume(i,j,l) * LE_conc(i,j,l,SA_ispec)*SA_frac(i,j,l,SA_ispec,:)+ int_m_SA(l-1,:)
      endif
    endif

    ! ok
    status = 0

  end subroutine SA_Adjust_Setup


! ---------------------------------------------------------------------------------------

  subroutine SA_Adjust_Conc(i,j,l,k,ix,dpvratio,volume,status)

    ! --- in/out -----------------------
    integer, intent(in)         ::  i, j, l,k, ix
    real, intent(in)            ::  dpvratio
    real, intent(in)            ::  volume
    integer, intent(out)        ::  status

    ! --- const ------------------------
    character(len=*), parameter :: rname = mname//'/SA_Adjust_Conc'

    ! --- local ------------------------
    integer                     ::  SA_ispec

    ! --- begin ------------------------

    SA_ispec = Labelled_specs(k)

    if (SA_ispec > 0) then
      ! Determine new masses (SA_conc_adj now in mass term)
      SA_conc_adj(l,:) = int_m_SA(ix-1,:) + ( int_m_SA(ix,:) - int_m_SA(ix-1,:) ) * dpvratio - SA_mass(:)

      ! update the cumulative mass assigned so far
      SA_mass(:) = SA_mass(:) + SA_conc_adj(l,:)

      ! divide by volume in order to get concentration: SA_conc_adj will now become a concentration term
      SA_conc_adj(l,:) = SA_conc_adj(l,:) / volume

    end if

    ! ok
    status = 0

  end subroutine SA_Adjust_Conc

! ---------------------------------------------------------------------------------------

  subroutine SA_Adjust_Fractions(i,j,k,status)

    ! --- in/out -----------------------
    integer, intent(in)         ::  i,j,k
    integer, intent(out)        ::  status

    ! --- const ------------------------
    character(len=*), parameter :: rname = mname//'/SA_Adjust_Fractions'

    ! --- local ------------------------
    integer                     ::  SA_ispec

    ! --- begin ------------------------

    integer                   ::  iz
    real                      ::  denum, inv_denum

    ! --- begin -------------------------

    SA_ispec = Labelled_specs(k)

    if (SA_ispec > 0 ) then
        ! calculate the fractions of each label by dividing concentration per label through tot concentration
      do iz = 1, nz
        denum = sum(SA_conc_adj(iz,:))
        if ( denum > numlim_conc ) then
          inv_denum = 1.0/denum
          SA_frac(i,j,iz,SA_ispec,:) = SA_conc_adj(iz,:) * inv_denum
        else
          SA_frac(i,j,iz,SA_ispec,:) = 0.0
        endif
      enddo
    endif

    ! ok
    status = 0

  end subroutine SA_Adjust_Fractions

!
! ---------------- VERTICAL DIFFUSION ---------------------------------------------------
!
  subroutine SA_Vdiff_Init(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Vdiff_Init'

    ! --- begin -------------------------

    ! allocate help arrays for this process
    allocate( SA_conc_vdiff(nz,nspec_labelled,SA_nlabel), stat = status )
    IF_NOTOK_RETURN(status=1)

    SA_conc_vdiff = 0.0

    ! ok
    status = 0

  end subroutine SA_Vdiff_Init

! ---------------------------------------------------------------------------------------

  subroutine SA_Vdiff_Done(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Vdiff_Done'

    ! --- begin -------------------------

    ! deallocate help arrays for this process
    deallocate( SA_conc_vdiff, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Vdiff_Done

! ---------------------------------------------------------------------------------------

  subroutine SA_Vdiff_Setup(i,j,status )

    ! --- in/out ------------------------
    integer, intent(in)       ::  i, j
    integer, intent(out)      ::  status

    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Vdiff_Setup'

    ! --- local -------------------------
    integer                   ::  SA_ispec, ispec
    integer                   ::  iz

    ! --- begin -------------------------
    
    
    do ispec = 1, nspec 
      ! help concentrations
      SA_ispec = Labelled_specs(ispec)

      if (SA_ispec > 0) then
        do iz = 1, nz
          ! fill help concentration
          SA_conc_vdiff(iz,SA_ispec,:) = LE_conc(i,j,iz,SA_ispec) * SA_frac(i,j,iz,SA_ispec,:)
        end do
    
      end if
    end do
  
    ! ok
    status = 0

  end subroutine SA_Vdiff_Setup

! ---------------------------------------------------------------------------------------

  subroutine SA_Vdiff_Conc(i,j,k,dtvd,dh,hh,caloft,diffmat, status)

    ! --- in/out ------------------------
    integer, intent(in)       ::  i, j, k
    real, intent(in)          ::  dtvd, dh, hh
    real, intent(in)          ::  caloft
    real, intent(inout)       ::  diffmat(nz,nz)
    integer, intent(out)      ::  status

    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Vdiff_Conc'

    ! --- local -------------------------
    real, allocatable         ::  rhs(:,:)
    integer                   ::  iz, SA_ilabel
    integer                   ::  SA_ispec
    real                      ::  pivot

    ! --- begin -------------------------

    SA_ispec = Labelled_specs(k)

    if (SA_ispec > 0 ) then
      allocate( rhs(nz,SA_nlabel),stat=status)
      rhs = 0.0

      ! implicit solve VDM * c_new = c_old (old concentration in rhs)
      do iz = 1, nz
        do SA_ilabel = 1, SA_nlabel
          rhs(iz,SA_ilabel) = SA_conc_vdiff(iz,SA_ispec,SA_ilabel)
        enddo
      enddo

      rhs(nz,Label_Aloft) = SA_conc_vdiff(nz,SA_ispec,Label_Aloft) + kz(i,j,nz)/dh/hh*dtvd*caloft
      ! Transform Diffmat into bi-diagonal form
      do iz=2,nz
        if (diffmat(iz-1,iz-1) == 0.0) then
          write (*,*) iz,diffmat(iz-1,iz-1)
          write (*,*) 'singular diffusion matrix!'
        endif
        pivot  = diffmat(iz,iz-1)/diffmat(iz-1,iz-1)
        diffmat(iz,iz) = diffmat(iz,iz)   - pivot*diffmat(iz-1,iz)
        if (iz /= nz) diffmat(iz,iz+1) = diffmat(iz,iz+1) - pivot*diffmat(iz-1,iz+1)
        rhs(iz,:) = rhs(iz,:) - pivot*rhs(iz-1,:)
      enddo

      ! performe backsolve
      SA_conc_vdiff(nz,SA_ispec,:) = rhs(nz,:)/diffmat(nz,nz)
      do iz=nz-1,1,-1
         SA_conc_vdiff(iz,SA_ispec,:) = &
              (rhs(iz,:)-diffmat(iz,iz+1)*SA_conc_vdiff(iz+1,SA_ispec,:))/diffmat(iz,iz)
      enddo

      ! clear
      deallocate( rhs, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! ok
    status = 0

  end subroutine SA_Vdiff_Conc

! ---------------------------------------------------------------------------------------

  subroutine SA_Vdiff_Fractions(i,j,status)

    ! --- in/out ------------------------
    integer, intent(in)       ::  i, j
    integer, intent(out)      ::  status

    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Vdiff_Conc'

    ! --- local -------------------------
    integer                   ::  SA_ispec, ispec
    integer                   ::  iz
    real, allocatable         ::  denum(:)
    real                      ::  inv_denum

    ! --- begin -------------------------

    
    do ispec = 1, nspec 
      ! total concentration in the denumerator
      SA_ispec = Labelled_specs(ispec)

      if ( SA_ispec > 0 ) then
        allocate( denum(nz), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! calculate fractions for each cell as the conc per label divide d by total conc
        denum = sum(SA_conc_vdiff(:,SA_ispec,:),dim=2)

        do iz = 1, nz
          if ( denum(iz) > numlim_conc ) then
            inv_denum = 1.0 / denum(iz)
            SA_frac(i,j,iz,SA_ispec,:) = SA_conc_vdiff(iz,SA_ispec,:) * inv_denum
  
            where(SA_frac(i,j,iz,SA_ispec,:) < 0.0 )
              SA_frac(i,j,iz,SA_ispec,:) = 0.0
            endwhere
          else
            SA_frac(i,j,iz,SA_ispec,:) = 0.0
          endif
        enddo

        deallocate( denum, stat=status )
        IF_NOTOK_RETURN(status=1)
      end if
    end do
 
    ! ok
    status = 0

  end subroutine SA_Vdiff_Fractions

!
! ---------------- SEDIMENTATION --------------------------------------------------------
!

  subroutine SA_Sedim_Init(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Sedim_Init'

    ! --- begin ------------------

    ! allocate and initialize help variables for this routine
    allocate( SA_conc_sedim( nz), stat = status )
    IF_NOTOK_RETURN(status=1)
    allocate( SA_frac_sedim( nz, SA_nlabel), stat = status )
    IF_NOTOK_RETURN(status=1)

    SA_conc_sedim = 0.0
    SA_frac_sedim = 0.0

    ! ok
    status = 0

  end subroutine SA_Sedim_Init

! ---------------------------------------------------------------------------------------

  subroutine SA_Sedim_Done(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Sedim_Done'

    ! --- begin ------------------

    ! clear help variables:
    deallocate( SA_conc_sedim, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( SA_frac_sedim, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine SA_Sedim_Done

! ---------------------------------------------------------------------------------------

  subroutine SA_Sedim_Setup(c,i,j,ispec, status)

    ! --- in/out ------------------------
    real, intent(in)          ::  c(nz)
    integer, intent(in)       ::  i,j, ispec
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Sedim_Setup'

    ! --- local ------------------
    integer                   :: SA_ispec

    ! --- begin ------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)

    if (SA_ispec > 0) then
      ! state help concentration and help array with fractions
      SA_conc_sedim(:)   = c
      SA_frac_sedim(:,:) = SA_frac(i,j,:,SA_ispec,:)
    endif

    ! ok
    status = 0

  end subroutine SA_Sedim_Setup

! ---------------------------------------------------------------------------------------

  subroutine SA_Sedim_Conc(i,j,ispec,flux,c, dh,status)

    ! --- in/out ------------------------
    integer, intent(in)       ::  i, j
    integer, intent(in)       ::  ispec
    real, intent(in)          ::  flux(nz+1)
    real, intent(in)          ::  c(nz)
    real, intent(in)          ::  dh(nz)
    integer, intent(out)      ::  status

    ! --- const ------------------
    character(len=*), parameter :: rname = mname//'/SA_Sedim_Conc'

    ! --- local ------------------
    real                      ::  inv_dh
    real                      ::  inv_c
    integer                   ::  SA_ispec, SA_ispec_budget, SA_ispec_drydep, SA_ilabel
    integer                   ::  iz

    ! --- begin ------------------

    ! Search corresponding specie parameter in list of labelled species
    SA_ispec = Labelled_specs(ispec)
    
    ! Drydep of budget output needed?
    SA_ispec_budget = Labelled_specs_budget_output(ispec)
    SA_ispec_drydep = Labelled_specs_budget_output(ispec)
    
    ! Calculate fractions with old fractions and fluxes from layer to layer
    if (SA_ispec > 0) then

      do iz = 1, nz
        if (c(iz) > numlim_conc ) then
          ! Use inverse values (less calculation time??)
          inv_c  = 1.0 / c(iz)
          inv_dh = 1.0 / dh(iz)
          if (iz < nz ) then
            SA_frac_sedim(iz,:) = ( SA_frac(i,j,iz,SA_ispec,:) * SA_conc_sedim(iz) + &
                                         ( flux(iz+1) * SA_frac(i,j,iz+1,SA_ispec,:) - flux(iz) * SA_frac(i,j,iz,SA_ispec,:) ) * inv_dh ) * inv_c

            ! update labelled dry deposition budget for sedimentation process
            if ( iz == 1 ) then
              ! Daily budget output needed?
              if ( SA_ispec_budget > 0 ) then
                SA_budget_drydepos_day(i,j,SA_ispec_budget,:)  = SA_budget_drydepos_day(i,j,SA_ispec_budget,:) + &
                                                                  ( flux(iz) * inv_dh * SA_frac(i,j,iz,SA_ispec,:) )
              end if
              ! hourly drydep output needed?
              if ( SA_ispec_drydep > 0 ) then
                SA_budget_drydepos_hour(i,j,SA_ispec_drydep,:) = SA_budget_drydepos_hour(i,j,SA_ispec_drydep,:) + &
                                                                  ( flux(iz) * inv_dh * SA_frac(i,j,iz,SA_ispec,:) )
              end if
            end if
          else
            SA_frac_sedim(iz,:) = ( SA_frac(i,j,iz,SA_ispec,:) * SA_conc_sedim(iz) + &
                                         ( flux(iz+1) * SA_caloft_frac(:) - flux(iz) * SA_frac(i,j,iz,SA_ispec,:) ) * inv_dh ) * inv_c
          endif
        else
          SA_frac_sedim(iz,:) = 0.0
        endif
      enddo

      ! Restate help variables
      do iz = 1, nz

        ! Fractions
        do SA_ilabel = 1, SA_nlabel
          if ( SA_frac_sedim(iz,SA_ilabel) > 0.0 ) then
            SA_frac(i,j,iz,SA_ispec,SA_ilabel) = SA_frac_sedim(iz,SA_ilabel)
          else
            SA_frac(i,j,iz,SA_ispec,SA_ilabel) = 0.0
          end if
        end do  ! labels
        
        ! concentrations
        SA_conc_sedim(iz) = c(iz)
        
      end do  ! iz
            
    end if  ! Tracer labelled?

    ! ok
    status = 0

  end subroutine SA_Sedim_Conc

  
!
! ---------------- DEPOSITION ROUTINES ----------------------------------------------------
!  
  
  subroutine SA_Depos_Init(rcF, status)
    
    use GO, only : ReadRc, TrcFile
    use GO, only : GoSplitString, GoMatchValues
    use LE_Landuse_Data, only : nlu, lu_name_abbr

    ! --- in/out --------------------
    type(TrcFile), intent(in) ::  rcF
    integer, intent(out)      ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Depos_Init'

    ! --- local -------------------------
    character(len=64)         ::  budget_label_files
    character(len=512)        ::  drydep_tracer_line, wetdep_tracer_line
    character(len=512)        ::  budget_fields_line
    integer                   ::  ibudget_field
    character(len=512)        ::  budget_tracer_line, budget_lu_line    
    integer                   ::  ispec,ilu
    
    ! --- begin -------------------------
    
    ! Hourly budget fields (drydep and wetdep fluxes)
    call ReadRc( rcF, 'le.output.drydep_label.drydepo-label.fields', drydep_tracer_line, status )
    IF_NOTOK_RETURN(status=1)
    
    ! Allocate necessary arrays      
    allocate( labelled_specs_drydep_output_names(nspec), stat=status )
    allocate( labelled_specs_drydep_output_numbers(nspec), stat=status )
    allocate( labelled_specs_drydep_output(nspec), stat=status)
    labelled_specs_drydep_output = -1
    
    if ( len(trim(drydep_tracer_line) ) > 0 ) then
      ! match species for drydep-label-output within list of all species
      call goMatchValues(drydep_tracer_line, specname,  &
                         ndrydep_label_output_tracer, labelled_specs_drydep_output_names, labelled_specs_drydep_output_numbers, &
                         status)
      IF_NOTOK_RETURN(status=1)
      do ispec = 1, ndrydep_label_output_tracer
        Labelled_specs_drydep_output(labelled_specs_drydep_output_numbers(ispec)) = ispec
      enddo    
      allocate( SA_budget_drydepos_hour(nx,ny,ndrydep_label_output_tracer,SA_nlabel), stat = status)
      IF_NOTOK_RETURN(status=1)
      SA_budget_drydepos_hour = 0.0
    else 
      ndrydep_label_output_tracer = 0
    end if
    
    ! Hourly budget fields (drydep and wetdep fluxes)
    call ReadRc( rcF, 'le.output.wetdep_label.wetdepo-label.fields', wetdep_tracer_line, status )
    IF_NOTOK_RETURN(status=1)
    
    ! Allocate necessary arrays      
    allocate( labelled_specs_wetdep_output_names(nspec), stat=status )
    allocate( labelled_specs_wetdep_output_numbers(nspec), stat=status )
    allocate( labelled_specs_wetdep_output(nspec), stat=status)
    labelled_specs_wetdep_output = -1
    
    if ( len(trim(wetdep_tracer_line)) > 0 ) then
      ! match species for wetdep-label-output within list of all species
      call goMatchValues(wetdep_tracer_line, specname,  &
                         nwetdep_label_output_tracer, labelled_specs_wetdep_output_names, labelled_specs_wetdep_output_numbers, &
                         status)
      IF_NOTOK_RETURN(status=1)
      do ispec = 1, nwetdep_label_output_tracer
        Labelled_specs_wetdep_output(labelled_specs_wetdep_output_numbers(ispec)) = ispec
      enddo

      allocate( SA_budget_wetdepos_hour(nx,ny,nwetdep_label_output_tracer,SA_nlabel), stat = status)
      IF_NOTOK_RETURN(status=1)
      SA_budget_wetdepos_hour = 0.0
    else
      nwetdep_label_output_tracer = 0
    endif
        
    ! Check which tracers should be in budget output
    ! Save memory by only saving those
    ! save calculation time by onlt calculating those 
    ! Allocate necessary arrays      
    allocate( labelled_specs_budget_output_names(nspec), stat=status )
    allocate( labelled_specs_budget_output_numbers(nspec), stat=status )
    allocate( labelled_specs_budget_output(nspec), stat=status)
    labelled_specs_budget_output = -1

    allocate( landuse_budget_output_names(nlu), stat=status )
    allocate( landuse_budget_output_numbers(nlu), stat=status )
    allocate( landuse_budget_output(nlu), stat=status )
    landuse_budget_output = -1
    
    ! Check if we need budget label output
    call ReadRc( rcF, 'le.output.budget_label.files', budget_label_files, status )
    IF_NOTOK_RETURN(status=1)
    
    if ( len_trim(budget_label_files) > 0 ) then

      ! daily budget fields (supported: drydep, wetdep, drydep_lu)
      call ReadRc( rcF, 'le.output.budget_label.budget-label.fields', budget_fields_line, status )
      IF_NOTOK_RETURN(status=1)
      call GoSplitString( budget_fields_line, nbudget_fields, budget_fields, status )
      IF_NOTOK_RETURN(status=1)


      do ibudget_field = 1, nbudget_fields

        ! Dry deposition budgets, 
        if ( trim(budget_fields(ibudget_field)) == 'drydep' ) then
          ! Find which tracers should be in output
          call ReadRc( rcF, 'le.output.budget_label.budget-label.tracers', budget_tracer_line, status )
          IF_NOTOK_RETURN(status=1)
          ! match species for budget-label-output within list of all species
          call goMatchValues(budget_tracer_line, specname,  &
                             nbudget_label_output_tracer, labelled_specs_budget_output_names, labelled_specs_budget_output_numbers, &
                             status)
          IF_NOTOK_RETURN(status=1)
          do ispec = 1, nbudget_label_output_tracer
            Labelled_specs_budget_output(labelled_specs_budget_output_numbers(ispec)) = ispec
          enddo

          ! Allocate an initialize
          allocate( SA_budget_drydepos_day(nx,ny,nbudget_label_output_tracer,SA_nlabel), stat = status)
          IF_NOTOK_RETURN(status=1)
          SA_budget_drydepos_day = 0.0

        ! Dry deposition budgets per landuse
        else if ( trim(budget_fields(ibudget_field)) == 'drydep_lu' ) then

          ! Find which tracers should be in output
          call ReadRc( rcF, 'le.output.budget_label.budget-label.tracers', budget_tracer_line, status )
          IF_NOTOK_RETURN(status=1)          
          ! match species for budget-label-output within list of all species
          call goMatchValues(budget_tracer_line, specname,  &
                             nbudget_label_output_tracer, labelled_specs_budget_output_names, labelled_specs_budget_output_numbers, &
                             status)
          IF_NOTOK_RETURN(status=1)
          do ispec = 1, nbudget_label_output_tracer
            Labelled_specs_budget_output(labelled_specs_budget_output_numbers(ispec)) = ispec
          enddo

          ! Find which landuses should be in output
          call ReadRc( rcF, 'le.output.budget_label.budget-label.landuses', budget_lu_line, status )
          IF_NOTOK_RETURN(status=1)          

          ! match landuse for budget-label-output within list of all landuses
          call goMatchValues(budget_lu_line, lu_name_abbr,  &
                             nbudget_label_output_lu, landuse_budget_output_names, landuse_budget_output_numbers, &
                             status)
          IF_NOTOK_RETURN(status=1)          
          do ilu = 1, nbudget_label_output_lu
            landuse_budget_output(landuse_budget_output_numbers(ilu)) = ilu
          enddo                   

          ! Allocate and initialize
          allocate( SA_budget_drydepos_lu_day(nx,ny,nbudget_label_output_tracer,nbudget_label_output_lu,SA_nlabel), stat=status)
          IF_NOTOK_RETURN(status=1)
          SA_budget_drydepos_lu_day = 0.0

        ! Wet deposition budgets
        else if ( trim( budget_fields(ibudget_field)) == 'wetdep' ) then
          ! Find which tracers should be in output
          call ReadRc( rcF, 'le.output.budget_label.budget-label.tracers', budget_tracer_line, status )
          IF_NOTOK_RETURN(status=1)

          ! match species for budget-label-output within list of all species
          call goMatchValues(budget_tracer_line, specname,  &
                             nbudget_label_output_tracer, labelled_specs_budget_output_names, labelled_specs_budget_output_numbers, &
                             status)
          IF_NOTOK_RETURN(status=1)
          do ispec = 1, nbudget_label_output_tracer
            Labelled_specs_budget_output(labelled_specs_budget_output_numbers(ispec)) = ispec
          enddo

          ! Allocate and initialize
          allocate( SA_budget_wetdepos_day(nx,ny,nbudget_label_output_tracer,SA_nlabel), stat = status)
          IF_NOTOK_RETURN(status=1)
          SA_budget_wetdepos_day = 0.0

        ! unsupported field types
        else 
          write( gol, '("Unknown data type for budget-label output: ", a)' ) trim(budget_fields(ibudget_field)) ; call goErr
          TRACEBACK;status=1;return
        end if
      end do ! budget fields         
    end if
        
    ! ok
    status = 0

  end subroutine SA_Depos_Init
  
! ----------------------------------------------

  subroutine SA_Depos_Done(status)
  
    ! --- in/out --------------------
    integer, intent(out)      ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_DryDepos_Done'

    ! --- begin -------------------------
    
    ! Done
    deallocate( labelled_specs_drydep_output_names, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( labelled_specs_drydep_output_numbers, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( labelled_specs_drydep_output, stat=status)
    IF_NOTOK_RETURN(status=1)

    deallocate( labelled_specs_wetdep_output_names, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( labelled_specs_wetdep_output_numbers, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( labelled_specs_wetdep_output, stat=status)
    IF_NOTOK_RETURN(status=1)

    deallocate( labelled_specs_budget_output_names, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( labelled_specs_budget_output_numbers, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( labelled_specs_budget_output, stat=status)
    IF_NOTOK_RETURN(status=1)

    deallocate( landuse_budget_output_names, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( landuse_budget_output_numbers, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( landuse_budget_output, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Depos_Done

! ----------------------------------------------

  subroutine SA_DryDepos_Fractions( budget_depos, status)
    
    use Binas, only : xm_air
    use Dims, only : nx, ny
    use Indices, only : specunit, specname

    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out --------------------
    real, intent(in)          ::  budget_depos(nx,ny,nspec, 2) ! Deposition budget (unit is ppb m or ug m-2)
    integer, intent(out)      ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_DryDepos_Fractions'
    
    ! --- local ---------------------
    integer                   ::  ispec, ix,iy,iz
    integer                   ::  SA_ilabel
    integer                   ::  SA_ispec
    real                      ::  conv_fact(nx,ny)
    real                      ::  SA_conc_before_dryemis(nx,ny,SA_nlabel)
    integer                   ::  SA_ispec_budget, SA_ispec_drydep
    
    ! meteo data:
    real, pointer          ::  dens(:,:,:)   ! (lon,lat,alt)    
    real, pointer          ::  h_m(:,:,:)   ! (lon,lat,alt)

    ! --- begin ---------------------
            
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m')    
    IF_NOTOK_RETURN(status=1)

    ! surface for dry deposition
    iz = 1
    
    do ispec = 1, nspec     

      ! Search corresponding specie parameter in list of labelled species
      SA_ispec = Labelled_specs(ispec)            
      
      ! Search indices for output in budget/drydep
      SA_ispec_budget = Labelled_specs_budget_output(ispec)
      SA_ispec_drydep = Labelled_specs_drydep_output(ispec)
      
      if (SA_ispec > 0 ) then
        
        ! Conversion factor from ppb m to ug/m2? or not needed for aerosols
        select case ( trim(specunit(ispec)) )
          !~ other notations ...
          case ( 'ug/m3' )
            ! aerosol tracer (no conversion needed)
            conv_fact = 1.0
          case ( 'ppb' )
            ! Conversion from ppb m to ug/m2
            !           ((mole tr/mole air)/ppb) * (kg tr/mole tr)   * (kg air/m3 air) / (kg air/mole air) * ug/kg
            conv_fact =       1e-9               *  specmolm(ispec)  * dens(1:nx,1:ny,iz)  /     xm_air        * 1e9
          case default

            write (gol,'("do not know how to convert `",a,"` from `",a," m` to `ug m-2")') &
                  trim(specname(ispec)), trim(specunit(ispec)); call goErr
                  TRACEBACK; status=1; return
        end select ! spec unit

        ! output deposition budget needed?
        if ( SA_ispec_budget > 0 ) then
                  
          do SA_ilabel = 1, SA_nlabel
            SA_budget_drydepos_day(:,:,SA_ispec_budget,SA_ilabel)  = SA_budget_drydepos_day(:,:,SA_ispec_budget,SA_ilabel) + &
                                                            budget_depos(:,:,ispec,iexdry_depo) * conv_fact * SA_frac(:,:,iz,SA_ispec,SA_ilabel)
          end do
          
        end if ! SA_ispec_budget
        
        ! Output for dry deposition fluxes needed
        if ( SA_ispec_drydep > 0 ) then
        
          do SA_ilabel = 1, SA_nlabel
              SA_budget_drydepos_hour(:,:,SA_ispec_drydep,SA_ilabel) = SA_budget_drydepos_hour(:,:,SA_ispec_drydep,SA_ilabel) + &
                                                            budget_depos(:,:,ispec,iexdry_depo) * conv_fact * SA_frac(:,:,iz,SA_ispec,SA_ilabel)
          end do                          
        end if ! SA_ispec_drydep
        
        ! Re-emission? than update label fractions with compensation point
        ! For new label fractions assume that we first calculate the deposition and than the re-emission
        do iy = 1, ny
        do ix = 1, nx

          if ( budget_depos(ix,iy,ispec,iexdry_emis) > 0.0 ) then
          
            ! calculate absolute concentrations in lowest layer
            do SA_ilabel = 1, SA_nlabel                    
              SA_conc_before_dryemis(ix,iy,SA_ilabel) = ( LE_conc(ix,iy,iz,SA_ispec) + budget_depos(ix,iy,ispec,iexdry_depo)/h_m(ix,iy,iz) ) * SA_frac(ix,iy,iz,SA_ispec,SA_ilabel)
            end do

            ! update concentrations in lowest layer  ( divide with height to get ug/m3 or ppb)
            LE_conc(ix,iy,iz,SA_ispec) = LE_conc(ix,iy,iz,SA_ispec) + sum( budget_depos(ix,iy,ispec,:)) / h_m(ix,iy,iz)
          
            ! Calculate updated fractions after re-emission
            if ( LE_conc(ix,iy,iz,SA_ispec) < -1e-4 ) then
              write( gol, '("LE_conc became negative")' ) ; call GoErr
              write( gol, '(i0," ",i0)' ) ix,iy; call GoErr
              write( gol, '(e12.5)') budget_depos(ix,iy,ispec,1) ; call GoErr
              write( gol, '(e12.5)') budget_depos(ix,iy,ispec,2) ; call GoErr
              write( gol, '(e12.5)') LE_conc(ix,iy,iz,SA_ispec) ; call GoErr
              TRACEBACK;status=1;return
            else if ( abs(LE_conc(ix,iy,iz,SA_ispec)) < 1e-4 ) then
              SA_frac(ix,iy,iz,SA_ispec,:) = 0.0
            else 
              do SA_ilabel = 1, SA_nlabel
                SA_frac(ix,iy,iz,SA_ispec,SA_ilabel) = ( SA_conc_before_dryemis(ix,iy,SA_ilabel) + &
                                                         budget_depos(ix,iy,ispec,iexdry_emis)/ h_m(ix,iy,iz) * SA_frac_comp_point(ix,iy,SA_ispec, SA_ilabel)  ) /&
                                                           LE_conc(ix,iy,iz,SA_ispec)
              end do
            end if  !
          else 
            ! only update concentrations in lowest layer  ( divide with height to get ug/m3 or ppb)
            LE_conc(ix,iy,iz,SA_ispec) = LE_conc(ix,iy,iz,SA_ispec) + sum( budget_depos(ix,iy,ispec,:)) / h_m(ix,iy,iz)

          end if ! Re-emission?

        end do  ! iy
        end do  ! ix
          
      end if ! Labelled tracer?
      
    end do ! loop over tracers
        
    ! ok
    status = 0
  
  end subroutine SA_DryDepos_Fractions

! ----------------------------------------------

  subroutine SA_DryDepos_Reset( status)
  
    ! --- in/out --------------------
    integer, intent(out)      ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_DryDepos_Reset'
    
    ! --- begin ---------------------
        
    if ( allocated(SA_budget_drydepos_hour) ) then    
      SA_budget_drydepos_hour = 0.0
    end if
    
    ! ok
    status = 0
  
  end subroutine SA_DryDepos_Reset
  
! ----------------------------------------------

  subroutine SA_DryDepos_Fractions_lu_day( budget_depos, status)
    
    use Binas, only : xm_air
    use Indices, only : specunit, specname
    use LE_Landuse_Data, only : nlu
    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out --------------------
    real, intent(in)          ::  budget_depos(nx,ny,nspec,nlu)
    integer, intent(out)      ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_DryDepos_Fractions_lu_day'
    
    ! --- local ---------------------
    integer                   ::  ispec, iz, ilu, ilu_budget
    integer                   ::  SA_ilabel
    integer                   ::  SA_ispec, SA_ispec_budget
    real                      ::  conv_fact(nx,ny)
    ! meteo data:
    real, pointer             ::  dens(:,:,:)   ! (lon,lat,alt)    
    
    ! --- begin ---------------------
            
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOTOK_RETURN(status=1)

    ! surface for dry deposition
    iz = 1
    
    do ispec = 1, nspec     
      ! Search corresponding specie parameter in list of labelled species
      SA_ispec = Labelled_specs(ispec)
      SA_ispec_budget = Labelled_specs_budget_output(ispec)
      
      if (SA_ispec > 0 .and. SA_ispec_budget > 0 ) then

        select case ( trim(specunit(ispec)) )
          !~ other notations ...
          case ( 'ug/m3' )
            ! Aerosol tracer, no conversion needed
            conv_fact = 1.0
          case ( 'ppb' )
            ! gaseous tracer (conversion to ug/m2 needed)
            !           ((mole tr/mole air)/ppb) * (kg tr/mole tr)   * (kg air/m3 air) / (kg air/mole air) * ug/kg
            conv_fact =       1e-9               *  specmolm(ispec)  * dens(1:nx,1:ny,iz)  /     xm_air        * 1e9
          !~ unkown ...
          case default
            write (gol,'("do not know how to convert `",a,"` from `",a," m` to `ug m-2")') &
                    trim(specname(ispec)), trim(specunit(ispec)); call goErr
            TRACEBACK; status=1; return
        end select
        
        ! Add deposition budgets to the daily deposition budget        
        do ilu = 1, nlu
          ilu_budget = landuse_budget_output(ilu)
          if ( ilu_budget > 0 ) then
            do SA_ilabel = 1, SA_nlabel
              SA_budget_drydepos_lu_day(:,:,SA_ispec_budget,ilu_budget,SA_ilabel) = &
                  SA_budget_drydepos_lu_day(:,:,SA_ispec_budget,ilu_budget,SA_ilabel) + &
                        budget_depos(:,:,ispec,ilu) * conv_fact * SA_frac(:,:,1,SA_ispec,SA_ilabel)
            end do
          end if  
        end do ! landuses
        
      end if ! Labelled tracer?

    end do ! loop over tracers
        
    ! ok
    status = 0
  
  end subroutine SA_DryDepos_Fractions_lu_day  
  
  ! ----------------------------------------------
 
  subroutine SA_Comp_point_Init( exdryemis, exdrydepo, status )
    
    integer, intent(in)   ::  exdryemis, exdrydepo
    integer, intent(out)  ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Comp_Point_Init'
    
    ! --- begin ---------------------
    
    iexdry_emis = exdryemis
    iexdry_depo = exdrydepo
        
    ! ok
    status = 0
    
    
  end subroutine SA_Comp_point_Init
  
  ! ----------------------------------------------
 
  subroutine SA_Comp_point_Done( status )
    
    integer, intent(out)  ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Comp_Point_Done'
    
    ! --- begin ---------------------
    
    ! done
    deallocate( SA_frac_comp_point, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( SA_conc_comp_point, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine SA_Comp_point_Done
  
  
  ! ----------------------------------------------

  subroutine SA_Comp_point_budget( renew, status, conc_sfc )
  
    use Indices      , only : i_so2, i_nh3
    
    ! --- in/out --------------------
    
    logical, intent(in)         ::  renew
    integer, intent(out)        ::  status
    real, intent(in), optional  ::  conc_sfc(nx,ny,nspec)
    
    ! --- const -------------------------
    
    character(len=*), parameter :: rname = mname//'/SA_Comp_Point_Budget'
    
    ! --- local ---------------------
    
    integer                     ::  ix,iy, ispec
    integer                     ::  SA_ispec, SA_ilabel
    real                        ::  inv_denum
    
    ! --- begin ---------------------
        
    if (present( conc_sfc) ) then
      do ispec = 1, nspec
        ! Only budgets for NH3 and SO2 are stored and used
        if ( .not. ( ispec == i_nh3 .or. ispec == i_so2 ) ) cycle

        ! Tracer labelled?
        SA_ispec = Labelled_specs(ispec)         
        if ( SA_ispec > 0 ) then
          do SA_ilabel = 1, SA_nlabel
            SA_conc_comp_point(:,:,SA_ispec,SA_ilabel) = SA_conc_comp_point(:,:,SA_ispec,SA_ilabel)  + conc_sfc(:,:,ispec) * SA_frac(:,:,1,SA_ispec,SA_ilabel)
          end do
        end if
      end do
    end if
    
    ! if end of the month, calculate fractional contribution of past month and reset budget to zero
    if ( renew ) then
      
      do ispec = 1, nspec 
        ! Only budgets for NH3 and SO2 are stored and used
        if ( .not. ( ispec == i_nh3 .or. ispec == i_so2 ) ) cycle
        
        ! Tracer labelled?
        SA_ispec = Labelled_specs(ispec)
        if ( SA_ispec > 0 ) then
          
          ! Calculate fractional contribution of compensation point
          do ix = 1, nx
          do iy = 1, ny
            if ( sum(SA_conc_comp_point(ix,iy,SA_ispec,:)) > 0 ) then
              ! Check division by zero
              inv_denum = 1.0/sum(SA_conc_comp_point(ix,iy,SA_ispec,:))
              do SA_ilabel = 1, SA_nlabel
                SA_frac_comp_point(ix,iy,SA_ispec,SA_ilabel) = SA_conc_comp_point(ix,iy,SA_ispec,SA_ilabel) * inv_denum
              end do
            else
              SA_frac_comp_point(ix,iy,SA_ispec,:) = 0.0
            end if
          end do
          end do    
        end if
      end do
      
      ! reset budget
      SA_conc_comp_point = 0.0
      
    end if
    
    ! ok
    status = 0
  
  end subroutine SA_Comp_point_budget
  
! ----------------------------------------------

  subroutine SA_WetDepos_Fractions( budget_depos, status)
    
    use Binas, only : xm_air
    use Indices, only : specunit, specname
    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out --------------------
    real, intent(in)          ::  budget_depos(nx,ny,nz,nspec)
    integer, intent(out)      ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_WetDepos_Fractions'
    
    ! --- local ---------------------
    integer                   ::  ispec
    integer                   ::  SA_ilabel
    integer                   ::  SA_ispec, SA_ispec_budget, SA_ispec_wetdep 
    real                      ::  conv_fact(nx,ny,nz)
    ! meteo data:
    real, pointer             ::  dens(:,:,:)   ! (lon,lat,alt)    
    
    ! --- begin ---------------------
    
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOTOK_RETURN(status=1)

    do ispec = 1, nspec     
      ! Search corresponding specie parameter in list of labelled species
      SA_ispec = Labelled_specs(ispec)
      SA_ispec_budget = Labelled_specs_budget_output(ispec)
      SA_ispec_wetdep = Labelled_specs_wetdep_output(ispec)

      if (SA_ispec > 0 .and. ( SA_ispec_budget > 0 .or. SA_ispec_wetdep > 0 ) ) then
        
        select case ( trim(specunit(ispec)) )
          !~ other notations ...
          case ( 'ug/m3' )
            ! aerosol tracer (no conversion needed)
            conv_fact = 1.0
          case ( 'ppb' )
            ! gaseous tracer (conversion to ug/m2 needed)
            ! Conversion from ppb m to ug/m2
            !                  ((mole tr/mole air)/ppb) * (kg tr/mole tr)   * (kg air/m3 air) / (kg air/mole air) * ug/kg
            conv_fact(:,:,:) =       1e-9               *  specmolm(ispec)  * dens(1:nx,1:ny,1:nz)    /     xm_air        * 1e9
          !~ unkown ...
          case default
            write (gol,'("do not know how to convert `",a,"` from `",a," m` to `ug m-2")') &
                    trim(specname(ispec)), trim(specunit(ispec)); call goErr
            TRACEBACK; status=1; return
        end select
        
        ! Update budgets on daily basis
        if ( SA_ispec_budget > 0 ) then     
            do SA_ilabel = 1, SA_nlabel 
              SA_budget_wetdepos_day(:,:,SA_ispec_budget,SA_ilabel)  = SA_budget_wetdepos_day(:,:,SA_ispec_budget,SA_ilabel) + &
                                                                sum( budget_depos(:,:,:,ispec) * conv_fact * SA_frac(:,:,:,SA_ispec,SA_ilabel), dim=3)
            end do                         
        end if

        ! Update budgets on hourly basis
        if ( SA_ispec_wetdep > 0 ) then     
            do SA_ilabel = 1, SA_nlabel 
              SA_budget_wetdepos_hour(:,:,SA_ispec_wetdep,SA_ilabel) = SA_budget_wetdepos_hour(:,:,SA_ispec_wetdep,SA_ilabel) + &
                                                                sum( budget_depos(:,:,:,ispec) * conv_fact * SA_frac(:,:,:,SA_ispec,SA_ilabel), dim=3)
            end do 
        end if


      end if ! Labelled tracer?
    end do ! loop over tracers
        
    ! ok
    status = 0
  
  end subroutine SA_WetDepos_Fractions

! ----------------------------------------------

  subroutine SA_WetDepos_Reset( status)
  
    ! --- in/out --------------------
    integer, intent(out)      ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_WetDepos_Reset'
    
    ! --- begin ---------------------
        
    if ( allocated( SA_budget_wetdepos_hour) ) then    
      SA_budget_wetdepos_hour = 0.0
    end if
    
    ! ok
    status = 0
  
  end subroutine SA_WetDepos_Reset


! ----------------------------------------------


  subroutine SA_Depos_Reset_day( status)
  
    ! --- in/out --------------------
    integer, intent(out)      ::  status
    
    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_DryDepos_Reset_lu_day'
    
    ! --- begin ---------------------
        
    SA_budget_drydepos_lu_day = 0.0    
    SA_budget_drydepos_day = 0.0
    SA_budget_wetdepos_day = 0.0
    
    ! ok
    status = 0
  
  end subroutine SA_Depos_Reset_day


!
! ---------------- CHEMISTRY ROUTINES -----------------------------------------------------
!

  subroutine SA_Chem_Done(status)

    ! --- in/out ------------------------
    integer, intent(out)      ::  status

    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Chem_Done'

    ! --- begin -------------------------

    deallocate( array_rows, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( array_cols, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Chem_Done

! ---------------------------------------------------------------------------------------

  subroutine SA_Chem_Gas_Setup(array_col_indices, array_row_indices, n_nonzeros, status)

    ! --- in/out ------------------
    integer, intent(in)       ::  array_col_indices(n_nonzeros)
    integer, intent(in)       ::  array_row_indices(n_nonzeros)
    integer, intent(in)       ::  n_nonzeros
    integer, intent(out)      ::  status

    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Chem_Gas_Setup'

    ! --- local -------------------------
    integer                   ::  i_nonzero
    integer, allocatable      ::  array_rowsh(:)
    integer, allocatable      ::  array_colsh(:)
    ! --- begin -------------------------

    ! allocate indices of non-zero elements
    allocate(array_cols(n_nonzeros), stat=status)
    IF_NOTOK_RETURN(status=1)
    allocate(array_rows(n_nonzeros), stat=status)
    IF_NOTOK_RETURN(status=1)
    allocate(Labelled_Chem_specs(nspec), stat=status)
    IF_NOTOK_RETURN(status=1)


    ! set column indices (order of species is still order in LE not in SA)
    array_rows = array_row_indices
    array_cols = array_col_indices

    ! Find the tracers which are labelled and involved in the chemistry scheme
    Labelled_Chem_Specs = -1
    do i_nonzero = 1, n_nonzeros
      Labelled_Chem_Specs(array_rows(i_nonzero)) = 1
      Labelled_Chem_Specs(array_cols(i_nonzero)) = 1
    enddo

    ! state number of nonzeros
    SA_n_nonzeros = n_nonzeros

    ! set the order of row and column indices from LE order to SA order)
    allocate( array_rowsh(SA_n_nonzeros) )
    allocate( array_colsh(SA_n_nonzeros) )
    array_colsh = array_cols
    array_rowsh = array_rows

    do i_nonzero = 1, SA_n_nonzeros
      array_cols(i_nonzero) = Labelled_specs(array_colsh(i_nonzero))
      array_rows(i_nonzero) = Labelled_specs(array_rowsh(i_nonzero))
    enddo

    ! clear
    deallocate(array_rowsh)
    deallocate(array_colsh)

    ! ok
    status = 0

  end subroutine SA_Chem_Gas_Setup

! ---------------------------------------------------------------------------------------

  subroutine SA_Chem_Gas_Iter( array_prod, array_prod_total, array_prod_total_backup, dt, dt_total, dt_total_backup, keyword)

    ! --- in/out ------------------
    real, intent(in)        ::  array_prod(SA_n_nonzeros)
    real, intent(inout)     ::  array_prod_total(SA_n_nonzeros)
    real, intent(inout)     ::  array_prod_total_backup(SA_n_nonzeros)
    real, intent(in)        ::  dt
    real, intent(inout)     ::  dt_total
    real, intent(inout)     ::  dt_total_backup
    character(len=*)        ::  keyword

    ! --- const -------------------
    character(len=*), parameter     :: rname = mname//'/SA_Chem_Gas_Iter'

    ! --- begin -------------------

    if (keyword == 'Backup' ) then
      array_prod_total_backup = array_prod_total
      dt_total_backup = dt_total

      array_prod_total = array_prod_total + array_prod * dt
      dt_total = dt_total + dt

    else if (keyword == 'failer' ) then
      array_prod_total = array_prod_total_backup
      dt_total = dt_total_backup
    else
       write (gol, '("incorrect keyword")' ) ; call goErr
       write (gol, '("keyword: ", a)' ) keyword; call goErr
       write (gol, '("keyword must be Backup or failer " )' ) ; call goErr
       TRACEBACK; return
    endif

  end subroutine SA_Chem_Gas_Iter


! ---------------------------------------------------------------------------------------


  subroutine SA_Chem_Gas_Step( ix,iy,iz, y_new, array_prod_total, status)

#ifdef with_mkl
    use MKL, only : sp_scoomm => mkl_scoomm
    use MKL, only : sp_dcoomm => mkl_dcoomm
#endif
#ifdef with_spblas
    use spblas, only : sp_scoomm, sp_dcoomm
#endif
    use dims, only : nspec

    ! --- in/out ------------------

    integer, intent(in)     ::  ix, iy, iz
    real, intent(in)        ::  y_new(nspec)
    real, intent(in)        ::  array_prod_total(SA_n_nonzeros)
    integer, intent(out)    ::  status

    ! --- const -------------------------

    character(len=*), parameter     :: rname = mname//'/SA_Chem_Gas_Step'

    ! description for matrix multiplication
    character(len=1), dimension(4)  :: mat_descra = (/'G', 'L', 'N', 'F'/)

    ! --- local -------------------------

    real                    ::  inv_denum
    integer                 ::  ispec, SA_ispec
    real                    ::  prod_plus_old
    real                    ::  sum_matprod
    real                    ::  matprod(nspec_labelled, SA_nlabel)


    ! --- begin -------------------------


    ! calculate the produced specie per label, via a sparse matrix mulitplication
    ! sp_scoomm : Sparse COOrdinate Matrix Multiplication: C = alpha *  A * B + beta * C
    ! in this case this will be                          : matprod  = 1.0 * productionmatrix * Label_fractions + 0 * matprod
    ! input : 'N'            : Normal Matrix multiplication
    !         nspec_labelled : # rows of A
    !         SA_nlabel      : # columns of C
    !         nspec_labelled : # columns of A
    !         1.0            : alpha = 1.0
    !         mat_descra     : description of matrix multiplication see : http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/cpp/win/mkl/refman/bla/bla_IC.html#tbl2-6
    !         array_prods    : array with nonzero elements of A
    !         array_rows     : array with row indices of nonzero elements of A
    !         array_cols     : array with column indices of nonzero elements of A
    !         SA_n_nonzeros  : Number of non-zero elements of A
    !         SA_frac        : Matrix B with old label fractions
    !         nspec_labelled : # rows of B
    !         0.0            : beta = 0.0
    !         nspec_labelled : # rows of C

    ! output: matprod        : Output of this routine with the production for each specie specified per label
    matprod = 0.0
#if defined (with_mkl) || defined (with_spblas)
    if ( kind(matprod) == 4 ) then
      call sp_scoomm('N', nspec_labelled, SA_nlabel, nspec_labelled, 1.0, mat_descra, array_prod_total, &
                      array_rows, array_cols, SA_n_nonzeros, SA_frac(ix,iy,iz,:,:), nspec_labelled, &
                      0.0, matprod, nspec_labelled )    
    else
      write( gol, '("Unsupported data kind for sparse blas matrix multiplication")' ) ; call GoErr
      TRACEBACK;status=1;return
    end if
#else
    write (gol,'("labelling needs MKL (ifort) or SPBLAS (gfortran) library, define either `with_mkl` or `with_spblas` ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    do ispec = 1, nspec
      ! tracer reactive ? and Labelled ?
      if (Labelled_Chem_Specs(ispec) > 0 ) then
        
        ! find corresponding index number in Labelled species
        SA_ispec = Labelled_specs(ispec)

        ! calculate production plus old concentration
        sum_matprod = sum(matprod(SA_ispec,:))

        prod_plus_old = LE_conc(ix,iy,iz,SA_ispec) + sum_matprod
        
        ! new concentration larger than production and old concentration ?
        ! theoretically not possible, but die to linearizations and lumping time steps it can happen
        ! Rescale production (equally over labels) such that production + old conc = new_conc
        ! loss is now stated equal to zero
        if ( y_new(ispec) > prod_plus_old .and. sum_matprod > numlim_conc) then
          matprod(SA_ispec,:) = ( (y_new(ispec) - LE_conc(ix,iy,iz,SA_ispec) )/ sum_matprod ) * matprod(SA_ispec,:)
          prod_plus_old = y_new(ispec)
        endif

        ! Calculate new fractions (if y_new == 0 the fractions do not matter)
        ! Fractions are calculated from the old concentration and the production,
        ! the loss does not change the labels. In this application it is assumed that for one time step
        ! first the production is taken and after that the loss.
        ! new labels are calculated with production and old conc.
        ! loss after this calculation does not change fractions.
        if (y_new(ispec) > numlim_conc .and. prod_plus_old > numlim_conc) then
          inv_denum = 1.0 / prod_plus_old
          SA_frac(ix,iy,iz,SA_ispec,:) = ( matprod(SA_ispec,:) + SA_frac(ix,iy,iz,SA_ispec,:) * (LE_conc(ix,iy,iz,SA_ispec) )  )* inv_denum
        else
          SA_frac(ix,iy,iz,SA_ispec,:) = 0.0
        endif

        ! synchronize the help concentration to the LE conc
        LE_conc(ix,iy,iz,SA_ispec) = y_new(ispec)
      endif
    enddo

    ! ok
    status = 0

  end subroutine SA_Chem_Gas_Step

! ---------------------------------------------------------------------------------------

  subroutine SA_Chem_Cloud(c,status)

    use indices

    ! --- in/out ------------------
    real, intent(in)        ::  c(nx,ny,nz,nspec)
    integer, intent(out)    ::  status

    ! --- const -------------------
    character(len=*), parameter     :: rname = mname//'/SA_Chem_Cloud_Setup'

    ! --- local -------------------
    real                    ::  inv_denum
    integer                 ::  SA_i_so4a_f, SA_i_so2
    integer                 ::  ix,iy,iz

    ! --- begin -------------------

    ! Cloud chemistry only changes so2, so4a_f, h2o2 and o3 concentrations
    ! h2o2 and o3 are not labelled, and for so2 there is only loss, thus no change of fractions
    ! all produced so4a_f is set with the labels from so2
    ! F_new_so4 = ( (F_old_so4 * c_old_so4) + F_old_so2 * (c_new_so4 - c_old_so4) )/ c_new_so4
    SA_i_so4a_f = Labelled_specs(i_so4a_f)
    SA_i_so2    = Labelled_specs(i_so2)

    ! sulphurs labelled ?
    if (SA_i_so4a_f > 0) then
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx
        if (c(ix,iy,iz,i_so4a_f) > numlim_conc ) then
          inv_denum = 1/c(ix,iy,iz,i_so4a_f)
          SA_frac(ix,iy,iz,SA_i_so4a_f,:) = ( SA_frac(ix,iy,iz,SA_i_so4a_f,:) *  LE_conc(ix,iy,iz,SA_i_so4a_f) + &
                                              SA_frac(ix,iy,iz,SA_i_so2,   :) * (c(ix,iy,iz,i_so4a_f) - LE_conc(ix,iy,iz,SA_i_so4a_f)) ) *inv_denum
        endif
      enddo
      enddo
      enddo

      ! State the new so2 and so4a concentration
      LE_conc(:,:,:,SA_i_so4a_f) = c(:,:,:,i_so4a_f)
      LE_conc(:,:,:,SA_i_so2   ) = c(:,:,:,i_so2)

      ! ok
      status = 0

    endif

  end subroutine SA_Chem_Cloud

! ---------------------------------------------------------------------------------------

  subroutine SA_Chem_SIA(status)

    use indices

    ! --- in/out ------------------
    integer, intent(out)    ::  status

    ! --- const -------------------
    character(len=*), parameter     :: rname = mname//'/SA_Chem_SIA'

    ! --- local -------------------
    real                    ::  mw_NH3 = 17.            ! molecular weight of NH3
    real                    ::  mw_HNO3 = 63.           ! molecular weight of HNO3
    real                    ::  inv_mol_volume = 0.0409 ! inverse of the volume of 1 mol gas
    real                    ::  inv_denum_ammo
    real                    ::  inv_denum_nitr
    integer                 ::  ix, iy, iz

    ! --- begin -------------------

    ! SIA tracers enabled ?
    if ( (i_nh4a_f > 0) .and. (i_nh3 > 0) .and. (i_so2 > 0) .and. (i_so4a_f > 0) ) then    
      ! SIA labels enabled ?
      if ( Labelled_specs(i_nh3) > 0 .and. Labelled_specs(i_nh4a_f ) > 0.and. Labelled_specs(i_hno3) > 0 .and. Labelled_specs(i_no3a_f) > 0 ) then

        !$OMP parallel &
        !$OMP default( none ) &
        !$OMP shared( nx,ny,nz ) &
        !$OMP shared( LE_conc ) &
        !$OMP shared( SA_frac ) &
        !$OMP shared( Labelled_specs ) &
        !$OMP shared( i_nh4a_f, i_nh3, i_no3a_f, i_hno3 ) &
        !$OMP shared( mw_nh3, mw_hno3, inv_mol_volume ) &
        !$OMP private( ix,iy,iz) &
        !$OMP private( inv_denum_ammo) &
        !$OMP private( inv_denum_nitr)
        !$OMP do

        ! Only nh3, nh4a_f, no3a_f and hno3 are updated by this routine, denumerator is the sum of concentrations (nh4a_f + nh3) and (no3a_f and hno3)
        ! Due to fast equilibrium the labels of (NH3 and NH4a_f) and (HNO3 and NO3a_f) are assumed the same as a weighted average
        do iz = 1, nz
        do iy = 1, ny
        do ix = 1, nx

          if (LE_conc(ix,iy,iz,Labelled_specs(i_NH3) ) * mw_NH3 * inv_mol_volume + LE_conc(ix,iy,iz,Labelled_specs(i_nh4a_f )) > numlim_conc ) then

            inv_denum_ammo = 1.0 / ( LE_conc(ix,iy,iz,Labelled_specs(i_NH3 )) * mw_NH3 * inv_mol_volume + LE_conc(ix,iy,iz,Labelled_specs(i_NH4a_f )) )

            SA_Frac(ix,iy,iz,Labelled_specs(i_NH3),:) = &
              (LE_conc(ix,iy,iz,Labelled_specs(i_NH4a_f)) * SA_Frac(ix,iy,iz,Labelled_specs(i_NH4a_f),:) + &
               LE_conc(ix,iy,iz,Labelled_specs(i_NH3   )) * SA_Frac(ix,iy,iz,Labelled_specs(i_NH3)   ,:) * mw_NH3 * inv_mol_volume) * inv_denum_ammo

            SA_Frac(ix,iy,iz,Labelled_specs(i_NH4a_f),:) = SA_Frac(ix,iy,iz,Labelled_specs(i_NH3),:)

          else
            SA_frac(ix,iy,iz,Labelled_specs(i_NH4a_f),:) = 0.0
            SA_frac(ix,iy,iz,Labelled_specs(i_NH3)   ,:) = 0.0
          endif

          if (LE_conc(ix,iy,iz,Labelled_specs(i_HNO3)) * mw_HNO3 * inv_mol_volume + LE_conc(ix,iy,iz,Labelled_specs(i_NO3a_f) ) > numlim_conc ) then
            inv_denum_nitr = 1.0 / ( LE_conc(ix,iy,iz,labelled_specs(i_HNO3)) * mw_HNO3 * inv_mol_volume + LE_conc(ix,iy,iz,Labelled_specs(i_NO3a_f) ) )

            SA_Frac(ix,iy,iz,Labelled_specs(i_HNO3),:) = &
                (LE_conc(ix,iy,iz,Labelled_specs(i_NO3a_f)) * SA_Frac(ix,iy,iz,Labelled_specs(i_NO3a_f),:) + &
                 LE_conc(ix,iy,iz,Labelled_specs(i_HNO3)  ) * SA_Frac(ix,iy,iz,Labelled_specs(i_HNO3)  ,:) * mw_HNO3 * inv_mol_volume) * inv_denum_nitr

            SA_Frac(ix,iy,iz,Labelled_specs(i_NO3a_f),:) = SA_Frac(ix,iy,iz,Labelled_specs(i_HNO3),:)
          else
            SA_Frac(ix,iy,iz,Labelled_specs(i_NO3a_f),:) = 0.0
            SA_Frac(ix,iy,iz,Labelled_specs(i_HNO3)  ,:) = 0.0

          endif
        enddo
        enddo
        enddo
        !$OMP end do
        !$OMP end parallel
      end if  ! sia tracers labeled ?
    end if  ! sia tracers enabled ?
    
    ! ok
    status = 0

  end subroutine SA_Chem_SIA

!
! ---------------- RESTART ROUTINES -----------------------------------------------------
!

  subroutine SA_Restart_Restore(c,labels,comp_point_label,comp_point_label_frac,status)

    ! --- in/out ------------------
    real, intent(in)          ::  c(nx,ny,nz,nspec)
    real, intent(in)          ::  labels(nx,ny,nz,nspec_labelled,SA_nlabel)
    real, intent(in)          ::  comp_point_label(nx,ny,nspec_labelled,SA_nlabel)
    real, intent(in)          ::  comp_point_label_frac(nx,ny,nspec_labelled,SA_nlabel)
    integer, intent(out)      ::  status

    ! --- const -------------------
    character(len=*), parameter :: rname = mname//'/SA_Restart_Restore'

    ! --- local -------------------
    integer                   ::  ispec, SA_ispec

    ! --- begin -------------------
    do ispec = 1, nspec
      SA_ispec = Labelled_specs(ispec)
      if (SA_ispec > 0 ) then
        LE_conc(:,:,:,SA_ispec) = c(:,:,:,ispec)
      end if
    end do
    
    ! Fill in fractions and budgets for compensation point
    SA_frac = labels    
    SA_conc_comp_point = comp_point_label
    SA_frac_comp_point = comp_point_label_frac

    ! ok
    status = 0

  end subroutine SA_Restart_Restore

!
! ---------------- BOUNDARY ROUTINES ----------------------------------------------------
!

  subroutine SA_Bound_Setup(t, status)

    use GO     , only : Get, wrtgol
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get, ExpandTime
    use GO     , only : Precisely, MidNight, operator(-), operator(+), operator(*),&
                        operator(/=), operator(==), wrtgol
    use NetCDF , only : NF90_Open, NF90_Close
    use NetCDF , only : NF90_NOWRITE, NF90_ENOTVAR
    use NetCDF , only : NF90_Inq_Dimid, NF90_Inquire_Dimension
    use NetCDF , only : NF90_Get_Var, NF90_Inq_Varid, NF90_Get_Att

    use C3PO   , only : T_Grid_Ugg    
    use LE_Grid, only : ugg
    use LE_Grid, only : ugg_west, ugg_east, ugg_south, ugg_north

    ! --- in/out ------------------
    type(Tdate), intent(in)   ::  t
    integer, intent(out)      ::  status

    ! --- const -------------------
    character(len=*), parameter  :: rname = mname//'/SA_Bound_Setup'

    ! --- local -------------------
    integer                   ::  year, month, day
    character(len=512)        ::  fname
    integer                   ::  ncid, dimid, varid
    real, allocatable         ::  lons(:), lats(:)
    real, allocatable         ::  lons_bnds(:,:), lats_bnds(:,:)
    real, allocatable         ::  label_bounds(:,:,:,:)
    type(T_Grid_Ugg)          ::  ugg_bc
    real*8, allocatable       ::  times(:)
    type(TDate)               ::  t_in
    type(TDate)               ::  time_ref
    type(TIncrDate)           ::  time_step
    logical                   ::  found
    integer                   ::  it, nt
    integer                   ::  irec
    character(len=32)         ::  vname
    integer                   ::  ilon, ilat
    integer                   ::  ispec
    integer                   ::  iz
    character(len=32)         ::  unit
    

    ! --- begin -------------------

    ! extract time
    call Get(t, year = year, month = month, day = day)

    ! Boundary file name
    write (fname, '(a,"/",a,2("_",a),"_",i4.4,2i2.2,".nc")') &
            trim(arch_path), trim(arch_model), trim(arch_expid), trim(arch_name), &
            year, month, day

    ! open file (read-only):
    status = NF90_Open( trim(fname), NF90_NOWRITE, ncid )
    if (status/=NF90_NOERR) then
      write (gol,'("NF90 error: ",a)') trim(nf90_strerror(status)); call goErr
      write (gol,'("opening file : ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    endif

    ! setup boundary grid?
    if (.not. cache_setup ) then

      status = NF90_Inq_DimID( ncid, 'longitude', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Inquire_Dimension( ncid, dimid, len=cache_nlon )
      IF_NF90_NOTOK_RETURN(status=1)
      ! get dimension:
      status = NF90_Inq_DimID( ncid, 'latitude', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Inquire_Dimension( ncid, dimid, len=cache_nlat )
      IF_NF90_NOTOK_RETURN(status=1)

      ! storage:
      allocate( lons(cache_nlon), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( lats(cache_nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( lons_bnds(2,cache_nlon), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( lats_bnds(2,cache_nlat), stat=status )
      IF_NOTOK_RETURN(status=1)

      ! get axis:
      status = NF90_Inq_VarID( ncid, 'longitude', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, lons )
      IF_NF90_NOTOK_RETURN(status=1)
      ! boundary axis
      status = NF90_Inq_VarID( ncid, 'longitude_bnds', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, lons_bnds )
      IF_NF90_NOTOK_RETURN(status=1)
      ! get axis:
      status = NF90_Inq_VarID( ncid, 'latitude', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, lats )
      IF_NF90_NOTOK_RETURN(status=1)
      ! boundary axis
      status = NF90_Inq_VarID( ncid, 'latitude_bnds', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, lats_bnds )
      IF_NF90_NOTOK_RETURN(status=1)

      ! setup grid:
      call ugg_bc%Init( lons(:), lons_bnds(:,:), &
                        lats(:), lats_bnds(:,:), status )
      IF_NF90_NOTOK_RETURN(status=1)
     
      ! clear:
      deallocate( lons, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( lats, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( lons_bnds, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( lats_bnds, stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! simply assign concentrations from a source cell to a boundary condition cell;
      ! select and store source cell indices:
      allocate( cache_ij_west(2,ugg_west%nlon,ugg_west%nlat), stat = status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg_west%nlat
        do ilon = 1, ugg_west%nlon
          call ugg_bc%GetLocation( ugg_west%longitude_1d(ilon), ugg_west%latitude_1d(ilat), &
                              cache_ij_west(1,ilon,ilat), cache_ij_west(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do
      allocate( cache_ij_east(2,ugg_east%nlon,ugg_east%nlat), stat = status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg_east%nlat
        do ilon = 1, ugg_east%nlon
          call ugg_bc%GetLocation( ugg_east%longitude_1d(ilon), ugg_east%latitude_1d(ilat), &
                              cache_ij_east(1,ilon,ilat), cache_ij_east(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do
      allocate( cache_ij_south(2,ugg_south%nlon,ugg_south%nlat), stat = status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg_south%nlat
        do ilon = 1, ugg_south%nlon
          call ugg_bc%GetLocation( ugg_south%longitude_1d(ilon), ugg_south%latitude_1d(ilat), &
                              cache_ij_south(1,ilon,ilat), cache_ij_south(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do
      allocate( cache_ij_north(2,ugg_north%nlon,ugg_north%nlat), stat = status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg_north%nlat
        do ilon = 1, ugg_north%nlon
          call ugg_bc%GetLocation( ugg_north%longitude_1d(ilon), ugg_north%latitude_1d(ilat), &
                              cache_ij_north(1,ilon,ilat), cache_ij_north(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do

      ! get dimension:
      status = NF90_Inq_DimID( ncid, 'level', dimid )
      if ( status /= NF90_NOERR ) then
        status = NF90_Inq_VarID( ncid, 'lev', varid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
      status = NF90_Inquire_Dimension( ncid, dimid, len=cache_nlev )
      IF_NF90_NOTOK_RETURN(status=1)

      ! check ...
      if ( cache_nlev /= nz ) then
        write (gol,'("number of levels in label nc output does not match with current model:")'); call goErr
        write (gol,'("  nc file levels : ",i6)') cache_nlev; call goErr
        write (gol,'("  model levels   : ",i6)') nz; call goErr
        !TRACEBACK; status=1; return
      end if

      ! reset flag:
      cache_setup = .true.

    end if   ! setup cache

    ! get time dimension:
    status = NF90_Inq_DimID( ncid, 'time', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=nt )
    IF_NF90_NOTOK_RETURN(status=1)
    ! storage for dates:
    allocate( times(nt) )
    ! time variable:
    status = NF90_Inq_VarID( ncid, 'time', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, times )
    IF_NF90_NOTOK_RETURN(status=1)
    ! time units:
    status = NF90_Get_Att( ncid, varid, 'units', unit )
    IF_NF90_NOTOK_RETURN(status=1)
    ! search:
    found = .false.
    do it = 1, nt
      ! expand:
      call ExpandTime( real(times(it)), unit, 'standard', t_in, status )
      IF_NOTOK_RETURN(status=1)
      ! compare:
      if ( t == t_in ) then
        irec = it
        found = .true.
      end if
    end do
    call wrtgol( '  model time      : ', t )
    write (gol,'("  output times    ")')
    if ( .not. found ) then
      write (gol,'("times in label nc output do not match with requested time:")'); call goErr
      call wrtgol( '  model time      : ', t ); call goErr
      write (gol,'("  output times    ")'); call goErr
      do it = 1, nt
        call ExpandTime( real(times(it)), unit, 'standard', t_in, status )
        IF_NOTOK_RETURN(status=1)
        call wrtgol( '     ', t_in ); call goPr
      end do
      TRACEBACK; status=1; return
    end if
    ! clear:
    deallocate( times )


    do ispec = 1, nspec_labelled

      ! storage for single field:
      allocate( label_bounds(cache_nlon,cache_nlat,cache_nlev, SA_nlabel), stat = status )
      IF_NOTOK_RETURN(status=1)

      ! variable name
      vname = trim(labelled_specs_names(ispec))

      status = NF90_Inq_VarID( ncid, trim(vname), varid )

      if ( status == NF90_ENOTVAR ) then ! not a variable in nc-file
        ! info ...
        write (gol,'("        label nc field not found : ",a)') trim(vname); call goPr
        ! not available; try next:
        cycle
      end if
      ! other errors
      IF_NF90_NOTOK_RETURN(status=1)

      ! info ..
      write (gol,'("        read Label nc field : ",a)') trim(vname); call goPr

      ! read record:
      status = NF90_Get_Var( ncid, varid, label_bounds, start=(/1,1,1,1,irec/), &
                       count=(/cache_nlon,cache_nlat,cache_nlev,SA_nlabel,1/) )
      IF_NF90_NOTOK_RETURN(status=1)

      ! assume levels are the same; loop over levels:
      do iz = 1, nz
        ! copy labelbounds from coarse grid to boundary conditions for the known species from the boundary run
        do ilat = 1, ugg_west%nlat
          do ilon = 1, ugg_west%nlon
            SA_frac_adv(:,0,ilat,iz,ispec)    = label_bounds(cache_ij_west( 1,ilon,ilat),cache_ij_west( 2,ilon,ilat),iz,:)
          end do
        end do
        do ilat = 1, ugg_east%nlat
          do ilon = 1, ugg_east%nlon
            SA_frac_adv(:,nx+1,ilat,iz,ispec) = label_bounds(cache_ij_east( 1,ilon,ilat),cache_ij_east( 2,ilon,ilat),iz,:)
          end do
        end do
        do ilat = 1, ugg_south%nlat
          do ilon = 1, ugg_south%nlon
            SA_frac_adv(:,ilon,0,iz,ispec)    = label_bounds(cache_ij_south(1,ilon,ilat),cache_ij_south(2,ilon,ilat),iz,:)
          end do
        end do
        do ilat = 1, ugg_north%nlat
          do ilon = 1, ugg_north%nlon
            SA_frac_adv(:,ilon,ny+1,iz,ispec) = label_bounds(cache_ij_north(1,ilon,ilat),cache_ij_north(2,ilon,ilat),iz,:)
          end do
        end do
      end do  ! levels

      ! clear
      deallocate(label_bounds, stat=status)
      IF_NOTOK_RETURN(status=1)

    end do ! specs

    ! Done
    status = NF90_Close(ncid)
    IF_NF90_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Bound_Setup

! ---------------------------------------------------------------------------------------

  subroutine SA_Bound_Done(status)

    ! --- in/out ------------------------
    integer, intent(out)          ::  status

    ! --- const -------------------------
    character(len=*), parameter :: rname = mname//'/SA_Bound_Done'

    ! --- begin -------------------------

    deallocate(cache_ij_west, stat=status)
    IF_NOTOK_RETURN(status=1)
    deallocate(cache_ij_north, stat=status)
    IF_NOTOK_RETURN(status=1)
    deallocate(cache_ij_east, stat=status)
    IF_NOTOK_RETURN(status=1)
    deallocate(cache_ij_south, stat=status)
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SA_Bound_Done


end module SA_Labeling
