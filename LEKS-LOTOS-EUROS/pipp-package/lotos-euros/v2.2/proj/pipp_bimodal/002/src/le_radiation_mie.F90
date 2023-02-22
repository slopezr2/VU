!###############################################################################
!
! Mie code: optical properties of spherical particles.
! For non-spherical particals a "T-matrix" code should be used.
!
! History
!   2012, Astrid Manders, TNO
!     Original implementation in "le_radiation.F90"
!   2015-06, Arjo Segers, TNO
!     Splitted this module off from "le_radiation.F90".
!     Introduced interpolation in lookup table.
!     Fixed computation of volume mean diameter.
!   2022-03, Janot Tokaya, TNO
!     Perform optics calculations for dust and SS in bimodal way.
!     Used to be 5 modes, but not done with correct sigma's.
!     The 5 modes are usefull for deposition and advection, not optics.
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

module LE_Radiation_Mie

  use GO, only : gol, goPr, goErr

  use Indices, only : N_AEROSOL_MODES

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  T_AerosolMode, T_AerosolModes
  public  ::  T_RefracIndex
  public  ::  calc_properties_mie
#ifdef with_m7
  public  ::  calc_properties_mie_m7
#endif


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Radiation_Mie'

  real, parameter :: le_RH_cutoff = 0.95


  ! --- types ----------------------------------

  ! properties of aerosol mode:
  type T_AerosolMode
    ! numer density definition:
    real        ::  Dg        ! number geometric mean [m]
    real        ::  sigmag    ! geometric std.dev.
    ! computed diameters:
    real        ::  Da2     ! 2nd moment arithmetic mean diameter (surface, cross-section)
    real        ::  Da3     ! 3rd moment arithmetic mean diameter (volume)
  contains
    procedure   ::  Init      => AerosolMode_Init
    procedure   ::  Done      => AerosolMode_Done
  end type T_AerosolMode

  ! *

  ! collection of multiple aerosol modes:
  type T_AerosolModes
    ! number of modes:
    integer                             ::  n
    ! definitions:
    type(T_AerosolMode), allocatable    ::  aerosol(:)  ! (n)
  contains
    procedure   ::  Init      => AerosolModes_Init
    procedure   ::  Done      => AerosolModes_Done
  end type T_AerosolModes
  ! *



  type T_RefracIndex
    complex, allocatable  ::  PerSpecie(:,:)   ! (swbands%n,n_aerosol)
    complex, allocatable  ::  Water(:)         ! (swbands%n)
  contains
    procedure   ::  Init      => RefracIndex_Init
    procedure   ::  Done      => RefracIndex_Done
  end type T_RefracIndex


contains


  ! ====================================================================
  ! ===
  ! === aerosol modes
  ! ===
  ! ====================================================================


  subroutine AerosolMode_Init( self, Dg, sigmag, status )

    use JAQL, only : ParticleMode_AMeanDiameter_Moment

    ! --- in/out ---------------------------------

    class(T_AerosolMode), intent(out)     ::  self
    real, intent(in)                      ::  Dg, sigmag
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AerosolMode_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! store:
    self%Dg     = Dg
    self%sigmag = sigmag

    ! BUG fixed:
    !   Wrong computation of surface/volume mean diameters:
    !     Dv = exp( log(Dpg) + 3.0*(log(sigma)**2.0) ) * exp((log(sigma)**2.0)/2.0)
    !     Ds = exp( log(Dpg) + 2.0*(log(sigma)**2.0) ) * exp((log(sigma)**2.0)/2.0)
    !   This seems a combination of formula (8.53) and (8.50) for Dg,
    !   and (8.44) for arithmetic mean in Seinfeld&Pands (2nd Edition).
    !   However, arithmetic mean diameter of a moment should be computed differently,
    !   see new implementation in JAQL_Particles module.

    ! arithmetic mean diameters representative for
    ! mean 2nd moment (surface,cross-section) and 3rd moment (volume):
    self%Da2 = ParticleMode_AMeanDiameter_Moment( self%Dg, self%sigmag, 2 )
    self%Da3 = ParticleMode_AMeanDiameter_Moment( self%Dg, self%sigmag, 3 )

    ! ok
    status = 0

  end subroutine AerosolMode_Init


  ! ***


  subroutine AerosolMode_Done( self, status )

    ! --- in/out ---------------------------------

    class(T_AerosolMode), intent(inout)   ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AerosolMode_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine AerosolMode_Done


  ! ***


  subroutine AerosolModes_Init( self, status )

    use Indices         , only : specname
    use Indices         , only : n_aerosol, ispecs_aerosol
    use Indices         , only : aerosol_radius_g, aerosol_sigma_g ! m, 1
    use Indices         , only : specmode ! m
    use LE_Config       , only : rcF
    use GO              , only : ReadRc    
    ! --- in/out ---------------------------------

    class(T_AerosolModes), intent(out)    ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AerosolModes_Init'

    ! --- local ----------------------------------

    integer     ::  ispec
    integer     ::  i_aerosol
    integer     ::  i
    real        ::  R_g, sigma_g
    real        ::  R_Na_f, R_Na_c, R_Du_f,  R_Du_c

    ! --- begin ----------------------------------

    ! size:
    self%n = n_aerosol+4
    ! storage:
    allocate( self%aerosol(self%n), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! only loop over aerosol tracers:
    do i_aerosol = 1, n_aerosol
      ! global index:
      ispec = ispecs_aerosol(i_aerosol)
      
      ! aerosol geometric radius and std.dev.:
      R_g     = aerosol_radius_g(ispec)  ! m
      sigma_g = aerosol_sigma_g (ispec)  ! 1

      !write (gol,'("aerosol radius for aerosol (`",a,"`) is f8.3")') &
      !                     trim(specname(ispec)),R_g; call goErr
      write(*,*) trim(specname(ispec))
      write(*,*) R_g

      ! check ...
      if ( R_g < 0.0 ) then
             write (gol,'("no aerosol radius defined for aerosol ",i0," `",a,"` (",i0,")")') &
                           i_aerosol, ispec, trim(specname(ispec)); call goErr
             TRACEBACK; status=1; return
      end if

      ! initialize aerosol description:  diameter [m],  sigma [1]
      call self%aerosol(i_aerosol)%Init(   2 * R_g   ,  sigma_g  , status )
      IF_NOTOK_RETURN(status=1)

    end do ! aerosols

    call ReadRc( rcF, 'my.r_Na_f', r_Na_f, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'my.r_Na_c', r_Na_c, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'my.r_Du_f', r_Du_f, status )
    IF_NOTOK_RETURN(status=1)  
    
    call ReadRc( rcF, 'my.r_Du_c', r_Du_c, status )
    IF_NOTOK_RETURN(status=1)
            
    ! dust and seasalt combined modes -> UGLY and hardcoded for now
    ! ss fine
    i_aerosol = n_aerosol+1
    ! initialize aerosol description:  diameter [m],  sigma [1]
    !call self%aerosol(i_aerosol)%Init(  2*2.5e-6    ,  2.0  , status )  ! Quinn 2015
    call self%aerosol(i_aerosol)%Init(   2*r_Na_f*1.e-6  ,  1.59  , status ) ! Sakai 2010 (we divide by 2 because RH~80 means radius times 2) and Quinn 2015
    IF_NOTOK_RETURN(status=1)
    
    ! ss coarse
    i_aerosol = n_aerosol+2
    ! initialize aerosol description:  diameter [m],  sigma [1]
    call self%aerosol(i_aerosol)%Init(  2*r_Na_c*1.e-6    ,  2.0  , status )  ! For bimodal test
    !call self%aerosol(i_aerosol)%Init(  2.3e-6    ,  2.0  , status )  ! Quinn 2015
    IF_NOTOK_RETURN(status=1)

    ! dust fine
    i_aerosol = n_aerosol+3 ! dust fine mode
    ! initialize aerosol description:  diameter [m],  sigma [1]
    !call self%aerosol(i_aerosol)%Init(   0.12e-6   ,  1.59  , status ) ! mode 2 from denjean 2016
    call self%aerosol(i_aerosol)%Init(   2*r_Du_f*1.e-6   ,  1.59  , status ) ! For bimodal test
    IF_NOTOK_RETURN(status=1)

    ! dust coarse
    i_aerosol = n_aerosol+4 ! dust coarse mode
    ! initialize aerosol description:  diameter [m],  sigma [1]
    call self%aerosol(i_aerosol)%Init(  2*r_Du_c*1.e-6    ,  2.0  , status )  ! For bimodal test
    !call self%aerosol(i_aerosol)%Init(   1.3e-6    ,  2.0  , status ) ! mode 4 from denjean 2016
    IF_NOTOK_RETURN(status=1)
    ! ok
    status = 0

  end subroutine AerosolModes_Init





  ! ***


  subroutine AerosolModes_Done( self, status )

    ! --- in/out ---------------------------------

    class(T_AerosolModes), intent(inout)    ::  self
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/AerosolModes_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%aerosol, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine AerosolModes_Done


  ! ====================================================================
  ! ===
  ! === RefracIndex
  ! ===
  ! ====================================================================


  subroutine RefracIndex_Init( self, lut, swbands, status )

    use Indices
    use LE_Radiation_LUT    , only : T_Radiation_LUT
    use LE_Radiation_SWBands, only : T_SWBands

    ! --- in/out ---------------------------------

    class(T_RefracIndex), intent(out)     ::  self
    class(T_Radiation_LUT), intent(in)    ::  lut
    class(T_SWBands), intent(in)          ::  swbands
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/RefracIndex_Init'

    ! --- local ----------------------------------

    integer                         ::  i, j
    integer                         ::  bandwidth
    real, allocatable               ::  tmp1(:), tmp2(:)
    real                            ::  Re, Im
    integer                         ::  i_aerosol, ispec

    ! --- begin ----------------------------------

    ! storage:
    allocate( self%PerSpecie(swbands%n,n_aerosol+4), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%Water(swbands%n), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over wavelength bands:
    do j = 1, swbands%n

      ! bandwith:
      bandwidth = nint( swbands%width(j)*1000.0 ) ! nm
      ! check ...
      if ( bandwidth < 1 ) then
        write (gol,'("bandwidth ",i0," seem less than 1 nm :",f8.3)') &
                        j, swbands%width(j)*1000.0; call goErr
        TRACEBACK; status=1; return
      end if
      !print *, 'aaa0 j ', j, ' bandwidth ', bandwidth, ' nm'

      ! subgrid of wavelengths per nm:
      allocate( tmp1(bandwidth) )
      do i = 1, bandwidth
        tmp1(i) = swbands%bounds(j,1) + (i-1.0) / 1000.0
      end do
      ! storage for interpolated values:
      allocate( tmp2(bandwidth) )

      ! loop over aersol tracers:
      do i_aerosol = 1, n_aerosol+4
        if (i_aerosol <= n_aerosol) then
        ! global index:
        ispec = ispecs_aerosol(i_aerosol)
        else if (i_aerosol <= n_aerosol+2) then
          ispec = ispec_Na_f1
        else if (i_aerosol <= n_aerosol+4) then
          ispec = ispec_dust_f1
        else
          write (gol,'("JPT: i_aerosol shouldnt be this high")'); call goPr
        end if
        ! switch per aerosol:
        select case ( ispec )

          !
          ! opac table sets BC, SALT, SO4, NO3
          !

#ifdef with_m7
          case ( ispec_SO4a_f, ispec_SO4a_c,  ispec_NH4a_f, ispec_NO3a_f, ispec_NO3a_c, ispec_so4ns, ispec_so4ks, ispec_so4as, ispec_so4cs)

            call rad_interpol( lut%opac(2,:), lut%opac(1,:), tmp1, tmp2 )
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%opac(3,:), lut%opac(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth
            
        
          case ( ispec_EC_f, ispec_EC_c, ispec_bcks, ispec_bcas, ispec_bccs, ispec_bcki )

            call rad_interpol(lut%opac(4,:), lut%opac(1,:), tmp1, tmp2)
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%opac(5,:), lut%opac(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth
           
           
        
          case (ispec_Na_f1, ispec_Na_f2, ispec_Na_c1, ispec_Na_c2, ispec_Na_c3, ispec_ssas, ispec_sscs )

            call rad_interpol(lut%opac(6,:), lut%opac(1,:), tmp1, tmp2)
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%opac(7,:), lut%opac(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth
            
          !
          ! echam_ham table sets PPM25 and DUST refractive indexes
          !
          
          case ( ispec_ppm_f, ispec_ppm_c, ispec_pom_f, ispec_pom_c, ispec_ocks, ispec_ocas, ispec_occs, ispec_ocki)

            call rad_interpol(lut%echam(2,:), lut%echam(1,:), tmp1, tmp2)
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%echam(3,:), lut%echam(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth
            
        
          case ( ispec_dust_f1,ispec_dust_f2, ispec_dust_c1, ispec_dust_c2, ispec_dust_c3, ispec_duas, ispec_ducs, ispec_duai, ispec_duci )

            call rad_interpol(lut%echam(4,:), lut%echam(1,:), tmp1, tmp2)
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%echam(5,:), lut%echam(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth

          !Jianbing: this is the number concentration, ignore 
          case ( ispec_nnus, ispec_nais, ispec_nacs, ispec_ncos, ispec_naii, ispec_naci, ispec_ncoi)
            ! write (gol,'("Jianbing: here should be igored")'); call goPr


#else
          case ( ispec_SO4a_f, ispec_SO4a_c,  ispec_NH4a_f, ispec_NO3a_f, ispec_NO3a_c )

            call rad_interpol( lut%opac(2,:), lut%opac(1,:), tmp1, tmp2 )
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%opac(3,:), lut%opac(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth


          case ( ispec_EC_f, ispec_EC_c )

            call rad_interpol(lut%opac(4,:), lut%opac(1,:), tmp1, tmp2)
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%opac(5,:), lut%opac(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth



          case (ispec_Na_f1, ispec_Na_f2, ispec_Na_c1, ispec_Na_c2, ispec_Na_c3 )

            call rad_interpol(lut%opac(6,:), lut%opac(1,:), tmp1, tmp2)
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%opac(7,:), lut%opac(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth

          !
          ! echam_ham table sets PPM25 and DUST refractive indexes
          !

          case ( ispec_ppm_f, ispec_ppm_c, ispec_pom_f, ispec_pom_c )

            call rad_interpol(lut%echam(2,:), lut%echam(1,:), tmp1, tmp2)
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%echam(3,:), lut%echam(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth


          case ( ispec_dust_f1,ispec_dust_f2, ispec_dust_c1, ispec_dust_c2, ispec_dust_c3 )

            call rad_interpol(lut%echam(4,:), lut%echam(1,:), tmp1, tmp2)
            Re = sum(tmp2)/bandwidth
            call rad_interpol(lut%echam(5,:), lut%echam(1,:), tmp1, tmp2)
            Im = sum(tmp2)/bandwidth

#endif


          case default
            write (gol,'(a,": could not fill RefracIndex%PerSpecie for aerosol ",i6," (",a,")")') rname, ispec, trim(specname(ispec)); call goErr
            TRACEBACK; status=1; return
        end select

        ! store:
        self%PerSpecie(j,i_aerosol) = cmplx(Re,Im)

      end do  ! aerosols

      ! Segelstein is about the water
      call rad_interpol(lut%segelstein(2,:), lut%segelstein(1,:), tmp1, tmp2)
      Re = sum(tmp2)/bandwidth
      call rad_interpol(lut%segelstein(3,:), lut%segelstein(1,:), tmp1, tmp2)
      Im = sum(tmp2)/bandwidth
      self%Water(j)=cmplx(Re,Im)

      ! clear:
      deallocate(tmp1,tmp2)

    end do   ! bands

    ! ok
    status = 0

  end subroutine RefracIndex_Init


  ! ***


  subroutine RefracIndex_Done( self, status )

    ! --- in/out ---------------------------------

    class(T_RefracIndex), intent(inout)   ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/RefracIndex_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%PerSpecie, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%Water, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine RefracIndex_Done


  ! ====================================================================
  ! ===
  ! === Mie
  ! ===
  ! ====================================================================


  subroutine calc_properties_mie( c, rh, swbands, lookup, lut, aer, RefracIndex, &
                                  ext, scat, asym, status, &
                                  sizeparam, refr_spec, extf_spec )

    use Binas               , only : pi
    use Binas               , only : massfrac_Na_in_seasalt
    use JAQL                , only : Particle_Volume_D, Particle_CrossSection_D
    use Num                 , only : T_LookUp
    use Indices
    use LE_Radiation_SWBands, only : T_SWBands
    use LE_Radiation_LUT    , only : T_Radiation_LUT

    use Indices             , only : nspec, specname, specmode
    use Indices             , only : n_aerosol, ispecs_aerosol
    use Indices             , only : ispec_SO4a_f,  ispec_SO4a_c, ispec_NH4a_f, ispec_NO3a_f,   ispec_NO3a_c, ispec_EC_f, ispec_EC_c
    use Indices             , only : ispec_POM_f,   ispec_POM_c,  ispec_PPM_f,  ispec_PPM_c
    use Indices             , only : ispec_Na_f1,   ispec_Na_f2,   ispec_Na_c1,   ispec_Na_c2,    ispec_Na_c3
    use Indices             , only : ispec_dust_f1, ispec_dust_f2, ispec_dust_c1, ispec_dust_c2 , ispec_dust_c3


    ! --- in/out ---------------------------------

    real, intent(in)                    ::  c(:)   ! (nspec) gas (ppb) and aerosol concentrations [ug/m3]
    real, intent(in)                    ::  rh     ! relative humidity [%]
    class(T_SWBands), intent(in)        ::  swbands
    class(T_LookUp), intent(inout)      ::  lookup
    class(T_Radiation_LUT), intent(in)  ::  lut
    class(T_AerosolModes), intent(in)   ::  aer
    class(T_RefracIndex), intent(in)    ::  RefracIndex
    real, intent(out)                   ::  ext(:)     ! (bands%n) extinction coefficent
    real, intent(out)                   ::  scat(:)    ! (bands%n) single scattering albedo
    real, intent(out)                   ::  asym(:)    ! (bands%n) asymmetry factor
    integer, intent(out)                ::  status

    ! testing ...
    real, intent(out), optional         ::  sizeparam(:,:)  ! (bands%n,nspec)
    complex, intent(out), optional      ::  refr_spec(:,:)  ! (bands%n,nspec)
    real, intent(out), optional         ::  extf_spec(:,:)  ! (bands%n,nspec)

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/calc_properties_mie'

    real, parameter  ::  kg_to_g = 1.0e3
    real, parameter  ::  mug_to_g= 1.0e-6

    logical, parameter :: print_warning = .true.

    ! --- local ----------------------------------

    ! logical                ::  mode_is_used(N_AEROSOL_MODES)

    real                   ::  TolVolumePerSpecie (n_aerosol+4) ! the 4 additional entries have the dust and ss fine and coarse content
    real                   ::  DryVolPerSpecie    (n_aerosol+4)
    real                   ::  WaterPerSpecie     (n_aerosol+4)
    real                   ::  Growth             (n_aerosol+4)
    real                   ::  Kappa              (n_aerosol+4)
    real                   ::  NumberConcentration(n_aerosol+4)
    real                   ::  MedianWetRadius    (n_aerosol+4)
    real                   ::  MeanWetCrossSection
    complex, allocatable   ::  RefracInd(:,:)  ! (n_aerosol, RRTM_NO_SWBANDS)
    complex                ::  RInd
    real                   ::  x
    real                   ::  k_ext(n_aerosol+4), k_sca(n_aerosol+4), AF(n_aerosol+4)
    real                   ::  lambda
    real                   ::  lu_ext, lu_a, lu_g
    real                   ::  rhl
    real                   ::  scattering
    real                   ::  totvolume, totnumber
    integer                ::  ispec
    integer                ::  jband
    integer                ::  i_aerosol
    real                   ::  csa

    ! --- begin ----------------------------------

    ! storage:
    allocate( RefracInd(n_aerosol+4, swbands%n), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! Calculate volume concentrations of the species

    ! init sum:
    DryVolPerSpecie = 0.0
    TolVolumePerSpecie = 0.0

    do i_aerosol = 1, n_aerosol
      ! global index:
      ispec = ispecs_aerosol(i_aerosol)

      if ( c(ispec) < 0.0 ) cycle

      DryVolPerSpecie(i_aerosol) = mug_to_g * c(ispec) / ( aerdens(ispec) * kg_to_g )

      ! convert from Na to seasalt if necessary:
      if ( (ispec == ispec_Na_f2) .or. (ispec == ispec_Na_f1) &
           .or. (ispec == ispec_Na_c1).or. (ispec == ispec_Na_c2).or. (ispec == ispec_Na_c3)) then
        DryVolPerSpecie(i_aerosol) = DryVolPerSpecie(i_aerosol) / massfrac_Na_in_seasalt  ! Na -> SS
      end if

    end do ! aerosols

    !write(*,*) 'DryVolPerSpecie', DryVolPerSpecie

    ! add the fine and coarse modes together in the last 4 entries for dust and seasalt
    do i_aerosol = 1, n_aerosol
      ! global index:
      ispec = ispecs_aerosol(i_aerosol)
      if ( c(ispec) < 0.0 ) cycle

      ! put f1 and f2 in one mode for seasalt:
      if ( (ispec == ispec_Na_f2) .or. (ispec == ispec_Na_f1)) then
        DryVolPerSpecie(n_aerosol+1) = DryVolPerSpecie(n_aerosol+1) + DryVolPerSpecie(i_aerosol) !/ massfrac_Na_in_seasalt  ! Na -> SS conversion was already done in previous loop
      end if
      ! put c1, c2 and c3 in one mode for seasalt:
      if ( (ispec == ispec_Na_c1).or. (ispec == ispec_Na_c2).or. (ispec == ispec_Na_c3)) then
        DryVolPerSpecie(n_aerosol+2) = DryVolPerSpecie(n_aerosol+2) + DryVolPerSpecie(i_aerosol) !/ massfrac_Na_in_seasalt  ! Na -> SS conversion was already done in previous loop
      end if

      ! put f1 and f2 in one mode for dust:
      if ( (ispec == ispec_dust_f2) .or. (ispec == ispec_dust_f1)) then
        DryVolPerSpecie(n_aerosol+3) = DryVolPerSpecie(n_aerosol+3) + DryVolPerSpecie(i_aerosol)
      end if
      ! put c1, c2 and c3 in one mode for seasalt:
      if ( (ispec == ispec_dust_c1).or. (ispec == ispec_dust_c2).or. (ispec == ispec_dust_c3)) then
        DryVolPerSpecie(n_aerosol+4) = DryVolPerSpecie(n_aerosol+4) + DryVolPerSpecie(i_aerosol)
      end if

    end do ! aerosols

    !for checking the sequence
    !do i_aerosol = 1, n_aerosol
    !  ispec = ispecs_aerosol(i_aerosol)
    !  write(*,*) trim(specname(ispec))
    !end do

    !write(*,*) 'DryVolPerSpecie 4 added', DryVolPerSpecie


    ! Calculate average hygroscopic parameter kappa per species,
    ! and afterwards determine hygroscopic growth factor
    Kappa = 0.0
    ! loop over species:
    do i_aerosol = 1, n_aerosol

      ispec = ispecs_aerosol(i_aerosol)
      ! check values
      if ( DryVolPerSpecie(i_aerosol) == 0.0 ) cycle

      Kappa(i_aerosol) = aerhygro(ispec)

    end do ! aerosols

    ! kappas for dust and seasalt are the same in the combined modes as in the quintmodes
    !do i_aerosol = n_aerosol+1, n_aerosol + 3
    Kappa(n_aerosol+1) = aerhygro(ispec_Na_c1)
    Kappa(n_aerosol+2) = aerhygro(ispec_Na_c1)
    !end do
    !do i_aerosol = n_aerosol + 3, n_aerosol + 5
    !  Kappa(i_aerosol) = aerhygro(ispec_dust_c1)
    !end do
    Kappa(n_aerosol+3) = aerhygro(ispec_dust_c1)
    Kappa(n_aerosol+4) = aerhygro(ispec_dust_c1)
    !write(*,*) 'Kappa 4 added', Kappa

    ! init:
    Growth = 0.0
    WaterPerSpecie = 0.0
    TolVolumePerSpecie = 0.0
    NumberConcentration = 0.0

    !
    ! Calculate water uptake and sum to total aerosol volume
    !

    ! RH from % to fraction:
    rhl = rh / 100.0  ! 0-1
    ! truncate:
    if ( rhl > le_RH_cutoff ) rhl = le_rh_cutoff
    if ( rhl < 0            ) rhl = 0.0

    ! loop over aerosol modes:
    do i_aerosol = 1, n_aerosol+4
      ! skip if no Dry volume or mass
      if ( DryVolPerSpecie(i_aerosol) == 0.0 ) cycle

      ! compute growth:
      Growth(i_aerosol) = (1 + Kappa(i_aerosol) * (rhl / (1 - rhl))) ** 0.333 !check units RH

      ! Calculate water uptake per species:
      WaterPerSpecie(i_aerosol) = (Growth(i_aerosol)**3.0 -1 ) * DryVolPerSpecie(i_aerosol)

      ! Calculate sum to total aerosol volume
      TolVolumePerSpecie(i_aerosol) = DryVolPerSpecie(i_aerosol) + WaterPerSpecie(i_aerosol)

      ! calculate number concentration [count/m3]
      NumberConcentration(i_aerosol) = DryVolPerSpecie(i_aerosol) / Particle_Volume_D( aer%aerosol(i_aerosol)%Da3 )

      ! part 6, calculate median radius of the humified aerosols (per mode)
      MedianWetRadius(i_aerosol)      =      0.5 * aer%aerosol( i_aerosol )%Dg * Growth(i_aerosol)

    end do ! modes

    !write(*,*) 'NumberConcentration 4 added', NumberConcentration


    ! init:
    RefracInd = (0.0,0.0)

    do jband = 1, swbands%n

      ! loop over aersol tracers:
      do i_aerosol = 1, n_aerosol+4

        ! skip if no Dry volume or mass
        if ( DryVolPerSpecie(i_aerosol) == 0.0 ) cycle

        ! add contribution:
        RefracInd(i_aerosol,jband) = RefracInd(i_aerosol,jband) + &
                                      DryVolPerSpecie(i_aerosol) * RefracIndex%PerSpecie(jband,i_aerosol)

        RefracInd(i_aerosol,jband) = RefracInd(i_aerosol,jband) + WaterPerSpecie(i_aerosol) * RefracIndex%Water(jband)

        if ( TolVolumePerSpecie(i_aerosol) > 0.0 ) then
          RefracInd(i_aerosol,jband) = RefracInd(i_aerosol,jband) / TolVolumePerSpecie(i_aerosol)
        else
          RefracInd(i_aerosol,jband) = (0.0,0.0)
        end if

      end do ! aerosol

    end do ! bands

    ! part 7. extract information from look-up tables, interpolate
    ext    = 0.0
    asym   = 0.0
    scat   = 0.0

    ! testing ...
    if ( present(sizeparam) ) sizeparam = 0.0
    if ( present(refr_spec) ) refr_spec = 0.0
    if ( present(extf_spec) ) extf_spec = 0.0

    ! loop over bands:
    do jband = 1, swbands%n

      scattering=0.
      totvolume =0.
      totnumber =0.

      ! loop over aerosols:
      do i_aerosol = 1, n_aerosol+4
        ! skip if no Dry volume or mass
        if ( DryVolPerSpecie(i_aerosol) == 0.0 ) cycle

        ! we now skip over dust and ss modes in the 'quintmodal' distributions
        if (i_aerosol<=n_aerosol) then
          ispec = ispecs_aerosol(i_aerosol)
          if ( (ispec == ispec_Na_f2)  .or. (ispec == ispec_Na_f1) &
          .or. (ispec == ispec_Na_c1)  .or. (ispec == ispec_Na_c2) .or. (ispec == ispec_Na_c3)  &
          .or. (ispec == ispec_dust_f2).or. (ispec == ispec_dust_f1) &
          .or. (ispec == ispec_dust_c1).or. (ispec == ispec_dust_c2).or.(ispec == ispec_dust_c3)) cycle
        end if

        ! extract:
        RInd   = RefracInd(i_aerosol, jband)

        ! wavelength in m
        lambda = 1.0D-6 * swbands%lambda(jband)

        ! size parameter [JadB-2013, eq. (6.9) or (6.19)],
        ! based on the geometric mean (median) radius,
        ! used as index in the lookup tables:
        x = 2.0 * pi * MedianWetRadius(i_aerosol) / lambda

        ! set interpolation indices:
        call lookup%InterpolSet( (/x,real(RInd),aimag(RInd)/), status )
        IF_NOTOK_RETURN(status=1)

        ! switch:
        ! modes with sigma=1.59 (or smaller):
        if ( (abs(aer%aerosol(i_aerosol)%sigmag - 1.59) < 0.01) .or. &
             (abs(aer%aerosol(i_aerosol)%sigmag - 1.35) < 0.01)      ) then

          ! apply interpolations in lookup tables:
          call lookup%InterpolApply( lut%ext_159, lu_ext, status )
          IF_NOTOK_RETURN(status=1)
          call lookup%InterpolApply( lut%a_159  , lu_a  , status )
          IF_NOTOK_RETURN(status=1)
          call lookup%InterpolApply( lut%g_159  , lu_g  , status )
          IF_NOTOK_RETURN(status=1)

        ! modes with sigma=2.00
        else if ( (abs(aer%aerosol(i_aerosol)%sigmag - 2.00) < 0.01) ) then

          ! apply interpolations in lookup tables:
          call lookup%InterpolApply( lut%ext_200, lu_ext, status )
          IF_NOTOK_RETURN(status=1)
          call lookup%InterpolApply( lut%a_200  , lu_a  , status )
          IF_NOTOK_RETURN(status=1)
          call lookup%InterpolApply( lut%g_200  , lu_g  , status )
          IF_NOTOK_RETURN(status=1)

        else
          write (gol,'("unsupported sigma ",f0.2," for aerosol ",i6)') aer%aerosol(i_aerosol)%sigmag, i_aerosol; call goErr
          TRACEBACK; status=1; return
        end if

        ! TOO MUCH OUTPUT
        !write (gol,'("Number conc is ",f0.2," for aerosol ",i6)') NumberConcentration(i_aerosol), i_aerosol; call goErr
        !write (gol,'("Dry volume is ",f0.2," for aerosol ",i6)') DryVolPerSpecie(i_aerosol), i_aerosol; call goErr

        ! fill:
        k_ext(i_aerosol) = lu_ext * lambda * lambda * NumberConcentration(i_aerosol)
        k_sca(i_aerosol) = lu_a * k_ext(i_aerosol)
        AF(i_aerosol)    = lu_g

        ! update sums:
        ext(jband)  = ext(jband)  + k_ext(i_aerosol)
        scattering  = scattering  + k_sca(i_aerosol)
        asym(jband) = asym(jband) + k_sca(i_aerosol)*AF(i_aerosol)

        
        ! ! testing ...
        ! if ( present(sizeparam)  ) then
          ! ! loop over aerosols:
          ! do i_aerosol = 1, n_aerosol
          !   ! current mode?
          !   if ( (DryVolPerSpecie(i_aerosol) > 0.0) .and. ( NumberConcentration(i_aerosol) > 0.0) ) then
        if ( present(sizeparam) ) sizeparam(jband,i_aerosol) = x
          !   end if
          ! end do ! aerosols
        ! end if  ! store


      end do ! loop i_aerosol


      ! Finally determine extinction coefficient, AOD (tau),
      ! single scattering coefficient ssa and asymmetry factor asy
      ! should be altered when there are more than 2 modes
      if ( ext(jband) > 1.0d-7 ) then
        scat(jband)  = scattering/ ext(jband)
        if (scat(jband) > 0.000000001) then
          asym(jband)  = asym(jband) / scattering
        else
          asym(jband)  = 1.0
        end if
      else
        scat(jband)  = 0.0
        asym(jband)  = 1.0
      end if

    end do ! loop j over swbands

    ! clear:
    deallocate( RefracInd, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine calc_properties_mie


  ! ***

#ifdef with_m7
  subroutine calc_properties_mie_m7( c, rad_wet_m7, rad_dry_m7, water_m7, rh, swbands, lookup, lut, RefracIndex, &
                                  ext, scat, asym, &
                                  status)
                                  ! , &
                                  ! refr_f, refr_c, &
                                  ! volf, volc, crosfine, croscoarse, &
                                  ! Nfine, Ncoarse, status, &
                                  ! sizeparam, refr_spec, extf_spec )

    use Binas               , only : pi
    ! use Binas               , only : massfrac_Na_in_seasalt
    ! use JAQL                , only : Particle_Volume_D, Particle_CrossSection_D
    use Num                 , only : T_LookUp
    use Indices
    use LE_Radiation_SWBands, only : T_SWBands
    use LE_Radiation_LUT    , only : T_Radiation_LUT
    use mo_aero_m7,           only : n_m7numb => nmod     ! n_m7numb = 7    
    use indices,              only : n_m7  ! if n_m7 > 0, then m7 mode is enabled
    use indices,              only : m7_mode
    ! use indices,              only : m7_numberconc_ispecs
    use indices,              only : tracer_is_m7, tracer_is_numberconc
    use mo_aero_m7,           only : inucs, iaits, iaccs, icoas, iaiti, iacci, icoai

    use Indices,              only : specname ! this is only for the Refraction index check 

    ! --- in/out ---------------------------------

    real, intent(in)                    ::  c(:)    ! (nspec) gas (ppb) and aerosol concentrations [ug/m3], and number concentration !!!!
    real, intent(in)                    ::  rh      ! ralative humidity, but actually I don't use it 
    real, intent(in)                    ::  rad_wet_m7(:)  ! wet radius each mode in m7  
    real, intent(in)                    ::  rad_dry_m7(:)  ! dry radius each mode in m7, it will not be used, only to check   
    real, intent(in)                    ::  water_m7(:)  ! wet radius each mode  
    class(T_SWBands), intent(in)        ::  swbands ! Jianbing: in the test case, we have 14 swands
    class(T_LookUp), intent(inout)      ::  lookup  ! the 
    class(T_Radiation_LUT), intent(in)  ::  lut
    ! class(T_AerosolModes), intent(in)   ::  aer    ! jianbing: the aerosol model (two mode?)
    class(T_RefracIndex), intent(in)    ::  RefracIndex  ! jianbing: set up from the 
    real, intent(out)                   ::  ext(:)     ! (bands%n) extinction coefficent
    real, intent(out)                   ::  scat(:)    ! (bands%n) single scattering albedo
    real, intent(out)                   ::  asym(:)    ! (bands%n) asymmetry factor
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/calc_properties_mie_m7'
    
    real, parameter  ::  kg_to_g  = 1.0e3
    real, parameter  ::  mug_to_g  = 1.0e-6
    real, parameter  ::  cm3_to_m3 = 1.0e6       ! per cm3 to per m3 
    real, parameter  ::  cm_to_m  = 1.0e-2       ! jianbing: this is wrong!!!!!
    logical, parameter :: print_warning = .true.

    ! --- local ----------------------------------
    
    logical                ::  mode_is_used(n_m7numb)
    logical                ::  num_mode_is_used(n_m7numb)

    ! jianbing: for mass conversion 
    real                   ::  mass_c(nspec) 
    real                   ::  mass_mode(n_m7numb)
    real                   ::  VolumePerSpecie(nspec)
    real                   ::  WaterVolPerMode(n_m7numb)
    real                   ::  DryVolPerMode  (n_m7numb)
    real                   ::  rad_wet_m7_correct(n_m7numb)
    ! real                   ::  Growth         (n_m7numb)
    ! real                   ::  Kappa          (n_m7numb)
    real                   ::  TotVolPerMode  (n_m7numb)
    real                   ::  NumberConcentration(n_m7numb)
    ! real                   ::  MedianWetRadius(N_AEROSOL_MODES)
    !real                   ::  MedianWetCrossection(N_AEROSOL_MODES)
    ! real                   ::  MeanWetCrossSection

    complex, allocatable   ::  RefracInd(:,:)  ! (N_AEROSOL_MODES,RRTM_NO_SWBANDS)
    complex                ::  RInd
    real                   ::  x
    real                   ::  k_ext(n_m7numb), k_sca(n_m7numb), AF(n_m7numb)
    real                   ::  lambda
    real                   ::  lu_ext, lu_a, lu_g
    ! real                   ::  rhl
    real                   ::  scattering
    real                   ::  totvolume, totnumber
    integer                ::  ispec
    integer                ::  i_m7mode
    integer                ::  jband
    integer                ::  i_aerosol
    real                   ::  csa

    ! --- begin ----------------------------------
    
    ! storage:
    allocate( RefracInd(n_m7numb,swbands%n), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! write (gol,'("Jianbing to see here 1 ")'); call goPr

    ! calculate the mass concentration    
    mass_c = 0.
    do ispec=1,nspec

      if ( tracer_is_numberconc(ispec) ) cycle

      if ( c(ispec) < 0.0 ) then 
        mass_c(ispec) = 0.0
      end if 


      if ( (ispec == ispec_so4ns) .or. (ispec == ispec_so4ks) &
           .or. (ispec == ispec_so4as).or. (ispec == ispec_so4cs) ) then
          mass_c( ispec ) =  c( ispec ) * 1.6e-10 ! to check this with Astrid 
      else  
          mass_c(ispec) = c( ispec )
          ! mass_c(ispec) = 0.0
      end if

    end do !species

    ! write (gol,'("Jianbing to see here 2 ")'); call goPr

    ! ! read the number concentration from c(;)
    Numberconcentration = 0.
    ! num_mode_is_used = .false.
    do ispec=1,nspec
      
      if ( c(ispec) == 0.0 ) cycle

      if ( tracer_is_numberconc(ispec) ) then

        i_m7mode = m7_mode(ispec)      
        ! ! set flag:
        ! num_mode_is_used(i_m7mode) = .true. 
        !               count/m3        count/cm3   cm3>m3   
        Numberconcentration(i_m7mode) = c(ispec) * cm3_to_m3

        ! normal value for c(ispec) should be around several thousand to 50 k at most
        if (Numberconcentration(i_m7mode) < 0.0 ) then 
           write (gol,'("Jianbing to see num conc is ", es12.3, " for mode ", i2)') Numberconcentration(i_m7mode), i_m7mode; call goErr
           Numberconcentration(i_m7mode) = 0.0 
          ! TRACEBACK; status=1; return
        end if 
      end if ! if tracer is the num conc 
    end do ! species 

    ! write (gol,'("Jianbing to see here 3 ")'); call goPr

    ! Calculate volume concentrations of the species, and calculate total volume per mode
    ! init flags:
    mode_is_used = .false.
    ! init sum:
    DryVolPerMode = 0.
    VolumePerSpecie = 0.
    !scat=0.0
    !asym=0.0
    mass_mode = 0.
    ! loop over species:
    do ispec = 1, nspec
      ! only when we have some concentration 
      if ( mass_c(ispec) == 0.0 ) cycle
      ! only when it is the m7 tracer 
      if ( .not. tracer_is_m7(ispec) ) cycle
      ! only when it is not numberconc tracer 
      if ( tracer_is_numberconc(ispec) ) cycle

      i_m7mode = m7_mode(ispec)      

      ! set flag:
      mode_is_used(i_m7mode) = .true.
      
      ! ! jianbing's check: whether we have number concentration eabled tother with the mass concentation 
      ! if ( .not. num_mode_is_used(i_m7mode) ) then 
      !     write (gol,'(" we have mass concentration but no num concentration at mode ", i2, " mode ")') i_m7mode; call goErr
      !         TRACEBACK; status=1; return
      ! end if 
         !       m3/m3              g/ug       ug/m3         kg/m3            g/kg   
         !Jianbing   aerodens denotes the aerosol density, for example the dust is 2.65 10^3 kg/m3 
      VolumePerSpecie(ispec) = mug_to_g * mass_c(ispec) / ( aerdens(ispec) * kg_to_g )
 
      mass_mode(i_m7mode)    = mass_mode(i_m7mode) + mass_c(ispec)    

      ! add contribution per mode:
      DryVolPerMode(i_m7mode) = DryVolPerMode(i_m7mode) + VolumePerSpecie(ispec) 

      ! !!!!!!!!!!!!!!!!! 
      ! if ( Numberconcentration(i_m7mode) == 0.0 .and. mass_mode(i_m7mode) > 1.0E-2 ) then 
      !     write (gol,'("Jianbing to see the mismatch at mode", i2)')  i_m7mode; call goPr
      !     write (gol,'("Jianbing to see the Mass conc ",   es12.3)')  mass_mode(i_m7mode); call goPr
      !     write (gol,'("Jianbing to see the Num  conc",   es12.3)')   Numberconcentration(i_m7mode); call goPr
      ! end if 

    end do !species

    ! write (gol,'("Jianbing to see here 4 ")'); call goPr


    ! ! this is only the mass and num conc consistence test 
    ! do i_m7mode = 1, n_m7numb
    !         !!!!!!!!!!!!!!!!! 
    !   ! if ( Numberconcentration(i_m7mode) == 0.0 .and. mass_mode(i_m7mode) > 1.0E-2 ) then 
    !   if ( mass_mode(i_m7mode) > 1.0) then 
    !       write (gol,'("Jianbing to see the mismatch at mode                                              : ", i2)')  i_m7mode; call goPr
    !       write (gol,'("Jianbing to see the Mass conc ",   es12.3)')  mass_mode(i_m7mode); call goPr
    !       write (gol,'("Jianbing to see the Num  conc",   es12.3)')   Numberconcentration(i_m7mode); call goPr
    !       if (i_m7mode > 4) then 
    !         write (gol,'("Jianbing to see the dry radius",   es12.3, "m")')   rad_wet_m7(i_m7mode) * 0.01; call goPr
    !       else 
    !         write (gol,'("Jianbing to see the dry radius",   es12.3, "m")')   rad_dry_m7(i_m7mode) * 0.01; call goPr
    !       end if 
    !   end if 
    ! end do 





    ! calculate the water volume concentration of permode
    ! intial water volume
    WaterVolPerMode = 0.
    TotVolPerMode = 0.
    rad_wet_m7_correct = 0.
    
    do i_m7mode = 1, n_m7numb


      if ( DryVolPerMode(i_m7mode) == 0.0 ) cycle
      if ( Numberconcentration(i_m7mode) == 0.0 ) cycle

      rad_wet_m7_correct(i_m7mode) = rad_wet_m7(i_m7mode)  ! copy to correct version 

      select case (i_m7mode)       
      ! this is the soluable mode  
      case( inucs, iaits, iaccs, icoas ) 
          !       m3 / m3                  g/cm3               m3/cm3      / (  kg/m3 * g/kg  )        
          ! to check the unit convert here                               aerdens%water
          ! WaterVolPerMode(i_m7mode) =  water_m7(i_m7mode)  * cm3_to_m3   / ( 1.0e3 ) 
          ! I don't know why, but it is like this, to discuss with Bas
          WaterVolPerMode(i_m7mode) =  water_m7(i_m7mode) * cm3_to_m3 * 1.0e-3  / ( 1.0e3 * kg_to_g ) 

          ! write (gol,'("Jianbing to see the Water ", es12.3)')  WaterVolPerMode(i_m7mode); call goPr
          ! write (gol,'("Jianbing to see the Dry ",   es12.3)')  DryVolPerMode(i_m7mode); call goPr
          ! write (gol,'("Jianbing to see the Num ",   f10.0)')  Numberconcentration(i_m7mode); call goPr

          ! if ( (WaterVolPerMode(i_m7mode) > 0.0 ) .and. ( DryVolPerMode(i_m7mode) < 0.0 ) ) then 
          !     write (gol,'("unsupported waterm7mode in ", i2, " mode ")') i_m7mode; call goErr
          !         TRACEBACK; status=1; return
          ! end if 

          TotVolPerMode(i_m7mode) = WaterVolPerMode(i_m7mode) + DryVolPerMode(i_m7mode) 

          ! if ( rad_dry_m7(i_m7mode) > 0 ) then 
          if ( rad_wet_m7(i_m7mode) / rad_dry_m7(i_m7mode) > 3.8 ) then  
              ! write (gol,'("Jianbing to see the dry radius is ", es12.4, "um for mode", i2)')    rad_dry_m7(i_m7mode)*1.0e5,   i_m7mode; call goPr
              ! write (gol,'("Jianbing to see the wet radius is ", es12.4, "um for mode", i2)')    rad_wet_m7(i_m7mode)*1.0e5,   i_m7mode; call goPr
              rad_wet_m7_correct(i_m7mode) = 3.8 * rad_dry_m7(i_m7mode)  ! sometimes we have too severe data here        
              ! rad_wet_m7(i_m7mode) = rad_dry_m7(i_m7mode) * 2.8
              ! write (gol,'("Jianbing: input break the Growth RULE!")'); call goErr
              ! TRACEBACK; status=1; return
          end if
          ! end if  

          ! write (gol,'("****                                          ***********")'); call goPr  
          ! write (gol,'("Jianbing to see the dry radius is ", es12.4, "um for mode", i2)')    rad_dry_m7(i_m7mode),   i_m7mode; call goPr
          ! write (gol,'("Jianbing to see the wet radius is ", es12.4, "um for mode", i2)')    rad_wet_m7(i_m7mode),   i_m7mode; call goPr
          ! write (gol,'("Jianbing to see the m7 water input is", es12.4, "for mode", i2)')    water_m7(i_m7mode),     i_m7mode; call goPr 
          ! write (gol,'("Jianbing to see the dry volume ", es12.4, "m3 for mode", i2)')       DryVolPerMode(i_m7mode),    i_m7mode; call goPr
          ! write (gol,'("Jianbing to see the water volume is ", es12.4, "m3 for mode", i2)')  WaterVolPerMode(i_m7mode),  i_m7mode; call goPr
          ! write (gol,'("Jianbing to see the ratio between water and dry volume is ", f20.10, " for mode", i2)')  WaterVolPerMode(i_m7mode)/DryVolPerMode(i_m7mode), i_m7mode; call goPr

      ! this is the insoluable mode
      case( iaiti, iacci, icoai )
          WaterVolPerMode(i_m7mode) = 0.0 
          TotVolPerMode(i_m7mode) = WaterVolPerMode(i_m7mode) + DryVolPerMode(i_m7mode)
      case default 
          write (gol,'("Jianbing: Kidding?")'); call goErr
          TRACEBACK; status=1; return
      end select 
    
    end do

    ! write (gol,'("Jianbing to see here 5 ")'); call goPr

    ! to calculate the RefracInd(i_m7mode, i_swband)
    ! init: size of imode, jband
    RefracInd = (0.0,0.0) 

    do jband = 1, swbands%n
      ! do ispec = 1, nspec
      do i_aerosol = 1, n_aerosol! this is a bad loop, need to be updated!!!! 
        ispec = ispecs_aerosol(i_aerosol)
        
        ! only when we have some concentration 
        ! if ( mass_c(ispec) == 0.0 ) cycle
        ! only when it is the m7 tracer 
        ! if ( .not. tracer_is_m7(ispec) ) cycle
        ! only when it is not numberconc tracer 
        ! if ( tracer_is_numberconc(ispec) ) cycle
        
        i_m7mode = m7_mode(ispec)

        RefracInd(i_m7mode,jband) = RefracInd(i_m7mode,jband) + &
                                   VolumePerSpecie(ispec) * RefracIndex%PerSpecie(jband,i_aerosol)  ! this must be i_aerosol, very important!!!!! 
      end do ! aerosol 
      
      ! loop over m7 mode:
      do i_m7mode = 1, n_m7numb
        
        if ( DryVolPerMode(i_m7mode) == 0.0 ) cycle
        if ( Numberconcentration(i_m7mode) == 0.0 ) cycle

        ! add water contribution:
        select case (i_m7mode)       
        ! this is the soluable mode  
        case( inucs, iaits, iaccs, icoas ) 
            RefracInd(i_m7mode,jband) = RefracInd(i_m7mode,jband) + WaterVolPerMode(i_m7mode) * RefracIndex%Water(jband)
        case( iaiti, iacci, icoai )
            RefracInd(i_m7mode,jband) = RefracInd(i_m7mode,jband) + WaterVolPerMode(i_m7mode) * RefracIndex%Water(jband)
        case default 
          write (gol,'("Jianbing: Kidding 2?")'); call goErr
          TRACEBACK; status=1; return
        end select

        if ( TotVolPerMode(i_m7mode) > 0.0 ) then
          RefracInd(i_m7mode,jband) = RefracInd(i_m7mode,jband) / TotVolPerMode(i_m7mode)
        else
          RefracInd(i_m7mode,jband) = (0.0, 0.0)
        end if
      end do ! modes
    end do ! bands






    ! write (gol,'("Jianbing to see here 6 ")'); call goPr


    ! write (gol,'("Jianbing to see Refraction index check" )') ; call goPr
    ! ! loop over aersol tracers:
    ! do i_aerosol = 1, n_aerosol
    !     ! global index:
    !     ispec = ispecs_aerosol(i_aerosol)


    ! ! ! the test to output the Refractin index 
    ! ! ! init:
    ! ! do ispec = 1, nspec

    !     ! if ( .not. tracer_is_m7(ispec) ) cycle
    !     ! ! only when it is not numberconc tracer 
    !     ! if ( tracer_is_numberconc(ispec) ) cycle
        
    !     ! RefracIndex%PerSpecie(10,ispec) 
    !     write (gol,'("Jianbing to see tracer", a, " has the Refraction Real: ", es12.3, ", Ima:" es12.3)') trim(specname(ispec)), real( RefracIndex%PerSpecie(10,i_aerosol) ), aimag( RefracIndex%PerSpecie(10, i_aerosol) ); call goPr
    !     ! write (gol,'("Jianbing to see tracer", a, " has the Refraction Real: ", es12.3, ", Ima:" es12.3)') trim(specname(ispec)),  RefracIndex%PerSpecie(10,ispec)(0), RefracIndex%PerSpecie(10,ispec)(1) ; call goPr

    ! end do ! species or tracers
      
    ! write (gol,'("Jianbing to see water has the Refraction Real: ", es12.3, ", ma:" es12.3)') real( RefracIndex%Water(10) ), aimag( RefracIndex%Water(10) ); call goPr
    ! ! write (gol,'("Jianbing to see water has the Refraction Real: ", es12.3, ", ma:" es12.3)') RefracIndex%Water(10)(0), RefracIndex%Water(10)(1); call goPr

    ! write (gol,'("Jianbing to see Refraction index check end" )') ; call goPr
    
    








    ! part 7. extract three infos from the lookuptable 
    ext    = 0.0  
    asym   = 0.0
    scat   = 0.0


    ! loop over bands:
    do jband = 1, swbands%n 

      scattering=0.0
      ! totvolume =0.0   
      ! totnumber =0.0
      k_ext = 0.
      k_sca = 0.
      AF    = 0.

      ! loop over modes: 
      do i_m7mode = 1, n_m7numb
        ! skip if not used:
        if ( .not. mode_is_used(i_m7mode) ) cycle
        if ( DryVolPerMode(i_m7mode) == 0.0 ) cycle
        if ( Numberconcentration(i_m7mode) == 0.0 ) cycle

        ! extract:
        RInd   = RefracInd(i_m7mode, jband)

        ! ! ! wavelength in m
        ! lambda = 1.0D-6 * swbands%lambda(jband)

        ! amazing unit: in standard LE, it is um, in M7, it is
        lambda = 1.0D-6 * swbands%lambda(jband)

        x = 0.0 

        ! this is different from the normal mie_calc
        ! since we have the 
        ! used as index in the lookup tables:
        ! to check the unit of the rad_wet_m7 and rad_dry_m7
        select case (i_m7mode)       
        ! soluable bins
        case( inucs, iaits, iaccs, icoas ) 
            if (rad_wet_m7_correct(i_m7mode) < 0.0) then 
              write (gol,'("Jianbing: strange wet radius")'); call goErr
              TRACEBACK; status=1; return
            end if
            ! jianbing here we have wrong rad_wet_m7 read, 10 times, need to discuss with astrid 
            x = 2.0 * pi * rad_wet_m7_correct(i_m7mode) * cm_to_m  / lambda
            ! if(  (x > 10.0) .and. (jband == 10) .and. (i_m7mode == 3) ) then 
            !   write (gol,'("JB 2 see the x is ", es12.4, "m for mode", i2)') x, i_m7mode; call goPr
            !   write (gol,'("JB 2 see the dry radius is ", es12.4, "m for mode", i2)')   rad_dry_m7(i_m7mode)*cm_to_m, i_m7mode; call goPr
            !   write (gol,'("JB 2 see the corr wet radius is ", es12.4, "m for mode", i2)')   rad_wet_m7_correct(i_m7mode)*cm_to_m, i_m7mode; call goPr
            !   write (gol,'("JB 2 see the lamda  ", es12.4, "m" )') lambda; call goPr
            !   write (gol,'("JB 2 see the num conc is ", es12.4, " per m3 for mode", i2)') Numberconcentration(i_m7mode),   i_m7mode; call goPr
            !   write (gol,'("JB 2 see the mass conc is ", es12.4, " ug per m3 sea salt accumulation")')  mass_c(i_ssas); call goPr
            ! end if 
        ! insoluable bins 
        case( iaiti, iacci, icoai )
            ! 3 is the number that need to be explain to Bas, here we are talking about the the 4 dry mode 
            if (rad_dry_m7(i_m7mode-3) < 0.0) then 
              write (gol,'("Jianbing: strange dry radius")'); call goErr
              TRACEBACK; status=1; return
            end if
            ! jianbing here we have wrong rad_wet_m7 read, 10 times, need to discuss with astrid 
            x = 2.0 * pi * rad_dry_m7(i_m7mode-3) * cm_to_m / lambda
            ! if (rad_wet_m7_correct(i_m7mode) < 0.0) then 
            !   write (gol,'("Jianbing: strange dry radius")'); call goErr
            !   TRACEBACK; status=1; return
            ! end if
            ! x = 2.0 * pi * rad_wet_m7_correct(i_m7mode) * cm_to_m / lambda
        case default 
          write (gol,'("Jianbing: Kidding 3?")'); call goErr
          TRACEBACK; status=1; return
        end select

        
        ! set interpolation indices:
        call lookup%InterpolSet( (/x,real(RInd),aimag(RInd)/), status )
        IF_NOTOK_RETURN(status=1)

        ! switch:
        ! Jianbing: here we have spaces for improvement: to discuss with Bas
        select case (i_m7mode)       
        ! fine bins !!! 
        case( inucs, iaits, iaccs, iaiti, iacci ) 
          ! apply interpolations in lookup tables:
          call lookup%InterpolApply( lut%ext_159, lu_ext, status )
          IF_NOTOK_RETURN(status=1)
          call lookup%InterpolApply( lut%a_159  , lu_a  , status )
          IF_NOTOK_RETURN(status=1)
          call lookup%InterpolApply( lut%g_159  , lu_g  , status )
          IF_NOTOK_RETURN(status=1)

        ! coarse bins 
        case( icoas, icoai )  
          ! apply interpolations in lookup tables:
          call lookup%InterpolApply( lut%ext_200, lu_ext, status )
          IF_NOTOK_RETURN(status=1)
          call lookup%InterpolApply( lut%a_200  , lu_a  , status )
          IF_NOTOK_RETURN(status=1)
          call lookup%InterpolApply( lut%g_200  , lu_g  , status )
          IF_NOTOK_RETURN(status=1)
          
        case default 
          write (gol,'("Jianbing: new mode???")'); call goErr
          TRACEBACK; status=1; return
        end select 

        ! fill:
        k_ext(i_m7mode)  = lu_ext * lambda * lambda * NumberConcentration(i_m7mode)
        k_sca(i_m7mode)  = lu_a * k_ext(i_m7mode)
        AF(i_m7mode)     = lu_g


        
        ! update sums:
        ext(jband)  = ext(jband)  + k_ext(i_m7mode)
        scattering  = scattering  + k_sca(i_m7mode)
        asym(jband) = asym(jband) + k_sca(i_m7mode)*AF(i_m7mode)

        
      end do ! loop imode
      
      ! write (gol,'("Jianbing to see nan test ")'); call goPr      
 
      ! if( jband == 10 ) then 
      !   ! write (gol,'("Jianbing to see nan test ")'); call goPr      
      !   call test_nan_1( ext(10) )
      ! end if 

      ! Finally determine extinction coefficient, AOD (tau), 
      ! single scattering coefficient ssa and asymmetry factor asy      
      ! should be altered when there are more than 2 modes
      if ( ext(jband) > 1.0d-7 ) then
        scat(jband)  = scattering/ ext(jband)
        if (scat(jband) > 0.000000001) then
          asym(jband)  = asym(jband) / scattering
        else 
          asym(jband)  = 1.0
        end if
      else 
        scat(jband)  = 0.0
        asym(jband)  = 1.0
      end if

      ! if ( jband == 10 ) then 
      !    if ( ext(10) > 0.0 ) then 
      !       jband = jband
      !    else 
      !     write (gol,'("Jianbing to see strange value )'); call goErr
      !     TRACEBACK; status=1; return
      !   end if 
      ! end if        

      ! if ( jband == 10 ) then ! the output band
      !     if ( ext(10) > 1.544d-3 ) then 
      !       write (gol,'("Jianbing to see the large extinction factor   ........................ ", es12.4)') ext(10); call goPr          
      !       do i_m7mode = 3, 4
      !         write (gol,'("Jianbing to see the dry radius is ", es12.4, "m for mode", i2)')   rad_dry_m7(i_m7mode),              i_m7mode; call goPr
      !         write (gol,'("Jianbing to see the wet radius is ", es12.4, "m for mode", i2)')   rad_wet_m7(i_m7mode),              i_m7mode; call goPr
      !         write (gol,'("Jianbing to see the corrected wet radius is ", es12.4, "m for mode", i2)')   rad_wet_m7_correct(i_m7mode),  i_m7mode; call goPr
      !         write (gol,'("Jianbing to see the num conc is ", es12.4, " per m3 for mode", i2)') Numberconcentration(i_m7mode),   i_m7mode; call goPr
      !         write (gol,'("Jianbing to see the mass conc is ", es12.4, " ug per m3 sea salt accumulation")')  mass_c(i_ssas); call goPr
      !         write (gol,'("Jianbing to see the k_ext is ", es12.4 )')  k_ext(i_m7mode); call goPr              
      !         write (gol,'("Jianbing to see the Refrac real and img are ", es12.4, es12.3 )')  real(RInd), aimag(RInd); call goPr
      !         write (gol,'("Jianbing to see the x  ", es12.4 )') x; call goPr
      !         write (gol,'("Jianbing to see the lamda  ", es12.4, "m" )') lambda; call goPr
      !       end do
      !     end if  ! if we have large ext
      ! end if ! if jband == 10

    end do ! loop j over swbands    

    ! write (gol,'("Jianbing to see here 10 ")'); call goPr

    ! clear:
    deallocate( RefracInd, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

    ! write (gol,'("Jianbing test mie m7 ends")'); call goPr

    ! ! jianbing stop here
    ! status = 1

  end subroutine calc_properties_mie_m7
#endif



  !-----------------------------------------------------------------

  subroutine rad_interpol(y,x,u,v)

    real, intent(in)    :: x(:), y(:)   ! (M) values and data
    real, intent(in)    :: u(:)         ! (N) target values
    real, intent(out)   :: v(:)         ! (N) interpolated data

    integer :: N, M, i, j

    N = size(x)
    M = size(u)

    !print *, 'n,m,u',n,m,u
    !print *, 'x1,xn,', x(1),x(N)

    ! loop over output values:
    do j = 1, M

      ! set to first data if below minimum:
      if ( u(j) < x(1) ) then
        ! set to first value:
        v(j) = y(1)
        ! next target value:
        cycle
      end if

      ! set to last data if above maximum:
      if ( u(j) > x(N) ) then
        ! set to last value:
        v(j) = y(N)
        ! warning ...
        !print *,'piep', v(j), y(N), j,N
        ! next target value:
        cycle
      end if

      ! loop over begin values of intervals:
      do i = 1, N-1
        ! target value in interval?
        if (u(j)>=x(i) .and. u(j)<x(i+1)) then
          ! interplate:
          v(j) = y(i) + (u(j)-x(i)) * (y(i+1)-y(i))/(x(i+1)-x(i))
          ! found, next target value:
          exit
        end if
      end do ! intervals

    end do ! output values

  end subroutine rad_interpol

end module LE_Radiation_Mie





