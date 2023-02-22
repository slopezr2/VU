!###############################################################################
!
! Compute radiation properties
!
! NOTES AM:
!  - do something with naerspec, number of aerosol species
!  - RefracIndex%PerSpecie(j,LE_IND_SO4) = cmplx(Re,Im)
!  - only works for 2 modes, loops over number of modes but sum over modes only for mode 1 and 2
!
! History
!
!   2012, Astrid Manders, TNO
!     Original implementation.
!
!   2015-06, Arjo Segers, TNO
!     Splitted module into sub-modules:
!       - swbands  : wavelength bands
!       - lut      : lookup table
!       - mie      : Mie calculations
!     Introduced interpolation in lookup table.
!     Fixed computation of volume mean diameter (mie module).
!
!   2019-11, Jianbing Jin, TNO
!     Introduced extra Angstrom fields for comparison with AERONET and MODIS products.
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

module LE_Radiation

  use GO                  , only : gol, goPr, goErr
  use Indices             , only : N_AEROSOL_MODES
  use Indices             , only : n_aerosol
  use Num                 , only : T_LookUp
  use LE_Radiation_SWBands, only : T_SWBands
  use LE_Radiation_LUT    , only : T_Radiation_LUT
  use LE_Radiation_Mie    , only : T_AerosolModes
  use LE_Radiation_Mie    , only : T_RefracIndex
 
  implicit none

  ! --- in/out -----------------------------------

  private

  public  ::  LE_Radiation_Init
  public  ::  LE_Radiation_Done
  public  ::  LE_Radiation_Calc
  
  !public  :: iswband_440nm, iswband_675nm, iswband_870nm, iswband_1020nm
  public  :: swbands
  public  :: tau
  public  :: angstrom, i_ang_aeronet, i_ang_modis, i_ang_polder
  public  :: extinction
  public  :: ssa
  public  :: asy
  public  :: refr_fine, refr_coarse, refr_t
  public  :: refr_3d_fine, refr_3d_coarse
  public  :: reff_fine, reff_coarse
  public  :: Ncolumn_fine, Ncolumn_coarse
  public  :: Ndens_fine, Ndens_coarse
  public  :: Vdens_fine, Vdens_coarse


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Radiation'
  
  ! angstrom bands:
  integer, parameter      ::  n_ang = 3
  integer, parameter      ::  i_ang_aeronet = 1
  integer, parameter      ::  i_ang_modis   = 2
  integer,parameter       ::  i_ang_polder   = 3
  ! angstrom wavelengts:
  real, parameter         ::  ang_wl(2,n_ang) = reshape( &
                                (/ 0.440, 0.870,  &       ! [um] AERONET
                                   0.470, 0.650 , &.      ! [um] MODIS
                                   0.443,0.875  /), &             ! [um] POLDER
                                (/2,n_ang/) )


  ! --- var --------------------------------------

  ! wavelength bands:
  type(T_SWBands)             ::  swbands

  ! lookup table:
  type(T_Radiation_LUT)       ::  lut
  type(T_LookUp)              ::  lookup

  ! aerosol properties:
  type(T_AerosolModes)        ::  aer
  ! refractive indices lookup table:
  type(T_RefracIndex)         ::  RefracIndex
  
  ! number of layers to compute: nlev or nlev_top
  integer                     ::  nlev_rad

  ! results
  real, allocatable           ::  tau        (:,:,:,:)  ! (nx,ny,nlev_rad,swbands%n)
  real, allocatable           ::  angstrom   (:,:,:)    ! (nx,ny,n_ang)
  
  real, allocatable           ::  extinction (:,:,:,:)  ! (nx,ny,nlev_rad,swbands%n)
  real, allocatable           ::  ssa        (:,:,:,:)  ! (nx,ny,nlev_rad,swbands%n)
  real, allocatable           ::  asy        (:,:,:,:)  ! (nx,ny,nlev_rad,swbands%n)
  complex, allocatable        ::  refr_fine  (:,:,:)       ! (nx,ny,swbands%n)
  complex, allocatable        ::  refr_coarse(:,:,:)       ! (nx,ny,swbands%n)
  complex, allocatable        ::  refr_3d_fine  (:,:,:,:)  ! (nx,ny,nlev_rad,swbands%n)
  complex, allocatable        ::  refr_3d_coarse(:,:,:,:)  ! (nx,ny,nlev_rad,swbands%n)
  complex, allocatable        ::  refr_t     (:,:,:)    ! (nx,ny,swbands%n)
  real, allocatable           ::  reff_fine  (:,:)      ! (nx,ny)
  real, allocatable           ::  reff_coarse(:,:)      ! (nx,ny)
  real, allocatable           ::  Ncolumn_fine  (:,:)   ! (nx,ny)
  real, allocatable           ::  Ncolumn_coarse(:,:)   ! (nx,ny)
  real, allocatable           ::  Ndens_fine  (:,:,:)   ! (nx,ny,nlev_rad) aerosol number density [#/m3], fine modes
  real, allocatable           ::  Ndens_coarse(:,:,:)   ! (nx,ny,nlev_rad) aerosol number density [#/m3], coarse modes
  real, allocatable           ::  Vdens_fine  (:,:,:)   ! (nx,ny,nlev_rad) aerosol volume density [m3/m3], fine modes
  real, allocatable           ::  Vdens_coarse(:,:,:)   ! (nx,ny,nlev_rad) aerosol volume density [m3/m3], coarse modes

 
  
contains 


  ! ====================================================================
  ! ===
  ! === module init/done
  ! ===
  ! ====================================================================


  subroutine LE_Radiation_Init( rcF, status )

    use GO            , only : TrcFile
    use dims          , only : nx, ny
    use LE_Data_Common, only : nlev, nlev_top
    
    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Radiation_Init'

    ! --- local ----------------------------------    
    logical     ::  with_top
    
    ! --- begin ----------------------------------
    
    ! also comppute for top layers?
    call rcF%Get( 'le.radiation.with_top', with_top, status )
    IF_NOTOK_RETURN(status=1)
    
    ! set number of layers:
    if ( with_top ) then
      nlev_rad = nlev_top
    else
      nlev_rad = nlev
    end if

    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! aerosol modes
    ! ~~~~~~~~~~~~~~~~~~~~~~~~

    ! init aerosol properties:
    call aer%Init( status )
    IF_NOTOK_RETURN(status=1)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! lookup table
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! read radiance lookup table:
    call lut%Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! init 3D lookup:
    call lookup%Init( 3, status )
    IF_NOTOK_RETURN(status=1)
    ! define axis:
    call lookup%SetAx( 1, lut%TabInd_Rg, status )  
    IF_NOTOK_RETURN(status=1)
    call lookup%SetAx( 2, lut%TabInd_Re, status )  
    IF_NOTOK_RETURN(status=1)
    call lookup%SetAx( 3, lut%TabInd_Im, status )  
    IF_NOTOK_RETURN(status=1)


    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! demo
    ! ~~~~~~~~~~~~~~~~~~~~~~~~

    ! first compute demo extinctions, write to table ;
    ! module variables used:
    !  - already defined above : lut, lookup, aer
    !  - temporary defined     : swbands, RefracIndex
    call LE_Radiation_Demo( rcF, status )
    IF_NOTOK_RETURN(status=1)
    

    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! swbands
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! init from rcfile:
    call swbands%Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! refactive indices
    ! ~~~~~~~~~~~~~~~~~~~~~~~~

    ! refraction indices computed from radiance lut:
    call RefracIndex%Init( lut, swbands, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! storage
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! storage:
    allocate( extinction(nx,ny,nlev_rad,swbands%n), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( tau       (nx,ny,nlev_rad,swbands%n), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( ssa       (nx,ny,nlev_rad,swbands%n), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( asy       (nx,ny,nlev_rad,swbands%n), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( refr_fine  (nx,ny,swbands%n), source=(0.0,0.0), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( refr_coarse(nx,ny,swbands%n), source=(0.0,0.0), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( refr_3d_fine  (nx,ny,nlev_rad,swbands%n), source=(0.0,0.0), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( refr_3d_coarse(nx,ny,nlev_rad,swbands%n), source=(0.0,0.0), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( refr_t     (nx,ny,swbands%n), source=(0.0,0.0), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( reff_fine     (nx,ny), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( reff_coarse   (nx,ny), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Ncolumn_fine  (nx,ny), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Ncolumn_coarse(nx,ny), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Ndens_fine    (nx,ny,nlev_rad), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Ndens_coarse  (nx,ny,nlev_rad), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Vdens_fine    (nx,ny,nlev_rad), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Vdens_coarse  (nx,ny,nlev_rad), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( angstrom      (nx,ny,n_ang), source=0.0, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    ! end
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! ok
    status = 0

  end subroutine LE_Radiation_Init


  ! ***


  subroutine LE_Radiation_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Radiation_Done'

    ! --- begin ----------------------------------

    ! clear:
    deallocate( extinction, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( tau, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( angstrom, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( ssa, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( asy, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( refr_fine, stat=status   )
    IF_NOTOK_RETURN(status=1)
    deallocate( refr_coarse, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( refr_3d_fine  , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( refr_3d_coarse, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( refr_t, stat=status      )
    IF_NOTOK_RETURN(status=1)
    deallocate( reff_fine, stat=status)
    IF_NOTOK_RETURN(status=1)
    deallocate( reff_coarse, stat=status)
    IF_NOTOK_RETURN(status=1)
    deallocate( Ncolumn_fine, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Ncolumn_coarse, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Ndens_fine  , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Ndens_coarse, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Vdens_fine  , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Vdens_coarse, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! done with refract index:
    call RefracIndex%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with aerosol properties:
    call aer%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with lookup:
    call lookup%Done( status )
    IF_NOTOK_RETURN(status=1)
    ! done with lookup table:
    call lut%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    call swbands%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Radiation_Done


  !----------------------------------------------------------------------  


  !
  ! Write table with extinctions per wavelength
  !

  subroutine LE_Radiation_Demo( rcF, status )

    use GO              , only : TrcFile, ReadRc
    use GO              , only : goGetFU, pathsep
    use Indices         , only : nspec, specname
    use Indices         , only : n_aerosol, ispecs_aerosol, NO_AEROSOL_MODE
    use LE_Radiation_Mie, only : calc_properties_mie

    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Radiation_Demo'

    ! --- local ----------------------------------
    
    type(T_SWBands)         ::  swb
    integer                 ::  i, n
    real                    ::  bwidth
    integer                 ::  ispec, i_aerosol, irh
    real, allocatable       ::  conc(:)   ! (nspec)
    real, allocatable       ::  c_ext(:), c_ssa(:), c_asy(:)  ! (nband)
    complex, allocatable    ::  refr_f(:)       ! (nband)
    complex, allocatable    ::  refr_c(:)       ! (nband)
    real, allocatable       ::  sizeparam(:,:)  ! (nband,nspec)
    real, allocatable       ::  extf_spec(:,:)  ! (nband,nspec)
    complex, allocatable    ::  refr_spec(:,:)  ! (nband,nspec)
    real                    ::  volf, volc
    real                    ::  volfine, volcoarse
    real                    ::  crosf, crosc
    real                    ::  Nfine, Ncoarse
    real                    ::  rh0
    integer                 ::  fu
    character(len=1024)     ::  outdir, fname
    character(len=256)      ::  fmt

    ! --- begin ----------------------------------
    
    ! 
    !~ testing extinctions: 300 - 1800 nm
    !   (plus eps to avoid problems with search for lambda
    !   that happends to be exactly on the edge)
    !
    ! n = 1500
    ! bwidth = 0.001  ! um
    ! call swb%Init( n+1, status )
    ! IF_NOTOK_RETURN(status=1)
    ! do i = 1, n+1
    !   call swb%SetBand( i, status, lambda=0.300+(i-1)*bwidth, width=bwidth )
    !   IF_NOTOK_RETURN(status=1)
    ! end do
    
    !
    !~ with wavebands defined in rcfile:
    !
    call swb%Init( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! refraction indices computed from radiance lut:
    call RefracIndex%Init( lut, swb, status )
    IF_NOTOK_RETURN(status=1)

    ! testing ...
    allocate( conc(nspec) )
    allocate( c_ext(swb%n) )
    allocate( c_ssa(swb%n) )
    allocate( c_asy(swb%n) )
    allocate( refr_f(swb%n) )
    allocate( refr_c(swb%n) )
    allocate( sizeparam(swb%n,nspec) )
    allocate( extf_spec(swb%n,nspec) )
    allocate( refr_spec(swb%n,nspec) )
    
    ! output dir:
    call ReadRc( rcF, 'le.output.outdir', outdir, status )
    IF_NOTOK_RETURN(status=1)
    ! output file:
    write (fname,'(3a)') trim(outdir), pathsep, 'aerosol-extinction.txt'

    ! file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open text file:
    open( unit=fu, file=trim(fname), form='formatted', iostat=status )
    if ( status /= 0 ) then
      write (gol,'("could not open aerosol extinction file")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! format:
    fmt = '(a20,f8.1,2000es12.3)'
    ! headers:
    write (fu,fmt) 'aerosol', 0.0, swb%lambda
    ! loop over aersol tracers:
    do i_aerosol = 1, n_aerosol
      ! global index:
      ispec = ispecs_aerosol(i_aerosol)
      ! dummy concentrations, only one tracer filled:
      conc = 0.0
      conc(ispec) = 1.0e6  ! 1 g/m3 as test
      ! humidities:
      do irh = 0, 5
        rh0 = irh * 20.0  ! %
        ! calc properties:
        call calc_properties_mie( conc, rh0, swb, lookup, lut, aer, RefracIndex, &
                                   c_ext, c_ssa, c_asy, &
                                   refr_f, refr_c , &
                                   volf  , volc   , &
                                   crosf , crosc  , &
                                   Nfine , Ncoarse, status, &
                                   sizeparam=sizeparam, extf_spec=extf_spec, refr_spec=refr_spec )
        IF_NOTOK_RETURN(status=1)
        ! info ...
        write (fu,fmt) trim(specname(ispec)), rh0, c_ext(:)
        write (fu,fmt) trim(specname(ispec)), -1.0,       sizeparam(:,ispec)
        !write (fu,fmt) trim(specname(ispec)), -2.0,  real(refr_spec(:,ispec))
        !write (fu,fmt) trim(specname(ispec)), -3.0, aimag(refr_spec(:,ispec))
        write (fu,fmt) trim(specname(ispec)), -4.0,       extf_spec(:,ispec)
        !! testing ...
        !exit
      end do  ! rh
      !! testing ...
      !exit
    end do  ! aerosols
    ! close:
    close( fu, iostat=status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( conc )
    deallocate( c_ext )
    deallocate( c_ssa )
    deallocate( c_asy )
    deallocate( refr_f )
    deallocate( refr_c )
    deallocate( sizeparam )
    deallocate( extf_spec )
    deallocate( refr_spec )

    ! done with refract index:
    call RefracIndex%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    call swb%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Radiation_Demo


  !----------------------------------------------------------------------  

  !
  ! calculates radiative properties of bulk aerosol, externally mixed, using mie theory and lookup-table
  !
  ! reference: TNO report-060-UT_2012-00508 and Aan de Brugh
  ! 
  ! input: aerosol concentrations, relative humidity, layer height
  ! output extinction coeffic, asymmetry factor, single scattering albedo, aerosol optical depth
  !

  subroutine LE_Radiation_Calc( c, status )

    use dims            , only : nx, ny
    use LE_Data_Common  , only : nlev, nlev_top
    use LE_Data         , only : LE_Data_GetPointer
    use LE_Bound        , only : caloft
    use Indices         , only : nspec, specname, ispecs_aerosol
#ifdef with_m7
    ! jianbing: only wet radius is needed
    use LE_M7_Data      , only : rwetm7modes, rdrym7modes    ![cm]
    use LE_M7_Data      , only : waterm7modes   ![g/cm3] to check with astrid?
    use LE_Radiation_Mie, only : calc_properties_mie_m7
#else
    use LE_Radiation_Mie, only : calc_properties_mie
#endif  

    ! --- in/out ---------------------------------

    real, intent(in)                ::  c(:,:,:,:)  ! (nx,ny,nlev,nspec)
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Radiation_Calc'

    ! --- local ----------------------------------

    real, pointer         ::  rh     (:,:,:)   ! (lon,lat,nlev)
    real, pointer         ::  delta_h(:,:,:)   ! (lon,lat,nlev)

    real, allocatable     ::  conc(:)  ! (nspec)
    complex, allocatable  ::  refr_f(:)
    complex, allocatable  ::  refr_c(:)
    real                  ::  volf, volc
    real                  ::  volfine, volcoarse, voltotal
    real                  ::  crosf, crosc
    real                  ::  crossfine
    real                  ::  crosscoarse
    real                  ::  Nfine, Ncoarse
    integer               ::  ix, iy, iz

    integer               ::  i_ang
    integer               ::  j_angstrom1, j_angstrom2
    real                  ::  lambda1, lambda2
    real                  ::  tau1, tau2

    ! testing ...
    !integer            ::  ix0, iy0

    ! --- begin ----------------------------------
    
    ! access meteo:
    call LE_Data_GetPointer( 'rh', rh     , status, check_units ='%', &
                               check_lbo=(/1,1,1/), check_ubo=(/nx,ny,nlev_rad/) )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dh', delta_h, status, check_units ='m', &
                               check_lbo=(/1,1,1/), check_ubo=(/nx,ny,nlev_rad/) )
    IF_NOTOK_RETURN(status=1)

    ! storage:
    allocate( conc(nspec), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( refr_f(swbands%n), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( refr_c(swbands%n), stat=status )
    IF_NOTOK_RETURN(status=1)

    !
    !loops: i=modes, j=wavelengths, ix, iy,iz lon, lat, lev
    refr_fine      = 0.0
    refr_coarse    = 0.0
    refr_3d_fine   = 0.0
    refr_3d_coarse = 0.0
    refr_t         = 0.0
    Ncolumn_fine   = 0.0
    Ncolumn_coarse = 0.0
    Ndens_fine     = 0.0
    Ndens_coarse   = 0.0
    Vdens_fine     = 0.0
    Vdens_coarse   = 0.0
    reff_fine      = 0.0
    reff_coarse    = 0.0
        
    ! loop over grid cells:
    do ix = 1, nx
      do iy = 1, ny

        ! init sums:
        volfine=0.
        volcoarse=0.
        crossfine=0.
        crosscoarse=0.

        ! loop over layers:
        do iz = 1, nlev_rad
        
          ! concentrations:
          if ( iz <= nlev ) then
            conc = c     (ix,iy,iz,:)
          else
            conc = caloft(ix,iy,iz,:)
          end if


#ifdef with_m7
          ! calculate Mie properties:
          call calc_properties_mie_m7( conc, rwetm7modes(ix, iy, iz, :), rdrym7modes(ix, iy, iz, :), waterm7modes(ix, iy, iz, :), rh(ix,iy,iz), &
                                    swbands, lookup, lut, RefracIndex, &
                                    extinction(ix,iy,iz,:), ssa(ix,iy,iz,:), asy(ix,iy,iz,:), & 
                                    status)
          IF_NOTOK_RETURN(status=1)

          refr_f = 0.0 
          refr_c = 0.0
          volf = 0.0
          volc = 0.0
          crosf = 0.0
          crosc = 0.0
          Nfine = 0.0
          Ncoarse = 0.0         
#else 
          ! calculate Mie properties:
          call calc_properties_mie( conc, rh(ix,iy,iz), &
                                     swbands, lookup, lut, aer, RefracIndex, &
                                     extinction(ix,iy,iz,:), ssa(ix,iy,iz,:), asy(ix,iy,iz,:),&
                                     refr_f, refr_c , &  ! refractive index (not volume normalized!)
                                     volf  , volc   , &  ! aerosol volume density [m3/m3]
                                     crosf , crosc  , &
                                     Nfine , Ncoarse, &  ! aerosol number density [#/m3]
                                     status )
          IF_NOTOK_RETURN(status=1)
#endif 

          ! AOD sum with height:     1/m                  m
          tau(ix,iy,iz,:) = extinction(ix,iy,iz,:) * delta_h(ix,iy,iz)    ! 1

          ! height integral:
          refr_fine  (ix,iy,:) = refr_fine  (ix,iy,:) + delta_h(ix,iy,iz) * refr_f
          refr_coarse(ix,iy,:) = refr_coarse(ix,iy,:) + delta_h(ix,iy,iz) * refr_c

          ! 3D refractive indices:
          if ( volf > 0.0 ) then
            refr_3d_fine  (ix,iy,iz,:) = refr_f / volf
          end if
          if ( volc > 0.0 ) then
            refr_3d_coarse(ix,iy,iz,:) = refr_c / volc
          end if

          ! aerosol volumne:
          volfine     = volfine     + delta_h(ix,iy,iz) * volf
          volcoarse   = volcoarse   + delta_h(ix,iy,iz) * volc
          ! ...
          crossfine   = crossfine   + delta_h(ix,iy,iz) * crosf
          crosscoarse = crosscoarse + delta_h(ix,iy,iz) * crosc
          ! ...
          Ncolumn_fine  (ix,iy) = Ncolumn_fine  (ix,iy) + Nfine   * delta_h(ix,iy,iz)   ! #/m2
          Ncolumn_coarse(ix,iy) = Ncolumn_coarse(ix,iy) + Ncoarse * delta_h(ix,iy,iz)   ! #/m2
          
          ! 3D
          Ndens_fine  (ix,iy,iz) =  Nfine    ! #/m3
          Ndens_coarse(ix,iy,iz) =  Ncoarse  ! #/m3
          Vdens_fine  (ix,iy,iz) =  volf     ! #/m3
          Vdens_coarse(ix,iy,iz) =  volc     ! #/m3

        end do !loop over lev
        
        ! *
        
        ! angstrom parameters:
        do i_ang = 1, n_ang

          ! indices of bands used for calculation:
          call swbands%FindBand( ang_wl(1,i_ang), j_angstrom1, status )
          IF_NOTOK_RETURN(status=1)
          call swbands%FindBand( ang_wl(2,i_ang), j_angstrom2, status )
          IF_NOTOK_RETURN(status=1)

          ! total AOD's columns for selected wavelengths:
          tau1 = sum( tau(ix,iy,:,j_angstrom1) )
          tau2 = sum( tau(ix,iy,:,j_angstrom2) )
          ! any aerosol at all ?
          if ( tau2 > 0.0 ) then
            ! band mids:
            lambda1 = swbands%lambda(j_angstrom1)
            lambda2 = swbands%lambda(j_angstrom2)
            ! compute Angstrom exponent:
            angstrom(ix,iy,i_ang) = - log( tau1 / tau2 ) / log( lambda1 / lambda2 )
          else
            ! no data:
            angstrom(ix,iy,i_ang) = -999.9
          end if
          
        end do  ! angstroms
        
        ! *

        ! refractive index of fine aerosols for each waveband:
        if ( volfine > 0.0 ) then
          refr_fine(ix,iy,:) = refr_fine(ix,iy,:) / volfine
        end if

        ! refractive index of coarse aerosols for each waveband:
        if ( volcoarse > 0.0 ) then
          refr_coarse(ix,iy,:) = refr_coarse(ix,iy,:) / volcoarse
        end if

        ! refractive index of total aerosol is volume weighted average:
        voltotal = volfine + volcoarse
        if ( voltotal > 0.0 ) then
          refr_t(ix,iy,:) = ( refr_fine(ix,iy,:) * volfine + refr_coarse(ix,iy,:) * volcoarse )/voltotal
        end if

        !! effective radius:
        !OLD:
        !if ( crossfine > 0.0 ) then
        !  reff_fine  (ix,iy) = (1.0/6.0) * volfine   / crossfine
        !end if
        !if ( crosscoarse > 0.0 ) then
        !  reff_coarse(ix,iy) = (1.0/6.0) * volcoarse / crosscoarse
        !end if
        !
        ! NEW: Following Aeronet definition:
        !           int r^3 dN   V / (4/3 pi)   3 V
        !   r_eff = ---------- = ------------ = - -
        !           int r^2 dN      C / pi      4 C
        if ( crossfine > 0.0 ) then
          reff_fine  (ix,iy) = (3.0/4.0) * volfine   / crossfine
        end if
        if ( crosscoarse > 0.0 ) then
          reff_coarse(ix,iy) = (3.0/4.0) * volcoarse / crosscoarse
        end if

      end do ! loop over lat
    end do ! loop over lon

    ! clear:
    deallocate( conc, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( refr_f, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( refr_c, stat=status )
    IF_NOTOK_RETURN(status=1)

    !! testing ...
    !print *, 'break after first LE_Radiation_Calc'
    !status=1; return

    ! ok
    status = 0

  end subroutine LE_Radiation_Calc

end module LE_Radiation




