!###############################################################################
!
! This module contains several routines that set the indices number for the model species
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

module indices

  use GO     , only : gol, goErr, goPr

  use Binas, only : xm_H, xm_C, xm_N, xm_O, xm_S, xm_F
  use Binas, only : xm_seasalt, xm_Na, xm_Ca, xm_K, xm_Mg, xm_Pb, xm_Cd
  use Binas, only : xm_dummy

#ifdef with_m7
  use mo_aero_m7, only : m7_nmod => nmod
#endif

  implicit none


  ! --- in/out -----------------------------------

  public

  private  ::  mname


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'indices'

  !
  ! accumulated species
  !

  ! conversion factor from Na mass to seasalt:
  real,parameter  ::  Na_to_SeaSalt = xm_seasalt / xm_Na


  !
  ! Generated tracer info:
  !
  !  ! number of tracers, including special:
  !  integer, parameter  ::  nspec_all = 59
  !
  !  ! number of tracers, excluding special:
  !  integer, parameter  ::  nspec = 48
  !
  !  ! tracer indices;
  !  ! value -1 for disabled tracers:
  !  integer, parameter  ::  ispec_NO2      =   1  ! Nitrogen dioxide
  !  integer, parameter  ::  ispec_NO       =   2  ! Nitric oxide
  !    :
  !  integer, parameter  ::  ispec_CRO      =  33  ! Methylphenoxy radical
  !  integer, parameter  ::  ispec_BZA      =  -1  ! Benzaldehyde
  !  integer, parameter  ::  ispec_PHEN     =  -2  ! Phenol
  !    :
  !
  !  ! tracer indices as variables, for backwards compatibility:
  !  integer             ::  i_NO2      =   1  ! Nitrogen dioxide
  !  integer             ::  i_NO       =   2  ! Nitric oxide
  !    :
  !
  !  ! tracer names:
  !  character(len=*), parameter  ::  specname(nspec) = (/ &
  !          'no2   ', &    ! NO2
  !          'no    ', &    ! NO
  !            :
  !
  !  ! tracer units:
  !  character(len=*), parameter  ::  specname(nspec) = (/ &
  !          'ppb  ', &    ! NO2
  !          'ppb  ', &    ! NO
  !            :
  !
  !  ! tracer mole mass:
  !  real, parameter  ::  specmolm(nspec) = (/ &
  !       xm_O * 2 + xm_N, &     ! NO2
  !       xm_O + xm_N, &     ! NO
  !            :
  !
  !  ! aerosol densities:
  !  real, parameter  ::  aerdens(nspec_all) = (/ &
  !          -999.9e3, &    ! [kg/m3]  NO2
  !          -999.9e3, &    ! [kg/m3]  NO
  !            :
  !
  !  ! aerosol geometric radius [m]:
  !  real, parameter  ::  aerosol_radius_g(nspec_all) = (/ &
  !          -999.9   , &    ! [m]  NO2
  !          -999.9   , &    ! [m]  NO
  !            :
  !
  !  ! aerosol geometric std.dev. [1]:
  !  real, parameter  ::  aerosol_sigma_g(nspec_all) = (/ &
  !          -999.9, &    ! [1]  NO2
  !          -999.9, &    ! [1]  NO
  !            :
  !
  !  ! aerosol hygroscopicities:
  !  real, parameter  ::  aerhygro(nspec_all) = (/ &
  !          -999.9, &    ! NO2
  !          -999.9, &    ! NO
  !            :
  !
  !  ! aerosol modes:
  !  integer, parameter              ::  NO_AEROSOL          = 0
  !  integer, parameter              ::  AEROSOL_FINE_MODE   = 1
  !  integer, parameter              ::  AEROSOL_COARSE_MODE = 2
  !  integer, parameter              ::  AEROSOL_ALL_MODES   = 3
  !  ! per tracer:
  !  integer, parameter  ::  specmode(nspec) = (/ &
  !      NO_AEROSOL         , &    ! NO2
  !      NO_AEROSOL         , &    ! NO
  !          :
  !      AEROSOL_FINE_MODE  , &    ! SO4a_f
  !      AEROSOL_COARSE_MODE, &    ! SO4a_c
  !          :
  !
  !  ! number of advected tracers:
  !  integer, parameter  ::  n_advected = 40
  !  ! list of tracer indices:
  !  integer, parameter  ::  ispecs_advected(n_advected) = (/ &
  !      ispec_NO2     , &   ! NO2
  !      ispec_NO      , &   ! NO
  !      ispec_O3      , &   ! O3
  !          :
  !
  !  ! number of sia tracers:
  !  integer, parameter  ::  n_sia = 5
  !  ! list of tracer indices:
  !  integer, parameter  ::  ispecs_sia(n_sia) = (/ &
  !      ispec_SO4a_f  , &   ! SO4a_f
  !      ispec_SO4a_c  , &   ! SO4a_c
  !          :
  !
  !  ! number of dust tracers:
  !  integer, parameter  ::  n_dust = 2
  !  ! list of tracer indices:
  !  integer, parameter  ::  ispecs_dust(n_dust) = (/ &
  !      ispec_dust_f  , &   ! dust_f
  !      ispec_dust_c   /)   ! dust_c
  !
  !  ! number of seasalt tracers:
  !  integer, parameter  ::  n_seasalt = 2
  !  ! list of tracer indices:
  !  integer, parameter  ::  ispecs_seasalt(n_seasalt) = (/ &
  !      ispec_Na_f    , &   ! Na_f
  !      ispec_Na_c     /)   ! Na_c
  !
  !  ! number of basecation tracers:
  !  integer, parameter  ::  n_basecation = 2
  !  ! list of tracer indices:
  !  integer, parameter  ::  ispecs_basecation(n_basecation) = (/ &
  !      ispec_Na_f    , &   ! Na_f
  !      ispec_Na_c     /)   ! Na_c
  !
  !  ! number of hm tracers:
  !  integer, parameter  ::  n_hm = 0
  !  ! dummy list of tracer indices:
  !  integer, parameter  ::  ispecs_hm(1) = (/ -1 /)
  !
  !  ! number of cg tracers:
  !  integer, parameter  ::  n_cg = 0
  !  ! dummy list of tracer indices:
  !  integer, parameter  ::  ispecs_cg(1) = (/ -1 /)
  !
  !  ! number of soa tracers:
  !  integer, parameter  ::  n_soa = 0
  !  ! dummy list of tracer indices:
  !  integer, parameter  ::  ispecs_soa(1) = (/ -1 /)
  !
  !  ! number of soa_prec tracers:
  !  integer, parameter  ::  n_soa_prec = 0
  !  ! dummy list of tracer indices:
  !  integer, parameter  ::  ispecs_soa_prec(1) = (/ -1 /)
  !
  !  ! number of pops tracers:
  !  integer, parameter  ::  n_pops = 0
  !  ! dummy list of tracer indices:
  !  integer, parameter  ::  ispecs_pops(1) = (/ -1 /)
  !
  !
  include "le_indices.inc"

  !
  ! Number of components, original tracer index per component,
  ! and weights of components in total sum:
  integer, allocatable  ::  accum_n(:)      ! (ispec)
  integer, allocatable  ::  accum_ii(:,:)   ! (ispec,1:n)
  real, allocatable     ::  accum_ww(:,:)   ! (ispec,1:n)
  !
  ! Factor used for conversion from ppb to ug/m3 :
  !   ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) (mole air)/(kg air) (kg air)/(m3 air)
  !    c *         1e-9           *     xm_tracer   *     1e9       /       xm_air      *     dens
  ! thus:
  !    c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
  ! Define air molar concentration:
  !    air_molar_conc = air_mass_conc / air_molar_mass
  !      [mole/m3]         [kg/m3]    /   [kg/mole]
  ! In terms of LE variables for a single grid cell:
  !    air_molar_conc =  dens(i,j,k) / xm_air
  ! Define extra flag for conversion:
  logical, allocatable  ::  accum_ppb_to_ugm3(:,:)  ! (ispec,1:n)
  ! If this flag is true, the weight factor in accum_ww is set to the required tracer molar mass.
  !
  ! Example:
  !   c_tot = 0.0
  !   do k = 1, accum_n(i_tot)
  !     if ( accum_ppb_to_ugm3(i_tot,k) ) then
  !       c_tot = c_tot + c(accum_ii(i_tot,k)) * accum_ww(i_tot,k) * dens/xm_air
  !     else
  !       c_tot = c_tot + c(accum_ii(i_tot,k)) * accum_ww(i_tot,k)
  !     end if
  !   end do

#ifdef with_m7
  ! special M7 arrays:
  !  mode number following parameters in "mo_aero_m7" :
  !         inucs=1,  iaits=2,  iaccs=3,  icoas=4,  iaiti=5,  iacci=6,  icoai=7
  !         nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !         soluble | soluble | soluble | soluble | insol.  | insol.  | insol.  |
  integer, allocatable  ::  m7_mode(:)    ! (nspec_all)
  ! tracer index for number concentration for each M7 mode:
  integer, allocatable  ::  m7_numberconc_ispecs(:)    ! (m7_mod)
#endif
  
  ! number of vbs bins, used to define some tables ;
  ! the actual numbers 'n_vbs_cg' and 'n_vbs_soa' should be
  ! either both -1 or both equal to this value:
  ! not anymore; n_vbs_species > n_vbs_bins due to multiple precursor classes, rj
  integer, parameter    ::  n_vbs_bins = 9

  ! particle density used for sedimentation and M7 deposition,
  ! represents seasalt+water ;
  ! in future, this should be defined in the tracer table:
  real, parameter   ::  rhopart_default = 1.14e3 ! kg/m3

  ! dust not in water
  ! "Development of a global model of mineral dust aerosol microphysics", Lee et al, 2009
  real, parameter   ::  rhopart_dust    = 2.65e3 ! kg/m3
#ifdef with_pollen  
  ! Pollen are different
  real, parameter   ::  rhopart_pollen = 800.0 ! kg/m3
#endif  
  real, allocatable ::  rhopart(:)  ! nspec

contains


  ! ============================================================


  subroutine LE_Indices_Init( status )

    use Binas, only : xm_dummy

    ! --- in/out ------------------------------

    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Indices_Init'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! accumulated species; not in concentration array ...
    call Indices_Set_Accum( status )
    IF_NOTOK_RETURN(status=1)

#ifdef with_m7
    ! special M7 stuff ..
    call Indices_Set_M7( status )
    IF_NOTOK_RETURN(status=1)
#endif
    
    ! check VBS:
    if ( n_vbs_cg /= n_vbs_soa ) then
      write (gol,'("number of vbs cg tracers (",i6,") should be the same")') n_vbs_cg; call goErr
      write (gol,'("as number of vbs soa tracers (",i6,")")') n_vbs_soa; call goErr
      TRACEBACK; status=1; return
    end if
    !! VBS enabled ? !rj: no longer true when using specific cg's for precursor classes
    !if ( n_vbs_cg > 0 ) then
    !  ! check ...
    !  if ( n_vbs_cg /= n_vbs_bins ) then
    !    write (gol,'("number of vbs cg or soa tracers (",i6,") should be the same")') n_vbs_cg; call goErr
    !    write (gol,'("as number of vbs classes (",i6,")")') n_vbs_bins; call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !end if

    ! ok
    status = 0

  end subroutine LE_Indices_Init


  ! ***


  subroutine LE_Indices_Done( status )

    ! --- in/out ------------------------------

    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Indices_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! clear:
    deallocate( accum_n )
    deallocate( accum_ii )
    deallocate( accum_ww )
    deallocate( accum_ppb_to_ugm3 )
#ifdef with_m7    
    deallocate( m7_mode )
    deallocate( m7_numberconc_ispecs )
#endif
    ! ok
    status = 0

  end subroutine LE_Indices_Done


  ! ***


  subroutine Indices_Set_Accum( status )

    ! --- in/out ------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Indices_Set_Accum'

    ! --- local -------------------------------

    integer   ::  ispec, k, icomp
    integer   ::  ieqv

    ! --- begin -------------------------------

!#ifdef with_m7
!    ! warning ...
!    if ( n_m7 > 0 ) then
!      write (gol,'("Check definition of accumulated tracers with M7, some tracers might be double now.")'); call goErr
!      write (gol,'("Please do not add all kind of adhoc tests, just ensure that either the M7 or the fine/coarse tracers are used.")'); call goErr
!      TRACEBACK; status=1; return
!    end if
!#endif

    ! accumulated species:
    !   c_tot = 0.0
    !   do k = 1, accum_n(i_tot)
    !     c_tot = c_tot + c(accum_ii(i_tot,k)) * accum_ww(i_tot,k)
    !   end do
    ! number of components, original components, weights:
    allocate( accum_n (nspec_all) )
    allocate( accum_ii(nspec_all,nspec) )
    allocate( accum_ww(nspec_all,nspec) )
    allocate( accum_ppb_to_ugm3(nspec_all,nspec) )
    ! init to zero:
    accum_n = 0
    ! for safety, maybe program crashes if not correctly used ...
    accum_ii = -999
    accum_ww = -9.99e20
    ! by default no conversion:
    accum_ppb_to_ugm3 = .false.

    ! not-accumulated tracers:
    do ispec = 1, nspec
      accum_n (ispec  ) = 1
      accum_ii(ispec,1) = ispec
      accum_ww(ispec,1) = 1.0
    end do

    ! total pm [0,10um] :
    ispec = ispec_tpm10
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      !
      ! fine/coarse
      call Indices_Accum_Component( ispec, ispec_so4a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4a_c,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_nh4a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_no3a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_no3a_c,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_f  ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_c  ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_pom_f ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_pom_c ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_f ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_c ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f1  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f2  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c1  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c2  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c3  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f1,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f2,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c1,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c2,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c3,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
      ! ppm tracers used for OPS emissions
      call Indices_Accum_Component( ispec, ispec_ppm_1   ,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_1_25,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_25_4,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_4_10,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_1    ,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_1_25 ,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
#ifdef with_m7      
      ! M7 tracers
      call Indices_Accum_Component( ispec, ispec_so4ns, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4ks, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4as, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4cs, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bcks , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bcas , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bccs , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bcki , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocks , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocas , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_occs , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocki , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ssas , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_sscs , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_duas , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ducs , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_duai , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_duci , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
#endif
      !
      ! POA/SOA tracers
      call Indices_Accum_Component( ispec, ispec_vbs_poa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_poa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_poa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_poa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_poa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_poa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_poa7  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_poa8  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_poa9  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)

      call Indices_Accum_Component( ispec, ispec_vbs_sisoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa7  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa8  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      

      call Indices_Accum_Component( ispec, ispec_vbs_asoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_asoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_asoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_asoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_asoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_asoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      

      call Indices_Accum_Component( ispec, ispec_vbs_bsoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)                      
      !
    end if                      

    ! SIA PM2.5 fraction (CAMS50) :
    ispec = ispec_sia25
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      !
      ! fine/coarse
      call Indices_Accum_Component( ispec, ispec_so4a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_nh4a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_no3a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
    end if

#ifdef with_labeling
    ! total pm [0,10um] use for labeling to combine nh4no3 as 1 tracer :
    ispec = ispec_tpm10_comb
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      !
      ! fine/coarse
      call Indices_Accum_Component( ispec, ispec_so4a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4a_c,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_nh4a_f_no3a_f,    1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_nh4a_f_on_so4a_f, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_no3a_c,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_f  ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_c  ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_pom_f ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_pom_c ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_f ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_c ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f1  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f2  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c1  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c2  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c3  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f1,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f2,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c1,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c2,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c3,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
    end if
#endif

    ! total pm [0,25um]
    ispec = ispec_tpm25
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      !
      ! fine/coarse modes
      call Indices_Accum_Component( ispec, ispec_so4a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_nh4a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_no3a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_f  ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_pom_f ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_f ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f1  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f2  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f1,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f2,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
      ! ppm tracers used for OPS emissions
      call Indices_Accum_Component( ispec, ispec_ppm_1   ,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_1_25,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_1    ,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_1_25 ,         1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
#ifdef with_m7
      ! M7 tracers
      call Indices_Accum_Component( ispec, ispec_so4ns, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4ks, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4as, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bcks , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bcas , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bcki , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocks , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocas , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocki , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ssas , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_duas , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_duai , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
#endif      
      !
      ! POA/SOA tracers
      call Indices_Accum_Component( ispec, ispec_vbs_poa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_poa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_poa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_poa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_poa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_poa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_poa7  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_poa8  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_poa9  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)

      call Indices_Accum_Component( ispec, ispec_vbs_sisoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa7  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_sisoa8  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)

      call Indices_Accum_Component( ispec, ispec_vbs_asoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_asoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_asoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_asoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_asoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_asoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)

      call Indices_Accum_Component( ispec, ispec_vbs_bsoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_vbs_bsoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      !
    end if

#ifdef with_labeling
    ! total pm [0,25um] use for labeling to combine nh4no3 as 1 tracer :
    ispec = ispec_tpm25_comb
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      !
      ! fine/coarse modes
      call Indices_Accum_Component( ispec, ispec_so4a_f,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_nh4a_f_no3a_f,    1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_nh4a_f_on_so4a_f, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_f  ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_pom_f ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ppm_f ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f1  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f2  , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f1,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f2,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
    end if
#endif      
    ! total sea-salt:
    ispec = ispec_tss
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      !
      ! fine/coarse
      call Indices_Accum_Component( ispec, ispec_na_f1 , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_f2 , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c1 , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c2 , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_na_c3 , Na_to_SeaSalt, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
#ifdef with_m7      
      ! M7 aerosols
      call Indices_Accum_Component( ispec, ispec_ssas ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_sscs ,           1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
#endif    
    end if

    ! total dust:
    ispec = ispec_tdust
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      !
      ! fine/coarse
      call Indices_Accum_Component( ispec, ispec_dust_f1, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_f2, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c1, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c2, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_dust_c3, 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
#ifdef with_m7      
      ! M7 aerosols
      call Indices_Accum_Component( ispec, ispec_duas  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ducs  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_duai  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_duci  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
#endif
    end if

    ! total carbon [0,10um] :
    ispec = ispec_tc
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      !
      ! fine/coarse
      call Indices_Accum_Component( ispec, ispec_ec_f  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_c  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_pom_f , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_pom_c , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
      ! EC components used for coupling with OPS
      call Indices_Accum_Component( ispec, ispec_ec_1    , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ec_1_25 , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !
#ifdef with_m7      
      ! M7 aerosols
      call Indices_Accum_Component( ispec, ispec_bcks  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bcas  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bccs  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_bcki  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocks  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocas  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_occs  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_ocki  , 1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
#endif      
    end if

    ! total soa
    ispec = ispec_tsoa
    if ( ispec > 0 ) then
      accum_n(ispec) = 0

      call Indices_Accum_Component( ispec, i_vbs_sisoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa7  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa8  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)

      call Indices_Accum_Component( ispec, i_vbs_asoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)

      call Indices_Accum_Component( ispec, i_vbs_bsoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
    end if

    ! Total CG (SVOC) 
    ! ~ convert from ppb [(1e-9 mole tracer)/(mole air)] to ug/(mole air)
    !   using mole masses of references carbon bonds ;
    !   enable 'ppb_to_ugm3' to have multiplication with  dens/xm_air [(mole air)/m3] :
    !      output = sum  conc(ispc)  specmolm      dens/xm_air
    !       ug/m3           ppb     kg/(mole tr)  (mole air)/m3
    ispec = ispec_tcg
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      do k = 1, n_vbs_cg
        icomp = ispecs_vbs_cg(k)
        ! check to avoid bound error:
        if ( icomp > 0 ) then
          !                                           kg/(mole tracer), unit conversion flag
          call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
          IF_NOTOK_RETURN(status=1)
        end if
      end do
    end if

    ! Total POG
    ! ~ convert from ppb [(1e-9 mole tracer)/(mole air)] to ug/(mole air)
    !   using mole masses of references carbon bonds ;
    !   enable 'ppb_to_ugm3' to have multiplication with  dens/xm_air [(mole air)/m3] :
    !      output = sum  conc(ispc)  specmolm      dens/xm_air
    !       ug/m3           ppb     kg/(mole tr)  (mole air)/m3
    ispec = ispec_tpog
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      do k = 1, n_vbs_pog
        icomp = ispecs_vbs_pog(k)
        ! check to avoid bound error:
        if ( icomp > 0 ) then
          !                                           kg/(mole tracer), unit conversion flag
          call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
          IF_NOTOK_RETURN(status=1)
        end if
      end do
    end if

    ! Total SISOG
    ! ~ convert from ppb [(1e-9 mole tracer)/(mole air)] to ug/(mole air)
    !   using mole masses of references carbon bonds ;
    !   enable 'ppb_to_ugm3' to have multiplication with  dens/xm_air [(mole air)/m3] :
    !      output = sum  conc(ispc)  specmolm      dens/xm_air
    !       ug/m3           ppb     kg/(mole tr)  (mole air)/m3
    ispec = ispec_tsisog
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      do k = 1, n_vbs_sisog
        icomp = ispecs_vbs_sisog(k)
        write(*,*) 'accum sisog', ispec, icomp, specmolm(icomp)
        ! check to avoid bound error:
        if ( icomp > 0 ) then
          !                                           kg/(mole tracer), unit conversion flag
          call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
          IF_NOTOK_RETURN(status=1)
        end if
      end do
    end if

    ! Total ASOG
    ! ~ convert from ppb [(1e-9 mole tracer)/(mole air)] to ug/(mole air)
    !   using mole masses of references carbon bonds ;
    !   enable 'ppb_to_ugm3' to have multiplication with  dens/xm_air [(mole air)/m3] :
    !      output = sum  conc(ispc)  specmolm      dens/xm_air
    !       ug/m3           ppb     kg/(mole tr)  (mole air)/m3
    ispec = ispec_tasog
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      do k = 1, n_vbs_asog
        icomp = ispecs_vbs_asog(k)
        write(*,*) 'accum asog', ispec, icomp, specmolm(icomp)
        ! check to avoid bound error:
        if ( icomp > 0 ) then
          !                                           kg/(mole tracer), unit conversion flag
          call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
          IF_NOTOK_RETURN(status=1)
        end if
      end do
    end if

    ! Total BSOG
    ! ~ convert from ppb [(1e-9 mole tracer)/(mole air)] to ug/(mole air)
    !   using mole masses of references carbon bonds ;
    !   enable 'ppb_to_ugm3' to have multiplication with  dens/xm_air [(mole air)/m3] :
    !      output = sum  conc(ispc)  specmolm      dens/xm_air
    !       ug/m3           ppb     kg/(mole tr)  (mole air)/m3
    ispec = ispec_tbsog
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      do k = 1, n_vbs_bsog
        icomp = ispecs_vbs_bsog(k)
        write(*,*) 'accum bsog', ispec, icomp, specmolm(icomp)
        ! check to avoid bound error:
        if ( icomp > 0 ) then
          !                                           kg/(mole tracer), unit conversion flag
          call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
          IF_NOTOK_RETURN(status=1)
        end if
      end do
    end if

!    ! Total POA
!    ispec = ispec_tpoa
!    if ( ispec > 0 ) then
!      accum_n(ispec) = 0
!      do k = 1, n_vbs_poa
!        icomp = ispecs_vbs_poa(k)
!        write(*,*) 'accum poa', ispec, icomp
!        ! check to avoid bound error:
!        if ( icomp > 0 ) then
!          !                                           kg/(mole tracer), unit conversion flag
!          !call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
!          call Indices_Accum_Component( ispec, i_comp  , 1.0, .false., status   )
!          IF_NOTOK_RETURN(status=1)
!        end if
!      end do
!    end if
!
!    ! Total SISOA
!    ispec = ispec_tsisoa
!    if ( ispec > 0 ) then
!      accum_n(ispec) = 0
!      do k = 1, n_vbs_sisoa
!        icomp = ispecs_vbs_sisoa(k)
!        write(*,*) 'accum sisoa', ispec, icomp
!        ! check to avoid bound error:
!        if ( icomp > 0 ) then
!          !                                           kg/(mole tracer), unit conversion flag
!          !call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
!          call Indices_Accum_Component( ispec, i_comp  , 1.0, .false., status   )
!          IF_NOTOK_RETURN(status=1)
!        end if
!      end do
!    end if
!
!    ! Total ASOA 
!    ispec = ispec_tasoa
!    if ( ispec > 0 ) then
!      accum_n(ispec) = 0
!      do k = 1, n_vbs_asoa
!        icomp = ispecs_vbs_asoa(k)
!        write(*,*) 'accum asoa', ispec, icomp
!        ! check to avoid bound error:
!        if ( icomp > 0 ) then
!          !                                           kg/(mole tracer), unit conversion flag
!          !call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
!          call Indices_Accum_Component( ispec, i_comp  , 1.0, .false., status   )
!          IF_NOTOK_RETURN(status=1)
!        end if
!      end do
!    end if
!
!    ! Total BSOA 
!    ispec = ispec_tbsoa
!    if ( ispec > 0 ) then
!      accum_n(ispec) = 0
!      do k = 1, n_vbs_bsoa
!        icomp = ispecs_vbs_bsoa(k)
!        write(*,*) 'accum bsoa', ispec, icomp
!        ! check to avoid bound error:
!        if ( icomp > 0 ) then
!          !                                           kg/(mole tracer), unit conversion flag
!          !call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
!          call Indices_Accum_Component( ispec, i_comp  , 1.0, .false., status   )
!          IF_NOTOK_RETURN(status=1)
!        end if
!      end do
!    end if

    ! total poa
    ispec = ispec_tpoa
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      call Indices_Accum_Component( ispec, i_vbs_poa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_poa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_poa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_poa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_poa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_poa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_poa7  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_poa8  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_poa9  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
    end if

    ! total sisoa
    ispec = ispec_tsisoa
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      call Indices_Accum_Component( ispec, i_vbs_sisoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa7  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_sisoa8  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
    end if

    ! total asoa
    ispec = ispec_tasoa
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      call Indices_Accum_Component( ispec, i_vbs_asoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_asoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
    end if

    ! total bsoa
    ispec = ispec_tbsoa
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      call Indices_Accum_Component( ispec, i_vbs_bsoa1  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa2  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa3  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa4  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa5  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, i_vbs_bsoa6  , 1.0, .false., status   )
      IF_NOTOK_RETURN(status=1)
    end if

    ! SOx as SO4a equivalent
    ispec = ispec_SOxa
    if ( ispec > 0 ) then
      ! init counter:
      accum_n(ispec) = 0
      ! check here already, if not present then no need to check on SO4a ..
      if ( ispec_SO2 > 0 ) then
        ! equivalent tracer:
        ieqv = i_SO4a_f
        ! check ...
        if ( ieqv <= 0 ) then
           write (gol,'("SO2 enabled, but need SO4a too for weight in accumulated tracer")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! add gas phase, convert from ppb to ug/m3:
        call Indices_Accum_Component( ispec, ispec_SO2   , specmolm(ieqv),  .true., status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! fine/coarse
      call Indices_Accum_Component( ispec, ispec_SO4a_f,            1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_SO4a_c,            1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
#ifdef with_m7      
      ! M7 aerosols
      call Indices_Accum_Component( ispec, ispec_so4ns ,            1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4ks ,            1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4as ,            1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_so4cs ,            1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
#endif      
    end if

    ! NOy as NO3a equivalent
    ispec = ispec_NOya
    if ( ispec > 0 ) then
      ! init counter:
      accum_n(ispec) = 0
      ! check here already, if not present then no need to check on NO3a ..
      if ( any( (/ispec_HNO3,ispec_NO2,ispec_NO,ispec_NO3,ispec_N2O5/) > 0 ) ) then
        ! equivalent tracer:
        ieqv = i_NO3a_f
        ! check ...
        if ( ieqv <= 0 ) then
          write (gol,'("one of HNO3/NO2/NO/NO3/N2O5 enabled, but need NO3a too for weight in accumulated tracer")'); call goErr
        TRACEBACK; status=1; return
        end if
        ! add gas phase, convert from ppb to ug/m3:
        call Indices_Accum_Component( ispec, ispec_HNO3  ,       specmolm(ieqv),  .true., status )
        IF_NOTOK_RETURN(status=1)
        call Indices_Accum_Component( ispec, ispec_NO2   ,       specmolm(ieqv),  .true., status )
        IF_NOTOK_RETURN(status=1)
        call Indices_Accum_Component( ispec, ispec_NO    ,       specmolm(ieqv),  .true., status )
        IF_NOTOK_RETURN(status=1)
        call Indices_Accum_Component( ispec, ispec_NO3   ,       specmolm(ieqv),  .true., status )
        IF_NOTOK_RETURN(status=1)
        call Indices_Accum_Component( ispec, ispec_N2O5  , 2.0 * specmolm(ieqv),  .true., status ) ! factor 2 because 2 N-molecules in N2O5
        IF_NOTOK_RETURN(status=1)
      end if
      ! add aerosols:
      call Indices_Accum_Component( ispec, ispec_NO3a_f,                  1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_NO3a_c,                  1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! NHx as NH4a equivalent
    ispec = ispec_NHxa
    if ( ispec > 0 ) then
      ! set components:
      accum_n(ispec) = 0
      ! check here already, if not present then no need to check on SO4a ..
      if ( ispec_NH3 > 0 ) then
        ! equivalent tracer:
        ieqv = i_NH4a_f
        ! check ...
        if ( ieqv <= 0 ) then
          write (gol,'("NH3 enabled, but need NH4a too for weight in accumulated tracer")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! add gas phase, convert from ppb to ug/m3:
        call Indices_Accum_Component( ispec, ispec_NH3   , specmolm(ieqv),  .true., status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! add aerosols:
      call Indices_Accum_Component( ispec, ispec_NH4a_f,            1.0, .false., status )
      IF_NOTOK_RETURN(status=1)
      !call Indices_Accum_Component( ispec, ispec_NH4a_c,            1.0, .false., status )
      !IF_NOTOK_RETURN(status=1)
    end if

    ! Total NMVOC
    ! ~ convert from ppb [(1e-9 mole tracer)/(mole air)] to ug/(mole air)
    !   using mole masses of references carbon bonds ;
    !   enable 'ppb_to_ugm3' to have multiplication with  dens/xm_air [(mole air)/m3] :
    !      output = sum  conc(ispc)  specmolm      dens/xm_air
    !       ug/m3           ppb     kg/(mole tr)  (mole air)/m3
    ispec = ispec_tnmvoc
    if ( ispec > 0 ) then
      accum_n(ispec) = 0
      do k = 1, n_nmvoc
        icomp = ispecs_nmvoc(k)
        ! check to avoid bound error:
        if ( icomp > 0 ) then
          !                                           kg/(mole tracer)
          call Indices_Accum_Component( ispec, icomp, specmolm(icomp), .true., status )
          IF_NOTOK_RETURN(status=1)
        end if
      end do
    end if

    !
    ! NO2 as observed by common used NOx measuring devices.
    ! What these equipment reports as "NO2" actually includes
    ! also other NOy components.
    ! In (Lamsal et al., 2008) the following is proposed (eq. 1):
    !
    !   NO2_obs ~ NO2 + Alkyl-Nitrates + 0.95 PAN + 0.35 HNO3
    !
    ! We assume that this is in mole.
    ! No Alkyl-Nitrates in model yet, or hidden in XO2N.
    !
    ! Reference:
    !   Lamsal, L. N., et al, (2008), 
    !   Ground-level nitrogen dioxide concentrations inferred from the 
    !   satellite-borne Ozone Monitoring Instrument, 
    !   J. Geophys. Res., 113, D16308, doi:10.1029/2007JD009235.
    !
    ispec = ispec_NO2_obs
    if ( ispec > 0 ) then
      ! init counter:
      accum_n(ispec) = 0
      ! add contributions:
      call Indices_Accum_Component( ispec, ispec_NO2   , 1.00,  .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_PAN   , 0.95,  .false., status )
      IF_NOTOK_RETURN(status=1)
      call Indices_Accum_Component( ispec, ispec_HNO3  , 0.35,  .false., status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! info ...
    write (gol,'("---[accumulated species]---------------------------------------")'); call goPr
    do ispec = 1, nspec_all
      if ( accum_n(ispec) < 1 ) cycle
      if ( ispec <= nspec ) then
        write (gol,'("tracer  ",i4," ",a," ",i4)') ispec, trim(specname(ispec)), accum_n(ispec); call goPr
      else
        write (gol,'("special ",i4," ",a," ",i4)') ispec, trim(specname(ispec)), accum_n(ispec); call goPr
      end if
      if ( accum_n(ispec) > 1 ) then
        do k = 1, accum_n(ispec)
          icomp = accum_ii(ispec,k)
          write (gol,'("  comp ",2i4," ",a," ",f5.2)') k, icomp, trim(specname(icomp)), accum_ww(ispec,k); call goPr
        end do
      end if
    end do
    write (gol,'("---------------------------------------------------------------")'); call goPr

    ! ok
    status = 0

  end subroutine Indices_Set_Accum


  ! ***


  ! try if component should be added to accumulated species ;
  ! weight in sum is 'w'

  subroutine Indices_Accum_Component( ispec, icomp, w, ppb_to_ugm3, status )

    ! --- in/out ------------------------------

    integer, intent(in)     ::  ispec
    integer, intent(in)     ::  icomp
    real, intent(in)        ::  w
    logical, intent(in)     ::  ppb_to_ugm3
    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Indices_Accum_Component'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! component defined ?
    if ( icomp > 0 ) then
      ! increase counter:
      accum_n(ispec) = accum_n(ispec) + 1
      ! store index of component:
      accum_ii(ispec,accum_n(ispec)) = icomp
      ! store weight:
      accum_ww(ispec,accum_n(ispec)) = w
      ! store conversion flag:
      accum_ppb_to_ugm3(ispec,accum_n(ispec)) = ppb_to_ugm3
    end if

    ! ok
    status = 0

  end subroutine Indices_Accum_Component


#ifdef with_m7
  ! ***


  subroutine Indices_Set_M7( status )

    use mo_aero_m7, only : inucs, iaits, iaccs, icoas, iaiti, iacci, icoai

    ! --- in/out ------------------------------

    integer, intent(out)    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Indices_Set_M7'

    ! --- local -------------------------------

    integer   ::  ispec
    integer   ::  imode

    ! --- begin -------------------------------

    !  inucs=1,  iaits=2,  iaccs=3,  icoas=4,  iaiti=5,  iacci=6,  icoai=7
    !  nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
    !  soluble | soluble | soluble | soluble | insol.  | insol.  | insol.  |

    ! mode number of all species:
    allocate( m7_mode(nspec_all) )
    ! defaults:
    m7_mode = -1   ! no M7 mode number
    ! set mode number:
    do ispec = 1, nspec_all
      ! should be m7 ...
      if ( .not. tracer_is_m7(ispec) ) cycle
      ! should be aerosol:
      if ( specmode(ispec) == NO_AEROSOL_MODE ) cycle
      ! soluble or insoluble ?
      if ( tracer_is_soluble(ispec) ) then
        ! which mode ?
        if ( specmode(ispec) == AEROSOL_NUCL_MODE ) then
          imode = inucs
        else if ( specmode(ispec) == AEROSOL_AITKEN_MODE ) then
          imode = iaits
        else if ( specmode(ispec) == AEROSOL_ACCUM_MODE ) then
          imode = iaccs
        else if ( specmode(ispec) == AEROSOL_COARSE_MODE ) then
          imode = icoas
        else
          write (gol,'("unsupported soluble specmode ",i6," for ispec ",i6)') specmode(ispec), ispec; call goErr
          TRACEBACK; status=1; return
        end if
      else
        ! which mode ?
        if ( specmode(ispec) == AEROSOL_AITKEN_MODE ) then
          imode = iaiti
        else if ( specmode(ispec) == AEROSOL_ACCUM_MODE ) then
          imode = iacci
        else if ( specmode(ispec) == AEROSOL_COARSE_MODE ) then
          imode = icoai
        else
          write (gol,'("unsupported soluble specmode ",i6," for ispec ",i6)') specmode(ispec), ispec; call goErr
          TRACEBACK; status=1; return
        end if
      end if
      ! store:
      m7_mode(ispec) = imode
    end do

    ! number concentration tracers per M7 mode:
    allocate( m7_numberconc_ispecs(m7_nmod) )
    ! defaults:
    m7_numberconc_ispecs = -1   ! no M7 mode number
    ! set mode number:
    do ispec = 1, nspec_all
      ! should be m7 ...
      if ( .not. tracer_is_m7(ispec) ) cycle
      ! should be number concentration:
      if ( .not. tracer_is_numberconc(ispec) ) cycle
      ! get mode number:
      imode = m7_mode(ispec)
      ! set tracer index:
      m7_numberconc_ispecs(imode) = ispec
    end do

    ! ok
    status = 0

  end subroutine Indices_Set_M7
#endif

end module
