!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Chem

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  LE_Chem_Init, LE_Chem_Done
  public  ::  LE_Chem_Step
  public  ::  Calc_r1r3


  ! --- const ------------------------------------

  ! module name:
  character(len=*), parameter  ::  mname = 'LE_Chem'


  ! --- var --------------------------------------

  ! correction factors for JNO2
  integer, parameter   ::  jno2_corr_nh = 10
  integer, parameter   ::  jno2_corr_na = 12
  real    ::      h_JNO2_corr(jno2_corr_nh)
  real    ::  angle_JNO2_corr(jno2_corr_na) = 0.0   ! <--- read if 0.0
  real    ::   corr_jno2(jno2_corr_nh,jno2_corr_na)
  character(len=512)   ::  JNO2_corr_file

  ! aerosol code:
  character(len=32)     ::  aero_code
  ! cloud code enabled ?
  logical               ::  with_chem_cloud

  ! minimum/maximum number of iterations for the chemistry
  integer   ::  miniter, maxiter

  ! timers:
  integer   ::  itim_chem_gas, itim_chem_cloud, itim_chem_aero


contains


  ! ====================================================================


  subroutine LE_Chem_Init( rcF, status )

    use GO               , only : TRcFile, ReadRc
    use GO               , only : GO_Timer_Def
    use LE_Chem_Work     , only : LE_Chem_Work_Init
    use LE_Chem_Hetro    , only : LE_Chem_Hetro_Init
    use LE_Chem_Cloud    , only : LE_Chem_Cloud_Init
    use LE_SIA_Isorropia2, only : LE_SIA_Isorropia2_Init
    use LE_Data          , only : LE_Data_Enable

    ! --- in/out ---------------------------------

    type(TRcFile), intent(in)           ::  rcF
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Chem_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! JNO2 correction table:
    call ReadRc( rcF, 'le.chem.JNO2_corr.file', JNO2_corr_file, status )
    IF_NOTOK_RETURN(status=1)

    ! SIA chemistry routine:
    call ReadRc( rcF, 'le.chem.sia.scheme', aero_code, status )
    IF_NOTOK_RETURN(status=1)

    ! with cloud chemistry ?
    call ReadRc( rcF, 'le.chem.cloud', with_chem_cloud, status )
    IF_NOTOK_RETURN(status=1)

    ! chemistry solver settings:
    call ReadRc( rcF, 'le.chem.solver.iter.min', miniter, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'le.chem.solver.iter.max', maxiter, status )
    IF_NOTOK_RETURN(status=1)

    ! define timers:
    call GO_Timer_Def( itim_chem_gas  , 'gas-phase chemistry', status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_chem_cloud, 'cloud chemistry'    , status )
    IF_NOTOK_RETURN(status=1)
    call GO_Timer_Def( itim_chem_aero , 'aerosol chemistry'  , status )
    IF_NOTOK_RETURN(status=1)

    ! setup gas-phase chemistry:
    call LE_Chem_Work_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! setup hetro chemistry:
    call LE_Chem_Hetro_Init( status )
    IF_NOTOK_RETURN(status=1)

    ! init cloud chemistry ?
    if ( with_chem_cloud ) then
      call LE_Chem_Cloud_Init( status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! switch:
    select case ( trim(aero_code) )
      case ( 'isorropia2' )
        call LE_SIA_Isorropia2_Init( status )
        IF_NOTOK_RETURN(status=1)
      case ( 'equisam' )
        ! nothing to be done
      case default
        write (gol,'("unsupported aerosol code : ",a)') trim(aero_code); call goErr
        TRACEBACK; status=1; return
    end select

    ! enable data:
    call LE_Data_Enable( 'p', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'q', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'occ', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Chem_Init


  ! ***


  subroutine LE_Chem_Done( status )

    use LE_Chem_Work     , only : LE_Chem_Work_Done
    use LE_Chem_Hetro    , only : LE_Chem_Hetro_Done
    use LE_Chem_Cloud    , only : LE_Chem_Cloud_Done
    use LE_SIA_Isorropia2, only : LE_SIA_Isorropia2_Done

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Chem_Done'

    ! --- begin ----------------------------------

    ! done with gas-phase chemistry:
    call LE_Chem_Work_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with hetro chemistry:
    call LE_Chem_Hetro_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with cloud chemistry:
    if ( with_chem_cloud ) then
      call LE_Chem_Cloud_Done( status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! switch:
    select case ( trim(aero_code) )
      case ( 'isorropia2' )
        call LE_SIA_Isorropia2_Done( status )
        IF_NOTOK_RETURN(status=1)
      case ( 'equisam' )
        ! nothing to be done
      case default
        write (gol,'("unsupported aerosol code : ",a)') trim(aero_code); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine LE_Chem_Done


  ! ***


  subroutine LE_Chem_Step( c, aerh2o, t ,deltat, lrestart, status )

    use Binas, only : h2odens => rol
    use Binas, only : xm_NH4HSO4, xm_SO4, rho_NH4HSO4a
    use Binas, only : xm_seasalt, rho_seasalt
    use Binas, only : xm_air,xm_h2o
    use Binas, only : Rgas, Avog
    use GO   , only : TDate
    use GO   , only : GO_Timer_Start, GO_Timer_End
    use LE_Logging, only : ident2
    !use LE_Logging, only : u_log_c
    use dims , only : zenith
    !use dims , only : h
    use dims , only : runF
    use indices
    use LE_SIA_Isorropia2   , only : LE_SIA_Isorropia2_Step
    use LE_SIA_Equisam      , only : LE_SIA_Equisam_Step
#ifdef with_vbs
    use LE_VBS              , only : i_soa_prec_par, i_soa_prec_ole, i_soa_prec_tol
    use LE_VBS              , only : i_soa_prec_xyl, i_soa_prec_iso, i_soa_prec_terp
    use LE_VBS              , only : with_vbs_soa_chemistry
    use LE_VBS              , only : LE_VBS_Apply
#endif

    use dims, only: outF, runF
    use dims, only: nx, ny, nz, nspec
    use dims, only: rk1rk3

    use LE_Chem_Hetro, only : nreac_het
    use LE_Chem_Hetro, only : ireac_N2O5_NH4HSO4a_f
    use LE_Chem_Hetro, only : ireac_N2O5_ss_f
    use LE_Chem_Hetro, only : ireac_N2O5_ss_c
    use LE_Chem_Hetro, only : ireac_HNO3_ss_f
    use LE_Chem_Hetro, only : ireac_HNO3_ss_c
    use LE_Chem_Hetro, only : LUT_Rk_Het
    use LE_Chem_Hetro, only : LUT_Rk_Het_Get
    use LE_Chem_Hetro, only : increment_factor_aerh2o
    use LE_Chem_Hetro, only : increment_factor_rh
    use LE_Chem_Work          , only : nreac
    use LE_Chem_Work          , only : naux
    use LE_Chem_Work          , only : iaux_air
    use LE_Chem_Work          , only : iaux_T, iaux_H2O, iaux_zen, iaux_Rh, iaux_p
    use LE_Chem_Work          , only : iaux_cloud, iaux_cldfac, iaux_cldsulf
    use LE_Chem_Work          , only : iaux_hetn2o5
    use LE_Chem_Work          , only : iaux_nox0
    use LE_Chem_Work          , only : iaux_ppm_to_mlccm3, iaux_ppb_to_mlccm3
    use LE_Chem_Work          , only : cldadj
#ifdef with_m7    
    use LE_Chem_Work          , only : ireac_M7_R71, ireac_M7_R72
#endif
    use LE_Chem_Work          , only : LE_Chem_Work_Rates
    use LE_Chem_Work          , only : LE_Chem_Work_Iter
    use LE_Chem_Solver_TwoStep, only : TwoStep
    use LE_Chem_Cloud         , only : wetS
    use LE_Data      , only : LE_Data_GetPointer

    ! --- in/out ---------------------------------

    real, intent(inout)       ::  c(nx,ny,nz,nspec)
    real, intent(inout)       ::  aerh2o(nx,ny,nz)
    type(TDate), intent(in)   ::  t
    real, intent(in)          ::  deltat
    logical, intent(in)       ::  lrestart
    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Chem_Step'

    ! air constant:
    real, parameter   ::  M = 2.5E19

    ! this is only used for chem_solve_nox (solver bdf1)
    !     -> set nox0 to high value
    real, parameter   ::  nox0 = 1.0e30

    ! conversion:
    real, parameter   ::  min2sec = 60.0  ! sec/min

    ! --- local ----------------------------------

    ! loop variables:
    integer     ::  i, j, k
    integer     ::  ispec
    integer     ::  status_par

    ! the reaction rates: gas phase
    real(8)     ::  rk(nreac)
    real(8)     ::  rk_het(nreac_het)
    real(8)     ::  aux(naux)
    integer     ::  icld
    real        ::  rho_air

    ! conversion mask:
    logical, allocatable  ::  ppbmask(:)

    ! hetrogeneous reacations:
    real        ::  finc
    real        ::  conc_aer

    ! concentration arrays:
    real        ::  ch(nspec)
    real        ::  Q(nspec)

    ! twostep solver:
    real        ::  ts
    real        ::  dtmin, dtmax, dtstart
    real        ::  atol(nspec), rtol(nspec)
    integer     ::  nfcn, naccpt, nrejec,nstart
    integer     ::  method
    ! iteration info:
    integer  ::  niter, niter_tot(nz), nitermin(nz), nitermax(nz)


#ifdef with_m7
    ! M7 variables:
    real        ::  mass, number, radavmass
    ! M7 cloud processing, should be moved elsewhere?
    real, parameter :: dens_so4   = 1.841    ! g/cm3
    real, parameter :: rad_so4ait = 0.03e-4  ! (cm)
    real, parameter :: rad_so4acc = 0.075e-4 ! (cm)
    real, parameter :: sigma      = 1.59
#endif

    ! to store amount of soa precursors before reaction and that reacts (in part giving rise to soa)
    ! soa precursors: isoprene (iso), terpenes (terp), alkanes (alk), olefines (ole), aromatics (toluene + xylene)
    !real        ::  soa_prec_before (n_soa_prec)
    real        ::  soa_prec_reacted(n_soa_prec)
    !! organic yields
    real        ::  drog(nspec)

    ! meteo data:
    real, pointer        ::   temp(:,:,:)   ! (lon,lat,nz)
    real, pointer        ::     rh(:,:,:)   ! (lon,lat,nz)
    real, pointer        ::   pres(:,:,:)   ! (lon,lat,nz)    
    real, pointer        ::   dens(:,:,:)   ! (lon,lat,nz)    
    real, pointer        :: shumid(:,:,:)   ! (lon,lat,nz)
    real, pointer        ::    occ(:,:,:)   ! (lon,lat,nz)

    ! --- begin ---------------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 't', temp, status, check_units='K' )    
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_GetPointer( 'rh', rh, status, check_units='%' )    
    IF_NOTOK_RETURN(status=1)   
    call LE_Data_GetPointer( 'p', pres, status, check_units='Pa' )        
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3' )     
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'q', shumid, status, check_units ='kg/kg' )            
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'occ', occ, status, check_units ='1')     
    IF_NOTOK_RETURN(status=1)

    ! emissions not inputted here so put to zero
    Q = 0.0

    !! info:
    !if ( .not. outF%suppress ) then
    !  write (gol,'(a,"<chemistry>")') ident2; call goPr
    !end if

#ifdef with_m7
    ! safety check ...
    if ( any((/i_so4ns,i_so4ks,i_so4as,i_so4cs/) > 0) .and. (i_so4a_f > 0) ) then
      write (gol,'("M7 in combination with SO4a_f tracer should not occure ...")'); call goErr
      TRACEBACK; status=1; return
    end if
#endif

    ! flag for tracers in ppb:
    allocate( ppbmask(nspec) )
    ! set:
    do ispec = 1, nspec
      ppbmask(ispec) = trim(specunit(ispec)) == 'ppb'
    end do

    ! bookkeeping of number of iterations
    niter_tot = 0

    ! counters for number of maximum / minimum iterations reached
    nitermin = 0
    nitermax = 0

    !! trick to get terrain elevation
    !maxvalue = maxval(h(:,:,nz))
    !do l = 1, nz
    !  if ( l == 1 ) then
    !    ! full level elevation  =  orography                     + full level height
    !    heff_3d(:,:,l) = max( 0.0, maxvalue - h(:,:,nz) )*1000.0 + 0.5 * h(:,:,l)*1000.0
    !  else
    !    ! full level elevation  =  orography                     + full level height
    !    heff_3d(:,:,l) = max( 0.0, maxvalue - h(:,:,nz) )*1000.0 + 0.5 * ( h(:,:,l-1)*1000.0 + h(:,:,l)*1000.0 )
    !  end if
    !end do

    ! Compiler internal error if the following OpenMP directive is added:
    !OMP     reduction(+:niter_tot,nitermin,nitermax)
    ! Adhoc: set to actual value:
    niter_tot = 3 * nx*ny

    ! init parallel status:
    status_par = 0

    ! start timing:
    call GO_Timer_Start( itim_chem_gas, status )
    IF_NOTOK_RETURN(status=1)
    do k = 1, nz
      !$OMP PARALLEL &
#ifndef __GFORTRAN__
      !$OMP   default ( none ) &
      !$OMP   shared ( cldadj ) &
      !$OMP   shared ( specmolm ) &
#endif
      !$OMP   shared ( k ) &
      !$OMP   shared ( nx, ny ) &
      !$OMP   shared ( c,aerh2o ) &
      !$OMP   shared ( Q ) &
      !$OMP   shared ( temp, occ, shumid, zenith, Rh, dens, pres ) &
      !$OMP   shared ( i_SO4a_f, i_SO4a_c, i_Na_ff, i_Na_f, i_Na_c, i_Na_cc, i_Na_ccc ) &
      !$OMP   shared ( i_no3a_f, i_no3a_c ) &
      !$OMP   shared ( i_so2, i_oh ) &
#ifdef with_m7      
      !$OMP   shared ( i_h2so4 ) &
      !$OMP   shared ( i_nais ) &
      !$OMP   shared ( i_so4ns, i_so4ks, i_so4as, i_so4cs ) &
      !$OMP   shared ( ireac_M7_R71, ireac_M7_R72 ) &
      !$OMP   private( mass, radavmass, number ) &
#endif
      !$OMP   shared ( LUT_Rk_Het ) &
      !$OMP   shared ( deltat ) &
      !$OMP   shared ( ppbmask ) &
      !$OMP   private( i, j ) &
      !$OMP   private( rk, rk_het ) &
      !$OMP   private( aux, icld ) &
      !$OMP   private( rho_air, conc_aer, finc ) &
#ifdef with_vbs      
      !$OMP   private( with_vbs_soa_chemistry ) &
#endif    
      !$OMP   private( soa_prec_reacted ) &
      !$OMP   private( drog ) &
      !$OMP   private( ts, dtmin, dtmax ) &
      !$OMP   private( ch ) &
      !$OMP   private( atol, rtol, method ) &
      !$OMP   private( nfcn, naccpt, nrejec, nstart, dtstart ) &
      !$OMP   private( niter ) &
      !$OMP   private( status ) &
      !$OMP   reduction( + : status_par )
      !$OMP   DO SCHEDULE( DYNAMIC )
      do j = 1, ny
        do i = 1, nx

          ! fill auxilary concentrations:
          aux(iaux_air) = M

          ! fill auxilary meteo variables:
          aux(iaux_T)       = temp(i,j,k)
          ! for backwards compatability, use total cloud cover field;
          ! to be changed to 3D value after improving cloud computations:
          aux(iaux_cloud)   = occ(i,j,1)   ! total value at surface!
          aux(iaux_cldsulf) = occ(i,j,1)   ! total value at surface!
          ! unit conversion to ppm
          !               (kg H2O)/(kg air)   (kg air)/(mole air) / (kg H2O)/(mole H2O)
          aux(iaux_H2O)     = shumid(i,j,k) *       xm_air        / xm_h2o              * 1e6 ! (mole H2O)/(mole air) ppb
          aux(iaux_zen)     = zenith(i,j)
          aux(iaux_Rh)      = Rh  (i,j,k) ! %
          aux(iaux_p)       = pres(i,j,k) ! Pa

          ! index in cloud factor array:
          icld = int( occ(i,j,k)*10.0 + 1.5)
          if ( icld > 11 ) icld = 11
          if  (icld <  1 ) icld = 1
          aux(iaux_cldfac) = cldadj(icld-1)

          ! some meteo parameters:
          rho_air = dens(i,j,k)

          ! air number density (mlc/cm3) :
          !   n/V      Avog     1e-6   =  p / R / T * Avog * 1e-6
          ! mole/m3  mlc/mole  m3/cm3  =       mlc/cm3
          aux(iaux_air) = aux(iaux_p) / Rgas / aux(iaux_T) * Avog * 1e-6  ! (mlc/cm3)

          ! conversion factor from volume mixing ratios (ppm,ppb)
          ! to number density (mlc/cm3) :
          aux(iaux_ppm_to_mlccm3) = 1e-6 * aux(iaux_air)
          aux(iaux_ppb_to_mlccm3) = 1e-9 * aux(iaux_air)

          aux(iaux_H2O)   = aux(iaux_ppm_to_mlccm3) * aux(iaux_H2O)

          ! * hetrogenious rates

          ! ~~ (NH4)SO4a  (fine only )

          ! enabled ?
          if ( i_SO4a_f > 0 ) then
            ! fine nh4hso4 concentrations:
            conc_aer = max( c(i,j,k,i_SO4a_f) * xm_NH4HSO4/xm_SO4, 1.e-6 )
            ! relative increase of aerosol radius due to water uptake ;
            ! use 'aerH2O' array as used by SIA routine (isorropia etc),
            ! this is the water around SO4 aerosols:
            finc = increment_factor_aerh2o( conc_aer, rho_NH4HSO4a, aerh2o(i,j,k), h2odens )
            ! extract rate:
            call LUT_Rk_Het_Get( LUT_Rk_Het(ireac_N2O5_NH4HSO4a_f), &
                                   finc, real(aux(iaux_p)), rk_het(ireac_N2O5_NH4HSO4a_f) )
            ! convert the heterogeneous reaction rate to the proper units:
            !!         1/ppb/min                       1/ppb/s               sec/min
            !rk_het(ireac_N2O5_NH4HSO4a_f) = rk_het(ireac_N2O5_NH4HSO4a_f) * min2sec
            !        1 / mlc/cm3 /s =                1/ppb/s              / (mlc/cm3) /ppb
            rk_het(ireac_N2O5_NH4HSO4a_f) = rk_het(ireac_N2O5_NH4HSO4a_f) / aux(iaux_ppb_to_mlccm3)
          else
            ! dummy:
            rk_het(ireac_N2O5_NH4HSO4a_f) = 0.0
          end if

          ! ~~ sea-salt (fine mode)

          ! enabled ?
          if ( i_Na_f > 0 ) then
          
            ! In look up table, we do not make a difference between the different
            ! classes of fine aerosols. In le_chem_hetro.F90, increment factors are
            ! calculated with one radius for fine aerosols.
            
            ! fine sea salt concentrations:
            conc_aer = max( c(i,j,k,i_Na_f) * Na_to_seasalt, 1.e-6 )
            ! relative increase of aerosol radius due to water uptake ;
            ! based on relative humidity around the seasalt particles:
            finc = increment_factor_rh( conc_aer, rho_seasalt, Rh(i,j,k) )
            !
            ! extract rate:
            call LUT_Rk_Het_Get( LUT_Rk_Het(ireac_N2O5_ss_f), &
                                   finc, real(aux(iaux_p)), rk_het(ireac_N2O5_ss_f) )
            ! convert the heterogeneous reaction rate to the proper units:
            !!         1/ppb/min               1/ppb/s           sec/min
            !rk_het(ireac_N2O5_ss_f) = rk_het(ireac_N2O5_ss_f) * min2sec
            !        1 / mlc/cm3 /s =    1/ppb/s              / (mlc/cm3) /ppb
            rk_het(ireac_N2O5_ss_f) = rk_het(ireac_N2O5_ss_f) / aux(iaux_ppb_to_mlccm3)
            !
            ! extract rate:
            call LUT_Rk_Het_Get( LUT_Rk_Het(ireac_HNO3_ss_f), &
                                   finc, real(aux(iaux_p)), rk_het(ireac_HNO3_ss_f) )
            ! convert the heterogeneous reaction rate to the proper units:
            !!         1/ppb/min               1/ppb/s           sec/min
            !rk_het(ireac_HNO3_ss_f) = rk_het(ireac_HNO3_ss_f) * min2sec
            !        1 / mlc/cm3 /s =    1/ppb/s              / (mlc/cm3) /ppb
            rk_het(ireac_HNO3_ss_f) = rk_het(ireac_HNO3_ss_f) / aux(iaux_ppb_to_mlccm3)
          else
            ! dummy:
            rk_het(ireac_N2O5_ss_f) = 0.0
            rk_het(ireac_HNO3_ss_f) = 0.0
          end if

          ! ~~ sea-salt (coarse mode)

          ! enabled ?
          if ( i_Na_c > 0 ) then

            ! In look up table, we do not make a difference between the different
            ! classes of coarse aerosols. In le_chem_hetro.F90, increment factors are
            ! calculated with one radius for coarse aerosols.
            
            ! coarse sea salt concentrations:
            conc_aer = max( c(i,j,k,i_Na_c) * Na_to_seasalt, 1.e-6 )
            ! relative increase of aerosol radius due to water uptake ;
            ! based on relative humidity around the seasalt particles:
            finc = increment_factor_rh( conc_aer, rho_seasalt, rh(i,j,k) )
            !
            ! extract rate:
            call LUT_Rk_Het_Get( LUT_Rk_Het(ireac_N2O5_ss_c), &
                                   finc, real(aux(iaux_p)), rk_het(ireac_N2O5_ss_c) )
            ! convert the heterogeneous reaction rate to the proper units:
            !!         1/ppb/min               1/ppb/s           sec/min
            !rk_het(ireac_N2O5_ss_c) = rk_het(ireac_N2O5_ss_c) * min2sec
            !        1 / mlc/cm3 /s =    1/ppb/s              / (mlc/cm3) /ppb
            rk_het(ireac_N2O5_ss_c) = rk_het(ireac_N2O5_ss_c) / aux(iaux_ppb_to_mlccm3)
            !
            ! extract rate:
            call LUT_Rk_Het_Get( LUT_Rk_Het(ireac_HNO3_ss_c), &
                                   finc, real(aux(iaux_p)), rk_het(ireac_HNO3_ss_c) )
            ! convert the heterogeneous reaction rate to the proper units:
            !!         1/ppb/min               1/ppb/s           sec/min
            !rk_het(ireac_HNO3_ss_c) = rk_het(ireac_HNO3_ss_c) * min2sec
            !        1 / mlc/cm3 /s =    1/ppb/s              / (mlc/cm3) /ppb
            rk_het(ireac_HNO3_ss_c) = rk_het(ireac_HNO3_ss_c) / aux(iaux_ppb_to_mlccm3)
          else
            ! dummy:
            rk_het(ireac_N2O5_ss_c) = 0.0
            rk_het(ireac_HNO3_ss_c) = 0.0
          end if

          ! * all reaction rates

          ! compute reaction rates, fill some aux values:
          call LE_Chem_Work_Rates( rk, aux, rk_het )

          ! * solve

          ! other aux values:
          aux(iaux_nox0) = nox0

          ! Set conc's and emissions in help array:
          ch = max( 0.0, c(i,j,k,:) )  ! ppb, ug/m3, mlc/cm3, #/cm3

          ! convert units of gas phase tracers from ppb to mlc/cm3
          where ( ppbmask )
            ch = ch * aux(iaux_ppb_to_mlccm3)
          endwhere
          ! convert the units of selected aerosols from ug/m3 to mlc/cm3 :
          !   (ug tr)/(m3 air) * (kg tr)/(ug tr) / ((kg tr)/(mole tr)) * mlc/mole * m3/cm3
          !         c          *      1e-9       /        xm           *   Avog   *  1e-6
          if ( i_so4a_f > 0 ) ch(i_so4a_f) = ch(i_so4a_f) / specmolm(i_so4a_f) * Avog * 1.0e-15  ! mlc/cm3
          if ( i_no3a_f > 0 ) ch(i_no3a_f) = ch(i_no3a_f) / specmolm(i_no3a_f) * Avog * 1.0e-15  ! mlc/cm3
          if ( i_no3a_c > 0 ) ch(i_no3a_c) = ch(i_no3a_c) / specmolm(i_no3a_c) * Avog * 1.0e-15  ! mlc/cm3
          if ( i_na_ff  > 0 ) ch(i_na_ff ) = ch(i_na_ff ) / specmolm(i_na_ff ) * Avog * 1.0e-15  ! mlc/cm3
          if ( i_na_f   > 0 ) ch(i_na_f  ) = ch(i_na_f  ) / specmolm(i_na_f  ) * Avog * 1.0e-15  ! mlc/cm3
          if ( i_na_c   > 0 ) ch(i_na_c  ) = ch(i_na_c  ) / specmolm(i_na_c  ) * Avog * 1.0e-15  ! mlc/cm3
          if ( i_na_cc  > 0 ) ch(i_na_cc ) = ch(i_na_cc ) / specmolm(i_na_cc ) * Avog * 1.0e-15  ! mlc/cm3
          if ( i_na_ccc > 0 ) ch(i_na_ccc) = ch(i_na_ccc) / specmolm(i_na_ccc) * Avog * 1.0e-15  ! mlc/cm3

#ifdef with_m7
          !>>> Combination of 'so4a_f' and M7 tracers should not occure anymore,
          !    this is checked at the top of this routine ...
          !! contribute all sulfate from m7 particles to bulk so4 mass for chemistry;
          !! use that the M7 sulphate tracers are already in mlc/cm3 ...
          !!                     mlc/cm3                      mlc/cm3
          !!should be replaced by decent SIA chemistry in M7
          !if (i_so4a_f > 0) then
          !   if ( i_so4ns > 0 ) ch(i_so4a_f) = ch(i_so4a_f) + ch(i_so4ns)  ! mlc/cm3
          !   if ( i_so4ks > 0 ) ch(i_so4a_f) = ch(i_so4a_f) + ch(i_so4ks)  ! mlc/cm3
          !   if ( i_so4as > 0 ) ch(i_so4a_f) = ch(i_so4a_f) + ch(i_so4as)  ! mlc/cm3
          !   if ( i_so4cs > 0 ) ch(i_so4a_c) = ch(i_so4a_c) + ch(i_so4cs)  ! mlc/cm3
          !end if
          !<<<
#endif          
          
          ! VBS code left here for sed for inspiration.
          ! First solve bug in computation of 'drog', see comment in top of 'LE_VBS'.
          ! VBS enabled ?
          !if ( n_vbs > 0 ) then
          !  ! soa chemistry enabled ?
          !  if ( with_vbs_soa_chemistry ) then
          !    ! store concentrations before chemistry step (not for alkanes = PAR)
          !    soa_prec_before(i_soa_prec_par ) = 0.0        ! alkanes
          !    soa_prec_before(i_soa_prec_ole ) = ch(i_ole ) ! olefines
          !    soa_prec_before(i_soa_prec_tol ) = ch(i_tol ) ! toluene
          !    soa_prec_before(i_soa_prec_xyl ) = ch(i_xyl ) ! xylene
          !    soa_prec_before(i_soa_prec_iso ) = ch(i_iso ) ! isoprene
          !    soa_prec_before(i_soa_prec_terp) = ch(i_terp) ! terpene
          !  end if
          !end if

          ! solver input:
          ts     = min2sec * 0.0
          dtmin  = min2sec * 1.0
          dtmax  = min2sec * deltat/4.0
          atol   = aux(iaux_ppb_to_mlccm3)*1e-3
          rtol   = 1e-2
          method = 0  ! not used in twostep
          niter  = 3  ! fixed value, not determined online

          ! check ...
          if ( min2sec * deltat < dtmin ) then
            write (gol,'("chemistry timestep over ",f7.2," sec while minium step is ",f7.2," sec")') &
                            min2sec * deltat, dtmin; call goErr
            !TRACEBACK; status=1; return
            status_par = status_par + 1
          
          else

            ! solve:
            call TWOSTEP( nspec, ts, min2sec * deltat, dtmin, dtmax, &
                          ch, Q, &
#ifdef with_vbs
                          drog, &
#endif                        
                          nreac, rk, naux, aux, &
                          atol, rtol, method, niter, &
                          nfcn, naccpt, nrejec, nstart, dtstart, &
                          LE_Chem_Work_Iter, status )
            ! add status:
            status_par = status_par + status
            ! add to number of iterations
            !niter_tot(l) = niter_tot(l) + niter       ! <-- compiler internal error on reduction statement

          end if   ! check timestep

#ifdef with_m7
          ! M7 aerosol scheme, needs some improvement ...
          ! - replace this by reaction in "reactions.csv" file ?
          ! - there is three times the same factor "deltat * rk(ireac_M7_R72) * ch(i_so2)"
          ! - sulphate masses are now reset, should they be updated ?
          if ( n_m7 > 0 ) then
            ! production of H2SO4, in mlc/cm3 ;
            ! reaction rate is taken from:
            !  R71  :  SO2 + OH -> SO4a_f + HO2
            ! Originally direct " nucleation"  to so4a_f assumed
            ch(i_h2so4) = ch(i_h2so4) + min2sec*deltat * rk(ireac_M7_R71) * ch(i_so2) * ch(i_oh)
            ! rk(72) represents multi-phase reactions that would yield
            ! soluble accumulation mode sulphate
            ! now neglected as source for sulphuric acid
            ! should be incorporated somewhere else into the right bin: XXX ACTION
            ch(i_so4ks) = ch(i_so4ks) + min2sec*deltat * rk(ireac_M7_R72) * ch(i_so2)                ! (molec/cm3)
            mass        =               min2sec*deltat * rk(ireac_M7_R72) * ch(i_so2) * 96.0/6.0e23  ! (g/cm3)
            ! update number concentration;
            ! these should be implemented with functions from 'jaql_particles' ...
            radavmass = rad_so4ait * exp( 1.5*(log(sigma))**2.0 ) !(cm)
            number=mass*3./(4.*3.14*(radavmass**3.)*dens_so4) !#/cm3
            ch(i_nais)=  ch(i_nais)+number !#/cm3
            ! reset bulk mass:
            if ( i_so4a_f > 0 ) ch(i_so4a_f) = ch(i_so4ns) + ch(i_so4ks) + ch(i_so4as) ! mlc/cm3
            if ( i_so4a_c > 0 ) ch(i_so4a_f) = ch(i_so4cs)  ! mlc/cm3
          end if
#endif
          ! convert units of gas phase tracers from mlc/cm3 back to ppb:
          where ( ppbmask )
            ch = ch / aux(iaux_ppb_to_mlccm3)
          endwhere
          ! convert aerosol concentrations from mlc/cm3 back to ug/m3 :
          if ( i_so4a_f > 0 ) ch(i_so4a_f) = ch(i_so4a_f) * specmolm(i_so4a_f) / Avog * 1.0e15  ! ug/m3
          if ( i_no3a_f > 0 ) ch(i_no3a_f) = ch(i_no3a_f) * specmolm(i_no3a_f) / Avog * 1.0e15  ! ug/m3
          if ( i_no3a_c > 0 ) ch(i_no3a_c) = ch(i_no3a_c) * specmolm(i_no3a_c) / Avog * 1.0e15  ! ug/m3
          if ( i_na_f   > 0 ) ch(i_na_f  ) = ch(i_na_f  ) * specmolm(i_na_f  ) / Avog * 1.0e15  ! ug/m3
          if ( i_na_ff  > 0 ) ch(i_na_ff ) = ch(i_na_ff ) * specmolm(i_na_ff ) / Avog * 1.0e15  ! ug/m3
          if ( i_na_c   > 0 ) ch(i_na_c  ) = ch(i_na_c  ) * specmolm(i_na_c  ) / Avog * 1.0e15  ! ug/m3
          if ( i_na_cc  > 0 ) ch(i_na_cc ) = ch(i_na_cc ) * specmolm(i_na_cc ) / Avog * 1.0e15  ! ug/m3
          if ( i_na_ccc > 0 ) ch(i_na_ccc) = ch(i_na_ccc) * specmolm(i_na_ccc) / Avog * 1.0e15  ! ug/m3

#ifdef with_vbs
          ! VBS enabled ?
          if ( n_vbs > 0 ) then

            ! In the "v1.7/proj/VBS/v1.7.005" version now a call was made to the TwoStep
            ! solver using a special iteration routine named "iter_vbs".
            ! The only reactions implemented there were the reactions R77-R79
            ! describing the loss of TERP due to presence of OH, O3, and NO3 respectively.
            ! Note that the result of these reactions on the concentrations of OH/O3/NO3
            ! were omitted, since these losses are assumed to be implicitly defined
            ! already in the reaction rates.
            ! The same result is now achieved with the automatically genereated code
            ! based on the "reactions.csv" by describing R77-R79 as catalyst reactions.

            ! soa chemistry enabled ?
            if ( with_vbs_soa_chemistry ) then
              ! after previous call of twostep (with function iter_cbm4 as argument)
              ! drog(1) contains the amount of PAR reacted per minute.
              ! Store this in soa_prec_reacted (and correct for time step deltat instead of minute)
              soa_prec_reacted(i_soa_prec_par)  = drog( ispecs_soa_prec(i_soa_prec_par) )  ! paraffines
              soa_prec_reacted(i_soa_prec_ole ) = drog( ispecs_soa_prec(i_soa_prec_ole) )  ! olefines
              soa_prec_reacted(i_soa_prec_tol ) = drog( ispecs_soa_prec(i_soa_prec_tol) )  ! toluene
              soa_prec_reacted(i_soa_prec_xyl ) = drog( ispecs_soa_prec(i_soa_prec_xyl) )  ! xylene
              soa_prec_reacted(i_soa_prec_iso ) = drog( ispecs_soa_prec(i_soa_prec_iso) )  ! isoprene
              soa_prec_reacted(i_soa_prec_terp) = drog( ispecs_soa_prec(i_soa_prec_terp) ) ! terpene
              
              ! convert back from mlc/cm3 to ppb
              soa_prec_reacted = soa_prec_reacted / aux(iaux_ppb_to_mlccm3)
            else
             ! turn off production of condensables from PAR, OLE, TOL, XYL, ISO and TERP chemistry
             soa_prec_reacted = 0.0
            end if

            ! soa_prec_reacted is in ppb (per timestep deltat) ;
            ! routines seem to assume that ch is in ppb for gasses and ug/m3 for aerosol:
            call LE_VBS_Apply( ch, soa_prec_reacted, nreac, rk, &
                                real(aux(iaux_p)), real(aux(iaux_T)), deltat )

          endif
#endif


          ! restore concentrations:
          c(i,j,k,:) = ch(:)

          ! *

        end do ! i
      end do ! j
      !$OMP   END DO
      !$OMP END PARALLEL

      ! check status ...
      if ( status_par /= 0 ) then
        write (gol,'("error status returned from OpenMP loop : ",i6)') status_par; call goErr
        TRACEBACK; status=1; return
      end if

    end do ! k

    ! stop timing:
    call GO_Timer_End( itim_chem_gas, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( ppbmask )

    !! integration statistics
    !! divide niter_tot by the number of grid cells
    !niter_tot = niter_tot/(nx*ny)
    !nitermin  = 100*nitermin/(nx*ny)
    !nitermax  = 100*nitermax/(nx*ny)

    !!print message
    !if (.NOT.outF%suppress) then
    !  print *,ident2,'  <updated gas phase concentrations>'
    !  print *,ident2,'    <average #iterations per layer>', niter_tot
    !endif
    !write (u_log_c,*) ' # iterations per layer: ',niter_tot
    !write (u_log_c,*) ' percentage miniter:     ',nitermin
    !write (u_log_c,*) ' percentage maxiter:     ',nitermax

    ! cloud chemistry enabled ?
    if ( with_chem_cloud ) then
      ! start timing:
      call GO_Timer_Start( itim_chem_cloud, status )
      IF_NOTOK_RETURN(status=1)
      ! perform cloud chemistry ; timestep should be in seconds:
      call wetS( 60*deltat, c, status )
      IF_NOTOK_RETURN(status=1)

      ! end timing:
      call GO_Timer_End( itim_chem_cloud, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! aerosols >>>
    call GO_Timer_Start( itim_chem_aero, status )
    IF_NOTOK_RETURN(status=1)

    ! solve aerosol equilibrium chemistry
    ! skip aerosol computations if first time step; model needs to be filled ....
    !if (runF%do_aerosol .AND. (.NOT.runF%first)) then
    if ( n_sia > 0 ) then
      !<--
      ! prevent aerosol concentrations  to be too small
      where (c(:,:,:,i_so4a_f) < 1.0e-4) c(:,:,:,i_so4a_f) = 1.0e-4
      where (c(:,:,:,i_nh4a_f) < 1.0e-4) c(:,:,:,i_nh4a_f) = 1.0e-4
      where (c(:,:,:,i_no3a_f) < 1.0e-4) c(:,:,:,i_no3a_f) = 1.0e-4
      !where (c(:,:,:,i_no3a_c) < 1.0e-4) c(:,:,:,i_no3a_c) = 1.0e-4
      !<-- not yet, testing backwards compatibility
      !! prevent aerosol concentrations  to be too small:
      !do itr = 1, n_sia
      !  ispec = ispec_sia(itr)
      !  where (c(:,:,:,ispec) < 1.0e-4) c(:,:,:,ispec) = 1.0e-4
      !end do
      !<--
      !if (.NOT.outF%suppress) print *,ident2,'  <performing aerosol equilibrium calculations>'
      !! info ...
      !write (gol,'("         chemistry - calling ",a," ...")') trim(aero_code); call goPr
      ! different codes ...
      select case ( trim(aero_code) )
        case ( 'isorropia2' )
          call LE_SIA_Isorropia2_Step( c, aerh2o, status )
          IF_NOTOK_RETURN(status=1)
        case ( 'equisam' )
          call LE_SIA_Equisam_Step( c, aerh2o, status )
          IF_NOTOK_RETURN(status=1)
        case default
          write (gol,'("unsupported aerosol code : ",a)') trim(aero_code); call goErr
          TRACEBACK; status=1; return
      end select
      !if (.NOT.outF%suppress) print *,ident2,'  <finished aerosol equilibrium calculations>'
    endif

    call GO_Timer_End( itim_chem_aero, status )
    IF_NOTOK_RETURN(status=1)
    ! <<< aerosols

    ! ok
    status = 0

  end subroutine LE_Chem_Step


  ! ***


  subroutine calc_r1r3( cloud,coszen, tk, status )

    use dims , only : zenith
    !use dims , only : h
    use dims, only : nx, ny
    use dims, only : rk1rk3
    use LE_Chem_Work, only : cldadj

    real, intent(in)          ::  cloud(nx, ny)
    real, intent(in)          ::  coszen(nx,ny)
    real, intent(in)          ::  tk(nx,ny)
    integer, intent(out)      ::  status

    character(len=*), parameter ::  rname = mname//'/calc_r1r3'

    real :: factor(nx, ny)
    integer :: index
    real    :: heff
    integer :: i, j

    ! correction factors for JNO2
    real :: corr_phot
    ! for the whole grid

    ! computes the ratio r1/r3 between two reactions
    !  r1 = NO2 --> O3 + NO
    !  r3 = O3 + NO --> NO2
    ! compute factor for cloud cover
    ! compute cloud cover correction factor * 60

    do i=1,nx
       do j=1,ny

          ! Compute height (middle of lower layer)
          !h1   = h(i,j,1)*1000.0
          !elev = max(0.0, hmax - h(i,j,nz) ) *1000 ! trick to get terrain elevation
          !heff = elev + 0.5*h1
          heff = 0

          call get_JNO2_corr( trim(JNO2_corr_file), heff, zenith(i,j), corr_phot, status )
          IF_NOTOK_RETURN(status=1)
          index = int(cloud(i,j)*10+1.5) - 1
          index = max(0, min(index, 10))
          factor(i,j)=60.0*cldadj(index)*corr_phot
       enddo
    enddo

    ! choose the cut-off such that the argument of
    ! the exponent is > -10.0 to prevent underflow
    where (coszen > 0.04)
      rk1rk3 = 1.45e-2*exp(-0.4/coszen)*factor/( 2.952*exp(-1450.0/tk) )
    elsewhere
      rk1rk3 = 0.0
    endwhere

    ! ok
    status = 0

  end subroutine calc_r1r3


  ! ***


  subroutine get_JNO2_corr( fname, heff, zen, corr_phot, status )

  use GO, only : pathsep
  use GO, only : goGetFU
  use dims, only : runF
  use LE_Logging, only : ident2

  character(len=*), intent(in)    ::  fname
  real, intent(in)                ::  heff
  real, intent(in)                ::  zen
  real, intent(out)               ::  corr_phot
  integer, intent(out)            ::  status

  character(len=*), parameter ::  rname = mname//'/get_JNO2_corr'

  integer             ::  u_tmp
  integer :: i, k, iz, izen
  real    :: fac, fac_zen, zen2

  !add by Ger, aug 2004
  if (angle_JNO2_corr(10) == 0.0) then
    !print *,ident2,'  <reading JNO2 correction factors>'
    ! new file unit:
    call goGetFU( u_tmp, status )
    IF_NOTOK_RETURN(status=1)
    open (u_tmp, file=trim(fname), status='old')
    read (u_tmp,*) angle_JNO2_corr(:)
    do i=1,10
      read (u_tmp,*) h_JNO2_corr(i), corr_jno2(i,:)
    enddo
    close (u_tmp)
  endif

  zen2 = zen*57.29577951
  ! effective height for which we need the photolyis rate correction
  if (zen2 >= 0.0 .AND. zen2 <= 90.0) then
     iz = -1
     do k=1,9
       if (heff >= h_JNO2_corr(k) .AND. heff <= h_JNO2_corr(k+1) ) then
           iz = k
           exit
       endif
     enddo
     !write(*,*)'k, h_JNO2_corr ', k, h_JNO2_corr(k), heff, h_JNO2_corr(k+1),iz

     ! if iz wasn't found: assume we are above the highest level
     if (iz == -1) iz=9
     fac = (heff-h_JNO2_corr(iz))/(h_JNO2_corr(iz+1)-h_JNO2_corr(iz))
     !if (fac.gt.0.0) write( *,*)' zen2,fac,iz,heff ', zen2,fac,iz,heff
     do k=1, 11
       !write(*,*)' k,izen, angle-no2 ',k, izen, angle_JNO2_corr(k)
       if (zen2 >= angle_JNO2_corr(k) .AND. zen2 <= angle_JNO2_corr(k+1) ) izen = k
     enddo
     fac_zen   = (zen2-angle_JNO2_corr(izen))/(angle_JNO2_corr(izen+1)-angle_JNO2_corr(izen))
     corr_phot = (1.0-fac)*( (1.0-fac_zen)*corr_jno2(iz  ,izen) + fac_zen*corr_jno2(iz  ,izen+1) ) + &
                      fac *( (1.0-fac_zen)*corr_jno2(iz+1,izen) + fac_zen*corr_jno2(iz+1,izen+1) )
  else
     corr_phot=1.0
  endif
  ! added by Ger, aug 2004
  !if (corr_phot <> 1.0) write(*,*)' izen, fac_zen, cor_phot', izen, fac_zen, corr_phot

    ! ok
    status = 0

  end subroutine get_JNO2_corr


end module LE_Chem
