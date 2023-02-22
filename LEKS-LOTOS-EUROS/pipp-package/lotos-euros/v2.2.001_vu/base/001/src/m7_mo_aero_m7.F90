MODULE mo_aero_m7
  
  ! *mo_aero_m7*  Contains parameters, switches and initialization
  !               routines for the m7 aerosol scheme.
  !
  ! Authors:
  ! --------
  ! E. Vignati (JRC/IES)						   2005
  ! P. Stier (MPI)                                 2001/2002 
  ! E. Vignati and J. Wilson (JRC / EI)            2000

  IMPLICIT NONE

  SAVE

  !--- 1) Define and pre-set switches for the processes of M7: -----------------------

  !--- Physical:

  LOGICAL :: lsnucl     = .TRUE., & ! nucleation
             lscoag     = .TRUE., & ! coagulation
             lscond     = .TRUE.    ! condensation of H2SO4
  
  INTEGER :: nnucl      = 3         ! Choice of the nucleation scheme:
                                    ! 
                                    !    nnucl = 1  Vehkamaeki (2002)
                                    !          = 2  Kulmala (1998)
                                    !          = 3  Activation (Kulmala 2006)

  !--- Technical:

  LOGICAL :: lmass_diag = .FALSE.   ! mass balance check in m7_interface

  !--- 2) Numbers of compounds and modes of m7: --------------------------------------

  INTEGER, PARAMETER :: naermod=18,        & !number of all compounds
                        nmod=7,            & !number of modes
                        nss=2,             & !number of sea salt compounds 
                        nsol=4,            & !number of soluble  compounds 
                        ngas=3,            & !number of gaseous  compounds 
                        nsulf=4              !number of sulfate  compounds 
                        
  !--- 3) List of indexes corresponding to the compound masses and mode numbers:------

  !--- 3.1) Mass index (in array aerml and ttn): 
  !
  !         Attention:
  !         The mass of sulfate compounds is always given in [molec. cm-3] 
  !         whilst the mass of other compounds is given in [ug cm-3].
  !
  !         Compounds:
  !
  !           so4 = sulphate
  !           bc  = black carbon
  !           oc  = organic carbon, 
  !           ss  = sea salt
  !           du  = dust 
  !
  !         Modes:
  !
  !           n   = nucleation mode
  !           k   = aitken mode 
  !           a   = accumulation mode
  !           c   = coarse mode
  !
  !         Type:
  !
  !           s   = soluble mode
  !           i   = insoluble mode
                                                                                  !  COMPOUND:
  INTEGER, PARAMETER ::                                                         &
           iso4ns=1, iso4ks=2, iso4as=3, iso4cs=4,                              & !- Sulfate
                     ibcks =5, ibcas =6, ibccs =7, ibcki =8,                    & !- Black Carbon
                     iocks =9, iocas=10, ioccs=11, iocki=12,                    & !- Organic Carbon
                               issas=13, isscs=14,                              & !- Sea Salt
                               iduas=15, iducs=16,           iduai=17, iduci=18   !- Dust  
  ! MODE:           |         |         |         |         |
  !         nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !         soluble | soluble | soluble | soluble | insol.  | insol.  | insol.  |

   
  !--- 3.2) Number index (in array aernl):
  !

  INTEGER, PARAMETER ::                                      &
           inucs=1,  iaits=2,  iaccs=3,  icoas=4,  iaiti=5,  iacci=6,  icoai=7    
  ! MODE:           |         |         |         |         |
  !         nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !         soluble | soluble | soluble | soluble | insol.  | insol.  | insol.  |


  !--- 4) Definition of the modes of M7: ------------------------------------------------------

  !--- 4.1) Threshold radii between the different modes [cm]:
  !         Used for the repartititioning in m7_dconc.
  !         crdiv(jmod) is the lower bound and crdiv(jmod+1) is 
  !         the upper bound of the respective geometric mode
  !         Default value for nucleation mode is modified by the 
  !         choice of the nuclation scheme.

  REAL :: crdiv(4)=(/ 0.0005E-4, 0.005E-4, 0.05E-4, 0.5E-4 /)    
  !                             |         |        |      
  !                             |         |        |
  !                 nucleation -- aitken  -  accum -- coarse mode

  !--- 4.2) Standard deviation for the modes:

  REAL, PARAMETER :: sigma(nmod)=(/ 1.59, 1.59, 1.59, 2.00, 1.59, 1.59, 2.00 /)

  !--- Natural logarithm of the standard deviation of each mode:
  !    Calulated in m7_initialize. 

  REAL            :: sigmaln(nmod)

  !--- 5) Conversion factors for lognormal particle size distributions: -------------
  !       Calulated in m7_initialize. 

  REAL            :: cmr2ras(nmod) ! Conversion factor: count median radius to radius of average surface

  REAL            :: cmr2mmr(nmod) ! Conversion factor: count median radius to mass mean radius

  REAL            :: cmedr2mmedr(nmod) ! Conversion factor: count median radius to mass median radius

  REAL            :: cmr2ram(nmod) ! Conversion factor: count median radius to radius of average mass

  REAL            :: ram2cmr(nmod) ! Conversion factor: radius of average mass to count median radius


  !--- 6) Assumed thresholds for occurence of specific quantities: -------------
  !@@@    To be done!
  
  !  REAL, PARAMETER :: cmin_aerml     = 1.E-15 , ! Aerosol mass
  !                     cmin_aernl     = 1.E-10 , ! Aerosol number
  !                     
  
  !--- 7) Chemical constants: ----------------------------------------------------
  !
  !--- Accomodation coefficient of H2SO4 on aerosols:
  !    (reduced for insoluble modes)

  REAL, PARAMETER :: caccso4(nmod) = (/ 1.0, 1.0, 1.0, 1.0, 0.3, 0.3, 0.3 /)

  !--- Critical relative humidity:

  REAL, PARAMETER :: crh    = 0.45           ! Assumed relative humidity for the 
                                             ! Na2SO4 / NaCl system below which 
                                             ! crystalization occurs.
                                             ! (estimated from Tang, I.N.; JGR 102, D2 1883-1893)

  !--- 8) Physical constants: ----------------------------------------------------
  !
  !--- 8.1) General physical constants: 

  REAL, PARAMETER :: bk      = 1.38e-16,   & ! Bolzman constant []
                     avo     = 6.02217E+23,& ! Avogadro number [mol-1]
                     rerg    = 8.314E+7,   & ! Ideal gas constant [erg.K-1.mole-1]
                     r_kcal  = 1.986E-3      ! Ideal gas constant [kcal K-1.mole-1]
  
  !--- 8.2) Type specific physical constants:
  !
  REAL, PARAMETER :: dh2so4  = 1.841,      & ! Density          H2SO4  [g cm-3]
                     ddust   = 2.650,      & ! Density          du     [g cm-3]
                     dbc     = 2.,         & ! Density          bc     [g cm-3]
                     doc     = 2.,         & ! Density          oc     [g cm-3]
                     dnacl   = 2.165,      & ! Density          NaCl   [g cm-3]
                     dna2so4 = 2.68,       & ! Density          Na2SO4 [g cm-3]
                     dnahso4 = 2.435,      & ! Density          NaHSO4 [g cm-3]
                     dh2o    = 1.0,        & ! Density          H2O    [g cm-3]

                     wh2so4  = 98.0734,    & ! Molecular weight H2SO4  [g mol-1]
                     wh2o    = 18.0,       & ! Molecular weight H2O    [g mol-1]
                     wso4    = 96.0576,    & ! Molecular weight SO4    [g mol-1]
                     wso2    = 64.0,       & ! Molecular weight SO2    [g mol-1]
                     wna     = 22.99,      & ! Atomic    weight Na     [g mol-1]
                     wcl     = 35.453,     & ! Atomic    weight Cl     [g mol-1]
                     wnacl   = 58.443,     & ! Molecular weight NaCl   [g mol-1]
                     wna2so4 = 142.0376,   & ! Molecular weight Na2SO4 [g mol-1]
                     wnahso4 = 120.0555      ! Molecular weight NaHSO4 [g mol-1]


  !--- 9) Assumed parameters: ------------------------------------------------------

  REAL, PARAMETER :: critn=100.,           & ! Assumed mass of an nucleated sulfate 
                                             ! particle for the Kulmala scheme [molecules]
                     fmax=0.95,            & ! Factor that limits the condensation 
                                             ! of sulfate to fmax times the available
                                             ! sulfate in the gas phase [1].
                                             ! (m7_dgas)
                     cLayerThickness = 1.0   ! Assumed required layer thickness of
                                             ! sulfate to transfer an insoluble 
                                             ! particle to a soluble mode. It is 
                                             ! given in units of layers of 
                                             ! monomolecular sulfate. Determines the
                                             ! transfer rate from insoluble to 
                                             ! soluble modes. 

  !--- 10) Computational constants: ------------------------------------------------

  REAL, PARAMETER :: sqrt2=1.4142136,    & 
                     pi=3.141592654

  !--- 11) Data used for the calculation of the aerosol properties -----------------
  !       under ambient conditions:
  !       (Included the conversion from Pa to hPa in the first parameter.)

  REAL, PARAMETER :: wvb(17)=                                                   &
                     (/   95.80188,     -28.5257,     -1.082153,     0.1466501, &
                         -20627.51,    0.0461242,     -0.003935,      -3.36115, &
                       -0.00024137,  0.067938345, 0.00000649899,   8616124.373, &
                       1.168155578, -0.021317481,   0.000270358, -1353332314.0, &
                      -0.002403805                                              /)

  REAL, PARAMETER :: gmb(9)=                                                 &
                     (/ 1.036391467, 0.00728531, -0.011013887, -0.068887407, &
                        0.001047842, 0.001049607, 0.000740534, -1.081202685, &
                       -0.0000029113                                         /)

  !--- 4) Logical mask for coagulation kernel: -------------------------------------
  !       (The coagulation kernel mask is symmetric and not all 
  !        values are used for physical considerations. As its 
  !        calculation is very expensive, a mask is used to 
  !        calculate only the necessarey elements.)

  LOGICAL :: locoagmask(nmod,nmod)

  DATA locoagmask(1:nmod,1) / .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.  /

  DATA locoagmask(1:nmod,2) / .FALSE., .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.  /

  DATA locoagmask(1:nmod,3) / .FALSE., .FALSE., .TRUE.,  .FALSE., .TRUE.,  .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,4) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,5) / .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.,  .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,6) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,7) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /


  !--- 12) Service routines for initialization and auxiliary computations ----------

CONTAINS

  SUBROUTINE m7_initialize

    ! Purpose:
    ! ---------
    ! Initializes constants and parameters 
    ! used in the m7 aerosol model.
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                          may 2001
    !
    ! Interface:
    ! ---------
    ! *m7_initialize*  is called from *call_init_submodels* in *call_submodels*
    !

    IMPLICIT NONE

    INTEGER :: jmod
      
    DO jmod=1, nmod

       !--- 1) Calculate conversion factors for lognormal distributions:----
       !       Radius of average mass (ram) to count median radius (cmr) and 
       !       vice versa. Count median radius to radius of average 
       !       mass (ram).
       !       These factors depend on the standard deviation (sigma)
       !       of the lognormal distribution.
       !       (Based on the Hatch-Choate Conversins Equations;
       !        see Hinds, Chapter 4.5, 4.6 for more details.
       !        In particular equation 4.53.)

       !--- Count Median Radius to Mass Median Radius:

       cmedr2mmedr(jmod) = EXP(3.0*(LOG(sigma(jmod)))**2)

       !--- Count Median Radius to Mass Mean Radius:

       cmr2mmr(jmod) = EXP(3.5*(LOG(sigma(jmod)))**2)

       !--- Count Median Radius to Radius of Average Mass:

       cmr2ram(jmod) = EXP(1.5*(LOG(sigma(jmod)))**2)

       !--- Radius of Average Mass to Count Median Radius:

       ram2cmr(jmod) = 1. / cmr2ram(jmod)

       !--- Count Median Radius to Radius of Average Surface:

       cmr2ras(jmod) = EXP(1.0*(LOG(sigma(jmod)))**2)


       !--- 2) Calculate the natural logarithm of the standard deviation:

       sigmaln(jmod) = LOG(sigma(jmod))

    END DO

    !--- 3) Set the lower mode boundary for the nucleation mode:
    !       (Depends on the choice of the nucleation scheme.)

    SELECT CASE (nnucl) 
    CASE(1)
       crdiv(1)=0.0005E-4

    CASE(2)
       crdiv(1)=( critn * wh2so4/avo/dh2so4*0.75/pi )**(1./3.)
    END SELECT

    !--------------------------------------------------------------------

  END SUBROUTINE m7_initialize

  SUBROUTINE m7_coat(pm6rp_lon_lev, pcrtcst)
        
    ! Purpose:
    ! ---------
    ! *m7_coat* calculates the number of sulfate 
    !           molecules required to coat a particle
    !           with cLayerThickness of sulfate
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                          2001
    !
    ! Interface:
    ! ---------
    ! *m7_coat* is called from *m7_concoag*
    !

    IMPLICIT NONE

    INTEGER         :: jmod 

    REAL            :: pm6rp_lon_lev(nmod)    ! Ambient radii for current
                                              ! longitude and level [cm]
    REAL            :: pcrtcst(nmod)          ! Critical constant, i.e. number of
                                              ! sulfate to cover an average particle
                                              ! of the mode with a layer of the 
                                              ! thickness determined by cLayerThickness.
    REAL            :: zras(nmod)             ! Radius of average surface 
                                              ! for a single particle [cm]
    REAL            :: zas(nmod)              ! Average surface 
                                              ! for single particle [cm+2]
        
    REAL, PARAMETER :: csurf_molec = 2.39E-15 ! Average cross-section 
                                              ! of a single H2SO4 molecule [cm+2]

    !--- 1) Calculate the radii of average surface for modes 5-7:

    zras(5) = pm6rp_lon_lev(5) * cmr2ras(5)
    zras(6) = pm6rp_lon_lev(6) * cmr2ras(6)
    zras(7) = pm6rp_lon_lev(7) * cmr2ras(7)

    DO jmod=5, 7
           
       !--- 2) Calculate the average surface of an particle for modes 5-7:

       zas(jmod)    = 4 * zras(jmod)**2 * pi

       !--- 3) Determine the number of sulfate molecules needed to form
       !       a cLayerThickness thick layer of sulfate on the particles
       !       in modes 5-7:
           
       pcrtcst(jmod) = (zas(jmod) / csurf_molec) * cLayerThickness

    END DO

  END SUBROUTINE m7_coat
  
!!$  SUBROUTINE setaeroM7
!!$
!!$    ! *setaeroM7* modifies pre-set switches of the aeroM7ctl
!!$    !             namelist for the configuration of the 
!!$    !             M7 component of the ECHAM/HAM aerosol model
!!$    ! 
!!$    ! Authors:
!!$    ! --------
!!$    ! Philip Stier, MPI-MET                        12/2002
!!$    !
!!$    ! *setaeroM7* is called from *call_init_submodels* in *call_submodels*
!!$    !
!!$
!!$    USE mo_mpi,       ONLY: p_parallel, p_parallel_io, p_bcast, p_io
!!$    USE mo_namelist,  ONLY: position_nml, nnml, POSITIONED
!!$    USE mo_doctor,    ONLY: nout
!!$    USE mo_exception, ONLY: finish
!!$
!!$    IMPLICIT NONE
!!$
!!$    INCLUDE 'aerom7ctl.inc'
!!$
!!$    !--- Local variables:
!!$
!!$    INTEGER :: ierr
!!$
!!$    !--- 1) Read namelist:
!!$
!!$    IF (p_parallel_io) THEN
!!$       CALL position_nml ('AEROM7CTL', status=ierr)
!!$       SELECT CASE (ierr)
!!$       CASE (POSITIONED)
!!$          READ (nnml, aerom7ctl)
!!$       END SELECT
!!$    ENDIF
!!$
!!$    !--- 2) Broadcast over processors:
!!$
!!$    IF (p_parallel) THEN
!!$       CALL p_bcast (lsnucl,     p_io)
!!$       CALL p_bcast (lscoag,     p_io)
!!$       CALL p_bcast (lscond,     p_io)
!!$       CALL p_bcast (nnucl,      p_io)
!!$       CALL p_bcast (lmass_diag, p_io)
!!$    END IF
!!$
!!$    IF (.NOT. p_parallel) THEN
!!$       WRITE(nout,*) ''
!!$       WRITE(nout,*) ''
!!$       WRITE(nout,*) '----------------------------------------------------------'
!!$       WRITE(nout,*) '----------------------------------------------------------'
!!$       WRITE(nout,*) '--- Initialization of settings for aerosol module M7   ---'
!!$       WRITE(nout,*) '---'
!!$       WRITE(nout,*) '---    Default values of aeroctl modified by setaero:'
!!$       WRITE(nout,*) '---'
!!$       WRITE(nout,*) '---    New settings: lsnucl = ', lsnucl
!!$       WRITE(nout,*) '---                  lscoag = ', lscoag
!!$       WRITE(nout,*) '---                  lscond = ', lscond
!!$       IF (nnucl==1) THEN
!!$       WRITE(nout,*) '---                  nnucl  = ', nnucl
!!$       WRITE(nout,*) '---                            => Vehkamaeki et al., 2002'
!!$       ELSE IF (nnucl==2)  THEN
!!$       WRITE(nout,*) '---                  nnucl  = ', nnucl
!!$       WRITE(nout,*) '---                            => Kulmala et al., 1998'
!!$       ELSE IF (lsnucl .AND. (nnucl/=1 .OR. nnucl/=2)) THEN
!!$       WRITE(nout,*) '---    Error:'
!!$       CALL finish('setaerom7',  'nucleation requested but no scheme selected')
!!$       END IF
!!$       WRITE(nout,*) '---'
!!$       IF (lmass_diag) THEN
!!$       WRITE(nout,*) '---    Mass balance check in m7_interface activated'
!!$       ELSE
!!$       WRITE(nout,*) '---    Mass balance check in m7_interface deactivated'
!!$       END IF
!!$       WRITE(nout,*) '---'
!!$       WRITE(nout,*) '----------------------------------------------------------'
!!$       WRITE(nout,*) '----------------------------------------------------------'
!!$       WRITE(nout,*) ''
!!$       WRITE(nout,*) ''
!!$    END IF
!!$
!!$  END SUBROUTINE setaeroM7


END MODULE mo_aero_m7
