!#######################################################################
!
! Gas-phase chemistry routines.
!
! Do not edit, automatically generated ...
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

module LE_Chem_Work

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  nreac

  public  ::  naux
  public  ::  iaux_air
  public  ::  iaux_H2O
  public  ::  iaux_p
  public  ::  iaux_T
  public  ::  iaux_zen
  public  ::  iaux_cloud
  public  ::  iaux_cldfac
  public  ::  iaux_rh
  public  ::  iaux_cldsulf
  public  ::  iaux_hetn2o5
  public  ::  iaux_nox0
  public  ::  iaux_rat75
  public  ::  iaux_rat10
  public  ::  iaux_alpha8
  public  ::  iaux_alpha10
  public  ::  iaux_ppm_to_mlccm3
  public  ::  iaux_ppb_to_mlccm3

  public  ::  cldadj

  public  ::  LE_Chem_Work_Init, LE_Chem_Work_Done
  public  ::  LE_Chem_Work_Rates
  public  ::  LE_Chem_Work_Iter

  ! reaction indices needed for M7:
  public  ::  ireac_M7_R71, ireac_M7_R72

  ! reaction indices needed for VBS:
  public  ::  ireac_VBS_R63, ireac_VBS_R65, ireac_VBS_R66

  ! Santiago, Reaction indice needed for SO2 perturbation
  public    ::  ireac_R71

  ! --- const ------------------------------------

  ! module name:
  character(len=*), parameter  ::  mname = 'LE_Chem_Work'

  ! number of gas-phase reactions:
  integer, parameter  ::  nreac = 85

  ! auxilary reaction rate parameters :
  integer, parameter  ::  naux = 17
  integer, parameter  ::  iaux_air      =  1  ! air (concentration)
  integer, parameter  ::  iaux_H2O      =  2  ! water vapour (concentration)
  integer, parameter  ::  iaux_p        =  3  ! air pressure (Pa)
  integer, parameter  ::  iaux_T        =  4  ! temperature (K)
  integer, parameter  ::  iaux_zen      =  5  ! solar zenith angle ?
  integer, parameter  ::  iaux_cloud    =  6  ! cloud fraction [0-1]
  integer, parameter  ::  iaux_cldfac   =  7  ! cloud fraction applied to photolysis rates [0-1]
  integer, parameter  ::  iaux_rh       =  8  ! relative humidity at 2m
  integer, parameter  ::  iaux_cldsulf  =  9  ! cloud factor for sulf reaction
  integer, parameter  ::  iaux_hetn2o5  = 10  ! heterogene chemistry
  integer, parameter  ::  iaux_nox0     = 11  ! chem_solve_nox_value
  integer, parameter  ::  iaux_rat75    = 12  ! ratio between reaction rates for lumping reactions
  integer, parameter  ::  iaux_rat10    = 13  ! ratio between reaction rates for lumping reactions
  integer, parameter  ::  iaux_alpha8   = 14  ! factor used for reverse reactions of ozone photolysis
  integer, parameter  ::  iaux_alpha10  = 15  ! factor used for lumping reactions
  integer, parameter  ::  iaux_ppm_to_mlccm3 = 16  ! factor used for conversion from ppm to mlc/cm3
  integer, parameter  ::  iaux_ppb_to_mlccm3 = 17  ! factor used for conversion from ppb to mlc/cm3

  ! cldadj(11) array of k1 adjustments for cloud cover (0-10 clear to overcast)
  ! source: maul,1980
  real, parameter :: cldadj(0:10) = &
        (/1.0,1.00,1.00,0.79,0.75,0.72,0.68,0.62,0.53,0.41,0.35/)

  ! reaction indices:
  integer, parameter  ::  ireac_R1       =  1  !  NO2             ->  NO + O3
  integer, parameter  ::  ireac_R3       =  2  !  O3+NO           ->  NO2
  integer, parameter  ::  ireac_R7       =  3  !  NO2+O3          ->  NO3
  integer, parameter  ::  ireac_R8       =  4  !  O3              ->  ftmp_oh  * OH + ftmp_o3  * O3
  integer, parameter  ::  ireac_R10      =  5  !  O3+OH           ->  HO2
  integer, parameter  ::  ireac_R11      =  6  !  O3+HO2          ->  OH
  integer, parameter  ::  ireac_R12      =  7  !  NO3+NO          ->  2*NO2
  integer, parameter  ::  ireac_R13      =  8  !  NO3+NO2         ->  NO + NO2
  integer, parameter  ::  ireac_R14      =  9  !  NO3+NO2         ->  N2O5
  integer, parameter  ::  ireac_R16      = 10  !  N2O5            ->  NO3 + NO2
  integer, parameter  ::  ireac_R17      = 11  !  NO+NO2+H2O      ->  2*HNO2
  integer, parameter  ::  ireac_R18      = 12  !  HNO2+HNO2       ->  NO + NO2
  integer, parameter  ::  ireac_R19      = 13  !  HNO2            ->  NO + OH
  integer, parameter  ::  ireac_R20      = 14  !  NO2+OH          ->  HNO3
  integer, parameter  ::  ireac_R21      = 15  !  NO+OH           ->  HNO2
  integer, parameter  ::  ireac_R22      = 16  !  HO2+NO          ->  OH + NO2
  integer, parameter  ::  ireac_R23      = 17  !  NO+NO           ->  2*NO2
  integer, parameter  ::  ireac_R26      = 18  !  OH+HNO2         ->  NO2
  integer, parameter  ::  ireac_R27      = 19  !  NO3             ->  NO2 + O3
  integer, parameter  ::  ireac_R28      = 20  !  NO3             ->  NO
  integer, parameter  ::  ireac_R29      = 21  !  HO2+HO2         ->  H2O2
  integer, parameter  ::  ireac_R30      = 22  !  HO2+HO2+H2O     ->  H2O2
  integer, parameter  ::  ireac_R31      = 23  !  OH+CO           ->  HO2 + CO2
  integer, parameter  ::  ireac_R32      = 24  !  FORM+OH         ->  HO2 + CO
  integer, parameter  ::  ireac_R33      = 25  !  FORM            ->  2*HO2 + CO
  integer, parameter  ::  ireac_R34      = 26  !  FORM            ->  CO
  integer, parameter  ::  ireac_R36      = 27  !  FORM+NO3        ->  HNO3 + HO2 + CO
  integer, parameter  ::  ireac_R38      = 28  !  ALD+OH          ->  C2O3
  integer, parameter  ::  ireac_R39      = 29  !  ALD+NO3         ->  C2O3 + HNO3
  integer, parameter  ::  ireac_R40      = 30  !  ALD             ->  CO + FORM + 2*HO2 + XO2
  integer, parameter  ::  ireac_R42      = 31  !  C2O3+NO         ->  NO2 + XO2 + FORM + HO2
  integer, parameter  ::  ireac_R43      = 32  !  C2O3+NO2        ->  PAN
  integer, parameter  ::  ireac_R44      = 33  !  PAN             ->  C2O3 + NO2
  integer, parameter  ::  ireac_R45      = 34  !  C2O3+C2O3       ->  2*FORM + 2*XO2 + 2*HO2
  integer, parameter  ::  ireac_R46      = 35  !  C2O3+HO2        ->  0.79 * FORM + 0.79 * HO2 + 0.79 * XO2 + 0.79 *OH
  integer, parameter  ::  ireac_R47      = 36  !  MGLY            ->  C2O3 + HO2 + CO
  integer, parameter  ::  ireac_R48      = 37  !  OH+MGLY         ->  XO2 + C2O3
  integer, parameter  ::  ireac_R49      = 38  !  CH4+OH          ->  XO2 + FORM + HO2
  integer, parameter  ::  ireac_R50      = 39  !  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
  integer, parameter  ::  ireac_R52      = 40  !  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
  integer, parameter  ::  ireac_R53      = 41  !  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
  integer, parameter  ::  ireac_R54      = 42  !  NO3+OLE         ->  0.910*HO2 + 0.910*XO2 + 0.090*XO2N - PAR
  integer, parameter  ::  ireac_R56      = 43  !  OH+ETH          ->  XO2 + 2*FORM + HO2
  integer, parameter  ::  ireac_R57      = 44  !  O3+ETH          ->  FORM + 0.370*CO + 0.130*HO2
  integer, parameter  ::  ireac_R58      = 45  !  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
  integer, parameter  ::  ireac_R59      = 46  !  CRES+NO3        ->  CRO + HNO3
  integer, parameter  ::  ireac_R60      = 47  !  CRO+NO2         ->  
  integer, parameter  ::  ireac_R61      = 48  !  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
  integer, parameter  ::  ireac_R62      = 49  !  OH+CRES         ->  0.400*CRO + 0.600*XO2N + 0.600*HO2
  integer, parameter  ::  ireac_R63      = 50  !  XO2+NO          ->  NO2
  integer, parameter  ::  ireac_R64      = 51  !  XO2N+NO         ->  
  integer, parameter  ::  ireac_R65      = 52  !  XO2+XO2         ->  
  integer, parameter  ::  ireac_R66      = 53  !  XO2+HO2         ->  
  integer, parameter  ::  ireac_R67      = 54  !  XO2N+HO2        ->  
  integer, parameter  ::  ireac_R68      = 55  !  XO2N+XO2N       ->  
  integer, parameter  ::  ireac_R69      = 56  !  XO2+XO2N        ->  
  integer, parameter  ::  ireac_R71      = 57  !  SO2+OH          ->  SO4a_f + HO2
  integer, parameter  ::  ireac_R73      = 58  !  OH+H2O2         ->  HO2
  integer, parameter  ::  ireac_R74      = 59  !  H2O2            ->  2*OH
  integer, parameter  ::  ireac_R75      = 60  !  HNO3            ->  OH + NO2
  integer, parameter  ::  ireac_R76      = 61  !  OH+HNO3         ->  NO3
  integer, parameter  ::  ireac_R80      = 62  !  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N
  integer, parameter  ::  ireac_R81      = 63  !  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
  integer, parameter  ::  ireac_R82      = 64  !  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
  integer, parameter  ::  ireac_R95      = 65  !  TO2+NO        ->  NO2 + 0.900*HO2 + 0.900*OPEN
  integer, parameter  ::  ireac_R96      = 66  !  TO2           ->  CRES + HO2
  integer, parameter  ::  ireac_R97      = 67  !  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
  integer, parameter  ::  ireac_R98      = 68  !  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
  integer, parameter  ::  ireac_R99      = 69  !  OPEN          ->  C2O3 + HO2 + CO
  integer, parameter  ::  ireac_R100     = 70  !  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR
  integer, parameter  ::  ireac_R101     = 71  !  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
  integer, parameter  ::  ireac_R102     = 72  !  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
  integer, parameter  ::  ireac_R103     = 73  !  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
  integer, parameter  ::  ireac_R104     = 74  !  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3
  integer, parameter  ::  ireac_RH1f     = 75  !  SO4a_f + N2O5    ->  SO4a_f + 2*HNO3
  integer, parameter  ::  ireac_RH2f2    = 76  !  Na_f2 + N2O5      ->  Na_f2 + 2*NO3a_f
  integer, parameter  ::  ireac_RH2f1    = 77  !  Na_f1 + N2O5      ->  Na_f1 + 2*NO3a_f
  integer, parameter  ::  ireac_RH2c3    = 78  !  Na_c3 + N2O5      ->  Na_c3 + 2*NO3a_c
  integer, parameter  ::  ireac_RH2c2    = 79  !  Na_c2 + N2O5      ->  Na_c2 + 2*NO3a_c
  integer, parameter  ::  ireac_RH2c1    = 80  !  Na_c1 + N2O5      ->  Na_c1 + 2*NO3a_c
  integer, parameter  ::  ireac_RH3f1    = 81  !  Na_f1 + HNO3      ->  Na_f1 + NO3a_f
  integer, parameter  ::  ireac_RH3f2    = 82  !  Na_f2 + HNO3      ->  Na_f2 + NO3a_f
  integer, parameter  ::  ireac_RH3c3    = 83  !  Na_c3 + HNO3      ->  Na_c3 + NO3a_c
  integer, parameter  ::  ireac_RH3c2    = 84  !  Na_c2 + HNO3      ->  Na_c2 + NO3a_c
  integer, parameter  ::  ireac_RH3c1    = 85  !  Na_c1 + HNO3      ->  Na_c1 + NO3a_c
  integer             ::  ireac_VBS_R63  = -1  ! dummy for hardcoded chemistry
  integer             ::  ireac_VBS_R65  = -1  ! dummy for hardcoded chemistry
  integer             ::  ireac_VBS_R66  = -1  ! dummy for hardcoded chemistry
  integer             ::  ireac_M7_R71   = -1  ! dummy for hardcoded chemistry
  integer             ::  ireac_M7_R72   = -1  ! dummy for hardcoded chemistry



contains


  ! ====================================================================


  subroutine LE_Chem_Work_Init( status )

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Chem_Work_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine LE_Chem_Work_Init


  ! ***


  subroutine LE_Chem_Work_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Chem_Work_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine LE_Chem_Work_Done


  ! ***


  !
  ! Reaction rate expressions.
  !
  !  A x <J_label>  k = A * j
  !
  !  A  @E          k = A * exp(-E/T)
  !
  !  A^B@E          k = A x (T/300)B x exp(-E/T)
  !
  !  k1 & k2 & F & n
  !                         k1 [M]                          log10(k1[M]/k2) +2  -1
  !                 k = [------------] F^G   ,   G = [ 1 + (---------------)  ]
  !                      1 + k1[M]/k2                              n
  !
  !     NOTE: in original document, the formula for G included a power '-2' ;
  !     Comparison with the TROE function used in CB99, and other references
  !     (e.g. http://jpldataeval.jpl. nasa.gov/pdf/JPL_15_AllInOne.pdf ),
  !     show that this power should be '+2' .
  !     Ferd Sauter, RIVM, july 2009
  !
  !     NOTE: the 'log' in the original document should be interpreted
  !     as the 10-log, not the e-log .
  !       Ferd Sauter, RIVM, july 2009
  !     However, this seems only true if 'n=1' is used.
  !     With 'n=log(10)' as seen in the original LOTOS-EUROS code
  !     it should be 'log' .
  !       Arjo Segers, TNO, sept 2011
  !
  !  %2 k1 & k2 & k3
  !                             k3[M]
  !                 k = k1 + ------------
  !                          1 + k3[M]/k2
  !
  !  %3 k1 & k2     k = k1 + k2[M]
  !

  elemental function rate1( k1, k2, F, n, M )

    ! --- in/out ---------------------------------

    real                      ::  rate1
    real, intent(in)          ::  k1
    real, intent(in)          ::  k2
    real, intent(in)          ::  F
    real, intent(in)          ::  n
    real, intent(in)          ::  M

    ! --- local ----------------------------------

    real         ::  G

    ! --- begin ----------------------------------

    !! CB05 definition, usually n=1 :
    !G = 1.0 / ( 1.0 + ( log10(k1*M/k2)/n )**2 )

    ! LOTOS-EUROS definition, usually n=log(10) :
    G = 1.0 / ( 1.0 + ( log(k1*M/k2)/n )**2 )

    rate1 = ( (k1*M)/(1.0+k1*M/k2) ) * F**G

  end function rate1

  ! *

  elemental function rate2( k1, k2, k3, M )

    ! --- in/out ---------------------------------

    real                      ::  rate2
    real, intent(in)          ::  k1
    real, intent(in)          ::  k2
    real, intent(in)          ::  k3
    real, intent(in)          ::  M

    ! --- begin ----------------------------------

    rate2 = k1 + (k3*M)/(1.0+k3*M/k2)

  end function rate2

  ! *

  elemental function rate3( k1, k2, M )

    ! --- in/out ---------------------------------

    real                      ::  rate3
    real, intent(in)          ::  k1
    real, intent(in)          ::  k2
    real, intent(in)          ::  M

    ! --- begin ----------------------------------

    rate3 = k1 + k2*M

  end function rate3

  ! *

  elemental function phux( x,y,z, zen )

    ! --- in/out ---------------------------------

    real                      ::  phux
    real, intent(in)          ::  x, y, z
    real, intent(in)          ::  zen

    ! --- local ----------------------------------

    real    ::  arg
    real    ::  fac

    ! --- begin ----------------------------------

    arg = zen*z
    if ( arg < 1.57 ) then
      fac = max( -10.0, y*(1.0-1.0/cos(arg)) )
      phux = x*exp(fac)
    else
      phux = 0.0
    end if

  end function phux


  ! ====================================================================


  ! computation of the photochemical reaction rates
  ! and the temperature dependent rates
  ! unit: ppb**(-n) min**(-1)
  !   with 'n' the order of reaction, thus n=1,2, ..

  subroutine LE_Chem_Work_Rates( rk, aux, rk_het )

    use LE_Chem_Hetro, only : nreac_het
    use LE_Chem_Hetro, only : ireac_N2O5_NH4HSO4a_f
    use LE_Chem_Hetro, only : ireac_N2O5_ss_f
    use LE_Chem_Hetro, only : ireac_N2O5_ss_c
    use LE_Chem_Hetro, only : ireac_HNO3_ss_f
    use LE_Chem_Hetro, only : ireac_HNO3_ss_c

    ! --- in/out ---------------------------------

    real(8), intent(out)       ::  rk(nreac)         ! reaction rates (1/ppb**n/min)
    real(8), intent(inout)     ::  aux(naux)         ! auxilary values
    real(8), intent(in)        ::  rk_het(nreac_het) ! hetr. reac. rates  (1/ppb**n/min)

    ! --- const ----------------------------------

    ! Constant used for 10-based logarithms in reaction rates.
    ! Should be replaced by the 10-log available in Fortran:
    !    log(x)/ln10 = log(x)/log(10) = log10(x)
    real, parameter   ::  ln10 = 2.3025851

    ! --- local ----------------------------------

    ! auxilary values
    real              ::  air
    real              ::  ch2o
    real              ::  T
    real              ::  zen
    real              ::  cldfac
    real              ::  cldsulf
    real              ::  rh
    real              ::  invT
    real              ::  Tdiv300

    ! --- begin ----------------------------------

    ! auxilary concentrations:
    air     = aux(iaux_air)      ! (mlc/cm3)
    ch2o    = aux(iaux_H2O)      ! water concentration

    ! extract auxilary values:
    T       = aux(iaux_T     )   ! temperature (K)
    zen     = aux(iaux_zen   )   ! solar zenith angle ?
    cldfac  = aux(iaux_cldfac)   ! cloud cover factor applied to photolysis rates
    cldsulf = aux(iaux_cldsulf)  ! cloud cover factor for sulf reaction
    rh      = aux(iaux_rh    )   ! relative humidity at 2m

    ! precompute inversions as a much faster multiplication ...
    invT    = 1.0/T
    Tdiv300 = T/300.0

    ! R1 :  NO2             ->  NO + O3
    ! rate expr. : 1.0 x <NO2_SAPRC99>
    rk(ireac_R1) = 1.0 * phux( 1.07e-2, 1.01319, 0.83330, zen ) * cldfac

    ! R3 :  O3+NO           ->  NO2
    ! rate expr. : 1.8E-12 @ 1450
    rk(ireac_R3) = 1.8E-12 * exp(-1450/T)

    ! R7 :  NO2+O3          ->  NO3
    ! rate expr. : 1.2E-13 @ 2450
    rk(ireac_R7) = 1.2E-13 * exp(-2450/T)

    ! R8 :  O3              ->  ftmp_oh  * OH + ftmp_o3  * O3
    ! rate expr. : 1.0 x <O3_O1D_IUPAC05>
    rk(ireac_R8) = phux( 3.22e-5, 4.45037, 0.78028, zen ) * cldfac

    ! R10 :  O3+OH           ->  HO2
    ! rate expr. : 1.6E-12 @ 940
    rk(ireac_R10) = 1.6E-12 * exp(-940/T)

    ! R11 :  O3+HO2          ->  OH
    ! rate expr. : 1.4E-14 @ 580
    rk(ireac_R11) = 1.4E-14 * exp(-580/T)

    ! R12 :  NO3+NO          ->  2*NO2
    ! rate expr. : 1.3E-11 @ -250
    rk(ireac_R12) = 1.3E-11 * exp(250/T)

    ! R13 :  NO3+NO2         ->  NO + NO2
    ! rate expr. : 2.5E-14 @ 1230
    rk(ireac_R13) = 2.5E-14 * exp(-1230/T)

    ! R14 :  NO3+NO2         ->  N2O5
    ! rate expr. : 5.3E-13 @ -256
    rk(ireac_R14) = 5.3E-13 * exp(256/T)

    ! R16 :  N2O5            ->  NO3 + NO2
    ! rate expr. : 3.5E14  @ 10897 & 1.0 x <N2O5_IUPAC05>
    rk(ireac_R16) = 3.5E14 * exp(-10897/T) + phux( 3.79e-5, 1.70537, 0.80153, zen ) * cldfac

    ! R17 :  NO+NO2+H2O      ->  2*HNO2
    ! rate expr. : 4.6D-46 @ -6348
    rk(ireac_R17) = 4.6D-46 * exp(6348/T)

    ! R18 :  HNO2+HNO2       ->  NO + NO2
    ! rate expr. : 1.0E-20
    rk(ireac_R18) = 1.0E-20

    ! R19 :  HNO2            ->  NO + OH
    ! rate expr. : 1.0 x <HNO2_IUPAC05>
    rk(ireac_R19) = phux( 8.96e-4, 0.99438, 0.83295, zen ) * cldfac

    ! R20 :  NO2+OH          ->  HNO3
    ! rate expr. : 1.0E-12 @ -713
    rk(ireac_R20) = 1.0E-12 * exp(713/T)

    ! R21 :  NO+OH           ->  HNO2
    ! rate expr. : 4.5E-13 @ -806
    rk(ireac_R21) = 4.5E-13 * exp(806/T)

    ! R22 :  HO2+NO          ->  OH + NO2
    ! rate expr. : 3.7E-12 @ -240
    rk(ireac_R22) = 3.7E-12 * exp(240/T)

    ! R23 :  NO+NO           ->  2*NO2
    ! rate expr. : 1.8E-20 @ -530
    rk(ireac_R23) = 1.8E-20 * exp(530/T)

    ! R26 :  OH+HNO2         ->  NO2
    ! rate expr. : 6.6E-12
    rk(ireac_R26) = 6.6E-12

    ! R27 :  NO3             ->  NO2 + O3
    ! rate expr. : 1.0 x <NO3NO2_SAPRC99>
    rk(ireac_R27) = phux( 2.73e-1, 0.29327, 0.92401, zen ) * cldfac

    ! R28 :  NO3             ->  NO
    ! rate expr. : 1.0 x <NO3NO_SAPRC99>
    rk(ireac_R28) = phux( 2.74e-2, 0.26226, 0.92849, zen ) * cldfac

    ! R29 :  HO2+HO2         ->  H2O2
    ! rate expr. : 5.9E-14 @ -1150
    rk(ireac_R29) = 5.9E-14 * exp(1150/T)

    ! R30 :  HO2+HO2+H2O     ->  H2O2
    ! rate expr. : 2.2E-38 @ -5800
    rk(ireac_R30) = 2.2E-38 * exp(5800/T)

    ! R31 :  OH+CO           ->  HO2 + CO2
    ! rate expr. : 2.2E-13
    rk(ireac_R31) = 2.2E-13

    ! R32 :  FORM+OH         ->  HO2 + CO
    ! rate expr. : 1.6E-11 @ 110
    rk(ireac_R32) = 1.6E-11 * exp(-110/T)

    ! R33 :  FORM            ->  2*HO2 + CO
    ! rate expr. : 1.0 x <HCHO_R_SAPRC99>
    rk(ireac_R33) = 1.0 * phux( 4.05e-5, 2.06917, 0.80267, zen ) * cldfac

    ! R34 :  FORM            ->  CO
    ! rate expr. : 1.0 x <HCHO_M_SAPRC99>
    rk(ireac_R34) = phux( 4.92e-5, 1.60973, 0.80184, zen ) * cldfac

    ! R36 :  FORM+NO3        ->  HNO3 + HO2 + CO
    ! rate expr. : 6.3E-16
    rk(ireac_R36) = 6.3E-16

    ! R38 :  ALD+OH          ->  C2O3
    ! rate expr. : 7.0E-12 @ -250
    rk(ireac_R38) = 7.0E-12 * exp(250/T)

    ! R39 :  ALD+NO3         ->  C2O3 + HNO3
    ! rate expr. : 1.4E-12 @ 1900
    rk(ireac_R39) = 1.4E-12 * exp(-1900/T)

    ! R40 :  ALD             ->  CO + FORM + 2*HO2 + XO2
    ! rate expr. : 1.0 x <CCHO_R_SAPRC99>
    rk(ireac_R40) = phux( 5.40e-6, 2.52915, 0.79722, zen ) * cldfac

    ! R42 :  C2O3+NO         ->  NO2 + XO2 + FORM + HO2
    ! rate expr. : 5.4E-12 @ -250
    rk(ireac_R42) = 5.4E-12 * exp(250/T)

    ! R43 :  C2O3+NO2        ->  PAN
    ! rate expr. : 9.7E-29^-5.6 & 9.3E-12^1.5 & 0.6 & ln10
    rk(ireac_R43) = rate1( (9.7E-29 * (Tdiv300)**(-5.6)), (9.3E-12 * (Tdiv300)**(1.5)), 0.6, ln10, air )

    ! R44 :  PAN             ->  C2O3 + NO2
    ! rate expr. : 9.7E-29^-5.6 & 9.3E-12^1.5 & 0.6 & ln10 / 9.0E-29 @ -14000
    rk(ireac_R44) = rate1( (9.7E-29 * (Tdiv300)**(-5.6)), (9.3E-12 * (Tdiv300)**(1.5)), 0.6, ln10, air ) / &
             (9.0E-29 * exp(14000/T))

    ! R45 :  C2O3+C2O3       ->  2*FORM + 2*XO2 + 2*HO2
    ! rate expr. : 2.5E-12
    rk(ireac_R45) = 2.5E-12

    ! R46 :  C2O3+HO2        ->  0.79 * FORM + 0.79 * HO2 + 0.79 * XO2 + 0.79 *OH
    ! rate expr. : 4.3E-13 @ -1040
    rk(ireac_R46) = 4.3E-13 * exp(1040/T)

    ! R47 :  MGLY            ->  C2O3 + HO2 + CO
    ! rate expr. : 0.02 x <NO2_SAPRC99>
    rk(ireac_R47) = 0.02 * phux( 1.07e-2, 1.01319, 0.83330, zen ) * cldfac

    ! R48 :  OH+MGLY         ->  XO2 + C2O3
    ! rate expr. : 1.7E-11
    rk(ireac_R48) = 1.7E-11

    ! R49 :  CH4+OH          ->  XO2 + FORM + HO2
    ! rate expr. : 2.6E-12 @ 1800
    rk(ireac_R49) = 2.6E-12 * exp(-1800/T)

    ! R50 :  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
    ! rate expr. : 8.1E-13
    rk(ireac_R50) = 8.1E-13

    ! R52 :  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
    ! rate expr. : 5.2E-12 @ -504
    rk(ireac_R52) = 5.2E-12 * exp(504/T)

    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! rate expr. : 1.4E-14 @ 2105
    rk(ireac_R53) = 1.4E-14 * exp(-2105/T)

    ! R54 :  NO3+OLE         ->  0.910*HO2 + 0.910*XO2 + 0.090*XO2N - PAR
    ! rate expr. : 7.7E-15
    rk(ireac_R54) = 7.7E-15

    ! R56 :  OH+ETH          ->  XO2 + 2*FORM + HO2
    ! rate expr. : 2.0E-12 @ -411
    rk(ireac_R56) = 2.0E-12 * exp(411/T)

    ! R57 :  O3+ETH          ->  FORM + 0.370*CO + 0.130*HO2
    ! rate expr. : 1.3E-14 @ 2633
    rk(ireac_R57) = 1.3E-14 * exp(-2633/T)

    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! rate expr. : 2.1E-12 @ -322
    rk(ireac_R58) = 2.1E-12 * exp(322/T)

    ! R59 :  CRES+NO3        ->  CRO + HNO3
    ! rate expr. : 2.2E-11
    rk(ireac_R59) = 2.2E-11

    ! R60 :  CRO+NO2         ->  
    ! rate expr. : 1.4E-11
    rk(ireac_R60) = 1.4E-11

    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
    ! rate expr. : 1.7E-11 @ -116
    rk(ireac_R61) = 1.7E-11 * exp(116/T)

    ! R62 :  OH+CRES         ->  0.400*CRO + 0.600*XO2N + 0.600*HO2
    ! rate expr. : 4.1E-11
    rk(ireac_R62) = 4.1E-11

    ! R63 :  XO2+NO          ->  NO2
    ! rate expr. : 8.1E-12
    rk(ireac_R63) = 8.1E-12

    ! R64 :  XO2N+NO         ->  
    ! rate expr. : 8.1E-13
    rk(ireac_R64) = 8.1E-13

    ! R65 :  XO2+XO2         ->  
    ! rate expr. : 8.5E-15 @ -1300
    rk(ireac_R65) = 8.5E-15 * exp(1300/T)

    ! R66 :  XO2+HO2         ->  
    ! rate expr. : 7.6E-14 @ -1300
    rk(ireac_R66) = 7.6E-14 * exp(1300/T)

    ! R67 :  XO2N+HO2        ->  
    ! rate expr. : 7.6E-14 @ -1300
    rk(ireac_R67) = 7.6E-14 * exp(1300/T)

    ! R68 :  XO2N+XO2N       ->  
    ! rate expr. : 8.5E-15 @ -1300
    rk(ireac_R68) = 8.5E-15 * exp(1300/T)

    ! R69 :  XO2+XO2N        ->  
    ! rate expr. : 3.4E-14 @ -1300
    rk(ireac_R69) = 3.4E-14 * exp(1300/T)

    ! R71 :  SO2+OH          ->  SO4a_f + HO2
    ! rate expr. : 1.0E-12
    rk(ireac_R71) = 1.0E-12

    ! R73 :  OH+H2O2         ->  HO2
    ! rate expr. : 3.1E-12 @ 187
    rk(ireac_R73) = 3.1E-12 * exp(-187/T)

    ! R74 :  H2O2            ->  2*OH
    ! rate expr. : 1.0 x <H2O2_SAPRC99>
    rk(ireac_R74) = phux( 7.78e-6, 1.91463, 0.79810, zen ) * cldfac

    ! R75 :  HNO3            ->  OH + NO2
    ! rate expr. : 1.0 x <HNO3_IUPAC05>
    rk(ireac_R75) = phux( 5.48e-7, 2.86922, 0.79561, zen ) * cldfac

    ! R76 :  OH+HNO3         ->  NO3
    ! rate expr. : 5.1E-15 @ -1000
    rk(ireac_R76) = 5.1E-15 * exp(1000/T)

    ! R80 :  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N
    ! rate expr. : 2.7E-11 @ -407.6
    rk(ireac_R80) = 2.7E-11 * exp(407.6/T)

    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! rate expr. : 8.5E-15 @ 1912
    rk(ireac_R81) = 8.5E-15 * exp(-1912/T)

    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! rate expr. : 3.3E-12 @ 448
    rk(ireac_R82) = 3.3E-12 * exp(-448/T)

    ! R95 :  TO2+NO        ->  NO2 + 0.900*HO2 + 0.900*OPEN
    ! rate expr. : 8.1E-12
    rk(ireac_R95) = 8.1E-12

    ! R96 :  TO2           ->  CRES + HO2
    ! rate expr. : 4.2
    rk(ireac_R96) = 4.2

    ! R97 :  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
    ! rate expr. : 3.0E-11
    rk(ireac_R97) = 3.0E-11

    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! rate expr. : 5.4E-17 @ 500
    rk(ireac_R98) = 5.4E-17 * exp(-500/T)

    ! R99 :  OPEN          ->  C2O3 + HO2 + CO
    ! rate expr. : 6.0 x <HCHO_R_SAPRC99>
    rk(ireac_R99) = 6.0 * phux( 4.05e-5, 2.06917, 0.80267, zen ) * cldfac

    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR
    ! rate expr. : 1.5E-19
    rk(ireac_R100) = 1.5E-19

    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! rate expr. : 3.4E-11
    rk(ireac_R101) = 3.4E-11

    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! rate expr. : 7.1E-18
    rk(ireac_R102) = 7.1E-18

    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! rate expr. : 1.0E-15
    rk(ireac_R103) = 1.0E-15

    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3
    ! rate expr. : 1.7E-4 x <NO2_SAPRC99>
    rk(ireac_R104) = 1.7E-4 * phux( 1.07e-2, 1.01319, 0.83330, zen ) * cldfac

    ! RH1f :  SO4a_f + N2O5    ->  SO4a_f + 2*HNO3
    ! rate expr. : rk_het(ireac_N2O5_NH4HSO4a_f)
    rk(ireac_RH1f) = rk_het(ireac_N2O5_NH4HSO4a_f)

    ! RH2f2 :  Na_f2 + N2O5      ->  Na_f2 + 2*NO3a_f
    ! rate expr. : rk_het(ireac_N2O5_ss_f)
    rk(ireac_RH2f2) = rk_het(ireac_N2O5_ss_f)

    ! RH2f1 :  Na_f1 + N2O5      ->  Na_f1 + 2*NO3a_f
    ! rate expr. : rk_het(ireac_N2O5_ss_f)
    rk(ireac_RH2f1) = rk_het(ireac_N2O5_ss_f)

    ! RH2c3 :   Na_c3 + N2O5      ->  Na_c3 + 2*NO3a_c
    ! rate expr. : rk_het(ireac_N2O5_ss_c)
    rk(ireac_RH2c3) = rk_het(ireac_N2O5_ss_c)

    ! RH2c2 :  Na_c2 + N2O5      ->  Na_c2 + 2*NO3a_c
    ! rate expr. : rk_het(ireac_N2O5_ss_c)
    rk(ireac_RH2c2) = rk_het(ireac_N2O5_ss_c)

    ! RH2c1 : Na_c1 + N2O5      ->  Na_c1 + 2*NO3a_c
    ! rate expr. : rk_het(ireac_N2O5_ss_c)
    rk(ireac_RH2c1) = rk_het(ireac_N2O5_ss_c)

    ! RH3f1 :   Na_f1 + HNO3     ->  Na_f1 + NO3a_f
    ! rate expr. : rk_het(ireac_HNO3_ss_f)
    rk(ireac_RH3f1) = rk_het(ireac_HNO3_ss_f)

    ! RH3f2 :    Na_f2 + HNO3      ->  Na_f2 + NO3a_f
    ! rate expr. : rk_het(ireac_HNO3_ss_f)
    rk(ireac_RH3f2) = rk_het(ireac_HNO3_ss_f)

    ! RH3c3 :    Na_c3 + HNO3      ->  Na_c3 + NO3a_c
    ! rate expr. : rk_het(ireac_HNO3_ss_c)
    rk(ireac_RH3c3) = rk_het(ireac_HNO3_ss_c)

    ! RH3c2 :   Na_c2 + HNO3     ->  Na_c2 + NO3a_c
    ! rate expr. : rk_het(ireac_HNO3_ss_c)
    rk(ireac_RH3c2) = rk_het(ireac_HNO3_ss_c)

    ! RH3c1 :  Na_c1 + HNO3    ->  Na_c1 + NO3a_c
    ! rate expr. : rk_het(ireac_HNO3_ss_c)
    rk(ireac_RH3c1) = rk_het(ireac_HNO3_ss_c)

    !
    ! Set some extra auxilary values, these are results of lumping reactions
    !
    ! For example the reaction:
    !       O + H2O   --> 2 OH
    ! is involved in the photolysis of ozone:
    !       O3        --> O2 + O
    !
    ! Further the reactions:
    !       ROR       --> 1.10 ALD + 0.96 XO2 + 0.94 HO2 +
    !                     0.04 XO2N + 0.02 ROR - 2.10 PAR
    !
    !       ROR       --> HO2
    !
    !       ROR + NO2 -->
    !
    ! are involved in the reaction:
    !       PAR + OH  --> 0.87 XO2 + 0.13 XO2N + 0.11 HO2 +
    !                     0.11 ALD + 0.76 ROR  - 0.11 PAR
    !
    ! The fractions below are the result of the different reaction rates
    ! (ratios depends on temperature):
    !
    ! rat 10 must be in order of mlc/cm3 thus multiplied with factor ppm_to_mlccm3  
    ! This factor is 2.4621e13 with st temp (298K ) and st pres (1.013hPa)
    aux(iaux_rat10  ) = 42900.0 * 1.e+06 * aux(iaux_ppm_to_mlccm3)
    aux(iaux_alpha8 ) = cH2O * 326000.0 / ( cH2O * 326000.0 + aux(iaux_rat10) )
    aux(iaux_rat75  ) = 90000.0 * exp( 7000.0*(1/298.0-1.0/T) )
    aux(iaux_alpha10) = aux(iaux_rat75) / ( aux(iaux_rat75) + 390000.0 )

  end subroutine LE_Chem_Work_Rates


  ! ***


  pure subroutine LE_Chem_Work_Iter( n, y, ysum, gdt, yp, yl, Q, &
                                       nreac, rk, naux, aux )

    use Indices

    ! --- in/out ---------------------------------

    integer, intent(in)     ::  n             ! number of concentrations
    real, intent(inout)     ::  y(1:n)        ! concentrations  (conc)
    real, intent(in)        ::  ysum(1:n)     ! work array
    real, intent(in)        ::  gdt           ! sec
    real, intent(out)       ::  yp(1:n)       ! production rate (conc/s)
    real, intent(out)       ::  yl(1:n)       ! loss rate (conc/s)
    real, intent(in)        ::  Q(1:n)        ! source term (conc/s)
    integer, intent(in)     ::  nreac         ! number of reactions
    real(8), intent(in)     ::  rk(nreac)     ! reaction rate coeff.
    integer, intent(in)     ::  naux          ! number of auxilary values
    real(8), intent(in)     ::  aux(naux)     ! auxilary values

    ! --- local ----------------------------------

    real                    :: alpha8
    real                    :: alpha9
    real                    :: alpha10
    real                    :: rat75

    real                    ::  ftmp_NO2 
    real                    ::  ftmp_O3 
    real                    ::  ftmp_ALD 
    real                    ::  ftmp_PAR 
    real                    ::  ftmp_OH 
    real                    ::  ftmp_HO2 
    real                    ::  ftmp_XO2 

    ! --- begin ----------------------------------


    !
    ! NO2
    !

    ! R1 :  NO2             ->  NO + O3
    ! R3 :  O3+NO           ->  NO2
    ! R7 :  NO2+O3          ->  NO3
    ! R12 :  NO3+NO          ->  2*NO2
    ! R13 :  NO3+NO2         ->  NO + NO2
    ! R14 :  NO3+NO2         ->  N2O5
    ! R16 :  N2O5            ->  NO3 + NO2
    ! R17 :  NO+NO2+H2O      ->  2*HNO2
    ! R18 :  HNO2+HNO2       ->  NO + NO2
    ! R20 :  NO2+OH          ->  HNO3
    ! R22 :  HO2+NO          ->  OH + NO2
    ! R23 :  NO+NO           ->  2*NO2
    ! R26 :  OH+HNO2         ->  NO2
    ! R27 :  NO3             ->  NO2 + O3
    ! R42 :  C2O3+NO         ->  NO2 + XO2 + FORM + HO2
    ! R43 :  C2O3+NO2        ->  PAN
    ! R44 :  PAN             ->  C2O3 + NO2
    ! R50 :  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
    ! R60 :  CRO+NO2         ->  
    ! R63 :  XO2+NO          ->  NO2
    ! R75 :  HNO3            ->  OH + NO2
    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! R95 :  TO2+NO        ->  NO2 + 0.900*HO2 + 0.900*OPEN
    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR

    rat75 = aux(iaux_rat75)
    alpha9 = 22.*aux(iaux_ppb_to_mlccm3)/(y(ispec_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))
    ftmp_NO2 = 0.8*alpha9 

    yl(ispec_NO2) = 0.0  &
        + rk(ireac_R1) &
        + rk(ireac_R7) * y(ispec_O3) &
        + rk(ireac_R14) * y(ispec_NO3) &
        + rk(ireac_R17) * y(ispec_NO) * aux(iaux_H2O) &
        + rk(ireac_R20) * y(ispec_OH) &
        + rk(ireac_R43) * y(ispec_C2O3) &
        + rk(ireac_R60) * y(ispec_CRO) &
        + rk(ireac_R100) * y(ispec_ISO)

    yp(ispec_NO2) = Q(ispec_NO2)  &
        + rk(ireac_R3) * y(ispec_O3) * y(ispec_NO) &
        + rk(ireac_R12) * 2 * y(ispec_NO3) * y(ispec_NO) &
        + rk(ireac_R16) * y(ispec_N2O5) &
        + rk(ireac_R18) * y(ispec_HNO2) * y(ispec_HNO2) &
        + rk(ireac_R22) * y(ispec_HO2) * y(ispec_NO) &
        + rk(ireac_R23) * 2 * y(ispec_NO) * y(ispec_NO) &
        + rk(ireac_R26) * y(ispec_OH) * y(ispec_HNO2) &
        + rk(ireac_R27) * y(ispec_NO3) &
        + rk(ireac_R42) * y(ispec_C2O3) * y(ispec_NO) &
        + rk(ireac_R44) * y(ispec_PAN) &
        + rk(ireac_R50) * - ftmp_no2 * y(ispec_PAR) * y(ispec_OH) &
        + rk(ireac_R63) * y(ispec_XO2) * y(ispec_NO) &
        + rk(ireac_R75) * y(ispec_HNO3) &
        + rk(ireac_R82) * y(ispec_NO3) * y(ispec_ISO) &
        + rk(ireac_R95) * y(ispec_TO2) * y(ispec_NO)

    y(ispec_NO2) = max( 0.0, ( ysum(ispec_NO2) + gdt*yp(ispec_NO2) ) / ( 1.0 + gdt*yl(ispec_NO2) ) )

    !
    ! NO
    !

    ! R1 :  NO2             ->  NO + O3
    ! R3 :  O3+NO           ->  NO2
    ! R12 :  NO3+NO          ->  2*NO2
    ! R13 :  NO3+NO2         ->  NO + NO2
    ! R17 :  NO+NO2+H2O      ->  2*HNO2
    ! R18 :  HNO2+HNO2       ->  NO + NO2
    ! R19 :  HNO2            ->  NO + OH
    ! R21 :  NO+OH           ->  HNO2
    ! R22 :  HO2+NO          ->  OH + NO2
    ! R23 :  NO+NO           ->  2*NO2
    ! R28 :  NO3             ->  NO
    ! R42 :  C2O3+NO         ->  NO2 + XO2 + FORM + HO2
    ! R63 :  XO2+NO          ->  NO2
    ! R64 :  XO2N+NO         ->  
    ! R95 :  TO2+NO        ->  NO2 + 0.900*HO2 + 0.900*OPEN
    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR

    yl(ispec_NO) = 0.0  &
        + rk(ireac_R3) * y(ispec_O3) &
        + rk(ireac_R12) * y(ispec_NO3) &
        + rk(ireac_R17) * y(ispec_NO2) * aux(iaux_H2O) &
        + rk(ireac_R21) * y(ispec_OH) &
        + rk(ireac_R22) * y(ispec_HO2) &
        + rk(ireac_R23) * y(ispec_NO) &
        + rk(ireac_R23) * y(ispec_NO) &
        + rk(ireac_R42) * y(ispec_C2O3) &
        + rk(ireac_R63) * y(ispec_XO2) &
        + rk(ireac_R64) * y(ispec_XO2N) &
        + rk(ireac_R95) * y(ispec_TO2)

    yp(ispec_NO) = Q(ispec_NO)  &
        + rk(ireac_R1) * y(ispec_NO2) &
        + rk(ireac_R13) * y(ispec_NO3) * y(ispec_NO2) &
        + rk(ireac_R18) * y(ispec_HNO2) * y(ispec_HNO2) &
        + rk(ireac_R19) * y(ispec_HNO2) &
        + rk(ireac_R28) * y(ispec_NO3) &
        + rk(ireac_R100) * 0.200 * y(ispec_NO2) * y(ispec_ISO)

    y(ispec_NO) = max( 0.0, ( ysum(ispec_NO) + gdt*yp(ispec_NO) ) / ( 1.0 + gdt*yl(ispec_NO) ) )

    !
    ! O3
    !

    ! R1 :  NO2             ->  NO + O3
    ! R3 :  O3+NO           ->  NO2
    ! R7 :  NO2+O3          ->  NO3
    ! R8 :  O3              ->  ftmp_oh  * OH + ftmp_o3  * O3
    ! R10 :  O3+OH           ->  HO2
    ! R11 :  O3+HO2          ->  OH
    ! R27 :  NO3             ->  NO2 + O3
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R57 :  O3+ETH          ->  FORM + 0.370*CO + 0.130*HO2
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD

    alpha8 = aux(iaux_alpha8)
    ftmp_O3 = 1.-alpha8 

    yl(ispec_O3) = 0.0  &
        + rk(ireac_R3) * y(ispec_NO) &
        + rk(ireac_R7) * y(ispec_NO2) &
        + rk(ireac_R8) &
        + rk(ireac_R10) * y(ispec_OH) &
        + rk(ireac_R11) * y(ispec_HO2) &
        + rk(ireac_R53) * y(ispec_OLE) &
        + rk(ireac_R57) * y(ispec_ETH) &
        + rk(ireac_R81) * y(ispec_ISO) &
        + rk(ireac_R98) * y(ispec_OPEN) &
        + rk(ireac_R102) * y(ispec_ISPD)

    yp(ispec_O3) = Q(ispec_O3)  &
        + rk(ireac_R1) * y(ispec_NO2) &
        + rk(ireac_R8) * ftmp_o3 * y(ispec_O3) &
        + rk(ireac_R27) * y(ispec_NO3)

    y(ispec_O3) = max( 0.0, ( ysum(ispec_O3) + gdt*yp(ispec_O3) ) / ( 1.0 + gdt*yl(ispec_O3) ) )

    !
    ! CO
    !

    ! R31 :  OH+CO           ->  HO2 + CO2
    ! R32 :  FORM+OH         ->  HO2 + CO
    ! R33 :  FORM            ->  2*HO2 + CO
    ! R34 :  FORM            ->  CO
    ! R36 :  FORM+NO3        ->  HNO3 + HO2 + CO
    ! R40 :  ALD             ->  CO + FORM + 2*HO2 + XO2
    ! R47 :  MGLY            ->  C2O3 + HO2 + CO
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R57 :  O3+ETH          ->  FORM + 0.370*CO + 0.130*HO2
    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R97 :  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R99 :  OPEN          ->  C2O3 + HO2 + CO
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3

    yl(ispec_CO) = 0.0  &
        + rk(ireac_R31) * y(ispec_OH)

    yp(ispec_CO) = Q(ispec_CO)  &
        + rk(ireac_R32) * y(ispec_FORM) * y(ispec_OH) &
        + rk(ireac_R33) * y(ispec_FORM) &
        + rk(ireac_R34) * y(ispec_FORM) &
        + rk(ireac_R36) * y(ispec_FORM) * y(ispec_NO3) &
        + rk(ireac_R40) * y(ispec_ALD) &
        + rk(ireac_R47) * y(ispec_MGLY) &
        + rk(ireac_R53) * 0.212 * y(ispec_O3) * y(ispec_OLE) &
        + rk(ireac_R57) * 0.370 * y(ispec_O3) * y(ispec_ETH) &
        + rk(ireac_R58) * 1.130 * y(ispec_TOL) * y(ispec_OH) &
        + rk(ireac_R81) * 0.066 * y(ispec_O3) * y(ispec_ISO) &
        + rk(ireac_R97) * 2 * y(ispec_OPEN) * y(ispec_OH) &
        + rk(ireac_R98) * 0.690 * y(ispec_OPEN) * y(ispec_O3) &
        + rk(ireac_R99) * y(ispec_OPEN) &
        + rk(ireac_R101) * 0.334 * y(ispec_OH) * y(ispec_ISPD) &
        + rk(ireac_R102) * 0.225 * y(ispec_O3) * y(ispec_ISPD) &
        + rk(ireac_R103) * 0.643 * y(ispec_NO3) * y(ispec_ISPD) &
        + rk(ireac_R104) * 0.33 * y(ispec_ISPD)

    y(ispec_CO) = max( 0.0, ( ysum(ispec_CO) + gdt*yp(ispec_CO) ) / ( 1.0 + gdt*yl(ispec_CO) ) )

    !
    ! FORM
    !

    ! R32 :  FORM+OH         ->  HO2 + CO
    ! R33 :  FORM            ->  2*HO2 + CO
    ! R34 :  FORM            ->  CO
    ! R36 :  FORM+NO3        ->  HNO3 + HO2 + CO
    ! R40 :  ALD             ->  CO + FORM + 2*HO2 + XO2
    ! R42 :  C2O3+NO         ->  NO2 + XO2 + FORM + HO2
    ! R45 :  C2O3+C2O3       ->  2*FORM + 2*XO2 + 2*HO2
    ! R46 :  C2O3+HO2        ->  0.79 * FORM + 0.79 * HO2 + 0.79 * XO2 + 0.79 *OH
    ! R49 :  CH4+OH          ->  XO2 + FORM + HO2
    ! R52 :  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R56 :  OH+ETH          ->  XO2 + 2*FORM + HO2
    ! R57 :  O3+ETH          ->  FORM + 0.370*CO + 0.130*HO2
    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R80 :  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R97 :  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3

    yl(ispec_FORM) = 0.0  &
        + rk(ireac_R32) * y(ispec_OH) &
        + rk(ireac_R33) &
        + rk(ireac_R34) &
        + rk(ireac_R36) * y(ispec_NO3)

    yp(ispec_FORM) = Q(ispec_FORM)  &
        + rk(ireac_R40) * y(ispec_ALD) &
        + rk(ireac_R42) * y(ispec_C2O3) * y(ispec_NO) &
        + rk(ireac_R45) * 2 * y(ispec_C2O3) * y(ispec_C2O3) &
        + rk(ireac_R46) * 0.79 * y(ispec_C2O3) * y(ispec_HO2) &
        + rk(ireac_R49) * y(ispec_CH4) * y(ispec_OH) &
        + rk(ireac_R52) * y(ispec_OH) * y(ispec_OLE) &
        + rk(ireac_R53) * 0.660 * y(ispec_O3) * y(ispec_OLE) &
        + rk(ireac_R56) * 2 * y(ispec_OH) * y(ispec_ETH) &
        + rk(ireac_R57) * y(ispec_O3) * y(ispec_ETH) &
        + rk(ireac_R58) * 1.130 * y(ispec_TOL) * y(ispec_OH) &
        + rk(ireac_R80) * 0.629 * y(ispec_OH) * y(ispec_ISO) &
        + rk(ireac_R81) * 0.600 * y(ispec_O3) * y(ispec_ISO) &
        + rk(ireac_R97) * y(ispec_OPEN) * y(ispec_OH) &
        + rk(ireac_R98) * 0.700 * y(ispec_OPEN) * y(ispec_O3) &
        + rk(ireac_R101) * 0.167 * y(ispec_OH) * y(ispec_ISPD) &
        + rk(ireac_R102) * 0.150 * y(ispec_O3) * y(ispec_ISPD) &
        + rk(ireac_R103) * 0.282 * y(ispec_NO3) * y(ispec_ISPD) &
        + rk(ireac_R104) * 0.900 * y(ispec_ISPD)

    y(ispec_FORM) = max( 0.0, ( ysum(ispec_FORM) + gdt*yp(ispec_FORM) ) / ( 1.0 + gdt*yl(ispec_FORM) ) )

    !
    ! ALD
    !

    ! R38 :  ALD+OH          ->  C2O3
    ! R39 :  ALD+NO3         ->  C2O3 + HNO3
    ! R40 :  ALD             ->  CO + FORM + 2*HO2 + XO2
    ! R50 :  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
    ! R52 :  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3

    rat75 = aux(iaux_rat75)
    alpha9 = y(i_NO2)*22./(y(i_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))
    alpha10 = aux(iaux_alpha10)
    ftmp_ALD = 0.13 + (0.8*(1.-alpha9)*(0.185+alpha10*1.115)) 

    yl(ispec_ALD) = 0.0  &
        + rk(ireac_R38) * y(ispec_OH) &
        + rk(ireac_R39) * y(ispec_NO3) &
        + rk(ireac_R40)

    yp(ispec_ALD) = Q(ispec_ALD)  &
        + rk(ireac_R50) * ftmp_ald * y(ispec_PAR) * y(ispec_OH) &
        + rk(ireac_R52) * y(ispec_OH) * y(ispec_OLE) &
        + rk(ireac_R53) * 0.500 * y(ispec_O3) * y(ispec_OLE) &
        + rk(ireac_R81) * 0.150 * y(ispec_O3) * y(ispec_ISO) &
        + rk(ireac_R82) * 0.800 * y(ispec_NO3) * y(ispec_ISO) &
        + rk(ireac_R98) * 0.030 * y(ispec_OPEN) * y(ispec_O3) &
        + rk(ireac_R100) * 0.800 * y(ispec_NO2) * y(ispec_ISO) &
        + rk(ireac_R101) * 0.273 * y(ispec_OH) * y(ispec_ISPD) &
        + rk(ireac_R102) * 0.020 * y(ispec_O3) * y(ispec_ISPD) &
        + rk(ireac_R103) * 0.357 * y(ispec_NO3) * y(ispec_ISPD) &
        + rk(ireac_R104) * 0.067 * y(ispec_ISPD)

    y(ispec_ALD) = max( 0.0, ( ysum(ispec_ALD) + gdt*yp(ispec_ALD) ) / ( 1.0 + gdt*yl(ispec_ALD) ) )

    !
    ! PAN
    !

    ! R43 :  C2O3+NO2        ->  PAN
    ! R44 :  PAN             ->  C2O3 + NO2

    yl(ispec_PAN) = 0.0  &
        + rk(ireac_R44)

    yp(ispec_PAN) = Q(ispec_PAN)  &
        + rk(ireac_R43) * y(ispec_C2O3) * y(ispec_NO2)

    y(ispec_PAN) = max( 0.0, ( ysum(ispec_PAN) + gdt*yp(ispec_PAN) ) / ( 1.0 + gdt*yl(ispec_PAN) ) )

    !
    ! MGLY
    !

    ! R47 :  MGLY            ->  C2O3 + HO2 + CO
    ! R48 :  OH+MGLY         ->  XO2 + C2O3
    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD

    yl(ispec_MGLY) = 0.0  &
        + rk(ireac_R47) &
        + rk(ireac_R48) * y(ispec_OH)

    yp(ispec_MGLY) = Q(ispec_MGLY)  &
        + rk(ireac_R58) * 0.560 * y(ispec_TOL) * y(ispec_OH) &
        + rk(ireac_R61) * 0.800 * y(ispec_OH) * y(ispec_XYL) &
        + rk(ireac_R98) * 0.200 * y(ispec_OPEN) * y(ispec_O3) &
        + rk(ireac_R101) * 0.168 * y(ispec_OH) * y(ispec_ISPD) &
        + rk(ireac_R102) * 0.850 * y(ispec_O3) * y(ispec_ISPD)

    y(ispec_MGLY) = max( 0.0, ( ysum(ispec_MGLY) + gdt*yp(ispec_MGLY) ) / ( 1.0 + gdt*yl(ispec_MGLY) ) )

    !
    ! PAR
    !

    ! R50 :  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
    ! R52 :  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R54 :  NO3+OLE         ->  0.910*HO2 + 0.910*XO2 + 0.090*XO2N - PAR
    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3

    rat75 = aux(iaux_rat75)
    alpha9 = y(i_NO2)*22./(y(i_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))
    alpha10 = aux(iaux_alpha10)
    ftmp_PAR = 0.13+(0.8*(1.-alpha9)*(0.415+alpha10*1.885)) 

    yl(ispec_PAR) = 0.0  &
        + rk(ireac_R50) * y(ispec_OH)

    yp(ispec_PAR) = Q(ispec_PAR)  &
        + rk(ireac_R50) * - ftmp_par * y(ispec_PAR) * y(ispec_OH) &
        + rk(ireac_R52) * -1.0 * y(ispec_OH) * y(ispec_OLE) &
        + rk(ireac_R53) * -1.0 * y(ispec_O3) * y(ispec_OLE) &
        + rk(ireac_R54) * -1.0 * y(ispec_NO3) * y(ispec_OLE) &
        + rk(ireac_R58) * 0.360 * y(ispec_TOL) * y(ispec_OH) &
        + rk(ireac_R61) * 1.100 * y(ispec_OH) * y(ispec_XYL) &
        + rk(ireac_R81) * 0.350 * y(ispec_O3) * y(ispec_ISO) &
        + rk(ireac_R82) * 2.400 * y(ispec_NO3) * y(ispec_ISO) &
        + rk(ireac_R100) * 2.400 * y(ispec_NO2) * y(ispec_ISO) &
        + rk(ireac_R101) * 1.565 * y(ispec_OH) * y(ispec_ISPD) &
        + rk(ireac_R102) * 0.360 * y(ispec_O3) * y(ispec_ISPD) &
        + rk(ireac_R103) * 1.282 * y(ispec_NO3) * y(ispec_ISPD) &
        + rk(ireac_R104) * 0.832 * y(ispec_ISPD)

    y(ispec_PAR) = max( 0.0, ( ysum(ispec_PAR) + gdt*yp(ispec_PAR) ) / ( 1.0 + gdt*yl(ispec_PAR) ) )

    !
    ! OLE
    !

    ! R52 :  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R54 :  NO3+OLE         ->  0.910*HO2 + 0.910*XO2 + 0.090*XO2N - PAR

    yl(ispec_OLE) = 0.0  &
        + rk(ireac_R52) * y(ispec_OH) &
        + rk(ireac_R53) * y(ispec_O3) &
        + rk(ireac_R54) * y(ispec_NO3)

    yp(ispec_OLE) = Q(ispec_OLE) 

    y(ispec_OLE) = max( 0.0, ( ysum(ispec_OLE) + gdt*yp(ispec_OLE) ) / ( 1.0 + gdt*yl(ispec_OLE) ) )

    !
    ! ETH
    !

    ! R56 :  OH+ETH          ->  XO2 + 2*FORM + HO2
    ! R57 :  O3+ETH          ->  FORM + 0.370*CO + 0.130*HO2

    yl(ispec_ETH) = 0.0  &
        + rk(ireac_R56) * y(ispec_OH) &
        + rk(ireac_R57) * y(ispec_O3)

    yp(ispec_ETH) = Q(ispec_ETH) 

    y(ispec_ETH) = max( 0.0, ( ysum(ispec_ETH) + gdt*yp(ispec_ETH) ) / ( 1.0 + gdt*yl(ispec_ETH) ) )

    !
    ! TOL
    !

    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2

    yl(ispec_TOL) = 0.0  &
        + rk(ireac_R58) * y(ispec_OH)

    yp(ispec_TOL) = Q(ispec_TOL) 

    y(ispec_TOL) = max( 0.0, ( ysum(ispec_TOL) + gdt*yp(ispec_TOL) ) / ( 1.0 + gdt*yl(ispec_TOL) ) )

    !
    ! CRES
    !

    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R59 :  CRES+NO3        ->  CRO + HNO3
    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
    ! R62 :  OH+CRES         ->  0.400*CRO + 0.600*XO2N + 0.600*HO2
    ! R96 :  TO2           ->  CRES + HO2

    yl(ispec_CRES) = 0.0  &
        + rk(ireac_R59) * y(ispec_NO3) &
        + rk(ireac_R62) * y(ispec_OH)

    yp(ispec_CRES) = Q(ispec_CRES)  &
        + rk(ireac_R58) * 0.360 * y(ispec_TOL) * y(ispec_OH) &
        + rk(ireac_R61) * 0.200 * y(ispec_OH) * y(ispec_XYL) &
        + rk(ireac_R96) * y(ispec_TO2)

    y(ispec_CRES) = max( 0.0, ( ysum(ispec_CRES) + gdt*yp(ispec_CRES) ) / ( 1.0 + gdt*yl(ispec_CRES) ) )

    !
    ! XYL
    !

    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2

    yl(ispec_XYL) = 0.0  &
        + rk(ireac_R61) * y(ispec_OH)

    yp(ispec_XYL) = Q(ispec_XYL) 

    y(ispec_XYL) = max( 0.0, ( ysum(ispec_XYL) + gdt*yp(ispec_XYL) ) / ( 1.0 + gdt*yl(ispec_XYL) ) )

    !
    ! SO4a_f
    !

    ! R71 :  SO2+OH          ->  SO4a_f + HO2
    ! RH1f :  SO4a_f + N2O5    ->  SO4a_f + 2*HNO3

    yl(ispec_SO4a_f) = 0.0 

    yp(ispec_SO4a_f) = Q(ispec_SO4a_f)  &
        + rk(ireac_R71) * y(ispec_SO2) * y(ispec_OH)

    y(ispec_SO4a_f) = max( 0.0, ( ysum(ispec_SO4a_f) + gdt*yp(ispec_SO4a_f) ) / ( 1.0 + gdt*yl(ispec_SO4a_f) ) )

    !
    ! SO2
    !

    ! R71 :  SO2+OH          ->  SO4a_f + HO2

    yl(ispec_SO2) = 0.0  &
        + rk(ireac_R71) * y(ispec_OH)

    yp(ispec_SO2) = Q(ispec_SO2) 

    y(ispec_SO2) = max( 0.0, ( ysum(ispec_SO2) + gdt*yp(ispec_SO2) ) / ( 1.0 + gdt*yl(ispec_SO2) ) )

    !
    ! HNO2
    !

    ! R17 :  NO+NO2+H2O      ->  2*HNO2
    ! R18 :  HNO2+HNO2       ->  NO + NO2
    ! R19 :  HNO2            ->  NO + OH
    ! R21 :  NO+OH           ->  HNO2
    ! R26 :  OH+HNO2         ->  NO2

    yl(ispec_HNO2) = 0.0  &
        + rk(ireac_R18) * y(ispec_HNO2) &
        + rk(ireac_R18) * y(ispec_HNO2) &
        + rk(ireac_R19) &
        + rk(ireac_R26) * y(ispec_OH)

    yp(ispec_HNO2) = Q(ispec_HNO2)  &
        + rk(ireac_R17) * 2 * y(ispec_NO) * y(ispec_NO2) * aux(iaux_H2O) &
        + rk(ireac_R21) * y(ispec_NO) * y(ispec_OH)

    y(ispec_HNO2) = max( 0.0, ( ysum(ispec_HNO2) + gdt*yp(ispec_HNO2) ) / ( 1.0 + gdt*yl(ispec_HNO2) ) )

    !
    ! CH4
    !

    ! R49 :  CH4+OH          ->  XO2 + FORM + HO2

    yl(ispec_CH4) = 0.0  &
        + rk(ireac_R49) * y(ispec_OH)

    yp(ispec_CH4) = Q(ispec_CH4) 

    y(ispec_CH4) = max( 0.0, ( ysum(ispec_CH4) + gdt*yp(ispec_CH4) ) / ( 1.0 + gdt*yl(ispec_CH4) ) )

    !
    ! NO3
    !

    ! R7 :  NO2+O3          ->  NO3
    ! R12 :  NO3+NO          ->  2*NO2
    ! R13 :  NO3+NO2         ->  NO + NO2
    ! R14 :  NO3+NO2         ->  N2O5
    ! R16 :  N2O5            ->  NO3 + NO2
    ! R27 :  NO3             ->  NO2 + O3
    ! R28 :  NO3             ->  NO
    ! R36 :  FORM+NO3        ->  HNO3 + HO2 + CO
    ! R39 :  ALD+NO3         ->  C2O3 + HNO3
    ! R54 :  NO3+OLE         ->  0.910*HO2 + 0.910*XO2 + 0.090*XO2N - PAR
    ! R59 :  CRES+NO3        ->  CRO + HNO3
    ! R76 :  OH+HNO3         ->  NO3
    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3

    yl(ispec_NO3) = 0.0  &
        + rk(ireac_R12) * y(ispec_NO) &
        + rk(ireac_R13) * y(ispec_NO2) &
        + rk(ireac_R14) * y(ispec_NO2) &
        + rk(ireac_R27) &
        + rk(ireac_R28) &
        + rk(ireac_R36) * y(ispec_FORM) &
        + rk(ireac_R39) * y(ispec_ALD) &
        + rk(ireac_R54) * y(ispec_OLE) &
        + rk(ireac_R59) * y(ispec_CRES) &
        + rk(ireac_R82) * y(ispec_ISO) &
        + rk(ireac_R103) * y(ispec_ISPD)

    yp(ispec_NO3) = Q(ispec_NO3)  &
        + rk(ireac_R7) * y(ispec_NO2) * y(ispec_O3) &
        + rk(ireac_R16) * y(ispec_N2O5) &
        + rk(ireac_R76) * y(ispec_OH) * y(ispec_HNO3)

    y(ispec_NO3) = max( 0.0, ( ysum(ispec_NO3) + gdt*yp(ispec_NO3) ) / ( 1.0 + gdt*yl(ispec_NO3) ) )

    !
    ! OH
    !

    ! R8 :  O3              ->  ftmp_oh  * OH + ftmp_o3  * O3
    ! R10 :  O3+OH           ->  HO2
    ! R11 :  O3+HO2          ->  OH
    ! R19 :  HNO2            ->  NO + OH
    ! R20 :  NO2+OH          ->  HNO3
    ! R21 :  NO+OH           ->  HNO2
    ! R22 :  HO2+NO          ->  OH + NO2
    ! R26 :  OH+HNO2         ->  NO2
    ! R31 :  OH+CO           ->  HO2 + CO2
    ! R32 :  FORM+OH         ->  HO2 + CO
    ! R38 :  ALD+OH          ->  C2O3
    ! R46 :  C2O3+HO2        ->  0.79 * FORM + 0.79 * HO2 + 0.79 * XO2 + 0.79 *OH
    ! R48 :  OH+MGLY         ->  XO2 + C2O3
    ! R49 :  CH4+OH          ->  XO2 + FORM + HO2
    ! R50 :  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
    ! R52 :  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R56 :  OH+ETH          ->  XO2 + 2*FORM + HO2
    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
    ! R62 :  OH+CRES         ->  0.400*CRO + 0.600*XO2N + 0.600*HO2
    ! R71 :  SO2+OH          ->  SO4a_f + HO2
    ! R73 :  OH+H2O2         ->  HO2
    ! R74 :  H2O2            ->  2*OH
    ! R75 :  HNO3            ->  OH + NO2
    ! R76 :  OH+HNO3         ->  NO3
    ! R80 :  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R97 :  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD

    alpha8 = aux(iaux_alpha8)
    ftmp_OH = 2.*alpha8 

    yl(ispec_OH) = 0.0  &
        + rk(ireac_R10) * y(ispec_O3) &
        + rk(ireac_R20) * y(ispec_NO2) &
        + rk(ireac_R21) * y(ispec_NO) &
        + rk(ireac_R26) * y(ispec_HNO2) &
        + rk(ireac_R31) * y(ispec_CO) &
        + rk(ireac_R32) * y(ispec_FORM) &
        + rk(ireac_R38) * y(ispec_ALD) &
        + rk(ireac_R48) * y(ispec_MGLY) &
        + rk(ireac_R49) * y(ispec_CH4) &
        + rk(ireac_R50) * y(ispec_PAR) &
        + rk(ireac_R52) * y(ispec_OLE) &
        + rk(ireac_R56) * y(ispec_ETH) &
        + rk(ireac_R58) * y(ispec_TOL) &
        + rk(ireac_R61) * y(ispec_XYL) &
        + rk(ireac_R62) * y(ispec_CRES) &
        + rk(ireac_R71) * y(ispec_SO2) &
        + rk(ireac_R73) * y(ispec_H2O2) &
        + rk(ireac_R76) * y(ispec_HNO3) &
        + rk(ireac_R80) * y(ispec_ISO) &
        + rk(ireac_R97) * y(ispec_OPEN) &
        + rk(ireac_R101) * y(ispec_ISPD)

    yp(ispec_OH) = Q(ispec_OH)  &
        + rk(ireac_R8) * ftmp_oh * y(ispec_O3) &
        + rk(ireac_R11) * y(ispec_O3) * y(ispec_HO2) &
        + rk(ireac_R19) * y(ispec_HNO2) &
        + rk(ireac_R22) * y(ispec_HO2) * y(ispec_NO) &
        + rk(ireac_R46) * 0.79 * y(ispec_C2O3) * y(ispec_HO2) &
        + rk(ireac_R53) * 0.080 * y(ispec_O3) * y(ispec_OLE) &
        + rk(ireac_R74) * 2 * y(ispec_H2O2) &
        + rk(ireac_R75) * y(ispec_HNO3) &
        + rk(ireac_R81) * 0.266 * y(ispec_O3) * y(ispec_ISO) &
        + rk(ireac_R98) * 0.080 * y(ispec_OPEN) * y(ispec_O3) &
        + rk(ireac_R102) * 0.266 * y(ispec_O3) * y(ispec_ISPD)

    y(ispec_OH) = max( 0.0, ( ysum(ispec_OH) + gdt*yp(ispec_OH) ) / ( 1.0 + gdt*yl(ispec_OH) ) )

    !
    ! HO2
    !

    ! R10 :  O3+OH           ->  HO2
    ! R11 :  O3+HO2          ->  OH
    ! R22 :  HO2+NO          ->  OH + NO2
    ! R29 :  HO2+HO2         ->  H2O2
    ! R30 :  HO2+HO2+H2O     ->  H2O2
    ! R31 :  OH+CO           ->  HO2 + CO2
    ! R32 :  FORM+OH         ->  HO2 + CO
    ! R33 :  FORM            ->  2*HO2 + CO
    ! R36 :  FORM+NO3        ->  HNO3 + HO2 + CO
    ! R40 :  ALD             ->  CO + FORM + 2*HO2 + XO2
    ! R42 :  C2O3+NO         ->  NO2 + XO2 + FORM + HO2
    ! R45 :  C2O3+C2O3       ->  2*FORM + 2*XO2 + 2*HO2
    ! R46 :  C2O3+HO2        ->  0.79 * FORM + 0.79 * HO2 + 0.79 * XO2 + 0.79 *OH
    ! R47 :  MGLY            ->  C2O3 + HO2 + CO
    ! R49 :  CH4+OH          ->  XO2 + FORM + HO2
    ! R50 :  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
    ! R52 :  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R54 :  NO3+OLE         ->  0.910*HO2 + 0.910*XO2 + 0.090*XO2N - PAR
    ! R56 :  OH+ETH          ->  XO2 + 2*FORM + HO2
    ! R57 :  O3+ETH          ->  FORM + 0.370*CO + 0.130*HO2
    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
    ! R62 :  OH+CRES         ->  0.400*CRO + 0.600*XO2N + 0.600*HO2
    ! R66 :  XO2+HO2         ->  
    ! R67 :  XO2N+HO2        ->  
    ! R71 :  SO2+OH          ->  SO4a_f + HO2
    ! R73 :  OH+H2O2         ->  HO2
    ! R80 :  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! R95 :  TO2+NO        ->  NO2 + 0.900*HO2 + 0.900*OPEN
    ! R96 :  TO2           ->  CRES + HO2
    ! R97 :  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R99 :  OPEN          ->  C2O3 + HO2 + CO
    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3

    alpha9 = y(i_NO2)*22./(y(i_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))
    ftmp_HO2 = 0.13 + (0.8*(1.-alpha9)) 

    yl(ispec_HO2) = 0.0  &
        + rk(ireac_R11) * y(ispec_O3) &
        + rk(ireac_R22) * y(ispec_NO) &
        + rk(ireac_R29) * y(ispec_HO2) &
        + rk(ireac_R29) * y(ispec_HO2) &
        + rk(ireac_R30) * y(ispec_HO2) * aux(iaux_H2O) &
        + rk(ireac_R30) * y(ispec_HO2) * aux(iaux_H2O) &
        + rk(ireac_R46) * y(ispec_C2O3) &
        + rk(ireac_R66) * y(ispec_XO2) &
        + rk(ireac_R67) * y(ispec_XO2N)

    yp(ispec_HO2) = Q(ispec_HO2)  &
        + rk(ireac_R10) * y(ispec_O3) * y(ispec_OH) &
        + rk(ireac_R31) * y(ispec_OH) * y(ispec_CO) &
        + rk(ireac_R32) * y(ispec_FORM) * y(ispec_OH) &
        + rk(ireac_R33) * 2 * y(ispec_FORM) &
        + rk(ireac_R36) * y(ispec_FORM) * y(ispec_NO3) &
        + rk(ireac_R40) * 2 * y(ispec_ALD) &
        + rk(ireac_R42) * y(ispec_C2O3) * y(ispec_NO) &
        + rk(ireac_R45) * 2 * y(ispec_C2O3) * y(ispec_C2O3) &
        + rk(ireac_R46) * 0.79 * y(ispec_C2O3) * y(ispec_HO2) &
        + rk(ireac_R47) * y(ispec_MGLY) &
        + rk(ireac_R49) * y(ispec_CH4) * y(ispec_OH) &
        + rk(ireac_R50) * ftmp_ho2 * y(ispec_PAR) * y(ispec_OH) &
        + rk(ireac_R52) * y(ispec_OH) * y(ispec_OLE) &
        + rk(ireac_R53) * 0.280 * y(ispec_O3) * y(ispec_OLE) &
        + rk(ireac_R54) * 0.910 * y(ispec_NO3) * y(ispec_OLE) &
        + rk(ireac_R56) * y(ispec_OH) * y(ispec_ETH) &
        + rk(ireac_R57) * 0.130 * y(ispec_O3) * y(ispec_ETH) &
        + rk(ireac_R58) * y(ispec_TOL) * y(ispec_OH) &
        + rk(ireac_R61) * 0.700 * y(ispec_OH) * y(ispec_XYL) &
        + rk(ireac_R62) * 0.600 * y(ispec_OH) * y(ispec_CRES) &
        + rk(ireac_R71) * y(ispec_SO2) * y(ispec_OH) &
        + rk(ireac_R73) * y(ispec_OH) * y(ispec_H2O2) &
        + rk(ireac_R80) * 0.912 * y(ispec_OH) * y(ispec_ISO) &
        + rk(ireac_R81) * 0.066 * y(ispec_O3) * y(ispec_ISO) &
        + rk(ireac_R82) * 0.800 * y(ispec_NO3) * y(ispec_ISO) &
        + rk(ireac_R95) * 0.900 * y(ispec_TO2) * y(ispec_NO) &
        + rk(ireac_R96) * y(ispec_TO2) &
        + rk(ireac_R97) * 2 * y(ispec_OPEN) * y(ispec_OH) &
        + rk(ireac_R98) * 0.760 * y(ispec_OPEN) * y(ispec_O3) &
        + rk(ireac_R99) * y(ispec_OPEN) &
        + rk(ireac_R100) * 0.800 * y(ispec_NO2) * y(ispec_ISO) &
        + rk(ireac_R101) * 0.500 * y(ispec_OH) * y(ispec_ISPD) &
        + rk(ireac_R102) * 0.154 * y(ispec_O3) * y(ispec_ISPD) &
        + rk(ireac_R103) * 0.925 * y(ispec_NO3) * y(ispec_ISPD) &
        + rk(ireac_R104) * 1.03 * y(ispec_ISPD)

    y(ispec_HO2) = max( 0.0, ( ysum(ispec_HO2) + gdt*yp(ispec_HO2) ) / ( 1.0 + gdt*yl(ispec_HO2) ) )

    !
    ! N2O5
    !

    ! R14 :  NO3+NO2         ->  N2O5
    ! R16 :  N2O5            ->  NO3 + NO2
    ! RH1f :  SO4a_f + N2O5    ->  SO4a_f + 2*HNO3
    ! RH2f2 :  Na_f2 + N2O5      ->  Na_f2 + 2*NO3a_f
    ! RH2f1 :  Na_f1 + N2O5      ->  Na_f1 + 2*NO3a_f
    ! RH2c3 :  Na_c3 + N2O5      ->  Na_c3 + 2*NO3a_c
    ! RH2c2 :  Na_c2 + N2O5      ->  Na_c2 + 2*NO3a_c
    ! RH2c1 :  Na_c1 + N2O5      ->  Na_c1 + 2*NO3a_c

    yl(ispec_N2O5) = 0.0  &
        + rk(ireac_R16) &
        + rk(ireac_RH1f) * y(ispec_SO4a_f) &
        + rk(ireac_RH2f2) * y(ispec_Na_f2) &
        + rk(ireac_RH2f1) * y(ispec_Na_f1) &
        + rk(ireac_RH2c3) * y(ispec_Na_c3) &
        + rk(ireac_RH2c2) * y(ispec_Na_c2) &
        + rk(ireac_RH2c1) * y(ispec_Na_c1)

    yp(ispec_N2O5) = Q(ispec_N2O5)  &
        + rk(ireac_R14) * y(ispec_NO3) * y(ispec_NO2)

    y(ispec_N2O5) = max( 0.0, ( ysum(ispec_N2O5) + gdt*yp(ispec_N2O5) ) / ( 1.0 + gdt*yl(ispec_N2O5) ) )

    !
    ! NO3a_f
    !

    ! RH2f2 :  Na_f1 + N2O5      ->  Na_f1 + 2*NO3a_f
    ! RH2f1 :  Na_f2 + N2O5      ->  Na_f2 + 2*NO3a_f
    ! RH3f1 :   Na_f2 + HNO3     ->  Na_f2 + NO3a_f
    ! RH3f2 :    Na_f1 + HNO3      ->  Na_f1 + NO3a_f

    yl(ispec_NO3a_f) = 0.0 

    yp(ispec_NO3a_f) = Q(ispec_NO3a_f)  &
        + rk(ireac_RH2f2) * 2 * y(ispec_Na_f1) * y(ispec_N2O5) &
        + rk(ireac_RH2f1) * 2 * y(ispec_Na_f2) * y(ispec_N2O5) &
        + rk(ireac_RH3f1) * y(ispec_Na_f2) * y(ispec_HNO3) &
        + rk(ireac_RH3f2) * y(ispec_Na_f1) * y(ispec_HNO3)

    y(ispec_NO3a_f) = max( 0.0, ( ysum(ispec_NO3a_f) + gdt*yp(ispec_NO3a_f) ) / ( 1.0 + gdt*yl(ispec_NO3a_f) ) )

    !
    ! NO3a_c
    !

    ! RH2c3 :   Na_c3 + N2O5     ->  Na_c3 + 2*NO3a_c
    ! RH2c2 :   Na_c2 + N2O5     ->  Na_c2 + 2*NO3a_c
    ! RH2c1 :   Na_c1 + N2O5     ->  Na_c1 + 2*NO3a_c
    ! RH3c3 :   Na_c3 + HNO3     ->  Na_c3 + NO3a_c
    ! RH3c2 :   Na_c2 + HNO3     ->  Na_c2 + NO3a_c
    ! RH3c1 :   Na_c1 + HNO3     ->  Na_c1 + NO3a_c

    yl(ispec_NO3a_c) = 0.0 

    yp(ispec_NO3a_c) = Q(ispec_NO3a_c)  &
        + rk(ireac_RH2c3) * 2 * y(ispec_Na_c3) * y(ispec_N2O5) &
        + rk(ireac_RH2c2) * 2 * y(ispec_Na_c2) * y(ispec_N2O5) &
        + rk(ireac_RH2c1) * 2 * y(ispec_Na_c1) * y(ispec_N2O5) &
        + rk(ireac_RH3c3) * y(ispec_Na_c3) * y(ispec_HNO3) &
        + rk(ireac_RH3c2) * y(ispec_Na_c2) * y(ispec_HNO3) &
        + rk(ireac_RH3c1) * y(ispec_Na_c1) * y(ispec_HNO3)

    y(ispec_NO3a_c) = max( 0.0, ( ysum(ispec_NO3a_c) + gdt*yp(ispec_NO3a_c) ) / ( 1.0 + gdt*yl(ispec_NO3a_c) ) )

    !
    ! C2O3
    !

    ! R38 :  ALD+OH          ->  C2O3
    ! R39 :  ALD+NO3         ->  C2O3 + HNO3
    ! R42 :  C2O3+NO         ->  NO2 + XO2 + FORM + HO2
    ! R43 :  C2O3+NO2        ->  PAN
    ! R44 :  PAN             ->  C2O3 + NO2
    ! R45 :  C2O3+C2O3       ->  2*FORM + 2*XO2 + 2*HO2
    ! R46 :  C2O3+HO2        ->  0.79 * FORM + 0.79 * HO2 + 0.79 * XO2 + 0.79 *OH
    ! R47 :  MGLY            ->  C2O3 + HO2 + CO
    ! R48 :  OH+MGLY         ->  XO2 + C2O3
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R97 :  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R99 :  OPEN          ->  C2O3 + HO2 + CO
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3

    yl(ispec_C2O3) = 0.0  &
        + rk(ireac_R42) * y(ispec_NO) &
        + rk(ireac_R43) * y(ispec_NO2) &
        + rk(ireac_R45) * y(ispec_C2O3) &
        + rk(ireac_R45) * y(ispec_C2O3) &
        + rk(ireac_R46) * y(ispec_HO2)

    yp(ispec_C2O3) = Q(ispec_C2O3)  &
        + rk(ireac_R38) * y(ispec_ALD) * y(ispec_OH) &
        + rk(ireac_R39) * y(ispec_ALD) * y(ispec_NO3) &
        + rk(ireac_R44) * y(ispec_PAN) &
        + rk(ireac_R47) * y(ispec_MGLY) &
        + rk(ireac_R48) * y(ispec_OH) * y(ispec_MGLY) &
        + rk(ireac_R81) * 0.200 * y(ispec_O3) * y(ispec_ISO) &
        + rk(ireac_R97) * y(ispec_OPEN) * y(ispec_OH) &
        + rk(ireac_R98) * 0.620 * y(ispec_OPEN) * y(ispec_O3) &
        + rk(ireac_R99) * y(ispec_OPEN) &
        + rk(ireac_R101) * 0.498 * y(ispec_OH) * y(ispec_ISPD) &
        + rk(ireac_R102) * 0.114 * y(ispec_O3) * y(ispec_ISPD) &
        + rk(ireac_R103) * 0.075 * y(ispec_NO3) * y(ispec_ISPD) &
        + rk(ireac_R104) * 1.667 * y(ispec_ISPD)

    y(ispec_C2O3) = max( 0.0, ( ysum(ispec_C2O3) + gdt*yp(ispec_C2O3) ) / ( 1.0 + gdt*yl(ispec_C2O3) ) )

    !
    ! XO2
    !

    ! R40 :  ALD             ->  CO + FORM + 2*HO2 + XO2
    ! R42 :  C2O3+NO         ->  NO2 + XO2 + FORM + HO2
    ! R45 :  C2O3+C2O3       ->  2*FORM + 2*XO2 + 2*HO2
    ! R46 :  C2O3+HO2        ->  0.79 * FORM + 0.79 * HO2 + 0.79 * XO2 + 0.79 *OH
    ! R48 :  OH+MGLY         ->  XO2 + C2O3
    ! R49 :  CH4+OH          ->  XO2 + FORM + HO2
    ! R50 :  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
    ! R52 :  OH+OLE          ->  FORM + ALD + XO2 + HO2 - PAR
    ! R53 :  O3+OLE          ->  0.500*ALD + 0.660*FORM  + 0.212*CO + 0.280*HO2 + 0.080*OH + 0.144*XO2 - PAR
    ! R54 :  NO3+OLE         ->  0.910*HO2 + 0.910*XO2 + 0.090*XO2N - PAR
    ! R56 :  OH+ETH          ->  XO2 + 2*FORM + HO2
    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
    ! R63 :  XO2+NO          ->  NO2
    ! R65 :  XO2+XO2         ->  
    ! R66 :  XO2+HO2         ->  
    ! R69 :  XO2+XO2N        ->  
    ! R80 :  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N
    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! R97 :  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3

    rat75 = aux(iaux_rat75)
    alpha9 = y(i_NO2)*22./(y(i_NO2)*22.+390000.*aux(iaux_ppb_to_mlccm3)+rat75*aux(iaux_ppb_to_mlccm3))
    alpha10 = aux(iaux_alpha10)
    ftmp_XO2 = 0.93 + (0.8*(1.-alpha9)*(0.615+alpha10*0.385)) 

    yl(ispec_XO2) = 0.0  &
        + rk(ireac_R63) * y(ispec_NO) &
        + rk(ireac_R65) * y(ispec_XO2) &
        + rk(ireac_R65) * y(ispec_XO2) &
        + rk(ireac_R66) * y(ispec_HO2) &
        + rk(ireac_R69) * y(ispec_XO2N)

    yp(ispec_XO2) = Q(ispec_XO2)  &
        + rk(ireac_R40) * y(ispec_ALD) &
        + rk(ireac_R42) * y(ispec_C2O3) * y(ispec_NO) &
        + rk(ireac_R45) * 2 * y(ispec_C2O3) * y(ispec_C2O3) &
        + rk(ireac_R46) * 0.79 * y(ispec_C2O3) * y(ispec_HO2) &
        + rk(ireac_R48) * y(ispec_OH) * y(ispec_MGLY) &
        + rk(ireac_R49) * y(ispec_CH4) * y(ispec_OH) &
        + rk(ireac_R50) * ftmp_xo2 * y(ispec_PAR) * y(ispec_OH) &
        + rk(ireac_R52) * y(ispec_OH) * y(ispec_OLE) &
        + rk(ireac_R53) * 0.144 * y(ispec_O3) * y(ispec_OLE) &
        + rk(ireac_R54) * 0.910 * y(ispec_NO3) * y(ispec_OLE) &
        + rk(ireac_R56) * y(ispec_OH) * y(ispec_ETH) &
        + rk(ireac_R58) * 0.640 * y(ispec_TOL) * y(ispec_OH) &
        + rk(ireac_R61) * 0.100 * y(ispec_OH) * y(ispec_XYL) &
        + rk(ireac_R80) * 0.991 * y(ispec_OH) * y(ispec_ISO) &
        + rk(ireac_R82) * y(ispec_NO3) * y(ispec_ISO) &
        + rk(ireac_R97) * y(ispec_OPEN) * y(ispec_OH) &
        + rk(ireac_R98) * 0.030 * y(ispec_OPEN) * y(ispec_O3) &
        + rk(ireac_R100) * y(ispec_NO2) * y(ispec_ISO) &
        + rk(ireac_R101) * 0.713 * y(ispec_OH) * y(ispec_ISPD) &
        + rk(ireac_R102) * 0.064 * y(ispec_O3) * y(ispec_ISPD) &
        + rk(ireac_R103) * 0.074 * y(ispec_NO3) * y(ispec_ISPD) &
        + rk(ireac_R104) * 0.700 * y(ispec_ISPD)

    y(ispec_XO2) = max( 0.0, ( ysum(ispec_XO2) + gdt*yp(ispec_XO2) ) / ( 1.0 + gdt*yl(ispec_XO2) ) )

    !
    ! XO2N
    !

    ! R50 :  PAR+OH          ->  ftmp_xo2*XO2 + 0.067*XO2N + ftmp_ho2 *HO2 + ftmp_ald*ALD - ftmp_no2*NO2 - ftmp_par*PAR
    ! R54 :  NO3+OLE         ->  0.910*HO2 + 0.910*XO2 + 0.090*XO2N - PAR
    ! R62 :  OH+CRES         ->  0.400*CRO + 0.600*XO2N + 0.600*HO2
    ! R64 :  XO2N+NO         ->  
    ! R67 :  XO2N+HO2        ->  
    ! R68 :  XO2N+XO2N       ->  
    ! R69 :  XO2+XO2N        ->  
    ! R80 :  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N

    yl(ispec_XO2N) = 0.0  &
        + rk(ireac_R64) * y(ispec_NO) &
        + rk(ireac_R67) * y(ispec_HO2) &
        + rk(ireac_R68) * y(ispec_XO2N) &
        + rk(ireac_R68) * y(ispec_XO2N) &
        + rk(ireac_R69) * y(ispec_XO2)

    yp(ispec_XO2N) = Q(ispec_XO2N)  &
        + rk(ireac_R50) * 0.067 * y(ispec_PAR) * y(ispec_OH) &
        + rk(ireac_R54) * 0.090 * y(ispec_NO3) * y(ispec_OLE) &
        + rk(ireac_R62) * 0.600 * y(ispec_OH) * y(ispec_CRES) &
        + rk(ireac_R80) * 0.088 * y(ispec_OH) * y(ispec_ISO)

    y(ispec_XO2N) = max( 0.0, ( ysum(ispec_XO2N) + gdt*yp(ispec_XO2N) ) / ( 1.0 + gdt*yl(ispec_XO2N) ) )

    !
    ! CRO
    !

    ! R59 :  CRES+NO3        ->  CRO + HNO3
    ! R60 :  CRO+NO2         ->  
    ! R62 :  OH+CRES         ->  0.400*CRO + 0.600*XO2N + 0.600*HO2

    yl(ispec_CRO) = 0.0  &
        + rk(ireac_R60) * y(ispec_NO2)

    yp(ispec_CRO) = Q(ispec_CRO)  &
        + rk(ireac_R59) * y(ispec_CRES) * y(ispec_NO3) &
        + rk(ireac_R62) * 0.400 * y(ispec_OH) * y(ispec_CRES)

    y(ispec_CRO) = max( 0.0, ( ysum(ispec_CRO) + gdt*yp(ispec_CRO) ) / ( 1.0 + gdt*yl(ispec_CRO) ) )

    !
    ! HNO3
    !

    ! R20 :  NO2+OH          ->  HNO3
    ! R36 :  FORM+NO3        ->  HNO3 + HO2 + CO
    ! R39 :  ALD+NO3         ->  C2O3 + HNO3
    ! R59 :  CRES+NO3        ->  CRO + HNO3
    ! R75 :  HNO3            ->  OH + NO2
    ! R76 :  OH+HNO3         ->  NO3
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! RH1f :  SO4a_f + N2O5    ->  SO4a_f + 2*HNO3
    ! RH3f1 :   Na_f1 + HNO3     ->  Na_f1 + NO3a_f
    ! RH3f2 :    Na_f2 + HNO3      ->  Na_f2 + NO3a_f
    ! RH3c3 :    Na_c3 + HNO3      ->  Na_c3 + NO3a_c
    ! RH3c2 :   Na_c2 + HNO3     ->  Na_c2 + NO3a_c
    ! RH3c1 :  Na_c1 + HNO3    ->  Na_c1 + NO3a_c

    yl(ispec_HNO3) = 0.0  &
        + rk(ireac_R75) &
        + rk(ireac_R76) * y(ispec_OH) &
        + rk(ireac_RH3f1) * y(ispec_Na_f1) &
        + rk(ireac_RH3f2) * y(ispec_Na_f2) &
        + rk(ireac_RH3c3) * y(ispec_Na_c3) &
        + rk(ireac_RH3c2) * y(ispec_Na_c2) &
        + rk(ireac_RH3c1) * y(ispec_Na_c1)

    yp(ispec_HNO3) = Q(ispec_HNO3)  &
        + rk(ireac_R20) * y(ispec_NO2) * y(ispec_OH) &
        + rk(ireac_R36) * y(ispec_FORM) * y(ispec_NO3) &
        + rk(ireac_R39) * y(ispec_ALD) * y(ispec_NO3) &
        + rk(ireac_R59) * y(ispec_CRES) * y(ispec_NO3) &
        + rk(ireac_R103) * 0.075 * y(ispec_NO3) * y(ispec_ISPD) &
        + rk(ireac_RH1f) * 2 * y(ispec_SO4a_f) * y(ispec_N2O5)

    y(ispec_HNO3) = max( 0.0, ( ysum(ispec_HNO3) + gdt*yp(ispec_HNO3) ) / ( 1.0 + gdt*yl(ispec_HNO3) ) )

    !
    ! H2O2
    !

    ! R29 :  HO2+HO2         ->  H2O2
    ! R30 :  HO2+HO2+H2O     ->  H2O2
    ! R73 :  OH+H2O2         ->  HO2
    ! R74 :  H2O2            ->  2*OH

    yl(ispec_H2O2) = 0.0  &
        + rk(ireac_R73) * y(ispec_OH) &
        + rk(ireac_R74)

    yp(ispec_H2O2) = Q(ispec_H2O2)  &
        + rk(ireac_R29) * y(ispec_HO2) * y(ispec_HO2) &
        + rk(ireac_R30) * y(ispec_HO2) * y(ispec_HO2) * aux(iaux_H2O)

    y(ispec_H2O2) = max( 0.0, ( ysum(ispec_H2O2) + gdt*yp(ispec_H2O2) ) / ( 1.0 + gdt*yl(ispec_H2O2) ) )

    !
    ! ISO
    !

    ! R80 :  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR

    yl(ispec_ISO) = 0.0  &
        + rk(ireac_R80) * y(ispec_OH) &
        + rk(ireac_R81) * y(ispec_O3) &
        + rk(ireac_R82) * y(ispec_NO3) &
        + rk(ireac_R100) * y(ispec_NO2)

    yp(ispec_ISO) = Q(ispec_ISO) 

    y(ispec_ISO) = max( 0.0, ( ysum(ispec_ISO) + gdt*yp(ispec_ISO) ) / ( 1.0 + gdt*yl(ispec_ISO) ) )

    !
    ! ISPD
    !

    ! R80 :  OH+ISO          ->  ISPD + 0.629*FORM + 0.991*XO2 + 0.912*HO2 + 0.088*XO2N
    ! R81 :  O3+ISO          ->  ISPD + 0.600*FORM + 0.066*HO2 + 0.266*OH + 0.200*C2O3 + 0.150*ALD + 0.350*PAR + 0.066*CO
    ! R82 :  NO3+ISO         ->  ISPD + XO2 + 0.800*HO2 + NO2 + 0.800*ALD + 2.400*PAR
    ! R100 :  NO2+ISO       ->  ISPD + XO2 + 0.800*HO2 + 0.200*NO + 0.800*ALD  + 2.400*PAR
    ! R101 :  OH+ISPD       ->  1.565*PAR + 0.167*FORM + 0.713*XO2 + 0.500*HO2 + 0.334*CO + 0.168*MGLY + 0.498*C2O3 + 0.273*ALD
    ! R102 :  O3+ISPD       ->  0.114*C2O3 + 0.150*FORM + 0.850*MGLY + 0.154*HO2 + 0.266*OH + 0.064*XO2 + 0.360*PAR + 0.225*CO + 0.020*ALD
    ! R103 :  NO3+ISPD      ->  0.357*ALD + 0.282*FORM + 1.282*PAR + 0.925*HO2 + 0.643*CO + 0.075*C2O3 + 0.074*XO2 + 0.075*HNO3
    ! R104 :  ISPD          ->  0.33*CO + 0.067*ALD + 0.900*FORM + 0.832*PAR + 1.03*HO2 + 0.700*XO2 + 1.667*C2O3

    yl(ispec_ISPD) = 0.0  &
        + rk(ireac_R101) * y(ispec_OH) &
        + rk(ireac_R102) * y(ispec_O3) &
        + rk(ireac_R103) * y(ispec_NO3) &
        + rk(ireac_R104)

    yp(ispec_ISPD) = Q(ispec_ISPD)  &
        + rk(ireac_R80) * y(ispec_OH) * y(ispec_ISO) &
        + rk(ireac_R81) * y(ispec_O3) * y(ispec_ISO) &
        + rk(ireac_R82) * y(ispec_NO3) * y(ispec_ISO) &
        + rk(ireac_R100) * y(ispec_NO2) * y(ispec_ISO)

    y(ispec_ISPD) = max( 0.0, ( ysum(ispec_ISPD) + gdt*yp(ispec_ISPD) ) / ( 1.0 + gdt*yl(ispec_ISPD) ) )

    !
    ! TO2
    !

    ! R58 :  TOL+OH          ->  1.130*CO +  1.130*FORM + 0.560*MGLY + 0.360*PAR + 0.360*CRES + HO2 + 0.640*XO2 + 0.560*TO2
    ! R61 :  OH+XYL          ->  0.700*HO2 + 1.100*PAR + 0.800*MGLY +  0.200*CRES + 0.300*TO2 + 0.100*XO2
    ! R95 :  TO2+NO        ->  NO2 + 0.900*HO2 + 0.900*OPEN
    ! R96 :  TO2           ->  CRES + HO2

    yl(ispec_TO2) = 0.0  &
        + rk(ireac_R95) * y(ispec_NO) &
        + rk(ireac_R96)

    yp(ispec_TO2) = Q(ispec_TO2)  &
        + rk(ireac_R58) * 0.560 * y(ispec_TOL) * y(ispec_OH) &
        + rk(ireac_R61) * 0.300 * y(ispec_OH) * y(ispec_XYL)

    y(ispec_TO2) = max( 0.0, ( ysum(ispec_TO2) + gdt*yp(ispec_TO2) ) / ( 1.0 + gdt*yl(ispec_TO2) ) )

    !
    ! OPEN
    !

    ! R95 :  TO2+NO        ->  NO2 + 0.900*HO2 + 0.900*OPEN
    ! R97 :  OPEN+OH       ->  XO2 + C2O3 + 2*HO2 + 2*CO + FORM
    ! R98 :  OPEN+O3       ->  0.030*ALD + 0.620*C2O3 + 0.700*FORM + 0.030*XO2 + 0.690*CO + 0.080*OH + 0.760*HO2 + 0.200*MGLY
    ! R99 :  OPEN          ->  C2O3 + HO2 + CO

    yl(ispec_OPEN) = 0.0  &
        + rk(ireac_R97) * y(ispec_OH) &
        + rk(ireac_R98) * y(ispec_O3) &
        + rk(ireac_R99)

    yp(ispec_OPEN) = Q(ispec_OPEN)  &
        + rk(ireac_R95) * 0.900 * y(ispec_TO2) * y(ispec_NO)

    y(ispec_OPEN) = max( 0.0, ( ysum(ispec_OPEN) + gdt*yp(ispec_OPEN) ) / ( 1.0 + gdt*yl(ispec_OPEN) ) )

    !
    ! Na_f1
    !

    ! RH2f1 :  Na_f1 + N2O5      ->  Na_f1 + 2*NO3a_f
    ! RH3f1 :  Na_f1 + HNO3      ->  Na_f1 + NO3a_f

    !
    ! Na_f2
    !

    ! RH2f2 :  Na_f2 + N2O5      ->  Na_f2 + 2*NO3a_f
    ! RH3f2 :  Na_f2 + HNO3      ->  Na_f2 + NO3a_f

    !
    ! Na_c3
    !

    ! RH2c3 :   Na_c3 + N2O5      ->  Na_c3 + 2*NO3a_c
    ! RH3c3 :   Na_c3 + HNO3      ->  Na_c3 + NO3a_c

    !
    ! Na_c2
    !

    ! RH2c2 :  Na_c2 + N2O5      ->  Na_c2 + 2*NO3a_c
    ! RH3c2 :  Na_c2 + HNO3     ->  Na_c2 + NO3a_c

    !
    ! Na_c1
    !

    ! RH2c1 : Na_c1 + N2O5      ->  Na_c1 + 2*NO3a_c
    ! RH3c1 :  Na_c1 + HNO3    ->  Na_c1 + NO3a_c

  end subroutine LE_Chem_Work_Iter


end module LE_Chem_Work
