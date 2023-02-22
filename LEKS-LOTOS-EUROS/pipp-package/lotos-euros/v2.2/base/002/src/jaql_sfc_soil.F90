!#######################################################################
!
! JAQL_SFC_Dust
!   Joint Air Quality Library
!   Surface modelling
!   Soil parameterizations
!
!
!#######################################################################


module JAQL_SFC_Soil

  implicit none
  

  ! --- in/out --------------------------------
  
  private
  
  public  ::  nssd, nmodes
  
  public  ::  issd_S, issd_LS, issd_SL, issd_SiL, issd_Si,issd_L, issd_SCL, &
      issd_SiCL, issd_CL, issd_SC, issd_SiC, issd_C
  ! S = sand, Si = silt, C = clay
  
  public  ::  ssd_gmean, ssd_gsigma, ssd_mfrac
  public  ::  rho_mineral_soil
  

  ! --- const --------------------------------
  
  ! * soil size modes

  ! number soil size modes (modes after Mokhtari et al., 2012)
  integer, parameter :: nssd = 12
  integer, parameter :: nmodes = 3

  integer, parameter  ::  issd_S    = 1   ! Sand
  integer, parameter  ::  issd_LS   = 2   ! Loamy sand
  integer, parameter  ::  issd_SL   = 3   ! Sandy loam
  integer, parameter  ::  issd_SiL  = 4   ! Silt loam
  integer, parameter  ::  issd_Si   = 5   ! Silt
  integer, parameter  ::  issd_L    = 6   ! Loam
  integer, parameter  ::  issd_SCL  = 7   ! Sandy clay loam
  integer, parameter  ::  issd_SiCL = 8   ! Silty clay loam
  integer, parameter  ::  issd_CL   = 9   ! Clay loam
  integer, parameter  ::  issd_SC   = 10  ! Sandy clay
  integer, parameter  ::  issd_SiC  = 11  ! Silty clay
  integer, parameter  ::  issd_C    = 12  ! Clay

  real, parameter    :: ssd_gmean(nssd, nmodes) = reshape( &
    (/  1000.0e-6,  690.0e-6, 520.0e-6, 520.0e-6, 520.0e-6, 520.0e-6, 210.0e-6, 210.0e-6, 125.0e-6, 100.0e-6, 100.0e-6, 100.0e-6,   &
      100.0e-6, 100.0e-6, 100.0e-6, 100.0e-6, 75.0e-6,  75.0e-6,  75.0e-6,  50.0e-6,  50.0e-6,  10.0e-6,  10.0e-6,  10.0e-6,    &
      10.0e-6,  10.0e-6,  5.0e-6,   5.0e-6,   2.5e-6,   2.5e-6,   2.5e-6,   2.5e-6,   1.0e-6,   1.0e-6,   0.5e-6,   0.5e-6    /), &
    (/nssd,nmodes/) )! geometric mean mass diameter [m]
  real, parameter    :: ssd_gsigma(nssd, nmodes) = reshape( &
    (/  1.6,    1.6,    1.6,    1.6,    1.6,    1.6,    1.7,    1.7,    1.7,    1.8,    1.8,    1.8,  &
      1.7,    1.7,    1.7,    1.7,    1.7,    1.7,    1.7,    1.7,    1.7,    1.8,    1.8,    1.8,  &
      1.8,    1.8,    1.8,    1.8,    1.8,    1.8,    1.8,    1.8,    1.8,    1.8,    1.8,    1.8   /), &
    (/nssd,nmodes/) )! geometric standard deviation [-]
  real, parameter    :: ssd_mfrac(nssd, nmodes) = reshape( &
    (/  0.9,    0.6,    0.6,    0.5,    0.45,   0.35,   0.3,    0.3,    0.2,    0.65,   0.6,    0.5,  &
        0.1,    0.3,    0.3,    0.35,   0.4,    0.5,    0.5,    0.5,    0.5,    0.0,    0.0,    0.0,  &
        0.0,    0.1,    0.1,    0.15,   0.15,   0.15,   0.2,    0.2,    0.3,    0.35,   0.4,    0.5   /), &
    (/nssd,nmodes/) )! geometric standard deviation [-]
                                                
  ! * particle density
  
  ! mass density of minerals present in arid soils   (Tegen et al., 2002) :
  real, parameter :: rho_mineral_soil = 2.65e3    ! kg/m3

end module JAQL_SFC_Soil
