!###############################################################################
!
! LandUse_Soil  -  Read soil texture database, convert
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_LandUse_Soil

  use GO, only : gol, goPr, goErr

  use JAQL_SFC_Soil, only : nssd

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  n_soiltype
  public  ::  i_soiltype_urban_and_industrial_area
  public  ::  i_soiltype_mineral_extraction_site
  public  ::  i_soiltype_arable_land
  public  ::  i_soiltype_rice_fields
  public  ::  i_soiltype_permanent_crops
  public  ::  i_soiltype_pastures
  public  ::  i_soiltype_complex_cultivation_patterns
  public  ::  i_soiltype_mixed_agriculture_and_natural
  public  ::  i_soiltype_mixed_forest
  public  ::  i_soiltype_broad_leaved_forest
  public  ::  i_soiltype_coniferous_forest
  public  ::  i_soiltype_natural_grassland
  public  ::  i_soiltype_shrub_area
  public  ::  i_soiltype_beaches_dunes_sands
  public  ::  i_soiltype_bare_rock
  public  ::  i_soiltype_sparsely_vegetated_area
  public  ::  i_soiltype_ice
  public  ::  i_soiltype_inland_wetland
  public  ::  i_soiltype_coastal_wetland
  public  ::  i_soiltype_continental_water
  public  ::  i_soiltype_marine_water

  public  ::  n_soiltexture

  !!!!!!!!!!!!!!!!!!!! old style !!!!!!!!!!!!!!!!!!!!
  !public  ::  i_soiltexture_coarse
  !public  ::  i_soiltexture_medium
  !public  ::  i_soiltexture_medium_fine
  !public  ::  i_soiltexture_fine
  !public  ::  i_soiltexture_very_fine
  !!!!!!!!!!!!!!!!!!!! old style !!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!! Ina Tegen style !!!!!!!!!!!!!!!!!!!!
  !public ::  i_soiltexture_coarse
  !public ::  i_soiltexture_medium
  !public ::  i_soiltexture_fine
  !public ::  i_soiltexture_coarse_medium
  !public ::  i_soiltexture_coarse_fine
  !public ::  i_soiltexture_medium_fine
  !public ::  i_soiltexture_coarse_medium_fine
  !!!!!!!!!!!!!!!!!!!! Ina Tegen style !!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! new style !!!!!!!!!!!!!!!!!!!!
  public    ::  i_soiltexture_coarse
  public    ::  i_soiltexture_medium
  public    ::  i_soiltexture_medium_fine
  public    ::  i_soiltexture_fine
  public    ::  i_soiltexture_very_fine
  
  public    ::  n_soiltexture_orig
  
  public  ::  i_soiltexture_sand_orig
  public  ::  i_soiltexture_loamy_sand_orig
  public  ::  i_soiltexture_sandy_loam_orig
  public  ::  i_soiltexture_silt_loam_orig
  public  ::  i_soiltexture_silt_orig
  public  ::  i_soiltexture_loam_orig
  public  ::  i_soiltexture_sandy_clay_loam_orig
  public  ::  i_soiltexture_silty_clay_loam_orig
  public  ::  i_soiltexture_clay_loam_orig
  public  ::  i_soiltexture_sandy_clay_orig
  public  ::  i_soiltexture_silty_clay_orig
  public  ::  i_soiltexture_clay_orig
  !!!!!!!!!!!!!!!!!!!! new style !!!!!!!!!!!!!!!!!!!!
  
  public    ::  clay_frac
  public    ::  ssd_frac

  public  ::  n_erodible
  public  ::  i_erodible_arable, i_erodible_bare
  public  ::  erodible_soil, erodible_soil_texture

  public  ::  rho_soil

  public  ::  LandUse_Soil_Init, LandUse_Soil_Done


  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LandUse_Soil'


  ! --- soil choice -----------------------
  
  !character(len=180), parameter  ::  soil_choice = 'Tegen'
  !character(len=180), parameter  ::  soil_choice = 'LE_old'
  character(len=180), parameter ::  soil_choice = 'LE_new'
  
  
  !
  ! Soil textures:
  !  --------------------------------------------------
  !  code label
  !  ---- ---------------------------------------------
  !   0   No information  clay     sand
  !   1   Coarse          0-18%   65-100%
  !   2   Medium         18-35%   15-100%  or
  !                       0-18%   15- 65%
  !   3   Medium fine     0-35%    0- 15%
  !   4   Fine           35-60%
  !   5   Very fine      60-100%
  !   6   No texture    (other cases)
  !   7   No texture    (because of rock outcrop)
  !   8   No texture    (because of organic layer)
  !  --------------------------------------------------
  !
  !      100% +-----------------------+
  !           |    |    |      |      |
  !           |  1 |    |      |      |
  !   s   65% +----+    |      |      |
  !   a       |         |      |      |
  !   n       |    2    |  4   |  5   |
  !   d       |         |      |      |
  !       15% +---------+      |      |
  !           |    3    |      |      |
  !        0% +----+----+------+------+
  !          0%   18%  35%    60%    100%
  !                        clay
  !
  ! upper value of soil texture labels:
  
  !!!!!!!!!!!!!!!!!!!! old style !!!!!!!!!!!!!!!!!!!!

  !integer, parameter  ::  n_soiltexture = 5   ! 1-5, skip the no-info and no-texture parts
  !integer, parameter  ::  n_soiltexture_max = 6 ! for checking the input file
  !integer, parameter  ::  n_soiltexture_min = 0 ! for checking the input file
  !character(len=64), parameter ::  soil_name = 'soil_measure'
  !! indices:
  !integer, parameter  ::  i_soiltexture_coarse      = 1
  !integer, parameter  ::  i_soiltexture_medium      = 2
  !integer, parameter  ::  i_soiltexture_medium_fine = 3
  !integer, parameter  ::  i_soiltexture_fine        = 4
  !integer, parameter  ::  i_soiltexture_very_fine   = 5
  !!
  !! associated clay fractions; soiltext class :   0    1     2    3     4     5     6    7    8
  !real, parameter :: clay_frac(n_soiltexture) = (/      0.1, 0.25, 0.25, 0.50, 0.65              /)
  !! soil density: estimate for a soil with equal fractions
  !! of sand, silt and clay particles:
  !real, parameter :: rho_soil = 1.3e3     ! kg/m3

  !! Mapping from soil texture database to soil size modes
  !! defined in JAQL_SFC_Soil :
  !!  real, parameter :: ssd_frac(n_soiltexture,nssd) = &
  !!              !   1     2    3     4     5
  !!           (/ (/ 0.1, 0.25, 0.25, 0.5,  0.65 /), &  ! FFS  "fery fine soil"
  !!              (/ 0.1, 0.25, 0.4,  0.4,  0.35 /), &  ! FS   "fine soil"
  !!              (/ 0.0, 0.40, 0.25, 0.1,  0.0  /), &  ! MS   "medium soil"
  !!              (/ 0.8, 0.10, 0.1,  0.0 , 0.0  /) /)  ! CS   "coarse soil"
  !real, parameter :: ssd_frac(n_soiltexture,nssd) =  reshape( &
  !!          !   1     2    3     4     5
  !    (/ 0.1, 0.25, 0.25, 0.5,  0.65, &     ! FFS  "fery fine soil"
  !     0.1, 0.25, 0.4,  0.4,  0.35, &     ! FS   "fine soil"
  !     0.0, 0.40, 0.25, 0.1,  0.0 , &     ! MS   "medium soil"
  !     0.8, 0.10, 0.1,  0.0 , 0.0  /), &  ! CS   "coarse soil"
  !    (/n_soiltexture,nssd/) )
  !!!!!!!!!!!!!!!!!!! old style !!!!!!!!!!!!!!!!!!!!
  

  !!!!!!!!!!!!!!!!!!!! Ina Tegen style !!!!!!!!!!!!!!!!!!!!
  !integer, parameter  ::  n_soiltexture = 7   ! 1-7, skip the no-info and no-texture parts
  !integer , parameter ::  n_soiltexture_max = 14 ! for checking the input file
  !integer, parameter  ::  n_soiltexture_min = 1 ! for checking the input file
  !character(len=64), parameter ::  soil_name = 'type'
  ! indices:
  !integer, parameter ::  i_soiltexture_coarse        = 1
  !integer, parameter ::  i_soiltexture_medium        = 2
  !integer, parameter ::  i_soiltexture_fine          = 3
  !integer, parameter ::  i_soiltexture_coarse_medium     = 4
  !integer, parameter ::  i_soiltexture_coarse_fine     = 5
  !integer, parameter ::  i_soiltexture_medium_fine     = 6
  !integer, parameter ::  i_soiltexture_coarse_medium_fine  = 7
  !! associated clay fractions; soiltext class :  1 2 3 4 5 6 7 8 9 ''
  !real, parameter ::   clay_frac(n_soiltexture) = (/0.0, 0.3, 0.67, 0.2, 0.38, 0.48, 0.35/)
  !! soil density: estimate for a soil with equal fractions
  !! of sand, silt and clay particles:
  !real, parameter :: rho_soil = 1.3e3     ! kg/m3
  !real, parameter :: ssd_frac(n_soiltexture,nssd) =  reshape( &
  !   !1    2   3   4   5   6   7
  !       (/0.43, 0.0,  0.0,  0.1,  0.0,  0.0,  0.23,   & ! Coarse sand
  !         0.4,  0.37, 0.0,  0.5,  0.5,  0.27, 0.23,   & ! Medium/fine sand
  !         0.17, 0.33, 0.33, 0.2,  0.12, 0.25, 0.19,   & ! Silt
  !         0.0,  0.3,  0.67, 0.2,  0.38, 0.48, 0.35/), & ! Clay
  !       (/n_soiltexture,nssd/) )
  !!!!!!!!!!!!!!!!!!!! Ina Tegen style !!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! new soil, new style !!!!!!!!!!!!!!!!!!!!
  !integer, parameter  ::  n_soiltexture = 12   ! 1-12, skip the no-info and no-texture parts
  !integer, parameter  ::  n_soiltexture_max = 16  ! for checking the input file ...
  !integer, parameter  ::  n_soiltexture_min = 1 ! for checking the input file
  !character(len=64), parameter ::  soil_name = 'soil'
  !! indices:
  !integer, parameter ::  i_soiltexture_sand          = 1
  !integer, parameter ::  i_soiltexture_loamy_sand      = 2
  !integer, parameter ::  i_soiltexture_sandy_loam      = 3
  !integer, parameter ::  i_soiltexture_silt_loam       = 4
  !integer, parameter ::  i_soiltexture_silt          = 5
  !integer, parameter ::  i_soiltexture_loam          = 6
  !integer, parameter ::  i_soiltexture_sandy_clay_loam   = 7
  !integer, parameter ::  i_soiltexture_silty_clay_loam   = 8
  !integer, parameter ::  i_soiltexture_clay_loam       = 9
  !integer, parameter ::  i_soiltexture_sandy_clay      = 10
  !integer, parameter ::  i_soiltexture_silty_clay      = 11
  !integer, parameter ::  i_soiltexture_clay          = 12
  !! the next items are not irodible
  !!integer, parameter  ::  i_soiltexture_organic       = 13
  !!integer, parameter  ::  i_soiltexture_water         = 14
  !!integer, parameter  ::  i_soiltexture_bedrock       = 15
  !!integer, parameter  ::  i_soiltexture_other         = 16
  
  ! !associated clay fractions; soiltext class :  1   2   3   4   5   6       
  !real, parameter ::   clay_frac(n_soiltexture) = (/0.03,  0.06, 0.1,  0.13, 0.05, 0.18, &
  !                   0.27, 0.34, 0.34, 0.42, 0.47, 0.58/)
  !! soil density: estimate for a soil with equal fractions
  !! of sand, silt and clay particles:
  !real, parameter :: rho_soil = 1.3e3     ! kg/m3
  !real, parameter :: ssd_frac(n_soiltexture,nssd) =  reshape( &
  !   !1    2   3   4   5   6   7   8   9   10    11    12
  !      (/ 0.92, 0.82, 0.58, 0.17, 0.1,  0.43, 0.58, 0.1,  0.32, 0.52, 0.06, 0.22, &   ! Sand
  !         0.05, 0.12, 0.32, 0.7,  0.85, 0.39, 0.15, 0.56, 0.34, 0.06, 0.47, 0.2,  & ! Silt
  !         0.03, 0.06, 0.1,  0.13, 0.05, 0.18, 0.27, 0.34, 0.34, 0.42, 0.47, 0.58 /),& ! Clay
  !       (/n_soiltexture,nssd/) )
  !!!!!!!!!!!!!!!!!!!! new soil, new style !!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!! new soil, old style !!!!!!!!!!!!!!!!!!!!
  integer, parameter  ::  n_soiltexture_orig    = 12   ! 1-12, skip the no-info and no-texture parts
  !integer, parameter ::  n_soiltexture       = 5
  integer, parameter  ::  n_soiltexture     = 12
  integer, parameter  ::  n_soiltexture_max_orig  = 16  ! for checking the input file ...
  !integer, parameter ::  n_soiltexture_max   = 6
  integer, parameter  ::  n_soiltexture_max   = 16
  integer, parameter  ::  n_soiltexture_min     = 1 ! for checking the input file
  character(len=64), parameter  ::  soil_name = 'soil'
  ! indices:
  integer, parameter  ::  i_soiltexture_coarse              = 1
  integer, parameter  ::  i_soiltexture_medium              = 2
  integer, parameter  ::  i_soiltexture_medium_fine         = 3
  integer, parameter  ::  i_soiltexture_fine                = 4
  integer, parameter  ::  i_soiltexture_very_fine           = 5
  
  integer, parameter  ::  i_soiltexture_sand_orig         = 1
  integer, parameter  ::  i_soiltexture_loamy_sand_orig     = 2
  integer, parameter  ::  i_soiltexture_sandy_loam_orig     = 3
  integer, parameter  ::  i_soiltexture_silt_loam_orig      = 4
  integer, parameter  ::  i_soiltexture_silt_orig         = 5
  integer, parameter  ::  i_soiltexture_loam_orig         = 6
  integer, parameter  ::  i_soiltexture_sandy_clay_loam_orig    = 7
  integer, parameter  ::  i_soiltexture_silty_clay_loam_orig    = 8
  integer, parameter  ::  i_soiltexture_clay_loam_orig      = 9
  integer, parameter  ::  i_soiltexture_sandy_clay_orig     = 10
  integer, parameter  ::  i_soiltexture_silty_clay_orig     = 11
  integer, parameter  ::  i_soiltexture_clay_orig         = 12
  !! the next items are not irodible
  !integer, parameter ::  i_soiltexture_organic       = 13
  !integer, parameter ::  i_soiltexture_water         = 14
  !integer, parameter ::  i_soiltexture_bedrock       = 15
  !integer, parameter ::  i_soiltexture_other         = 16
      
  real, parameter ::  clay_frac(n_soiltexture) = (/ 0.03, 0.06, 0.1,  0.13, 0.05, 0.18, &
                      0.27, 0.34, 0.34, 0.42, 0.47, 0.58/)
                      
  !real, parameter :: clay_frac(n_soiltexture) = (/0.1, 0.25, 0.25, 0.50, 0.65 /)
  !real, parameter :: clay_frac(n_soiltexture) = (/0.0, 0.0,  0.0,  0.0,  0.0 /)

  ! soil density: estimate for a soil with equal fractions
  ! of sand, silt and clay particles:
  !real, parameter :: rho_soil = 1.3e3     ! kg/m3
  real, parameter :: rho_soil = 2.65e3  ! kg/m3 [value obtained from T Cheng et al. Improved dust emission....]

  ! real, parameter :: ssd_frac(n_soiltexture,nssd) =  reshape( &
    !     !   1   2     3   4   5 
  !    (/0.1, 0.25, 0.25, 0.5,  0.65,   &     ! FFS  "fery fine soil" 
  !    0.1, 0.25, 0.4,  0.4,  0.35,     &     ! FS   "fine soil"
  !    0.0, 0.4,  0.25, 0.1,  0.0,    &       ! MS   "medium soil"
  !    0.8, 0.1,  0.1,  0.0,  0.0   /), &     ! CS   "coarse soil"
  !    (/n_soiltexture,nssd/) )
  !real, parameter ::   clay_frac(n_soiltexture) = (/ 0.03, 0.06, 0.10, 0.13, 0.05, 0.18, &
  !                           0.27, 0.34, 0.34, 0.42, 0.47, 0.58  /)
  real, parameter :: ssd_frac(n_soiltexture,nssd) =  reshape( &
         !   1    2     3   4   5   6   7   8   9   10    11    12
      (/1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &     ! Sand
      0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  &     ! 
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0 /), &     ! Clay
       (/n_soiltexture,nssd/) )
  
  !real, parameter :: ssd_frac(n_soiltexture,nssd) =  reshape( &
  !   !1    2   3   4   5   6   7   8   9   10    11    12
  !      (/ 0.92, 0.82, 0.58, 0.17, 0.1,  0.43, 0.58, 0.1,  0.32, 0.52, 0.06, 0.22, &   ! Sand
  !         0.05, 0.12, 0.32, 0.7,  0.85, 0.39, 0.15, 0.56, 0.34, 0.06, 0.47, 0.2,  & ! Silt
  !         0.03, 0.06, 0.1,  0.13, 0.05, 0.18, 0.27, 0.34, 0.34, 0.42, 0.47, 0.58 /),& ! Clay
  !       (/n_soiltexture,nssd/) )
  
  !
  ! The data file provides for each location the "soil texture"
  ! following the above definition, but also a "soil type".
  ! Heare a table with the distribution over soil texture classes
  ! per soil type, generated with the test code in LandUse_Soil_Init.
  ! The dominant soil texture is marked with a '*',
  ! which seems to be "medium" for almost all types.
  !
  !  ---------------------------------- ----- ----- ----- ----- ----- ----- ----- ----- -----
  !                                      soil texture class %
  !  soil type                              0     1     2     3     4     5     6     7     8
  !                                         -     C     M    MF     F    VF     -     -     -
  !  ---------------------------------- ----- ----- ----- ----- ----- ----- ----- ----- -----
  !   1  Urban and Industrial Area         19    20    38*   15     7     0     0     0     1
  !   2  Mineral Extraction Site            6    28    41*   15     6     3     0     0     1
  !   3  Arable Land                        1    18    40*   25    13     1     0     0     2
  !   4  Rice Fields                        4     7    67*   10     8     1     0     0     2
  !   5  Permanent Crops                    1    10    53*   17    18     0     0     0     0
  !   6  Pastures                           1    17    49*   16    10     0     0     0     6
  !   7  Complex Cultivation Patterns       1    20    51*   17    10     0     0     0     1
  !   8  Mixed Agriculture and Natural      1    22    50*   14    10     1     0     0     2
  !   9  Mixed Forest                       1    35    42*    8     4     0     0     0    10
  !  10  Broad Leaved Forest                1    17    54*   18     8     0     0     0     1
  !  11  Coniferous Forest                  1    33    52*    5     3     0     0     0     6
  !  12  Natural Grassland                  2    17    61*    8     6     1     0     0     7
  !  13  Shrub Area                         3    20    55*    6     6     0     0     0    10
  !  14  Beaches, Dunes, Sands             21    26    36*   11     5     0     0     0     2
  !  15  Bare Rock                          9    16    65*    4     2     0     1     0     3
  !  16  Sparsely Vegetated Area            3    17    62*    5     5     0     0     0     8
  !  17  Ice                               35    20    36*    0     0     0     8     0     0
  !  18  Inland Wetland                     3    22    39*    3     2     0     0     0    32
  !  19  Coastal Wetland                   45*   13    32     4     5     0     0     0     1
  !  20  Continental Water                 27    22    37*    5     3     0     0     0     5
  !  21  Marine Water                      62*    9    16     5     5     0     0     0     2
  !  ---------------------------------- ----- ----- ----- ----- ----- ----- ----- ----- -----
  !
  ! number of soil types:
  integer, parameter  ::  n_soiltype = 21
  ! index parameters:
  integer, parameter  ::  i_soiltype_urban_and_industrial_area     =  1
  integer, parameter  ::  i_soiltype_mineral_extraction_site       =  2
  integer, parameter  ::  i_soiltype_arable_land                   =  3
  integer, parameter  ::  i_soiltype_rice_fields                   =  4
  integer, parameter  ::  i_soiltype_permanent_crops               =  5
  integer, parameter  ::  i_soiltype_pastures                      =  6
  integer, parameter  ::  i_soiltype_complex_cultivation_patterns  =  7
  integer, parameter  ::  i_soiltype_mixed_agriculture_and_natural =  8
  integer, parameter  ::  i_soiltype_mixed_forest                  =  9
  integer, parameter  ::  i_soiltype_broad_leaved_forest           = 10
  integer, parameter  ::  i_soiltype_coniferous_forest             = 11
  integer, parameter  ::  i_soiltype_natural_grassland             = 12
  integer, parameter  ::  i_soiltype_shrub_area                    = 13
  integer, parameter  ::  i_soiltype_beaches_dunes_sands           = 14
  integer, parameter  ::  i_soiltype_bare_rock                     = 15
  integer, parameter  ::  i_soiltype_sparsely_vegetated_area       = 16
  integer, parameter  ::  i_soiltype_ice                           = 17
  integer, parameter  ::  i_soiltype_inland_wetland                = 18
  integer, parameter  ::  i_soiltype_coastal_wetland               = 19
  integer, parameter  ::  i_soiltype_continental_water             = 20
  integer, parameter  ::  i_soiltype_marine_water                  = 21

  ! selected erodible soils for use in wind-blown dust module:
  integer, parameter  ::  n_erodible = 2
  integer, parameter  ::  i_erodible_arable = 1
  integer, parameter  ::  i_erodible_bare   = 2



  ! --- var -------------------------------

  ! maps on LOTOS-EUROS resolution:
  logical, allocatable  ::  erodible_soil        (:,:)      ! (nx,nx)
  real, allocatable     ::  erodible_soil_texture(:,:,:,:)  ! (nx,ny,n_erodible,n_soiltexture)


contains


  ! ========================================================================


  subroutine LandUse_Soil_Init( rcF, status, arable_soil, bare_soil )

    use GO            , only : TrcFile, ReadRc
    use LE_Grid       , only : ugg
    use C3PO          , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors
  
    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status
    real, intent(in), optional  ::  arable_soil(:,:)   ! (nx,ny) fracion with arable soil [0-1]
    real, intent(in), optional  ::  bare_soil(:,:)     ! (nx,ny) fracion with bare  soil [0-1]

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_LandUse_Soil_Init'

    ! --- local ------------------------------------------

    integer                 ::  i
  
    character(len=512)      ::  fname
    integer                 ::  fu
    character(len=256)      ::  line
    integer                 ::  iline
    real                    ::  lon, lat, dlon, dlat
    integer                 ::  ix, iy
    integer                 ::  i_soiltype, i_soiltext
    real                    ::  soiltexture_r
    integer                 ::  i_soiltexture
    integer                 ::  i_erodible
    integer                 ::  i_erodible_soil
    integer                 ::  cmode
    integer                 ::  ncid, dimid, varid
    real, allocatable       ::  erodible_soil_fraction(:,:)    ! (nx,ny)
  
    integer                 ::  i_filetype
    integer, parameter      ::  i_filetype_dominant   = 1
    integer, parameter      ::  i_filetype_fractional   = 2
    integer                 ::  level
    integer                 ::  j

    ! soilmeasure
    real, allocatable       ::  soil_measureFX(:,:,:)         !(nlon, nlat, n_soiltexture)
    real, allocatable       ::  soil_measureDX(:,:)           !(nlon, nlat)
    integer, allocatable    ::  soil_measureX(:,:)    ! (nlon,nlat) on input grid
    real, allocatable       ::  soil_measure(:,:,:)   ! (nx,ny,n_soiltexture)
    real, allocatable       ::  soil_measureXX(:,:,:) !(nlon, nlat,n_soiltexture+1)

    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    character(len=64)     ::  units_in
    character(len=64)     ::  description
    integer               ::  start_ind(2)
    integer               ::  count_ind(2)
    
    !! testing: table to find relation between soil type and soil texture:
    !real   ::  soiltypetexture(nsoil,0:usoiltext)

    ! --- begin ------------------------------------------

    ! info ...
    write( gol, '("read landuse data base with soil texture information ...")' ) ; call goPr

    ! map with erodible soil texture
    allocate( erodible_soil         (ugg%nlon,ugg%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( erodible_soil_texture (ugg%nlon,ugg%nlat,n_erodible,n_soiltexture), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! init to none:
    erodible_soil         = .false.
    erodible_soil_texture = 0.0

    ! read filename:
    call ReadRc( rcF, 'landuse.soil_texture.file', fname, status )
    IF_NOTOK_RETURN(status=1)
    ! variable description (should be in settings in future ...)
    description = 'var_name=soil'

    ! info ...
    write( gol, '("  input file: ",a)') trim(fname); call goPr
  
    ! determine file type
    if ( index(fname, 'dominant') > 0 ) then
      i_filetype = i_filetype_dominant
      write( gol, '("  file type: dominant")'); call goPr
    else if ( index(fname, 'fractional') > 0 ) then
      i_filetype = i_filetype_fractional
      write( gol, '("  file type: fractional")'); call goPr
    else
      write (gol, '("Soil data file is neither dominant nor fractional!")') ; call goErr
      TRACEBACK;status=1;return
    end if

    ! initialize soil measure on target grid (from LE-ugg)
    allocate( soil_measure( ugg%nlon,ugg%nlat,n_soiltexture), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( erodible_soil_fraction(ugg%nlon,ugg%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    call file_in%Open( trim(fname), status )
    IF_NOTOK_RETURN(status=1)

    ! variable id:
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! init grid definition, clip to target grid and return input range:
    call file_in%Get_Grid( varid, grid_in, status, &
                            ugg_to=ugg, start_ind=start_ind, count_ind=count_ind )
    IF_NOTOK_RETURN(status=1)

    ! read data values
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( i_filetype == i_filetype_dominant ) then
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! storage:
      allocate( soil_measureDX(grid_in%nlon,grid_in%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)       
      allocate( soil_measureX(grid_in%nlon,grid_in%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( soil_measureXX( ugg%nlon,ugg%nlat,n_soiltexture+1), stat=status )
      IF_NOTOK_RETURN(status=1)

      ! read in data
      call file_in%Get_Var( trim(description), soil_measureDX(:,:), units_in, status, &
                             start=start_ind, count=count_ind )
      IF_NOTOK_RETURN(status=1)

      ! redirect original soiltexture types to 5 classes (coarse, medium, medium-fine, fine, very-fine)
      do i = 1,n_soiltexture_max_orig
        if ( i == 1 .or. i == 2 )  i_soiltext = i_soiltexture_coarse
        if ( i == 3 .or. i == 6  .or. i == 7 )  i_soiltext = i_soiltexture_medium
        if ( i == 4 .or. i == 5  .or. i == 8 )  i_soiltext = i_soiltexture_medium_fine
        if ( i == 9 .or. i == 10 .or. i == 11 ) i_soiltext = i_soiltexture_fine
        if ( i == 12 ) i_soiltext = i_soiltexture_very_fine
        if ( i == 13 .or. i == 14 .or. i == 15 .or. i == 16 ) i_soiltext = n_soiltexture + 1  ! rest category          
        where ( soil_measureDX == i )  soil_measureX = i_soiltext
      end do ! original soiltextures
      ! regrid:
      call Grid_Convertors%Ugg_IndexFractions( grid_in, soil_measureX, ugg, soil_measureXX, status )
      IF_NOTOK_RETURN(status=1)
      ! pick only first 12 classes:
      soil_measure(:,:,1:n_soiltexture) = soil_measureXX(:,:,1:n_soiltexture)

      ! clear:
      deallocate( soil_measureDX, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( soil_measureX, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( soil_measureXX, stat=status )
      IF_NOTOK_RETURN(status=1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if ( i_filetype == i_filetype_fractional ) then
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! storage:
      allocate( soil_measureFX(grid_in%nlon,grid_in%nlat,n_soiltexture_max), stat=status )
      IF_NOTOK_RETURN(status=1)

      ! read:
      call file_in%Get_Var( trim(description), soil_measureFX(:,:,:), units_in, status, &
                             start=(/start_ind(1),start_ind(2),1/), &
                             count=(/count_ind(1),count_ind(2),n_soiltexture_max/) )
      ! scale:
      soil_measureFX = 0.01 * soil_measureFX

      ! regrid:                            
      do i_soiltexture = 1, n_soiltexture
        call Grid_Convertors%Ugg_AreaAver( grid_in, soil_measureFX(:,:,i_soiltexture), ugg, soil_measure(:,:,i_soiltexture), status )
        IF_NOTOK_RETURN(status=1)
      end do

      ! clear:
      deallocate( soil_measureFX, stat=status )
      IF_NOTOK_RETURN(status=1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
      write( gol, '("Unknown filetype for soiltexture file: ",i0)' ) i_filetype ; call GoErr
      TRACEBACK;status=1;return
    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! clear:
    call grid_in%Done( status )
    IF_NOTOK_RETURN(status=1)
    ! close:
    call file_in%Close( status )
    IF_NOTOK_RETURN(status=1)

    ! extend to fraction per erodible soil and texture class:
    do i_erodible_soil = 1, n_erodible
    ! check ...
      select case ( i_erodible_soil )
        case ( i_erodible_arable )
          if ( .not. present(arable_soil) ) then
            write (gol,'("argument arable_soil should be present")'); call goErr
            TRACEBACK; status=1; return
          end if
          erodible_soil_fraction = arable_soil
        case ( i_erodible_bare )
          if ( .not. present(bare_soil) ) then
            write (gol,'("argument bare_soil should be present")'); call goErr
            TRACEBACK; status=1; return
          end if
          erodible_soil_fraction = bare_soil
        case default
          write (gol,'("unsupported erodible soil index ",i4)'); call goErr
          TRACEBACK; status=1; return
      end select

      ! loop over supported soil textures:
      do i_soiltexture = 1, n_soiltexture
        ! fill fractions:
        erodible_soil_texture(:,:,i_erodible_soil,i_soiltexture) = soil_measure(:,:,i_soiltexture) * erodible_soil_fraction(:,:)
      end do ! textures
    end do ! erodible soils

    ! set flags if grid cell could have any erosion at all:
    erodible_soil = .false.
    do i_erodible_soil = 1, n_erodible
      do i_soiltexture = 1, n_soiltexture
        erodible_soil = erodible_soil .or. ( erodible_soil_texture(:,:,i_erodible_soil,i_soiltexture) > 0 )
      end do        
    end do     

    ! clear:
    deallocate( soil_measure, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( erodible_soil_fraction, stat=status )
    IF_NOTOK_RETURN(status=1)
        
    ! ok
    status = 0

  end subroutine LandUse_Soil_Init


  ! ***


  subroutine LandUse_Soil_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_LandUse_Soil_Done'

    ! --- begin ----------------------------------

    ! clear:
    deallocate( erodible_soil, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( erodible_soil_texture, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

    end subroutine LandUse_Soil_Done

end module LE_LandUse_Soil
