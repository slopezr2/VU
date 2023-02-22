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

    use GO     , only : TrcFile, ReadRc
    use GO     , only : goGetFU
    !use Grid   , only : GetLocation
    use Grid        , only : TllGridInfo, Init, Done, IndexFractions, FillGrid_AreaAverage

    use Dims   , only : nx, ny
    use LE_Grid, only : lli

#ifdef with_netcdf    
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att, NF90_NOWRITE, NF90_NOERR
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_StrError
#endif
  
    use file_nc

    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status
    real, intent(in), optional  ::  arable_soil(:,:)   ! (nx,ny) fracion with arable soil [0-1]
    real, intent(in), optional  ::  bare_soil(:,:)     ! (nx,ny) fracion with bare  soil [0-1]

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_LandUse_Soil_Init'

    ! --- local ------------------------------------------

  integer                 ::  l, k, i
    character(len=3)        ::  ext
  
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
    integer, allocatable    ::  nval(:,:)
    integer                 ::  cmode
    integer                 ::  ncid, dimid, varid
    integer                 ::  nlon, nlat, ilon, ilat
    real, allocatable       ::  lons(:), lats(:)
    real, allocatable       ::  lon_bounds(:,:), lat_bounds(:,:) ! (2, nlon)
    integer, allocatable    ::  soil_measureX(:,:)    ! (nlon,nlat) on input grid
    real, allocatable       ::  soil_measure(:,:,:)   ! (nx,ny,n_soiltexture)
    real, allocatable       ::  erodible_soil_fraction(:,:)    ! (nx,ny)
    !integer, allocatable    ::  erodible_soil_texture_counter(:,:)
  
    integer                 ::  i_filetype
    integer, parameter      ::  i_filetype_dominant   = 1
    integer, parameter      ::  i_filetype_fractional   = 2
    integer                 ::  level
    integer                 ::  j
    character(len=64)       ::  dumpname

    !real, allocatable    ::  soil_measure_fractionalX(:,:,:)     !(nlon, nlat, level)
    byte, allocatable   ::  soil_measure_fractionalX(:,:,:)     !(nlon, nlat, level)
    real, allocatable   ::  soil_measureFX(:,:,:)         !(nlon, nlat, n_soiltexture)
    real, allocatable   ::  soil_measureDX(:,:)           !(nlon, nlat)

    type(TllGridInfo)       ::  lliX

    !! testing: table to find relation between soil type and soil texture:
    !real   ::  soiltypetexture(nsoil,0:usoiltext)

    ! --- begin ------------------------------------------

    ! map with erodible soil texture
    allocate( erodible_soil        (nx,ny)                          )
    allocate( erodible_soil_texture(nx,ny,n_erodible,n_soiltexture) )
    ! init to none:
    erodible_soil         = .false.
    erodible_soil_texture = 0.0

    ! counter:
    allocate( nval(nx,ny) )

    !! testing: init table to zero:
    !soiltypetexture = 0.0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Read new landuse data base (based on Corine 2000) with soil texture information
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! info ...
    write (gol,'("read landuse data base with soil texture information ...")'); call goPr

    ! read filename:
    call ReadRc( rcF, 'landuse.soil_texture.file', fname, status )
    IF_NOTOK_RETURN(status=1)
  
    ! determine file type
    if ( index(fname, 'dominant') > 0 ) then
      i_filetype = i_filetype_dominant
    else if ( index(fname, 'fractional') > 0 ) then
      i_filetype = i_filetype_fractional
    else
      write (gol, '("Soil data file is neither dominant nor fractional!")') ; call goErr
      TRACEBACK;status=1;return
    end if
  
    ! extract extension:
    l = len_trim(fname)
    k = index( fname(1:l), '.', back=.true. )
    if ( (k == 0) .or. (l-k < 1) .or. (l-k > 3) ) then
      write (gol,'("could not extract 3-character extension of : ",a)') trim(fname); call goErr
      IF_NOTOK_RETURN(status=1)
    end if
    ext = fname(k+1:l)
  
    select case( trim(ext) )
      
      case ( 'csv' )
      
        write (gol,'("regridding of soil texture with csv format not ok")'); call goErr
        TRACEBACK; status=1; return

    
      case ( 'nc' )   
#ifdef with_netcdf
        ! set open mode flag:
        cmode = NF90_NOWRITE   ! read-only

        ! open file:
        status = NF90_Open( trim(fname), cmode, ncid )
        if (status/=NF90_NOERR) then
          write (gol,'("NF90 error: ",a)') trim(NF90_StrError(status)); call goErr
          write (gol,'("opening file : ",a)') trim(fname); call goErr
          TRACEBACK; status=1; return
        end if

        ! extract lon/lat
        status = NF90_Inq_DimID( ncid, 'lon', dimid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! extract dimension length:
        status = NF90_Inquire_Dimension( ncid, dimid, len=nlon )
        IF_NF90_NOTOK_RETURN(status=1)

        status = NF90_Inq_DimID( ncid, 'lat', dimid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! extract dimension length:
        status = NF90_Inquire_Dimension( ncid, dimid, len=nlat )
        IF_NF90_NOTOK_RETURN(status=1)

        ! storage:
        allocate ( lons( nlon ) )
        allocate ( lats( nlat ) ) 
        allocate( soil_measure(nx,ny,n_soiltexture_min:n_soiltexture) )
        allocate( erodible_soil_fraction(nx,ny) )


        status = NF90_Inq_VarID( ncid, 'lon', varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Get_Var( ncid, varid, lons )
        IF_NF90_NOTOK_RETURN(status=1)

        status = NF90_Inq_VarID( ncid, 'lat', varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Get_Var( ncid, varid, lats )
        IF_NF90_NOTOK_RETURN(status=1)

        if ( i_filetype == i_filetype_dominant ) then
          allocate( soil_measureDX(nlon,nlat) )
          allocate( lon_bounds(2,nlon) )
          allocate( lat_bounds(2,nlat) )

          status = NF90_Inq_VarID( ncid, 'lon_bounds', varid )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Get_Var( ncid, varid, lon_bounds )
          IF_NF90_NOTOK_RETURN(status=1)

          status = NF90_Inq_VarID( ncid, 'lat_bounds', varid )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Get_Var( ncid, varid, lat_bounds )
          IF_NF90_NOTOK_RETURN(status=1)

          ! regrid
          call Init( lliX, lon_bounds, lat_bounds, status )
          IF_NOTOK_RETURN(status=1)

          ! extract soiltexture types
          status = NF90_Inq_VarID( ncid, trim(soil_name), varid )
          IF_NF90_NOTOK_RETURN(status=1)

          status = NF90_Get_Var( ncid, varid, soil_measureDX )
          IF_NF90_NOTOK_RETURN(status=1)

          ! storage:
          allocate( soil_measureX(nlon,nlat), stat=status )
          IF_NOTOK_RETURN(status=1)
          ! init with dummy:
          soil_measureX = -999
          ! assign texture classes used in model to original classes;
          ! loop over original soil texture classes:
          do i = 1,n_soiltexture_max_orig
            ! assign model class to original:
            select case( i )
              case ( 1 )
                i_soiltext = i_soiltexture_coarse
              case ( 2 )
                i_soiltext = i_soiltexture_coarse
              case ( 3 )
                i_soiltext = i_soiltexture_medium
              case ( 4 )
                i_soiltext = i_soiltexture_medium_fine
              case ( 5 )
                i_soiltext = i_soiltexture_medium_fine
              case ( 6 )
                i_soiltext = i_soiltexture_medium
              case ( 7 )
                i_soiltext = i_soiltexture_medium
              case ( 8 )
                i_soiltext = i_soiltexture_medium_fine
              case ( 9 )
                i_soiltext = i_soiltexture_fine
              case ( 10 )
                i_soiltext = i_soiltexture_fine
              case ( 11 )
                i_soiltext = i_soiltexture_fine
              case ( 12 )
                i_soiltext = i_soiltexture_very_fine
              case default
                i_soiltext = n_soiltexture_max
            end select
            ! fill soil texture index in corresponding cells:
            where ( soil_measureDX == i ) 
              soil_measureX = i_soiltext
            endwhere
          end do
          
          ! clear:
          deallocate( soil_measureDX, stat=status )
          IF_NOTOK_RETURN(status=1)

        else
          ! extract levels
          status = NF90_Inq_DimID( ncid, 'level', dimid )
          IF_NF90_NOTOK_RETURN(status=1)
          ! extract dimension length:
          status = NF90_Inquire_Dimension( ncid, dimid, len=level )
          IF_NF90_NOTOK_RETURN(status=1)

          allocate( soil_measure_fractionalX(nlon,nlat,level) )

          call Init( lliX, lons(1), lons(2)-lons(1), nlon, lats(1), lats(2)-lats(1), nlat, status)
          IF_NOTOK_RETURN(status=1)

          ! extract soiltexture types
          status = NF90_Inq_VarID( ncid, trim(soil_name), varid )
          IF_NF90_NOTOK_RETURN(status=1)

          status = NF90_Get_Var( ncid, varid, soil_measure_fractionalX )
          IF_NF90_NOTOK_RETURN(status=1)


          allocate( soil_measureFX(nlon,nlat,n_soiltexture_max) )
          soil_measureFX = 0.0
          !do i = 1,n_soiltexture_max_orig
            !select case( i )
              !case ( 1 )
              ! soil_measureFX(:,:,1) = soil_measureFX(:,:,1) + soil_measure_fractionalX(:,:,i)
              !case ( 2 )
              ! soil_measureFX(:,:,1) = soil_measureFX(:,:,1) + soil_measure_fractionalX(:,:,i)
              !case ( 3 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + soil_measure_fractionalX(:,:,i) 
              !case ( 4 )
              ! soil_measureFX(:,:,3) = soil_measureFX(:,:,3) + soil_measure_fractionalX(:,:,i)
              !case ( 5 )
              ! soil_measureFX(:,:,3) = soil_measureFX(:,:,3) + soil_measure_fractionalX(:,:,i)
              !case ( 6 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + soil_measure_fractionalX(:,:,i)
              !case ( 7 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + soil_measure_fractionalX(:,:,i)
              !case ( 8 )
              ! soil_measureFX(:,:,3) = soil_measureFX(:,:,3) + soil_measure_fractionalX(:,:,i)
              !case ( 9 )
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + soil_measure_fractionalX(:,:,i)
              !case ( 10 )
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + soil_measure_fractionalX(:,:,i)
              !case ( 11 )
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + soil_measure_fractionalX(:,:,i)
              !case ( 12 )            
              ! soil_measureFX(:,:,5) = soil_measureFX(:,:,5) + soil_measure_fractionalX(:,:,i)
              !case default
              ! ! not eroble soil textures

              !!! new method !!!
              !case ( 1 )
              ! soil_measureFX(:,:,1) = soil_measureFX(:,:,1) + soil_measure_fractionalX(:,:,i)
              !case ( 2 )
              ! soil_measureFX(:,:,1) = soil_measureFX(:,:,1) + soil_measure_fractionalX(:,:,i)
              !case ( 3 )
              ! soil_measureFX(:,:,1) = soil_measureFX(:,:,1) + 0.5 * soil_measure_fractionalX(:,:,i)
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + 0.5 * soil_measure_fractionalX(:,:,i)
              !case ( 4 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + 0.66 * soil_measure_fractionalX(:,:,i)
              ! soil_measureFX(:,:,3) = soil_measureFX(:,:,3) + 0.34 * soil_measure_fractionalX(:,:,i)
              !case ( 5 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + 0.1 * soil_measure_fractionalX(:,:,i)
              ! soil_measureFX(:,:,3) = soil_measureFX(:,:,3) + 0.9 * soil_measure_fractionalX(:,:,i)
              !case ( 6 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + soil_measure_fractionalX(:,:,i)
              !case ( 7 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + 0.9 * soil_measure_fractionalX(:,:,i)
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + 0.1 * soil_measure_fractionalX(:,:,i)
              !case ( 8 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + 0.1 * soil_measure_fractionalX(:,:,i)
              ! soil_measureFX(:,:,3) = soil_measureFX(:,:,3) + 0.3 * soil_measure_fractionalX(:,:,i)
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + 0.6 * soil_measure_fractionalX(:,:,i)
              !case ( 9 )
              ! soil_measureFX(:,:,2) = soil_measureFX(:,:,2) + 0.4 * soil_measure_fractionalX(:,:,i)
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + 0.6 * soil_measure_fractionalX(:,:,i)
              !case ( 10 )
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + soil_measure_fractionalX(:,:,i)
              !case ( 11 )
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + soil_measure_fractionalX(:,:,i)
              !case ( 12 )            
              ! soil_measureFX(:,:,4) = soil_measureFX(:,:,4) + 0.5 * soil_measure_fractionalX(:,:,i)
              ! soil_measureFX(:,:,5) = soil_measureFX(:,:,5) + 0.5 * soil_measure_fractionalX(:,:,i)
              !case default
                ! not eroble soil textures
            !end select
          !end do
          ! input is percentage, output is fraction
          !soil_measureFX = soil_measureFX * 0.01
          soil_measureFX(:,:,:) = soil_measure_fractionalX(:,:,1:n_soiltexture_max) * 0.01
          !deallocate( soil_measure_fractionalX )

        end if
        ! close:
        status = NF90_Close( ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      
#else
        write (gol,'("not linked with NetCDF library")'); call goErr
        TRACEBACK; status=1; return
#endif

        ! temporary storage:

        if ( i_filetype == i_filetype_dominant ) then
          if ( (minval(soil_measureX) < n_soiltexture_min) .or. (maxval(soil_measureX) > n_soiltexture_max) ) then
            write (gol,'("found soil measures out of range 0 to ",i4)') n_soiltexture_max; call goErr
            TRACEBACK; status=1; return
          end if

          ! convert from 2D field with indices to 3D field with fraction per cell for each index:
          call IndexFractions( lliX, soil_measureX, lli, n_soiltexture_min, &
                    n_soiltexture, soil_measure, status, n_soiltexture_max )
          IF_NOTOK_RETURN(status=1)

        else
          do i_soiltexture = 1, n_soiltexture
            call FillGrid_AreaAverage( lli, soil_measure(:,:,i_soiltexture), lliX, soil_measureFX(:,:,i_soiltexture), status )
            IF_NOTOK_RETURN(status=1)

            !write (dumpname,'("soiltext_",i,".nc")') i_soiltexture
            !call nc_dump( trim(dumpname), soil_measure(:,:,i_soiltexture), 'var', (/'nx','ny'/), status )
            !IF_NOTOK_RETURN(status=1)
          !end do


          end do        
          !stop('dumped')
        end if

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

          ! check
          !do i = 1,nx
          ! do j = 1, ny
          !   if ( erodible_soil_fraction(i,j) < 0.0 ) then
          !     if ( erodible_soil_fraction(i,j) < -0.01 ) then
          !       ! too small
          !       stop 'temporary stop in landuse_soil, fraction too small'
          !     else
          !       erodible_soil_fraction(i,j) = 0.0
          !     end if
          !   else if ( erodible_soil_fraction(i,j) > 1.0 ) then
          !     if ( erodible_soil_fraction(i,j) < 1.01 ) then
          !       ! too large
          !       stop 'temporary stop in landuse_soil, fraction too large'
          !     else
          !       erodible_soil_fraction(i,j) = 1.0
          !     end if
          !   end if
          ! end do
          !end do

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

        ! clear 

        !deallocate( soil_measure )
        !deallocate( erodible_soil_fraction ) 

        if ( i_filetype == i_filetype_dominant ) then
          ! clear:
          deallocate( soil_measureX )
        else
          deallocate( soil_measureFX )
        end if  
        ! done with griddef:
        call Done( lliX, status )

      case default

        write (gol,'("unsupported extension for landuse soil texture file : ",a)') trim(ext); call goErr
        IF_NOTOK_RETURN(status=1)

    end select
    
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
    deallocate( erodible_soil )
    deallocate( erodible_soil_texture )

    ! ok
    status = 0

    end subroutine LandUse_Soil_Done

end module LE_LandUse_Soil
