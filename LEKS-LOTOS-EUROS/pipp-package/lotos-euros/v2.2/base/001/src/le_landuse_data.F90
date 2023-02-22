!###############################################################################
!
! landuse  - landuse database
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!
#include "le.inc"
!
!###############################################################################

!
!  Proposal for new convention
!  ---------------------------
!
!    2012-03 Roy Wichink Kruit, Arjo Segers (TNO)
!
!    New convention for canopy .
!    Consider a vegetation of trees:
!       vegetation_height = 20  ! m
!    Assume that the 'orography' in the model is at 2/3 of this;
!    the distance from the actual soil to the orogrpahy is called:
!       displacement_height = 2.0/3.0 * vegetation_height
!    The "canopy" (="cover") top is the height of the remaining cover.
!
!    A rule-of-thumb for the roughness lenght is:
!        z0 = 10% veg_height
!
!        -----------------------------------------------------------------+-------
!                                                                         |
!                                                                         |
!                                                                         |
!        --------------------------------------------                     |
!        ^             ( )( )( )( )    ^                                  |
!        |             ( )( )( )( )    | zcanopy_top = 1/3 veg_height     |
!        |             ( )( )( )( )    |                                _ _ _ _ _ "surface" concentration model (z0 + 2.5 m)
!        |             ( )( )( )( )    |                                  |
!        |             ( )( )( )( )    |                                _ _ _ _ _
!        |             ( )( )( )( )    v                                  |      }  z0
!        | veg_height  ( )( )( )( ) ---------------- model surface -------+------         
!        |             ( )( )( )( )    ^                                        
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    | displacement_height = 2/3 veg_height   
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        |             ( )( )( )( )    |
!        v             || |} || ||     v                                         
!       =========================================================================  soil                     
!                                                  ^                               
!                                                  v  orography
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  sea level
!
!    All parameters summarized:
!
!        ------------------  ----------------  ----------  ------- ------   --------
!        landuse class       includes           veg.h.     disp.h. zcanopy   z0 2)
!                                                          2/3 vh  1/3 vh
!        ------------------  ----------------  ----------  ------- ------   --------
!        grass                                  0.30   1)   0.20    0.10      0.03
!        arable              corn               2.500  1)   1.667   0.833     0.25
!        permanent_crops     fruit trees        2.500  1)   1.667   0.833     0.25
!        conifereous_forest                    17.4    2)  11.600   5.800     2.00
!        deciduous_forest                      19.4    2)  12.933   6.467     2.00
!        water                                         0)                     0.002
!        urban                                         0)                     3.00
!        other               tundra, savanne?   2.5    1)   1.667   0.833     0.25
!        desert                                        0)                     0.10
!        ice                                           0)                     0.0005
!        wheat                                         0)                     0.25
!        beech                                         0)                     2.00
!        spruce                                        0)                     2.00
!        semi_natural_grass                            0)                     0.25
!        ------------------  ----------------  ----------  ------- ------   --------
!
!        Notes:
!          0) no vegetation
!          1) use veg.height = 10 * z0
!          2) from exisiting code ...
!

module LE_Landuse_Data

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- const -----------------------------

  character(len=*), parameter, private   ::  mname = 'LE_Landuse_Data'


  integer, parameter    ::  ntree_type            = 200
  integer, parameter    ::  nwaterbody_type       = 4

!  integer, parameter    ::  nlu_ozone_specials    = 18
!  integer, parameter    ::  nlu_normal            = 15
!
!  integer, parameter    ::  ilu_grass             = 1
!  integer, parameter    ::  ilu_arable            = 2
!  integer, parameter    ::  ilu_permanent_crops   = 3
!  integer, parameter    ::  ilu_coniferous_forest = 4
!  integer, parameter    ::  ilu_deciduous_forest  = 5
!  integer, parameter    ::  ilu_water_sea         = 6
!  integer, parameter    ::  ilu_urban             = 7
!  integer, parameter    ::  ilu_other             = 8
!  integer, parameter    ::  ilu_desert            = 9
!  integer, parameter    ::  ilu_ice               = 10
!  integer, parameter    ::  ilu_savanna           = 11
!  integer, parameter    ::  ilu_tropical_forest   = 12
!  integer, parameter    ::  ilu_water_inland      = 13
!  integer, parameter    ::  ilu_mediterrean_scrub = 14
!  integer, parameter    ::  ilu_semi_natural_veg  = 15
!  integer, parameter    ::  ilu_wheat             = 16
!  integer, parameter    ::  ilu_beech             = 17
!  integer, parameter    ::  ilu_spruce            = 18
!
! Landuse definitions matched with depac classes (numbering is as above)
#include "depac_lu.inc"

  ! roughness length for momentum per landuse:
  real, parameter :: z0m_depac(nlu_ozone_specials) = &
            (/0.03, 0.10,  0.25, 2.0,  2.0,   0.002, 2.0,  0.05, 0.013,  0.002, 0.03, 2.0, 0.002, 0.25, 0.02, 0.10, 2.0,  2.0/)
  !           grass arabl  crops conif decid  water  urba  othr  desr      ice   sav   trf  wai    med   sem   wheat beech spruce

  ! roughness length for heat and trace gases per landuse:
  real, parameter :: z0h_depac(nlu_ozone_specials) = &
            (/0.03, 0.10,  0.25, 2.0,  2.0,   0.002, 2.0,  0.05, 0.0013,  0.002,  0.03, 2.0, 0.002, 0.25, 0.02, 0.10, 2.0,  2.0/)
  !           grass  arabl   crops conif decid    water  urba   othr  desr         ice   sav   trf   wai      med   sem    wheat beech spruce

  ! roughness length for dust emissions per landuse
  ! Parameters on local scale instead of regional scale for deposition
  ! (only desert area)
  real, parameter :: z0dust_emis(nlu_ozone_specials) = &
            (/-999., -999.,  -999., -999., -999., -999., -999., -999., 0.0008, -999., -999., -999., -999., -999., -999., -999., -999., -999./)
  !           grass  arabl  crops  conif  decid   water  urba   othr    desr    ce    sav     trf    wai    med   sem    wheat  beech spruce
  
  real, parameter :: z0h_snow_surface = 0.01
  real, parameter :: z0m_snow_surface = 0.01
  
  ! canopy top
  ! as the concentrations are calculated relative to the displacement height, 
  ! the canopy top is 1/3 of the vegetation height, which should always be larger than z0h_depac
  real, parameter :: zcanopytop(nlu_ozone_specials) = &
            (/0.10, 0.33,  0.83, 6.67, 6.67,  0.007, 6.67, 0.17, 0.043, 0.007, 0.10, 6.67, 0.007, 0.83, 0.07, 0.33, 6.67, 6.67/)
  !           grass arabl  crops conif decid  water  urba  othr  desr    ice      sav    trf    wai    med   sem   wheat beech spruce

  !deposition Zhang
  real, parameter :: alfa_lu(nlu_ozone_specials) = &
            (/1.2,  1.2,   1.2,  1.0,  1.0,   100.0, 1.5,  1.2, 50.0, 100.0, 1.2, 1.0, 100.0, 1.2, 50.0, 1.2,  1.0,  1.0/)
  !           grass arabl  crops conif decid  water  urba  othr  desr  ice   sav   trf  wai   med  sem   wheat beech spruce
  real, parameter :: gamma_lu(nlu_ozone_specials) = &
            (/0.54, 0.54,  0.54, 0.56, 0.56,  0.50,  0.56, 0.54, 0.58, 0.50, 0.54, 0.56, 0.50, 0.54, 0.54, 0.54, 0.56, 0.54/)
  !           grass arabl  crops conif decid  water  urba  othr  desr  ice   sav   trf   wai   med   sem   wheat beech spruce
  real, parameter ::A_lu(nlu_ozone_specials) =&
            (/3.0,  3.0,   2.0,  2.0,  7.0,  -99.,   10.0, 3.0, -99., -99.,  3.0, 7.0, -99., 2.0, -99., 3.0,  7.0,  2.0  /)
  !           grass arabl  crops conif decid  water  urba  othr  desr  ice   sav  trf   wai  med   sem wheat beech spruce


  !--------gas deposition------------------
  !

                                 !                             grass    arabl    crops   conif   decid    water     urba     othr  desr    ice     sav    trf    wai     med    wheat   beech  spruce  semi
  real, dimension(nlu_ozone_specials), parameter:: rsmin=   (/  100.0,  100.0,   100.0,  145.0,  232.0,   -999.0,  -999.0,  100.0, -999.0, -999.0, 100.0, 232.0, -999.0, 100.0, 100.0, 100.0,  232.0,  145.0/)  ! used for Baldocchi only???

                                                 !              grass   arabl    crops   conif   decid     water    urba    othr    desr     ice    sav   trf    wai   med   semi   wheat  beech  spruce
  !-- begin Baldocci specific declarations
  real, dimension(nlu_ozone_specials), parameter:: brs      = (/ 66.0,   66.0,    66.0,   22.0,   25.0,    -999.0,  -999.0, 66.0,   -999.0, -999.0, 66.0, 25.0, -999.0, 66.0,  66.0, 66.0,  25.0,   22.0   /)
  real, dimension(nlu_ozone_specials), parameter:: be       = (/  0.0,    0.0,     0.0,    0.0,    0.0,    -999.0,  -999.0,  0.0,   -999.0, -999.0,  0.0,  0.0, -999.0,  0.0,   0.0,  0.0,   0.0,    0.0   /)
  real, dimension(nlu_ozone_specials), parameter:: Tmin_bal = (/  5.0,    5.0,     5.0,    -5.0,  10.0,    -999.0,  -999.0,  5.0,   -999.0, -999.0,  5.0, 10.0, -999.0,  5.0,   5.0,  5.0,  10.0,   -5.0   /)
  real, dimension(nlu_ozone_specials), parameter:: Topt_bal = (/ 25.0,   25.0,    25.0,   9.0,    25.0,    -999.0,  -999.0, 25.0,   -999.0, -999.0, 25.0, 25.0, -999.0, 25.0,  25.0, 25.0,  25.0,    9.0   /)
  real, dimension(nlu_ozone_specials), parameter:: Tmax_bal = (/ 45.0,   45.0,    45.0,   35.0,   45.0,    -999.0,  -999.0, 45.0,   -999.0, -999.0, 45.0, 45.0, -999.0, 45.0,  45.0, 45.0,  45.0,   35.0   /)
  real, dimension(nlu_ozone_specials), parameter ::laimx    = (/  6.0,    5.0,     5.0,    5.0,    5.0,    -999.0,  -999.0,  6.0,   -999.0, -999.0,  5.0,  5.0, -999.0,  5.0,   6.0,  6.0,   5.0,    5.0   /)
  !--end Baldocci specific declarations


  !--end Emberson specific declarations------------

  ! for get_Rb_mcnaughton
  ! cross wind leaf dimension, taken from Simpson et al 2007.
  !                                             grass arabl  crops conif decid  water  urba  othr  desr  ice   sav   trf   wai   med   semi wheat beech spruce
  real, parameter ::  Ld(nlu_ozone_specials) =(/0.00, 0.02,  0.02, 0.04, 0.04,  0.00,  0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 0.00, 0.02, 0.05, 0.02, 0.07, 0.008/)

  ! Wesely
  ! NOTE: gcc does not allow a 2D parameter (/(/../),(/../),../),
  ! therefore use a reshape of a 1D array:
  real, parameter :: riwes(nlu_ozone_specials,4) = reshape( &
                             (/   60,   60,    60,  130,   70,  -999,  -999,   60, -999, -999,  60,   70, -999,   60,   60,   60,   70, 130,     &  ! 1=Summer
                                -999, -999,  -999,  250, -999,  -999,  -999, -999, -999, -999,-999, -999, -999, -999, -999, -999, -999, 250,     &  ! 2=Autumn
                                -999, -999,  -999,  400, -999,  -999,  -999, -999, -999, -999,-999, -999, -999, -999, -999, -999, -999, 400,     &  ! 3=Winter
                                 120,  120,   120,  250,  140,  -999,  -999,  120, -999, -999, 120,  140, -999,  120,  120,  120,  140, 250 /), &  ! 4=Spring
                             (/nlu_ozone_specials,4/) )

  ! o3 stomatal flux threshold (Y) in nmol.m-2 PLA,s-1
  ! (Mills et al. 2004, ICP Mappping Manual, Ch3: Mapping Critical Levels for Vegetation)
  ! note: for o3 threshold stomatal flux Y in a single time step use:
  !   60.0*deltat*fst_th_Y(k) [nmol.m-2 per time step in seconds]
  real, parameter   ::  fst_th_Y(nlu_ozone_specials) = &
                                     (/0.0,  6.0,   6.0,  1.6,  1.6,   0.0,   0.0,  0.0,  0.0,  0.0, 0.0, 1.6, 0.0, 6.0, 1.0, 6.0,  1.0,  1.0/)  ! nmol.m-2 PLA.s-1
  !                                    grass arabl  crops conif decid  water  urba  othr  desr  ice  sav  trf  wai  med  sem wheat beech spruce


  ! maximum number of growing season dependent landuses
  integer, parameter            ::  max_names_gs = 10

  ! factors for temperature dependency no emissions
  real, parameter  :: A_bio_no(nlu_ozone_specials) = &
                                (/ 0.9,  0.9,  0.9,  0.07, 0.07,  0.0,  0.0,  0.0,  0.0,  0.0, 0.9, 0.07, 0.0, 0.07, 0.9,  0.9,  0.07,  0.07/)
  !                              grass arabl  crops conif decid  water  urba  othr  desr  ice  sav   trf  wai   med  sem  wheat beech spruce

  real, parameter  :: temp_coeff_no(nlu_ozone_specials) = &
                                (/0.67, 0.67, 0.67, 0.84,  0.84,  0.0,  0.0,  0.0,  0.0,  0.0, 0.67, 0.84, 0.0, 0.84, 0.67, 0.67,  0.84, 0.84/)
  !                              grass arabl  crops conif decid  water  urba  othr  desr  ice  sav   trf   wai   med  sem   wheat beech  spruce
  real, parameter  :: base_coeff_no(nlu_ozone_specials) = &
                                (/8.80, 8.80, 8.80, 3.60,  3.60,  0.0,  0.0,  0.0,  0.0,  0.0, 8.80, 3.60, 0.0, 3.60, 8.80, 8.80,  3.60, 3.60/)


  ! ---data---

  ! fractions of each landuse type
  real, allocatable     ::  lu_fracs(:,:,:)    ! (nx,ny,nlu)

  ! vegetation fractions (tree types, common lu-types)
  real, allocatable     ::  veg_fracs(:,:,:)   ! (nx,ny,nveg)

  ! aerodynamic resistances
  ! Ra for heat (htop-->zobs, htop-->z0, htop-->zcanopy, zcanopy-->z0)
  real, allocatable   ::  Ra_h_htop_zobs_lu(:,:,:)       ! (nx,ny,nlu)
  real, allocatable   ::  Ra_h_htop_z0_lu(:,:,:)         ! (nx,ny,nlu)
  real, allocatable   ::  Ra_h_htop_zcanopytop_lu(:,:,:) ! (nx,ny,nlu)
  real, allocatable   ::  Ra_h_zcanopytop_z0_lu(:,:,:)   ! (nx,ny,nlu)

  ! friction velocities
  real, allocatable   ::  ustar_lu(:,:,:)                ! (nx,ny,nlu)
  real, allocatable   ::  u_zcanopytop_lu(:,:,:)         ! (nx,ny,nlu)
  ! For dust emissions roughness parameters on very local scale are used!!
  real, allocatable   ::  ustar_lu_dust_emis(:,:,:)      ! (nx,ny,nlu)


  ! growing season info
  real, allocatable   ::  start_gs(:,:,:)                ! (nx,ny,nlu)
  real, allocatable   ::  end_gs(:,:,:)                  ! (nx,ny,nlu)

  ! number of vegatation types  (200 tree-types + number of landusetypes)
  integer               ::  nveg
  integer               ::  nlu

  ! Seafraction
  real, allocatable     ::  waterbodies(:,:,:)

  ! use CLC 2006 landuse cover if available? (for deposition)
  logical               ::  use_clc2006_landuse
  ! soil texture data for dust
  logical               ::  with_soil_texture
  ! special ozone flux landuses
  logical               ::  with_ozone_specials
  ! deposition output on whole grid
  logical,allocatable   ::  ludep_output_whole_grid(:)

  ! growing season files
  character(len=256), allocatable ::  filenames_gs(:)    ! filenames for growing seasons
  integer                         ::  nlu_names_gs
  character(len=32)               ::  lu_names_gs(max_names_gs)
  ! current year (need for growing seaseon files)
  integer                         ::  current_year

contains

  subroutine LE_Landuse_Data_Alloc( status )

    use dims, only : nx,ny

    ! --- in/out ------------------------------

    integer, intent(out)    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Landuse_Data_Alloc'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! nuber of landuse categories
    if ( with_ozone_specials ) then
      nlu = nlu_ozone_specials
    else
      nlu = nlu_normal
    end if

    ! lu dependent output on whole grid
    allocate( ludep_output_whole_grid( nlu ) )
    ludep_output_whole_grid = .false.

    ! current year for grwoing season files
    current_year = -9999

    ! number of vegetations
    nveg = ntree_type + nlu

    ! landuse fractions
    allocate( lu_fracs(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( veg_fracs(nx,ny,nveg), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! Atmospheric resistance
    allocate( Ra_h_htop_zobs_lu(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Ra_h_htop_z0_lu(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Ra_h_htop_zcanopytop_lu(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( Ra_h_zcanopytop_z0_lu(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! friction velocities:
    allocate( ustar_lu(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( u_zcanopytop_lu(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( ustar_lu_dust_emis(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! start and end date of growing seasons
    allocate(  start_gs(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate(    end_gs(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! initialize with dummy values:
    start_gs = -999.
    end_gs   = -999.

    allocate( waterbodies(nx,ny,nwaterbody_type), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Landuse_Data_Alloc


  subroutine LE_Landuse_Data_Dealloc( status )

    ! --- in/out ------------------------------

    integer, intent(out)    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Landuse_Data_Dealloc'

    ! --- begin ----

    deallocate( lu_fracs, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( veg_fracs, stat=status )
    IF_NOTOK_RETURN(status=1)

    deallocate( Ra_h_htop_zobs_lu, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Ra_h_htop_z0_lu, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Ra_h_htop_zcanopytop_lu, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( Ra_h_zcanopytop_z0_lu, stat=status )
    IF_NOTOK_RETURN(status=1)

    deallocate( ustar_lu, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( u_zcanopytop_lu, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( ustar_lu_dust_emis, stat=status )
    IF_NOTOK_RETURN(status=1)

    deallocate(  start_gs, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(    end_gs, stat=status )
    IF_NOTOK_RETURN(status=1)

    deallocate( waterbodies, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Landuse_Data_Dealloc

end module LE_Landuse_Data
