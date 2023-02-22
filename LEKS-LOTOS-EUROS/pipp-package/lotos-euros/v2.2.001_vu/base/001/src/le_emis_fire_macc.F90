!###############################################################################
!
! LE_Emis_Fire_MACC
!   routines to read and work with daily and hourly fire emissions
!
!                                          Henk Eskes, KNMI, March 2010
!
! Emissions from the MACC project (jan 2010)
!   The MACC Global Fire Assimilation System: First Emission Products (GFASv0)
!   J.W. Kaiser, J. Flemming, M.G. Schultz, M. Suttie and M.J. Wooster
!   ECMWF Tech. Memo. 596, August 2009
!
!   Experiment "fagg"
!   http://www.gmes-atmosphere.eu/about/project_structure/input_data/d_fire/ProductsInMARS/
!   Available for 2003-2008, daily, 2D maps, 0.1 deg resol, redistributed GFED3
!   Species groups:
!   a = BC, CO2, CO, CH4, OC, PM2.5, SO2, TPM, C
!   b = H, NOx, N2O, NMHC
!   c = dry matter combustion rate, total carbon emission, C2H4, C2H4O (CH3CHO), C2H5OH,
!       C2H6, C2H6S (DMS), C3H6, C3H6O (CH3COCH3), C3H8, C5H8 (ISOP), CH2O, CH3OH, Higher_Alkanes,
!       Higher_Alkenes, NH3, Terpenes, Toluene_lump
!
! The original "fagg" data has been extracted from MARS,
!   regridded to 0.5x0.25 degree on the European GEMS/MACC RAQ domain
!   and stored in netcdf files as 2D fields
!   one field per day, 365 days in a file, one file per species
!   File name convention = <species id>_yyyy_le.nc
!      Example: "CO_2007_le.nc"
!   <species_id> =
!   CO  NO  NH3  SO2  DMS  CH2O CH4
!   C2H4  C10H16  C2H5OH  C2H6  C3H6  C3H8  CH3CHO  CH3COCH3  CH3OH
!      HIGH_ALKANES  HIGH_ALKENES  ISOP  TERP  TOLUENE
!   NMHC  TCFIRE
!   BCFIRE  OCFIRE  PM  PM2.5
!
! This module includes some code from TM5, Vincent Huijnen, KNMI:
! vh According to Angelika Heil (FZ Juelich) the mol. weights for lumped species
! vh    in biomass burning are different from the anthropogenic emissions, as follows:
! vh Ankelika Heil:
! We updated these molecular weights for biomass burning for the MOZART
! lumped groups as follows:

! Higher_Alkanes:  58 g/mole            ==> 78.8 g/mole
! Higher_Alkenes:  56 g/mole            ==> 64.0 g/mole
! Toluene_lump: 92 g/mole               ==> 85.7 g/mole

! The MOZART lumped groups were formed as follows:

! Toluene_lump (C7H8+ C6H6 + C8H10): MOZART lumped toluene species,
! incorporating benzene, toluene, xylene

! Higher Alkenes (CnH2n, C>=4): all alkenes (C>=4) specified in Andreae
! and Merlet (2001) which are not contained in the Toluene_lump group.
! These are: Butenes (1-butene + i-butene + tr-2-butene + cis-2-butene)
! (C4H8), Pentenes  (1-pentene + 2-pentene) (C5H10), Hexene  (C6H12),
! Octene  (C8H16)

! Higher Alkanes (CnH2n+2, C>=4): all alkanes (C>=4) specified in Andreae
! and Merlet (2001). These are: Butanes (n-butane + i-butane) (C4H10),
! Pentanes (n-pentane + i-pentane) (C5H12), Hexane (n-hexane + i-hexane)
! (C6H14), Heptane (C7H16)

! We calculated the molecular weight of these lumped groups from weighting
! the molecular mass of the individual species with their relative mass
! contribution to the total fire emissions of each group during 2001-2006
! (calculated from GFEDv2 data).
!
!
! 2011-11-03, Richard Kranenburg (TNO)
!   Changed routines such that fire emissions can be used for grids
!   other than MACC Emission grid. Emissions in grid cell from fire-files
!   distributed over underlying grid cell from run-domain
!
! 2012-04-11, Henk Eskes (KNMI)
!   Support files with list of fire sources as produced by IDL script.
!
! 2012-07-02, Arjo Segers (TNO)
!   Support files with compressed coordinates following CF-conventions
!   as produced by LEIP scripts.
!
! 2015-10, Arjo Segers (TNO)
!   Support CO2 emissions.
!   Bugs fixed:
!   - Only latest MACC/FIRE species that contributes to a
!     model CBM4 lumped species was taken into account, now changed.
!   - In "voc_to_cbm4_ole" an input molemass for "C3H6" of "56.0e-3"
!     was used, now changed to correct value.
!
! 2018-04, John Douros (KNMI)
!   Read and use hourly emissions files. Implementation is still a bit slopy, as 
!   it is based on a number of if-checks for the tres key.
!
!   NOTE: current implementation with hardcoded indices is tricky ..
!     In future better consider:
!      - read macc-to-model species mapping from csv file
!      - use loops over all species instead of subsets
!      - use specunits to convert to target units
!     See template settings in:
!       proj/macc3-co2/005/rc/lotos-euros-emissions-gfas.rc
!
! 2019-12, AMM
!   Fixed bug: emitted total PM is not just ppm, but also includes EC and POM;
!   corrected assignment to LE tracers.
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

module LE_Emis_Fire_MACC

  use GO, only : gol, goPr, goErr
  use GO, only : TDate
  use binas, only : xm_C, xm_O, xm_H
  use Indices, only : ispec_co, ispec_no, ispec_no2, ispec_nh3, ispec_so2, ispec_ch4, ispec_par, &
                      ispec_eth, ispec_ole, ispec_ald, ispec_ald2, ispec_form, ispec_iso, ispec_terp, &
                      ispec_co2

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  LE_Emis_Fire_MACC_Init, LE_Emis_Fire_MACC_Done
  public  ::  LE_Emis_Fire_MACC_Setup

  ! --- const ------------------------------------

  ! module name:
  character(len=*), parameter ::  mname = 'LE_Emis_Fire_MACC'

  ! Number of species names / emission files
  integer, parameter          ::  nspec_macc = 23
  ! which part of this is gas
  integer, parameter          ::  nspec_macc_gas = 19  ! CO .. CH3COCH3, CO2
  ! how many NMVOC species
  integer, parameter          ::  nvoc = 13
  ! MACC/GFED3 species names (used to construct file names)
  character(len=12), parameter :: speciesNCFNames(nspec_macc) =  (/ &
       'CO          ', 'NO          ', 'NH3         ', 'SO2         ', 'CH4         ', &  ! 1-5  , inorganic+ch4
       'C2H5OH      ', 'C2H6        ', 'C3H8        ', 'HIGH_ALKANES',                 &  ! 6-9  , nmvoc
       'C2H4        ', 'C3H6        ', 'ISOP        ', 'TERP        ',                 &  ! 10-13, nmvoc
       'HIGH_ALKENES', 'TOL         ', 'CH2O        ', 'CH3CHO      ',                 &  ! 14-17, nmvoc
       'CH3COCH3    ',                                                                 &  ! 18   , nmvoc
       'CO2         ',                                                                 &  ! 19     CO2
       'BCFIRE      ', 'OCFIRE      ', 'PM          ', 'PM2.5       '               /)    ! 20-23  aerosol

  ! CBM4 gas species for which fire emissions are defined
  integer, parameter             ::  nspec_cbm4 = 14
  ! indices of target trace gas species;
  ! switch from ald (cbm4) to ald2 (cb99) if necessary:
  integer, parameter  ::  cb4_spec(nspec_cbm4) = &
    (/  ispec_co  , ispec_no  , ispec_no2 , ispec_nh3 , ispec_so2 , ispec_ch4 ,  ispec_par , &
        ispec_eth , ispec_ole , max(ispec_ald,ispec_ald2) , ispec_form, ispec_iso, ispec_terp,&
        ispec_co2   /)
  character(4), parameter :: cb4_spec_name(nspec_cbm4) = (/ &
       'CO  ', 'NO  ', 'NO2 ', 'NH3 ', 'SO2 ', 'CH4 ',         & !  1 - 6
       'PAR ', 'ETH ', 'OLE ', 'ALD ', 'FORM', 'ISOP', 'TERP', & !  7 - 13
       'CO2 '    /)                                              ! 14

  ! mol mass of several hydrocarbons
  real, parameter     ::  xm_C2H6O  = xm_C*2  + xm_H*6  + xm_O
  real, parameter     ::  xm_C2H6   = xm_C*2  + xm_H*6
  real, parameter     ::  xm_C3H8   = xm_C*3  + xm_H*8
  real, parameter     ::  xm_C3H6   = xm_C*3  + xm_H*6
  real, parameter     ::  xm_C2H4   = xm_C*2  + xm_H*4
  real, parameter     ::  xm_ISOP   = xm_C*5  + xm_H*8        ! C5H8
  real, parameter     ::  xm_TERP   = xm_C*10 + xm_H*18       ! C10H18
  real, parameter     ::  xm_CH2O   = xm_C    + xm_H*2  + xm_O
  real, parameter     ::  xm_CH3COCH3 = xm_C + xm_H*3 + xm_C + xm_O + xm_C + xm_H*3

  ! PAR
  !                           mol/kg    (ispec_macc)
  real, parameter     ::  voc_to_cbm4_par(6:nvoc+5) = (/ &
                             1.0/xm_C2H6O   , &  !  6 - C2H5OH
                             1.0/xm_C2H6    , &  !  7 - C2H6
                             1.5/xm_C3H8    , &  !  8 - C3H8
                             4.0/78.8E-3    , &  !  9 - HIGH_ALKANES  (Heil MOZART)
                             0.0            , &  ! 10 - C2H4
                             1.0/xm_C3H6    , &  ! 11 - C3H6
                             0.0            , &  ! 12 - ISOP
                             0.0            , &  ! 13 - TERP
                             0.0            , &  ! 14 - HIGH_ALKENES
                             1.0/85.7e-3    , &  ! 15 - TOL            (Heil MOZART)
                             0.0            , &  ! 16 - CH2O
                             0.0            , &  ! 17 - CH3CHO
                             3.0/xm_CH3COCH3 /)  ! 18 - CH3COCH3
  ! OLE
  !                            mol/kg    (ispec_macc)
  real, parameter     ::  voc_to_cbm4_ole(6:nvoc+5) = (/ &
                             0.0          , &  !  6 - C2H5OH
                             0.0          , &  !  7 - C2H6
                             0.0          , &  !  8 - C3H8
                             0.0          , &  !  9 - HIGH_ALKANES
                             0.0          , &  ! 10 - C2H4
                             1.0/xm_c3h6  , &  ! 11 - C3H6             BUG FIXED: was 1.0/56.0e-3
                             0.0          , &  ! 12 - ISOP
                             0.0          , &  ! 13 - TERP
                             1.0/64.0e-3  , &  ! 14 - HIGH_ALKENES     (Heil MOZART)
                             0.0          , &  ! 15 - TOL
                             0.0          , &  ! 16 - CH2O
                             0.0          , &  ! 17 - CH3CHO
                             0.0           /)  ! 18 - CH3COCH3
  ! ALD (Aceetaldehyde of ethanal)
  !                            mol/kg    (ispec_macc)
  real, parameter     ::  voc_to_cbm4_ald(6:nvoc+5) = (/ &
                             0.0          , &   !  6 - C2H5OH
                             0.0          , &   !  7 - C2H6
                             0.0          , &   !  8 - C3H8
                             0.0          , &   !  9 - HIGH_ALKANES
                             0.0          , &   ! 10 - C2H4
                             0.0          , &   ! 11 - C3H6
                             0.0          , &   ! 12 - ISOP
                             0.0          , &   ! 13 - TERP
                             0.0          , &   ! 14 - HIGH_ALKENES
                             0.0          , &   ! 15 - TOL
                             0.0          , &   ! 16 - CH2O
                             1.0/44.0e-3  , &   ! 17 - CH3CHO   = 2 C + 4 H + O
                             0.0           /)   ! 18 - CH3COCH3       

  ! Add fire emissions to the boundary layer, level 2
  integer, parameter ::  emis_level = 2

  ! fraction NO from NOx   from  Andrea and Merlet 2001
  real, parameter    ::  NO_rfire = 0.9

  ! conversion factors:
  real, parameter    ::  ug_per_kg   = 1.0e9   ! ug/kg
  real, parameter    ::  sec_per_min = 60.0    ! sec/min

  ! --- var --------------------------------------

  ! total emission cocktail:
  ! real, allocatable  ::  em_fire(:,:)

  ! location of emissions files:
  character(len=1024)  ::  emispath_fire

  ! format of emissions files:
  character(len=32)    ::  emis_maccfire_tres

  ! store emissions:
  real, allocatable    ::  emis_maccfire(:,:,:)
  real, allocatable    ::  emis_maccfire_hourly(:,:,:,:)

  logical              ::  firstcall = .true.

  ! allow missing files ?
  logical              ::  arch_allow_missing


contains


  ! ===============================================================


  subroutine LE_Emis_Fire_MACC_Init( rcF, rckey, status )

    use GO, only : TrcFile, ReadRc
    use Dims, only : nx, ny, nspec

    use LE_Data, only : LE_Data_Enable

    implicit none

    ! --- in/out ------------------------------

    type(TrcFile), intent(in)       ::  rcF
    character(len=*), intent(in)    ::  rckey
    integer, intent(out)            ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Fire_MACC_Init'

    ! --- local ------------------------------

    ! --- begin -------------------------------

    ! read settings:
    call ReadRc( rcF, trim(rckey)//'.files', emispath_fire, status )
    IF_NOTOK_RETURN(status=1)

    ! Temporal resolution, hourly or daily file
    call ReadRc( rcF, trim(rckey)//'.tres', emis_maccfire_tres, status )
    IF_NOTOK_RETURN(status=1)

    ! allow missing files ?
    call ReadRc( rcF, trim(rckey)//'.allow_missing', arch_allow_missing, status )
    IF_NOTOK_RETURN(status=1)

    ! storage for emissions:
    if ( emis_maccfire_tres == 'hourly' ) then
      allocate(  emis_maccfire_hourly(nx,ny,0:23,nspec) )
    elseif ( emis_maccfire_tres == 'daily' ) then
      allocate(  emis_maccfire(nx,ny,nspec) )
    else
       write (gol,'("unsupported macc temporal resolution format `",a,"`")') trim(emis_maccfire_tres); call goErr
       TRACEBACK; status=1; return
    endif

    ! enable data
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Fire_MACC_Init


  ! ***


  subroutine LE_Emis_Fire_MACC_Done( status )

    implicit none

    ! --- in/out ------------------------------

    integer, intent(out)  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Fire_MACC_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! clear:
    if (allocated( emis_maccfire_hourly ) ) deallocate(  emis_maccfire_hourly )
    if (allocated( emis_maccfire ) ) deallocate(  emis_maccfire )

    ! ok
    status = 0

  end subroutine LE_Emis_Fire_MACC_Done


  ! ***


  subroutine LE_Emis_Fire_MACC_Setup( emis_a, t, status )

    use GO     , only : TDate, Midnight
    use Grid   , only : AreaOper
    use LE_Grid, only : ugg
    use Dims   , only : runF
    use Dims   , only : nx, ny, nz, nspec
    use indices

    ! --- in/out ---------------------------

    real, intent(inout)       ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter :: rname = mname//'/LE_Emis_Fire_MACC_Setup'
    real, parameter             :: tiny = 1.0e-30

    ! multiplication factors GFED3/MACC -> CBM4
    real                        ::  macc_to_cbm4(nspec_macc,nspec_cbm4)

    ! --- local ----------------------------

    real, dimension(:,:), allocatable :: emis_ncf
    real, dimension(:,:,:), allocatable :: emis_ncf_hourly

    integer                           :: i, j, hh, i_macc, i_cbm4, i_max, j_max, hour_max
    real                              :: weight, max_emis
    character(20)                     :: speciesIdString
    character(20)                     :: units
    integer                           :: ispec, ispec1, ispec2, ispec3

    ! --- begin ----------------------------

    ! info ...
     write (gol,'("          MACC forest fires ...")'); call goPr


    ! read new emission data from file
    if ( firstcall .or. Midnight( t ) ) then

       ! create storage
       ! initialise fire emissions
       if ( emis_maccfire_tres == 'hourly' ) then
          emis_maccfire_hourly(:,:,:,:) = 0.0
          allocate ( emis_ncf_hourly(ugg%nlon,ugg%nlat,0:23) )
       elseif ( emis_maccfire_tres == 'daily' ) then
          emis_maccfire(:,:,:) = 0.0
          allocate ( emis_ncf(ugg%nlon,ugg%nlat) )
       else
          write (gol,'("unsupported macc temporal resolution format `",a,"`")') trim(emis_maccfire_tres); call goErr
          TRACEBACK; status=1; return
       endif

       ! set up the MACC -> CBM4 conversion table
       ! macc_to_cbm4(nspec_macc,nspec_cbm4)
       !   rows 1-19   :  mol/kg
       !   other       :  undefined
       macc_to_cbm4(:,:) = 0.0

       ! fill matrix with 1/molarweight = mol/kg factors
       do ispec = 1, nspec
         select case ( ispec )
           !                     mol/kg          =   1 /  (kg/mol)                       ! macc   -> cb4
           case (ispec_co ); macc_to_cbm4( 1, 1) = 1.0 / specmolm(ispec)                 !  1 CO  ->  1 CO
           case (ispec_no ); macc_to_cbm4( 2, 2) =        NO_rfire  / specmolm(ispec)    !  2 NO  ->  2 NO
           case (ispec_no2); macc_to_cbm4( 2, 3) = (1.0 - NO_rfire) / specmolm(ispec)    !  2 NO  ->  3 NO2
           case (ispec_nh3); macc_to_cbm4( 3, 4) = 1.0 / specmolm(ispec)                 !  3 NH3 ->  4 NH3
           case (ispec_so2); macc_to_cbm4( 4, 5) = 1.0 / specmolm(ispec)                 !  4 SO2 ->  5 SO2
           case (ispec_ch4); macc_to_cbm4( 5, 6) = 1.0 / specmolm(ispec)                 !  5 CH4 ->  6 CH4
           case (ispec_co2); macc_to_cbm4(19,14) = 1.0 / specmolm(ispec)                 ! 19 CO2 -> 14 CO2
         end select
       end do
       ! VOC part                                                  macc                 -> cb4
       do i = 6, nvoc+5
          macc_to_cbm4(i,7) = voc_to_cbm4_par(i)              ! 6 C2H5OH .. 18 CH3COCH3 ->  7 PAR
       end do
       macc_to_cbm4(10,8) = 1.0 / xm_c2h4                     ! 10 C2H4                 ->  8 ETH
       do i = 6, nvoc+5
          macc_to_cbm4(i,9) = voc_to_cbm4_ole(i)              ! 6 C2H5OH .. 18 CH3COCH3 ->  9 OLE
       end do
       do i = 6, nvoc+5
          macc_to_cbm4(i,10) = voc_to_cbm4_ald(i)             ! 6 C2H5OH .. 18 CH3COCH3 -> 10 ALD
       end do
       macc_to_cbm4(16,11) = 1.0 / xm_ch2o                    ! 16 Ch2O                 -> 11 FORM
       macc_to_cbm4(12,12) = 1.0 / xm_isop                    ! 12 ISOP                 -> 12 ISO
       macc_to_cbm4(13,13) = 1.0 / xm_terp                    ! 13 TERP                 -> 13 TERP

       ! initialise fire emissions
       emis_maccfire(:,:,:) = 0.0

       ! read the fire emissions from netcdf file
       ! write (gol,'("LE:       Read MACC fire emissions for new day ...")'); call goPr

       ! loop over species / emission files
       spec_loop: do i_macc = 1, nspec_macc

          ! current species:
          speciesIdString = trim(speciesNCFNames(i_macc))
          ! target units:
          units = 'kg m**-2 s**-1'

          ! read from file with list of fire locations and strengths
          !
          if ( emis_maccfire_tres == 'hourly' ) then

             !! info ...
             !write (gol,'("          Hourly MACC forest fires ...")'); call goPr
             ! read:
             call maccFire_nc_readCompress_hourly( emispath_fire, speciesIdString, units, t, &
                                              i_macc==1, ugg%nlon,ugg%nlat, emis_ncf_hourly, status )
             IF_NOTOK_RETURN(status=1)

          else if ( emis_maccfire_tres == 'daily' ) then

             !! info ...
             !write (gol,'("          Daily MACC forest fires ...")'); call goPr
             ! read:
             call maccFire_nc_readCompress( emispath_fire, speciesIdString, units, t, &
                                              i_macc==1, ugg%nlon,ugg%nlat, emis_ncf, status )
             IF_NOTOK_RETURN(status=1)

          else
             write (gol,'("unsupported macc temporal resolution format `",a,"`")') trim(emis_maccfire_tres); call goErr
             TRACEBACK; status=1; return
          endif

          ! find location of max emissions
          if ( emis_maccfire_tres == 'hourly' ) then
             max_emis = -1.0
             i_max = 1
             j_max = 1
             hour_max =0
             do i = 1, ugg%nlon
                do j = 1, ugg%nlat
                   do hh = 0, 23
                      if ( emis_ncf_hourly(i,j,hh) > max_emis ) then
                         max_emis = emis_ncf_hourly(i,j,hh)
                         i_max = i
                         j_max = j
                         hour_max = hh
                      end if
                   end do !hh
                end do !j
             end do !i
             if ( i_macc == 1 ) then
               write (gol,'("            Maximum fire emissions found in grid cell  i, j =",2i4)') i_max, j_max; call goPr
               write (gol,'("            at hour = ",i4)') hour_max; call goPr
             end if
          elseif ( emis_maccfire_tres == 'daily' ) then
             max_emis = -1.0
             i_max = 1
             j_max = 1
             do i = 1, ugg%nlon
                do j = 1, ugg%nlat
                   if ( emis_ncf(i,j) > max_emis ) then
                      max_emis = emis_ncf(i,j)
                      i_max = i
                      j_max = j
                   end if
                end do
             end do
             if ( i_macc == 1 ) then
               write (gol,'("            Maximum fire emissions found in grid cell  i, j =",2i4)') i_max, j_max; call goPr
             end if
          else
             write (gol,'("unsupported macc temporal resolution format `",a,"`")') trim(emis_maccfire_tres); call goErr
             TRACEBACK; status=1; return
          endif


          ! store emissions

          if ( i_macc <= nspec_macc_gas ) then

             ! Trace gas fire emissions - store in "emis_maccfire"
             ! gas species: convert from (kg m-2 s-1) to (mol min-1)
             !       Multiply by: gridCellArea(1:nx,1:ny) * sec_per_min / specmolw
             ! (area multiplication after all contribution have been added)
             
             ! loop over cbm4 species created from macc/fire species:
             do i_cbm4 = 1, nspec_cbm4
                ! target tracer; not enabled ?
                ispec = cb4_spec(i_cbm4)
                if ( ispec <= 0 ) cycle
                ! does current macc/fire emission contribut to this cbm4 species in model?
                if ( macc_to_cbm4(i_macc,i_cbm4) > tiny ) then
                   ! contribution weight:
                   !            s/min          mol/kg
                   weight = sec_per_min * macc_to_cbm4(i_macc,i_cbm4)    ! mol/kg s/min
                   ! info ...
                   !print '(11x,a12,a,a4,a,e11.4,a,e11.4)', &
                   !     speciesNCFNames(i_macc),'-> ',specname(cb4_spec(i_cbm4)),', weight = ',weight,', max emis = ',emis_ncf(i_max,j_max) * gridCellArea(i_max,j_max)
                   !>>> BUG: This ignores multiple contributions to same ispec ...
                   !    BUG: AJS, 2015-10
                   !emis_maccfire(:,:,ispec) = emis_ncf(:,:) * weight
                   !! total in grid cell:
                   !call AreaOper( lli, emis_maccfire(:,:,ispec), '*', 'm2', status )
                   !IF_NOTOK_RETURN(status=1)
                   !<<<
                   ! add contribution:
                   !     mol/m2/min                                      kg/m2/s  mol/kg s/min
                   !
                   if ( emis_maccfire_tres == 'hourly' ) then
                      emis_maccfire_hourly(:,:,:,ispec) = emis_maccfire_hourly(:,:,:,ispec) + emis_ncf_hourly * weight
                   elseif ( emis_maccfire_tres == 'daily' ) then
                      emis_maccfire(:,:,ispec) = emis_maccfire(:,:,ispec) + emis_ncf * weight
                   else
                      write (gol,'("unsupported macc temporal resolution format `",a,"`")') trim(emis_maccfire_tres); call goErr
                      TRACEBACK; status=1; return
                   endif
                   !
                   !
                   !
                end if
             end do

          else

             ! aerosol fire emissions
             ! PM species: convert from (kg m-2 s-1) to (ug min-1)
             !       Multiply by: gridCellArea(1:ny) * ug_per_kg * sec_per_min
             ! (area multiplication after all contribution have been added)

             ! conversion factor:
             !          ug/kg         s/min
             weight = ug_per_kg * sec_per_min   ! ug/kg s/min
             ! model tracer index:
             select case ( trim(speciesNCFNames(i_macc)) )
               case ( 'BCFIRE' ); ispec = i_ec_f
               case ( 'OCFIRE' ); ispec = i_pom_f
               case ( 'PM'     ); ispec = i_ppm_c
               case ( 'PM2.5'  ); ispec = i_ppm_f
               case default
                write (gol,'("could not assign MACC fire variable `",a,"` to tracer")') trim(speciesNCFNames(i_macc)); call goErr
                TRACEBACK; status=1; return
             end select

             ! enabled ?
             if ( ispec > 0 ) then
               ! info ...
               !print '(11x,a12,a,a4,a,e11.4,a,e11.4)', &
               !     speciesNCFNames(i_macc),'-> ',specname(ispec),', weight = ',weight,', max emis = ',emis_ncf(i_max,j_max) * gridCellArea(i_max,j_max)
               !! info ...
               !write (gol,'("       MACC/Fire: add `",a,"` to `",a,"`")') &
               !        trim(speciesNCFNames(i_macc)), trim(specname(ispec)); call goPr
               !>>> BUG: This ignores multiple contributions to same ispec ...
               !    BUG: AJS, 2015-10
               !! add:
               !emis_maccfire(:,:,ispec) = emis_ncf(:,:) * weight
               !! total in grid cell:
               !call AreaOper( lli, emis_maccfire(:,:,ispec), '*', 'm2', status )
               !IF_NOTOK_RETURN(status=1)
               !<<<
               ! add:
               !      ug/m2/min                                      kg/m2/s    ug/kg s/min
               if ( emis_maccfire_tres == 'hourly' ) then
                  emis_maccfire_hourly(:,:,:,ispec) = emis_maccfire_hourly(:,:,:,ispec) + emis_ncf_hourly * weight
               elseif ( emis_maccfire_tres == 'daily' ) then
                  emis_maccfire(:,:,ispec) = emis_maccfire(:,:,ispec) + emis_ncf * weight
               else
                  write (gol,'("unsupported macc temporal resolution format `",a,"`")') trim(emis_maccfire_tres); call goErr
                  TRACEBACK; status=1; return
               endif
             end if

          end if   ! gas, else aerosol

       end do spec_loop

       ! correct for pm2.5 included in pm10
       ! PM read from file is total PM, including OC and BC

       if ( emis_maccfire_tres == 'hourly' ) then
          if ( (i_ppm_c > 0) .and. (i_ppm_f > 0) ) then
            ispec  = i_ppm_c
            ispec1 = i_ppm_f
            do hh=0,23
               ! "ppm coarse" = "pm total" - "pm fine"
               emis_maccfire_hourly(:,:,hh,ispec) = emis_maccfire_hourly(:,:,hh,ispec) - emis_maccfire_hourly(:,:,hh,ispec1)
               ! ec and pom tracers in model ?
               if ( (i_ec_f > 0) .and. (i_pom_f > 0) ) then
                  ispec2=i_ec_f
                  ispec3=i_pom_f
                  ! "ppm fine" := "ppm fine" - "ec fine" - "pom fine"
                  emis_maccfire_hourly(:,:,hh,ispec1) = emis_maccfire_hourly(:,:,hh,ispec1) - emis_maccfire_hourly(:,:,hh,ispec2)- emis_maccfire_hourly(:,:,hh,ispec3)
               !else
               !  ! only total pm relevant, do nothing to separate
               end if 
            end do
          end if
       elseif ( emis_maccfire_tres == 'daily' ) then
          if ( (i_ppm_c > 0) .and. (i_ppm_f > 0) ) then
            ispec  = i_ppm_c
            ispec1 = i_ppm_f
            emis_maccfire(:,:,ispec) = emis_maccfire(:,:,ispec) - emis_maccfire(:,:,ispec1)
            ! "ppm coarse" = "pm total" - "pm fine"
            if ( (i_ec_f > 0) .and. (i_pom_f > 0) ) then
                ispec2=i_ec_f
                ispec3=i_pom_f
                  ! "ppm fine" := "ppm fine" - "ec fine" - "pom fine"
                emis_maccfire(:,:,ispec1) = emis_maccfire(:,:,ispec1) - emis_maccfire(:,:,ispec2)- emis_maccfire(:,:,ispec3)
            !else
            !  ! only total pm relevant, do nothing to separate
            end if
          end if
       else
          write (gol,'("unsupported macc temporal resolution format `",a,"`")') trim(emis_maccfire_tres); call goErr
          TRACEBACK; status=1; return
       endif

       ! convert from  (mol|ug)/m2/min  to  (mol|ug)/min  in grid cell:
       call ugg%AreaOper( emis_maccfire, '*', 'm2', status )
       IF_NOTOK_RETURN(status=1)

       ! clear:
       if (allocated( emis_ncf ) ) deallocate ( emis_ncf )
       if (allocated( emis_ncf_hourly ) ) deallocate ( emis_ncf_hourly )

       ! reset flag:
       firstcall = .false.

    end if   ! firstcall or midnight - read new emission data from file

    !-----------------------------------------------------------------
    ! add emissions to 'emis_a' for hour 00..23 with a diurnal profile
    !-----------------------------------------------------------------

    if ( emis_maccfire_tres == 'hourly' ) then
        !emissions level is equal to 2 
        !TODO: use injection height to calculate injection model level
        !emis_a(:,:,2,:) = emis_maccfire_hourly(:,:,t%hour+1,:) 
        emis_a(:,:,2,:) = emis_maccfire_hourly(:,:,t%hour,:) 
    elseif ( emis_maccfire_tres == 'daily' ) then
       call make_fire_emis( emis_a, t%hour+1, status )
       IF_NOTOK_RETURN(status=1)
    else
       write (gol,'("unsupported macc temporal resolution format `",a,"`")') trim(emis_maccfire_tres); call goErr
       TRACEBACK; status=1; return
    endif

    ! ok
    status = 0

  end subroutine LE_Emis_Fire_MACC_Setup


  ! ***

  subroutine make_fire_emis( emis_a, hour, status )
    ! add emissions to 'emis_a' for hour 01..24

    use dims,    only : nx, ny, nz, nspec
    use indices, only : i_ec_f, i_ppm_c, i_ppm_f, i_pom_f

    implicit none

    ! --- in/out ---------------------------------

    integer, intent(in)         ::  hour
    real, intent(inout)         ::  emis_a(nx,ny,nz,nspec)
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter :: rname = mname//'/make_fire_emis'

    ! add to level 2 only ( = boundary layer )
    integer, parameter          ::  emis_level = 2

    ! Daily cycle of fire emissions:
    !
    ! Profile introduced for SILAM fire emissions (Arjo Segers)
    !real, dimension(1:24), parameter ::  hour_frac_ffires = &
    !     (/ 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, &
    !        1.25, 1.25, 1.25, 1.25, 1.25, 1.25, &
    !        1.25, 1.25, 1.25, 1.25, 1.25, 1.25, &
    !        0.75, 0.75, 0.75, 0.75, 0.75, 0.75 /)
    !
    ! Profile from Johannes Kaiser
    ! ECMWF Tech. Memo. 596, August 2009, p4
    ! ( May be suitable for the tropics, but probably too little
    !   night fires for Europe )
    !real, dimension(1:24), parameter ::  hour_frac_Kaiser = &
    !     (/ 0.2000, 0.2000, 0.2000, 0.2000, 0.2005, 0.2034, &
    !        0.2195, 0.2873, 0.5047, 1.0283, 1.9534, 3.0909, &
    !        3.9120, 3.9120, 3.0909, 1.9534, 1.0283, 0.5047, &
    !        0.2873, 0.2195, 0.2034, 0.2005, 0.2000, 0.2000 /)
    !
    ! The following form corresponds roughly to the results found for the US by
    ! Zhang et al, RSE 2008, doi: 10.1016/j.rse.2008.02.006, Fig.10
    real, dimension(1:24), parameter ::  hour_frac_Zhang = &
         (/ 0.4000, 0.4000, 0.4000, 0.4000, 0.4003, 0.4025, &
            0.4146, 0.4655, 0.6285, 1.0212, 1.7151, 2.5682, &
            3.1840, 3.1840, 2.5682, 1.7151, 1.0212, 0.6285, &
            0.4655, 0.4146, 0.4025, 0.4003, 0.4000, 0.4000 /)

    ! --- local ------------------------------------

    integer  :: ispec
!    integer  :: i_cbm4
    real     ::  delta_emis_a(nx,ny)

    ! --- start ------------------------------------

    ! loop over target species:
    do ispec = 1, nspec
      ! current emissions:
      delta_emis_a = emis_maccfire(:,:,ispec) * hour_frac_Zhang(hour)
      ! any emissions?
      if ( any( delta_emis_a > 0.0 ) ) then
        ! add:
        emis_a(:,:,emis_level,ispec) = emis_a(:,:,emis_level,ispec) + delta_emis_a
      end if  ! any emis?
    end do  ! spec


    ! ok
    status = 0

  end subroutine make_fire_emis


  ! ***




  ! ---------------------------------------------------------------------------------------------
  !
  ! Read fire emissions as a list of point sources.
  ! From NetCDF files with the following content:
  !
  !  netcdf FIRE_20091201_1d {
  !  dimensions:
  !          lon = 1300 ;
  !          lat = 520 ;
  !          time = UNLIMITED ; // (1 currently)
  !          nv = 2 ;
  !          point = 45 ;
  !
  !  variables:
  !
  !          int point(point) ;
  !                  point:long_name = "zero-based indices in 1D array" ;
  !                  point:compress = "lat lon" ;
  !                  point:description = "original zero-based indices: lat_index = point/len(lon), lon_index = point mod len(lon)" ;
  !
  !          double lon(lon) ;
  !                  lon:long_name = "longitude" ;
  !                  lon:units = "degrees_east" ;
  !                  lon:standard_name = "longitude" ;
  !                  lon:axis = "X" ;
  !                  lon:bounds = "lon_bounds" ;
  !          double lon_bounds(lon, nv) ;
  !                  lon_bounds:long_name = "longitude" ;
  !                  lon_bounds:units = "degrees_east" ;
  !                  lon_bounds:standard_name = "longitude" ;
  !
  !          double lat(lat) ;
  !                  lat:long_name = "latitude" ;
  !                  lat:units = "degrees_north" ;
  !                  lat:standard_name = "latitude" ;
  !                  lat:axis = "Y" ;
  !                  lat:bounds = "lat_bounds" ;
  !          double lat_bounds(lat, nv) ;
  !                  lat_bounds:long_name = "latitude" ;
  !                  lat_bounds:units = "degrees_north" ;
  !                  lat_bounds:standard_name = "latitude" ;
  !
  !          double time(time) ;
  !                  time:units = "hours since 2009-12-01 12:00:00" ;
  !                  time:calendar = "proleptic_gregorian" ;
  !                  time:standard_name = "time" ;
  !                  time:axis = "T" ;
  !          float co2fire(time, point) ;
  !                  co2fire:table = 210 ;
  !                  co2fire:long_name = "Wildfire flux of Carbon Dioxide" ;
  !                  co2fire:units = "kg m**-2 s**-1" ;
  !                  co2fire:code = 80 ;
  !          ...
  !
  !  // global attributes:
  !                  :CDI = "Climate Data Interface version 1.5.0 (http://code.zmaw.de/projects/cdi)" ;
  !                  :Conventions = "CF-1.0" ;
  !                  :history = "Mon Jul 02 09:12:06 2012: gathered non-zero points ..."
  !                  :institution = "European Centre for Medium-Range Weather Forecasts" ;
  !                  :mars_keywords = "retrieve;class=rd;expver=fl6z;stream=oper;type=fc;..."
  !                  :CDO = "Climate Data Operators version 1.5.0 (http://code.zmaw.de/projects/cdo)" ;
  !  }
  !
  ! status on return:
  !    0 : OK, field read from file
  !   -1 : Emissions for species "speciesIdString" not found
  !    1 : Error
  !
  ! Arjo Segers, 2012-07-02
  ! ---------------------------------------------------------------------------------------------

  subroutine maccFire_nc_ReadCompress( path, speciesIdString, units, t, &
                                         verbose, nx, ny, fire_emis, status )

    use NetCDF, only : NF90_NoWrite, NF90_NOERR
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_StrError

    use GO     , only : TDate, NewDate, Get, Compare_Date_Num
    use GO     , only : operator(-), operator(+), operator(<), operator(<=)
    use GO     , only : IsAnyDate, IncrDate, wrtgol
    use GO     , only : goReplace, goLoCase, goSplitString
    use Grid   , only : ll_area_deg_m2
    use LE_Grid, only : ugg
    use Dims   , only : runF

    ! --- in/out --------------------------------

    character(len=*), intent(in)            ::  path
    character(len=*), intent(in)            ::  speciesIdString
    character(len=*), intent(in)            ::  units
    type(TDate), intent(in)                 ::  t
    logical, intent(in)                     ::  verbose
    integer, intent(in)                     ::  nx
    integer, intent(in)                     ::  ny
    real, dimension(nx,ny), intent(out)     ::  fire_emis
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/maccFire_nc_ReadCompress'

    ! max number of file templates:
    integer, parameter    ::  maxfile = 10

    ! --- local ---------------------------------

    character(len=1024)              ::  fnames(maxfile)
    integer                          ::  nfile, ifile
    character(len=1024)              ::  fname
    logical                          ::  exist
    integer                          ::  ncid, dimid, varid
    character(len=32)                ::  vname
    integer                          ::  l
    integer                          ::  year, month, day
    integer                          ::  itime
    integer                          ::  ntime
    real, allocatable                ::  times(:)
    character(len=32)                ::  time_units
    character(len=32)                ::  time_calendar
    type(TDate)                      ::  time_curr
    integer                          ::  npoint
    integer                          ::  ipoint
    integer, allocatable             ::  points(:)
    integer                          ::  nlon, nlat
    integer                          ::  ilon, ilat
    real, allocatable                ::  lon_bounds(:,:), lat_bounds(:,:)
    real, allocatable                ::  emis(:)
    character(len=32)                ::  emis_units
    integer                          ::  n_in_domain
    integer                          ::  i1, i2, j1, j2
    integer                          ::  i, j
    real, pointer                    ::  fracs(:,:)
    real                             ::  west, east, south, north
    real                             ::  area

    ! --- begin -------------------------------

    ! Initialise
    fire_emis(:,:) = 0.0

    ! forecast base specified?
    if ( IsAnyDate(runF%t_base) ) then
      ! target time:
      call Get( t, year=year, month=month, day=day )
    else
      ! before forecast base ?
      if ( t < runF%t_base ) then
        ! as usual:
        call Get( t, year=year, month=month, day=day )
        !
      ! within 2 days after forecast base?
      else if ( t <= runF%t_base+IncrDate(day=2) ) then
        ! info ...
        write (gol,'("WARNING - no fire emission forecast, use data from day before forecast base!")'); call goPr
        ! use analysis day:
        call Get( runF%t_base-IncrDate(day=1), year=year, month=month, day=day )
        !
      ! no data ..
      else
        ! info ...
        write (gol,'("WARNING - no fire emission forecast, keep zero at this day!")'); call goPr
        ! ok
        status=0; return
      end if
    end if

    ! path might consist of several values,
    ! should be tested until the first match:
    call goSplitString( trim(path), nfile, fnames, status )
    IF_NOTOK_RETURN(status=1)
    ! dummy:
    fname = '/no/file/found/yet'
    ! loop over directories:
    do ifile = 1, nfile
      ! replace some values:
      call goReplace( fnames(ifile), '%{year}', '(i4.4)', year, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( fnames(ifile), '%{month}', '(i2.2)', month, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( fnames(ifile), '%{day}', '(i2.2)', day, status )
      IF_NOTOK_RETURN(status=1)
      ! check, if present then store match and leave:
      inquire( file=trim(fnames(ifile)), exist=exist )
      if ( exist ) then
        fname = trim(fnames(ifile))
        exit
      end if
    end do
    ! check ..
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      if ( arch_allow_missing) then
        write (gol,'("WARNING - missing MACC fire emissions file; continue ...")'); call goPr
        status=0; return
      else
        write (gol,'("None of the specified MACC fire emission files found:")'); call goErr
        do ifile = 1, nfile
          write (gol,'("  ",a)') trim(fnames(ifile)); call goErr
        end do
        TRACEBACK; status=1; return
      end if
    end if

    ! info:
    if ( verbose ) then
      write (gol,'("         Read MACC fire emissions for ",i4, 2i3)') year, month, day; call goPr
      write (gol,'("           from file : ",a)') trim(fname); call goPr
    end if

    ! open netcdf file:
    status = NF90_Open( trim(fname), NF90_NoWrite, ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! variable name; first translate species name:
    select case ( trim(speciesIdString) )
      case ( 'NO'           ) ; vname = 'nox'
      case ( 'HIGH_ALKANES' ) ; vname = 'hialkanes'
      case ( 'ISOP'         ) ; vname = 'c5h8'
      case ( 'TERP'         ) ; vname = 'terpenes'
      case ( 'HIGH_ALKENES' ) ; vname = 'hialkenes'
      case ( 'TOL'          ) ; vname = 'toluene'
      case ( 'CH3CHO'       ) ; vname = 'c2h4o'
      case ( 'CH3COCH3'     ) ; vname = 'c3h6o'
      case ( 'BCFIRE'       ) ; vname = 'bc'
      case ( 'OCFIRE'       ) ; vname = 'oc'
      case ( 'PM'           ) ; vname = 'tpm'
      case ( 'PM2.5'        ) ; vname = 'pm2p5'
      case default
        ! just lower case:
        l = len_trim(speciesIdString)
        vname = goLoCase(speciesIdString(1:l))
    end select
    ! add extension:
    vname = trim(vname)//'fire'

    ! check for missing variable immediatelly:
    status = NF90_Inq_VarID( ncid, trim(vname), varid )
    if ( status /= 0 ) then
       ! info ...
       write (gol,'("         WARNING missing emission for species ",a)') trim(speciesIdString) ; call goPr
       ! close:
       status = NF90_Close( ncid )
       IF_NF90_NOTOK_RETURN(status=1)
       ! leave with warning status:
       status = -1; return
    end if

    ! * time

    ! single time per file:
    itime = 1
    ! number of time value:
    status = NF90_Inq_DimID( ncid, 'time', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=ntime )
    IF_NF90_NOTOK_RETURN(status=1)
    ! storage:
    allocate( times(ntime), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! read:
    status = NF90_Inq_VarID( ncid, 'time', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, times )
    IF_NF90_NOTOK_RETURN(status=1)
    ! time values:
    !   units    = 'days since 0-01-00 00:00:00' ;
    !   calendar = '366_day' ;
    status = NF90_Get_Att( ncid, varid, 'units', time_units )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Att( ncid, varid, 'calendar', time_calendar )
    IF_NF90_NOTOK_RETURN(status=1)

    
    ! expect field to be valid for 12:00 or 24:00 :
    !~ 12:00
    time_curr = NewDate( year=year, month=month, day=day, hour=12, calendar=time_calendar )
    ! check:
    call Compare_Date_Num( time_curr, times(itime), time_units, time_calendar, status )
    IF_ERROR_RETURN(status=1)
    !~ try other?
    if ( status < 0 ) then
      !~ 24:00
      time_curr = NewDate( year=year, month=month, day=day, hour=24, calendar=time_calendar )
      ! check:
      call Compare_Date_Num( time_curr, times(itime), time_units, time_calendar, status )
      IF_ERROR_RETURN(status=1)
      !~ try other?
      if ( status < 0 ) then
        !~ 00:00
        time_curr = NewDate( year=year, month=month, day=day, hour=0, calendar=time_calendar )
        ! check:
        call Compare_Date_Num( time_curr, times(itime), time_units, time_calendar, status )
        IF_ERROR_RETURN(status=1)
        if ( status < 0 ) then
          !~ 15:00
          time_curr = NewDate( year=year, month=month, day=day, hour=15, calendar=time_calendar )
          time_curr = time_curr - IncrDate(day=1)
          ! check:
          call Compare_Date_Num( time_curr, times(itime), time_units, time_calendar, status )
          IF_ERROR_RETURN(status=1)
          if ( status < 0 ) then
            write (gol,'("time in fire emis file is not one of 12:00, 15:00, 24:00, or 00:00")'); call goErr
            call wrtgol( '  current time : ', time_curr ); call goErr
            write (gol,'("  filename     : ",a)') trim(fname); call goErr
            TRACEBACK; status=1; return
          end if
        end if
      end if
    end if
    ! clear:
    deallocate( times, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! * compressed coordinate

    ! number of fires:
    status = NF90_Inq_DimID( ncid, 'point', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=npoint )
    IF_NF90_NOTOK_RETURN(status=1)
    ! check ...
    if ( npoint < 1 ) then
       write (gol,'("MACC fire emission file: nr of fires < 1 ")'); call goErr
       TRACEBACK; status=1; return
    end if
    ! compressed cell indices:
    allocate( points(npoint), stat=status )
    IF_NOTOK_RETURN(status=1)
    status = NF90_Inq_VarID( ncid, 'point', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, points )
    IF_NF90_NOTOK_RETURN(status=1)

    ! * grid

    ! lon axis:
    status = NF90_Inq_DimID( ncid, 'lon', dimid )
    if ( status /= NF90_NOERR ) then
      status = NF90_Inq_DimID( ncid, 'longitude', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    status = NF90_Inquire_Dimension( ncid, dimid, len=nlon )
    IF_NF90_NOTOK_RETURN(status=1)
    ! lon axis bounds:
    allocate( lon_bounds(2,nlon), stat=status )
    IF_NOTOK_RETURN(status=1)
    status = NF90_Inq_VarID( ncid, 'lon_bounds', varid )
    if ( status /= NF90_NOERR ) then
      status = NF90_Inq_VarID( ncid, 'longitude_bounds', varid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    status = NF90_Get_Var( ncid, varid, lon_bounds )
    IF_NF90_NOTOK_RETURN(status=1)

    ! lat axis:
    status = NF90_Inq_DimID( ncid, 'lat', dimid )
    if ( status /= NF90_NOERR ) then
      status = NF90_Inq_DimID( ncid, 'latitude', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    status = NF90_Inquire_Dimension( ncid, dimid, len=nlat )
    IF_NF90_NOTOK_RETURN(status=1)
    ! lat axis bounds:
    allocate( lat_bounds(2,nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    status = NF90_Inq_VarID( ncid, 'lat_bounds', varid )
    if ( status /= NF90_NOERR ) then
      status = NF90_Inq_VarID( ncid, 'latitude_bounds', varid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    status = NF90_Get_Var( ncid, varid, lat_bounds )
    IF_NF90_NOTOK_RETURN(status=1)

    ! * data

    ! storage:
    allocate( emis(npoint), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! read the emission time slab:
    status = NF90_Inq_VarID( ncid, trim(vname), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, emis, start=(/1,itime/), count=(/npoint,1/) )
    IF_NF90_NOTOK_RETURN(status=1)
    ! units:
    status = NF90_Get_Att( ncid, varid, 'units', emis_units )
    IF_NF90_NOTOK_RETURN(status=1)
    ! check  ...
    if ( trim(emis_units) /= 'kg m**-2 s**-1' ) then
      write (gol,'("macc fire emis only implemented for units `kg m**-2 s**-1`, found `",a,"`")') trim(emis_units); call goErr
      TRACEBACK; status=1; return
    end if


    ! *

    ! close netcdf file
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    !
    ! Reading done,
    ! information stored in "points", "lons", "lats", "fire_emis_list" (unit kg/s)
    !

    ! counter:
    n_in_domain = 0
    ! init fractions:
    nullify( fracs )
    ! loop over points:
    do ipoint = 1, npoint
      ! point gives the zero-based 1D index of the compressed coordinates,
      ! expand to 1-based:
      ilat = points(ipoint) / nlon + 1
      ilon = modulo( points(ipoint), nlon ) + 1
      ! boundingbox of emission:
      west  = lon_bounds(1,ilon)
      east  = lon_bounds(2,ilon)
      south = min( lat_bounds(1,ilat), lat_bounds(2,ilat) )
      north = max( lat_bounds(1,ilat), lat_bounds(2,ilat) )
      ! emission area:
      area = ll_area_deg_m2( west, east, south, north ) ! m2
      ! Given a box with sides [west,east]x[south,north] in degrees;
      ! determine how to distribute it's area over a regular lat/lon grid.
      ! Returns indices ranges i1..i2 and j1..j2 and fractions
      !   fracs(i1:i2,j1:j2)
      ! Return status:
      !   -1          : location outside domain
      !   0           : ok
      !   >0          : error
      call ugg%GetDistribution( west, east, south, north, &
                                 i1, i2, j1, j2, fracs, status )
      if ( status == -1 ) cycle  ! cell outside of target domain
      IF_NOTOK_RETURN(status=1)
      ! increase counter:
      n_in_domain = n_in_domain + 1
      ! loop over covered cells:
      do i = i1, i2
        do j = j1, j2
          ! add:              kg/s            kg/m2/s       m2
          fire_emis(i,j) = fire_emis(i,j) + emis(ipoint) * area * fracs(i,j)
        end do
      end do
    end do  ! points
    if ( verbose ) then
       write (gol,'("         nr of fire locations = ",i5,",  inside domain = ",i5)' ) &
                                      npoint, n_in_domain ; call goPr
    end if

    ! devide by grid cell area to convert to kg/m2/s :
    call ugg%AreaOper( fire_emis, '/', 'm2', status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    if ( trim(units) /= 'kg m**-2 s**-1' ) then
      write (gol,'("macc fire emis only implemented for target units `kg m**-2 s**-1`, requested `",a,"`")') trim(units); call goErr
      TRACEBACK; status=1; return
    end if

    ! clear:
    deallocate( points, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( lon_bounds, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( lat_bounds, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( emis, stat=status )
    IF_NOTOK_RETURN(status=1)
    if ( associated(fracs) ) then
      deallocate( fracs, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! info ...
    write (gol,'("         ",a," fire emissions added")' ) trim(speciesIdString) ; call goPr

    ! ok
    status = 0

  end subroutine maccFire_nc_ReadCompress

  ! *

  subroutine maccFire_nc_ReadCompress_hourly( path, speciesIdString, units, t, &
                                         verbose, nx, ny, fire_emis, status )

    use NetCDF, only : NF90_NoWrite, NF90_NOERR
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_StrError

    use GO     , only : TDate, NewDate, Get, Compare_Date_Num
    use GO     , only : operator(-), operator(+), operator(<), operator(<=)
    use GO     , only : IsAnyDate, IncrDate, wrtgol
    use GO     , only : goReplace, goLoCase, goSplitString
    use Grid   , only : ll_area_deg_m2
    use LE_Grid, only : ugg
    use Dims   , only : runF

    ! --- in/out --------------------------------

    character(len=*), intent(in)            ::  path
    character(len=*), intent(in)            ::  speciesIdString
    character(len=*), intent(in)            ::  units
    type(TDate), intent(in)                 ::  t
    logical, intent(in)                     ::  verbose
    integer, intent(in)                     ::  nx
    integer, intent(in)                     ::  ny
    real, intent(out)                       ::  fire_emis(nx,ny,0:23)
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/maccFire_nc_ReadCompress_hourly'

    ! max number of file templates:
    integer, parameter    ::  maxfile = 10

    ! --- local ---------------------------------

    character(len=1024)              ::  fnames(maxfile)
    integer                          ::  nfile, ifile
    character(len=1024)              ::  fname
    logical                          ::  exist
    integer                          ::  ncid, dimid, varid
    character(len=32)                ::  vname
    integer                          ::  l
    integer                          ::  year, month, day
    integer                          ::  itime
    integer                          ::  ntime
    real, allocatable                ::  times(:)
    character(len=32)                ::  time_units
    character(len=32)                ::  time_calendar
    type(TDate)                      ::  time_curr
    integer                          ::  npoint
    integer                          ::  ipoint
    integer, allocatable             ::  points(:)
    integer                          ::  nlon, nlat
    integer                          ::  ilon, ilat
    real, allocatable                ::  lon_bounds(:,:), lat_bounds(:,:)
    real, allocatable                ::  emis(:,:)
    character(len=32)                ::  emis_units
    integer                          ::  n_in_domain
    integer                          ::  i1, i2, j1, j2
    integer                          ::  i, j, hh
    real, pointer                    ::  fracs(:,:)
    real                             ::  west, east, south, north
    real                             ::  area

    ! --- begin -------------------------------

    ! Initialise
    fire_emis(:,:,:) = 0.0

    ! forecast base specified? 
    if ( IsAnyDate(runF%t_base) ) then
      ! not a forecast run i.e. hindcast or NRT analysis for CAMS;
      ! target time:
      call Get( t, year=year, month=month, day=day )
    else
      ! forecast run;
      ! if before 2 days after forecast base use fires
      ! from day before forecast day:
      if ( t <= runF%t_base+IncrDate(day=2) ) then
        ! use analysis day:
        call Get( runF%t_base-IncrDate(day=1), year=year, month=month, day=day )
        ! info ...
        if ( verbose ) then
          write (gol,'("WARNING - no fire emission forecast, use data from day before forecast base!")'); call goPr
          write (gol,'("        - i.e. ",i4, 2i3)') year, month, day; call goPr
        end if
        !
      ! no data ..
      else
        ! info ...
        if ( verbose ) then
          write (gol,'("WARNING - no fire emission forecast, keep zero at this day!")'); call goPr
        end if
        ! ok
        status=0; return
      end if
    end if

    ! path might consist of several values,
    ! should be tested until the first match:
    call goSplitString( trim(path), nfile, fnames, status )
    IF_NOTOK_RETURN(status=1)
    ! dummy:
    fname = '/no/file/found/yet'
    ! loop over directories:
    do ifile = 1, nfile
      ! replace some values:
      call goReplace( fnames(ifile), '%{year}', '(i4.4)', year, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( fnames(ifile), '%{month}', '(i2.2)', month, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( fnames(ifile), '%{day}', '(i2.2)', day, status )
      IF_NOTOK_RETURN(status=1)
      ! check, if present then store match and leave:
      inquire( file=trim(fnames(ifile)), exist=exist )
      if ( exist ) then
        fname = trim(fnames(ifile))
        exit
      end if
    end do
    ! check ..
    inquire( file=trim(fname), exist=exist )

    ! if file does not exist and it is not a forecast run, try the previous day.
    ! Useful for NRT analyses 
    if ( .not. exist .and. IsAnyDate(runF%t_base) ) then

      ! info ..
      if ( verbose ) then
        write (gol,'("WARNING - no fire emission found for the day, will try to use the one from the day previous to forecast base!")'); call goPr
      end if
      ! time values:
      call Get( t-IncrDate(day=1), year=year, month=month, day=day )

      ! path might consist of several values,
      ! should be tested until the first match:
      call goSplitString( trim(path), nfile, fnames, status )
      IF_NOTOK_RETURN(status=1)
      ! dummy:
      fname = '/no/file/found/yet'
      ! loop over directories:
      do ifile = 1, nfile
        ! replace some values:
        call goReplace( fnames(ifile), '%{year}', '(i4.4)', year, status )
        IF_NOTOK_RETURN(status=1)
        call goReplace( fnames(ifile), '%{month}', '(i2.2)', month, status )
        IF_NOTOK_RETURN(status=1)
        call goReplace( fnames(ifile), '%{day}', '(i2.2)', day, status )
        IF_NOTOK_RETURN(status=1)
        ! check, if present then store match and leave:
        inquire( file=trim(fnames(ifile)), exist=exist )
        if ( exist ) then
          fname = trim(fnames(ifile))
          exit
        end if
      end do
      ! check ..
      inquire( file=trim(fname), exist=exist )

    end if

    ! if file still can't be found
    if ( .not. exist ) then
      if ( arch_allow_missing) then
        write (gol,'("WARNING - missing MACC fire emissions file; continue ...")'); call goPr
        status=0; return
      else
        write (gol,'("None of the specified MACC fire emission files found:")'); call goErr
        do ifile = 1, nfile
          write (gol,'("  ",a)') trim(fnames(ifile)); call goErr
        end do
        TRACEBACK; status=1; return
      end if
    end if

    ! info:
    if ( verbose ) then
      write (gol,'("         Read MACC fire emissions for ",i4, 2i3)') year, month, day; call goPr
      write (gol,'("           from file : ",a)') trim(fname); call goPr
    end if

    ! open netcdf file:
    status = NF90_Open( trim(fname), NF90_NoWrite, ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! variable name; first translate species name:
    select case ( trim(speciesIdString) )
      case ( 'NO'           ) ; vname = 'nox'
      case ( 'HIGH_ALKANES' ) ; vname = 'hialkanes'
      case ( 'ISOP'         ) ; vname = 'c5h8'
      case ( 'TERP'         ) ; vname = 'terpenes'
      case ( 'HIGH_ALKENES' ) ; vname = 'hialkenes'
      case ( 'TOL'          ) ; vname = 'toluene'
      case ( 'CH3CHO'       ) ; vname = 'c2h4o'
      case ( 'CH3COCH3'     ) ; vname = 'c3h6o'
      case ( 'BCFIRE'       ) ; vname = 'bc'
      case ( 'OCFIRE'       ) ; vname = 'oc'
      case ( 'PM'           ) ; vname = 'tpm'
      case ( 'PM2.5'        ) ; vname = 'pm2p5'
      case default
        ! just lower case:
        l = len_trim(speciesIdString)
        vname = goLoCase(speciesIdString(1:l))
    end select
    ! add extension:
    vname = trim(vname)//'fire'

    ! check for missing variable immediatelly:
    status = NF90_Inq_VarID( ncid, trim(vname), varid )
    if ( status /= 0 ) then
       ! info ...
       write (gol,'("         WARNING missing emission for species ",a)') trim(speciesIdString) ; call goPr
       ! close:
       status = NF90_Close( ncid )
       IF_NF90_NOTOK_RETURN(status=1)
       ! leave with warning status:
       status = -1; return
    end if

    ! * time

    ! single time per file:
    itime = 24
    ! number of time value:
    status = NF90_Inq_DimID( ncid, 'time', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=ntime )
    IF_NF90_NOTOK_RETURN(status=1)
    ! storage:
    allocate( times(ntime), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! read:
    status = NF90_Inq_VarID( ncid, 'time', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, times )
    IF_NF90_NOTOK_RETURN(status=1)
    ! time values:
    !   units    = 'days since 0-01-00 00:00:00' ;
    !   calendar = '366_day' ;
    status = NF90_Get_Att( ncid, varid, 'units', time_units )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Att( ncid, varid, 'calendar', time_calendar )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( times, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! * compressed coordinate

    ! number of fires:
    status = NF90_Inq_DimID( ncid, 'point', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=npoint )
    IF_NF90_NOTOK_RETURN(status=1)
    ! check ...
    if ( npoint < 1 ) then
       write (gol,'("MACC fire emission file: nr of fires < 1 ")'); call goErr
       TRACEBACK; status=1; return
    end if
    ! compressed cell indices:
    allocate( points(npoint), stat=status )
    IF_NOTOK_RETURN(status=1)
    status = NF90_Inq_VarID( ncid, 'point', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, points )
    IF_NF90_NOTOK_RETURN(status=1)

    ! * grid

    ! lon axis:
    status = NF90_Inq_DimID( ncid, 'lon', dimid )
    if ( status /= NF90_NOERR ) then
      status = NF90_Inq_DimID( ncid, 'longitude', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    status = NF90_Inquire_Dimension( ncid, dimid, len=nlon )
    IF_NF90_NOTOK_RETURN(status=1)
    ! lon axis bounds:
    allocate( lon_bounds(2,nlon), stat=status )
    IF_NOTOK_RETURN(status=1)
    status = NF90_Inq_VarID( ncid, 'lon_bounds', varid )
    if ( status /= NF90_NOERR ) then
      status = NF90_Inq_VarID( ncid, 'longitude_bounds', varid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    status = NF90_Get_Var( ncid, varid, lon_bounds )
    IF_NF90_NOTOK_RETURN(status=1)

    ! lat axis:
    status = NF90_Inq_DimID( ncid, 'lat', dimid )
    if ( status /= NF90_NOERR ) then
      status = NF90_Inq_DimID( ncid, 'latitude', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    status = NF90_Inquire_Dimension( ncid, dimid, len=nlat )
    IF_NF90_NOTOK_RETURN(status=1)
    ! lat axis bounds:
    allocate( lat_bounds(2,nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    status = NF90_Inq_VarID( ncid, 'lat_bounds', varid )
    if ( status /= NF90_NOERR ) then
      status = NF90_Inq_VarID( ncid, 'latitude_bounds', varid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    status = NF90_Get_Var( ncid, varid, lat_bounds )
    IF_NF90_NOTOK_RETURN(status=1)

    ! * data

    ! storage:
    allocate( emis(npoint,24), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! read the emission time slab:
    status = NF90_Inq_VarID( ncid, trim(vname), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, emis, start=(/1,1/), count=(/npoint,itime/) )
    IF_NF90_NOTOK_RETURN(status=1)
    ! units:
    status = NF90_Get_Att( ncid, varid, 'units', emis_units )
    IF_NF90_NOTOK_RETURN(status=1)
    ! check  ...
    if ( trim(emis_units) /= 'kg m**-2 s**-1' .and. trim(emis_units) /= 'kg m-2 s-1') then
      write (gol,'("macc fire emis only implemented for units `kg m**-2 s**-1`, found `",a,"`")') trim(emis_units); call goErr
      TRACEBACK; status=1; return
    end if

    ! *

    ! close netcdf file
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    !
    ! Reading done,
    ! information stored in "points", "lons", "lats", "fire_emis_list" (unit kg/s)
    !

    ! counter:
    n_in_domain = 0
    ! init fractions:
    nullify( fracs )
    ! loop over points:
    do ipoint = 1, npoint
      ! point gives the zero-based 1D index of the compressed coordinates,
      ! expand to 1-based:
      ilat = points(ipoint) / nlon + 1
      ilon = modulo( points(ipoint), nlon ) + 1
      ! boundingbox of emission:
      west  = lon_bounds(1,ilon)
      east  = lon_bounds(2,ilon)
      south = min( lat_bounds(1,ilat), lat_bounds(2,ilat) )
      north = max( lat_bounds(1,ilat), lat_bounds(2,ilat) )
      ! emission area:
      area = ll_area_deg_m2( west, east, south, north ) ! m2
      ! Given a box with sides [west,east]x[south,north] in degrees;
      ! determine how to distribute it's area over a regular lat/lon grid.
      ! Returns indices ranges i1..i2 and j1..j2 and fractions
      !   fracs(i1:i2,j1:j2)
      ! Return status:
      !   -1          : location outside domain
      !   0           : ok
      !   >0          : error
      call ugg%GetDistribution( west, east, south, north, &
                                 i1, i2, j1, j2, fracs, status )
      if ( status == -1 ) cycle  ! cell outside of target domain
      IF_NOTOK_RETURN(status=1)
      ! increase counter:
      n_in_domain = n_in_domain + 1
      ! loop over covered cells:
      do i = i1, i2
        do j = j1, j2
          do hh = 0, 23
            ! add:              kg/s            kg/m2/s       m2
            fire_emis(i,j,hh) = fire_emis(i,j,hh) + emis(ipoint,hh+1) * area * fracs(i,j)
            ! hh+1 because emis is read in the fortran way
          end do
        end do
      end do
    end do  ! points
    if ( verbose ) then
       write (gol,'("         nr of fire locations = ",i5,",  inside domain = ",i5)' ) &
                                      npoint, n_in_domain ; call goPr
    end if

    ! devide by grid cell area to convert to kg/m2/s :
    call ugg%AreaOper( fire_emis, '/', 'm2', status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    if ( trim(units) /= 'kg m**-2 s**-1' .and. trim(emis_units) /= 'kg m-2 s-1' ) then
      write (gol,'("macc fire emis only implemented for target units `kg m**-2 s**-1`, requested `",a,"`")') trim(units); call goErr
      TRACEBACK; status=1; return
    end if

    ! clear:
    deallocate( points, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( lon_bounds, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( lat_bounds, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( emis, stat=status )
    IF_NOTOK_RETURN(status=1)
    if ( associated(fracs) ) then
      deallocate( fracs, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! info ...
    write (gol,'("         ",a," fire emissions added")' ) trim(speciesIdString) ; call goPr

    ! ok
    status = 0

  end subroutine maccFire_nc_ReadCompress_hourly


end module LE_Emis_Fire_MACC
