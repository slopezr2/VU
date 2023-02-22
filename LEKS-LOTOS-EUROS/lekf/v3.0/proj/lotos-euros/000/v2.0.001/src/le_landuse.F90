!###############################################################################
!
! landuse  !
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

module LE_LandUse

  use GO, only : gol, goPr, goErr
!  use LE_Landuse_Data
  
  implicit none
  
  ! --- in/out -----------
  
  private
  
  public    ::  LE_Landuse_Init
  public    ::  LE_Landuse_Done
  public    ::  LE_Landuse_Update  
  
  ! --- const -------------------

  character(len=*), parameter   ::  mname = 'LE_Landuse'

  ! --- var -------------------
  
contains

  subroutine LE_Landuse_Init( rcF, status )
    
    use GO, only : ReadRc, TrcFile
    use GO, only : GoVarValue
    use GO, only : goSplitString
    
    use indices, only : n_dust
    
    use Dims, only : nx, ny

    use LE_Landuse_Data   , only : lu_fracs
    use LE_Landuse_Data   , only : ilu_water_sea, ilu_desert, ilu_arable, ilu_water_inland
    use LE_Landuse_Data   , only : with_soil_texture 
    use LE_Landuse_Data   , only : with_ozone_specials 
    use LE_Landuse_Data   , only : use_clc2006_landuse
    use LE_Landuse_Data   , only : LE_Landuse_Data_Alloc
    use LE_Landuse_Data   , only : nlu_names_gs, lu_names_gs
    use LE_Landuse_Data   , only : filenames_gs
    use LE_Landuse_Data   , only : ludep_output_whole_grid
    use LE_Landuse_Data   , only : waterbodies
     
    use LE_Landuse_Smiatek, only : LE_Landuse_Smiatek_Init
    use LE_Landuse_GLC2000, only : LE_Landuse_GLC2000_Init
    use LE_Landuse_CLC2006, only : LE_Landuse_CLC2006_Init
    use LE_Landuse_z0     , only : LE_Landuse_z0_Init
    use LE_LandUse_LSM    , only : LandUse_LSM_Init
    use LE_LandUse_Soil   , only : LandUse_Soil_Init
    use LE_LandUse_Traffic, only : LandUse_Traffic_Init
    use LE_LandUse_BCatIon, only : LandUse_BCatIon_Init
    
    ! ---in/out --------------
    type(TrcFile), intent(in) ::  rcF
    integer, intent(out)      ::  status
    
    ! --- const --------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Landuse_Init'
        
    ! --- local --------------
    
    character(len=512)    ::  query
    character(len=32)     ::  landuse_type
    character(len=512)    ::  list
    integer               ::  ilu_name_gs
    integer               ::  ix, iy
    
    ! --- begin --------------
    
    ! Use special landuse categories for ozone flux?
    call ReadRC( rcf, 'deposition.ozone_flux_specials', with_ozone_specials, status)
    IF_NOTOK_RETURN(status=1)
    
    ! allocate variables
    call LE_Landuse_Data_Alloc( status )
    IF_NOTOK_RETURN(status=1)
    
    ! read in landuse fractions
    call ReadRc( rcF, 'landuse.file', query, status )
    IF_NOTOK_RETURN(status=1)
    
    call GoVarValue( trim(query), ';', 'ftype', '=', landuse_type, status )
    IF_NOTOK_RETURN(status=1)
    
    ! which type of landuse database
    select case ( trim(landuse_type) )
      
      case ( 'smiatek' )
        
        ! initialize Smiatek landuse
        call LE_Landuse_Smiatek_Init( trim(query), status )
        IF_NOTOK_RETURN(status=1)
        
      case ( 'glc2000' )

        ! initialize GLC2000 landuse
        call LE_Landuse_GLC2000_Init( trim(query), status )
        IF_NOTOK_RETURN(status=1)
        
      case default
        
        ! unknown landcover
        write( gol, '("Unknown Landuse database: ", a)') trim(landuse_type) ; call GoErr
        TRACEBACK;status=1;return

    end select
    
    ! Read in Corine 2006 landuse database ?
    call ReadRc( rcF, 'my.landuse.clc2006.enabled', use_clc2006_landuse, status )
    IF_NOTOK_RETURN(status=1)

    if ( use_clc2006_landuse ) then

      ! read in clc2006 file
      call LE_Landuse_CLC2006_Init( rcF, status )
      IF_NOTOK_RETURN(status=1)
      
    endif

    ! * land/sea mask
    call Landuse_LSM_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
           
    ! Divide waterfraction in original file in Sea and Inland water.
    do ix = 1, nx
    do iy = 1, ny
      
      if ( lu_fracs(ix,iy,ilu_water_sea) > 0 ) then
      
        ! check waterbody 
        
        ! all water in cell is sea
        if ( waterbodies(ix,iy, 2) >= lu_fracs(ix,iy,ilu_water_sea) ) then
          cycle
          
        ! part of the water in the cell is inland  
        else 
          lu_fracs(ix,iy,ilu_water_inland) = lu_fracs(ix,iy,ilu_water_sea) - waterbodies(ix,iy, 2)
          lu_fracs(ix,iy,ilu_water_sea   ) = waterbodies(ix,iy, 2)
        end if
      
      end if ! water in cell?
      
    end do ! iy
    end do ! ix

    ! initialize z0 parameters
    call LE_Landuse_z0_Init( status )
    IF_NOTOK_RETURN(status=1)
       
    !  o dust emissions (wind blown dust)
    with_soil_texture = n_dust > 0
    if ( with_soil_texture ) then
      ! initialize soil-texture map
      call Landuse_Soil_Init( rcF, status, &
                               bare_soil = lu_fracs(:,:,ilu_desert), &
                               arable_soil = lu_fracs(:,:,ilu_arable) )
      IF_NOTOK_RETURN(status=1)
    end if
    

    ! * base-cat-ion composition

    ! read soil composition:
    call Landuse_BCatIon_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    

    ! * traffic

    ! read traffic maps:
    call Landuse_Traffic_Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
       
    ! Read in info on growing seasons
    call ReadRc( rcF, 'le.landuse.accumulation.types', list, status)
    IF_NOTOK_RETURN(status=1)
    
    ! which classes
    call goSplitString( list, nlu_names_gs, lu_names_gs, status )
    IF_NOTOK_RETURN(status=1)
  
    ! read in filenames  
    allocate( filenames_gs( nlu_names_gs ) )    
    do ilu_name_gs = 1, nlu_names_gs
      call ReadRc( rcF, 'le.landuse.accumulation.'//trim(lu_names_gs(ilu_name_gs))//'.file', filenames_gs(ilu_name_gs), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! ok
    status = 0
    
  end subroutine LE_Landuse_Init


  ! ***
  
  
  subroutine LE_Landuse_Done( status )
      
    use LE_Landuse_Data   , only : with_soil_texture  
    use LE_Landuse_Data   , only : LE_Landuse_Data_Dealloc
    use LE_LandUse_LSM    , only : LandUse_LSM_Done
    use LE_LandUse_Soil   , only : LandUse_Soil_Done
    use LE_LandUse_Traffic, only : LandUse_Traffic_Done
    use LE_LandUse_BCatIon, only : LandUse_BCatIon_Done

    ! --- in/out ----------
    integer, intent(out)    ::  status
    
    ! --- const -----------

    character(len=*), parameter ::  rname = mname//'/LE_Landuse_Done'
    
    ! --- local -----------
    
    ! deallocate
    call LE_Landuse_Data_Dealloc( status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with land/sea mask
    call Landuse_LSM_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! for dust purposes:
    if ( with_soil_texture ) then
      call LandUse_Soil_Done( status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! done with traffic maps:
    call Landuse_Traffic_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with base-cat-ion composition:
    call Landuse_BCatIon_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok 
    status = 0
    
  end subroutine LE_Landuse_Done


  ! ***
  
  
  subroutine LE_Landuse_Update( t, lfirst, status )
    
    use GO              , only : TDate, MidNight
   
    use LE_Landuse_z0   , only : LE_Landuse_z0_Setup

    
    ! --- in/out ----------
    type(TDate)             ::  t
    logical, intent(in)     ::  lfirst
    integer, intent(out)    ::  status
    
    ! --- const -----------

    character(len=*), parameter ::  rname = mname//'/LE_Landuse_Update'
    
    ! --- local -----------

    ! --- begin -----------    
    
    if ( Midnight(t) .or. lfirst ) then
      ! update information on growing season files
      call LE_Landuse_Update_GS( t, status )
      IF_NOTOK_RETURN(status=1)
    
      ! Update Leaf area index and surface area index (lai/sai)    
      call LE_Landuse_Update_LAI_SAI( t, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! update z0 values
    call LE_Landuse_z0_Setup( status )
    IF_NOTOK_RETURN(status=1)

    
    ! ok 
    status = 0
    
  end subroutine LE_Landuse_Update

  
  ! Update information about growing seasons
  
  subroutine LE_Landuse_Update_GS( t, status )
  
    use GO, only : TDate
    use GO, only : goReplace
    use Dims, only : nx, ny
    use LE_Landuse_Data, only : start_gs, end_gs
    use LE_Landuse_Data, only : nlu_names_gs, lu_names_gs
    use LE_Landuse_Data, only : filenames_gs
    use LE_Landuse_Data, only : current_year
    use LE_Landuse_Data, only : ilu_semi_natural_veg, ilu_wheat
    use LE_Landuse_Data, only : ilu_beech, ilu_spruce
    use Grid, only : TllGridInfo, Init, Done,FillGrid_AreaAverage
    use LE_Grid, only : lli
    use MDF, only : MDF_Open, MDF_NETCDF, MDF_READ, MDF_Close
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var
    
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)     ::  t
    integer, intent(out)        ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_LandUse_Update_GS'
        
    ! --- local -------------------------------

    character(len=512)  ::  fname
    
    integer             ::  i, j, ilu_name_gs
    logical             ::  exist
    
    
    integer             ::  ncid, dimid, varid
    integer             ::  nlon, nlat
    real, allocatable   ::  lons(:), lats(:)
    real, allocatable   ::  sgs(:,:), egs(:,:)

    type(TllGridInfo)   ::  lli_lu
    
    ! --- begin -------
    
    
    if ( .not. (t%year == current_year) ) then
  
      ! loop over landuses with growing seasons    
      do ilu_name_gs = 1, nlu_names_gs
        
        ! pick filenamebase
        fname = filenames_gs(ilu_name_gs)
              
        call goReplace( fname, '%{year}', '(i4.4)', t%year, status )
        IF_NOTOK_RETURN(status=1)

        inquire( file=trim(fname), exist=exist)
        if ( .not. exist ) then
          write (gol,'("file not found:")'); call goErr
          write (gol,'("  ",a)') trim(fname); call goErr
          TRACEBACK; status=1; return
        end if

        write ( gol, '("Open accumulation file for: ",a," -- ", a)') trim(lu_names_gs(ilu_name_gs)), trim(fname); call goPr

        ! open file:
        call MDF_Open( trim(fname), MDF_NETCDF, MDF_READ, ncid, status )
        IF_NOTOK_RETURN(status=1)

        !get dimid: longitudes
        call MDF_Inq_DimID( ncid, 'longitude', dimid, status )
        IF_NOTOK_RETURN(status=1)
        call MDF_Inquire_Dimension( ncid, dimid, status, length=nlon )
        IF_NOTOK_RETURN(status=1)

        allocate(lons (nlon) )
        ! get variable id longitudes
        call MDF_Inq_VarID( ncid, 'grid_longitude', varid, status )
        IF_NOTOK_RETURN(status=1)
        call MDF_Get_Var( ncid, varid, lons, status )
        IF_NOTOK_RETURN(status=1)

        !get dimid: latitudes
        call MDF_Inq_DimID( ncid, 'latitude', dimid, status )
        IF_NOTOK_RETURN(status=1)
        call MDF_Inquire_Dimension( ncid, dimid, status, length=nlat )
        IF_NOTOK_RETURN(status=1)

        allocate( lats(nlat) )
        ! get variable id longitudes
        call MDF_Inq_VarID( ncid, 'grid_latitude', varid, status )
        IF_NOTOK_RETURN(status=1)
        call MDF_Get_Var( ncid, varid, lats, status )
        IF_NOTOK_RETURN(status=1)

        ! setup input grid:
        call Init( lli_lu, lons(1), lons(2)-lons(1), nlon, &
                        lats(1), lats(2)-lats(1), nlat, status )
        IF_NOTOK_RETURN(status=1)

        !  read in start and end of growing season
        allocate( sgs(nlon,nlat) )
        call MDF_Inq_VarID( ncid, 'start_days', varid, status )
        IF_NOTOK_RETURN(status=1)
        call MDF_Get_Var( ncid, varid, sgs, status)
        IF_NOTOK_RETURN(status=1)

        allocate( egs(nlon,nlat) )
        call MDF_Inq_VarID( ncid, 'end_days', varid, status )
        IF_NOTOK_RETURN(status=1)
        call MDF_Get_Var( ncid, varid, egs, status)
        IF_NOTOK_RETURN(status=1)

        ! Convert start/end growing season to landuse dependent variables
        if ( trim(lu_names_gs(ilu_name_gs)) == 'Beech' ) then
          call FillGrid_AreaAverage( lli, start_gs(:,:,ilu_beech), lli_lu, sgs, status)
          IF_NOTOK_RETURN(status=1)
          call FillGrid_AreaAverage( lli, end_gs(:,:,ilu_beech), lli_lu, egs, status)
          IF_NOTOK_RETURN(status=1)

        else if ( trim(lu_names_gs(ilu_name_gs)) == 'Clover' ) then 
          call FillGrid_AreaAverage( lli, start_gs(:,:,ilu_semi_natural_veg), lli_lu, sgs, status)
          IF_NOTOK_RETURN(status=1)
          call FillGrid_AreaAverage( lli, end_gs(:,:,ilu_semi_natural_veg), lli_lu, egs, status)
          IF_NOTOK_RETURN(status=1)

        else if ( trim(lu_names_gs(ilu_name_gs)) == 'Wheat' ) then 
          call FillGrid_AreaAverage( lli, start_gs(:,:,ilu_wheat), lli_lu, sgs, status)
          IF_NOTOK_RETURN(status=1)
          call FillGrid_AreaAverage( lli, end_gs(:,:,ilu_wheat), lli_lu, egs, status)
          IF_NOTOK_RETURN(status=1)

        else if ( trim(lu_names_gs(ilu_name_gs)) == 'Spruce' ) then 
          call FillGrid_AreaAverage( lli, start_gs(:,:,ilu_spruce), lli_lu, sgs, status)
          IF_NOTOK_RETURN(status=1)
          call FillGrid_AreaAverage( lli, end_gs(:,:,ilu_spruce), lli_lu, egs, status)
          IF_NOTOK_RETURN(status=1)

        end if
        
        ! filled
        current_year = t%year

        ! done with inut grid
        call Done( lli_lu, status)
        IF_NOTOK_RETURN(status=1)

        call MDF_Close( ncid, status ) 
        IF_NOTOK_RETURN(status=1)

        deallocate( lons )
        deallocate( lats )
        deallocate( sgs )
        deallocate( egs )
      end do
  
    end if ! check current year
    
    ! ok
    status = 0
  
  end subroutine LE_Landuse_Update_GS
  
  ! Update information about growing seasons
  
  subroutine LE_Landuse_Update_LAI_SAI( t, status )
    
    use GO, only : TDate, DayNumber
    use Dims            , only : nx, ny
    
    use LE_grid         , only : lli 
    use LE_Landuse_Data , only : nlu, lu_fracs
    use LE_Landuse_Data , only : ilu_grass, ilu_arable, ilu_permanent_crops   
    use LE_Landuse_Data , only : ilu_coniferous_forest, ilu_deciduous_forest  
    use LE_Landuse_Data , only : ilu_water_sea, ilu_urban, ilu_other             
    use LE_Landuse_Data , only : ilu_desert, ilu_ice, ilu_savanna            
    use LE_Landuse_Data , only : ilu_tropical_forest, ilu_water_inland   
    use LE_Landuse_Data , only : ilu_mediterrean_scrub, ilu_semi_natural_veg
    use LE_Landuse_Data , only : ilu_wheat, ilu_beech, ilu_spruce
    use LE_Landuse_Data , only : with_ozone_specials
    use LE_Landuse_Data , only : lai, sai
    use LE_Landuse_Data , only : laitype, lai_par
    use LE_Landuse_Data , only : start_gs, end_gs
    use LE_Landuse_Data , only : ludep_output_whole_grid
  
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)     ::  t
    integer, intent(out)        ::  status
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_LandUse_Update_GS'
    
    ! --- local -------------------------------
    
    type(laitype) :: lai_par1 ! lai parameters 
    
    real    ::  lat   ! latitude of the grid cell 
    real    ::  lai_loc, sai_loc  ! lai/sai index of a grid cell
    integer ::  ix, iy, ilu
    integer ::  day_of_year
    integer ::  sgs   ! start growing season at certain latitude (days)
    integer ::  egs   ! end growing season at certain latitude (days)
    
    
    ! --- begin ----
    
    ! variables from module
    !       LAI
    !        |
    !  laimax|               ---------
    !        |              /         \
    !        |             /           \
    !        |            /             \
    !        |           /               \
    !        |          /                 \
    !  laimin|         /                   \
    !        |        |                     |
    !        |        |                     |
    !        ---------|------|--------|-----|-----------
    !                sgs    sgs+     egs-   egs
    !                     s_lai_len  e_lai_len
    !
    ! contents of lai_par:
    ! sgs50     ! start growing season at 50 degrees latitude (days)
    ! dsgs      ! shift in start growing season (days/degree latitude)
    !             dsgs > 0 -> start-growing-season earlier in the south, later in the north
    ! egs50     ! end growing season at 50 degrees latitude (days)
    ! degs      ! shift in end growing season (days/degree latitude)
    !             degs < 0 -> end-growing-season later in the south, earlier in the north
    ! laimin    ! leaf area index at start and end of growing season (m2 leaf/m2 ground surface);
                ! outside growing season LAI = 0. Note that the SAI can be > 0 outside the
                ! growing season.
    ! laimax    ! maximal leaf area index (m2 leaf/m2 ground surface)
    ! s_lai_len ! length of starting phase of LAI (days)
    ! e_lai_len ! length of end phase of LAI (days)

    ! Compute start and end of growing season for current latitude;
    ! dsgs > 0 -> sgs earlier in the south, later in the north
    ! degs < 0 -> egs later in the south, earlier in the north
    
    day_of_year = DayNumber(t)
    ! fill lai and sai for domain and landuse class
    do ix = 1, nx
    do iy = 1, ny
      ! ~ latitude :
      lat = lli%lat_deg(iy)
      do ilu = 1, nlu
        if ( lu_fracs(ix,iy,ilu) <= 0.0 .and. ( .not.(ludep_output_whole_grid(ilu)) )) cycle
        
        ! parameters for this landuse
        lai_par1 = lai_par(ilu)
        
        ! set start/end grwoing season from file for extra landuse classes. 
        ! Data based on Northern Hemisphere, 20<lat<80, centered around 50
        if ( with_ozone_specials .and. &
              ( ilu .eq. ilu_wheat  .or. &
                ilu .eq. ilu_beech  .or. &
                ilu .eq. ilu_spruce .or. & 
                ilu .eq. ilu_semi_natural_veg) ) then
          sgs = start_gs(ix,iy,ilu)
          egs = end_gs(ix,iy,ilu)
        else     ! use original setting for 'normal' landuse classes 
          sgs = int ( 0.5 +  lai_par1%sgs50 + lai_par1%dsgs * (abs(lat)-50.0) )
          egs = int ( 0.5 +  lai_par1%egs50 + lai_par1%degs * (abs(lat)-50.0) )
        end if

        if (missing(lai_par1%laimax)) then
          lai_loc = -999.0
          sai_loc = -999.0
        else
          ! calculation of lai
          if (lat > 20) then
            !Northern hemisphere
       
            if ( day_of_year < sgs .or. day_of_year > egs ) then
              lai_loc = 0.0;
            else if (day_of_year <=  sgs + lai_par1%s_lai_len ) then

              ! (lai    - laimin)   (day_of_year - sgs)
              ! ----------------- = -------------------
              ! (laimax - laimin)        s_lai_len

              lai_loc = lai_par1%laimin + (lai_par1%laimax-lai_par1%laimin)*(day_of_year-sgs)/lai_par1%s_lai_len

            else if (day_of_year >=  egs - lai_par1%e_lai_len) then

              ! (lai    - laimin)   (egs - day_of_year)
              ! ----------------- = -------------------
              ! (laimax - laimin)        e_lai_len

              lai_loc = lai_par1%laimin + (lai_par1%laimax-lai_par1%laimin)*(egs-day_of_year)/lai_par1%e_lai_len
            else
              lai_loc = lai_par1%laimax
            endif
         else if (lat <-20) then
            !Southern hemisphere
            !shift growing seasons
            !growing seasons defined for Norhtern hemisphere 
            !sgs=sgs+183, assume that start of growing season always in second half of yearsgs <183
            !egs=egs+183-365, assume that end of growing season always in first half of year egs>183
            if ( day_of_year < sgs+183 .and. day_of_year > egs+183-365 ) then
                lai_loc = 0.0;
            else if (day_of_year >=  sgs+183 + lai_par1%s_lai_len ) then

              ! (lai    - laimin)   (day_of_year - (sgs+183))
              ! ----------------- = -------------------
              ! (laimax - laimin)        s_lai_len

                lai_loc = lai_par1%laimin + (lai_par1%laimax-lai_par1%laimin)*(day_of_year-sgs-183)/lai_par1%s_lai_len

             else if (day_of_year <=  (egs+183-365) - lai_par1%e_lai_len) then

              ! (lai    - laimin)   (egs+183-365 - day_of_year)
              ! ----------------- = -------------------
              ! (laimax - laimin)        e_lai_len

              lai_loc = lai_par1%laimin + (lai_par1%laimax-lai_par1%laimin)*(egs+183-385-day_of_year)/lai_par1%e_lai_len
            else
              lai_loc = lai_par1%laimax
            endif
          else
            !Equatorial band
            !tropical region, no seasons, always maximum leaf area index
            lai_loc = lai_par1%laimax
          endif

          ! calculation of sai following EMEP report 1/2003 Simpson et al.
          if (ilu .eq. ilu_coniferous_forest .or. ilu .eq. ilu_deciduous_forest .or. ilu .eq. ilu_tropical_forest ) then ! forest
            sai_loc = lai_loc + 1.                                     
          else if (ilu .eq. ilu_permanent_crops) then  ! permanent crop
            sai_loc = lai_loc + .5                             ! estimate based on pers. comm. Roy Wichink Kruit

          else if (ilu .eq. ilu_arable ) then          ! arable land or wheat
            if (lat >20) then
				! Northern Hemisphere
                if (day_of_year < sgs .or. day_of_year > egs) then            
                  sai_loc = lai_loc
                else if (day_of_year <=  sgs + lai_par1%s_lai_len) then
                  sai_loc = lai_loc + lai_loc/lai_par1%laimax *1.5    ! avoid overshoot for sai (different from Simpson et al.(2003) due to bug in EMEP description)
                else if (day_of_year >=  egs - lai_par1%e_lai_len) then
                  sai_loc = lai_loc + 1.5
                else
                 sai_loc = lai_loc + 1.5
                endif
            else if (lat < -20) then
				! Southern Hemisphere
                if (day_of_year < sgs+183 .and. day_of_year > egs+183-365) then            
                  sai_loc = lai_loc
                else if (day_of_year <=  sgs+183 + lai_par1%s_lai_len) then
                  sai_loc = lai_loc + lai_loc/lai_par1%laimax *1.5    ! avoid overshoot for sai (different from Simpson et al.(2003) due to bug in EMEP description)
                else if (day_of_year >=  egs+183-365 - lai_par1%e_lai_len) then
                  sai_loc = lai_loc + 1.5
                else
                 sai_loc = lai_loc + 1.5
                endif
            else !equatorial band
                sai_loc = lai_loc + 1.5
            end if
          else if (ilu .eq. ilu_grass             .or. &
                   ilu .eq. ilu_water_sea         .or. &
                   ilu .eq. ilu_water_inland      .or. &
                   ilu .eq. ilu_urban             .or. &
                   ilu .eq. ilu_other             .or. &
                   ilu .eq. ilu_desert            .or. &
                   ilu .eq. ilu_ice               .or. &
                   ilu .eq. ilu_savanna           .or. &
                   ilu .eq. ilu_mediterrean_scrub .or. &
                   ilu .eq. ilu_semi_natural_veg ) then ! rest
            sai_loc = lai_loc
          else
            write( gol, '(" Error win calculation of LAI and SAI")' ) ; call goErr
            write( gol, '(" Landuse index: ", i0, " not supported ")' ) ilu ; call goErr
            TRACEBACK; status=1; return
          end if
        end if

        lai(ix,iy,ilu) = lai_loc
        sai(ix,iy,ilu) = sai_loc
      end do ! landuse classes
    
    end do ! ny
    end do ! nx

    ! ok
    status = 0
    
    
  end subroutine LE_Landuse_Update_LAI_SAI
  
  !-------------------------------------------------------------------
  ! missing: check for data that correspond with a missing deposition path
  !          this data is represented by -999
  !-------------------------------------------------------------------

  logical function missing(x)

    real, intent(in) :: x

    ! bandwidth for checking (in)equalities of floats
    real, parameter :: EPS = 1.0e-5

    missing = (abs(x + 999.) .le. EPS)

  end function missing
  
end module LE_Landuse
