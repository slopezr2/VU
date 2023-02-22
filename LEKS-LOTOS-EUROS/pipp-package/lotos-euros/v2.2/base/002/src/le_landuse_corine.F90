!###############################################################################
!
! le_landuse_Corine  -  Read NetCDF landuse file from COrine (2006,2012,2018,...)
!
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

module LE_LandUse_Corine

  use GO, only : gol, goPr, goErr
  
  implicit none


  ! --- in/out -------------------------------------

  private
  
  public  ::  LE_Landuse_Corine_Init
  public  ::  T_Corine_Indices
  
  
  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_Landuse_Corine'
  
  ! --- types -----------------------------
  type T_Corine_Indices
    integer     ::  ilu_cont_urban_fabric
    integer     ::  ilu_disc_urban_fabric
    integer     ::  ilu_ind_units
    integer     ::  ilu_road_rail
    integer     ::  ilu_port
    integer     ::  ilu_airport
    integer     ::  ilu_min_extr
    integer     ::  ilu_dump
    integer     ::  ilu_constr
    integer     ::  ilu_green_urb
    integer     ::  ilu_sport
    integer     ::  ilu_non_irr_arable
    integer     ::  ilu_perm_irr_arable
    integer     ::  ilu_rice
    integer     ::  ilu_vine
    integer     ::  ilu_fruit_trees
    integer     ::  ilu_olives
    integer     ::  ilu_pastures
    integer     ::  ilu_annual_crops
    integer     ::  ilu_compl_cult
    integer     ::  ilu_princ_agri
    integer     ::  ilu_agro_forest
    integer     ::  ilu_broad_leaved_forest
    integer     ::  ilu_conif_forest
    integer     ::  ilu_mixed_forest
    integer     ::  ilu_nat_grass
    integer     ::  ilu_moors_heath
    integer     ::  ilu_sclero_veg
    integer     ::  ilu_trans_woodland
    integer     ::  ilu_beach
    integer     ::  ilu_bare_rocks
    integer     ::  ilu_sparse_veg_area
    integer     ::  ilu_burnt_area
    integer     ::  ilu_glac_snow
    integer     ::  ilu_inland_marsh
    integer     ::  ilu_peat_bogs
    integer     ::  ilu_salt_marsh
    integer     ::  ilu_salines
    integer     ::  ilu_inter_flats
    integer     ::  ilu_water_courses
    integer     ::  ilu_water_bodies
    integer     ::  ilu_coastal_lagoons
    integer     ::  ilu_estuaries
    integer     ::  ilu_sea_ocean
    integer     ::  ilu_NODATA
    integer     ::  ilu_UNCLASS_LAND
    integer     ::  ilu_UNCLASS_WATER
    integer     ::  ilu_UNCLASS
  end type
   
  ! --- interfaces ---------------------------------    

contains
  
  
  ! ==================================================================
  
  
  subroutine LE_Landuse_Corine_Init( rcF, status )
    
    use GO, only : TrcFile, ReadRC
    use MDF, only : MDF_Open, MDF_NETCDF, MDF_READ, MDF_Close
    use MDF, only : MDF_Inq_Dimid, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var

    use LE_grid, only : ugg
    use C3PO   , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors
    
   ! --- in/out -----------------------------
        
    type(TrcFile), intent(in)         ::  rcF
    integer, intent(out)              ::  status
    
    ! --- const ------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Landuse_Corine_Init'
    
    ! --- local ------------------------------
    
    character(len=256)  ::  fname
    logical             ::  exist
    integer             ::  ncid, dimid
    integer             ::  nlu_corine, len_lunames 
    real, allocatable   ::  lufrac_corine_le(:,:,:)
    type(T_Corine_Indices) ::  corine_ind  
    character(len=128), allocatable ::  lunames_corine(:)
    integer             ::  ilu_corine, ilon, ilat
        
    integer               ::  varid
    character(len=128)    ::  description
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    real, allocatable     ::  values_in(:,:,:)
    character(len=64)     ::  units_in
    real                  ::  missing_value

    ! --- begin ------------------------------
    
    ! clc file
    call ReadRc( rcF, 'landuse.corine.file', fname, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found:")'); call goErr
      write (gol,'("  ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! open file:
    call MDF_Open( trim(fname), MDF_NETCDF, MDF_READ, ncid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! number of landuse types, and max length of the names
    call MDF_Inq_DimID( ncid, 'lu_types', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( ncid, dimid, status, length=nlu_corine )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inq_DimID( ncid, 'lu_name_len', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( ncid, dimid, status, length=len_lunames )
    IF_NOTOK_RETURN(status=1)
    ! read Landuse names
    allocate (lunames_corine(nlu_corine) )
    call MDF_Get_StrArr( ncid, 'landuse_types', len_lunames, nlu_corine, lunames_corine, status )
    IF_NOTOK_RETURN(status=1) 
    ! fill in LU indices
    call Fill_indices (lunames_corine, nlu_corine, corine_ind, status )
    IF_NOTOK_RETURN(status=1)
    ! close file
    call MDF_Close( ncid, status )
    IF_NOTOK_RETURN(status=1)
    
     
    ! open file:
    call file_in%Open( trim(fname), status )
    IF_NOTOK_RETURN(status=1)

    ! variable id:
    description='var_name=lu_fraction'
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! init grid definition
    call file_in%Get_Grid( varid, grid_in, status )
    IF_NOTOK_RETURN(status=1)

    ! storage
    allocate( values_in(grid_in%nlon, grid_in%nlat,nlu_corine), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! read:
    call file_in%Get_Var( trim(description), values_in, units_in, status, &
                           start=(/1,1,1/), count=(/grid_in%nlon,grid_in%nlat,nlu_corine/) )
    IF_NOTOK_RETURN(status=1)

    ! target array    
    allocate( lufrac_corine_le(ugg%nlon,ugg%nlat,nlu_corine) )

    ! regrid to LE-domain
    do ilu_corine = 1, nlu_corine
      call Grid_Convertors%Ugg_AreaAver( grid_in, values_in(:,:,ilu_corine), ugg, lufrac_corine_le(:,:,ilu_corine), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! translate to LE-DEPAC classes
    call Translate_LU_to_DEPAC( nlu_corine,lufrac_corine_le(:,:,:), corine_ind, status )
    IF_NOTOK_RETURN(status=1)
    
    deallocate( lufrac_corine_le )
    deallocate( lunames_corine )
    
    ! clear:
    deallocate( values_in )
    ! clear:
    call grid_in%Done(status)
    IF_NOTOK_RETURN(status=1)
    ! close:
    call file_in%Close( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine LE_Landuse_Corine_Init
  
  ! ***    
  
  ! Read 2d character variable in a variable with oversized length

  subroutine MDF_Get_StrArr( hid, varname, nchar, nval, values, status )
  
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var

    ! --- in/out -----------------------

    integer, intent(in)             ::  hid
    character(len=*), intent(in)    ::  varname
    integer, intent(in)             ::  nchar, nval
    character(len=*), intent(out)   ::  values(nval)
    integer, intent(out)            ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MDF_Get_StrArr'
    
    ! --- local -------------------------

    integer                 ::  varid
    character(len=nchar)    ::  tmp_values(nval)
    integer                 ::  ival
    integer                 ::  l

    ! --- begin -------------------------

    ! get variable id given name:      
    call MDF_Inq_VarID( hid, varname, varid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read data:
    call MDF_Get_Var( hid, varid, tmp_values, status )
    IF_NOTOK_RETURN(status=1)
    
    ! copy:
    do ival = 1, nval
      ! seach <nul> characters:
      l = index( tmp_values(ival), char(0) )
      ! found a <nul> character ?
      if ( l < 1 ) then
        ! no <nul>, just copy:
        values(ival) = trim(tmp_values(ival))
      else if ( l == 1 ) then
        ! first character is <nul>, thus empty:
        values(ival) = ''
      else
        ! found a <nul> somewhere in the string; copy part before:
        values(ival) = trim(tmp_values(ival)(1:l-1))
      end if
    end do

    ! ok
    status = 0

  end subroutine MDF_Get_StrArr
  
  ! Fill in in landuse indices for corine
  
  subroutine Fill_Indices( lunames, nlu_names, corine_ind, status )
  
    ! --- in/out ------
    character(len=*), intent(in)      ::  lunames(:)
    integer, intent(in)               ::  nlu_names
    type(T_Corine_Indices), intent(out)  ::  corine_ind
    integer, intent(out)              ::  status
    
    ! --- const -----
    
    character(len=*), parameter ::  rname = mname//'/Fill_Indices'
        
    ! --- local -----
    integer           ::  ilu

    ! --- begin -----
    
    do ilu = 1, nlu_names
      write( gol, '("ilu: ", i0, " Name: ", a)') ilu, trim(lunames(ilu)) ; call goPr
      
      if ( trim(lunames(ilu) ) == 'Continuous urban fabric' ) then
        corine_ind%ilu_cont_urban_fabric = ilu
      else if ( trim(lunames(ilu) ) == 'Discontinuous urban fabric') then
        corine_ind%ilu_disc_urban_fabric = ilu
      else if ( trim(lunames(ilu) ) == 'Industrial or commercial units') then
        corine_ind%ilu_ind_units = ilu
      else if ( trim(lunames(ilu) ) == 'Road and rail networks and associated land' ) then
        corine_ind%ilu_road_rail = ilu
      else if ( trim(lunames(ilu) ) == 'Port areas' ) then
        corine_ind%ilu_port = ilu
      else if ( trim(lunames(ilu) ) == 'Airports' ) then
        corine_ind%ilu_airport = ilu
      else if ( trim(lunames(ilu) ) == 'Mineral extraction sites' ) then
        corine_ind%ilu_min_extr = ilu
      else if ( trim(lunames(ilu) ) == 'Dump sites' ) then
        corine_ind%ilu_dump = ilu
      else if ( trim(lunames(ilu) ) == 'Construction sites' ) then
        corine_ind%ilu_constr = ilu
      else if ( trim(lunames(ilu) ) == 'Green urban areas' ) then
        corine_ind%ilu_green_urb = ilu
      else if ( trim(lunames(ilu) ) == 'Sport and leisure facilities' ) then
        corine_ind%ilu_sport = ilu
      else if ( trim(lunames(ilu) ) == 'Non-irrigated arable land' ) then
        corine_ind%ilu_non_irr_arable = ilu
      else if ( trim(lunames(ilu) ) == 'Permanently irrigated land' ) then
        corine_ind%ilu_perm_irr_arable = ilu
      else if ( trim(lunames(ilu) ) == 'Rice fields' ) then
        corine_ind%ilu_rice = ilu
      else if ( trim(lunames(ilu) ) == 'Vineyards' ) then
        corine_ind%ilu_vine = ilu
      else if ( trim(lunames(ilu) ) == 'Fruit trees and berry plantations' ) then
        corine_ind%ilu_fruit_trees = ilu
      else if ( trim(lunames(ilu) ) == 'Olive groves' ) then
        corine_ind%ilu_olives = ilu
      else if ( trim(lunames(ilu) ) == 'Pastures' ) then
        corine_ind%ilu_pastures = ilu
      else if ( trim(lunames(ilu) ) == 'Annual crops associated with permanent crops' ) then
        corine_ind%ilu_annual_crops = ilu
      else if ( trim(lunames(ilu) ) == 'Complex cultivation patterns' ) then
        corine_ind%ilu_compl_cult = ilu
      else if ( trim(lunames(ilu) ) == 'Land principally occupied by agriculture with significant areas of natural vegetation' ) then
        corine_ind%ilu_princ_agri = ilu
      else if ( trim(lunames(ilu) ) == 'Land principally occupied by agriculture, with significant areas of natural vegetation' ) then
        corine_ind%ilu_princ_agri = ilu
      else if ( trim(lunames(ilu) ) == 'Agro-forestry areas' ) then
        corine_ind%ilu_agro_forest = ilu
      else if ( trim(lunames(ilu) ) == 'Broad-leaved forest' ) then
        corine_ind%ilu_broad_leaved_forest = ilu
      else if ( trim(lunames(ilu) ) == 'Coniferous forest' ) then
        corine_ind%ilu_conif_forest = ilu
      else if ( trim(lunames(ilu) ) == 'Mixed forest' ) then
        corine_ind%ilu_mixed_forest = ilu
      else if ( trim(lunames(ilu) ) == 'Natural grasslands' ) then
        corine_ind%ilu_nat_grass = ilu
      else if ( trim(lunames(ilu) ) == 'Moors and heathland' ) then
        corine_ind%ilu_moors_heath = ilu
      else if ( trim(lunames(ilu) ) == 'Sclerophyllous vegetation' ) then
        corine_ind%ilu_sclero_veg = ilu
      else if ( trim(lunames(ilu) ) == 'Transitional woodland-shrub' ) then
        corine_ind%ilu_trans_woodland = ilu
      else if ( trim(lunames(ilu) ) == 'Beaches - dunes - sands' ) then
        corine_ind%ilu_beach = ilu
      else if ( trim(lunames(ilu) ) == 'Beaches, dunes, sands' ) then
        corine_ind%ilu_beach = ilu
      else if ( trim(lunames(ilu) ) == 'Bare rocks' ) then
        corine_ind%ilu_bare_rocks = ilu
      else if ( trim(lunames(ilu) ) == 'Sparsely vegetated areas' ) then
        corine_ind%ilu_sparse_veg_area = ilu
      else if ( trim(lunames(ilu) ) == 'Burnt areas' ) then
        corine_ind%ilu_burnt_area = ilu
      else if ( trim(lunames(ilu) ) == 'Glaciers and perpetual snow' ) then
        corine_ind%ilu_glac_snow = ilu
      else if ( trim(lunames(ilu) ) == 'Inland marshes' ) then
        corine_ind%ilu_inland_marsh = ilu
      else if ( trim(lunames(ilu) ) == 'Peat bogs' ) then
        corine_ind%ilu_peat_bogs = ilu
      else if ( trim(lunames(ilu) ) == 'Salt marshes' ) then
        corine_ind%ilu_salt_marsh = ilu
      else if ( trim(lunames(ilu) ) == 'Salines' ) then
        corine_ind%ilu_salines = ilu
      else if ( trim(lunames(ilu) ) == 'Intertidal flats' ) then
        corine_ind%ilu_inter_flats = ilu
      else if ( trim(lunames(ilu) ) == 'Water courses' ) then
        corine_ind%ilu_water_courses = ilu
      else if ( trim(lunames(ilu) ) == 'Water bodies' ) then
        corine_ind%ilu_water_bodies = ilu
      else if ( trim(lunames(ilu) ) == 'Coastal lagoons' ) then
        corine_ind%ilu_coastal_lagoons = ilu
      else if ( trim(lunames(ilu) ) == 'Estuaries' ) then
        corine_ind%ilu_estuaries = ilu
      else if ( trim(lunames(ilu) ) == 'Sea and ocean' ) then
        corine_ind%ilu_sea_ocean = ilu
      else if ( trim(lunames(ilu) ) == 'NODATA' ) then
        corine_ind%ilu_NODATA = ilu
      else if ( trim(lunames(ilu) ) == 'UNCLASSIFIED LAND SURFACE' ) then
        corine_ind%ilu_UNCLASS_LAND = ilu
      else if ( trim(lunames(ilu) ) == 'UNCLASSIFIED WATER BODIES' ) then
        corine_ind%ilu_UNCLASS_WATER = ilu
      else if ( trim(lunames(ilu) ) == 'UNCLASSIFIED' ) then
        corine_ind%ilu_UNCLASS = ilu
      else
        write( gol, '("unknown Index for Landuse type: ", a) ' ) trim(lunames(ilu)) ; call goErr
        TRACEBACK; status=1; return
      
      end if
    end do

    ! ok 
    status = 0
    
  end subroutine Fill_Indices
  
  
  subroutine Translate_LU_to_DEPAC( nlu_corine, lufrac_corine, corine_ind, status )
    
    use dims, only : nx, ny
    
    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_Data, only : ilu_grass, ilu_arable,ilu_permanent_crops
    use LE_Landuse_Data, only : ilu_coniferous_forest, ilu_deciduous_forest
    use LE_Landuse_Data, only : ilu_water_sea, ilu_urban, ilu_other, ilu_desert
    use LE_Landuse_Data, only : ilu_ice, ilu_savanna, ilu_tropical_forest
    use LE_Landuse_Data, only : ilu_water_inland, ilu_mediterrean_scrub
    use LE_Landuse_Data, only : ilu_semi_natural_veg, ilu_wheat, ilu_beech
    use LE_Landuse_Data, only : ilu_spruce
    
    use LE_Landuse_Data, only : with_ozone_specials
    use file_nc, only : nc_dump
    ! --- in/out ---
    integer, intent(in)   ::  nlu_corine
    real, intent(in)      ::  lufrac_corine(nx,ny,nlu_corine)
    type(T_Corine_Indices), intent(in) :: corine_ind
    integer, intent(out)  ::  status
    
    ! --- const -----

    character(len=*), parameter ::  rname = mname//'/Translate_LU_to_DEPAC'

    ! --- local ---
    integer       ::  ix,iy
    
    ! --- begin ---
    
    do ix = 1, nx
    do iy = 1, ny
      if ( abs( sum(lufrac_corine(ix,iy,:)) - 1.0 ) < 0.001 ) then 
        if ( lufrac_corine(ix,iy,corine_ind%ilu_NODATA) > 0.0 .or. &
             lufrac_corine(ix,iy,corine_ind%ilu_UNCLASS_LAND) > 0.0 .or. &
             lufrac_corine(ix,iy,corine_ind%ilu_UNCLASS_WATER) > 0.0 ) then
             ! not all data in this cell of corine land cover are valid,
             ! For example in Corine 2012 land area of undefined countries (RUS,UKR,GRC,...),
             ! are defined as Unclassified Water, so for all Undefined classes:
             !  --- DO not change depac fractions ---
             !
             cycle
        else 
          lu_fracs(ix,iy,ilu_grass)              = lufrac_corine(ix,iy,corine_ind%ilu_green_urb) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_sport ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_pastures ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_peat_bogs )

          lu_fracs(ix,iy,ilu_arable)             = lufrac_corine(ix,iy,corine_ind%ilu_non_irr_arable ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_perm_irr_arable ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_rice )+ &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_annual_crops )+ &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_compl_cult ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_princ_agri ) 


          lu_fracs(ix,iy,ilu_permanent_crops)    = lufrac_corine(ix,iy,corine_ind%ilu_vine ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_fruit_trees ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_olives ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_annual_crops ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_agro_forest ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_princ_agri )

          lu_fracs(ix,iy,ilu_coniferous_forest)  = lufrac_corine(ix,iy,corine_ind%ilu_conif_forest ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_mixed_forest ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_trans_woodland )

          lu_fracs(ix,iy,ilu_deciduous_forest)   = lufrac_corine(ix,iy,corine_ind%ilu_broad_leaved_forest ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_mixed_forest ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_inland_marsh )+ &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_trans_woodland )

          lu_fracs(ix,iy,ilu_water_sea)          = lufrac_corine(ix,iy,corine_ind%ilu_glac_snow ) + &                                             
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_salt_marsh ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_salines ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_inter_flats ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_coastal_lagoons ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_sea_ocean )

          lu_fracs(ix,iy,ilu_urban)              = lufrac_corine(ix,iy,corine_ind%ilu_cont_urban_fabric ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_disc_urban_fabric ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_ind_units ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_road_rail ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_port ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_airport ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_min_extr ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_dump ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_constr )

          lu_fracs(ix,iy,ilu_other)              = lufrac_corine(ix,iy,corine_ind%ilu_moors_heath ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_sclero_veg ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_nat_grass ) 

          lu_fracs(ix,iy,ilu_desert)             = lufrac_corine(ix,iy,corine_ind%ilu_beach ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_bare_rocks ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_sparse_veg_area ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_burnt_area ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_salt_marsh ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_salines ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_inter_flats )

          lu_fracs(ix,iy,ilu_ice)                = 0.0 
          lu_fracs(ix,iy,ilu_savanna)            = 0.0
          lu_fracs(ix,iy,ilu_tropical_forest)    = 0.0
          lu_fracs(ix,iy,ilu_water_inland)       = lufrac_corine(ix,iy,corine_ind%ilu_water_courses ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_inland_marsh ) + &
                                             0.5 * lufrac_corine(ix,iy,corine_ind%ilu_peat_bogs ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_water_bodies ) + &
                                                   lufrac_corine(ix,iy,corine_ind%ilu_estuaries )
          lu_fracs(ix,iy,ilu_mediterrean_scrub)  = 0.0                                       
          lu_fracs(ix,iy,ilu_semi_natural_veg)   = 0.0
          if ( with_ozone_specials ) then           
            lu_fracs(ix,iy,ilu_wheat)            = 0.0
            lu_fracs(ix,iy,ilu_beech)            = 0.0
            lu_fracs(ix,iy,ilu_spruce)           = 0.0
          end if
        end if ! all data valid?
      end if ! data available for corine?
               
    end do ! ny    
    end do ! nx
    
    call nc_dump( 'lufracs_corine.nc', lufrac_corine, 'lus', (/'nx ','ny ','nlu'/), status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Translate_LU_to_DEPAC
 
end module LE_LandUse_Corine

