!###############################################################################
!
! landuse  -  Read landuse database, convert
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

module LE_LandUse_GLC2000

  use GO, only : gol, goPr, goErr
  use LE_Landuse_Data  

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  LE_LandUse_GLC2000_Init

  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_Landuse_GLC2000'

  ! number of landuse classes:
  integer, parameter  ::  nlu_glc2000 = 23

  integer, parameter  ::  ilu_glc2000_Tree_Cover_broadleaved_evergreen                  =  1
  integer, parameter  ::  ilu_glc2000_Tree_Cover_broadleaved_deciduous_closed           =  2
  integer, parameter  ::  ilu_glc2000_Tree_Cover_broadleaved_deciduous_open             =  3
  integer, parameter  ::  ilu_glc2000_Tree_Cover_needle_leaved_evergreen                =  4
  integer, parameter  ::  ilu_glc2000_Tree_Cover_needle_leaved_deciduous                =  5
  integer, parameter  ::  ilu_glc2000_Tree_Cover_mixed_leaf_type                        =  6
  integer, parameter  ::  ilu_glc2000_Tree_Cover_regularly_flooded_fresh                =  7
  integer, parameter  ::  ilu_glc2000_Tree_Cover_regularly_flooded_saline               =  8
  integer, parameter  ::  ilu_glc2000_Mosaic_Other_natural_vegetation                   =  9
  integer, parameter  ::  ilu_glc2000_Tree_Cover_burnt                                  = 10
  integer, parameter  ::  ilu_glc2000_Shrub_Cover_closed_open_evergreen                 = 11
  integer, parameter  ::  ilu_glc2000_Shrub_Cover_closed_open_deciduous                 = 12
  integer, parameter  ::  ilu_glc2000_Herbaceous_Cover_closed_open                      = 13
  integer, parameter  ::  ilu_glc2000_Sparse_Herbaceous_or_sparse_shrub_cover           = 14
  integer, parameter  ::  ilu_glc2000_Regularly_flooded_shrub_and_herbaceous_cover      = 15
  integer, parameter  ::  ilu_glc2000_Cultivated_and_managed_areas                      = 16
  integer, parameter  ::  ilu_glc2000_Mosaic_Cropland_Tree_Cover_Other_Nat_Veg          = 17
  integer, parameter  ::  ilu_glc2000_Mosaic_Cropland_and_Shrub_and_Herbaceous_cover    = 18
  integer, parameter  ::  ilu_glc2000_Bare_Areas                                        = 19
  integer, parameter  ::  ilu_glc2000_Water_Bodies                                      = 20
  integer, parameter  ::  ilu_glc2000_Snow_and_Ice                                      = 21
  integer, parameter  ::  ilu_glc2000_Artificial_surfaces_and_associated_areas          = 22
  integer, parameter  ::  ilu_glc2000_No_data                                           = 23

  ! --- var --------------------------------------

contains


  ! ========================================================================


  subroutine LE_Landuse_GLC2000_Init( query, status )

    use GO, only : goVarValue
    use dims, only : nx, ny
    
    ! --- in/out ---------------------------------

    character(len=*), intent(in)  ::  query
    integer, intent(out)          ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_LandUse_GLC2000_Init'

    ! --- local ----------------------------------

    character(len=512)  ::  fname
    real, allocatable   ::  lu_glc2000(:,:,:)  ! (nx,ny,nlu_smtk)

    ! --- begin ----------------------------------

    ! * landuse map, vegetation map:    
    allocate( lu_glc2000(nx,ny,nlu_glc2000) )
    
    ! name of landuse file:
    call goVarValue( trim(query), ';', 'file', '=', fname, status )
    IF_NOTOK_RETURN(status=1)

    ! Read land use database:
    call get_lu( trim(fname), lu_glc2000, status )
    IF_NOTOK_RETURN(status=1)

    ! setup DEPAC landuse classes for deposition ;
    ! fill "lu_fracs(:,:,1:nlu)" from "lu_glc200(:,:,nlu_glc2000)" :
    call Translate_LU_to_Depac(lu_glc2000, status)
    IF_NOTOK_RETURN(status=1)
    
    ! fill in vegetations ( no specific information of trees, so all tree valus zero)
    ! common classes identical to landuse type
    veg_fracs(:,:,1:ntree_type)      = 0.0
    ! use depac classes for vegetation parameters:
    veg_fracs(:,:,ntree_type+1:nveg) = lu_fracs(:,:,:)
        
    ! done
    deallocate( lu_glc2000 )       
    ! ok
    status = 0

  end subroutine LE_Landuse_GLC2000_Init


  ! ***


  subroutine get_lu( fname, lu_glc2000, status )

    use GO          , only : pathsep
    use GO          , only : goGetFU
    use Grid        , only : TllGridInfo, Init, Done, IndexFractions

    use dims        , only : runF
    use dims        , only : nx, ny

    use LE_Grid     , only : lli
    use LE_Grid     , only : indomain
    use LE_LandUse_File_nc, only : T_Landuse_File_nc, Init, Done, Get

    ! --- in/out ---------------------------------

    character(len=*), intent(in)    ::  fname
    real, intent(inout)             ::  lu_glc2000(nx,ny,nlu_glc2000)
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/get_lu'

    ! --- local ------------------------------------------

    integer                 ::  l, k
    character(len=3)        ::  ext

    ! landuse and grid:
    real                    ::  west_bound, east_bound, dlon
    real                    ::  south_bound, north_bound, dlat
    integer                 ::  nlon, nlat
    integer, allocatable    ::  landuseX(:,:)
    integer, allocatable    ::  vegetationX(:,:)
    type(TllGridInfo)       ::  lliX

    ! read from nc file:
    type(T_Landuse_File_nc) ::  luf

    ! read from .txt file:
    integer                 ::  fu
    character(len=256)      ::  line
    integer                 ::  iline, ilon, ilat
    real                    ::  cell_lon, cell_lat
    integer                 ::  cell_landuse, cell_forest

    ! --- begin ------------------------------------------

    ! info ...
    write (gol,'("read landuse data base ...")'); call goPr

    ! extract extension:
    l = len_trim(fname)
    k = index( fname(1:l), '.', back=.true. )
    if ( (k == 0) .or. (l-k < 1) .or. (l-k > 3) ) then
      write (gol,'("could not extract 3-character extension of : ",a)') trim(fname); call goErr
      IF_NOTOK_RETURN(status=1)
    end if
    ext = fname(k+1:l)

    ! open using different routines given extension:
    select case ( ext )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'nc' )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! open database:
        call Init( luf, trim(fname), status )
        IF_NOTOK_RETURN(status=1)

        ! extract grid parameters:
        call Get( luf, status, &
                    west_bound=west_bound, east_bound=east_bound, nlon=nlon, &
                    south_bound=south_bound, north_bound=north_bound, nlat=nlat )
        IF_NOTOK_STOP

        ! spacing:
        dlon = (  east_bound -  west_bound ) / nlon
        dlat = ( north_bound - south_bound ) / nlat

        ! setup storage:
        allocate( landuseX(nlon,nlat) )
        
        ! read:
        call Get( luf, status, glc=landuseX )
        IF_NOTOK_STOP

        ! close:
        call Done( luf, status )
        IF_NOTOK_STOP

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        write (gol,'("unsupported extension for landuse file : ",a)') trim(ext); call goErr
        IF_NOTOK_RETURN(status=1)

    end select

         !! dump
         !call nc_dump( 'lux-'//trim(ext)//'.nc', landuseX, 'landuseX', (/'x','y'/), status )
         !IF_NOTOK_STOP

    ! info ...
    write (gol,'("  convert to LE grid ...")'); call goPr
    
    ! setup input grid:
    call Init( lliX, west_bound+0.5*dlon, dlon, nlon, &
                    south_bound+0.5*dlat, dlat, nlat, status )
    IF_NOTOK_STOP

    ! convert from 2D field with indices to 3D field with fraction per cell for each index:
    call IndexFractions( lliX, landuseX, lli, 1, nlu_glc2000, lu_glc2000, status )
    IF_NOTOK_RETURN(status=1)

    ! done:
    call Done( lliX, status )
    IF_NOTOK_STOP

    ! clear:
    if ( allocated( landuseX   ) ) deallocate( landuseX    )

    ! ok
    status = 0

  end subroutine get_lu

  ! ***

  subroutine Translate_LU_to_Depac(lu_glc2000, status )
  
    use dims            , only :  nx, ny

    ! --- in/out ---------------------------------
    
    real, intent(in)              ::  lu_glc2000(nx,ny,nlu_glc2000)
    integer, intent(out)          ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   :: rname = mname//'/Translate_LU_to_Depac'
    
    ! --- local ----------------------------------
    
    integer                       ::  i, j
    
    ! --- begin -------------------------
    ! storage:       
  
    ! Now appoint the classes used for DEPAC.
    do i =1,nx
      do j = 1, ny
        ! assign depac classes from GLC classes 
        lu_fracs(i,j,ilu_grass)              = (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Shrub_Cover_closed_open_evergreen) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Shrub_Cover_closed_open_deciduous) + &
                                                           lu_glc2000(i,j,ilu_glc2000_Herbaceous_Cover_closed_open ) + &
                                               (1.0/2.0) * lu_glc2000(i,j,ilu_glc2000_Sparse_Herbaceous_or_sparse_shrub_cover  ) + &
                                               (1.0/4.0) * lu_glc2000(i,j,ilu_glc2000_Regularly_flooded_shrub_and_herbaceous_cover ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Mosaic_Cropland_and_Shrub_and_Herbaceous_cover )
        
        lu_fracs(i,j,ilu_arable)             = (1.0/2.0) * lu_glc2000(i,j,ilu_glc2000_Cultivated_and_managed_areas ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Mosaic_Cropland_and_Shrub_and_Herbaceous_cover )
        
        lu_fracs(i,j,ilu_permanent_crops)    = (1.0/2.0) * lu_glc2000(i,j,ilu_glc2000_Cultivated_and_managed_areas ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Mosaic_Cropland_Tree_Cover_Other_Nat_Veg ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Mosaic_Cropland_and_Shrub_and_Herbaceous_cover )
        
        lu_fracs(i,j,ilu_coniferous_forest)  =             lu_glc2000(i,j,ilu_glc2000_Tree_Cover_needle_leaved_evergreen ) + &
                                                           lu_glc2000(i,j,ilu_glc2000_Tree_Cover_needle_leaved_deciduous ) + &
                                               (1.0/2.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_mixed_leaf_type ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_regularly_flooded_fresh ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_regularly_flooded_saline ) + &
                                               (1.0/2.0) * lu_glc2000(i,j,ilu_glc2000_Mosaic_Other_natural_vegetation ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_burnt ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Shrub_Cover_closed_open_evergreen) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Shrub_Cover_closed_open_deciduous) + &
                                               (1.0/4.0) * lu_glc2000(i,j,ilu_glc2000_Regularly_flooded_shrub_and_herbaceous_cover ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Mosaic_Cropland_Tree_Cover_Other_Nat_Veg )
        
        lu_fracs(i,j,ilu_deciduous_forest)   =             lu_glc2000(i,j,ilu_glc2000_Tree_Cover_broadleaved_evergreen ) + &
                                                           lu_glc2000(i,j,ilu_glc2000_Tree_Cover_broadleaved_deciduous_closed ) + &
                                                           lu_glc2000(i,j,ilu_glc2000_Tree_Cover_broadleaved_deciduous_open ) + &
                                               (1.0/2.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_mixed_leaf_type ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_regularly_flooded_fresh ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_regularly_flooded_saline ) + &
                                               (1.0/2.0) * lu_glc2000(i,j,ilu_glc2000_Mosaic_Other_natural_vegetation ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_burnt ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Shrub_Cover_closed_open_evergreen) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Shrub_Cover_closed_open_deciduous) + &
                                               (1.0/2.0) * lu_glc2000(i,j,ilu_glc2000_Sparse_Herbaceous_or_sparse_shrub_cover  ) + &
                                               (1.0/4.0) * lu_glc2000(i,j,ilu_glc2000_Regularly_flooded_shrub_and_herbaceous_cover ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Mosaic_Cropland_Tree_Cover_Other_Nat_Veg )
                                               
        lu_fracs(i,j,ilu_water_sea)          = (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_regularly_flooded_fresh ) + &
                                               (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_regularly_flooded_saline ) + &
                                               (1.0/4.0) * lu_glc2000(i,j,ilu_glc2000_Regularly_flooded_shrub_and_herbaceous_cover ) + &
                                                           lu_glc2000(i,j,ilu_glc2000_Water_Bodies)
                                                           
                                                           
        lu_fracs(i,j,ilu_urban)              =             lu_glc2000(i,j,ilu_glc2000_Artificial_surfaces_and_associated_areas)
         
        lu_fracs(i,j,ilu_other)              =             lu_glc2000(i,j,ilu_glc2000_No_data )
         
        lu_fracs(i,j,ilu_desert)             = (1.0/3.0) * lu_glc2000(i,j,ilu_glc2000_Tree_Cover_burnt ) + &
                                                           lu_glc2000(i,j,ilu_glc2000_Bare_Areas )                                                            
      
        lu_fracs(i,j,ilu_ice)                =             lu_glc2000(i,j,ilu_glc2000_Snow_and_Ice)
        lu_fracs(i,j,ilu_savanna)            = 0.0
        lu_fracs(i,j,ilu_tropical_forest)    = 0.0
        lu_fracs(i,j,ilu_water_inland)       = 0.0
        lu_fracs(i,j,ilu_mediterrean_scrub)  = 0.0
        lu_fracs(i,j,ilu_semi_natural_veg)   = 0.0
        if ( ilu_wheat  <= nlu ) lu_fracs(i,j,ilu_wheat ) = 0.0
        if ( ilu_beech  <= nlu ) lu_fracs(i,j,ilu_beech ) = 0.0
        if ( ilu_spruce <= nlu ) lu_fracs(i,j,ilu_spruce) = 0.0
  	    ! check
      	if (abs(sum(lu_fracs(i,j,:))-1.0) > 0.001) then
              print *, 'lu_depac = ', lu_fracs(i, j, :)
              print *, 'abs(sum(lu_depac(i,j,:))-1.0) = ', abs(sum(lu_fracs(i,j,:))-1.0)
              print *, 'landuse fractions do not count up to 1 for cell (', i, j, ')'
              stop
        end if
      
      end do
    end do                                                            
           
    ! ok
    status = 0
           
  end subroutine Translate_LU_to_Depac 


end module LE_Landuse_GLC2000
