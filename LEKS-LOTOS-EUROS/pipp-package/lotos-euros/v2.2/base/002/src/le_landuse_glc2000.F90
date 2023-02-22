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
    use LE_Grid     , only : ugg
    
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
    allocate( lu_glc2000(ugg%nlon,ugg%nlat,nlu_glc2000) )
    
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

    use GO            , only : pathsep
    use GO            , only : goGetFU
    use LE_Grid       , only : ugg
    use C3PO          , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors


    ! --- in/out ---------------------------------

    character(len=*), intent(in)    ::  fname
    real, intent(inout)             ::  lu_glc2000(ugg%nlon,ugg%nlat,nlu_glc2000)
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/get_lu'

    ! --- local ------------------------------------------

    integer               ::  l, k
    character(len=3)      ::  ext

    integer               ::  varid
    character(len=128)    ::  description
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    integer, allocatable  ::  values_in(:,:)
    character(len=64)     ::  units_in
    
    integer               ::  start_ind(2)
    integer               ::  count_ind(2)

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

        description = 'var_name=glc;units=1'
        ! open file:
        call file_in%Open( trim(fname), status )
        IF_NOTOK_RETURN(status=1)

        ! variable id:
        call file_in%Inq_VarID( trim(description), varid, status )
        IF_NOTOK_RETURN(status=1)
        ! init grid definition
        call file_in%Get_Grid( varid, grid_in, status, &
                                  ugg_to=ugg, start_ind=start_ind, count_ind=count_ind )
        IF_NOTOK_RETURN(status=1)

        ! storage
        allocate( values_in(grid_in%nlon, grid_in%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! read:
        call file_in%Get_Var( trim(description), values_in, units_in, status, &
                               start=start_ind, count=count_ind )
        IF_NOTOK_RETURN(status=1)

        ! convert
        call Grid_Convertors%Ugg_IndexFractions( grid_in, values_in, ugg, lu_glc2000, status )
        IF_NOTOK_RETURN(status=1)

        ! clear:
        deallocate( values_in )
        call grid_in%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! close:
        call file_in%Close( status )
        IF_NOTOK_RETURN(status=1)
       
       
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        write (gol,'("unsupported extension for landuse file : ",a)') trim(ext); call goErr
        IF_NOTOK_RETURN(status=1)

    end select

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
          write (gol,*) 'lu_depac = ', lu_fracs(i, j, :); call goErr
          write (gol,*) 'abs(sum(lu_depac(i,j,:))-1.0) = ', abs(sum(lu_fracs(i,j,:))-1.0); call goErr
          write (gol,*) 'landuse fractions do not count up to 1 for cell (', i, j, ')'; call goErr
          TRACEBACK; status=1; return
        end if
      
      end do   ! j
    end do   ! i
           
    ! ok
    status = 0
           
  end subroutine Translate_LU_to_Depac 


end module LE_Landuse_GLC2000
