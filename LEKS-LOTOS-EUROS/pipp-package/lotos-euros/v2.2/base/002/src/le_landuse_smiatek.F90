!###############################################################################
!
! landuse  -  Read landuse database from smiatek
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

module LE_LandUse_Smiatek

  use GO, only : gol, goPr, goErr
  
  use LE_Landuse_Data

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  LE_LandUse_Smiatek_Init
  
  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_Landuse_Smiatek'

  ! number of landuse classes:
  integer, parameter  ::  nlu_smtk = 13
  !
  ! the  lotos landuse database consists of the following "Smiatek" clases:
  integer, parameter  ::  ilu_smtk_urban_areas        =  1
  integer, parameter  ::  ilu_smtk_agriculture        =  2
  integer, parameter  ::  ilu_smtk_grassland          =  3
  integer, parameter  ::  ilu_smtk_deciduous_forest   =  4
  integer, parameter  ::  ilu_smtk_coniferous_forest  =  5
  integer, parameter  ::  ilu_smtk_mixed_forest       =  6
  integer, parameter  ::  ilu_smtk_water              =  7
  integer, parameter  ::  ilu_smtk_marsh_or_wetland   =  8
  integer, parameter  ::  ilu_smtk_sand_or_bare_rocks =  9
  integer, parameter  ::  ilu_smtk_tundra             = 10
  integer, parameter  ::  ilu_smtk_permanent_ice      = 11
  integer, parameter  ::  ilu_smtk_tropical_forest    = 12
  integer, parameter  ::  ilu_smtk_woodland_scrub     = 13
  
  ! --- var --------------------------------------

contains


  ! ========================================================================


  subroutine LE_Landuse_Smiatek_Init( query, status )

    use GO, only : goVarValue
    use dims, only : nx, ny
    
    ! --- in/out ---------------------------------

    character(len=*), intent(in)  ::  query
    integer, intent(out)          ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Landuse_Smiatek_Init'
    ! --- local ----------------------------------

    character(len=512)  ::  fname
    real, allocatable     ::  lu_smtk(:,:,:)              ! (nx,ny,nlu_smtk)
    real, allocatable     ::  veg_smtk(:,:,:)             ! (nx,ny,nveg)
    real, allocatable     ::  veg_smtk_to_depac(:,:,:)    ! (nx,ny,nlu)

    ! --- begin ----------------------------------

    ! storage:
    allocate( lu_smtk(nx,ny,nlu_smtk), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( veg_smtk(nx,ny,nveg), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( veg_smtk_to_depac(nx,ny,nlu), stat=status )
    IF_NOTOK_RETURN(status=1)
       
    ! * landuse map, vegetation map:    
    lu_smtk  = 0.0
    veg_smtk = 0.0
    veg_smtk_to_depac = 0.0
    
    ! name of landuse file:
    call goVarValue( trim(query), ';', 'file', '=', fname, status )
    IF_NOTOK_RETURN(status=1)

    ! Read land use database:
    call get_lu( trim(fname), lu_smtk, veg_smtk, status )
    IF_NOTOK_RETURN(status=1)
   
    ! setup DEPAC landuse classes for deposition:
    call Translate_LU_to_Depac( lu_smtk, veg_smtk, veg_smtk_to_depac, status)
    IF_NOTOK_RETURN(status=1)
    
    ! fill in vegetation data in common array
    ! For Smiatek this is identical
    veg_fracs(:,:,1:ntree_type)      = veg_smtk(:,:,1:ntree_type)
    veg_fracs(:,:,ntree_type+1:nveg) = veg_smtk_to_depac(:,:,:)

    ! clear:
    deallocate( lu_smtk, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( veg_smtk, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( veg_smtk_to_depac, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Landuse_Smiatek_Init


  ! ***


  subroutine get_lu( fname, lu_smtk, veg_smtk, status )

    use GO          , only : pathsep
    use GO          , only : goGetFU
    use LE_Grid     , only : ugg
    use C3PO   , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors
    
    use file_nc
    
    ! --- in/out ---------------------------------

    character(len=*), intent(in)    ::  fname
    real, intent(inout)             ::  lu_smtk(:,:,:)  ! (ugg%nlon,ugg%nlat,nlu_smtk)
    real, intent(inout)             ::  veg_smtk(:,:,:) ! (ugg%nlon,ugg%nlat,nveg)
    integer, intent(out)            ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/get_lu'

    ! --- local ------------------------------------------

    ! landuse and grid:
    integer, allocatable    ::  landuseX(:,:)
    integer, allocatable    ::  vegetationX(:,:)

    integer               ::  varid
    character(len=128)    ::  description
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    character(len=64)     ::  units_in
    
    integer               ::  start_ind(2)  ! start indices of input file to match target grid
    integer               ::  count_ind(2)  ! number of cells of input file to match target grid
    
    ! --- begin ------------------------------------------

    ! info ...
    write (gol,'("read landuse data base ...")'); call goPr
    
    ! check ...
    if ( any( shape(lu_smtk) /= (/ugg%nlon,ugg%nlat,nlu_smtk/) ) ) then
      write (gol,'("shape of lu_smtk is (",i0,2(",",i0),") while ugg has shape (",i0,",",i0,") and nlu_smtk=",i0)') &
                      shape(lu_smtk), ugg%nlon, ugg%nlat, nlu_smtk; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( any( shape(veg_smtk) /= (/ugg%nlon,ugg%nlat,nveg/) ) ) then
      write (gol,'("shape of lu_smtk is (",i0,2(",",i0),") while ugg has shape (",i0,",",i0,") and nveg=",i0)') &
                      shape(veg_smtk), ugg%nlon, ugg%nlat, nveg; call goErr
      TRACEBACK; status=1; return
    end if

    ! info ...
    write (gol,'("  open ",a," ...")') trim(fname); call goPr
    
    ! open file:
    call file_in%Open( trim(fname), status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    write (gol,'("  read grid definition ...")'); call goPr

    ! variable id:
    description = 'var_name=landuse'
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! init grid definition for source data,
    ! limit to target region to save memory:
    call file_in%Get_Grid( varid, grid_in, status, &
                            ugg_to=ugg, start_ind=start_ind, count_ind=count_ind )
    IF_NOTOK_RETURN(status=1)

    ! storage
    allocate( landuseX(grid_in%nlon, grid_in%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( vegetationX(grid_in%nlon, grid_in%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    write (gol,'("  read landuse ...")'); call goPr
    ! read:
    description = 'var_name=landuse'
    call file_in%Get_Var( trim(description), landuseX, units_in, status, &
                           start=start_ind, count=count_ind )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("    convert ...")'); call goPr
    ! convert from 2D field with indices to 3D field with fraction per cell for each index:
    call Grid_Convertors%UGG_IndexFractions( grid_in, landuseX, ugg, lu_smtk, status, &
                                                         clear=.true. )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    write (gol,'("  read vegetation ...")'); call goPr
    ! read:
    description = 'var_name=vegetation'
    call file_in%Get_Var( trim(description), vegetationX, units_in, status, &
                           start=start_ind, count=count_ind )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    write (gol,'("    convert ...")'); call goPr
    ! convert from 2D field with indices to 3D field with fraction per cell for each index:
    call Grid_Convertors%UGG_IndexFractions( grid_in, vegetationX, ugg, veg_smtk, status, &
                                                         clear=.true. )
    IF_NOTOK_RETURN(status=1)
    
    !! testing ...
    !call nc_dump( 'landuseX.nc', landuseX, 'landuseX', (/'longitude','latitude '/), status )
    !IF_NOTOK_RETURN(status=1)
    !call nc_dump( 'landuse.nc', lu_smtk, 'landuse', (/'longitude','latitude ','lu       '/), status)
    !IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( landuseX, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( vegetationX, stat=status )
    IF_NOTOK_RETURN(status=1)
    call grid_in%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! close:
    call file_in%Close( status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    write (gol,'("  ok.")'); call goPr

    ! ok
    status = 0

  end subroutine get_lu


  ! ***


  subroutine Translate_LU_to_Depac(lu_smtk, veg_smtk, veg_smtk_to_depac, status )
  
    use dims            , only :  nx, ny

    ! --- in/out ---------------------------------
    
    real, intent(in)              ::  lu_smtk(nx,ny,nlu_smtk)
    real, intent(in)              ::  veg_smtk(nx,ny,nveg)
    real, intent(inout)           ::  veg_smtk_to_depac(nx,ny,nlu)
    integer, intent(out)          ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   :: rname = mname//'/Translate_LU_to_Depac'
    
    ! --- local ----------------------------------
    
    integer                       ::  i, j, k
    
    ! --- begin ------------
    
    ! storage:        
    do i =1,nx
      do j = 1, ny
      
        ! assign depac classes from smiatek classes
        lu_fracs(i,j,ilu_grass)              =       lu_smtk(i,j,ilu_smtk_grassland)
        lu_fracs(i,j,ilu_arable)             = 0.9 * lu_smtk(i,j,ilu_smtk_agriculture)
        lu_fracs(i,j,ilu_permanent_crops)    = 0.1 * lu_smtk(i,j,ilu_smtk_agriculture)
        lu_fracs(i,j,ilu_coniferous_forest)  =       lu_smtk(i,j,ilu_smtk_coniferous_forest) + &
                                               0.5 * lu_smtk(i,j,ilu_smtk_mixed_forest) 
        lu_fracs(i,j,ilu_deciduous_forest)   =       lu_smtk(i,j,ilu_smtk_deciduous_forest) + &
                                               0.5 * lu_smtk(i,j,ilu_smtk_mixed_forest) + &
                                               0.5 * lu_smtk(i,j,ilu_smtk_marsh_or_wetland) + &
                                                     lu_smtk(i,j,ilu_smtk_woodland_scrub) + &
                                                     lu_smtk(i,j,ilu_smtk_tropical_forest)                                                                                               
        lu_fracs(i,j,ilu_water_sea)          =       lu_smtk(i,j,ilu_smtk_water) + &
                                               0.5 * lu_smtk(i,j,ilu_smtk_marsh_or_wetland) + &
                                                     lu_smtk(i,j,ilu_smtk_permanent_ice)   
        lu_fracs(i,j,ilu_urban)              =       lu_smtk(i,j,ilu_smtk_urban_areas)
        lu_fracs(i,j,ilu_other)              =       lu_smtk(i,j,ilu_smtk_tundra)
        lu_fracs(i,j,ilu_desert)             =       lu_smtk(i,j,ilu_smtk_sand_or_bare_rocks)
        lu_fracs(i,j,ilu_ice)                =       0.0  
!        lu_fracs(i,j,ilu_ice)                =       lu_smtk(i,j,ilu_smtk_permanent_ice)   
        lu_fracs(i,j,ilu_savanna)            =       0.0
        lu_fracs(i,j,ilu_tropical_forest)    =       0.0
!        lu_fracs(i,j,ilu_tropical_forest)    =       lu_smtk(i,j,ilu_smtk_tropical_forest)
        lu_fracs(i,j,ilu_water_inland)       =       0.0
        lu_fracs(i,j,ilu_mediterrean_scrub)  =       0.0
        lu_fracs(i,j,ilu_semi_natural_veg)   =       0.0
        if ( with_ozone_specials ) then
          lu_fracs(i,j,ilu_wheat)            =       0.0
          lu_fracs(i,j,ilu_beech)            =       0.0
          lu_fracs(i,j,ilu_spruce)           =       0.0
        end if
        
        ! check 
        if ( abs(sum(lu_fracs(i,j,:))-1.0) > 0.001 ) then
          write (gol,'("landuse fractions do not count up to 1 for cell (",i0,",",i0,")")') i, j; call goErr
          write (gol,'("  sum(lu_smtk (:)) = ",f10.4)') sum(lu_smtk (i,j,:)); call goErr
          write (gol,'("  sum(lu_fracs(:)) = ",f10.4)') sum(lu_fracs(i,j,:)); call goErr
          do k = lbound(lu_fracs,3), ubound(lu_fracs,3)
            write (gol,'("  lu_fracs(",i3,") = ",f10.4)') k, lu_fracs(i, j, k); call goErr
          end do
          TRACEBACK; status=1; return
        end if      
        
        ! translation of common vegations classes
        veg_smtk_to_depac(i,j,ilu_grass)              =       veg_smtk(i,j,ntree_type+ilu_smtk_grassland)
        veg_smtk_to_depac(i,j,ilu_arable)             = 0.9 * veg_smtk(i,j,ntree_type+ilu_smtk_agriculture)
        veg_smtk_to_depac(i,j,ilu_permanent_crops)    = 0.1 * veg_smtk(i,j,ntree_type+ilu_smtk_agriculture)
        veg_smtk_to_depac(i,j,ilu_coniferous_forest)  =       veg_smtk(i,j,ntree_type+ilu_smtk_coniferous_forest) + &
                                                        0.5 * veg_smtk(i,j,ntree_type+ilu_smtk_mixed_forest) 
        veg_smtk_to_depac(i,j,ilu_deciduous_forest)   =       veg_smtk(i,j,ntree_type+ilu_smtk_deciduous_forest) + &
                                                        0.5 * veg_smtk(i,j,ntree_type+ilu_smtk_mixed_forest) + &
                                                        0.5 * veg_smtk(i,j,ntree_type+ilu_smtk_marsh_or_wetland) + &
                                                              veg_smtk(i,j,ntree_type+ilu_smtk_woodland_scrub) + &
                                                              veg_smtk(i,j,ntree_type+ilu_smtk_tropical_forest)                                                                                                   
        veg_smtk_to_depac(i,j,ilu_water_sea)          =       veg_smtk(i,j,ntree_type+ilu_smtk_water) + &
                                                        0.5 * veg_smtk(i,j,ntree_type+ilu_smtk_marsh_or_wetland) + &
                                                              veg_smtk(i,j,ntree_type+ilu_smtk_permanent_ice)
        veg_smtk_to_depac(i,j,ilu_urban)              =       veg_smtk(i,j,ntree_type+ilu_smtk_urban_areas)
        veg_smtk_to_depac(i,j,ilu_other)              =       veg_smtk(i,j,ntree_type+ilu_smtk_tundra)
        veg_smtk_to_depac(i,j,ilu_desert)             =       veg_smtk(i,j,ntree_type+ilu_smtk_sand_or_bare_rocks)
        veg_smtk_to_depac(i,j,ilu_ice)                =       0.0
!        veg_smtk_to_depac(i,j,ilu_ice)                =       veg_smtk(i,j,200+ilu_smtk_permanent_ice)   
        veg_smtk_to_depac(i,j,ilu_savanna)            =       0.0
        veg_smtk_to_depac(i,j,ilu_tropical_forest)    =       0.0
!        veg_smtk_to_depac(i,j,ilu_tropical_forest)    =       veg_smtk(i,j,200+ilu_smtk_tropical_forest)
        veg_smtk_to_depac(i,j,ilu_water_inland)       =       0.0
        veg_smtk_to_depac(i,j,ilu_mediterrean_scrub)  =       0.0
        veg_smtk_to_depac(i,j,ilu_semi_natural_veg)   =       0.0
        if ( with_ozone_specials ) then
          veg_smtk_to_depac(i,j,ilu_wheat)              =       0.0
          veg_smtk_to_depac(i,j,ilu_beech)              =       0.0
          veg_smtk_to_depac(i,j,ilu_spruce)             =       0.0
        end if

      end do
    end do
    
    ! ok
    status = 0
           
  end subroutine Translate_LU_to_Depac 


end module LE_Landuse_Smiatek
