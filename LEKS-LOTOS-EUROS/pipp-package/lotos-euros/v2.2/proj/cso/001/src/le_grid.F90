!###############################################################################
!
! NAME
!   LE_Grid  -  LOTOS-EUROS grid stuff
!
! DESCRIPTION
!
!   Grid definitions.
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

module le_grid

  use GO  , only : gol, goPr, goErr
  use GO  , only : T_Domains
  
  use C3PO, only : T_Grid_Ugg
  use C3PO, only : CRNR_LL, CRNR_UL, CRNR_UR, CRNR_LR
  use C3PO, only : EDGE_LEFT, EDGE_UPPER, EDGE_RIGHT, EDGE_LOWER
  
  implicit none
  
  
  ! --- in/out ------------------------------------
  
  private
  
  public  ::  LE_Grid_Init, LE_Grid_Done
  
  public  ::  dom
  
  public  ::  glb_ugg
  public  ::  glb_ugg_bnds
  public  ::  glb_ugg_crnr  
  public  ::  glb_ugg_crnr_bnds
  public  ::  glb_ugg_ustag
  public  ::  glb_ugg_vstag
  
  public  ::  ugg
  public  ::  ugg_bnds
  public  ::  ugg_crnr  
  public  ::  ugg_crnr_bnds
  public  ::  ugg_west, ugg_east
  public  ::  ugg_north, ugg_south
  public  ::  ugg_ustag
  public  ::  ugg_vstag
  
  public  ::  CRNR_LL, CRNR_UL, CRNR_UR, CRNR_LR
  public  ::  EDGE_LEFT, EDGE_UPPER, EDGE_RIGHT, EDGE_LOWER
  
  !public  :: debug_i, debug_j

  
  ! --- const -------------------------------------
    
  character(len=*), parameter   ::  mname = 'le_grid'
  
  
  ! --- var ---------------------------------------

  ! global grids:
  type(T_Grid_Ugg), target      ::  glb_ugg
  type(T_Grid_Ugg), target      ::  glb_ugg_bnds
  type(T_Grid_Ugg), target      ::  glb_ugg_crnr
  type(T_Grid_Ugg), target      ::  glb_ugg_crnr_bnds
  type(T_Grid_Ugg), target      ::  glb_ugg_ustag
  type(T_Grid_Ugg), target      ::  glb_ugg_vstag
  
  ! local grids:
  type(T_Grid_Ugg), target      ::  ugg
  type(T_Grid_Ugg), target      ::  ugg_bnds
  type(T_Grid_Ugg), target      ::  ugg_crnr
  type(T_Grid_Ugg), target      ::  ugg_crnr_bnds
  type(T_Grid_Ugg), target      ::  ugg_west, ugg_east
  type(T_Grid_Ugg), target      ::  ugg_south, ugg_north
  type(T_Grid_Ugg), target      ::  ugg_ustag
  type(T_Grid_Ugg), target      ::  ugg_vstag
  
  ! local domain definition:
  type(T_Domains)               ::  dom
  
  !! testing ...
  !integer         ::  debug_i, debug_j
  
  
contains


  ! ===================================================================
  ! ===
  ! === universial grids
  ! ===
  ! ===================================================================
  
  
  subroutine LE_Grid_Init( rcF, status )
  
    use GO  , only : TrcFile, ReadRc
    use GO  , only : goc
    use C3PO, only : T_File_Ugg
    use Dims, only : nx,ny,nz
    
    ! in/out  ------------------------------------
    
    type(TrcFile), intent(in)         ::  rcF
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Grid_Init'
    
    ! --- local ----------------------------------

    character(len=32)    ::  grid_type
    
    integer              ::  gnx, gny
    integer              ::  ndomx, ndomy
    integer              ::  shp(2)
    integer              ::  off(2)
    integer              ::  i1, i2, j1, j2
    
    real, allocatable    ::  lons(:)
    real, allocatable    ::  lons_bnds(:,:)
    real, allocatable    ::  lats(:)
    real, allocatable    ::  lats_bnds(:,:)
    
    integer              ::  ix, iy
    integer              ::  nh
    real                 ::  dlon, dlat
    real                 ::  westb, southb

    type(T_File_Ugg)     ::  file_in
    character(len=1024)  ::  filename
    character(len=128)   ::  description
    integer              ::  subset(4)
    integer              ::  varid
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (gol,'(a,": define grids")') rname; call goPr
    
    ! number of levels:
    call ReadRc( rcF, 'levels.nz', nz, status )
    IF_NOTOK_RETURN(status=1)

    ! grid type:
    call ReadRc( rcF, 'grid.type', grid_type, status )
    IF_NOTOK_RETURN(status=1)
    
    ! info ..
    write (gol,'(a,":   grid type: ",a)') rname, trim(grid_type); call goPr
    
    ! switch:
    select case ( trim(grid_type) )
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! regular cartesian grid:
      case ( 'cartesian' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        ! read corner:
        call ReadRc( rcF, 'grid.west' , westb , status )
        IF_NOTOK_RETURN(status=1)
        call ReadRc( rcF, 'grid.south', southb, status )
        IF_NOTOK_RETURN(status=1)

        ! resolutions:
        call ReadRc( rcF, 'grid.dlon', dlon, status )
        IF_NOTOK_RETURN(status=1)
        call ReadRc( rcF, 'grid.dlat', dlat, status )
        IF_NOTOK_RETURN(status=1)

        ! global grid size:
        call ReadRc( rcF, 'grid.nx', gnx, status )
        IF_NOTOK_RETURN(status=1)
        call ReadRc( rcF, 'grid.ny', gny, status )
        IF_NOTOK_RETURN(status=1)
        
        ! decomposition:
        call ReadRc( rcF, 'domains.x', ndomx, status )
        IF_NOTOK_RETURN(status=1)
        call ReadRc( rcF, 'domains.y', ndomy, status )
        IF_NOTOK_RETURN(status=1)
        ! define data bounds for local domain ; provide:
        !  - decomposition shape
        !  - global lower and upper bounds ;
        call dom%Init( (/ndomx,ndomy/), (/1,1/), (/gnx,gny/), status )
        IF_NOTOK_RETURN(status=1)
        !! info ...
        !call dom%Show( status )
        !IF_NOTOK_RETURN(status=1)
        
        ! number of halo cells
        !call ReadRc( rcF, 'grid.halocells', nh, status )
        !IF_NOTOK_RETURN(status=1)
        !~ default ...
        nh = 2

        ! * global grids
        
        ! get local shape:
        call dom%Get( status, glb_shp=shp )
        IF_NOTOK_RETURN(status=1)
        ! copy:
        nx = shp(1)
        ny = shp(2)
        ! info ...
        write (gol,'("global shape : (",i0,",",i0,")")') nx, ny; call goPr

        ! storage for cell centers:
        allocate( lons(1-nh:nx+nh), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( lats(1-nh:ny+nh), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! storage for cell bounds:
        allocate( lons_bnds(2,1-nh:nx+nh), stat=status )
        IF_NOTOK_RETURN(status=1)
        allocate( lats_bnds(2,1-nh:ny+nh), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! loop over lon cells:
        do ix = 1-nh, nx+nh
          ! fill center:
          lons(ix) = westb+0.5*dlon + (ix-1)*dlon
          ! fill bounds:
          lons_bnds(1,ix) = westb + (ix-1)*dlon
          lons_bnds(2,ix) = westb + (ix  )*dlon
        end do

        ! loop over lat cells:    
        do iy = 1-nh, ny+nh 
          ! fill center:
          lats(iy) = southb+0.5*dlat + (iy-1)*dlat
          ! fill bounds:
          lats_bnds(1,iy) = southb + (iy-1)*dlat
          lats_bnds(2,iy) = southb + (iy  )*dlat
        end do    

        ! init grid structure:
        call glb_ugg%Init( lons(1:nx), lons_bnds(:,1:nx), &
                           lats(1:ny), lats_bnds(:,1:ny), status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure:
        call glb_ugg_crnr%Init( lons(1:nx+1)-0.5*dlon, lons_bnds(:,1:nx+1)-0.5*dlon, &
                                lats(1:ny+1)-0.5*dlat, lats_bnds(:,1:ny+1)-0.5*dlat, status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure:
        call glb_ugg_bnds%Init( lons(0:nx+1), lons_bnds(:,0:nx+1), &
                                lats(0:ny+1), lats_bnds(:,0:ny+1), status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure:
        call glb_ugg_crnr_bnds%Init( lons(0:nx+2)-0.5*dlon, lons_bnds(:,0:nx+2)-0.5*dlon, &
                                     lats(0:ny+2)-0.5*dlat, lats_bnds(:,0:ny+2)-0.5*dlat, status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure; lon edges, lat centers:
        call glb_ugg_ustag%Init( lons(1:nx+1)-0.5*dlon, lons_bnds(:,1:nx+1)-0.5*dlon, &
                                 lats(1:ny  )         , lats_bnds(:,1:ny  )         , status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure; lon centers, lat edges:
        call glb_ugg_vstag%Init( lons(1:nx)           , lons_bnds(:,1:nx  )         , &
                                 lats(1:ny+1)-0.5*dlat, lats_bnds(:,1:ny+1)-0.5*dlat, status )
        IF_NOTOK_RETURN(status=1)
        
        ! * local grids
        !   copy from global grids to ensure same results for different decompositions
        
        ! get local shape:
        call dom%Get( status, shp=shp, off=off )
        IF_NOTOK_RETURN(status=1)
        ! copy:
        nx = shp(1)
        ny = shp(2)
        ! info ...
        write (gol,'("local  shape : (",i0,",",i0,")")') nx, ny; call goPr
        
        ! range:
        i1 = off(1) + 1
        i2 = off(1) + shp(1)
        j1 = off(2) + 1
        j2 = off(2) + shp(2)

        ! init grid structure; centers:
        call ugg%Init( lons(i1:i2), lons_bnds(:,i1:i2), &
                       lats(j1:j2), lats_bnds(:,j1:j2), status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure; corners:
        call ugg_crnr%Init( lons(i1:i2+1)-0.5*dlon, lons_bnds(:,i1:i2+1)-0.5*dlon, &
                            lats(j1:j2+1)-0.5*dlat, lats_bnds(:,j1:j2+1)-0.5*dlat, status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure; centers, with halo:
        call ugg_bnds%Init( lons(i1-1:i2+1), lons_bnds(:,i1-1:i2+1), &
                            lats(j1-1:j2+1), lats_bnds(:,j1-1:j2+1), status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure; corners, with halo:
        call ugg_crnr_bnds%Init( lons(i1-1:i2+2)-0.5*dlon, lons_bnds(:,i1-1:i2+2)-0.5*dlon, &
                                 lats(j1-1:j2+2)-0.5*dlat, lats_bnds(:,j1-1:j2+2)-0.5*dlat, status )
        IF_NOTOK_RETURN(status=1)
        
        ! init grid structure; centers, west bound:
        call ugg_west%Init( lons(i1-1:i1-1), lons_bnds(:,i1-1:i1-1), &
                            lats(  j1:j2  ), lats_bnds(:,  j1:j2  ), status )
        IF_NOTOK_RETURN(status=1)
                                    
        ! init grid structure; centers, east bound:
        call ugg_east%Init( lons(i2+1:i2+1), lons_bnds(:,i2+1:i2+1), &
                            lats(  j1:j2  ), lats_bnds(:,  j1:j2  ), status )
        IF_NOTOK_RETURN(status=1)
                                    
        ! init grid structure; centers, south bound:
        call ugg_south%Init( lons(  i1:i2  ), lons_bnds(:,  i1:i2  ), &
                             lats(j1-1:j1-1), lats_bnds(:,j1-1:j1-1), status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure; centers, north bound:
        call ugg_north%Init( lons(  i1:i2  ), lons_bnds(:,  i1:i2  ), &
                             lats(j2+1:j2+1), lats_bnds(:,j2+1:j2+1), status )
        IF_NOTOK_RETURN(status=1)
                             
        ! init grid structure; lon edges, lat centers:
        call ugg_ustag%Init( lons(i1:i2+1)-0.5*dlon, lons_bnds(:,i1:i2+1)-0.5*dlon, &
                             lats(j1:j2  )         , lats_bnds(:,j1:j2  )         , status )
        IF_NOTOK_RETURN(status=1)

        ! init grid structure; lon centers, lat edges:
        call ugg_vstag%Init( lons(i1:i2  )         , lons_bnds(:,i1:i2  )         , &
                             lats(j1:j2+1)-0.5*dlat, lats_bnds(:,j1:j2+1)-0.5*dlat, status )
        IF_NOTOK_RETURN(status=1)
        
        ! * clear:

        ! clear:
        deallocate( lons, stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( lats, stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( lons_bnds, stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( lats_bnds, stat=status )
        IF_NOTOK_RETURN(status=1)
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! non-cartesian grid, cell definition from sample file
      case ( 'non-cartesian' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! sample file:
        call ReadRc( rcF, 'grid.file.name', filename, status )
        IF_NOTOK_RETURN(status=1)
        ! sample variable:
        call ReadRc( rcF, 'grid.file.var', description, status )
        IF_NOTOK_RETURN(status=1)
        ! index space:
        call ReadRc( rcF, 'grid.file.subset', subset, status )
        IF_NOTOK_RETURN(status=1)
    
        ! info ..
        write (gol,'(a,":   read grid definition from   : ",a)') rname, trim(filename); call goPr
        write (gol,'(a,":   sample variable description : ",a)') rname, trim(description); call goPr
        write (gol,'(a,":   index range                 : ",4i4)') rname, subset; call goPr

        ! open file:
        call file_in%Open( trim(filename), status )
        IF_NOTOK_RETURN(status=1)
        ! access variable:
        call file_in%Inq_VarID( trim(description), varid, status )
        IF_NOTOK_RETURN(status=1)
        ! init grid definition of core region:
        call file_in%Get_Grid( varid, ugg, status, subset=subset )
        IF_NOTOK_RETURN(status=1)
        ! init grid definition of core+bounds region:
        call file_in%Get_Grid( varid, ugg_bnds, status, subset=subset+(/-1,1,-1,1/) )
        IF_NOTOK_RETURN(status=1)
        ! close:
        call file_in%Close( status )
        IF_NOTOK_RETURN(status=1)
        
        ! info ...
        write (gol,'(a,":   define corner grids ...")') rname; call goPr
        
        ! init grid on corner points:
        call ugg_crnr%InitFromCorners( ugg, status )
        IF_NOTOK_RETURN(status=1)
        ! init grid on corner points:
        call ugg_crnr_bnds%InitFromCorners( ugg_bnds, status )
        IF_NOTOK_RETURN(status=1)
        
        ! info ...
        write (gol,'(a,":   define staggered grids ...")') rname; call goPr
        
        ! init staggered grids:
        call ugg_ustag%InitFromUEdges( ugg, status )
        IF_NOTOK_RETURN(status=1)
        call ugg_vstag%InitFromVEdges( ugg, status )
        IF_NOTOK_RETURN(status=1)
        
        ! copy shape ...
        nx = ugg%nlon
        ny = ugg%nlat
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! unknown ...
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        write (gol,'("unsupported grid type `",a,"`")') trim(grid_type); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! info ...
    write (gol,'(a,": end")') rname; call goPr
    
    ! * testing ...
    
    !! check ..
    !call ReadRc( rcF, 'test.debug.i', debug_i, status )
    !IF_NOTOK_RETURN(status=1)
    !call ReadRc( rcF, 'test.debug.j', debug_j, status )
    !IF_NOTOK_RETURN(status=1)
    !! get local shape:
    !call dom%Get( status, shp=shp, off=off )
    !IF_NOTOK_RETURN(status=1)
    !! remove offset:
    !debug_i = debug_i - off(1)
    !debug_j = debug_j - off(2)
    !! reset if not local ...
    !if ( any( (/debug_i,debug_j/) > shp ) ) then
    !  debug_i = -999
    !  debug_j = -999
    !end if
    !! info ..
    !write (gol,'("debug local cell : (",i0,",",i0,")")') debug_i, debug_j; call goPr

    ! ok
    status = 0
    
  end subroutine LE_Grid_Init
  
  ! ***
  
  subroutine LE_Grid_Done( status )
    
    ! in/out  ------------------------------------
    
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Grid_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! done with grid on core region:
    call ugg%Done( status )
    IF_NOTOK_RETURN(status=1)
    call ugg_bnds%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! done with grids on corner points:
    call ugg_crnr%Done( status )
    IF_NOTOK_RETURN(status=1)
    call ugg_crnr_bnds%Done( status )
    IF_NOTOK_RETURN(status=1)


    ! Done with bounday grids
    call ugg_west%Done(  status )
    IF_NOTOK_RETURN(status=1)
    call ugg_east%Done( status )
    IF_NOTOK_RETURN(status=1)
    call ugg_south%Done( status )
    IF_NOTOK_RETURN(status=1)
    call ugg_north%Done( status )
    IF_NOTOK_RETURN(status=1)
                             
    ! done with staggered grids:
    call ugg_ustag%Done( status )
    IF_NOTOK_RETURN(status=1)
    call ugg_vstag%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Grid_Done

  

  ! ***  



end module LE_Grid
