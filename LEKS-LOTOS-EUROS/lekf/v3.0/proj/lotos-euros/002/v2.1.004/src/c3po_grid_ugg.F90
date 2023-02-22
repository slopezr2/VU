!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
! DESCRIPTION
!
!    Type to define universal gaussian grid.
!
! EXAMPLE
!
!    ! variables:
!    type(T_Grid_Ugg)      ::  ugg
!    integer               ::  ilon, ilat
!    integer               ::  n       ! number of points in interpolation (1..4)
!    integer               ::  ii(4)   ! indices of points in grid
!    real                  ::  ww(4)   ! interpolation weights
!    real                  ::  value
!    real                  ::  data(:,:) ! input field of size (nlon,nlat)
!    integer               ::  k
!
!    ! Define cartesian grid using lons/lats and cell boundaries.
!    !   real, intent(in)                ::  lons(:)          ! (  nlon)
!    !   real, intent(in)                ::  lons_bnds(:,:)   ! (2,nlon)
!    !   real, intent(in)                ::  lats(:)          ! (  nlat)
!    !   real, intent(in)                ::  lats_bnds(:,:)   ! (2,nlat)
!    !
!    call ugg%Init( lons, lons_bnds, lats, lats_bnds, status )
!    if (status/=0) stop
!
!    ! indices and weights for interpolation to (lon,lat) location:
!    call ugg%GetInterpol( 5.2, 52.4, n, ii, ww, status )
!    if (status/=0) stop
!    ! apply interpolation to data field:
!    value = 0.0
!    do k = 1, n
!      value = value + data(ii(k)) * ww(k)
!    end do
!
!    ! done with grid:
!    call ugg%Done( status )
!    if (status/=0) stop
!
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "c3po.inc"
!
!######################################################################

module C3PO_Grid_Ugg

  use GO    , only : gol, goPr, goErr
  use GO    , only : T_Polygon
  use NetCDF, only : NF90_StrError, NF90_NOERR

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Grid_Ugg
  public  ::  T_Grid_NcDef
  
  public  ::  CRNR_LL, CRNR_UL, CRNR_UR, CRNR_LR
  public  ::  EDGE_LEFT, EDGE_UPPER, EDGE_RIGHT, EDGE_LOWER
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Grid_Ugg'
  
  !            3
  !     UL   upper   UR
  !      4 o<------o 3
  !        |       ^
  ! 4 left |       | right 2
  !        v       | 
  !      1 o------>o 2
  !     LL   lower   LR
  !            1
  
  ! corner indices for (rotated) lon/lat cells, 
  ! should be countr-clock wise for polygon definition:
  integer, parameter    ::  CRNR_LL = 1
  integer, parameter    ::  CRNR_LR = 2
  integer, parameter    ::  CRNR_UR = 3
  integer, parameter    ::  CRNR_UL = 4
  
  ! edge indices for (rotated) lon/lat cells,
  ! use order of corner definitions:
  integer, parameter    ::  EDGE_LOWER = 1
  integer, parameter    ::  EDGE_RIGHT = 2
  integer, parameter    ::  EDGE_UPPER = 3
  integer, parameter    ::  EDGE_LEFT  = 4
  

  ! --- types ----------------------------------------
  
  ! ~ GrADS grid definitions
  
  ! cylindrical (lon,lat)
  type T_Cylindrical
    real                    ::  west
    real                    ::  south
    integer                 ::  nlon
    integer                 ::  nlat
    real                    ::  dlon
    real                    ::  dlat
  contains
    procedure     :: Init               => ProjInit_Cylindrical
    procedure     :: Done               => ProjDone_Cylindrical
    procedure     :: WriteCtlProjection => WriteCtlProjection_Cylindrical
  end type
  
  ! Lambert project (=rotated pole?)
  type T_Lambert
    integer                 ::  nx
    integer                 ::  ny
    real                    ::  lonref
    real                    ::  latref
    integer                 ::  iref
    integer                 ::  jref
    real                    ::  STrueLat
    real                    ::  NTrueLat
    real                    ::  slon
    real                    ::  dx
    real                    ::  dy
  contains
    procedure     :: Init               => ProjInit_Lambert
    procedure     :: Done               => ProjDone_Lambert
    procedure     :: WriteCtlProjection => WriteCtlProjection_Lambert
  end type
  
  ! *
  
  type T_Grid_Ugg
    
    ! type of grid (cartesian, cartesian regular, rotated pole, ...)
    character(len=64)       ::  type
    ! flags:
    logical                 ::  is_cartesian
    
    ! definitions for all gridtypes
    integer                 ::  nlon                     ! number of lons
    integer                 ::  nlat                     ! number of lats
    integer                 ::  nv                       ! number of boundary values (4)

    real, allocatable       ::  longitude(:,:)           ! (nlon,nlat) lons
    real, allocatable       ::  longitude_bnds(:,:,:)    ! (nv,nlon,nlat) lon bounds
    real, allocatable       ::  latitude(:,:)            ! (nlon,nlat) lats
    real, allocatable       ::  latitude_bnds(:,:,:)     ! (nv,nlon,nlat) lat bounds

    real, allocatable       ::  area(:,:)                ! (nlon,nlat) area [m2]
    
    ! corner locations in 2D grids:
    real, allocatable       ::  longitude_corners(:,:)   ! (0:nlon,0:nlat)
    real, allocatable       ::  latitude_corners(:,:)    ! (0:nlon,0:nlat)
    
    ! components of normal vectors:
    real, allocatable       ::  unormal(:,:,:)           ! (2,0:nlon,1:nlat)
    real, allocatable       ::  vnormal(:,:,:)           ! (2,1:nlon,0:nlat)
    
    ! polygons connecting mid points:
    type(T_Polygon), allocatable    ::  mid_pg(:,:)      ! (nlon-1,nlat-1)
    integer, allocatable            ::  mid_ii(:,:,:)    ! (nlon-1,nlat-1,4)
    integer, allocatable            ::  mid_jj(:,:,:)    ! (nlon-1,nlat-1,4)
    
    ! cell polygons
    type(T_Polygon), allocatable    ::  pg(:,:)          ! (nlon,nlat)
    
    ! specials for cartesian gridtypes:
    integer                 ::  nv_1d                    ! number of boundary values in cartesian grid (2)
    real, allocatable       ::  longitude_1d(:)          ! (nlon) lons
    real, allocatable       ::  longitude_bnds_1d(:,:)   ! (nv_1d,nlon) lon bounds
    real, allocatable       ::  latitude_1d(:)           ! (nlon) lats
    real, allocatable       ::  latitude_bnds_1d(:,:)    ! (nv_1d,nlat) lat bounds
    
    real                    ::  dlon                     ! used for cartesian regular
    real                    ::  dlat                     ! used for cartesian regular 
    
    ! specials for rotated pole grid types:
    !real                    ::  rpol_lon
    !real                    ::  rpol_lat
    ! ...
    
    ! projection types (cylindrical, Lambert, ... )
    type(T_Cylindrical), pointer  ::  Cylindrical
    type(T_Lambert), pointer      ::  Lambert

  contains
    procedure   ::                     Grid_Ugg_Init_Cartesian 
    procedure   ::                     Grid_Ugg_Init_RotatedPole 
    generic     ::  Init            => Grid_Ugg_Init_Cartesian, &
                                       Grid_Ugg_Init_RotatedPole
    procedure   ::  InitFromCorners => Grid_Ugg_InitFromCorners
    procedure   ::  InitFromUEdges  => Grid_Ugg_InitFromUEdges
    procedure   ::  InitFromVEdges  => Grid_Ugg_InitFromVEdges
    procedure   ::  InitFromStags   => Grid_Ugg_InitFromStags
    procedure   ::  Done            => Grid_Ugg_Done
    
    procedure   ::  GetIndexDomain  => Grid_Ugg_GetIndexDomain
    procedure   ::  GetBoundingBox  => Grid_Ugg_GetBoundingBox
    procedure   ::  GetEdgeLength   => Grid_Ugg_GetEdgeLength
    
    procedure   ::  GetLonIndex        => Grid_Ugg_GetLonIndex
    procedure   ::  GetLatIndex        => Grid_Ugg_GetLatIndex
    procedure   ::  GetLocation        => Grid_Ugg_GetLocation
    procedure   ::  GetInterpol        => Grid_Ugg_GetInterpol
    procedure   ::  GetDistribution    => Grid_Ugg_GetDistribution_Cartesian
    
    procedure   ::  GetCellDistance    => Grid_Ugg_GetCellDistance
    procedure   ::  WithinDistance     => Grid_Ugg_WithinDistance
    procedure   ::  DistanceGrid       => Grid_Ugg_DistanceGrid
    procedure   ::  RoundToResolution  => Grid_Ugg_RoundToResolution
                                              
    procedure   ::                        Grid_Ugg_InDomain_Point
    procedure   ::                        Grid_Ugg_InDomain_Cell
    generic     ::  InDomain           => Grid_Ugg_InDomain_Point, &
                                          Grid_Ugg_InDomain_Cell

    procedure   ::                        Grid_Ugg_AreaOper_2d
    procedure   ::                        Grid_Ugg_AreaOper_3d
    generic     ::  AreaOper           => Grid_Ugg_AreaOper_2d, &
                                          Grid_Ugg_AreaOper_3d
                                    
    procedure   ::  DefGrid_NetCDF     => Grid_Ugg_DefGrid_NetCDF
    procedure   ::  DefCoor_NetCDF     => Grid_Ugg_DefCoor_NetCDF
    procedure   ::  PutGrid_NetCDF     => Grid_Ugg_PutGrid_NetCDF
    
    procedure   ::  WriteCtlProjection => Grid_Ugg_WriteCtlProjection
                                        
  end type T_Grid_Ugg
  
  ! *
  
  ! stuff to write grid definition to netCDF file:
  type T_Grid_NcDef
    ! cell range:
    integer                ::  subset(4)   ! (/i1,i2,j1,j2/)
    ! output shape:
    integer                ::  nx
    integer                ::  ny
    ! file id:
    integer                ::  ncid
    ! dimensions:
    integer                ::  dimid_x
    integer                ::  dimid_y
    integer                ::  dimid_xc
    integer                ::  dimid_yc
    ! variables:
    integer                ::  varid_x
    integer                ::  varid_y
    integer                ::  varid_lon
    integer                ::  varid_lat
    integer                ::  varid_lon_bnds
    integer                ::  varid_lat_bnds
    integer                ::  varid_lon_crnr
    integer                ::  varid_lat_crnr
    ! attribue:
    character(len=128)     ::  coordinates
  end type T_Grid_NcDef

  

contains


  ! ********************************************************************
  ! ***
  ! *** ugg
  ! ***
  ! ********************************************************************


  subroutine Grid_Ugg_Init_Cartesian( self, lons, lons_bnds, lats, lats_bnds, status )
  
    use GO             , only : T_Vector
    use C3PO_Grid_Tools, only : ll_area_deg_m2, ll_distance
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(out)  ::  self
    real, intent(in)                ::  lats(:)          ! (   nlat)
    real, intent(in)                ::  lats_bnds(:,:)   ! (nv,nlat)
    real, intent(in)                ::  lons(:)          ! (   nlon)
    real, intent(in)                ::  lons_bnds(:,:)   ! (nv,nlon)
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_Init_Cartesian'
    
    ! --- local ----------------------------------
    
    integer                         ::  i, j
    integer                         ::  k
    real                            ::  d
    
    ! --- begin ----------------------------------
    
    ! size:
    self%nlon  = size(lons)
    self%nlat  = size(lats)
    self%nv    = 4
    self%nv_1d = 2
    
    ! check ...
    if ( any( shape(lons_bnds) /=  (/self%nv_1d,self%nlon/) ) ) then
      write (gol,'("shape of argument lons_bnds is ",2i6," not ",2i6)') shape(lons_bnds), self%nv,self%nlon; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( any( shape(lats_bnds) /=  (/self%nv_1d,self%nlat/) ) ) then
      write (gol,'("shape of argument lats_bnds is ",2i6," not ",2i6)') shape(lats_bnds), self%nv,self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! default type, assume regular:
    self%type  = 'cartesian-regular'
    ! init step sizes:
    self%dlon = -999.9
    ! loop over lon cells:
    do i = 1, self%nlon
      ! current spacing:
      d = abs( lons_bnds(2,i) - lons_bnds(1,i) )
      ! first?
      if ( self%dlon <= 0.0 ) then
        self%dlon = d
      else
        ! compare:
        if ( abs( d - self%dlon ) > 1.0e-4 ) then
          ! irregular spacing ..
          self%type = 'cartesian'
          ! leave:
          exit
        end if ! compare
      end if ! first
    end do ! i
    ! init step sizes:
    self%dlat = -999.9
    ! loop over lat cells:
    do j = 1, self%nlat
      ! current spacing:
      d = abs( lats_bnds(2,j) - lats_bnds(1,j) )
      ! first?
      if ( self%dlat <= 0.0 ) then
        self%dlat = d
      else
        ! compare:
        if ( abs( d - self%dlat ) > 1.0e-4 ) then
          ! irregular spacing ..
          self%type = 'cartesian'
          ! leave:
          exit
        end if ! compare
      end if ! first
    end do ! j
    
    ! set flag:
    self%is_cartesian = .true.
            
    ! storage for cell centers:
    allocate( self%longitude(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for cell bounds:
    allocate( self%longitude_bnds(self%nv,self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude_bnds(self%nv,self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for area:
    allocate( self%area(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! no GrADS grid definitions yet:
    nullify( self%Cylindrical )
    nullify( self%Lambert )
    
    ! ~ bounds

    ! loop over cells:
    do j = 1, self%nlat
      do i = 1, self%nlon

        ! fill center:
        self%longitude(i,j) = lons(i)
        self%latitude(i,j)  = lats(j)

        ! fill corners, support reversed orders:
        self%latitude_bnds (CRNR_LL,i,j) = minval( lats_bnds(:,j) )   ! lower
        self%longitude_bnds(CRNR_LL,i,j) = minval( lons_bnds(:,i) )   ! left
        self%latitude_bnds (CRNR_LR,i,j) = minval( lats_bnds(:,j) )   ! lower
        self%longitude_bnds(CRNR_LR,i,j) = maxval( lons_bnds(:,i) )   ! right
        self%latitude_bnds (CRNR_UR,i,j) = maxval( lats_bnds(:,j) )   ! upper
        self%longitude_bnds(CRNR_UR,i,j) = maxval( lons_bnds(:,i) )   ! right
        self%latitude_bnds (CRNR_UL,i,j) = maxval( lats_bnds(:,j) )   ! upper
        self%longitude_bnds(CRNR_UL,i,j) = minval( lons_bnds(:,i) )   ! left
        
        !if (i==368 .and. j==221) then
        !  print *, 'uuu1 longitude_bnds = ', self%longitude_bnds(:,i,j)
        !  print *, '  u1 latitude_bnds  = ', self%latitude_bnds(:,i,j)
        !  if ( self%latitude_bnds(3,i,j) < self%latitude_bnds(2,i,j) ) then
        !    print *, '  u1 wrong lat order'
        !    TRACEBACK; status=1; return
        !  end if
        !end if

      end do ! i
    end do ! j
    
    ! ~ cell corners

    ! storage for cell corners:
    allocate( self%longitude_corners(0:self%nlon,0:self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude_corners(0:self%nlon,0:self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! lon bounds:
    self%longitude_corners(0,:) = lons_bnds(1,1)  ! left
    do i = 1, self%nlon
      self%longitude_corners(i,:) = lons_bnds(2,i)  ! right
    end do
    
    ! lat bounds:
    self%latitude_corners(:,0) = lats_bnds(1,1)   ! lower
    do j = 1, self%nlat
      self%latitude_corners(:,j) = lats_bnds(2,j)  ! upper
    end do
    
    ! ~ 1D axes
    
    ! storage for 1-dimensional arrays:
    allocate( self%longitude_1d(self%nlon), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%longitude_bnds_1d(self%nv_1d,self%nlon), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude_1d(self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude_bnds_1d(self%nv_1d,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    self%longitude_1d      = lons
    self%longitude_bnds_1d = lons_bnds
    self%latitude_1d       = lats
    self%latitude_bnds_1d  = lats_bnds
  
    ! loop over cells:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! area in m2:
        self%area(i,j) = ll_area_deg_m2( self%longitude_bnds_1d(1,i), self%longitude_bnds_1d(2,i), &
                                         self%latitude_bnds_1d(1,j),  self%latitude_bnds_1d(2,j) )
      end do
    end do
    
    ! ~ cell polygons

    ! storage for polygons:
    allocate( self%pg(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over cells:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! define polygon by corners, order is counter clockwise:
        call self%pg(i,j)%Init( real(self%longitude_bnds(:,i,j),8), &
                                real(self%latitude_bnds (:,i,j),8), status )
        IF_NOT_OK_RETURN(status=1)
      end do ! i
    end do ! j
    
    ! ~ normals
    
    ! fill components of normal vectors:
    call Grid_Ugg_FillNormals( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ GrADS projections

    ! Storage for projection:
    allocate( self%Cylindrical, stat=status)
    IF_NOT_OK_RETURN(status=1)
    ! Fill cartesian projection properties:
    call self%Cylindrical%Init( self%longitude_bnds_1d(1,1), self%latitude_bnds_1d(1,1), &
                                 self%nlon, self%nlat, self%dlon, self%dlat, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine Grid_Ugg_Init_Cartesian


  ! ***


  subroutine Grid_Ugg_Init_RotatedPole( self, xlons, xlats, status )
  
    use GO             , only : T_Vector
    use Num            , only : CenterToCorners
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(out)  ::  self
    real, intent(in)                ::  xlons(:,:)
    real, intent(in)                ::  xlats(:,:)
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_Init_RotatedPole'
    
    ! --- local ----------------------------------
    
    integer                         ::  i, j
    integer                         ::  k
    integer                         ::  mi, mj
    real(8)                         ::  xx(4), yy(4)
    real(8)                         ::  area
    
    ! --- begin ----------------------------------
    
    ! set type:
    self%type = 'rotated-pole'
    ! set flag:
    self%is_cartesian = .false.
    
    ! size:
    self%nlon  = size(xlons,1)
    self%nlat  = size(xlons,2)
    self%nv    = 4
    
    ! no 1D yet ...
    self%nv_1d = -999
    
    ! no regular spacing:
    self%dlon = -999.9
    self%dlat = -999.9
    
    ! ~ centers
            
    ! storage for cell centers:
    allocate( self%longitude(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    self%longitude = xlons
    self%latitude  = xlats
    
    ! ~ polygons connecting mid points:
    
    ! storage for polygons between corners:
    allocate( self%mid_pg(self%nlon-1,self%nlat-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%mid_ii(self%nlon-1,self%nlat-1,self%nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%mid_jj(self%nlon-1,self%nlat-1,self%nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over polygons between mid points:
    do j = 1, self%nlat-1
      do i = 1, self%nlon-1
        ! define center indices:
        !     4 o--o 3  j+1
        !     1 o--o 2  j
        !       i i+1
        self%mid_ii(i,j,:) = (/ i  , i+1, i+1, i   /)
        if ( self%latitude(i,j) < self%latitude(i,j+1) ) then
          self%mid_jj(i,j,:) = (/ j  , j  , j+1, j+1 /)
        else
          self%mid_jj(i,j,:) = (/ j+1, j+1, j  , j   /)
        end if
        ! collect center locations:
        do k = 1, self%nv ! 4
          ! center index:
          mi = self%mid_ii(i,j,k)
          mj = self%mid_jj(i,j,k)
          ! define point:
          xx(k) = real(self%longitude(mi,mj),8)
          yy(k) = real(self%latitude (mi,mj),8)
        end do ! cornters
        ! define polygon:
        call self%mid_pg(i,j)%Init( xx, yy, status )
        IF_NOT_OK_RETURN(status=1)
      end do ! i
    end do ! j
    
    ! ~ cell corners

    ! storage for cell corners:
    allocate( self%longitude_corners(0:self%nlon,0:self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude_corners(0:self%nlon,0:self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! fill corners, just interpolate/extrapolate from centers:
    call CenterToCorners( self%longitude, self%longitude_corners, status )
    IF_NOT_OK_RETURN(status=1)
    call CenterToCorners( self%latitude, self%latitude_corners, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ cell bounds

    ! storage for cell bounds:
    allocate( self%longitude_bnds(self%nv,self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude_bnds(self%nv,self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy from corners:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! fill, take care of reversed orders:
        self%longitude_bnds(CRNR_LL,i,j) = minval( self%longitude_corners(i-1:i,j-1) )
        self%longitude_bnds(CRNR_LR,i,j) = maxval( self%longitude_corners(i-1:i,j-1) )
        self%longitude_bnds(CRNR_UR,i,j) = maxval( self%longitude_corners(i-1:i,j  ) )
        self%longitude_bnds(CRNR_UL,i,j) = minval( self%longitude_corners(i-1:i,j  ) )
        ! fill, take care of reversed orders:
        self%latitude_bnds (CRNR_LL,i,j) = minval( self%latitude_corners (i-1,j-1:j) )
        self%latitude_bnds (CRNR_LR,i,j) = minval( self%latitude_corners (  i,j-1:j) )
        self%latitude_bnds (CRNR_UR,i,j) = maxval( self%latitude_corners (  i,j-1:j) )
        self%latitude_bnds (CRNR_UL,i,j) = maxval( self%latitude_corners (i-1,j-1:j) )
        
        !if (i==368 .and. j==221) then
        !  print *, 'uuu2 longitude_bnds = ', self%longitude_bnds(:,i,j)
        !  print *, '  u2 latitude_bnds  = ', self%latitude_bnds(:,i,j)
        !  if ( self%latitude_bnds(3,i,j) < self%latitude_bnds(2,i,j) ) then
        !    print *, '  u2 wrong lat order'
        !    TRACEBACK; status=1; return
        !  end if
        !end if

      end do ! i
    end do ! j
    
    ! ~ cell polygons
    
    ! storage for polygons:
    allocate( self%pg(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for area:
    allocate( self%area(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over cells:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! define polygon:
        call self%pg(i,j)%Init( real(self%longitude_bnds(:,i,j),8), &
                                real(self%latitude_bnds (:,i,j),8), status )
        IF_NOT_OK_RETURN(status=1)
        ! fill area:
        call self%pg(i,j)%LonLat_Area( area, status )
        IF_NOT_OK_RETURN(status=1)
        ! copy:
        self%area(i,j) = real(area)
      end do ! i
    end do ! j

    ! ~ normals
    
    ! fill components of normal vectors:
    call Grid_Ugg_FillNormals( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ GrADS projection info
    
    ! Storage for projection:
    allocate( self%Lambert, stat=status)
    IF_NOT_OK_RETURN(status=1)

    ! Fill projection properties        
    call self%Lambert%Init( self%nlon, self%nlat, &
                            self%longitude(1,1), self%latitude(1,1), &
                            1, 1, self%latitude(1,1), self%latitude(1,self%nlat), self%longitude(1,1), &
                            1.0, 1.0, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine Grid_Ugg_Init_RotatedPole


  ! ***

  
  subroutine Grid_Ugg_FillNormals( self, status )
  
    use Num            , only : CenterToCorners
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(inout)  ::  self
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_FillNormals'
    
    ! --- local ----------------------------------
    
    integer                    ::  i, j
    real(8)                    ::  xx(4), yy(4)
    
    ! --- begin ----------------------------------
    
    ! storage:
    allocate( self%unormal(2,0:self%nlon,1:self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%vnormal(2,1:self%nlon,0:self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over cells
    do j = 1, self%nlat
      do i = 1, self%nlon
      
        ! get normal vectors:
        call self%pg(i,j)%GetNormals( xx, yy, status )
        IF_NOT_OK_RETURN(status=1)
        
        !print *, 'nnn1 xx = ', xx
        !print *, '  n1 yy = ', yy
        !stop 'break'
        
        ! fill left-most u-edge from first polygon:
        if ( i == 1 ) then
          ! swap from leftward to rightward:
          self%unormal(1,0,j) = -1.0 * xx(EDGE_LEFT)
          self%unormal(2,0,j) = -1.0 * yy(EDGE_LEFT)
        end if
        ! u-edges at right side of cell:
        self%unormal(1,i,j) = xx(EDGE_RIGHT)
        self%unormal(2,i,j) = yy(EDGE_RIGHT)
        
        ! fill lower-most v-edge from first polygon:
        if ( j == 1 ) then
          ! swap from downward to upward:
          self%vnormal(1,i,0) = -1.0 * xx(EDGE_LOWER)
          self%vnormal(2,i,0) = -1.0 * yy(EDGE_LOWER)
        end if
        ! u-edges at right side of cell:
        self%vnormal(1,i,j) = xx(EDGE_UPPER)
        self%vnormal(2,i,j) = yy(EDGE_UPPER)
        
      end do ! i
    end do ! j

    ! ok
    status = 0
  
  end subroutine Grid_Ugg_FillNormals


  ! ***

  
  ! 
  ! Init grid from the corners of other grid.
  !

  subroutine Grid_Ugg_InitFromCorners( self, ugg, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(out)  ::  self
    class(T_Grid_Ugg), intent(in)   ::  ugg
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_InitFromCorners'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init from corners:
    call self%Init( ugg%longitude_corners, ugg%latitude_corners, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Grid_Ugg_InitFromCorners


  ! ***

  
  ! 
  ! Init grid on mid of the u-edges of other grid.
  !

  subroutine Grid_Ugg_InitFromUEdges( self, ugg, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(out)  ::  self
    class(T_Grid_Ugg), intent(in)   ::  ugg
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_InitFromUEdges'
    
    ! --- local ----------------------------------
    
    real, allocatable     ::  lons(:,:)
    real, allocatable     ::  lats(:,:)
    integer               ::  i, j
    
    ! --- begin ----------------------------------
    
    ! storage for u-edges:
    allocate( lons(0:ugg%nlon,1:ugg%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( lats(0:ugg%nlon,1:ugg%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over edges:
    do j = 1, ugg%nlat
      do i = 0, ugg%nlon
        ! fill average of lower and upper corner:
        lons(i,j) = 0.5*( ugg%longitude_corners(i,j-1) + ugg%longitude_corners(i,j) )
        lats(i,j) = 0.5*( ugg%latitude_corners (i,j-1) + ugg%latitude_corners (i,j) )
      end do  ! i
    end do ! j
    
    ! init from arrays:
    call self%Init( lons, lats, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( lons, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( lats, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Grid_Ugg_InitFromUEdges


  ! ***

  
  ! 
  ! Init grid on mid of the v-edges of other grid.
  !

  subroutine Grid_Ugg_InitFromVEdges( self, ugg, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(out)  ::  self
    class(T_Grid_Ugg), intent(in)   ::  ugg
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_InitFromVEdges'
    
    ! --- local ----------------------------------
    
    real, allocatable     ::  lons(:,:)
    real, allocatable     ::  lats(:,:)
    integer               ::  i, j
    
    ! --- begin ----------------------------------
    
    ! storage for v-edges:
    allocate( lons(1:ugg%nlon,0:ugg%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( lats(1:ugg%nlon,0:ugg%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over edges:
    do j = 0, ugg%nlat
      do i = 1, ugg%nlon
        ! fill average of left and right corner:
        lons(i,j) = 0.5*( ugg%longitude_corners(i-1,j) + ugg%longitude_corners(i,j) )
        lats(i,j) = 0.5*( ugg%latitude_corners (i-1,j) + ugg%latitude_corners (i,j) )
      end do  ! i
    end do ! j
    
    ! init from arrays:
    call self%Init( lons, lats, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( lons, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( lats, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Grid_Ugg_InitFromVEdges


  ! ***

  
  ! 
  ! Init grid on centers between staggered grids.
  ! If padding=.true., the staggered grids are missing the first edge
  ! and the first grid points are extrapolated.
  !

  subroutine Grid_Ugg_InitFromStags( self, ugg_u, ugg_v, status, &
                                       padding )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(out)  ::  self
    class(T_Grid_Ugg), intent(in)   ::  ugg_u
    class(T_Grid_Ugg), intent(in)   ::  ugg_v
    integer, intent(out)            ::  status

    logical, intent(in), optional   ::  padding

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_InitFromStags'
    
    ! --- local ----------------------------------
    
    logical               ::  without_first
    integer               ::  nlon, nlat
    real, allocatable     ::  lons(:,:)
    real, allocatable     ::  lats(:,:)
    integer               ::  i, j
    
    ! --- begin ----------------------------------
    
    ! set flag:
    without_first = .false.
    if ( present(padding) ) without_first = padding
    
    ! target shape:
    nlon = ugg_v%nlon
    nlat = ugg_u%nlat

    ! storage for centers:
    allocate( lons(nlon,nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( lats(nlon,nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! switch ...
    if ( without_first ) then
    
      ! check ...
      if ( any( (/ugg_u%nlon,ugg_u%nlat/) /= (/ugg_v%nlon,ugg_v%nlat/) ) ) then
        write (gol,'("shapes of input grids do not seem to be staggered:")'); call goErr
        write (gol,'("  u-grid shape : ",2i6)') ugg_u%nlon, ugg_u%nlat; call goErr
        write (gol,'("  v-grid shape : ",2i6)') ugg_v%nlon, ugg_v%nlat; call goErr
        write (gol,'("  target shape : ",2i6)') nlon, nlat; call goErr
        TRACEBACK; status=1; return
      end if

      ! fill with average over all 4 edges if present:
      do j = 1, nlat
        do i = 1, nlon
          ! switch:
          if ( (i == 1) .and. (j == 1) ) then
            lons(i,j) =          ugg_v%longitude(i  ,j  )
            lats(i,j) =          ugg_u%latitude (i  ,j  )
          else if ( i == 1 ) then
            lons(i,j) = 0.50 * ( ugg_v%longitude(i  ,j-1) + &
                                 ugg_v%longitude(i  ,j  ) )
            lats(i,j) = 0.50 * ( ugg_v%latitude (i  ,j-1) + &
                                 ugg_v%latitude (i  ,j  ) )
          else if ( j == 1 ) then
            lons(i,j) = 0.50 * ( ugg_u%longitude(i-1,j  ) + &
                                 ugg_u%longitude(i  ,j  ) )
            lats(i,j) = 0.50 * ( ugg_u%latitude (i-1,j  ) + &
                                 ugg_u%latitude (i  ,j  ) )
          else
            lons(i,j) = 0.25 * ( ugg_u%longitude(i-1,j  ) + &
                                 ugg_u%longitude(i  ,j  ) + &
                                 ugg_v%longitude(i  ,j-1) + &
                                 ugg_v%longitude(i  ,j  ) )
            lats(i,j) = 0.25 * ( ugg_u%latitude (i-1,j  ) + &
                                 ugg_u%latitude (i  ,j  ) + &
                                 ugg_v%latitude (i  ,j-1) + &
                                 ugg_v%latitude (i  ,j  ) )
          end if
        end do  ! i
      end do  ! j

    else
    
      ! check ...
      if ( any( (/ugg_u%nlon,ugg_u%nlat/) /= (/ugg_v%nlon+1,ugg_v%nlat-1/) ) ) then
        write (gol,'("shapes of input grids do not seem to be staggered:")'); call goErr
        write (gol,'("  u-grid shape : ",2i6)') ugg_u%nlon, ugg_u%nlat; call goErr
        write (gol,'("  v-grid shape : ",2i6)') ugg_v%nlon, ugg_v%nlat; call goErr
        write (gol,'("  target shape : ",2i6)') nlon, nlat; call goErr
        TRACEBACK; status=1; return
      end if

      ! fill with average over all 4 edges:
      do j = 1, nlat
        do i = 1, nlon
          lons(i,j) = 0.25 * ( ugg_u%longitude(i  ,j  ) + &
                               ugg_u%longitude(i+1,j  ) + &
                               ugg_v%longitude(i  ,j  ) + &
                               ugg_v%longitude(i  ,j+1) )
          lats(i,j) = 0.25 * ( ugg_u%latitude (i  ,j  ) + &
                               ugg_u%latitude (i+1,j  ) + &
                               ugg_v%latitude (i  ,j  ) + &
                               ugg_v%latitude (i  ,j+1) )
        end do  ! i
      end do  ! j
                      
    end if
    
    ! init from arrays:
    call self%Init( lons, lats, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( lons, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( lats, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Grid_Ugg_InitFromStags


  ! ***
  
  
  subroutine Grid_Ugg_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(inout)    ::  self
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_Done'
    
    ! --- local ----------------------------------
    
    integer   ::  i, j
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%longitude, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%latitude, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! cell bounds:
    deallocate( self%longitude_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%latitude_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! cell properties:
    deallocate( self%area, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! corners:
    deallocate( self%longitude_corners, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%latitude_corners, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! normals:
    deallocate( self%unormal, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%vnormal, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! cartesian Grid?
    if ( allocated(self%longitude_1d) ) then
      deallocate( self%longitude_1d, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%longitude_bnds_1d, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%latitude_1d, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%latitude_bnds_1d, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! polygons?
    if ( allocated(self%mid_pg) ) then
      ! loop:
      do j = 1, self%nlat-1
        do i = 1, self%nlon-1
          ! done:
          call self%mid_pg(i,j)%Done( status )
          IF_NOT_OK_RETURN(status=1)
        end do
      end do
      ! clear:
      deallocate( self%mid_pg, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%mid_ii, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%mid_jj, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! cell polygons?
    if ( allocated(self%pg) ) then
      ! loop:
      do j = 1, self%nlat
        do i = 1, self%nlon
          ! done:
          call self%pg(i,j)%Done( status )
          IF_NOT_OK_RETURN(status=1)
        end do
      end do
      ! clear:
      deallocate( self%pg, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! Clear Cylindrical projection properties
    if ( associated( self%Cylindrical ) ) then
      call self%Cylindrical%Done( status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%Cylindrical, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! Clear Lambert projection properties
    if ( associated( self%Lambert ) ) then
      call self%Lambert%Done( status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%Lambert, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
  
  end subroutine Grid_Ugg_Done
  

  ! ***

  
  !
  ! Return index range  (/i1,i2,j1,j2/) for which the cell centers
  ! enclose the target box (/x1,x2,y1,y2/).
  !
  

  subroutine Grid_Ugg_GetIndexDomain( self, bbox, ibox, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(in)                ::  bbox(4)  ! (/west,east,south,north/)
    integer, intent(out)            ::  ibox(4)
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_GetIndexDomain'
    
    ! --- local ----------------------------------
    
    integer             ::  i, j
    integer             ::  idir, jdir
    
    ! --- begin ----------------------------------
    
    ! init box to maximum:
    ibox = (/1,self%nlon,1,self%nlat/)
    
    ! quick scan to select a small i-range
    ! that contains the target point ;
    !~ direction:
    if ( self%longitude(1,1) < self%longitude(self%nlon,1) ) then
      idir = 1
    else
      idir = -1
    end if
    
    ! switch:
    if ( self%is_cartesian ) then

      ! search for i with west bound:
      do i = 1, self%nlon
        ! in box?
        if ( (idir*self%longitude_bnds_1d(1,i) <= idir*bbox(1)) .and. &
             (idir*self%longitude_bnds_1d(2,i) >= idir*bbox(1)) ) then
          ! previous:
          ibox(1) = i
          ! leave:
          exit
        end if
      end do    
      ! search for i with east bound:
      do i = 1, self%nlon
        ! in box?
        if ( (idir*self%longitude_bnds_1d(1,i) <= idir*bbox(2)) .and. &
             (idir*self%longitude_bnds_1d(2,i) >= idir*bbox(2)) ) then
          ! next:
          ibox(2) = i
          ! leave:
          exit
        end if
      end do
      
    else  ! non-cartesian
      
      !~ loop over columns:
      do i = 1, self%nlon
        ! all left of target ?
        if ( all( idir*self%longitude(i,:) <= idir*bbox(1) ) ) then
          ! reset:
          ibox(1) = max( i, ibox(1) )
        end if
        ! all right of traget?
        if ( all( idir*self%longitude(i,:) >= idir*bbox(2) ) ) then
          ! reset:
          ibox(2) = min( i, ibox(2) )
          ! found end, leave:
          exit
        end if
      end do ! i
      
    end if  ! cartesian?

    ! quick scan to select a small j-range
    ! that contains the target point ;
    !~ direction:
    if ( self%latitude(1,1) < self%latitude(1,self%nlat) ) then
      jdir = 1
    else
      jdir = -1
    end if
    ! switch:
    if ( self%is_cartesian ) then

      ! search for j with south bound:
      do j = 1, self%nlat
        ! in box?
        if ( (jdir*self%latitude_bnds_1d(1,j) <= jdir*bbox(3)) .and. &
             (jdir*self%latitude_bnds_1d(2,j) >= jdir*bbox(3)) ) then
          ! previous:
          ibox(3) = j
          ! leave:
          exit
        end if
      end do    
      ! search for j with north bound:
      do j = 1, self%nlat
        ! in box?
        if ( (jdir*self%latitude_bnds_1d(1,j) <= jdir*bbox(4)) .and. &
             (jdir*self%latitude_bnds_1d(2,j) >= jdir*bbox(4)) ) then
          ! next:
          ibox(4) = j
          ! leave:
          exit
        end if
      end do
      
    else  ! non-cartesian

      !~ loop over columns:
      do j = 1, self%nlat
        ! all below target ?
        if ( all( jdir*self%latitude(ibox(1):ibox(2),j) <= jdir*bbox(3) ) ) then
          ! reset:
          ibox(3) = max( j, ibox(3) )
        end if
        ! all above traget?
        if ( all( jdir*self%latitude(ibox(1):ibox(2),j) >= jdir*bbox(4) ) ) then
          ! reset:
          ibox(4) = min( j, ibox(4) )
          ! found end, leave:
          exit
        end if
      end do ! j
      
    end if  ! cartesian?
    
    ! check ...
    if ( (ibox(2) - ibox(1) + 1 >= 100) .or. (ibox(4) - ibox(3) + 1 >= 100) ) then
      write (gol,'("found very large index range; something wrong?")'); call goErr
      write (gol,'("  target box   : ",4f11.5)') bbox; call goErr
      write (gol,'("  index range  : ",4i8)') ibox; call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine Grid_Ugg_GetIndexDomain
  

  ! ***

  
  !
  ! Return bounding box enclosing domain.
  !
  

  subroutine Grid_Ugg_GetBoundingBox( self, bbox, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(out)               ::  bbox(4)  ! (/west,east,south,north/)
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_GetBoundingBox'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! bounds enclosing corners:
    bbox(1) = minval( self%longitude_corners )
    bbox(2) = maxval( self%longitude_corners )
    bbox(3) = minval( self%latitude_corners  )
    bbox(4) = maxval( self%latitude_corners  )
   
    ! ok
    status = 0

  end subroutine Grid_Ugg_GetBoundingBox
  

  ! ***

  
  !
  ! Return edge length of specified cell
  !
  

  subroutine Grid_Ugg_GetEdgeLength( self, i, j, edge, length, status )
  
    use C3PO_Grid_Tools, only : ll_distance

    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    integer, intent(in)             ::  i, j
    integer, intent(in)             ::  edge
    real, intent(out)               ::  length
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_GetEdgeLength'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (i < 1) .or. (i > self%nlon) .or. &
         (j < 1) .or. (j > self%nlat)      ) then
      write (gol,'("cell (",i0,",",i0,") outside grid shape (",i0,",",i0,")")') i,j, self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! switch:
    select case ( edge )
      case ( EDGE_LOWER )
        length = ll_distance( self%longitude_bnds(CRNR_LR,i,j),self%latitude_bnds(CRNR_LR,i,j), &
                              self%longitude_bnds(CRNR_LL,i,j),self%latitude_bnds(CRNR_LL,i,j) )
      case ( EDGE_RIGHT )
        length = ll_distance( self%longitude_bnds(CRNR_UR,i,j),self%latitude_bnds(CRNR_UR,i,j), &
                              self%longitude_bnds(CRNR_LR,i,j),self%latitude_bnds(CRNR_LR,i,j) )
      case ( EDGE_UPPER )
        length = ll_distance( self%longitude_bnds(CRNR_UL,i,j),self%latitude_bnds(CRNR_UL,i,j), &
                              self%longitude_bnds(CRNR_UR,i,j),self%latitude_bnds(CRNR_UR,i,j) )
      case ( EDGE_LEFT )
        length = ll_distance( self%longitude_bnds(CRNR_LL,i,j),self%latitude_bnds(CRNR_LL,i,j), &
                              self%longitude_bnds(CRNR_UL,i,j),self%latitude_bnds(CRNR_UL,i,j) )
      case default
        write (gol,'("unsupported edge id: ",i0)') edge; call goErr
        TRACEBACK; status=1; return
    end select
   
    ! ok
    status = 0

  end subroutine Grid_Ugg_GetEdgeLength
  
  
  ! ***

  
  ! index of band including longitude

  subroutine Grid_Ugg_GetLonIndex( self, lon, ilon, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(in)                ::  lon  ! [deg]
    integer, intent(out)            ::  ilon
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_GetLonIndex'
    
    ! --- local ----------------------------------
    
    integer                 ::  i
    
    ! --- begin ----------------------------------
    
    ! init result:
    ilon = -999
    
    ! check grid type:
    select case ( trim(self%type) )
      case ( 'cartesian', 'cartesian-regular' )
        ! loop over  bands:
        do i = 1, self%nlon
          ! within bounds ?
          if ( ((lon       >= minval(self%longitude_bnds_1d(:,i))) .and. (lon       <= maxval(self%longitude_bnds_1d(:,i)))) .or. &
               ((lon-360.0 >= minval(self%longitude_bnds_1d(:,i))) .and. (lon-360.0 <= maxval(self%longitude_bnds_1d(:,i))))      ) then
            ! fill result:
            ilon = i
            ! leave:
            exit
          end if
        end do
        ! check ..
        if ( ilon < 0 ) then
          write (gol,'("longitude ",f12.6," not covered by grid")') lon; call goErr
          write (gol,'("gridtype: ",a )' ) trim(self%type) ; call goErr
          write (gol,'("grid longitude range: [",f7.2,",",f7.2,"]")' ) minval(self%longitude_bnds_1d), maxval(self%longitude_bnds_1d); call goErr
          TRACEBACK; status=1; return
        end if
      
      case default
        write (gol,'("Function GetLonIndex not defined for grid type: ",a)') trim(self%type) ; call goErr
        TRACEBACK; status=1; return
        
        ! NB. For 2d defined lon/lat grids, ilon ant ilat must be found at the same time
        ! Functions getLocation needed
          
    end select

    ! ok
    status = 0

  end subroutine Grid_Ugg_GetLonIndex


  ! *
  
  ! index of band including latitude

  subroutine Grid_Ugg_GetLatIndex( self, lat, ilat, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(in)                ::  lat  ! [deg]
    integer, intent(out)            ::  ilat
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_GetLatIndex'
    
    ! --- local ----------------------------------
    
    integer                 ::  j
    
    ! --- begin ----------------------------------
    
    ! init result:
    ilat = -999
    ! check grid type:
    select case ( trim(self%type) )
      case ( 'cartesian', 'cartesian-regular' )
        ! loop over  bands:
        do j = 1, self%nlat
          ! within bouonds ?
          if ( (lat >= minval(self%latitude_bnds_1d(:,j))) .and. (lat <= maxval(self%latitude_bnds_1d(:,j))) ) then
            ! fill result:
            ilat = j
            ! leave:
            exit
          end if
        end do
        ! check ..
        if ( ilat < 0 ) then
          write (gol,'("latitude ",f12.6," not covered by grid")') lat; call goErr
          write (gol,'("Gridtype: ",a )' ) trim(self%type) ; call goErr
          TRACEBACK; status=1; return
        end if
        
      case default
        write (gol,'("Function GetLonIndex not defined for grid type: ",a)') trim(self%type) ; call goErr
        TRACEBACK; status=1; return
        
        ! NB. For 2d defined lon/lat grids, ilon ant ilat must be found at the same time
        ! Functions getLocation Needed
          
    end select

    ! ok
    status = 0

  end subroutine Grid_Ugg_GetLatIndex

  
  ! *  


  !
  ! Fill indices (ilon,ilat) of grid cell including point (lon,lat).
  ! Return status -1 if location is outside domain.
  !

  subroutine Grid_Ugg_GetLocation( self, lon, lat, ilon, ilat, status, quiet )
  
    use GO, only : T_Vector
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(in)                ::  lon  ! [deg]
    real, intent(in)                ::  lat  ! [deg]
    integer, intent(out)            ::  ilon
    integer, intent(out)            ::  ilat
    integer, intent(out)            ::  status
    logical, intent(in),optional    ::  quiet

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_GetLocation'
    
    ! --- local ----------------------------------
    
    logical                 ::  verbose
    integer                 ::  i, j
    type(T_Vector)          ::  v
    integer                 ::  ibox(4)
    logical                 ::  inside
    
    ! --- begin ----------------------------------
    
    ! error output?
    verbose = .false.
    if ( present(quiet) ) verbose = .not. quiet
    
    ! init result:
    ilat = -999
    ilon = -999
    ! check grid type:
    select case ( trim(self%type) )

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'cartesian', 'cartesian-regular' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        ! Longitude index:
        ! loop over  bands:
        do i = 1, self%nlon
          ! within bouonds ?
          if ( ((lon       >= minval(self%longitude_bnds_1d(:,i))) .and. (lon       < maxval(self%longitude_bnds_1d(:,i)))) .or. &
               ((lon-360.0 >= minval(self%longitude_bnds_1d(:,i))) .and. (lon-360.0 < maxval(self%longitude_bnds_1d(:,i))))      ) then
            ! fill result:
            ilon = i
            ! leave:
            exit
          end if
        end do
        ! check ..
        if ( ilon < 0 ) then
          if ( verbose ) then
            write (gol,'("longitude ",f12.6," not covered by grid")') lon; call goErr
            write (gol,'("Gridtype: ",a )' ) trim(self%type) ; call goErr
            TRACEBACK
          endif
          status=-1;return
        end if
        
        ! Latitude index:
        ! loop over  bands:
        do j = 1, self%nlat
          ! within bouonds ?
          if ( (lat >= minval(self%latitude_bnds_1d(:,j))) .and. (lat < maxval(self%latitude_bnds_1d(:,j))) ) then
            ! fill result:
            ilat = j
            ! leave:
            exit
          end if
        end do
        ! check ..
        if ( ilat < 0 ) then
          if ( verbose ) then
            write (gol,'("latitude ",f12.6," not covered by grid")') lat; call goErr
            write (gol,'("Gridtype: ",a )' ) trim(self%type) ; call goErr
            TRACEBACK
          end if
          status=-1;return
        end if

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'rotated-pole')
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! define location as vector:
        call v%Init( real(lon,8), real(lat,8), status )
        IF_NOT_OK_RETURN(status=1)
      
        ! quick scan to select a small range of cells:
        call self%GetIndexDomain( (/lon,lon,lat,lat/), ibox, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! loop over cells:
        do j = ibox(3), ibox(4)
          do i = ibox(1), ibox(2)
            ! check if inside:
            call self%pg(i,j)%InsideConvex( v, inside, status )
            IF_NOT_OK_RETURN(status=1)
            ! found?
            if ( inside ) then
              ! store indices:
              ilon = i
              ilat = j
              ! leave:
              exit
            end if ! inside
          end do ! i
          ! leave?
          if ( inside ) exit
        end do
        
        ! check ...
        if ( .not. inside ) then
          if ( verbose ) then
            write (gol,'("location (",f12.6,",",f12.6,") not covered by grid")') lon, lat; call goErr
            write (gol,'("Gridtype: ",a )' ) trim(self%type) ; call goErr
            TRACEBACK
          end if
          status=-1;return
        end if
      
        ! clear:
        call v%Done( status )
        IF_NOT_OK_RETURN(status=1)
        
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        write (gol,'("Function GetLocation not yet defined for grid type: ",a)') trim(self%type) ; call goErr
        TRACEBACK; status=1; return
                  
    end select

    ! ok
    status = 0

  end subroutine Grid_Ugg_GetLocation

  
  ! *  


  !
  ! Return distance between centers of cell (i1,j1) and (i2,j2) in m.
  !

  subroutine Grid_Ugg_GetCellDistance( self, i1,j1, i2,j2, dist, status )
  
    use C3PO_Grid_Tools, only : ll_distance
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    integer, intent(in)             ::  i1, j1
    integer, intent(in)             ::  i2, j2
    real, intent(out)               ::  dist
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_GetCellDistance'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! distance in m:
    dist = ll_distance( self%longitude(i1,j1), self%latitude(i1,j1), &
                        self%longitude(i2,j2), self%latitude(i2,j2) )
    ! ok
    status = 0

  end subroutine Grid_Ugg_GetCellDistance

  
  ! *  


  !
  ! Check if point is within distance (m) of domain;
  ! (lon,lat) in deg.
  !

  subroutine Grid_Ugg_WithinDistance( self, lon, lat, dist, nearby, status )
  
    use C3PO_Grid_Tools, only : ll_distance
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(in)                ::  lon, lat    ! deg
    real, intent(in)                ::  dist        ! m
    logical, intent(out)            ::  nearby
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_WithinDistance'
    
    ! --- local ----------------------------------

    real      ::  bbox(4)  ! (/west,east,south,north/)
    real      ::  lon0, lat0
    real      ::  d
    
    ! --- begin ----------------------------------
    
    ! check grid type:
    select case ( trim(self%type) )
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'cartesian', 'cartesian-regular' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        ! first check if is in domain:
        call self%InDomain( lon, lat, nearby, status )
        IF_NOT_OK_RETURN(status=1)
        ! not in domain? then might be nearby:
        if ( .not. nearby ) then
          ! edges:
          call self%GetBoundingBox( bbox, status )
          IF_NOT_OK_RETURN(status=1)
          ! find coordinates of nearby point on edge:
          lon0 = min( max( bbox(1), lon ), bbox(2)  )
          lat0 = min( max( bbox(3), lat ), bbox(4) )
          ! distance:
          d = ll_distance( lon0, lat0, lon, lat )
          ! within range?
          nearby = d <= dist
        end if

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        write (gol,'("grid type `",a,"` not supported yet")') trim(self%type); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0

  end subroutine Grid_Ugg_WithinDistance


  ! ***
  
  
  !
  ! Return 2D field with 'distance' to (lon,lat) position.
  ! The target poisition might be outside the grid.
  !
  
  subroutine Grid_Ugg_DistanceGrid( self, lon, lat, dist, status )
  
    use Grid_Tools, only : ll_distance
    
    ! --- in/out ---------------------------
    
    class(T_Grid_Ugg), intent(in)    ::  self
    real, intent(in)                 ::  lon, lat    ! deg
    real, intent(out)                ::  dist(:,:)   ! m
    integer, intent(out)             ::  status

    ! --- const -----------------------------

    character(len=*), parameter   ::  rname = mname//'/Grid_Ugg_DistanceGrid'
    
    ! --- local -----------------------------
    
    integer   ::  i, j
    real      ::  d
    
    ! --- begin -----------------------------
    
    ! loop over cells:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! distance:
        dist(i,j) = ll_distance( lon, lat, self%longitude(i,j), self%latitude(i,j) ) ! m
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Grid_Ugg_DistanceGrid

  
  ! *  


  !
  ! Round location to cell resolution;
  ! cell center if inside grid.
  !

  subroutine Grid_Ugg_RoundToResolution( self, lon, lat, lon0, lat0, status )
  
    use C3PO_Grid_Tools, only : ll_distance
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(in)                ::  lon, lat    ! deg
    real, intent(out)               ::  lon0, lat0  ! deg
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_RoundToResolution'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! check grid type:
    select case ( trim(self%type) )
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'cartesian-regular' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        ! round to nearby multiple of dlon:
        lon0 = self%longitude_1d(1) + nint((lon-self%longitude_1d(1))/self%dlon) * self%dlon
        lat0 = self%latitude_1d (1) + nint((lat-self%latitude_1d (1))/self%dlat) * self%dlat

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        write (gol,'("grid type `",a,"` not supported yet")') trim(self%type); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0

  end subroutine Grid_Ugg_RoundToResolution

  
  ! *
  
  
  !
  ! Mapping from self to specified domain [west,east,south,north] ;
  ! returns the grid indices that specify the box: i1:i2,j1:j2
  ! and fracs on same index space:  fracs(i1:i2,j1:j2)
  !  
    
  subroutine Grid_Ugg_GetDistribution_Cartesian( self, west, east, south, north, &
                                                        i1,i2,j1,j2, fracs, &
                                                        status )
  
    use GO             , only : T_Polygon
    use C3PO_Grid_Tools, only : ll_area_frac_deg
    use Num            , only : Interval
    
    ! --- in/out --------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(in)                ::  west, east, south, north  ! [deg]
    integer, intent(out)            ::  i1,i2,j1,j2
    real, pointer                   ::  fracs(:,:)
    integer, intent(out)            ::  status  

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Grid_Ugg_GetDistribution_Cartesian'
    
    ! --- local -----------------------------------
    
    integer               ::  i, j
    integer               ::  ibox(4)
    type(T_Polygon)       ::  pbox

    ! --- begin -----------------------------------
    
    ! check grid type:
    select case ( trim(self%type) )
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'cartesian', 'cartesian-regular' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! Check if cell is in target domain
        if ( west  > self%longitude_bnds_1d(2,self%nlon) .or. &
             east  < self%longitude_bnds_1d(1,1) .or. &
             south > self%latitude_bnds_1d(2,self%nlat) .or. &
             north < self%latitude_bnds_1d(1,1) ) then
          ! outside target domain
          status = -1; return
        endif

        ! init indices at the boundaries of target domain
        i1 = 1
        i2 = self%nlon
        j1 = 1
        j2 = self%nlat

        ! Check west boundary
        do i = 1, self%nlon
          if ( self%longitude_bnds_1d(1,i) <= west .and. self%longitude_bnds_1d(2,i) >= west ) then
            i1 = i
            exit
          end if
        end do    
        ! Check east boundary
        do i = 1, self%nlon
          if ( self%longitude_bnds_1d(1,i) <= east .and. self%longitude_bnds_1d(2,i) >= east ) then
            i2 = i
            exit
          end if
        end do

        ! Check south boundary
        do j = 1, self%nlat
          if ( self%latitude_bnds_1d(1,j) <= south .and. self%latitude_bnds_1d(2,j) >= south ) then
            j1 = j
            exit
          end if
        end do    
        ! Check north boundary
        do j = 1, self%nlat
          if ( self%latitude_bnds_1d(1,j) <= north .and. self%latitude_bnds_1d(2,j) >= north ) then
            j2 = j
            exit
          end if
        end do
      
        ! TESTING ...................
        !
        !! quick scan to select a small range of cells:
        !call self%GetIndexDomain( (/west,east,south,north/), ibox, status )
        !IF_NOT_OK_RETURN(status=1)
        !! compare ...
        !if ( (i1 < ibox(1)) .or. (i2 > ibox(2)) .or. (j1 < ibox(3)) .or. (j2 > ibox(4)) ) then
        !  write (gol,'("TESTING: index range (",i0,":",i0,",",i0,":",i0,") outside ibox (",i0,3(",",i0),")")') &
        !                  i1, i2, j1, j2, ibox; call goErr
        !  TRACEBACK; status=1; return
        !end if
        !
        ! ...........................

        ! storage for fractions; clear current:
        if ( associated(fracs) ) then
          deallocate( fracs, stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if
        ! new:
        allocate( fracs(i1:i2,j1:j2), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! fill area fractions:
        do j = j1, j2
          do i = i1, i2
            ! fraction of box that covers grid cell:
            fracs(i,j) = ll_area_frac_deg( west, east, south, north, &
                                           self%longitude_bnds_1d(1,i), self%longitude_bnds_1d(2,i), &
                                           self%latitude_bnds_1d(1,j), self%latitude_bnds_1d(2,j) )
          end do  ! j
        end do  ! i

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'rotated-pole')
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! polygon from box, init corners counter-clockwise:
        call pbox%Init( (/real(west ,8),real(east ,8),real(east ,8),real(west ,8)/), &
                        (/real(south,8),real(south,8),real(north,8),real(north,8)/), status )
        IF_NOT_OK_RETURN(status=1)       

        ! quick scan to select a small range of cells:
        call self%GetIndexDomain( (/west,east,south,north/), ibox, status )
        IF_NOT_OK_RETURN(status=1)
        ! extract elements:
        i1 = ibox(1)
        i2 = ibox(2)
        j1 = ibox(3)
        j2 = ibox(4)

        ! storage for fractions; clear current:
        if ( associated(fracs) ) then
          deallocate( fracs, stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if
        ! new:
        allocate( fracs(i1:i2,j1:j2), stat=status )
        IF_NOT_OK_RETURN(status=1)
        
        ! fill area fractions:
        do j = j1, j2
          do i = i1, i2
            ! fraction of box that is covered by grid cell:
            call self%pg(i,j)%LonLat_AreaFraction( pbox, fracs(i,j), status )
            IF_NOT_OK_RETURN(status=1)
          end do  ! j
        end do  ! i
        
        ! clear:
        call pbox%Done( status )
        IF_NOT_OK_RETURN(status=1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        write (gol,'("grid type `",a,"` not supported yet")') trim(self%type) ; call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
    
  end subroutine Grid_Ugg_GetDistribution_Cartesian
  

  ! ***

  
  ! indices and weights for bi-linear interpolation

  subroutine Grid_Ugg_GetInterpol( self, lon, lat, n, ii, jj, ww, status )
  
    use GO , only : T_Vector, T_PlotFile
    use Num, only : GetInterpolWeights
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)   ::  self
    real, intent(in)                ::  lon  ! [deg]
    real, intent(in)                ::  lat  ! [deg]
    integer, intent(out)            ::  n       ! number of points in interpolation (1..4)
    integer, intent(out)            ::  ii(4)   ! indices of lon cells
    integer, intent(out)            ::  jj(4)   ! indices of lat cells
    real, intent(out)               ::  ww(4)   ! interpolation weights
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_GetInterpol'
    
    ! --- local ----------------------------------
    
    integer             ::  lon_n
    integer             ::  lon_ii(2)
    real                ::  lon_ww(2)
    integer             ::  lat_n
    integer             ::  lat_jj(2)
    real                ::  lat_ww(2)
    integer             ::  i, j
    
    integer             ::  ibox(4)
    type(T_Vector)      ::  p
    integer             ::  k
    real(8)             ::  cw(4)
    type(T_PlotFile)    ::  pf
    
    ! --- begin ----------------------------------
    
    ! no points selected yet:
    n = 0
    
    ! check gridtype
    select case (trim(self%type) )

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'cartesian', 'cartesian-regular')
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! interpolation weights in lon direction,
        ! might require shift over 360 degrees;
        ! return number of points involved (1 or 2), which points, and weights:
        call GetInterpolWeights( self%longitude_1d, lon, lon_n, lon_ii, lon_ww, status, period=360.0 )
        IF_NOT_OK_RETURN(status=1)

        ! interpolation weights in lat direction;
        ! return number of points involved (1 or 2), which points, and weights:
        call GetInterpolWeights( self%latitude_1d, lat, lat_n, lat_jj, lat_ww, status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over lat bands involved:
        do j = 1, lat_n
          ! loop:
          do i = 1, lon_n
            ! next contribution:
            n = n + 1
            ! copy indices:
            ii(n) = lon_ii(i)
            jj(n) = lat_jj(j)
            ! combined weight:
            ww(n) = lon_ww(i) * lat_ww(j)
          end do ! i
        end do ! j

        !! testing ...
        !print *, 'interpol to ', lon, lat, ' from ', n, ' points:'
        !do i = 1, n
        !  print *, i, ii(i), jj(i), self%longitude(ii(i)), self%latitude(jj(i)), ww(i)
        !end do
        !print *, '  sum of weights : ', sum(ww(1:n))
        !stop 'break after interpol weights'

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'rotated-pole')
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! quick scan to select a small range of cells:
        call self%GetIndexDomain( (/lon,lon,lat,lat/), ibox, status )
        IF_NOT_OK_RETURN(status=1)

        ! create target vector:
        call p%Init( real(lon,8), real(lat,8), status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over boxes between centers:
        do j = ibox(3), max(ibox(3),ibox(4)-1)
          do i = ibox(1), max(ibox(1),ibox(2)-1)
            ! interpolation weights ;
            ! status=-1 if not inside:
            call self%mid_pg(i,j)%CornerWeights( p, cw, status )
            if ( status == 0 ) then
              ! weights defined, copy indices if non zero:
              n = 0
              do k = 1, 4
                if ( cw(k) > 0.0 ) then
                  n = n + 1
                  ii(n) = self%mid_ii(i,j,k)
                  jj(n) = self%mid_jj(i,j,k)
                  ww(n) = cw(k)
                end if
              end do
              ! leave:
              exit
            else if ( status > 0 ) then
              ! error ....
              TRACEBACK; status=1; return
            end if
          end do ! i
          ! found?
          if ( status == 0 ) exit
        end do ! j
        ! check ...
        if ( status /= 0 ) then
          ! write plot file:
          call pf%Init( 'debug.py', status )
          IF_NOT_OK_RETURN(status=1)
          !~ source cells:
          do j = ibox(3), ibox(4)
            do i = ibox(1), ibox(2)
              ! add to plot:
              call self%pg(i,j)%Plot( 'edges', pf, status, kwargs='color="red"' )
              IF_NOT_OK_RETURN(status=1)
              !call self%pg(i,j)%Plot( 'edge-arrows', pf, status, kwargs='color="red"' )
              !IF_NOT_OK_RETURN(status=1)
            end do ! i
          end do ! j
          !~ polygons between source cell centers:
          do j = ibox(3), max(ibox(3),ibox(4)-1)
            do i = ibox(1), max(ibox(1),ibox(2)-1)
              ! add to plot:
              call self%mid_pg(i,j)%Plot( 'edges', pf, status, kwargs='color="green"' )
              IF_NOT_OK_RETURN(status=1)
              !call self%mid_pg(i,j)%Plot( 'edge-arrows', pf, status, kwargs='color="green"' )
              !IF_NOT_OK_RETURN(status=1)
              !! again for info:
              !call self%mid_pg(i,j)%CornerWeights( p, cw, status, debug=.true. )
              !print *, 'xxx1 ', i, j, status
              !print *, 'xxx2 longitude_bnds = ', self%longitude_bnds(:,i,j)
              !print *, '  x2 latitude_bnds  = ', self%latitude_bnds(:,i,j)
            end do ! i
          end do ! j
          call p%Plot( 'point', pf, status, kwargs='marker="o", color="blue"' )
          IF_NOT_OK_RETURN(status=1)
          !~ write:
          call pf%Done( status )
          IF_NOT_OK_RETURN(status=1)
          ! info ...
          write (gol,'("no polygon between centers found holding location")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! clear:
        call p%Done( status )
        IF_NOT_OK_RETURN(status=1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        write (gol, '("Interpolation not defined for gridtype: ",a )' ) trim(self%type) ; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0

  end subroutine Grid_Ugg_GetInterpol

  
  ! *
  
  
  subroutine Grid_Ugg_InDomain_Point( self, lon, lat, indomain, status )
  
    ! --- in/out ---------------------------------
  
    class(T_Grid_Ugg), intent(in)  ::  self
    real, intent(in)               ::  lon
    real, intent(in)               ::  lat
    logical, intent(out)           ::  indomain
    integer, intent(out)           ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_InDomain_Point'
    
    ! --- local ----------------------------------
    
    real          ::  lonX

    ! --- begin -----------------------------
    
    
    ! bring  lon  in  [-180,180.0]
    lonX = modulo(lon,360.0)
    if ( lonX > 180.0 ) lonX = lonX - 360.0
    
    select case ( trim(self%type) )
      case ( 'cartesian', 'cartesian-regular' )
        indomain = ( lonX >= self%longitude_bnds_1d(1,1) .and. (lonX <= self%longitude_bnds_1d(2,self%nlon)) .and. &
                     lat  >= self%latitude_bnds_1d(1,1)  .and. (lat  <= self%latitude_bnds_1d(2,self%nlat)) )
      case default
        write (gol,'("Function GetLonIndex not defined for grid type: ",a)') trim(self%type) ; call goErr
        TRACEBACK; status=1; return
        
    end select

    ! ok
    status = 0

  end subroutine Grid_Ugg_InDomain_Point
  
  
  ! *
  
  
  subroutine Grid_Ugg_InDomain_Cell( self, west,east,south,north, indomain, status )
  
    ! in/out
    class(T_Grid_Ugg), intent(in)  ::  self
    real, intent(in)               ::  west, east
    real, intent(in)               ::  south, north
    logical, intent(out)           ::  indomain
    integer, intent(out)           ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_InDomain_Cell'
    
    ! --- local ----------------------------------
    
    real          ::  westX, eastX

    ! --- begin -----------------------------
    
    
    ! bring  lon  in  [-180,180.0]
    westX = modulo(west,360.0)
    if ( westX > 180.0 ) westX = westX - 360.0
    eastX = modulo(east,360.0)
    if ( eastX > 180.0 ) eastX = eastX - 360.0
    
    ! check ...
    select case ( trim(self%type) )
      case ( 'cartesian', 'cartesian-regular' )
        indomain = ( westX >= self%longitude_bnds_1d(1,1) .and. (eastX <= self%longitude_bnds_1d(2,self%nlon)) .and. &
                     south >= self%latitude_bnds_1d(1,1)  .and. (north <= self%latitude_bnds_1d(2,self%nlat)) )
      case default
        write (gol,'("Function GetLonIndex not defined for grid type: ",a)') trim(self%type) ; call goErr
        TRACEBACK; status=1; return
        
    end select

    ! ok
    status = 0

  end subroutine Grid_Ugg_InDomain_Cell
  
  
  ! =====================================================
  
  !  call AreaOper( ugg, ll, '/' | '*' | '=', 'rad2' | 'm2' )
  
  subroutine Grid_Ugg_AreaOper_2d( self, ll, oper, unit, status )
  
    ! --- in/out ----------------------------------
    
    class(T_Grid_Ugg), intent(in)           ::  self
    real, intent(inout)                     ::  ll(:,:)
    character(len=*), intent(in)            ::  unit, oper    
    integer, intent(out)                    ::  status
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_AreaOper_2d'
    
    ! --- local --------------------------------
    
    integer   ::  j,i
    real      ::  cell_area
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( any(shape(ll) /= (/self%nlon,self%nlat/)) ) then
      write (gol,'("unexpected size of ll grid:")'); call goErr
      write (gol,'("  shape(target) : ",i4," x ",i4)') shape(ll); call goErr
      write (gol,'("  shape grid    : ",i4," x ",i4)') self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if

    ! loop over cells:    
    do j = 1, self%nlat
    do i = 1, self%nlon  
      ! select correct area for cells in this row:
      select case ( unit )
        case ( 'm2' )
          cell_area = self%area(i,j)
        ! case ( 'rad' )
          !cell_area = self%area_rad(i,j)
        case default
          write (gol,'("unknown unit : ",a)') trim(unit); call goErr
          TRACEBACK; status=1; return
      end select
      
      ! assign/mult/div by cell area:
      select case ( oper )
        case ( '=' )
          ll(i,j) = cell_area
        case ( '/' )
          ll(i,j) = ll(i,j) / cell_area
        case ( '*' )
          ll(i,j) = ll(i,j) * cell_area
        case default
          write (gol,'("unknown operation : ",a)') trim(oper); call goErr
          TRACEBACK; status=1; return
      end select
    
    end do ! i
    end do ! j
    
    ! ok
    status = 0
    
  end subroutine Grid_Ugg_AreaOper_2d
  
  ! ***
  
  subroutine Grid_Ugg_AreaOper_3d( self, ll, oper, unit, status )
  
    ! --- in/out ----------------------------------
    
    class(T_Grid_Ugg), intent(in)           ::  self
    real, intent(inout)                     ::  ll(:,:,:)
    character(len=*), intent(in)            ::  oper
    character(len=*), intent(in)            ::  unit
    integer, intent(out)                    ::  status
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/llgrid_AreaOper_3d'
    
    ! --- local --------------------------------
    
    integer   ::  l
    
    ! --- begin ----------------------------------
    
    ! loop over layers
    do l = 1, size(ll,3)
      ! apply 2d operator:
      call self%AreaOper( ll(:,:,l), oper, unit, status )
      IF_NOT_OK_RETURN(status=1)
    end do  ! layers
    
    ! ok
    status = 0
    
  end subroutine Grid_Ugg_AreaOper_3d
  
  ! *
  
  !
  ! Define coordinate variables for grid in NetCDF file.
  ! Actual data will be written by 'Grid_Ugg_PutGrid_NetCDF'
  !
  ! Input:
  !   ncid   :  access to netcdf file
  ! Optional:
  !   subset = (/i1,i2,j1,j2/)   : range of output cells, used for limitted output
  !
  ! Output:
  !   ncd    :  T_Grid_NcDef with dimension ids, variable ids, etc.
  ! Optional:
  !   dimid_lon, dimid_lat  : dimension ids of 1D array, 
  !                           should be used to define output variables
  !
  
  subroutine Grid_Ugg_DefGrid_NetCDF( self, ncd, ncid, status, &
                                      dimid_lon, dimid_lat, subset )
    
    use NetCDF, only : NF90_Def_Dim, NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_REAL, NF90_INT, NF90_DOUBLE
    
    ! --- in/out -------------
    
    class(T_Grid_Ugg), intent(in)           ::  self
    type(T_Grid_NcDef), intent(out)         ::  ncd
    integer, intent(in)                     ::  ncid
    integer, intent(out)                    ::  status
    
    integer, intent(out), optional          ::  dimid_lon
    integer, intent(out), optional          ::  dimid_lat
    integer, intent(in), optional           ::  subset(4)

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_DefGrid_NetCDF'
    
    ! --- local ----------------------------------
    
    integer      ::  dimid_nv
    integer      ::  varid
    
    ! --- begin ----------------------------------
    
    ! subset or full grid?
    if ( present(subset) ) then
      ! check ...
      if ( (subset(1) < 1) .or. (subset(2) > self%nlon) .or. &
           (subset(3) < 1) .or. (subset(4) > self%nlat) ) then
        write (gol,'("subset (",i0,3(",",i0),") not in range (",i0,3(",",i0),")")') &
                 subset, 1, self%nlon, 1, self%nlat; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      ncd%subset = subset
    else
      ! by default all:
      ncd%subset = (/ 1, self%nlon, 1, self%nlat /)
    end if
    ! shape:
    ncd%nx = ncd%subset(2) - ncd%subset(1) + 1
    ncd%ny = ncd%subset(4) - ncd%subset(3) + 1
    
    ! store file unit:
    ncd%ncid = ncid
        
    ! switch:
    select case ( trim(self%type) )
      
      ! cartesian grids, defined by 1D lon/lat axes:
      case ( 'cartesian', 'cartesian-regular' )
        
        ! Define dimensions (longitude, latitude)
        status = NF90_Def_Dim( ncid, 'longitude', ncd%nx, ncd%dimid_x )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'latitude' , ncd%ny, ncd%dimid_y )
        IF_NF90_NOTOK_RETURN(status=1)
        ! extra corners:
        status = NF90_Def_Dim( ncid, 'longitude_crnr', ncd%nx+1, ncd%dimid_xc )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'latitude_crnr' , ncd%ny+1, ncd%dimid_yc )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! longitude     
        status = NF90_Def_Var( ncid, 'longitude', NF90_REAL, (/ncd%dimid_x/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'longitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Lon' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'bounds', 'longitude_bnds' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lon = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'latitude', NF90_REAL, (/ncd%dimid_y/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'latitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Lat' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'bounds', 'latitude_bnds' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lat = varid
        
        ! vertices:
        status = NF90_Def_Dim( ncid, 'nv', self%nv_1d, dimid_nv )
        IF_NF90_NOTOK_RETURN(status=1)

        ! longitude     
        status = NF90_Def_Var( ncid, 'longitude_bnds', NF90_REAL, (/dimid_nv,ncd%dimid_x/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'longitude bounds' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lon_bnds = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'latitude_bnds', NF90_REAL, (/dimid_nv,ncd%dimid_y/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'latitude bounds' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lat_bnds = varid

        ! longitude     
        status = NF90_Def_Var( ncid, 'longitude_crnr', NF90_REAL, (/ncd%dimid_xc/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'longitude corners' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lon_crnr = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'latitude_crnr', NF90_REAL, (/ncd%dimid_yc/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'latitude corners' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lat_crnr = varid
        
        ! no coordinates:
        ncd%coordinates = ''
        
      ! 2D grids
      case ( 'rotated-pole' )
        
        ! Define 1D dimensions:
        status = NF90_Def_Dim( ncid, 'x', ncd%nx, ncd%dimid_x )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'y', ncd%ny, ncd%dimid_y )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'xc', ncd%nx+1, ncd%dimid_xc )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'yc', ncd%ny+1, ncd%dimid_yc )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! ~ 1D variables
        
        ! longitude     
        status = NF90_Def_Var( ncid, 'x', NF90_INT, (/ncd%dimid_x/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'x-coordinate in Cartesian system' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', '1' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'X' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_x = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'y', NF90_INT, (/ncd%dimid_y/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'y-coordinate in Cartesian system' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', '1' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, '_CoordinateAxisType', 'Y' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_y = varid
        
        ! ~ 2D variables
        
        ! auxilary coordinates:
        ncd%coordinates = 'longitude latitude'
        
        ! longitude     
        status = NF90_Def_Var( ncid, 'longitude', NF90_REAL, (/ncd%dimid_x,ncd%dimid_y/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'longitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'bounds', 'longitude_bnds' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lon = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'latitude', NF90_REAL, (/ncd%dimid_x,ncd%dimid_y/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'standard_name', 'latitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'bounds', 'latitude_bnds' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lat = varid
        
        ! vertices:
        status = NF90_Def_Dim( ncid, 'nv', self%nv, dimid_nv )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! longitude     
        status = NF90_Def_Var( ncid, 'longitude_bnds', NF90_REAL, (/dimid_nv,ncd%dimid_x,ncd%dimid_y/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'longitude bounds' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lon_bnds = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'latitude_bnds', NF90_REAL, (/dimid_nv,ncd%dimid_x,ncd%dimid_y/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'latitude bounds' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lat_bnds = varid
        
        ! longitude     
        status = NF90_Def_Var( ncid, 'longitude_crnr', NF90_REAL, (/ncd%dimid_xc,ncd%dimid_yc/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'longitude corners' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lon_crnr = varid

        ! latitude
        status = NF90_Def_Var( ncid, 'latitude_crnr', NF90_REAL, (/ncd%dimid_xc,ncd%dimid_yc/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'long_name', 'latitude corners' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        ! store:
        ncd%varid_lat_crnr = varid

      case default 
        write( gol, '("unsupported grid type `",a,"`")') trim(self%type); call GoErr
        TRACEBACK;status=1;return
    end select
    
    ! return 1D dimensions?
    if ( present(dimid_lon) ) dimid_lon = ncd%dimid_x
    if ( present(dimid_lat) ) dimid_lat = ncd%dimid_y

    ! ok
    status = 0    
  
  end subroutine Grid_Ugg_DefGrid_NetCDF
  
  ! *
  
  ! Add "coordinates" attribute to variable if necessary.
  
  subroutine Grid_Ugg_DefCoor_NetCDF( self, ncd, varid, status )
    
    use NetCDF, only : NF90_Put_Att
    
    ! --- in/out -------------
    
    class(T_Grid_Ugg), intent(in)           ::  self
    type(T_Grid_NcDef), intent(in)          ::  ncd
    integer, intent(inout)                  ::  varid
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_DefCoor_NetCDF'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! coordinates defined?
    if ( len_trim(ncd%coordinates) > 0 ) then
      ! add attribute:
      status = NF90_Put_Att( ncd%ncid, varid, 'coordinates', trim(ncd%coordinates) )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0    
  
  end subroutine Grid_Ugg_DefCoor_NetCDF

  ! *
  
  ! Write coordinate data.

  subroutine Grid_Ugg_PutGrid_NetCDF( self, ncd, status )
    
    use NetCDF, only : NF90_Put_Var
    
    ! --- in/out ---------------------------------
    
    class(T_Grid_Ugg), intent(in)           ::  self
    type(T_Grid_NcDef), intent(in)          ::  ncd
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_PutGrid_NetCDF'
    
    ! --- local ----------------------------------
    
    integer       ::  i1, i2, j1, j2
    integer       ::  i, j
    
    ! --- begin ----------------------------------
    
    ! short:
    i1 = ncd%subset(1)
    i2 = ncd%subset(2)
    j1 = ncd%subset(3)
    j2 = ncd%subset(4)
    
    ! switch:
    select case ( trim(self%type) )

      ! regular
      case ( 'cartesian', 'cartesian-regular' )
        
        ! Fill lon/lat array
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lon, self%longitude_1d(i1:i2) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lat, self%latitude_1d (j1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! Fill lon/lat array
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lon_bnds, self%longitude_bnds_1d(:,i1:i2) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lat_bnds, self%latitude_bnds_1d (:,j1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! Fill lon/lat array
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lon_crnr, self%longitude_corners(i1-1:i2,0) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lat_crnr, self%latitude_corners (0,j1-1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)

      ! 2D
      case ( 'rotated-pole' )
        
        ! Fill 1D cartesian lists:
        status = NF90_Put_Var( ncd%ncid, ncd%varid_x, (/ (i, i=i1, i2) /) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( ncd%ncid, ncd%varid_y, (/ (j, j=j1, j2) /) )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! Fill lon/lat array
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lon, self%longitude(i1:i2,j1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lat, self%latitude (i1:i2,j1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! Fill lon/lat array
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lon_bnds, self%longitude_bnds(:,i1:i2,j1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lat_bnds, self%latitude_bnds (:,i1:i2,j1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! Fill lon/lat array
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lon_crnr, self%longitude_corners(i1-1:i2,j1-1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( ncd%ncid, ncd%varid_lat_crnr, self%latitude_corners (i1-1:i2,j1-1:j2) )
        IF_NF90_NOTOK_RETURN(status=1)

      case default 
        write( gol, '("Define grid structure for type: ", a)') trim(self%type); call GoErr
        TRACEBACK;status=1;return

    end select
    
    ! ok
    status = 0    
  
  end subroutine Grid_Ugg_PutGrid_NetCDF


  ! CTL header
  
  subroutine Grid_Ugg_WriteCtlProjection( self, ctl, status, nhalo_x, nhalo_y, westb, ni, southb, nj )
    
    use Grads_Ctl, only : T_Grads_Ctl
    
    ! --- in/out -------------
    
    class(T_Grid_Ugg), intent(in)           ::  self
    type(T_Grads_Ctl), intent(inout)        ::  ctl
    integer, intent(out)                    ::  status 
    integer, intent(in), optional           ::  nhalo_x
    integer, intent(in), optional           ::  nhalo_y
    real, intent(in), optional              ::  westb
    integer, intent(in), optional           ::  ni
    real, intent(in), optional              ::  southb
    integer, intent(in), optional           ::  nj
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Grid_Ugg_WriteCtlProjection'
    
    ! --- local --------------------------------

    integer                   ::  nhx, nhy
    character(len=256)        ::  line
    
    ! --- begin --------------
    
    ! number of halo cells written?
    nhx = 0
    nhy = 0    
    if (present(nhalo_x) ) nhx = nhalo_x      
    if (present(nhalo_y) ) nhy = nhalo_y
    
    ! Cylindrical projection ?
    if ( associated(self%Cylindrical) ) then
      call self%Cylindrical%WriteCtlProjection(ctl, nhx, nhy, status, westb, ni, southb, nj)
      IF_NOT_OK_RETURN(status=1)
  
    ! Lambert projection ?
    else if ( associated(self%Lambert) ) then
      call self%Lambert%WriteCtlProjection(ctl, status)
      IF_NOT_OK_RETURN(status=1)
    
    ! Other:
    else
      write( gol, '("Unknown projection type to write CTL header: ")') ; call goErr
      TRACEBACK;status=1;return
    end if       

  end subroutine Grid_Ugg_WriteCtlProjection
  

  ! *** Cylindrical Projection functions.
  
  
  subroutine ProjInit_Cylindrical( self, west, south, nlon, nlat, dlon, dlat, status )
    
    
    ! --- in/out -------------
    
    class(T_Cylindrical), intent(out)       ::  self
    real, intent(in)                        ::  west
    real, intent(in)                        ::  south
    integer, intent(in)                     ::  nlon
    integer, intent(in)                     ::  nlat
    real, intent(in)                        ::  dlon
    real, intent(in)                        ::  dlat
    integer, intent(out)                    ::  status 
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ProjInit_Cylindrical'
    
    ! --- local --------------------------------
    
    ! --- begin --------------
    
    ! fill
    self%west  = west
    self%south = south  
    self%nlon  = nlon
    self%nlat  = nlat
    self%dlon  = dlon
    self%dlat  = dlat

    ! ok
    status = 0
    
  end subroutine ProjInit_Cylindrical
  
  ! *
  
  subroutine ProjDone_Cylindrical( self, status )
    
    
    ! --- in/out -------------
    
    class(T_Cylindrical), intent(in)        ::  self
    integer, intent(out)                    ::  status 
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ProjDone_Cylindrical'
    
    ! --- local --------------------------------
    
    ! --- begin --------------
    
    ! ok
    status = 0
    
  end subroutine ProjDone_Cylindrical
  
  ! *
  
  subroutine WriteCtlProjection_Cylindrical( self, ctl, nhx, nhy, status, westb, ni, southb, nj )
    
    use Grads_Ctl, only : T_Grads_ctl
    use Grads_Ctl, only : GrADS_Ctl_XYDef
    
    ! --- in/out -------------
    
    class(T_Cylindrical), intent(in)        ::  self
    type(T_Grads_ctl), intent(inout)        ::  ctl
    integer, intent(in)                     ::  nhx
    integer, intent(in)                     ::  nhy
    real, intent(in), optional              ::  westb
    integer, intent(in), optional           ::  ni
    real, intent(in), optional              ::  southb
    integer, intent(in), optional           ::  nj
    integer, intent(out)                    ::  status 
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/WriteCtlProjection_Cylindrical'
    
    ! --- local --------------------------------
    
    ! --- begin --------------
    
    
    if ( present( westb ) ) then
      ! specific domain for ctl fils is geive
      ! xdef:    
      call GrADS_Ctl_XYDef( ctl, 'XDEF', ni+2*nhx, westb+(0.5-nhx)*self%dlon, self%dlon, status )
      IF_NOT_OK_RETURN(status=1)
      ! ydef:
      call GrADS_Ctl_XYDef( ctl, 'YDEF', nj+2*nhy, southb+(0.5-nhy)*self%dlat, self%dlat, status )
      IF_NOT_OK_RETURN(status=1)
    else 
      ! xdef:    
      call GrADS_Ctl_XYDef( ctl, 'XDEF', self%nlon+2*nhx, self%west +(0.5-nhx)*self%dlon, self%dlon, status )
      IF_NOT_OK_RETURN(status=1)
      ! ydef:
      call GrADS_Ctl_XYDef( ctl, 'YDEF', self%nlat+2*nhy, self%south+(0.5-nhy)*self%dlat, self%dlat, status )
      IF_NOT_OK_RETURN(status=1)
    end if  
    
    ! ok
    status = 0
    
  end subroutine WriteCtlProjection_Cylindrical
  

  ! *** Lambert Projection functions.
  

  subroutine ProjInit_Lambert( self, nx, ny, lonref, latref, iref, jref, &
                                     STrueLat, NTrueLat, slon, dx, dy, status )
    
    
    ! --- in/out -------------
    
    class(T_Lambert), intent(out)       ::  self
    
    integer, intent(in)                 ::  nx
    integer, intent(in)                 ::  ny
    
    real, intent(in)                    ::  lonref
    real, intent(in)                    ::  latref
    
    integer, intent(in)                 ::  iref
    integer, intent(in)                 ::  jref
    
    real, intent(in)                    ::  STrueLat
    real, intent(in)                    ::  NTrueLat
    real, intent(in)                    ::  slon
    
    real, intent(in)                    ::  dx
    real, intent(in)                    ::  dy
    
    integer, intent(out)                ::  status 
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ProjInit_Lambert'
    
    ! --- local --------------------------------
    
    ! --- begin --------------
    
    ! fill
    self%nx     = nx
    self%ny     = ny
    
    self%lonref = lonref
    self%latref = latref  
    
    self%iref   = iref
    self%jref   = jref

    self%STrueLat = STrueLat
    self%NTrueLat = NTrueLat
    self%slon   = slon

    self%dx     = dx
    self%dy     = dy

    ! ok
    status = 0
    
  end subroutine ProjInit_Lambert
  
  ! *
  
  subroutine ProjDone_Lambert( self, status )
    
    
    ! --- in/out -------------
    
    class(T_Lambert), intent(in)            ::  self
    integer, intent(out)                    ::  status 
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ProjDone_Lambert'
    
    ! --- local --------------------------------
    
    ! --- begin --------------
    
    ! ok
    status = 0
    
  end subroutine ProjDone_Lambert
  
  ! *

  subroutine WriteCtlProjection_Lambert( self, ctl, status )

    use Grads_Ctl, only : T_Grads_ctl
    use Grads_Ctl, only : GrADS_Ctl_PDef
    
    ! --- in/out -------------
    
    class(T_Lambert), intent(in)           ::  self
    type(T_Grads_Ctl), intent(inout)       :: ctl
    integer, intent(out)                   ::  status 
    
    ! --- const -----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/WriteCtlProjection_Lambert'
    
    ! --- local --------------------------------
    
    ! --- begin --------------
          
    ! pdef:    
    call GrADS_Ctl_PDef( ctl, 'PDEF', 'LAMBERT', &
                              self%nx, self%ny, self%lonref, self%latref, &
                              self%iref, self%jref, self%STrueLat, self%NTrueLat, &
                              self%slon, self%dx, self%dy, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine WriteCtlProjection_Lambert
  
  
  
end module C3PO_Grid_Ugg
