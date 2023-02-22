!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
! DESCRIPTION
!
!    Type to define reduced gaussian grid as used by ECMWF models.
!    A reduced grid consists of a number of latitude bands,
!    with on each band a number of longitude points that get less
!    towards the poles (reduced at higher latitudes). The result
!    is a grid that has more or less the same distance between
!    grid points.
!
! EXAMPLE
!
!    ! variables:
!    type(T_Grid_Rgg)      ::  rgg
!    integer               ::  ilat
!    integer               ::  n       ! number of points in interpolation (1..4)
!    integer               ::  ii(4)   ! indices of points in grid
!    real                  ::  ww(4)   ! interpolation weights
!    real                  ::  value
!    real                  ::  data(:) ! input field of size (npoint)
!    integer               ::  k
!
!    ! Define reduced gaussian grid by specifing properties of latitiude bands
!    ! (latitude of center and boundaries, and number of longitudes per band),
!    ! and properties for all points (longitude values and cell bounds).
!    ! Arguments, the size 'nulat' is the number of latitude bands,
!    ! while 'npoints' is the total number of grid points ;
!    ! consistency of argument size is checked:
!    !   real, intent(in)                ::  band_lats(:)          ! (  nulat)
!    !   real, intent(in)                ::  band_lats_bnds(:,:)   ! (2,nulat)
!    !   integer, intent(in)             ::  band_nlon(:)          ! (  nulat)
!    !   real, intent(in)                ::  lons(:)               ! (  npoint)
!    !   real, intent(in)                ::  lons_bnds(:,:)        ! (2,npoint)
!    !
!    call rgg%Init( band_lats, band_lats_bnds, band_nlon, &
!                   lons, lons_bnds, status )
!    if (status/=0) stop
!
!    ! return index of latitude band holding the specified latitude:
!    call rgg%GetBandIndex( 52.5, ilat, status )
!    if (status/=0) stop
!
!    ! indices and weights for interpolation to (lon,lat) location:
!    call rgg%GetInterpol( 5.2, 52.4, n, ii, ww, status )
!    if (status/=0) stop
!    ! apply interpolation to data field:
!    value = 0.0
!    do k = 1, n
!      value = value + data(ii(k)) * ww(k)
!    end do
!
!    ! done with grid:
!    call rgg%Done( status )
!    if (status/=0) stop
!
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#include "c3po.inc"
!
!######################################################################

module C3PO_Grid_Rgg

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Grid_Rgg
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Grid_Rgg'
  

  ! --- types ----------------------------------------
  
  type T_Grid_Rgg
    integer                 ::  nlat                  ! number of lat bands
    integer                 ::  nulon                 ! number of lons in unreduced grid
    integer                 ::  nv                    ! number of boundary values (2)
    real, allocatable       ::  band_lats(:)          ! (nlat) latitude per band
    real, allocatable       ::  band_lats_bnds(:,:)   ! (nv,nlat) latitude bounds per band
    integer, allocatable    ::  band_i0(:)            ! (nlat) zero-based index of first
                                                      ! point on this band in (npoint) arrays
    integer, allocatable    ::  band_nlon(:)          ! (nlat) number of lons per band
    integer                 ::  npoint                ! number of grid points
    real, allocatable       ::  longitude(:)          ! (npoint) lons per point
    real, allocatable       ::  longitude_bnds(:,:)   ! (nv,npoint) lon bounds per point
    real, allocatable       ::  latitude(:)           ! (npoint) lons per point
    real, allocatable       ::  latitude_bnds(:,:)    ! (nv,npoint) lat bounds per point
    integer, allocatable    ::  irgrid(:)             ! (npoint) zero based flat indices ;
                                      !  values provide location in 2D unreduced array:
                                      !    j = irgrid(n) / unlon
                                      !    i = irgrid(n) - unlon * j
    real, allocatable       ::  area(:)          ! (npoint) area [m2] per point
  contains
    procedure   ::  Init         => Grid_Rgg_Init
    procedure   ::  Done         => Grid_Rgg_Done
    procedure   ::  GetBandIndex => Grid_Rgg_GetBandIndex
    procedure   ::  GetInterpol  => Grid_Rgg_GetInterpol
    procedure   ::  GetRegion    => Grid_Rgg_GetRegion
  end type T_Grid_Rgg


contains


  ! ********************************************************************
  ! ***
  ! *** rgg
  ! ***
  ! ********************************************************************


  subroutine Grid_Rgg_Init( self, band_lats, band_lats_bnds, band_nlon, &
                             lons, lons_bnds, status )
  
    use C3PO_Grid_Tools, only : ll_area_deg_m2
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Rgg), intent(out)  ::  self
    real, intent(in)                ::  band_lats(:)          ! (   nulat)
    real, intent(in)                ::  band_lats_bnds(:,:)   ! (nv,nulat)
    integer, intent(in)             ::  band_nlon(:)          ! (   nulat)
    real, intent(in)                ::  lons(:)               ! (   npoint)
    real, intent(in)                ::  lons_bnds(:,:)        ! (nv,npoint)
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Rgg_Init'
    
    ! --- local ----------------------------------
    
    integer   ::  ilon, ilat
    integer   ::  ipoint
    
    ! --- begin ----------------------------------
    
    ! size:
    self%nlat = size(band_lats)
    self%nv   = 2
    self%npoint = size(lons)
    ! maximum unreduced size:
    self%nulon = maxval( band_nlon )
    
    ! check ...
    if ( any( shape(band_lats_bnds) /=  (/self%nv,self%nlat/) ) ) then
      write (gol,'("shape of argument band_lats_bnds is ",2i6," not ",2i6)') shape(band_lats_bnds), self%nv,self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( size(band_nlon) /= self%nlat ) then
      write (gol,'("size of argument band_nlon is ",i6," not ",i6)') size(band_nlon), self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( any( shape(lons_bnds) /=  (/self%nv,self%npoint/) ) ) then
      write (gol,'("shape of argument lons_bnds is ",2i6," not ",2i6)') shape(lons_bnds), self%nv,self%npoint; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( sum(band_nlon) /= self%npoint ) then
      write (gol,'("sum of nlon values is ",i8," while number of grid points is ",i8)') sum(band_nlon), self%npoint; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( self%band_lats(self%nlat) )
    allocate( self%band_lats_bnds(self%nv,self%nlat) )
    allocate( self%band_i0(self%nlat) )
    allocate( self%band_nlon(self%nlat) )
    allocate( self%longitude(self%npoint) )
    allocate( self%longitude_bnds(self%nv,self%npoint) )
    allocate( self%latitude(self%npoint) )
    allocate( self%latitude_bnds(self%nv,self%npoint) )
    allocate( self%irgrid(self%npoint) )
    allocate( self%area(self%npoint) )
    
    ! copy:
    self%band_lats = band_lats
    self%band_lats_bnds = band_lats_bnds
    self%band_nlon = band_nlon
    ! start indices:
    self%band_i0(1) = 0
    do ilat = 2, self%nlat
      self%band_i0(ilat) = self%band_i0(ilat-1) + self%band_nlon(ilat-1)
    end do
    
    ! copy:
    self%longitude = lons
    self%longitude_bnds = lons_bnds
    
    ! loop over points; update point index every step:
    ipoint = 0
    ! loop over lat bands:
    do ilat = 1, self%nlat
      ! loop over longitudes in this band:
      do ilon = 1, self%band_nlon(ilat)
        ! increase counter:
        ipoint = ipoint + 1
        ! fill coordinates:
        self%latitude(ipoint) = self%band_lats(ilat)
        self%latitude_bnds(:,ipoint) = self%band_lats_bnds(:,ilat)
        ! unpacking index in unreduced 2D array:
        self%irgrid(ipoint) = self%nulon * (ilat-1) + ilon-1
      end do  ! lons in lat band
    end do ! lat bands
    
    ! fill area:
    do ipoint = 1, self%npoint
      ! area in m2:
      self%area(ipoint) = ll_area_deg_m2( self%longitude_bnds(1,ipoint), self%longitude_bnds(2,ipoint), &
                                          self%latitude_bnds (1,ipoint), self%latitude_bnds (2,ipoint) )
    end do

    ! ok
    status = 0
  
  end subroutine Grid_Rgg_Init


  ! ***
  
  
  subroutine Grid_Rgg_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Rgg), intent(inout)    ::  self
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Rgg_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%band_lats )
    deallocate( self%band_lats_bnds )
    deallocate( self%band_i0 )
    deallocate( self%band_nlon )
    deallocate( self%longitude )
    deallocate( self%longitude_bnds )
    deallocate( self%latitude )
    deallocate( self%latitude_bnds )
    deallocate( self%irgrid )
    deallocate( self%area )
    
    ! ok
    status = 0
  
  end subroutine Grid_Rgg_Done
  
  
  ! ***

  
  ! index of band including latitude

  subroutine Grid_Rgg_GetBandIndex( self, lat, ilat, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Rgg), intent(in)   ::  self
    real, intent(in)                ::  lat  ! [deg]
    integer, intent(out)            ::  ilat
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Rgg_GetBandIndex'
    
    ! --- local ----------------------------------
    
    integer                 ::  j
    
    ! --- begin ----------------------------------
    
    ! init result:
    ilat = -999
    ! loop over  bands:
    do j = 1, self%nlat
      ! within bouonds ?
      if ( (lat >= minval(self%band_lats_bnds(:,j))) .and. (lat <= maxval(self%band_lats_bnds(:,j))) ) then
        ! fill result:
        ilat = j
        ! leave:
        exit
      end if
    end do
    ! check ..
    if ( ilat < 0 ) then
      write (gol,'("latitude ",f12.6," not covered by grid")') lat; call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine Grid_Rgg_GetBandIndex


  ! ***

  
  ! indices and weights for bi-linear interpolation

  subroutine Grid_Rgg_GetInterpol( self, lon, lat, n, ii, ww, status )
  
    use Num, only : GetInterpolWeights
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Rgg), intent(in)   ::  self
    real, intent(in)                ::  lon  ! [deg]
    real, intent(in)                ::  lat  ! [deg]
    integer, intent(out)            ::  n       ! number of points in interpolation (1..4)
    integer, intent(out)            ::  ii(4)   ! indices of points in grid
    real, intent(out)               ::  ww(4)   ! interpolation weights
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Rgg_GetInterpol'
    
    ! --- local ----------------------------------
    
    integer             ::  lat_n
    integer             ::  lat_jj(2)
    real                ::  lat_ww(2)
    integer             ::  lon_n
    integer             ::  lon_ii(2)
    real                ::  lon_ww(2)
    integer             ::  j, jlat
    integer             ::  i, ip
    integer             ::  ip1, ip2
    
    ! --- begin ----------------------------------
    
    ! no points selected yet:
    n = 0

    ! interpolation weights in lat direction:
    call GetInterpolWeights( self%band_lats, lat, lat_n, lat_jj, lat_ww, status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over lat bands involved:
    do j = 1, lat_n
      ! band index:
      jlat = lat_jj(j)
      ! cells on band:
      ip1 = self%band_i0(jlat) + 1
      ip2 = self%band_i0(jlat) + self%band_nlon(jlat)
      ! get weights in lon direction:
      call GetInterpolWeights( self%longitude(ip1:ip2), lon, lon_n, lon_ii, lon_ww, status )
      IF_NOT_OK_RETURN(status=1)
      ! loop:
      do i = 1, lon_n
        ! point index:
        ip = ip1-1 + lon_ii(i)
        ! next contribution:
        n = n + 1
        ii(n) = ip
        ww(n) = lat_ww(j) * lon_ww(i)
      end do
    end do

    !! testing ...
    !print *, 'interpol to ', lon, lat, ' from ', n, ' points:'
    !do i = 1, n
    !  print *, '  point ', i, ii(i), self%longitude(ii(i)), self%latitude(ii(i)), ww(i)
    !end do
    !print *, '  sum of weights : ', sum(ww(1:n))
    !stop 'break after interpol weights'

    ! ok
    status = 0

  end subroutine Grid_Rgg_GetInterpol



  ! ***


  subroutine Grid_Rgg_GetRegion( self, domain, rgg, nselected, selection, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Rgg), intent(in)   ::  self
    real, intent(in)                ::  domain(4)  ! (w,e,s,n) [deg]
    class(T_Grid_Rgg), intent(out)  ::  rgg
    integer, intent(out)            ::  nselected       ! number of seleceted points
    integer, intent(out)            ::  selection(:)    ! (npoint) 1-based indices
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Rgg_GetRegion'
    
    ! --- local ----------------------------------
    
    real                    ::  west, east, south, north
    real, allocatable       ::  lons(:), lons_bnds(:,:)
    integer, allocatable    ::  selected_band_nlon(:)
    real, allocatable       ::  selected_lons(:), selected_lons_bnds(:,:)
    integer                 ::  ilon, ilat
    integer                 ::  ilat1, ilat2
    integer                 ::  ipoint
    
    ! --- begin ----------------------------------
    
    ! extract:
    west  = domain(1)
    east  = domain(2)
    south = domain(3)
    north = domain(4)

    ! storage for longitudes harmonized to [0,360) :
    allocate( lons(self%npoint) )
    allocate( lons_bnds(self%nv,self%npoint) )
    ! copy:
    lons      = self%longitude
    lons_bnds = self%longitude_bnds
    ! harmonize to [0,360.0) :
    where ( self%longitude < 0.0 )
      lons           = lons           + 360.0
      lons_bnds(1,:) = lons_bnds(1,:) + 360.0
      lons_bnds(2,:) = lons_bnds(2,:) + 360.0
    endwhere
    ! idem for requested domain:
    if ( west < 0.0 ) west = west + 360.0
    if ( east < 0.0 ) east = east + 360.0
    
    ! temporary storage, maximum size is full grid:
    allocate( selected_band_nlon(self%nlat) )
    allocate( selected_lons(self%npoint) )
    allocate( selected_lons_bnds(self%nv,self%npoint) )
    ! init current size:
    selected_band_nlon = 0
    
    ! init selected indices:
    selection = -999
    nselected = 0
    
    ! init range of lat bands:
    ilat1 = self%nlat + 1000
    ilat2 = 1         - 1000
    ! loop over current bands:
    do ilat = 1, self%nlat
      ! check if band is completely outside domain;
      ! also take bands with edge on border of domain to avoid problems ...
      if ( (maxval(self%band_lats_bnds(:,ilat)) < south) .or. &
           (minval(self%band_lats_bnds(:,ilat)) > north)      ) cycle
      ! band (partly) overlaps domain;
      ! update range:
      if ( ilat < ilat1 ) ilat1 = ilat
      if ( ilat > ilat2 ) ilat2 = ilat;
      ! requested domain lon range might be in [0,360] or [-180,180] ; 
      ! for the later, west > east after modulo 360 :
      if ( west < east ) then
        ! loop over longitudes on this band:
        do ilon = 1, self%band_nlon(ilat)
          ! current original point:
          ipoint = self%band_i0(ilat) + ilon
          ! skip if not covers (edge of) domain:
          if ( (maxval(lons_bnds(:,ipoint)) < west) .or. &
               (minval(lons_bnds(:,ipoint)) > east)      ) cycle
          ! cell (partly) overlaps domain;
          ! update selection index:
          nselected = nselected + 1
          ! update counter:
          selected_band_nlon(ilat) = selected_band_nlon(ilat) + 1
          ! copy values:
          selected_lons     (  nselected) = lons     (  ipoint)
          selected_lons_bnds(:,nselected) = lons_bnds(:,ipoint)
          ! store original index:
          selection(nselected) = ipoint
        end do
      else
        ! select longitudes in [west,360.0) :
        do ilon = 1, self%band_nlon(ilat)
          ! current original point:
          ipoint = self%band_i0(ilat) + ilon
          ! skip if not covers (edge of) domain:
          if ( maxval(lons_bnds(:,ipoint)) < west ) cycle
          ! cell (partly) overlaps domain;
          ! update selection index:
          nselected = nselected + 1
          ! update counter:
          selected_band_nlon(ilat) = selected_band_nlon(ilat) + 1
          ! copy values:
          selected_lons     (  nselected) = lons     (  ipoint) - 360.0
          selected_lons_bnds(:,nselected) = lons_bnds(:,ipoint) - 360.0
          ! store original index:
          selection(nselected) = ipoint
        end do
        ! select longitudes in [0,east] :
        do ilon = 1, self%band_nlon(ilat)
          ! current original point:
          ipoint = self%band_i0(ilat) + ilon
          ! skip if not covers (edge of) domain:
          if ( minval(lons_bnds(:,ipoint)) > east ) cycle
          ! cell (partly) overlaps domain;
          ! update selection index:
          nselected = nselected + 1
          ! update counter:
          selected_band_nlon(ilat) = selected_band_nlon(ilat) + 1
          ! copy values:
          selected_lons     (  nselected) = lons     (  ipoint)
          selected_lons_bnds(:,nselected) = lons_bnds(:,ipoint)
          ! store original index:
          selection(nselected) = ipoint
        end do
      end if  ! lon range accros lon=0
    end do  ! lat bands
    
    ! check ...
    if ( nselected == 0 ) then
      write (gol,'("no grid points selected for domain")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init new region:
    call Grid_Rgg_Init( rgg, self%band_lats(ilat1:ilat2), self%band_lats_bnds(:,ilat1:ilat2), &
                     selected_band_nlon(ilat1:ilat2), &
                     selected_lons(1:nselected), selected_lons_bnds(:,1:nselected), &
                     status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( lons )
    deallocate( lons_bnds )
    deallocate( selected_band_nlon )
    deallocate( selected_lons )
    deallocate( selected_lons_bnds )
    
    ! ok
    status = 0

  end subroutine Grid_Rgg_GetRegion


end module C3PO_Grid_Rgg
