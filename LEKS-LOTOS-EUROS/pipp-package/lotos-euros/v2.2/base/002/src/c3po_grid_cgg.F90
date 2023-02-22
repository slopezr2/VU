!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
! DESCRIPTION
!
!    Type to define carthesian gaussian grid.
!
! EXAMPLE
!
!    ! variables:
!    type(T_Grid_Cgg)      ::  cgg
!    integer               ::  ilon, ilat
!    integer               ::  n       ! number of points in interpolation (1..4)
!    integer               ::  ii(4)   ! indices of points in grid
!    real                  ::  ww(4)   ! interpolation weights
!    real                  ::  value
!    real                  ::  data(:,:) ! input field of size (nlon,nlat)
!    integer               ::  k
!
!    ! Define carthesian grid using lons/lats and cell boundaries.
!    !   real, intent(in)                ::  lons(:)          ! (  nlon)
!    !   real, intent(in)                ::  lons_bnds(:,:)   ! (2,nlon)
!    !   real, intent(in)                ::  lats(:)          ! (  nlat)
!    !   real, intent(in)                ::  lats_bnds(:,:)   ! (2,nlat)
!    !
!    call cgg%Init( lons, lons_bnds, lats, lats_bnds, status )
!    if (status/=0) stop
!
!    ! indices and weights for interpolation to (lon,lat) location:
!    call cgg%GetInterpol( 5.2, 52.4, n, ii, ww, status )
!    if (status/=0) stop
!    ! apply interpolation to data field:
!    value = 0.0
!    do k = 1, n
!      value = value + data(ii(k)) * ww(k)
!    end do
!
!    ! done with grid:
!    call cgg%Done( status )
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

module C3PO_Grid_Cgg

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Grid_Cgg
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Grid_Cgg'
  

  ! --- types ----------------------------------------
  
  type T_Grid_Cgg
    integer                 ::  nlon                  ! number of lons
    integer                 ::  nlat                  ! number of lats
    integer                 ::  nv                    ! number of boundary values (2)
    real, allocatable       ::  longitude(:)          ! (nlon) lons
    real, allocatable       ::  longitude_bnds(:,:)   ! (nv,nlon) lon bounds
    real, allocatable       ::  latitude(:)           ! (nlon) lats
    real, allocatable       ::  latitude_bnds(:,:)    ! (nv,nlat) lat bounds
    real, allocatable       ::  area(:,:)             ! (nlon,nlat) area [m2]
  contains
    procedure   ::  Init         => Grid_Cgg_Init
    procedure   ::  Done         => Grid_Cgg_Done
    procedure   ::  GetLonIndex  => Grid_Cgg_GetLonIndex
    procedure   ::  GetLatIndex  => Grid_Cgg_GetLatIndex
    procedure   ::  GetInterpol  => Grid_Cgg_GetInterpol
  end type T_Grid_Cgg

  

contains


  ! ********************************************************************
  ! ***
  ! *** cgg
  ! ***
  ! ********************************************************************


  subroutine Grid_Cgg_Init( self, lons, lons_bnds, lats, lats_bnds, status )
  
    use C3PO_Grid_Tools, only : ll_area_deg_m2
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Cgg), intent(out)  ::  self
    real, intent(in)                ::  lats(:)          ! (   nlat)
    real, intent(in)                ::  lats_bnds(:,:)   ! (nv,nlat)
    real, intent(in)                ::  lons(:)          ! (   nlon)
    real, intent(in)                ::  lons_bnds(:,:)   ! (nv,nlon)
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Cgg_Init'
    
    ! --- local ----------------------------------
    
    integer   ::  i, j
    
    ! --- begin ----------------------------------
    
    ! size:
    self%nlon = size(lons)
    self%nlat = size(lats)
    self%nv   = 2
    
    ! check ...
    if ( any( shape(lons_bnds) /=  (/self%nv,self%nlon/) ) ) then
      write (gol,'("shape of argument lons_bnds is ",2i6," not ",2i6)') shape(lons_bnds), self%nv,self%nlon; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( any( shape(lats_bnds) /=  (/self%nv,self%nlat/) ) ) then
      write (gol,'("shape of argument lats_bnds is ",2i6," not ",2i6)') shape(lats_bnds), self%nv,self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( self%longitude(self%nlon), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%longitude_bnds(self%nv,self%nlon), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude(self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude_bnds(self%nv,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%area(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy:
    self%longitude      = lons
    self%longitude_bnds = lons_bnds
    self%latitude       = lats
    self%latitude_bnds  = lats_bnds
    
    ! fill area:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! area in m2:
        self%area(i,j) = ll_area_deg_m2( self%longitude_bnds(1,i), self%longitude_bnds(2,i), &
                                          self%latitude_bnds(1,j),  self%latitude_bnds(2,j) )
      end do ! i
    end do ! j

    ! ok
    status = 0
  
  end subroutine Grid_Cgg_Init


  ! ***
  
  
  subroutine Grid_Cgg_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Cgg), intent(inout)    ::  self
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Cgg_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%longitude, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%longitude_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%latitude, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%latitude_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%area, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Grid_Cgg_Done
  
  
  ! ***

  
  ! index of band including longitude

  subroutine Grid_Cgg_GetLonIndex( self, lon, ilon, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Cgg), intent(in)   ::  self
    real, intent(in)                ::  lon  ! [deg]
    integer, intent(out)            ::  ilon
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Cgg_GetLonIndex'
    
    ! --- local ----------------------------------
    
    integer                 ::  i
    
    ! --- begin ----------------------------------
    
    ! init result:
    ilon = -999
    ! loop over  bands:
    do i = 1, self%nlon
      ! within bouonds ?
      if ( ((lon       >= minval(self%longitude_bnds(:,i))) .and. (lon       <= maxval(self%longitude_bnds(:,i)))) .or. &
           ((lon-360.0 >= minval(self%longitude_bnds(:,i))) .and. (lon-360.0 <= maxval(self%longitude_bnds(:,i))))      ) then
        ! fill result:
        ilon = i
        ! leave:
        exit
      end if
    end do
    ! check ..
    if ( ilon < 0 ) then
      write (gol,'("longitude ",f12.6," not covered by grid; longitude bounds:")') lon; call goErr
      do i = 1, self%nlon
        write (gol,'(i6,2f12.4)') i, self%longitude_bnds(:,i); call goErr
      end do
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine Grid_Cgg_GetLonIndex


  ! *
  
  ! index of band including latitude

  subroutine Grid_Cgg_GetLatIndex( self, lat, ilat, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Cgg), intent(in)   ::  self
    real, intent(in)                ::  lat  ! [deg]
    integer, intent(out)            ::  ilat
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Cgg_GetLatIndex'
    
    ! --- local ----------------------------------
    
    integer                 ::  j
    
    ! --- begin ----------------------------------
    
    ! init result:
    ilat = -999
    ! loop over  bands:
    do j = 1, self%nlat
      ! within bouonds ?
      if ( (lat >= minval(self%latitude_bnds(:,j))) .and. (lat <= maxval(self%latitude_bnds(:,j))) ) then
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

  end subroutine Grid_Cgg_GetLatIndex


  ! ***

  
  ! indices and weights for bi-linear interpolation

  subroutine Grid_Cgg_GetInterpol( self, lon, lat, n, ii, jj, ww, status )
  
    use Num, only : GetInterpolWeights
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Cgg), intent(in)   ::  self
    real, intent(in)                ::  lon  ! [deg]
    real, intent(in)                ::  lat  ! [deg]
    integer, intent(out)            ::  n       ! number of points in interpolation (1..4)
    integer, intent(out)            ::  ii(4)   ! indices of lon cells
    integer, intent(out)            ::  jj(4)   ! indices of lat cells
    real, intent(out)               ::  ww(4)   ! interpolation weights
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Cgg_GetInterpol'
    
    ! --- local ----------------------------------
    
    integer             ::  lon_n
    integer             ::  lon_ii(2)
    real                ::  lon_ww(2)
    integer             ::  lat_n
    integer             ::  lat_jj(2)
    real                ::  lat_ww(2)
    integer             ::  i, j
    
    ! --- begin ----------------------------------
    
    ! no points selected yet:
    n = 0

    ! interpolation weights in lon direction,
    ! might require shift over 360 degrees:
    call GetInterpolWeights( self%longitude, lon, lon_n, lon_ii, lon_ww, status, period=360.0 )
    IF_NOT_OK_RETURN(status=1)
    ! interpolation weights in lat direction:
    call GetInterpolWeights( self%latitude, lat, lat_n, lat_jj, lat_ww, status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over lat bands involved:
    do j = 1, lat_n
      ! loop:
      do i = 1, lon_n
        ! next contribution:
        n = n + 1
        ii(n) = lon_ii(i)
        jj(n) = lat_jj(j)
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

    ! ok
    status = 0

  end subroutine Grid_Cgg_GetInterpol


end module C3PO_Grid_Cgg
