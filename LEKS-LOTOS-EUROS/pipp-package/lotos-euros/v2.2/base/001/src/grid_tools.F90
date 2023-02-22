!
! Grid tools.
!
! may 2002, Arjo Segers
!

module grid_tools


  use Binas, only : deg2rad, pi, ae
  
  implicit none
  
  
  ! --- in/out --------------------------------------
  
  private
  
  public  ::  deg2rad, pi, ae
  
  public  ::  ll_area
  public  ::  ll_area_deg_m2
  public  ::  ll_area_frac
  public  ::  ll_area_frac_deg
  public  ::  ll_distance
  
  

contains


  ! ===================================================================
    
  
  ! given rectangle [west,east]x[south,north] in radians,
  ! compute area in rad^2
  
  real function ll_area( west, east, south, north )
  
    ! --- in/out -------------------------------------
    
    real, intent(in)  ::  west, east, south, north    ! rad
    
    ! --- begin ------------------------------------

    ll_area = ( sin(max(north,south)) - sin(min(north,south)) ) * abs(east-west)
    
  end function ll_area
  
  
  ! *
    
  
  ! given rectangle [west,east]x[south,north] in degrees,
  ! compute area in m^2
  
  real function ll_area_deg_m2( west, east, south, north )
  
    ! --- in/out -------------------------------------
    
    real, intent(in)  ::  west, east, south, north    ! deg
    
    ! --- begin ------------------------------------
    
    ll_area_deg_m2 = ll_area( west*deg2rad, east*deg2rad, south*deg2rad, north*deg2rad ) * ae**2
    
  end function ll_area_deg_m2
  
  
  ! ===
  
  
  ! Compute fraction of rectangle 1 that covers rectangle 2,
  ! both defined in radians.
  ! Fraction is equal to ratio of shaded area and area 1.
  !
  !   +-----------------+         |
  !   |           +-----+---------+--
  !   |    1      |/////|         |
  !   +-----------------+   2     |
  !               |               |
  !             --+---------------+--
  !               |               |
  !
  
  real function ll_area_frac( west1, east1, south1, north1, &
                               west2, east2, south2, north2 )
    
    ! --- in/out ----------------------------
    
    real, intent(in)      ::  west1, east1, south1, north1   ! rad
    real, intent(in)      ::  west2, east2, south2, north2   ! rad

    ! --- local ---------------------------

    real        ::  xwest2, xeast2
    
    ! --- begin -----------------------------

    ! check ..    
    if ( east1 < west1 .or. east2 < west2 .or. &
         north1 < south1 .or. north2 < south2 ) then
      print *, 'found strange area:'
      print *, '  1:', west1, east1, south1, north1
      print *, '  2:', west2, east2, south2, north2
      stop 'FATAL ERROR IN ll_area_frac'
    end if
    
    ! shift rect 2 over 360 deg if it does not cover rect 1 at all;
    ! if they still not cover, the fraction will be set to zero later on:
    if ( west2 > east1 ) then
      ! cell 2 is east of cell 1; try to shift 360.0 westwards
      xwest2 = west2 - 2*pi
      xeast2 = east2 - 2*pi
    else if ( east2 < west1 ) then
      ! cell 2 is west of cell 1; try to shift 360.0 eastwards
      xwest2 = west2 + 2*pi
      xeast2 = east2 + 2*pi
    else
      ! just copy ...
      xwest2 = west2
      xeast2 = east2
    end if

    ! compute fraction; zero if rectangles do not cover:
    if ( (xeast2 <= west1 ) .or. (xwest2 >= east1 ) .or. &
         (north2 <= south1) .or. (south2 >= north1) ) then
      ll_area_frac = 0.0
    else
      ll_area_frac = &
         ll_area( max(west1 ,xwest2), min(east1 ,xeast2), &
                  max(south1,south2), min(north1,north2) ) &
          / ll_area( west1, east1, south1, north1 )
    end if

  end function ll_area_frac
  
  ! *
  
  subroutine ll_area_frac_deg( west1, east1, south1, north1, &
                               west2, east2, south2, north2, &
                               frac, status )
    
    ! --- in/out ----------------------------
    
    real, intent(in)      ::  west1, east1, south1, north1   ! deg
    real, intent(in)      ::  west2, east2, south2, north2   ! deg
    real, intent(out)     ::  frac
    integer, intent(out)  ::  status

    ! --- local ---------------------------

    ! --- begin -----------------------------
    
    ! call function with converted arguments:
    frac = ll_area_frac( west1*deg2rad, east1*deg2rad, south1*deg2rad, north1*deg2rad, &
                         west2*deg2rad, east2*deg2rad, south2*deg2rad, north2*deg2rad )

    ! ok
    status = 0

  end subroutine ll_area_frac_deg
    
  
  ! *


  ! distance in m between two locations in degrees
  
  real function ll_distance( lon0, lat0, lon, lat )
  
    use Binas, only : ae
    use Binas, only : deg2rad
  
    ! --- in/out -------------------------------------
    
    real, intent(in)  ::  lon0, lat0, lon, lat   ! deg
    
    ! --- local ------------------------------------
    
    real      ::  lam0, phi0, lam, phi
    real      ::  v0(3), v(3)
    real      ::  cos_alpha
    
    ! --- begin ------------------------------------
    
    ! location in polar coordinates:
    lam0 = lon0 * deg2rad  ! rad
    phi0 = lat0 * deg2rad  ! rad
    ! vector position:
    !  x-axis from origin to (lon= 0,lat= 0)
    !  y-axis from origin to (lon=90,lat= 0)
    !  z-axis from origin to (lon= 0,lat=90
    v0(1) = cos(lam0) * cos(phi0)
    v0(2) = sin(lam0) * cos(phi0)
    v0(3) =             sin(phi0)

    ! location in polar coordinates:
    lam = lon * deg2rad  ! rad
    phi = lat * deg2rad  ! rad
    ! vector position:
    !  x-axis from origin to (lon= 0,lat= 0)
    !  y-axis from origin to (lon=90,lat= 0)
    !  z-axis from origin to (lon= 0,lat=90
    v(1) = cos(lam) * cos(phi)
    v(2) = sin(lam) * cos(phi)
    v(3) =            sin(phi)
    
    ! angle between two positions:
    !   cos(alpha) = v1.v2 / (|v1||v2|)
    ! here v1 and v2 have both unit length:
    cos_alpha = dot_product( v0, v )
    ! truncate to avoid rounding errors:
    cos_alpha = min( max( -1.0, cos_alpha ), 1.0 )
    
    ! distance over globe:
    !             m         
    ll_distance = ae * acos(cos_alpha)  ! m
    
  end function ll_distance
  

end module grid_tools
  
