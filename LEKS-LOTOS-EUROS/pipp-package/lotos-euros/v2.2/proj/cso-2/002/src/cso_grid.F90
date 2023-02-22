!###############################################################################
!
! CSO_Grid - tools for regular lon/lat grid
!
!###############################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; csol=nf90_strerror(status); call csoErr; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!###############################################################################


!! testing ...
!module CSO_Logging
!  character(len=1024)   ::  csol
!contains
!  subroutine csoPr()
!    write (*,'(a)') trim(csol)
!  end subroutine csoPr
!  subroutine csoErr()
!    write (*,'("ERROR - ",a)') trim(csol)
!  end subroutine csoErr
!end module CSO_Logging


  
module CSO_Grid

  use CSO_Logging     , only : csol, csoPr, csoErr
  use NetCDF          , only : NF90_StrError, NF90_NOERR

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public    ::  TriangleArea
  public    ::  LonLatTriangleArea
  public    ::  LonLatRectangleArea
  
  public    ::  GetPolygonPoints
  public    ::  T_GridMapping
  public    ::  GridSampleFile
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Grid'
  
  ! earth radius used in area compuations:
  real, parameter         ::  EarthRadius = 6.371e6     ! m
  
  ! something with circles ...
  real, parameter         ::  pi = 3.14159265358979
  ! conversion from degrees to radians:
  real, parameter         ::  deg2rad = pi/180.0     ! rad/deg

  
  ! --- types ------------------------------
  
  type :: T_GridMapping
    ! grid description:
    real                      ::  west, south
    real                      ::  dlon, dlat
    integer                   ::  ilon0, ilat0
    integer                   ::  nlon, nlat
    ! mapping type:
    integer                   ::  levels
    ! sampling points:
    real, allocatable         ::  xxp(:), yyp(:)
    real, allocatable         ::  wwp(:)
    ! source cells, weights:
    integer, pointer          ::  ii(:), jj(:)
    real, pointer             ::  ww(:)
    ! work array:
    real, allocatable         ::  wmap(:,:)  ! (nlon,nlat)
    !
    ! idem for arrays for pixel arrays:
    integer                   ::  nall
    integer                   ::  mxall
    integer, pointer          ::  all_iw0(:), all_nw(:)   ! (npix)
    integer, pointer          ::  all_ii(:), all_jj(:)   ! (mxall)
    real, pointer             ::  all_ww(:)    ! (mxall)
    real, pointer             ::  all_area(:)  ! (npix)
    !
  contains
    procedure ::  Init         =>  GridMapping_Init
    procedure ::  Done         =>  GridMapping_Done
    procedure ::                   GridMapping_GetWeights_0d
    procedure ::                   GridMapping_GetWeights_1d
    generic   ::  GetWeights   =>  GridMapping_GetWeights_0d, &
                                   GridMapping_GetWeights_1d
  end type T_GridMapping

  
contains



  !
  ! Use Heron's formula for area of triangle given corners:
  !
  !   https://en.wikipedia.org/wiki/Heron%27s_formula
  !
  ! Formula below is the standard form:
  ! - sides of triangle: a, b, c
  ! - semi-perimeter:   s = 0.5 * ( a + b + c )
  ! - area:  A^2 =  s * (s-a) * (s-b) * (s-c) 
  ! 
  ! Alternative is to use the "numerical stable" form?
  ! - order: a >= b >= c
  ! - area: (4*A)^2 = (a+(b+c)) * (c - (a-b)) * (c+(a-b)) * (a+(b-c))
  !
  ! However, for problematic "needle" shaped triangles also the second
  ! formula could give slightly negative area^2. 
  ! Because it requires less if-statements the starndard form is used,
  ! and negative area's are truncated to zero.
  !
  
  subroutine TriangleArea( xx, yy, area, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(3), yy(3)
    real, intent(out)                    ::  area
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/TriangleArea'
    
    ! --- local ----------------------------------
    
    real                    ::  a, b, c
    real                    ::  s
    real                    ::  area2

    !real                    ::  qa, qb, qc

    ! --- begin ----------------------------------
    
    ! side lengths:
    a = sqrt( (xx(2) - xx(1))**2 + (yy(2) - yy(1))**2 )
    b = sqrt( (xx(3) - xx(2))**2 + (yy(3) - yy(2))**2 )
    c = sqrt( (xx(1) - xx(3))**2 + (yy(1) - yy(3))**2 )
    ! semi-perimeter:
    s = 0.5 * ( a + b + c )
    ! squared area:
    area2 = s * (s-a) * (s-b) * (s-c)
    ! check ...
    if ( area2 < 0.0 ) then
      ! truncate ...
      area = 0.0
      !! testing ...
      !write (csol,*) 'could not calculate area of triangle:'; call csoErr
      !write (csol,*) '  xx     = ', xx; call csoErr
      !write (csol,*) '  yy     = ', yy; call csoErr
      !write (csol,*) '  area^2 = ', area2; call csoErr
      !TRACEBACK; status=1; return
    else
      ! area:
      area = sqrt( area2 )
    end if
    
    !! sort such that qa >= qb >= qc:
    !!~ initial copy: qa, qb, qc
    !qa = a
    !qb = b
    !qc = c
    !! swap first pair if needed:
    !if ( qa < qb ) then
    !  qa = b
    !  qb = a
    !end if
    !! now:  qa >= qb, qc
    !! need to swap second pair?
    !if ( qb < qc ) then
    !  qc = qb
    !  qb = c
    !  ! now:  qa, qb >= qc
    !  ! need to swap firt pair again?
    !  if ( qa < qb ) then
    !    qb = qa
    !    qa = c
    !  end if
    !  ! now: qa >= qb >= qc
    !end if
    !! compute area:
    !TriangleArea = 0.25*sqrt( (qa+(qb+qc)) * (qc - (qa-qb)) * (qc+(qa-qb)) * (qa+(qb-qc)) )

    ! ok
    status = 0

  end subroutine TriangleArea
  
  !
  ! Area of triangle (xx(:),yy(:)) as (lon,lat) in degrees
  ! with R the earth radius:
  !
  !    area =  R**2 * int cos(y) dy dx
  !                   x,y
  !
  ! Approximate this by first computing the area in degrees**2
  ! and use the average latitude:
  !
  !    area ~  R**2 * int dy dx * cos(y_aver)
  !                   x,y
  !
  
  subroutine LonLatTriangleArea( xx, yy, A, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(3), yy(3)   ! degrees
    real, intent(out)                    ::  A
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LonLatTriangleArea'
    
    ! --- local ----------------------------------
    
    real                    ::  lambda
    real                    ::  Adeg2

    ! --- begin ----------------------------------
    
    ! average latitude:
    lambda = sum(yy)/size(yy) * deg2rad  ! rad
    ! area in deg2:
    call TriangleArea( xx, yy, Adeg2, status )
    IF_NOT_OK_RETURN(status=1)
    ! use triangle area:
    !                            deg2              rad2/deg2                      m2
    A = Adeg2 * deg2rad**2 * cos(lambda) * EarthRadius**2  ! m2

    ! ok
    status = 0

  end subroutine LonLatTriangleArea
  
  !
  ! Approximate area of rectangle in (lon,lat) in degrees using 2 triangles:
  !
  !      4 o---o 3
  !        | / | 
  !      1 o--o 2
  !
  ! Althouth the exact formula is not very difficult, this routine is used
  ! to ensure that area's are approximated in the same way ...
  !
  
  subroutine LonLatRectangleArea( xx, yy, A, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(4), yy(4)   ! degrees
    real, intent(out)                    ::  A
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LonLatRectangleArea'
    
    ! --- local ----------------------------------

    real    ::  A1, A2

    ! --- begin ----------------------------------
    
    ! triangle:
    call LonLatTriangleArea( xx(1:3), yy(1:3), A1, status )
    IF_NOT_OK_RETURN(status=1)
    ! triangle:
    call LonLatTriangleArea( (/xx(3),xx(4),xx(1)/), (/yy(3),yy(4),yy(1)/), A2, status )
    IF_NOT_OK_RETURN(status=1)
    ! combine:
    A = A1 + A2
    
    ! ok
    status = 0

  end subroutine LonLatRectangleArea



  
  !
  ! Arrays xx(:) and yy(:) define corners of a polygon.
  ! If the polygon has 4 or more sides, first devide into triangles:
  !
  !       o-----o
  !       |\  / |
  !       |  o  |
  !       |/  \ |
  !       o-----o
  !
  ! Recursively divide triangles into 4 new triangles 
  ! until level==maxlevel :
  !
  !          o
  !         /_\
  !        /\/\
  !       o----o
  !
  ! Return values:
  !   xxp(:), yyp(:) : centroids (mass middle points) of triangles
  !   wwp(:)         : triangel areas in [(units xx and yy)**2]
  !
  ! Level 0 returns one point.
  !

  recursive subroutine GetTriangleCentroids( xx, yy, level, maxlevel, xxp, yyp, wwp, np, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(:), yy(:)   ! [degree]
    integer, intent(in)                  ::  level
    integer, intent(in)                  ::  maxlevel
    real, intent(inout)                  ::  xxp(:), yyp(:)   ! (np+) centroids   [degree]
    real, intent(inout)                  ::  wwp(:)           ! (np+) triangle area's [m2]
    integer, intent(out)                 ::  np
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GetTriangleCentroids'
    
    ! --- local ----------------------------------
    
    integer                 ::  nc
    integer                 ::  nmax
    integer                 ::  i
    real                    ::  xc, yc
    real                    ::  x1, x2, xm12, xm23, xm31
    real                    ::  y1, y2, ym12, ym23, ym31
    real, allocatable       ::  xxt(:), yyt(:)
    real, allocatable       ::  wwt(:)
    integer                 ::  nt

    ! --- begin ----------------------------------
    
    ! number of corners:
    nc = size(xx)
    
    ! centroid:
    xc = sum(xx)/nc
    yc = sum(yy)/nc
    
    ! reached end?
    if ( level == maxlevel ) then
    
      ! return center:
      xxp(1) = xc
      yyp(1) = yc
      ! single centroid requested (in polygon) ?
      if ( maxlevel == 0 ) then
        ! single centroid, assign dummy weight:
        wwp(1) = 1.0
      else
        ! weight by area:
        call LonLatTriangleArea( xx, yy, wwp(1), status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! set counter:
      np = 1
      
    else
    
      ! storage:
      allocate( xxt(size(xxp)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( yyt(size(yyp)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( wwt(size(wwp)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! no data yet:
      np = 0

      ! devide in triangles?
      if ( nc /= 3 ) then
      
        ! loop over sides:
        do i = 1, nc
          ! start and end point:
          x1 = xx(i)
          y1 = yy(i)
          if ( i == nc ) then
            x2 = xx(1)
            y2 = yy(1)
          else
            x2 = xx(i+1)
            y2 = yy(i+1)
          end if
          ! new triangle:
          call GetTriangleCentroids( (/x1,x2,xc/), (/y1,y2,yc/), level+1, maxlevel, &
                                       xxt, yyt, wwt, nt, status )
          IF_NOT_OK_RETURN(status=1)
          ! extend:
          xxp(np+1:np+nt) = xxt(1:nt)
          yyp(np+1:np+nt) = yyt(1:nt)
          wwp(np+1:np+nt) = wwt(1:nt)
          np = np+nt
        end do ! sides
      
      else
      
        !           3
        !           o
        !         /  \
        !    31 *-----* 23
        !     /  \  /  \
        !  1 o----*-----o 2
        !         12
      
        ! mids:
        xm12 = ( xx(1) + xx(2) )/2
        ym12 = ( yy(1) + yy(2) )/2
        xm23 = ( xx(2) + xx(3) )/2
        ym23 = ( yy(2) + yy(3) )/2
        xm31 = ( xx(3) + xx(1) )/2
        ym31 = ( yy(3) + yy(1) )/2
        
        ! new triangle:
        call GetTriangleCentroids( (/xx(1),xm12,xm31/), (/yy(1),ym12,ym31/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOT_OK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
        ! new triangle:
        call GetTriangleCentroids( (/xx(2),xm23,xm12/), (/yy(2),ym23,ym12/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOT_OK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
        ! new triangle:
        call GetTriangleCentroids( (/xx(3),xm31,xm23/), (/yy(3),ym31,ym23/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOT_OK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
        ! new triangle:
        call GetTriangleCentroids( (/xm12,xm23,xm31/), (/ym12,ym23,ym31/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOT_OK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
      end if

      ! clear:
      deallocate( xxt, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( yyt, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( wwt, stat=status )
      IF_NOT_OK_RETURN(status=1)
    
    end if  ! end recursion?
    
    !! top level?
    !if ( level == 0 ) then
    !  ! nornmalize weights:
    !  wwp(1:np) = wwp(1:np) / sum(wwp(1:np))
    !end if
    
    ! ok
    status = 0
    
  end subroutine GetTriangleCentroids
  
  ! *
  
  !
  ! Given polygon defined by corners [(xx(:),yy(:)],
  ! return arrays [xxp(:),yyp(:)] with locations of points distributed within the polygon;
  ! also return wwp(:) with area corresponding to each point
  !
  ! The polygon is recursively devided in triangles, the points are their centroids.
  ! Each triangle is devided into 4 new triangles:
  !          o
  !         /_\
  !        /\/\
  !       o----o
  ! Recursion is repeated up to "levels" deep;
  ! for a 6-sided polygon the number of points is
  ! 1 (levels==0), 6 (levels=1), 24 (levels==2), 96 (levels=3, default)
  ! With "levels==0" a single point is returned.
  ! Output arrays are allocated at return.
  !
  
  subroutine GetPolygonPoints( xx, yy, xxp, yyp, wwp, status, levels )

    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(:), yy(:)     ! [degree]
    real, allocatable, intent(inout)     ::  xxp(:), yyp(:)   ! [degree]
    real, allocatable, intent(inout)     ::  wwp(:)           ! [m2]
    integer, intent(out)                 ::  status
    integer, intent(in), optional        ::  levels

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GetPolygonPoints'
    
    ! --- local ----------------------------------
    
    integer                 ::  nc
    integer                 ::  ns
    integer                 ::  maxp
    integer                 ::  np
    integer                 ::  i, j
    integer                 ::  maxlevel

    ! --- begin ----------------------------------
    
    ! number of corners:
    nc = size(xx)
    
    ! default:
    maxlevel = 3
    if ( present(levels) ) maxlevel = levels

    ! number of centroids expected:
    do i = 0, maxlevel
      if ( i == 0 ) then
        maxp = 1
      else if ( i == 1 ) then
        maxp = nc
      else
        maxp = maxp * 4
      end if
    end do

    ! already allocated?
    if ( allocated(xxp) ) then
      ! check size:
      if ( size(xxp) /= maxp ) then
        ! clear:
        deallocate( xxp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( yyp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( wwp, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
    end if
    ! not allocated ?
    if ( .not. allocated(xxp) ) then
      ! storage:
      allocate( xxp(maxp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( yyp(maxp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( wwp(maxp), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! recursive call:
    call GetTriangleCentroids( xx, yy, 0, maxlevel, xxp, yyp, wwp, np, status )
    IF_NOT_OK_RETURN(status=1)

    ! check ...
    if ( np /= maxp ) then
      write (csol,'("filled ",i0," centroids while expected ",i0)') np, maxp; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine GetPolygonPoints
  
  
  ! ====================================================================
  ! ===
  ! === GridMapping
  ! ===
  ! ====================================================================


  !
  ! Regular grid defined by:
  ! - south west corner (corner of first grid cell) in global domain
  ! - grid resolution
  ! - grid shape
  ! - cell offset in global domain
  !

  subroutine GridMapping_Init( self, west , dlon, ilon0, nlon, &
                                     south, dlat, ilat0, nlat, &
                                     levels, status )

    use CSO_PArray, only : CSO_PArray_Init
  
    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(out)    ::  self
    real, intent(in)                     ::  west, south
    real, intent(in)                     ::  dlon, dlat
    integer, intent(in)                  ::  ilon0, ilat0
    integer, intent(in)                  ::  nlon, nlat
    integer, intent(in)                  ::  levels
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%west  = west
    self%dlon  = dlon
    self%ilon0 = ilon0
    self%nlon  = nlon
    self%south = south
    self%dlat  = dlat
    self%ilat0 = ilat0
    self%nlat  = nlat
    
    ! store:
    self%levels = levels
    
    ! initialize with maximum storage for single mapping: all cells needed
    allocate( self%ii(self%nlon*self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj(self%nlon*self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ww(self%nlon*self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
              
    ! allocate 2D array with weight sum:
    allocate( self%wmap(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! no arrays for "all" pixels yet:
    self%nall = 0
    self%mxall = 0
    call CSO_PArray_Init( self%all_ii, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%all_jj, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%all_ww, status )
    IF_NOT_OK_RETURN(status=1)
    ! no area per pixel yet:
    call CSO_PArray_Init( self%all_iw0, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%all_nw, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%all_area, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine GridMapping_Init
  
  
  ! *
  

  subroutine GridMapping_Done( self, status )

    use CSO_PArray, only : CSO_PArray_Done
  
    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! sampling points?
    if ( allocated(self%xxp) ) then
      deallocate( self%xxp, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%yyp, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! source info:
    deallocate( self%ii, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%jj, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%ww, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%wmap, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear arrays per pixel:
    call CSO_PArray_Done( self%all_iw0, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Done( self%all_nw, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Done( self%all_area, status )
    IF_NOT_OK_RETURN(status=1)
    ! clear mapping arrays for all pixels:
    call CSO_PArray_Done( self%all_ii, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Done( self%all_jj, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Done( self%all_ww, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine GridMapping_Done
  
  ! *
  
  !
  ! Devide polygon defined by corners (xx(:),yy(:)) in triangles or quadrangles,
  ! and store centroids (xxp(:),yyp(:)) such that these can be used as sampling points.
  ! Assign centroids (xxp(:),yyp(:)) to grid cells (ii(1:n),jj(1:n)),
  ! the weights ww(1:n) have approximately the part of the polygon area
  ! that covers a grid cell; the total polygon area is returned in "area"
  ! Note that the polygon might only partly overlap with the domain, thus:
  !   sum(ww) <= area
  ! On input the arrays ii/jj/ww should have sufficient size,
  ! number of elements filled on exit is n.
  !
  
  subroutine GridMapping_GetWeights_0d( self, xx, yy, area, n, ii, jj, ww, status )

    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    real, intent(in)                     ::  xx(:), yy(:)  ! [degree]
    real, intent(out)                    ::  area          ! [m2]
    integer, intent(out)                 ::  n
    integer, pointer                     ::  ii(:), jj(:)  ! (n)
    real, pointer                        ::  ww(:)         ! (n)  [m2]
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_GetWeights_0d'
    
    ! --- local ----------------------------------
    
    integer                 ::  i1, i2, j1, j2
    integer                 ::  ip
    integer                 ::  i, j
    integer                 ::  k
    integer                 ::  np

    ! --- begin ----------------------------------
        
    ! devide polygon in triangles or quadrangles,
    ! return centroids that can be used as sampling points,
    ! arrays xxp and yyp are allocated on output to maximum size needed:
    call GetPolygonPoints( xx, yy, self%xxp, self%yyp, self%wwp, status, levels=self%levels )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of sampling points:
    np = size(self%xxp)
    
    ! total area of polygon:
    area = sum(self%wwp)
 
    ! index box:
    i1 = max(         1, ceiling( (minval(self%xxp) - self%west )/self%dlon ) - self%ilon0 )
    i2 = min( self%nlon, ceiling( (maxval(self%xxp) - self%west )/self%dlon ) - self%ilon0 )
    j1 = max(         1, ceiling( (minval(self%yyp) - self%south)/self%dlat ) - self%ilat0 )
    j2 = min( self%nlat, ceiling( (maxval(self%yyp) - self%south)/self%dlat ) - self%ilat0 )
    
    ! reset sum:
    self%wmap = 0.0
    ! loop over points:
    do ip = 1, np
      ! target cell:
      i = ceiling( (self%xxp(ip)-self%west )/self%dlon ) - self%ilon0
      j = ceiling( (self%yyp(ip)-self%south)/self%dlat ) - self%ilat0
      ! in range?
      if ( (i >= 1) .and. (i <= self%nlon) .and. (j >= 1) .and. (j <= self%nlat) ) then
        ! increase area on map:
        self%wmap(i,j) = self%wmap(i,j) + self%wwp(ip)
      end if
    end do

    ! number of cells with contributions:
    n = count( self%wmap > 0.0 )
    ! loop over grid cells:
    k = 0
    do i = i1, i2
      do j = j1, j2
        ! any contribution?
        if ( self%wmap(i,j) > 0.0 ) then
          ! increase counter:
          k = k + 1
          ! store location:
          self%ii(k) = i
          self%jj(k) = j
          ! copy total area over this cell:
          self%ww(k) = self%wmap(i,j)
        end if
      end do
    end do
    
    ! assign pointers:
    ii => self%ii
    jj => self%jj
    ww => self%ww
    
    ! ok
    status = 0
    
  end subroutine GridMapping_GetWeights_0d

  ! *
  
  !
  ! Idem for arrays with pixels.
  ! Input:
  !   xx, yy  : pixel footprints
  ! Ouptut:
  !   area       : pixel area
  !   iw0, nw    : per pixel the offset and number of elements in ii/jj/ww arrays
  !   ii, jj     : source cell indices
  !   ww         : source cell weights
  !
  
  subroutine GridMapping_GetWeights_1d( self, xx, yy, &
                                         area, iw0, nw, ii, jj, ww, status )

    use CSO_PArray, only : CSO_PArray_Reshape

    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    real, intent(in)                     ::  xx(:,:), yy(:,:)  ! (ncorner,npix) [degree]
    real, pointer                        ::  area(:)           ! (npix) [m2]
    integer, pointer                     ::  iw0(:)            ! (npix)
    integer, pointer                     ::  nw(:)             ! (npix)
    integer, pointer                     ::  ii(:), jj(:)      ! (nw)
    real, pointer                        ::  ww(:)             ! (nw) [m2]
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_GetWeights_1d'
    
    ! --- local ----------------------------------
    
    integer                   ::  npix
    integer                   ::  ipix
    real                      ::  pix_area
    integer, pointer          ::  pix_ii(:), pix_jj(:)
    real, pointer             ::  pix_ww(:)
    integer                   ::  pix_nw
    integer                   ::  nnew

    ! --- begin ----------------------------------

    ! number of pixels:
    npix = size(xx,2)
    ! check ...
    if ( size(yy,2) /= npix ) then
      write (csol,'("arrays xx (",i0,",",i0,") and yy (",i0,",",i0,") should have same shape")') &
                size(xx,1), size(xx,2), size(yy,1), size(yy,2); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! storage for pixel area; check current storage:
    call CSO_PArray_Reshape( self%all_area, npix, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%all_iw0 , npix, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%all_nw  , npix, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! reset counter:
    self%nall = 0
    self%all_nw = 0
    
    ! loop over pixels:
    do ipix = 1, npix
      ! pixel mapping weights:
      call self%GetWeights( xx(:,ipix), yy(:,ipix), &
                             pix_area, pix_nw, pix_ii, pix_jj, pix_ww, status )
      if ( status /= 0 ) then
        write (csol,'("could not compute mapping weights for pixel ",i0)') ipix; call csoErr
        write (csol,*) '  xx = ', xx(:,ipix); call csoErr
        write (csol,*) '  yy = ', yy(:,ipix); call csoErr
        TRACEBACK; status=ipix; return
      end if

      ! store pixel area:
      self%all_area(ipix) = pix_area

      ! offset and number of overlapping cells (might be 0 ...):
      self%all_iw0 (ipix) = self%nall
      self%all_nw  (ipix) = pix_nw

      ! any overlap? some pixels might not overlap with domain ...
      if ( pix_nw > 0 ) then
        ! exceeds maximum storage?
        if ( self%nall + pix_nw > self%mxall ) then
          ! new size, extend with 1 value extra per cell until it fits ...
          do
            self%mxall = self%mxall + self%nlon*self%nlat
            if ( self%nall + pix_nw <= self%mxall ) exit
          end do
          ! extend arrays, copy current:
          call CSO_PArray_Reshape( self%all_ii , self%mxall, status )
          IF_NOT_OK_RETURN(status=1)
          call CSO_PArray_Reshape( self%all_jj , self%mxall, status )
          IF_NOT_OK_RETURN(status=1)
          call CSO_PArray_Reshape( self%all_ww , self%mxall, status )
          IF_NOT_OK_RETURN(status=1)
        end if
        ! store pixel mapping:
        self%all_ii  (self%nall+1:self%nall+pix_nw) = pix_ii(1:pix_nw)
        self%all_jj  (self%nall+1:self%nall+pix_nw) = pix_jj(1:pix_nw)
        self%all_ww  (self%nall+1:self%nall+pix_nw) = pix_ww(1:pix_nw)
        ! increase counter:
        self%nall = self%nall + pix_nw
      end if ! nw > 0
      
    end do ! ipix
    
    ! truncate to length that is actually used:
    call CSO_PArray_Reshape( self%all_ii , self%nall, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%all_jj , self%nall, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%all_ww , self%nall, status )
    IF_NOT_OK_RETURN(status=1)
    ! reset maximum size:
    self%mxall = self%nall

    ! set pointers:
    area => self%all_area
    iw0  => self%all_iw0
    nw   => self%all_nw
    ii   => self%all_ii
    jj   => self%all_jj
    ww   => self%all_ww
    
    ! ok
    status = 0
    
  end subroutine GridMapping_GetWeights_1d
  
  
  ! ====================================================================
  ! ===
  ! === grid defintion file
  ! ===
  ! ====================================================================
  
  !
  ! Create netcdf file with grid definition.
  ! Used by postprocessing for averaging pixels over grid.
  !
  
  subroutine GridSampleFile( filename, lons, lats, area, status, &
                                 ilon0, ilat0 )
  
    use CSO_NcFile, only : T_NcFile

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)          ::  filename
    real, intent(in)                      ::  lons(:)
    real, intent(in)                      ::  lats(:)
    real, intent(in)                      ::  area(:,:)
    integer, intent(out)                  ::  status
    
    ! offsets in global grid:
    integer, intent(in)                   ::  ilon0, ilat0

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridSampleFile'
    
    ! --- local ----------------------------------
  
    ! grid definition file:
    type(T_NcFile)                    ::  GridFile
    integer                           ::  nlon, nlat
    integer                           ::  ivar_lon, ivar_lat, ivar_area

    ! --- begin ----------------------------------
  
    ! create file:
    call GridFile%Init( trim(filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! global attributes:
    call GridFile%Set_Attr( 0, 'conventions', 'CF-1.7', status )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Set_Attr( 0, 'title', 'CSO grid description', status )
    IF_NOT_OK_RETURN(status=1)
    
    ! shape:
    nlon = size(lons)
    nlat = size(lats)

    ! add dimensions, provide local size and offset in global grid:
    call GridFile%Def_Dim( 'longitude', nlon, status, offset=ilon0 )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Def_Dim( 'latitude', nlat, status, offset=ilat0 )
    IF_NOT_OK_RETURN(status=1)

    ! variable, get internal index back:
    call GridFile%Def_Var( 'longitude', (/'longitude'/), status, ivar=ivar_lon )
    IF_NOT_OK_RETURN(status=1)
    !~ add attributes:
    call GridFile%Set_Attr( ivar_lon, 'standard_name', 'longitude', status )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Set_Attr( ivar_lon, 'units', 'degrees_east', status )
    IF_NOT_OK_RETURN(status=1)

    ! variable, get internal index back:
    call GridFile%Def_Var( 'latitude', (/'latitude'/), status, ivar=ivar_lat )
    IF_NOT_OK_RETURN(status=1)
    !~ add attributes:
    call GridFile%Set_Attr( ivar_lat, 'standard_name', 'latitude', status )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Set_Attr( ivar_lat, 'units', 'degrees_north', status )
    IF_NOT_OK_RETURN(status=1)

    ! variable, get internal index back:
    call GridFile%Def_Var( 'cell_area', (/'longitude','latitude '/), status, ivar=ivar_area )
    IF_NOT_OK_RETURN(status=1)
    !~ add attributes:
    call GridFile%Set_Attr( ivar_area, 'standard_name', 'area', status )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Set_Attr( ivar_area, 'units', 'm2', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call GridFile%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write 1D array, gathered on root from first processor row:
    call GridFile%Put_Var( ivar_lon, lons, status, empty=(ilon0 > 0) )
    IF_NOT_OK_RETURN(status=1)
    ! write 1D array, gathered on root from first processor column:
    call GridFile%Put_Var( ivar_lat, lats, status, empty=(ilat0 > 0) )
    IF_NOT_OK_RETURN(status=1)
    ! write 2D array, gathered on root:
    call GridFile%Put_Var2D( ivar_area, area, status )
    IF_NOT_OK_RETURN(status=1)

    ! close:
    call GridFile%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine GridSampleFile


end module CSO_Grid


! ***


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Testing generation of sample points in footprints
!!!
!!! - Uncomment in this file:
!!!   - the dummy "CSO_Logging" module in the top
!!!   - the "Test" program below
!!!
!!! - In subroutine "GetTriangleCentroids" enable the print statement:
!!!     print *, xc, yc, xx, yy
!!!   and similar in "GetPolygonPoints":
!!!     print *, xxp(np), yyp(np), ...
!!!
!!! - Compile with:
!!!     f90 -o test.x cso_tools.F90
!!!   or after uncommenting some lines in "Makefile":
!!!     make test.x
!!!
!!! - Run with redirection of output:
!!!     ./test.x > test.out
!!!
!!! - Create a python script from the code at the bottom of this file
!!!   and try to make a plot of the footprint.
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!program Test
!
!  use CSO_Logging, only : csol, csoPr, csoErr  
!  use CSO_Tools  , only : GetPolygonPoints
!  
!  implicit None
!    
!  character(len=*), parameter   :: rname = 'Test'
!  
!  integer, parameter    ::  nc = 4
!  real, parameter       ::  xx(nc) = (/0.5,4.5,5.5,1.5/), yy(nc) = (/2.5,1.5,3.5,4.5/)
!  
!  !integer, parameter    ::  nc = 6
!  !real, parameter       ::  xx(nc) = (/0.5,2.5,4.5,5.5,3.5,1.5/), yy(nc) = (/2.5,0.5,1.5,3.5,5.5,4.5/) 
!  
!  integer, parameter    ::  maxlevel = 2
!
!  real, allocatable     ::  xxp(:), yyp(:), wwp(:)
!  integer               ::  status
!
!  ! get points:
!  call GetPolygonPoints( xx, yy, xxp, yyp, wwp, status, levels=maxlevel )
!  IF_NOT_OK_RETURN(status=1)
!
!end program Test
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Test program for plotting.
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  #! /usr/bin/env python
!
!  # read output from 'test.x' program
!  with open('test.out','r') as f : lines = f.readlines()
!
!  # new:
!  fig = plt.figure()
!  ax = fig.add_axes([0.1,0.1,0.8,0.8])
!
!  # rows with:
!  #   xc, yc, xx, yy
!  for line in lines :
!      # convert:
!      values = list(map(float,line.split()))
!      # extract:
!      xc = values[0]
!      yc = values[1]
!      nc = int((len(values)-2)/2)
!      xx = values[2:2+nc]
!      yy = values[2+nc:2+2*nc]
!      # add:
!      ax.plot( xx, yy, linestyle='-', color='k' )
!      ax.plot( xc, yc, marker='o', color='k' )
!  #endfor
!
!  # axes:
!  ax.set_xlim([0,6])
!  ax.set_ylim([0,6])
!
!  # show:
!  plt.show()
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! End
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
