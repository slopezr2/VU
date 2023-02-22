!###############################################################################
!
! Tools for measurement update.
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
!###############################################################################


module LEKF_Meas_Tools

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------

  private
  
  public    ::  GetPolygonPoints
  public    ::  T_GridMapping
  
  public    ::  GetSpecApply
  
  public    ::  TGainBoundInfo
  public    ::  T_XCorr
  

  ! --- const --------------------------------

  character(len=*), parameter   ::  mname = 'LEKF_Meas_Tools'
    

  ! --- types --------------------------------
  
  ! chem and other correlation info:
  type T_XCorr
    ! decorrelation between obs and species:
    logical, allocatable              ::  analyse_spec(:)  ! (nspec_all)
    logical                           ::  analyse_sia
    ! decorrelation between obs and noise factors:
    logical, allocatable              ::  analyse_noise(:)  ! (nnoise)
  end type T_XCorr

  ! *

  ! gain bounds
  
  type TGainBoundInfo
    !  o list off all cells correlated with a cell (ix,iy)
    integer, allocatable    ::  nc_corr(:,:)
    integer, allocatable    ::  ix_corr(:,:,:)
    integer, allocatable    ::  iy_corr(:,:,:)
    real, allocatable       ::  fc_corr(:,:,:)
  contains
    procedure ::                      gainbound_Init_0d
    procedure ::                      gainbound_Init_2d
    generic   ::  Init            =>  gainbound_Init_0d, &
                                      gainbound_Init_2d
    procedure ::  Done            =>  gainbound_Done
  end type TGainBoundInfo
  

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
  contains
    procedure ::  Init         =>  GridMapping_Init
    procedure ::  Done         =>  GridMapping_Done
    procedure ::  GetWeights   =>  GridMapping_GetWeights
  end type T_GridMapping

contains



  !
  ! Use Heron's formula for aera of triangle given corners.
  !
  
  real function TriangleArea( xx, yy )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(3), yy(3)

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/TriangleArea'
    
    ! --- local ----------------------------------
    
    real                    ::  a, b, c
    real                    ::  s

    ! --- begin ----------------------------------
    
    ! side lengths:
    a = sqrt( (xx(2) - xx(1))**2 + (yy(2) - yy(1))**2 )
    b = sqrt( (xx(3) - xx(2))**2 + (yy(3) - yy(2))**2 )
    c = sqrt( (xx(1) - xx(3))**2 + (yy(1) - yy(3))**2 )
    ! semi-perimeter:
    s = 0.5 * ( a + b + c )
    ! area:
    TriangleArea = sqrt( s * (s-a) * (s-b) * (s-c) )

  end function TriangleArea



  
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
  ! Return the centroids (mass middle point) and relative weight based on area.
  !
  ! Level 0 returns one point.
  !

  recursive subroutine GetTriangleCentroids( xx, yy, level, maxlevel, xxp, yyp, wwp, np, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(:), yy(:)
    integer, intent(in)                  ::  level
    integer, intent(in)                  ::  maxlevel
    real, intent(inout)                  ::  xxp(:), yyp(:)   ! (np+) centroids
    real, intent(inout)                  ::  wwp(:)
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
        wwp(1) = TriangleArea( xx, yy )
      end if
      ! set counter:
      np = 1
      
    else
    
      ! storage:
      allocate( xxt(size(xxp)), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( yyt(size(yyp)), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( wwt(size(wwp)), stat=status )
      IF_NOTOK_RETURN(status=1)
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
          IF_NOTOK_RETURN(status=1)
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
        IF_NOTOK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
        ! new triangle:
        call GetTriangleCentroids( (/xx(2),xm23,xm12/), (/yy(2),ym23,ym12/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOTOK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
        ! new triangle:
        call GetTriangleCentroids( (/xx(3),xm31,xm23/), (/yy(3),ym31,ym23/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOTOK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
        ! new triangle:
        call GetTriangleCentroids( (/xm12,xm23,xm31/), (/ym12,ym23,ym31/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOTOK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
      end if

      ! clear:
      deallocate( xxt, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( yyt, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( wwt, stat=status )
      IF_NOTOK_RETURN(status=1)
    
    end if  ! end recursion?
    
    ! top level?
    if ( level == 0 ) then
      ! nornmalize weights:
      wwp(1:np) = wwp(1:np) / sum(wwp(1:np))
    end if
    
    ! ok
    status = 0
    
  end subroutine GetTriangleCentroids
  
  ! *
  
  !
  ! Allocate arrays (xxp(:),yyp(:)) with locations of points distributed 
  ! in polygon defined by corners ((xx(:),yy(:)).
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
    
    real, intent(in)                     ::  xx(:), yy(:)
    real, allocatable, intent(inout)     ::  xxp(:), yyp(:)
    real, allocatable, intent(inout)     ::  wwp(:)
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
        IF_NOTOK_RETURN(status=1)
        deallocate( yyp, stat=status )
        IF_NOTOK_RETURN(status=1)
        deallocate( wwp, stat=status )
        IF_NOTOK_RETURN(status=1)
      end if
    end if
    ! not allocated ?
    if ( .not. allocated(xxp) ) then
      ! storage:
      allocate( xxp(maxp), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( yyp(maxp), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( wwp(maxp), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! recursive call:
    call GetTriangleCentroids( xx, yy, 0, maxlevel, xxp, yyp, wwp, np, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    if ( np /= maxp ) then
      write (gol,'("filled ",i0," centroids while expected ",i0)') np, maxp; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine GetPolygonPoints
  
  
  ! ====================================================================
  ! ===
  ! === Pixels
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
    
    ! maximum storage:
    allocate( self%ii(self%nlon*self%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%jj(self%nlon*self%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%ww(self%nlon*self%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
              
    ! allocate 2D array with weight sum:
    allocate( self%wmap(self%nlon,self%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine GridMapping_Init
  
  
  ! *
  

  subroutine GridMapping_Done( self, status )
  
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
      IF_NOTOK_RETURN(status=1)
      deallocate( self%yyp, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! source info:
    deallocate( self%ii, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%jj, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%ww, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( self%wmap, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine GridMapping_Done
  
  ! *
  
  !
  ! Assign locations (xxp,yyp) to grid cells (ii,jj),
  ! the weights ww has the fraction of the total number
  ! of points assigned to a cell.
  ! On input the arrays ii/jj/ww should have sufficient size,
  ! number of elements filled on exit is n.
  !
  
  subroutine GridMapping_GetWeights( self, xx, yy, ii, jj, ww, n, status )

    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    real, intent(in)                     ::  xx(:), yy(:)
    integer, pointer                     ::  ii(:), jj(:)
    real, pointer                        ::  ww(:)
    integer, intent(out)                 ::  n
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GetGridWeights'
    
    ! --- local ----------------------------------
    
    integer                 ::  i1, i2, j1, j2
    integer                 ::  ip
    integer                 ::  i, j
    integer                 ::  k
    integer                 ::  np
    real                    ::  wsum

    ! --- begin ----------------------------------
        
    ! devide footprint in triangles or quadrangles,
    ! return centroids that can be used as sampling points,
    ! arrays xxp and yyp are allocated on output to maximum size needed:
    call GetPolygonPoints( xx, yy, self%xxp, self%yyp, self%wwp, status, levels=self%levels )
    IF_NOTOK_RETURN(status=1)
    
    ! number of sampling points:
    np = size(self%xxp)
 
    ! index box:
    i1 = max(         1, ceiling( (minval(self%xxp) - self%west )/self%dlon ) - self%ilon0 )
    i2 = min( self%nlon, ceiling( (maxval(self%xxp) - self%west )/self%dlon ) - self%ilon0 )
    j1 = max(         1, ceiling( (minval(self%yyp) - self%south)/self%dlat ) - self%ilat0 )
    j2 = min( self%nlat, ceiling( (maxval(self%yyp) - self%south)/self%dlat ) - self%ilat0 )
    
    ! reset sum:
    self%wmap = 0    
    ! loop over points:
    do ip = 1, np
      ! target cell:
      i = ceiling( (self%xxp(ip)-self%west )/self%dlon ) - self%ilon0
      j = ceiling( (self%yyp(ip)-self%south)/self%dlat ) - self%ilat0
      ! in range?
      if ( (i >= 1) .and. (i <= self%nlon) .and. (j >= 1) .and. (j <= self%nlat) ) then
        ! increase weights on map:
        self%wmap(i,j) = self%wmap(i,j) + self%wwp(ip)
      end if
    end do
    
    ! total weight sum:
    wsum = sum(self%wwp)

    ! number of cells with contributins:
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
          ! weight:
          self%ww(k) = real(self%wmap(i,j)) / real(wsum)
        end if
      end do
    end do
    
    ! assign pointers:
    ii => self%ii
    jj => self%jj
    ww => self%ww
    
    ! ok
    status = 0
    
  end subroutine GridMapping_GetWeights

  ! =========================================================================
  ! ===
  ! === tracer selection
  ! ===
  ! =========================================================================

  
  ! Given a keyword, determine analysis weights for species.
  ! This is used to analyze selected species only for certain observations.
  ! Specify one of the keywords: 'gasses','aerosols', 'both'
  ! Return value is logical array of length nspec which is .true. if a specie 
  ! is analyzed and .false. if not.
  
  subroutine GetSpecApply( line, apply, apply_to_sia, status )
  
    use Indices, only : nspec_all, specname
    use Indices, only : specmode, NO_AEROSOL_MODE
    use Indices, only : n_sia, ispecs_sia

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)     ::  line
    logical, intent(out)             ::  apply(:)   ! (nspec_all)
    logical, intent(out)             ::  apply_to_sia
    integer, intent(out)             ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GetSpecapply'

    ! --- local ----------------------------------
    
    integer                        ::  ispec
    integer                        ::  k
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(apply) /= nspec_all ) then
      write (gol,'("output apply should have length ",i6," not ",i6)') nspec_all, size(apply); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! switch:
    select case ( trim(line) )
      case ( 'gasses'   ) ; apply = specmode == NO_AEROSOL_MODE
      case ( 'aerosols' ) ; apply = specmode /= NO_AEROSOL_MODE
      case ( 'both'     ) ; apply = .true.
      case default
        write (gol,'("unsupported spec apply keyword `",a,"`")') trim(line); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! apply to sia? used to enable update of aerh2o:
    apply_to_sia = .false.
    ! loop over sia tracers:
    do k = 1, n_sia
      ! global index:
      ispec = ispecs_sia(k)
      ! enabled?
      if ( apply(ispec) ) then
        apply_to_sia = .true.
        exit
      end if
    end do  ! sia specs

    ! info ...
    write (gol,'("LEKF:       spec weigts:")'); call goPr
    write (gol,'("LEKF:         line  : ",a)') trim(line); call goPr
    write (gol,'("LEKF:         number name (apply)")'); call goPr
    do ispec = 1, nspec_all
      write (gol,'("LEKF:         ",i2," ",a," (",l1,")")') ispec, trim(specname(ispec)), apply(ispec); call goPr
    end do
    write (gol,'("LEKF:         apply to sia : ",l1)') apply_to_sia; call goPr
    
    ! ok
    status = 0
    
  end subroutine GetSpecApply
  

  ! =========================================================================
  ! ===
  ! === gainbound
  ! ===
  ! =========================================================================

  
  subroutine gainbound_Init_0d( self, rho, status )
  
    use dims, only : nx, ny
    
    ! --- in/out ----------------------------
    
    class(TGainBoundInfo), intent(out)    ::  self
    real, intent(in)                      ::  rho
    integer, intent(out)                  ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'gainbound_Init_2d'
    
    ! --- local -----------------------------
    
    real, allocatable   ::  rho_2d(:,:)
    
    ! --- begin -----------------------------
    
    ! storage:
    allocate( rho_2d(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! fill
    rho_2d = rho
    
    ! init:
    call self%Init( rho_2d, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( rho_2d, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine gainbound_Init_0d
  
  ! *  
  
  subroutine gainbound_Init_2d( self, rho, status )
  
    use dims   , only : nx, ny
    use LE_Grid, only : ugg
    
    ! --- in/out ----------------------------
    
    class(TGainBoundInfo), intent(out)    ::  self
    real, intent(in)                      ::  rho(nx,ny)
    integer, intent(out)                  ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'gainbound_Init_2d'
    
    ! --- local -----------------------------
    
    integer                 ::  ix, iy
    integer                 ::  i, j
    integer                 ::  mx_corr
    integer                 ::  icell
    real                    ::  dist
    
    ! --- begin -----------------------------
    
    ! setup counter array:
    allocate( self%nc_corr(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! no cells counted yet:    
    self%nc_corr = 0

    ! loop over all cells:
    do iy = 1, ny
      do ix = 1, nx
        ! loop over all possibly correlated cells:
        do j = 1, ny
          do i = 1, nx
            ! distance in m:
            call ugg%GetCellDistance( ix, iy, i, j, dist, status )  ! m
            IF_NOTOK_RETURN(status=1)
            ! convert to km:
            dist = dist * 1.0e-3  ! km
            ! within range (rho in km)
            if ( dist <= 3.5*rho(ix,iy) ) self%nc_corr(ix,iy) = self%nc_corr(ix,iy) + 1
          end do
        end do
        ! end for this cell
      end do
    end do
    
    ! maximum number of correlated cells:
    mx_corr = maxval(self%nc_corr)

    ! indices of correlated cells and weight for all cells:
    allocate( self%ix_corr(nx,ny,mx_corr), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%iy_corr(nx,ny,mx_corr), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%fc_corr(nx,ny,mx_corr), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over all cells:
    do iy = 1, ny
      do ix = 1, nx
        ! no cells filled yet:
        icell = 0
        ! loop over all possibly correlated cells:
        do j = 1, ny
          do i = 1, nx
            ! distance in m:
            call ugg%GetCellDistance( ix, iy, i, j, dist, status )  ! m
            IF_NOTOK_RETURN(status=1)
            ! convert to km:
            dist = dist * 1.0e-3  ! km
            ! within range ?
            if ( dist <= 3.5*rho(ix,iy) ) then
              ! next index:
              icell = icell + 1
              ! check ...
              if ( icell > mx_corr ) stop 'BUG - icell > mx_corr'
              ! store location and weight:
              self%ix_corr(ix,iy,icell) = i
              self%iy_corr(ix,iy,icell) = j
              ! compute correlation factor:
              if ( rho(ix,iy) == 0.0 ) then
                self%fc_corr(ix,iy,icell) = 1.0
              else
                self%fc_corr(ix,iy,icell) = exp( -0.5*(dist/rho(ix,iy))**2 )
              end if
            end if
          end do
        end do
        ! end for this cell
      end do
    end do
    
    ! ok
    status = 0
  
  end subroutine gainbound_Init_2d


  ! ***
  
  
  subroutine gainbound_Done( self, status )
  
    ! --- in/out ----------------------------
    
    class(TGainBoundInfo), intent(inout)    ::  self
    integer, intent(out)                    ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'gainbound_Done'
    
    ! --- local -----------------------------
    
    ! --- begin -----------------------------
    
    ! clear:
    deallocate( self%nc_corr, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%ix_corr, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%iy_corr, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%fc_corr, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine gainbound_Done
  


end module LEKF_Meas_Tools
