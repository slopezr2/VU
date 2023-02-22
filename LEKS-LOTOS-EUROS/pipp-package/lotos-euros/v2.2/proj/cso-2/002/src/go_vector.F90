!###############################################################################
!
! NAME
!   GO_Vector - vector and polygon
!
! REFERENCE
!   Sutherland-Hodgman polygon clipping on Rosetta Code:
!     https://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Fortran
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action)  if (status> 0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#################################################################


module GO_Vector

  use GO_Print, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public    ::  T_Vector
  public    ::  T_Polygon
  

  ! --- const ------------------------------------
  
  character(len=*), parameter   ::  mname = 'GO_Vector'
  
  ! working precission for internal data:
  integer, parameter     ::  wpr = 8
  
    
  ! --- types ------------------------------------
  
  ! 2D vector
  type T_Vector
    ! compontents:
    real(wpr)            ::  x, y
    ! flag to define if valid:
    logical              ::  filled
  contains
    procedure   ::  Init                => Vector_Init
    procedure   ::  Done                => Vector_Done
    procedure   ::  IsEmpty             => Vector_IsEmpty
    procedure   ::  Scale               => Vector_Scale
    procedure   ::  Normalize           => Vector_Normalize
    procedure   ::  LeftFromLine        => Vector_LeftFromLine
    procedure   ::  NormalVector        => Vector_NormalVector
    procedure   ::  String              => Vector_String
    procedure   ::  Plot                => Vector_Plot
  end type T_Vector
  
  !
  ! Polygon defined by vertices.
  ! Fill vertices counter clock-wise to define the polygon to be the interior
  ! enclosed by the edges.
  !
  
  type T_Polygon
    ! number of vertices (= number of edges):
    integer                          ::  n
    ! vertices, one extra for end=start vertex:
    type(T_Vector), allocatable           ::  vertex(:)   ! (n)
  contains
    procedure                                   Polygon_Init_n
    procedure                                   Polygon_Init_xy
    generic     ::  Init                     => Polygon_Init_n, &
                                                Polygon_Init_xy
    procedure   ::  InitClippedPolygon       => Polygon_InitClippedPolygon
    procedure   ::  Done                     => Polygon_Done
    procedure   ::  IsEmpty                  => Polygon_IsEmpty
    procedure   ::  IsTriangle               => Polygon_IsTriangle
    procedure   ::                              Polygon_AddVertex_v
    procedure   ::                              Polygon_AddVertex_xy
    generic     ::  AddVertex                => Polygon_AddVertex_v, &
                                                Polygon_AddVertex_xy
    procedure   ::  EdgeLength               => Polygon_EdgeLength
    procedure   ::  GetNormals               => Polygon_GetNormals
    procedure   ::  Copy                     => Polygon_Copy
    procedure   ::  RemoveVertex             => Polygon_RemoveVertex
    procedure   ::  RemoveTriangle           => Polygon_RemoveTriangle
    procedure   ::  CleanUp                  => Polygon_CleanUp
    procedure   ::  InsideConvex             => Polygon_InsideConvex
    procedure   ::  BoundBox                 => Polygon_BoundBox
    procedure   ::  BoundBoxDiameter         => Polygon_BoundBoxDiameter
    procedure   ::  CopyClippedByLine        => Polygon_CopyClippedByLine
    procedure   ::  Triangles                => Polygon_Triangles
    procedure   ::  LonLat_TriangleArea      => Polygon_LonLat_TriangleArea
    procedure   ::  LonLat_Area              => Polygon_LonLat_Area
    procedure   ::  LonLat_AreaFraction      => Polygon_LonLat_AreaFraction
    procedure   ::  LonLat_ArcLengths        => Polygon_LonLat_ArcLengths
    procedure   ::  CornerWeights            => Polygon_CornerWeights
    procedure   ::  Plot                     => Polygon_Plot
  end type T_Polygon
  
  
  ! --- interfaces -------------------------------------
  
  interface operator(+)
    module procedure vector_plus_vector
  end interface
  
  interface operator(-)
    module procedure vector_minus_vector
  end interface
  
  interface operator(*)
    module procedure scalar_times_vector
  end interface
 

contains


  ! ********************************************************************
  ! ***
  ! *** Vector
  ! ***
  ! ********************************************************************

  !
  ! Init vector using (x,y) coordinates.
  !

  subroutine Vector_Init( self, x, y, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(out)    ::  self
    real(wpr), intent(in)           ::  x, y
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Vector_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%x = x
    self%y = y
    ! set flag:
    self%filled = .true.
    
    ! ok
    status = 0
  
  end subroutine Vector_Init


  ! ***
  
  
  subroutine Vector_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(inout)      ::  self
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Vector_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! reset flag:
    self%filled = .false.
    
    ! ok
    status = 0
  
  end subroutine Vector_Done


  ! ***
  
  
  logical function Vector_IsEmpty( self )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(in)      ::  self
    
    ! --- begin ----------------------------------
    
    ! return flag:
    Vector_IsEmpty = .not. self%filled
  
  end function Vector_IsEmpty


  ! ***


  type(T_Vector) function vector_plus_vector( a, b )
  
    ! --- in/out ---------------------------------
    
    type(T_Vector), intent(in)        ::  a
    type(T_Vector), intent(in)        ::  b
    
    ! --- begin ----------------------------------
    
    ! fill with sum:
    vector_plus_vector%x = a%x + b%x
    vector_plus_vector%y = a%y + b%y
    ! set flag:
    vector_plus_vector%filled = .true.
    
  end function vector_plus_vector


  ! ***


  type(T_Vector) function vector_minus_vector( a, b )
  
    ! --- in/out ---------------------------------
    
    type(T_Vector), intent(in)        ::  a
    type(T_Vector), intent(in)        ::  b
    
    ! --- begin ----------------------------------
    
    ! fill with difference:
    vector_minus_vector%x = a%x - b%x
    vector_minus_vector%y = a%y - b%y
    ! set flag:
    vector_minus_vector%filled = .true.
    
  end function vector_minus_vector


  ! ***


  type(T_Vector) function scalar_times_vector( a, v )
  
    ! --- in/out ---------------------------------
    
    real(wpr), intent(in)             ::  a
    type(T_Vector), intent(in)        ::  v
    
    ! --- begin ----------------------------------
    
    ! fill with multiplication:
    scalar_times_vector%x = a * v%x
    scalar_times_vector%y = a * v%y
    ! set flag:
    scalar_times_vector%filled = .true.
    
  end function scalar_times_vector


  ! ***
  

  !
  ! Vector length in l2 norm:
  !                       x.y
  !   length = sqrt( x**2 + y**2 )
  !
  
  real(wpr) function Length( self )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(in)      ::  self
    
    ! --- begin ----------------------------------
    
    ! compute length:
    Length = sqrt( self%x**2 + self%y**2 )
  
  end function Length


  ! ***
  
  
  !
  ! Scale vector with supplied factor.
  !
  
  subroutine Vector_Scale( self, factor, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(inout)      ::  self
    real(wpr), intent(in)               ::  factor
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Vector_Scale'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! scale:
    self%x = self%x * factor
    self%y = self%y * factor
    
    ! ok
    status = 0
  
  end subroutine Vector_Scale


  ! ***
  
  
  !
  ! Normalize vector to length 1
  !
  
  subroutine Vector_Normalize( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(inout)      ::  self
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Vector_Normalize'
    
    ! --- local ----------------------------------
    
    real(wpr)       ::  l
    
    ! --- begin ----------------------------------
    
    ! current length:
    l = Length( self )
    
    ! check length, should be non-zero:
    if ( l == 0.0 ) then
      write (gol,'("could not normize length-zero vector")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! scale:
    call self%Scale( 1.0/l, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Vector_Normalize
  
  
  ! ***


  !
  ! dot product between vectors
  !
  
  real(wpr) function DotProduct( a, b ) result( c )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(in)      ::  a
    class(T_Vector), intent(in)      ::  b

    ! --- begin ----------------------------------
    
    ! compute cross product:
    c = a%x * b%x + a%y * b%y
  
  end function DotProduct
  
  
  ! ***


  !
  ! cross product between vectors
  !
  
  real(wpr) function CrossProduct( a, b ) result( c )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(in)      ::  a
    class(T_Vector), intent(in)      ::  b

    ! --- begin ----------------------------------
    
    ! compute cross product:
    c = a%x * b%y - a%y * b%x
  
  end function CrossProduct


  ! ***
  

  !
  ! Angle in radians between vectors:
  !                       a.b
  !        cos(alpha) = ------
  !                     |a |b|
  !
  
  real(wpr) function Angle( a, b )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(in)         ::  a
    class(T_Vector), intent(in)         ::  b
    
    ! --- begin ----------------------------------
    
    ! return angle:
    Angle = acos( min( max( -1.0, DotProduct(a,b) / Length(a) / Length(b) ), 1.0 ) )
  
  end function Angle
  
  
  ! ***
  

  !
  ! Return inside=.true. if self point is in half-plane
  ! left of the line (y1 -> y2)
  !
  !                o y2
  !      self *   ^
  !           ^  /
  !        v2 | / v1
  !           |/
  !           o  y1
  !
  
  subroutine Vector_LeftFromLine( self, y1, y2, inside, status )!, debug )
  
    ! --- in/out ---------------------------------
 
    class(T_Vector), intent(in)    ::  self
    class(T_Vector), intent(in)    ::  y1
    class(T_Vector), intent(in)    ::  y2
    logical, intent(out)           ::  inside
    integer, intent(out)           ::  status
    
    !logical, intent(in), optional  ::  debug

    ! --- const --------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Vector_LeftFromLine'
    
    ! --- local ----------------------------------
    
    !logical           ::  dbg
    type(T_Vector)    ::  v1, v2
    
    ! --- begin ----------------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! vector from y1 to y2:
    v1 = y2 - y1
    ! vector from y1 to self point:
    v2 = self - y1
    
    !! tesing ...
    !if ( dbg) then
    !  print *, '  lll1 self = ', self%x, self%y
    !  print *, '    l1 y1   = ', y1%x, y1%y
    !  print *, '    l1 y2   = ', y2%x, y2%y
    !end if
    
    ! cross product (v1 x v2) is positive if v1 should be rotated
    ! counter-clock wise to be aligned with v2 ;
    ! this is true for all points right of the v1,
    ! thus v1 is left from v2 if the crossproduct is positive:
    inside = CrossProduct( v1, v2 ) >= 0.0
    
    ! clear:
    call v1%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call v2%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
 
  end subroutine Vector_LeftFromLine
  
  
  ! ***


  !
  ! Return intersection between finite segment [x1,x2] 
  ! and infinite line through (y1,y2)
  !
  
  type(T_Vector) function Intersection_Segment_Line( x1, x2, y1, y2 ) result ( p )
  
    ! --- in/out ---------------------------------

    type(T_Vector), intent(in)      ::  x1, x2
    type(T_Vector), intent(in)      ::  y1, y2
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Intersection_Segment_Line'
    
    ! --- local ----------------------------------
  
    type(T_Vector)          ::  vx, vy
    type(T_Vector)          ::  x1y1
    real(wpr)               ::  a

    ! --- begin ----------------------------------
    
    ! init result as empty:
    p%filled = .false.
      
    ! vectors from start to end:
    vx = x2 - x1
    vy = y2 - y1
    
    ! vector between start points:
    x1y1 = y1 - x1
 
    ! vectors are colinear ?
    if ( CrossProduct(vx,vy) == 0.0 ) then
      
      ! if the the segment [x1x2] is included in the line (y1y2)
      if ( CrossProduct(x1y1,vx) == 0.d0 ) then
        ! the intersection is the last point of the segment
        p = x2
      end if
      
    else 
    
      ! the vectors are not colinear
      ! we want to find the inersection between [x1x2] and (y1,y2).
      ! mathematically, we want to find a in [0;1] such that :
      !     x1 + a vx    = y1   + b vy        
      ! <=>      a vx    = x1y1 + b vy
      ! <=>      a vx^vy = x1y1^vy      , ^ is cross product
      ! <=>      a       = x1y1^vy / vx^vy
      ! 
      a = CrossProduct(x1y1,vy) / CrossProduct(vx,vy)
      
      ! if a is not in [0,1]
      if ( (a > 1.0) .or. (a < 0.0) ) then
        ! no intersection
      else
        p = x1 + a*vx
      end if

    end if

  end function Intersection_Segment_Line


  ! ***
  
  
  !
  ! Return normal vector with length 1.
  !
  !        ^ self
  !        |
  !        | 
  !        o--> n
  
  subroutine Vector_NormalVector( self, nvec, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(in)         ::  self
    class(T_Vector), intent(out)        ::  nvec
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Vector_NormalVector'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
        
    ! init result:
    call nvec%Init( self%y, -self%x, status )
    IF_NOT_OK_RETURN(status=1)
    ! normalize:
    call nvec%Normalize( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Vector_NormalVector
  
  
  ! ***
  
  
  function Vector_String( self, fmt )
  
    ! --- in/out ---------------------------------
    
    character(len=128)                  ::  Vector_String
    class(T_Vector), intent(in)         ::  self
    character(len=*), intent(in), optional  ::  fmt

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Vector_String'
    
    ! --- local ----------------------------------
    
    character(len=64)     ::  form

    ! --- begin ----------------------------------

    ! fill format:
    if ( present(fmt) ) then
      form = '("(",'//fmt//'",",'//fmt//',")")'
    else
      form = '("(",f10.4,",",f10.4,")")'
    end if
    
    ! formatted string:
    write (Vector_String,form) self%x, self%y
  
  end function Vector_String


  ! ***
  
  
  !
  ! Add python matplotlib lines to text file with unit fu.
  ! Command defines what is plotted from the vector:
  !  'point'    : mark end point
  !  'line'     : line  from origin to end
  !  'arrow'    : arrow from origin to end
  !
  
  subroutine Vector_Plot( self, command, pf, status, kwargs, origin )
  
    use GO_Plot, only : T_PlotFile

    ! --- in/out ---------------------------------
    
    class(T_Vector), intent(in)         ::  self
    character(len=*), intent(in)        ::  command
    type(T_PlotFile), intent(in)        ::  pf
    integer, intent(out)                ::  status
    
    character(len=*), intent(in), optional  ::  kwargs
    class(T_Vector), intent(in), optional   ::  origin

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Vector_Plot'
    
    ! --- local ----------------------------------
    
    character(len=1024)     ::  line
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(command) )
    
      ! end point:
      case ( 'point' )
    
        ! plot command:
        write (line,*) 'ax.plot( [', self%x, '], [', self%y, ']'
        if ( present(kwargs) ) line = trim(line)//', '//trim(kwargs)
        line = trim(line)//' )'
        write (pf%fu,'(a)') '# add vector location:'
        write (pf%fu,'(a)') trim(adjustl(line))
        
      ! line segment:
      case ( 'line' )
    
        ! plot command:
        write (pf%fu,'(a)') '# add line:'
        if ( present(origin) ) then
          write (line,*) 'ax.plot( [', origin%x, ', ', self%x, '], [', origin%y, ', ', self%y, ']'
        else
          write (line,*) 'ax.plot( [0.0, ', self%x, '], [0.0, ', self%y, ']'
        end if
        if ( present(kwargs) ) line = trim(line)//', '//trim(kwargs)
        line = trim(line)//' )'
        write (pf%fu,'(a)') trim(adjustl(line))
        
      ! arrow:
      case ( 'arrow' )
    
        ! plot command:
        write (pf%fu,'(a)') '# add vector arrow:'
        if ( present(origin) ) then
          write (line,*) 'ax.arrow( ', origin%x, ', ', origin%y, ', ', self%x-origin%x, ', ', self%y-origin%y
        else
          write (line,*) 'ax.arrow( 0.0, 0.0, ', self%x, ', ', self%y
        end if
        if ( index(kwargs,'head_width') < 0 ) line = trim(line)//', head_width=0.1'
        line = trim(line)//', length_includes_head=True'
        if ( present(kwargs) ) line = trim(line)//', '//trim(kwargs)
        line = trim(line)//' )'
        write (pf%fu,'(a)') trim(adjustl(line))
        
      ! unknown ...
      case default
        write (gol,'("unsupported command: ",a)') trim(command); call goErr
        TRACEBACK; status=1; return
        
    end select
    
    ! ok
    status = 0
  
  end subroutine Vector_Plot


  ! ********************************************************************
  ! ***
  ! *** Polygon
  ! ***
  ! ********************************************************************
  
  
  !
  ! Init empty polygon, provide maximum number of vertices
  ! (max number of edges plus one for end=start vertex)
  !

  subroutine Polygon_Init_n( self, maxvertex, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(out)             ::  self
    integer, intent(in)                       ::  maxvertex
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_Init_n'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! check ...
    if ( maxvertex < 3 ) then
      write (gol,'("polygon needs storage for at least 3 vertices, requested ",i0)') maxvertex; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! allocate storage:
    allocate( self%vertex(maxvertex), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Polygon_Init_n
  
  
  ! ***
  
  !
  ! Init polygon by two arrays of (x,y) coordinates.
  !

  subroutine Polygon_Init_xy( self, xx, yy, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(out)             ::  self
    real(wpr), intent(in)                     ::  xx(:)
    real(wpr), intent(in)                     ::  yy(:)
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_Init_xy'
    
    ! --- local ----------------------------------
    
    integer         ::  n
    integer         ::  i

    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(xx) /= size(yy) ) then
      write (gol,'("size of xx ",i0," does not match with size of yy ",i0)') size(xx), size(yy); call goErr
      TRACEBACK; status=1; return
    end if
    ! size:
    n = size(xx)
    
    ! init empty storage, one extra for end=start:
    call self%Init( n+1, status )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do i = 1, n
      ! add vertex:
      call self%AddVertex( xx(i), yy(i), status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! ok
    status = 0

  end subroutine Polygon_Init_xy
  
  
  ! ***
  
  
  !
  ! Return copy of input polygon clipped by other polygon.
  ! Sutherland Hodgman algorithm for 2d polygons.
  !

  subroutine Polygon_InitClippedPolygon( self, poly, clip, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(out)             ::  self
    class(T_Polygon), intent(in)              ::  poly
    class(T_Polygon), intent(in)              ::  clip
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_InitClippedPolygon'
    
    ! --- local ----------------------------------
    
    type(T_Polygon)         ::  workPolygon   ! polygon clipped step by step 
    integer                 ::  i1, i2

    ! --- begin ----------------------------------
    
    ! allocate with the maximal possible size,
    ! which is the sum of the size of input polygon and clip:
    call self%Init( poly%n+clip%n, status )
    IF_NOT_OK_RETURN(status=1)
 
    ! allocate workPolygon with the maximal possible size:
    call workPolygon%Init( poly%n+clip%n, status )
    IF_NOT_OK_RETURN(status=1)
 
    !  initialise the work polygon with clip:
    call workPolygon%Copy( clip, status )
    IF_NOT_OK_RETURN(status=1)
 
    ! loop over input polgon edges:
    do i1 = 1, poly%n
      ! end point:
      i2 = modulo(i1,poly%n) + 1
    
      ! fill ouptut polygon with copy of work polygon
      ! clipped by edge i (y1 -> y2):
      call self%CopyClippedByLine( workPolygon, poly%vertex(i1), poly%vertex(i2), status )
      IF_NOT_OK_RETURN(status=1)

      ! workPolygon <= outputPolygon
      call workPolygon%Copy( self, status )
      IF_NOT_OK_RETURN(status=1)
 
    end do ! edges of input polygon
    
    ! clear:
    call workPolygon%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! remove vertices that are very close to others:
    call self%Cleanup( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Polygon_InitClippedPolygon
  
  
  ! ***
  
  
  !
  ! Clear storage
  !

  subroutine Polygon_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_Done'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%vertex, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Polygon_Done
  
  
  ! ***
  
  
  logical function Polygon_IsEmpty( self )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)           ::  self
    
    ! --- begin ----------------------------------
    
    ! at least 3 vertices needed:
    Polygon_IsEmpty = self%n < 3

  end function Polygon_IsEmpty
  
  
  ! ***
  
  
  logical function Polygon_IsTriangle( self )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)           ::  self
    
    ! --- begin ----------------------------------
    
    ! exactly 3+1 vertices needed:
    Polygon_IsTriangle = self%n == 3

  end function Polygon_IsTriangle
  
  
  ! ***
  
  
  !
  ! Add vertex defined by vector.
  !
  
  subroutine Polygon_AddVertex_v( self, v, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    type(T_Vector), intent(in)                ::  v
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_AddVertex_v'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! check ...
    if ( v%IsEmpty() ) then
      write (gol,'("could not add empty vertex")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! add:
    call self%AddVertex( v%x, v%y, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Polygon_AddVertex_v
  
  
  ! ***
  
  
  !
  ! Add vertex defined by (x,y) coordinates
  !
  
  subroutine Polygon_AddVertex_xy( self, x, y, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    real(wpr), intent(in)                     ::  x, y
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_AddVertex_xy'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! increase counter:
    self%n = self%n + 1

    ! check ...
    if ( .not. allocated(self%vertex) ) then
      write (gol,'("no storage allocated in polygon")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( self%n > size(self%vertex) ) then
      write (gol,'("insufficient storage in output polygone: ",i0)') size(self%vertex); call goErr
      TRACEBACK; status=1; return
    end if

    ! store:
    call self%vertex(self%n)%Init( x, y, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Polygon_AddVertex_xy
  
  
  ! ***
  
  
  !
  ! Return length of edge i between vertices i and i+1.
  !
  
  subroutine Polygon_EdgeLength( self, i, length, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    integer, intent(in)                       ::  i
    real(wpr), intent(out)                    ::  length
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_EdgeLength'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
  
    ! check ...
    if ( self%IsEmpty() ) then
      write (gol,'("could not get vertex from empty polygon")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! switch:
    if ( (i < 1) .or. (i > self%n) ) then
      write (gol,'("vertex index ",i0," not in range 1,..,",i0)') i; call goErr
      TRACEBACK; status=1; return
    else if ( i < self%n ) then
      ! distance between vertices i and i+1:
      length = sqrt( ( self%vertex(i+1)%x - self%vertex(i)%x )**2 + &
                     ( self%vertex(i+1)%y - self%vertex(i)%y )**2   )
    else
      ! distance beteen last vertex and first vertex:
      length = sqrt( ( self%vertex(1)%x - self%vertex(i)%x )**2 + &
                     ( self%vertex(1)%y - self%vertex(i)%y )**2   )
    end if
    
    ! ok
    status = 0
    
  end subroutine Polygon_EdgeLength
  
  
  ! ***


  !
  ! Return vector for numbered vertex.
  !
  
  subroutine Polygon_GetNormals( self, xx, yy, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)              ::  self
    real(wpr), intent(out)                    ::  xx(4), yy(4)
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_GetNormals'
    
    ! --- local ----------------------------------
    
    integer         ::  k
    type(T_Vector)  ::  edge
    type(T_Vector)  ::  nvec

    ! --- begin ----------------------------------
  
    ! check ...
    if ( self%IsEmpty() ) then
      write (gol,'("could not get vertex from empty polygon")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop over edges:
    do k = 1, self%n
    
      ! edge vector:
      if ( k < self%n ) then
        edge = self%vertex(k+1) - self%vertex(k)
      else
        edge = self%vertex(1) - self%vertex(k)
      end if
      
      ! normal vector:
      call edge%NormalVector( nvec, status )
      IF_NOT_OK_RETURN(status=1)
      ! copy elements:
      xx(k) = nvec%x
      yy(k) = nvec%y
      
      !print *, '  vvv1 edge ', k, ' = ', edge%x, edge%y
      !print *, '    v1 nvec ', nvec%x, nvec%y
      
    end do  ! edges

    ! ok
    status = 0
    
  end subroutine Polygon_GetNormals
  
  
  ! ***
  
  
  !
  ! Copy content of polygon, storage should be sufficient.
  !

  subroutine Polygon_Copy( self, poly, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    class(T_Polygon), intent(in)              ::  poly
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_Copy'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! check ...
    if ( poly%n > size(self%vertex) ) then
      write (gol,'("insufficient storage in output polygon (",i0,") for input (",i0,")")') &
                         size(self%vertex), poly%n; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! copy target size:
    self%n = poly%n
    ! copy vertices:
    self%vertex(1:self%n) = poly%vertex(1:self%n)
    
    ! ok
    status = 0

  end subroutine Polygon_Copy
  
  
  ! ***
  
  
  !
  ! Remove vertex from polygon
  !

  subroutine Polygon_RemoveVertex( self, i, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    integer, intent(in)                       ::  i
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_RemoveVertex'
    
    ! --- local ----------------------------------
    
    integer        ::  k

    ! --- begin ----------------------------------
  
    ! check ...
    if ( self%IsEmpty() ) then
      write (gol,'("could not get vertex from empty polygon")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! vertex number in 1,..,n:
    k = modulo( i-1, self%n ) + 1
    
    ! remove vertex before last one ?
    if ( k < self%n ) then
      ! shift vertices after k:
      self%vertex(k:self%n-1) = self%vertex(k+1:self%n)
    end if
    ! reset size:
    self%n = self%n - 1

    ! ok
    status = 0

  end subroutine Polygon_RemoveVertex
  
  
  ! ***
  
  
  !
  ! Remove vertices that are very close to another.
  ! First compute a length scale as the maximum length of the bounding box.
  ! Vertices within a distance that is a small fraction 
  ! of this length scale are removed.
  !

  subroutine Polygon_CleanUp( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_CleanUp'
    
    ! threshold fraction for distant
    real, parameter   ::  thrf = 1.0e-4
    
    ! --- local ----------------------------------
    
    real(wpr)   ::  L
    real(wpr)   ::  edgelength
    integer     ::  i

    ! --- begin ----------------------------------
  
    ! filled?
    if ( .not. self%IsEmpty() ) then
    
      ! get length scale:
      call self%BoundBoxDiameter( L, status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over vertices, backwards to first remove the latest:
      do i = self%n, 1, -1
        ! length of edge towards next vertex:
        call self%EdgeLength( i, edgelength, status )
        IF_NOT_OK_RETURN(status=1)
        ! too small?
        if ( edgelength < thrf * L ) then
          ! remove vertex:
          call self%RemoveVertex( i, status )
          IF_NOT_OK_RETURN(status=1)
          ! leave if empty:
          if ( self%IsEmpty() ) exit
        end if  ! close to next?
      end do  ! vertices
    
    end if  ! not empty
    
    ! ok
    status = 0

  end subroutine Polygon_CleanUp


  ! ***
  
  
  !
  ! Cut off triangle with tip at specified vertex.
  ! On input, triangle should have enough storage to store 3 vertices.
  ! Empty triangle is returned if part of it is outside the polygon,
  ! which could happen if it is not convex.
  !
  
  subroutine Polygon_RemoveTriangle( self, i0, tri, status )
  
    !use GO_Plot, only : T_PlotFile

    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    integer, intent(in)                       ::  i0
    class(T_Polygon), intent(inout)           ::  tri
    integer, intent(out)                      ::  status
    
    !type(T_PlotFile), intent(inout)           ::  pf
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_RemoveTriangle'
    
    ! --- local ----------------------------------
    
    !logical           ::  dbg
    integer           ::  i1, i2, i
    logical           ::  vertex_inside

    ! --- begin ----------------------------------
    
    ! reset to empty:
    tri%n = 0

    ! previous vertex:
    i1 = i0 - 1
    if ( i1 == 0 ) i1 = self%n
    ! add to triangle:
    call tri%AddVertex( self%vertex(i1), status )
    IF_NOT_OK_RETURN(status=1)

    ! add current vertex:
    call tri%AddVertex( self%vertex(i0), status )
    IF_NOT_OK_RETURN(status=1)

    ! next vertex:
    i2 = i0 + 1
    if ( i2 == self%n+1 ) i2 = 1
    ! add to triangle:
    call tri%AddVertex( self%vertex(i2), status )
    IF_NOT_OK_RETURN(status=1)

    ! init flag as if none of the other vertices is inside triangle:
    vertex_inside = .false.
    ! loop over all vertices:
    do i = 1, self%n
      ! skip current triangle:
      if ( i == i0 ) cycle
      if ( i == i1 ) cycle
      if ( i == i2 ) cycle
      ! check if inside triangle:
      call tri%InsideConvex( self%vertex(i), vertex_inside, status )
      IF_NOT_OK_RETURN(status=1)
      ! leave if inside ..
      if ( vertex_inside ) exit
    end do  ! other vertices
    
    ! some other vertices inside?
    if ( vertex_inside ) then
      ! return empty result:
      tri%n = 0
    else
      ! remove corner with polygon from triangle:
      call self%RemoveVertex( i0, status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine Polygon_RemoveTriangle
  
  
  ! ***
  
  
  !
  ! Return triangle area in m2 if (x,y) coordinates are to
  ! be interpreted as (lon,lat) in degrees.
  ! Input polygon is first split into a lower and upper triangle
  ! for easier integration in y direction.
  !
  !        y^
  !         |
  !    ymax o
  !         |\
  !         |\\    line  Dx/(ymid-ymax) * (y-ymax)     
  !         |\\\         \_____a______/
  !    ymid +---o
  !         |///   
  !         |//    line  Dx/(ymid-ymin) * (y-ymin) 
  !         |/           \_____a______/
  !    ymin o
  !         |
  !    -----+---+--> x
  !         |   Dx
  !     
  !  Use that:
  !
  !    [(y-ya) sin(y) + cos(y)]' = sin(y) + (y-ya) cos(y) - sin(y) = (y-ya) cos(y)
  !  
  !  Area in m2, weight with cos(y) to account for smaller cells towards poles:
  !
  !         y2                                                  y2
  !    A  =  | a (y-ya) cos(y) dy = a [ (y-ya) sin(y) + cos(y) ]
  !         y=y1                                                y=y1
  !
  !  For upper triangle:
  !                                                      ymax
  !   A  =  Dx/(ymid-ymax) * [ (y-ymax) sin(y) + cos(y) ]  
  !                                                      y=ymid
  !
  !      =  Dx/(ymid-ymax) * ( { cos(ymax) } - { (ymid-ymax) sin(ymid) + cos(ymid) } )
  !
  !      =  Dx/(ymid-ymax) * ( (ymax-ymid) sin(ymid) + cos(ymax) - cos(ymid) )
  !
  !  For lower triangle:
  !                                                      ymid
  !   A  =  Dx/(ymid-ymin) * [ (y-ymin) sin(y) + cos(y) ]
  !                                                      y=ymin
  !
  !      =  Dx/(ymid-ymin) * ( { (ymid-ymin) sin(ymid) + cos(ymid) } - { cos{ymin} } )
  !
  !      =  Dx/(ymid-ymin) * ( (ymid-ymin) sin(ymid) + cos(ymid) - cos{ymin} )
  !

  subroutine Polygon_LonLat_TriangleArea( self, area, status )!, debug )
  
    use Binas, only : deg2rad, ae

    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    real(wpr), intent(out)                    ::  area
    integer, intent(out)                      ::  status
    
    !logical, intent(in), optional             ::  debug
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_LonLat_TriangleArea'
    
    ! --- local ----------------------------------
    
    !logical           ::  dbg
    integer           ::  imin, imax, imid
    real(wpr)         ::  xmin, xmax, xmid
    real(wpr)         ::  ymin, ymax, ymid
    integer           ::  i
    real(wpr)         ::  xmin_rad, xmid_rad, xmax_rad
    real(wpr)         ::  ymin_rad, ymid_rad, ymax_rad
    type(T_Vector)    ::  vmin, vmax, vmid
    type(T_Vector)    ::  vmid2
    type(T_Vector)    ::  p
    real(wpr)         ::  xp_rad
    real(wpr)         ::  Dx

    ! --- begin ----------------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! check ...
    if ( .not. self%IsTriangle() ) then
      write (gol,'("input polygon is not a triangle")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! find bottom and top vertices:
    imin = 1
    imax = 1
    ymin = self%vertex(imin)%y
    ymax = self%vertex(imax)%y
    do i = 2, self%n
      if ( self%vertex(i)%y < ymin ) then
        imin = i
        ymin = self%vertex(imin)%y
      end if
      if ( self%vertex(i)%y > ymax ) then
        imax = i
        ymax = self%vertex(imax)%y
      end if
    end do
    
    ! only non-empty results if polygon is not a line ...
    if ( ymin < ymax ) then
    
      ! find mid point:
      imid = 1
      do i = 1, 3
        if ( i == imin ) cycle
        if ( i == imax ) cycle
        imid = i
        exit
      end do
      ! mid value:
      ymid = self%vertex(imid)%y
      
      ! x-coordinates:
      xmin = self%vertex(imin)%x
      xmid = self%vertex(imid)%x
      xmax = self%vertex(imax)%x

      ! coordinates in rad:
      xmin_rad = xmin * deg2rad
      xmid_rad = xmid * deg2rad
      xmax_rad = xmax * deg2rad
      ymin_rad = ymin * deg2rad
      ymid_rad = ymid * deg2rad
      ymax_rad = ymax * deg2rad
      
      !! testing ...
      !if (dbg) print *, '    ttt1 ', ymin, ymid, ymax

      ! only lower triangle ?
      if ( ymid == ymax ) then
        !
        !  ymax o---o ymid o---o vmax
        !        \ /        \ /
        !         o   ymin   o
        !
        !  A  =  Dx/(ymid-ymin) * ( (ymid-ymin) sin(ymid) + cos(ymid) - cos{ymin} )
        !
        ! width of triangle base:
        Dx = abs( xmid_rad - xmax_rad )
        ! area:
        area = Dx/(ymid_rad-ymin_rad) * ( (ymid_rad-ymin_rad) * sin(ymid_rad) + cos(ymid_rad) - cos(ymin_rad) )
      
        !! testing ...
        !if (dbg) print *, '    ttt1 lower ', area
        
      ! only upper triangle ?
      else if ( ymid == ymin ) then
        !
        !         o   ymax   o
        !        / \        / \
        !  ymin o---o ymid o---o ymin
        !
        !  A  =  Dx/(ymid-ymax) * ( (ymax-ymid) sin(ymid) + cos(ymax) - cos(ymid) )
        !
        ! width of triangle base:
        Dx = abs( xmid_rad - xmin_rad )
        ! area:
        area = Dx/(ymid_rad-ymax_rad) * ( (ymax_rad-ymid_rad) * sin(ymid_rad) + cos(ymax_rad) - cos(ymid_rad) )
      
        !! testing ...
        !if (dbg) print *, '    ttt1 upper ', area
      
      else
      
        !
        !         o ymax o
        !        /|      |\
        !  ymid o-* p  p *-o ymid
        !        \|      |/
        !         o ymin o
        !
        
        ! define vertices as vectors:
        call vmin%Init( xmin, ymin, status )
        IF_NOT_OK_RETURN(status=1)
        call vmid%Init( xmid, ymid, status )
        IF_NOT_OK_RETURN(status=1)
        call vmax%Init( xmax, ymax, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! just another point right of vmid ...
        call vmid2%Init( xmid+100.0, ymid, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! intersection point between edge [vmin,vmax]
        ! and horizontal line through vmid:
        p = Intersection_Segment_Line( vmin, vmax, vmid, vmid2 )
        
        ! x-coordinate in radians:
        xp_rad = p%x * deg2rad
      
        !
        !  A  =  Dx/(ymid-ymin) * ( (ymid-ymin) sin(ymid) + cos(ymid) - cos{ymin} )
        !
        ! width of triangle base:
        Dx = abs( xp_rad - xmid_rad )
        ! area:
        area = Dx/(ymid_rad-ymin_rad) * ( (ymid_rad-ymin_rad) * sin(ymid_rad) + cos(ymid_rad) - cos(ymin_rad) )
      
        !! testing ...
        !if (dbg) print *, '    ttt1 lower ', area
        
        !
        !  A  =  Dx/(ymid-ymax) * ( (ymax-ymid) sin(ymid) + cos(ymax) - cos(ymid) )
        !
        ! width of triangle base:
        Dx = abs( xp_rad - xmid_rad )
        ! area:
        area = area + Dx/(ymid_rad-ymax_rad) * ( (ymax_rad-ymid_rad) * sin(ymid_rad) + cos(ymax_rad) - cos(ymid_rad) )
      
        !! testing ...
        !if (dbg) print *, '    ttt1 +upper ', area
      
      end if
      
      ! trap tiny negatives, result of subtracting sin and cos:
      if ( (-1.0e-6 < area) .and. (area < 0.0) ) area = 0.0
      
      ! convert from rad2 to m2:
      area = area * ae**2
      
    else
    
      ! empty:
      area = 0.0
    
    end if ! non-line
    
    ! ok
    status = 0
    
  end subroutine Polygon_LonLat_TriangleArea



  ! ***
  
  
  !
  ! Compute polygon area in m2 if coordinates represent (lon,lat).
  ! Algorithm devides the polygon in triangles and sums their areas.
  !
  
  subroutine Polygon_LonLat_Area( self, area, status )!, debug )
  
    !use GO_Plot, only : T_PlotFile

    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)              ::  self
    real(wpr), intent(out)                    ::  area   ! m2
    integer, intent(out)                      ::  status
    
    !logical, intent(in), optional             ::  debug
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_LonLat_Area'
    
    ! --- local ----------------------------------
    
    !logical           ::  dbg
    type(T_Polygon)   ::  wp
    type(T_Polygon)   ::  tri
    real(wpr)         ::  tri_area
    integer           ::  i0
    logical           ::  found
    
    !integer           ::  k
    !type(T_PlotFile)  ::  pf

    ! --- begin ----------------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! init work polygon as copy:
    call wp%Init( self%n, status )
    IF_NOT_OK_RETURN(status=1)
    ! init as copy:
    call wp%Copy( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! init empty triangle:
    call tri%Init( 3, status )
    IF_NOT_OK_RETURN(status=1)
 
    ! init result:
    area = 0.0
    ! loop unit all found defined:
    do

      !! testing ..
      !if (dbg) then
      !  print *, '  aaa1 current work polygon '
      !  do k = 1, wp%n
      !    print *, '    a1 vertex ', k, ' = ', wp%vertex(k)%x, wp%vertex(k)%y
      !  end do
      !  print *, '    a1 is triangle: ', wp%IsTriangle()
      !end if

      ! end if work polygon is trangle:
      if ( wp%IsTriangle() ) then
        ! compute area:
        call wp%LonLat_TriangleArea( tri_area, status )!, debug=debug )
        IF_NOT_OK_RETURN(status=1)
        ! add to sum:
        area = area + tri_area
        !! testing ..
        !if (dbg) then
        !  print *, '  aaa1 add final ', tri_area, area
        !  do k = 1, wp%n
        !    print *, '    a1 vertex ', k, ' = ', wp%vertex(k)%x, wp%vertex(k)%y
        !  end do
        !end if
        ! leave:
        exit
      end if  ! wp is triangle

      ! init flag:
      found = .false.
      ! loop over vertices that serve as tip of triangle:
      do i0 = 1, wp%n

        !! testing ..
        !if (dbg) then
        !  print *, '  aaa1 try to form triangle with tip ', i0
        !end if
        
        ! try to remove triangle, result could be empty:
        call wp%RemoveTriangle( i0, tri, status )
        IF_NOT_OK_RETURN(status=1)

        ! result filled?
        found = .not. tri%IsEmpty()
        ! found?
        if ( found ) exit

      end do  ! vertices i0 as tip of triangle

      ! new triangle found?
      if ( found ) then

        ! compute area:
        call tri%LonLat_TriangleArea( tri_area, status )!, debug=debug )
        IF_NOT_OK_RETURN(status=1)
        ! add to sum:
        area = area + tri_area

        !! testing ..
        !if (dbg) then
        !  print *, '  aaa1 add tria  ', tri_area, area
        !  do k = 1, tri%n
        !    print *, '    a1 vertex ', k, ' = ', tri%vertex(k)%x, tri%vertex(k)%y
        !  end do
        !end if
        
      else

        !! testing ...
        !if (dbg) then
        !  write (gol,'("could not define new triangle")'); call goErr
        !  ! ... plot ...
        !  call pf%Init( 'debug-lonlat-area.py', status )
        !  IF_NOT_OK_RETURN(status=1)
        !  call wp%Plot( 'edges', pf, status, kwargs='color="red"' )
        !  IF_NOT_OK_RETURN(status=1)
        !  call pf%Done( status )
        !  IF_NOT_OK_RETURN(status=1)
        !  ! ..............
        !  TRACEBACK; status=1; return
        !end if
        
        ! probably workpolygon is just a line now; leave:
        exit

      end if
      
    end do ! until wp is triangle or has zero area
    
    ! clear:
    call wp%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call tri%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Polygon_LonLat_Area



  ! ***
  
  
  !
  ! Compute fraction of input polygon 'pg' that is covered by self
  ! if coordinates represent (lon,lat).
  !
  
  subroutine Polygon_LonLat_AreaFraction( self, pg, frac, status )!, debug )
  
    use GO_Plot, only : T_PlotFile

    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)              ::  self
    class(T_Polygon), intent(in)              ::  pg
    real, intent(out)                         ::  frac   ! [0,1]
    integer, intent(out)                      ::  status
    
    !logical, intent(in), optional             ::  debug
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_LonLat_AreaFraction'
    
    ! --- local ----------------------------------
    
    !logical               ::  dbg
    real(wpr)             ::  pg_area
    type(T_Polygon)       ::  wp
    real(wpr)             ::  wp_area
    type(T_PlotFile)      ::  pf

    ! --- begin ----------------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug

    ! init work polygon as copy of self clipped by input polygon:
    call wp%InitClippedPolygon( self, pg, status )
    IF_NOT_OK_RETURN(status=1)
    ! no overlap?
    if ( wp%IsEmpty() ) then
      ! no overlap fraction:
      frac = 0.0
    else
      ! intersection area:
      call wp%LonLat_Area( wp_area, status )
      IF_NOT_OK_RETURN(status=1)
      ! source area:
      call pg%LonLat_Area( pg_area, status )
      IF_NOT_OK_RETURN(status=1)
      ! fraction of source
      frac = real( wp_area / pg_area )
      ! trap strange fractions ...
      if ( (frac < 0.0) .or. (1.0 < frac) ) then
        ! round if just above 1.0 ...
        if ( abs(frac-1.0) < 1.0e-6 ) then
          frac = 1.0
        else
          write (gol,'("strange area fraction  : ",es12.4)') frac; call goErr
          write (gol,'("  intersection area    : ",es12.4)') wp_area; call goErr
          write (gol,'("  source area          : ",es12.4)') pg_area; call goErr
          ! ... plot ...
          call pf%Init( 'debug.py', status )
          IF_NOT_OK_RETURN(status=1)
          call self%Plot( 'fill', pf, status, kwargs='color="red", alpha=0.1' )
          IF_NOT_OK_RETURN(status=1)
          call self%Plot( 'edges', pf, status, kwargs='color="red"' )
          IF_NOT_OK_RETURN(status=1)
          call pg%Plot( 'fill', pf, status, kwargs='color="blue", alpha=0.1' )
          IF_NOT_OK_RETURN(status=1)
          call pg%Plot( 'edges', pf, status, kwargs='color="blue"' )
          IF_NOT_OK_RETURN(status=1)
          call wp%Plot( 'fill', pf, status, kwargs='color="green", alpha=0.1' )
          IF_NOT_OK_RETURN(status=1)
          call wp%Plot( 'edges', pf, status, kwargs='color="green"' )
          IF_NOT_OK_RETURN(status=1)
          call pf%Done( status )
          IF_NOT_OK_RETURN(status=1)
          ! ..............
          TRACEBACK; status=1; return
        end if  ! no roundoff
      end if  ! strange fraction
    end if  ! overlap
    
    ! clear:
    call wp%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Polygon_LonLat_AreaFraction


  ! ***
  
  
  !
  ! Compute archlengths in m for polygon sides
  ! if coordinates represent (lon,lat).
  !
  
  subroutine Polygon_LonLat_ArcLengths( self, lengths, status )
  
    use Binas, only : deg2rad
    use Binas, only : ae   ! m
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)              ::  self
    real, intent(out)                         ::  lengths(:)
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_LonLat_Area'
    
    ! --- local ----------------------------------
    
    integer      ::  i
    real(8)      ::  lam0, phi0, lam1, phi1
    real(8)      ::  v0(3), v1(3)
    real(8)      ::  cos_alpha
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(lengths) < self%n ) then
      write (gol,'("storage for arc lengths is ",i0," while polygon as ",i0," sides")') &
              size(lengths), self%n; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop over edges:
    do i = 1, self%n
    
      ! start location in polar coordinates:
      lam0 = self%vertex(i)%x * deg2rad  ! rad
      phi0 = self%vertex(i)%y * deg2rad  ! rad
      ! vector position:
      !  x-axis from origin to (lon= 0,lat= 0)
      !  y-axis from origin to (lon=90,lat= 0)
      !  z-axis from origin to (lon= 0,lat=90)
      v0(1) = cos(lam0) * cos(phi0)
      v0(2) = sin(lam0) * cos(phi0)
      v0(3) =             sin(phi0)

      ! end location in polar coordinates:
      if ( i < self%n ) then
        lam1 = self%vertex(i+1)%x * deg2rad  ! rad
        phi1 = self%vertex(i+1)%y * deg2rad  ! rad
      else
        lam1 = self%vertex(  1)%x * deg2rad  ! rad
        phi1 = self%vertex(  1)%y * deg2rad  ! rad
      end if
      ! vector position:
      !  x-axis from origin to (lon= 0,lat= 0)
      !  y-axis from origin to (lon=90,lat= 0)
      !  z-axis from origin to (lon= 0,lat=90)
      v1(1) = cos(lam1) * cos(phi1)
      v1(2) = sin(lam1) * cos(phi1)
      v1(3) =             sin(phi1)

      ! angle between two positions:
      !   cos(alpha) = v1.v2 / (|v1||v2|)
      ! here v1 and v2 have both unit length:
      cos_alpha = dot_product( v0, v1 )
      ! truncate to avoid rounding errors:
      cos_alpha = min( max( -1.0, cos_alpha ), 1.0 )

      ! distance over globe:
      lengths(i) = real( ae * acos(cos_alpha) )  ! m

    end do
    
    ! ok
    status = 0

  end subroutine Polygon_LonLat_ArcLengths

  ! ***
  
  
  !
  ! Initialize a polygon as copy of the input polygon 
  ! clipped by the line (y1,y2).
  ! Everything right of the line is removed.
  !
  
  subroutine Polygon_CopyClippedByLine( self, poly, y1, y2, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(inout)           ::  self
    class(T_Polygon), intent(in)              ::  poly
    type(T_Vector), intent(in)                ::  y1, y2
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_CopyClippedByLine'
    
    ! --- local ----------------------------------
    
    type(T_Vector)       ::  p
    integer              ::  i1, i2
    logical              ::  inside1, inside2

    ! --- begin ----------------------------------
 
    ! reset counter for the output polygon:
    self%n = 0
 
    ! loop over edges:
    do i1 = 1, poly%n
      ! end vertex of edge:
      i2 = modulo(i1,poly%n) + 1

      ! test if vertices are inside half-plane left from (y1,y2):
      call poly%vertex(i1)%LeftFromLine( y1, y2, inside1, status )!, debug=.true. )
      IF_NOT_OK_RETURN(status=1)
      call poly%vertex(i2)%LeftFromLine( y1, y2, inside2, status )!, debug=.true. )
      IF_NOT_OK_RETURN(status=1)
      
      ! vertex 1 in inside clipping region ?
      if ( inside1 ) then

        ! vertex 2 in inside clipping region?
        if ( inside2 ) then
        
          ! add the vertex 2 to the output polygon:
          call self%AddVertex( poly%vertex(i2), status )
          IF_NOT_OK_RETURN(status=1)
          
        else
        
          ! vertex i+1 is outside
          p = Intersection_Segment_Line( poly%vertex(i1), poly%vertex(i2), y1, y2 )
          ! add intersection point as new vertex:
          call self%AddVertex( p, status )
          IF_NOT_OK_RETURN(status=1)
          
        end if
        
      else
      
        ! vertex i is outside
        if ( inside2 ) then
        
          ! intersection between edge and line:
          p = Intersection_Segment_Line( poly%vertex(i1), poly%vertex(i2), y1, y2 )
          ! defined?
          if ( .not. p%IsEmpty() ) then
            ! add intersection point as new vertex:
            call self%AddVertex( p, status )
            IF_NOT_OK_RETURN(status=1)
          end if
          
          ! add end point new vertex:
          call self%AddVertex( poly%vertex(i2), status )
          IF_NOT_OK_RETURN(status=1)
          
        end if ! second vertex of edge inside halfplane
        
      end if ! first vertex of edge inside halfplane

    end do  ! edges of poly
    
    ! ok
    status = 0

  end subroutine Polygon_CopyClippedByLine
  
  
  ! ***
  
  
  !
  ! Return .true. if vector is inside polygon,
  ! where the polygon is assumed to be convex.
  !

  subroutine Polygon_InsideConvex( self, v, inside, status )!, debug )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)              ::  self
    type(T_Vector), intent(in)                ::  v
    logical, intent(out)                      ::  inside
    integer, intent(out)                      ::  status
    
    !logical, intent(in), optional             ::  debug
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_InsideConvex'
    
    ! --- local ----------------------------------
    
    !logical           ::  dbg
    integer           ::  i1, i2
    
    ! --- begin ----------------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! loop over edges:
    do i1 = 1, self%n
      ! end point of edge:
      i2 = modulo(i1,self%n) + 1
      ! inside half plane defined by edge?
      call v%LeftFromLine( self%vertex(i1), self%vertex(i2), inside, status )!, debug=debug )
      IF_NOT_OK_RETURN(status=1)
      !! testing ...
      !if ( dbg) then
      !  print *, '      insideconvex p ', v%x, v%y, inside
      !  print *, '                   v1', self%vertex(i1)%x, self%vertex(i1)%y
      !  print *, '                   v2', self%vertex(i2)%x, self%vertex(i2)%y
      !end if
      ! leave?
      if ( .not. inside ) exit
    end do
    
    ! ok
    status = 0

  end subroutine Polygon_InsideConvex


  ! ***
  
  
  !
  ! Return bounding box [x0,x1,y0,y1].
  !
  
  subroutine Polygon_BoundBox( self, box, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)        ::  self
    real(wpr), intent(out)              ::  box(4)
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_BoundBox'
    
    ! --- local ----------------------------------
    
    integer     ::  i

    ! --- begin ----------------------------------
    
    ! init from first:
    box(1) = self%vertex(1)%x
    box(2) = self%vertex(1)%x
    box(3) = self%vertex(1)%y
    box(4) = self%vertex(1)%y
    ! loop over other vertices:
    do i = 2, self%n
      ! update:
      box(1) = min( box(1), self%vertex(i)%x )
      box(2) = max( box(2), self%vertex(i)%x )
      box(3) = min( box(3), self%vertex(i)%y )
      box(4) = max( box(4), self%vertex(i)%y )
    end do

    ! ok
    status = 0
  
  end subroutine Polygon_BoundBox


  ! ***
  
  
  !
  ! Return diameter of bounding box:
  !    [x0,x1,y0,y1]
  !
  
  subroutine Polygon_BoundBoxDiameter( self, L, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)        ::  self
    real(wpr), intent(out)              ::  L
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_BoundBoxDiameter'
    
    ! --- local ----------------------------------
    
    real(wpr)   ::  box(4)

    ! --- begin ----------------------------------
    
    ! get bounding box:
    call self%BoundBox( box, status )
    IF_NOT_OK_RETURN(status=1)
      
    ! diagonal length:
    L = sqrt( (box(2)-box(1))**2 + (box(4)-box(3))**2 )

    ! ok
    status = 0
  
  end subroutine Polygon_BoundBoxDiameter


  ! ***
  
  
  !
  ! Return list of triangular polgons covering input.
  !
  
  subroutine Polygon_Triangles( self, ntri, tris, status )
  
    !use GO_Plot, only : T_PlotFile

    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)              ::  self
    integer, intent(out)                      ::  ntri
    class(T_Polygon), intent(inout)           ::  tris(:)
    integer, intent(out)                      ::  status
    
    !type(T_PlotFile), intent(inout)           ::  pf
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_Triangles'
    
    ! --- local ----------------------------------
    
    !logical           ::  dbg
    integer           ::  i0
    type(T_Polygon)   ::  wp
    logical           ::  found

    ! --- begin ----------------------------------
    
    ! init work polygon as copy:
    call wp%Init( self%n, status )
    IF_NOT_OK_RETURN(status=1)
    ! init as copy:
    call wp%Copy( self, status )
    IF_NOT_OK_RETURN(status=1)
 
    ! init counter:
    ntri = 0
    ! loop unit all triangles defined:
    do
      ! increase counter:
      ntri = ntri + 1
      ! check ..
      if ( ntri > size(tris) ) then
        write (gol,'("need storage for more then ",i0," triangles")') size(tris); call goErr
        TRACEBACK; status=1; return
      end if

      !! testing ...
      !if (dbg) call wp%Plot( 'fill', pf, status, kwargs='color="green", alpha=0.1' )

      ! end if work polygon is trangle:
      if ( wp%IsTriangle() ) then
        ! add copy:
        call tris(ntri)%Copy( wp, status )
        IF_NOT_OK_RETURN(status=1)
        ! leave:
        exit
      end if

      ! init flag:
      found = .false.
      ! loop over vertices that serve as tip of triangle:
      do i0 = 1, wp%n
      
        ! try to remove triangle, result could be empty:
        call wp%RemoveTriangle( i0, tris(ntri), status )
        IF_NOT_OK_RETURN(status=1)
        
        ! result filled?
        found = .not. tris(ntri)%IsEmpty()
        ! found?
        if ( found ) exit
      
      end do  ! vertices i0 as tip of triangle
      
      ! check ...
      if ( .not. found ) then
        write (gol,'("could not define triange ",i0)') ntri; call goErr
        TRACEBACK; status=1; return
      end if

    end do ! until wp is triangle
    
    ! clear:
    call wp%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Polygon_Triangles


  ! ***
  
  
  !
  ! Fill interpolation weights based on distance to corners:
  !
  !     di  = distance to corner i
  !
  !            1/di
  !     wi = --------
  !          sum 1/dj
  !
  ! In 1D this is equivalent to linear interpolation.
  ! For example, target point p at 1/4 of interval [x1,x2] with length 12:
  !
  !           * _
  !               o _
  !                   - _
  !                       - _
  !                           *
  !   --------+---|-----------+----------------
  !           x1  p           x2
  !
  !    w1 = 1/3 / (1/3 + 1/9) = 1/3 / 4/9 = 3/4
  !    w2 = 1/9 / (1/3 + 1/9) = 1/9 / 4/9 = 1/4
  !
  ! Return status:
  !   -1  : p not inside 
  !
  
  subroutine Polygon_CornerWeights( self, p, w, status )!, debug )
  
    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)        ::  self
    class(T_Vector), intent(in)         ::  p
    real(wpr), intent(out)              ::  w(:)  ! (n)
    integer, intent(out)                ::  status

    !logical, intent(in), optional       ::  debug

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_Inside'
    
    ! --- local ----------------------------------
    
    !logical               ::  dbg
    logical               ::  inside
    type(T_Vector)        ::  line
    real(wpr)             ::  d
    integer               ::  i
    
    ! --- begin ----------------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! inside convex around vertices?
    call self%InsideConvex( p, inside, status )!, debug=debug )
    IF_NOT_OK_RETURN(status=1)
    ! check ...
    if ( .not. inside ) then
      ! leave with warning status:
      status=-1; return
    end if
    
    ! check ...
    if ( size(w) /= self%n ) then
      write (gol,'("output array w (",i0,") should have size equal to number of corners ",i0)') &
                       size(w), self%n; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop over corners:
    do i = 1, self%n
      ! vector from p to corner:
      line = self%vertex(i) - p
      ! distance:
      d = Length(line)
      ! zero length is exactly on corner:
      if ( d == 0.0 ) then
        ! reset all weights to zero:
        w = 0.0
        ! only this one needed:
        w(i) = 1.0
        ! leave:
        exit
      else
        ! distance weighted:
        w(i) = 1.0 / d
      end if
    end do
    ! normalize:
    w = w / sum(w)
    
    ! ok
    status = 0
    
  end subroutine Polygon_CornerWeights


  ! ***
  
  
  !
  ! Add python matplotlib lines to text file with unit fu.
  ! Command defines what is plotted from the polygon:
  !  'fill'         : filled area
  !  'center'       : mark center location
  !  'edges'        : lines at edges
  !  'edge-arrows'  : edge arrows
  !  'normals'      : normal vectors at edges
  !
  
  subroutine Polygon_Plot( self, command, pf, status, kwargs )
  
    use GO_Plot, only : T_PlotFile

    ! --- in/out ---------------------------------
    
    class(T_Polygon), intent(in)              ::  self
    character(len=*), intent(in)              ::  command
    type(T_PLotFile), intent(in)              ::  pf
    integer, intent(out)                      ::  status
    
    character(len=*), intent(in), optional  ::  kwargs

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Polygon_Plot'
    
    ! --- local ----------------------------------
    
    character(len=1024)     ::  line
    integer                 ::  i, i1, i2
    character(len=1024)     ::  xx, yy
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(command) )
    
      ! area fill:
      case ( 'fill' )
    
        ! collect corners:
        xx = '['
        yy = '['
        do i = 1, self%n
          if ( i > 1 ) then
            xx = trim(xx)//','
            yy = trim(yy)//','
          end if
          write (xx,'(a,e16.8)') trim(xx), self%vertex(i)%x
          write (yy,'(a,e16.8)') trim(yy), self%vertex(i)%y
        end do
        xx = trim(xx)//']'
        yy = trim(yy)//']'

        ! plot command:
        write (pf%fu,'(a)') '# add filled polygon:'
        write (line,'("ax.fill( ",a,", ",a)') trim(xx), trim(yy)
        if ( present(kwargs) ) line = trim(line)//', '//trim(kwargs)
        line = trim(line)//' )'
        write (pf%fu,'(a)') trim(line)
    
      ! edge lines:
      case ( 'edges' )
    
        ! loop over edges:
        do i1 = 1, self%n
          ! end point:
          i2 = modulo(i1,self%n) + 1
          ! add plot commands:
          call self%vertex(i2)%Plot( 'line', pf, status, kwargs=kwargs, origin=self%vertex(i1) )
          IF_NOT_OK_RETURN(status=1)
        end do
    
      ! edge arrows:
      case ( 'edge-arrows' )
    
        ! loop over edges:
        do i1 = 1, self%n
          ! end point:
          i2 = modulo(i1,self%n) + 1
          ! add plot commands:
          call self%vertex(i2)%Plot( 'arrow', pf, status, kwargs=kwargs, origin=self%vertex(i1) )
          IF_NOT_OK_RETURN(status=1)
        end do
    
      ! unknown ...
      case default
        write (gol,'("unsupported command: ",a)') trim(command); call goErr
        TRACEBACK; status=1; return
        
    end select

    ! ok
    status = 0
  
  end subroutine Polygon_Plot
 

end module GO_Vector




!#######################################################################
!###
!### end
!###
!#######################################################################
  
