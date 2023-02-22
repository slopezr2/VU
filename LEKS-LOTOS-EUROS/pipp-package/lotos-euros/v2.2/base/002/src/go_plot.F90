!###############################################################################
!
! NAME
!   GO_Plot - general objects to write python plotting file
!
! USAGE
!
!   type(T_PlotFile)    ::  pf
!
!   ! init new figure:
!   call pf%Init( 'debug.py', status )
!   IF_NOTOK_RETURN(status=1)
!
!   ! add lines for new figure:
!   call pf%NewFigure( status )
!   IF_NOTOK_RETURN(status=1)
!
!   ! add line:
!   call pf%Add( '# add line:' )
!   IF_NOTOK_RETURN(status=1)
!   call pf%Add( 'ax.plot( [0,1], [1,2], color="red" )' )
!   IF_NOTOK_RETURN(status=1)
!   call pf%Add( '' )
!   IF_NOTOK_RETURN(status=1)
!
!   ! write and close:
!   call pf%Done( status )
!   IF_NOTOK_RETURN(status=1)
!
! EXAMPLE FILE
!
!   # modules:
!   import matplotlib.pyplot as plt
!   
!   # new figure:
!   fig = plt.figure()'
!   ax = fig.add_axes([0.1,0.1,0.8,0.8])
!   #ax.set_aspect("equal")
!
!   # add line:
!   ax.plot( [0,1], [1,2], color="red" )
!   
!   # show:
!   plt.show()
!   
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#################################################################

module GO_Plot

  use GO_Print, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------------

  private

  public  ::  T_PlotFile
  

  ! --- const ------------------------------------
  
  character(len=*), parameter  ::  mname = 'GO_Plot'
  
  
  ! --- types ------------------------------------
  
  ! tools
  type T_PlotFile
    ! file unit:
    integer               ::  fu
    ! target filename:
    character(len=1024)   ::  filename
    ! temporary content:
    character(len=1024)   ::  line
    !
  contains
    procedure   ::  Init                     => PlotFile_Init
    procedure   ::  Done                     => PlotFile_Done
    procedure   ::  NewFigure                => PlotFile_NewFigure
    procedure   ::  NextColor                => PlotFile_NextColor
    procedure   ::  Add                      => PlotFile_Add
    procedure   ::  Line_Write               => PlotFile_Line_Write
    procedure   ::  Line_AddKwargs           => PlotFile_Line_AddKwargs
    procedure   ::  Line_Add                 => PlotFile_Line_Add
    procedure   ::  Line_AddList             => PlotFile_Line_AddList
    procedure   ::  AddLines                 => PlotFile_AddLines
    procedure   ::  AddPolygon               => PlotFile_AddPolygon
    procedure   ::  AddRectangle             => PlotFile_AddRectangle
    procedure   ::  AddArrow                 => PlotFile_AddArrow
  end type T_PlotFile


contains


  ! ********************************************************************
  ! ***
  ! *** plot tools
  ! ***
  ! ********************************************************************
  
  subroutine PlotFile_Init( self, filename, status )
  
    use GO_FU, only : goGetFU
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(out)            ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
  
    ! store:
    self%filename = trim(filename)
    
    ! new file unit:
    call goGetFU( self%fu, status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    write (gol,'("open file for plotting: ",a)') trim(filename); call goPr

    ! open file:
    open( self%fu, file=trim(filename), form='formatted', status='replace', iostat=status )
    IF_NOTOK_RETURN(status=1)
  
    ! header:
    write (self%fu,'(a)') '# modules:'
    write (self%fu,'(a)') 'import matplotlib.pyplot as plt'
    write (self%fu,'(a)') ''
    write (self%fu,'(a)') '# list of standard colors:'
    write (self%fu,'(a)') 'colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]'
    write (self%fu,'(a)') '# init color index:'
    write (self%fu,'(a)') 'icolor = 0'
    write (self%fu,'(a)') ''

    ! ok:
    status = 0

  end subroutine PlotFile_Init
  
  ! *
  
  subroutine PlotFile_NewFigure( self, status, kwargs )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(in)           ::  self
    integer, intent(out)                    ::  status
    character(len=*), intent(in), optional  ::  kwargs

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_NewFigure'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! add:
    write (self%fu,'(a)') '# new figure:'
    if ( present(kwargs) ) then
      write (self%fu,'(a)') 'fig = plt.figure('//trim(kwargs)//')'
    else
      write (self%fu,'(a)') 'fig = plt.figure()'
    end if
    write (self%fu,'(a)') 'ax = fig.add_axes([0.1,0.1,0.8,0.8])'
    write (self%fu,'(a)') '#ax.set_aspect("equal")'
    write (self%fu,'(a)') ''
    
    ! ok
    status = 0

  end subroutine PlotFile_NewFigure
  
  ! *
  
  subroutine PlotFile_Add( self, line, status, &
                            comment )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(in)             ::  self
    character(len=*), intent(in)              ::  line
    integer, intent(out)                      ::  status
    
    character(len=*), intent(in), optional    ::  comment

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_Add'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! add comment?
    if ( present(comment) ) then
      write (self%fu,'(a)') trim(comment)
    end if
    ! add new line:
    write (self%fu,'(a)') trim(line)
    
    ! ok
    status = 0

  end subroutine PlotFile_Add

  
  ! *
  
  
  !
  ! Add line to plot script that defines next color:
  !   color = colors[icolor] ; ..
  !
  
  subroutine PlotFile_NextColor( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(in)           ::  self
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_NextColor'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! extract current color, increase index:
    call self%Add( 'color = colors[icolor]; icolor = (icolor+1) % len(colors)', status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine PlotFile_NextColor
  

  ! *
  
  
  !
  ! Write current "line" attribute.
  !
  
  subroutine PlotFile_Line_Write( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)          ::  self
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_Line__Write'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! add:
    write (self%fu,'(a)') trim(self%line)
    
    ! ok
    status = 0

  end subroutine PlotFile_Line_Write


  ! *

  
  ! Add text to internal line buffer.
  ! If clear=.true, first make the buffer empty.
  
  subroutine PlotFile_Line_Add( self, text, status, clear )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)        ::  self
    character(len=*), intent(in)            ::  text
    integer, intent(out)                    ::  status
    logical, intent(in), optional           ::  clear

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_Line_Add'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! new?
    if ( present(clear) ) then
      if ( clear ) self%line = ''
    end if
    ! add content:
    self%line = trim(self%line)//trim(adjustl(text))
        
    ! ok
    status = 0
  
  end subroutine PlotFile_Line_Add


  ! *

  
  !
  ! Add list to self%line filled from values.
  ! If "closed=.true." the first element is added again to have a closed line.
  !
  
  subroutine PlotFile_Line_AddList( self, x, status, &
                                      closed )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)        ::  self
    real, intent(in)                        ::  x(:)
    integer, intent(out)                    ::  status
    
    logical, intent(in), optional           ::  closed

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_Line_AddList'
    
    ! --- local ----------------------------------
    
    character(len=32)       ::  s
    integer                 ::  i
    
    ! --- begin ----------------------------------
    
    ! init list:
    write (s,*) '[', x(1)
    self%line = trim(self%line)//trim(adjustl(s))
    ! add other elements:
    do i = 2, size(x)
      write (s,*) ',', x(i)
      self%line = trim(self%line)//trim(adjustl(s))
    end do
    ! add start again?
    if ( present(closed) ) then
      if ( closed ) then
        write (s,*) ',', x(1)
        self%line = trim(self%line)//trim(adjustl(s))
      end if
    end if
    ! close:
    self%line = trim(self%line)//']'
        
    ! ok
    status = 0
  
  end subroutine PlotFile_Line_AddList


  ! *

  
  ! Add optional elements
  
  subroutine PlotFile_Line_AddKwargs( self, status, &
                                       kwargs, alpha )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)        ::  self
    integer, intent(out)                    ::  status
    
    character(len=*), intent(in), optional  ::  kwargs
    real, intent(in), optional              ::  alpha  

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_Line_AddKwargs'
    
    ! --- local ----------------------------------
    
    character(len=32)     ::  arg
    
    ! --- begin ----------------------------------
    
    ! extra keyword arguments?
    if ( present(kwargs) ) self%line = trim(self%line)//', '//trim(kwargs)
    
    ! extra alpha:
    if ( present(alpha ) ) then
      write (arg,'("alpha=",f4.2)') alpha
      self%line = trim(self%line)//', '//trim(arg)
    end if
    
    ! ok
    status = 0
  
  end subroutine PlotFile_Line_AddKwargs


  ! *

  
  subroutine PlotFile_AddLines( self, x, y, status, &
                                  kwargs, xmap, closed )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)        ::  self
    real, intent(in)                        ::  x(:), y(:)
    integer, intent(out)                    ::  status
    
    real, intent(in), optional              ::  xmap(2)   ! (scale,center)
    character(len=*), intent(in), optional  ::  kwargs
    logical, intent(in), optional           ::  closed

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_AddLines'
    
    ! --- local ----------------------------------
    
    character(len=1024)     ::  line
    character(len=1024)     ::  xx, yy
    integer                 ::  i
    real                    ::  xscale, xcenter
    
    ! --- begin ----------------------------------
    
    ! scale in x-direction around center?
    xscale = 1.0
    xcenter = 0.0
    if ( present(xmap) ) then
      xscale  = xmap(1)
      xcenter = xmap(2)
    endif

    ! check ..
    if ( size(x) /= size(y) ) then
      write (gol,'("x (",i0,") and y (",i0,") should have same lenght")') size(x), size(y); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init line:
    call self%Line_Add( 'ax.plot( ', status, clear=.true. )
    IF_NOTOK_RETURN(status=1)
    ! add x list:
    call self%Line_AddList( xcenter+(x-xcenter)*xscale, status, closed=closed )
    IF_NOTOK_RETURN(status=1)
    ! sep:
    self%line = trim(self%line)//','
    ! add y list:
    call self%Line_AddList( y, status, closed=closed )
    IF_NOTOK_RETURN(status=1)
    ! extra arguments?
    call self%Line_AddKwargs( status, kwargs=kwargs )
    IF_NOTOK_RETURN(status=1)
    ! close:
    call self%Line_Add( ')', status )
    IF_NOTOK_RETURN(status=1)
    ! write:
    call self%Line_Write( status )
    IF_NOTOK_RETURN(status=1)
        
    ! ok
    status = 0
  
  end subroutine PlotFile_AddLines


  ! *
  
  
  !
  ! Add line with filled polygon.
  ! The "alpha" in [0,1] defines the transparancy.
  ! With "edge=.true." also a line around the polygon is drawn.
  !
  
  subroutine PlotFile_AddPolygon( self, x, y, status, &
                                     xmap, kwargs, alpha, edge )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)        ::  self
    real, intent(in)                        ::  x(:), y(:)
    integer, intent(out)                    ::  status
    
    real, intent(in), optional              ::  xmap(2)   ! (scale,center)
    character(len=*), intent(in), optional  ::  kwargs
    real, intent(in), optional              ::  alpha
    logical, intent(in), optional           ::  edge

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_AddPolygon'
    
    ! --- local ----------------------------------
    
    character(len=1024)     ::  line
    character(len=1024)     ::  xx, yy
    integer                 ::  i
    real                    ::  xscale, xcenter
    
    ! --- begin ----------------------------------
    
    ! scale in x-direction around center?
    xscale = 1.0
    xcenter = 0.0
    if ( present(xmap) ) then
      xscale  = xmap(1)
      xcenter = xmap(2)
    endif
    
    ! check ..
    if ( size(x) /= size(y) ) then
      write (gol,'("x (",i0,") and y (",i0,") should have same lenght")') size(x), size(y); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init line:
    call self%Line_Add( 'ax.fill( ', status, clear=.true. )
    IF_NOTOK_RETURN(status=1)
    ! add x list:
    call self%Line_AddList( xcenter+(x-xcenter)*xscale, status )
    IF_NOTOK_RETURN(status=1)
    ! sep:
    self%line = trim(self%line)//','
    ! add y list:
    call self%Line_AddList( y, status )
    IF_NOTOK_RETURN(status=1)
    ! extra arguments?
    call self%Line_AddKwargs( status, kwargs=kwargs, alpha=alpha )
    IF_NOTOK_RETURN(status=1)
    ! close:
    call self%Line_Add( ')', status )
    IF_NOTOK_RETURN(status=1)
    ! write:
    call self%Line_Write( status )
    IF_NOTOK_RETURN(status=1)
    
    ! edge?
    if ( present(edge) ) then
      if ( edge ) then
        ! closed line, no alpha argument:
        call self%AddLines( x, y, status, xmap=xmap, kwargs=kwargs, closed=.true. )
        IF_NOTOK_RETURN(status=1)
      end if
    end if
        
    ! ok
    status = 0
  
  end subroutine PlotFile_AddPolygon


  ! *

  
  !
  ! Add plotting commands for rectangle (x1,x2) x (y1,y2).
  ! Keywords the same as for "AddPolygon".
  !
  
  subroutine PlotFile_AddRectangle( self, x,y, status, &
                                      xmap, kwargs, alpha, edge )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)        ::  self
    real, intent(in)                        ::  x(2), y(2)
    integer, intent(out)                    ::  status
    
    real, intent(in), optional              ::  xmap(2)   ! (scale,center)
    character(len=*), intent(in), optional  ::  kwargs
    real, intent(in), optional              ::  alpha
    logical, intent(in), optional           ::  edge

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_AddRectangle'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! add solid line:
    call self%AddPolygon( (/x(1),x(2),x(2),x(1)/), &
                          (/y(1),y(1),y(2),y(2)/), &
                          status, &
                          xmap=xmap, &
                          kwargs=kwargs, alpha=alpha, edge=edge )
    IF_NOTOK_RETURN(status=1)
        
    ! ok
    status = 0
  
  end subroutine PlotFile_AddRectangle


  ! *
  
  
  !
  ! Add plotting commands for arrow from (x,y) to (x2,y2).
  !
  
  subroutine PlotFile_AddArrow( self, x, y, x2, y2, status, &
                                  xmap, xmap2, kwargs )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)        ::  self
    real, intent(in)                        ::  x, y
    real, intent(in)                        ::  x2, y2
    integer, intent(out)                    ::  status
    
    real, intent(in), optional              ::  xmap(2)   ! (scale,center)
    real, intent(in), optional              ::  xmap2(2)   ! (scale,center)
    character(len=*), intent(in), optional  ::  kwargs

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_AddArrow'
    
    ! --- local ----------------------------------
    
    character(len=1024)     ::  line
    real                    ::  xscale, xcenter
    real                    ::  xscale2, xcenter2
    
    ! --- begin ----------------------------------
    
    ! scale in x-direction around center?
    xscale = 1.0
    xcenter = 0.0
    if ( present(xmap) ) then
      xscale  = xmap(1)
      xcenter = xmap(2)
    endif
    
    ! scale end point in x-direction around center?
    xscale2 = 1.0
    xcenter2 = 0.0
    if ( present(xmap2) ) then
      xscale2  = xmap2(1)
      xcenter2 = xmap2(2)
    endif
    
    ! use the annotation function to have proper arrows 
    ! also in an ax with non-unity aspect ratio:
    write (line,*) 'ax.annotate( "", (',xcenter2+(x2-xcenter2)*xscale2,',',y2,'), ', &
                    'xytext=(',xcenter+(x-xcenter)*xscale,',',y,'), ', &
                    'arrowprops=dict(arrowstyle="-|>"'
    ! store:
    call self%Line_Add( line, status, clear=.true. )
    IF_NOTOK_RETURN(status=1)
    ! extra arguments?
    call self%Line_AddKwargs( status, kwargs=kwargs )
    IF_NOTOK_RETURN(status=1)
    ! close:
    call self%Line_Add( '))', status )
    IF_NOTOK_RETURN(status=1)
    ! write:
    call self%Line_Write( status )
    IF_NOTOK_RETURN(status=1)
        
    ! ok
    status = 0
  
  end subroutine PlotFile_AddArrow
  
  ! *
  
  subroutine PlotFile_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PlotFile), intent(inout)          ::  self
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/PlotFile_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! add:
    write (self%fu,'(a)') '# show:'
    write (self%fu,'(a)') 'plt.show()'
    write (self%fu,'(a)') ''
    
    ! close:
    close( self%fu, iostat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine PlotFile_Done

  
end module GO_Plot
