!######################################################################
!
! Num - numerical tools
!
!######################################################################
!
#define TRACEBACK write (*,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action)  if (status> 0) then; TRACEBACK; action; return; end if
!
!######################################################################

module Num_Tools

  implicit none
  
  ! --- in/out -------------------
  
  private
  
  public   ::  Interval
  public   ::  Interval_modulo
  public   ::  Swap
  public   ::  GetInterpolWeights
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'Num_Tools' 

  
contains


  !------------------------------------------------------------------
  !
  ! Computes ileft = max( i, 1 <= i <= n, .and. x(i) < a )
  !
  ! Input:
  !       x:      a real sequence of length n, assumed nondecreasing.
  !       a:      the point whose location with tespect to the sequence
  !               x is to be determined.
  !       ileft:  a 'first guess' estimate of the index to be found.
  !               this variable should be set to 1 for the first call to
  !               INTERVAL.
  !               Interpretate 'ileft' as the interval number {1,..,n-1}
  !
  ! Output:
  !       ileft   iflag     meaning
  !       ------  --------  -------------------------------------
  !        1       -1        a  <  x(1)
  !        i        0              x(i)  <=  a  <  x(i+1)
  !        n-1      0                        a  =  x( n )
  !        n        1                              x( n )  <  a
  !
  !                   +-----+--- ... ---+-----------+---- ... ---+
  !
  !     boundary      1     2           i          i+1           n
  !     interval  1      1                    i              n-1     n
  !     status   -1      0                    0               0      1
  !
  !           (-inf,x1)
  !                   [x1,x2)
  !                                     [x(i),x(i+1))
  !                                                              [x(n),inf)
  !
  ! Method:
  !  Same as de Boor's INTERV sr, but not using a local index variable to
  !  store the last index found. Instead the 'first guess' index LEFT is
  !  used to initiate a search for the 'true' LEFT.
  !
  !------------------------------------------------------------------

  subroutine Interval( x, a, ileft, iflag )

    ! --- in/out ---------------------------------

    real, intent(in)        ::  x(:)
    real, intent(in)        ::  a
    integer, intent(inout)  ::  ileft
    integer, intent(out)    ::  iflag

    ! --- local ----------------------------------

    integer     ::  iright, istep, middle
    integer     ::  n

    ! --- begin -----------------------------------

    n = size(x)

    ! check correctness of first guess:
    if ( (ileft < 1) .or. (ileft>n) ) ileft = 1

    iright = ileft + 1
    if ( iright >= n ) then
      if ( a == x(n) ) then
        iflag = 0
        ileft = n-1
        return
      else if ( a > x(n) ) then
        iflag = 1
        ileft = n
        return
      end if
      if ( n <= 1 ) then
        iflag = -1
        ileft = 1
        return
      end if
      ileft = n-1
      iright = n
    end if

    if ( a >= x(iright) ) then

      ! now a >= x(iright). increase iright to capture a.
      istep = 1
      do
        ileft = iright
        iright = ileft + istep
        if ( iright >= n ) then
          if ( a == x(n) ) then
            iflag = 0
            ileft = n-1
            return
          else if ( a > x(n) ) then
            iflag = 1
            ileft = n
            return
          end if
          iright = n
          exit
        end if
        if ( a < x(iright) ) exit
        istep = istep*2
      end do

    else if ( a >= x(ileft) ) then

      iflag = 0
      return

    else

      ! now a  <  x(ileft).  decrease ileft to capture a

      istep = 1
      do
        iright = ileft
        ileft = iright - istep
        if ( ileft <= 1 ) then
          ileft = 1
          if ( a < x(1) ) then
            iflag = -1
            ileft = 1
            return
          end if
          exit
        end if
        if ( a >= x(ileft) ) exit
        istep = istep*2
      end do

    end if

    ! now x(ileft) <= x(iright). narrow the interval.

    do
      middle = ( ileft + iright ) / 2
      if ( middle == ileft ) then
        iflag = 0
        return
      end if
      if ( a < x(middle) ) then
        iright = middle
      else
        ileft = middle
      end if
    end do

  end subroutine Interval

! --------------------------------------------------

! Same procedure as Interval, but then such that the routine can handle with
! intervals which are shifted by a certain number of degrees.
! For example longitude -20 degree, must handle with a file with boundary conditions
! which has longitude 340 degree in the definition.

  subroutine Interval_modulo( x, ain, period, ileft, iflag )

    ! --- in/out ---------------------------------

    real, intent(in)        ::  x(:)
    real, intent(in)        ::  ain
    real, intent(in)        ::  period
    integer, intent(inout)  ::  ileft
    integer, intent(out)    ::  iflag

    ! --- local ----------------------------------

    integer     ::  iright, istep, middle
    integer     ::  n
    real        ::  a

    ! --- begin -----------------------------------

    a = x(1) + modulo(ain-x(1),period)
    n = size(x)

    ! check correctness of first guess:
    if ( (ileft < 1) .or. (ileft>n) ) ileft = 1

    iright = ileft + 1
    if ( iright >= n ) then
      if ( a == x(n) ) then
        iflag = 0
        ileft = n-1
        return
      else if ( a > x(n) ) then
        iflag = 1
        ileft = n
        return
      end if
      if ( n <= 1 ) then
        iflag = -1
        ileft = 1
        return
      end if
      ileft = n-1
      iright = n
    end if

    if ( a >= x(iright) ) then

      ! now a >= x(iright). increase iright to capture a.
      istep = 1
      do
        ileft = iright
        iright = ileft + istep
        if ( iright >= n ) then
          if ( a == x(n) ) then
            iflag = 0
            ileft = n-1
            return
          else if ( a > x(n) ) then
            iflag = 1
            ileft = n
            return
          end if
          iright = n
          exit
        end if
        if ( a < x(iright) ) exit
        istep = istep*2
      end do

    else if ( a >= x(ileft) ) then

      iflag = 0
      return

    else

      ! now a  <  x(ileft).  decrease ileft to capture a

      istep = 1
      do
        iright = ileft
        ileft = iright - istep
        if ( ileft <= 1 ) then
          ileft = 1
          if ( a < x(1) ) then
            iflag = -1
            ileft = 1
            return
          end if
          exit
        end if
        if ( a >= x(ileft) ) exit
        istep = istep*2
      end do

    end if

    ! now x(ileft) <= x(iright). narrow the interval.

    do
      middle = ( ileft + iright ) / 2
      if ( middle == ileft ) then
        iflag = 0
        return
      end if
      if ( a < x(middle) ) then
        iright = middle
      else
        ileft = middle
      end if
    end do

  end subroutine Interval_modulo

  ! ==========================================================
  
  !
  ! call Swap( lx, xx )
  !
  !   Swaps the elements in array xx with length lx .
  !  

  subroutine Swap( lx, xx )

    ! --- in/out -------------------------------
    
    integer, intent(in)       ::  lx
    real, intent(inout)       ::  xx(lx)
    
    ! --- local --------------------------------
    
    real        ::  swp
    integer     ::  l
    
    ! --- begin --------------------------------

    do l = 1, lx/2
      swp = xx(l)
      xx(l) = xx(lx-l+1)
      xx(lx-l+1) = swp
    end do

  end subroutine Swap
  
  
  ! *
  
  ! linear interpolation weights:
  
  subroutine GetInterpolWeights( xx, x0, n, ii, ww, status, period )
  
    ! --- in/out -------------------------------
  
    real, intent(in)        ::  xx(:)
    real, intent(in)        ::  x0
    integer, intent(out)    ::  n     ! involved points (1 or 2)
    integer, intent(out)    ::  ii(2)
    real, intent(out)       ::  ww(2)
    integer, intent(out)    ::  status

    real, intent(in), optional  ::  period

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GetInterpolWeights'
    
    ! --- local --------------------------------
    
    integer   ::  nx
    logical   ::  reverse
    real      ::  x
    integer   ::  i1, i2
    integer   ::  iflag

    ! --- begin --------------------------------
    
    ! size:
    nx = size(xx)
    
    ! normal order or reverse ?
    reverse = xx(2) < xx(1)
    
    ! shift if necessary:
    if ( present(period) ) then
      ! start at period base of first point in xx,
      ! add normalized target:
      if ( reverse ) then
        if ( x0 >= 0 .and. x0 < period ) then
          x = x0
        else 
          x = floor(xx(nx)/period)*period + modulo(x0,period)
        end if  
      else
        if ( x0 >= 0 .and. x0 < period ) then
          x = x0
        else 
          x = floor(xx( 1)/period)*period + modulo(x0,period)
        end if
      end if
    else
      ! copy:
      x = x0
    end if

    ! no points found yet:
    n = 0
    
    ! find interval with value:
    if ( reverse ) then
      call Interval( -xx, -x, i1, iflag )
    else
      call Interval(  xx,  x, i1, iflag )
    end if
    if ( iflag /= 0 ) then
      print *, 'ERROR - value ', x0 
      print *, 'ERROR - interval including ', x, ' not found in values:'      
      print *, 'ERROR - ', xx
      TRACEBACK; status=1; return
    end if
    ! exact ?
    if ( x == xx(i1) ) then
      ! add point:
      n = n + 1
      ii(n) = i1
      ww(n) = 1.0
    else
      ! select second point:
      if ( ((.not. reverse) .and. (x < xx(i1))) .or. &
           (       reverse  .and. (x > xx(i1)))      ) then
        i2 = i1 - 1
      else
        i2 = i1 + 1
      end if
      ! add first point:
      n     = n + 1
      ii(n) = i1
      ww(n) = ( xx(i2) - x )/( xx(i2) - xx(i1) )
      ! add second point:
      n     = n + 1
      ii(n) = i2
      ww(n) = 1.0 - ww(n-1)
    end if
    
    ! ok
    status = 0
    
  end subroutine GetInterpolWeights
  

end module num_tools
