!######################################################################
!
! Num - numerical tools - statistics
!
!######################################################################
!
#define TRACEBACK write (*,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action)  if (status> 0) then; TRACEBACK; action; return; end if
!
!######################################################################

module Num_Stat

  use GO, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out -------------------
  
  private
  
  public   ::  Mean, StdDev
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'Num_Stat'
  

  ! --- interfaces ---------------------------------------- 
  
  interface Mean
    module procedure Mean_2d
  end interface Mean
  
  interface StdDev
    module procedure StdDev_2d
  end interface StdDev

  
contains
  
  
  ! ====================================================================
  
  ! mean value
  
  subroutine Mean_2d( x, mu, status )
  
    ! --- in/out -------------------------------
  
    real, intent(in)        ::  x(:,:)
    real, intent(out)       ::  mu
    integer, intent(out)    ::  status

    ! --- const --------------------------------

    character(len=*), parameter  ::  rname = mname//'/Mean_2d'

    ! --- local --------------------------------
    
    integer   ::  n

    ! --- begin --------------------------------
    
    ! total size:
    n = size(x)
    
    ! mean value:
    mu = sum(x) / real(n)
    
    ! ok
    status = 0
    
  end subroutine Mean_2d
  
  ! *
  
  ! standard deviation
  
  subroutine StdDev_2d( x, sigma, status )
  
    ! --- in/out -------------------------------
  
    real, intent(in)        ::  x(:,:)
    real, intent(out)       ::  sigma
    integer, intent(out)    ::  status

    ! --- const --------------------------------

    character(len=*), parameter  ::  rname = mname//'/StdDev_2d'

    ! --- local --------------------------------
    
    real      ::  mu
    integer   ::  n

    ! --- begin --------------------------------
    
    ! total size:
    n = size(x)
    
    ! check ...
    if ( n < 2 ) then
      write (gol,'("standard deviation requires at least 2 elements, found ",i0)') n; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! mean:
    call Mean( x, mu, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! std.dev. value:
    sigma = sqrt( sum( ( x - mu )**2 ) / real(n-1) )
    
    ! ok
    status = 0
    
  end subroutine StdDev_2d


end module Num_Stat

