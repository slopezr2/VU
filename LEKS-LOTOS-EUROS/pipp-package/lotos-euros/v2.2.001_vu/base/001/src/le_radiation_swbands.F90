!###############################################################################
!
! Define radiation shortwave bands
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Radiation_SWBands
 
  use GO, only : gol, goPr, goErr
 
  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  T_SWBands


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Radiation_SWBands'


  ! --- types ------------------------------------

  type T_SWBands
    ! number of bands:
    integer               ::  n
    ! band mid and borders:
    real, allocatable     ::  lambda(:)    ! (n)     ! um
    real, allocatable     ::  bounds(:,:)  ! (n,2)   ! um
    real, allocatable     ::  width(:)     ! (n)     ! um
  contains
    procedure   ::               SWBands_Init_n
    procedure   ::               SWBands_Init_rc
    generic     ::  Init      => SWBands_Init_n, &
                                 SWBands_Init_rc
    procedure   ::  Done      => SWBands_Done
    procedure   ::  SetBand   => SWBands_SetBand
    procedure   ::  FindBand  => SWBands_FindBand
  end type T_SWBands


contains

  ! ====================================================================
  ! ===
  ! === SW bands
  ! ===
  ! ====================================================================


  subroutine SWBands_Init_n( self, n, status )

    ! --- in/out ---------------------------------

    class(T_SWBands), intent(out)   ::  self
    integer, intent(in)             ::  n
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/SWBand_Init_n'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%n = n
    
    ! storage:
    allocate( self%lambda(self%n), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%bounds(self%n,2), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%width(self%n), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! undefined yet:
    self%lambda = -999.9
    self%bounds = -999.9
    self%width  = -999.9
    
    ! ok
    status = 0
    
  end subroutine SWBands_Init_n
  
  
  ! *


  subroutine SWBands_Init_rc( self, rcF, status )

    use GO, only : TrcFile, ReadRc

    ! --- in/out ---------------------------------

    class(T_SWBands), intent(out)   ::  self
    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/SWBand_Init_rc'

    ! --- local ----------------------------------
    
    character(len=32)   ::  name
    character(len=256)  ::  rcbase
    character(len=256)  ::  rckey
    character(len=256)  ::  line
    integer             ::  n
    integer             ::  i
    real                ::  bounds(2)

    ! --- begin ----------------------------------
    
    ! name of swband definition:
    call ReadRc( rcF, 'le.radiation.swbands.name', name, status )
    IF_NOTOK_RETURN(status=1)
    
    ! first part of rcfile keys:
    rcbase = 'le.radiation.swbands.'//trim(name)
    
    ! number of bands:
    call ReadRc( rcF, trim(rcbase)//'.n', n, status )
    IF_NOTOK_RETURN(status=1)
    
    ! store:
    self%n = n
    
    ! storage:
    allocate( self%lambda(self%n), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%bounds(self%n,2), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%width(self%n), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! undefined yet:
    self%lambda = -999.9
    self%bounds = -999.9
    self%width  = -999.9
    
    ! loop:
    do i = 1, n
      ! definition key for this band:
      write (rckey,'(a,".",i2.2)') trim(rcbase), i
      ! read line:
      call ReadRc( rcF, rckey, line, status )
      IF_NOTOK_RETURN(status=1)
      ! read bounds:
      read (line,*) bounds
      ! fill:
      call self%SetBand( i, status, bounds=bounds )
      IF_NOTOK_RETURN(status=1)
    end do ! i

    ! ok
    status = 0
    
  end subroutine SWBands_Init_rc


  ! *
  

  subroutine SWBands_Done( self, status )

    ! --- in/out ---------------------------------

    class(T_SWBands), intent(inout)     ::  self
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/SWBand_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%lambda, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%bounds, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%width, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! undefined yet:
    self%n = -999
  
    ! ok
    status = 0
    
  end subroutine SWBands_Done


  ! *
  

  subroutine SWBands_SetBand( self, i, status, &
                               lambda, bounds, width )

    ! --- in/out ---------------------------------

    class(T_SWBands), intent(inout)       ::  self
    integer, intent(in)                 ::  i
    integer, intent(out)                ::  status
    real, intent(in), optional          ::  lambda
    real, intent(in), optional          ::  bounds(2)
    real, intent(in), optional          ::  width

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/SWBand_SetBand'

    ! --- local ----------------------------------
    
    real          ::  bnds(2)
    integer       ::  k

    ! --- begin ----------------------------------
    
    ! check ...
    if ( (i < 1) .or. (i > self%n) ) then
      write (gol,'("index ",i0," not in range 1,..,",i0)') i, self%n; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! fill bounds:
    !
    if ( present(bounds) ) then
      ! check ...
      if ( present(lambda) .or. present(width) ) then
        write (gol,'("argument `bounds` not in combination with `lambda` or `width`")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      bnds = bounds
    !
    else if ( present(lambda) .and. present(width) ) then 
      ! check ...
      if ( present(bounds) ) then
        write (gol,'("arguments `lambdsa` and `width` not in combination with `bounds`")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! set bounds:
      bnds(1) = lambda - 0.5*width
      bnds(2) = lambda + 0.5*width
    !
    else
      write (gol,'("provide argument `bounds`, or `lambda` or `width`")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check for overlap:
    do k = 1, self%n
      ! defined ?
      if ( self%lambda(k) > 0.0 ) then
        ! overlap?
        if ( ( (self%bounds(k,1)+0.00001 < bnds(1)) .and. (bnds(1) < self%bounds(k,2)-0.00001) ) .or. &
             ( (self%bounds(k,1)+0.00001 < bnds(2)) .and. (bnds(2) < self%bounds(k,2)-0.00001) )      ) then
          write (gol,'("provided band ",i0," (",f8.3,",",f8.3,") overlaps with existing ",i0," (",f8.3,",",f8.3,")")') &
                  i, bnds, k, self%bounds(k,:); call goErr
          TRACEBACK; status=1; return
        end if ! overlap?
      end if ! defined
    end do ! bands
    
    ! store:
    self%lambda(i) = 0.5 * ( bnds(1) + bnds(2) )  ! mid value
    self%bounds(i,:) = bnds  ! boundaries
    self%width(i) = self%bounds(i,2) - self%bounds(i,1) ! width

    !! info ...
    !write (gol,'("defined band ",i4," ",f8.3," (",f8.3,",",f8.3,") ",f8.3)') &
    !              i, self%lambda(i), self%bounds(i,:), self%width(i); call goPr

    ! ok
    status = 0
    
  end subroutine SWBands_SetBand


  ! *
  

  subroutine SWBands_FindBand( self, lambda, i, status )

    ! --- in/out ---------------------------------

    class(T_SWBands), intent(in)        ::  self
    real, intent(in)                    ::  lambda
    integer, intent(out)                ::  i
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/SWBand_FindBand'

    ! --- local ----------------------------------
    
    integer   ::  k
    integer   ::  nmatch

    ! --- begin ----------------------------------
    
    ! init counter:
    nmatch = 0
    ! loop over bands:
    do k = 1, self%n
      ! defined ?
      if ( self%lambda(k) > 0.0 ) then
        ! compare:
        if ( (self%bounds(k,1) <= lambda) .and. (lambda <= self%bounds(k,2)) ) then
          ! store:
          i = k
          ! increase counter:
          nmatch = nmatch + 1
        end if
      end if ! defined?
    end do  ! bands
    ! check ...
    if ( nmatch /= 1 ) then
      write (gol,'(i0," matches found for lambda = ",f8.3)') nmatch, lambda; call goErr
      do k = 1, self%n
        if ( self%lambda(k) > 0.0 ) then
          if ( (self%bounds(k,1) <= lambda) .and. (lambda <= self%bounds(k,2)) ) then
            write (gol,'("  band ",i3," (",f8.3,",",f8.3,")  <--")') k, self%bounds(k,:); call goErr
          else
            write (gol,'("  band ",i3," (",f8.3,",",f8.3,")")') k, self%bounds(k,:); call goErr
          end if
        end if
      end do
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine SWBands_FindBand


end module LE_Radiation_SWBands



