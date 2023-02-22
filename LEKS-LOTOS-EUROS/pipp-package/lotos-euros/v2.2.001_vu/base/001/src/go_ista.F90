!#######################################################################
!
! Incremental statistics: temporal averages.
!
! USAGE
!
!   use GO_ISTA, only : T_ISTA, ISTA_Init, ISTA_Done, ISTA_Add, ISTA_Next, ISTA_Shift
!
!   ! data size:
!   integer, parameter  ::  n = 5
!   ! storage:
!   real                ::  sample(n)
!   real                ::  sample_aver(n)
!   real                ::  sample_prev(n)
!   ! averaging info:
!   type(T_ISTA)        ::  ista
!   ! time resolution:
!   type(TDate)         ::  t0
!   type(TIncrDate)     ::  dt
!   logical             ::  finished_aver
!   
!   ! timing:
!   t0 = NewDate( 2006, 01, 01, 00, 00 )
!   dt = IncrDate( day=1 )
!
!   ! setup temporal averages starting at 't0' with step 'dt'
!   call ISTA_Init( ista, t0, dt, status )
!   if (status/=0) stop
!   sample_aver = 0.0
!
!   ! init samples:
!   sample = 0.0
!
!   ! time loop:
!   do
!
!     ! next sample:
!     sample = sample + 1.0
!
!     ! add sample:
!     call ISTA_Add( ista, t, w_aver, w_prev, w_samp, finished_aver, status )
!     if (status/=0) stop
!     sample_aver = sample_aver * w_aver + values_prev * w_prev + sample * w_samp
!
!     ! finished_aver ? write sample:
!     if ( finished_aver ) then
!       ! write:
!       call putout( sample_aver )
!       ! init new interval:
!       call ISTA_Next( ista, w_prev, w_samp, status )
!       if (status/=0) stop
!       sample_aver = sample_prev * w_prev + sample * w_samp
!     end if
!
!     ! store current sample:
!     call ISTA_Shift( ista, status )
!     if (status/=0) stop
!     values_prev = sample
!
!   end do
!
!   ! clear:
!   call ISTA_Done( ista, status )
!   if (status/=0) stop
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#######################################################################


module GO_ISTA

  use GO_Print, only : gol, goPr, goErr
  use GO_Date , only : TDate, TIncrDate
  
  implicit none
  
  
  ! --- in/out -----------------------------------

  private
  
  public  ::  T_ISTA
  public  ::  ISTA_Init, ISTA_Done
  public  ::  ISTA_Add, ISTA_Shift, ISTA_Next


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'GO_ISTA'
  
  
  ! --- types ------------------------------------
  
  type T_ISTA
    ! flags:
    logical               ::  isfirst
    ! length of average interval:
    type(TIncrDate)       ::  dt
    ! current target interval for averaging:
    type(TDate)           ::  ival_t1, ival_t2
    ! actual average interval:
    type(TDate)           ::  aver_t1, aver_t2
    ! previous sample:
    type(TDate)           ::  prev_t
    ! added sample:
    type(TDate)           ::  samp_t
  end type T_ISTA
  
  

contains


  ! ====================================================================
  
  
  subroutine ISTA_Init( ista, t0, dt, status )
  
    use GO_Date, only : TDate, TIncrDate, AnyDate, operator(+)
    
    ! --- in/out ---------------------------------
    
    type(T_ISTA), intent(out)               ::  ista
    type(TDate), intent(in)                 ::  t0
    type(TIncrDate), intent(in)             ::  dt
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ISTA_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store interval length:
    ista%dt = dt
    
    ! setup first interval:
    ista%ival_t1 = t0
    ista%ival_t2 = ista%ival_t1 + dt
    
    ! dummy vales:
    ista%aver_t1 = AnyDate()
    ista%aver_t2 = AnyDate()
    ista%prev_t  = AnyDate()
    ista%samp_t  = AnyDate()
    
    ! set flag:
    ista%isfirst = .true.
    
    ! ok
    status = 0
    
  end subroutine ISTA_Init


  ! ***
  
  
  subroutine ISTA_Done( ista, status )
  
    ! --- in/out ---------------------------------
    
    type(T_ISTA), intent(inout)             ::  ista
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ISTA_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! dummy to avoid warnings about unused variables:
    ista%isfirst = .true.
    
    ! ok
    status = 0
    
  end subroutine ISTA_Done


  ! ***
  
  
  subroutine ISTA_Add( ista, t, w_aver, w_prev, w_samp, finished_aver, status )
  
    use GO_Date, only : wrtgol, operator(<), operator(<=), operator(>=)
    use GO_Date, only : operator(-), operator(+), rTotal
    
    ! --- in/out ---------------------------------
    
    type(T_ISTA), intent(inout)             ::  ista
    type(TDate), intent(in)                 ::  t
    real, intent(out)                       ::  w_aver, w_prev, w_samp
    logical, intent(out)                    ::  finished_aver
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ISTA_Add'
    
    ! --- local ----------------------------------
    
    type(TDate)   ::  t_add
    real          ::  dsec_aver, dsec_samp, dsec_add
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( t < ista%ival_t1 ) then
      write (gol,'("sample time prior to current average interval:")'); call goErr
      call wrtgol( '  sample time      : ', t ); call goErr
      call wrtgol( '  average interval : ', ista%ival_t1, ' - ', ista%ival_t2 ); call goErr
      TRACEBACK; status=1; return
    end if
    
    !
    !                                                t   current sample
    !                                      t_add         interpolation between latest and current
    !                           prev_t                   latest sample
    !            aver_t1       aver_t2                   previously avaraged
    !   ival_t1                          ival_t2         target interval for average
    !  ----+--------+-------------+---------+--------+-----------
    !
    !                             |<------->|
    !                               dsec_add
    !               |<----------->|<---------------->|
    !                   dsec_aver        dsec_samp
    !
    
    ! store sample time:
    ista%samp_t = t
    
    ! first sample ?
    if ( ista%isfirst ) then
    
      ! check ...
      if ( ista%ival_t2 < t ) then
        write (gol,'("first sample time after current average interval:")'); call goErr
        call wrtgol( '  sample time      : ', t ); call goErr
        call wrtgol( '  average interval : ', ista%ival_t1, ' - ', ista%ival_t2 ); call goErr
        TRACEBACK; status=1; return
      end if

      ! set weights:
      w_aver = 0.0    ! no previous average yet
      w_prev = 0.0    ! no previous sample yet
      w_samp = 1.0    ! set average equal to sample
      
      ! 'previous' sample:
      ista%prev_t = t

      ! update average info:
      ista%aver_t1  = t
      ista%aver_t2  = t

      ! reset flag:
      ista%isfirst = .false.
      
    else
    
      ! check ...
      if ( t <= ista%prev_t ) then
        write (gol,'("sample time prior or equal to previous sample:")'); call goErr
        call wrtgol( '  sample time      : ', t ); call goErr
        call wrtgol( '  previous sample  : ', ista%prev_t ); call goErr
        TRACEBACK; status=1; return
      end if

      ! add (interpolated) sample at time t_add:
      if ( t >= ista%ival_t2 ) then
        t_add = ista%ival_t2
      else
        t_add = t
      end if

      ! lenght of added interval:
      dsec_aver = rTotal( ista%aver_t2 - ista%aver_t1, 'sec' )
      dsec_add  = rTotal( t_add - ista%prev_t, 'sec' )
      dsec_samp = rTotal( t     - ista%prev_t, 'sec' )
      
      ! set weights:
      w_aver =                                     dsec_aver / (dsec_aver + dsec_add)
      w_prev = (dsec_samp-0.5*dsec_add)/dsec_samp * dsec_add / (dsec_aver + dsec_add)
      w_samp = (          0.5*dsec_add)/dsec_samp * dsec_add / (dsec_aver + dsec_add)
      
      ! extend average interval:
      ista%aver_t2  = t_add
    
    end if
    
    !! info ...
    !call wrtgol( '  ista: ival ', ista%ival_t1, ' - ', ista%ival_t2 ); call goPr
    !call wrtgol( '  ista: aver ', ista%aver_t1, ' - ', ista%aver_t2 ); call goPr
    !call wrtgol( '  ista: samp ', t ); call goPr

    ! reached end of target interval ?
    finished_aver = t >= ista%ival_t2
      
    ! ok
    status = 0
    
  end subroutine ISTA_Add


  ! ***
  
  
  subroutine ISTA_Next( ista, w_prev, w_samp, status )
  
    use GO_Date, only : rTotal, operator(-), operator(+), operator(>)
  
    ! --- in/out ---------------------------------
    
    type(T_ISTA), intent(inout)             ::  ista
    real, intent(out)                       ::  w_prev, w_samp
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ISTA_Next'
    
    ! --- local ----------------------------------
    
    real    ::  prev_sec, add_sec
    
    ! --- begin ----------------------------------
    
    ! fill weights:
    if ( ista%samp_t > ista%ival_t2 ) then
      ! interpolate over interval border:
      prev_sec = rTotal( ista%samp_t - ista%prev_t , 'sec' )
      add_sec  = rTotal( ista%samp_t - ista%ival_t2, 'sec' )
      w_prev = 0.5*add_sec/prev_sec
      w_samp = 1.0 - w_prev
      ista%aver_t1 = ista%ival_t2
      ista%aver_t2 = ista%samp_t
    else
      ! new sample valid at start of new interval, thus no previous info needed:
      w_prev = 0.0
      w_samp = 1.0
      ista%aver_t1 = ista%samp_t
      ista%aver_t2 = ista%samp_t
    end if
    
    ! new interval:
    ista%ival_t1 = ista%ival_t1 + ista%dt
    ista%ival_t2 = ista%ival_t1 + ista%dt
    
    ! ok
    status = 0
    
  end subroutine ISTA_Next


  ! ***
  
  
  subroutine ISTA_Shift( ista, status )
  
    ! --- in/out ---------------------------------
    
    type(T_ISTA), intent(inout)             ::  ista
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ISTA_Shift'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! current sample becomes previous:
    ista%prev_t = ista%samp_t
    
    ! ok
    status = 0
    
  end subroutine ISTA_Shift


end module GO_ISTA


! ######################################################################
! ###
! ### test
! ###
! ######################################################################
!
! gfortran -o test.x -ffree-line-length-none -fbounds-check go_fu.F90 go_print.F90 go_string.F90 go_date.F90 go_file.F90 go_rc.F90 go_system.F90 go.F90 ista.F90
!
!#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!
!program test
!
!  use GO, only : gol, goPr, goErr
!  use GO, only : TDate, TIncrDate, NewDate, IncrDate
!  use GO, only : wrtgol, operator(+), operator(>=), rTotal
!  use ISTA
!  
!  implicit none
!  
!  ! --- const ------------------------------------
!
!  character(len=*), parameter   ::  rname = 'test'
!  
!  ! number of parameters:
!  integer, parameter    ::  n = 1
!  
!  ! --- local ------------------------------------
!  
!  integer           ::  status
!  type(T_ISTA)      ::  ista
!  type(TDate)       ::  t0, t
!  type(TIncrDate)   ::  dt
!  real              ::  w_aver, w_prev, w_samp
!  logical           ::  finished_aver
!  
!  ! storage:
!  real              ::  sample_aver(n)
!  real              ::  sample_prev(n)
!  real              ::  sample(n)
!  
!  ! --- begin ------------------------------------
!  
!  print *, 'testing ista ...'
!  
!  ! timing:
!  t0 = NewDate( 2006, 01, 01, 00, 00 )
!  dt = IncrDate( day=1 )
!
!  ! setup temporal averages starting at 't' with step 'dt'
!  call ISTA_Init( ista, t0, dt, status )
!  IF_NOTOK_STOP
!  sample_aver = 0.0
!  
!  ! start time:
!  t = t0
!  sample = 0.0
!  ! add sample:
!  call ISTA_Add( ista, t, w_aver, w_prev, w_samp, finished_aver, status )
!  IF_NOTOK_STOP
!  sample_aver = sample_aver * w_aver + sample_prev * w_prev + sample * w_samp
!  print *, '    weights : ', w_aver, w_prev, w_samp, ';', sample_aver(1)
!
!  ! start time:
!  t = t0 + IncrDate(hour=3)
!  sample = t%hour - 6.0
!
!  ! time loop:
!  do
!  
!    ! info ...
!    call wrtgol( '  time : ', t ); call goPr
!    
!    ! fill sample:
!    sample = sample + 6.0
!
!    ! add sample:
!    call ISTA_Add( ista, t, w_aver, w_prev, w_samp, finished_aver, status )
!    IF_NOTOK_STOP
!    sample_aver = sample_aver * w_aver + sample_prev * w_prev + sample * w_samp
!    print *, '    weights : ', w_aver, w_prev, w_samp, ';', sample_aver(1)
!
!    ! finished_aver ? write sample:
!    if ( finished_aver ) then
!      ! write:
!      print *, '      ---> write ', sample_aver(1)
!      ! init new interval:
!      call ISTA_Next( ista, w_prev, w_samp, status )
!      IF_NOTOK_STOP
!      sample_aver = sample_prev * w_prev + sample * w_samp
!      print *, '    next : ', 0.0, w_prev, w_samp, ';', sample_aver(1)
!    end if
!
!    ! save current sample:
!    call ISTA_Shift( ista, status )
!    IF_NOTOK_STOP
!    sample_prev = sample
!
!    ! end ?
!    if ( t >= NewDate(2006,01,03,00,00) ) exit
!    ! next sample:
!    t = t + IncrDate(hour=6)
!
!  end do
!
!  ! clear:
!  call ISTA_Done( ista, status )
!  IF_NOTOK_STOP
!
!  print *, 'ok'
!  
!end program test
