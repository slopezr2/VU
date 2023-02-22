!###############################################################################
!
! NAME
!   LE_Install - communication between domains I/O task
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action)  if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Install

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Install_Init, LE_Install_Done
  public  ::  LE_Install_Step


  ! --- const ------------------------------------
    
  character(len=*), parameter   ::  mname = 'LE_Install'
  

  ! --- var --------------------------------------
  
  ! script and standard arguments:
  character(len=1024)       ::  script
  ! flag:
  logical                   ::  apply_install
  logical                   ::  is_first



contains


  ! ====================================================================
  
  
  subroutine LE_Install_Init( rcF, status )
  
    use GO     , only : TrcFile
      
    ! --- in/out ---------------------------------
    
    integer, intent(out)                      ::  status
    type(TrcFile), intent(in)                 ::  rcF
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Install_Init'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! call install script?
    call rcF%Get( 'le.io.apply', apply_install, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! enabled?
    if ( apply_install ) then

      ! script and arguments:
      call rcF%Get( 'le.io.command', script, status )
      IF_NOT_OK_RETURN(status=1)
    
      ! set flag:
      is_first = .true.
      
    end if

    ! ok
    status = 0

  end subroutine LE_Install_Init
  
  
  ! *
  
  
  subroutine LE_Install_Step( t, status )
  
    use GO     , only : TDate, IncrDate, IsAnyDate
    use GO     , only : operator(<), operator(-)
    use GO     , only : MidNight, Get_End_Of
    use GO     , only : goSystem
    use LE_Comm, only : goc_io
    use Dims   , only : runF
      
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Install_Step'
    
    ! --- local ----------------------------------
  
    character(len=1024)       ::  command
    character(len=64)         ::  tstamp
    character(len=64)         ::  fcstamp
    type(TDate)               ::  t24
    type(TDate)               ::  tprev

    ! --- begin ----------------------------------

    ! call script?
    if ( apply_install ) then

      ! first install?
      if ( is_first ) then
        ! processing task?
        if ( goc_io%last ) then
          ! also yesterday? current meteo files have data (00,24] ..
          if ( MidNight(t) ) then
            ! yesterday:
            tprev = t - IncrDate(day=1)
            ! timerange of 2 days:
            write (tstamp,'(i4.4,2("-",i2.2)," upto ",i4.4,2("-",i2.2)," by 1 day")') &
                       tprev%year, tprev%month, tprev%day, t%year, t%month, t%day
          else
            ! just today:
            write (tstamp,'(i4.4,2("-",i2.2))') t%year, t%month, t%day
          end if
          ! Forecast made:         
          if ( .not. IsAnyDate(runF%t_base) ) then
            write( fcstamp,'(i4.4,2("-",i2.2))') runF%t_base%year, runF%t_base%month, runF%t_base%day
          else 
            write( fcstamp,'("None")' )
          end if
          
          ! install today:
          write (command,'(a," --install-timerange=""",a,""""," --forecast-base=""",a,"""")') &
                             trim(script), trim(tstamp), trim(fcstamp)     
          ! testing ...
          write (gol,'(a,":       ",a)') rname, trim(command); call goPr
          ! run:
          call goSystem( command, status )
          IF_NOT_OK_RETURN(status=1)
          ! info ..
          write (gol,'(a,":       ready")') rname; call goPr
        else
          ! info ...
          write (tstamp,'(i4.4,2("-",i2.2))') t%year, t%month, t%day
          write (gol,'(a,":     wait for input for ",a," to be ready ..")') rname, trim(tstamp); call goPr
        end if
        ! wait:
        call goc_io%Barrier( status )
        IF_NOT_OK_RETURN(status=1)
        ! info ...
        write (tstamp,'(i4.4,2("-",i2.2))') t%year, t%month, t%day
        write (gol,'(a,":     input for start day ",a," ready")') rname, trim(tstamp); call goPr
      end if ! first install
      
      ! at midnight (or at start), install for next day (24-48 hours ahead):
      if ( MidNight(t) .or. is_first ) then

        ! processing task?
        if ( goc_io%last ) then
          write (gol,'(a,":     wait for other processors ...")') rname; call goPr
        end if
        ! wait here to ensure that all data for 00-24 is ready ..
        call goc_io%Barrier( status )
        IF_NOT_OK_RETURN(status=1)
        ! today:
        write (tstamp,'(i4.4,2("-",i2.2))') t%year, t%month, t%day
        write (gol,'(a,":     input for ",a," ready")') rname, trim(tstamp); call goPr
      
        ! processing task?
        if ( goc_io%last ) then
          ! next midnight:
          t24 = Get_End_Of( t, 'day' )
          ! prepare next?
          if ( t24 < runF%t_end ) then
            ! next day:
            write (tstamp,'(i4.4,2("-",i2.2))') t24%year, t24%month, t24%day
            ! install today:
            ! Forecast made:
            if ( .not. IsAnyDate(runF%t_base) ) then
              write( fcstamp,'(i4.4,2("-",i2.2))') runF%t_base%year, runF%t_base%month, runF%t_base%day
            else 
              write( fcstamp,'("None")' )
            end if          
            ! install today:
            write (command,'(a," --install-timerange=""",a,""""," --forecast-base=""",a,"""")') &
                             trim(script), trim(tstamp), trim(fcstamp)     
    
            ! cleanup day-before-yesterday:
            tprev = t - IncrDate(day=4)
            write (tstamp,'(i4.4,2("-",i2.2))') tprev%year, tprev%month, tprev%day
            write (command,'(a," --cleanup-timerange=""",a,""""," --forecast-base=""",a,"""")') &
                              trim(command), trim(tstamp), trim(fcstamp)
            ! testing ...
            write (gol,'(a,":       ",a)') rname, trim(command); call goPr
            ! run:
            call goSystem( command, status )
            IF_NOT_OK_RETURN(status=1)
            ! info ..
            write (gol,'(a,":       ready")') rname; call goPr
          end if  ! prepare next day
        end if ! i/o processor

      end if ! midnight or first

      ! reset flag:
      is_first = .false.

    end if ! call install script

    ! ok
    status = 0

  end subroutine LE_Install_Step
  
  
  ! ***
  
  
  subroutine LE_Install_Done( t, status )

    use GO     , only : TDate, IncrDate, IsAnyDate, MidNight, operator(-), max
    use GO     , only : goSystem
    use LE_Comm, only : goc_io
    use Dims   , only : runF
  
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)           ::  t
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Install_Done'
    
    ! --- local ----------------------------------
  
    character(len=1024)       ::  command
    character(len=64)         ::  tstamp
    character(len=64)         ::  fcstamp
    type(TDate)               ::  tm1
    
    ! --- begin ----------------------------------
    
    ! enabled?
    if ( apply_install ) then
    
      ! part of i/o network?
      if ( goc_io%enabled ) then

        ! wait here to remove data that is still being used ..
        call goc_io%Barrier( status )
        IF_NOT_OK_RETURN(status=1)

        ! i/o processing task?
        if ( goc_io%last ) then
          ! info ..
          write (gol,'(a,": cleanup ..")') rname; call goPr
          ! fill base command:
          command = '../build/bin/le-copy-input lotos-euros.rc'
          ! cleanup from day-before-yesterday to today:
          if ( MidNight(t) ) then
            ! Run finished at midnight, so t is already one day later
            tm1 = max( runF%t_start, t - IncrDate(day=5) )
          else 
            tm1 = max( runF%t_start, t - IncrDate(day=4) )
          endif
          write (tstamp,'(i4.4,2("-",i2.2)," upto ",i4.4,2("-",i2.2)," by 1 day")') &
                   tm1%year, tm1%month, tm1%day, t%year, t%month, t%day
          
          ! Forecast made:
          if ( .not. IsAnyDate(runF%t_base) ) then
            write( fcstamp,'(i4.4,2("-",i2.2))') runF%t_base%year, runF%t_base%month, runF%t_base%day
          else 
            write( fcstamp,'("None")' )
          end if          
          ! Command line    
          write (command,'(a," --cleanup-timerange=""",a,""""," --forecast-base=""",a,"""")') &
                             trim(command), trim(tstamp), trim(fcstamp)
          ! testing ...
          write (gol,'(a,":   ",a)') rname, trim(command); call goPr
          ! run:
          call goSystem( command, status )
          IF_NOT_OK_RETURN(status=1)
          ! info ..
          write (gol,'(a,":   ready")') rname; call goPr
        else
          ! info ..
          write (gol,'(a,": wait for cleanup ..")') rname; call goPr
        end if
        ! info ..
        write (gol,'(a,":   cleaned up")') rname; call goPr

      end if ! i/o network

    end if ! call install script
    
    ! ok
    status = 0

  end subroutine LE_Install_Done



end module LE_Install


