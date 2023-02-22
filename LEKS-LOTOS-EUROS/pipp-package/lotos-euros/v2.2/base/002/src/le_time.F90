!###############################################################################
!
! NAME
!   LE_Time  -  LOTOS-EUROS time stuff
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


module LE_Time

  use GO  , only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Time_Init, LE_Time_Done
  
  public  ::  date_le_date
  public  ::  local_time2
  public  ::  day_of_the_week

  
  ! --- const ------------------------------------
    
  character(len=*), parameter   ::  mname = 'LE_Time'
  
  
  ! --- var --------------------------------------
  
contains


  ! ====================================================================
  
  
  subroutine LE_Time_Init( rcF, status )
  
    use GO   , only : TrcFile, ReadRc
    use GO   , only : NewDate, AnyDate
    use GO   , only : goTranslate
    use Dims , only : runF
    use go   , only : operator(<=)
      
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)         ::  rcF
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Time_Init'
    
    ! --- local ----------------------------------

    character(len=3)    ::  yesno
    character(len=32)   ::  stime
    integer             ::  minu, seco

    ! --- begin ----------------------------------
    
    !
    ! time range
    !
    
    ! start time:
    call ReadRc( rcF, 'timerange.start' , stime, status )
    IF_NOTOK_RETURN(status=1)
    call goTranslate( stime, '/-:', ' ', status )
    IF_NOTOK_RETURN(status=1)
    read (stime,*,iostat=status) runF%yy_s, runF%mm_s, runF%dd_s, runF%hh_s, minu, seco
    if ( status /= 0 ) then
      write (gol,'("could not read start time from : ",a)') trim(stime); call goErr
      TRACEBACK; status=1; return
    end if
    ! info ...
    write (gol,'("simulation will start at: ",i4.4,2("-",i2.2)," ",i2.2,":00")') runF%yy_s, runF%mm_s, runF%dd_s, runF%hh_s; call goPr
    ! store in structure:
    runF%t_start = NewDate( year=runF%yy_s, month=runF%mm_s, day=runF%dd_s, hour=runF%hh_s )

    ! start time:
    call ReadRc( rcF, 'timerange.end' , stime, status )
    IF_NOTOK_RETURN(status=1)
    call goTranslate( stime, '/-:', ' ', status )
    IF_NOTOK_RETURN(status=1)
    read (stime,*,iostat=status) runF%yy_e, runF%mm_e, runF%dd_e, runF%hh_e, minu, seco
    if ( status /= 0 ) then
      write (gol,'("could not read end time from : ",a)') trim(stime); call goErr
      TRACEBACK; status=1; return
    end if
    ! info ...
    write (gol,'("simulation will end at:      ",i4.4,2("-",i2.2)," ",i2.2,":00")') runF%yy_e, runF%mm_e, runF%dd_e, runF%hh_e; call goPr
    !
    ! store in structure:
    runF%t_end = NewDate( year=runF%yy_e, month=runF%mm_e, day=runF%dd_e, hour=runF%hh_e )
    
    if ( runF%t_end <= runF%t_start ) then
      write(gol, '("end time is before start time") '); call goErr
      TRACEBACK; status=1; return
    end if

    !
    ! forecast mode
    !

    ! forecast mode ?
    call ReadRc( rcF, 'forecast.mode', runF%forecast_mode, status )
    IF_NOTOK_RETURN(status=1)

    ! forecast mode ?
    if ( runF%forecast_mode ) then
      ! read forecast base time:
      call ReadRc( rcF, 'forecast.timerange.base' , stime, status )
      IF_NOTOK_RETURN(status=1)
      call goTranslate( stime, '/-:', ' ', status )
      IF_NOTOK_RETURN(status=1)
      read (stime,*,iostat=status) runF%yy_b, runF%mm_b, runF%dd_b, runF%hh_b, minu, seco
      if ( status /= 0 ) then
        write (gol,'("could not read forecast base time from : ",a)') trim(stime); call goErr
        TRACEBACK; status=1; return
      end if
      ! store in structure:
      runF%t_base = NewDate( year=runF%yy_b, month=runF%mm_b, day=runF%dd_b, hour=runF%hh_b )
    !
    else
      ! dummy:
      runF%t_base = AnyDate()
    end if
    
    !
    ! logging
    !
    
    ! write details to the log file
    write (gol,'("start date of the simulation at: ",i4,2("-",i2.2)," ",i2.2,":00")') runF%yy_s, runF%mm_s, runF%dd_s, runF%hh_s; call goPr
    write (gol,'("end date of the simulation at:   ",i4,2("-",i2.2)," ",i2.2,":00")') runF%yy_e, runF%mm_e, runF%dd_e, runF%hh_e; call goPr
    write (gol,'("")'); call goPr

    !
    ! time structures
    !

    ! store in structure:
    runF%t_start = NewDate( year=runF%yy_s, month=runF%mm_s, day=runF%dd_s, hour=runF%hh_s )
    runF%t_end   = NewDate( year=runF%yy_e, month=runF%mm_e, day=runF%dd_e, hour=runF%hh_e )
    
    !
    ! done
    !
    
    ! ok
    status = 0

  end subroutine LE_Time_Init
  
  ! ***
  
  
  subroutine LE_Time_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Time_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0

  end subroutine LE_Time_Done

 
  ! ***
  

  subroutine local_time2(country_code, yy, mm, dd, hh, yyh, mmh, ddh, hhh, iday, status)

    !use LE_Logging, only : days
    use GO, only : calc_days_in_month

    ! --- const ------------
    character(len=*), parameter   ::  rname = mname//'/local_time2'

    ! --- local ------------
    ! computes local time, based on the country table
    ! hour range from 0 - 23
    character(len=3), intent(in) :: country_code
    integer, intent(in)          :: yy, mm, dd, hh
    integer, intent(out)         :: yyh, mmh, ddh, hhh, iday    
    integer, intent(out)         :: status

    ! set help date equal to input
    yyh = yy
    ddh   = dd
    mmh   = mm

    ! determine the shift from GMT base on the country code
    select case(country_code)
      
      case( 'GRL' )
        hhh = hh - 3
      case( 'ATL','ATX','GBR','INT','IRL','ISL','PRT', 'MAR', &
            'GGY','FRO','JEY','IMN','IRC','ENC','GRS')
        hhh = hh
      case( 'ALB','AUT','BAS','BEL','BIH','CHE','CZE','DEU','DNK','ESP','FRA','GIB','HRV','HUN','ITA', &
            'LIE','LUX','MED','MKD','MLT','NLD','NOA','NOR','NOS','POL','SUN','SVK','SVN','SWE','YUG', &
            'DZA', 'TUN', 'SMR', 'MCO', 'SJM', 'AND','NWS','KOS','MNE', 'SRB' )
        hhh = hh + 1
      case( 'BGR','BLS','CYP','EST','FIN','GRC','LTU','LVA','MDA','ROU','TUR','UKR', &
            'EGY', 'ISR', 'JOR', 'LBN', 'LBY', 'PSE', 'SYR', 'ALA' )
        hhh = hh + 2
      case( 'BLR','RUR','RUS','RUX', 'IRQ', 'KWT', 'SAU','BAR','KAR', 'PSG' )
        hhh = hh + 3
      case( 'ASI','ARM','AZE','GEO', 'IRN', 'CAS' )
        hhh = hh + 4
      case( 'KAZ','TJK','TKM','UZB' )
        hhh = hh + 5        
      case( 'KGZ' )
        hhh = hh + 6
      case default
        write (gol,'("Time zone not defined for country:", a)') country_code; call goErr       
        TRACEBACK; status=1; return        
    end select

    ! correct for summer time: one hour extra shift from GMT
    if (mm > 3 .AND. mm < 11) then
      select case(country_code)
        case( 'ALB','ASI','ARM','AUT','AZE','BAS','BEL','BGR','BIH','BLS','CHE','CYP','CZE','DEU', &
              'DNK','ESP','EST','FIN','FRA','GBR','GIB','GRC','HRV','HUN','IRL','ITA','LIE','LTU','LUX', &
              'LVA','MDA','MED','MKD','MLT','NLD','NOA','NOR','NOS','POL','PRT','ROU','RUR','RUS','RUX', &
              'SUN','SVK','SVN','SWE','TUR','UKR','YUG', &
              'GRL', 'IRN', 'ISR', 'KWT', 'LBN', 'MAR', 'PSE', 'SYR', &
              'GGY', 'FRO', 'JEY', 'IMN', 'SMR', 'MCO', 'ALA', &
              'SJM', 'IRC', 'AND', 'ENC', 'NWS', 'KOS', 'MNE', 'SRB' )
          hhh = hhh + 1
        case( 'ATL','ATX','BLR', 'GEO','INT','ISL','KAZ','KGZ','TJK','TKM','UZB', &
              'DZA', 'EGY', 'IRQ', 'JOR', 'LBY', 'SAU', 'TUN','CAS','BAR', 'GRS','PSG','KAR' )
          hhh = hhh
        case default
          write (gol,'("Daylight saving time defined for country:", a)') country_code; call goErr         
          TRACEBACK; status=1; return
      end select
    end if
    
    if (hhh > 23) then
       ! entered a new day
       hhh=hhh-24
       ddh=ddh+1
       ! also entered new month?
       !if (ddh > days(mmh)) then
       if ( ddh > calc_days_in_month('gregorian',yyh,mmh) ) then
          mmh = mmh + 1
          ddh = 1
          if (mmh > 12) then
            mmh=1
            yyh = yyh+1
          endif
       endif
    endif

    call day_of_the_week(yyh,mmh,ddh,iday)
    
    ! ok
    status = 0
    
  end subroutine

 
  ! ***
  

  subroutine day_of_the_week(yy,mm,dd,iday)

    !use LE_Logging, only : days
    use GO, only : calc_DayNumber

    ! input: year, month, day. Input/output: iday, the day of the week
    ! 1=monday, 2=tuesday, ......, 7=sunday (should actually be 1!)
    integer :: yy, mm, dd, iday, iday2

    !local
    integer :: i, nd

    ! now compute day-of-the-week
    ! used: 1-1-1990 is a monday (day-of-week=1)
    iday = calc_DayNumber( 'gregorian', yy, mm, dd )

    if (yy < 1990) then
    !new, AM 9-6-2010: count backwards from 1 jan 1990
        if (yy.eq.1989) then
         nd=(365-iday)/7
         iday=7-(365-iday-nd*7)
        else
           iday2=0
           do i=1989,yy,-1
            if (mod(i,4)==0) then
              iday2 = iday2 + 366
            else
              iday2 = iday2 + 365
            endif
           enddo
         nd=(iday2-iday)/7
         iday=7-(iday2-iday-nd*7)
        endif
         
         
    
       
   
    else
    !count foreward from 1990, original code
    do i=1990,yy-1
       if (mod(i,4)==0) then
          iday = iday + 366
       else
          iday = iday + 365
       endif
    enddo

    ! now we can compute the day number
    nd=iday/7
    iday=iday-nd*7
    endif 
    
    if (iday == 0) iday = 7
    if (iday < 1 .OR. iday > 7)  print *, 'wrong day number ',iday
  end subroutine

 
  ! ***
  

!  subroutine local_time(ix, iy, yy, mm, dd, hh, yyh, mmh, ddh, hhh, iday)
!
!    !use LE_Logging, only : days
!    use GO, only : calc_days_in_month
!
!    ! computes local time, based on the country table
!    ! hour range from 1 - 24
!    use dims, only : lst
!    integer :: ix, iy, yy, mm, dd, hh, yyh, mmh, ddh, hhh, iday
!
!    ! set years equal
!    yyh = yy
!
!    ! compute local date and time
!    hhh   = hh + lst(ix,iy)
!    ! correct for summer time: one hour extra shift from GMT
!    if (mm > 3 .AND. mm < 10) hhh = hhh + 1
!    ddh   = dd
!    mmh   = mm
!
!    if (hhh > 24) then
!       ! entered a new day
!       hhh=hhh-24
!       mmh=mmh+1
!       ! also entered new month?
!       !if (mmh > days(mmh)) then
!       if ( mmh > calc_days_in_month('gregorian',yyh,mmh) ) then
!          mmh = mmh + 1
!          if (mmh > 12) then
!             mmh=1
!             yyh=yyh+1
!          endif
!       endif
!    endif
!
!    ! now compute day-of-the-week
!    call day_of_the_week(yyh,mmh,ddh,iday)
!
!  end subroutine

 
  ! ***
  

  logical function date_le_date(iy1,im1,id1,ih1,iy2,im2,id2,ih2)

    ! true if date1 is smaller than or equal to date 2
    integer :: date1, date2
    integer :: iy1,im1,id1,ih1,iy2,im2,id2,ih2

    date_le_date = .false.
    if (iy1 < iy2) then
      date_le_date = .true.
    else if (iy1 == iy2) then
      date1 = ih1 + id1*100 + im1*10000
      date2 = ih2 + id2*100 + im2*10000
      if (date1 <= date2) date_le_date = .true.
    endif

  end function


end module LE_Time
