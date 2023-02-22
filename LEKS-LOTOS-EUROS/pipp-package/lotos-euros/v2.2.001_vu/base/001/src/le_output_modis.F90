!###############################################################################
!
! NAME
!
!   LE_Output_MODIS  -  put out simulations on MODIS pixels
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################
  
module le_output_modis

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif
  
  use LE_Output_Common, only : T_LE_Output_Common
  
  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
!  public  ::  LE_Output_MODIS_Init, LE_Output_MODIS_Done
!  public  ::  LE_Output_MODIS_GetMeas
!  public  ::  MeasUpdate_MODIS
!  public  ::  CalcMeas_MODIS
!  
  public  ::  maxp
  !public  ::  status_default, status_nodata, status_screened
  public  ::  bpos_nodata, bpos_validation, bpos_screened, bpos_analysed
  
  public  ::  T_Data_MODIS
  public  ::  T_Simulation_MODIS
  public  ::  T_LE_Output_MODIS
  public  ::  Init, Done
  public  ::  Setup
  public  ::  Simulate
  public  ::  PutOut
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'le_output_modis'
  
  ! maximum number of satellite platforms:
  integer, parameter    ::  maxplatform = 2
  
  ! character lines:
  integer, parameter    ::  len_line      = 512
  integer, parameter    ::  len_units     =  64
  integer, parameter    ::  len_long_name = 512
  
  ! maximum number of modis measurements within time step:
  integer, parameter    ::  maxp = 60000 ! 60e3

  !   
  ! Status number is sum of one or more numbers 2^n,
  ! for example   1+4+8
  ! To test if it is composed of a certain number, use:
  !   if ( iand(status,4) /= 0 )  ..
  !
  !integer, parameter            ::  status_default    = 0
  !integer, parameter            ::  status_outside    = 1
  !integer, parameter            ::  status_nodata     = 2
  !integer, parameter            ::  status_validation = 4
  !integer, parameter            ::  status_screened   = 8
  !character(len=*), parameter   ::  status_description = &
  !    'status flag, 0=default, +1=outside-domain, +2=no-data, +4=validation, +8=screened'
  !
  ! use bit-by-bit functions using positions:
  !
  !   integer, parameter  ::  bpos_outside     = 1
  !   integer, parameter  ::  bpos_nodata      = 2
  !   integer, parameter  ::  bpos_validation  = 3
  !
  ! Set and clear the bits using:
  !   status = 0
  !   status = IBSet( status, bpos_nodata )
  !   status = IBClr( status, bpos_nodata )
  ! Test if a bit is set:
  !   if ( BTest( status, bpos_nodata ) ) ...
  !
  integer, parameter            ::  bpos_outside    = 1
  integer, parameter            ::  bpos_nodata     = 2
  integer, parameter            ::  bpos_validation = 3
  integer, parameter            ::  bpos_screened   = 4
  integer, parameter            ::  bpos_analysed   = 5
  character(len=*), parameter   ::  status_description = &
      'status flag, 0=default, +1=outside-domain, +2=no-data, +4=validation, +8=screened, +16=analysed'

  
  ! --- types --------------------------------
    
  ! storage of MODIS data:
  type T_Data_MODIS
    ! archive directory:
    character(len=512)            ::  arch_dir
    ! platforms:
    integer                       ::  nplatform
    character(len=16)             ::  platforms(maxplatform)
    ! product names:
    character(len=16)             ::  products(maxplatform)
    ! list files
    character(len=len_line)       ::  list_file  (maxplatform)
    integer                       ::  list_fu    (maxplatform)
    logical                       ::  list_opened(maxplatform)
    integer                       ::  list_iline (maxplatform)
    character(len=len_line)       ::  list_fname (maxplatform)
    type(TDate)                   ::  list_t     (maxplatform)
    integer                       ::  list_jday  (maxplatform)
    ! data:
    integer                       ::  npix
    integer                       ::  pix_ipl(maxp)
    integer                       ::  pix_jday(maxp)
    real                          ::  pix_time(maxp)
    integer                       ::  pix_date(6,maxp)
    real                          ::  pix_Longitude(maxp)
    real                          ::  pix_Latitude (maxp)
    real                          ::  pix_Cloud_Fraction_Land(maxp)
    real                          ::  pix_Cloud_Fraction_Ocean(maxp)
    real                          ::  pix_Optical_Depth_Land_And_Ocean(maxp)
    real                          ::  pix_Optical_Depth_Ratio_Small_Land_And_Ocean(maxp)
    character(len=len_units)      ::  pix_Cloud_Fraction_Land__units
    character(len=len_units)      ::  pix_Cloud_Fraction_Ocean__units
    character(len=len_units)      ::  pix_Optical_Depth_Land_And_Ocean__units
    character(len=len_units)      ::  pix_Optical_Depth_Ratio_Small_Land_And_Ocean__units
    character(len=len_long_name)  ::  pix_Cloud_Fraction_Land__long_name
    character(len=len_long_name)  ::  pix_Cloud_Fraction_Ocean__long_name
    character(len=len_long_name)  ::  pix_Optical_Depth_Land_And_Ocean__long_name
    character(len=len_long_name)  ::  pix_Optical_Depth_Ratio_Small_Land_And_Ocean__long_name
    ! cell indices in model grid:
    integer                       ::  pix_ix(maxp)
    integer                       ::  pix_iy(maxp)
    ! assimilation flag:
    integer                       ::  pix_assim_status(maxp)
    ! times are computed as days since t0
    type(TDate)                   ::  t0
    ! error std.dev. fraction:
    real                          ::  sigma_frac
    ! analyse this set ?
    logical                       ::  analyse
  end type T_Data_MODIS
  
  
  ! ***
  
  
  ! simulated measurements:
  
  type T_Simulation_MODIS
    ! number of pixels:
    integer                 ::  npix
    ! simulated AOD:
    real                    ::  aod(maxp)
  end type T_Simulation_MODIS
    

  ! **
  

  ! output stuff
  
  type T_LE_Output_MODIS
    ! name for this file:
    character(len=16)           ::  typ
    character(len=16)           ::  name
    ! state name:
    character(len=16)           ::  state
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
    ! file opened ?
    logical                     ::  opened
    ! current time range:
    type(TDate)                 ::  tr(2)
    ! latest time:
    type(TDate)                 ::  tprev
    ! record counter:
    integer                     ::  irec
    ! file name:
    character(len=512)          ::  fname
    ! file handle:
    integer                     ::  ncid
    ! dimension handles:
    integer                     ::  dimid_pixel
    ! meta variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_time
    integer                     ::  varid_jday
    integer                     ::  varid_Cloud_Fraction_Land
    integer                     ::  varid_Cloud_Fraction_Ocean
    integer                     ::  varid_Optical_Depth_Land_And_Ocean
    integer                     ::  varid_Optical_Depth_Ratio_Small_Land_And_Ocean
    ! simulation:
    integer                     ::  varid_sim
  end type T_LE_Output_MODIS
  
  
  ! --- interfaces -------------------------
  
  interface Init
     module procedure dat_Init
     module procedure leo_Init
  end interface
  
  interface Done
     module procedure dat_Done
     module procedure leo_Done
  end interface
  
  interface Setup
     module procedure dat_Setup
     module procedure leo_Done
  end interface
  
  interface Simulate
     module procedure dat_Simulate
  end interface
  
  interface PutOut
     module procedure leo_PutOut
  end interface

  
  ! --- var ----------------------------------

!  

contains


  ! ======================================================================
  ! ===
  ! === input, update
  ! ===
  ! ======================================================================


  subroutine dat_Init( dat, rcfile, rckey, typ, status )

    use GO, only : TrcFile, Init, Done, ReadRc
    use GO, only : NewDate, AnyDate
    use GO, only : goSplitString
    use GO, only : goGetFU

    ! --- in/out --------------------------------

    type(T_Data_MODIS), intent(out)   ::  dat
    character(len=*), intent(in)      ::  rcfile
    character(len=*), intent(in)      ::  rckey
    character(len=*), intent(in)      ::  typ
    integer, intent(out)              ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/modis_Init'

    ! --- local ----------------------------------

    type(TRcFile)         ::  rcF
    character(len=64)     ::  key
    integer               ::  i
    logical               ::  exist

    ! --- begin ----------------------------------

    write (gol,'("init MODIS data ...")') ; call goPr

    ! ~~ settings from .rc file

    ! open rcfile:
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)

    ! archive location:
    call ReadRc( rcF, trim(rckey)//'.'//trim(typ)//'.archive', dat%arch_dir, status )
    IF_NOTOK_RETURN(status=1)

    ! space seperated list with platform names:
    call ReadRc( rcF, trim(rckey)//'.'//trim(typ)//'.platforms', key, status )
    IF_NOTOK_RETURN(status=1)
    ! split into platforms:
    call goSplitString( trim(key), dat%nplatform, dat%platforms, status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over platforms:
    do i = 1, dat%nplatform
      ! read product name:
      call ReadRc( rcF, trim(rckey)//'.'//trim(typ)//'.'//trim(dat%platforms(i))//'.product', dat%products(i), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! error std.dev. fraction:
    call ReadRc( rcF, trim(rckey)//'.'//trim(typ)//'.sigma_frac', dat%sigma_frac, status )
    IF_NOTOK_RETURN(status=1)
    
    ! analyse this set ?
    call ReadRc( rcF, trim(rckey)//'.'//trim(typ)//'.analyse', dat%analyse, status )
    IF_NOTOK_RETURN(status=1)
    
    ! close:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! ~~ init list files:
    
    ! loop over platforms:
    do i = 1, dat%nplatform
      ! full path:
      write (dat%list_file(i),'(a,"/",a,"/",a,".list")') &
               trim(dat%arch_dir), trim(dat%platforms(i)), trim(dat%products(i))
      ! check ...
      inquire( file=trim(dat%list_file(i)), exist=exist )
      if ( .not. exist ) then
        write (gol,'("list file not found: ",a)') trim(dat%list_file(i)); call goPr
        TRACEBACK; status=1; return
      end if
      ! free file unit:
      call goGetFU( dat%list_fu(i), status )
      IF_NOTOK_RETURN(status=1)
      ! open :
      open( dat%list_fu(i), file=trim(dat%list_file(i)), &
              status='old', form='formatted', iostat=status )
      if (status/=0) then
        write (gol,'("opening modis list file:")'); call goErr
        write (gol,'("  platform    : ",a," (",i2,")")') trim(dat%platforms(i)), i; call goErr
        write (gol,'("  list file   : ",a)') trim(dat%list_file(i)); call goErr
        TRACEBACK; status=1; return
      end if
      ! set flag:
      dat%list_opened(i) = .true.
      ! nothing read yet:
      dat%list_iline(i) = 0
      dat%list_t(i) = AnyDate()
    end do

    ! times are computed as days since 2000-1-1 00:00
    dat%t0 = NewDate( time6=(/2000,01,01,00,00,00/) )

    ! nothing stored yet ...
    dat%npix = 0

    ! ok
    status = 0

  end subroutine dat_Init


  ! ***


  ! read station locations etc

  subroutine dat_Done( dat, status )

    ! --- in/out --------------------------------

    type(T_Data_MODIS), intent(inout)   ::  dat
    integer, intent(out)                ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/dat_Done'

    ! --- local ----------------------------------

    integer                   ::  i

    ! --- begin ------------------------------
    
    ! loop over list files:
    do i = 1, dat%nplatform
      ! list file open ?
      if ( dat%list_opened(i) ) then
        ! close ...
        close( dat%list_fu(i), iostat=status )
        IF_NOTOK_RETURN(status=1)
        ! reset flag:
        dat%list_opened(i) = .false.
      end if
    end do

    ! ok
    status = 0

  end subroutine dat_Done


  ! ***


  ! read measurments for time interval (t1,t2] :

  subroutine dat_Setup( dat, t1, t2, status )

    use GO, only : TDate, NewDate, IncrDate, AnyDate, IsAnyDate, Get
    use GO, only : operator(>), operator(<=), operator(+), operator(-)
    use GO, only : rTotal, wrtgol
    use GO, only : goReadFromLine
    use Grid   , only : GetLocation
    use MODIS  , only : T_MODIS, Init, Done
    use LE_Grid, only : ugg
    
    ! --- in/out --------------------------------

    type(T_Data_MODIS), intent(inout)   ::  dat
    type(TDate), intent(in)             ::  t1, t2
    integer, intent(out)                ::  status

    ! --- const --------------------------------

    character(len=*), parameter   ::  rname = mname//'/dat_Setup'

    ! --- local ----------------------------------

    character(len=32)         ::  prompt
    integer                   ::  ipl
    integer                   ::  i, j
    character(len=len_line)   ::  line
    character(len=32)         ::  key
    integer                   ::  year, jday, hour, minu
    type(T_MODIS)             ::  tmodis
    character(len=512)        ::  fname
    logical                   ::  indomain
    logical                   ::  isfirst
    real(8)                   ::  sst_first
    type(TDate)               ::  t_first
    type(TDate)               ::  pix_t0, pix_t
    integer                   ::  time6(6)
    integer                   ::  ipix

    ! --- begin ----------------------------------
    
    ! clear buffer:
    dat%npix = 0
    
    ! loop over platforms:
    do ipl = 1, dat%nplatform
    
      ! set prompt:
      prompt = '    MODIS/'//trim(dat%platforms(ipl))//': '
      ! info:
      call wrtgol( trim(prompt)//'setup for ', t1, ' - ', t2 ); call goPr

      ! if list file is closed, no new data is available anymore ...
      if ( .not. dat%list_opened(ipl) ) cycle
      
      ! loop until time fits:
      do
      
        ! item available ?
        if ( .not. IsAnyDate(dat%list_t(ipl)) ) then

          ! too new ? then leave:
          if ( dat%list_t(ipl) > t2 ) then
            call wrtgol( trim(prompt)//'  keep orbit in buffer; valid for : ', dat%list_t(ipl) ); call goPr
            exit
          end if
      
          ! in interval (t1,t2] ?
          if ( (dat%list_t(ipl) > t1) .and. (dat%list_t(ipl) <= t2) ) then
            ! info ...
            write (gol,'(a,"  read orbit ...")') trim(prompt); call goPr
            ! create file name: <archive>/<platform>/<fname>
            write (fname,'(a,"/",a,"/",a)')  trim(dat%arch_dir), trim(dat%platforms(ipl)), trim(dat%list_fname(ipl))
            !! info ...
            !write (gol,'(a,"    filename     : ",a)') trim(prompt), trim(fname); call goPr
            ! read input file:
            call Init( tmodis, trim(fname), status )
            IF_NOTOK_RETURN(status=1)
            ! handle time:
            if ( trim(tmodis%Scan_Start_Time__units) == 'Seconds since 1993-1-1 00:00:00.0 0' ) then
              pix_t0 = NewDate( 1993, 01, 01, 00, 00, 00 )
            else
              write (gol,'("unsupported MODIS time unit : ",a)') trim(tmodis%Scan_Start_Time__units); call goErr
              TRACEBACK; status=1; return
            end if
            ! copy units:
            dat%pix_Cloud_Fraction_Land__units                      = tmodis%Cloud_Fraction_Land__units
            dat%pix_Cloud_Fraction_Ocean__units                     = tmodis%Cloud_Fraction_Ocean__units
            dat%pix_Optical_Depth_Land_And_Ocean__units             = tmodis%Optical_Depth_Land_And_Ocean__units
            dat%pix_Optical_Depth_Ratio_Small_Land_And_Ocean__units = tmodis%Optical_Depth_Ratio_Small_Land_And_Ocean__units
            ! copy long_name:
            dat%pix_Cloud_Fraction_Land__long_name                      = tmodis%Cloud_Fraction_Land__long_name
            dat%pix_Cloud_Fraction_Ocean__long_name                     = tmodis%Cloud_Fraction_Ocean__long_name
            dat%pix_Optical_Depth_Land_And_Ocean__long_name             = tmodis%Optical_Depth_Land_And_Ocean__long_name
            dat%pix_Optical_Depth_Ratio_Small_Land_And_Ocean__long_name = tmodis%Optical_Depth_Ratio_Small_Land_And_Ocean__long_name
            ! info ...
            write (gol,'(a,"    checking new pixels : ",i6)') trim(prompt), tmodis%nx*tmodis%ny; call goPr
            ! loop over pixels:
            isfirst = .true.
            do j = 1, tmodis%ny
              do i = 1, tmodis%nx
                ! no-data ? then skip:
                if ( tmodis%Optical_Depth_Land_And_Ocean(i,j) < 0.0 ) cycle
                ! check location:
                call ugg%InDomain( real(tmodis%Longitude(i,j)), real(tmodis%Latitude(i,j)), indomain, status )
                IF_NOTOK_RETURN(status=1)
                ! not in domain ? then skip:
                if ( .not. indomain ) cycle
                ! full ?
                if ( dat%npix == maxp ) then
                  write (gol,'("reached maximum number of pixels : ",i6)') maxp; call goErr
                  write (gol,'("please increase parameter `maxp` ...")'); call goErr
                  TRACEBACK; status=1; return
                end if
                ! increase counter:
                ipix = dat%npix + 1
                ! store platform index:
                dat%pix_ipl (ipix) = ipl
                ! store julian day:
                dat%pix_jday(ipix) = dat%list_jday(ipl)
                ! pre-compute expensive time stuf:
                if ( isfirst ) then
                  ! store:
                  sst_first = tmodis%Scan_Start_Time(i,j)
                  ! full time:
                  t_first = pix_t0 + IncrDate( sec=nint(sst_first) )
                  ! reset flag:
                  isfirst = .false.
                end if
                ! convert from 'seconds since ...' to year,month,etc
                pix_t = t_first + IncrDate( sec=nint(tmodis%Scan_Start_Time(i,j)-sst_first) )
                call Get( pix_t, time6=time6 )
                dat%pix_date(:,ipix) = time6
                ! convert to 'days since t0' :
                dat%pix_time(ipix) = rTotal( pix_t - dat%t0, 'day' )
                ! copy data:
                dat%pix_Longitude(ipix)                                = tmodis%Longitude(i,j)                               
                dat%pix_Latitude (ipix)                                = tmodis%Latitude (i,j)                               
                dat%pix_Cloud_Fraction_Land (ipix)                     = tmodis%Cloud_Fraction_Land (i,j)                    
                dat%pix_Cloud_Fraction_Ocean(ipix)                     = tmodis%Cloud_Fraction_Ocean(i,j)                    
                dat%pix_Optical_Depth_Land_And_Ocean(ipix)             = tmodis%Optical_Depth_Land_And_Ocean(i,j)            
                dat%pix_Optical_Depth_Ratio_Small_Land_And_Ocean(ipix) = tmodis%Optical_Depth_Ratio_Small_Land_And_Ocean(i,j)
                ! grid cell with pixel center:
                call ugg%GetLocation( dat%pix_Longitude(ipix), dat%pix_Latitude(ipix), &
                                      dat%pix_ix(ipix), dat%pix_iy(ipix), status )
                ! set flags:
                !dat%pix_assim_status(ipix) = status_default
                dat%pix_assim_status(ipix) = 0
                ! store counter:
                dat%npix = ipix
              end do
            end do
            ! info ...
            write (gol,'(a,"    last pixel   : ",i6)') trim(prompt), dat%npix; call goPr
            ! close:
            call Done( tmodis, status )
            IF_NOTOK_RETURN(status=1)
          else
            ! file too old ...
            call wrtgol( trim(prompt)//'  orbit in buffer too old : ', dat%list_t(ipl) ); call goPr
          end if

        end if
      
        ! increase counter:
        dat%list_iline(ipl) = dat%list_iline(ipl) + 1

        ! read next line:
        read (dat%list_fu(ipl),'(a)',iostat=status) line
        ! eof ?
        if ( status < 0 ) then
          ! close:
          close( unit=dat%list_fu(ipl), iostat=status )
          IF_NOTOK_RETURN(status=1)
          ! reset flag:
          dat%list_opened(ipl) = .false.
          ! dummy time ..
          dat%list_t(ipl) = AnyDate()
          ! leave list loop:
          exit
        end if
        ! error ...
        if (status/=0) then
          write (gol,'("reading line from list file : ")'); call goErr
          write (gol,'("  file   : ",a)') trim(dat%list_file(ipl)); call goErr
          write (gol,'("  line   : ",i6)') trim(line); call goErr
          TRACEBACK; status=1; return
        end if
        
        ! store for later usage:
        dat%list_fname(ipl) = trim(line)

        ! example file name:
        !   MYD04_L2.A2006001.0850.005.2007128013650.hdf
        
        ! extract fields:
        ! MYD04_L2
        call goReadFromLine( line, key, status, sep='.' )
        IF_NOTOK_RETURN(status=1)
        ! A2006001
        call goReadFromLine( line, key, status, sep='.' )
        IF_NOTOK_RETURN(status=1)
        read (key(2:8),'(i4,i3)',iostat=status) year, jday
        if (status/=0 ) then
          write (gol,'("reading year and jday from line : ",a)') trim(line); call goErr
          TRACEBACK; status=1; return
        end if
        ! 0850
        call goReadFromLine( line, key, status, sep='.' )
        IF_NOTOK_RETURN(status=1)
        read (key,'(2i2)',iostat=status) hour, minu
        if (status/=0 ) then
          write (gol,'("reading hour and minu from line : ",a)') trim(line); call goErr
          TRACEBACK; status=1; return
        end if
        
        ! set time; jday 1..365 should be converted to correct month/day :
        dat%list_jday(ipl) = jday
        dat%list_t   (ipl) = NewDate( year=year, month=1, day=jday, hour=hour, min=minu )

      end do  ! lines in file list

    end do  ! platforms

    ! ok
    status = 0

  end subroutine dat_Setup


  ! ======================================================================
  ! ===
  ! === simulation
  ! ===
  ! ======================================================================
  
  
  subroutine dat_Simulate( dat, lli, aod, sim, status )
  
    use Grid, only : TllGridInfo, Check
  
    ! --- in/out ---------------------------------
    
    type(T_Data_MODIS), intent(in)          ::  dat
    type(TllGridInfo), intent(in)           ::  lli
    real, intent(in)                        ::  aod(:,:)
    type(T_Simulation_MODIS), intent(out)   ::  sim
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/dat_Simulate'
    
    ! --- local ----------------------------------
    
    integer         ::  ipix
    
    ! --- begin ----------------------------------
    
    ! check ...
    call Check( lli, 'n', aod, status )
    IF_NOTOK_RETURN(status=1)
    
    ! reset counter:
    sim%npix = dat%npix

    ! loop over pixels:
    do ipix = 1, dat%npix
    
      ! copy aod value:
      sim%aod(ipix) = aod(dat%pix_ix(ipix),dat%pix_iy(ipix))

    end do  ! pixels
    
    ! ok
    status = 0
    
  end subroutine dat_Simulate


  ! ======================================================================
  ! ===
  ! === output
  ! ===
  ! ======================================================================


  subroutine leo_Init( leo, rcfile, rckey, typ, name, state, status )

    use GO     , only : AnyDate
    use LE_Output_Common, only : Init

    ! --- in/out --------------------------------

    type(T_LE_Output_MODIS), intent(out)    ::  leo
    character(len=*), intent(in)            ::  rcfile
    character(len=*), intent(in)            ::  rckey
    character(len=*), intent(in)            ::  typ
    character(len=*), intent(in)            ::  name
    character(len=*), intent(in)            ::  state
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_Init'

    ! --- local ---------------------------------

    ! --- begin ---------------------------------

    ! store file name keys:
    leo%typ   = typ
    leo%name  = name
    ! store state name:
    leo%state = trim(state)
    
    ! init common stuff:
    call Init( leo%com, rcfile, rckey, status )
    IF_NOTOK_RETURN(status=1)

    ! files not open yet:
    leo%opened = .false.

    ! no time range set yet:
    leo%tr(1) = AnyDate()
    leo%tr(2) = AnyDate()

    ! ok
    status = 0

  end subroutine leo_Init


  ! ***


  subroutine leo_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif
    use LE_Output_Common, only : Done

    ! --- in/out --------------------------------

    type(T_LE_Output_MODIS), intent(inout)      ::  leo
    integer, intent(out)                        ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_Done'

    ! --- begin ---------------------------------

    ! file opened ?
    if ( leo%opened ) then
       ! close:
#ifdef with_netcdf
       status = NF90_Close( leo%ncid )
       IF_NF90_NOTOK_RETURN(status=1)
#endif
       ! reset flag:
       leo%opened = .false.
    end if

    ! done with common stuff ...
    call Done( leo%com, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine leo_Done


  ! ***


  subroutine leo_PutOut( leo, t, dat, sim, status, without_data )

    use GO, only : TDate, IncrDate, Get, NewDate, AnyDate
    use GO, only : operator(+), operator(-), operator(<), operator(>), rTotal
    use GO, only : wrtgol, Pretty

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE

    use LE_Output_Common, only : PutOut_GlobalAttributes

    ! --- in/out --------------------------------

    type(T_LE_Output_MODIS), intent(inout)    ::  leo
    type(TDate), intent(in)                   ::  t
    type(T_Data_MODIS)                        ::  dat
    type(T_Simulation_MODIS)                  ::  sim
    integer, intent(out)                      ::  status
    
    logical, intent(in), optional             ::  without_data

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_PutOut'

    ! --- local ---------------------------------

    integer               ::  time6(6)

    integer               ::  cmode
    integer               ::  varid
    character(len=512)    ::  units
    
    logical               ::  with_data

    ! --- begin ---------------------------------
    
    ! put out standard data too ?
    with_data = .true.
    if ( present(without_data) ) with_data = .not. without_data

    ! check ...
    if ( dat%npix /= sim%npix ) then
      write (gol,'("found different number of pixels in data and simulation:")'); call goErr
      write (gol,'("  data       : ",i6)') dat%npix; call goErr
      write (gol,'("  simulated  : ",i6)') sim%npix; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no measurements ? then leave
    if ( dat%npix < 1 ) then
       call wrtgol('    no MODIS data for time : ',t); call goPr
       status=0; return
    end if

#ifndef with_netcdf
    write (gol,'("LE Output MODIS module requires compilation with NetCDF")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! ** open new file ?

    ! extract time fields:
    call Get( t, time6=time6 )
    
    ! current time not in time range ?
    if ( (t < leo%tr(1)) .or. (leo%tr(2) < t) ) then

      ! file opened ?
      if ( leo%opened ) then
        ! close:
#ifdef with_netcdf    
        status = NF90_Close( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif
        ! reset flag:
        leo%opened = .false.
      end if

      ! set time range for this day:
      leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
      leo%tr(2) = leo%tr(1) + IncrDate( day=1 )

      ! new file name:
      write (leo%fname,'(a,a,"_",a,"_",a,"_",i4.4,2i2.2)') &
           trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
           trim(leo%name), time6(1:3)
      if ( len_trim(leo%state) > 0 ) write (leo%fname,'(a,"_",a)') trim(leo%fname), trim(leo%state)
      write (leo%fname,'(a,".nc")') trim(leo%fname)

      ! set creation mode flag:
      cmode = NF90_NOCLOBBER     ! do not overwrite existing files

      ! create file:
      status = NF90_Create( leo%fname, cmode, leo%ncid )
      if ( status /= 0 ) then
         write (gol,'("creating file :")'); call goErr
         write (gol,'("  ",a)') trim(leo%fname); call goErr
         TRACEBACK; status=1; return
      end if

      ! reset flag:
      leo%opened = .true.

      ! write global attributes:
      call PutOut_GlobalAttributes( leo%com, leo%ncid, status )
      IF_NOTOK_RETURN(status=1)

      ! define dimensions:
      status = NF90_Def_Dim( leo%ncid, 'pixel', NF90_UNLIMITED, leo%dimid_pixel )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variables:

      if ( with_data ) then
        status = NF90_Def_Var( leo%ncid, 'lon', NF90_REAL, (/leo%dimid_pixel/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'standard_name', 'longitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'longitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_lon = varid
      end if

      if ( with_data ) then
        status = NF90_Def_Var( leo%ncid, 'lat', NF90_REAL, (/leo%dimid_pixel/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'standard_name', 'latitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'latitude' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_lat = varid
      end if

      if ( with_data ) then
        status = NF90_Def_Var( leo%ncid, 'time', NF90_REAL, (/leo%dimid_pixel/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'time' )
        IF_NF90_NOTOK_RETURN(status=1)
        write (units,'("days since ",a)') trim(Pretty(dat%t0))
        status = NF90_Put_Att( leo%ncid, varid, 'units', trim(units) )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'calendar', 'standard' )
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_time = varid
      end if

      ! julian day:
      if ( with_data ) then
        status = NF90_Def_Var( leo%ncid, 'jday', NF90_INT, (/leo%dimid_pixel/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'julian day' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'units', '1=jan 1, 2=jan 2, ...' )
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_jday = varid
      end if

      ! data:
      if ( with_data ) then
        status = NF90_Def_Var( leo%ncid, 'Cloud_Fraction_Land', NF90_REAL, (/leo%dimid_pixel/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'long_name', dat%pix_Cloud_Fraction_Land__long_name )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'units', dat%pix_Cloud_Fraction_Land__units )
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_Cloud_Fraction_Land = varid
      end if

      ! data:
      if ( with_data ) then
        status = NF90_Def_Var( leo%ncid, 'Cloud_Fraction_Ocean', NF90_REAL, (/leo%dimid_pixel/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'long_name', dat%pix_Cloud_Fraction_Ocean__long_name )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'units', dat%pix_Cloud_Fraction_Ocean__units )
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_Cloud_Fraction_Ocean = varid
      end if

      ! data:
      if ( with_data ) then
        status = NF90_Def_Var( leo%ncid, 'Optical_Depth_Land_And_Ocean', NF90_REAL, (/leo%dimid_pixel/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'long_name', dat%pix_Optical_Depth_Land_And_Ocean__long_name )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'units', dat%pix_Optical_Depth_Land_And_Ocean__units )
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_Optical_Depth_Land_And_Ocean = varid
      end if

      ! data:
      if ( with_data ) then
        status = NF90_Def_Var( leo%ncid, 'Optical_Depth_Ratio_Small_Land_And_Ocean', NF90_REAL, (/leo%dimid_pixel/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'long_name', dat%pix_Optical_Depth_Ratio_Small_Land_And_Ocean__long_name )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Att( leo%ncid, varid, 'units', dat%pix_Optical_Depth_Ratio_Small_Land_And_Ocean__units )
        IF_NF90_NOTOK_RETURN(status=1)
        leo%varid_Optical_Depth_Ratio_Small_Land_And_Ocean = varid
      end if

      ! simulation:
      status = NF90_Def_Var( leo%ncid, 'aod', NF90_REAL, (/leo%dimid_pixel/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( leo%ncid, varid, 'long_name', 'aerosol optical depth' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( leo%ncid, varid, 'units', 'm' )
      IF_NF90_NOTOK_RETURN(status=1)
      leo%varid_sim = varid

      ! end defintion mode:

#ifdef with_netcdf    
      status = NF90_EndDef( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#endif

      ! no records written yet:
      leo%irec = 0
      leo%tprev = AnyDate()
    
    end if   ! open new file
    
    ! **
    
    ! add pixels if not done before:
    if ( t > leo%tprev ) then
    
      ! put out data ?
      if ( with_data ) then

        ! add longitudes:
        status = NF90_Put_Var( leo%ncid, leo%varid_lon, dat%pix_longitude(1:dat%npix), &
                                 start=(/leo%irec+1/), count=(/dat%npix/) )
        IF_NF90_NOTOK_RETURN(status=1)

        ! add latitudes:
        status = NF90_Put_Var( leo%ncid, leo%varid_lat, dat%pix_latitude(1:dat%npix), &
                                 start=(/leo%irec+1/), count=(/dat%npix/) )
        IF_NF90_NOTOK_RETURN(status=1)

        ! add time since t0
        status = NF90_Put_Var( leo%ncid, leo%varid_time, dat%pix_time(1:dat%npix), &
                                 start=(/leo%irec+1/), count=(/dat%npix/) )
        IF_NF90_NOTOK_RETURN(status=1)

        ! julian day (used in file name)
        status = NF90_Put_Var( leo%ncid, leo%varid_jday, dat%pix_jday(1:dat%npix), &
                                 start=(/leo%irec+1/), count=(/dat%npix/) )
        IF_NF90_NOTOK_RETURN(status=1)

        ! add data:
        status = NF90_Put_Var( leo%ncid, leo%varid_Cloud_Fraction_Land, &
                                 dat%pix_Cloud_Fraction_Land(1:dat%npix), &
                                 start=(/leo%irec+1/), count=(/dat%npix/) )
        IF_NF90_NOTOK_RETURN(status=1)

        ! add data:
        status = NF90_Put_Var( leo%ncid, leo%varid_Cloud_Fraction_Ocean, &
                                 dat%pix_Cloud_Fraction_Ocean(1:dat%npix), &
                                 start=(/leo%irec+1/), count=(/dat%npix/) )
        IF_NF90_NOTOK_RETURN(status=1)

        ! add data:
        status = NF90_Put_Var( leo%ncid, leo%varid_Optical_Depth_Land_And_Ocean, &
                                 dat%pix_Optical_Depth_Land_And_Ocean(1:dat%npix), &
                                 start=(/leo%irec+1/), count=(/dat%npix/) )
        IF_NF90_NOTOK_RETURN(status=1)

        ! add data:
        status = NF90_Put_Var( leo%ncid, leo%varid_Optical_Depth_Ratio_Small_Land_And_Ocean, &
                                 dat%pix_Optical_Depth_Ratio_Small_Land_And_Ocean(1:dat%npix), &
                                 start=(/leo%irec+1/), count=(/dat%npix/) )
        IF_NF90_NOTOK_RETURN(status=1)

      end if  ! with data
      
      ! add simulation:
      status = NF90_Put_Var( leo%ncid, leo%varid_sim, sim%aod(1:sim%npix), &
                               start=(/leo%irec+1/), count=(/sim%npix/) )
      IF_NF90_NOTOK_RETURN(status=1)

      ! increase counter:
      leo%irec = leo%irec + dat%npix

      ! store this time:
      leo%tprev = t
      
    end if  ! new time

    ! ok
    status = 0

  end subroutine leo_PutOut


end module LE_Output_MODIS
