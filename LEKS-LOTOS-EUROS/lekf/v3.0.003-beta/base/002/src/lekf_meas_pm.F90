!###############################################################################
!
! NAME
!
!   kf_meas_pm  -  interface to daily pm10 measurements
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
!###############################################################################


module LEKF_meas_pm

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  use NetCDF, only : NF90_StrError, NF90_NOERR
  
  use LE_Output_Common, only : T_LE_Output_Common
  
  use LEKF_State, only : maxmeas
  use LEKF_Meas , only : TGainBoundInfo
  
  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  pm25_rho
  public  ::  InitMeas_PM, DoneMeas_PM
  public  ::  GetMeas_PM, MeasUpdate_PM
  public  ::  UpdateState_PM
 
  public  ::  meas_output_pm
  public  ::  Init, Done, PutOut
  
  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'kf_meas_pm'
  
  ! the number of measurements in the EMEP file
  integer, parameter :: emepmeas = 84

  ! name lengths:
  integer, parameter            ::  codelen = 7
  integer, parameter            ::  namelen = 20
  
  ! filter states to be written:
  integer, parameter            ::  nfs = 6
  character(len=6), parameter   ::  fsname(nfs) = &
        (/'y     ','v     ','status','xb    ','x     ','s     '/)
  
  ! filter moments:
  integer, parameter            ::  nfm = 2
  
  ! status flags:
  integer, parameter            ::  status_validation = 1
  integer, parameter            ::  status_nan        = 2
  integer, parameter            ::  status_screened   = 4
  character(len=*), parameter   ::  status_description = &
      'status flag, 0=default, +1=validation, +2=no-data, +4=screened'  

  ! --- types ----------------------------------
   
  type measurement_pm
    ! the name and code of the station
    character(len=codelen)  :: code
    character(len=namelen)  :: name
    ! its coordinates
    real              :: lon, lat, alt
    ! its position in the hourly arrays:
    integer           :: irec
    ! flag if station is in current domain
    logical           :: indomain
    ! the grid indices in the 4D array
    integer           :: ix, iy, iz, is
    ! its position in the statevector
    integer           :: inds(0:7)
    ! the height of the inlet
    real              :: href
    ! radius of influence to be applied (km.)7
    real              :: rho
    ! analyse or for validation only ?
    ! data accepted or nan etc ?
    logical           ::  analyse
    !logical           ::  accepted
    integer           ::  status
    ! about the measurement: the value + standard deviation
    real              ::  c, sd
  end type measurement_pm

  ! output stuff
  
  type meas_output_pm
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
    ! file opened ?
    logical                     ::  opened
    ! current time range:
    type(TDate)                 ::  tr(2)
    ! latest time:
    type(TDate)                 ::  tprev
    ! time record counter:
    integer                     ::  itrec
    ! file name:
    character(len=512)          ::  fname
    ! file handle:
    integer                     ::  ncid
    ! dimension handles:
    integer                     ::  dimid_namelen
    integer                     ::  dimid_codelen
    integer                     ::  dimid_meas
!    integer                     ::  dimid_lon
!    integer                     ::  dimid_lat
    integer                     ::  dimid_time
    integer                     ::  dimid_datelen
    ! dimension variables:
    integer                     ::  varid_time
    integer                     ::  varid_date
    ! station variables:
    integer                     ::  varid_stn_name
    integer                     ::  varid_stn_code
    integer                     ::  varid_stn_lon
    integer                     ::  varid_stn_lat
    integer                     ::  varid_stn_ix
    integer                     ::  varid_stn_iy
    integer                     ::  varid_stn_alt
    integer                     ::  varid_stn_ana
!    integer                     ::  varid_lon
!    integer                     ::  varid_lat
    ! tracer variables:
    real                        ::  unitconv
    integer                     ::  varid(nfs,nfm)
  end type meas_output_pm


  ! --- interfaces -------------------------
  
  interface Init
    module procedure meas_output_pm_Init
  end interface
  
  interface Done
    module procedure meas_output_pm_Done
  end interface

  interface PutOut
    module procedure meas_output_pm_PutOut
  end interface

  
  ! --- var --------------------------------
  
  !
  ! pm25s = so4a + no3a + nh4a + pm25 + bc + dust_f + na_f
  !         \_sec.inorg.aer._/   \_______prim.aer._______/
  !
  integer, parameter   ::  n_pm25s = 7
  integer              ::  ii_pm25s(n_pm25s) 

  ! analyse this set ?
  logical                   ::  analyse

  ! radius of influence:
  real                      ::  pm25_rho
  type(TGainBoundInfo)      ::  gbi
  
  ! declaration of the meas set
  type(measurement_pm)  ::  meas_pm(maxmeas)
  
  ! actual number read:
  integer                   ::  nmeas
  
!  integer   ::  meas_record
!  integer   ::  out_unit
!  integer   ::  out_record

  ! single tracer only, thus define globaly
  integer     ::  icomp  

  ! time range:
  type(TDate)           ::  tr(2)

  ! file name and unit:
  character(len=256)    ::  fname
  integer               ::  fu


contains


  ! ======================================================================
  ! ===
  ! === input, update
  ! ===
  ! ======================================================================
  
  
  ! opens the measurement file (if present)
  ! end sets some parameters

  subroutine InitMeas_PM( rcfile, t, status )

    use GO, only : TrcFile, Init, Done, ReadRc
    use GO, only : goGetFU
    use GO, only : TDate, NewDate, Get
    
    use dims      , only : nx, ny, runF, locF
    use constants , only : cum_days
    use units     , only : u_tmp
    !use sysdep
    use utils     , only : indomain
    use indices   , only : i_so4a,i_no3a,i_nh4a,i_pm25,i_bc,i_dust_f,i_na_f

!    use LEKF_state  , only :  m_s
    use LEKF_state  , only :  substate_m
    !use LEKF_meas, only : meas, nmeas, emepmeas
    !use LEKF_meas, only : meas_unit, u_meas_min
    use LEKF_meas, only : Init

    ! --- in/out --------------------------------
    
    character(len=*), intent(in)    ::  rcfile
    type(TDate), intent(in)         ::  t
    integer, intent(out)            ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/initmeas_pm'
    
    ! --- local ----------------------------------
    
    type(TRcFile)         ::  rcF
    character(len=256)    ::  fdir
    integer               ::  yy, mm
    character(len=6)      ::  period
    logical               ::  exist
    integer               ::  istn, nstn
    character(len=16)     ::  name
    real                  ::  lon, lat, alt
    integer               ::  irec
    integer               ::  ix, iy
    logical               ::  in_domain
    integer               ::  recl
    real(4)               ::  cc(emepmeas)
    
    ! --- begin ----------------------------------
    
    write (gol,'("KF:     setup pm measurements ...")'); call goPr
    
    !
    ! ~~ settings
    !
    
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.pm25.analyse', analyse, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.pm25.dir', fdir, status )
    IF_NOTOK_RETURN(status=1)

    call ReadRc( rcF, 'kf.meas.pm25.rho', pm25_rho, status )
    IF_NOTOK_RETURN(status=1)

    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! set name and start date for binary files with measurements:
    call Get( t, year=yy, month=mm )
    if ( mm < 7 ) then
      period = 'winter'
      tr(1) = NewDate( yy, 02, 15, 01, 00 )  ! feb  15  01:00
      tr(2) = NewDate( yy, 03, 16, 00, 00 )  !  plus 29 days
      ! NOTE: at least one or maybe 2 hours extra in file!
    else
      period = 'summer'
      tr(1) = NewDate( yy, 07, 15, 01, 00 )  ! july 15  01:00
      tr(2) = NewDate( yy, 08, 16, 00, 00 )  !  plus 32 days
    end if

    ! gain bounds:
    if ( analyse ) then
      call Init( gbi, pm25_rho, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
   
    !
    ! ~~ tracers
    !
    
    ! set tracer to be pm25, note that it is in fact a sum of several components:
    icomp = i_pm25
    
    ! components of pm25 sum:    
    ii_pm25s = (/i_so4a,i_no3a,i_nh4a,i_pm25,i_bc,i_dust_f,i_na_f/)
    

    !
    ! ~~ station locations
    !
    
    ! nothing stored yet ...
    nmeas = 0
    
    write (gol,'("KF:     stn name lon lat alt irec indomain")'); call goPr

    ! name with station locations:
    write (fname,'(a,"/statlist.txt")') trim(fdir)
    
    ! check ...
    inquire( file=fname, exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found:")'); call goErr
      write (gol,'("  ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! select file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! adhoc, to avoid interfernece with LE units ..
    fu = fu + 456
    
    ! open:
    open( unit=fu, file=trim(fname), status='old', &
            form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("from opening file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if

    ! read the number of measurements
    read (fu,*,iostat=status) nstn
    if (status/=0) then
      write (gol,'("from reading number of stations from file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if

    ! read all locations
    do istn = 1, nstn
    
      ! read line:
      read(fu,*,iostat=status) name, lon, lat, alt, irec
      if (status/=0) then
        write (gol,'("from reading station line:")'); call goErr
        write (gol,'("  file name : ",a)') trim(fname); call goErr
        write (gol,'("  station   : ",i6)') istn; call goErr
        write (gol,'("  iostat    : ",i6)') status; call goErr
        TRACEBACK; status=1; return
      end if

      ! check if meas is in domain;
      ! return cell indices (ix,iy) and flag:
      call InDomain( lon, lat, runF%dlon, runF%dlat, runF%southb, runF%westb, &
                      ix, iy, in_domain )

      ! info ...
      write (gol,'("KF:     ",i3," ",a,3f8.2,i4," ",l2)') &
              istn, trim(name), lon, lat, alt, irec,&
              in_domain; call goPr

      ! not in domain ? then skip ...
      if ( .not. in_domain ) cycle
      
      ! above 700 m ? then skip ...
      if ( alt > 700.0 ) cycle

      ! increase counter:
      nmeas = nmeas + 1

      ! check ...
      if ( nmeas > maxmeas ) then
        write (gol,'("total number of stations in data sets exceeds maxmeas = ",i6)') maxmeas; call goErr
        TRACEBACK; status=1; return
      end if

      ! store:
      ! the name and code of the station
      meas_pm(nmeas)%code      = trim(name)
      meas_pm(nmeas)%name      = trim(name)    ! no long name available, use code ...
      ! its coordinates
      meas_pm(nmeas)%lon       = lon
      meas_pm(nmeas)%lat       = lat
      meas_pm(nmeas)%alt       = alt
      ! position in data record:
      meas_pm(nmeas)%irec      = irec
      ! flag if station is in current domain
      meas_pm(nmeas)%indomain  = in_domain
      ! the grid indices in the 4D array
      meas_pm(nmeas)%ix        = ix
      meas_pm(nmeas)%iy        = iy
      meas_pm(nmeas)%iz        = 0  ! surface concentrations
      meas_pm(nmeas)%is        = icomp
      
      ! the index of the measurement in the measurement part 
      ! of the state vector:
!      meas_pm(nmeas)%ind        = m_s-1 + nmeas
      meas_pm(nmeas)%inds(0:1)   = (/ substate_m, nmeas /)

      ! the height of the inlet
      meas_pm(nmeas)%href      = 0.0            ! dummy

      ! radius of influence to be applied (km.)
      meas_pm(nmeas)%rho       = pm25_rho   ! km

      ! analyse this site ?
      meas_pm(nmeas)%analyse   = analyse

      !! not accepted by default ...
      !meas_pm(nmeas)%accepted  = .false.
      ! accepted by default ...
      meas_pm(nmeas)%status = 0

      ! no data yet ..
      meas_pm(nmeas)%c         = 0.0
      meas_pm(nmeas)%sd        = 0.0

    end do

    ! close:
    close( unit=fu, iostat=status )
    if (status/=0) then
      write (gol,'("from closing file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    

    !
    ! ~~ data file
    !
    
    ! name with data:
    write (fname,'(a,"/",a,"_hourly_pm25.dat")') trim(fdir), trim(period)
    
    ! check ...
    inquire( file=fname, exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found:")'); call goErr
      write (gol,'("  ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! get record length 
    inquire( iolength=recl ) cc

    ! open the meas file
    !open(meas_unit,file='/user4/lotos/measurements/meas/grnd/emep/emep97.dat', &
    !     form='unformatted', &
    !     recl=emepmeas*reclen_fac_bin_io ,access='direct',status='old')
    !open(meas_unit,file='/linuxMA5a/lotos/input/eumetsat2/summer_hourly_pm25.dat', &
    !     form='unformatted', recl=emepmeas*reclen_fac_bin_io ,access='direct',status='old') 

    ! open file:
    open( unit=fu, file=trim(fname), status='old', &
            form='unformatted', access='direct', recl=recl, &
            iostat=status )
    if (status/=0) then
      write (gol,'("from opening file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
    !out_unit = u_meas_min+2
    !open(out_unit,file=trim(locF%outpath)//'pm_stat.dat', &
    !     form='unformatted', recl=1*reclen_fac_bin_io ,access='direct',status='new') 

    !meas_record = (cum_days(mm) + dd -1)*24 +hh-1- (cumdays(2)+14)*24
    !out_record = 0

  end subroutine InitMeas_PM


  ! ***
  

  ! read station locations etc

  subroutine DoneMeas_PM( status )

    use LEKF_meas, only : Done

    ! --- in/out --------------------------------
    
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/DoneMeas_PM'
    
    ! --- begin ------------------------------
    
    ! close data file:
    close( unit=fu, iostat=status )
    if (status/=0) then
      write (gol,'("from closing file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! done with gain bounds:
    if ( analyse ) then
      call Done( gbi, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0

  end subroutine DoneMeas_PM
  
  
  ! ***
  
  
  subroutine GetMeas_pm( t, status )

#ifdef pgf90
    use GO
#else  
    use GO, only : TDate, Get, operator(<), operator(-), rTotal, wrtgol
#endif

    !use LEKF_meas, only : meas, nmeas, emepmeas
    !use LEKF_meas, only : meas_unit
    
    ! --- in/out ----------------------------
    
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GetMeas_pm'
    
    ! --- local -----------------------------
    
    real      ::  dhr
    integer   ::  meas_rec
    real(4)   ::  cc(emepmeas)
    integer   ::  i, ind
    
    ! --- begin -----------------------------

    ! fill with nan:
    do i = 1, nmeas
      meas_pm(i)%c  = -999.9
      meas_pm(i)%sd = -999.9
      !meas_pm(i)%accepted = .false.
      ! no data by default ...
      meas_pm(i)%status = status_nan
    end do

    ! used for validation only ?
    do i = 1, nmeas
      if ( .not. meas_pm(i)%analyse ) then
        meas_pm(i)%status = meas_pm(i)%status + status_validation
      end if
    end do

    ! check ...
    if ( t < tr(1) ) then
      call wrtgol( 'KF:     pm25 not available before ', tr(1) ); call goPr
      status=0; return
    end if 
    if ( tr(2) < t ) then
      call wrtgol( 'KF:     pm25 not available after ', tr(2) ); call goPr
      status=0; return
    end if 

    ! hourly records:
    dhr = rTotal( t - tr(1), 'hour' )
    if ( dhr /= nint(dhr) ) then
      call wrtgol( 'KF:     pm25 only hourly, not for ', t ); call goPr
      status=0; return
    end if
    
    ! hourly record number:
    meas_rec = nint(dhr) + 1
    
    ! info ...    
    call wrtgol( 'KF:     load pm25 for ', t ); call goPr
    write (gol,'("KF:       record ",i6)') meas_rec; call goPr

    ! read record:
    read (fu,rec=meas_rec,iostat=status) cc
    if (status/=0) then
      write (gol,'("from reading measurement record:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname); call goErr
      write (gol,'("  file unit : ",i6)') fu; call goErr
      write (gol,'("  record    : ",i6)') meas_rec; call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if

    ! set meas:
    do i = 1, nmeas

      ! index of this station in measurement array:
      ind = meas_pm(i)%irec

      ! copy concentration:
      meas_pm(i)%c = cc(ind)
      
      ! valid number ?
      !meas_pm(i)%accepted = meas_pm(i)%c > 0.0
      if ( meas_pm(i)%c > 0.0 ) then
        ! reset status flag:
        meas_pm(i)%status = meas_pm(i)%status - status_nan
        ! compute the s.d.:
        !  20% of the concentration,
        !   minimum of 1,
        !   no maximum to prevent amplification of outlyers:
        !meas_pm(i)%sd = min(max( 1.0, 0.2*meas_pm(i)%c ), 7.5 )
        meas_pm(i)%sd = max( 1.0, 0.20 * meas_pm(i)%c )
      else
        ! dummy error standard deviation:
        meas_pm(i)%sd = -999.0
      end if

      !! convert to ug/m3
      !meas_pm(i)%sd = meas_pm(i)%sd 
      !meas_pm(i)%c  = meas_pm(i)%c 

      ! print *, 'read meas #',i,' value = ', c(ind)

    end do

    ! ok
    status = 0

  end subroutine GetMeas_pm


  ! ***
  

  subroutine UpdateState_PM( xm, cg, status )

    use LEKF_state , only : nx, ny, nspec, maxmeas

    ! --- in/out --------------------------------

    real, intent(inout)             ::  xm(maxmeas)
    real, intent(in)                ::  cg(nx,ny,nspec)
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/UpdateStates_PM'
    
    ! --- local ----------------------------------
    
    integer         ::  i, k
    integer         ::  ix, iy, iz, itr
    
    ! --- begin ------------------------------
    
    write (gol,'("KF:     fill measurements in state")'); call goPr
    
    ! loop over measurements:
    do i = 1, nmeas
    
      ! cell indices:
      ix = meas_pm(i)%ix
      iy = meas_pm(i)%iy
      
      ! lowest layer:
      iz = 1
      
      ! init to zero:
      xm(i) = 0.0

      ! add pm25 components:
      do k = 1, n_pm25s
        ! tracer index:
        itr = ii_pm25s(k)
        ! add:
        xm(i) = xm(i) + cg(ix,iy,itr)
      end do
    
    end do  ! nmeas
    
    ! ok
    status = 0

  end subroutine UpdateState_PM
  
  
  ! ***
  
  
  ! implementation of measurement update
  ! routine expects a measurement list +
  ! a list specifying the corresponding element nr.
  ! of the statevector.

  subroutine MeasUpdate_PM( t, nmodes, status )

    use GO, only : TDate, Midnight
    use LEKF_state, only : maxmodes
!    use LEKF_state, only : nstate
    use LEKF_meas , only : meas_update
    use LEKF_meas , only : meas_update_scalar
!    use LEKF_meas , only : meas_update_profile

    ! --- in/out -----------------------------

    type(TDate), intent(in)   ::  t    
    integer, intent(in)       ::  nmodes     ! actual number of modes
    integer, intent(out)      ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MeasUpdate_PM'

    ! measurements are 1D    
    integer, parameter  ::  prf_nz = 1

    ! --- local ------------------------------
    
    integer       ::  i
!    integer       ::  ind
    !real          ::  z1(nmeas)
    real          ::  c, sd
    integer       ::  ix, iy
!    integer       ::  ihist
!    real          ::  rho

    real          ::  prf_yr(prf_nz)
!    real          ::  prf_ya(prf_nz)
!    real          ::  prf_A(prf_nz,prf_nz)
    real          ::  prf_R(prf_nz,prf_nz)
    integer       ::  inds(0:7,prf_nz)
          
    ! --- begin ------------------------------
    
    ! check ...
    if ( nmodes > maxmodes ) then
      write (gol,'("KF:   nmodes > maxmodes : ",2i6)') nmodes, maxmodes; call goPr
      TRACEBACK; status=1; return
    end if

    ! no measurements ? then leave
    if ( nmeas < 1 ) then
      write (gol,'("KF:   no measurement updating: no measurements available")'); call goPr
      status=0; return
    end if

    ! info ...
    write (gol,'("KF:   measurement update PM25")'); call goPr
    
    ! loop over measurements:
    do i = 1, nmeas

             !write (gol,'("KF:   aaa ",2i6)') i, meas_pm(i)%status; call goPr

      ! skip if not accepted (validation, no data):
      !if ( .not. meas_pm(i)%accepted ) cycle
      if ( meas_pm(i)%status /= 0 ) cycle

      ! hoezee!
      write (gol,'("KF:   analyse ",i3)') i; call goPr

      ! extract:
      c    = meas_pm(i)%c
      sd   = meas_pm(i)%sd
      ix   = meas_pm(i)%ix
      iy   = meas_pm(i)%iy
!
!      ! ihist       time level under consideration 
!      !             (determines together with rho the influence area)
!      ! rho         influence radius
!      ihist = 4  ! <-- why 4 ?
!      rho   = meas_pm(i)%rho
!
!      ! check ...
!      if ( ind > nstate ) then
!        write (gol,'("KF:   ind > nstate : ",2i6)') ind, nstate; call goPr
!        TRACEBACK; status=1; return
!      end if

      ! dummy profile ...
      prf_yr = c      ! retrieved
!      prf_ya = 0.0    ! a-priori
!      prf_A  = 1.0    ! kernel
      prf_R  = sd**2  ! observation error covariance

!      ! interpolation from state:
!      inds = ind
      inds(:,1) = meas_pm(i)%inds

      ! analyse profile:
!      if ( prf_nz == 1 ) then
!        call meas_update_scalar( prf_yr(1), prf_ya(1), prf_A(1,1), prf_R(1,1), &
!                                  rho, ix, iy, inds(:,1), nmodes, ihist, &
!                                  status )
!        IF_NOTOK_RETURN(status=1)
!      else
!        call meas_update_profile( prf_nz, prf_yr, prf_ya, prf_A, prf_R, &
!                                  rho, ix, iy, inds, nmodes, ihist, &
!                                  status )
!        IF_NOTOK_RETURN(status=1)
!        stop 'meas_update_profile not supported'
!      end if
!
!      call meas_update( prf_nz, prf_yr, prf_R, inds, &
!                          ix, iy, gbi, &
!                          nmodes, status )
!      IF_NOTOK_RETURN(status=1)
!
      call meas_update_scalar( prf_yr(1), prf_R(1,1), inds(:,1), &
                                ix, iy, gbi, &
                                nmodes, status )
      if ( status == -1 ) then
        meas_pm(i)%status = meas_pm(i)%status + status_screened
      else if ( status /= 0 ) then
        TRACEBACK; status=1; return
      end if

    end do
    
    ! ok
    status = 0

  end subroutine measupdate_pm


  ! ======================================================================
  ! ===
  ! === output
  ! ===
  ! ======================================================================
  
  
  subroutine meas_output_pm_Init( kfmo, rcfile, rckey, status )

    use GO     , only : AnyDate
    use LE_Output_Common, only : Init
  
    ! --- in/out --------------------------------
    
    type(meas_output_pm), intent(out)   ::  kfmo
    character(len=*), intent(in)            ::  rcfile
    character(len=*), intent(in)            ::  rckey
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/meas_output_pm_Init'
    
    ! --- local ---------------------------------
    
    ! --- begin ---------------------------------
    
    ! init common stuff:
    call Init( kfmo%com, rcfile, rckey, status )
    IF_NOTOK_RETURN(status=1)
    
    ! files not open yet:
    kfmo%opened = .false.
    
    ! no time range set yet:
    kfmo%tr(1) = AnyDate()
    kfmo%tr(2) = AnyDate()
    
    ! ok
    status = 0
    
  end subroutine meas_output_pm_Init
  
  
  ! ***
  

  subroutine meas_output_pm_Done( kfmo, status )
  
    use NetCDF          , only : NF90_Close
    use LE_Output_Common, only : Done
  
    ! --- in/out --------------------------------
    
    type(meas_output_pm), intent(inout)   ::  kfmo
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/meas_output_pm_Done'
    
    ! --- begin ---------------------------------
    
    ! file opened ?
    if ( kfmo%opened ) then
      ! close:
      status = NF90_Close( kfmo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! reset flag:
      kfmo%opened = .true.
    end if
    
    ! done with common stuff ...
    call Done( kfmo%com, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine meas_output_pm_Done
  
  
  ! ***
  
  
  subroutine meas_output_pm_PutOut( kfmo, key, t, status )

#ifdef pgf90
    use GO_Date
#else  
    use GO, only : TDate, NewDate, IncrDate, AnyDate, Get
    use GO, only : operator(+), operator(-), operator(<), operator(>)
    use GO, only : MidNight, rTotal
    use GO, only : wrtgol
#endif

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE
    
#ifdef pgf90
    use JAQL_Output_CF_Conventions
#else  
    use JAQL   , only : JAQL_Output_CF_names
#endif

    use Dims   , only : runF
    use Indices, only : specname, specunit
    use LE_Output_Common, only : PutOut_GlobalAttributes
    
    use LEKF_State, only : kf_with_xb, kf_with_xm
    use LEKF_State, only : xb, x, sigma
!    use LEKF_State, only : m_s, m_e
 
    ! --- in/out --------------------------------
    
    type(meas_output_pm), intent(inout)   ::  kfmo
    character(len=*), intent(in)              ::  key
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/kfmo_meas_PutOut'
    
    ! --- local ---------------------------------

    logical               ::  newfile    
    type(TDate)           ::  tday
    integer               ::  time6(6)
    real                  ::  time

    integer               ::  cmode

    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units
    character(len=512)    ::  comment
    
    integer               ::  itr !, ilev
    integer               ::  varid
    type(TDate)           ::  t0

    integer               ::  ifs
    character(len=16)     ::  varname
    character(len=512)    ::  descr

    integer               ::  i
    real                  ::  xm(maxmeas)
    real                  ::  cc(maxmeas)
    
    integer(1)            ::  bytes(maxmeas)
    
    integer               ::  ifm
    
    ! --- begin ---------------------------------
    
    ! new file if current time is not in timerange:
    newfile = (t < kfmo%tr(1)) .or. (kfmo%tr(2) < t)

    ! close old file if necessary:
    if ( newfile .and. kfmo%opened ) then
      ! close:
      status = NF90_Close( kfmo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! reset flag:
      kfmo%opened = .false.
    end if
    
    !! at midnight only ...
    !if ( .not. MidNight(t) ) then
    !  status=0; return
    !end if
    
    ! info ...
    call wrtgol('KF: put out meas_pm for : ',t); call goPr
    
    ! day is defined for (00,24]
    tday = t
    if ( MidNight(t) ) tday = tday - IncrDate(day=1)
    
    ! create new file ?
    if ( newfile ) then

      ! extract time fields for day assigned to file:
      call Get( tday, time6=time6 )

      ! time range for this file is (00,24]
      kfmo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
      kfmo%tr(2) = kfmo%tr(1) + IncrDate( day=1 )

      ! new file name:
      write (kfmo%fname,'(a,a,"_",a,"_",a,"_",i4.4,2i2.2,".nc")') &
                trim(kfmo%com%outdir), trim(kfmo%com%model), trim(kfmo%com%expid), &
                'meas-pm', time6(1:3)

      ! set creation mode flag:
      cmode = NF90_NOCLOBBER     ! do not overwrite existing files

      ! create file:
      status = NF90_Create( kfmo%fname, cmode, kfmo%ncid )
      if ( status /= 0 ) then
        write (gol,'("creating file :")'); call goErr
        write (gol,'("  ",a)') trim(kfmo%fname); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! reset flag:
      kfmo%opened = .true.

      ! write global attributes:
      call PutOut_GlobalAttributes( kfmo%com, kfmo%ncid, status )
      IF_NOTOK_RETURN(status=1)
      
      ! define dimensions:

      status = NF90_Def_Dim( kfmo%ncid, 'namelen', namelen, kfmo%dimid_namelen )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'codelen', codelen, kfmo%dimid_codelen )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'meas', nmeas, kfmo%dimid_meas )
      IF_NF90_NOTOK_RETURN(status=1)
!      status = NF90_Def_Dim( kfmo%ncid, 'longitude', nx, kfmo%dimid_lon )
!      IF_NF90_NOTOK_RETURN(status=1)
!      status = NF90_Def_Dim( kfmo%ncid, 'latitude' , ny, kfmo%dimid_lat )
!      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'time', NF90_UNLIMITED, kfmo%dimid_time )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'datelen', 6, kfmo%dimid_datelen )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variables:

!      status = NF90_Def_Var( kfmo%ncid, 'lon', NF90_REAL, kfmo%dimid_lon, varid )
!      IF_NF90_NOTOK_RETURN(status=1)
!      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'longitude' )
!      IF_NF90_NOTOK_RETURN(status=1)
!      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'longitude' )
!      IF_NF90_NOTOK_RETURN(status=1)
!      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'degrees_east' )
!      IF_NF90_NOTOK_RETURN(status=1)
!      kfmo%varid_lon = varid
!
!      status = NF90_Def_Var( kfmo%ncid, 'lat', NF90_REAL, kfmo%dimid_lat, varid )
!      IF_NF90_NOTOK_RETURN(status=1)
!      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'latitude' )
!      IF_NF90_NOTOK_RETURN(status=1)
!      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'latitude' )
!      IF_NF90_NOTOK_RETURN(status=1)
!      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'degrees_north' )
!      IF_NF90_NOTOK_RETURN(status=1)
!      kfmo%varid_lat = varid

      status = NF90_Def_Var( kfmo%ncid, 'time', NF90_REAL, (/kfmo%dimid_time/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'time' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'time' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'days since 2000-01-01 00:00:00' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'calender', 'gregorian' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_time = varid

      status = NF90_Def_Var( kfmo%ncid, 'date', NF90_INT, (/kfmo%dimid_datelen,kfmo%dimid_time/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'date and time' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'year, month, day' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_date = varid
      
      ! define station description data:
      
      status = NF90_Def_Var( kfmo%ncid, 'station_name', NF90_CHAR, (/kfmo%dimid_namelen,kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_name = varid

      status = NF90_Def_Var( kfmo%ncid, 'station_code', NF90_CHAR, (/kfmo%dimid_codelen,kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_code = varid

      status = NF90_Def_Var( kfmo%ncid, 'station_lon', NF90_REAL, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'longitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'longitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'degrees_east' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_lon = varid

      status = NF90_Def_Var( kfmo%ncid, 'station_lat', NF90_REAL, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'latitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'latitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'degrees_north' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_lat = varid

      status = NF90_Def_Var( kfmo%ncid, 'station_alt', NF90_REAL, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'altitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'altitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'm' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_alt = varid

      status = NF90_Def_Var( kfmo%ncid, 'station_analyse', NF90_BYTE, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', '1=analyse, 0=validation' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_ana = varid

      status = NF90_Def_Var( kfmo%ncid, 'station_ix', NF90_INT, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_ix = varid

      status = NF90_Def_Var( kfmo%ncid, 'station_iy', NF90_INT, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_iy = varid

      ! global tracer index
      itr = icomp

      ! CF standard name for concentration/mixing ratio/column:

      ! no comment yet
      comment = ''

      ! get names following CF conventions;
      ! store conversion factor for later usage:
      call JAQL_Output_CF_names( &
                   specname(itr), specunit(itr), &
                   cf_standard_name, cf_long_name, cf_units, &
                   kfmo%unitconv, comment, &
                   status )
      IF_NOTOK_RETURN(status=1)
      
      ! loop over filter states:
      do ifs = 1, nfs
      
        ! skip ?
        select case ( trim(fsname(ifs)) )
          case ( 'xb'     ) ; if ( .not. kf_with_xb ) cycle
          case ( 'x', 's' ) ; if ( .not. kf_with_xm ) cycle
        end select
        
        ! filter state description:
        select case ( trim(fsname(ifs)) )
          case ( 'y'      ) ; descr = 'measurement'
          case ( 'v'      ) ; descr = 'measurement error standard deviation'
          case ( 'status' ) ; descr = 'status flag'
          case ( 'xb'     ) ; descr = 'model background run'
          case ( 'x'      ) ; descr = 'ensemble mean'
          case ( 's'      ) ; descr = 'ensemble standard deviation'
          case default
            write (*,'("unsupported var name : ",a)') trim(fsname(ifs)); call goErr
            TRACEBACK; status=1; return
        end select
        
        ! loop over filter moments (1=analysis,2=forecast)
        do ifm = 1, nfm
        
          ! skip ?
          select case ( trim(fsname(ifs)) )
            case ( 'y', 'v', 'status' ) ; if ( ifm == 2 ) cycle
          end select

          ! variable name:
          write (varname,'(a,"_",a)') trim(specname(itr)), trim(fsname(ifs))
          
          ! extend for forecast moment:
          if ( ifm == 2 ) then
            varname = trim(varname)//'_f'
            descr = trim(descr)//' (forecast, prior to analysis)'
          end if

          select case ( trim(fsname(ifs)) )
            case ( 'status' )
              ! define variable:
              status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_BYTE, &
                                       (/kfmo%dimid_meas,kfmo%dimid_time/), varid )
              IF_NF90_NOTOK_RETURN(status=1)
              ! write attributes:
              status = nf90_put_att( kfmo%ncid, varid, 'description', status_description )
              IF_NF90_NOTOK_RETURN(status=1)
            case default
              ! define variable:
              status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_REAL, &
                                       (/kfmo%dimid_meas,kfmo%dimid_time/), varid )
              IF_NF90_NOTOK_RETURN(status=1)
              ! write attributes:
              status = nf90_put_att( kfmo%ncid, varid, 'standard_name', trim(cf_standard_name) )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( kfmo%ncid, varid, 'long_name', trim(cf_long_name) )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( kfmo%ncid, varid, 'units', trim(cf_units) )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( kfmo%ncid, varid, 'description', trim(descr) )
              IF_NF90_NOTOK_RETURN(status=1)
              if ( len_trim(comment) > 0 ) then
                status = nf90_put_att( kfmo%ncid, varid, 'comment', trim(comment) )
                IF_NF90_NOTOK_RETURN(status=1)
              end if
          end select

          ! store variable id:
          kfmo%varid(ifs,ifm) = varid
          
        end do  ! filter moments
        
      end do  ! filter states

      ! end defintion mode:

      status = NF90_EndDef( kfmo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! no records written yet:
      kfmo%itrec = 0
      kfmo%tprev = AnyDate()
    
    end if
    
    ! next time record ?
    if ( t > kfmo%tprev ) kfmo%itrec = kfmo%itrec + 1

    ! store this time:
    kfmo%tprev = t
    
    ! write station data only once ...    
    if ( kfmo%itrec == 1 ) then
      ! write station info; 
      ! for some reason, character arrays need to be written one-by-one:
      do i = 1, nmeas
        status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_name, trim(meas_pm(i)%name), &
                                 start=(/1,i/), count=(/len_trim(meas_pm(i)%name),1/) )
        IF_NF90_NOTOK_RETURN(status=1)
      end do
      do i = 1, nmeas
        status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_code, meas_pm(i)%code, &
                                 start=(/1,i/), count=(/codelen,1/) )
        IF_NF90_NOTOK_RETURN(status=1)
      end do
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_lon, meas_pm(1:nmeas)%lon )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_lat, meas_pm(1:nmeas)%lat )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_alt, meas_pm(1:nmeas)%alt )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_ix, meas_pm(1:nmeas)%ix )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_iy, meas_pm(1:nmeas)%iy )
      IF_NF90_NOTOK_RETURN(status=1)

      ! write analysis flags:      
      where ( meas_pm(1:nmeas)%analyse )
        bytes(1:nmeas) = 1
      else where
        bytes(1:nmeas) = 0
      end where
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_ana, bytes(1:nmeas) )
      IF_NF90_NOTOK_RETURN(status=1)

!      ! write longitudes:
!      do i = 1, nx
!        lons(i) = runF%westb + (i-0.5)*runF%dlon
!      end do
!      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_lon, lons )
!      IF_NF90_NOTOK_RETURN(status=1)
!
!      ! write latitudes:
!      do j = 1, ny
!        lats(j) = runF%southb + (j-0.5)*runF%dlat
!      end do
!      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_lat, lats )
!      IF_NF90_NOTOK_RETURN(status=1)

    end if  ! first record

    ! extract time fields for instant time:
    call Get( t, time6=time6 )

    ! time since 2000-1-1 00:00
    t0 = NewDate( time6=(/2000,01,01,00,00,00/) )
    time = rTotal( t - t0, 'day' )
    
    ! write time record:
    status = NF90_Put_Var( kfmo%ncid, kfmo%varid_time, (/time/), start=(/kfmo%itrec/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! write date record:
    status = NF90_Put_Var( kfmo%ncid, kfmo%varid_date, reshape(time6,(/6,1/)), &
                                           start=(/1,kfmo%itrec/), count=(/6,1/) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! loop over filter states:
    do ifs = 1, nfs
    
      ! loop over filter moments:
      do ifm = 1, nfm
      
        ! skip ?
        select case ( key )
          ! forecast written to ifm=2
          case ( 'forecast' )
            if ( ifm /= 2 ) cycle
          ! analysis written to ifm=1
          case ( 'analysis' )
            if ( ifm /= 1 ) cycle
          ! error ...
          case default
            write (gol,'("unsupported key : ",a)') trim(key); call goErr
            TRACEBACK; status=1; return
        end select
    
        ! init all concentrations to zero:
        cc = 0.0

        ! extract simulated concentrations:
        select case ( trim(fsname(ifs)) )
          case ( 'y' )
            if ( ifm == 2 ) cycle  ! no measurement forecast
            do i = 1, nmeas
              cc(i) = meas_pm(i)%c
            end do
          case ( 'v' )
            if ( ifm == 2 ) cycle  ! no measurement forecast
            do i = 1, nmeas
              cc(i) = meas_pm(i)%sd
            end do
          case ( 'status' )
            if ( ifm == 2 ) cycle  ! no measurement forecast
            do i = 1, nmeas
              cc(i) = meas_pm(i)%status
            end do
          case ( 'xb' )
            if ( .not. kf_with_xb ) cycle
!            xm = reshape( xb(m_s:m_e), (/maxmeas/) )
!            cc(1:nmeas) = xm(1:nmeas)
            cc(1:nmeas) = xb%m(1:nmeas)
          case ( 'x' )
            if ( .not. kf_with_xm ) cycle
!            xm = reshape( x(m_s:m_e), (/maxmeas/) )
!            cc(1:nmeas) = xm(1:nmeas)
            cc(1:nmeas) = x%m(1:nmeas)
          case ( 's' )
            if ( .not. kf_with_xm ) cycle 
!            xm = reshape( sigma(m_s:m_e), (/maxmeas/) )
!            cc(1:nmeas) = xm(1:nmeas)
            cc(1:nmeas) = sigma%m(1:nmeas)
          case default
            write (*,'("unsupported filter state name : ",a)') trim(fsname(ifs)); call goErr
            TRACEBACK; status=1; return
        end select

        ! write record:
        select case ( trim(fsname(ifs)) )
          case ( 'status' )
            status = NF90_Put_Var( kfmo%ncid, kfmo%varid(ifs,ifm), int(cc(1:nmeas),kind=1), &
                           start=(/1,kfmo%itrec/), count=(/nmeas,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            ! unit conversion:
            cc = cc * kfmo%unitconv
            ! write concentrations:
            status = NF90_Put_Var( kfmo%ncid, kfmo%varid(ifs,ifm), cc(1:nmeas), &
                           start=(/1,kfmo%itrec/), count=(/nmeas,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
        end select
        
      end do  ! filter momeents
      
    end do  ! filter states
    
    ! ok
    status = 0
    
  end subroutine meas_output_pm_PutOut


end module
