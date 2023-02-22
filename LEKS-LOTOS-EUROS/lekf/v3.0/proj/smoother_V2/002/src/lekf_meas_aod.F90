!###############################################################################
!
! NAME
!
!   kf_meas_aod  -  interface to satelite AOD measurements
!
!
! AOD DATA
!
!  File names:  aod_sx_<day>_<hour>_<instr>_hr6_summer_v2_ret.dat
!
!        day number  (00=15 July 2003)  ;  00 .. 30
!        hour number (00=midnight)      ;  05 .. 20
!
!        instruments  :
!          sr1000  :  imager
!          sr1     :  sounder (A-band instrument)
!
!        _ret :  retrieved optical depths in each layer
!        _err :  estimated random errors on the retrieval at each box/level
!
!  Data sizes:
!    _ret   :  [140,144,3]
!    _err   :  [140,144,3,3]
!
!  Kernel:
!    yr =  ya +  A(H'x-ya) + v   ,   v  ~ N(o, R  )
!   Cyr = Cya + CA(H'x-ya) + vc  ,   vc ~ N(o,CRC')   total column, C=[1,..,1]
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
!###############################################################################


module LEKF_meas_aod

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  use NetCDF, only : NF90_StrError, NF90_NOERR
  
  use LE_Output_Common, only : T_LE_Output_Common
  
  use LEKF_State, only : nz_aod
  use LEKF_Meas , only : TGainBoundInfo
  
  use bAOD, only : RAL_nx, RAL_ny, RAL_nlay
  
  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  aod_rho
  public  ::  Initmeas_AOD, Donemeas_AOD
  public  ::  Getmeas_AOD, MeasUpdate_AOD
 
  public  ::  Meas_Output_AOD
  public  ::  Init, Done, PutOut
  
  
  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'kf_meas_aod'
  
  ! maximum number of simulated AOD profiles:
  integer, parameter            ::   maxaod    = RAL_nx * RAL_ny

  ! dummy index ...
  integer, parameter      ::  i_aod = -10
  
  ! filter states to be written:
  integer, parameter            ::  nfs = 7
  character(len=6), parameter   ::  fsname(nfs) = &
             (/'y     ','R     ','status','xb    ','x     ','P     ','cfr   '/)
  
  ! filter moments: 1=analysis, 2=forecast
  integer, parameter            ::  nfm = 2
  
  ! status flags:
  integer, parameter            ::  status_validation = 1
  integer, parameter            ::  status_nan        = 2
  integer, parameter            ::  status_screened   = 4
  character(len=*), parameter   ::  status_description = &
      'status flag, 0=default, +1=validation, +2=no-data, +4=screened'  
  
  
  ! --- types --------------------------------

  ! define data type measurement
  type Measurement_AOD
    ! longitude and latitude
    real              ::  lon, lat
    ! flag if location is in current domain
    logical           ::  indomain
    !! the index numbers of simulated aod profile in the 1D state vector:
    !integer           ::  indices(:)
    integer, pointer  ::  inds(:,:)
    ! the grid indices in the 3D array
    integer           ::  ix, iy
    !! AOD layer boundaries (m):
    !real              ::  hb(0::)
    ! location in aod grid:
    integer           ::  ix_aod, iy_aod
    ! analyse or for validation only ?
    ! data accepted or nan etc ?
    logical           ::  analyse
    !logical           ::  accepted
    integer           ::  status
    ! cloud fraction:
    real              ::  cfr
    ! the AOD value:
    real, pointer     ::  yr(:)
    ! apriori, kernel, error covariance:
    !real, pointer     ::  ya(:)  ! implicit 0
    !real, pointer     ::  A(:,:) ! implicit I
    real, pointer     ::  R(:,:)
    ! horizontal radius of influence to be applied (km.)
    real              ::  rho
  end type Measurement_AOD
  
  ! output stuff
  
  type Meas_Output_AOD
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
    integer                     ::  dimid_lon
    integer                     ::  dimid_lat
    integer                     ::  dimid_lev
    integer                     ::  dimid_time
    integer                     ::  dimid_datelen
    ! dimension variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_time
    integer                     ::  varid_date
    ! tracer variables:
    real                        ::  unitconv
    integer                     ::  varid(nfs,nfm)
  end type Meas_Output_aod


  ! --- interfaces -------------------------
  
  interface Init
    module procedure meas_output_aod_Init
  end interface
  
  interface Done
    module procedure meas_output_aod_Done
  end interface

  interface PutOut
    module procedure meas_output_aod_PutOut
  end interface


  ! --- var ----------------------------------
  
  ! archive description:
  character(len=256)        ::  arch_dir
  character(len=64)         ::  arch_subdir
  character(len=64)         ::  arch_filebase
  !
  character(len=16)         ::  arch_instrument
  character(len=16)         ::  arch_scenario
  character(len=16)         ::  arch_sr
  character(len=16)         ::  arch_hr
  character(len=16)         ::  arch_period
  character(len=16)         ::  arch_retr_version
  character(len=16)         ::  arch_file_version

  ! declaration of the meas set
  type (measurement_aod)    ::  meas_aod(maxaod)

  ! actual number read:
  integer                   ::  nmeas

  ! single tracer only, thus define globaly
  integer                   ::  icomp  

  ! analyse this set ?
  logical                   ::  analyse

  ! the radius of influence for AOD values
  ! "wild guess" to be sophisticated in future
  real                      ::  aod_rho
  type(TGainBoundInfo)      ::  gbi
  
  ! profile or total column ?
  logical                   ::  aod_column
  
  ! minimum value for noise/signal ratio:
  real                      ::  aod_column_nsr_min
  
  ! adhoc error std.dev. scaling:
  real                      ::  aod_err_scale
  
  ! temporal resolution in hours:
  real                      ::  aod_dhour
  
  ! time range:
  type(TDate)               ::  ral_t0     ! time of first retrieval
  type(TDate)               ::  ral_tr(2)  ! time range with valid data


contains


  ! ======================================================================
  ! ===
  ! === input, update
  ! ===
  ! ======================================================================
  
  
  ! opens the measurement file (if present)
  ! end sets some parameters

  subroutine InitMeas_AOD( rcfile, t, status )

    use GO, only : TrcFile, Init, Done, ReadRc
    use GO, only : goGetFU
    use GO, only : TDate, NewDate, Get
    
    use dims      , only : nx, ny, runF, locF
    use constants , only : cum_days
    use units     , only : u_tmp
    !use sysdep
    use utils     , only : indomain
!    use indices   , only : i_so4a,i_no3a,i_nh4a,i_pm25,i_bc,i_dust_f,i_na_f
    use bAOD      , only : RAL_westb, RAL_southb, RAL_dlon, RAL_dlat
    use bAOD      , only : RAL_nx, RAL_ny, RAL_nlay

!    use LEKF_state  , only :  aod_s
    use LEKF_state  , only :  substate_aod
    !use LEKF_meas, only : meas, nmeas, emepmeas
    !use LEKF_meas, only : meas_unit, u_meas_min
    use LEKF_meas, only : Init

    ! --- in/out --------------------------------
    
    character(len=*), intent(in)    ::  rcfile
    type(TDate), intent(in)         ::  t
    integer, intent(out)            ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/initmeas_AOD'
    
    ! --- local ----------------------------------
    
    type(TRcFile)         ::  rcF
    integer               ::  yy, mm
    logical               ::  exist
    integer               ::  ix_aod, iy_aod
    real                  ::  lon, lat
    logical               ::  in_domain
    integer               ::  ix, iy, iz
    
    ! --- begin ----------------------------------
    
    write (gol,'("KF:     setup AOD measurements ...")'); call goPr
    
    !
    ! ~~ settings
    !
    
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.analyse', analyse, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.adir', arch_dir, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.instrument', arch_instrument, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.scenario', arch_scenario, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.retr.version', arch_retr_version, status )
    IF_NOTOK_RETURN(status=1)

    call ReadRc( rcF, 'kf.meas.aod.file.version', arch_file_version, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.rho', aod_rho, status )
    IF_NOTOK_RETURN(status=1)

    call ReadRc( rcF, 'kf.meas.aod.column', aod_column, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.column.nsr.min', aod_column_nsr_min, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.err.scale', aod_err_scale, status )
    IF_NOTOK_RETURN(status=1)
    
    call ReadRc( rcF, 'kf.meas.aod.dhour', aod_dhour, status )
    IF_NOTOK_RETURN(status=1)

    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! set name and start date for binary files with measurements:
    call Get( t, year=yy, month=mm )
    if ( mm < 7 ) then
      arch_period = 'winter'
      ! day 00 = feb 15 ??  , day 13 = aug 28 ??
      ! day 27 = aug 14 ??
      ral_t0    = NewDate( yy, 02, 15, 00, 00 )  ! feb 15 00:00
      ! all data valid:
      ral_tr(1) = NewDate( yy, 02, 15, 00, 00 )  ! feb 15 00:00
      ral_tr(2) = NewDate( yy, 03, 14, 23, 30 )  ! t0 + 27 days
    else
      arch_period = 'summer'
      ! day 00 = jul 10  ,  day 21 = jul 31
      ! day 39 = aug 18
      ral_t0    = NewDate( yy, 07, 10, 00, 00 )
      ! no cloud mask first 5 days
      ral_tr(1) = NewDate( yy, 07, 15, 00, 00 )  ! jul 15 00:00
      ral_tr(2) = NewDate( yy, 08, 18, 23, 30 )  ! t0 + 39 days
    end if
    
    ! set resolutions
    select case ( arch_instrument )
      case ( 'imager' )
        arch_sr = 'sr1000'
        arch_hr = 'hr2.5'
      case ( 'sounder' )
        arch_sr = 'sr1'
        arch_hr = 'hr6'
      case default
        write (gol,'("unsupported aod instrument : ",a)') trim(arch_instrument); call goErr
        TRACEBACK; status=1; return
    end select

    ! subdirectory: simulated_rets_v3p1_sr1_summer[<scenario>]_v2
    write (arch_subdir,'("simulated_rets",4("_",a))') &
              trim(arch_retr_version), trim(arch_sr), &
               trim(arch_period)//trim(arch_scenario), &
               trim(arch_file_version)
    
    ! file base: aod_sx_00_00_sr1_hr6_summer[<scenario>]_v2
    write (arch_filebase,'("aod_sx",2("_",i2.2),4("_",a))') &
              00, 00, trim(arch_sr), trim(arch_hr), &
               trim(arch_period)//trim(arch_scenario), &
               trim(arch_file_version)

    ! set number of measured layers:
    if ( aod_column ) then
      nz_aod = 1
    else
      nz_aod = RAL_nlay
    end if
    
    ! gain bounds:
    if ( analyse ) then
      call Init( gbi, aod_rho, status )
      IF_NOTOK_RETURN(status=1)
    end if
    

    !
    ! ~~ tracers
    !
    
    ! set tracer to be the dummy tracer 'aod'
    icomp = i_aod
    
    ! components of pm25 sum:    
    !ii_pm25s = (/i_so4a,i_no3a,i_nh4a,i_pm25,i_bc,i_dust_f,i_na_f/)
    

    !
    ! ~~ measurement locations
    !
    
    ! nothing stored yet ...
    nmeas = 0
    
    ! loop over AOD grid cells:
    do iy_aod = 1, RAL_ny
      do ix_aod = 1, RAL_nx
      
        ! center of grid cell:
        lon =  RAL_westb + (ix_aod-0.5)*RAL_dlon
        lat = RAL_southb + (iy_aod-0.5)*RAL_dlat
      
        ! check if meas is in domain;
        ! return cell indices (ix,iy) and flag:
        call InDomain( lon, lat, runF%dlon, runF%dlat, runF%southb, runF%westb, &
                       ix, iy, in_domain )

                !if ( ix_aod <= 10 .and. (iy_aod <= 10) ) then
                !print *, 'xxx ', ix_aod, iy_aod, lon, lat, in_domain, ix, iy
                !end if

        ! skip ?        
        if ( .not. in_domain ) cycle

        ! increase measuremnt counter:
        nmeas = nmeas + 1
        
        ! check ...
        if ( nmeas > maxaod ) then
          write (gol,'("total number of grid cells exceeds maxaod = ",i6)') maxaod; call goErr
          TRACEBACK; status=1; return
        end if

        ! storage:
        allocate( meas_aod(nmeas)%inds(0:7,RAL_nlay)   )
        allocate( meas_aod(nmeas)%yr(nz_aod)         )
        !allocate( meas_aod(nmeas)%ya(nz_aod)         )
        !allocate( meas_aod(nmeas)%A(nz_aod,RAL_nlay) )
        allocate( meas_aod(nmeas)%R(nz_aod,nz_aod) )

        ! store location:
        meas_aod(nmeas)%ix_aod   = ix_aod
        meas_aod(nmeas)%iy_aod   = iy_aod
        meas_aod(nmeas)%lon      = lon
        meas_aod(nmeas)%lat      = lon
        meas_aod(nmeas)%indomain = .true.
        meas_aod(nmeas)%ix       = ix
        meas_aod(nmeas)%iy       = iy
        
        ! indices in kf state:
        meas_aod(nmeas)%inds = 0
        do iz = 1, RAL_nlay
          ! element 0 defines substate, rest the index in the substate array
          meas_aod(nmeas)%inds(0:3,iz) = (/substate_aod,ix,iy,iz/)
        end do
      
        ! radius of influence to be applied (km.)
        meas_aod(nmeas)%rho       = aod_rho   ! km

        ! analyse ?
        meas_aod(nmeas)%analyse   = analyse

        !! not accepted by default ...
        !meas_aod(nmeas)%accepted  = .false.
        ! accepted by default ...
        meas_aod(nmeas)%status = 0

        ! no data yet ..
        meas_aod(nmeas)%yr = 0.0
        !meas_aod(nmeas)%ya = 0.0
        !meas_aod(nmeas)%A  = 0.0
        meas_aod(nmeas)%R  = 0.0
        
      end do  ! ix
    end do  ! iy
    
    
    !
    ! done
    !
    
    ! ok
    status = 0

  end subroutine InitMeas_AOD


  ! ***
  

  ! read station locations etc

  subroutine DoneMeas_AOD( status )

    use LEKF_meas, only : Done

    ! --- in/out --------------------------------
    
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/DoneMeas_AOD'
    
    ! --- local ------------------------------
    
    integer     :: imeas
    
    ! --- begin ------------------------------

    do imeas = 1, nmeas
      deallocate( meas_aod(imeas)%inds )
      deallocate( meas_aod(imeas)%yr )
      !deallocate( meas_aod(imeas)%ya )
      !deallocate( meas_aod(imeas)%A  )
      deallocate( meas_aod(imeas)%R  )
    end do

    ! done with gain bounds:
    if ( analyse ) then
      call Done( gbi, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0

  end subroutine DoneMeas_AOD
  
  
  ! ***
  
  
  subroutine GetMeas_AOD( t, status )

    use GO, only : TDate, Get, Precisely, wrtgol
    use GO, only : operator(<), operator(>), operator(-), rTotal
    use GO, only : goGetFU
    
    ! --- in/out ----------------------------
    
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Getmeas_aod'
    
    ! --- local -----------------------------
    
    integer             ::  day_nr, halfhour_nr
    integer             ::  hour, minu
    character(len=64)   ::  bname
    character(len=256)  ::  fname_ret, fname_err, fname_cfr
    logical             ::  exist_ret, exist_err
    real(4)             ::  pat(RAL_nx,RAL_ny)
    integer             ::  recl
    integer             ::  irec
    integer             ::  fu
    integer             ::  ix, iy, iz, k, l
    integer             ::  imeas
    
    ! --- begin -----------------------------
    
    ! no data by default ...
    do imeas = 1, nmeas
      meas_aod(imeas)%yr = -999.0
      meas_aod(imeas)%R  = -999.0
      ! no data by default:
      meas_aod(imeas)%status = status_nan
    end do

    ! validation only ?
    do imeas = 1, nmeas
      if ( .not. meas_aod(imeas)%analyse ) then
        meas_aod(imeas)%status = meas_aod(imeas)%status + status_validation
      end if
    end do

    ! check ...
    if ( t < ral_t0 ) then
      call wrtgol( 'KF:     AOD not available before ', ral_t0 ); call goPr
      status=0; return
    end if 
    if ( t < ral_tr(1) ) then
      call wrtgol( 'KF:     AOD not valid before ', ral_tr(1) ); call goPr
      status=0; return
    end if 
    if ( t > ral_tr(2) ) then
      call wrtgol( 'KF:     AOD not valid after ', ral_tr(2) ); call goPr
      status=0; return
    end if 

    ! check time:
    if ( .not. Precisely(t,aod_dhour,'hour')  ) then
      write (gol,'("KF:     AOD available every ",f4.1," hour, ")') aod_dhour; call goPr
      call wrtgol( 'KF:     not for ', t ); call goPr
      status=0; return
    end if
    
    ! days since start of data set (00,01,...)
    day_nr = floor(rTotal(t-ral_t0,'day'))
    
    ! half hour in day: 
    !   00:00 -> 00
    !   00:30 -> 01
    !   01:00 -> 02
    !      ..
    call Get(t,hour=hour,min=minu)
    halfhour_nr = floor( ( hour*60.0 + minu*1.0 )/30.0 )
    
    ! subdirectory: simulated_rets_v3p1_sr1_summer_v2
    write (arch_subdir,'("simulated_rets",4("_",a))') &
              trim(arch_retr_version), trim(arch_sr), &
              trim(arch_period), trim(arch_file_version)
    
    ! basename of data files:
    !   aod_sx_00_00_sr1_hr6_summer_v2
    write (bname,'("aod_sx",2("_",i2.2),4("_",a))') &
                     day_nr, halfhour_nr, &
                     trim(arch_sr), trim(arch_hr), &
                     trim(arch_period), trim(arch_file_version)
    
    ! data files:
    write (fname_ret,'(a,"/",a,"/",a,"_",a,".dat")') trim(arch_dir), trim(arch_subdir), trim(bname), 'ret'
    write (fname_err,'(a,"/",a,"/",a,"_",a,".dat")') trim(arch_dir), trim(arch_subdir), trim(bname), 'err'
    write (fname_cfr,'(a,"/",a,"/",a,"_",a,".dat")') trim(arch_dir), trim(arch_subdir), trim(bname), 'cfr'
    
    ! available ?
    inquire( file=trim(fname_ret), exist=exist_ret )
    inquire( file=trim(fname_err), exist=exist_err )
    ! both missing ? then no data availble for this hour:
    if ( (.not. exist_ret) .and. (.not. exist_err) ) then
      write (gol,'("KF:     no retrievals available; files not found:")'); call goPr
      write (gol,'("KF:       archive  : ",a)') trim(arch_dir); call goPr
      write (gol,'("KF:       subdir   : ",a)') trim(arch_subdir); call goPr
      write (gol,'("KF:       files    : ",a,"_*.dat")') trim(bname); call goPr
      status=0; return
    end if
    ! one missing ? strange ...
    if ( exist_ret .neqv. exist_err ) then
      write (gol,'("missing one of the two aod data files: ")'); call goErr
      write (gol,'("  ret  : ",a)') trim(fname_ret); call goErr
      write (gol,'("  err  : ",a)') trim(fname_err); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! info ...    
    call wrtgol( 'KF:     load AOD for ', t ); call goPr
    write (gol,'("KF:       files    : ",a,"_*.dat")') trim(bname); call goPr
    
    ! record length:
    inquire( iolength=recl ) pat
    
    !
    ! ret
    !
    
    ! free file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( unit=fu, file=trim(fname_ret), status='old', &
            form='unformatted', access='direct', recl=recl, &
            iostat=status )
    if (status/=0) then
      write (gol,'("from opening file:")'); call goErr
      write (gol,'("  ret name  : ",a)') trim(fname_ret); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! read records:
    irec = 0
    do iz = 1, RAL_nlay
      ! next record number:
      irec = irec + 1
      ! read 2D field
      read (fu,rec=irec,iostat=status) pat
      if (status/=0) then
        write (gol,'("from reading measurement record:")'); call goErr
        write (gol,'("  file name : ",a)') trim(fname_ret); call goErr
        write (gol,'("  record    : ",i6)') irec; call goErr
        write (gol,'("  iostat    : ",i6)') status; call goErr
        TRACEBACK; status=1; return
      end if
      ! store:
      do imeas = 1, nmeas
        ! aod retrieval:
        ix = meas_aod(imeas)%ix
        iy = meas_aod(imeas)%iy
        if ( aod_column ) then
          ! total column ...
          if ( iz == 1 ) meas_aod(imeas)%yr(iz) = 0.0
          meas_aod(imeas)%yr(1) = meas_aod(imeas)%yr(1) + pat(ix,iy)
        else
          ! store point in profile:
          meas_aod(imeas)%yr(iz) = pat(ix,iy)
        end if
      end do
    end do
    
    ! valid ?
    do imeas = 1, nmeas
      !meas_aod(imeas)%accepted = all( meas_aod(imeas)%yr > 0.0 )
      if ( all( meas_aod(imeas)%yr > 0.0 ) ) then
        meas_aod(imeas)%status = meas_aod(imeas)%status - status_nan
      end if
    end do

    ! close data file:
    close( unit=fu, iostat=status )
    if (status/=0) then
      write (gol,'("from closing file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname_ret); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
    !
    ! err
    !

    ! free file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( unit=fu, file=trim(fname_err), status='old', &
            form='unformatted', access='direct', recl=recl, &
            iostat=status )
    if (status/=0) then
      write (gol,'("from opening file:")'); call goErr
      write (gol,'("  ret name  : ",a)') trim(fname_err); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! read records:
    irec = 0
    do k = 1, RAL_nlay
      do l = 1, RAL_nlay
        ! next record number:
        irec = irec + 1
        ! read 2D field
        read (fu,rec=irec,iostat=status) pat
        if (status/=0) then
          write (gol,'("from reading measurement record:")'); call goErr
          write (gol,'("  file name : ",a)') trim(fname_ret); call goErr
          write (gol,'("  record    : ",i6)') irec; call goErr
          write (gol,'("  iostat    : ",i6)') status; call goErr
          TRACEBACK; status=1; return
        end if
        ! store:
        do imeas = 1, nmeas
          ix = meas_aod(imeas)%ix
          iy = meas_aod(imeas)%iy
          if ( aod_column ) then
            ! error variance of total column ...
            if ( (k==1) .and. (l==1) ) meas_aod(imeas)%R(1,1) = 0.0
            meas_aod(imeas)%R(1,1) = meas_aod(imeas)%R(1,1) + pat(ix,iy)
          else
            ! store point in profile:
            meas_aod(imeas)%R(k,l) = pat(ix,iy)
          end if
          !
        end do
      end do  ! l
    end do  ! k
    
    ! apply extra scaling:
    meas_aod(imeas)%R = meas_aod(imeas)%R * (aod_err_scale**2)

    ! minimum of standard deviation in column ?
    if ( aod_column .and. (aod_column_nsr_min > 0.0) ) then
      ! loop over all measurements:
      do imeas = 1, nmeas
        ! valid data ?
        if ( meas_aod(imeas)%yr(1) > 0.0 ) then
          ! variance should exceed square of fraction of retrieval:
          meas_aod(imeas)%R(1,1) = max( (meas_aod(imeas)%yr(1) * aod_column_nsr_min)**2, meas_aod(imeas)%R(1,1) )
        end if
      end do  ! measurements
    end if  ! columns only
    
    ! close data file:
    close( unit=fu, iostat=status )
    if (status/=0) then
      write (gol,'("from closing file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname_err); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
                !! adhoc: R in aod files too small ...
                !write (gol,'("WARNING - set aod std.dev. to 20% ...")'); call goPr
                !if ( .not. aod_column ) stop 'TO BE DONE: aod std.dev. for profile'
                !do imeas = 1, nmeas
                !  meas_aod(imeas)%R(1,1) = ( 0.20 * meas_aod(imeas)%yr(1) )**2
                !end do
                
    
    !
    ! cloud fraction
    !
    
    ! free file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( unit=fu, file=trim(fname_cfr), status='old', &
            form='unformatted', access='direct', recl=recl, &
            iostat=status )
    if (status/=0) then
      write (gol,'("from opening file:")'); call goErr
      write (gol,'("  cfr name  : ",a)') trim(fname_cfr); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! single record only:
    irec = 1
    ! read 2D field
    read (fu,rec=irec,iostat=status) pat
    if (status/=0) then
      write (gol,'("from reading measurement record:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname_cfr); call goErr
      write (gol,'("  record    : ",i6)') irec; call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    ! store:
    do imeas = 1, nmeas
      ! aod retrieval:
      ix = meas_aod(imeas)%ix
      iy = meas_aod(imeas)%iy
      meas_aod(imeas)%cfr = pat(ix,iy)
    end do
    
    ! close data file:
    close( unit=fu, iostat=status )
    if (status/=0) then
      write (gol,'("from closing file:")'); call goErr
      write (gol,'("  file name : ",a)') trim(fname_cfr); call goErr
      write (gol,'("  iostat    : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
    !
    ! apriori and kernel
    !
    
    !! loop over measurements:
    !do imeas = 1, nmeas
    !
    !  ! dummy apriori:
    !  meas_aod(imeas)%ya = 0.0
    !
    !  ! kernel
    !  if ( aod_column ) then
    !    ! A = [1,..,1] * I = [1,..,1]
    !    meas_aod(imeas)%A = 1
    !  else
    !    if ( nz_aod /= RAL_nlay ) stop 'kernel: nz_aod /= RAL_nlay'
    !    ! A = I
    !    meas_aod(imeas)%A = 0.0
    !    do k = 1, nz_aod
    !      meas_aod(imeas)%A(k,k) = 1.0
    !    end do
    !  end if
    !
    !end do  ! measurements    

    !
    ! done
    !

    ! ok
    status = 0

  end subroutine GetMeas_AOD


  ! ***
  
  
  ! implementation of measurement update
  ! routine expects a measurement list +
  ! a list specifying the corresponding element nr.
  ! of the statevector.

  subroutine MeasUpdate_AOD( t, nmodes, status )

!    use GO, only : TDate, Midnight
    use LEKF_state, only : maxmodes
!    use LEKF_state, only : x, modes
    use LEKF_state, only : x, Ens
!    use LEKF_state, only : nstate
!    use LEKF_meas , only : meas_update_profile
!    use LEKF_meas , only : meas_update_scalar
!    use LEKF_meas , only : meas_update
    use LEKF_meas , only : meas_update
    use LEKF_meas , only : meas_update_scalar

    ! --- in/out -----------------------------

    type(TDate), intent(in)   ::  t    
    integer, intent(in)       ::  nmodes     ! actual number of modes
    integer, intent(out)      ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MeasUpdate_AOD'

    ! --- local ------------------------------
    
    integer       ::  i, j
!    integer       ::  ind
!    !real          ::  z1(nmeas)
!    real          ::  c, sd
!    integer       ::  ix, iy
!    integer       ::  ihist
!    real          ::  rho
!
!    real          ::  prf_yr(nz_aod)
!    real          ::  prf_ya(nz_aod)
!    real          ::  prf_A(nz_aod,nz_aod)
!    real          ::  prf_R(nz_aod,nz_aod)
!    integer       ::  inds(nz_aod)

    integer       ::  nana, iana

!    integer       ::  ip
          
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
    write (gol,'("KF:   measurement update AOD")'); call goPr

!    ! ensemble members:
!    write (gol,'("KF:     form ensemble ...")'); call goPr
!    do j = 1, nmodes
!      modes(:,j) = x + modes(:,j) * sqrt(nmodes-1.0)
!    end do

    write (gol,'("KF:     nmeas : ",i6)') nmeas; call goPr
    
    !nana = count( meas_aod(1:nmeas)%analyse .and. meas_aod(1:nmeas)%accepted )
    nana = count( meas_aod(1:nmeas)%analyse .and. (meas_aod(1:nmeas)%status == 0) )
    iana = 0
    ! loop over measurements:
    do i = 1, nmeas
     
      !write (gol,'("KF:   analyse, accepted : ",i6,2l3)') i, meas_aod(i)%analyse, meas_aod(i)%accepted; call goPr

      ! skip if not accepted (validation, no data?):
      !if ( .not. meas_aod(i)%accepted ) cycle
      if ( meas_aod(i)%status /= 0 ) cycle

      ! hoezee!
      iana = iana + 1
      if ( modulo(iana,1000) == 0 ) then
        write (gol,'("KF:     analyse ",2i6)') iana, nana; call goPr
      end if

      ! ihist       time level under consideration 
      !             (determines together with rho the influence area)
      ! rho         influence radius
      !ihist = 0  ! <-- dummy; no idea what to do with this ...

      ! analyse profile;
      ! only modes are adjust, thus update mean later on:
!      call meas_update_profile_fast( &
!                   nz_aod, meas_aod(i)%c, meas_aod(i)%c_apriori, meas_aod(i)%A, meas_aod(i)%R, &
!                   meas_aod(i)%rho, meas_aod(i)%ix, meas_aod(i)%iy, meas_aod(i)%indices, nmodes, ihist, &
!                   status )
      if ( nz_aod == 1 ) then
        call meas_update_scalar( meas_aod(i)%yr(1), meas_aod(i)%R(1,1), meas_aod(i)%inds(:,1), &
                            meas_aod(i)%ix, meas_aod(i)%iy, gbi, &
                            nmodes, status )
        if ( status == -1 ) then
          meas_aod(i)%status = meas_aod(i)%status + status_screened
        else if ( status /= 0 ) then
          TRACEBACK; status=1; return
        end if
      else
        call meas_update( nz_aod, meas_aod(i)%yr, meas_aod(i)%R, meas_aod(i)%inds, &
                            meas_aod(i)%ix, meas_aod(i)%iy, gbi, &
                            nmodes, status )
        if ( status == -1 ) then
          meas_aod(i)%status = meas_aod(i)%status + status_screened
        else if ( status /= 0 ) then
          TRACEBACK; status=1; return
        end if
      end if
    end do
    
!     ! modes are deviation of ensemble from mean: 
!     write (gol,'("KF:     form modes ...")'); call goPr
!     do j = 1, nmodes
!       modes(:,j) = ( modes(:,j) - x(:) )/sqrt(nmodes-1.0)
!     enddo

    ! ok
    status = 0

  end subroutine MeasUpdate_AOD


  ! ======================================================================
  ! ===
  ! === output
  ! ===
  ! ======================================================================
  
  subroutine meas_output_aod_Init( kfmo, rcfile, rckey, status )

    use GO     , only : AnyDate
    use LE_Output_Common, only : Init
  
    ! --- in/out --------------------------------
    
    type(meas_output_aod), intent(out)      ::  kfmo
    character(len=*), intent(in)            ::  rcfile
    character(len=*), intent(in)            ::  rckey
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/meas_output_aod_Init'
    
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
    
  end subroutine meas_output_aod_Init
  
  
  ! ***
  

  subroutine meas_output_aod_Done( kfmo, status )
  
    use NetCDF          , only : NF90_Close
    use LE_Output_Common, only : Done
  
    ! --- in/out --------------------------------
    
    type(meas_output_aod), intent(inout)      ::  kfmo
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/meas_output_aod_Done'
    
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
    
  end subroutine meas_output_aod_Done
  
  
  ! ***
  
  
  subroutine meas_output_aod_PutOut( kfmo, key, t, status )

    use GO, only : TDate, IncrDate, Get, NewDate, AnyDate
    use GO, only : operator(+), operator(-), operator(<), operator(>), rTotal
    use GO, only : Precisely, MidNight, wrtgol

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE
    
    use JAQL   , only : JAQL_Output_CF_names

    use Dims   , only : runF
    use Indices, only : specname, specunit
    use LE_Output_Common, only : PutOut_GlobalAttributes
    
    use LEKF_State, only : kf_with_xb, kf_with_xm
    use LEKF_State, only : xb, x, sigma, P_aod
    use LEKF_State, only : nx, ny
!    use LEKF_State, only : aod_s, aod_e
 
    ! --- in/out --------------------------------
    
    type(meas_output_aod), intent(inout)      ::  kfmo
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

    integer               ::  i, j, k, l
    integer               ::  ix, iy
    real                  ::  lons(nx), lats(ny)
    real                  ::  aod(nx,ny,nz_aod)
    real                  ::  cc(nx,ny,nz_aod,nz_aod)
    logical               ::  cov
    
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
    
    ! half hourly only ...
    if ( .not. Precisely(t,0.5,'hour') ) then
      status=0; return
    end if
    
    ! info ...
    call wrtgol('KF: put out AOD for : ',t); call goPr
    
    ! day is defined for (00,24]
    tday = t
    if ( MidNight(t) ) tday = tday - IncrDate(day=1)
    
    ! extract time fields:
    call Get( tday, time6=time6 )

    ! create new file ?
    if ( newfile ) then

      ! time range for this file is (00,24]
      kfmo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
      kfmo%tr(2) = kfmo%tr(1) + IncrDate( day=1 )

      ! new file name:
      write (kfmo%fname,'(a,a,"_",a,"_",a,"_",i4.4,2i2.2,".nc")') &
                trim(kfmo%com%outdir), trim(kfmo%com%model), trim(kfmo%com%expid), &
                'meas-aod', time6(1:3)

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
      
      ! add archive description
      status = NF90_Put_Att( kfmo%ncid, NF90_GLOBAL, 'ral_instrument', trim(arch_instrument) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, NF90_GLOBAL, 'ral_sr', trim(arch_sr) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, NF90_GLOBAL, 'ral_hr', trim(arch_hr) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, NF90_GLOBAL, 'ral_retr_version', trim(arch_retr_version) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, NF90_GLOBAL, 'ral_file_version', trim(arch_file_version) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, NF90_GLOBAL, 'ral_files', trim(arch_subdir)//'/'//trim(arch_filebase) )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! define dimensions:

      status = NF90_Def_Dim( kfmo%ncid, 'longitude', nx, kfmo%dimid_lon )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'latitude' , ny, kfmo%dimid_lat )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'level' , nz_aod, kfmo%dimid_lev )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'time', NF90_UNLIMITED, kfmo%dimid_time )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'datelen', 6, kfmo%dimid_datelen )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variables:

      status = NF90_Def_Var( kfmo%ncid, 'lon', NF90_REAL, kfmo%dimid_lon, varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'longitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'longitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'degrees_east' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_lon = varid

      status = NF90_Def_Var( kfmo%ncid, 'lat', NF90_REAL, kfmo%dimid_lat, varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'latitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'latitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'degrees_north' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_lat = varid

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
      
      ! global tracer index
      itr = icomp

      ! CF standard name for concentration/mixing ratio/column:

      ! no comment yet
      comment = ''

      ! get names following CF conventions;
      ! store conversion factor for later usage:
      call JAQL_Output_CF_names( &
                   'AOD', 'm', &
                   cf_standard_name, cf_long_name, cf_units, &
                   kfmo%unitconv, comment, &
                   status )
      IF_NOTOK_RETURN(status=1)
      
      ! loop over filter states:
      do ifs = 1, nfs
      
        ! skip ?
        select case ( trim(fsname(ifs)) )
          case ( 'xb'     ) ; if ( .not. kf_with_xb ) cycle
          case ( 'x', 'P' ) ; if ( .not. kf_with_xm ) cycle
        end select
        
        ! filter state description:
        select case ( trim(fsname(ifs)) )
          case ( 'y'      ) ; descr = 'measurement'
          case ( 'R'      ) ; descr = 'measurement error covariance'
          case ( 'status' ) ; descr = 'error flag'
          case ( 'xb'     ) ; descr = 'model background run'
          case ( 'x'      ) ; descr = 'ensemble mean'
          case ( 'P'      ) ; descr = 'ensemble covariance'
          case ( 'cfr'    ) ; descr = 'cloud fraction'
          case default
            write (gol,'("unsupported var name : ",a)') trim(fsname(ifs)); call goErr
            TRACEBACK; status=1; return
        end select
        
        ! loop over filter moments (1=analysis,2=forecast)
        do ifm = 1, nfm
        
          ! skip ?
          !select case ( trim(fsname(ifs)) )
          !  case ( 'y', 'R' ) ; if ( ifm == 2 ) cycle
          !end select
          ! always skip the forecast, otherwise files become too large ...
          if ( ifm == 2 ) cycle

          ! covariance ?
          select case ( trim(fsname(ifs)) )
            case ( 'R', 'P' ) ; cov = .true.
            case default      ; cov = .false.
          end select

          ! variable name:
          select case ( itr )
            case ( i_aod ) ;  write (varname,'(a,"_",a)') 'AOD', trim(fsname(ifs))
            case default   ;  write (varname,'(a,"_",a)') trim(specname(itr)), trim(fsname(ifs))
          end select
          
          ! extend for forecast moment:
          if ( ifm == 2 ) varname = trim(varname)//'_f'

          select case ( trim(fsname(ifs)) )
            case ( 'status' )
              ! define variable:
              status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_BYTE, &
                                       (/kfmo%dimid_lon,kfmo%dimid_lat,kfmo%dimid_time/), varid )
              IF_NF90_NOTOK_RETURN(status=1)
              ! write attributes:
              status = nf90_put_att( kfmo%ncid, varid, 'description', status_description )
              IF_NF90_NOTOK_RETURN(status=1)
            case ( 'cfr' )
              ! define variable:
              status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_REAL, &
                                       (/kfmo%dimid_lon,kfmo%dimid_lat,kfmo%dimid_time/), varid )
              IF_NF90_NOTOK_RETURN(status=1)
              ! write attributes:
              status = nf90_put_att( kfmo%ncid, varid, 'description', 'cloud fraction' )
              IF_NF90_NOTOK_RETURN(status=1)
            case default
              ! define variable:
              if ( cov ) then
                status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_REAL, &
                                         (/kfmo%dimid_lon,kfmo%dimid_lat,kfmo%dimid_lev,kfmo%dimid_lev,kfmo%dimid_time/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
              else
                status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_REAL, &
                                         (/kfmo%dimid_lon,kfmo%dimid_lat,kfmo%dimid_lev,kfmo%dimid_time/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
              end if
              ! write attributes:
              status = nf90_put_att( kfmo%ncid, varid, 'standard_name', trim(cf_standard_name) )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( kfmo%ncid, varid, 'long_name', trim(cf_long_name) )
              IF_NF90_NOTOK_RETURN(status=1)
              if ( cov ) then
                status = nf90_put_att( kfmo%ncid, varid, 'units', '('//trim(cf_units)//')**2' )
                IF_NF90_NOTOK_RETURN(status=1)
              else
                status = nf90_put_att( kfmo%ncid, varid, 'units', trim(cf_units) )
                IF_NF90_NOTOK_RETURN(status=1)
              end if
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
    
      ! write longitudes:
      do i = 1, nx
        lons(i) = runF%westb + (i-0.5)*runF%dlon
      end do
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_lon, lons )
      IF_NF90_NOTOK_RETURN(status=1)

      ! write latitudes:
      do j = 1, ny
        lats(j) = runF%southb + (j-0.5)*runF%dlat
      end do
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_lat, lats )
      IF_NF90_NOTOK_RETURN(status=1)

    end if  ! first record

    ! time since 2000-1-1 00:00
    t0 = NewDate( time6=(/2000,01,01,00,00,00/) )
    time = rTotal( t - t0, 'day' )
    
    ! write time record:
    status = NF90_Put_Var( kfmo%ncid, kfmo%varid_time, (/time/), start=(/kfmo%itrec/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! current time up to seconds:
    call Get( t, time6=time6 )
    
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
            !if ( ifm /= 2 ) cycle
            cycle   ! do no write forecast fields, too large ...
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
        
        ! no covariance by default ...
        cov = .false.

        ! extract simulated concentrations:
        select case ( trim(fsname(ifs)) )
          case ( 'y' )
            if ( ifm == 2 ) cycle  ! no measurement forecast
            do i = 1, nmeas
              ix = meas_aod(i)%ix
              iy = meas_aod(i)%iy
              cc(ix,iy,:,1) = meas_aod(i)%yr
            end do
          case ( 'R' )
            if ( ifm == 2 ) cycle  ! no measurement forecast
            do i = 1, nmeas
              ix = meas_aod(i)%ix
              iy = meas_aod(i)%iy
              cc(ix,iy,:,:) = meas_aod(i)%R
            end do
            cov = .true.
          case ( 'cfr' )
            if ( ifm == 2 ) cycle  ! no measurement forecast
            do i = 1, nmeas
              ix = meas_aod(i)%ix
              iy = meas_aod(i)%iy
              cc(ix,iy,1,1) = meas_aod(i)%cfr
            end do
          case ( 'status' )
            if ( ifm == 2 ) cycle  ! no measurement forecast
            do i = 1, nmeas
              ix = meas_aod(i)%ix
              iy = meas_aod(i)%iy
              cc(ix,iy,1,1) = meas_aod(i)%status
            end do
          case ( 'xb' )
            if ( .not. kf_with_xb ) cycle
!            aod = reshape( xb(aod_s:aod_e), (/nx,ny,nz_aod/) )
            do i = 1, nmeas
              ix = meas_aod(i)%ix
              iy = meas_aod(i)%iy
              do k = 1, nz_aod
                cc(ix,iy,k,1) = xb%aod(ix,iy,k)
              end do
            end do
          case ( 'x' )
            if ( .not. kf_with_xm ) cycle
!            aod = reshape( x(aod_s:aod_e), (/nx,ny,nz_aod/) )
            do i = 1, nmeas
              ix = meas_aod(i)%ix
              iy = meas_aod(i)%iy
              do k = 1, nz_aod
                cc(ix,iy,k,1) = x%aod(ix,iy,k)
              end do
            end do
          case ( 'P' )
            if ( .not. kf_with_xm ) cycle
!            aod = reshape( sigma(aod_s:aod_e), (/nx,ny,nz_aod/) )
            do i = 1, nmeas
              ix = meas_aod(i)%ix
              iy = meas_aod(i)%iy
              do k = 1, nz_aod
                do l = 1, nz_aod
                  cc(ix,iy,k,l) = P_aod(ix,iy,k,l)
                end do
              end do
            end do
            cov = .true.
          case default
            write (*,'("unsupported filter state name : ",a)') trim(fsname(ifs)); call goErr
            TRACEBACK; status=1; return
        end select

        ! write record:
        select case ( trim(fsname(ifs)) )
          case ( 'status' )
            ! write 2D error flags:
            status = NF90_Put_Var( kfmo%ncid, kfmo%varid(ifs,ifm), int(cc(:,:,1,1),kind=1), &
                           start=(/1,1,kfmo%itrec/), count=(/nx,ny,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 'cfr' )
            ! write 2D cloud fractions:
            status = NF90_Put_Var( kfmo%ncid, kfmo%varid(ifs,ifm), cc(:,:,1,1), &
                           start=(/1,1,kfmo%itrec/), count=(/nx,ny,1/) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            ! write concentrations:
            if ( cov ) then
              ! unit conversion:
              cc = cc * kfmo%unitconv**2
              ! write:
              status = NF90_Put_Var( kfmo%ncid, kfmo%varid(ifs,ifm), cc, &
                             start=(/1,1,1,1,kfmo%itrec/), count=(/nx,ny,nz_aod,nz_aod/) )
              IF_NF90_NOTOK_RETURN(status=1)
            else
              ! unit conversion:
              cc = cc * kfmo%unitconv
              ! write:
              status = NF90_Put_Var( kfmo%ncid, kfmo%varid(ifs,ifm), cc(:,:,:,1), &
                             start=(/1,1,1,kfmo%itrec/), count=(/nx,ny,nz_aod,1/) )
              IF_NF90_NOTOK_RETURN(status=1)
            end if
        end select
        
      end do  ! filter momeents
      
    end do  ! filter states
    
    ! ok
    status = 0
    
  end subroutine meas_output_aod_PutOut


end module



