!###############################################################################
!
! NAME
!
!   kf_meas_ground  -  interface to measurements from LML network
!
! SETS OF GROUND OBSERVATIONS
!
!   Ground observations are read from a number of data sets.
!   The rcfile contains a list of data set names:
!
!      kf.meas.ground.sets     :   lml-ana lml-val
!
!   Each set is described by:
!    o data set type : lml, ...
!    o list of stations; for each station (depending on data set type)
!      the following info could be provided:
!         * station name
!         * eventually station code
!         * longitude, latitude, altitude
!         * measured components: O3, NO, NO2, SO2, SO4a, ...
!           and unit(s)
!    o flags for analsysis or validation all components and exceptions
!    o correlation cut-off scale for analysis gain (km)
!    o representation error stuff
!
!
! LML DATA
!
!  Example of rcfile settings:
!
!      ! type of data file:
!      kf.meas.ground.lml-ana.type        :  lml
!      ! name of station list file:
!      kf.meas.ground.lml-ana.statlist    :  /data/lml/lml-ana-statlist
!      ! location of data files:
!      kf.meas.ground.lml-ana.datadir     :  /data/lml/ascii/
!      ! observed component:
!      kf.meas.ground.lml-ana.specname    :  O3
!      ! analyse measurements (T|F) ?
!      kf.meas.ground.lml-ana.analyse     :  T
!      ! correlation cut-off scale (km) :
!      kf.meas.ground.lml-ana.rho         :  100.0
!      ! representation error fraction of concentration, minimu, maximum
!      kf.meas.ground.lml-ana.fmm         :  0.10, 0.0, 999.9
!
!  The station list file contains the names of sites.
!  Location is read from the header in the data file.
!  For each measured spec a data file should be available.
!  Example:
!
!    Posterholt-Vlodropperweg                
!    Vredepeel-Vredeweg                      
!    Wijnandsrade-Opfergeltstraat            
!      :
!
!  Example of data file:
!
!    ---[o3ned12m2006o.spr]---------------------------------------------
!
!    Component: O3               
!
!     unit        ;ug/m3     ;ug/m3     ;...
!     latitude (y);        51.12;        51.54;...
!     longitude(x);         6.04;         5.85;...
!     station numb;          107;          131;...
!     Height/type ;observation;observation;...
!     Origin      ;RILPLUS   ;RILPLUS   ;...
!    Hour\Station ;Posterholt-Vlodropperweg                ;Vredepeel-Vredeweg                      ;...
!    -------------;     ----------;     ----------;...
!    01/01/2006 01; .4815E+02; .5583E+02; ...
!    01/01/2006 02; .5089E+02; .5480E+02; ...
!    01/01/2006 03; .5189E+02; .5378E+02; ...
!    :
!    ---------------------------------------------------------------------
!
!  NOTE: time in CET, not UTC
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


module LEKF_meas_ground

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  use NetCDF, only : NF90_StrError, NF90_NOERR
  
  use LE_Output_Common, only : T_LE_Output_Common
  
  use LEKF_Dims, only : maxmeas
  use LEKF_Meas, only : TGainBoundInfo
  
  use LML_Data_File, only : T_LML_Data_File
  
  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  InitMeas_Ground, DoneMeas_Ground
  public  ::  GetMeas_Ground, MeasUpdate_Ground
  public  ::  CalcMeas_Ground
 
  public  ::  meas_output_ground
  public  ::  Init, Done, PutOut
  
  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'kf_meas_ground'
  
  ! screening factor: measurements are rejected
  ! if square of departure exceeds factor times variance:
  real, parameter   ::  screening_factor = 5.0
    
  ! data sets:
  integer, parameter            ::  max_set = 2
  
  !! the number of measurements in the EMEP file
  !integer, parameter :: emepmeas = 84

  ! name lengths:
  integer, parameter            ::  codelen = 7
  integer, parameter            ::  namelen = 40
  integer, parameter            ::  complen = 8
  integer, parameter            ::  unitlen = 8
  
  !! maximum number of observed components:
  !integer, parameter            ::  maxcomp
  
  ! filter states to be written:
  integer, parameter            ::  nfs = 3
  character(len=*), parameter   ::  fsname(nfs) = (/'xb','x ','s '/)
  
  ! filter moments:
  integer, parameter            ::  nfm = 2
  integer, parameter            ::  ifm_forecast = 1
  integer, parameter            ::  ifm_analysis = 2
  character(len=*), parameter   ::  fm_descr(nfm) = &
        (/ 'after filter forecast', &
           'after filter analysis' /)
  
  ! Status number is sum of one or more numbers 2^n,
  ! for example   1+4+8
  ! To test if it is composed of a certain number, use:
  !   if ( iand(status,4) /= 0 )  ..
  !   
  integer, parameter            ::  status_default    = 0
  integer, parameter            ::  status_outside    = 1
  integer, parameter            ::  status_nodata     = 2
  integer, parameter            ::  status_validation = 4
  integer, parameter            ::  status_screened   = 8
  character(len=*), parameter   ::  status_description = &
      'status flag, 0=default, +1=outside-domain, +2=no-data, +4=validation, +8=screened'

  
  ! --- types ----------------------------------
  
  type T_Ground_Set
    ! set name:
    character(len=16)       ::  name
    ! type of data set: lml, ...
    character(len=4)        ::  typ
    ! range of observations:
    integer                 ::  meas1
    integer                 ::  meas2
    ! interface to LML set:
    character(len=256)      ::  lml_datadir
    character(len=256)      ::  lml_datafile
    type(T_LML_Data_File)   ::  lml_data_file
    character(len=8)        ::  lml_specname
    character(len=8)        ::  lml_specunit
    logical                 ::  lml_analyse
    type(TGainBoundInfo)    ::  lml_gbi
    real                    ::  lml_r_frac
    real                    ::  lml_r_min
    real                    ::  lml_r_max
  end type T_Ground_Set
  
  ! *
   
  type T_Ground_Meas
    ! data set and observation number within set:
    integer                   ::  iset
    integer                   ::  iobs
    !integer                   ::  iobs
    ! the name and code of the station
    character(len=namelen)    ::  name
    character(len=codelen)    ::  code
    ! its coordinates
    real                      ::  lon, lat, alt
    ! the height of the inlet
    real                      ::  href
    ! le spec:
    character(len=complen)    ::  specname
    character(len=unitlen)    ::  specunit
    integer                   ::  ix, iy
    integer                   ::  is
    ! its position in the statevector
    integer                   ::  inds(0:7)
    ! radius of influence to be applied (km.)7
    real                      ::  rho
    ! analyse or for validation only ?
    logical                   ::  analyse
    ! data accepted or nan etc ?
    integer                   ::  status
    ! about the measurement: the value + standard deviation
    real                      ::  c, sd
  end type T_Ground_Meas

  ! output stuff
  
  type meas_output_ground
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
    integer                     ::  dimid_complen
    integer                     ::  dimid_unitlen
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
    integer                     ::  varid_stn_comp
    integer                     ::  varid_stn_unit
!    integer                     ::  varid_stn_ana
!    integer                     ::  varid_lon
!    integer                     ::  varid_lat
    ! tracer variables:
!    real                        ::  unitconv
    integer                     ::  varid_status
    integer                     ::  varid_y
    integer                     ::  varid_r
    integer                     ::  varid_c(nfs,nfm)
  end type meas_output_ground


  ! --- interfaces -------------------------
  
  interface Init
    module procedure meas_output_ground_Init
  end interface
  
  interface Done
    module procedure meas_output_ground_Done
  end interface

  interface PutOut
    module procedure meas_output_ground_PutOut
  end interface

  
  ! --- var --------------------------------
  
  ! data set description:
  type(T_Ground_Set)          ::  dset(max_set)
  integer                     ::  n_set
  
  !! analyse this set ?
  !logical                   ::  analyse

!  ! correlation cut-off:
!  real                        ::  lml_rho
!  type(TGainBoundInfo)        ::  gbi
  
  ! declaration of the meas set
  type(T_Ground_Meas)         ::  meas(maxmeas)
  
  ! actual number read:
  integer                     ::  nmeas
  
!  integer   ::  meas_record
!  integer   ::  out_unit
!  integer   ::  out_record

  ! single tracer only, thus define globaly
  integer                     ::  icomp  

!  ! file name and unit:
!  character(len=256)          ::  fname
!  integer                     ::  fu


contains


  ! ======================================================================
  ! ===
  ! === input, update
  ! ===
  ! ======================================================================
  
  
  ! opens the measurement file (if present)
  ! end sets some parameters

  subroutine InitMeas_Ground( rcfile, t, status )

    use GO, only : TrcFile, Init, Done, ReadRc
    use GO, only : goGetFU
    use GO, only : TDate, NewDate, Get
    use GO, only : goReadFromLine, goMatchValue
    
!    use dims      , only : nx, ny, runF, locF
    use dims      , only : nspec
!    use constants , only : cum_days
!    use units     , only : u_tmp
    !use sysdep
!    use utils     , only : indomain
    use indices   , only : specname

!    use LEKF_state  , only :  m_s
    use LEKF_state  , only :  substate_m
    !use LEKF_meas, only : meas, nmeas, emepmeas
    !use LEKF_meas, only : meas_unit, u_meas_min
    use LEKF_meas, only : Init
    
    use lml_statlist_file, only : T_LML_Statlist_File, Init, Done, ReadRecord

    ! --- in/out --------------------------------
    
    character(len=*), intent(in)    ::  rcfile
    type(TDate), intent(in)         ::  t
    integer, intent(out)            ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/initmeas_ground'
    
    ! --- local ----------------------------------
    
    type(TRcFile)               ::  rcF
    character(len=256)          ::  setnames
    character(len=16)           ::  setname
    integer                     ::  iset
    
    character(len=256)          ::  fname
    
    integer                     ::  icomps
    integer                     ::  stn_nspec
    character(len=8)            ::  stn_specnames(nspec)
    integer                     ::  stn_ispecs(nspec)
    
    type(T_LML_Statlist_File)   ::  lml_statlist
    character(len=40)           ::  lml_station

    integer                     ::  ispec
    logical                     ::  analyse
    real                        ::  rho
    real                        ::  fmm(3)

!    character(len=256)    ::  fdir
!    integer               ::  yy, mm
!    character(len=6)      ::  period
!    logical               ::  exist
!    integer               ::  istn, nstn
!    character(len=16)     ::  name
!    real                  ::  lon, lat, alt
!    integer               ::  irec
!    integer               ::  ix, iy
!    logical               ::  in_domain
!    integer               ::  recl
!    real(4)               ::  cc(emepmeas)
    
    ! --- begin ----------------------------------
    
    write (gol,'("KF:     setup LML measurements ...")'); call goPr

    ! open rcfile:    
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)
    
    ! data set names:
    call ReadRc( rcF, 'kf.meas.ground.sets', setnames, status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over data sets:
    iset = 0
    do

      ! end ?
      if ( len_trim(setnames) == 0 ) then
        n_set = iset
        exit
      end if

      ! next number:
      iset = iset + 1

      ! extract name:
      call goReadFromLine( setnames, setname, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      
      ! store:
      dset(iset)%name = trim(setname)
      
      ! set type:
      call ReadRc( rcF, 'kf.meas.ground.'//trim(setname)//'.type', dset(iset)%typ, status )
      IF_NOTOK_RETURN(status=1)
      
      ! read other settings based on type:
      select case ( dset(iset)%typ )
      
        case ( 'lml' )

          ! name of file with station list:
          call ReadRc( rcF, 'kf.meas.ground.'//trim(setname)//'.statlist', fname, status )
          IF_NOTOK_RETURN(status=1)
          ! measured spec:
          call ReadRc( rcF, 'kf.meas.ground.'//trim(setname)//'.specname', dset(iset)%lml_specname, status )
          IF_NOTOK_RETURN(status=1)
          ! data directory:
          call ReadRc( rcF, 'kf.meas.ground.'//trim(setname)//'.datadir', dset(iset)%lml_datadir, status )
          IF_NOTOK_RETURN(status=1)
          ! analyse ?
          call ReadRc( rcF, 'kf.meas.ground.'//trim(setname)//'.analyse', analyse, status, default=.false. )
          IF_NOTOK_RETURN(status=1)
          ! cut-off radius; by default infinite:
          call ReadRc( rcF, 'kf.meas.ground.'//trim(setname)//'.rho', rho, status, default=huge(1.0) )
          IF_NOTOK_RETURN(status=1)
          
          ! representation error frac/min/max :
          call ReadRc( rcF, 'kf.meas.ground.'//trim(setname)//'.fmm', fmm, status )
          IF_NOTOK_RETURN(status=1)
          ! store:
          dset(iset)%lml_r_frac = fmm(1)
          dset(iset)%lml_r_min  = fmm(2)
          dset(iset)%lml_r_max  = fmm(3)
          
          ! compare component names with tracer names:
          call goMatchValue( dset(iset)%lml_specname, specname, ispec, status )
          IF_NOTOK_RETURN(status=1)
          
          ! store ...
          dset(iset)%lml_analyse = analyse

          ! start of measurements for this set:
          dset(iset)%meas1 = nmeas+1
          dset(iset)%meas2 = -1  ! dummy
          
          ! open station list file:
          call Init( lml_statlist, fname, status )
          IF_NOTOK_RETURN(status=1)
          ! loop over records:
          do
            ! read next record:
            call ReadRecord( lml_statlist, lml_station, status )
            if (status<0) exit  ! eof
            IF_NOTOK_RETURN(status=1)
            ! check ...
            if ( nmeas >= maxmeas ) then
              write (gol,'("all meas filled; please increase maxmeas")'); call goErr
              TRACEBACK; status=1; return
            end if
            ! next meas:
            nmeas = nmeas + 1
            ! setup observation
            meas(nmeas)%iset     = iset
            meas(nmeas)%name     = trim(lml_station)
            meas(nmeas)%code     = 'LML'   ! dummy
            meas(nmeas)%lon      = 0.0     ! dummy
            meas(nmeas)%lat      = 0.0     ! dummy
            meas(nmeas)%alt      = 0.0     ! dummy
            meas(nmeas)%href     = 2.5     ! dummy
            meas(nmeas)%specname = trim(dset(iset)%lml_specname)
            meas(nmeas)%specunit = '-'     ! dummy
            meas(nmeas)%is       = ispec
            meas(nmeas)%inds     = -1      ! dummy
            meas(nmeas)%analyse  = analyse
            meas(nmeas)%status   = status_default + status_nodata + status_outside
            meas(nmeas)%c        = 0.0     ! dummy
            meas(nmeas)%sd       = 0.0     ! dummy
          end do ! stations
          ! close:
          call Done( lml_statlist, status )
          IF_NOTOK_RETURN(status=1)
          
          ! set index of final measurement:
          dset(iset)%meas2 = nmeas
          
          ! pre-compute gain boounds if necessary:
          if ( dset(iset)%lml_analyse ) then
            write (gol,'("KF:       pre-compute correlation ranges ...")'); call goPr
            call Init( dset(iset)%lml_gbi, rho, status )
            IF_NOTOK_RETURN(status=1)
          end if
          
          ! no data file opened yet:
          dset(iset)%lml_datafile = 'none'

        case default
        
          write (gol,'("unsupported ground set type : ",a)') trim(dset(iset)%typ); call goErr
          TRACEBACK; status=1; return
      
      end select

    end do  ! data sets
   
!    !
!    ! ~~ tracers
!    !
!    
!    ! set tracer to be ozone
!    icomp = i_o3
!    
!    
!
!    !
!    ! ~~ station locations
!    !
!    
!    ! nothing stored yet ...
!    nmeas = 0
!    
!    write (gol,'("KF:     stn name lon lat alt irec indomain")'); call goPr
!
!    ! name with station locations:   
!    write (fname,'(a,"/statlist.txt")') trim(fdir)
!    
!    ! check ...
!    inquire( file=fname, exist=exist )
!    if ( .not. exist ) then
!      write (gol,'("file not found:")'); call goErr
!      write (gol,'("  ",a)') trim(fname); call goErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! select file unit:
!    call goGetFU( fu, status )
!    IF_NOTOK_RETURN(status=1)
!    ! adhoc, to avoid interfernece with LE units ..
!    fu = fu + 456
!    
!    ! open:
!    open( unit=fu, file=trim(fname), status='old', &
!            form='formatted', iostat=status )
!    if (status/=0) then
!      write (gol,'("from opening file:")'); call goErr
!      write (gol,'("  file name : ",a)') trim(fname); call goErr
!      write (gol,'("  iostat    : ",i6)') status; call goErr
!      TRACEBACK; status=1; return
!    end if
!
!    ! read the number of measurements
!    read (fu,*,iostat=status) nstn
!    if (status/=0) then
!      write (gol,'("from reading number of stations from file:")'); call goErr
!      write (gol,'("  file name : ",a)') trim(fname); call goErr
!      write (gol,'("  iostat    : ",i6)') status; call goErr
!      TRACEBACK; status=1; return
!    end if
!
!    ! read all locations
!    do istn = 1, nstn
!    
!      ! read line:
!      read(fu,*,iostat=status) name, lon, lat, alt, irec
!      if (status/=0) then
!        write (gol,'("from reading station line:")'); call goErr
!        write (gol,'("  file name : ",a)') trim(fname); call goErr
!        write (gol,'("  station   : ",i6)') istn; call goErr
!        write (gol,'("  iostat    : ",i6)') status; call goErr
!        TRACEBACK; status=1; return
!      end if
!
!      ! check if meas is in domain;
!      ! return cell indices (ix,iy) and flag:
!      call InDomain( lon, lat, runF%dlon, runF%dlat, runF%southb, runF%westb, &
!                      ix, iy, in_domain )
!
!      ! info ...
!      write (gol,'("KF:     ",i3," ",a,3f8.2,i4," ",l2)') &
!              istn, trim(name), lon, lat, alt, irec,&
!              in_domain; call goPr
!
!      ! not in domain ? then skip ...
!      if ( .not. in_domain ) cycle
!      
!      ! above 700 m ? then skip ...
!      if ( alt > 700.0 ) cycle
!
!      ! increase counter:
!      nmeas = nmeas + 1
!
!      ! check ...
!      if ( nmeas > maxmeas ) then
!        write (gol,'("total number of stations in data sets exceeds maxmeas = ",i6)') maxmeas; call goErr
!        TRACEBACK; status=1; return
!      end if
!
!      ! store:
!      ! the name and code of the station
!      meas(nmeas)%code      = trim(name)
!      meas(nmeas)%name      = trim(name)    ! no long name available, use code ...
!      ! its coordinates
!      meas(nmeas)%lon       = lon
!      meas(nmeas)%lat       = lat
!      meas(nmeas)%alt       = alt
!      ! position in data record:
!      meas(nmeas)%irec      = irec
!      ! flag if station is in current domain
!      meas(nmeas)%indomain  = in_domain
!      ! the grid indices in the 4D array
!      meas(nmeas)%ix        = ix
!      meas(nmeas)%iy        = iy
!      meas(nmeas)%iz        = 0  ! surface concentrations
!      meas(nmeas)%is        = icomp
!      
!      ! the index of the measurement in the measurement part 
!      ! of the state vector:
!!      meas(nmeas)%ind        = m_s-1 + nmeas
!      meas(nmeas)%inds(0:1)   = (/ substate_m, nmeas /)
!
!      ! the height of the inlet
!      meas(nmeas)%href      = 0.0            ! dummy
!
!      ! radius of influence to be applied (km.)
!      meas(nmeas)%rho       = lml_rho   ! km
!
!      ! analyse this site ?
!      meas(nmeas)%analyse   = analyse
!
!      !! not accepted by default ...
!      !meas(nmeas)%accepted  = .false.
!      ! accepted by default ...
!      meas(nmeas)%status = 0
!
!      ! no data yet ..
!      meas(nmeas)%c         = 0.0
!      meas(nmeas)%sd        = 0.0
!
!    end do
!
!    ! close:
!    close( unit=fu, iostat=status )
!    if (status/=0) then
!      write (gol,'("from closing file:")'); call goErr
!      write (gol,'("  file name : ",a)') trim(fname); call goErr
!      write (gol,'("  iostat    : ",i6)') status; call goErr
!      TRACEBACK; status=1; return
!    end if
!    
!
!    !
!    ! ~~ data file
!    !
!    
!    ! name with data:
!    write (fname,'(a,"/",a,"_hourly_ground.dat")') trim(fdir), trim(period)
!    
!    ! check ...
!    inquire( file=fname, exist=exist )
!    if ( .not. exist ) then
!      write (gol,'("file not found:")'); call goErr
!      write (gol,'("  ",a)') trim(fname); call goErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! get record length 
!    inquire( iolength=recl ) cc
!
!    ! open the meas file
!    !open(meas_unit,file='/user4/lotos/measurements/meas/grnd/emep/emep97.dat', &
!    !     form='unformatted', &
!    !     recl=emepmeas*reclen_fac_bin_io ,access='direct',status='old')
!    !open(meas_unit,file='/linuxMA5a/lotos/input/eumetsat2/summer_hourly_ground.dat', &
!    !     form='unformatted', recl=emepmeas*reclen_fac_bin_io ,access='direct',status='old') 
!
!    ! open file:
!    open( unit=fu, file=trim(fname), status='old', &
!            form='unformatted', access='direct', recl=recl, &
!            iostat=status )
!    if (status/=0) then
!      write (gol,'("from opening file:")'); call goErr
!      write (gol,'("  file name : ",a)') trim(fname); call goErr
!      write (gol,'("  iostat    : ",i6)') status; call goErr
!      TRACEBACK; status=1; return
!    end if
!    
!    !out_unit = u_meas_min+2
!    !open(out_unit,file=trim(locF%outpath)//'pm_stat.dat', &
!    !     form='unformatted', recl=1*reclen_fac_bin_io ,access='direct',status='new') 
!
!    !meas_record = (cum_days(mm) + dd -1)*24 +hh-1- (cumdays(2)+14)*24
!    !out_record = 0

    ! ok
    status = 0

  end subroutine InitMeas_Ground


  ! ***
  

  ! read station locations etc

  subroutine DoneMeas_Ground( status )

    use LML_Data_File, only : Done
    use LEKF_Meas      , only : Done

    ! --- in/out --------------------------------
    
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/DoneMeas_Ground'
    
    ! --- local -----------------------------
    
    integer               ::  iset

    ! --- begin ------------------------------
    
    ! loop over data sets:
    do iset = 1, n_set

      ! use appropriate routine:
      select case ( dset(iset)%typ )
      
        case ( 'lml' )
        
          ! dat file opened ? then close:
          if ( trim(dset(iset)%lml_datafile) /= 'none' ) then
            ! close:
            call Done( dset(iset)%lml_data_file, status )
            IF_NOTOK_RETURN(status=1)
          end if

          ! done with gain bounds:
          if ( dset(iset)%lml_analyse ) then
            call Done( dset(iset)%lml_gbi, status )
            IF_NOTOK_RETURN(status=1)
          end if
    
        case default
        
          write (gol,'("unsupported ground set type : ",a)') trim(dset(iset)%typ); call goErr
          TRACEBACK; status=1; return
      
      end select

    end do  ! data sets
   
    ! ok
    status = 0

  end subroutine DoneMeas_Ground
  
  
  ! ***
  
  
  subroutine GetMeas_Ground( t, status )

#ifdef pgf90
    use GO
#else  
    use GO, only : TDate, Get, operator(<), operator(-), rTotal, wrtgol, Precisely
    use GO, only : goLoCase
#endif

    use LML_Data_File, only : Init, Done, GetObservation, SearchStation, FindRecord
    
    use dims , only : runF
    use utils, only : InDomain

    use LEKF_State, only : substate_m
    use LEKF_Meas , only : Done
    
    ! --- in/out ----------------------------
    
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GetMeas_Ground'
    
    ! --- local -----------------------------
    
    integer               ::  i
    integer               ::  iset
    integer               ::  imeas
    character(len=256)    ::  datafile
    integer               ::  lml_num
    character(len=8)      ::  lml_comp
    integer               ::  ix, iy
    logical               ::  in_domain
   
    ! --- begin -----------------------------

    ! fill with nan:
    do i = 1, nmeas
      ! dummy values:
      meas(i)%c  = -999.9
      meas(i)%sd = -999.9
      ! no data by default ...
      meas(i)%status = status_default + status_nodata
    end do

    ! used for validation only ?
    do i = 1, nmeas
      if ( .not. meas(i)%analyse ) then
        meas(i)%status = meas(i)%status + status_validation
      end if
    end do

    ! info ...    
    call wrtgol( 'KF:     load ground measurements for ', t ); call goPr
    
    ! loop over data sets:
    do iset = 1, n_set

      ! info ...    
      write (gol,'("KF:       data set : ",a)') trim(dset(iset)%name); call goPr
      
      ! use appropriate read routine:
      select case ( dset(iset)%typ )
      
        case ( 'lml' )
        
          ! hourly records only:
          if ( .not. Precisely(t,1.0,'hour') ) then
            call wrtgol( 'KF:         only hourly, not for ', t ); call goPr
            status=0; return
          end if

          ! name of data file:
          write (datafile,'(a,"/",a,"ned12m",i4,"o.spr")') &
              trim(dset(iset)%lml_datadir), goLoCase(trim(dset(iset)%lml_specname)), t%year
              
          ! different from current name ? then open:
          if ( trim(datafile) /= trim(dset(iset)%lml_datafile) ) then

            ! already a file opened ? then close:
            if ( trim(dset(iset)%lml_datafile) /= 'none' ) then
              ! close:
              call Done( dset(iset)%lml_data_file, status )
              IF_NOTOK_RETURN(status=1)
            end if

            ! store new name:
            dset(iset)%lml_datafile = trim(datafile)

            ! open new file:
            call Init( dset(iset)%lml_data_file, trim(dset(iset)%lml_datafile), status )
            IF_NOTOK_RETURN(status=1)

            ! loop over associated measurements to extract station info:
            do imeas = dset(iset)%meas1, dset(iset)%meas2

              ! column in lml files:
              call SearchStation( dset(iset)%lml_data_file, trim(meas(imeas)%name), &
                                    meas(imeas)%iobs, status )
              IF_NOTOK_RETURN(status=1)

              ! extract station stuff:
              call GetObservation( dset(iset)%lml_data_file, meas(imeas)%iobs, status, &
                                     comp     =lml_comp, &
                                     unit     =meas(imeas)%specunit, &
                                     longitude=meas(imeas)%lon, &
                                     latitude =meas(imeas)%lat, &
                                     station_numb=lml_num )
              IF_NOTOK_RETURN(status=1)
              
              ! fill code:
              write (meas(imeas)%code,'("LML",i3)') lml_num

              ! check ...
              if ( goLoCase(trim(lml_comp)) /= goLoCase(trim(meas(imeas)%specname)) ) then
                write (gol,'("component in lml file is not what is expected ...")'); call goErr
                write (gol,'("  expected   : ",a)') trim(meas(imeas)%specname); call goErr
                write (gol,'("  found      : ",a)') trim(lml_comp); call goErr
                TRACEBACK; status=1; return
              end if

              ! check if meas is in domain;
              ! return cell indices (ix,iy) and flag:
              call InDomain( meas(imeas)%lon, meas(imeas)%lat, &
                             runF%dlon, runF%dlat, runF%southb, runF%westb, &
                              ix, iy, in_domain )

              ! info ...
              write (gol,'("KF:     ",i3," ",a,3f8.2,i4," ",l2)') &
                      imeas, trim(meas(imeas)%name), &
                      meas(imeas)%lon, meas(imeas)%lat, 0.0, &
                      meas(imeas)%iobs, in_domain; call goPr

              ! fill location:
              if ( in_domain ) then
                ! fill dummy cell location:
                meas(imeas)%ix = ix
                meas(imeas)%iy = iy
              else
                ! fill dummy cell location:
                meas(imeas)%ix = -999
                meas(imeas)%iy = -999
                ! set outside-domain bit:
                meas(imeas)%status = meas(imeas)%status + status_outside
              end if
              
             ! the index of the measurement in the measurement 
             ! part of the state vector:
             meas(imeas)%inds(0:1) = (/ substate_m, imeas /)

            end do

          end if  ! measurements
    
          ! read record for requested time (if not done yet):
          call FindRecord( dset(iset)%lml_data_file, t, status )
          IF_NOTOK_RETURN(status=1)
          
          ! loop over associated measurements to extract station info:
          do imeas = dset(iset)%meas1, dset(iset)%meas2
            ! extract measured value:
            call GetObservation( dset(iset)%lml_data_file, meas(imeas)%iobs, status, value=meas(imeas)%c )
            IF_NOTOK_RETURN(status=1)
            ! valid data ?
            if ( meas(imeas)%c >= 0.0 ) then
              ! remove no-data bit:
              meas(imeas)%status = meas(imeas)%status - status_nodata
              ! fill representation error as bounded fraction of measured value:
              meas(imeas)%sd = min(max( dset(iset)%lml_r_min, dset(iset)%lml_r_frac*meas(imeas)%c ), dset(iset)%lml_r_max )
            else
              ! fill dummy representation error:
              meas(imeas)%sd = -9999.9
            end if
          end do  ! measurements
    
        case default
        
          write (gol,'("unsupported ground set type : ",a)') trim(dset(iset)%typ); call goErr
          TRACEBACK; status=1; return
      
      end select

    end do  ! data sets
   
    ! ok
    status = 0

  end subroutine GetMeas_Ground


  ! ***
  

  subroutine CalcMeas_Ground( cg, xm, status )

    use binas  , only : XM_O, XM_air
    use indices, only : specname, specunit
    use dims   , only : dens
    use LEKF_Dims, only : nx, ny, nspec, maxmeas

    ! --- in/out --------------------------------

    real, intent(in)                ::  cg(nx,ny,nspec)
    real, intent(out)               ::  xm(maxmeas)
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/CalcMeas_Ground'
    
    ! --- local ----------------------------------
    
    integer             ::  i
    integer             ::  ix, iy, iz, is
    character(len=32)   ::  conversion_name
    real                ::  XM_tracer
    
    ! --- begin ------------------------------
    
    write (gol,'("KF:     fill measurements ...")'); call goPr
    
    ! loop over measurements:
    do i = 1, nmeas

      ! outside domain ? skip ..
      if ( iand(meas(i)%status,status_outside) /= 0 ) then

        xm(i) = -999.9
      
      else
      
        ! shorthands ..
        ix = meas(i)%ix
        iy = meas(i)%iy
        iz = 1  ! first layer if density etc is required
        is = meas(i)%is
      
        ! extract:
        if ( trim(specunit(meas(i)%is)) == trim(meas(i)%specunit) ) then

          ! same units; copy:
          xm(i) = cg(ix,iy,is)

        else

          ! combine units in a name:
          conversion_name = trim(specunit(meas(i)%is))//' -> '//trim(meas(i)%specunit)
          ! convert
          select case ( conversion_name )

            case ( 'ppb -> ug/m3' )

              !      (mol tracer)/(mol air) (kg tracer)/(mol tracer)       ug/kg         (ug tracer)
              ! ppb ----------------------- ------------------------ ----------------- = -----------
              !                 ppb            (kg air)/(mol air)    (m3 air)/(kg air)    (m3 air)
              
              ! mole mass:
              select case ( specname(meas(i)%is) )
                case ( 'o3', 'O3' ) ; XM_tracer = XM_O * 3
                case default
                  write (gol,'("mole mass not implemented for tracer : ",a)') trim(specname(meas(i)%is)); call goErr
                  TRACEBACK; status=1; return
              end select

              ! convert:              
              !xm(i) = cg(ix,iy,is) * 1e-9 * XM_tracer / XM_air * 1e9 * dens(ix,iy,1)
              xm(i) = cg(ix,iy,is)        * XM_tracer / XM_air       * dens(ix,iy,1)

            case default
              write (gol,'("convserion not implemented yet : ",a)') trim(conversion_name); call goErr
              TRACEBACK; status=1; return

          end select

        end if
        
      end if
    
    end do  ! nmeas
    
    ! ok
    status = 0

  end subroutine CalcMeas_Ground
  
  
  ! ***
  
  
  ! implementation of measurement update
  ! routine expects a measurement list +
  ! a list specifying the corresponding element nr.
  ! of the statevector.

  subroutine MeasUpdate_Ground( t, nmodes, status )

    use GO     , only : TDate, Midnight
    use LEKF_Dims, only : maxmodes
    use LEKF_meas, only : meas_update_scalar

    ! --- in/out -----------------------------

    type(TDate), intent(in)   ::  t    
    integer, intent(in)       ::  nmodes     ! actual number of modes
    integer, intent(out)      ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MeasUpdate_Ground'

    ! measurements are 1D    
    integer, parameter  ::  prf_nz = 1

    ! --- local ------------------------------
    
    integer       ::  i
    integer       ::  iset
    real          ::  prf_yr(prf_nz)
    !real          ::  prf_ya(prf_nz)
    !real          ::  prf_A(prf_nz,prf_nz)
    real          ::  prf_R(prf_nz,prf_nz)
          
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
    write (gol,'("KF:   measurement update ground observations ...")'); call goPr
    
    ! loop over measurements:
    do i = 1, nmeas

             !write (gol,'("KF:   aaa ",2i6)') i, meas(i)%status; call goPr

      ! skip if not accepted (out of domain, no data, validation, ..):
      if ( meas(i)%status /= status_default ) cycle

      ! hoezee!
      write (gol,'("KF:   analyse ",i3)') i; call goPr
      
      ! data set number:
      iset = meas(i)%iset
      
      ! analyse given data set type:
      ! use appropriate read routine:
      select case ( dset(iset)%typ )
      
        case ( 'lml' )
        
          ! dummy profile ...
          prf_yr = meas(i)%c      ! retrieved
          !prf_ya = 0.0            ! a-priori
          !prf_A  = 1.0            ! kernel
          prf_R  = meas(i)%sd**2  ! observation error covariance

          ! analyse scalar measurment:
          call meas_update_scalar( prf_yr(1), prf_R(1,1), &
                                    meas(i)%inds, meas(i)%ix, meas(i)%iy, &
                                    dset(iset)%lml_gbi, screening_factor, &
                                    nmodes, status )
          ! screening: innovation too large compared with expected variance:
          if ( status == -1 ) then
            meas(i)%status = meas(i)%status + status_screened
          else if ( status /= 0 ) then
            TRACEBACK; status=1; return
          end if

        case default
        
          write (gol,'("unsupported ground set type : ",a)') trim(dset(iset)%typ); call goErr
          TRACEBACK; status=1; return
      
      end select

    end do
    
    ! ok
    status = 0

  end subroutine measupdate_ground


  ! ======================================================================
  ! ===
  ! === output
  ! ===
  ! ======================================================================
  
  
  subroutine meas_output_ground_Init( kfmo, rcfile, rckey, status )

    use GO     , only : AnyDate
    use LE_Output_Common, only : Init
  
    ! --- in/out --------------------------------
    
    type(meas_output_ground), intent(out)   ::  kfmo
    character(len=*), intent(in)            ::  rcfile
    character(len=*), intent(in)            ::  rckey
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/meas_output_ground_Init'
    
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
    
  end subroutine meas_output_ground_Init
  
  
  ! ***
  

  subroutine meas_output_ground_Done( kfmo, status )
  
    use NetCDF          , only : NF90_Close
    use LE_Output_Common, only : Done
  
    ! --- in/out --------------------------------
    
    type(meas_output_ground), intent(inout)   ::  kfmo
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/meas_output_ground_Done'
    
    ! --- begin ---------------------------------
    
    ! file opened ?
    if ( kfmo%opened ) then
      ! close:
      status = NF90_Close( kfmo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! reset flag:
      kfmo%opened = .false.
    end if
    
    ! done with common stuff ...
    call Done( kfmo%com, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine meas_output_ground_Done
  
  
  ! ***
  
  
  subroutine meas_output_ground_PutOut( kfmo, key, t, status )

#ifdef pgf90
    use GO_Date
#else  
    use GO, only : TDate, NewDate, IncrDate, AnyDate, Get
    use GO, only : operator(+), operator(-), operator(<), operator(>)
    use GO, only : rTotal, Precisely
    use GO, only : wrtgol
#endif

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE
    
!#ifdef pgf90
!    use JAQL_Output_CF_Conventions
!#else  
!    use JAQL   , only : JAQL_Output_CF_names
!#endif

    use Dims   , only : runF
    use Indices, only : specname, specunit
    use LE_Output_Common, only : PutOut_GlobalAttributes
    
    use LEKF_State, only : kf_with_xb, kf_with_xm
    use LEKF_State, only : xb, x, sigma
!    use LEKF_State, only : m_s, m_e
 
    ! --- in/out --------------------------------
    
    type(meas_output_ground), intent(inout)   ::  kfmo
    character(len=*), intent(in)              ::  key
    type(TDate), intent(in)                   ::  t
    integer, intent(out)                      ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/kfmo_meas_PutOut'
    
    ! --- local ---------------------------------

    logical               ::  newfile    
    integer               ::  time6(6)
    real                  ::  time

    integer               ::  cmode

!    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units
!    character(len=512)    ::  comment
    
    integer               ::  varid
    type(TDate)           ::  t0

    integer               ::  ifs
    character(len=16)     ::  varname
    character(len=512)    ::  descr

    integer               ::  ifm
    
    integer               ::  i
    real                  ::  cc(maxmeas)
    
    ! --- begin ---------------------------------
    
    ! check ...
    if ( nmeas <= 0 ) then
      write (gol,'(a," - WARNING - no measurements available ...")') trim(rname); call goPr
      status=0; return
    end if
    
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
    
    ! hourly only ...
    if ( .not. Precisely(t,1.0,'hour') ) then
      status=0; return
    end if
    
    ! info ...
    call wrtgol('KF: put out ground observations for : ',t); call goPr
    
    ! create new file ?
    if ( newfile ) then

      ! extract time fields for day assigned to file:
      call Get( t, time6=time6 )

      ! set time range [00,24) for this day:
      ! NOTE: measurements are probably averaged over the past interval,
      !       but this is now ignored until the model could produce time
      !       averaged concentrations too ...
      kfmo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
      kfmo%tr(2) = kfmo%tr(1) + IncrDate( day=1 ) - IncrDate(mili=1)

      ! new file name:
      write (kfmo%fname,'(a,a,"_",a,"_",a,"_",i4.4,2i2.2,".nc")') &
                trim(kfmo%com%outdir), trim(kfmo%com%model), trim(kfmo%com%expid), &
                'meas-ground', time6(1:3)

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
      status = NF90_Def_Dim( kfmo%ncid, 'complen', complen, kfmo%dimid_complen )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'unitlen', unitlen, kfmo%dimid_unitlen )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'meas', nmeas, kfmo%dimid_meas )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'time', NF90_UNLIMITED, kfmo%dimid_time )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfmo%ncid, 'datelen', 6, kfmo%dimid_datelen )
      IF_NF90_NOTOK_RETURN(status=1)

      ! define variables:

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

      status = NF90_Def_Var( kfmo%ncid, 'station_ix', NF90_INT, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'description', 'lotos-euros grid cell in x' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_ix = varid

      status = NF90_Def_Var( kfmo%ncid, 'station_iy', NF90_INT, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'description', 'lotos-euros grid cell in y' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_iy = varid
      
      status = NF90_Def_Var( kfmo%ncid, 'station_alt', NF90_REAL, (/kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'standard_name', 'altitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'long_name', 'altitude' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'units', 'm' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_alt = varid

      varname = 'station_comp'
      status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_CHAR, (/kfmo%dimid_complen,kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'description', 'observed component' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_comp = varid

      varname = 'station_unit'
      status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_CHAR, (/kfmo%dimid_unitlen,kfmo%dimid_meas/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( kfmo%ncid, varid, 'description', 'observation unit' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_stn_unit = varid

      ! observations
      
      varname = 'status'
      status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_INT, (/kfmo%dimid_meas,kfmo%dimid_time/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( kfmo%ncid, varid, 'description', trim(status_description) )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_status = varid

      varname = 'conc_y'
      status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_REAL, (/kfmo%dimid_meas,kfmo%dimid_time/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( kfmo%ncid, varid, 'description', 'observed' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_y = varid

      varname = 'conc_r'
      status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_REAL, (/kfmo%dimid_meas,kfmo%dimid_time/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( kfmo%ncid, varid, 'description', 'representation error std.dev.' )
      IF_NF90_NOTOK_RETURN(status=1)
      kfmo%varid_r = varid

!      ! global tracer index
!      itr = icomp
!
!      ! CF standard name for concentration/mixing ratio/column:
!
!      ! no comment yet
!      comment = ''
!
!      ! get names following CF conventions;
!      ! store conversion factor for later usage:
!      call JAQL_Output_CF_names( &
!                   specname(itr), specunit(itr), &
!                   cf_standard_name, cf_long_name, cf_units, &
!                   kfmo%unitconv, comment, &
!                   status )
!      IF_NOTOK_RETURN(status=1)
      
      ! reset id's for safety:
      kfmo%varid_c = 0
      
      ! loop over filter moments (forecast, analysis)
      do ifm = 1, nfm

        ! loop over filter states (xb, x, s)
        do ifs = 1, nfs

          ! filter state description, skip if possible
          select case ( trim(fsname(ifs)) )
            case ( 'xb' )
              ! skip ?
              if ( .not. kf_with_xb ) cycle
              if ( ifm /= ifm_analysis ) cycle   ! xb is not changed by analyis, write for analysis
              ! description:
              write (descr,'("model background run")')
            case ( 'x' )
              ! skip ?
              if ( .not. kf_with_xm ) cycle
              ! description:
              write (descr,'("ensemble mean; ",a)') trim(fm_descr(ifm))
            case ( 's' )
              ! skip ?
              if ( .not. kf_with_xm ) cycle
              ! description:
              write (descr,'("ensemble std.dev.; ",a)') trim(fm_descr(ifm))
            case default
              write (gol,'("unsupported filter state : ",a)') trim(fsname(ifs)); call goErr
              TRACEBACK; status=1; return
          end select

          ! variable name:
          write (varname,'("conc_",a)') trim(fsname(ifs))

          ! extend for forecast moment:
          if ( ifm == ifm_forecast ) varname = trim(varname)//'_f'

          ! define veriable:
          status = NF90_Def_Var( kfmo%ncid, trim(varname), NF90_REAL, (/kfmo%dimid_meas,kfmo%dimid_time/), varid )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Put_Att( kfmo%ncid, varid, 'description', &
                     'LOTOS-EUROS simulated concentration; '//trim(descr) )
          IF_NF90_NOTOK_RETURN(status=1)
          kfmo%varid_c(ifs,ifm) = varid

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
        status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_name, trim(meas(i)%name), &
                                 start=(/1,i/), count=(/len_trim(meas(i)%name),1/) )
        IF_NF90_NOTOK_RETURN(status=1)
      end do
      do i = 1, nmeas
        status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_code, trim(meas(i)%code), &
                                 start=(/1,i/), count=(/len_trim(meas(i)%code),1/) )
        IF_NF90_NOTOK_RETURN(status=1)
      end do
      do i = 1, nmeas
        status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_comp, trim(meas(i)%specname), &
                                 start=(/1,i/), count=(/len_trim(meas(i)%specname),1/) )
        IF_NF90_NOTOK_RETURN(status=1)
      end do
      do i = 1, nmeas
        status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_unit, trim(meas(i)%specunit), &
                                 start=(/1,i/), count=(/len_trim(meas(i)%specunit),1/) )
        IF_NF90_NOTOK_RETURN(status=1)
      end do
      ! write non-character fields:
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_lon, meas(1:nmeas)%lon )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_lat, meas(1:nmeas)%lat )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_alt, meas(1:nmeas)%alt )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_ix, meas(1:nmeas)%ix )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_stn_iy, meas(1:nmeas)%iy )
      IF_NF90_NOTOK_RETURN(status=1)

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
    
    ! write status:
    status = NF90_Put_Var( kfmo%ncid, kfmo%varid_status, meas(1:nmeas)%status, &
                     start=(/1,kfmo%itrec/), count=(/nmeas,1/) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write observations:
    status = NF90_Put_Var( kfmo%ncid, kfmo%varid_y, meas(1:nmeas)%c, &
                     start=(/1,kfmo%itrec/), count=(/nmeas,1/) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Var( kfmo%ncid, kfmo%varid_r, meas(1:nmeas)%sd, &
                     start=(/1,kfmo%itrec/), count=(/nmeas,1/) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! which filter moment ?
    select case ( key )
      case ( 'forecast' ) ; ifm = 1
      case ( 'analysis' ) ; ifm = 2
      case default
        write (gol,'("unsupported key : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! loop over filter states:
    do ifs = 1, nfs
    
      ! init all concentrations to zero:
      cc = 0.0

      ! extract simulated concentrations:
      select case ( trim(fsname(ifs)) )
        case ( 'xb' )
          if ( .not. kf_with_xb ) cycle
          if ( ifm /= ifm_analysis ) cycle   ! xb is not changed by analyis, write for analysis
          cc(1:nmeas) = xb%m(1:nmeas)
        case ( 'x' )
          if ( .not. kf_with_xm ) cycle
          cc(1:nmeas) = x%m(1:nmeas)
        case ( 's' )
          if ( .not. kf_with_xm ) cycle 
          cc(1:nmeas) = sigma%m(1:nmeas)
        case default
          write (*,'("unsupported filter state name : ",a)') trim(fsname(ifs)); call goErr
          TRACEBACK; status=1; return
      end select

      ! write concentrations:
      status = NF90_Put_Var( kfmo%ncid, kfmo%varid_c(ifs,ifm), cc(1:nmeas), &
                     start=(/1,kfmo%itrec/), count=(/nmeas,1/) )
      IF_NF90_NOTOK_RETURN(status=1)
      
    end do  ! filter states
    
    ! ok
    status = 0
    
  end subroutine meas_output_ground_PutOut


end module
