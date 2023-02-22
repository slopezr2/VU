!#######################################################################
!
! NAME
!
!   MAORI - Model And Output Routine Interface
!
! DESCRIPTION
!
!   Mode independent output routines.
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_Set_Sample

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate
  use GO, only : T_ISTA

  use MAORI_Param, only : MAORI_LEN_NAME, MAORI_LEN_FILE, MAORI_LEN_LINE
  use MAORI_LocList, only : T_MAORI_LocList
  use MAORI_VarList, only : T_MAORI_VarList

  use MAORI_DataFile_EMEP, only : T_MAORI_DataFile_EMEP
  
#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  implicit none

  ! --- in/out ---------------------------

  private

  public  ::  T_MAORI_Set_Sample_Data
  public  ::  MAORI_Set_Sample_Data_Init, MAORI_Set_Sample_Data_Done
  public  ::  MAORI_Set_Sample_Data_Start, MAORI_Set_Sample_Data_Setup
  public  ::  MAORI_Set_Sample_Data_Inq
  public  ::  MAORI_Set_Sample_Data_Param_Inq, MAORI_Set_Sample_Data_Param_Put
  public  ::  MAORI_Set_Sample_Data_Obs_Get, MAORI_Set_Sample_Data_Obs_Put
  public  ::  MAORI_Set_Sample_Data_Var_Inq
  public  ::  MAORI_Set_Sample_Data_Loc_Inq

  public  ::  T_MAORI_Set_Sample_State, T_MAORI_Set_Output_State
  public  ::  MAORI_Set_Sample_State_Init, MAORI_Set_Sample_State_Done
  public  ::  MAORI_Set_Sample_State_Start, MAORI_Set_Sample_State_Setup
  public  ::  MAORI_Set_Sample_State_Put, MAORI_Set_Sample_State_Post
  public  ::  MAORI_Set_Sample_State_Values_Get, MAORI_Set_Sample_State_Values_Put
  public  ::  MAORI_Set_Sample_State_LocValues_Get, MAORI_Set_Sample_State_LocValues_Put
  public  ::  MAORI_Set_Sample_State_Obs_Get

  public  ::  MAORI_Set_Output_State_Init, MAORI_Set_Output_State_Done
  public  ::  MAORI_Set_Output_State_Start
  public  ::  MAORI_Set_Output_State_Write


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI_Set_Sample'

  ! parameter indices:
  integer, parameter  ::  ipar_bound_west      = 1
  integer, parameter  ::  ipar_bound_east      = 2
  integer, parameter  ::  ipar_bound_south     = 3
  integer, parameter  ::  ipar_bound_north     = 4
  integer, parameter  ::  ipar_nlayer          = 5
  integer, parameter  ::  ipar_bias_correction = 6
  integer, parameter  ::  npar = 6

  ! --- types ------------------------------

  type T_MAORI_Set_Sample_Data
    ! name:
    character(len=MAORI_LEN_NAME)         ::  name
    ! flags:
    logical                               ::  started
    ! parameter flags:
    logical                               ::  param_set(npar)
    ! model parameters:
    real                                  ::  bounds(4)
    character(len=MAORI_LEN_LINE)         ::  bias_correction
    ! time resolution:
    character(len=5)                      ::  tres_period    ! 'day'  (in future: 'month', ...)
    !character(len=4)                      ::  tres_type      ! 'aver' (in future: 'inst', ..)
    character(len=6)                      ::  tres_unit       ! 'hour'
    real                                  ::  tres_start
    real                                  ::  tres_step
    type(TDate)                           ::  tres_t0
    type(TIncrDate)                       ::  tres_dt
    logical                               ::  tres_hist
    ! time step stuff:
    logical                               ::  the_end
    type(TDate)                           ::  setup_tr(2)
    type(TDate)                           ::  period_tr(2)
    logical                               ::  end_of_period
    integer                               ::  nstep
    integer                               ::  istep
    ! ~ history
    integer                               ::  nhist
    integer                               ::  ihist
    logical                               ::  new_period
    ! ~ averaging
    type(T_ISTA)                          ::  ista
    real                                  ::  w_aver, w_prev, w_curr
    logical                               ::  finished_aver
    ! vertical resolution:
    logical                               ::  profile
    integer                               ::  nlayer   ! in the model
    integer                               ::  maxlev   ! 1 or nlayer+1
    ! list with locations:
    character(len=MAORI_LEN_FILE)         ::  loc_query
    type(T_MAORI_LocList)                 ::  loclist
    integer                               ::  nloc
    ! list with input variables:
    type(T_MAORI_VarList)                 ::  obs_varlist
    integer                               ::  obs_nvar
    character(len=MAORI_LEN_NAME), pointer ::  obs_vname(:)   ! (obs_nvar)
    integer                               ::  obs_nlev
    character(len=MAORI_LEN_LINE)         ::  obs_query
    character(len=MAORI_LEN_FILE)         ::  obs_file
    logical                               ::  obs_opened
    type(T_MAORI_DataFile_EMEP), pointer  ::  obs_emep
    real, pointer                         ::  obs_value (:,:,:,:) ! nloc,obs_nlev,obs_nvar,nhist
    real, pointer                         ::  obs_errstd(:,:,:,:) ! nloc,obs_nlev,obs_nvar,nhist
    real                                  ::  obs_err_frac
    real                                  ::  obs_err_min_value
    character(len=MAORI_LEN_NAME)         ::  obs_err_min_units
    integer, pointer                      ::  obs_astat(:,:,:,:) ! nloc,obs_nlev,obs_nvar,nhist
    integer, pointer                      ::  obs_iloc2iobs(:)
    integer, pointer                      ::  obs_ivar(:)
    type(TDate)                           ::  obs_t
    ! list with variables:
    type(T_MAORI_VarList)                 ::  varlist
    integer                               ::  nvar
    integer, pointer                      ::  nlev(:)  ! nvar ; 1, nlayer, or nlayer+1
    ! assimilation flag:
    logical                               ::  assim_analyse
    real                                  ::  assim_rho
    real                                  ::  assim_alfa
    ! number of 'state' values:
    integer                               ::  nvalue
    integer                               ::  nlocvalue
  end type T_MAORI_Set_Sample_Data
  
  ! simulation from a state:
  type T_MAORI_Set_Sample_State
    ! id
    character(len=MAORI_LEN_NAME)         ::  name
    ! flags:
    logical                               ::  started
    ! simulated values:
    logical, pointer                      ::  sample_set (:,:,:)    ! nloc,maxlev,nvar
    real, pointer                         ::  sample_curr(:,:,:)    ! nloc,maxlev,nvar
    real, pointer                         ::  sample_prev(:,:,:)    ! nloc,maxlev,nvar
    real, pointer                         ::  sample_aver(:,:,:,:)  ! nloc,maxlev,nvar,nhist
  end type T_MAORI_Set_Sample_State
  
  ! output: file units etc
  type T_MAORI_Set_Output_State
    ! id
    character(len=MAORI_LEN_NAME)         ::  name
    ! file opened ?
    logical                               ::  opened
    ! time record counter:
    integer                               ::  itrec
#ifdef with_netcdf
    ! file handle:
    integer                               ::  ncid
    ! dimension handles:
    integer                               ::  dimid_namelen
    integer                               ::  dimid_codelen
    integer                               ::  dimid_loc
    integer                               ::  dimid_slev
    integer                               ::  dimid_flev
    integer                               ::  dimid_hlev
    integer                               ::  dimid_time
    integer                               ::  dimid_datelen
    ! dimension variables:
    integer                               ::  varid_time, varid_time2
    integer                               ::  varid_date, varid_date2
    ! station variables:
    integer                               ::  varid_stn_name
    integer                               ::  varid_stn_code
    integer                               ::  varid_stn_lon
    integer                               ::  varid_stn_lat
    integer                               ::  varid_stn_alt
    ! simulated variables:
    integer, pointer                      ::  varid(:)
    ! observations:
    integer, pointer                      ::  varid_y(:)
    integer, pointer                      ::  varid_r(:)
    integer, pointer                      ::  varid_astat(:)
#endif
  end type T_MAORI_Set_Output_State


  ! --- interfaces -------------------------------
  
  interface MAORI_Set_Sample_Data_Param_Put
    module procedure MAORI_Set_Sample_Data_Param_Put_i
    module procedure MAORI_Set_Sample_Data_Param_Put_r
    module procedure MAORI_Set_Sample_Data_Param_Put_s
  end interface


contains


  ! ********************************************************************
  ! ***
  ! *** data
  ! ***
  ! ********************************************************************


  subroutine MAORI_Set_Sample_Data_Init( mdata, name, rcF, rckey, t0, status )

    use GO, only : TRcFile, ReadRc
    use GO, only : goVarValue
    use GO, only : TDate, NewDate, AnyDate, IncrDate, operator(*)
    use GO, only : ISTA_Init

    use MAORI_Param  , only : MAORI_LEN_LINE
    use MAORI_VarList, only : MAORI_VarList_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(out)        ::  mdata
    character(len=*), intent(in)                  ::  name
    type(TRcFile), intent(in)                     ::  rcF
    character(len=*), intent(in)                  ::  rckey
    type(TDate), intent(in)                       ::  t0
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_Data_Init'

    ! --- local -----------------------------------

    character(len=MAORI_LEN_LINE)       ::  key
    integer                             ::  maxstep
    integer                             ::  ivar
    character(len=MAORI_LEN_NAME)       ::  varname

    ! --- begin -----------------------------------

    ! store name:
    mdata%name = trim(name)

    ! read key that describes station list file:
    call ReadRc( rcF, trim(rckey)//'.query', mdata%loc_query, status )
    IF_NOTOK_RETURN(status=1)
    
    ! time resolution key:
    call ReadRc( rcF, trim(rckey)//'.tres', key, status )
    IF_NOTOK_RETURN(status=1)
    ! extract values:
    mdata%tres_period = 'day'
      call goVarValue( key, ';', 'period', '=', mdata%tres_period, status )
      IF_ERROR_RETURN(status=1)
    !mdata%tres_type = 'aver'
    !  call goVarValue( key, ';', 'type', '=', mdata%tres_type, status )
    !  IF_ERROR_RETURN(status=1)
    mdata%tres_unit = 'hour'
      call goVarValue( key, ';', 'unit', '=', mdata%tres_unit, status )
      IF_ERROR_RETURN(status=1)
    mdata%tres_step = 1.0
      call goVarValue( key, ';', 'step', '=', mdata%tres_step, status )
      IF_ERROR_RETURN(status=1)
    mdata%tres_hist = .false.
      call goVarValue( key, ';', 'hist', '=', mdata%tres_hist, status )
      IF_ERROR_RETURN(status=1)

    ! set time resolution start and increment:
    select case ( trim(mdata%tres_unit) )
      !~
      case ( 'hour' )
        ! maximum number of steps in a period:
        select case ( trim(mdata%tres_period) )
          case ( 'day' ) ; maxstep = 24
          case default
            write (gol,'("unsupported period to determine nstep : ",a)') trim(mdata%tres_period); call goErr
            TRACEBACK; status=1; return
        end select
        ! first guess of start hour:
        mdata%tres_start = t0%hour
        if ( mdata%tres_step == 24.0 ) mdata%tres_start = 0
          ! eventually overwrite:
          call goVarValue( key, ';', 'start', '=', mdata%tres_start, status )
          IF_ERROR_RETURN(status=1)
        ! fill start time:
        mdata%tres_t0 = NewDate( year=t0%year, month=t0%month, day=t0%day, hour=nint(mdata%tres_start) )
        ! fill time step:
        if ( modulo(real(maxstep),mdata%tres_step) /= 0.0 ) then
          write (gol,'("tres hourly only support integer number of steps within ",i4," hour:")') mdata%nstep; call goErr
          write (gol,'("  tres (hours)    : ",f8.2)') mdata%tres_step; call goErr
          TRACEBACK; status=1; return
        end if
        mdata%tres_dt = IncrDate(hour=1) * mdata%tres_step
        ! number of steps within a period:
        mdata%nstep = nint(real(maxstep)/mdata%tres_step)
      !~
      case ( 'day'  )
        ! maximum number of steps in a period:
        select case ( trim(mdata%tres_period) )
          case ( 'day' ) ; maxstep = 1
          case default
            write (gol,'("unsupported period to determine nstep : ",a)') trim(mdata%tres_period); call goErr
            TRACEBACK; status=1; return
        end select
        ! fill start time:
        mdata%tres_t0 = NewDate( year=t0%year, month=t0%month, day=t0%day, hour=0 )
        ! fill time step:
        mdata%tres_dt = IncrDate(day=1) * mdata%tres_step
        ! number of steps within a period:
        mdata%nstep = nint(real(maxstep)/mdata%tres_step)
      !~
      case default
        write (gol,'("unsupported tres unit : ",a)') trim(mdata%tres_unit); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! store history ?
    if ( mdata%tres_hist ) then
      ! store history for whole file:
      mdata%nhist = mdata%nstep
      ! replace if explicitly defined in settings:
      call goVarValue( key, ';', 'nhist', '=', mdata%nhist, status )
      IF_ERROR_RETURN(status=1)
      ! init current value:
      mdata%ihist = 0
    else
      ! no history, only one value:
      mdata%nhist = 1
    end if

    ! setup sample averaging:
    call ISTA_Init( mdata%ista, mdata%tres_t0, mdata%tres_dt, status )
    IF_NOTOK_RETURN(status=1)
    ! init flags:
    mdata%finished_aver = .false.
    ! init time range:
    mdata%setup_tr(1) = t0
    mdata%setup_tr(2) = t0
    ! init period time range:
    mdata%period_tr(1) = AnyDate()
    mdata%period_tr(2) = AnyDate()
    
    ! profile ? otherwise surface only:
    call ReadRc( rcF, trim(rckey)//'.profile', mdata%profile, status )
    IF_NOTOK_RETURN(status=1)

    ! read space seperated list with variables:
    call ReadRc( rcF, trim(rckey)//'.var', key, status )
    IF_NOTOK_RETURN(status=1)
    ! store list:
    call MAORI_VarList_Init( mdata%varlist, key, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read space seperated list with input variables:
    call ReadRc( rcF, trim(rckey)//'.obs.var', key, status )
    IF_NOTOK_RETURN(status=1)
    ! store list:
    call MAORI_VarList_Init( mdata%obs_varlist, key, status )
    IF_NOTOK_RETURN(status=1)
    ! store number:
    mdata%obs_nvar = mdata%obs_varlist%n
    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
      ! data query:
      call ReadRc( rcF, trim(rckey)//'.obs.query', mdata%obs_query, status )
      IF_NOTOK_RETURN(status=1)
      ! variable names in observation files:
      allocate( mdata%obs_vname(mdata%obs_nvar), stat=status )
      IF_NOTOK_RETURN(status=1)
      do ivar = 1, mdata%obs_nvar
        ! current:
        varname = trim(mdata%obs_varlist%var(ivar)%name)
        ! read from settings:
        call ReadRc( rcF, trim(rckey)//'.obs.vname.'//trim(varname), &
                            mdata%obs_vname(ivar), status )
        IF_ERROR_RETURN(status=1)
      end do
      ! single level for the moment ...
      mdata%obs_nlev = 1
      ! error std.dev. fraction:
      call ReadRc( rcF, trim(rckey)//'.obs.err.frac', mdata%obs_err_frac, status )
      IF_NOTOK_RETURN(status=1)
      ! error std.dev. minimum:
      call ReadRc( rcF, trim(rckey)//'.obs.err.min', key, status, default='value=0.0;units=ppb' )
      IF_NOTOK_RETURN(status=1)
      call goVarValue( key, ';', 'value', '=', mdata%obs_err_min_value, status )
      IF_NOTOK_RETURN(status=1)
      call goVarValue( key, ';', 'units', '=', mdata%obs_err_min_units, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! analyse in assimilation ?
    call ReadRc( rcF, trim(rckey)//'.assim.analyse', mdata%assim_analyse, status )
    IF_NOTOK_RETURN(status=1)
    ! assimilation settings:
    if ( mdata%assim_analyse ) then
      ! check ...
      if ( mdata%obs_nvar < 1 ) then
        write (gol,'("could not analyse `",a,"` if no observations are available ...")') trim(name); call goErr
        TRACEBACK; status=1; return
      end if
      ! spatial correlation scale:
      call ReadRc( rcF, trim(rckey)//'.assim.rho', mdata%assim_rho, status )
      IF_NOTOK_RETURN(status=1)
      ! screening factor:
      call ReadRc( rcF, trim(rckey)//'.assim.alfa', mdata%assim_alfa, status )
      IF_NOTOK_RETURN(status=1)
    else
      ! dummy ...
      mdata%assim_rho  = 0.0   ! no correlation
      mdata%assim_alfa = 999.9   ! no screening ...
    end if

    ! no parameters set yet:
    mdata%param_set = .false.
    
    ! not started yet:
    mdata%started = .false.

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Init


  ! ***


  subroutine MAORI_Set_Sample_Data_Start( mdata, status )

    use GO, only : AnyDate
    use GO, only : goMatchValue, goVarValue
    use MAORI_Param  , only : MAORI_ASTAT_VALIDATION
    use MAORI_LocList, only : MAORI_LocList_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(inout)  ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_Data_Start'
    
    ! --- local ----------------------------------
    
    integer                         ::  ivar
    character(len=MAORI_LEN_NAME)   ::  obs_type

    ! --- begin -----------------------------------
    
    ! all parameters set ?
    if ( .not. all(mdata%param_set) ) then
      write (gol,'("not all parameters set: ")'); call goErr
      write (gol,'("  ",1000l2)') mdata%param_set
      TRACEBACK; status=1; return
    end if

    ! read locations; skip all outside bound box:
    call MAORI_LocList_Init( mdata%loclist, mdata%loc_query, status, &
                               bounds=mdata%bounds )
    IF_NOTOK_RETURN(status=1)
    
    ! extract dimensions:
    mdata%nloc = mdata%loclist%n
    mdata%nvar = mdata%varlist%n
    
    ! setup storage for number of output levels:
    allocate( mdata%nlev(mdata%nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! set maximum and actual number of output levels:
    do ivar = 1, mdata%nvar
      ! profile output or surface only ?
      if ( mdata%profile ) then
        ! all layers plus 1 extra for half levels
        mdata%maxlev = mdata%nlayer+1
        ! half level field ?
        if ( mdata%varlist%var(ivar)%halflevel ) then
          mdata%nlev(ivar) = mdata%nlayer + 1
        else
          mdata%nlev(ivar) = mdata%nlayer
        end if
      else
        ! single level output only:
        mdata%maxlev     = 1
        mdata%nlev(ivar) = 1
      end if
    end do
    
    ! number of state values on a single location:
    mdata%nlocvalue = 0
    do ivar = 1, mdata%nvar
      ! add contribution of this variable (times 2+nhist for 'curr', 'prev', and 'aver'):
      mdata%nlocvalue = mdata%nlocvalue + mdata%nlev(ivar) * (2+mdata%nhist)
    end do
    
    ! total number of state values:
    mdata%nvalue = mdata%nlocvalue * mdata%nloc

    ! info ...
    write (gol,'("  MAORI model variables:")'); call goPr
    write (gol,'("    -------------------------")'); call goPr
    write (gol,'("      nr name                ")'); call goPr
    write (gol,'("    ---- --------------------")'); call goPr
    do ivar = 1, mdata%nvar
      write (gol,'("    ",i4," ",a20)') ivar, trim(mdata%varlist%var(ivar)%name); call goPr
    end do
    write (gol,'("    -------------------------")'); call goPr
    
    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
      ! storage:
      allocate( mdata%obs_iloc2iobs(mdata%nloc), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( mdata%obs_value (mdata%nloc,mdata%obs_nlev,mdata%obs_nvar,mdata%nhist), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( mdata%obs_errstd(mdata%nloc,mdata%obs_nlev,mdata%obs_nvar,mdata%nhist), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! no file openend yet:
      mdata%obs_file = 'no-data-file-yet.txt'
      mdata%obs_opened = .false.
      mdata%obs_iloc2iobs = 0
      mdata%obs_t = AnyDate()
      ! extract from query:
      obs_type='none'
        call goVarValue( mdata%obs_query, ';', 'type', '=', obs_type, status )
        IF_ERROR_RETURN(status=1)
      ! storage:
      select case ( trim(obs_type) )
        case ( 'dummy' )
          ! nothing extra required ...
        case ( 'emep-daily-csv' )
          ! setup access to observation data:
          allocate( mdata%obs_emep, stat=status )
          IF_NOTOK_RETURN(status=1)
        case ( 'csv' )
          ! nothing to be init ...
        case ( 'obsground-nc' )
          ! nothing to be init ...
        case ( 'lml-nc' )
          ! nothing to be init ...
        case default
          write (gol,'("unsupported input type : ",a)') trim(obs_type); call goErr
          TRACEBACK; status=1; return
      end select

      ! output variable corresponding to observed variable:    
      allocate( mdata%obs_ivar(mdata%obs_nvar), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over observed variables:
      do ivar = 1, mdata%obs_nvar
        ! search in list of input variables:
        call goMatchValue( trim(mdata%obs_varlist%var(ivar)%name), &
                mdata%varlist%var(1:mdata%nvar)%name, mdata%obs_ivar(ivar), status )
        IF_NOTOK_RETURN(status=1)
      end do

      ! assimilation flags:
      allocate( mdata%obs_astat(mdata%nloc,mdata%obs_nlev,mdata%obs_nvar,mdata%nhist), stat=status )
      IF_NOTOK_RETURN(status=1)
      mdata%obs_astat = 0
      ! set some bits ...
      if ( .not. mdata%assim_analyse ) mdata%obs_astat = IBSet( mdata%obs_astat, MAORI_ASTAT_VALIDATION )

      ! info ...
      write (gol,'("  MAORI observed variables:")'); call goPr
      write (gol,'("    ----------------------------------------------")'); call goPr
      write (gol,'("      nr name                 model variable      ")'); call goPr
      write (gol,'("    ---- -------------------- --------------------")'); call goPr
      do ivar = 1, mdata%obs_nvar
        write (gol,'("    ",i4," ",a20," ",a20)') ivar, trim(mdata%obs_varlist%var(ivar)%name), trim(mdata%varlist%var(mdata%obs_ivar(ivar))%name); call goPr
      end do
      write (gol,'("    ----------------------------------------------")'); call goPr

    end if  ! obs_nvar > 0

    ! reset flag:
    mdata%started = .true.

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Start


  ! ***


  subroutine MAORI_Set_Sample_Data_Setup( mdata, t1, t2, the_end, status )

    use Binas        , only : Rgas, xm_O, xm_N, xm_S, p0
    use GO           , only : NewDate, IncrDate, IsAnyDate, Get, rTotal, wrtgol, Get_Begin_Of
    use GO           , only : operator(/=), operator(==), operator(>)
    use GO           , only : operator(+), operator(-), operator(*)
    use GO           , only : goReplace, goVarValue
    use GO           , only : ISTA_Add, ISTA_Next, ISTA_Shift

    use MAORI_LocList       , only : MAORI_LocList_Find
    use MAORI_DataFile_EMEP , only : Init, Done, FindRecord, Get, GetObservation
    use MAORI_DataFile_CSV  , only : T_MAORI_DataFile_CSV
    use MAORI_DataFile_AQORD, only : MAORI_DataFile_AQORD_ReadRecord
    use MAORI_DataFile_LML  , only : MAORI_DataFile_LML_ReadRecord

    use MAORI_Param  , only : MAORI_LEN_FILE, MAORI_LEN_NAME
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(inout)  ::  mdata
    type(TDate), intent(in)                   ::  t1, t2
    logical, intent(in)                       ::  the_end
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_Data_Setup'

    ! --- local -----------------------------------
    
    integer                         ::  date1(6)
    real                            ::  w_aver, w_prev, w_curr
    real                            ::  frac
    
    type(T_MAORI_DataFile_CSV)      ::  obs_CSV

    character(len=MAORI_LEN_NAME)   ::  obs_type
    character(len=MAORI_LEN_FILE)   ::  fname, fname_tmp
    integer                         ::  nobs, iobs
    character(len=MAORI_LEN_NAME)   ::  vname, obsvname
    character(len=MAORI_LEN_NAME)   ::  station_code
    integer                         ::  iloc, ivar
    integer                         ::  nloc_found
    integer                         ::  k
    type(TDate)                     ::  tcur
    character(len=MAORI_LEN_NAME)   ::  comp, unit
    real                            ::  val
    real                            ::  unitfac

    logical                         ::  exist
    logical                         ::  obs_allow_missing
    
    type(TDate)                     ::  read_tr(2)
    
    real                            ::  molevolume, molemass
    
    ! --- begin -----------------------------------
    
    !
    ! time step settings
    !
    
    ! setup at begin of time step over [t1,t2] :
    !  o a sample will be received at t2
    
    ! info ...
    write (gol,'("MAORI Setup:  setup for ",a)') trim(mdata%name); call goPr
    call wrtgol( 'MAORI Setup:    time step       : ', t1, ' - ', t2 ); call goPr
    
    ! check ...
    if ( t1 /= mdata%setup_tr(2) ) then
      write (gol,'("setup interval not connected with previous interval:")'); call goErr
      call wrtgol( '  previous setup interval : ', mdata%setup_tr(1), ' - ', mdata%setup_tr(2) ); call goErr
      call wrtgol( '  t1, t2                  : ', t1, ' - ', t2 ); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract date arrays:
    call Get( t1, time6=date1 )
    
    ! new period ?
    if ( IsAnyDate(mdata%period_tr(1)) ) then
      mdata%new_period = .true.
    else
      mdata%new_period =  t2 > mdata%period_tr(2)
    end if
    ! info ...
    call wrtgol( 'MAORI Setup:    (prev) period   : ', mdata%period_tr(1), ' - ', mdata%period_tr(2) ); call goPr
    write (gol,'("MAORI Setup:    new period      : ",l2)') mdata%new_period; call goPr
    
    ! set time period including t1:
    select case ( trim(mdata%tres_period) )
      case ( 'day' )
        mdata%period_tr(1) = NewDate( year=date1(1), month=date1(2), day=date1(3), hour=00 )
        mdata%period_tr(2) = mdata%period_tr(1) + IncrDate(hour=24)
      case default
        write (gol,'("could not determine time range period for : ",a)') trim(mdata%tres_period); call goErr
        TRACEBACK; status=1; return
    end select
    ! check ...
    if ( t2 > mdata%period_tr(2) ) then
      write (gol,'("setup period partly outside periode defined by t1:")'); call goErr
      call wrtgol( '  t1, t2     : ', t1, ' - ', t2 ); call goErr
      call wrtgol( '  period     : ', mdata%period_tr(1), ' - ', mdata%period_tr(2) ); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! end of period ?
    mdata%end_of_period = t2 == mdata%period_tr(2)
    
    ! info ...
    call wrtgol( 'MAORI Setup:    current period  : ', mdata%period_tr(1), ' - ', mdata%period_tr(2) ); call goPr
    write (gol,'("MAORI Setup:    end of period   : ",l2)') mdata%end_of_period; call goPr
    write (gol,'("MAORI Setup:    prev. fin. aver : ",l2)') mdata%finished_aver; call goPr

    ! previously finished an averaging interval ? then step to new inteval:
    if ( mdata%finished_aver ) then
      ! info ...
      write (gol,'("MAORI Setup:    next aver interval ...")'); call goPr
      ! init new interval; return weights to build new average:
      !    aver :=   aver * w_aver1 + prev * w_prev1 + curr * w_curr1
      call ISTA_Next( mdata%ista, w_prev, w_curr, status )
      IF_NOTOK_RETURN(status=1)
      ! store:
      mdata%w_aver = 0.0
      mdata%w_prev = w_prev
      mdata%w_curr = w_curr
      ! reset flag:
      mdata%finished_aver = .false.
      ! increase step counter:
      mdata%istep = mdata%istep + 1
      if ( mdata%istep > mdata%nstep ) mdata%istep = 1
      ! info ...
      write (gol,'("MAORI Setup:      weights       : ",3f6.2)') 0.0, w_prev, w_curr; call goPr
      write (gol,'("MAORI Setup:      istep         : ",i6)') mdata%istep; call goPr
    else
      ! reset weights:
      !    aver :=   aver * w_aver1 + prev * w_prev1 + curr * w_curr1
      mdata%w_aver = 1.0
      mdata%w_prev = 0.0
      mdata%w_curr = 0.0
    end if

    ! obtain weights and flags for samples added at t2:
    call ISTA_Add( mdata%ista, t2, w_aver, w_prev, w_curr, mdata%finished_aver, status )
    IF_NOTOK_RETURN(status=1)
    ! add weights:
    !    aver  ~   aver * w_aver1 + prev * w_prev1 + curr * w_curr1
    !    aver :=   aver * w_aver2 + prev * w_prev2 + curr * w_curr2
    !          =   aver * (w_aver1 * w_aver2          ) + 
    !              prev * (w_prev1 * w_aver2 + w_prev2) + 
    !              curr * (w_curr1 * w_aver2 + w_curr2)
    mdata%w_aver = mdata%w_aver * w_aver
    mdata%w_prev = mdata%w_prev * w_aver + w_prev
    mdata%w_curr = mdata%w_curr * w_aver + w_curr
    ! info ...
    write (gol,'("MAORI Setup:    add weights     : ",3f6.2)') w_aver, w_prev, w_curr; call goPr
    write (gol,'("MAORI Setup:    tot weights     : ",3f6.2)') mdata%w_aver, mdata%w_prev, mdata%w_curr; call goPr
    ! check ...
    if ( abs( mdata%w_aver + mdata%w_prev + mdata%w_curr - 1.0 ) > 0.01 ) then
      write (gol,'("weights do not sum to 1 : ")'); call goErr
      write (gol,'("  w_aver, w_prev, w_curr  : ",3f8.4)') mdata%w_aver, mdata%w_prev, mdata%w_curr; call goErr
      write (gol,'("  sum                     : ",3f8.4)') mdata%w_aver+mdata%w_prev+mdata%w_curr; call goErr
      TRACEBACK; status=1; return
    end if
    write (gol,'("MAORI Setup:    is finish aver  : ",l2)') mdata%finished_aver; call goPr

    ! determine step within period; trap initial setup for instant time:
    frac = rTotal( t2 - mdata%period_tr(1), mdata%tres_unit ) / mdata%tres_step
    mdata%istep = ceiling(frac)
    if ( (t1 == t2) .and. (frac == ceiling(frac)) ) mdata%istep = mdata%istep + 1
    ! check ...
    if ( (mdata%istep < 1) .or. (mdata%istep > mdata%nstep) ) then
      write (gol,'("step index ",i6," out of range 1 .. ",i6)') mdata%istep, mdata%nstep; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! history:
    if ( mdata%tres_hist ) then
      mdata%ihist = modulo( mdata%istep-1, mdata%nhist ) + 1
    else
      mdata%ihist = 1
    end if

    ! info ...
    write (gol,'("MAORI Setup:    current step    : ",i6," / ",i6)') mdata%istep, mdata%nstep; call goPr
    write (gol,'("MAORI Setup:    current hist    : ",i6," / ",i6)') mdata%ihist, mdata%nhist; call goPr

    ! shift current into previous:
    call ISTA_Shift( mdata%ista, status )
    IF_NOTOK_RETURN(status=1)
    
    ! store:
    mdata%setup_tr(1) = t1
    mdata%setup_tr(2) = t2
    mdata%the_end     = the_end
    
    !
    ! observations
    !
    
    ! read data ?
    if ( mdata%obs_nvar > 0 ) then

      ! extract from query:
      obs_type='none'
        call goVarValue( mdata%obs_query, ';', 'type', '=', obs_type, status )
        IF_ERROR_RETURN(status=1)
    
      ! split by type:
      select case ( trim(obs_type) )
      
        ! for testing ...
        case ( 'dummy' )
        
          ! loop over observed variables:
          do ivar = 1, mdata%obs_nvar
            ! set conversion factor:
            select case ( trim(mdata%obs_varlist%var(ivar)%cf_units) )
              case ( 'mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
              case ( 'kg m-3'      ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
              case default
                write (gol,'("do not know how to set dummy unitfac for  `",a,"`")') &
                        trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                TRACEBACK; status=1; return
            end select
            ! just some value ...
            mdata%obs_value (:,:,ivar,mdata%ihist) = 1.5 * unitfac
            ! set error fraction
            mdata%obs_errstd(:,:,ivar,mdata%ihist) = mdata%obs_value(:,:,ivar,mdata%ihist) * mdata%obs_err_frac
            ! apply minimum
            select case ( trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
              case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
              case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
              case default
                write (gol,'("do not know how to set unitfac for  `",a,"`")') &
                        trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                TRACEBACK; status=1; return
            end select
            mdata%obs_errstd(:,:,ivar,mdata%ihist) = max( mdata%obs_errstd(:,:,ivar,mdata%ihist), mdata%obs_err_min_value*unitfac )
          end do
              
        ! daily averages
        case ( 'emep-daily-csv' )

          ! extract filename from query, might contain %{.} keys:
          fname='none'
            call goVarValue( mdata%obs_query, ';', 'file', '=', fname, status )
            IF_ERROR_RETURN(status=1)
          ! set variable name used in filename and header:
          vname = trim(mdata%obs_varlist%var(1)%name)
          obsvname = trim(mdata%obs_vname(1))
          ! replace some %{.} keys in the filename:
          call goReplace( fname, '%{var}', trim(vname), status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( fname, '%{vname}', trim(obsvname), status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( fname, '%{year4}', '(i4.4)', t1%year, status )
          IF_NOTOK_RETURN(status=1)

          ! new file ?
          if ( trim(fname) /= trim(mdata%obs_file) ) then
            ! already opened ?
            if ( mdata%obs_opened ) then
              ! close current file:
              call Done( mdata%obs_emep, status )
              IF_NOTOK_RETURN(status=1)
              ! reset flag:
              mdata%obs_opened = .false.
            end if
            ! open new data file:
            call Init( mdata%obs_emep, trim(fname), status )
            IF_NOTOK_RETURN(status=1)
            ! reset flags:
            mdata%obs_opened = .true.
            mdata%obs_file   = trim(fname)
          end if

          ! search corresponding record in observations:
          call FindRecord( mdata%obs_emep, t1, t2, status )
          IF_NOTOK_RETURN(status=1)
          
          ! fill mapping from locations to observations ?
          if ( any(mdata%obs_iloc2iobs == 0) ) then
            ! number of observations:
            call Get( mdata%obs_emep, status, nobs=nobs )
            IF_NOTOK_RETURN(status=1)
            ! no iobs set yet:
            mdata%obs_iloc2iobs = -1
            ! loop over observations:
            do iobs = 1, nobs
              ! extract station code:
              call GetObservation( mdata%obs_emep, iobs, status, station_code=station_code )
              IF_NOTOK_RETURN(status=1)
              ! compare with required location code's:
              iloc = -1
              do k = 1, mdata%nloc
                if ( trim(station_code) == mdata%loclist%loc(k)%code ) then
                  iloc = k
                  exit
                end if
              end do
              if ( iloc > 0 ) then
                ! found; store observation index:
                mdata%obs_iloc2iobs(iloc) = iobs
              end if
            end do
            ! check ...
            if ( any(mdata%obs_iloc2iobs < 0) ) then
              write (gol,'("WARNING - some station codes not found in observation file:")'); call goPr
              do iloc = 1, mdata%nloc
                if ( mdata%obs_iloc2iobs(iloc) < 0 ) then
                  write (gol,'("WARNING -   ",i4," ",a)') iloc, trim(mdata%loclist%loc(iloc)%code); call goPr
                end if
              end do
              write (gol,'("WARNING - fill with no-data value")'); call goPr
              write (gol,'("WARNING - in ",a)') rname; call goPr
            end if
          end if  ! initialise iloc2iobs
          
          ! extract observation time:
          call Get( mdata%obs_emep, status, t=tcur )
          IF_ERROR_RETURN(status=1)
          ! different ?
          if ( IsAnyDate(mdata%obs_t) .or. (tcur /= mdata%obs_t) ) then
            ! only one variable ...
            if ( mdata%obs_nvar /= 1 ) then
              write (gol,'("emep csv files for single variable only yet ...")'); call goErr
              TRACEBACK; status=1; return
            end if
            ivar = 1
            ! loop over locations:
            do iloc = 1, mdata%nloc
              ! observation index:
              iobs = mdata%obs_iloc2iobs(iloc)
              ! no observations for this location ?
              if ( iobs < 0 ) then
                ! no data:
                val = -999.9
              else
                ! extract values:
                call GetObservation( mdata%obs_emep, iobs, status, comp=comp, unit=unit, value=val )
                IF_NOTOK_RETURN(status=1)
                ! check ...
                if ( (trim(comp) == trim(mdata%obs_varlist%var(ivar)%name)) .or. &
                     ((trim(comp) == 'PM10') .and. (trim(mdata%obs_varlist%var(ivar)%name)=='tpm10')) ) then
                  ! match ...
                else
                  ! no match ....
                  write (gol,'("component mismatch; required `",a,"`, found `",a,"`")') &
                           trim(mdata%obs_varlist%var(ivar)%name), trim(comp); call goErr
                  TRACEBACK; status=1; return
                end if
                ! unit conversion:
                select case ( trim(unit)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
                  case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
                  case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
                  case default
                    write (gol,'("do not know how to convert from `",a,"` to `",a,"`")') &
                            trim(unit), trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                    TRACEBACK; status=1; return
                end select
                ! convert:
                val = val * unitfac
              end if
              ! store in lowest level:
              mdata%obs_value(iloc,1,ivar,mdata%ihist) = val
              ! error standard deviation fraction:
              mdata%obs_errstd(iloc,1,ivar,mdata%ihist) = val * mdata%obs_err_frac
              ! apply minimum
              select case ( trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
                case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
                case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
                case default
                  write (gol,'("do not know how to set unitfac for  `",a,"`")') &
                          trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                  TRACEBACK; status=1; return
              end select
              mdata%obs_errstd(iloc,:,ivar,mdata%ihist) = max( mdata%obs_errstd(iloc,:,ivar,mdata%ihist), mdata%obs_err_min_value*unitfac )
              ! info ...
              !write (gol,'("  aaa3 --> read observation ",i4," ",f12.4)') iloc, val; call goPr
            end do  ! locations
            ! reset time:
            mdata%obs_t = tcur
          end if  ! new time
          
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              
        ! csv files, one per hour
        case ( 'csv' )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              

          ! extract filename from query, might contain %{.} keys:
          fname='none'
            call goVarValue( mdata%obs_query, ';', 'file', '=', fname, status )
            IF_ERROR_RETURN(status=1)
          obs_allow_missing = .false.
            call goVarValue( mdata%obs_query, ';', 'allow_missing', '=', obs_allow_missing, status )
            IF_ERROR_RETURN(status=1)

          ! set variable name used in filename and header:
          vname = trim(mdata%obs_varlist%var(1)%name)
          obsvname = trim(mdata%obs_vname(1))
          ! replace some %{.} keys in the filename:
          call goReplace( fname, '%{var}', trim(vname), status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( fname, '%{vname}', trim(obsvname), status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( fname, '%{year4}', '(i4.4)', t2%year, status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( fname, '%{month2}', '(i2.2)', t2%month, status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( fname, '%{day2}', '(i2.2)', t2%day, status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( fname, '%{hour2}', '(i2.2)', t2%hour, status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( fname, '%{min2}', '(i2.2)', t2%min, status )
          IF_NOTOK_RETURN(status=1)

          ! only one variable ...
          if ( mdata%obs_nvar /= 1 ) then
            write (gol,'("emep csv files for single variable only yet ...")'); call goErr
            TRACEBACK; status=1; return
          end if
          ivar = 1

          ! init with no-data, fill into lowest level:
          mdata%obs_value (:,1,ivar,mdata%ihist) = -999.9
          mdata%obs_errstd(:,1,ivar,mdata%ihist) = -999.9

          ! check ...
          inquire( file=trim(fname), exist=exist )
          if ( .not. exist ) then
            ! allowed ?
            if ( .not. obs_allow_missing ) then
              ! error ...
              write (gol,'("file not found : ",a)') trim(fname); call goErr
              TRACEBACK; status=1; return
            else
              ! warning:
              write (gol,'("WARNING - file not found : ",a)') trim(fname); call goPr
              ! no data value:
              val = -9999.9
            end if
          
          else

            ! open new data file:
            call obs_CSV%Init( trim(fname), mdata%obs_query, status )
            IF_NOTOK_RETURN(status=1)
            ! reset flags:
            mdata%obs_opened = .true.
            mdata%obs_file   = trim(fname)

            ! loop over data records:
            do
              ! read record, eventually filter records on end time:
              call obs_CSV%ReadRecord( station_code, obsvname, val, unit, status, t=t2 )
              if ( status < 0 ) then
                ! end of file:
                exit
              else if ( status > 0 ) then
                TRACEBACK; status=1; return
              end if

              ! find matching location:
              call MAORI_LocList_Find( mdata%loclist, iloc, status, &
                                          code=trim(station_code) )
              IF_NOTOK_RETURN(status=1)

              ! fix if file does not have component in columns:
              if ( trim(obsvname) == 'None' ) then
                ! copy:
                obsvname = trim(mdata%obs_vname(ivar))
              end if

              ! check ...
              if ( trim(obsvname) /= trim(mdata%obs_vname(ivar)) ) then
                write (gol,'("component name in data file `",a,"` does not match with `",a,"`")') &
                           trim(obsvname), trim(mdata%obs_vname(ivar)); call goErr
                TRACEBACK; status=1; return
              end if

              ! unit conversion:
              select case ( trim(unit)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
                case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
                case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
                case ( 'kg/m3 -> kg m-3'    ) ; unitfac = 1.0
                case ( 'kg m-3 -> kg m-3'    ) ; unitfac = 1.0
                !>>> use actual air density instead of reference ...
                !----> The problem is that in 'maori_varlist' the 'LE_CF_Standard'
                !      is called to obtain the prefered units ;
                !      should be changed to a 'MAORI_CF_Standard' to get prefered
                !      units 'ug/m3' for tracers!
                case ( 'kg/m3 -> mole mole-1' )
                  ! mole volume: V = RT/p  [m3/mole]
                  molevolume = Rgas * 293.0 / p0
                  ! set molemass:
                  select case ( trim(vname) )
                    case ( 'o3', 'o3_biascorr'  ) ; molemass = xm_O * 3
                    case ( 'no2', 'no2_obs'     ) ; molemass = xm_N + xm_O * 2
                    case ( 'so2'                ) ; molemass = xm_S + xm_O * 2
                    case default
                      write (gol,'("no molemass defined for ",a)') trim(vname); call goErr
                      TRACEBACK; status=1; return
                  end select
                  ! conc_in * molevolume /         molemass         
                  !  ug/m3  * (m3/mole)  / (kg tracer)/(mol tracer) 
                  unitfac   = molevolume / molemass
                !<<<
                case default
                  write (gol,'("do not know how to convert from `",a,"` to `",a,"`")') &
                          trim(unit), trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                  TRACEBACK; status=1; return
              end select
              ! convert:
              val = val * unitfac
              
              ! store in lowest level:
              mdata%obs_value(iloc,1,ivar,mdata%ihist) = val
              ! error standard deviation fraction:
              mdata%obs_errstd(iloc,1,ivar,mdata%ihist) = val * mdata%obs_err_frac
              ! apply minimum
              select case ( trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
                case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
                case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
                case default
                  write (gol,'("do not know how to set unitfac for  `",a,"`")') &
                          trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                  TRACEBACK; status=1; return
              end select
              mdata%obs_errstd(iloc,:,ivar,mdata%ihist) = max( mdata%obs_errstd(iloc,:,ivar,mdata%ihist), mdata%obs_err_min_value*unitfac )

            end do ! records

            ! close:
            call obs_CSV%Done( status )
            IF_NOTOK_RETURN(status=1)
            ! reset flag:
            mdata%obs_opened = .false.

            ! reset time:
            mdata%obs_t = t2

          end if ! file present

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              
        ! 'standard' nc files
        case ( 'obsground-nc' )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              
        
          !! read at new period:
          !if ( mdata%new_period ) then
          ! <-- was only useful when reading daily averages ..
          ! now read every interval again, might be the same data as before

            ! extract filename from query, might contain %{.} keys:
            fname='none'
              call goVarValue( mdata%obs_query, ';', 'file', '=', fname_tmp, status )
              IF_ERROR_RETURN(status=1)
            obs_allow_missing = .false.
              call goVarValue( mdata%obs_query, ';', 'allow_missing', '=', obs_allow_missing, status )
              IF_ERROR_RETURN(status=1)

            ! loop over observed variables:
            do ivar = 1, mdata%obs_nvar

              ! set variable name used in filename and header:
              vname = trim(mdata%obs_varlist%var(ivar)%name)
              ! name of concentration record:
              obsvname = trim(mdata%obs_vname(ivar))

              ! init counter:
              nloc_found = 0
              ! loop over stations:
              do iloc = 1, mdata%nloc

                ! station code:
                station_code = trim(mdata%loclist%loc(iloc)%code)
                
                ! filename template:
                fname = trim(fname_tmp)

                ! replace some %{.} keys in the filename:
                call goReplace( fname, '%{station_code}', trim(station_code), status )
                IF_NOTOK_RETURN(status=1)
                call goReplace( fname, '%{country_code}', station_code(1:2), status )
                IF_NOTOK_RETURN(status=1)
                call goReplace( fname, '%{year4}', '(i4.4)', t1%year, status )
                IF_NOTOK_RETURN(status=1)
                call goReplace( fname, '%{var}', trim(vname), status )
                IF_NOTOK_RETURN(status=1)
                call goReplace( fname, '%{vname}', trim(obsvname), status )
                IF_NOTOK_RETURN(status=1)

                ! check ...
                inquire( file=trim(fname), exist=exist )
                if ( .not. exist ) then
                  ! allowed ?
                  if ( .not. obs_allow_missing ) then
                    ! error ...
                    write (gol,'("file not found : ",a)') trim(fname); call goErr
                    TRACEBACK; status=1; return
                  else
                    ! warning:
                    write (gol,'("WARNING - file not found : ",a)') trim(fname); call goPr
                    ! no data value:
                    val = -9999.9
                  end if
                !~ no aver inteval yet ?
                else if ( mdata%ista%ival_t1 == mdata%ista%ival_t2 ) then
                  ! skip ...
                  val = -9999.9
                !~ data should be available 
                else
                  !! info ...
                  !write (gol,'("MAORI/Sample: file found : ",a)') trim(fname); call goPr
                  !! set input time range:
                  !if ( trim(mdata%tres_unit) == 'hour' ) then
                  !  read_tr(1) = mdata%setup_tr(1)
                  !  read_tr(2) = read_tr(1) + IncrDate(hour=1) * mdata%tres_step
                  !else if ( trim(mdata%tres_unit) == 'day' ) then
                  !  read_tr(1) = Get_Begin_Of( mdata%setup_tr(1), 'day' )
                  !  read_tr(2) = read_tr(1) + IncrDate(day=1) * mdata%tres_step
                  !else
                  !  write (gol,'("do not know how to set time range for reading observed samples")'); call goErr
                  !  write (gol,'("for tres with step ",f6.2," and unit `",a,"`")') mdata%tres_step, trim(mdata%tres_unit); call goErr
                  !  TRACEBACK; status=1; return
                  !end if
                  ! setup time range, same as current target interval:
                  read_tr(1) = mdata%ista%ival_t1
                  read_tr(2) = mdata%ista%ival_t2
                  ! read record ...
                  call MAORI_DataFile_AQORD_ReadRecord( trim(fname), read_tr, &
                                                  val, unit, status, &
                                                  concname=trim(obsvname) )
                  IF_NOTOK_RETURN(status=1)
                  ! different units ?
                  if ( trim(unit) /= trim(mdata%obs_varlist%var(ivar)%cf_units) ) then
                    ! unit conversion:
                    select case ( trim(unit)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
                      case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
                      case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
                      !>>> use actual air density instead of reference ...
                      !----> The problem is that in 'maori_varlist' the 'LE_CF_Standard'
                      !      is called to obtain the prefered units ;
                      !      should be changed to a 'MAORI_CF_Standard' to get prefered
                      !      units 'ug/m3' for tracers!
                      case ( 'ug/m3 -> mole mole-1' )
                        ! mole volume: V = RT/p  [m3/mole]
                        molevolume = Rgas * 293.0 / p0
                        ! set molemass:
                        select case ( trim(vname) )
                          case ( 'o3', 'o3_biascorr'  ) ; molemass = xm_O * 3
                          case ( 'no2', 'no2_obs'     ) ; molemass = xm_N + xm_O * 2
                          case ( 'so2'                ) ; molemass = xm_S + xm_O * 2
                          case default
                            write (gol,'("no molemass defined for ",a)') trim(vname); call goErr
                            TRACEBACK; status=1; return
                        end select
                        ! conc_in * 1e-9  * molevolume *         molemass         
                        !  ug/m3  * kg/ug * (m3/mole)  / (kg tracer)/(mol tracer) 
                        unitfac = 1e-9 * molevolume / molemass
                      !<<<
                      case default
                        write (gol,'("do not know how to convert from `",a,"` to `",a,"`")') &
                                trim(unit), trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                        TRACEBACK; status=1; return
                    end select
                    ! convert:
                    val = val * unitfac
                  end if  ! different untis
                  ! counter:
                  nloc_found = nloc_found + 1
                end if  ! obs file present or not
                ! store in lowest level:
                mdata%obs_value(iloc,1,ivar,mdata%ihist) = val
                ! guess error standard deviation ...
                mdata%obs_errstd(iloc,1,ivar,mdata%ihist) = val * mdata%obs_err_frac
                ! apply minimum
                select case ( trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
                  case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
                  case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
                  case default
                    write (gol,'("do not know how to set unitfac for  `",a,"`")') &
                            trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                    write (gol,'("   filename : ",a)') trim(fname); call goErr
                    write (gol,'("   variable : ",a)') trim(vname); call goErr
                    TRACEBACK; status=1; return
                end select
                mdata%obs_errstd(iloc,1,ivar,mdata%ihist) = max( mdata%obs_errstd(iloc,1,ivar,mdata%ihist), mdata%obs_err_min_value*unitfac )
                
                !! info ...
                !write (*,*) '     ---> ', trim(station_code), trim(aname), t1%year, mdata%obs_value(iloc,1,ivar,mdata%ihist), mdata%obs_errstd(iloc,1,ivar,mdata%ihist)

              end do  ! stations
              
              ! info ...
              write (gol,'("    --> found ",i4," of ",i4," locations")') nloc_found, mdata%nloc; call goPr
              
              ! safety check ...
              if ( nloc_found == 0 ) then
                write (gol,'("no valid locations found, probably something wrong with filename expansion ...")'); call goErr
                write (gol,'("add `allow_missing=F` to the query to list the files that could not be found")'); call goErr
                TRACEBACK; status=1; return
              end if

            end do  ! variables
            
          !end if  ! new period
          
        ! Dutch LML nc files
        case ( 'lml-nc' )

          !! read at new period:
          !if ( mdata%new_period ) then
          ! <-- was only useful when reading daily averages ..
          ! now read every interval again, might be the same data as before

          ! extract filename from query, might contain %{.} keys:
          fname='none'
          call goVarValue( mdata%obs_query, ';', 'file', '=', fname_tmp, status )
          IF_ERROR_RETURN(status=1)
          obs_allow_missing = .false.
          call goVarValue( mdata%obs_query, ';', 'allow_missing', '=', obs_allow_missing, status )
          IF_ERROR_RETURN(status=1)

          ! loop over observed variables:
          do ivar = 1, mdata%obs_nvar

            ! filename template:
            fname = trim(fname_tmp)

            ! replace some %{.} keys in the filename:
            call goReplace( fname, '%{year4}', '(i4.4)', t1%year, status )
            IF_NOTOK_RETURN(status=1)

            ! set variable name used in filename and header:
            vname = trim(mdata%obs_varlist%var(ivar)%name)

            ! name of concentration record:
            obsvname = trim(mdata%obs_vname(ivar))

            ! init counter:
            nloc_found = 0

            ! loop over stations:
            do iloc = 1, mdata%nloc

              ! station code:
              station_code = trim(mdata%loclist%loc(iloc)%code)

              ! check ...
              inquire( file=trim(fname), exist=exist )
              if ( .not. exist ) then
                ! allowed ?
                if ( .not. obs_allow_missing ) then
                   ! error ...
                   write (gol,'("file not found : ",a)') trim(fname); call goErr
                   TRACEBACK; status=1; return
                else
                   !! warning:
                   !write (gol,'("WARNING - file not found : ",a)') trim(fname); call goPr
                   ! no data value:
                   val = -9999.9
                end if
                !~ no aver inteval yet ?
              else if ( mdata%ista%ival_t1 == mdata%ista%ival_t2 ) then
                ! skip ...
                val = -9999.9
                !~ data should be available 
              else
                ! setup time range, same as current target interval:
                read_tr(1) = mdata%ista%ival_t1
                read_tr(2) = mdata%ista%ival_t2
                ! read record ...
                call MAORI_DataFile_LML_ReadRecord( trim(fname), read_tr, &
                     station_code, val, unit, status)
                IF_NOTOK_RETURN(status=1)
                ! different units ?
                if ( trim(unit) /= trim(mdata%obs_varlist%var(ivar)%cf_units) ) then
                  ! unit conversion:
                  select case ( trim(unit)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
                    case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
                    case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
                    !>>> use actual air density instead of reference ...
                    !case ( 'ug/m3 -> mole mole-1' )
                    !  ! mole volume: V = RT/p  [m3/mole]
                    !  molevolume = Rgas * 293.0 / p0
                    !  ! set molemass:
                    !  select case ( trim(vname) )
                    !  case ( 'o3'  ) ; molemass = xm_O * 3
                    !  case ( 'no2' ) ; molemass = xm_N + xm_O * 2
                    !  case default
                    !     write (gol,'("no molemass defined for ",a)') trim(vname); call goErr
                    !     TRACEBACK; status=1; return
                    !  end select
                    !  ! conc_in * 1e-9  * molevolume *         molemass         
                    !  !  ug/m3  * kg/ug * (m3/mole)  / (kg tracer)/(mol tracer) 
                    !  unitfac = 1e-9 * molevolume / molemass
                    !<<<
                    case default
                     write (gol,'("do not know how to convert from `",a,"` to `",a,"`")') &
                          trim(unit), trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                     TRACEBACK; status=1; return
                  end select
                  ! convert:
                  val = val * unitfac
                end if  ! different untis
                ! counter:
                nloc_found = nloc_found + 1
              end if  ! obs file present or not
              ! store in lowest level:
              mdata%obs_value(iloc,1,ivar,mdata%ihist) = val
              ! guess error standard deviation ...
              mdata%obs_errstd(iloc,1,ivar,mdata%ihist) = val * mdata%obs_err_frac
              ! apply minimum
              select case ( trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units) )
              case ( 'ppb -> mole mole-1' ) ; unitfac = 1.0e-9   ! ppb -> mole/mole
              case ( 'ug/m3 -> kg m-3'    ) ; unitfac = 1.0e-9   ! ug/m3 -> kg/m3
              case default
                 write (gol,'("do not know how to set unitfac for  `",a,"`")') &
                      trim(mdata%obs_err_min_units)//' -> '//trim(mdata%obs_varlist%var(ivar)%cf_units); call goErr
                 TRACEBACK; status=1; return
              end select
              mdata%obs_errstd(iloc,1,ivar,mdata%ihist) = max( mdata%obs_errstd(iloc,1,ivar,mdata%ihist), mdata%obs_err_min_value*unitfac )

              !! info ...
              !write (*,*) '     ---> ', trim(station_code), trim(aname), t1%year, mdata%obs_value(iloc,1,ivar,mdata%ihist), mdata%obs_errstd(iloc,1,ivar,mdata%ihist)

            end do  ! stations

            ! info ...
            write (gol,'("    --> found ",i4," of ",i4," locations")') nloc_found, mdata%nloc; call goPr

          end do  ! variables

          !end if  ! new period

        ! unknown ...
        case default
          write (gol,'("unsupported input type : ",a)') trim(obs_type); call goErr
          TRACEBACK; status=1; return
      end select
      
    end if  ! input

    !
    ! done
    !
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Setup


  ! ***


  subroutine MAORI_Set_Sample_Data_Done( mdata, status )

    use GO, only : goVarValue
    use GO, only : ISTA_Done
    
    use MAORI_LocList, only : MAORI_LocList_Done
    use MAORI_VarList, only : MAORI_VarList_Done
#ifdef with_netcdf
    use NetCDF       , only : NF90_Close
#endif

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(inout)       ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_Data_Done'
    
    ! --- local ----------------------------------
    
    character(len=MAORI_LEN_NAME)   ::  obs_type

    ! --- begin -----------------------------------

    ! done with sample averaging:
    call ISTA_Done( mdata%ista, status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with locations:
    call MAORI_LocList_Done( mdata%loclist, status )
    IF_NOTOK_RETURN(status=1)

    ! done with variables:
    call MAORI_VarList_Done( mdata%varlist, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear levels:
    deallocate( mdata%nlev, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
      ! extract from query:
      obs_type='none'
        call goVarValue( mdata%obs_query, ';', 'type', '=', obs_type, status )
        IF_ERROR_RETURN(status=1)
      ! file input ?
      select case ( trim(obs_type) )
        case ( 'dummy' )
          ! nothing extra required ...
        case ( 'emep-daily-csv' )
          ! setup access to observation data:
          deallocate( mdata%obs_emep, stat=status )
          IF_NOTOK_RETURN(status=1)
        case ( 'csv' )
          ! nothing extra required ...
        case ( 'obsground-nc' )
          ! nothing extra required ...
       case ( 'lml-nc' )
          ! nothing extra required ...
        case default
          write (gol,'("unsupported input type : ",a)') trim(obs_type); call goErr
          TRACEBACK; status=1; return
      end select
      ! clear:
      deallocate( mdata%obs_vname, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mdata%obs_iloc2iobs, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mdata%obs_value, stat=status  )
      IF_NOTOK_RETURN(status=1)
      deallocate( mdata%obs_errstd, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mdata%obs_astat, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mdata%obs_ivar, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Done


  ! ***


  subroutine MAORI_Set_Sample_Data_Inq( mdata, status, &
                              nparam, nloc, nvar, obs_nvar, &
                              assim_analyse, assim_analyse_now, assim_rho, &
                              nvalue, nlocvalue, varname )
                              
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(in)          ::  mdata
    integer, intent(out)                      ::  status

    integer, intent(out), optional            ::  nparam
    integer, intent(out), optional            ::  nloc
    integer, intent(out), optional            ::  nvar
    integer, intent(out), optional            ::  obs_nvar
    logical, intent(out), optional            ::  assim_analyse
    logical, intent(out), optional            ::  assim_analyse_now
    real, intent(out), optional               ::  assim_rho
    integer, intent(out), optional            ::  nvalue
    integer, intent(out), optional            ::  nlocvalue
    character(len=*), intent(out), optional   ::  varname

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_Data_Inq'

    ! --- begin -----------------------------------

    ! number of required parameters:
    if ( present(nparam) ) nparam = npar

    ! return number of locations:
    if ( present(nloc) ) nloc = mdata%nloc

    ! return number of variables to be put out:
    if ( present(nvar) ) nvar = mdata%nvar

    ! return number of observation variables:
    if ( present(obs_nvar) ) obs_nvar = mdata%obs_nvar

    ! analyse during assimilation ?
    if ( present(assim_analyse) ) assim_analyse = mdata%assim_analyse
    
    ! something to analyse now ?
    if ( present(assim_analyse_now) ) then
      ! true if this set is to be analysed, and averaging of last element of history is finished:
      assim_analyse_now = mdata%assim_analyse .and. (mdata%finished_aver .and. (mdata%ihist == mdata%nhist))
    end if

    ! assimilation correlation scale:
    if ( present(assim_rho) ) assim_rho = mdata%assim_rho

    ! return number of elements in the state:
    if ( present(nvalue) ) nvalue = mdata%nvalue

    ! return number of elements in the state for a single location:
    if ( present(nlocvalue) ) nlocvalue = mdata%nlocvalue
    
    ! variable name ?
    if ( present(varname) ) then
      ! only possible for single var:
      if ( mdata%nvar /= 1 ) then
        write (gol,'("could only return varname for single variable set")'); call goErr
        write (gol,'("  set name : ",a)') trim(mdata%name); call goErr
        write (gol,'("  nvar     : ",i6)') mdata%nvar; call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      varname = trim(mdata%varlist%var(1)%name)
    end if

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Inq


  ! ********************************************************************
  ! ***
  ! *** parameters
  ! ***
  ! ********************************************************************


  subroutine MAORI_Set_Sample_Data_Param_Inq( mdata, ipar, type, name, unit, status )

    use MAORI_Param, only : MAORI_INT, MAORI_REAL, MAORI_CHAR
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(inout)       ::  mdata
    integer, intent(in)                       ::  ipar
    integer, intent(out)                      ::  type
    character(len=*), intent(out)             ::  name
    character(len=*), intent(out)             ::  unit
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_Data_Param_Inq'

    ! --- begin -----------------------------------

    ! inquire parameter:
    select case ( ipar )
      case ( ipar_bound_west      ) ; type = MAORI_REAL ; name = 'bound_west'      ; unit = 'degree_east'
      case ( ipar_bound_east      ) ; type = MAORI_REAL ; name = 'bound_east'      ; unit = 'degree_east'
      case ( ipar_bound_south     ) ; type = MAORI_REAL ; name = 'bound_south'     ; unit = 'degree_north'
      case ( ipar_bound_north     ) ; type = MAORI_REAL ; name = 'bound_north'     ; unit = 'degree_north'
      case ( ipar_nlayer          ) ; type = MAORI_INT  ; name = 'nlayer'          ; unit = '1'
      case ( ipar_bias_correction ) ; type = MAORI_CHAR ; name = 'bias_correction' ; unit = 'free'
      case default
        write (gol,'("unsupported parameter index : ",i6)') ipar; call goErr
        write (gol,'("  set sample name : ",a)') trim(mdata%name); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Param_Inq


  ! ***
  

  subroutine MAORI_Set_Sample_Data_Param_Put_i( mdata, ipar, ival, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(inout)       ::  mdata
    integer, intent(in)                       ::  ipar
    integer, intent(in)                       ::  ival
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_Data_Param_Put_i'

    ! --- begin -----------------------------------

    ! inquire parameter:
    select case ( ipar )
      case ( ipar_nlayer ) ; mdata%nlayer = ival
      case default
        write (gol,'("unsupported parameter index : ",i6)') ipar; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set flag:
    mdata%param_set(ipar) = .true.

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Param_Put_i


  ! ***
  

  subroutine MAORI_Set_Sample_Data_Param_Put_r( mdata, ipar, rval, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(inout)       ::  mdata
    integer, intent(in)                       ::  ipar
    real, intent(in)                          ::  rval
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_Data_Param_Put_r'

    ! --- begin -----------------------------------

    ! inquire parameter:
    select case ( ipar )
      case ( ipar_bound_west  ) ; mdata%bounds(1) = rval
      case ( ipar_bound_east  ) ; mdata%bounds(2) = rval
      case ( ipar_bound_south ) ; mdata%bounds(3) = rval
      case ( ipar_bound_north ) ; mdata%bounds(4) = rval
      case default
        write (gol,'("unsupported parameter index : ",i6)') ipar; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set flag:
    mdata%param_set(ipar) = .true.

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Param_Put_r


  ! ***
  

  subroutine MAORI_Set_Sample_Data_Param_Put_s( mdata, ipar, sval, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(inout)        ::  mdata
    integer, intent(in)                             ::  ipar
    character(len=*), intent(in)                    ::  sval
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_Data_Param_Put_s'

    ! --- begin -----------------------------------

    ! inquire parameter:
    select case ( ipar )
      case ( ipar_bias_correction ) ; mdata%bias_correction = sval
      case default
        write (gol,'("unsupported parameter index : ",i6)') ipar; call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set flag:
    mdata%param_set(ipar) = .true.

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Param_Put_s


  ! ********************************************************************
  ! ***
  ! *** locations
  ! ***
  ! ********************************************************************


  subroutine MAORI_Set_Sample_Data_Loc_Inq( mdata, iloc, status, &
                                      name, code, lon, lat, alt, &
                                      varname )

    use MAORI_Param  , only : MAORI_SAMPLE
    use MAORI_Param  , only : MAORI_TYPE_NAME
    use MAORI_LocList, only : MAORI_Loc_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(in)       ::  mdata
    integer, intent(in)                         ::  iloc
    integer, intent(out)                        ::  status

    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  code
    real, intent(out), optional                 ::  lon
    real, intent(out), optional                 ::  lat
    real, intent(out), optional                 ::  alt
    character(len=*), intent(out), optional     ::  varname

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_Data_Loc_Inq'

    ! --- begin ----------------------------------

    ! inquire standaard location structure:
    call MAORI_Loc_Inq( mdata%loclist, iloc, status, &
                           name=name, code=code, lon=lon, lat=lat, alt=alt )
    IF_NOTOK_RETURN(status=1)
    
    ! variable name ?
    if ( present(varname) ) then
      ! only possible for single var:
      if ( mdata%nvar /= 1 ) then
        write (gol,'("could only return varname for single variable set")'); call goErr
        write (gol,'("  set name : ",a)') trim(mdata%name); call goErr
        write (gol,'("  nvar     : ",i6)') mdata%nvar; call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      varname = trim(mdata%varlist%var(1)%name)
    end if

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Loc_Inq


  ! ********************************************************************
  ! ***
  ! *** variables
  ! ***
  ! ********************************************************************


  subroutine MAORI_Set_Sample_Data_Var_Inq( mdata, ivar, status, &
                                          name, unit, nlev )

    use MAORI_VarList, only : MAORI_Var_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(in)            ::  mdata
    integer, intent(in)                         ::  ivar
    integer, intent(out)                        ::  status
    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  unit
    integer, intent(out), optional              ::  nlev

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_Data_Var_Inq'

    ! --- begin -----------------------------------

    ! inquire output set:
    call MAORI_Var_Inq( mdata%varlist, ivar, status, &
                             name=name, unit=unit )
    IF_NOTOK_RETURN(status=1)

    ! return number of levels to be put out:
    if ( present(nlev) ) nlev = mdata%nlev(ivar)

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Var_Inq


  ! ********************************************************************
  ! ***
  ! *** observed values
  ! ***
  ! ********************************************************************
  

  subroutine MAORI_Set_Sample_Data_Obs_Get( mdata, obs_ivar, status, iloc, y, r, alfa )

    use MAORI_LocList, only : MAORI_Loc_Inq
    use MAORI_VarList, only : MAORI_Var_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(in)           ::  mdata
    integer, intent(in)                             ::  obs_ivar
    integer, intent(out)                            ::  status
    
    integer, intent(in), optional                   ::  iloc
    real, intent(out), optional                     ::  y       ! measured value
    real, intent(out), optional                     ::  r       ! error std.dev.
    real, intent(out), optional                     ::  alfa    ! screening factor

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_Data_Obs_Get'

    ! --- local ----------------------------------
    
    character(len=MAORI_LEN_LINE)         ::  code, name
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( mdata%obs_nlev /= 1  ) then
      write (gol,'("getobs not implemented for obs_nlev ",i4)') mdata%obs_nlev; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! copy measured value:
    if ( present(y) ) then
      ! check ...
      if ( .not. present(iloc) ) then
        write (gol,'("get `y` but `iloc` not present.")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! fill:
      y = mdata%obs_value(iloc,1,obs_ivar,mdata%ihist)
      ! debug ...
      if ( y > 1e10 ) then
        call MAORI_Loc_Inq( mdata%loclist, iloc, status, code=code )
        IF_NOTOK_RETURN(status=1)
        call MAORI_Var_Inq( mdata%obs_varlist, obs_ivar, status, name=name )
        IF_NOTOK_RETURN(status=1)
        write (gol,*) 'BUG: extreme observation value: ', y; call goErr
        write (gol,'("BUG: location ",i6," code `",a,"`")') iloc, trim(code); call goErr
        write (gol,'("BUG: variable ",i6," name `",a,"`")') obs_ivar, trim(name); call goErr
        TRACEBACK; status=1; return
      end if        
    end if

    ! copy error std.dev.
    if ( present(r) ) then
      ! check ...
      if ( .not. present(iloc) ) then
        write (gol,'("get `r` but `iloc` not present.")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! fill:
      r = mdata%obs_errstd(iloc,1,obs_ivar,mdata%ihist)
    end if

    ! screening factor:
    if ( present(alfa) ) then
      ! fill:
      alfa = mdata%assim_alfa
    end if

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Obs_Get


  ! ***


  subroutine MAORI_Set_Sample_Data_Obs_Put( mdata, obs_ivar, status, iloc, astat_ibset )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_Data), intent(inout)        ::  mdata
    integer, intent(in)                             ::  obs_ivar
    integer, intent(out)                            ::  status
    
    integer, intent(in), optional                   ::  iloc
    integer, intent(in), optional                   ::  astat_ibset  ! MAORI_ASTAT_*

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_Data_Obs_Get'

    ! --- begin ----------------------------------
    
    ! set assimilation flag:
    if ( present(astat_ibset) ) then
      ! check ...
      if ( .not. present(iloc) ) then
        write (gol,'("put `astat_ibset` but `iloc` not present.")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! fill:
      mdata%obs_astat(iloc,:,obs_ivar,mdata%ihist) = IBSet( mdata%obs_astat(iloc,:,obs_ivar,mdata%ihist), astat_ibset )
    end if

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_Data_Obs_Put
  

  ! ********************************************************************
  ! ***
  ! *** simulated values
  ! ***
  ! ********************************************************************
  

  subroutine MAORI_Set_Sample_State_Init( mstate, mdata, name, status )

    use GO  , only : AnyDate

    use MAORI_Param, only : MAORI_LEN_LINE

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(out)       ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)         ::  mdata
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_State_Init'

    ! --- local -----------------------------------
    
    logical     ::  ldum

    ! --- begin -----------------------------------
    
    !write (gol,'("MAORI:     init sampling ...")'); call goPr
    
    ! dummy to avoid warnings about unused variables:
    ldum = mdata%started

    ! store name:
    mstate%name = trim(name)
    
    ! no values stored yet:
    nullify( mstate%sample_set )
    nullify( mstate%sample_curr )
    nullify( mstate%sample_prev )
    nullify( mstate%sample_aver )
    
    ! not started yet:
    mstate%started = .false.

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_Init


  ! ***


  subroutine MAORI_Set_Sample_State_Done( mstate, mdata, status )

    use MAORI_LocList, only : MAORI_LocList_Done
    use MAORI_VarList, only : MAORI_VarList_Done

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(inout)       ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)           ::  mdata
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_State_Done'

    ! --- local -----------------------------------
    
    logical     ::  ldum

    ! --- begin ----------------------------------
    
    ! dummy to avoid warnings about unused variables:
    ldum = mdata%started

    ! either none or all were allocated ...
    if ( .not. associated(mstate%sample_set) ) then
      ! clear values:
      deallocate( mstate%sample_set, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mstate%sample_curr, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mstate%sample_prev, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mstate%sample_aver, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_Done


  ! ***


  subroutine MAORI_Set_Sample_State_Start( mstate, mdata, status )

    use GO, only : AnyDate
    use MAORI_LocList, only : MAORI_LocList_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(inout)       ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)           ::  mdata
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_State_Start'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. mdata%started ) then
      write (gol,'("could not start state while sample not started ...")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! setup storage for simulated values:
    allocate( mstate%sample_set(mdata%nloc,mdata%maxlev,mdata%nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( mstate%sample_curr(mdata%nloc,mdata%maxlev,mdata%nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( mstate%sample_prev(mdata%nloc,mdata%maxlev,mdata%nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( mstate%sample_aver(mdata%nloc,mdata%maxlev,mdata%nvar,mdata%nhist), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! initial average:
    mstate%sample_curr = 0.0
    mstate%sample_prev = 0.0
    mstate%sample_aver = 0.0
    
    ! reset flag:
    mstate%started = .true.
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_Start


  ! ***


  subroutine MAORI_Set_Sample_State_Setup( mstate, mdata, status )
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(inout)       ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)           ::  mdata
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_State_Setup'

    ! --- local ----------------------------------
    
    logical   ::  ldum
    
    ! --- begin ----------------------------------

    ! dummy to avoid warnings about unused variables:
    ldum = mdata%started

    ! shift current values into previous:
    mstate%sample_prev = mstate%sample_curr

    ! reset current simulated values:
    mstate%sample_curr = 0.0
    mstate%sample_set  = .false.
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_Setup


  ! ***


  subroutine MAORI_Set_Sample_State_Put( mstate, mdata, ivar, status, iloc, values )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(inout)       ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)           ::  mdata
    integer, intent(in)                             ::  ivar
    integer, intent(out)                            ::  status

    integer, intent(in), optional                   ::  iloc
    real, intent(in), optional                      ::  values(:)
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_State_Put'

    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. all( (/present(iloc),present(values)/) ) ) then
      write (gol,'("optional arguments `iloc` and `values` should be provided ...")'); call goErr
      TRACEBACK; status=1; return 
    end if
    ! check ...
    if ( size(values) /= mdata%nlev(ivar) ) then
      write (gol,'("put values with ",i4," levels, while expected ",i4)') &
               size(values), mdata%nlev(ivar); call goErr
      TRACEBACK; status=1; return
    end if

    ! store value in appropriate levels:
    mstate%sample_curr(iloc,1:mdata%nlev(ivar),ivar) = values
    
    ! set flag at all levels:
    mstate%sample_set(iloc,:,ivar) = .true.

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_Put
  

  ! ***
  
  
  subroutine MAORI_Set_Sample_State_Post( mstate, mdata, status )

    use GO, only : wrtgol
    use MAORI_Param, only : MAORI_NAN

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(inout)       ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)           ::  mdata
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_State_Post'

    ! --- local ----------------------------------

    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. all(mstate%sample_set) ) then
      write (gol,'("not all simulated values seem to be set")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    !! info ...
    !call wrtgol( rname//': post process `'//trim(mdata%name)//'` samples for : ',mdata%setup_tr(2)); call goPr
    !write (gol,'(a,":   flags and weights : ",l2,3f8.4)') rname, mdata%new_period, mdata%w_aver,mdata%w_prev,mdata%w_curr ; call goPr
    !write (gol,*) rname//':   aver,prev,curr    : ', mdata%ihist, sum(mstate%sample_aver(:,:,:,mdata%ihist))/size(mstate%sample_aver(:,:,:,mdata%ihist)), sum(mstate%sample_prev)/size(mstate%sample_prev), sum(mstate%sample_curr)/size(mstate%sample_curr) ; call goPr
    
    ! clean history ?
    if ( mdata%new_period ) then
      !! info ...
      !write (gol,'(a,":   clean history ...")') rname; call goPr
      ! check ...
      if ( mdata%w_aver /= 0.0 ) then
        write (gol,'(a,": new period but non-zero w_aver : ",f12.6)') rname, mdata%w_aver; call goErr
        TRACEBACK; status=1; return
      end if
      ! clear historical data:
      mstate%sample_aver = MAORI_NAN
    end if
    
    ! check ...
    if ( (mdata%ihist < 1) .or. (mdata%ihist > mdata%nhist) ) then
      write (gol,'("found ihist ",i0," out of range 1,..,",i0)') mdata%ihist, mdata%nhist; call goErr
      TRACEBACK; status=1; return
    end if

    ! apply weights:
    mstate%sample_aver(:,:,:,mdata%ihist) = mstate%sample_aver(:,:,:,mdata%ihist) * mdata%w_aver + &
                                            mstate%sample_prev                    * mdata%w_prev + &
                                            mstate%sample_curr                    * mdata%w_curr
                         
    ! ok
    status = 0
    
  end subroutine MAORI_Set_Sample_State_Post
  

  ! ***


  subroutine MAORI_Set_Sample_State_Obs_Get( mstate, mdata, obs_ivar, status, iloc, value )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(in)          ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)           ::  mdata
    integer, intent(in)                             ::  obs_ivar
    integer, intent(out)                            ::  status

    integer, intent(in), optional                   ::  iloc
    real, intent(out), optional                     ::  value

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Sample_State_Obs_Get'

    ! --- begin -----------------------------------
    
    ! return simulated value:
    if ( present(value) ) then
      ! check ...
      if ( .not. present(iloc) ) then
        write (gol,'("get `value` but `iloc` not present.")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( mdata%obs_nlev /= 1  ) then
        write (gol,'("getobs not implemented for obs_nlev ",i4)') mdata%obs_nlev; call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      value = mstate%sample_aver(iloc,1,mdata%obs_ivar(obs_ivar),mdata%ihist)
    end if

    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_Obs_Get


  ! ***


  subroutine MAORI_Set_Sample_State_Values_Get( mstate, mdata, values, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(in)              ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)               ::  mdata
    real, intent(out)                                   ::  values(:)
    integer, intent(out)                                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_State_Values_Get'

    ! --- local ----------------------------------
    
    integer       ::  ivar
    integer       ::  i
    integer       ::  n
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(values) /= mdata%nvalue ) then
      write (gol,'("size of output array (",i6,") not equal to expected number (",i6,")")') &
                               size(values), mdata%nvalue; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! copy values:
      n = mdata%nloc * mdata%nlev(ivar)
      values(i+1:i+n) = reshape(mstate%sample_curr(1:mdata%nloc,1:mdata%nlev(ivar),ivar),(/n/)) ; i = i + n
      values(i+1:i+n) = reshape(mstate%sample_prev(1:mdata%nloc,1:mdata%nlev(ivar),ivar),(/n/)) ; i = i + n
      n = mdata%nloc * mdata%nlev(ivar) * mdata%nhist
      values(i+1:i+n) = reshape(mstate%sample_aver(1:mdata%nloc,1:mdata%nlev(ivar),ivar,1:mdata%nhist),(/n/)) ; i = i + n
    end do
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_Values_Get


  ! ***


  subroutine MAORI_Set_Sample_State_Values_Put( mstate, mdata, values, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(inout)           ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)               ::  mdata
    real, intent(in)                                    ::  values(:)
    integer, intent(out)                                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_State_Values_Put'

    ! --- local ----------------------------------
    
    integer       ::  ivar
    integer       ::  i
    integer       ::  n
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(values) /= mdata%nvalue ) then
      write (gol,'("size of input array (",i6,") not equal to expected number (",i6,")")') &
                               size(values), mdata%nvalue; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! copy values:
      n = mdata%nloc * mdata%nlev(ivar)
      mstate%sample_curr(1:mdata%nloc,1:mdata%nlev(ivar),ivar) = reshape(values(i+1:i+n),(/mdata%nloc,mdata%nlev(ivar)/)) ; i = i + n
      mstate%sample_prev(1:mdata%nloc,1:mdata%nlev(ivar),ivar) = reshape(values(i+1:i+n),(/mdata%nloc,mdata%nlev(ivar)/)) ; i = i + n
      n = mdata%nloc * mdata%nlev(ivar) * mdata%nhist
      mstate%sample_aver(1:mdata%nloc,1:mdata%nlev(ivar),ivar,1:mdata%nhist) = reshape(values(i+1:i+n),(/mdata%nloc,mdata%nlev(ivar),mdata%nhist/)) ; i = i + n
    end do
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_Values_Put


  ! ***


  subroutine MAORI_Set_Sample_State_LocValues_Get( mstate, mdata, iloc, values, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(in)              ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)               ::  mdata
    integer, intent(in)                                 ::  iloc
    real, intent(out)                                   ::  values(:)
    integer, intent(out)                                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_State_LocValues_Get'

    ! --- local ----------------------------------
    
    integer       ::  ivar
    integer       ::  i
    integer       ::  n
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(values) /= mdata%nlocvalue ) then
      write (gol,'("size of output array (",i6,") not equal to expected number (",i6,")")') &
                               size(values), mdata%nlocvalue; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! copy values:
      n = mdata%nlev(ivar)
      values(i+1:i+n) = mstate%sample_curr(iloc,1:mdata%nlev(ivar),ivar) ; i = i + n
      values(i+1:i+n) = mstate%sample_prev(iloc,1:mdata%nlev(ivar),ivar) ; i = i + n
      n = mdata%nlev(ivar) * mdata%nhist
      values(i+1:i+n) = reshape(mstate%sample_aver(iloc,1:mdata%nlev(ivar),ivar,1:mdata%nhist),(/n/)) ; i = i + n
    end do
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_LocValues_Get


  ! ***


  subroutine MAORI_Set_Sample_State_LocValues_Put( mstate, mdata, iloc, values, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Sample_State), intent(inout)           ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)               ::  mdata
    integer, intent(in)                                 ::  iloc
    real, intent(in)                                    ::  values(:)
    integer, intent(out)                                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Sample_State_LocValues_Put'

    ! --- local ----------------------------------
    
    integer       ::  ivar
    integer       ::  i
    integer       ::  n
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(values) /= mdata%nlocvalue ) then
      write (gol,'("size of input array (",i6,") not equal to expected number (",i6,")")') &
                               size(values), mdata%nlocvalue; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over variables:
    do ivar = 1, mdata%nvar
      ! copy values:
      n = mdata%nlev(ivar)
      mstate%sample_curr(iloc,1:mdata%nlev(ivar),ivar) = values(i+1:i+n) ; i = i + n
      mstate%sample_prev(iloc,1:mdata%nlev(ivar),ivar) = values(i+1:i+n) ; i = i + n
      n = mdata%nlev(ivar) * mdata%nhist
      mstate%sample_aver(iloc,1:mdata%nlev(ivar),ivar,1:mdata%nhist) = reshape(values(i+1:i+n),(/mdata%nlev(ivar),mdata%nhist/)) ; i = i + n
    end do
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Sample_State_LocValues_Put


  ! ********************************************************************
  ! ***
  ! *** output
  ! ***
  ! ********************************************************************
  

  subroutine MAORI_Set_Output_State_Init( moutput, mdata, name, status )

    use GO, only : AnyDate

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output_State), intent(out)          ::  moutput
    type(T_MAORI_Set_Sample_Data), intent(in)             ::  mdata
    character(len=*), intent(in)                      ::  name
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Output_State_Init'

    ! --- local -----------------------------------
    
    logical     ::  ldum
    
    ! --- begin -----------------------------------
    
    ! dummy to avoid warnings about unused variables:
    ldum = mdata%started
    
    ! store:
    moutput%name = trim(name)
    
    ! no file opened yet:
    moutput%opened = .false.

    ! ok
    status = 0

  end subroutine MAORI_Set_Output_State_Init
  

  ! ***
  
  
  subroutine MAORI_Set_Output_State_Done( moutput, mdata, status )

#ifdef with_netcdf
    use NetCDF, only : NF90_Close
#endif

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output_State), intent(inout)        ::  moutput
    type(T_MAORI_Set_Sample_Data), intent(in)             ::  mdata
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Output_State_Done'

    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------

    ! file still open ?
    if ( moutput%opened ) then
#ifdef with_netcdf
      status = NF90_Close( moutput%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#endif
      ! reset flag:
      moutput%opened = .false.
    end if
    
    ! clear storage for netcdf id's :
    deallocate( moutput%varid, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
      deallocate( moutput%varid_y    , stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( moutput%varid_r    , stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( moutput%varid_astat, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0

  end subroutine MAORI_Set_Output_State_Done
  

  ! ***
  
  
  subroutine MAORI_Set_Output_State_Start( moutput, mdata, status )

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output_State), intent(inout)        ::  moutput
    type(T_MAORI_Set_Sample_Data), intent(in)             ::  mdata
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Output_State_Start'

    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! setup storage for netcdf id's :
    allocate( moutput%varid(mdata%nvar) )

    ! read observations ?
    if ( mdata%obs_nvar > 0 ) then
      ! setup storage for netcdf id's :
      allocate( moutput%varid_y    (mdata%obs_nvar), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( moutput%varid_r    (mdata%obs_nvar), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( moutput%varid_astat(mdata%obs_nvar), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! no file opened yet:
    moutput%opened = .false.

    ! ok
    status = 0

  end subroutine MAORI_Set_Output_State_Start
  

  ! ***
  
  
  subroutine MAORI_Set_Output_State_Create( moutput, mdata, t, status )

    use GO, only : TDate, NewDate, IncrDate, AnyDate, Get, Pretty
    use GO, only : operator(>), operator(+), operator(-)
    
    use MAORI_Param , only : MAORI_ASTAT_DESCRIPTION
    use MAORI_Common, only : com_path, com_model_name, com_exp_id, com_t0
    use MAORI_Common, only : MAORI_Common_NF90_GlobalAttributes
    
#ifdef with_netcdf
    use NetCDF, only : NF90_Create, NF90_Close
    use NetCDF, only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF, only : NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF, only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE
#endif

    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output_State), intent(inout)        ::  moutput
    type(T_MAORI_Set_Sample_Data), intent(in)             ::  mdata
    type(TDate), intent(in)                           ::  t
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Output_State_Create'

    ! --- local -----------------------------------
    
    integer                         ::  time6(6)
    character(len=MAORI_LEN_FILE)   ::  fname
    integer                         ::  cmode
    integer                         ::  ivar
    integer                         ::  dimid_z
    integer                         ::  varid
    character(len=MAORI_LEN_LINE)   ::  units
    
    ! --- begin -----------------------------------
    
    ! extract time fields from time assigned to current record:
    call Get( t, time6=time6 )

    ! new file name:
    write (fname,'(a,"_",a,"_",a,"_",i4.4,2i2.2)') &
              trim(com_model_name), trim(com_exp_id), &
              trim(mdata%name), time6(1:3)
    if ( len_trim(moutput%name) > 0 ) fname = trim(fname)//'_'//trim(moutput%name)
    fname = trim(fname)//'.nc'

    ! info ...
    write (gol,'("MAORI:     file name : ",a)') trim(fname); call goPr
    
    ! add path:
    fname = trim(com_path)//'/'//trim(fname)

#ifdef with_netcdf
    ! set creation mode flag:
    cmode = NF90_CLOBBER       ! overwrite existing files
    !cmode = NF90_NOCLOBBER     ! do not overwrite existing files

    ! create file:
    status = NF90_Create( fname, cmode, moutput%ncid )
    if ( status /= 0 ) then
      write (gol,'("creating file :")'); call goErr
      write (gol,'("  ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
#endif

    ! reset flag:
    moutput%opened = .true.

    ! write global attributes:
    call MAORI_Common_NF90_GlobalAttributes( moutput%ncid, status )
    IF_NOTOK_RETURN(status=1)

    ! extra global attributes:
    status = NF90_Put_Att( moutput%ncid, NF90_GLOBAL, 'bias_correction', trim(mdata%bias_correction) )
    IF_NF90_NOTOK_RETURN(status=1)

#ifdef with_netcdf
    ! define dimensions:

    status = NF90_Def_Dim( moutput%ncid, 'namelen', MAORI_LEN_NAME, moutput%dimid_namelen )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( moutput%ncid, 'codelen', MAORI_LEN_NAME, moutput%dimid_codelen )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( moutput%ncid, 'station', mdata%nloc, moutput%dimid_loc )
    IF_NF90_NOTOK_RETURN(status=1)
    ! output levels:
    status = NF90_Def_Dim( moutput%ncid, 'surf_level', 1, moutput%dimid_slev )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( moutput%ncid, 'full_level', mdata%nlayer, moutput%dimid_flev )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( moutput%ncid, 'half_level', mdata%nlayer+1, moutput%dimid_hlev )
    IF_NF90_NOTOK_RETURN(status=1)
    ! timing:
    status = NF90_Def_Dim( moutput%ncid, 'time', NF90_UNLIMITED, moutput%dimid_time )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( moutput%ncid, 'datelen', 6, moutput%dimid_datelen )
    IF_NF90_NOTOK_RETURN(status=1)

    ! define time variables:

    status = NF90_Def_Var( moutput%ncid, 'time', NF90_REAL, (/moutput%dimid_time/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'time' )
    IF_NF90_NOTOK_RETURN(status=1)
    write (units,'("days since ",a)') trim(Pretty(com_t0))
    status = NF90_Put_Att( moutput%ncid, varid, 'units', trim(units) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'calendar', 'standard' )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_time = varid

    status = NF90_Def_Var( moutput%ncid, 'date', NF90_INT, (/moutput%dimid_datelen,moutput%dimid_time/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'date and time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'units', 'year, month, day' )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_date = varid

    ! define time variables (end time) :

    status = NF90_Def_Var( moutput%ncid, 'time2', NF90_REAL, (/moutput%dimid_time/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'time' )
    IF_NF90_NOTOK_RETURN(status=1)
    write (units,'("days since ",a)') trim(Pretty(com_t0))
    status = NF90_Put_Att( moutput%ncid, varid, 'units', trim(units) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'calendar', 'standard' )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_time2 = varid

    status = NF90_Def_Var( moutput%ncid, 'date2', NF90_INT, (/moutput%dimid_datelen,moutput%dimid_time/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'date and time' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'units', 'year, month, day' )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_date2 = varid

    ! define station description data:

    status = NF90_Def_Var( moutput%ncid, 'station_name', NF90_CHAR, (/moutput%dimid_namelen,moutput%dimid_loc/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_stn_name = varid

    status = NF90_Def_Var( moutput%ncid, 'station_code', NF90_CHAR, (/moutput%dimid_codelen,moutput%dimid_loc/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_stn_code = varid

    status = NF90_Def_Var( moutput%ncid, 'station_lon', NF90_REAL, (/moutput%dimid_loc/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'units', 'degrees_east' )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_stn_lon = varid

    status = NF90_Def_Var( moutput%ncid, 'station_lat', NF90_REAL, (/moutput%dimid_loc/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'units', 'degrees_north' )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_stn_lat = varid

    status = NF90_Def_Var( moutput%ncid, 'station_alt', NF90_REAL, (/moutput%dimid_loc/), varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'standard_name', 'altitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'altitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( moutput%ncid, varid, 'units', 'm' )
    IF_NF90_NOTOK_RETURN(status=1)
    moutput%varid_stn_alt = varid

    ! observations

    ! loop over variables:
    do ivar = 1, mdata%nvar

      ! select level id:
      if ( mdata%nlev(ivar) == 1 ) then
        dimid_z = moutput%dimid_slev
      else if ( mdata%nlev(ivar) == mdata%nlayer ) then
        dimid_z = moutput%dimid_flev
      else if ( mdata%nlev(ivar) == mdata%nlayer+1 ) then
        dimid_z = moutput%dimid_hlev
      else
        write (gol,'("unsupported number of levels : ")'); call goErr
        write (gol,'("  variable    : ",i6," (",a,")")') ivar, trim(mdata%varlist%var(ivar)%name); call goErr
        write (gol,'("  nlev        : ",i6)') mdata%nlev(ivar); call goErr
        write (gol,'("  nlayer      : ",i6)') mdata%nlayer; call goErr
        TRACEBACK; status=1; return
      end if

      ! define veriable:
      status = NF90_Def_Var( moutput%ncid, trim(mdata%varlist%var(ivar)%name), &
                   NF90_REAL, (/moutput%dimid_loc,dimid_z,moutput%dimid_time/), &
                   varid )
      IF_NF90_NOTOK_RETURN(status=1)

      ! add attributes:
      status = NF90_Put_Att( moutput%ncid, varid, 'description', 'model simulation' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%varlist%var(ivar)%cf_standard_name) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%varlist%var(ivar)%cf_long_name) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%varlist%var(ivar)%cf_units) )
      IF_NF90_NOTOK_RETURN(status=1)
      ! store:
      moutput%varid(ivar) = varid

    end do

    ! loop over observerd variables:
    do ivar = 1, mdata%obs_nvar

      ! select level id:
      if ( mdata%obs_nlev == 1 ) then
        dimid_z = moutput%dimid_slev
      else
        write (gol,'("unsupported number of obs levels : ",i6)') mdata%obs_nlev; call goErr
        TRACEBACK; status=1; return
      end if

      ! define veriable:
      status = NF90_Def_Var( moutput%ncid, trim(mdata%obs_varlist%var(ivar)%name)//'_y', &
                      NF90_REAL, (/moutput%dimid_loc,dimid_z,moutput%dimid_time/), &
                      varid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! add attributes:
      status = NF90_Put_Att( moutput%ncid, varid, 'description', 'observation' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%obs_varlist%var(ivar)%cf_standard_name) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%obs_varlist%var(ivar)%cf_long_name) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%obs_varlist%var(ivar)%cf_units) )
      IF_NF90_NOTOK_RETURN(status=1)
      ! store:
      moutput%varid_y(ivar) = varid

      ! define veriable:
      status = NF90_Def_Var( moutput%ncid, trim(mdata%obs_varlist%var(ivar)%name)//'_r', &
                      NF90_REAL, (/moutput%dimid_loc,dimid_z,moutput%dimid_time/), &
                      varid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! add attributes:
      status = NF90_Put_Att( moutput%ncid, varid, 'description', 'observation error std.dev.' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'standard_name', trim(mdata%obs_varlist%var(ivar)%cf_standard_name) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'long_name', trim(mdata%obs_varlist%var(ivar)%cf_long_name) )
      IF_NF90_NOTOK_RETURN(status=1)
      status = nf90_put_att( moutput%ncid, varid, 'units', trim(mdata%obs_varlist%var(ivar)%cf_units) )
      IF_NF90_NOTOK_RETURN(status=1)
      ! store:
      moutput%varid_r(ivar) = varid

      ! define veriable:
      status = NF90_Def_Var( moutput%ncid, trim(mdata%obs_varlist%var(ivar)%name)//'_astat', &
                      NF90_REAL, (/moutput%dimid_loc,dimid_z,moutput%dimid_time/), &
                      varid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! add attributes:
      status = NF90_Put_Att( moutput%ncid, varid, 'long_name', 'assimilation status' )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Put_Att( moutput%ncid, varid, 'units', trim(MAORI_ASTAT_DESCRIPTION) )
      IF_NF90_NOTOK_RETURN(status=1)
      ! store:
      moutput%varid_astat(ivar) = varid
      
    end do

    ! end defintion mode:

    status = NF90_EndDef( moutput%ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif

    ! no records written yet:
    moutput%itrec = 0

    ! ok
    status = 0

  end subroutine MAORI_Set_Output_State_Create
  

  ! ***
  
  
  subroutine MAORI_Set_Output_State_Write( moutput, mstate, mdata, status )

    use GO, only : TDate, NewDate, AnyDate, Get
    use GO, only : operator(>), operator(-), operator(+), operator(*)
    use GO, only : rTotal
use GO, only : wrtgol
    
    use MAORI_Common, only : com_t0

#ifdef with_netcdf
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Close
#endif
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output_State), intent(inout)        ::  moutput
    type(T_MAORI_Set_Sample_State), intent(in)            ::  mstate
    type(T_MAORI_Set_Sample_Data), intent(in)             ::  mdata
    integer, intent(out)                              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Output_State_Write'

    ! --- local -----------------------------------
    
    integer                         ::  i
    character(len=MAORI_LEN_NAME)   ::  sval
    integer                         ::  ivar
    type(TDate)                     ::  t
    integer                         ::  time6(6)
    real                            ::  time
    integer                         ::  ihist
    
    ! --- begin -----------------------------------
    
    !! debug ...
    !write (gol,'(a,": finished_aver2 = ",l1)') rname, mdata%finished_aver2; call goPr
    ! info ...
    write (gol,'("MAORI: write set ",a," `",a,"` ...")') trim(mdata%name), trim(moutput%name); call goPr
      
    ! write if either:
    !  o averaging of last element of history is finished;
    !  o end of the run is reached:
    if ( (mdata%finished_aver .and. (mdata%ihist == mdata%nhist)) .or. mdata%the_end ) then
      ! info ...
      write (gol,'("MAORI:   write now ...")'); call goPr
    else
      ! info ...
      write (gol,'("MAORI:   do not write yet: finished aver ",l2,"; ihist ",i4,"/",i4)') mdata%finished_aver, mdata%ihist, mdata%nhist; call goPr
      ! ok
      status=0; return
    end if

    ! new file ?
    if ( .not. moutput%opened ) then
      ! info ...
      write (gol,'("MAORI:   create new file ...")'); call goPr
      ! new file:
      call MAORI_Set_Output_State_Create( moutput, mdata, mdata%period_tr(1), status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! info ...
    write (gol,'("MAORI:   write records; ihist is ",i4," of ",i4)') mdata%ihist, mdata%nhist; call goPr

#ifndef with_netcdf
    ! check libraries ...
    write (gol,'("not compiled with netcdf support")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! loop over historical records:
    do ihist = 1, mdata%nhist

      ! next time record:
      moutput%itrec = moutput%itrec + 1

      ! write station info only once ...    
      if ( moutput%itrec == 1 ) then
#ifdef with_netcdf
        ! write station info; 
        ! for some reason, character arrays need to be written one-by-one:
        do i = 1, mdata%nloc
          sval = trim(mdata%loclist%loc(i)%name)
          status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_name, trim(sval), &
                                   start=(/1,i/), count=(/len_trim(sval),1/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end do
        do i = 1, mdata%nloc
          sval = trim(mdata%loclist%loc(i)%code)
          status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_code, trim(sval), &
                                   start=(/1,i/), count=(/len_trim(sval),1/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end do
        ! write non-character fields:
        status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_lon, mdata%loclist%loc(1:mdata%nloc)%lon )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_lat, mdata%loclist%loc(1:mdata%nloc)%lat )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Put_Var( moutput%ncid, moutput%varid_stn_alt, mdata%loclist%loc(1:mdata%nloc)%alt )
        IF_NF90_NOTOK_RETURN(status=1)
#endif
      end if  ! first record

#ifdef with_netcdf
      ! averaging interval is split into nhist history steps;
      ! set start of interval for current history:
      t = mdata%ista%aver_t1 - mdata%ista%dt * (mdata%nhist-ihist)
      ! extract 6-value array:
      call Get( t, time6=time6 )
      ! write time record:
      time = rTotal( t - com_t0, 'day' )
      status = NF90_Put_Var( moutput%ncid, moutput%varid_time, (/time/), start=(/moutput%itrec/) )
      IF_NF90_NOTOK_RETURN(status=1)
      ! write date record:
      status = NF90_Put_Var( moutput%ncid, moutput%varid_date, reshape(time6,(/6,1/)), &
                                             start=(/1,moutput%itrec/), count=(/6,1/) )
      IF_NF90_NOTOK_RETURN(status=1)

      ! averaging interval is split into nhist history steps;
      ! set end of interval for current history:
      t = mdata%ista%aver_t2 - mdata%ista%dt * (mdata%nhist-ihist)
      ! averaging interval is split into nhist history steps;
      ! add offset for history:
      t = t + (ihist-1)*mdata%ista%dt
      ! extract 6-value array:
      call Get( t, time6=time6 )
      ! write time record:
      time = rTotal( t - com_t0, 'day' )
      status = NF90_Put_Var( moutput%ncid, moutput%varid_time2, (/time/), start=(/moutput%itrec/) )
      IF_NF90_NOTOK_RETURN(status=1)
      ! write date record:
      status = NF90_Put_Var( moutput%ncid, moutput%varid_date2, reshape(time6,(/6,1/)), &
                                             start=(/1,moutput%itrec/), count=(/6,1/) )
      IF_NF90_NOTOK_RETURN(status=1)


      ! loop over variables:
      do ivar = 1, mdata%nvar
        ! write simulations:
        status = NF90_Put_Var( moutput%ncid, moutput%varid(ivar), &
                       mstate%sample_aver(:,1:mdata%nlev(ivar),ivar,ihist), &
                       start=(/1,1,moutput%itrec/), &
                       count=(/mdata%nloc,mdata%nlev(ivar),1/) )
        IF_NF90_NOTOK_RETURN(status=1)
      end do

      ! write observations ?
      if ( mdata%obs_nvar > 0 ) then
        ! loop over observerd variables:
        do ivar = 1, mdata%obs_nvar
          ! write observations:
          status = NF90_Put_Var( moutput%ncid, moutput%varid_y(ivar), &
                         mdata%obs_value(:,1:mdata%obs_nlev,ivar,ihist), &
                         start=(/1,1,moutput%itrec/), &
                         count=(/mdata%nloc,mdata%obs_nlev,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
          ! write error standard deviations:
          status = NF90_Put_Var( moutput%ncid, moutput%varid_r(ivar), &
                         mdata%obs_errstd(:,1:mdata%obs_nlev,ivar,ihist), &
                         start=(/1,1,moutput%itrec/), &
                         count=(/mdata%nloc,mdata%obs_nlev,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
          ! write assimilation status:
          status = NF90_Put_Var( moutput%ncid, moutput%varid_astat(ivar), &
                         mdata%obs_astat(:,1:mdata%obs_nlev,ivar,ihist), &
                         start=(/1,1,moutput%itrec/), &
                         count=(/mdata%nloc,mdata%obs_nlev,1/) )
          IF_NF90_NOTOK_RETURN(status=1)
        end do
      end if
#endif

    end do  ! hist

    ! close if end of period or end of run is reached:
    if ( mdata%end_of_period .or. mdata%the_end ) then
      ! info ...
      write (gol,'("MAORI:   close file ...")'); call goPr
      ! close:
#ifdef with_netcdf
      status = NF90_Close( moutput%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#endif
      ! reset flag:
      moutput%opened = .false.
    end if

    ! ok
    status = 0

  end subroutine MAORI_Set_Output_State_Write



end module MAORI_Set_Sample
