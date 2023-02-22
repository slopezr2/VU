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
#include "mdf.inc"
!
!#######################################################################

module MAORI_Set

  use GO, only : gol, goPr, goErr

  use MAORI_Param    , only : MAORI_LEN_NAME
  use MAORI_Set_Sample   , only : T_MAORI_Set_Sample_Data   , T_MAORI_Set_Sample_State   , T_MAORI_Set_Output_State
  use MAORI_Set_Satellite, only : T_MAORI_Data_Satellite, T_MAORI_State_Satellite, T_MAORI_Output_Satellite
  
  implicit none


  ! --- in/out ---------------------------

  private

  public  ::  T_MAORI_Set_Data
  public  ::  MAORI_Set_Data_Init, MAORI_Set_Data_Done
  public  ::  MAORI_Set_Data_Start, MAORI_Set_Data_Setup
  public  ::  MAORI_Set_Data_Check
  public  ::  MAORI_Set_Data_Inq
  public  ::  MAORI_Set_Data_Param_Inq, MAORI_Set_Data_Param_Put
  public  ::  MAORI_Set_Data_Loc_Inq
  public  ::  MAORI_Set_Data_Var_Inq
  public  ::  MAORI_Set_Data_Obs_Get, MAORI_Set_Data_Obs_Put
  
  public  ::  T_MAORI_Set_State
  public  ::  MAORI_Set_State_Init, MAORI_Set_State_Done
  public  ::  MAORI_Set_State_Start, MAORI_Set_State_Setup
  public  ::  MAORI_Set_State_Put, MAORI_Set_State_Post
  public  ::  MAORI_Set_State_Obs_Get
  public  ::  MAORI_Set_State_Values_Get, MAORI_Set_State_Values_Put
  public  ::  MAORI_Set_State_LocValues_Get, MAORI_Set_State_LocValues_Put
  
  public  ::  T_MAORI_Set_Output
  public  ::  MAORI_Output_Set_Init, MAORI_Output_Set_Done
  public  ::  MAORI_Output_Set_Start, MAORI_Output_Set_Write


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI_Set'


  ! --- types ------------------------------

  type T_MAORI_Set_Data
    ! name:
    character(len=MAORI_LEN_NAME)             ::  name
    ! type (model field, station)
    integer                                   ::  type
    ! output sets:
    type(T_MAORI_Set_Sample_Data   ), pointer     ::  sample
    type(T_MAORI_Data_Satellite), pointer     ::  satellite
  end type T_MAORI_Set_Data

  type T_MAORI_Set_State
    ! name:
    character(len=MAORI_LEN_NAME)             ::  name
    ! type (model field, station)
    integer                                   ::  type
    ! output sets:
    type(T_MAORI_Set_Sample_State   ), pointer    ::  sample
    type(T_MAORI_State_Satellite), pointer    ::  satellite
  end type T_MAORI_Set_State

  type T_MAORI_Set_Output
    ! name:
    character(len=MAORI_LEN_NAME)             ::  name
    ! type (model field, station)
    integer                                   ::  type
    ! output sets:
    type(T_MAORI_Set_Output_State   ), pointer   ::  sample
    type(T_MAORI_Output_Satellite), pointer   ::  satellite
  end type T_MAORI_Set_Output


  ! --- interfaces -------------------------------
  
  interface MAORI_Set_Data_Param_Put
    module procedure MAORI_Set_Data_Param_Put_i
    module procedure MAORI_Set_Data_Param_Put_r
    module procedure MAORI_Set_Data_Param_Put_s
  end interface


contains


  ! ====================================================================
  ! ===
  ! === sets
  ! ===
  ! ====================================================================


  subroutine MAORI_Set_Data_Init( mdata, rcF, name, t, status )

    use GO, only : TRcFile, Init, Done, ReadRc
    use GO, only : goLoCase
    use GO, only : TDate
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Init
    use MAORI_Set_Satellite, only : MAORI_Data_Satellite_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(out)   ::  mdata
    type(TRcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  name
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Data_Init'

    ! --- local -----------------------------------

    character(len=6+MAORI_LEN_NAME)       ::  rckey
    character(len=MAORI_LEN_NAME)         ::  ctype

    ! --- begin -----------------------------------

    ! store name:
    mdata%name = trim(name)
    
    write (gol,'("MAORI:   init data for set `",a,"` ...")') trim(mdata%name); call goPr

    ! init empty:
    nullify( mdata%sample  )
    !nullify( mdata%profile )
    !nullify( mdata%field   )

    ! first part of keys: maori.<name>
    write (rckey,'("maori.",a)') trim(name)

    ! read type:
    call ReadRc( rcF, trim(rckey)//'.type', ctype, status )
    IF_NOTOK_RETURN(status=1)

    ! convert to lower case:
    ctype = goLoCase( trim(ctype) )

    ! type specific stuff:
    select case ( trim(ctype) )
      ! ~ sample
      case ( 'sample'  )
        ! store type:
        mdata%type = MAORI_SAMPLE
        ! storage:
        allocate( mdata%sample )
        ! init type:
        call MAORI_Set_Sample_Data_Init( mdata%sample, trim(name), rcF, trim(rckey), t, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( 'satellite'  )
        ! store type:
        mdata%type = MAORI_SATELLITE
        ! storage:
        allocate( mdata%satellite )
        ! init type:
        call MAORI_Data_satellite_Init( mdata%satellite, trim(name), rcF, trim(rckey), t, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type : `",a,"`")') trim(ctype); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Init


  ! ***


  subroutine MAORI_Set_Data_Start( mdata, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Start
    use MAORI_Set_Satellite, only : MAORI_Data_satellite_Start

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(inout)  ::  mdata
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Data_Start'

    ! --- begin -----------------------------------

    write (gol,'("MAORI:    start data for set `",a,"` ...")') trim(mdata%name); call goPr

    ! type specific stuff:
    select case ( mdata%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Sample_Data_Start( mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! start satellite output:
        call MAORI_Data_satellite_Start( mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Start


  ! ***


  subroutine MAORI_Set_Data_Setup( mdata, t1, t2, the_end, status )
  
    use GO             , only : TDate
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Setup
    use MAORI_Set_Satellite, only : MAORI_Data_satellite_Setup

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(inout)   ::  mdata
    type(TDate), intent(in)                 ::  t1, t2
    logical, intent(in)                     ::  the_end
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Data_Setup'

    ! --- begin -----------------------------------

    ! type specific stuff:
    select case ( mdata%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Sample_Data_Setup( mdata%sample, t1, t2, the_end, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! start satellite output:
        call MAORI_Data_satellite_Setup( mdata%satellite, t1, t2, the_end, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Setup


  ! ***


  subroutine MAORI_Set_Data_Done( mdata, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Done
    use MAORI_Set_Satellite, only : MAORI_Data_satellite_Done

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(inout)  ::  mdata
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Data_Done'

    ! --- begin -----------------------------------

    ! type specific stuff:
    select case ( mdata%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! done with samples:
        call MAORI_Set_Sample_Data_Done( mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( mdata%sample )
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! done with satellites:
        call MAORI_Data_satellite_Done( mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( mdata%satellite )
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Done


  ! ***


  subroutine MAORI_Set_Data_Check( mdata, status )

    use MAORI_Param, only : MAORI_TYPE_MIN, MAORI_TYPE_MAX, MAORI_TYPE_NAME
    use MAORI_Param, only : MAORI_SAMPLE, MAORI_satellite

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(in)     ::  mdata
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Data_Check'

    ! --- begin -----------------------------------

    ! check range:
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("type id ",i6," out of range ",i6," .. ",i6)') &
                      mdata%type, MAORI_TYPE_MIN, MAORI_TYPE_MAX; call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific checks:
    select case ( mdata%type )
      case ( MAORI_SAMPLE )
        ! set allocated ?
        if ( .not. associated(mdata%sample) ) then
          write (gol,'("set `sample` not allocated")'); call goErr
          TRACEBACK; status=1; return
        end if
      case ( MAORI_satellite )
        ! set allocated ?
        if ( .not. associated(mdata%satellite) ) then
          write (gol,'("set `satellite` not allocated")'); call goErr
          TRACEBACK; status=1; return
        end if
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Check



  ! ***


  subroutine MAORI_Set_Data_Inq( mdata, status, name, type, &
                                    nparam, nloc, loc_id_range, nlon, nlat, &
                                    nvar, obs_nvar, &
                                    assim_analyse, assim_analyse_now, assim_rho, &
                                    nvalue, nlocvalue, varname )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Inq
    use MAORI_Set_Satellite, only : MAORI_Data_satellite_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(in)        ::  mdata
    integer, intent(out)                      ::  status

    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  type
    
    integer, intent(out), optional            ::  nparam
    integer, intent(out), optional            ::  nloc
    integer, intent(out), optional            ::  loc_id_range(2)
    integer, intent(out), optional            ::  nlon, nlat
    integer, intent(out), optional            ::  nvar
    integer, intent(out), optional            ::  obs_nvar
    logical, intent(out), optional            ::  assim_analyse
    logical, intent(out), optional            ::  assim_analyse_now
    real, intent(out), optional               ::  assim_rho
    integer, intent(out), optional            ::  nvalue
    integer, intent(out), optional            ::  nlocvalue
    character(len=*), intent(out), optional   ::  varname

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_Data_Inq'

    ! --- begin -----------------------------------

    ! return set name:
    if ( present(name) ) name = trim(mdata%name)

    ! return ouptut set type:
    if ( present(type) ) type = mdata%type

    ! inquire per type:
    select case ( mdata%type )
      ! sample
      case ( MAORI_SAMPLE )
        ! check ...
        if ( any((/present(nlon),present(nlat)/)) ) then
          write (gol,'("optional arguments `nlon` and `nlat` not defined for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! inquire ...
        call MAORI_Set_Sample_Data_Inq( mdata%sample, status, &
                               nparam=nparam, nloc=nloc, loc_id_range=loc_id_range, &
                                nvar=nvar, obs_nvar=obs_nvar, &
                               assim_analyse=assim_analyse, assim_analyse_now=assim_analyse_now, &
                               assim_rho=assim_rho, &
                               nvalue=nvalue, nlocvalue=nlocvalue, varname=varname )
        IF_NOTOK_RETURN(status=1)
      ! satellite
      case ( MAORI_satellite )
        if (present(loc_id_range)) then
          write (gol,'("loc_id_range not supported yet")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! inquire ...
        call MAORI_Data_Satellite_Inq( mdata%satellite, status, &
                               nparam=nparam, nlon=nlon, nlat=nlat, nloc=nloc, nvar=nvar, obs_nvar=obs_nvar, &
                               assim_analyse=assim_analyse, assim_analyse_now=assim_analyse_now, &
                               assim_rho=assim_rho, &
                               nvalue=nvalue, nlocvalue=nlocvalue )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Inq


  ! ********************************************************************
  ! ***
  ! *** locations
  ! ***
  ! ********************************************************************


  subroutine MAORI_Set_Data_Loc_Inq( mdata, iloc, status, &
                                      loc_id, name, code, lon, lat, alt, &
                                      varname )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Loc_Inq
    use MAORI_Set_Satellite, only : MAORI_Data_Satellite_Loc_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    integer, intent(in)                         ::  iloc
    integer, intent(out)                        ::  status

    integer, intent(out), optional              ::  loc_id
    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  code
    real, intent(out), optional                 ::  lon
    real, intent(out), optional                 ::  lat
    real, intent(out), optional                 ::  alt
    character(len=*), intent(out), optional     ::  varname

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Data_Loc_Inq'

    ! --- begin ----------------------------------

    ! inquire per type:
    select case ( mdata%type )
      ! ~ sample output:
      case ( MAORI_SAMPLE )
        ! inquire:
        call MAORI_Set_Sample_Data_Loc_Inq( mdata%sample, iloc, status, &
                               loc_id=loc_id, name=name, code=code, lon=lon, lat=lat, alt=alt, &
                               varname=varname )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite output:
      case ( MAORI_SATELLITE )
        ! check ...
        if ( any((/present(loc_id),present(name),present(code),present(alt)/)) ) then
          write (gol,'("optional arguments `name`, `code`, and `alt` not defined for type `satellite`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! inquire:
        call MAORI_Data_Satellite_Loc_Inq( mdata%satellite, iloc, status, &
                                             lon=lon, lat=lat )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Loc_Inq


  ! ********************************************************************
  ! ***
  ! *** variables
  ! ***
  ! ********************************************************************
  
  
  subroutine MAORI_Set_Data_Var_Inq( mdata, ivar, status, name, unit, nlev )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Var_Inq
    use MAORI_Set_Satellite, only : MAORI_Data_Satellite_Var_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    integer, intent(in)                         ::  ivar
    integer, intent(out)                        ::  status
    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  unit
    integer, intent(out), optional              ::  nlev

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Data_Var_Inq'

    ! --- begin -----------------------------------

    ! inquire per type:
    select case ( mdata%type )
      ! ~ sample output:
      case ( MAORI_SAMPLE )
        ! inquire output set:
        call MAORI_Set_Sample_Data_Var_Inq( mdata%sample, ivar, status, &
                                 name=name, unit=unit, nlev=nlev )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite output:
      case ( MAORI_SATELLITE )
        ! inquire output set:
        call MAORI_Data_Satellite_Var_Inq( mdata%satellite, ivar, status, &
                                 name=name, unit=unit, nlev=nlev )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Var_Inq


  ! ********************************************************************
  ! ***
  ! *** observatoins
  ! ***
  ! ********************************************************************
  
  
  subroutine MAORI_Set_Data_Obs_Get( mdata, obs_ivar, status, iloc, ilon, ilat, y, r, alfa )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Obs_Get
    use MAORI_Set_Satellite, only : MAORI_Data_Satellite_Obs_Get
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    integer, intent(in)                         ::  obs_ivar
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  iloc
    integer, intent(in), optional               ::  ilon, ilat
    real, intent(out), optional                 ::  y       ! measured value
    real, intent(out), optional                 ::  r       ! error std.dev.
    real, intent(out), optional                 ::  alfa    ! screening factor

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Data_Obs_Get'

    ! --- begin -----------------------------------

    ! inquire per type:
    select case ( mdata%type )
      ! ~ sample output:
      case ( MAORI_SAMPLE )
        ! check ...
        if ( any((/present(ilon),present(ilat)/)) ) then
          write (gol,'("arguments ilon,ilat not supported for samples")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! inquire output set:
        call MAORI_Set_Sample_Data_Obs_Get( mdata%sample, obs_ivar, status, &
                                          iloc=iloc, y=y, r=r, alfa=alfa )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite output:
      case ( MAORI_SATELLITE )
        ! inquire output set:
        call MAORI_Data_satellite_Obs_Get( mdata%satellite, obs_ivar, status, &
                                          iloc=iloc, ilon=ilon, ilat=ilat, y=y, r=r, alfa=alfa )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Obs_Get


  ! ***
  
  
  subroutine MAORI_Set_Data_Obs_Put( mdata, obs_ivar, status, iloc, ilon, ilat, astat_ibset )

    use MAORI_Param , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample, only : MAORI_Set_Sample_Data_Obs_Put
    use MAORI_Set_Satellite, only : MAORI_Data_Satellite_Obs_Put

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(inout)       ::  mdata
    integer, intent(in)                         ::  obs_ivar
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  iloc
    integer, intent(in), optional               ::  ilon, ilat
    integer, intent(in), optional               ::  astat_ibset

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Data_Obs_Put'

    ! --- begin -----------------------------------

    ! inquire per type:
    select case ( mdata%type )
      ! ~ sample output:
      case ( MAORI_SAMPLE )
        ! check ...
        if ( any((/present(ilon),present(ilat)/)) ) then
          write (gol,'("arguments ilon,ilat not supported for samples")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! put observation entities:
        call MAORI_Set_Sample_Data_Obs_Put( mdata%sample, obs_ivar, status, &
                                           iloc=iloc, astat_ibset=astat_ibset )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite output:
      case ( MAORI_SATELLITE )
        ! put observation entities:
        call MAORI_Data_Satellite_Obs_Put( mdata%satellite, obs_ivar, status, &
                                           iloc=iloc, ilon=ilon, ilat=ilat, astat_ibset=astat_ibset )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Obs_Put


  ! ********************************************************************
  ! ***
  ! *** parameters
  ! ***
  ! ********************************************************************


  subroutine MAORI_Set_Data_Param_Inq( mdata, ipar, type, name, unit, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Param_Inq
    use MAORI_Set_Satellite, only : MAORI_Data_satellite_Param_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(inout)     ::  mdata
    integer, intent(in)                       ::  ipar
    integer, intent(out)                      ::  type
    character(len=*), intent(out)             ::  name
    character(len=*), intent(out)             ::  unit
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Data_Param_Inq'

    ! --- begin -----------------------------------

    ! inquire per type:
    select case ( mdata%type )
      ! ~ sample output:
      case ( MAORI_SAMPLE )
        ! inquire ...
        call MAORI_Set_Sample_Data_Param_Inq( mdata%sample, ipar, type, name, unit, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite output:
      case ( MAORI_satellite )
        ! inquire ...
        call MAORI_Data_satellite_Param_Inq( mdata%satellite, ipar, type, name, unit, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Param_Inq


  ! ***
  
  
  subroutine MAORI_Set_Data_Param_Put_i( mdata, ipar, ival, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Param_Put
    use MAORI_Set_Satellite, only : MAORI_Data_satellite_Param_Put
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(inout)     ::  mdata
    integer, intent(in)                       ::  ipar
    integer, intent(in)                       ::  ival
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Data_Param_Put_r'

    ! --- begin -----------------------------------

    ! put per type:
    select case ( mdata%type )
      ! ~ sample output:
      case ( MAORI_SAMPLE )
        ! put ...
        call MAORI_Set_Sample_Data_Param_Put( mdata%sample, ipar, ival, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite output:
      case ( MAORI_satellite )
        ! put ...
        call MAORI_Data_satellite_Param_Put( mdata%satellite, ipar, ival, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Param_Put_i


  ! ***
  
  
  subroutine MAORI_Set_Data_Param_Put_r( mdata, ipar, rval, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Param_Put
    use MAORI_Set_Satellite, only : MAORI_Data_satellite_Param_Put
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(inout)     ::  mdata
    integer, intent(in)                       ::  ipar
    real, intent(in)                          ::  rval
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Data_Param_Put_r'

    ! --- begin -----------------------------------

    ! put per type:
    select case ( mdata%type )
      ! ~ sample output:
      case ( MAORI_SAMPLE )
        ! put ...
        call MAORI_Set_Sample_Data_Param_Put( mdata%sample, ipar, rval, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite output:
      case ( MAORI_satellite )
        ! put ...
        call MAORI_Data_satellite_Param_Put( mdata%satellite, ipar, rval, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Param_Put_r


  ! ***
  
  
  subroutine MAORI_Set_Data_Param_Put_s( mdata, ipar, sval, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_Data_Param_Put
    use MAORI_Set_Satellite, only : MAORI_Data_satellite_Param_Put
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Data), intent(inout)     ::  mdata
    integer, intent(in)                       ::  ipar
    character(len=*), intent(in)              ::  sval
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Set_Data_Param_Put_s'

    ! --- begin -----------------------------------

    ! put per type:
    select case ( mdata%type )
      ! ~ sample output:
      case ( MAORI_SAMPLE )
        ! put ...
        call MAORI_Set_Sample_Data_Param_Put( mdata%sample, ipar, sval, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite output:
      case ( MAORI_satellite )
        ! put ...
        call MAORI_Data_satellite_Param_Put( mdata%satellite, ipar, sval, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("not implemented for type ",i6," (",a,")")') &
                            mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_Data_Param_Put_s


  ! ====================================================================
  ! ===
  ! === state sets
  ! ===
  ! ====================================================================


  subroutine MAORI_Set_State_Init( mstate, mdata, name, status )

    use GO, only : TRcFile, Init, Done, ReadRc
    use GO, only : goLoCase
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Init
    use MAORI_Set_Satellite, only : MAORI_State_satellite_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(out)      ::  mstate
    type(T_MAORI_Set_Data), intent(in)        ::  mdata
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Init'

    ! --- local -----------------------------------

    ! --- begin -----------------------------------

    write (gol,'("MAORI:   init state for set `",a,"` ...")') trim(mdata%name); call goPr

    ! store name:
    mstate%name = trim(name)
    
    ! copy type:
    mstate%type = mdata%type
    
    ! init empty:
    nullify( mstate%sample     )
    nullify( mstate%satellite  )

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! storage:
        allocate( mstate%sample )
        ! init type:
        call MAORI_Set_Sample_State_Init( mstate%sample, mdata%sample, trim(mstate%name), status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! storage:
        allocate( mstate%satellite )
        ! init type:
        call MAORI_State_satellite_Init( mstate%satellite, mdata%satellite, trim(mstate%name), status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Init


  ! ***


  subroutine MAORI_Set_State_Done( mstate, mdata, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Done
    use MAORI_Set_Satellite, only : MAORI_State_satellite_Done

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(inout)    ::  mstate
    type(T_MAORI_Set_Data), intent(in)        ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Done'

    ! --- begin ----------------------------------
    
    ! check ...
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! done with samples:
        call MAORI_Set_Sample_State_Done( mstate%sample, mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( mstate%sample )
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! done with satellites:
        call MAORI_State_satellite_Done( mstate%satellite, mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( mstate%satellite )
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Done


  ! ***


  subroutine MAORI_Set_State_Start( mstate, mdata, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_satellite
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Start
    use MAORI_Set_Satellite, only : MAORI_State_satellite_Start

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(inout)    ::  mstate
    type(T_MAORI_Set_Data), intent(in)        ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Start'

    ! --- begin ----------------------------------
    
    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Sample_State_Start( mstate%sample, mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! start satellite output:
        call MAORI_State_satellite_Start( mstate%satellite, mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Start


  ! ***


  subroutine MAORI_Set_State_Setup( mstate, mdata, status )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_satellite
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Setup
    use MAORI_Set_Satellite, only : MAORI_State_satellite_Setup

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(inout)      ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Setup'

    ! --- begin -----------------------------------

    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Sample_State_Setup( mstate%sample, mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! start satellite output:
        call MAORI_State_satellite_Setup( mstate%satellite, mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Setup


  ! ***


  subroutine MAORI_Set_State_Put( mstate, mdata, ivar, status, ilon, ilat, iloc, values )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Put
    use MAORI_Set_Satellite, only : MAORI_State_Satellite_Put

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(inout)      ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    integer, intent(in)                         ::  ivar
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  ilon, ilat, iloc
    real, intent(in), optional                  ::  values(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Put'

    ! --- begin -----------------------------------

    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! check ...
        if ( any((/present(ilon),present(ilat)/)) ) then
          write (gol,'("optional arguments `ilon` and `ilat` not supported for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! put sample output:
        call MAORI_Set_Sample_State_Put( mstate%sample, mdata%sample, ivar, status, &
                                       iloc=iloc, values=values )
        IF_NOTOK_RETURN(status=1)
      ! ~ Satellite:
      case ( MAORI_SATELLITE  )
        ! put Satellite output:
        call MAORI_State_Satellite_Put( mstate%satellite, mdata%satellite, ivar, status, &
                                       ilon=ilon, ilat=ilat, iloc=iloc, values=values )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Put


  ! ***


  subroutine MAORI_Set_State_Post( mstate, mdata, status )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Post
    use MAORI_Set_Satellite, only : MAORI_State_satellite_Post

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(inout)      ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Post'

    ! --- begin -----------------------------------

    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Sample_State_Post( mstate%sample, mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! start satellite output:
        call MAORI_State_satellite_Post( mstate%satellite, mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Post
  

  ! ***


  subroutine MAORI_Set_State_Obs_Get( mstate, mdata, obs_ivar, status, &
                                         iloc, ilon, ilat, value )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Obs_Get
    use MAORI_Set_Satellite, only : MAORI_State_Satellite_Obs_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(in)         ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    integer, intent(in)                         ::  obs_ivar
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  iloc
    integer, intent(in), optional               ::  ilon, ilat
    real, intent(out), optional                 ::  value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Obs_Get'

    ! --- begin ----------------------------------

    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! check ...
        if ( any((/present(ilon),present(ilat)/)) ) then
          write (gol,'("optional arguments `ilon` and `ilat` not supported for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! start sample output:
        call MAORI_Set_Sample_State_Obs_Get( mstate%sample, mdata%sample, obs_ivar, status, &
                                          iloc=iloc, value=value )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_SATELLITE  )
        ! start satellite output:
        call MAORI_State_Satellite_Obs_Get( mstate%satellite, mdata%satellite, obs_ivar, status, &
                                          iloc=iloc, ilon=ilon, ilat=ilat, value=value )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Obs_Get
  

  ! ***


  subroutine MAORI_Set_State_Values_Get( mstate, mdata, values, status )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Values_Get
    use MAORI_Set_Satellite, only : MAORI_State_Satellite_Values_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(in)         ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    real, intent(out)                           ::  values(:)
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Values_Get'

    ! --- begin ----------------------------------

    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Sample_State_Values_Get( mstate%sample, mdata%sample, values, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_SATELLITE  )
        ! start satellite output:
        call MAORI_State_Satellite_Values_Get( mstate%satellite, mdata%satellite, values, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Values_Get


  ! ***


  subroutine MAORI_Set_State_Values_Put( mstate, mdata, values, status )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_Values_Put
    use MAORI_Set_Satellite, only : MAORI_State_Satellite_Values_Put

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(inout)      ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    real, intent(in)                            ::  values(:)
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_Values_Put'

    ! --- begin ----------------------------------

    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Sample_State_Values_Put( mstate%sample, mdata%sample, values, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ Satellite:
      case ( MAORI_SATELLITE  )
        ! start Satellite output:
        call MAORI_State_Satellite_Values_Put( mstate%satellite, mdata%satellite, values, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_Values_Put
  

  ! ***


  subroutine MAORI_Set_State_LocValues_Get( mstate, mdata, values, status, iloc, ilon, ilat )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_LocValues_Get
    use MAORI_Set_Satellite, only : MAORI_State_Satellite_LocValues_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(in)         ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    real, intent(out)                           ::  values(:)
    integer, intent(out)                        ::  status
    integer, intent(in), optional               ::  iloc
    integer, intent(in), optional               ::  ilon, ilat

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_LocValues_Get'

    ! --- begin ----------------------------------

    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! check ...
        if ( .not. all((/present(iloc)/)) ) then
          write (gol,'("optional arguments `iloc` should be provided for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! check ...
        if ( any((/present(ilon),present(ilat)/)) ) then
          write (gol,'("optional arguments `ilon` and `ilat` not supported for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! start sample output:
        call MAORI_Set_Sample_State_LocValues_Get( mstate%sample, mdata%sample, iloc, values, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_SATELLITE  )
        ! check ...
        if ( .not. all((/present(ilon),present(ilat)/)) ) then
          write (gol,'("optional arguments `ilon` and `ilat` should be provided for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! check ...
        if ( any((/present(iloc)/)) ) then
          write (gol,'("optional argument `iloc` not supported for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! start satellite output:
        call MAORI_State_Satellite_LocValues_Get( mstate%satellite, mdata%satellite, ilon, ilat, values, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_LocValues_Get


  ! ***


  subroutine MAORI_Set_State_LocValues_Put( mstate, mdata, values, status, iloc, ilon, ilat )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME, MAORI_TYPE_MIN, MAORI_TYPE_MAX
    use MAORI_Set_Sample   , only : MAORI_Set_Sample_State_LocValues_Put
    use MAORI_Set_Satellite, only : MAORI_State_Satellite_LocValues_Put

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_State), intent(inout)      ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    real, intent(in)                            ::  values(:)
    integer, intent(out)                        ::  status
    integer, intent(in), optional               ::  iloc
    integer, intent(in), optional               ::  ilon, ilat

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Set_State_LocValues_Put'

    ! --- begin ----------------------------------

    ! check ...
    if ( (mdata%type < MAORI_TYPE_MIN) .or. (mdata%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mdata out of range : ",i8)') mdata%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( (mstate%type < MAORI_TYPE_MIN) .or. (mstate%type > MAORI_TYPE_MAX) ) then
      write (gol,'("MAORI type of mstate out of range : ",i8)') mstate%type; call goErr
      TRACEBACK; status=1; return
    end if
    if ( mdata%type /= mstate%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  state type : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type   : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( mstate%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! check ...
        if ( .not. all((/present(iloc)/)) ) then
          write (gol,'("optional arguments `iloc` should be provided for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! check ...
        if ( any((/present(ilon),present(ilat)/)) ) then
          write (gol,'("optional arguments `ilon` and `ilat` not supported for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! start sample output:
        call MAORI_Set_Sample_State_LocValues_Put( mstate%sample, mdata%sample, iloc, values, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_SATELLITE  )
        ! check ...
        if ( .not. all((/present(ilon),present(ilat)/)) ) then
          write (gol,'("optional arguments `ilon` and `ilat` should be provided for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! check ...
        if ( any((/present(iloc)/)) ) then
          write (gol,'("optional argument `iloc` not supported for type `sample`")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! start satellite output:
        call MAORI_State_Satellite_LocValues_Put( mstate%satellite, mdata%satellite, ilon, ilat, values, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Set_State_LocValues_Put


  ! ====================================================================
  ! ===
  ! === state sets
  ! ===
  ! ====================================================================


  subroutine MAORI_Output_Set_Init( moutput, mdata, name, status )

    use GO, only : TRcFile, Init, Done, ReadRc
    use GO, only : goLoCase
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Output_State_Init
    use MAORI_Set_Satellite, only : MAORI_Output_satellite_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output), intent(out)     ::  moutput
    type(T_MAORI_Set_Data), intent(in)        ::  mdata
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Set_Init'

    ! --- local -----------------------------------

    ! --- begin -----------------------------------

    write (gol,'("MAORI:   init output for set `",a,"` ...")') trim(mdata%name); call goPr

    ! store name:
    moutput%name = trim(name)
    
    ! copy type:
    moutput%type = mdata%type
    
    ! init empty:
    nullify( moutput%sample  )

    ! type specific stuff:
    select case ( moutput%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! storage:
        allocate( moutput%sample )
        ! init type:
        call MAORI_Set_Output_State_Init( moutput%sample, mdata%sample, trim(moutput%name), status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! storage:
        allocate( moutput%satellite )
        ! init type:
        call MAORI_Output_satellite_Init( moutput%satellite, mdata%satellite, trim(moutput%name), status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            moutput%type, trim(MAORI_TYPE_NAME(moutput%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Output_Set_Init


  ! ***


  subroutine MAORI_Output_Set_Done( moutput, mdata, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Output_State_Done
    use MAORI_Set_Satellite, only : MAORI_Output_satellite_Done

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output), intent(inout)   ::  moutput
    type(T_MAORI_Set_Data), intent(in)        ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Set_Done'

    ! --- begin ----------------------------------

    ! check ...
    if ( mdata%type /= moutput%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  output type : ",i6," (",a,")")') moutput%type, trim(MAORI_TYPE_NAME(moutput%type)); call goErr
      write (gol,'("  set type    : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( moutput%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! done with samples:
        call MAORI_Set_Output_State_Done( moutput%sample, mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( moutput%sample )
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! done with satellites:
        call MAORI_Output_satellite_Done( moutput%satellite, mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
        ! clear:
        deallocate( moutput%satellite )
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            moutput%type, trim(MAORI_TYPE_NAME(moutput%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Output_Set_Done


  ! ***


  subroutine MAORI_Output_Set_Start( moutput, mdata, status )

    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Output_State_Start
    use MAORI_Set_Satellite, only : MAORI_Output_satellite_Start

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output), intent(inout)    ::  moutput
    type(T_MAORI_Set_Data), intent(in)        ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Set_Start'

    ! --- begin ----------------------------------
    
    ! check ...
    if ( mdata%type /= moutput%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  output type : ",i6," (",a,")")') moutput%type, trim(MAORI_TYPE_NAME(moutput%type)); call goErr
      write (gol,'("  set type    : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( moutput%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Output_State_Start( moutput%sample, mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! start satellite output:
        call MAORI_Output_satellite_Start( moutput%satellite, mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            moutput%type, trim(MAORI_TYPE_NAME(moutput%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Output_Set_Start


  ! ***


  subroutine MAORI_Output_Set_Write( moutput, mstate, mdata, status )
  
    use MAORI_Param    , only : MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI_Param    , only : MAORI_TYPE_NAME
    use MAORI_Set_Sample   , only : MAORI_Set_Output_State_Write
    use MAORI_Set_Satellite, only : MAORI_Output_satellite_Write

    ! --- in/out ---------------------------------

    type(T_MAORI_Set_Output), intent(inout)     ::  moutput
    type(T_MAORI_Set_State), intent(in)         ::  mstate
    type(T_MAORI_Set_Data), intent(in)          ::  mdata
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Set_Write'

    ! --- begin -----------------------------------

    ! check ...
    if ( mdata%type /= moutput%type ) then
      write (gol,'("MAORI types of arguments should match:")'); call goErr
      write (gol,'("  output type : ",i6," (",a,")")') moutput%type, trim(MAORI_TYPE_NAME(moutput%type)); call goErr
      write (gol,'("  state type  : ",i6," (",a,")")') mstate%type, trim(MAORI_TYPE_NAME(mstate%type)); call goErr
      write (gol,'("  set type    : ",i6," (",a,")")') mdata%type, trim(MAORI_TYPE_NAME(mdata%type)); call goErr
      TRACEBACK; status=1; return
    end if

    ! type specific stuff:
    select case ( moutput%type )
      ! ~ sample:
      case ( MAORI_SAMPLE  )
        ! start sample output:
        call MAORI_Set_Output_State_Write( moutput%sample, mstate%sample, mdata%sample, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ satellite:
      case ( MAORI_satellite  )
        ! start satellite output:
        call MAORI_Output_satellite_Write( moutput%satellite, mstate%satellite, mdata%satellite, status )
        IF_NOTOK_RETURN(status=1)
      ! ~ unknown ...
      case default
        write (gol,'("unsupported set type ",i6," (",a,")")') &
                            moutput%type, trim(MAORI_TYPE_NAME(moutput%type)); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine MAORI_Output_Set_Write


end module MAORI_Set
