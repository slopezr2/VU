!###############################################################################
!
! LOTOS-EUROS model data
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

module LE_Data_VData

  use GO  , only : gol, goPr, goErr, goLabel
  use LE_Data_Variables, only : T_Variables
  
  implicit none


  ! --- in/out -----------------------------------
  
  private

  public  ::  LE_Data_VData_Init, LE_Data_VData_Done  
  public  ::  LE_Data_VData_Inquire
  public  ::  LE_Data_VData_InqVar
  public  ::  LE_Data_VData_Enable
  public  ::  LE_Data_VData_GetPointer
  public  ::  LE_Data_VData_Setup


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Data_VData'
  

  ! --- interfaces -------------------------------
  
  interface LE_Data_VData_InqVar
    module procedure  ::  LE_Data_VData_InqVar_ivar
    module procedure  ::  LE_Data_VData_InqVar_name
  end interface LE_Data_VData_InqVar
  

  ! --- var --------------------------------------
  
  ! variables:
  type(T_Variables)             ::  vdata  
  

contains


  ! ====================================================================
  ! ===
  ! === module
  ! ===
  ! ====================================================================
  

  subroutine LE_Data_VData_Init( rcF, status )
  
    use GO     , only : TrcFile, ReadRc
    use GO     , only : goSplitString
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), pointer                  ::  rcF
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Data_VData_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------
    
    ! init variables:
    call vdata%Init( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Data_VData_Init
  
  
  ! ***
  
  
  subroutine LE_Data_VData_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_VData_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! done with variables:
    call vdata%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Data_VData_Done


  ! ====================================================================
  ! ===
  ! === data access
  ! ===
  ! ====================================================================
  
  
  subroutine LE_Data_VData_Enable( name, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)            ::  name
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_VData_Enable'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
  
    ! enable:
    call vdata%SetVariable( trim(name), status, enabled=.true. )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Data_VData_Enable
  
  
  ! ***
    
  
  subroutine LE_Data_VData_Inquire( status, nvar )

    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    integer, intent(out), optional          ::  nvar
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_VData_Inquire'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! number of variables?
    if ( present(nvar) ) nvar = vdata%nvar
    
    ! ok
    status = 0
    
  end subroutine LE_Data_VData_Inquire
    
  
  ! *
  
  
  !
  ! Arguments:
  !  - ivar   :  index 1,..,nvar, could return name
  !  - name 
  !
  
  subroutine LE_Data_VData_InqVar_ivar( ivar, status, &
                                          name, long_name, units, &
                                          rank, const, levtype )

    ! --- in/out ---------------------------------
    
    integer, intent(in)                       ::  ivar
    integer, intent(out)                      ::  status
    
    character(len=*), intent(out), optional   ::  name
    character(len=*), intent(out), optional   ::  long_name
    character(len=*), intent(out), optional   ::  units
    integer, intent(out), optional            ::  rank
    logical, intent(out), optional            ::  const
    character(len=*), intent(out), optional   ::  levtype
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_VData_InqVar_ivar'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( (ivar < 1) .or. (ivar > vdata%nvar) ) then
      write (gol,'("variable index ",i0," out of range 1,..,",i0)') ivar, vdata%nvar; call goErr
      TRACEBACK; status=1; return
    end if

    ! copy:
    if ( present(name     ) ) name      = vdata%var(ivar)%name
    if ( present(long_name) ) long_name = vdata%var(ivar)%long_name
    if ( present(units    ) ) units     = vdata%var(ivar)%units
    if ( present(rank     ) ) rank      = vdata%var(ivar)%rank
    if ( present(const    ) ) const     = vdata%var(ivar)%const
    if ( present(levtype  ) ) levtype   = vdata%var(ivar)%levtype
    
    ! ok
    status = 0
    
  end subroutine LE_Data_VData_InqVar_ivar
  
  ! *
  
  subroutine LE_Data_VData_InqVar_name( name, status, &
                                          long_name, units, &
                                          rank, const, levtype )
  
    use LE_Data_Variable, only : T_Variable

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  status
    
    character(len=*), intent(out), optional   ::  long_name
    character(len=*), intent(out), optional   ::  units
    integer, intent(out), optional            ::  rank
    logical, intent(out), optional            ::  const
    character(len=*), intent(out), optional   ::  levtype
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_VData_InqVar_name'
    
    ! --- local ----------------------------------
    
    class(T_Variable), pointer    ::  V
    logical                       ::  enabled
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call vdata%GetVariable( name, status, V=V )
    IF_NOTOK_RETURN(status=1)

    ! get flag:
    call V%Get( status, enabled=enabled )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    if ( .not. enabled ) then
      write (gol,'("requested properties from variable `",a,"` which is not enabled")') trim(name); call goErr
      TRACEBACK; status=1; return
    end if

    ! copy:
    if ( present(long_name) ) long_name = trim(V%long_name)
    if ( present(units    ) ) units     = trim(V%units)
    if ( present(rank     ) ) rank      = V%rank
    if ( present(const    ) ) const     = V%const
    if ( present(levtype  ) ) levtype   = V%levtype
    
    ! ok
    status = 0
    
  end subroutine LE_Data_VData_InqVar_name
  
  
  ! ***
    
  
  subroutine LE_Data_VData_GetPointer( name, pdata, status, component, &
                                         units, check_units, check_lbo, check_ubo )

    use LE_Data_Common  , only : LE_Data_CompareUnits
    use LE_Data_Variable, only : T_Variable

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)              ::  name
    real, pointer                             ::  pdata(:,:,:)
    integer, intent(out)                      ::  status
    character(len=1), intent(in), optional    ::  component
    character(len=*), intent(out), optional   ::  units
    character(len=*), intent(in), optional    ::  check_units
    integer, intent(in), optional             ::  check_lbo(3), check_ubo(3)
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_VData_GetPointer'
    
    ! --- local ----------------------------------
    
    class(T_Variable), pointer    ::  V
    logical                       ::  enabled
    character(len=1)              ::  comp
    
    ! --- begin ----------------------------------
    
    ! vector component?
    comp = '-'
    if ( present(component) ) comp = trim(component)
    
    ! get variable by name:
    call vdata%GetVariable( name, status, V=V )
    IF_NOTOK_RETURN(status=1)

    ! get flag:
    call V%Get( status, enabled=enabled )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    if ( .not. enabled ) then
      write (gol,'("requested pointer from variable `",a,"` which is not enabled")') trim(name); call goErr
      TRACEBACK; status=1; return
    end if

    ! standard data or vector component?
    select case ( trim(comp) )
      !~ standard data:
      case ( '-' )
        ! get pointer to data (and units), perform requested checks:
        call V%Get( status, pdata=pdata, units=units, &
                      check_units=check_units, &
                      check_lbo=check_lbo, check_ubo=check_ubo )
        IF_NOTOK_RETURN(status=1)
      !~ vector component:
      case ( 'u' )
        ! get pointer to data (and units), perform requested checks:
        call V%Get( status, pudata=pdata, units=units, &
                      check_units=check_units, &
                      check_lbo=check_lbo, check_ubo=check_ubo )
        IF_NOTOK_RETURN(status=1)
      !~ vector component:
      case ( 'v' )
        ! get pointer to data (and units), perform requested checks:
        call V%Get( status, pvdata=pdata, units=units, &
                      check_units=check_units, &
                      check_lbo=check_lbo, check_ubo=check_ubo )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unusupported component `",a,"`")') trim(comp); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! check ...
    if ( .not. associated(pdata) ) then
      write (gol,'("variable `",a,"` pdata not associated")') trim(name); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine LE_Data_VData_GetPointer
  
  
  ! ***
  
  
  subroutine LE_Data_VData_Setup( tref, t1, t2, status, var_name )
  
    use GO  , only : TDate

    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)                 ::  tref
    type(TDate), intent(in)                 ::  t1, t2
    integer, intent(out)                    ::  status
    character(len=*), intent(in), optional  ::  var_name   ! only setup/update for specific varname
        
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_VData_Setup'

    ! order of time stamps:
    integer, parameter           ::  ntime = 4
    character(len=1), parameter  ::  times(ntime) = (/ 'p', 'c', 'n', 'a' /)
    
    !! testing ...
    !character(len=*), parameter  ::  vnames(3) = (/ 'oro', 'lsm', 'vol' /)
    
    ! --- local ----------------------------------
    
    integer               ::  itime
    integer               ::  ivar
    logical               ::  enabled
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a,": setup data for current time interval ...")') rname; call goPr

    ! loop over time stamps:
    do itime = 1, ntime

      !! info ...
      !write (gol,'(a,":  setup for time `",a,"`")') rname, times(itime); call goPr

      ! loop over all variables:
      do ivar = 1, vdata%nvar

        ! filter on time stamp:
        if ( trim(vdata%var(ivar)%time) /= times(itime)  ) cycle
        
        ! filter on specific var_name, not present continue for all variables
        if ( present(var_name) ) then
          ! check var name 
          if ( .not. trim(vdata%var(ivar)%name) == trim(var_name) ) cycle
        end if
        
        ! get flag:
        call vdata%var(ivar)%Get( status, enabled=enabled )
        IF_NOTOK_RETURN(status=1)
        ! skip if not enabled:
        if ( .not. enabled ) cycle

        ! switch on time stamp:
        select case ( times(itime) )

          !~ previous:
          case ( 'p' )
            ! setup for start of time interval:
            call vdata%Setup( trim(vdata%var(ivar)%name), tref, t1, t1, status )
            IF_NOTOK_RETURN(status=1)

          !~ current:
          case ( 'c' )
            ! setup for time interval (or mid):
            call vdata%Setup( trim(vdata%var(ivar)%name), tref, t1, t2, status )
            IF_NOTOK_RETURN(status=1)

          !~ next:
          case ( 'n' )
            ! setup for end of time interval:
            call vdata%Setup( trim(vdata%var(ivar)%name), tref, t2, t2, status )
            IF_NOTOK_RETURN(status=1)

          !~ all (p, c, and n should have been setup):
          case ( 'a' )
            ! setup for interval (or mid);
            ! do not update dependencies because fields now setup 
            ! for 'previous' time should remain valid for that time:
            call vdata%Setup( trim(vdata%var(ivar)%name), tref, t1, t2, status, nodeps=.true. )
            IF_NOTOK_RETURN(status=1)

          !~
          case default
            write (gol,'("unsupported time stamp `",a,"`")') times(itime); call goErr
            TRACEBACK; status=1; return
        end select

      end do ! variables
    
    end do ! time stamps
    
    ! ok
    status = 0
    
  end subroutine LE_Data_VData_Setup
  
  
end module LE_Data_VData

