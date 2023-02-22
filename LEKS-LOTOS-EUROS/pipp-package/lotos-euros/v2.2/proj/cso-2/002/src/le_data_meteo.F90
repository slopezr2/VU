!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Data_Meteo

  use GO     , only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Data_Setup_cFS_sfc
  public  ::  LE_Data_Setup_cFS_3D
  public  ::  LE_Data_Setup_aFS_sfc
  public  ::  LE_Data_Setup_iFS_sfc
  public  ::  LE_Data_Setup_iFS_3D
  public  ::  LE_Data_Setup_iVFS_sfc
  public  ::  LE_Data_Setup_iVFS_3D


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Data_Meteo'


contains


  ! ====================================================================
  ! ===
  ! === module
  ! ===
  ! ====================================================================
  

  subroutine LE_Data_Setup_cFS_sfc( mf, tref, tt, rcF, gridtype, mapping, status )
  
    use GO     , only : TDate, wrtgol, IsAnyDate
    use GO     , only : TDate, operator(+), operator(-), operator(*)
    use GO     , only : TrcFile
    use GO     , only : T_Constant_Field_Series
    use LE_Grid, only : ugg, ugg_bnds, ugg_crnr, ugg_crnr_bnds
    use C3PO   , only : T_Grid_Ugg
    
    ! --- in/out ---------------------------------
    
    type(T_Constant_Field_Series), intent(inout)    ::  mf
    type(TDate), intent(in)                         ::  tref
    type(TDate), intent(in)                         ::  tt(2)
    type(TrcFile), intent(in)                       ::  rcF
    character(len=*), intent(in)                    ::  gridtype
    character(len=*), intent(in)                    ::  mapping
    integer, intent(out)                            ::  status
     
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Setup_cFS_sfc'
    
    ! --- local ----------------------------------
    
    type(T_Grid_Ugg), pointer   ::  gri
    character(len=64)           ::  gmapping
    character(len=1024)         ::  inputs, input
    type(TDate)                 ::  tmid
    type(TDate)                 ::  tt_input(2)
    type(TDate)                 ::  tt_rec(2)
    character(len=64)           ::  name, units
    real, allocatable           ::  mdata(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'corner' )
        gri => ugg_crnr
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case ( 'corner_bnds' )
        gri => ugg_crnr_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! set grid mapping:
    if ( trim(mapping) == 'None' ) then
      ! reset default based on grid type:
      select case ( trim(gridtype) )
        case ( 'cell', 'cell_bnds' )
          gmapping = 'area-aver'
        case ( 'corner', 'corner_bnds' )
          gmapping = 'interpol'
        case default
          write (gol,'("could not set grid mapping for grid type: ",a)') trim(gridtype); call goErr
          TRACEBACK; status=1; return
      end select
    else
      ! trap incorrect use of area-aver for corner variables
      select case ( trim(gridtype) )
        case ( 'corner', 'corner_bnds' )
          if ( trim(mapping) == 'area-aver' ) then
            write (gol,'("Do not use area-average mapping for variables with gridtype: ",a)') trim(gridtype) ; call goErr
            TRACEBACK;status=1;return
          end if
        case default
          ! ~ no action needed
      end select
      
      ! copy:
      gmapping = trim(mapping)
    end if
    
    ! setup for requested time interval; return status -1 if new data should be put:
    call mf%Setup( tt, status )
    if ( status > 0 ) then
      ! error status ...
      TRACEBACK; status=1; return
    else if ( status < 0 ) then

      ! new data needed; get input description and other:
      call mf%Get( status, input=inputs, name=name, units=units )
      IF_NOTOK_RETURN(status=1)
      ! info ...
      call wrtgol( '    read field `'//trim(name)//'` valid for ', tt(1), ' - ', tt(2) ); call goPr
      ! mid time:
      tmid = tt(1) + (tt(2)-tt(1))*0.5
      ! extract input description valid for requested time,
      ! returns also time range for which description is valid:
      call GetInput( inputs, tmid, tt_input, input, status )
      IF_NOTOK_RETURN(status=1)

      ! obtain data from somewhere, allocate mdata, and store into there;
      ! time interval is used to form filename;
      ! output interval 'tt_rec' defines for which time range record is valid,
      ! equal to (/AnyDate,AnyDate/) for time independed variable:
      call MeteoField_ReadRecord_sfc( tref, tt, rcF, input, units, gri, gmapping, &
                                       mdata, tt_rec, status )
      IF_NOTOK_RETURN(status=1)
      ! no interval defined?
      if ( IsAnyDate(tt_rec(1)) ) then
        ! valid for the interval specified in the input description:
        tt_rec = tt_input
      end if
      ! store field:
      call mf%Setup_Put( mdata, tt_rec, status )
      IF_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( mdata, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! new data neeeded

    ! ok
    status = 0
    
  end subroutine LE_Data_Setup_cFS_sfc


  ! ***
  
  
  subroutine LE_Data_Setup_cFS_3D( mf, tref, tt, rcF, gridtype, mapping, status, halflevels )
  
    use GO     , only : TDate, wrtgol, IsAnyDate
    use GO     , only : TDate, operator(+), operator(-), operator(*)
    use GO     , only : TrcFile
    use GO     , only : T_Constant_Field_Series
    use GO     , only : goReadFromLine
    use LE_Grid, only : ugg, ugg_bnds, ugg_crnr, ugg_crnr_bnds
    use C3PO   , only : T_Grid_Ugg
    
    ! --- in/out ---------------------------------
    
    type(T_Constant_Field_Series), intent(inout)    ::  mf
    type(TDate), intent(in)                         ::  tref
    type(TDate), intent(in)                         ::  tt(2)
    type(TrcFile), intent(in)                       ::  rcF
    character(len=*), intent(in)                    ::  gridtype
    character(len=*), intent(in)                    ::  mapping
    integer, intent(out)                            ::  status
    logical, intent(in), optional                   ::  halflevels
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Setup_cFS_3D'
    
    ! --- local ----------------------------------
    
    type(T_Grid_Ugg), pointer   ::  gri
    character(len=64)           ::  gmapping
    character(len=1024)         ::  inputs, input
    type(TDate)                 ::  tmid
    type(TDate)                 ::  tt_input(2)
    type(TDate)                 ::  tt_rec(2)
    character(len=64)           ::  name, units
    real, allocatable           ::  mdata(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'corner' )
        gri => ugg_crnr
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case ( 'corner_bnds' )
        gri => ugg_crnr_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set grid mapping:
    if ( trim(mapping) == 'None' ) then
      ! reset default based on grid type:
      select case ( trim(gridtype) )
        case ( 'cell', 'cell_bnds' )
          gmapping = 'area-aver'
        case ( 'corner', 'corner_bnds' )
          gmapping = 'interpol'
        case default
          write (gol,'("could not set grid mapping for grid type: ",a)') trim(gridtype); call goErr
          TRACEBACK; status=1; return
      end select
    else
      ! trap incorrect use of area-aver for corner variables
      select case ( trim(gridtype) )
        case ( 'corner', 'corner_bnds' )
          if ( trim(mapping) == 'area-aver' ) then
            write (gol,'("Do not use area-average mapping for variables with gridtype: ",a)') trim(gridtype) ; call goErr
            TRACEBACK;status=1;return
          end if
        case default
          ! ~ no action needed
      end select
      ! copy:
      gmapping = trim(mapping)
    end if
    
    ! setup for requested time interval; return status -1 if new data should be put:
    call mf%Setup( tt, status )
    if ( status > 0 ) then
      TRACEBACK; status=1; return
    else if ( status < 0 ) then

      ! new data needed; get input description and other:
      call mf%Get( status, input=inputs, name=name, units=units )
      IF_NOTOK_RETURN(status=1)
      ! info ...
      call wrtgol( '    read field `'//trim(name)//'` valid for ', tt(1), ' - ', tt(2) ); call goPr
      ! mid time:
      tmid = tt(1) + (tt(2)-tt(1))*0.5
      ! extract input description valid for requested time,
      ! returns also time range for which description is valid:
      call GetInput( inputs, tmid, tt_input, input, status )
      IF_NOTOK_RETURN(status=1)

      ! obtain data from somewhere, allocate mdata, and store into there;
      ! time interval is used to form filename;
      ! output interval 'tt_rec' defines for which time range record is valid,
      ! equal to (/AnyDate,AnyDate/) for time independed variable:
      call MeteoField_ReadRecord_3D( tref, tt, rcF, input, units, &
                                       gri, gmapping, &
                                       mdata, tt_rec, status, &
                                       halflevels=halflevels )
      IF_NOTOK_RETURN(status=1)
      ! no interval defined?
      if ( IsAnyDate(tt_rec(1)) ) then
        ! valid for the interval specified in the input description:
        tt_rec = tt_input
      end if
      ! store field:
      call mf%Setup_Put( mdata, tt_rec, status, &
                          lbo=lbound(mdata), ubo=ubound(mdata) )
      IF_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( mdata, stat=status )
      IF_NOTOK_RETURN(status=1)
      
    end if ! new data neeeded

    ! ok
    status = 0
    
  end subroutine LE_Data_Setup_cFS_3D
  
  
  ! ***
  
  
  subroutine LE_Data_Setup_aFS_sfc( mf, tref, tt, rcF, gridtype, mapping, status )
  
    use GO     , only : TDate, wrtgol, IsAnyDate, operator(-), IncrDate
    use GO     , only : TrcFile
    use GO     , only : T_Accumulated_Field_Series
    use LE_Grid, only : ugg, ugg_bnds, ugg_crnr, ugg_crnr_bnds
    use C3PO   , only : T_Grid_Ugg

    ! --- in/out ---------------------------------
    
    type(T_Accumulated_Field_Series), intent(inout) ::  mf
    type(TDate), intent(in)                         ::  tref
    type(TDate), intent(in)                         ::  tt(2)
    type(TrcFile), intent(in)                       ::  rcF
    character(len=*), intent(in)                    ::  gridtype
    character(len=*), intent(in)                    ::  mapping
    integer, intent(out)                            ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Setup_aFS_sfc'
    
    ! --- local ----------------------------------
    
    type(T_Grid_Ugg), pointer   ::  gri
    character(len=64)           ::  gmapping
    character(len=1024)         ::  inputs, input
    type(TDate)                 ::  tt_input(2)
    type(TDate)                 ::  tt_rec(2)
    character(len=64)           ::  name, units
    integer                     ::  nreceive, ireceive
    type(TDate)                 ::  receive_t
    real, allocatable           ::  mdata(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case ( 'corner' )
        gri => ugg_crnr
      case ( 'corner_bnds' )
        gri => ugg_crnr_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set grid mapping:
    if ( trim(mapping) == 'None' ) then
      ! reset default based on grid type:
      select case ( trim(gridtype) )
        case ( 'cell', 'cell_bnds' )
          gmapping = 'area-aver'
        case ( 'corner', 'corner_bnds' )
          gmapping = 'interpol'
        case default
          write (gol,'("could not set grid mapping for grid type: ",a)') trim(gridtype); call goErr
          TRACEBACK; status=1; return
      end select
    else
      ! trap incorrect use of area-aver for corner variables
      select case ( trim(gridtype) )
        case ( 'corner', 'corner_bnds' )
          if ( trim(mapping) == 'area-aver' ) then
            write (gol,'("Do not use area-average mapping for variables with gridtype: ",a)') trim(gridtype) ; call goErr
            TRACEBACK;status=1;return
          end if
        case default
          ! ~ no action needed
      end select
      ! copy:
      gmapping = trim(mapping)
    end if
    
    ! setup for time t; return status -1 if new data should be put:
    call mf%Setup( tt, status )
    if ( status > 0 ) then
      ! error status ...
      TRACEBACK; status=1; return
    else if ( status < 0 ) then

      ! prepare to receive new data:
      call mf%Setup_Prepare( tref, nreceive, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over fields to be received:
      do ireceive = 1, nreceive

        ! get input description:
        call mf%Get( status, input=inputs, name=name, units=units )
        IF_NOTOK_RETURN(status=1)
        ! get time value:
        call mf%Setup_InqTime( ireceive, receive_t, status )
        IF_NOTOK_RETURN(status=1)
        
        ! previous field might have been skipped by Setup_Put for the current field:
        if ( IsAnyDate(receive_t) ) then
          ! info ...
          write (gol,'("    skip `",a,"`, no previous data needed ...")') trim(name); call goPr
          ! leave:
          exit
        end if
        
        ! info ...
        call wrtgol( '    read `'//trim(name)//'` valid for ', receive_t ); call goPr
        ! extract input description valid for requested time:
        call GetInput( inputs, receive_t, tt_input, input, status )
        IF_NOTOK_RETURN(status=1)

        ! obtain data for interval ending at target time, just subtract one minute to have an interval;
        ! allocate mdata, and store into there;
        ! return interval 'tt_rec' with accumulation period:
        call MeteoField_ReadRecord_sfc( tref, (/receive_t-IncrDate(min=1),receive_t/), &
                                         rcF, input, units, &
                                         gri, gmapping, &
                                         mdata, tt_rec, status )
        IF_NOTOK_RETURN(status=1)
        ! store:
        call mf%Setup_Put( ireceive, mdata, tt_rec, status )
        IF_NOTOK_RETURN(status=1)
        
      end do  ! fields to receive

      ! clear:
      deallocate( mdata, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! new data neeeded

    ! ok
    status = 0
    
  end subroutine LE_Data_Setup_aFS_sfc
  
  
  ! ***
  
  
  subroutine LE_Data_Setup_iFS_sfc( mf, tref, t, rcF, gridtype, mapping, status )
  
    use GO     , only : TDate, wrtgol
    use GO     , only : TrcFile
    use GO     , only : T_Instant_Field_Series
    use LE_Grid, only : ugg, ugg_bnds, ugg_crnr, ugg_crnr_bnds
    use C3PO   , only : T_Grid_Ugg

    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(inout)     ::  mf
    type(TDate), intent(in)                         ::  tref
    type(TDate), intent(in)                         ::  t
    type(TrcFile), intent(in)                       ::  rcF
    character(len=*), intent(in)                    ::  gridtype
    character(len=*), intent(in)                    ::  mapping
    integer, intent(out)                            ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Setup_iFS_sfc'
    
    ! --- local ----------------------------------
    
    type(T_Grid_Ugg), pointer  ::  gri
    character(len=64)           ::  gmapping
    character(len=1024)         ::  inputs, input
    type(TDate)                 ::  tt(2)
    type(TDate)                 ::  tt_rec(2)
    character(len=64)           ::  name, units
    integer                     ::  nreceive, ireceive
    type(TDate)                 ::  receive_t
    real, allocatable           ::  mdata(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case ( 'corner' )
        gri => ugg_crnr
      case ( 'corner_bnds' )
        gri => ugg_crnr_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set grid mapping:
    if ( trim(mapping) == 'None' ) then
      ! reset default based on grid type:
      select case ( trim(gridtype) )
        case ( 'cell', 'cell_bnds' )
          gmapping = 'area-aver'
        case ( 'corner', 'corner_bnds' )
          gmapping = 'interpol'
        case default
          write (gol,'("could not set grid mapping for grid type: ",a)') trim(gridtype); call goErr
          TRACEBACK; status=1; return
      end select
    else
      ! trap incorrect use of area-aver for corner variables
      select case ( trim(gridtype) )
        case ( 'corner', 'corner_bnds' )
          if ( trim(mapping) == 'area-aver' ) then
            write (gol,'("Do not use area-average mapping for variables with gridtype: ",a)') trim(gridtype) ; call goErr
            TRACEBACK;status=1;return
          end if
        case default
          ! ~ no action needed
      end select
      ! copy:
      gmapping = trim(mapping)
    end if
    
    ! setup for time t; return status -1 if new data should be put:
    call mf%Setup( t, status )
    if ( status > 0 ) then
      ! error status ...
      TRACEBACK; status=1; return
    else if ( status < 0 ) then

      ! prepare to receive new data:
      call mf%Setup_Prepare( tref, nreceive, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over fields to be received:
      do ireceive = 1, nreceive

        ! get input description:
        call mf%Get( status, input=inputs, name=name, units=units )
        IF_NOTOK_RETURN(status=1)
        ! get time value:
        call mf%Setup_InqTime( ireceive, receive_t, status )
        IF_NOTOK_RETURN(status=1)
        ! info ...
        call wrtgol( '    read `'//trim(name)//'` valid for ', receive_t ); call goPr
        ! extract input description valid for requested time:
        call GetInput( inputs, receive_t, tt, input, status )
        IF_NOTOK_RETURN(status=1)

        ! obtain data from somewhere, allocate mdata, and store into there;
        ! time interval 'tt_rec' for which field is valid is not used:
        call MeteoField_ReadRecord_sfc( tref, (/receive_t,receive_t/), rcF, input, units, &
                                         gri, gmapping, &
                                         mdata, tt_rec, status )
        IF_NOTOK_RETURN(status=1)
        ! store; optionally define the lower and upper bounds:
        call mf%Setup_Put( ireceive, mdata, receive_t, status )
        IF_NOTOK_RETURN(status=1)
        
      end do  ! fields to receive

      ! clear:
      deallocate( mdata, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! new data neeeded

    ! ok
    status = 0
    
  end subroutine LE_Data_Setup_iFS_sfc


  ! ***
  
  
  subroutine LE_Data_Setup_iFS_3D( mf, tref, t, rcF, gridtype, mapping, status, halflevels )
  
    use GO            , only : TDate, wrtgol
    use GO            , only : TrcFile
    use GO            , only : T_Instant_Field_Series
    use GO            , only : goReadFromLine
    use C3PO          , only : T_Grid_Ugg
    use LE_Grid       , only : ugg, ugg_bnds, ugg_crnr, ugg_crnr_bnds
    
    ! --- in/out ---------------------------------
    
    type(T_Instant_Field_Series), intent(inout)     ::  mf
    type(TDate), intent(in)                         ::  tref
    type(TDate), intent(in)                         ::  t
    type(TrcFile), intent(in)                       ::  rcF
    character(len=*), intent(in)                    ::  gridtype
    character(len=*), intent(in)                    ::  mapping
    integer, intent(out)                            ::  status
    logical, intent(in), optional                   ::  halflevels
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Setup_iFS_3D'
    
    ! --- local ----------------------------------
    
    type(T_Grid_Ugg), pointer   ::  gri
    character(len=64)           ::  gmapping
    character(len=1024)         ::  inputs, input
    character(len=1024)         ::  input_element
    type(TDate)                 ::  tt(2)
    type(TDate)                 ::  tt_rec(2)
    character(len=64)           ::  units
    integer                     ::  nreceive, ireceive
    type(TDate)                 ::  receive_t
    integer                     ::  nreset
    logical                     ::  reset
    real, allocatable           ::  mdata(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'corner' )
        gri => ugg_crnr
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case ( 'corner_bnds' )
        gri => ugg_crnr_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set grid mapping:
    if ( trim(mapping) == 'None' ) then
      ! reset default based on grid type:
      select case ( trim(gridtype) )
        case ( 'cell', 'cell_bnds' )
          gmapping = 'area-aver'
        case ( 'corner', 'corner_bnds' )
          gmapping = 'interpol'
        case default
          write (gol,'("could not set grid mapping for grid type: ",a)') trim(gridtype); call goErr
          TRACEBACK; status=1; return
      end select
    else
      ! trap incorrect use of area-aver for corner variables
      select case ( trim(gridtype) )
        case ( 'corner', 'corner_bnds' )
          if ( trim(mapping) == 'area-aver' ) then
            write (gol,'("Do not use area-average mapping for variables with gridtype: ",a)') trim(gridtype) ; call goErr
            TRACEBACK;status=1;return
          end if
        case default
          ! ~ no action needed
      end select

      ! copy:
      gmapping = trim(mapping)
    end if
    
    ! setup for time t; return status -1 if new data should be put:
    call mf%Setup( t, status )
    if ( status > 0 ) then
      TRACEBACK; status=1; return
    else if ( status < 0 ) then

      ! first try with complete reset ;
      ! this assumes that series are continous (without change in layers etc):
      nreset = 0
      ! loop unit ok:
      do

        ! increase counter:
        nreset = nreset + 1
        ! check ...
        if ( nreset > 2 ) then
          write (gol,'("still errors after reset, something wrong ...")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! prepare to receive new data; eventually total reset in case of discontinous series:
        call mf%Setup_Prepare( tref, nreceive, status, reset_levels=(nreset>1) )
        IF_NOTOK_RETURN(status=1)
        ! loop over fields to be received:
        do ireceive = 1, nreceive
          ! get time value:
          call mf%Setup_InqTime( ireceive, receive_t, status )
          IF_NOTOK_RETURN(status=1)
          ! info ...
          call wrtgol( '    read "'//trim(mf%name)//'" valid for ', receive_t ); call goPr
          ! get input description, per time interval a key for the further input description:
          !   __input_cams_nrt_VAR.input   :  2010-01-01 03:00, 2015-01-01 00:00, inputA1, inputA2 | \
          !                                   2015-01-01 03:00, 2020-01-01 00:00, inputB
          call mf%Get( status, input=inputs, units=units )
          IF_NOTOK_RETURN(status=1)
          ! extract input for requested time, e.g. "inputA1, inputA2", "inputB", ..
          ! set optional flag to ensure that interval [receive_t,..] 
          ! is selected instead of [...,receive_t] :
          call GetInput( inputs, receive_t, tt, input, status, &
                          isfirst=((ireceive==1) .and. (nreceive>1) ) )
          IF_NOTOK_RETURN(status=1)

          ! loop over input elements:
          !   inputA1, inputA2, ...
          do
            ! extract leading element:
            call goReadFromLine( input, input_element, status, sep=',' )
            IF_NOTOK_RETURN(status=1)
            ! obtain data from somewhere, allocate mdata, and store into there ;
            ! trap warning status (missing file?), then try next element:
            call MeteoField_ReadRecord_3d( tref, (/receive_t,receive_t/), rcF, trim(input_element), units, &
                                           gri, gmapping, &
                                           mdata, tt_rec, status, &
                                           halflevels=halflevels )
            if ( status < 0 ) then
              ! info ...
              write (gol,'("WARNING - could not read record, try next input description")'); call goPr
              ! check ..
              if ( len_trim(input) == 0 ) then
                write (gol,'("no new file input descriptions available anymore")'); call goErr
                TRACEBACK; status=1; return
              end if
              ! next:
              cycle
            end if
            IF_NOTOK_RETURN(status=1)
            ! store; use index space of allocated field if necessary;
            ! return status -1 if shape changed (number of layers):
            call mf%Setup_Put( ireceive, mdata, receive_t, status, &
                                 lbo=lbound(mdata), ubo=ubound(mdata) )
            if ( status == -1 ) then
              ! info ...
              write (gol,'("WARNING - probably number of layers changed, try again with full reset of temporal interpolations ...")'); call goPr
              ! request reset:
              reset = .true.
              ! leave loop:
              exit
            else if ( status == 0 ) then
              ! ok, no reset needed:
              reset = .false.
            else
              TRACEBACK; status=1; return
            endif
            ! data found, no new input element needed:
            exit
          end do ! input elements

          ! reset needed?
          if ( reset ) exit
        end do  ! fields to receive
        
        ! leave if no reset is needed:
        if ( .not. reset ) exit
      end do ! reset loop

      ! clear:
      deallocate( mdata, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! new data neeeded

    ! ok
    status = 0
    
  end subroutine LE_Data_Setup_iFS_3D
  
  
  ! ***
  
  
  subroutine LE_Data_Setup_iVFS_sfc( mf, tref, t, rcF, gridtype, mapping, status )
  
    use GO     , only : TDate, wrtgol
    use GO     , only : TrcFile
    use GO     , only : T_Instant_VectorField_Series
    use C3PO   , only : T_Grid_Ugg
    use LE_Grid, only : ugg, ugg_crnr, ugg_crnr_bnds

    ! --- in/out ---------------------------------
    
    type(T_Instant_VectorField_Series), intent(inout)       ::  mf
    type(TDate), intent(in)                                 ::  tref
    type(TDate), intent(in)                                 ::  t
    type(TrcFile), intent(in)                               ::  rcF
    character(len=*), intent(in)                            ::  gridtype
    character(len=*), intent(in)                            ::  mapping
    integer, intent(out)                                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Setup_iVFS_sfc'
    
    ! --- local ----------------------------------
    
    type(T_Grid_Ugg), pointer   ::  gri
    character(len=64)           ::  gmapping
    character(len=1024)         ::  inputs
    character(len=1024)         ::  input
    type(TDate)                 ::  tt(2)
    type(TDate)                 ::  tt_rec(2)
    character(len=64)           ::  name, units
    integer                     ::  nreceive, ireceive
    type(TDate)                 ::  receive_t
    real, allocatable           ::  mdata_u(:,:,:)
    real, allocatable           ::  mdata_v(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'corner' )
        gri => ugg_crnr
      case ( 'corner_bnds' )
        gri => ugg_crnr_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set grid mapping:
    if ( trim(mapping) == 'None' ) then
      ! reset default based on grid type:
      select case ( trim(gridtype) )
        case ( 'cell', 'cell_bnds' )
          gmapping = 'area-aver'
        case ( 'corner', 'corner_bnds' )
          gmapping = 'interpol'
        case default
          write (gol,'("could not set grid mapping for grid type: ",a)') trim(gridtype); call goErr
          TRACEBACK; status=1; return
      end select
    else
      ! trap incorrect use of area-aver for corner variables
      select case ( trim(gridtype) )
        case ( 'corner', 'corner_bnds' )
          if ( trim(mapping) == 'area-aver' ) then
            write (gol,'("Do not use area-average mapping for variables with gridtype: ",a)') trim(gridtype) ; call goErr
            TRACEBACK;status=1;return
          end if
        case default
          ! ~ no action needed
      end select

      ! copy:
      gmapping = trim(mapping)
    end if
    
    ! setup for time t; return status -1 if new data should be put:
    call mf%Setup( t, status )
    if ( status > 0 ) then
      ! error status ...
      TRACEBACK; status=1; return
    else if ( status < 0 ) then

      ! prepare to receive new data:
      call mf%Setup_Prepare( tref, nreceive, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over fields to be received:
      do ireceive = 1, nreceive

        ! get input description:
        call mf%Get( status, input=inputs, name=name, units=units )
        IF_NOTOK_RETURN(status=1)
        ! get time value:
        call mf%Setup_InqTime( ireceive, receive_t, status )
        IF_NOTOK_RETURN(status=1)

        ! info ...
        call wrtgol( '    read `'//trim(name)//'` (u,v) valid for ', receive_t ); call goPr
        ! extract input description valid for requested time:
        call GetInput( inputs, receive_t, tt, input, status )
        IF_NOTOK_RETURN(status=1)
        ! obtain data from somewhere, allocate mdata, and store into there;
        ! time interval 'tt_rec' for which field is valid is not used:
        call MeteoField_ReadRecordUV_sfc( tref, (/receive_t,receive_t/), rcF, input, units, &
                                           gri, gridtype, gmapping, &
                                           mdata_u, mdata_v, tt_rec, status )
        IF_NOTOK_RETURN(status=1)

        ! store; optionally define the lower and upper bounds:
        call mf%Setup_Put( ireceive, mdata_u, mdata_v, receive_t, status )
        IF_NOTOK_RETURN(status=1)
        
      end do  ! fields to receive

      ! clear:
      deallocate( mdata_u, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mdata_v, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! new data neeeded

    ! ok
    status = 0
    
  end subroutine LE_Data_Setup_iVFS_sfc


  ! ***
  
  
  subroutine LE_Data_Setup_iVFS_3D( mf, tref, t, rcF, gridtype, mapping, status, halflevels )
  
    use GO     , only : TDate, wrtgol
    use GO     , only : TrcFile
    use GO     , only : T_Instant_VectorField_Series
    use GO     , only : goReadFromLine
    use C3PO   , only : T_Grid_Ugg
    use LE_Grid, only : ugg, ugg_crnr, ugg_bnds, ugg_crnr_bnds, ugg_ustag, ugg_vstag

    ! --- in/out ---------------------------------
    
    type(T_Instant_VectorField_Series), intent(inout)       ::  mf
    type(TDate), intent(in)                                 ::  tref
    type(TDate), intent(in)                                 ::  t
    type(TrcFile), intent(in)                               ::  rcF
    character(len=*), intent(in)                            ::  gridtype
    character(len=*), intent(in)                            ::  mapping
    integer, intent(out)                                    ::  status
    logical, intent(in), optional                           ::  halflevels
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Setup_iVFS_3D'
    
    ! --- local ----------------------------------
    
    type(T_Grid_Ugg), pointer   ::  gri
    character(len=64)           ::  gmapping
    character(len=1024)         ::  inputs, input
    character(len=1024)         ::  input_element
    type(TDate)                 ::  tt(2)
    character(len=64)           ::  units
    integer                     ::  nreceive, ireceive
    type(TDate)                 ::  receive_t
    integer                     ::  nreset
    logical                     ::  reset
    real, allocatable           ::  mdata_u(:,:,:)
    real, allocatable           ::  mdata_v(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell', 'center' )
        gri => ugg
      case ( 'corner' )
        gri => ugg_crnr
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case ( 'corner_bnds' )
        gri => ugg_crnr_bnds
      case ( 'u-edge' )
        gri => ugg_ustag
      case ( 'v-edge' )
        gri => ugg_vstag
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! set grid mapping:
    if ( trim(mapping) == 'None' ) then
      ! reset default based on grid type:
      select case ( trim(gridtype) )
        case ( 'cell', 'cell_bnds' )
          gmapping = 'area-aver'
        case ( 'corner', 'corner_bnds', 'u-edge', 'v-edge', 'center' )
          gmapping = 'interpol'
        case default
          write (gol,'("could not set grid mapping for grid type: ",a)') trim(gridtype); call goErr
          TRACEBACK; status=1; return
      end select
    else
      ! trap incorrect use of area-aver for corner,edge and center variables
      select case ( trim(gridtype) )
        case ( 'corner', 'corner_bnds', 'u-edge', 'v-edge', 'center' )
          if ( trim(mapping) == 'area-aver' ) then
            write (gol,'("Do not use area-average mapping for variables with gridtype: ",a)') trim(gridtype) ; call goErr
            TRACEBACK;status=1;return
          end if
        case default
          ! ~ no action needed
      end select

      ! copy:
      gmapping = trim(mapping)
    end if
    
    ! setup for time t; return status -1 if new data should be put:
    call mf%Setup( t, status )
    if ( status > 0 ) then
      TRACEBACK; status=1; return
    else if ( status < 0 ) then

      ! first try with complete reset ;
      ! this assumes that series are continous (without change in layers etc):
      nreset = 0
      ! loop unit ok:
      do

        ! increase counter:
        nreset = nreset + 1
        ! check ...
        if ( nreset > 2 ) then
          write (gol,'("still errors after reset, something wrong ...")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! prepare to receive new data:
        call mf%Setup_Prepare( tref, nreceive, status, reset_levels=(nreset>1) )
        IF_NOTOK_RETURN(status=1)
        ! loop over fields to be received:
        do ireceive = 1, nreceive
          ! get time value:
          call mf%Setup_InqTime( ireceive, receive_t, status )
          IF_NOTOK_RETURN(status=1)
          ! info ...
          call wrtgol( '    read "'//trim(mf%name)//'" valid for ', receive_t ); call goPr
          ! get input description:
          call mf%Get( status, input=inputs, units=units )
          IF_NOTOK_RETURN(status=1)
          ! extract for requested time:
          call GetInput( inputs, receive_t, tt, input, status, &
                          isfirst=((ireceive==1) .and. (nreceive>1) ) )
          IF_NOTOK_RETURN(status=1)
          ! loop over input elements:
          do
            ! extract leading element:
            call goReadFromLine( input, input_element, status, sep=',' )
            IF_NOTOK_RETURN(status=1)
            ! obtain data from somewhere, allocate mdata, and store into there ;
            ! trap warning status (missing file?), then try next element:
            call MeteoField_ReadRecordUV_3d( tref, receive_t, rcF, trim(input_element), units, &
                                             gri, gmapping, &
                                             mdata_u, mdata_v, status, &
                                             halflevels=halflevels )
            if ( status < 0 ) then
              ! info ...
              write (gol,'("WARNING - could not read record, try next input description")'); call goPr
              ! check ..
              if ( len_trim(input) == 0 ) then
                write (gol,'("no new file input descriptions available anymore")'); call goErr
                TRACEBACK; status=1; return
              end if
              ! next:
              cycle
            end if
            IF_NOTOK_RETURN(status=1)
            ! store; use index space of allocated field if necessary:
            call mf%Setup_Put( ireceive, mdata_u, mdata_v, receive_t, status, &
                                 lbo=lbound(mdata_u), ubo=ubound(mdata_v) )
            if ( status == -1 ) then
              ! info ...
              write (gol,'("WARNING - probably number of layers changed, try again with full reset of temporal interpolations ...")'); call goPr
              ! request reset:
              reset = .true.
              ! leave loop:
              exit
            else if ( status == 0 ) then
              ! ok, no reset needed:
              reset = .false.
            else
              TRACEBACK; status=1; return
            endif
            IF_NOTOK_RETURN(status=1)
            ! leave:
            exit
          end do ! input elements
        end do  ! fields to receive
        ! leave if no reset is needed:
        if ( .not. reset ) exit
      end do ! reset loop

      ! clear:
      deallocate( mdata_u, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( mdata_v, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if ! new data neeeded

    ! ok
    status = 0
    
  end subroutine LE_Data_Setup_iVFS_3D



  ! ====================================================================
  ! ===
  ! === meteo files
  ! ===
  ! ====================================================================
  
  
  !
  ! Description of input series:
  !    t1, t2, description | t1, t2, description | ...
  !
  ! Return first description with t in [t1,t2],
  ! or second in case start time is the same as first end time:
  !    ...,  t, description1 |\
  !      t, .., description2
  !

  subroutine GetInput( line, t, tt, input, status, isfirst )

    use GO, only : TDate, wrtgol, operator(<=)
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)            ::  line
    type(TDate), intent(in)                 ::  t
    type(TDate), intent(out)                ::  tt(2)
    character(len=*), intent(out)           ::  input
    integer, intent(out)                    ::  status
    logical, intent(in), optional           ::  isfirst
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/GetInput'
    
    ! --- local ----------------------------------
    
    character(len=len(line))    ::  inputs
    logical                     ::  found
    logical                     ::  check_next_start
    type(TDate)                 ::  next_tt(2)
    character(len=len(input))   ::  next_input
    
    ! --- begin ----------------------------------
    
    ! set flag:
    check_next_start = .false.
    if ( present(isfirst) ) check_next_start = isfirst
    
    ! copy content:
    inputs = line
    ! init flag:
    found = .false.
    ! loop over (t1,t2,description) parts:
    do
      ! empty ?
      if ( len_trim(inputs) == 0 ) exit
      ! up to next seperation:
      call goReadFromLine( inputs, input, status, sep='|' )
      IF_NOTOK_RETURN(status=1)
      ! extract time range, remainder is input description:
      call goReadFromLine( input, tt(1), status, sep=',' )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( input, tt(2), status, sep=',' )
      IF_NOTOK_RETURN(status=1)
      ! valid for requested time?
      if ( (tt(1) <= t) .and. (t <= tt(2)) ) then
        found = .true.
        exit
      end if
    end do
    
    ! check next start?
    if ( found .and. (len_trim(inputs) > 0) .and. check_next_start ) then
      ! up to next seperation:
      call goReadFromLine( inputs, next_input, status, sep='|' )
      IF_NOTOK_RETURN(status=1)
      ! extract time range, remainder is input description:
      call goReadFromLine( next_input, next_tt(1), status, sep=',' )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( next_input, next_tt(2), status, sep=',' )
      IF_NOTOK_RETURN(status=1)
      ! valid for requested time?
      if ( (next_tt(1) <= t) .and. (t <= next_tt(2)) ) then
        ! reset output:
        input = trim(next_input)
        tt(1) = next_tt(1)
        tt(2) = next_tt(2)
        ! set flag:
        found = .true.
      end if
    end if  ! try next
    
    ! check ...
    if ( .not. found ) then
      call wrtgol( 'no description found for time: ', t ); call goErr
      write (gol,'("input line:")'); call goErr
      write (gol,'(a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine GetInput
  
  !
  ! Form filename from input line, expand time placeholders:
  !

  subroutine GetFilename( rcF, rcbase, tref, t, filename, status )

    use GO, only : TrcFile, ReadRc
    use GO, only : TDate, NewDate, IncrDate, IsAnyDate, rTotal
    use GO, only : operator(>=), operator(<=), operator(+), operator(-)
    use GO, only : goReplace
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)               ::  rcF
    character(len=*), intent(in)            ::  rcbase
    type(TDate), intent(in)                 ::  tref
    type(TDate), intent(in)                 ::  t
    character(len=*), intent(out)           ::  filename
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/GetFilename'
    
    ! --- local ----------------------------------

    integer                ::  dtsec
    character(len=32)      ::  fc00ext, fc12ext, fc72ext, fcext
    character(len=32)      ::  fc00, fc12
    type(TDate)            ::  tf
    real                   ::  rday
    integer                ::  fcday
    type(TDate)            ::  t12, t72
    type(TDate)            ::  tref00
    
    ! --- begin ----------------------------------
    
    ! extract filename template:
    call ReadRc( rcF, trim(rcbase)//'.name', filename, status )
    IF_NOTOK_RETURN(status=1)
    
    ! optional offset for time values in filename template:
    call ReadRc( rcF, trim(rcbase)//'.name_dtsec', dtsec, status, default=0 )
    IF_ERROR_RETURN(status=1)
    
    ! forecast extenstion:
    call ReadRc( rcF, trim(rcbase)//'.name_fcext', fcext, status, default='' )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rcbase)//'.name_fc00ext', fc00ext, status, default='' )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rcbase)//'.name_fc12ext', fc12ext, status, default='' )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rcbase)//'.name_fc72ext', fc72ext, status, default='' )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rcbase)//'.name_fc00', fc00, status, default='' )
    IF_ERROR_RETURN(status=1)
    call ReadRc( rcF, trim(rcbase)//'.name_fc12', fc12, status, default='' )
    IF_ERROR_RETURN(status=1)
    
    ! reference time of forecast?
    if ( .not. IsAnyDate(tref) ) then

      ! form time value in filename:
      tf = t + IncrDate(sec=dtsec)
      ! after forecast base ? then use base:
      if ( tf >= tref ) tf = tref
    
      ! replace placeholders for time values:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tf%year , status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tf%month, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tf%day  , status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{hh}'  , '(i2.2)', tf%hour , status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{yyyymmdd}'  , '(i8.8)'  , tf%year*10000  +tf%month*100  +tf%day            , status )
      IF_NOTOK_RETURN(status=1)
      
      ! after forecast base ?
      if ( tf >= tref ) then
        !~ before or after 12:00 ?
        t12 = NewDate( time6=(/tf%year,tf%month,tf%day,12,00,00/) )
        t72 = NewDate( time6=(/tf%year,tf%month,tf%day,72,00,00/) )
        if ( t <= t12 ) then
          ! insert extension for forecast files if present:
          call goReplace( filename, '%{fcext}', trim(fc00ext), status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( filename, '%{fchhext}', trim(fc00ext), status )
          IF_NOTOK_RETURN(status=1)
        else if ( t <= t72 ) then
          ! insert extension for forecast files if present:
          call goReplace( filename, '%{fchhext}', trim(fc12ext), status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( filename, '%{fcext}', trim(fc12ext), status )
          IF_NOTOK_RETURN(status=1)
        else
          ! insert extension for forecast files if present:
          call goReplace( filename, '%{fchhext}', trim(fc72ext), status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( filename, '%{fcext}', trim(fc72ext), status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! remove extension if present:
        call goReplace( filename, '%{fchh}', '', status )
        IF_NOTOK_RETURN(status=1)
      else
        ! remove extension if present:
        call goReplace( filename, '%{fchhext}', '', status )
        IF_NOTOK_RETURN(status=1)
        call goReplace( filename, '%{fcext}', '', status )
        IF_NOTOK_RETURN(status=1)
        ! on 'analysis' day eventually different files
        ! for (00,12] and (12,24] ; insert extensions if present:
        !~ before or after 12:00 ?
        t12 = NewDate( time6=(/tf%year,tf%month,tf%day,12,00,00/) )
        if ( t <= t12 ) then
          ! insert extension for forecast files if present:
          call goReplace( filename, '%{fchh}', trim(fc00), status )
          IF_NOTOK_RETURN(status=1)
        else
          ! insert extension for forecast files if present:
          call goReplace( filename, '%{fchh}', trim(fc12), status )
          IF_NOTOK_RETURN(status=1)
        end if
      end if

      ! Special: %{Dfcday}
      !   Assign labels given 'tf' (>= tref, see above):
      !     tref00+[00,24)  : D1
      !     tref00+[24,48)  : D2
      !     :
      !     tref00+[96,120) : D5
      !   Use the key '%{Dfcday5}' for files D1..D5 with in D5
      !   an extra record at the end for 24:00.
      if ( index( filename, '%{Dfcday' ) > 0 ) then
        ! start of day holding tref:
        tref00 = NewDate( time6=(/tref%year,tref%month,tref%day,0,0,0/) )
        ! days since start of forecast as real number:
        rday = rTotal( t+IncrDate(sec=dtsec) - tref00, 'day' )
        ! round to floor, add one to have 1,2,3,..
        fcday = int(floor(rday)) + 1
        ! minimum value of one in case previous days (spinup) is included in the run;
        ! for these days, use D1 files:
        fcday = max(fcday,1)
        ! replace by 'D1',... etc:
        call goReplace( filename, '%{Dfcday}', '("D",i1)', fcday, status )
        IF_NOTOK_RETURN(status=1)
        ! replace by 'D1',..,'D5':
        call goReplace( filename, '%{Dfcday5}', '("D",i1)', min(fcday,5), status )
        IF_NOTOK_RETURN(status=1)
      end if ! special key?
      
    else
    
      ! form time value:
      tf = t + IncrDate(sec=dtsec)

      ! replace placeholders for time values:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tf%year , status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tf%month, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tf%day  , status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{hh}'  , '(i2.2)', tf%hour , status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{yyyymmdd}'  , '(i8.8)'  , tf%year*10000  +tf%month*100  +tf%day            , status )
      IF_NOTOK_RETURN(status=1)

      ! replace by empty string:
      call goReplace( filename, '%{fcext}', '', status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filename, '%{fchhext}', '', status )
      IF_NOTOK_RETURN(status=1)
      ! on 'analysis' day eventually different files
      ! for (00,12] and (12,24] ; insert extensions if present:
      !~ before or after 12:00 ?
      t12 = NewDate( time6=(/tf%year,tf%month,tf%day,12,00,00/) )
      if ( t <= t12 ) then
        ! insert extension for forecast files if present:
        call goReplace( filename, '%{fchh}', trim(fc00), status )
        IF_NOTOK_RETURN(status=1)
      else
        ! insert extension for forecast files if present:
        call goReplace( filename, '%{fchh}', trim(fc12), status )
        IF_NOTOK_RETURN(status=1)
      end if

    end if ! tref present?
    
    ! ok
    status = 0
    
  end subroutine GetFilename
  

  ! ***
  
  
  subroutine MeteoField_ReadRecord_sfc( tref, tt, rcF, input, units, &
                                          ugg, mapping, &
                                          values, tt_rec, status )
  
    use GO            , only : TDate
    use GO            , only : operator(+), operator(-), operator(*)
    use Go            , only : TrcFile, ReadRc
    use C3PO          , only : T_File_Nc
    use C3PO          , only : T_Grid_Ugg

    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)            ::  tref
    type(TDate), intent(in)            ::  tt(2)
    type(TrcFile), intent(in)          ::  rcF
    character(len=*), intent(in)       ::  input
    character(len=*), intent(in)       ::  units
    type(T_Grid_Ugg), intent(in)       ::  ugg
    character(len=*), intent(in)       ::  mapping
    real, allocatable, intent(out)     ::  values(:,:,:)  ! (nlon,nlat,1)
    type(TDate), intent(out)           ::  tt_rec(2)
    integer, intent(out)               ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MeteoField_ReadRecord_sfc'
    
    ! --- local ----------------------------------
    
    !type(TDate)           ::  tmid
    character(len=512)    ::  filename
    character(len=512)    ::  description
    character(len=512)    ::  varkeys
    character(len=512)    ::  climat
    integer               ::  varid
    integer               ::  shp(3)
    type(T_File_Nc)       ::  file_in
    
    ! --- begin ----------------------------------

    !! mid time:
    !tmid = tt(1) + (tt(2)-tt(1))*0.5
    !! form filename using mid of interval:
    !call GetFilename( rcF, input, tref, tmid, filename, status )
    !IF_NOTOK_RETURN(status=1)
    
    ! form filename using end of interval:
    call GetFilename( rcF, input, tref, tt(2), filename, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read variable description line:
    call ReadRc( rcF, trim(input)//'.var', description, status )
    IF_NOTOK_RETURN(status=1)
    ! special instructions, for example for 3D pressure:
    call ReadRc( rcF, trim(input)//'.special', varkeys, status, default='' )
    IF_ERROR_RETURN(status=1)
    ! climatological time elements:
    call ReadRc( rcF, trim(input)//'.climat', climat, status, default='' )
    IF_ERROR_RETURN(status=1)

    ! open file:
    call file_in%Open( trim(filename), status )
    IF_NOTOK_RETURN(status=1)
    ! variable id:
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! close:
    call file_in%Close( status )
    IF_NOTOK_RETURN(status=1)

    ! read record:    
    call MeteoField_ReadRecord_sfc_ugg( tt, filename, trim(description), units, &
                                          ugg, mapping, varkeys, climat, &
                                          values, tt_rec, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine MeteoField_ReadRecord_sfc


  ! *** 

  
  !
  ! Read surface field on universial grid.
  !
  ! Return values:
  !   values(:,:,:)      : 3D values
  !   tt_rec(2)          : time range for which values are valid
  !
  ! Special keys:
  !    constant=T        :  read first value from time series, store as if constant
  !    time_name=XTIME   :  time is stored in  variable with name different from dimension
  !    start_time_attr=SIMULATION_START_DATE :
  !                         reset the field to zero if (equal) time values match attribute value;
  !                         this is used to reset accumulated fields to zero
  !
  
  subroutine MeteoField_ReadRecord_sfc_ugg( tt, filename, description, units, &
                                             ugg, mapping, keys, climat, &
                                             values, tt_rec, status )
  
    use Binas         , only : grav
    use GO            , only : goVarValue
    use GO            , only : TDate, AnyDate, operator(==), wrtgol
    use C3PO          , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : LE_Data_ClimatMask
    use LE_Data_Common, only : Grid_Convertors
    use LE_Data_Common, only : LE_Data_ConversionFactor

    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)            ::  tt(2)
    character(len=*), intent(in)       ::  filename
    character(len=*), intent(in)       ::  description
    character(len=*), intent(in)       ::  units
    type(T_Grid_Ugg), intent(in)       ::  ugg
    character(len=*), intent(in)       ::  mapping
    character(len=*), intent(in)       ::  keys
    character(len=*), intent(in)       ::  climat
    real, allocatable, intent(out)     ::  values(:,:,:)  ! (nlon,nlat,1)
    type(TDate), intent(out)           ::  tt_rec(2)
    integer, intent(out)               ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MeteoField_ReadRecord_sfc_ugg'
    
    ! --- local ----------------------------------
    
    logical               ::  constant
    character(len=64)     ::  start_time_attr
    type(TDate)           ::  t_attr
    integer               ::  varid
    character(len=64)     ::  time_dim_name
    logical               ::  cmask(6)
    integer               ::  irec
    character(len=64)     ::  units_in
    real                  ::  factor
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    integer               ::  start_ind(3)
    integer               ::  count_ind(3)
    real, allocatable     ::  values_in(:,:)   ! (lon,lat)

    ! --- begin ----------------------------------

    ! info ...
    write (gol,'("      open file: ",a)') trim(filename); call goPr
    if ( len_trim(keys) > 0 ) then
      write (gol,'("        specials: ",a)') trim(keys); call goPr
    end if

    ! init climat mask:
    call LE_Data_ClimatMask( climat, cmask, status )
    IF_NOTOK_RETURN(status=1)

    ! specials:
    constant = .false.
    start_time_attr = ''
    ! reset from keys if present:
    call goVarValue( keys, ';', 'constant', '=', constant, status )
    IF_ERROR_RETURN(status=1)
    call goVarValue( keys, ';', 'start_time_attr', '=', start_time_attr, status )
    IF_ERROR_RETURN(status=1)
    
    ! info ...
    if ( len_trim(start_time_attr) > 0 ) then
      write (gol,'("        start_time_attr: ",a)') trim(start_time_attr); call goPr
    end if
    
    ! open file:
    call file_in%Open( trim(filename), status )
    IF_NOTOK_RETURN(status=1)

    ! variable id:
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! get name of time dimension (or 'None' if not present):
    call file_in%Inq_Var( varid, status, time_dim_name=time_dim_name )
    IF_NOTOK_RETURN(status=1)

    ! init grid definition from file:
    call file_in%Get_Grid( varid, grid_in, status, &
                            ugg_to=ugg, start_ind=start_ind(1:2), count_ind=count_ind(1:2) )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( values_in(grid_in%nlon,grid_in%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! read data, and set time range for which data is valid;
    ! contant or time dependend ?
    if ( trim(time_dim_name) == 'None' ) then
    
      !! info ...
      !write (gol,'("      read constant field from file ...")'); call goPr
      
      ! read data:
      call file_in%Get_Var( trim(description), values_in, units_in, status, &
                                start=start_ind(1:2), count=count_ind(1:2) )
      IF_NOTOK_RETURN(status=1)
      
      ! valid for a long time:
      tt_rec = (/ AnyDate(), AnyDate() /)
      
    else
    
      ! index of  time record, and timerange for which record is valid
      !~ constant data, just read first record:
      if ( constant ) then
        !! info ...
        !write (gol,'("      read first time record (constant field)")'); call goPr
        ! set record index to first:
        irec = 1
        ! valid for a long time:
        tt_rec = (/ AnyDate(), AnyDate() /)
      !
      !~ instant time
      else if ( tt(1) == tt(2) ) then
        !! info ...
        !call wrtgol( '      read time record valid for ', tt(1) ); call goPr
        ! instant time, record time should match exactly:
        call file_in%Inq_TimeRecord( trim(time_dim_name), tt(1), irec, status, &
                                        keys=keys, climat=cmask )
        IF_NOTOK_RETURN(status=1)
        ! fill 'interval' for which input record is valid with single time:
        tt_rec = (/ tt(1), tt(2) /)
      !
      !~ time interval, start time not defined by time bounds but attribute:
      else if ( len_trim(start_time_attr) > 0 ) then
        !! info ...
        !call wrtgol( '      read time record ending at ', tt(2) ); call goPr
        ! instant time, record time should match exactly:
        call file_in%Inq_TimeRecord( trim(time_dim_name), tt(2), irec, status, &
                                       keys=keys, climat=cmask )
        IF_NOTOK_RETURN(status=1)
        !! info ...
        !write (gol,'("      check attribute `",a,"` for start time ...")') trim(start_time_attr); call goPr
        ! read time value from global attribute:
        call file_in%Get_Var_Attr( '', trim(start_time_attr), t_attr, status )
        IF_NOTOK_RETURN(status=1)
        !! info ...
        !call wrtgol( '        attribute value : ', t_attr ); call goPr
        ! reset start time for record:
        tt_rec = (/ t_attr, tt(2) /)
      !
      !~ time interval, time bounds are available:
      else
        !! info ...
        !call wrtgol( '      read time record enclosing ', tt(1), ' - ', tt(2) ); call goPr
        ! time interval; time variable might have associated 'bounds',
        ! of which one record should enclose the interval tt;
        ! this record interval is returned as 'tt_rec':
        call file_in%Inq_TimeRecord( trim(time_dim_name), tt, irec, tt_rec, status, &
                                        keys=keys, climat=cmask )
        IF_NOTOK_RETURN(status=1)
      end if
      ! fill into into space:
      start_ind(3) = irec
      count_ind(3) = 1
    
      ! read data slab:
      call file_in%Get_Var( trim(description), values_in, units_in, status, &
                                start=start_ind, count=count_ind )
      IF_NOTOK_RETURN(status=1)
      
    end if
    
    ! close:
    call file_in%Close( status )
    IF_NOTOK_RETURN(status=1)

    !! info ...
    !write (gol,*) '       value range : ', minval(values_in), maxval(values_in); call goPr
    !call wrtgol( '        valid for ', tt_rec ); call goPr
    
    ! get conversion factor:
    call LE_Data_ConversionFactor( trim(units_in), trim(units), factor, status )
    IF_NOTOK_RETURN(status=1)
    ! apply ?
    if ( factor /= 1.0 ) then
      values_in = values_in * factor
    end if ! conversion needed ?
    
    ! storage:
    allocate( values(ugg%nlon,ugg%nlat,1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! regrid:
    select case ( mapping )
      !~ interpolation
      case ( 'interpol' )
        ! interpolate to grid points:
        call Grid_Convertors%Ugg_Interpol( grid_in, values_in, ugg, values(:,:,1), status )
        IF_NOTOK_RETURN(status=1)
      !~ nearest neighbour:
      case ( 'nearest' )
        ! interpolate to grid points:
        call Grid_Convertors%Ugg_Nearest( grid_in, values_in, ugg, values(:,:,1), status )
        IF_NOTOK_RETURN(status=1)
      !~ area average:
      case ( 'area-aver' )
        ! grid cell averages:
        call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, values(:,:,1), status )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unsupported mapping `",a,"`")') trim(mapping); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! clear:
    deallocate( values_in, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    call grid_in%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine MeteoField_ReadRecord_sfc_ugg
  

  ! ***
  
  
  subroutine MeteoField_ReadRecordUV_sfc( tref, tt, rcF, input, units, &
                                          ugg, gridtype, mapping, &
                                          values_u, values_v, tt_rec, status )
  
    use GO            , only : TDate
    use GO            , only : operator(+), operator(-), operator(*)
    use Go            , only : TrcFile, ReadRc
    use C3PO          , only : T_File_Ugg, T_Grid_Ugg

    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)            ::  tref
    type(TDate), intent(in)            ::  tt(2)
    type(TrcFile), intent(in)          ::  rcF
    character(len=*), intent(in)       ::  input
    character(len=*), intent(in)       ::  units
    type(T_Grid_Ugg), intent(in)       ::  ugg
    character(len=*), intent(in)       ::  gridtype
    character(len=*), intent(in)       ::  mapping
    real, allocatable, intent(out)     ::  values_u(:,:,:)  ! (nlon,nlat,1)
    real, allocatable, intent(out)     ::  values_v(:,:,:)  ! (nlon,nlat,1)
    type(TDate), intent(out)           ::  tt_rec(2)
    integer, intent(out)               ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MeteoField_ReadRecordUV_sfc'
    
    ! --- local ----------------------------------
    
    type(TDate)           ::  tmid
    integer               ::  k
    character(len=1)      ::  keys(2)
    character(len=512)    ::  filenames(2)
    character(len=512)    ::  descriptions(2)
    character(len=512)    ::  varkeys
    character(len=512)    ::  climat
    integer               ::  varid
    integer               ::  shp(3)
    type(T_File_Ugg)      ::  file_in
    
    ! --- begin ----------------------------------

    ! mid time:
    tmid = tt(1) + (tt(2)-tt(1))*0.5
    
    ! uv keys:
    keys = (/ 'u', 'v' /)
    ! loop over uv:
    do k = 1, 2
      ! form filename using mid of interval:
      call GetFilename( rcF, trim(input)//'.'//keys(k), tref, tmid, filenames(k), status )
      IF_NOTOK_RETURN(status=1)
      ! read variable description line:
      call ReadRc( rcF, trim(input)//'.'//keys(k)//'.var', descriptions(k), status )
      IF_NOTOK_RETURN(status=1)
    end do ! k
    
    ! special instructions, for example for 3D pressure:
    call ReadRc( rcF, trim(input)//'.special', varkeys, status, default='' )
    IF_ERROR_RETURN(status=1)
    ! climatological time elements:
    call ReadRc( rcF, trim(input)//'.climat', climat, status, default='' )
    IF_ERROR_RETURN(status=1)
    
    ! read record:
    call MeteoField_ReadRecordUV_sfc_ugg( tt, filenames, descriptions, units, &
                                          ugg, mapping, varkeys, climat, &
                                          values_u, values_v, tt_rec, status )
    IF_NOTOK_RETURN(status=1)
        
    ! ok
    status = 0
    
  end subroutine MeteoField_ReadRecordUV_sfc
  
  
  ! *
  
  
  subroutine MeteoField_ReadRecordUV_sfc_ugg( tt, filenames, descriptions, &
                                              units, ugg, mapping, keys, climat, &
                                              values_u, values_v, tt_rec, status )
  
    use Binas         , only : grav
    use GO            , only : TDate, AnyDate, operator(==), wrtgol
    use GO            , only : goVarValue
    use C3PO          , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : LE_Data_ClimatMask
    use LE_Data_Common, only : Grid_Convertors
    use LE_Data_Common, only : LE_Data_ConversionFactor

    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)            ::  tt(2)
    character(len=*), intent(in)       ::  filenames(2)
    character(len=*), intent(in)       ::  descriptions(2)
    character(len=*), intent(in)       ::  units
    type(T_Grid_Ugg), intent(in)       ::  ugg
    character(len=*), intent(in)       ::  mapping
    character(len=*), intent(in)       ::  keys
    character(len=*), intent(in)       ::  climat
    real, allocatable, intent(out)     ::  values_u(:,:,:)  ! (nlon,nlat,1)
    real, allocatable, intent(out)     ::  values_v(:,:,:)  ! (nlon,nlat,1)
    type(TDate), intent(out)           ::  tt_rec(2)
    integer, intent(out)               ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MeteoField_ReadRecordUV_sfc_ugg'
    
    ! --- local ----------------------------------
    
    character(len=64)     ::  compute
    integer               ::  varid
    character(len=64)     ::  time_dim_name
    logical               ::  cmask(6)
    integer               ::  irec
    character(len=64)     ::  units_in
    real                  ::  factor
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  ugg_in
    integer               ::  start_ind(3)
    integer               ::  count_ind(3)
    real, allocatable     ::  values_in(:,:,:)   ! (lon,lat,2)
    integer               ::  k

    ! --- begin ----------------------------------
    
    ! compute task, default none:
    compute = 'None'
    call goVarValue( keys, ';', 'compute', '=', compute, status )

    ! init climat mask:
    call LE_Data_ClimatMask( climat, cmask, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over u,v
    do k = 1, 2

      ! open file:
      call file_in%Open( trim(filenames(k)), status )
      IF_NOTOK_RETURN(status=1)

      ! variable id:
      call file_in%Inq_VarID( trim(descriptions(k)), varid, status )
      IF_NOTOK_RETURN(status=1)
      ! get name of time dimension (or 'None' if not present):
      call file_in%Inq_Var( varid, status, time_dim_name=time_dim_name )
      IF_NOTOK_RETURN(status=1)

      ! first? then setup grid and storage:
      if ( k == 1 ) then
        ! init grid definition from file:
        call file_in%Get_Grid( varid, ugg_in, status, &
                                ugg_to=ugg, start_ind=start_ind(1:2), count_ind=count_ind(1:2) )
        IF_NOTOK_RETURN(status=1)
        ! storage:
        allocate( values_in(ugg_in%nlon,ugg_in%nlat,2), stat=status )
        IF_NOTOK_RETURN(status=1)
      end if ! first of u,v

      ! index of  time record:
      if ( tt(1) == tt(2) ) then
        !! info ...
        !call wrtgol( '      read time record valid for ', tt(1) ); call goPr
        ! instant time, record tim should match exactly:
        call file_in%Inq_TimeRecord( trim(time_dim_name), tt(1), irec, status, &
                                           keys=keys, climat=cmask )
        IF_NOTOK_RETURN(status=1)
        ! fill 'interval' for which input record is valid, here singel time:
        tt_rec = (/ tt(1), tt(2) /)
      else
        !! info ...
        !call wrtgol( '      read time record enclosing ', tt(1), ' - ', tt(2) ); call goPr
        ! time interval; time variable should have associated 'bounds'
        ! of which one record should enclose the interval:
        call file_in%Inq_TimeRecord( trim(time_dim_name), tt, irec, tt_rec, status, &
                                           keys=keys )
        IF_NOTOK_RETURN(status=1)
      end if      
      ! fill into into space:
      start_ind(3) = irec
      count_ind(3) = 1

      ! read data slab:
      call file_in%Get_Var( trim(descriptions(k)), values_in(:,:,k), units_in, status, &
                                start=start_ind, count=count_ind )
      IF_NOTOK_RETURN(status=1)

      ! get conversion factor:
      call LE_Data_ConversionFactor( trim(units_in), trim(units), factor, status )
      IF_NOTOK_RETURN(status=1)
      ! apply ?
      if ( factor /= 1.0 ) then
        values_in(:,:,k) = values_in(:,:,k) * factor
      end if ! conversion needed ?

      ! close:
      call file_in%Close( status )
      IF_NOTOK_RETURN(status=1)

    end do ! u,v

    ! storage:
    allocate( values_u(ugg%nlon,ugg%nlat,1), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( values_v(ugg%nlon,ugg%nlat,1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! regrid:
    select case ( mapping )
      !~ interpolation
      case ( 'interpol' )
        ! interpolate to grid points:
        call Grid_Convertors%Ugg_Interpol( ugg_in, values_in(:,:,1), values_in(:,:,2), &
                                              ugg, values_u(:,:,1), values_v(:,:,1), status )
        IF_NOTOK_RETURN(status=1)
      
      !~ nearest neighbour:
      !case ( 'nearest' )
      !  ! interpolate to grid points:
      !  call Grid_Convertors%Ugg_Nearest( ugg_in, values_in, ugg, values(:,:,1), status )
      !  IF_NOTOK_RETURN(status=1)
      !~ area average:
      case ( 'area-aver' )
        ! grid cell averages for vector field:
        call Grid_Convertors%Ugg_AreaAver( ugg_in, values_in(:,:,1), values_in(:,:,2), &
                                              ugg, values_u(:,:,1), values_v(:,:,1), status )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unsupported mapping `",a,"`")') trim(mapping); call goErr
        TRACEBACK; status=1; return
    end select

    ! clear:
    deallocate( values_in, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    call ugg_in%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine MeteoField_ReadRecordUV_sfc_ugg
  
  
  ! ***
  
  
  subroutine MeteoField_ReadRecord_3d( tref, tt, rcF, input, units, &
                                        ugg, mapping, &
                                        values, tt_rec, status, &
                                        halflevels )
  
    use GO            , only : TDate, IncrDate, wrtgol  !, Set
    use GO            , only : operator(+), operator(-), operator(*)
    use GO            , only : TrcFile, ReadRc
    use C3PO          , only : T_File_nc, T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)            ::  tref
    type(TDate), intent(in)            ::  tt(2)
    type(TrcFile), intent(in)          ::  rcF
    character(len=*), intent(in)       ::  input
    character(len=*), intent(in)       ::  units
    type(T_Grid_Ugg), intent(in)       ::  ugg
    character(len=*), intent(in)       ::  mapping
    real, allocatable, intent(out)     ::  values(:,:,:)
    type(TDate), intent(out)           ::  tt_rec(2)
    integer, intent(out)               ::  status
    logical, intent(in), optional      ::  halflevels
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MeteoField_ReadRecord_3d'
    
    ! --- local ----------------------------------
    
    type(TDate)           ::  tmid
    logical               ::  on_halflevels
    character(len=512)    ::  filename
    character(len=512)    ::  description
    character(len=512)    ::  varkeys
    character(len=512)    ::  climat
    logical               ::  exist
    logical               ::  fixed_leap
    type(TDate)           ::  ttx(2)
    
    ! --- begin ----------------------------------
    
    ! halflevel field ?
    on_halflevels = .false.
    if ( present(halflevels) ) on_halflevels = halflevels
    
    !! info ..
    !write (gol,'(a,": open file ",a," ...")') rname, trim(filename); call goPr

    ! read variable description line:
    call ReadRc( rcF, trim(input)//'.var', description, status )
    IF_NOTOK_RETURN(status=1)
    ! special instructions, for example for 3D pressure:
    call ReadRc( rcF, trim(input)//'.special', varkeys, status, default='' )
    IF_ERROR_RETURN(status=1)
    ! climatological time elements:
    call ReadRc( rcF, trim(input)//'.climat', climat, status, default='' )
    IF_ERROR_RETURN(status=1)

    ! mid time:
    tmid = tt(1) + (tt(2)-tt(1))*0.5
    ! copy:
    ttx = tt
    
    ! form filename:
    call GetFilename( rcF, input, tref, tmid, filename, status )
    IF_NOTOK_RETURN(status=1)
    
    ! init flag ..
    fixed_leap = .false.
    ! adhoc fix: climatology might be missing leap day ...
    if ( len_trim(climat) > 0 ) then
      ! check if file exists:
      inquire( file=filename, exist=exist, iostat=status )
      IF_NOTOK_RETURN(status=1)
      ! Feb 29 requested? then use Feb 28 instead ..
      ! Mar 01 requested and file not present (for Feb 29) ? then also use Feb 28 instead ..
      if ( ( (tmid%month == 2) .and. (tmid%day == 29)                     ) .or. &
           ( (tmid%month == 3) .and. (tmid%day ==  1) .and. (.not. exist) )      ) then
        ! info ...
        write (gol,'("WARNING - climatology request, but file for leap-day not found:")'); call goPr
        write (gol,'("WARNING -   ",a)') trim(filename); call goPr
        write (gol,'("WARNING - use previous day instead ...",a)') trim(filename); call goPr
        ! reset:
        if ( tmid%month == 2 ) then
          ttx(1) = ttx(1) - IncrDate(day=1)
          ttx(2) = ttx(2) - IncrDate(day=1)
        else
          ttx(1) = ttx(1) - IncrDate(day=2)
          ttx(2) = ttx(2) - IncrDate(day=2)
        end if
        tmid = ttx(1) + (ttx(2)-ttx(1))*0.5
        ! create new filename:
        call GetFilename( rcF, input, tref, tmid, filename, status )
        IF_NOTOK_RETURN(status=1)
        ! set flag ..
        fixed_leap = .true.
      end if ! file for leap-day not found
    end if   ! climat series
    
    ! read record:
    call MeteoField_ReadRecord_3d_ugg( ttx, filename, trim(description), units, &
                                         ugg, mapping, on_halflevels, &
                                         varkeys, climat, &
                                         values, tt_rec, status )
    IF_NOTOK_RETURN(status=1)
    
    ! testing ..
    if ( fixed_leap ) then
      call wrtgol( 'WARNING - requested timerange: ', ttx ); call goPr
      write (gol,'("WARNING -   source file  : ",a)') trim(filename); call goPr
      call wrtgol( 'WARNING -   source record: ', tt_rec ); call goPr
    end if

    ! ok
    status = 0

  end subroutine MeteoField_ReadRecord_3d
  
  
  ! *
  
  
  subroutine MeteoField_ReadRecord_3d_ugg( tt, filename, description, units, &
                                           ugg, mapping, on_halflevels, keys, climat, &
                                           values, tt_rec, status )
  
    use GO            , only : TDate, AnyDate, operator(==), wrtgol
    use GO            , only : goVarValue, goReadFromLine
    use C3PO          , only : T_File_Ugg
    use C3PO          , only : T_Grid_Ugg
    use C3PO          , only : T_Levs_Hyb
    use LE_Data_Common, only : LE_Data_ClimatMask
    use LE_Data_Common, only : Grid_Convertors
    use LE_Data_Common, only : LE_Data_ConversionFactor
  
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)            ::  tt(2)
    character(len=*), intent(in)       ::  filename
    character(len=*), intent(in)       ::  description
    character(len=*), intent(in)       ::  units
    type(T_Grid_Ugg), intent(in)       ::  ugg
    character(len=*), intent(in)       ::  mapping
    logical, intent(in)                ::  on_halflevels
    character(len=*), intent(in)       ::  keys
    character(len=*), intent(in)       ::  climat
    real, allocatable, intent(out)     ::  values(:,:,:)
    type(TDate), intent(out)           ::  tt_rec(2)
    integer, intent(out)               ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MeteoField_ReadRecord_3d_ugg'
    
    ! --- local ----------------------------------
    
    character(len=64)     ::  compute
    logical               ::  constant
    logical               ::  swap_levels
    integer               ::  varid
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    integer               ::  start_ind(4)
    integer               ::  count_ind(4)
    integer               ::  nlev
    type(T_Levs_Hyb)      ::  levs
    character(len=64)     ::  time_dim_name
    logical               ::  cmask(6)
    integer               ::  irec
    real, allocatable     ::  values_in(:,:,:)    ! (i,j,lev)
    real, allocatable     ::  ps_in(:,:)          ! (i,j)
    integer               ::  ndim
    integer               ::  lbo(3), ubo(3)
    real, allocatable     ::  values_hor(:,:,:)    ! (lon,lat,lev)
    character(len=64)     ::  units_in
    real                  ::  factor
    integer               ::  k
    real                  ::  newval
    
    ! --- begin ----------------------------------
    
    ! compute task, default none:
    compute = 'None'
    call goVarValue( keys, ';', 'compute', '=', compute, status )
    IF_ERROR_RETURN(status=1)

    ! by default not a constant field:
    constant = .false.
    ! reset from keys if present:
    call goVarValue( keys, ';', 'constant', '=', constant, status )
    IF_ERROR_RETURN(status=1)

    ! for safety, always need definition of this variable:
    call goVarValue( keys, ';', 'swap_levels', '=', swap_levels, status )
    if ( status < 0 ) then
      write (gol,'("for safety, define `swap_levels` flag in special settings for 3D file")'); call goErr
      write (gol,'("  filename             : ",a)') trim(filename); call goErr
      write (gol,'("  variable description : ",a)') trim(description); call goErr
      write (gol,'("  special keys         : ",a)') trim(keys); call goErr
      TRACEBACK; status=1; return
    else if ( status > 0 ) then
      TRACEBACK; status=1; return
    end if

    ! open file:
    call file_in%Open( trim(filename), status )
    IF_NOTOK_RETURN(status=1)  

    !! info ..
    !write (gol,'(a,": reading record ",i2," ...")') rname, irec; call goPr

    ! variable id:
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! init grid definition from file, use coordinates from variable:
    call file_in%Get_Grid( varid, grid_in, status, &
                            ugg_to=ugg, start_ind=start_ind(1:2), count_ind=count_ind(1:2) )
    IF_NOTOK_RETURN(status=1)

    ! init climat mask:
    call LE_Data_ClimatMask( climat, cmask, status )
    IF_NOTOK_RETURN(status=1)
    
    ! get name of time dimension:
    call file_in%Inq_Var( varid, status, time_dim_name=time_dim_name )
    IF_NOTOK_RETURN(status=1)
    !! info ...
    !call wrtgol( '      read time record valid for ', t ); call goPr

    ! determine time record index
    if ( constant ) then
      ! info ...
      !write (gol,'("      read first time record (constant field)")'); call goPr
      ! set record index to first:
      irec = 1
      ! valid for a long time:
      tt_rec = (/ AnyDate(), AnyDate() /)
    !
    else if ( tt(1) == tt(2) ) then
      !! info ...
      !call wrtgol( '      read time record valid for ', tt(1) ); call goPr
      ! index of  time record:
      call file_in%Inq_TimeRecord( trim(time_dim_name), tt(1), irec, status, &
                                     keys=keys, climat=cmask )
      IF_NOTOK_RETURN(status=1)
      ! fill 'interval' for which input record is valid with single time:
      tt_rec = (/ tt(1), tt(2) /)
    !
    else
      !! info ...
      !call wrtgol( '      read time record enclosing ', tt(1), ' - ', tt(2) ); call goPr
      ! time interval; time variable might have associated 'bounds'
      ! of which one record should enclose the interval tt (record interval is returned too):
      call file_in%Inq_TimeRecord( trim(time_dim_name), tt, irec, tt_rec, status, &
                                     keys=keys, climat=cmask )
      IF_NOTOK_RETURN(status=1)
    end if
    ! fill into into space:
    start_ind(4) = irec
    count_ind(4) = 1
    
    ! switch:
    select case ( trim(compute) )
    
      ! fill half-level pressure from surface pressure:
      case ( 'sp_to_hp', 'lnsp_to_hp' )
      
        ! get level definition, for (lon,lat,lev) this is the third dimension;
        ! alternatively extract dimension/variable names from key=value pairs in keys:
        call file_in%Get_Levs( varid, trim(keys), levs, status )
        IF_NOTOK_RETURN(status=1)
        
        ! count:
        nlev = levs%nlev
        
        ! storage for half level pressures:
        allocate( values_in(grid_in%nlon,grid_in%nlat,0:nlev), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! storage for surface pressure:
        allocate( ps_in(grid_in%nlon,grid_in%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! some files have (ln)sp with a level dimension; count number of dimensions:
        call file_in%Inq_Var( varid, status, ndim=ndim )
        IF_NOTOK_RETURN(status=1)
        ! switch:
        if ( ndim == 3 ) then
          ! reset 3rd dimension of index space to time: (no level dimension):
          start_ind(3) = irec
          count_ind(3) = 1
          ! read surface pressure:
          call file_in%Get_Var( trim(description), ps_in, units_in, status, &
                                  start=start_ind(1:3), count=count_ind(1:3) )
          IF_NOTOK_RETURN(status=1)
        else if ( ndim == 4 ) then
          ! fill 3rd dimension of index space:
          start_ind(3) = 1
          count_ind(3) = 1
          ! read surface pressure at first level:
          call file_in%Get_Var( trim(description), ps_in, units_in, status, &
                                  start=start_ind, count=count_ind )
          IF_NOTOK_RETURN(status=1)
        else
          write (gol,'("unsupported number of dimensions (",i0,") for (ln)sp in file:")') ndim; call goErr
          write (gol,'("  ",a)') trim(file_in%filename); call goPr
          TRACEBACK; status=1; return
        end if
        
        ! log pressure?
        if ( compute(1:4) == 'lnsp' ) then
          ! check:
          if ( .not. ( (trim(units_in) == '1') .or. (trim(units_in) == '~') ) ) then
            write (gol,'("unexpected log surface pressure units `",a,"`")') trim(units_in); call goErr
            TRACEBACK; status=1; return
          end if
          ! convert:
          ps_in = exp(ps_in)
          ! reset units:
          units_in = 'Pa'
        else
          ! adhoc fix: files from some providers are a bit sloppy in following cf-conventions,
          ! and contain surface pressure in 'Pa' in a variable named 'lnsp' with units '~' ...
          if ( (trim(units_in) == '~') .and. (minval(ps_in) > 500.0e2) .and. (maxval(ps_in) < 1100.0e2) ) then
            ! reset units:
            units_in = 'Pa'
          end if
          ! check:
          if ( trim(units_in) /= 'Pa' ) then
            write (gol,'("expected surface pressure in units `Pa`, got `",a,"`")') trim(units_in); call goErr
            TRACEBACK; status=1; return
          end if
        end if
        
        ! extend to half levels:
        call levs%Ps_to_Ph( ps_in, values_in, status )
        IF_NOTOK_RETURN(status=1)
        
        ! clear:
        deallocate( ps_in, stat=status )
        IF_NOTOK_RETURN(status=1)
        
        ! clear input level definition:
        call levs%Done( status )
        IF_NOTOK_RETURN(status=1)
        
      ! read field, but replace by constant value:
      case ( 'value' )
      
        ! number of levels:
        call file_in%Inq_Var( varid, status, dim=3, dim_length=nlev )
        IF_NOTOK_RETURN(status=1)
        ! fill into into space:
        start_ind(3) = 1
        count_ind(3) = nlev
      
        ! storage:
        allocate( values_in(grid_in%nlon,grid_in%nlat,nlev), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! read data slab:
        call file_in%Get_Var( trim(description), values_in, units_in, status, &
                                start=start_ind, count=count_ind )
        IF_NOTOK_RETURN(status=1)

        ! read value from keys, should be defined:
        call goVarValue( keys, ';', 'value', '=', newval, status )
        IF_NOTOK_RETURN(status=1)
        ! reset:
        values_in = newval
              
      ! no special computations needed:
      case ( 'None' )
      
        ! number of levels:
        call file_in%Inq_Var( varid, status, dim=3, dim_length=nlev )
        IF_NOTOK_RETURN(status=1)
        ! fill into into space:
        start_ind(3) = 1
        count_ind(3) = nlev
      
        ! storage:
        allocate( values_in(grid_in%nlon,grid_in%nlat,nlev), stat=status )
        IF_NOTOK_RETURN(status=1)
        
        ! read data slab:
        call file_in%Get_Var( trim(description), values_in, units_in, status, &
                                start=start_ind, count=count_ind )
        IF_NOTOK_RETURN(status=1)
        
      ! unknown ..
      case default
        write (gol,'("unsupported compute action `",a,"`")') trim(compute); call goErr
        TRACEBACK; status=1; return
        
    end select ! special

    ! close input file:
    call file_in%Close( status )
    IF_NOTOK_RETURN(status=1)

    ! reset units ? sometimes the actual units are different ...
    call goVarValue( keys, ';', 'units', '=', units_in, status )
    IF_ERROR_RETURN(status=1)
    
    ! get conversion factor:
    call LE_Data_ConversionFactor( trim(units_in), trim(units), factor, status, keys=keys )
    IF_NOTOK_RETURN(status=1)
    ! apply ?
    if ( factor /= 1.0 ) then
      values_in = values_in * factor
    end if ! conversion needed ?
    
    !! info ..
    !write (gol,*) rname//': value range : ', minval(values_in), maxval(values_in); call goPr
    
    ! bounds of target array:
    if ( on_halflevels ) then
      lbo = (/        1,        1,                   0 /)
      ubo = (/ ugg%nlon, ugg%nlat, size(values_in,3)-1 /)
    else
      lbo = (/        1,        1,                 1 /)
      ubo = (/ ugg%nlon, ugg%nlat, size(values_in,3) /)
    end if
    
    ! storage for horizontally interpolated field:
    allocate( values_hor(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! remap:
    select case ( trim(mapping) )
      !~
      case ( 'interpol' )
        ! interpolate to grid points (corners) :
        call Grid_Convertors%Ugg_Interpol( grid_in, values_in, ugg, values_hor, status )
        IF_NOTOK_RETURN(status=1)
      !~
      case ( 'area-aver' )
        ! averages over grid cell:
        call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, values_hor, status )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unsupported mapping `",a,"`")') trim(mapping); call goErr
        TRACEBACK; status=1; return
    end select

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! copy into target array:
    if ( swap_levels ) then
      ! swap levels:
      do k = lbo(3), ubo(3)
        values(:,:,k) = values_hor(:,:,ubo(3)+lbo(3)-k)
      end do
    else
      ! copy:
      values = values_hor
    end if
    
    ! clear input field:
    deallocate( values_in, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! clear temporary fields:
    deallocate( values_hor, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear input grid definition:
    call grid_in%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine MeteoField_ReadRecord_3d_ugg
  
  
  ! ***
  
  
  subroutine MeteoField_ReadRecordUV_3d( tref, t, rcF, input, units, &
                                        ugg, mapping, &
                                        values_u, values_v, status, halflevels )
  
    use GO            , only : TDate
    use GO            , only : TrcFile, ReadRc
    use C3PO          , only : T_File_nc, T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)            ::  tref
    type(TDate), intent(in)            ::  t
    type(TrcFile), intent(in)          ::  rcF
    character(len=*), intent(in)       ::  input
    character(len=*), intent(in)       ::  units
    type(T_Grid_Ugg), intent(in)       ::  ugg
    character(len=*), intent(in)       ::  mapping
    real, allocatable, intent(out)     ::  values_u(:,:,:)
    real, allocatable, intent(out)     ::  values_v(:,:,:)
    integer, intent(out)               ::  status
    logical, intent(in), optional      ::  halflevels
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MeteoField_ReadRecordUV_3d'
    
    ! --- local ----------------------------------
    
    logical               ::  on_halflevels
    integer               ::  k
    character(len=1)      ::  keys(2)
    character(len=512)    ::  filenames(2)
    character(len=512)    ::  descriptions(2)
    character(len=512)    ::  varkeys
    character(len=512)    ::  climat
    integer               ::  varid
    character(len=32)     ::  gridtype_in
    type(T_File_Nc)       ::  file_in
    
    ! --- begin ----------------------------------
    
    ! halflevel field ?
    on_halflevels = .false.
    if ( present(halflevels) ) on_halflevels = halflevels
    
    !! info ..
    !write (gol,'(a,": open file ",a," ...")') rname, trim(filename); call goPr

    ! uv keys:
    keys = (/ 'u', 'v' /)
    ! loop over uv:
    do k = 1, 2
      ! form filename:
      call GetFilename( rcF, trim(input)//'.'//keys(k), tref, t, filenames(k), status )
      IF_NOTOK_RETURN(status=1)
      ! read variable description line:
      call ReadRc( rcF, trim(input)//'.'//keys(k)//'.var', descriptions(k), status )
      IF_NOTOK_RETURN(status=1)
    end do ! k
    
    ! special instructions, for example for 3D pressure:
    call ReadRc( rcF, trim(input)//'.special', varkeys, status, default='' )
    IF_ERROR_RETURN(status=1)
    ! climatological time elements:
    call ReadRc( rcF, trim(input)//'.climat', climat, status, default='' )
    IF_ERROR_RETURN(status=1)
    
    ! read record:
    call MeteoField_ReadRecordUV_3d_ugg( t, filenames, descriptions, units, &
                                         ugg, mapping, on_halflevels, &
                                         varkeys, climat, values_u, values_v, status )
    IF_NOTOK_RETURN(status=1)
        
    ! ok
    status = 0

  end subroutine MeteoField_ReadRecordUV_3d
  
  
  ! *
  
  
  subroutine MeteoField_ReadRecordUV_3d_ugg( t, filenames, descriptions, units, &
                                           ugg, mapping, on_halflevels, &
                                           keys, climat, values_u, values_v, status )
  
    use GO            , only : TDate
    use GO            , only : goVarValue, goReadFromLine
    use C3PO          , only : T_File_ugg
    use C3PO          , only : T_Grid_ugg
    use C3PO          , only : T_Levs_Hyb
    use LE_Data_Common, only : LE_Data_ClimatMask
    use LE_Data_Common, only : Grid_Convertors
    use LE_Data_Common, only : LE_Data_ConversionFactor
    
    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)            ::  t
    character(len=*), intent(in)       ::  filenames(2)
    character(len=*), intent(in)       ::  descriptions(2)
    character(len=*), intent(in)       ::  units
    type(T_Grid_Ugg), intent(in)       ::  ugg
    character(len=*), intent(in)       ::  mapping
    logical, intent(in)                ::  on_halflevels
    character(len=*), intent(in)       ::  keys
    character(len=*), intent(in)       ::  climat
    real, allocatable, intent(out)     ::  values_u(:,:,:)
    real, allocatable, intent(out)     ::  values_v(:,:,:)
    integer, intent(out)               ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MeteoField_ReadRecordUV_3d_ugg'
    
    ! --- local ----------------------------------
    
    character(len=64)     ::  compute
    logical               ::  swap_levels
    integer               ::  varid
    type(T_File_ugg)      ::  file_in
    type(T_Grid_ugg)      ::  ugg_in_uv(2)
    type(T_Grid_ugg)      ::  ugg_in
    integer               ::  start_ind(4)
    integer               ::  count_ind(4)
    integer               ::  subset(4)
    integer               ::  nlev
    character(len=64)     ::  time_dim_name
    logical               ::  cmask(6)
    integer               ::  irec
    real, allocatable     ::  values_in_u(:,:,:)    ! (i,j,lev)
    real, allocatable     ::  values_in_v(:,:,:)    ! (i,j,lev)
    real, allocatable     ::  values_in(:,:,:,:)    ! (i,j,lev,uv)
    integer               ::  lbo(3), ubo(3)
    real, allocatable     ::  values_hor(:,:,:,:)    ! (lon,lat,lev,uv)
    character(len=64)     ::  units_in
    real                  ::  factor
    integer               ::  i, j, k
    real                  ::  newval
    
    ! --- begin ----------------------------------
    
    ! compute task, default none:
    compute = 'None'
    call goVarValue( keys, ';', 'compute', '=', compute, status )
    IF_ERROR_RETURN(status=1)

    ! for safety, always need definition of this variable:
    call goVarValue( keys, ';', 'swap_levels', '=', swap_levels, status )
    if ( status < 0 ) then
      write (gol,'("for safety, define `swap_levels` flag in special settings for 3D file")'); call goErr
      do k = 1, 2
        write (gol,'("  filename(",i0,")    : ",a)') k, trim(filenames(k)); call goErr
        write (gol,'("  description(",i0,") : ",a)') k, trim(descriptions(k)); call goErr
      end do
      TRACEBACK; status=1; return
    else if ( status > 0 ) then
      TRACEBACK; status=1; return
    end if
    
    ! init climat mask:
    call LE_Data_ClimatMask( climat, cmask, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over u,v
    do k = 1, 2

      ! open file:
      call file_in%Open( trim(filenames(k)), status )
      IF_NOTOK_RETURN(status=1)
      
      ! variable id:
      call file_in%Inq_VarID( trim(descriptions(k)), varid, status )
      IF_NOTOK_RETURN(status=1)

      ! init grid definition from file, use coordinates from variable;
      if ( k == 1 ) then
        ! clip to target domain:
        call file_in%Get_Grid( varid, ugg_in_uv(k), status, &
                              ugg_to=ugg, start_ind=start_ind(1:2), count_ind=count_ind(1:2) )
        IF_NOTOK_RETURN(status=1)
      else
        ! define index space for "v", different from current if staggered:
        if ( trim(compute) == 'uv_stag_to_center' ) then
          count_ind(1) = count_ind(1) - 1  ! from nlon+1 to nlon
          count_ind(2) = count_ind(2) + 1  ! from nlat   to nlat+1
        end if
        ! fill:
        subset = (/ start_ind(1), start_ind(1)+count_ind(1)-1, &
                    start_ind(2), start_ind(2)+count_ind(2)-1 /)
        ! extract subset:
        call file_in%Get_Grid( varid, ugg_in_uv(k), status, subset=subset )
        IF_NOTOK_RETURN(status=1)
      end if

      ! number of levels:
      call file_in%Inq_Var( varid, status, dim=3, dim_length=nlev )
      IF_NOTOK_RETURN(status=1)
      ! fill into into space:
      start_ind(3) = 1
      count_ind(3) = nlev

      ! get name of time dimension:
      call file_in%Inq_Var( varid, status, time_dim_name=time_dim_name )
      IF_NOTOK_RETURN(status=1)
      ! index of  time record:
      call file_in%Inq_TimeRecord( trim(time_dim_name), t, irec, status, &
                                         keys=keys, climat=cmask )
      IF_NOTOK_RETURN(status=1)
      ! fill into into space:
      start_ind(4) = irec
      count_ind(4) = 1

      ! switch:
      select case ( trim(compute) )

        ! index space (0:nlon,1:nlat) and (1:nlon,0:nlat)
        case ( 'None', 'uv_stag_to_center' )

          ! storage for u or v:
          if ( k == 1 ) then
            ! storage:
            allocate( values_in_u(0:ugg_in_uv(k)%nlon-1,1:ugg_in_uv(k)%nlat,nlev), stat=status )
            IF_NOTOK_RETURN(status=1)
            ! read data slab:
            call file_in%Get_Var( trim(descriptions(k)), values_in_u, units_in, status, &
                                    start=start_ind, count=count_ind )
            IF_NOTOK_RETURN(status=1)
          else
            ! storage:
            allocate( values_in_v(1:ugg_in_uv(k)%nlon,0:ugg_in_uv(k)%nlat-1,nlev), stat=status )
            IF_NOTOK_RETURN(status=1)
            ! read data slab:
            call file_in%Get_Var( trim(descriptions(k)), values_in_v, units_in, status, &
                                    start=start_ind, count=count_ind )
            IF_NOTOK_RETURN(status=1)
          end if

        ! index space (1:nlon,1:nlat)
        case ( 'uv_stag1_to_center' )

          ! storage for u or v:
          if ( k == 1 ) then
            ! storage:
            allocate( values_in_u(1:ugg_in_uv(k)%nlon,1:ugg_in_uv(k)%nlat,nlev), stat=status )
            IF_NOTOK_RETURN(status=1)
            ! read data slab:
            call file_in%Get_Var( trim(descriptions(k)), values_in_u, units_in, status, &
                                    start=start_ind, count=count_ind )
            IF_NOTOK_RETURN(status=1)
          else
            ! storage:
            allocate( values_in_v(1:ugg_in_uv(k)%nlon,1:ugg_in_uv(k)%nlat,nlev), stat=status )
            IF_NOTOK_RETURN(status=1)
            ! read data slab:
            call file_in%Get_Var( trim(descriptions(k)), values_in_v, units_in, status, &
                                    start=start_ind, count=count_ind )
            IF_NOTOK_RETURN(status=1)
          end if

        ! unknown ..
        case default
          write (gol,'("unsupported compute action `",a,"`")') trim(compute); call goErr
          TRACEBACK; status=1; return

      end select ! special

      ! reset units ? sometimes the actual units are different ...
      call goVarValue( keys, ';', 'units', '=', units_in, status )
      IF_ERROR_RETURN(status=1)

      ! get conversion factor:
      call LE_Data_ConversionFactor( trim(units_in), trim(units), factor, status )
      IF_NOTOK_RETURN(status=1)
      ! apply ?
      if ( factor /= 1.0 ) then
        values_in_u = values_in_u * factor
        values_in_v = values_in_v * factor
      end if ! conversion needed ?
      
      ! close input file:
      call file_in%Close( status )
      IF_NOTOK_RETURN(status=1)
    
    end do  ! u,v

    ! postproc?
    select case ( trim(compute) )
      ! staggered fields (1:nlon+1,1:nlat) and (1:nlon,1:nlat+1),
      ! average to center (1:nlon,1:nlat)
      case ( 'uv_stag_to_center' )
      
        ! init centered grid from staggered grid definitions:
        call ugg_in%InitFromStags( ugg_in_uv(1), ugg_in_uv(2), status )
        IF_NOTOK_RETURN(status=1)

        ! shared storage:
        allocate( values_in(ugg_in%nlon,ugg_in%nlat,nlev,2), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! loop over levels:
        do k = 1, nlev
          ! loop over centers:
          do j = 1, ugg_in%nlat
            do i = 1, ugg_in%nlon
              ! vector in center is computed as vector sum of 'u' and 'v' vectors;
              ! these two are computed as averages from the vectors at oposite sites;
              ! note that lower bounds are 0-based:
              !
              !                       i,j      
              !                        |
              !                    +---o---+      
              !                    |       |      
              !  unormal,u   i-1,j o--     o-- i,j  
              !                    |   |   |      
              !                    +---o---+      
              !                      i,j-1      
              !                   vnormal,v
              !
              values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i-1,j  ,k) + &
                                         ugg_in%unormal(:,i  ,j  )*values_in_u(i  ,j  ,k)   ) &
                                 + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j-1,k) + &
                                         ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j  ,k)   )
            end do ! i
          end do ! j
          ! loop over centers
        end do ! k

      ! staggered fields (1:nlon,1:nlat) and (1:nlon,1:nlat),
      ! average to center (1:nlon,1:nlat)
      case ( 'uv_stag1_to_center' )
      
        ! init centered grid from staggered grid definitions:
        call ugg_in%InitFromStags( ugg_in_uv(1), ugg_in_uv(2), status, padding=.true. )
        IF_NOTOK_RETURN(status=1)

        ! shared storage:
        allocate( values_in(ugg_in%nlon,ugg_in%nlat,nlev,2), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! loop over levels:
        do k = 1, nlev
          ! loop over centers:
          do j = 1, ugg_in%nlat
            do i = 1, ugg_in%nlon
              ! vector in center is computed as vector sum of 'u' and 'v' vectors;
              ! these two are computed as averages from the vectors at oposite sites;
              ! note that lower bounds of normal vectors in ugg_in include the initial 0,
              ! but values are available from index 1 onwards only:
              !
              !                       i,j      
              !                        |
              !                    +---o---+      
              !                    |       |      
              !  unormal,u   i-1,j o--     o-- i,j  
              !                    |   |   |      
              !                    +---o---+      
              !                      i,j-1      
              !                   vnormal,v
              !
              if ( (i == 1) .and. (j == 1) ) then
                ! use first values also for edge:
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i  ,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i  ,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j  ,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j  ,k)   )
              else if ( i == 1 .and. j == ugg_in%nlat ) then
                ! use first/last value also for edge:
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i  ,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i  ,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j-1,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j-1,k)   )
              else if ( i == ugg_in%nlon  .and. (j == 1) ) then
                ! use last/first value also for edge:
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i-1,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i-1,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j  ,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j  ,k)   )
              else if ( i == ugg_in%nlon .and. j == ugg_in%nlat  ) then
                ! use last values also for edge:
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i-1,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i-1,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j-1,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j-1,k)   )
              else if ( i == 1 ) then
                ! use first value also for edge:
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i  ,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i  ,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j-1,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j  ,k)   )
              else if ( i == ugg_in%nlon ) then
                ! use last value also for edge:
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i-1,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i-1,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j-1,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j  ,k)   )
              else if ( j == 1 ) then
                ! use first value also for edge:
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i-1,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i  ,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j  ,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j  ,k)   )
              else if ( j == ugg_in%nlat ) then
                ! use last value also for edge:
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i-1,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i  ,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j-1,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j-1,k)   )
              else
                values_in(i,j,k,:) = 0.5*( ugg_in%unormal(:,i-1,j  )*values_in_u(i-1,j  ,k) + &
                                           ugg_in%unormal(:,i  ,j  )*values_in_u(i  ,j  ,k)   ) &
                                   + 0.5*( ugg_in%vnormal(:,i  ,j-1)*values_in_v(i  ,j-1,k) + &
                                           ugg_in%vnormal(:,i  ,j  )*values_in_v(i  ,j  ,k)   )
              end if ! edge?
            end do ! i
          end do ! j
          ! loop over centers
        end do ! k

      ! no computations
      case default

        ! uv non-staggerd, should have same shape ...
        if ( any( shape(values_in_u) /= shape(values_in_v) ) ) then
          write (gol,'("non-staggered (u,v) arrays should have same shape:")'); call goErr
          write (gol,'("  u shape : ",3i6)') shape(values_in_u); call goErr
          write (gol,'("  v shape : ",3i6)') shape(values_in_v); call goErr
          TRACEBACK;status=1;return
        end if

        ! init shared source grid:
        call ugg_in%Init( ugg_in_uv(1)%longitude, ugg_in_uv(1)%latitude, status )
        IF_NOTOK_RETURN(status=1)

        ! shared storage:
        allocate( values_in(ugg_in%nlon,ugg_in%nlat,nlev,2), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! copy:
        values_in(:,:,:,1) = values_in_u
        values_in(:,:,:,2) = values_in_v

    end select
    
    ! bounds of target array:
    lbo = (/        1,        1, lbound(values_in,3) /)
    ubo = (/ ugg%nlon, ugg%nlat, ubound(values_in,3) /)
    
    ! storage for horizontally interpolated field:
    allocate( values_hor(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3),2), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! remap:
    select case ( trim(mapping) )
      !~
      case ( 'interpol' )
        ! interpolate to grid points (corners) :
        call Grid_Convertors%Ugg_Interpol( ugg_in, values_in(:,:,:,1), values_in(:,:,:,2), &
                                              ugg, values_hor(:,:,:,1), values_hor(:,:,:,2), status )
        IF_NOTOK_RETURN(status=1)
      !~
      case ( 'area-aver' )
        ! averages over grid cell:
        call Grid_Convertors%Ugg_AreaAver( ugg_in, values_in (:,:,:,1), values_in (:,:,:,2), &
                                              ugg, values_hor(:,:,:,1), values_hor(:,:,:,2), status )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unsupported mapping `",a,"`")') trim(mapping); call goErr
        TRACEBACK; status=1; return
    end select

    ! storage for target field:
    allocate( values_u(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( values_v(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! copy into target array:
    if ( swap_levels ) then
      ! swap levels:
      do k = lbo(3), ubo(3)
        values_u(:,:,k) = values_hor(:,:,ubo(3)+lbo(3)-k,1)
        values_v(:,:,k) = values_hor(:,:,ubo(3)+lbo(3)-k,2)
      end do
    else
      ! copy:
      values_u = values_hor(:,:,:,1)
      values_v = values_hor(:,:,:,2)
    end if

    ! clear input fields:
    deallocate( values_in_u, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( values_in_v, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( values_in, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! clear temporary fields:
    deallocate( values_hor, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear input grid definition:
    do k = 1, 2
      call ugg_in_uv(k)%Done( status )
      IF_NOTOK_RETURN(status=1)
    end do
    call ugg_in%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine MeteoField_ReadRecordUV_3d_ugg
  

end module LE_Data_Meteo
