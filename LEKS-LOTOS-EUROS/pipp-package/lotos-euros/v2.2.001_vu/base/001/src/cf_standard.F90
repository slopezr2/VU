!###############################################################################
!
! NAME
!   CF_Standard  -  tools for standard names and units following CF convention
!
! USAGE
!
!   use CF_Standard
!
!   ! initialize module by reading standard name table:
!   call CF_Standard_Init( 'cf-standard-name-table-v15.xml', status )
!     character(len=*_, intent(in)        ::  cf_standard_name_table
!     integer, intent(out)                ::  status
!
!   !
!   ! create standard name from components, return units too;
!   ! currently supported:
!   !
!   !    'mass_concentration_of_X_in_Y', X='ozone', Y='air'
!   !
!   call CF_Standard_Get_Name_and_Units( name, &
!                  cf_standard_name, cf_units, status [,X=] [,Y=] )
!     character(len=*), intent(in)                ::  name
!     character(len=*), intent(out)               ::  cf_standard_name
!     character(len=*), intent(out)               ::  cf_units
!     integer, intent(out)                        ::  status
!     character(len=*), intent(in), optional      ::  X, Y
!
!   ! get standard units for standard name:
!   call CF_Standard_Get_Units( cf_standard_name, cf_units, status )
!     character(len=*), intent(in)                ::  cf_standard_name
!     character(len=*), intent(out)               ::  cf_units
!     integer, intent(out)                        ::  status
!
!   ! get conversion factor from <units> to <cf_units>
!   call CF_Standard_Get_Conversion( units, cf_units, status )
!     character(len=*), intent(in)                ::  units
!     character(len=*), intent(in)                ::  cf_units
!     real, intent(out)                           ::  factor
!     integer, intent(out)                        ::  status
!
!   ! done with module:
!   call CF_Standard_Done( status )
!     integer, intent(out)                        ::  status
!
!
! EXTERNAL LIBRARIES
!
!   Module uses the 'UDUnits' interface to either
!   the udunits version 1 or 2 library.
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_UDUNITS_NOTOK_RETURN(action) if (status/=UDUNITS_NOERR) then; gol=trim(UDUnits_StrError(status)); call goErr; TRACEBACK; action; return; end if
!
! safety check ...
#define IF_NOT_INIT_RETURN(action) if (.not.module_init) then; write (gol,'("module not initialized")'); call goErr; TRACEBACK; action; return; end if
!
!###############################################################################

module CF_Standard

  use GO     , only : gol, goErr, goPr
  use UDUnits, only : UDUNITS_NOERR, UDUnits_StrError
  
  ! table read from xml file:
  use Standard_Name_Table, only : T_Standard_Name_Table
  
  implicit none
  

  ! --- in/out -----------------------------------
  
  private
  
  public  ::  CF_Standard_Init, CF_Standard_Done
  public  ::  CF_Standard_Get_Name_and_Units
  public  ::  CF_Standard_Get_Units
  public  ::  CF_Standard_Get_Conversion

  
  ! --- const ------------------------------------
  
  character(len=*), parameter  ::  mname = 'CF_Standard'

  
  ! --- local ------------------------------------

  ! flag ...
  logical                         ::  module_init = .false.
  
  ! storage for table:
  type(T_Standard_Name_Table)     ::  cf_table



contains


  ! ====================================================================
  ! ===
  ! === module
  ! ===
  ! ====================================================================


  subroutine CF_Standard_Init( cf_standard_name_table, status )
  
    use Standard_Name_Table, only : Standard_Name_Table_Init

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)        ::  cf_standard_name_table
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CF_Standard_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! standard names and units:
    call Standard_Name_Table_Init( cf_table, trim(cf_standard_name_table), status )
    IF_NOTOK_RETURN(status=1)

    ! reset flag:
    module_init = .true.

    ! ok
    status = 0
    
  end subroutine CF_Standard_Init


  ! ***


  subroutine CF_Standard_Done( status )
  
    use Standard_Name_Table, only : Standard_Name_Table_Done

    ! --- in/out ---------------------------------
    
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CF_Standard_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    IF_NOT_INIT_RETURN(status=1)

    ! done with standard names and units:
    call Standard_Name_Table_Done( cf_table, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine CF_Standard_Done


  ! ***


  subroutine CF_Standard_Get_Name_and_Units( &
                              name, &
                              cf_standard_name, cf_units, &
                              status, X, Y )

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)                ::  name
    character(len=*), intent(out)               ::  cf_standard_name
    character(len=*), intent(out)               ::  cf_units
    integer, intent(out)                        ::  status
    character(len=*), intent(in), optional      ::  X
    character(len=*), intent(in), optional      ::  Y

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CF_Standard_Get_Name_and_Units'
    
    ! --- local ----------------------------------
    
    character(len=256)  ::  cf_tracer
    integer             ::  id
    
    ! --- begin ----------------------------------
    
    ! check ...
    IF_NOT_INIT_RETURN(status=1)

    ! switch:
    select case ( trim(name) )

      ! concentrations:
      case ( 'mass_concentration_of_X_in_Y' )

        ! check ...
        if ( .not. present(X) ) then
          write (gol,'("argument X required for name : ",a)') trim(name); call goErr
          TRACEBACK; status=1; return
        end if
        ! translate:
        call CF_Standard_Get_Tracer( X, cf_tracer, status )
        IF_NOTOK_RETURN(status=1)

        ! check ...
        if ( .not. present(Y) ) then
          write (gol,'("argument Y required for name : ",a)') trim(name); call goErr
          TRACEBACK; status=1; return
        end if

        ! combine:
        write (cf_standard_name,'("mass_concentration_of_",a,"_in_",a)') trim(cf_tracer), trim(Y)

        ! search ...
        call CF_Standard_Get_ID( cf_standard_name, id, status )
        IF_ERROR_RETURN(status=1)
        if ( status == 0 ) then
          ! extract:
          cf_units = trim(cf_table%entry(id)%canonical_units)
        else
          ! gues something ..
          cf_units = 'kg m-3'
        end if

      ! problem ...
      case default
        write (gol,'("name not supported : ",a)') trim(name); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0
    
  end subroutine CF_Standard_Get_Name_and_Units


  ! ***


  subroutine CF_Standard_Get_Tracer( name, cf_name, status )

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)                ::  name
    character(len=*), intent(out)               ::  cf_name
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CF_Standard_Get_Tracer'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(name) )
      case ( 'ozone' )
        cf_name = 'ozone'
      case ( 'nitrogenDioxide' )
        cf_name = 'nitrogen_dioxide'
      case ( 'sulphurDioxide' )
        cf_name = 'sulfur_dioxide'
      case ( 'carbonMonoxide' )
        cf_name = 'carbon_monoxide'
      case ( 'pm10Aerosol' )
        cf_name = 'pm10_aerosol'
      case default
        ! assume it is ok already, will be checked later on:
        cf_name = trim(name)
    end select
    
    ! ok
    status = 0
    
  end subroutine CF_Standard_Get_Tracer


  ! ***


  ! return index number in table

  subroutine CF_Standard_Get_ID( cf_standard_name, id, status )

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)                ::  cf_standard_name
    integer, intent(out)                        ::  id
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CF_Standard_Get_ID'
    
    ! --- local ----------------------------------
    
    integer                         ::  i, j

    ! --- begin ----------------------------------
    
    ! check ...
    IF_NOT_INIT_RETURN(status=1)

    ! not found yet ..
    id = -1

    ! loop over all entries:
    do i = 1, size(cf_table%entry)
      ! compare ...
      if ( trim(cf_table%entry(i)%id) == trim(cf_standard_name) ) then
        ! store:
        id = i
        ! found!
        exit
      end if
    end do   ! CF table entries

    ! not found yet ? then check alia:
    if ( id < 0 ) then
      ! loop over all alia:
      do j = 1, size(cf_table%alias)
        ! compare ...
        if ( trim(cf_table%alias(j)%id) == trim(cf_standard_name) ) then
          ! match; now search real entries:
          do i = 1, size(cf_table%entry)
            ! compare ...
            if ( trim(cf_table%entry(i)%id) == trim(cf_table%alias(j)%entry_id) ) then
              ! store:
              id = i
              ! leave:
              exit
            end if
          end do   ! CF table entries
          ! found!
          exit
        end if
      end do   ! CF table alia
    end if
    
    ! set return status:
    if ( id < 0 ) then
      ! warning:
      status = -1
    else
      ! ok:
      status = 0
    end if
    
  end subroutine CF_Standard_Get_ID


  ! ***


  subroutine CF_Standard_Get_Units( cf_standard_name, cf_units, status )

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)                ::  cf_standard_name
    character(len=*), intent(out)               ::  cf_units
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CF_Standard_Get_Units'
    
    ! --- local ----------------------------------
    
    integer                         ::  id

    ! --- begin ----------------------------------
    
    ! check ...
    IF_NOT_INIT_RETURN(status=1)

    ! get index:
    call CF_Standard_Get_ID( cf_standard_name, id, status )
    IF_ERROR_RETURN(status=0)
    ! found ?
    if ( status == 0 ) then
      ! extract:
      cf_units = trim(cf_table%entry(id)%canonical_units)
    else
      write (gol,'("id not found in cf standard name table : ",a)') trim(cf_standard_name); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine CF_Standard_Get_Units


  ! ***


  subroutine CF_Standard_Get_Conversion( units, cf_units, factor, status )
  
    use UDUnits, only : UDUnits_ConversionFactor

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)                ::  units
    character(len=*), intent(in)                ::  cf_units
    real, intent(out)                           ::  factor
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CF_Standard_Get_Conversion'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    IF_NOT_INIT_RETURN(status=1)

    ! same already?
    if ( trim(units) == trim(cf_units) ) then

      ! same:
      factor = 1.0

    ! trap some common equivalents,
    ! this avoids a request for the udunits library in most cases:
    else if ( ( (trim(units) == '(0 - 1)'                    ) .and. (trim(cf_units) == '1'            ) ) .or. &
              ( (trim(units) == 'degrees'                    ) .and. (trim(cf_units) == 'degrees_east' ) ) .or. &
              ( (trim(units) == 'degrees'                    ) .and. (trim(cf_units) == 'degrees_north') ) .or. &
              ( (trim(units) == 'kg kg**-1'                  ) .and. (trim(cf_units) == 'kg/kg'        ) ) .or. &
              ( (trim(units) == 'm s**-1'                    ) .and. (trim(cf_units) == 'm/s'          ) ) .or. &
              ( (trim(units) == 'm**3 m**-3'                 ) .and. (trim(cf_units) == 'm3/m3'        ) ) .or. &
              ( (trim(units) == 'J m**-2 s**-1'              ) .and. (trim(cf_units) == 'J/m2/s'       ) ) .or. &
              ( (trim(units) == 'm of water equivalent'      ) .and. (trim(cf_units) == 'm'            ) ) .or. &
              ( (trim(units) == 'm of water equivalent s**-1') .and. (trim(cf_units) == 'm/s'          ) )      ) then

       ! same:
       factor = 1.0

    ! use udunits:
    else

      ! unit conversion factor:
      call UDUnits_ConversionFactor( trim(units), trim(cf_units), factor, status )
      if ( status /= UDUNITS_NOERR ) then
        write (gol,'("from conversion of model units to CF units:")'); call goErr
        write (gol,'("  model units  : ",a)') trim(units); call goErr
        write (gol,'("  CF units     : ",a)') trim(cf_units); call goErr
        TRACEBACK; status=1; return
      end if

    end if

    ! ok
    status = 0
    
  end subroutine CF_Standard_Get_Conversion


end module CF_Standard



