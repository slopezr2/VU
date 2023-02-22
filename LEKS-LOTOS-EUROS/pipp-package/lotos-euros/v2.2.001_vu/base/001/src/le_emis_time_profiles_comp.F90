!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################

module LE_Emis_Time_Prof_Comp

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Time_Prof_Comp
  public  ::  LE_Emis_Time_Prof_Comp_Init, LE_Emis_Time_Prof_Comp_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Time_Prof_Comp'
  
  
  ! --- types --------------------------------
  
  ! storage for emission composition:
  type T_Emis_Time_Prof_Comp
    ! dimensions:
    integer                   ::  ncat
    integer                   ::  ncountry
    integer                   ::  year
    integer                   ::  nhour  ! (8760 for a non leap year)
    integer                   ::  ncomp ! number of emitted tracers
    ! fraction assigned to component:
    real, allocatable         ::  profile(:,:,:,:)  ! (hour,sector,country,component)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_Time_Prof_Comp_Init( timeprof, query, year, &
                                          cat_codes, country_codes, component_names, &
                                          status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue
    use Go, only :  calc_DayNumber, days_in_year
    use MDF, only : MDF_Open, MDF_Close
    use MDF, only : MDF_Inquire
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Inquire_Variable, MDF_Get_Var
    use MDF, only : MDF_Get_Att
    use MDF, only : MDF_NETCDF, MDF_READ
    use LE_Emis_Tools, only : ShortSNAP_to_Code

    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof_Comp), intent(out)  ::  timeprof
    character(len=*), intent(in)              ::  query
    integer, intent(in)                       ::  year
    character(len=*), intent(in)              ::  cat_codes(:)
    character(len=3), intent(in)              ::  country_codes(:)
    character(len=*), intent(in)              ::  component_names(:)
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Comp_Init'
    integer, parameter          ::  max_years = 10
    
    ! --- local -------------------------------
        
    character(len=512)              ::  filename
    character(len=1)                ::  comment
    character(len=1)                ::  sep
    logical                         ::  exist
    
    integer                         ::  hid
    integer                         ::  dimid
    integer                         ::  varid
    
    real, allocatable               ::  profile_help(:,:,:,:)
    integer                         ::  ncountry
    integer                         ::  ncat
    integer                         ::  ncomp
    integer                         ::  len_country_code, len_component_name
    character(len=3), allocatable   ::  countries_time_prof(:)
    integer, allocatable            ::  cat_codes_time_prof(:)
    character(len=3)                ::  country_name
    integer                         ::  icountry, icat, icomp
    integer, allocatable            ::  country_match(:), cat_match(:), component_match(:)
    character(len=8)                ::  cat_code
    character(len=10), allocatable  ::  components_time_prof(:)
    character(len=10)               ::  component_name
  
    integer                         ::  nfile, ifile

    ! --- begin -------------------------------
    
    ! store:
    timeprof%year         = year
    
    ! store currently known dimensions:
    timeprof%ncat         = size(cat_codes)
    timeprof%ncountry     = size(country_codes(:))
    timeprof%ncomp        = size(component_names(:))
      
    ! extract filename, raise error message if not defined (status<0):
    call goVarValue( trim(query), ';', 'file', '=', filename, status )
    if ( status < 0 ) then
      write (gol,'("no keyword `file` in time profile query : ",a)') trim(query); call goErr
      TRACEBACK; status=1; return
    end if
    IF_NOTOK_RETURN(status=1)
    
    ! seperation character:
    sep = ';'
    call goVarValue( trim(query), ';', 'sep', '=', sep, status )
    IF_ERROR_RETURN(status=1)
    ! comment character:
    comment = '#'
    call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    IF_ERROR_RETURN(status=1)
    
    ! file should be present:
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! open file:
    call MDF_Open( trim(filename), MDF_NETCDF, MDF_READ, hid, status )
    IF_NOTOK_RETURN(status=1)

    ! number of countries in time profile file
    call MDF_Inq_DimID( hid, 'country', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ncountry )
    IF_NOTOK_RETURN(status=1)
    ! number of SNAP categories in time profile file
    call MDF_Inq_DimID( hid, 'emis_category', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ncat )
    IF_NOTOK_RETURN(status=1)
    ! number of components in time profile file
    call MDF_Inq_DimID( hid, 'tracers', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=ncomp )
    IF_NOTOK_RETURN(status=1)                    
    ! length of component name in time profile file
    call MDF_Inq_DimID( hid, 'tracer_name_len', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=len_component_name )
    IF_NOTOK_RETURN(status=1)                    
    ! number of hours in time profile file
    call MDF_Inq_DimID( hid, 'hours', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=timeprof%nhour )
    IF_NOTOK_RETURN(status=1)                    

    ! no strage yet ?
    if ( .not. allocated(timeprof%profile) ) then
      ! allocate timeprofile with dimension hour per year
      allocate( timeprof%profile(timeprof%nhour,timeprof%ncat,timeprof%ncountry,timeprof%ncomp) )
      ! by default fill with -999.0 
      timeprof%profile = -999.0
    end if

    ! read country codes:
    allocate( countries_time_prof(ncountry) )
    len_country_code = 3
    call MDF_Get_StrArr( hid, 'country_code'  , len_country_code, ncountry, countries_time_prof, status )
    IF_NOTOK_RETURN(status=1)          

    ! read emission category codes (LE SNAP):
    allocate( cat_codes_time_prof(ncat) )
    call MDF_Inq_VarID( hid, 'sectorID', varid, status )
    if ( status /= 0 ) then
      call MDF_Inq_VarID( hid, 'emis_cat', varid, status )
      IF_NOTOK_RETURN(status=1)
    end if
    call MDF_Get_Var( hid, varid, cat_codes_time_prof, status )
    IF_NOTOK_RETURN(status=1)

    ! read component names:
    allocate( components_time_prof(ncomp) )
    call MDF_Get_StrArr( hid, 'tracer_name', len_component_name, ncomp, components_time_prof, status )
    IF_NOTOK_RETURN(status=1)

    ! match countries in time profiles with countries in model
    allocate( country_match(ncountry))
    country_match = -1
    do icountry = 1, ncountry
      country_name = countries_time_prof(icountry)
      call goMatchValue(trim(country_name), country_codes, country_match(icountry), status)
      IF_NOTOK_RETURN(status=1)
    end do

    ! match categories in time profiles with categories in model
    allocate( cat_match(ncat) )
    cat_match = -1
    do icat = 1, ncat
      call ShortSNAP_to_Code( cat_codes_time_prof(icat), cat_code, status )
      IF_NOTOK_RETURN(status=1)
      call goMatchValue(cat_code, cat_codes, cat_match(icat), status)
      IF_NOTOK_RETURN(status=1)
    end do

    ! match components in time profiles with components in model
    allocate( component_match(ncomp))
    component_match = -1
    do icomp = 1, ncomp
      component_name = components_time_prof(icomp)
      call goMatchValue(trim(component_name), component_names, component_match(icomp), status)
      if ( status /= 0 .and. ( trim(component_name) == 'PM10' .or. trim(component_name) == 'pm10') ) then
        write( gol, '("Match pm10 profile with pm25_pm10")' ) ; call goPr
        call goMatchValue( 'pm25_pm10', component_names, component_match(icomp), status )
        IF_NOTOK_RETURN(status=1)
      else if ( status /= 0 .and. ( trim(component_name) == 'PM25' .or. trim(component_name) == 'pm25') ) then
        call goMatchValue( 'pm2_5', component_names, component_match(icomp), status )
        IF_NOTOK_RETURN(status=1)
      else if ( status /= 0 ) then
        write(gol, '("Cannot match component: ", a)' ) trim(component_name) ; call goErr
        TRACEBACK;status=1;return
      end if
    end do

    ! read variable
    allocate( profile_help(ncountry,ncat,ncomp,timeprof%nhour) )          
    call MDF_Inq_VarID( hid, 'time_factors', varid, status )
    IF_NOTOK_RETURN(status=1)           
    call MDF_Get_Var( hid, varid, profile_help, status )          
    IF_NOTOK_RETURN(status=1)

    ! loop over indices in file:
    do icat = 1, ncat
      do icountry = 1, ncountry
        do icomp = 1, ncomp
          ! store following order used in LE:
          timeprof%profile(:,cat_match(icat),country_match(icountry),component_match(icomp)) = profile_help(icountry,icat,icomp,:)
          ! check ...
          if ( abs(sum(timeprof%profile(:,cat_match(icat),country_match(icountry),component_match(icomp)))/timeprof%nhour -1 )  >= 1e-2)  then
            write(gol, '("Time profile do not have average 1 ")'); call goErr            
            write(gol, '("Country: ", a)') trim(countries_time_prof(icountry)); call goErr
            write(gol, '("Cat code: ", i0)' ) cat_codes_time_prof(icat); call goErr
            write(gol, '("Component_name: ", a)' ) trim(components_time_prof(icomp))  ; call goErr
            write(gol, '("Average = ", f10.5)' ) sum(timeprof%profile(:,cat_match(icat),country_match(icountry),component_match(icomp)))/timeprof%nhour   ; call goErr
            TRACEBACK; status=1; return
          end if 
        end do
      end do
    end do

    ! clear          
    deallocate( countries_time_prof )
    deallocate( cat_codes_time_prof )
    deallocate( components_time_prof )
    deallocate( country_match )
    deallocate( cat_match )
    deallocate( component_match )
    deallocate( profile_help )

    ! close file:
    call MDF_close( hid, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Time_Prof_Comp_Init
  
  
  ! ***
  

  subroutine LE_Emis_Time_Prof_Comp_Done( timeprof, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof_Comp), intent(inout) ::  timeprof
    integer, intent(out)                       ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Comp_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    if ( allocated(timeprof%profile     ) ) deallocate( timeprof%profile    )

    ! ok
    status = 0

  end subroutine LE_Emis_Time_Prof_Comp_Done
  
  
  ! ***
  

  ! Read 2d character variable in a variable with oversized length

  subroutine MDF_Get_StrArr( hid, varname, nchar, nval, values, status )
  
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var

    ! --- in/out -----------------------

    integer, intent(in)             ::  hid
    character(len=*), intent(in)    ::  varname
    integer, intent(in)             ::  nchar, nval
    character(len=*), intent(out)   ::  values(nval)
    integer, intent(out)            ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MDF_Get_StrArr'
    
    ! --- local -------------------------

    integer                 ::  varid
    character(len=nchar)    ::  tmp_values(nval)
    integer                 ::  ival
    integer                 ::  l

    ! --- begin -------------------------

    ! get variable id given name:      
    call MDF_Inq_VarID( hid, varname, varid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read data:
    call MDF_Get_Var( hid, varid, tmp_values, status )
    IF_NOTOK_RETURN(status=1)
    
    ! copy:
    do ival = 1, nval
      ! seach <nul> characters:
      l = index( tmp_values(ival), char(0) )
      ! found a <nul> character ?
      if ( l < 1 ) then
        ! no <nul>, just copy:
        values(ival) = trim(tmp_values(ival))
      else if ( l == 1 ) then
        ! first character is <nul>, thus empty:
        values(ival) = ''
      else
        ! found a <nul> somewhere in the string; copy part before:
        values(ival) = trim(tmp_values(ival)(1:l-1))
      end if
    end do

    ! ok
    status = 0

  end subroutine MDF_Get_StrArr


end module LE_Emis_Time_Prof_Comp

