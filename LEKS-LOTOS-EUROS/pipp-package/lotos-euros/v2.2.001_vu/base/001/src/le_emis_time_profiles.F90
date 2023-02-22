!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################

module LE_Emis_Time_Prof

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Time_Prof
  public  ::  LE_Emis_Time_Prof_Init, LE_Emis_Time_Prof_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Time_Prof'
  
  
  ! --- types --------------------------------
  
  ! storage for emission composition:
  type T_Emis_Time_Prof
    ! dimensions:
    integer                   ::  ncat
    integer                   ::  ncountry
    integer                   ::  year
    integer                   ::  nhour  ! (8760 for a non leap year)
    ! fraction assigned to component:
    real, allocatable         ::  profile(:,:,:)  ! (hour,sector,country)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_Time_Prof_Init( timeprof, query, year, &
                                     cat_codes, country_codes, status )

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
    use LE_Emis_Tools, only : MDF_Get_StrArr

    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof), intent(out)       ::  timeprof
    character(len=*), intent(in)              ::  query
    integer, intent(in)                       ::  year
    character(len=*), intent(in)              ::  cat_codes(:)
    character(len=3), intent(in)              ::  country_codes(:)
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Init'

    integer, parameter          ::  max_years = 10
    
    ! --- local -------------------------------
        
    character(len=512)             ::  filename
    character(len=1)               ::  comment
    character(len=1)               ::  sep
    logical                        ::  exist
    
    integer                        ::  hid
    integer                        ::  dimid
    integer                        ::  varid
    
    real, allocatable              ::  profile_help(:,:,:)
    integer                        ::  ncountry
    integer                        ::  ncat
    integer                        ::  len_country_code
    character(len=3), allocatable  ::  countries_time_prof(:)
    character(len=8), allocatable  ::  cat_codes_time_prof(:)
    integer                        ::  emis_cat_code_len
    character(len=3)               ::  country_name
    integer                        ::  icountry, icat
    integer, allocatable           ::  country_match(:), cat_match(:)
    character(len=8)               ::  cat_code
    integer, allocatable           ::  shortsnap(:)
  
    integer                        ::  nfile, ifile

    ! --- begin -------------------------------
    
    ! store:
    timeprof%year         = year
    
    ! store currently known dimensions:
    timeprof%ncat         = size(cat_codes)
    timeprof%ncountry     = size(country_codes(:))
      
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
    ! number of hours in time profile file
    call MDF_Inq_DimID( hid, 'hours', dimid, status )
    IF_NOTOK_RETURN(status=1)
    call MDF_Inquire_Dimension( hid, dimid, status, length=timeprof%nhour )
    IF_NOTOK_RETURN(status=1)                    

    ! no strage yet ?
    if ( .not. allocated(timeprof%profile) ) then
      ! allocate timeprofile with dimension hour per year
      allocate( timeprof%profile(timeprof%nhour,timeprof%ncat,timeprof%ncountry) )
      ! by default fill with -999.0 
      timeprof%profile = -999.0
    end if

    ! read country codes:
    allocate( countries_time_prof(ncountry) )
    len_country_code = 3
    call MDF_Get_StrArr( hid, 'country_code'  , len_country_code, ncountry, countries_time_prof, status )
    IF_NOTOK_RETURN(status=1)          

    ! storage for category codes:
    ! read emission category codes (LE SNAP):
    allocate( cat_codes_time_prof(ncat), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! try latest name first:
    call MDF_Inq_VarID( hid, 'emis_cat_code', varid, status )
    ! found?
    if ( status == 0 ) then
      ! maximum length of name:
      call MDF_Inq_DimID( hid, 'emis_cat_code_len', dimid, status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Inquire_Dimension( hid, dimid, status, length=emis_cat_code_len )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call MDF_Get_StrArr( hid, 'emis_cat_code', emis_cat_code_len, ncat, cat_codes_time_prof, status )
      IF_NOTOK_RETURN(status=1) 
    else
      ! try old variable with integer codes:
      call MDF_Inq_VarID( hid, 'emis_cat', varid, status )
      IF_NOTOK_RETURN(status=1)
      ! storage:
      allocate( shortsnap(ncat), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call MDF_Get_Var( hid, varid, shortsnap, status )
      IF_NOTOK_RETURN(status=1)
      ! covert:
      do icat = 1, ncat
        call ShortSNAP_to_Code( shortsnap(icat), cat_codes_time_prof(icat), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( shortsnap, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

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
      cat_code = cat_codes_time_prof(icat)
      call goMatchValue(cat_code, cat_codes, cat_match(icat), status)
      IF_NOTOK_RETURN(status=1)
    end do

    ! read variable
    allocate( profile_help(ncountry,ncat,timeprof%nhour) )          
    call MDF_Inq_VarID( hid, 'time_factors', varid, status )
    IF_NOTOK_RETURN(status=1)           
    call MDF_Get_Var( hid, varid, profile_help, status )          
    IF_NOTOK_RETURN(status=1)

    ! loop over indices in file:
    do icat = 1, ncat
      do icountry = 1, ncountry
        ! store following order used in LE:
        timeprof%profile(:,cat_match(icat),country_match(icountry)) = profile_help(icountry,icat,:)
        ! check ...
        if ( abs(sum(timeprof%profile(:,cat_match(icat),country_match(icountry)))/timeprof%nhour -1 )  >= 1e-2)  then
          write(gol, '("Time profile do not have average 1 ")'); call goErr            
          write(gol, '("Country: ", a)') trim(countries_time_prof(icountry)); call goErr
          write(gol, '("Cat code: ", a)' ) trim(cat_codes_time_prof(icat))  ; call goErr
          write(gol, '("Average = ", f10.5)' ) sum(timeprof%profile(:,cat_match(icat),country_match(icountry)))/timeprof%nhour   ; call goErr
          TRACEBACK; status=1; return
        end if 
      end do
    end do

    ! clear          
    deallocate( countries_time_prof )
    deallocate( cat_codes_time_prof )
    deallocate( country_match )
    deallocate( cat_match )
    deallocate( profile_help )

    ! close file:
    call MDF_close( hid, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Time_Prof_Init
  
  
  ! ***
  

  subroutine LE_Emis_Time_Prof_Done( timeprof, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof), intent(inout) ::  timeprof
    integer, intent(out)                  ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    if ( allocated(timeprof%profile     ) ) deallocate( timeprof%profile    )

    ! ok
    status = 0

  end subroutine LE_Emis_Time_Prof_Done


end module LE_Emis_Time_Prof

