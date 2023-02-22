!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################

module LE_Emis_EDGAR_Time_Prof

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
    integer                   ::  nhour  ! (8760 for a noN LEAP year)
    real                      ::  nyear
    character(len=64)         ::  years  ! which year are specified in the rc file
    ! fraction assigned to component:
    real, allocatable         ::  profile(:,:,:,:)  ! (year,hour,sector,zone)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_Time_Prof_Init( timeprof, query, cat_codes, &
                                     country_codes, status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue
    use Go, only :  calc_DayNumber, days_in_year
    use MDF, only : MDF_Open, MDF_Close
    use MDF, only : MDF_Inquire
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Inquire_Variable, MDF_Get_Var
    use MDF, only : MDF_Get_Att
    use MDF, only : MDF_NETCDF, MDF_READ

    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof), intent(out)       ::  timeprof
    character(len=*), intent(in)              ::  query
    character(len=*), intent(in)              ::  cat_codes(:)
    character(len=*), intent(in)              ::  country_codes(:)
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Init'
    integer, parameter          ::  max_years = 10
    
    ! --- local -------------------------------
        
    character(len=2048)     ::  filenames
    character(len=512)      ::  filename
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    character(len=32)       ::  fileformat
    logical                 ::  exist
    
    integer                 ::  hid
    integer                 ::  dimid
    integer                 ::  varid
    
    real, allocatable       ::  profile_help(:,:,:)
    integer                 ::  ncountry
    integer                 ::  ncat
    integer                 ::  len_country_code
    integer                 ::  len_emiscatname    
    character(len=4), allocatable  ::  countries_time_prof(:)
    character(len=44), allocatable    ::  cat_codes_time_prof(:)
    character(len=4)        ::  country_name
    integer                 ::  icountry, icat
    integer, allocatable    ::  country_match(:), cat_match(:)
    character(len=44)       ::  cat_code
  
    integer                 ::  nfile, ifile
    character(len=512)      ::  fnames(max_years)
    integer                 ::  years(max_years)
    integer                 ::  iyear
    integer                 ::  min_year, max_year
    

    ! --- begin -------------------------------
    
    ! store currently known dimensions:
    timeprof%ncat         = size(cat_codes)
    timeprof%ncountry     = size(country_codes(:))
      
    ! extract filename, by default empty:
    filenames = ''
    ! line with filenames:
    !   years=2000[,2001[,..]]
    call goVarValue( trim(query), ';', 'files', '=', filenames, status )
    IF_ERROR_RETURN(status=1)
    ! split into list (comma seperated)
    call goSplitString( filenames, nfile, fnames, status )
    IF_NOTOK_RETURN(status=1)
    
    ! extract description of years for which files are valid:
    !   years=2000[,2001[,..]]
    call goVarValue( trim(query), ';', 'years', '=', timeprof%years, status )
    IF_ERROR_RETURN(status=1)
    ! init list of integers with dummy values:
    years = -999
    ! split into list:
    call goSplitString( timeprof%years, nfile, years, status )
    IF_NOTOK_RETURN(status=1)
    
    ! year range of files
    min_year = 9999
    max_year = 0
    do ifile = 1, nfile
      min_year = min(min_year,years(ifile))
      max_year = max(max_year,years(ifile))
    end do
    timeprof%nyear = max_year - min_year

    ! seperation character:
    sep = ';'
    call goVarValue( trim(query), ';', 'sep', '=', sep, status )
    IF_ERROR_RETURN(status=1)
    ! comment character:
    comment = '#'
    call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    IF_ERROR_RETURN(status=1)
    ! format description:
    fileformat = 'csv'
    call goVarValue( trim(query), ';', 'format', '=', fileformat, status )
    IF_ERROR_RETURN(status=1)
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( len_trim(filenames) == 0 ) then   ! no file, fill with default value (-999.0)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ntime_vars = 24, one value for each hour
      timeprof%nhour = days_in_year( years(1) ) * 24
      
      allocate( timeprof%profile(1,timeprof%nhour,timeprof%ncat,timeprof%ncountry) )
      ! by default fill with -999.0
      timeprof%profile = -999.0
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else  ! read table from file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      ! file should be present:
      do ifile = 1, nfile

        filename = fnames(ifile)
        iyear = years(ifile)

        inquire( file=trim(filename), exist=exist )
        if ( .not. exist ) then
          write (gol,'("file not found : ",a)') trim(filename); call goErr
          TRACEBACK; status=1; return
        end if
            
        ! switch per format:
        select case ( trim(fileformat) )

          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          case ( 'nc' )
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
            ! open file:
            call MDF_Open( trim(filename), MDF_NETCDF, MDF_READ, hid, status )
            IF_NOTOK_RETURN(status=1)

            ! number of countries in time profile file
            call MDF_Inq_DimID( hid, 'zone', dimid, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Inquire_Dimension( hid, dimid, status, length=ncountry )
            IF_NOTOK_RETURN(status=1)
            ! length of country/zone names in time profile file
            call MDF_Inq_DimID( hid, 'zone_id', dimid, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Inquire_Dimension( hid, dimid, status, length=len_country_code )
            IF_NOTOK_RETURN(status=1)
            ! number of SNAP categories in time profile file
            call MDF_Inq_DimID( hid, 'emis_category', dimid, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Inquire_Dimension( hid, dimid, status, length=ncat )
            IF_NOTOK_RETURN(status=1)
            ! length of category names in time profile file
            call MDF_Inq_DimID( hid, 'emis_category_id', dimid, status )
            if (status /= 0 ) then
              call MDF_Inq_DimID( hid, 'emis_category_name_len', dimid, status )
              IF_NOTOK_RETURN(status=1)
            endif
            call MDF_Inquire_Dimension( hid, dimid, status, length=len_emiscatname )
            IF_NOTOK_RETURN(status=1)
            ! number of hours in time profile file
            call MDF_Inq_DimID( hid, 'hours', dimid, status )
            IF_NOTOK_RETURN(status=1)
            call MDF_Inquire_Dimension( hid, dimid, status, length=timeprof%nhour )
            IF_NOTOK_RETURN(status=1)                    

            if ( .not. allocated(timeprof%profile) ) then
              ! allocate timeprofile with dimension hour per year
              allocate( timeprof%profile(min_year:max_year,timeprof%nhour,timeprof%ncat,timeprof%ncountry) )
              ! by default fill with -999.0 
              timeprof%profile = -999.0
            end if
            

            ! 
            ! read country codes:
            allocate( countries_time_prof(ncountry) )
            call MDF_Get_StrArr( hid, 'zone'  , len_country_code, ncountry, countries_time_prof, status )
            IF_NOTOK_RETURN(status=1)          
            
            ! read emission category codes (LE SNAP):
            allocate( cat_codes_time_prof(ncat) )
            call MDF_Get_StrArr( hid, 'emis_cat'  , len_emiscatname, ncat, cat_codes_time_prof, status )
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

            ! loop over input categories:
            do icat = 1, ncat
              ! loop over input countries:
              do icountry = 1, ncountry
                ! store following order used in LE:
                timeprof%profile(iyear,:,cat_match(icat),country_match(icountry)) = profile_help(icountry,icat,:)
                ! check ...
                if ( abs( sum(profile_help(icountry,icat,:))/timeprof%nhour - 1.0 )  >= 1e-2)  then
                  write(gol, '("Time profile do not have average 1 ")'); call goErr            
                  write(gol, '("Country: ", a)') trim(countries_time_prof(icountry)); call goErr
                  !write(gol, '("Cat code: ", i2)' ) cat_codes_time_prof(icat)  ; call goErr
                  write(gol, '("Cat code: ",i0," `",a,"`")' ) icat, trim(cat_codes_time_prof(icat))  ; call goErr
                  write(gol, '("Average = ", f10.5)' ) sum(profile_help(icountry,icat,:))/timeprof%nhour   ; call goErr       
                  write(gol,*) 'test: ', sum(profile_help(icountry,icat,:))/timeprof%nhour , &
                                abs( sum(profile_help(icountry,icat,:))/timeprof%nhour - 1.0 ), 1e-2, &
                                abs( sum(profile_help(icountry,icat,:))/timeprof%nhour - 1.0 )  >= 1e-2; call goErr                
                  TRACEBACK; status=1; return
                 end if 
              end do
            end do
            
            ! check ...
            if ( any(timeprof%profile < 0.0) ) then
              write (gol,'("found undefined time profiles ...")'); call goErr
              do icat = 1, timeprof%ncat
                if ( any(timeprof%profile(:,:,icat,:) < 0.0) ) then
                  write (gol,'("  cat : ",i0," `",a,"`")') icat, trim(cat_codes(icat)); call goErr
                end if
              end do
              TRACEBACK; status=1; return
            end if
            
            ! clear          
            deallocate( countries_time_prof )
            deallocate( cat_codes_time_prof )
            deallocate( country_match )
            deallocate( cat_match )
            deallocate( profile_help )
            
            call MDF_close( hid, status )
            IF_NOTOK_RETURN(status=1)
            
          
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          case default
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            write (gol,'("unsupported format : ",a)') trim(fileformat); call goErr
            TRACEBACK; status=1; return

        end select

      end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

end module LE_Emis_EDGAR_Time_Prof

