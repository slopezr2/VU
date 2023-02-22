!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################

module LE_Emis_Time_Prof_Grid

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Time_Prof_Grid
  public  ::  LE_Emis_Time_Prof_Grid_Init, LE_Emis_Time_Prof_Grid_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Time_Prof_Grid'
  
  
  ! --- types --------------------------------
  
  ! storage for emission composition:
  type T_Emis_Time_Prof_Grid
    ! dimensions:
    integer                         ::  ncat
    character(len=8),allocatable    ::  icats(:)  ! for which categories is this valid??
    integer                         ::  nhour  ! (8760 for a NON LEAP year)
    ! fraction assigned to component:
    real, allocatable             ::  profile(:,:,:,:)  ! (hour,sector,nlon,nlat)
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_Time_Prof_Grid_Init( timeprof, query, cat_codes, &
                                     year, status )

    use GO,  only :  goVarValue
    use GO,  only :  days_in_year
    use MDF, only : MDF_Open, MDF_Close
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var
    use MDF, only : MDF_NETCDF, MDF_READ
    use LE_grid       , only : ugg
    use C3PO          , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors

    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof_Grid), intent(out)  ::  timeprof
    character(len=*), intent(in)              ::  query
    character(len=*), intent(in)              ::  cat_codes(:)
    integer, intent(in)                       ::  year
    integer, intent(out)                      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Grid_Init'
    
    ! --- local -------------------------------
        
    character(len=512)              ::  filename
    character(len=32)               ::  fileformat
    logical                         ::  exist
    
    integer                         ::  hid
    integer                         ::  dimid
    integer                         ::  varid
    
    integer                         ::  icat
    integer                        ::  ilon,ilat, ih
    character(len=8), allocatable   ::  cat_codes_time_prof(:)
    character(len=8)                ::  cat_code
  
    character(len=1024)            ::  description
    type(T_File_Ugg)               ::  file_in
    type(T_Grid_Ugg)               ::  grid_in
    integer                        ::  start_ind(2)
    integer                        ::  count_ind(2)
    integer                        ::  start_ind_4d(4)
    integer                        ::  count_ind_4d(4)
    real, allocatable              ::  values_in(:,:)
    character(len=32)              ::  units_in
    

    ! --- begin -------------------------------
    
    ! extract filename, raise error message if not defined (status<0):
    call goVarValue( trim(query), ';', 'file', '=', filename, status )
    if ( status < 0 ) then
      write (gol,'("no keyword `file` in time profile query : ",a)') trim(query); call goErr
      TRACEBACK; status=1; return
    else if ( status > 0 ) then
      TRACEBACK; status=1; return
    end if
        
    ! format description:
    fileformat = 'nc'
    call goVarValue( trim(query), ';', 'format', '=', fileformat, status )
    IF_ERROR_RETURN(status=1)
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( len_trim(filename) == 0 ) then   ! no file, fill with default value (-999.0)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ntime_vars number of days
      timeprof%nhour = days_in_year( year ) * 24
      
      allocate( timeprof%profile(timeprof%nhour,timeprof%ncat,ugg%nlon,ugg%nlat ) )
      ! by default fill with -999.0
      timeprof%profile = -999.0
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else  ! read table from file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

      ! Check existance
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

          ! Read in valid categories for this profile
          call MDF_Open( trim(filename), MDF_NETCDF, MDF_READ, hid, status )
          IF_NOTOK_RETURN(status=1)
          
          ! number of SNAP categories in time profile file
          call MDF_Inq_DimID( hid, 'emis_cat', dimid, status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Inquire_Dimension( hid, dimid, status, length=timeprof%ncat )
          IF_NOTOK_RETURN(status=1)
          ! number of hours in time profile file
          call MDF_Inq_DimID( hid, 'hour', dimid, status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Inquire_Dimension( hid, dimid, status, length=timeprof%nhour )
          IF_NOTOK_RETURN(status=1)                    

          ! allocate target array
          if ( .not. allocated(timeprof%profile) ) then
            ! allocate timeprofile with dimension hour per year
            allocate( timeprof%profile(timeprof%nhour,timeprof%ncat,ugg%nlon,ugg%nlat) )
            ! by default fill with -999.0 
            timeprof%profile = -999.0
          end if

          ! read emission category codes (LE SNAP):
          allocate( cat_codes_time_prof(timeprof%ncat), stat=status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Inq_VarID( hid, 'emis_cat', varid, status )
          IF_NOTOK_RETURN(status=1)
          call MDF_Get_Var( hid, varid, cat_codes_time_prof, status )
          IF_NOTOK_RETURN(status=1)

          ! category indices:
          allocate( timeprof%icats( timeprof%ncat), stat=status )
          IF_NOTOK_RETURN(status=1)
          ! fill:
          timeprof%icats=cat_codes_time_prof          

          ! close:
          call MDF_close( hid, status )
          IF_NOTOK_RETURN(status=1)           
          !
          ! Use c3po tools to read data and convert to LE-grid
          ! 
          
          ! open file:
          call file_in%Open( trim(filename), status )
          IF_NOTOK_RETURN(status=1)

          ! Hard-coded description --> in future bring to rc-file
          description = 'long_name=Hourly factor of emission based on heating data, grid cell dependent'
          
          ! variable id:
          call file_in%Inq_VarID( trim(description), varid, status )
          IF_NOTOK_RETURN(status=1)
          ! grid definition:
          call file_in%Get_Grid( varid, grid_in, status, &
                                 ugg_to=ugg, start_ind=start_ind, count_ind=count_ind )
          
          ! storage:
          allocate( values_in( grid_in%nlon, grid_in%nlat), stat=status )
          IF_NOTOK_RETURN(status=1)

          ! loop over categories and times in file
          do icat = 1, timeprof%ncat
            do ih = 1, timeprof%nhour
            
            ! indices in 4d array
            start_ind_4d = (/start_ind(1),start_ind(2),icat,ih/)
            count_ind_4d = (/count_ind(1),count_ind(2),1,1/)
            ! read:
            call file_in%Get_Var( trim(description), values_in, units_in, status, &
                                  start=start_ind_4d, count=count_ind_4d )

            ! regrid to LE-grid
            call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, timeprof%profile(ih,icat,:,:), status )
            IF_NOTOK_RETURN(status=1)
            end do
          end do
          
          ! clear:
          deallocate( values_in, stat=status )
          IF_NOTOK_RETURN(status=1)
          call grid_in%Done(status)
          IF_NOTOK_RETURN(status=1)
          ! close:
          call file_in%Close(status)
          IF_NOTOK_RETURN(status=1)

          ! Check profiles on average=1
          ! loop over indices in file:
          do icat = 1, timeprof%ncat
            do ilon = 1, ugg%nlon
            do ilat = 1, ugg%nlat
              ! check ...
              if ( abs(sum(timeprof%profile(:,icat,ilon,ilat))/timeprof%nhour -1 )  >= 1e-2)  then
                write(gol, '("Time profile do not have average 1 ")'); call goErr            
                write(gol, '("Cat code: ", a)' ) trim(cat_codes_time_prof(icat)); call goErr
                write(gol, '("lon/lat indices: ", 2i4)' ) ilon,ilat ; call goErr
                write(gol, '("Average = ", f10.5)' ) sum(timeprof%profile(:,icat,ilon,ilat))/timeprof%nhour   ; call goErr                  
                TRACEBACK; status=1; return
               end if 
            end do
            end do
          end do


        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          write (gol,'("unsupported format : ",a)') trim(fileformat); call goErr
          TRACEBACK; status=1; return

      end select

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! ok
    status = 0

  end subroutine LE_Emis_Time_Prof_Grid_Init
  
  
  ! ***
  

  subroutine LE_Emis_Time_Prof_Grid_Done( timeprof, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof_Grid), intent(inout) ::  timeprof
    integer, intent(out)                  ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Grid_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    if ( allocated(timeprof%profile     ) ) deallocate( timeprof%profile    )

    ! ok
    status = 0

  end subroutine LE_Emis_Time_Prof_Grid_Done

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

end module LE_Emis_Time_Prof_Grid

