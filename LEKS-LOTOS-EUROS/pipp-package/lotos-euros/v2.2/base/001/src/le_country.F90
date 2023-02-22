!#################################################################
!
! NAME
!   LE_Country  -  country map tools
!
! VARIABLES
!
!   integer             ::  country_map_n                ! number of countries
!   character(len=3)    ::  country_map_code(country_map_n)  ! 'ALB', 'BEL', ...
!   real                ::  country_map_frac(nx,ny,country_map_n)  ! [0-1]
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



module LE_Country

  use GO, only : gol, goPr, goErr
      
  implicit none
  
  
  ! --- in/out --------------------------
  
  private
  
  public  ::  country_map_n, country_map_code, country_map_frac
  public  ::  LE_Country_Init, LE_Country_Done
  public  ::  read_countries
  

  ! --- const ---------------------------

  character(len=*), parameter   ::  mname = 'LE_Country'
  
  
  ! --- variables -----------------------
   
  integer                         ::  country_map_n
  character(len=3), allocatable   ::  country_map_code(:)       ! (country_map_n)
  real, allocatable               ::  country_map_frac(:,:,:)   ! (nlon,nlat,country_map_n)
  
  ! read countries?
  logical                         ::  read_countries = .false.

contains


  ! =========================================================  

  
  subroutine LE_Country_Init( rcF, status )
  
    use GO            , only : TrcFile, ReadRc
    use MDF           , only : MDF_Open, MDF_Close
    use MDF           , only : MDF_NETCDF, MDF_READ
    use MDF           , only : MDF_Inq_Dimid, MDF_Inquire_Dimension
    use LE_grid       , only : ugg
    use C3PO          , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors
    
    ! --- in/out ----------------------------------
       
    type(TRcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status
       
    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Country_Init'
    
    ! --- local -----------------------------------

    character(len=512)        ::  fname
    logical                   ::  exist
    integer                   ::  ncid, dimid, varid
    integer                   ::  icountry, country_code_len
    character(len=32)         ::  country
        
    character(len=1024)                 ::  description
    type(T_File_Ugg)                    ::  file_in
    type(T_Grid_Ugg)                    ::  grid_in
    integer                             ::  start_ind(2)
    integer                             ::  count_ind(2)
    real, allocatable                   ::  values_in(:,:)
    character(len=32)                   ::  country_units
    
    ! --- begin -----------------------------------
    
    
    ! check:
    if ( read_countries ) then
      
      write( gol, '("Read general coutnry fractions")' ) ; call goPr
      
      ! Read country file
      call ReadRc( rcF, 'le.country.fraction_file', fname, status )
      IF_NOTOK_RETURN(status=1)

      ! check ...
      inquire( file=trim(fname), exist=exist )
      if ( .not. exist ) then
        write( gol, '("file not found:")') ; call goErr
        write( gol, '("  ", a)') trim(fname) ; call goErr
        TRACEBACK;status=1;return
      end if

      ! Read country codes
      ! open:
      call MDF_Open( trim(fname), MDF_NETCDF, MDF_READ, ncid, status )
      IF_NOTOK_RETURN(status=1)

      ! country list
      call MDF_Inq_Dimid( ncid, 'country', dimid, status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Inquire_Dimension( ncid, dimid, status, length=country_map_n )
      IF_NOTOK_RETURN(status=1)
      call MDF_Inq_Dimid( ncid, 'country_code_len', dimid, status )
      IF_NOTOK_RETURN(status=1)
      call MDF_Inquire_Dimension( ncid, dimid, status, length=country_code_len )
      IF_NOTOK_RETURN(status=1)

      ! allocate
      allocate( country_map_frac(ugg%nlon,ugg%nlat,country_map_n), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( country_map_code(country_map_n), stat=status )
      IF_NOTOK_RETURN(status=1)

      ! read country names:
      call MDF_Get_StrArr( ncid, 'country', country_code_len, country_map_n, country_map_code, status )
      IF_NOTOK_RETURN(status=1)

      ! close file:
      call MDF_Close( ncid,  status )
      IF_NOTOK_RETURN(status=1)


      ! Open file again in structure with quick interpolation to target grid:
      call file_in%Open( trim(fname), status )
      IF_NOTOK_RETURN(status=1)

      ! loop over countries to find fractions per country
      do icountry = 1, country_map_n

        country = trim(country_map_code(icountry) )

        ! description in netCDF file
        description = 'var_name='//trim(country)

        ! variable id:
        call file_in%Inq_VarID( trim(description), varid, status )
        IF_NOTOK_RETURN(status=1)
        ! grid definition:
        call file_in%Get_Grid( varid, grid_in, status, &
                              ugg_to=ugg, start_ind=start_ind, count_ind=count_ind )
        IF_NOTOK_RETURN(status=1)

        ! storage
        allocate( values_in(grid_in%nlon, grid_in%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! read :
        call file_in%Get_Var( trim(description), values_in, country_units, status, &
                              start=start_ind, count=count_ind )                              
        IF_NOTOK_RETURN(status=1)

        ! regrid to LE-grid
        call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, country_map_frac(:,:,icountry), status )
        IF_NOTOK_RETURN(status=1)

        ! clear:
        deallocate( values_in, stat=status )
        IF_NOTOK_RETURN(status=1)
        call grid_in%Done(status)
        IF_NOTOK_RETURN(status=1)

      end do  ! countries

      ! close file
      call file_in%Close(status)
      IF_NOTOK_RETURN(status=1)
  
    else
      write( gol, '("No need to read general country fractions")' ) ; call goPr

    end if  ! check
                         
    ! ok
    status = 0
      
  end subroutine LE_Country_Init
  

  ! =========================================================  
  

  subroutine LE_Country_Done( status )
  
    ! --- in/out --------------------------
    
    integer, intent(out)          ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Country_Done'
    
    ! --- local ----------------------------
    
    
    ! --- begin ----------------------------
    
    if ( read_countries ) then
      ! clear:
      deallocate( country_map_code, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( country_map_frac, stat=status )
      IF_NOTOK_RETURN(status=1)            

    end if
        
    ! ok
    status = 0
    
  end subroutine LE_country_Done


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

end module LE_Country
