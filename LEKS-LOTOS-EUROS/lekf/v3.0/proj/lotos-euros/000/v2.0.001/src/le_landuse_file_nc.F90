!###############################################################################
!
! Landuse_File_nc  -  Read NetCDF landuse/forest file
!
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

module LE_LandUse_File_nc

  use GO, only : gol, goPr, goErr
  
  implicit none


  ! --- in/out -------------------------------------

  private
  
  public  ::  T_Landuse_File_nc
  public  ::  Init, Done
  public  ::  Get
  
  
  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'Landuse_File_nc'
  

  ! --- types --------------------------------------
  
  type T_Landuse_File_nc
    ! file name:
    character(len=512)    ::  fname
    ! netcdf id
    integer               ::  ncid
  end type T_Landuse_File_nc
  
  
  ! --- interfaces ---------------------------------
  
  interface Init
    module procedure lu_Init
  end interface
    
  interface Done
    module procedure lu_Done
  end interface
  
  interface Get
    module procedure lu_Get
  end interface
  

contains
  
  
  ! ==================================================================
  
  
  subroutine lu_Init( lu, fname, status )
  
    use MDF, only : MDF_Open, MDF_NETCDF, MDF_READ

    ! --- in/out -----------------------------
    
    type(T_Landuse_File_nc), intent(out)   ::  lu
    character(len=*), intent(in)        ::  fname
    integer, intent(out)                ::  status
    
    ! --- const ------------------------------
    
    character(len=*), parameter ::  rname = mname//'/lu_Init'
    
    ! --- local ------------------------------
    
    logical       ::  exist
    
    ! --- begin ------------------------------
    
    ! store:
    lu%fname = trim(fname)
    
    ! check ...
    inquire( file=trim(lu%fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found:")'); call goErr
      write (gol,'("  ",a)') trim(lu%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! open file:
    call MDF_Open( trim(lu%fname), MDF_NETCDF, MDF_READ, lu%ncid, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine lu_Init
  
  
  ! ***
  
  
  subroutine lu_Done( lu, status )
  
    use MDF, only : MDF_Close

    ! --- in/out -----------------------------
    
    type(T_Landuse_File_nc), intent(inout)   ::  lu
    integer, intent(out)                  ::  status
    
    ! --- const ------------------------------
    
    character(len=*), parameter ::  rname = mname//'/lu_Done'
    
    ! --- local ------------------------------
    
    ! --- begin ------------------------------

    ! close file:
    call MDF_Close( lu%ncid, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine lu_Done
  
  
   ! ***
  
  
  subroutine lu_Get( lu, status, &
                     west_bound, east_bound, nlon, &
                     south_bound, north_bound, nlat, &
                     nlanduse, landuse, &
                     vegetation, glc, wwb )
  
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var
    use MDF, only : MDF_Get_Att

    ! --- in/out -----------------------------
    
    type(T_Landuse_File_nc), intent(inout)  ::  lu
    integer, intent(out)                    ::  status
    real, intent(out), optional             ::  west_bound, east_bound
    real, intent(out), optional             ::  south_bound, north_bound
    integer, intent(out), optional          ::  nlon, nlat
    integer, intent(out), optional          ::  nlanduse
    integer, intent(out), optional          ::  landuse(:,:)
    integer, intent(out), optional          ::  vegetation(:,:)
    integer, intent(out), optional          ::  glc(:,:)
    integer, intent(out), optional          ::  wwb(:,:)
   
    ! --- const ------------------------------
    
    character(len=*), parameter ::  rname = mname//'/lu_Get'
    
    ! --- local ------------------------------
    
    integer                   ::  dimid
    integer                   ::  varid
    integer(1), allocatable   ::  landuse_byte(:,:)
    integer(1), allocatable   ::  wwb_byte(:,:)
    integer(4), allocatable   ::  vegetation_in(:,:)
    
    ! --- begin ------------------------------

    ! extract data ?
    if ( present(west_bound) ) then
      ! get variable id:
      call MDF_Inq_VarID( lu%ncid, 'longitudes' , varid , status, quiet=.true. )
      if (status /= 0 ) then
        call MDF_Inq_VarID( lu%ncid, 'lon' , varid, status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! extract attribute value:
      call MDF_Get_Att( lu%ncid, varid, 'west_bound', west_bound, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! extract data ?
    if ( present(east_bound) ) then
      ! get variable id:
      call MDF_Inq_VarID( lu%ncid, 'longitudes' , varid , status, quiet=.true. )
      if ( status /= 0 ) then
        call MDF_Inq_VarID( lu%ncid, 'lon' , varid, status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! extract attribute value:
      call MDF_Get_Att( lu%ncid, varid, 'east_bound', east_bound, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! extract data ?
    if ( present(nlon) ) then
      ! get variable id:landuse
      call MDF_Inq_DimID( lu%ncid, 'longitude' , dimid , status, quiet=.true. )
      if ( status /= 0 ) then
        call MDF_Inq_DimID( lu%ncid, 'lon' , dimid, status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! extract dimension length:
      call MDF_Inquire_Dimension( lu%ncid, dimid, status, length=nlon )
      IF_NOTOK_RETURN(status=1)
    end if

    ! extract data ?
    if ( present(south_bound) ) then
      ! get variable id:
      call MDF_Inq_VarID( lu%ncid, 'latitudes' , varid , status, quiet=.true. )
      if ( status /= 0 ) then
        call MDF_Inq_VarID( lu%ncid, 'lat' , varid, status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! extract attribute value:
      call MDF_Get_Att( lu%ncid, varid, 'south_bound', south_bound, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! extract data ?
    if ( present(north_bound) ) then
      ! get variable id:
      call MDF_Inq_VarID( lu%ncid, 'latitudes' , varid , status, quiet=.true. )
      if ( status /= 0 ) then
        call MDF_Inq_VarID( lu%ncid, 'lat' , varid, status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! extract attribute value:
      call MDF_Get_Att( lu%ncid, varid, 'north_bound', north_bound, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! extract data ?
    if ( present(nlat) ) then
      ! get variable id:
      call MDF_Inq_DimID( lu%ncid, 'latitude' , dimid , status, quiet=.true. )
      if ( status /= 0 ) then
        call MDF_Inq_DimID( lu%ncid, 'lat' , dimid, status )
        if ( status /= 0 ) then
          call MDF_Inq_DimID( lu%ncid, 'latitudes' , dimid, status )
        endif
        IF_NOTOK_RETURN(status=1)
      end if
      ! extract dimension length:
      call MDF_Inquire_Dimension( lu%ncid, dimid, status, length=nlat )
      IF_NOTOK_RETURN(status=1)
    end if

    ! extract data ?
    if ( present(nlanduse) ) then
      ! get variable id:
      call MDF_Inq_DimID( lu%ncid, 'landuse_type' , dimid , status )
      IF_NOTOK_RETURN(status=1)
      ! extract dimension length:
      call MDF_Inquire_Dimension( lu%ncid, dimid, status, length=nlanduse )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! extract data ?
    if ( present(landuse) ) then

      ! get variable id:
      call MDF_Inq_VarID( lu%ncid, 'landuse' , varid , status )
      IF_NOTOK_RETURN(status=1)

      !! read data:
      !call MDF_Get_Var( lu%ncid, varid, landuse, status )
      !IF_NOTOK_RETURN(status=1)

      ! allocate byte array with same dimensions as "landuse" ;
      ! necessary, since directly reading into integer array
      ! failed on some platforms:
      allocate( landuse_byte(size(landuse,1),size(landuse,2)) )
      ! read data:
      call MDF_Get_Var( lu%ncid, varid, landuse_byte, status )
      IF_NOTOK_RETURN(status=1)
      ! convert:
      landuse = int(landuse_byte)
      ! clear:
      deallocate( landuse_byte )
      
    end if
    
    ! extract data ?
    if ( present(glc) ) then

      ! get variable id:
      call MDF_Inq_VarID( lu%ncid, 'glc' , varid, status )
      IF_NOTOK_RETURN(status=1)

      ! allocate byte array with same dimensions as "landuse" ;
      ! necessary, since directly reading into integer array
      ! failed on some platforms:
      allocate( landuse_byte(size(glc,1),size(glc,2)) )
      ! read data:
      call MDF_Get_Var( lu%ncid, varid, landuse_byte, status )
      IF_NOTOK_RETURN(status=1)
      ! convert:
      glc = int(landuse_byte)
      ! clear:
      deallocate( landuse_byte )
      
    end if
   
    ! extract data ?
    if ( present(vegetation) ) then

      ! get variable id:
      call MDF_Inq_VarID( lu%ncid, 'vegetation' , varid , status )
      IF_NOTOK_RETURN(status=1)

      !! read data:
      !call MDF_Get_Var( lu%ncid, varid, forest, status )
      !IF_NOTOK_RETURN(status=1)

      ! allocate byte array with same dimensions as "forest" ;
      ! necessary, since directly reading into integer array
      ! failed on some platforms:
      allocate( vegetation_in(size(vegetation,1),size(vegetation,2)), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! read data:
      call MDF_Get_Var( lu%ncid, varid, vegetation_in, status )
      IF_NOTOK_RETURN(status=1)
      ! convert:
      vegetation = int(vegetation_in)
      ! clear:
      deallocate( vegetation_in, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if
    
    ! extract data ?
    if ( present(wwb) ) then

      ! get variable id:
      call MDF_Inq_VarID( lu%ncid, 'wwb' , varid, status )
      IF_NOTOK_RETURN(status=1)

      ! allocate byte array with same dimensions as "landuse" ;
      ! necessary, since directly reading into integer array
      ! failed on some platforms:
      allocate( wwb_byte(size(wwb,1),size(wwb,2)) )
      ! read data:
      call MDF_Get_Var( lu%ncid, varid, wwb_byte, status )
      IF_NOTOK_RETURN(status=1)
      ! convert:
      wwb = int(wwb_byte)
      ! clear:
      deallocate( wwb_byte )
      
      
    end if
    
    ! ok
    status = 0
  
  end subroutine lu_Get
  
  
 
end module LE_LandUse_File_nc

