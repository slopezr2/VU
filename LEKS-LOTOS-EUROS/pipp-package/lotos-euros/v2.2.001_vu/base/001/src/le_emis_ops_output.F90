!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_OPS_Output

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -------------------------------

  private
  public     :: LE_Emis_OPS_Write_Annual_Total
  
  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_OPS_Output'
  
contains

  subroutine LE_Emis_OPS_Write_Annual_Total( npol, polnames, emis, emis_units, basename, status )
  
    !
    use NetCDF
    use dims, only : nx,ny
    use LE_grid, only : ugg
    
    ! ---in/out---
    integer, intent(in)             ::  npol
    character(len=*), intent(in)    ::  polnames(npol)
    real, intent(in)                ::  emis(npol,nx,ny)
    character(len=*), intent(in)    ::  emis_units
    character(len=*), intent(in)    ::  basename
    integer, intent(out)            ::  status
    
    ! --- const -----
    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_Write_Annual_Total'
    
    ! local
    integer                   ::  ipol
    integer                   ::  maxlen_polname, polname_len
    
    character(len=1024)       ::  fname
    
    integer                   ::  ncid
    integer                   ::  dimid_lon, dimid_lat, dimid_pol, dimid_polname_len
    integer                   ::  varid_lon, varid_lat, varid_pol, varid_emis

    ! begin
    
    ! maximum lenght of pollutant names
    maxlen_polname = 0
    do ipol = 1, npol
      maxlen_polname = max( maxlen_polname, len(trim(polnames(ipol))) )
    end do

    ! create file
    write (fname,'(a,"_gridded.nc")') trim(basename)

    fname = trim(basename)//'_gridded.nc'
    status = NF90_Create( trim(fname), NF90_CLOBBER, ncid )
    if ( status /= NF90_NOERR ) then
      write( gol, '("creating file : ")'); call goErr
      write (gol,'("  file name  : ",a)') trim(fname); call goErr
      write (gol,'("  nf90 error : ",a)') trim(nf90_strerror(status)); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! define dimension
    status = NF90_Def_Dim( ncid, 'longitude', nx, dimid_lon )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( ncid, 'latitude', ny, dimid_lat )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( ncid, 'pollutant', npol, dimid_pol )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( ncid, 'pollutant_name_len', maxlen_polname, dimid_polname_len )
    
    ! define variables
    status = NF90_Def_Var( ncid, 'longitude', NF90_REAL, dimid_lon, varid_lon )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_lon, 'standard_name', 'longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_lon, 'long_name', 'longitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_lon, 'units', 'degrees_east' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_lon, '_CoordinateAxisType', 'Lon' )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Def_Var( ncid, 'latitude', NF90_REAL, dimid_lat, varid_lat )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_lat, 'standard_name', 'latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_lat, 'long_name', 'latitude' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_lat, 'units', 'degrees_north' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_lat, '_CoordinateAxisType', 'Lat' )
    IF_NF90_NOTOK_RETURN(status=1)

    status = NF90_Def_Var( ncid, 'pollutants', NF90_CHAR, (/dimid_polname_len,dimid_pol/), varid_pol )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_pol, 'long_name', 'Emitted pollutants' )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Def_Var( ncid, 'emission', NF90_REAL, (/dimid_pol, dimid_lon, dimid_lat/), varid_emis )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_emis, 'long_name', 'Annual emissions' )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid_emis, 'units', emis_units )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_EndDef( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! fill variables
    status = NF90_Put_Var( ncid, varid_lon, ugg%longitude_1d )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Var( ncid, varid_lat, ugg%latitude_1d)
    IF_NF90_NOTOK_RETURN(status=1)
    do ipol = 1, npol 
      polname_len = len(trim(polnames(ipol)))
      status = NF90_Put_Var( ncid, varid_pol, trim(polnames(ipol)), start=(/1,ipol/), count=(/polname_len,1/) ) 
    end do    
    status = NF90_Put_Var( ncid, varid_emis, emis )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! close file
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_OPS_Write_Annual_Total  
   
end module LE_Emis_OPS_Output

