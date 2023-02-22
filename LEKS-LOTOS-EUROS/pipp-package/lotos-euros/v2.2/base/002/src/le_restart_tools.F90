!#######################################################################
!
! LE_Restart_Tools
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

module LE_Restart_Tools

  use GO, only : gol, goPr, goErr
#ifdef with_netcdf
  use NetCDF, only : NF90_NOERR, nf90_strerror
#endif

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Restart_FileName
  public  ::  LE_Restart_Restore


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Restart_Tools'


  ! --- interfaces -------------------------------
  
  interface LE_Restart_Restore
    module procedure LE_Restart_Restore_1D
    module procedure LE_Restart_Restore_2D
    module procedure LE_Restart_Restore_3D
    module procedure LE_Restart_Restore_4D
  end interface LE_Restart_Restore


contains



  ! ====================================================================
  

  subroutine LE_Restart_FileName( fname, key, t, status )
  
    use GO, only : TDate, Get
    use GO, only : goVarValue
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(out)     ::  fname
    character(len=*), intent(in)      ::  key
    type(TDate), intent(in)           ::  t
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_FileName'
    
    ! --- local ----------------------------------
    
    character(len=32)     ::  model
    character(len=32)     ::  expid
    character(len=32)     ::  name
    integer               ::  time6(6)
  
    ! --- begin ----------------------------------
    
    ! extract parts from key:
    !
    !   model=LE;expid=base;name=conc
    !
    model = 'LE'
      call goVarValue( key, ';', 'model', '=', model, status )
      IF_ERROR_RETURN(status=1)
    expid = 'base'
      call goVarValue( key, ';', 'expid', '=', expid, status )
      IF_ERROR_RETURN(status=1)
    name = 'conc'
      call goVarValue( key, ';', 'name' , '=', name , status )
      IF_ERROR_RETURN(status=1)

    ! extract time:
    call Get( t, time6=time6 )
    
    ! file name: <model>_<runid>_<name>_20070101_0000.nc
    write (fname,'(a,2("_",a),"_",i4,2i2.2,"_",2i2.2,".nc")') &
              trim(model), trim(expid), trim(name), time6(1:5)
              
    ! ok
    status = 0
    
  end subroutine LE_Restart_FileName
  
  

  
  
  ! *
  
 
  !
  ! return status:
  !   -1  : variable not found in file
  !    0  : ok
  !  other : error
  !
  
  subroutine LE_Restart_Restore_1D( varname, values, t, path, key, status, cnt )
  
    use GO  , only : TDate
#ifdef with_netcdf
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
#endif
    
    ! --- in/out ----------------------------------
    
    character(len=*), intent(in)    ::  varname
    real, intent(out)               ::  values(:)
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status
    
    integer, intent(out), optional  ::  cnt
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Restore_1D'
    
    ! --- local -----------------------------------
    
    character(len=256)    ::  fname
    logical               ::  exist
    
    integer       ::  ncid
    integer       ::  varid
    integer       ::  dimids(1)
    integer       ::  cmode
    
    ! --- begin -----------------------------------
    
    ! check implementation:    
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read boundary conditions from LE NetCDF output.")'); call goErr
    write (gol,'("Define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! format file name:
    call LE_Restart_FileName( fname, key, t, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    inquire( file=trim(path)//'/'//trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("restart file not found : ")'); call goErr
      write (gol,'("  ",a)') trim(path)//'/'//trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
#ifdef with_netcdf
    ! set open mode flag:
    cmode = NF90_NOWRITE   ! read-only

    ! open file:
    status = NF90_Open( trim(path)//'/'//trim(fname), cmode, ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! search variable, return id:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=-1)
    
    ! return size ?
    if ( present(cnt) ) then
      ! obtain dimension id:
      status = nf90_inquire_variable( ncid, varid, dimids=dimids )
      IF_NF90_NOTOK_RETURN(status=1)
      ! get length:
      status = nf90_inquire_dimension( ncid, dimids(1), len=cnt )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    
    ! read data:
    status = NF90_Get_Var( ncid, varid, values )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_1D
  

  ! *
  

  subroutine LE_Restart_Restore_2D( varname, values, t, path, key, status )
  
    use GO     , only : TDate
    use LE_Grid, only : dom
#ifdef with_netcdf
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
#endif
    
    ! --- in/out ----------------------------------
    
    character(len=*), intent(in)    ::  varname
    real, intent(out)               ::  values(:,:)
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Restore_2D'
    
    ! --- local -----------------------------------
    
    character(len=256)    ::  fname
    logical               ::  exist
    integer               ::  ncid
    integer               ::  varid
    integer               ::  cmode
    integer               ::  glbo(2)
    
    ! --- begin -----------------------------------
    
    ! check implementation:    
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read boundary conditions from LE NetCDF output.")'); call goErr
    write (gol,'("Define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! domain start in global index space:
    call dom%Get( status, glbo=glbo )
    IF_NOTOK_RETURN(status=1)

    ! format file name:
    call LE_Restart_FileName( fname, key, t, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    inquire( file=trim(path)//'/'//trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("restart file not found : ")'); call goErr
      write (gol,'("  ",a)') trim(path)//'/'//trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
#ifdef with_netcdf
    ! set open mode flag:
    cmode = NF90_NOWRITE   ! read-only

    ! open file:
    status = NF90_Open( trim(path)//'/'//trim(fname), cmode, ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read data:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, values, &
                             start=(/glbo(1),glbo(2)/), count=shape(values) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_2D
  

  ! *
  

  subroutine LE_Restart_Restore_3D( varname, values, t, path, key, status )
  
    use GO     , only : TDate
    use LE_Grid, only : dom
#ifdef with_netcdf
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
#endif
    
    ! --- in/out ----------------------------------
    
    character(len=*), intent(in)    ::  varname
    real, intent(out)               ::  values(:,:,:)
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Restore_3D'
    
    ! --- local -----------------------------------
    
    character(len=256)    ::  fname
    logical               ::  exist    
    integer               ::  ncid
    integer               ::  varid
    integer               ::  cmode
    integer               ::  glbo(2)
    
    ! --- begin -----------------------------------
    
    ! check implementation:    
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read boundary conditions from LE NetCDF output.")'); call goErr
    write (gol,'("Define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! domain start in global index space:
    call dom%Get( status, glbo=glbo )
    IF_NOTOK_RETURN(status=1)

    ! format file name:
    call LE_Restart_FileName( fname, key, t, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    inquire( file=trim(path)//'/'//trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("restart file not found : ")'); call goErr
      write (gol,'("  ",a)') trim(path)//'/'//trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
#ifdef with_netcdf
    ! set open mode flag:
    cmode = NF90_NOWRITE   ! read-only

    ! open file:
    status = NF90_Open( trim(path)//'/'//trim(fname), cmode, ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read data:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, values, &
                             start=(/glbo(1),glbo(2),1/), count=shape(values) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_3D
  
  
  ! *
  

  subroutine LE_Restart_Restore_4D( varname, values, t, path, key, status )
  
    use GO     , only : TDate
    use LE_Grid, only : dom
#ifdef with_netcdf
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
#endif
    
    ! --- in/out ----------------------------------
    
    character(len=*), intent(in)    ::  varname
    real, intent(out)               ::  values(:,:,:,:)
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Restore_4D'
    
    ! --- local -----------------------------------
    
    character(len=256)    ::  fname
    logical               ::  exist
    integer               ::  ncid
    integer               ::  varid
    integer               ::  cmode
    integer               ::  glbo(2)
    
    ! --- begin -----------------------------------
    
    ! check implementation:    
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read boundary conditions from LE NetCDF output.")'); call goErr
    write (gol,'("Define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! domain start in global index space:
    call dom%Get( status, glbo=glbo )
    IF_NOTOK_RETURN(status=1)

    ! format file name:
    call LE_Restart_FileName( fname, key, t, status )
    IF_NOTOK_RETURN(status=1)

    ! check ...
    inquire( file=trim(path)//'/'//trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("restart file not found : ")'); call goErr
      write (gol,'("  ",a)') trim(path)//'/'//trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
#ifdef with_netcdf
    ! set open mode flag:
    cmode = NF90_NOWRITE   ! read-only

    ! open file:
    status = NF90_Open( trim(path)//'/'//trim(fname), cmode, ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read data:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, values, &
                             start=(/glbo(1),glbo(2),1,1/), count=shape(values) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_4D
  
  
end module LE_Restart_Tools

