!#######################################################################
!
! LE_Restart - save and restore LOTOS-EUROS state
!
! EXAMPLE
!
!   integer, parameter           ::  nx = 100, ny = 140, nz = 4, nspec = 38
!   character(len=*), parameter  ::  path = '/scratch/output/'
!   character(len=*), parameter  ::  key = 'model=LE;expid=base;name=conc'
!   
!   real                  ::  c(nx,ny,nz,nspec)
!   type(TDate)           ::  t
!   integer               ::  status
!
!   type(T_LE_Restart_File)   ::  F
!   integer                   ::  dimid_nnoise
!   integer                   ::  varid_dc
!
!   ! *
!
!   ! dump concentrations and auxilary data:
!   call LE_Restart_Save( c, t, path, key, status )
!
!   ! ... or use expert routines to add extra fields:
!
!   ! create file, define standard dimensions:
!   call LE_Restart_Create( F, t, path, key, status )
!   ! extra dimensions:
!   status = NF90_Def_Dim( F%ncid, 'nnoise', nnoise, dimid_nnoise )
!   ! extra variables:
!   status = NF90_Def_Var( F%ncid, 'dc', NF90_FLOAT, &
!               (/F%dimid_nx,F%dimid_ny,F%dimid_nnoise/), varid_dc )
!   
!   ! end definition, write standard fields:
!   call LE_Restart_Write( F, c, status )
!   ! write extra fields:
!   status = NF90_Put_Var( F%ncid, varid_dc, dc )
! 
!   ! close file:
!   call LE_Restart_Close( F, status )
!
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

module LE_Restart

  use GO, only : gol, goPr, goErr
#ifdef with_netcdf
  use NetCDF, only : NF90_NOERR, nf90_strerror
#endif
  use LE_Restart_Tools, only : LE_Restart_Filename
  use LE_Restart_Tools, only : LE_Restart_Restore

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_LE_Restart_File

  public  ::  LE_Restart_Save
  public  ::  LE_Restart_Create, LE_Restart_Write, LE_Restart_Close

  public  ::  LE_Restart_Restore_State
  public  ::  LE_Restart_Restore_Data
  
  public  ::  LE_Restart_FileName
  public  ::  LE_Restart_Restore


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Restart'


  ! --- types ------------------------------------

  ! interface to restart file:
  type T_LE_Restart_File
    integer       ::  ncid
    integer       ::  dimid_nx, dimid_ny, dimid_nz, dimid_nz_bnd, dimid_nspec
    integer       ::  varid_volume
    integer       ::  varid_airmass
    integer       ::  varid_hpres
    integer       ::  varid_expcls
    integer       ::  varid_aerh2o
    integer       ::  varid_c
    integer       ::  varid_cg
    integer       ::  varid_cnh3_sum
    integer       ::  varid_cso2_sum
    integer       ::  varid_cnh3_ave_prev
    integer       ::  varid_cso2_ave_prev
  end type T_LE_Restart_File


contains


  ! ********************************************************************
  ! ***
  ! *** save
  ! ***
  ! ********************************************************************
  
  subroutine LE_Restart_Save( c, cg, aerh2o, &
                              bud, t, path, key, status )
  
    use GO, only : TDate
    use dims, only : nx, ny, nz, nspec
    use LE_Budget, only : T_Budget
    
    ! --- in/out ----------------------------------
    
    real, intent(in)                ::  c(nx,ny,nz,nspec)
    real, intent(in)                ::  cg(nx,ny,nspec)
    real, intent(in)                ::  aerh2o(nx,ny,nz)        
    type(T_Budget), intent(in)      ::  bud
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Save'
    
    ! --- local -----------------------------------
    
    type(T_LE_Restart_File)   ::  F

    ! --- begin -----------------------------------
    
    ! create file, define standard dimensions:
    call LE_Restart_Create( F, t, path, key, status )
    IF_NOTOK_RETURN(status=1)
    
    ! end definition, write standard fields:
    call LE_Restart_Write( F, c, cg, aerh2o, &
                           bud, status )
    IF_NOTOK_RETURN(status=1)
    
    ! close file:
    call LE_Restart_Close( F, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Restart_Save
    

  ! *
  

  subroutine LE_Restart_Create( F, t, path, key, status )
  
#ifdef with_netcdf
    use NetCDF , only : NF90_CLOBBER, NF90_NOCLOBBER
    use NetCDF , only : NF90_FLOAT
    use NetCDF , only : NF90_Create
    use NetCDF , only : NF90_Def_Dim
    use NetCDF , only : NF90_Def_Var
#endif

    use GO     , only : goc
    use GO     , only : TDate
    use GO     , only : CheckDir, pathsep
    use LE_Grid, only : glb_ugg
    use Dims   , only : nz
    use Indices, only : nspec
    use LE_Restart_Tools, only : LE_Restart_Filename

    ! --- in/out ----------------------------------
    
    type(T_LE_Restart_File), intent(out)    ::  F
    type(TDate), intent(in)                 ::  t
    character(len=*), intent(in)            ::  path
    character(len=*), intent(in)            ::  key
    integer, intent(out)                    ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Create'
    
    ! --- local -----------------------------------
    
    character(len=256)        ::  fname
    integer                   ::  cmode

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

    ! info ...
    write (gol,'("LE:       save ",a)') trim(fname); call goPr
    
    ! root only:
    if ( goc%root ) then
    
      ! create directory if necessary:
      call CheckDir( trim(path)//pathsep//'.', status )
      IF_NOTOK_RETURN(status=1)
    
#ifdef with_netcdf
      ! set creation mode flag:
      !cmode = NF90_NOCLOBBER     ! do not overwrite existing files, raise error instead
      cmode = NF90_CLOBBER       ! overwrite existing files if necessary

      ! create file:
      status = NF90_Create(  trim(path)//'/'//trim(fname), cmode, F%ncid )
      IF_NF90_NOTOK_RETURN(status=1)

      ! create dimensions:
      status = NF90_Def_Dim( F%ncid, 'nx', glb_ugg%nlon, F%dimid_nx )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( F%ncid, 'ny', glb_ugg%nlat, F%dimid_ny )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( F%ncid, 'nz', nz, F%dimid_nz )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( F%ncid, 'nz_bnd', nz+1, F%dimid_nz_bnd )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( F%ncid, 'nspec', nspec, F%dimid_nspec )
      IF_NF90_NOTOK_RETURN(status=1)
    
      ! create variable:
      status = NF90_Def_Var( F%ncid, 'volume', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny,F%dimid_nz/), F%varid_volume )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'airmass', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny,F%dimid_nz/), F%varid_airmass )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'hpres', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny,F%dimid_nz_bnd/), F%varid_hpres )
      IF_NF90_NOTOK_RETURN(status=1)

      ! create variable:
      status = NF90_Def_Var( F%ncid, 'expcls', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny/), F%varid_expcls )
      IF_NF90_NOTOK_RETURN(status=1)

      ! create variable:
      status = NF90_Def_Var( F%ncid, 'aerh2o', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny,F%dimid_nz/), F%varid_aerh2o )
      IF_NF90_NOTOK_RETURN(status=1)

      ! create variable:
      status = NF90_Def_Var( F%ncid, 'c', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny,F%dimid_nz,F%dimid_nspec/), &
                             F%varid_c )
      IF_NF90_NOTOK_RETURN(status=1)
      ! create variable:
      status = NF90_Def_Var( F%ncid, 'cg', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny,F%dimid_nspec/), &
                             F%varid_cg )
      IF_NF90_NOTOK_RETURN(status=1)

      ! create variable:
      status = NF90_Def_Var( F%ncid, 'cnh3_sum', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny/), &
                             F%varid_cnh3_sum )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'cso2_sum', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny/), &
                             F%varid_cso2_sum )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'cnh3_ave_prev', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny/), &
                             F%varid_cnh3_ave_prev )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'cso2_ave_prev', NF90_FLOAT, &
                             (/F%dimid_nx,F%dimid_ny/), &
                             F%varid_cso2_ave_prev )
      IF_NF90_NOTOK_RETURN(status=1)

                          
#endif
    end if ! root

    ! ok
    status = 0
    
  end subroutine LE_Restart_Create
    

  ! *
  
  subroutine LE_Restart_Write( F, c, cg, aerh2o, &
                               bud, status )
 
#ifdef with_netcdf
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_EndDef
#endif
    use GO           , only : goc
    use Dims         , only : nx, ny, nz
    use LE_Grid      , only : dom
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Meteo_Data, only : volume
    use LE_Meteo_Data, only : airmass
    use LE_Meteo_Data, only : hpres
    use Dims         , only : expcls
    use Indices      , only : nspec
    use Indices      , only : i_nh3, i_so2
    use LE_Budget    , only : T_Budget

    ! --- in/out ----------------------------------
    
    type(T_LE_Restart_File), intent(inout)  ::  F
    real, intent(in)                        ::  c(nx,ny,nz,nspec)
    real, intent(in)                        ::  cg(nx,ny,nspec)
    real, intent(in)                        ::  aerh2o(nx,ny,nz)
    type(T_Budget), intent(in)              ::  bud
    integer, intent(out)                    ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Write'
    
    ! --- local -----------------------------------
    
    integer               ::  ival
    integer               ::  ispec
    real, allocatable     ::  pat(:,:)

    !real, pointer          ::  volume(:,:,:)   ! (lon,lat,nz)
    
    ! --- begin -----------------------------------
    
    !call LE_Data_GetPointer( 'vol', volume, status, check_units ='m3' )
    !IF_NOTOK_RETURN(status=1)    
    
    ! check implementation:    
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read boundary conditions from LE NetCDF output.")'); call goErr
    write (gol,'("Define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! root only:
    if ( goc%root ) then

      ! NH3 compensation point
      ival = 0
      if ( i_nh3 > 0 ) ival = bud%drydepos%cnh3_nsum
      ! write counter as attribute:
#ifdef with_netcdf
      status = NF90_Put_Att( F%ncid, F%varid_cnh3_sum, 'nsum', ival )
      IF_NF90_NOTOK_RETURN(status=1)
#endif

      ! SO2 for co-deposition
      ival = 0
      if ( i_so2 > 0 ) ival = bud%drydepos%cso2_nsum
      ! write counter as attribute:
#ifdef with_netcdf
      status = NF90_Put_Att( F%ncid, F%varid_cso2_sum, 'nsum', ival )
      IF_NF90_NOTOK_RETURN(status=1)
#endif

      ! end definition phase:
#ifdef with_netcdf
      status = NF90_EndDef( F%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#endif
      
    end if ! root
    
    ! write data:
    call dom%Put_Var( F%ncid, F%varid_volume , -1, volume(1:nx,1:ny,1:nz), status )
    IF_NOTOK_RETURN(status=1)
    call dom%Put_Var( F%ncid, F%varid_airmass, -1, airmass, status )
    IF_NOTOK_RETURN(status=1)
    call dom%Put_Var( F%ncid, F%varid_hpres  , -1, hpres, status )
    IF_NOTOK_RETURN(status=1)
    
    ! write data:
    call dom%Put_Var( F%ncid, F%varid_expcls , -1, expcls, status )
    IF_NOTOK_RETURN(status=1)
    
    ! write data:
    call dom%Put_Var( F%ncid, F%varid_aerh2o , -1, aerh2o, status )
    IF_NOTOK_RETURN(status=1)
    
    ! write c per tracer to limit memory:
    do ispec = 1, nspec
      ! write slab:
      call dom%Put_Var( F%ncid, F%varid_c, ispec, c(:,:,:,ispec), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! write data:
    call dom%Put_Var( F%ncid, F%varid_cg , -1, cg, status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for fields that might not be present:
    allocate( pat(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! init with dummy value:
    pat = 0.0
    
    ! summed concentrations:
    if ( i_nh3 > 0 ) pat = bud%drydepos%cnh3_sum
    ! write data:
    call dom%Put_Var( F%ncid, F%varid_cnh3_sum , -1, pat, status )
    IF_NOTOK_RETURN(status=1)

    ! summed concentrations:
    if ( i_so2 > 0 ) pat = bud%drydepos%cso2_sum
    ! write data:
    call dom%Put_Var( F%ncid, F%varid_cso2_sum , -1, pat, status )
    IF_NOTOK_RETURN(status=1)

    ! average:
    if ( i_nh3 > 0 ) pat = bud%drydepos%cnh3_ave_prev
    ! write data:
    call dom%Put_Var( F%ncid, F%varid_cnh3_ave_prev , -1, pat, status )
    IF_NOTOK_RETURN(status=1)

    ! average:
    if ( i_so2 > 0 ) pat = bud%drydepos%cso2_ave_prev
    ! write data:
    call dom%Put_Var( F%ncid, F%varid_cso2_ave_prev , -1, pat, status )
    IF_NOTOK_RETURN(status=1)


    
    ! clear:
    deallocate( pat, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Restart_Write
    

  ! *
  

  subroutine LE_Restart_Close( F, status )
  
#ifdef with_netcdf
    use NetCDF, only : NF90_Close
#endif
    use GO, only : goc

    ! --- in/out ----------------------------------
    
    type(T_LE_Restart_File), intent(inout)    ::  F
    integer, intent(out)                      ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Close'
    
    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! check implementation:    
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read boundary conditions from LE NetCDF output.")'); call goErr
    write (gol,'("Define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

#ifdef with_netcdf
    ! root only:
    if ( goc%root ) then
      ! close file:
      status = NF90_Close( F%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Close
    

  ! ********************************************************************
  ! ***
  ! *** restore
  ! ***
  ! ********************************************************************


  subroutine LE_Restart_Restore_Data( t, path, key, status )
  
    use GO           , only : TDate
    use Dims         , only : nx, ny, nz
    use LE_Grid      , only : dom
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Meteo_Data, only : ovolume  ! <-- restore into old volume !
    use LE_Meteo_Data, only : oairmass  ! <-- restore into old airmass !
    use LE_Meteo_Data, only : ohpres  ! <-- restore into old hpres !
    use dims         , only : expcls
#ifdef with_netcdf
    use NetCDF       , only : NF90_NOWRITE
    use NetCDF       , only : NF90_Open, NF90_Close
    use NetCDF       , only : NF90_Inq_VarID, NF90_Get_Var
#endif
    use LE_Restart_Tools, only : LE_Restart_Filename
    
    ! --- in/out ----------------------------------
    
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Restore_Data'
    
    ! --- local -----------------------------------
    
    character(len=256)    ::  fname
    logical               ::  exist
    integer               ::  ncid
    integer               ::  varid
    integer               ::  cmode
    integer               ::  glbo(2)

    !real, pointer          ::  ovolume(:,:,:)   ! (lon,lat,nz)
    
    ! --- begin -----------------------------------
    
    !call LE_Data_GetPointer( 'pvol', ovolume, status, check_units ='m3' )
    !IF_NOTOK_RETURN(status=1)    
    
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

    ! info ...
    write (gol,'("LE:       restore data from ",a," ...")') trim(fname); call goPr
    
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

    ! search variable:
    status = NF90_Inq_Varid( ncid, 'volume', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    ! restore into old volume:
    status = NF90_Get_Var( ncid, varid, ovolume(1:nx,1:ny,1:nz), &
                             start=(/glbo(1),glbo(2),1/), count=(/nx,ny,nz/) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! search variable:
    status = NF90_Inq_Varid( ncid, 'airmass', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    ! restore into old airmass:
    status = NF90_Get_Var( ncid, varid, oairmass, &
                             start=(/glbo(1),glbo(2),1/), count=(/nx,ny,nz/) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! search variable:
    status = NF90_Inq_Varid( ncid, 'hpres', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    ! restore into old hpres:
    status = NF90_Get_Var( ncid, varid, ohpres, &
                             start=(/glbo(1),glbo(2),1/), count=(/nx,ny,nz+1/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read data:
    status = NF90_Inq_Varid( ncid, 'expcls', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, expcls, &
                             start=glbo, count=(/nx,ny/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_Data
    

  ! ***
  

  subroutine LE_Restart_Restore_State( c, cg, aerh2o, bud, t, path, key, status )
  
#ifdef with_netcdf
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
#endif

    use GO         , only : TDate
    use dims       , only : nx, ny, nz
    use LE_Grid    , only : dom
    use Indices    , only : nspec, i_nh3, i_so2
    use LE_Budget  , only : T_Budget
    use LE_Restart_Tools, only : LE_Restart_Filename
    
    ! --- in/out ----------------------------------
    
    real, intent(out)               ::  c(nx,ny,nz,nspec)
    real, intent(out)               ::  cg(nx,ny,nspec)
    real, intent(out)               ::  aerh2o(nx,ny,nz)
    type(T_Budget), intent(inout)   ::  bud
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  path
    character(len=*), intent(in)    ::  key
    integer, intent(out)            ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Restart_Restore_State'
    
    ! --- local -----------------------------------
    
    character(len=256)    ::  fname
    logical               ::  exist
    integer               ::  glbo(2)
    
    integer       ::  ncid
    integer       ::  varid
    integer       ::  cmode
    
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

    ! info ...
    write (gol,'("LE:       restore state from ",a," ...")') trim(fname); call goPr
    
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
    status = NF90_Inq_Varid( ncid, 'c', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, c, &
                             start=(/glbo(1),glbo(2),1,1/), &
                             count=(/nx,ny,nz,nspec/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read data:
    status = NF90_Inq_Varid( ncid, 'cg', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, cg, &
                             start=(/glbo(1),glbo(2),1/), &
                             count=(/nx,ny,nspec/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read data:
    status = NF90_Inq_Varid( ncid, 'aerh2o', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, aerh2o, &
                             start=(/glbo(1),glbo(2),1/), &
                             count=(/nx,ny,nz/) )
    IF_NF90_NOTOK_RETURN(status=1)


    ! read NH3 average field used for compensation point in depos;
    ! should be part of the state rather than the data,
    ! but since depos is already dependend on N/S ration this
    ! fact is already ignored ...
    if ( i_nh3 > 0 ) then
      ! summed concentrations:
      status = NF90_Inq_Varid( ncid, 'cnh3_sum', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, bud%drydepos%cnh3_sum, &
                             start=glbo, count=(/nx,ny/) )
      IF_NF90_NOTOK_RETURN(status=1)
      ! counter:
      status = NF90_Get_Att( ncid, varid, 'nsum', bud%drydepos%cnh3_nsum )
      IF_NF90_NOTOK_RETURN(status=1)
      ! average:
      status = NF90_Inq_Varid( ncid, 'cnh3_ave_prev', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, bud%drydepos%cnh3_ave_prev, &
                             start=glbo, count=(/nx,ny/) )
      IF_NF90_NOTOK_RETURN(status=1)
    end if

    ! read SO2 average field used for co-deposition;
    if ( i_so2 > 0 ) then
      ! summed concentrations:
      status = NF90_Inq_Varid( ncid, 'cso2_sum', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, bud%drydepos%cso2_sum, &
                             start=glbo, count=(/nx,ny/) )
      IF_NF90_NOTOK_RETURN(status=1)
      ! counter:
      status = NF90_Get_Att( ncid, varid, 'nsum', bud%drydepos%cso2_nsum )
      IF_NF90_NOTOK_RETURN(status=1)
      ! average:
      status = NF90_Inq_Varid( ncid, 'cso2_ave_prev', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, bud%drydepos%cso2_ave_prev, &
                             start=glbo, count=(/nx,ny/) )
      IF_NF90_NOTOK_RETURN(status=1)
    end if

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif

    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_State
  

end module LE_Restart
