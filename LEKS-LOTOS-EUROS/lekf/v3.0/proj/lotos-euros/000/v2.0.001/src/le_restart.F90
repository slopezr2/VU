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
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_LE_Restart_File

  public  ::  LE_Restart_Save
  public  ::  LE_Restart_Create, LE_Restart_Write, LE_Restart_Close

  public  ::  LE_Restart_Restore_Data
  public  ::  LE_Restart_Restore_State
  public  ::  LE_Restart_Restore
  public  ::  LE_Restart_FileName


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Restart'


  ! --- types ------------------------------------

  ! interface to restart file:
  type T_LE_Restart_File
    integer       ::  ncid
    integer       ::  dimid_nx, dimid_ny, dimid_nz, dimid_nspec
    integer       ::  varid_volume
    integer       ::  varid_expcls
    integer       ::  varid_aerh2o
    integer       ::  varid_c
    integer       ::  varid_cg
    integer       ::  varid_cnh3_sum
    integer       ::  varid_cso2_sum
    integer       ::  varid_cnh3_ave_prev
    integer       ::  varid_cso2_ave_prev
#ifdef with_pollen    
    integer       ::  varid_heatsum_polb
    integer       ::  varid_amt_polb_left
    integer       ::  varid_ripened_polb_left
    integer       ::  varid_heatsum_polo
    integer       ::  varid_amt_polo_left
    integer       ::  varid_ripened_polo_left
    integer       ::  varid_amt_polg_left
    integer       ::  varid_ripened_polg_left
#endif
  end type T_LE_Restart_File

  
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
  
  
  ! ********************************************************************
  ! ***
  ! *** save
  ! ***
  ! ********************************************************************
  
  subroutine LE_Restart_Save( c, cg, aerh2o, &
                              bud, t, path, key, status )
  
    use GO, only : TDate
    use dims, only : nx, ny, nz, nspec
#ifdef with_pollen    
    use indices, only : i_pol_b
#endif
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
  
    use GO, only : TDate
    use dims, only : nx, ny, nz, nspec
#ifdef with_pollen
    use Indices, only : i_pol_b, i_pol_g, i_pol_o
#endif
#ifdef with_netcdf
    use NetCDF, only : NF90_CLOBBER, NF90_NOCLOBBER
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Create
    use NetCDF, only : NF90_Def_Dim
    use NetCDF, only : NF90_Def_Var
#endif

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
    
#ifdef with_netcdf
    ! set creation mode flag:
    !cmode = NF90_NOCLOBBER     ! do not overwrite existing files, raise error instead
    cmode = NF90_CLOBBER       ! overwrite existing files if necessary

    ! create file:
    status = NF90_Create(  trim(path)//'/'//trim(fname), cmode, F%ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! create dimensions:
    status = NF90_Def_Dim( F%ncid, 'nx', nx, F%dimid_nx )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( F%ncid, 'ny', ny, F%dimid_ny )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( F%ncid, 'nz', nz, F%dimid_nz )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Def_Dim( F%ncid, 'nspec', nspec, F%dimid_nspec )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! create variable:
    status = NF90_Def_Var( F%ncid, 'volume', NF90_FLOAT, &
                           (/F%dimid_nx,F%dimid_ny,F%dimid_nz/), F%varid_volume )
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

#ifdef with_pollen
    if (i_pol_b > 0 ) then
      status = NF90_Def_Var( F%ncid, 'heatsum_polb', NF90_FLOAT, &
                            (/F%dimid_nx,F%dimid_ny/), F%varid_heatsum_polb )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'amt_polb_left', NF90_FLOAT, &
                            (/F%dimid_nx,F%dimid_ny/), F%varid_amt_polb_left )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'ripened_polb_left', NF90_FLOAT, &
                            (/F%dimid_nx,F%dimid_ny/), F%varid_ripened_polb_left )
      IF_NF90_NOTOK_RETURN(status=1)
    end if                              
    
    if (i_pol_g > 0 ) then
      status = NF90_Def_Var( F%ncid, 'amt_polg_left', NF90_FLOAT, &
                            (/F%dimid_nx,F%dimid_ny/), F%varid_amt_polg_left )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'ripened_polg_left', NF90_FLOAT, &
                            (/F%dimid_nx,F%dimid_ny/), F%varid_ripened_polg_left )
      IF_NF90_NOTOK_RETURN(status=1)
    end if                              
    
    if (i_pol_o > 0 ) then
      status = NF90_Def_Var( F%ncid, 'heatsum_polo', NF90_FLOAT, &
                            (/F%dimid_nx,F%dimid_ny/), F%varid_heatsum_polo )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'amt_polo_left', NF90_FLOAT, &
                            (/F%dimid_nx,F%dimid_ny/), F%varid_amt_polo_left )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Var( F%ncid, 'ripened_polo_left', NF90_FLOAT, &
                            (/F%dimid_nx,F%dimid_ny/), F%varid_ripened_polo_left )
      IF_NF90_NOTOK_RETURN(status=1)
    end if    
#endif
                          
#endif

    ! ok
    status = 0
    
  end subroutine LE_Restart_Create
    

  ! *
  
  subroutine LE_Restart_Write( F, c, cg, aerh2o, &
                               bud, status )
 
    use dims, only : nx, ny, nz, nspec
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Meteo_Data, only : volume
    use dims, only : expcls
    use Indices  , only : i_nh3, i_so2
#ifdef with_pollen
    use Indices, only : i_pol_b, i_pol_g, i_pol_o
#endif
    use LE_Emis, only : emis_set, max_emis
    use LE_Budget, only : T_Budget
#ifdef with_netcdf
    use NetCDF, only : NF90_EndDef
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Put_Att
#endif

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
    real, allocatable     ::  pat(:,:)
    
    integer               ::  iemis, itr

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

#ifdef with_netcdf
    ! NH3 compensation point
    ival = 0
    if ( i_nh3 > 0 ) ival = bud%drydepos%cnh3_nsum
    ! write counter as attribute:
    status = NF90_Put_Att( F%ncid, F%varid_cnh3_sum, 'nsum', ival )
    IF_NF90_NOTOK_RETURN(status=1)

    ! SO2 for co-deposition
    ival = 0
    if ( i_so2 > 0 ) ival = bud%drydepos%cso2_nsum
    ! write counter as attribute:
    status = NF90_Put_Att( F%ncid, F%varid_cso2_sum, 'nsum', ival )
    IF_NF90_NOTOK_RETURN(status=1)

    ! end definition phase:
    status = NF90_EndDef( F%ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    status = NF90_Put_Var( F%ncid, F%varid_volume, volume )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    status = NF90_Put_Var( F%ncid, F%varid_expcls, expcls )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    status = NF90_Put_Var( F%ncid, F%varid_aerh2o, aerh2o )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    status = NF90_Put_Var( F%ncid, F%varid_c, c )
    IF_NF90_NOTOK_RETURN(status=1)
    
    status = NF90_Put_Var( F%ncid, F%varid_cg, cg )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( pat(nx,ny) ) ; pat = 0.0
    ! summed concentrations:
    if ( i_nh3 > 0 ) pat = bud%drydepos%cnh3_sum
    status = NF90_Put_Var( F%ncid, F%varid_cnh3_sum, pat )
    IF_NF90_NOTOK_RETURN(status=1)
    ! summed concentrations:
    if ( i_so2 > 0 ) pat = bud%drydepos%cso2_sum
    status = NF90_Put_Var( F%ncid, F%varid_cso2_sum, pat )
    IF_NF90_NOTOK_RETURN(status=1)
    ! average:
    if ( i_nh3 > 0 ) pat = bud%drydepos%cnh3_ave_prev
    status = NF90_Put_Var( F%ncid, F%varid_cnh3_ave_prev, pat )
    IF_NF90_NOTOK_RETURN(status=1)
    ! average:
    if ( i_so2 > 0 ) pat = bud%drydepos%cso2_ave_prev
    status = NF90_Put_Var( F%ncid, F%varid_cso2_ave_prev, pat )
    IF_NF90_NOTOK_RETURN(status=1)


#ifdef with_pollen    
    ! accumulated heatsum for pollen ripening
    if ( i_pol_b > 0 .or. i_pol_o > 0 .or. i_pol_g > 0 ) then
      ! loop over emission input sets and find pollen set
      do iemis = 1, max_emis
        ! loop over emission sets
        if (trim(emis_set(iemis)%name) == 'silam-pollen' ) then
          ! loop over emitted tracers for this emission set
          do itr = 1, emis_set(iemis)%pollen%ntr
            ! if tracer is birch pollen
            if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer ) == 'pol_b' ) then
              ! extract heatsum from subset birch pollen and save to restart file
              status = NF90_Put_Var( F%ncid, F%varid_heatsum_polb, emis_set(iemis)%pollen%emp(itr)%Birch_Pollen%heatsum )
              IF_NF90_NOTOK_RETURN(status=1)
              ! extract amount of pollen left and save to restart file
              status = NF90_Put_Var( F%ncid, F%varid_amt_polb_left, emis_set(iemis)%pollen%emp(itr)%Birch_Pollen%rest_avail_grains )
              IF_NF90_NOTOK_RETURN(status=1)
              ! extract amount of ripened pollen left and save to restart file
              status = NF90_Put_Var( F%ncid, F%varid_ripened_polb_left, emis_set(iemis)%pollen%emp(itr)%Birch_Pollen%ripened_left )
              IF_NF90_NOTOK_RETURN(status=1)
 
            else if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer ) == 'pol_o' ) then
              ! extract heatsum from subset birch pollen and save to restart file
              status = NF90_Put_Var( F%ncid, F%varid_heatsum_polo, emis_set(iemis)%pollen%emp(itr)%Olive_Pollen%heatsum )
              IF_NF90_NOTOK_RETURN(status=1)
              ! extract amount of pollen left and save to restart file
              status = NF90_Put_Var( F%ncid, F%varid_amt_polo_left, emis_set(iemis)%pollen%emp(itr)%Olive_Pollen%rest_avail_grains )
              IF_NF90_NOTOK_RETURN(status=1)
              ! extract amount of ripened pollen left and save to restart file
              status = NF90_Put_Var( F%ncid, F%varid_ripened_polo_left, emis_set(iemis)%pollen%emp(itr)%Olive_Pollen%ripened_left )
              IF_NF90_NOTOK_RETURN(status=1)
              
            else if ( trim(emis_set(iemis)%pollen%emp(itr)%tracer ) == 'pol_g' ) then
              ! extract amount of pollen left and save to restart file
              status = NF90_Put_Var( F%ncid, F%varid_amt_polg_left, emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%rest_avail_grains )
              IF_NF90_NOTOK_RETURN(status=1)
              ! extract amount of ripened pollen left and save to restart file
              status = NF90_Put_Var( F%ncid, F%varid_ripened_polg_left, emis_set(iemis)%pollen%emp(itr)%Grass_Pollen%ripened_left )
              IF_NF90_NOTOK_RETURN(status=1)
              
            endif 
          end do
        end if
      end do
    end if
#endif
    
    ! clear:
    deallocate( pat )
#endif

    ! ok
    status = 0
    
  end subroutine LE_Restart_Write
    

  ! *
  

  subroutine LE_Restart_Close( F, status )
  
#ifdef with_netcdf
    use NetCDF, only : NF90_Close
#endif

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
    ! close file:
    status = NF90_Close( F%ncid )
    IF_NF90_NOTOK_RETURN(status=1)
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
  
    use GO     , only : TDate
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Meteo_Data, only : ovolume  ! <-- restore into old volume !
    use dims   , only : expcls
#ifdef with_netcdf
    use NetCDF , only : NF90_NOWRITE
    use NetCDF , only : NF90_Open, NF90_Close
    use NetCDF , only : NF90_Inq_VarID, NF90_Get_Var
#endif
    
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
    
    integer       ::  ncid
    integer       ::  varid
    integer       ::  cmode

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

    ! read data:
    status = NF90_Inq_Varid( ncid, 'volume', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, ovolume )   ! <-- restore into old volume !
    IF_NF90_NOTOK_RETURN(status=1)

    ! read data:
    status = NF90_Inq_Varid( ncid, 'expcls', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, expcls )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = nf90_close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_Data
    

  ! ***
  

  subroutine LE_Restart_Restore_State( c, cg, aerh2o, bud, t, path, key, status )
  
    use GO  , only : TDate
    use dims, only : nx, ny, nz, nspec
    use Indices  , only : i_nh3, i_so2
    use LE_Budget, only : T_Budget
#ifdef with_netcdf
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
#endif
    
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
    status = NF90_Get_Var( ncid, varid, c )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inq_Varid( ncid, 'cg', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, cg )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read data:
    status = NF90_Inq_Varid( ncid, 'aerh2o', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, aerh2o )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read NH3 average field used for compensation point in depos;
    ! should be part of the state rather than the data,
    ! but since depos is already dependend on N/S ration this
    ! fact is already ignored ...
    if ( i_nh3 > 0 ) then
      ! summed concentrations:
      status = NF90_Inq_Varid( ncid, 'cnh3_sum', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, bud%drydepos%cnh3_sum )
      IF_NF90_NOTOK_RETURN(status=1)
      ! counter:
      status = NF90_Get_Att( ncid, varid, 'nsum', bud%drydepos%cnh3_nsum )
      IF_NF90_NOTOK_RETURN(status=1)
      ! average:
      status = NF90_Inq_Varid( ncid, 'cnh3_ave_prev', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, bud%drydepos%cnh3_ave_prev )
      IF_NF90_NOTOK_RETURN(status=1)
    end if

    ! read SO2 average field used for co-deposition;
    if ( i_so2 > 0 ) then
      ! summed concentrations:
      status = NF90_Inq_Varid( ncid, 'cso2_sum', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, bud%drydepos%cso2_sum )
      IF_NF90_NOTOK_RETURN(status=1)
      ! counter:
      status = NF90_Get_Att( ncid, varid, 'nsum', bud%drydepos%cso2_nsum )
      IF_NF90_NOTOK_RETURN(status=1)
      ! average:
      status = NF90_Inq_Varid( ncid, 'cso2_ave_prev', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, bud%drydepos%cso2_ave_prev )
      IF_NF90_NOTOK_RETURN(status=1)
    end if

    ! close file:
    status = nf90_close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif

    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_State
  

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
    status = nf90_close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_1D
  

  ! *
  

  subroutine LE_Restart_Restore_2D( varname, values, t, path, key, status )
  
    use GO  , only : TDate
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
    status = NF90_Get_Var( ncid, varid, values )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = nf90_close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_2D
  

  ! *
  

  subroutine LE_Restart_Restore_3D( varname, values, t, path, key, status )
  
    use GO  , only : TDate
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
    status = NF90_Get_Var( ncid, varid, values )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = nf90_close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_3D
  
  
  ! *
  

  subroutine LE_Restart_Restore_4D( varname, values, t, path, key, status )
  
    use GO  , only : TDate
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
    status = NF90_Get_Var( ncid, varid, values )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = nf90_close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif
    
    ! ok
    status = 0
    
  end subroutine LE_Restart_Restore_4D
  
  
end module LE_Restart
