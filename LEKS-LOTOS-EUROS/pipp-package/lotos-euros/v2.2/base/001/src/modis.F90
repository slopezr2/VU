!#################################################################
!
! NAME
!
!   modis  -  read MODIS fields
!
!### GO routines #######################################################
!
!
!  ! substitutes for message routines from GO modules
!  
!  ! message line:
!  character(len=256)            ::  gol
!
!  ! display message:
!  subroutine goPr
!    write (*,'(a)') trim(gol)
!  end subroutine goPr
!
!  ! display error message:
!  subroutine goErr
!    write (*,'("ERROR - ",a)') trim(gol)
!  end subroutine goErr
!  
!  ! free file unit:
!  subroutine goGetFU( fu, status )
!    integer, intent(out)    ::  fu
!    integer, intent(out)    ::  status
!    logical                 ::  opened
!    fu = 456
!    do
!      inquire( unit=fu, opened=opened )
!      if ( .not. opened ) exit
!      fu = fu + 1
!    end do
!    status = 0
!  end subroutine goGetFU
!
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_HF90_NOTOK_RETURN(action) if (status/=HF90_SUCCEED) then; gol=trim(HF90_ERROR_MSG); call goErr; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#################################################################

module MODIS

  use GO, only : gol, goPr, goErr
  use GO, only : goLoCase

#ifdef with_hdf4
  use HDF4, only : HF90_SUCCEED, HF90_ERROR_MSG
#endif
  
  
  implicit none


  ! --- in/out -----------------------------

  private
  
  public  ::  T_MODIS
  public  ::  Init, Done

  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MODIS'

  ! expected maximum length of units:
  integer, parameter    ::  len_units     = 64
  integer, parameter    ::  len_long_name = 512

  ! --- types ------------------------------
  
  ! data file:
  type T_MODIS
    ! full path:
    character(len=512)              ::  hdf_name
    ! grid size:
    integer                         ::  nx, ny
    ! storage:
    real(8), pointer                ::  Scan_Start_Time(:,:)
    real(4), pointer                ::  Longitude(:,:)
    real(4), pointer                ::  Latitude (:,:)
    real(4), pointer                ::  Cloud_Fraction_Land (:,:)
    real(4), pointer                ::  Cloud_Fraction_Ocean(:,:)
    real(4), pointer                ::  Optical_Depth_Land_And_Ocean(:,:)
    real(4), pointer                ::  Optical_Depth_Ratio_Small_Land_And_Ocean(:,:)
    character(len=len_units)        ::  Scan_Start_Time__units
    character(len=len_units)        ::  Longitude__units
    character(len=len_units)        ::  Latitude__units
    character(len=len_units)        ::  Cloud_Fraction_Land__units
    character(len=len_units)        ::  Cloud_Fraction_Ocean__units
    character(len=len_units)        ::  Optical_Depth_Land_And_Ocean__units
    character(len=len_units)        ::  Optical_Depth_Ratio_Small_Land_And_Ocean__units
    character(len=len_long_name)    ::  Scan_Start_Time__long_name
    character(len=len_long_name)    ::  Longitude__long_name
    character(len=len_long_name)    ::  Latitude__long_name
    character(len=len_long_name)    ::  Cloud_Fraction_Land__long_name
    character(len=len_long_name)    ::  Cloud_Fraction_Ocean__long_name
    character(len=len_long_name)    ::  Optical_Depth_Land_And_Ocean__long_name
    character(len=len_long_name)    ::  Optical_Depth_Ratio_Small_Land_And_Ocean__long_name
  end type T_MODIS


  ! --- interfaces -------------------------

  interface Init
    module procedure modis_Init
  end interface

  interface Done
    module procedure modis_Done
  end interface

  
contains


  ! ====================================================================
  
  
  subroutine modis_Init( modis, hdf_name, status )
  
#ifdef with_hdf4
    use HDF4, only : HF90_DEBUG_LEVEL, HF90_DEBUG_NONE, HF90_DEBUG_ERRORS
    use HDF4, only : HF90_READ, HF90_GLOBAL
    use HDF4, only : HF90_Open, HF90_Close
    use HDF4, only : HF90_Inq_VarID, HF90_Get_Var
    use HDF4, only : HF90_Get_Att
#endif

    ! --- in/out --------------------------------

    type(T_MODIS), intent(out)              ::  modis
    character(len=*), intent(in)            ::  hdf_name
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/modis_Init'

    ! --- local ---------------------------------
    
    logical                       ::  exist
#ifdef with_hdf4
    integer                       ::  hdf_id
#endif
    
    ! --- begin -------------------------------
    
    ! check ...
    inquire( file=trim(hdf_name), exist=exist )
    if ( .not. exist ) then
      write (gol,'("MODIS file not found:")'); call goErr
      write (gol,'("  ",a)') trim(hdf_name); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! store:
    modis%hdf_name = trim(hdf_name)
    
#ifdef with_hdf4

    ! write extra messages in case of errors:
    HF90_DEBUG_LEVEL = HF90_DEBUG_ERRORS

    ! open hdf file:
    call HF90_Open( trim(hdf_name), HF90_READ, hdf_id, status )
    IF_HF90_NOTOK_RETURN(status=1)
    
    ! setup storage:
    nullify( modis%Scan_Start_Time )
    nullify( modis%Longitude       )
    nullify( modis%Latitude        )
    nullify( modis%cloud_fraction_land  )
    nullify( modis%cloud_fraction_ocean )
    nullify( modis%Optical_Depth_Land_And_Ocean             )
    nullify( modis%Optical_Depth_Ratio_Small_Land_And_Ocean )

    ! read times:
    call modis_Read_Record_r8( hdf_id, 'Scan_Start_Time', &
                                 modis%Scan_Start_Time, &
                                 modis%Scan_Start_Time__long_name, &
                                 modis%Scan_Start_Time__units, &
                                 status )
    IF_NOTOK_RETURN(status=1)
    
    ! store dimensions:
    modis%nx = size(modis%Scan_Start_Time,1)
    modis%ny = size(modis%Scan_Start_Time,2)

    ! read longitudes:
    call modis_Read_Record_r4( hdf_id, 'Longitude', &
                                modis%Longitude, &
                                modis%Longitude__long_name, &
                                modis%Longitude__units, &
                                status )
    IF_NOTOK_RETURN(status=1)

    ! read latitudes:
    call modis_Read_Record_r4( hdf_id, 'Latitude', &
                                modis%Latitude, &
                                modis%Latitude__long_name, &
                                modis%Latitude__units, &
                                status )
    IF_NOTOK_RETURN(status=1)
    
    ! read cloud fractions over land:
    call modis_Read_Record_r4( hdf_id, 'Cloud_Fraction_Land', &
                                modis%cloud_fraction_land, &
                                modis%cloud_fraction_land__long_name, &
                                modis%cloud_fraction_land__units, &
                                 status )
    IF_NOTOK_RETURN(status=1)
    
    ! read cloud fractions over land:
    call modis_Read_Record_r4( hdf_id, 'Cloud_Fraction_Ocean', &
                                modis%cloud_fraction_ocean, &
                                modis%cloud_fraction_ocean__long_name, &
                                modis%cloud_fraction_ocean__units, &
                                 status )
    IF_NOTOK_RETURN(status=1)
    
    ! read cloud fractions over land:
    call modis_Read_Record_r4( hdf_id, 'Optical_Depth_Land_And_Ocean', &
                                modis%Optical_Depth_Land_And_Ocean, &
                                modis%Optical_Depth_Land_And_Ocean__long_name, &
                                modis%Optical_Depth_Land_And_Ocean__units, &
                                 status )
    IF_NOTOK_RETURN(status=1)
    
    ! read cloud fractions over land:
    call modis_Read_Record_r4( hdf_id, 'Optical_Depth_Ratio_Small_Land_And_Ocean', &
                                modis%Optical_Depth_Ratio_Small_Land_And_Ocean, &
                                modis%Optical_Depth_Ratio_Small_Land_And_Ocean__long_name, &
                                modis%Optical_Depth_Ratio_Small_Land_And_Ocean__units, &
                                 status )
    IF_NOTOK_RETURN(status=1)

    ! open hdf file:
    call HF90_Close( hdf_id, status )
    IF_HF90_NOTOK_RETURN(status=1)

    ! no extra messages anymore:
    HF90_DEBUG_LEVEL = HF90_DEBUG_NONE

#else

    write (gol,'("MODIS module requires compilation with HDF4")'); call goErr
    TRACEBACK; status=1; return

#endif
    
    ! ok
    status = 0
    
  end subroutine modis_Init
  
  
  ! ***
  
  
  
  
  subroutine modis_Done( modis, status )
  
    ! --- in/out --------------------------------

    type(T_MODIS), intent(inout)            ::  modis
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/modis_Done'

    ! --- local ---------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate( modis%Scan_Start_Time )
    deallocate( modis%longitude       )
    deallocate( modis%latitude        )
    deallocate( modis%cloud_fraction_land  )
    deallocate( modis%cloud_fraction_ocean )
    deallocate( modis%Optical_Depth_Land_And_Ocean             )
    deallocate( modis%Optical_Depth_Ratio_Small_Land_And_Ocean )

    ! empty:
    modis%hdf_name = ''
    
    ! ok
    status = 0
    
  end subroutine modis_Done
  
  
  ! ***
  

  subroutine modis_Read_Record_r4( hdf_id, var_name, values, long_name, units, status )
  
#ifdef with_hdf4
    use HDF4, only : HF90_Inq_VarID, HF90_Inquire_Variable, HF90_Get_Var
    use HDF4, only : HF90_Inquire_Dimension
    use HDF4, only : HF90_Get_Att
#endif

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/modis_Read_Record_r4'
    
    ! variable kind and rank:
    integer, parameter    ::  vkind = 4
    integer, parameter    ::  vrank = 2

    ! --- in/out ---------------------------------
    
    integer, intent(in)             ::  hdf_id
    character(len=*), intent(in)    ::  var_name
    real(vkind), pointer            ::  values(:,:)    ! null on input
    character(len=*), intent(out)   ::  long_name
    character(len=*), intent(out)   ::  units
    integer, intent(out)            ::  status
    
    ! --- local ----------------------------------

#ifdef with_hdf4    
    integer         ::  var_id
    real(8)         ::  add_offset, scale_factor
    integer         ::  dimids(vrank)
    integer         ::  vshape(vrank)
    integer         ::  irank
#endif

    ! --- begin ----------------------------------

#ifdef with_hdf4
      
    ! get variable handle:
    call HF90_Inq_VarID( hdf_id, trim(var_name), var_id, status )
    IF_HF90_NOTOK_RETURN(status=1)

    ! get dimension id's :    
    call HF90_Inquire_Variable( hdf_id, var_id, status, dimids=dimids )
    IF_HF90_NOTOK_RETURN(status=1)
    ! get shape:
    do irank = 1, vrank
      call HF90_Inquire_Dimension( hdf_id, dimids(irank), status, len=vshape(irank) )
      IF_HF90_NOTOK_RETURN(status=1)
    end do
    ! setup output:
    allocate( values(vshape(1),vshape(2)) )
    
    ! read data:
    call HF90_Get_Var( hdf_id, var_id, values, status )
    IF_HF90_NOTOK_RETURN(status=1)
    
    ! get units:
    call HF90_Get_Att( hdf_id, var_id, 'long_name', long_name, status )
    IF_HF90_NOTOK_RETURN(status=1)
    
    ! get units:
    call HF90_Get_Att( hdf_id, var_id, 'units', units, status )
    IF_HF90_NOTOK_RETURN(status=1)
    
    ! offset and scale factor:
    call HF90_Get_Att( hdf_id, var_id, 'add_offset', add_offset, status )
    IF_HF90_NOTOK_RETURN(status=1)
    call HF90_Get_Att( hdf_id, var_id, 'scale_factor', scale_factor, status )
    IF_HF90_NOTOK_RETURN(status=1)
    ! apply:
    values = values * scale_factor + add_offset
    
    !! info ..
    !print *, 'modis_Read_Record: ', trim(var_name), ', ', trim(units), add_offset, scale_factor

#else

    ! dummy assignment to avoid warnings ...
    values    = -999.9
    long_name = '-'
    units     = '-'
    ! error ...
    write (gol,'("MODIS module requires compilation with HDF4")'); call goErr
    TRACEBACK; status=1; return

#endif
    
    ! ok
    status = 0
    
  end subroutine modis_Read_Record_r4
  
  
  ! ***
  

  subroutine modis_Read_Record_r8( hdf_id, var_name, values, long_name, units, status )
  
#ifdef with_hdf4
    use HDF4, only : HF90_Inq_VarID, HF90_Inquire_Variable, HF90_Get_Var
    use HDF4, only : HF90_Inquire_Dimension
    use HDF4, only : HF90_Get_Att
#endif

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/modis_Read_Record_r4'
    
    ! variable kind and rank:
    integer, parameter    ::  vkind = 8
    integer, parameter    ::  vrank = 2

    ! --- in/out ---------------------------------
    
    integer, intent(in)             ::  hdf_id
    character(len=*), intent(in)    ::  var_name
    real(vkind), pointer            ::  values(:,:)    ! null on input
    character(len=*), intent(out)   ::  long_name
    character(len=*), intent(out)   ::  units
    integer, intent(out)            ::  status
    
    ! --- local ----------------------------------

#ifdef with_hdf4
    integer         ::  var_id
    real(8)         ::  add_offset, scale_factor
    integer         ::  dimids(vrank)
    integer         ::  vshape(vrank)
    integer         ::  irank
#endif

    ! --- begin ----------------------------------

#ifdef with_hdf4
      
    ! get variable handle:
    call HF90_Inq_VarID( hdf_id, trim(var_name), var_id, status )
    IF_HF90_NOTOK_RETURN(status=1)

    ! get dimension id's :    
    call HF90_Inquire_Variable( hdf_id, var_id, status, dimids=dimids )
    IF_HF90_NOTOK_RETURN(status=1)
    ! get shape:
    do irank = 1, vrank
      call HF90_Inquire_Dimension( hdf_id, dimids(irank), status, len=vshape(irank) )
      IF_HF90_NOTOK_RETURN(status=1)
    end do
    ! setup output:
    allocate( values(vshape(1),vshape(2)) )
    
    ! read data:
    call HF90_Get_Var( hdf_id, var_id, values, status )
    IF_HF90_NOTOK_RETURN(status=1)
    
    ! get long_name:
    call HF90_Get_Att( hdf_id, var_id, 'long_name', long_name, status )
    IF_HF90_NOTOK_RETURN(status=1)
    
    ! get units:
    call HF90_Get_Att( hdf_id, var_id, 'units', units, status )
    IF_HF90_NOTOK_RETURN(status=1)
    
    ! offset and scale factor:
    call HF90_Get_Att( hdf_id, var_id, 'add_offset', add_offset, status )
    IF_HF90_NOTOK_RETURN(status=1)
    call HF90_Get_Att( hdf_id, var_id, 'scale_factor', scale_factor, status )
    IF_HF90_NOTOK_RETURN(status=1)
    ! apply:
    values = values * scale_factor + add_offset
    
    ! info ..
    print *, 'modis_Read_Record: ', trim(var_name), ', ', trim(units), add_offset, scale_factor

#else

    ! dummy assignment to avoid warnings:
    values    = -999.9
    long_name = '-'
    units     = '-'
    ! error ...
    write (gol,'("MODIS module requires compilation with HDF4")'); call goErr
    TRACEBACK; status=1; return

#endif
    
    ! ok
    status = 0
    
  end subroutine modis_Read_Record_r8
  
  
  ! ***
    
  

end module MODIS


! ######################################################################
! ###
! ### test
! ###
! ######################################################################
!
! gfortran -o test.x -ffree-line-length-none -Dwith_hdf4 -I${HDF4_HOME}/include hdf4.F90 modis.F90 -L${HDF4_HOME}/lib -l mfhdf -l df -l z -l jpeg
!
!program test
!
!  use MODIS
!  
!  character(len=*), parameter  ::  fname = '/modas/archive/projects/SATELLIETEN/MODIS_c5/2006/AQUA1/MYD04_L2.A2006001.1210.005.2007128013807.hdf'
!  
!  type(T_MODIS)     ::  modis
!  integer           ::  status
!  
!  print *, 'TEST: begin'
!
!  call Init( modis, fname, status )
!  if (status/=0) stop 'from init'
!  
!  call Done( modis, status )
!  if (status/=0) stop 'from done'
!
!  print *, 'TEST: end'
!
!end program test
!
