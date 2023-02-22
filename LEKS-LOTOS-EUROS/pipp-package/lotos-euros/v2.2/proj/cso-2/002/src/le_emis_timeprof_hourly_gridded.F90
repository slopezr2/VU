!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################

module LE_Emis_Timeprof_Hourly_Gridded

  use GO, only : gol, goErr, goPr

  use C3PO   , only : T_File_Ugg, T_Grid_Ugg
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Timeprof_Hourly_Gridded
  public  ::  LE_Emis_Timeprof_Hourly_Gridded_Init, LE_Emis_Timeprof_Hourly_Gridded_Done
  public  ::  LE_Emis_Timeprof_Hourly_Gridded_Setup
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Timeprof_Hourly_Gridded'
  
  
  ! --- types --------------------------------
  
  ! storage for emission timeprofile:
  type T_Emis_Timeprof_Hourly_Gridded
    
    ! dimensions:
    integer                   ::  year

    ! fraction assigned to hour:
    real, allocatable         ::  factor(:,:)  ! (nx,ny)
    real, allocatable         ::  factor_in(:,:)  ! ( grid_in%nlon,grid_in%nlat)

    integer                   ::  varid_fact
    type(T_File_Ugg)          ::  file_in
    type(T_Grid_Ugg)          ::  grid_in
    integer                   ::  start_ind(2)
    integer                   ::  count_ind(2)
    character(len=1024)       ::  description
  end type



contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_Timeprof_Hourly_Gridded_Init( tprof_hourly_Gridded, indir, fname, status )

    use LE_Grid, only : ugg
    
    ! --- in/out ------------------------------
    
    type(T_Emis_Timeprof_Hourly_Gridded), intent(out) ::  tprof_hourly_Gridded
    character(len=*), intent(in)                      ::  indir
    character(len=*), intent(in)                      ::  fname
    integer, intent(out)                              ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Timeprof_Hourly_Init'
    
    ! --- local -------------------------------
    
    character(len=512)             ::  filename        
    logical                        ::  exist    
    

    ! --- begin -------------------------------
    
    ! store:
    !tprof_hourly%year         = year
    
    !complete file name
    write( filename, '(a,"/",a)' ) trim(indir), trim(fname)             
    
    ! file should be present:
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! Open input file
    call tprof_hourly_gridded%file_in%Open( trim(filename), status )
    IF_NOTOK_RETURN(status=1)
        
    ! Variable id:
    tprof_hourly_gridded%description = 'var_name=time_factors'
    call tprof_hourly_gridded%file_in%Inq_VarID( trim(tprof_hourly_gridded%description), &
                                                   tprof_hourly_gridded%varid_fact, status )
    IF_NOTOK_RETURN(status=1)
        
    ! init grid definition
    call tprof_hourly_gridded%file_in%Get_Grid( tprof_hourly_gridded%varid_fact, tprof_hourly_gridded%grid_in, status, &
                            ugg_to=ugg, start_ind=tprof_hourly_gridded%start_ind, count_ind=tprof_hourly_gridded%count_ind )
    IF_NOTOK_RETURN(status=1)
        
    ! init:
    allocate( tprof_hourly_gridded%factor_in( tprof_hourly_gridded%grid_in%nlon,tprof_hourly_gridded%grid_in%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! redistributed factors
    allocate( tprof_hourly_gridded%factor( ugg%nlon, ugg%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Emis_Timeprof_Hourly_Gridded_Init
  
  
  ! ***
  

  subroutine LE_Emis_Timeprof_Hourly_Gridded_Done( tprof_hourly_gridded, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Timeprof_Hourly_Gridded), intent(inout) ::  tprof_hourly_gridded
    integer, intent(out)                                ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Gridded_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate(tprof_hourly_gridded%factor, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    
    ! ok
    status = 0

  end subroutine LE_Emis_Timeprof_Hourly_Gridded_Done
  
  
  ! ***

  
  subroutine LE_Emis_Timeprof_Hourly_Gridded_Setup( tprof_hourly_gridded, ihour, status )
  
    use LE_Grid, only : ugg
    use LE_Data_Common, only : Grid_Convertors

    ! --- in/out ------------------------------
    
    type(T_Emis_Timeprof_Hourly_Gridded), intent(inout) ::  tprof_hourly_gridded
    integer, intent(in)                                 ::  ihour                 ! Hour of the year (UTC)
    integer, intent(out)                                ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Time_Prof_Gridded_Done'
    
    ! --- local -------------------------------
    
    integer             ::  start_ind_3d(3)
    integer             ::  count_ind_3d(3)
    character(len=32)   ::  units_in
    
    ! --- begin -------------------------------
    
    ! Read factors for this hour        
    
    start_ind_3d = (/tprof_hourly_gridded%start_ind(1),tprof_hourly_gridded%start_ind(2),ihour/)
    count_ind_3d = (/tprof_hourly_gridded%count_ind(1),tprof_hourly_gridded%count_ind(2),1/)
    
    ! Read factors on input grid
    call tprof_hourly_gridded%file_in%Get_Var( trim(tprof_hourly_gridded%description),tprof_hourly_gridded%factor_in, units_in, status, &
                                        start=start_ind_3d,count=count_ind_3d)
    IF_NOTOK_RETURN(status=1)
    
    ! Regrid to LE-grid
    call Grid_Convertors%Ugg_AreaAver( tprof_hourly_gridded%grid_in, tprof_hourly_gridded%factor_in, &
                                        ugg, tprof_hourly_gridded%factor, status )
    IF_NOTOK_RETURN(status=1)
                                            
    ! ok
    status = 0

  end subroutine LE_Emis_Timeprof_Hourly_Gridded_Setup
  

end module LE_Emis_Timeprof_Hourly_Gridded

