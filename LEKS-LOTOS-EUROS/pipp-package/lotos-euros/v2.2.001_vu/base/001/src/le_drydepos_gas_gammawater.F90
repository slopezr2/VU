!###############################################################################
!
! Deposition
!
! Version including aerosol sedimentation following Zhang.
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

module LE_DryDepos_Gas_GammaWater

  use GO, only : gol, goPr, goErr

  use NetCDF, only : NF90_NOERR, nf90_strerror

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  gamma_water
  
  public  ::  LE_DryDepos_Gas_GammaWater_Init, LE_DryDepos_Gas_GammaWater_Done


  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_DryDepos_Gas_GammaWater' 
  

  ! --- var --------------------------------------

  ! gamma soil for water (from gammawater.nc file)
  real, allocatable   ::  gamma_water(:,:)



contains


  ! ========================================================================

  
  subroutine LE_DryDepos_Gas_GammaWater_Init( rcF, status )
  
    use GO, only : TrcFile, ReadRc

    use LE_Grid, only : ugg
    use C3PO   , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors

    !use file_nc
    
    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_DryDepos_Gas_GammaWater_Init'

    ! --- local ----------------------------------
    
    character(len=512)    ::  fname
    logical               ::  exist
    integer               ::  varid
    character(len=128)    ::  description
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    real, allocatable     ::  values_in(:,:)
    character(len=64)     ::  units_in
    real                  ::  missing_value
    
    ! --- begin ----------------------------------
    
    ! storage:
    allocate( gamma_water(ugg%nlon,ugg%nlat) )

    ! name of input file:
    call ReadRc( rcF, 'deposition.gammawater.file', fname, status )
    IF_NOTOK_RETURN(status=1)
    
    ! check for specials ...
    select case ( trim(fname) )
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'zero' )  ! set field to zero
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! for production of the file, first perform a run
        ! with gamma-water set to zero:
        gamma_water = 0.0
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default   ! filename provided
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        ! check ...
        inquire( file=trim(fname), exist=exist )
        if ( .not.exist ) then
          write (gol,'("file not found : ",a)') trim(fname); call goErr
          TRACEBACK; status=1; return
        end if

        description = 'var_name=gammawater'
        ! open file:
        call file_in%Open( trim(fname), status )
        IF_NOTOK_RETURN(status=1)

        ! variable id:
        call file_in%Inq_VarID( trim(description), varid, status )
        IF_NOTOK_RETURN(status=1)
        ! init grid definition
        call file_in%Get_Grid( varid, grid_in, status )
        IF_NOTOK_RETURN(status=1)

        ! storage
        allocate( values_in(grid_in%nlon, grid_in%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! read:
        call file_in%Get_Var( trim(description), values_in, units_in, status, &
                               start=(/1,1/), count=(/grid_in%nlon,grid_in%nlat/) )
        IF_NOTOK_RETURN(status=1)

        ! map from input to model grid:
        call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, gamma_water, status )
        IF_NOTOK_RETURN(status=1)

        ! clear:
        deallocate( values_in )
        call grid_in%Done(status)
        IF_NOTOK_RETURN(status=1)

        ! clear:
        call file_in%Close( status )
        IF_NOTOK_RETURN(status=1)
        
    end select
    
    !! testing ...
    !call nc_dump( 'gw.nc', gamma_water, 'gw', (/'x','y'/), status )
    !IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_DryDepos_Gas_GammaWater_Init


  ! ***
    

  subroutine LE_DryDepos_Gas_GammaWater_Done( status )
  
    ! --- in/out ---------------------------------

    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_DryDepos_Gas_GammaWater_Done'

    ! --- local ----------------------------------
        
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( gamma_water )
    
    ! ok
    status = 0
    
  end subroutine LE_DryDepos_Gas_GammaWater_Done


end module LE_DryDepos_Gas_GammaWater

