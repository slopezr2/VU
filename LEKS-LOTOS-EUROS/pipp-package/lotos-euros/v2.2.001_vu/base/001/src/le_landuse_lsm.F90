!###############################################################################
!
! Landuse_LSM  -  Make land sea mask
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!
#include "le.inc"
!
!###############################################################################


module LE_LandUse_LSM

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  Landuse_LSM_Init, Landuse_LSM_Done


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'Landuse_LSM'


  ! --- var --------------------------------------

contains


  ! ========================================================================


  subroutine Landuse_LSM_Init( rcF, status )

    use GO, only : TrcFile, ReadRc
    use GO, only : GoVarValue
    use LE_Grid     , only : ugg
    use C3PO   , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors

    use LE_Landuse_Data, only : waterbodies
    
    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_LandUse_LSM_Init'
  
    ! --- local ----------------------------------
    
    character(len=512)          ::  query

    character(len=512)    ::  fname
    logical               ::  exist
    integer               ::  varid
    character(len=128)    ::  description
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    integer, allocatable  ::  values_in(:,:)
    character(len=64)     ::  units_in
    
    integer               ::  start_ind(2)
    integer               ::  count_ind(2)
    
    ! --- begin ----------------------------------

    ! info ...
    write (gol,'("waterbody map ...")'); call goPr
    
    call ReadRc( RcF, 'my.landsea.waterbody.file', query, status )
    IF_NOTOK_RETURN(status=1)
    
    ! name of waterbody file:
    call goVarValue( trim(query), ';', 'file', '=', fname, status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    write (gol,'("  read ",a," ...")') trim(fname); call goPr
    
    description = 'var_name=wwb;units=-'
    ! open file:
    call file_in%Open( trim(fname), status )
    IF_NOTOK_RETURN(status=1)

    ! variable id:
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! init grid definition
    call file_in%Get_Grid( varid, grid_in, status, &
                             ugg_to=ugg, start_ind=start_ind, count_ind=count_ind  )
    IF_NOTOK_RETURN(status=1)

    ! storage
    allocate( values_in(grid_in%nlon, grid_in%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! read:
    call file_in%Get_Var( trim(description), values_in, units_in, status, &
                           start=start_ind, count=count_ind )
    IF_NOTOK_RETURN(status=1)

    ! close:
    call file_in%Close( status )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    write (gol,'("  regrid ...")'); call goPr
    
    ! convert, input variable is index 0,nwaterbody-1 --> change to 1, nwaterbody
    values_in = values_in + 1
    ! collect fractions per water type, don't keep regridding weights:
    call Grid_Convertors%Ugg_IndexFractions( grid_in, values_in, ugg, waterbodies, status, &
                                               clear=.true. )
    IF_NOTOK_RETURN(status=1)
        
    ! clear:
    deallocate( values_in, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! done:
    call grid_in%Done( status )
    IF_NOTOK_RETURN(status=1)
                
    ! ok
    status = 0

  end subroutine Landuse_LSM_Init


  ! ***


  subroutine Landuse_LSM_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_LandUse_LSM_Done'

    ! --- begin ----------------------------------

    ! ok
    status = 0

  end subroutine Landuse_LSM_Done


end module LE_LandUse_LSM

