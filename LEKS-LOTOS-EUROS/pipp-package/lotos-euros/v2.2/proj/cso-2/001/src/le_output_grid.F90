!###############################################################################
!
! NAME
!
!   LE_Output_grid  -  LOTOS-EUROS output of grid
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

module LE_Output_Grid

  use GO              , only : gol, goPr, goErr
#ifdef with_netcdf
  use NetCDF          , only : NF90_StrError, NF90_NOERR
#endif
  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_grid


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Grid'


  ! --- types ------------------------------

  type T_LE_Output_Grid
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
    ! file created ?
    logical                     ::  created
    ! file name:
    character(len=1024)         ::  fname
    ! file handle:
    integer                     ::  ncid
    ! dimension handles:
    integer                     ::  dimid_lon
    integer                     ::  dimid_lat
    !
  contains
    procedure :: Init            => LE_Output_Grid_Init
    procedure :: Done            => LE_Output_Grid_Done
    procedure :: PutOut          => LE_Output_Grid_PutOut
  end type T_LE_Output_Grid


contains


  ! ====================================================


  subroutine LE_Output_Grid_Init( self, rcF, rckey, status )

    use GO     , only : TrcFile
    use LE_Data, only : LE_Data_Enable
    use LE_Output_Common, only : Init

    ! --- in/out --------------------------------

    class(T_LE_Output_Grid), intent(out)    ::  self
    type(TrcFile), intent(in)               ::  rcF
    character(len=*), intent(in)            ::  rckey
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Grid_Init'

    ! --- local ---------------------------------

    ! --- begin ---------------------------------

    ! init common stuff:
    call Init( self%com, rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)
    
    ! enable meteo:
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)
    
    ! file no created yet, first meteo varialbes should be defined ...
    self%created = .false.

    ! ok
    status = 0

  end subroutine LE_Output_Grid_Init


  ! ***


  subroutine LE_Output_Grid_Done( self, status )

    use LE_Output_Common, only : Done

    ! --- in/out --------------------------------

    class(T_LE_Output_Grid), intent(inout)    ::  self
    integer, intent(out)                      ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Grid_Done'

    ! --- begin ---------------------------------

    ! done with common stuff ...
    call Done( self%com, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Grid_Done


  ! ***


  subroutine LE_Output_Grid_PutOut( self, status )

#ifdef with_netcdf
    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Var
    use NetCDF , only : NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF , only : NF90_NETCDF4
    use NetCDF , only : NF90_REAL
#endif

    use GO              , only : goc
    use C3PO            , only : T_Grid_NcDef
    use Dims            , only : nx, ny
    use LE_Grid         , only : glb_ugg
    use LE_Data         , only : LE_Data_GetPointer
    use LE_Output_Common, only : PutOut_GlobalAttributes
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains

    ! --- in/out --------------------------------

    class(T_LE_Output_Grid), intent(inout)    ::  self
    integer, intent(out)                      ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Grid_PutOut'

    ! --- local ---------------------------------

#ifdef with_netcdf
    integer               ::  cmode
#endif
    type(T_Grid_NcDef)    ::  gncd
    real, pointer         ::  area(:,:,:)      ! (lon,lat,1)
    character(len=32)     ::  units
    integer               ::  varid
    real                  ::  pat(nx,ny)

    ! --- begin ---------------------------------

    ! not created yet?
    if ( .not. self%created ) then
    
      ! point to meteo data, also obtain units:
      call LE_Data_GetPointer( 'area', area, status, units=units )
      IF_NOTOK_RETURN(status=1)

#ifdef with_netcdf
      ! write on root only:
      if ( goc%root ) then

        ! new file name:
        write (self%fname,'(a,"/",a,"_",a,"_",a,".nc")') &
                        trim(self%com%outdir), trim(self%com%model), trim(self%com%expid), 'grid'

        ! set creation mode flag:
        cmode = NF90_CLOBBER       ! overwrite existing files

        ! enable large file support:
        cmode = or( cmode, NF90_NETCDF4 )

        ! create file:
        status = NF90_Create( trim(self%fname), cmode, self%ncid )
        if ( status /= NF90_NOERR ) then
          write (gol,'("creating file : ")'); call goErr
          write (gol,'("  file name  : ",a)') trim(self%fname); call goErr
          write (gol,'("  nf90 error : ",a)') trim(nf90_strerror(status)); call goErr
          TRACEBACK; status=1; return
        end if

         ! write global attributes:
        call PutOut_GlobalAttributes( self%com, self%ncid, status )
        IF_NOTOK_RETURN(status=1)

        ! grid dimensions/variables
        call glb_ugg%DefGrid_NetCDF( gncd, self%ncid, status, &
                                      dimid_lon=self%dimid_lon, dimid_lat=self%dimid_lat )
        IF_NOTOK_RETURN(status=1)

        ! define variable:
        status = NF90_Def_Var( self%ncid, 'area', NF90_REAL, &
                                 (/self%dimid_lon,self%dimid_lat/), varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! attributes:
        status = nf90_put_att( self%ncid, varid, 'standard_name', 'cell_area' )
        IF_NF90_NOTOK_RETURN(status=1)
        status = nf90_put_att( self%ncid, varid, 'units', trim(units) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! grid attributes:
        call glb_ugg%DefCoor_NetCDF( gncd, varid, status )
        IF_NOTOK_RETURN(status=1)

        ! end defintion mode:
        status = NF90_EndDef( self%ncid )
        IF_NF90_NOTOK_RETURN(status=1)

        ! write grid to netCDF file
        call glb_ugg%PutGrid_NetCDF( gncd, status )
        IF_NOTOK_RETURN(status=1)

      end if ! root

      ! select:
      pat = area(1:nx,1:ny,1)
      ! write 2D field, without level, no time index:
      call LE_Output_Put_Var_Domains( self%ncid, varid, -1, -1, pat, status )
      IF_NOTOK_RETURN(status=1)

      ! writen on root only:
      if ( goc%root ) then
        ! close:
        status = NF90_Close( self%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if ! root
#endif

      ! reset flag:
      self%created = .true.

    end if  ! not created yet

    ! ok
    status = 0

  end subroutine LE_Output_Grid_PutOut


end module LE_Output_Grid
