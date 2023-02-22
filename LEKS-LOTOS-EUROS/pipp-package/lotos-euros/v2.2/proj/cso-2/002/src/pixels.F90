!###############################################################################
!
! Pixels - storage for satellite data
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
  
module Pixels

  use GO              , only : gol, goPr, goErr
  use NetCDF          , only : NF90_StrError, NF90_NOERR

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_Pixels_0D_i
  public  ::  T_Pixels_0D
  public  ::  T_Pixels_1D
  public  ::  T_Pixels_2D
  
  public  ::  T_Track_0D
  public  ::  T_Track_1D
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'Pixels'


  ! --- types --------------------------------
  
  ! 1 value per pixel:
  type :: T_Pixels
    ! size:
    integer               ::  npix
    ! attributes:
    character(len=64)     ::  units
    character(len=1024)   ::  long_name
    ! adhoc attributes ...
    character(len=64)     ::  compress
    character(len=1024)   ::  description
    ! netcdf output:
    integer               ::  varid
  contains
    procedure ::                     Pixels_Init
    procedure ::                     Pixels_Done
    procedure :: NcGetAttrs       => Pixels_NcGetAttrs
    procedure :: NcPutAttrs       => Pixels_NcPutAttrs
  end type T_Pixels
  
  ! *
  
  ! 1 integer value per pixel:
  type, extends(T_Pixels) :: T_Pixels_0D_i
    ! storage:
    integer, allocatable     ::  data(:)  ! (npix)
    ! netcdf output:
    integer                  ::  fill_value
  contains
    procedure :: Init            => Pixels_0D_i_Init
    procedure :: Done            => Pixels_0D_i_Done
    procedure :: NcInit          => Pixels_0D_i_NcInit
    procedure :: NcDef           => Pixels_0D_i_NcDef
    procedure :: NcPutGlbSelect  => Pixels_0D_i_NcPutGlbSelect
    procedure :: NcPutGather     => Pixels_0D_i_NcPutGather
  end type T_Pixels_0D_i
  
  ! *
  
  ! 1 real value per pixel:
  type, extends(T_Pixels) :: T_Pixels_0D
    ! storage:
    real, allocatable     ::  data(:)  ! (npix)
    ! netcdf output:
    real                  ::  fill_value
  contains
    procedure :: Init            => Pixels_0D_Init
    procedure :: Done            => Pixels_0D_Done
    procedure :: NcInit          => Pixels_0D_NcInit
    procedure :: NcDef           => Pixels_0D_NcDef
    procedure :: NcPutGlb        => Pixels_0D_NcPutGlb
    procedure :: NcPutGlbSelect  => Pixels_0D_NcPutGlbSelect
    procedure :: NcPutGather     => Pixels_0D_NcPutGather
  end type T_Pixels_0D
  
  ! *
  
  ! vector per pixel:
  type, extends(T_Pixels) :: T_Pixels_1D
    ! size:
    integer               ::  nv
    ! storage:
    real, allocatable     ::  data(:,:)  ! (nv,npix)
    ! netcdf output:
    real                  ::  fill_value
  contains
    procedure :: Init            => Pixels_1D_Init
    procedure :: Done            => Pixels_1D_Done
    procedure :: NcInit          => Pixels_1D_NcInit
    procedure :: NcDef           => Pixels_1D_NcDef
    procedure :: NcPutGlb        => Pixels_1D_NcPutGlb
    procedure :: NcPutGlbSelect  => Pixels_1D_NcPutGlbSelect
    procedure :: NcPutGather     => Pixels_1D_NcPutGather
  end type T_Pixels_1D
  
  ! *
  
  ! matrix per pixel:
  type, extends(T_Pixels) :: T_Pixels_2D
    ! size:
    integer               ::  mv, nv
    ! storage:
    real, allocatable     ::  data(:,:,:)  ! (mv,nv,npix)
    ! netcdf output:
    real                  ::  fill_value
  contains
    procedure :: Init            => Pixels_2D_Init
    procedure :: Done            => Pixels_2D_Done
    procedure :: NcInit          => Pixels_2D_NcInit
    procedure :: NcDef           => Pixels_2D_NcDef
    procedure :: NcPutGlb        => Pixels_2D_NcPutGlb
    procedure :: NcPutGlbSelect  => Pixels_2D_NcPutGlbSelect
    procedure :: NcPutGather     => Pixels_2D_NcPutGather
  end type T_Pixels_2D
  
  
  ! ***
  
  
  ! 2D track with values:
  type :: T_Track
    ! size:
    integer               ::  ntx, nty
    ! attributes:
    character(len=64)     ::  units
    character(len=1024)   ::  long_name
    ! netcdf output:
    integer               ::  varid
  contains
    procedure ::                     Track_Init
    procedure ::                     Track_Done
    procedure :: NcGetAttrs       => Track_NcGetAttrs
    procedure :: NcPutAttrs       => Track_NcPutAttrs
  end type T_Track
  
  ! *
  
  ! 1 real value per pixel:
  type, extends(T_Track) :: T_Track_0D
    ! storage:
    real, allocatable     ::  data(:,:)  ! (ntx,nty)
    ! netcdf output:
    real                  ::  fill_value
  contains
    procedure :: Init            => Track_0D_Init
    procedure :: Done            => Track_0D_Done
    procedure :: NcInit          => Track_0D_NcInit
    procedure :: NcDef           => Track_0D_NcDef
    procedure :: NcPutGlb        => Track_0D_NcPutGlb
  end type T_Track_0D
  
  ! *
  
  ! vector per pixel:
  type, extends(T_Track) :: T_Track_1D
    ! size:
    integer               ::  nv
    ! storage:
    real, allocatable     ::  data(:,:,:)  ! (nv,ntx,nty)
    ! netcdf output:
    real                  ::  fill_value
  contains
    procedure :: Init            => Track_1D_Init
    procedure :: Done            => Track_1D_Done
    procedure :: NcInit          => Track_1D_NcInit
    procedure :: NcDef           => Track_1D_NcDef
    procedure :: NcPutGlb        => Track_1D_NcPutGlb
  end type T_Track_1D
  

  
contains


  ! ====================================================================
  ! ===
  ! === Pixels
  ! ===
  ! ====================================================================


  subroutine Pixels_Init( self, npix, status, &
                             long_name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels), intent(out)      ::  self
    integer, intent(in)               ::  npix
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  long_name
    character(len=*), intent(in), optional    ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%npix = npix
    
    ! defaults:
    self%long_name = ''
    self%units = ''
    
    ! attributes?
    if ( present(long_name) ) self%long_name = long_name
    if ( present(units) ) self%units = units
    
    ! ok
    status = 0
    
  end subroutine Pixels_Init


  ! ***


  subroutine Pixels_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Pixels), intent(inout)      ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Pixels_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! empty:
    self%npix = 0

    ! ok
    status = 0
    
  end subroutine Pixels_Done
  
  
  ! ***


  subroutine Pixels_NcGetAttrs( self, ncid, varid, status )
  
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels), intent(inout)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  varid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_NcGetAttrs'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! read attributes:
    status = NF90_Get_Att( ncid, varid, 'units', self%units )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Att( ncid, varid, 'long_name', self%long_name )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! known other attribute, to be generalized ...
    status = NF90_Get_Att( ncid, varid, 'compress', self%compress )
    if ( status == NF90_ENOTATT ) then
      self%compress = ''
    else if ( status /= NF90_NOERR ) then
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! known other attribute, to be generalized ...
    status = NF90_Get_Att( ncid, varid, 'description', self%description )
    if ( status == NF90_ENOTATT ) then
      self%description = ''
    else if ( status /= NF90_NOERR ) then
      gol=NF90_StrError(status); call goErr
      TRACEBACK; status=1; return
    end if    
    
    ! ok
    status = 0
    
  end subroutine Pixels_NcGetAttrs
  
  
  ! ***


  subroutine Pixels_NcPutAttrs( self, ncid, varid, status )
  
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels), intent(in)       ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  varid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_NcPutAttrs'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! define attribute:
    status = NF90_Put_Att( ncid, varid, 'units', self%units )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid, 'long_name', self%long_name )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! optional ...
    if ( len_trim(self%compress) > 0 ) then
      status = NF90_Put_Att( ncid, varid, 'compress', self%compress )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    if ( len_trim(self%description) > 0 ) then
      status = NF90_Put_Att( ncid, varid, 'description', self%description )
      IF_NF90_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine Pixels_NcPutAttrs
  
  


  ! ====================================================================
  ! ===
  ! === Pixels (integer values)
  ! ===
  ! ====================================================================


  subroutine Pixels_0D_i_Init( self, npix, status, &
                                  source, long_name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D_i), intent(out)   ::  self
    integer, intent(in)                 ::  npix
    integer, intent(out)                ::  status
    
    integer, intent(in), optional             ::  source
    character(len=*), intent(in), optional    ::  long_name
    character(len=*), intent(in), optional    ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_i_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init base class:
    call self%Pixels_Init( npix, status, long_name=long_name, units=units )
    IF_NOTOK_RETURN(status=1)

    ! fill value:
    self%fill_value = - huge(1)
    ! storage:
    allocate( self%data(max(1,self%npix)), source=self%fill_value, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! replace if necessary:
    if ( present(source) ) self%data = source
    
    ! ok
    status = 0
    
  end subroutine Pixels_0D_i_Init


  ! ***


  subroutine Pixels_0D_i_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D_i), intent(inout)   ::  self
    integer, intent(out)                  ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Pixels_0D_i_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%data, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! empty:
    self%npix = 0

    ! ok
    status = 0
    
  end subroutine Pixels_0D_i_Done
  
  
  ! ***


  subroutine Pixels_0D_i_NcInit( self, ncid, varname, status )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D_i), intent(out)   ::  self
    integer, intent(in)                 ::  ncid
    character(len=*), intent(in)        ::  varname
    integer, intent(out)                ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_i_NcInit'
    
    integer, parameter     ::  ndim = 1
    
    ! --- local ----------------------------------
    
    integer     ::  varid
    integer     ::  dimids(ndim)
    integer     ::  shp(ndim)
    integer     ::  idim
    
    ! --- begin ----------------------------------

    ! get variable id:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, dimids=dimids )
    IF_NF90_NOTOK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, ndim
      ! get size:
      status = NF90_Inquire_Dimension( ncid, dimids(idim), len=shp(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do

    ! init storage:
    call self%Init( shp(1), status )
    IF_NOTOK_RETURN(status=1)
    
    ! read attributes:
    call self%NcGetAttrs( ncid, varid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read:
    status = NF90_Get_Var( ncid, varid, self%data )
    IF_NF90_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Pixels_0D_i_NcInit


  ! ***


  subroutine Pixels_0D_i_NcDef( self, ncid, varname, dimids, status )
  
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_INT
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D_i), intent(inout)   ::  self
    integer, intent(in)                   ::  ncid
    character(len=*), intent(in)          ::  varname
    integer, intent(in)                   ::  dimids(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_i_NcDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! define variable:
    status = NF90_Def_Var( ncid, trim(varname), NF90_INT, dimids, self%varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! add attributes:
    call self%NcPutAttrs( ncid, self%varid, status )
    IF_NOTOK_RETURN(status=1)
    ! define attribute:
    status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Pixels_0D_i_NcDef
  
  
  ! ***


  ! write from root, this is a global array
  
  subroutine Pixels_0D_i_NcPutGlb( self, ncid, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D_i), intent(in)    ::  self
    integer, intent(in)                 ::  ncid
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_i_NcPutGlb'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
    
      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, self%data )
      IF_NF90_NOTOK_RETURN(status=1)
        
    end if ! root

    ! ok
    status = 0
    
  end subroutine Pixels_0D_i_NcPutGlb
  
  
  ! ***


  ! write from root, selection from global array
  
  subroutine Pixels_0D_i_NcPutGlbSelect( self, ncid, mapping, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D_i), intent(in)    ::  self
    integer, intent(in)                 ::  ncid
    integer, intent(in)                 ::  mapping(:)  ! (nglb) target index, or <0 if not used
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_i_NcPutGlbSelect'
    
    ! --- local ----------------------------------

    integer                 ::  nglb
    integer                 ::  n
    integer, allocatable    ::  data(:)  ! (n)
    integer                 ::  k
    integer                 ::  iglb

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
    
      ! number of global pixels:
      nglb = size(mapping)
      ! number of selected pixels:
      n = maxval(mapping)
      ! target array:
      allocate( data(n), source=self%fill_value, stat=status )
      IF_NOTOK_RETURN(status=1)

      ! loop over global pixels:
      do iglb = 1, nglb
        ! target index:
        k = mapping(iglb)
        ! defined?
        if ( k > 0 ) then
          ! copy:
          data(k) = self%data(iglb)
        end if
      end do ! glb

      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, data )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( data, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! root

    ! ok
    status = 0
    
  end subroutine Pixels_0D_i_NcPutGlbSelect
  
  
  ! ***


  ! gather on root using mapping weights and write variable
  
  subroutine Pixels_0D_i_NcPutGather( self, ncid, mapping, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D_i), intent(in)  ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  mapping(:)   ! (npix_all) target index
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_i_NcPutGather'
    
    ! --- local ----------------------------------

    integer                 ::  npix_all
    integer, allocatable    ::  data_all(:)      ! (npix_all)
    integer                 ::  nout
    integer                 ::  iout
    integer, allocatable    ::  data(:)          ! (n)
    integer                 ::  iall
    
    ! --- begin ----------------------------------
    
    ! gather on root only ..
    if ( goc%root ) then
      ! count:
      npix_all = size(mapping)
      ! storage for all data:
      allocate( data_all(npix_all), stat=status )
      IF_NOTOK_RETURN(status=1)
    else
      ! dummy ...
      allocate( data_all(1), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if ! root

    ! gather data arrays on root, local size might be zero
    call goc%GatherV( self%data, data_all, status, nloc=self%npix )
    IF_NOTOK_RETURN(status=1)

    ! gather on root only ..
    if ( goc%root ) then

      ! number of selected pixels:
      nout = maxval(mapping)
      ! target array:
      allocate( data(nout), source=self%fill_value, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over collected pixels:
      do iall = 1, npix_all
        ! target index:
        iout = mapping(iall)
        !! add contribution (IN FUTURE: APPLY WEIGHTS!)
        !data(:,iout) = data(:,iout) + data_all(:,iall)
        ! just copy, might overwrite previous copy:
        data(iout) = data_all(iall)
      end do ! iall

      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, data )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( data, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! root

    ! clear:
    deallocate( data_all, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Pixels_0D_i_NcPutGather


  ! ====================================================================
  ! ===
  ! === Pixels (real values)
  ! ===
  ! ====================================================================


  subroutine Pixels_0D_Init( self, npix, status, &
                                  source, long_name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D), intent(out)   ::  self
    integer, intent(in)               ::  npix
    integer, intent(out)              ::  status
    
    real, intent(in), optional                ::  source
    character(len=*), intent(in), optional    ::  long_name
    character(len=*), intent(in), optional    ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init base class:
    call self%Pixels_Init( npix, status, long_name=long_name, units=units )
    IF_NOTOK_RETURN(status=1)

    ! fill value:
    self%fill_value = - huge(1.0)
    ! storage:
    allocate( self%data(max(1,self%npix)), source=self%fill_value, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! replace if necessary:
    if ( present(source) ) self%data = source
    
    ! ok
    status = 0
    
  end subroutine Pixels_0D_Init


  ! ***


  subroutine Pixels_0D_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D), intent(inout)   ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Pixels_0D_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%data, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! empty:
    self%npix = 0

    ! ok
    status = 0
    
  end subroutine Pixels_0D_Done
  
  
  ! ***


  subroutine Pixels_0D_NcInit( self, ncid, varname, status )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D), intent(out)   ::  self
    integer, intent(in)               ::  ncid
    character(len=*), intent(in)      ::  varname
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_NcInit'
    
    integer, parameter     ::  ndim = 1
    
    ! --- local ----------------------------------
    
    integer     ::  varid
    integer     ::  dimids(ndim)
    integer     ::  shp(ndim)
    integer     ::  idim
    
    ! --- begin ----------------------------------

    ! get variable id:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, dimids=dimids )
    IF_NF90_NOTOK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, ndim
      ! get size:
      status = NF90_Inquire_Dimension( ncid, dimids(idim), len=shp(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do

    ! init storage:
    call self%Init( shp(1), status )
    IF_NOTOK_RETURN(status=1)

    ! read attributes:
    call self%NcGetAttrs( ncid, varid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read:
    status = NF90_Get_Var( ncid, varid, self%data )
    IF_NF90_NOTOK_RETURN(status=1)
    
      !! testing ...
      !write (gol,*) rname//': long_name = ', trim(self%long_name), ' ; varname = ', trim(varname); call goPr
      !write (gol,*) rname//':   varid = ', varid; call goPr
      !write (gol,*) rname//':   shape = ', shape(self%data); call goPr
      !write (gol,*) rname//':   range = ', minval(self%data), maxval(self%data); call goPr
      
    ! ok
    status = 0
    
  end subroutine Pixels_0D_NcInit


  ! ***


  subroutine Pixels_0D_NcDef( self, ncid, varname, dimids, status )
  
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_FLOAT
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D), intent(inout)     ::  self
    integer, intent(in)                   ::  ncid
    character(len=*), intent(in)          ::  varname
    integer, intent(in)                   ::  dimids(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_NcDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! define variable:
    status = NF90_Def_Var( ncid, trim(varname), NF90_FLOAT, dimids, self%varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! add attributes:
    call self%NcPutAttrs( ncid, self%varid, status )
    IF_NOTOK_RETURN(status=1)
    ! define attribute:
    status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Pixels_0D_NcDef
  
  
  ! ***


  ! write from root, this is a global array
  
  subroutine Pixels_0D_NcPutGlb( self, ncid, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_NcPutGlb'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, self%data )
      IF_NF90_NOTOK_RETURN(status=1)
    end if ! root

    ! ok
    status = 0
    
  end subroutine Pixels_0D_NcPutGlb
  
  
  ! ***


  ! write from root, selection from global array
  
  subroutine Pixels_0D_NcPutGlbSelect( self, ncid, mapping, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D), intent(in)      ::  self
    integer, intent(in)                 ::  ncid
    integer, intent(in)                 ::  mapping(:)  ! (nglb) target index, or <0 if not used
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_NcPutGlbSelect'
    
    ! --- local ----------------------------------

    integer                 ::  nglb
    integer                 ::  n
    real, allocatable       ::  data(:)  ! (n)
    integer                 ::  k
    integer                 ::  iglb

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
    
      ! number of global pixels:
      nglb = size(mapping)
      ! number of selected pixels:
      n = maxval(mapping)
      ! target array:
      allocate( data(n), source=self%fill_value, stat=status )
      IF_NOTOK_RETURN(status=1)

      ! loop over global pixels:
      do iglb = 1, nglb
        ! target index:
        k = mapping(iglb)
        ! defined?
        if ( k > 0 ) then
          ! copy:
          data(k) = self%data(iglb)
        end if
      end do ! glb

      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, data )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( data, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! root

    ! ok
    status = 0
    
  end subroutine Pixels_0D_NcPutGlbSelect
  
  
  ! ***


  ! gather on root using mapping weights and write variable
  
  subroutine Pixels_0D_NcPutGather( self, ncid, mapping, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_0D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  mapping(:)   ! (npix_all) target index
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_0D_NcPutGather'
    
    ! --- local ----------------------------------

    integer                 ::  npix_all
    real, allocatable       ::  data_all(:)      ! (npix_all)
    integer                 ::  nout
    integer                 ::  iout
    real, allocatable       ::  data(:)          ! (n)
    integer                 ::  iall
    
    ! --- begin ----------------------------------
    
    ! gather on root only ..
    if ( goc%root ) then
      ! count:
      npix_all = size(mapping)
      ! storage for all data:
      allocate( data_all(npix_all), stat=status )
      IF_NOTOK_RETURN(status=1)
    else
      ! dummy ...
      allocate( data_all(1), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if ! root

    ! gather data arrays on root, local size might be zero
    call goc%GatherV( self%data, data_all, status, nloc=self%npix )
    IF_NOTOK_RETURN(status=1)

    ! gather on root only ..
    if ( goc%root ) then

      ! number of selected pixels:
      nout = maxval(mapping)
      ! target array:
      allocate( data(nout), source=self%fill_value, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over collected pixels:
      do iall = 1, npix_all
        ! target index:
        iout = mapping(iall)
        !! add contribution (IN FUTURE: APPLY WEIGHTS!)
        !data(:,iout) = data(:,iout) + data_all(:,iall)
        ! just copy, might overwrite previous copy:
        data(iout) = data_all(iall)
      end do ! iall

      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, data )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( data, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! root

    ! clear:
    deallocate( data_all, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Pixels_0D_NcPutGather


  ! ====================================================================
  ! ===
  ! === Pixels 1D
  ! ===
  ! ====================================================================


  subroutine Pixels_1D_Init( self, shp, status, &
                                  source, long_name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_1D), intent(out)   ::  self
    integer, intent(in)               ::  shp(2)  ! (/nv,npix/)
    integer, intent(out)              ::  status
    
    real, intent(in), optional                ::  source
    character(len=*), intent(in), optional    ::  long_name
    character(len=*), intent(in), optional    ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_1D_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init base class:
    call self%Pixels_Init( shp(2), status, long_name=long_name, units=units )
    IF_NOTOK_RETURN(status=1)
    
    ! store:
    self%nv   = shp(1)

    ! fill value:
    self%fill_value = - huge(1.0)
    ! storage:
    allocate( self%data(self%nv,max(1,self%npix)), source=self%fill_value, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! replace if necessary:
    if ( present(source) ) self%data = source
    
    ! ok
    status = 0
    
  end subroutine Pixels_1D_Init


  ! ***


  subroutine Pixels_1D_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Pixels_1D), intent(inout)   ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Pixels_1D_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%data, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! empty:
    self%npix = 0
    self%nv   = 0

    ! ok
    status = 0
    
  end subroutine Pixels_1D_Done
  
  
  ! ***


  subroutine Pixels_1D_NcInit( self, ncid, varname, status, &
                                  nselect, select )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_1D), intent(out)   ::  self
    integer, intent(in)               ::  ncid
    character(len=*), intent(in)      ::  varname
    integer, intent(out)              ::  status
    
    integer, intent(in), optional     ::  nselect    ! npix or 0
    integer, intent(in), optional     ::  select(:)  ! (npix) indices in global arrays

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_1D_NcInit'
    
    integer, parameter     ::  ndim = 2
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    integer             ::  nd
    integer             ::  dimids(ndim)
    integer             ::  shp(ndim)
    integer             ::  idim
    integer             ::  npix
    integer             ::  ipix
    real, allocatable   ::  data(:,:)  ! (nv,npix)
    
    ! --- begin ----------------------------------

    ! get variable id:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, ndims=nd, dimids=dimids )
    IF_NF90_NOTOK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, ndim
      ! scalar?
      if ( idim <= ndim-nd ) then
        shp(idim) = 1
      else
        ! get size:
        status = NF90_Inquire_Dimension( ncid, dimids(idim-(ndim-nd)), len=shp(idim) )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
    end do
    
    ! selection?
    if ( present(select) ) then
      ! number of selected pixels:
      if ( present(nselect) ) then
        ! check ...
        if ( (nselect < 0) .or. (nselect > size(select)) ) then
          write (gol,'("nselect=",i0," while size(select)=",i0)') nselect, size(select); call goErr
          TRACEBACK; status=1; return
        end if
        ! copy:
        npix = nselect
      else
        ! full size:
        npix = size(select)
      end if
      
      ! init storage, npix could be 0:
      call self%Init( (/shp(1),npix/), status )
      IF_NOTOK_RETURN(status=1)

      ! read attributes:
      call self%NcGetAttrs( ncid, varid, status )
      IF_NOTOK_RETURN(status=1)

      ! any local pixels?
      if ( npix > 0 ) then

        ! storage for all data:
        allocate( data(shp(1),shp(2)), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! read, seems needed to read exact number of dimensions ...
        if ( nd == 1 ) then
          status = NF90_Get_Var( ncid, varid, data(1,:) )
          IF_NF90_NOTOK_RETURN(status=1)
        else
          status = NF90_Get_Var( ncid, varid, data(:,:) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if

        !! testing...
        !write (gol,*) rname//': long_name = ', trim(self%long_name), ' ; varname = ', trim(varname); call goPr
        !write (gol,*) rname//':   range = ', minval(data), maxval(data); call goPr
        !write (gol,*) rname//':   shape = ', shape(data); call goPr
        !write (gol,*) rname//':   varid = ', varid; call goPr

        ! check ..
        if ( any(select < 1) .or. any(select > shp(ndim)) ) then
          write (gol,'("selection indices in range ",i0,",..,",i0," while number of pixels is ",i0)') &
                   minval(select), maxval(select), shp(ndim); call goErr
          TRACEBACK; status=1; return
        end if

        ! loop over selected pixels:
        do ipix = 1, npix
          ! copy value from global array:
          self%data(:,ipix) = data(:,select(ipix))
        end do

        ! clear:
        deallocate( data, stat=status )
        IF_NOTOK_RETURN(status=1)
      
      end if ! npix > 0
    
    else

      ! init storage for all pixels:
      call self%Init( shp, status )
      IF_NOTOK_RETURN(status=1)
      ! read attributes:
      call self%NcGetAttrs( ncid, varid, status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      status = NF90_Get_Var( ncid, varid, self%data )
      IF_NF90_NOTOK_RETURN(status=1)
          
    end if

    ! ok
    status = 0
    
  end subroutine Pixels_1D_NcInit


  ! ***


  subroutine Pixels_1D_NcDef( self, ncid, varname, dimids, status )
  
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_FLOAT
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_1D), intent(inout)     ::  self
    integer, intent(in)                   ::  ncid
    character(len=*), intent(in)          ::  varname
    integer, intent(in)                   ::  dimids(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_1D_NcDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! define variable:
    status = NF90_Def_Var( ncid, trim(varname), NF90_FLOAT, dimids, self%varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! add attributes:
    call self%NcPutAttrs( ncid, self%varid, status )
    IF_NOTOK_RETURN(status=1)
    ! define attribute:
    status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Pixels_1D_NcDef
  
  
  ! ***


  ! write global array from root, each pe has same copy
  
  subroutine Pixels_1D_NcPutGlb( self, ncid, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_1D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_1D_NcPutGlb'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, self%data )
      IF_NF90_NOTOK_RETURN(status=1)
    end if ! root
    
    ! ok
    status = 0
    
  end subroutine Pixels_1D_NcPutGlb
  
  
  ! ***


  ! write from root, selection from global array
  
  subroutine Pixels_1D_NcPutGlbSelect( self, ncid, mapping, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_1D), intent(in)      ::  self
    integer, intent(in)                 ::  ncid
    integer, intent(in)                 ::  mapping(:)  ! (nglb) target index, or <0 if not used
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_1D_NcPutGlbSelect'
    
    ! --- local ----------------------------------

    integer                 ::  nglb
    integer                 ::  n
    real, allocatable       ::  data(:,:)  ! (nv,n)
    integer                 ::  k
    integer                 ::  iglb

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
    
      ! number of global pixels:
      nglb = size(mapping)
      ! number of selected pixels:
      n = maxval(mapping)
      ! target array:
      allocate( data(self%nv,n), source=self%fill_value, stat=status )
      IF_NOTOK_RETURN(status=1)

      ! loop over global pixels:
      do iglb = 1, nglb
        ! target index:
        k = mapping(iglb)
        ! defined?
        if ( k > 0 ) then
          ! copy:
          data(:,k) = self%data(:,iglb)
        end if
      end do ! glb

      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, data )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( data, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! root

    ! ok
    status = 0
    
  end subroutine Pixels_1D_NcPutGlbSelect
  
  
  ! ***


  ! gather on root using mapping weights and write variable
  
  subroutine Pixels_1D_NcPutGather( self, ncid, mapping, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_1D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  mapping(:)   ! (npix_all) target index
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_1D_NcPutGather'
    
    ! --- local ----------------------------------

    integer                 ::  npix_all
    real, allocatable       ::  data_all(:,:)      ! (nv,npix_all)
    integer                 ::  nout
    integer                 ::  iout
    real, allocatable       ::  data(:,:)          ! (nv,n)
    integer                 ::  iall
    
    ! --- begin ----------------------------------
    
    ! gather on root only ..
    if ( goc%root ) then
      ! count:
      npix_all = size(mapping)
      ! storage for all data:
      allocate( data_all(self%nv,npix_all), stat=status )
      IF_NOTOK_RETURN(status=1)
    else
      ! dummy ...
      allocate( data_all(1,1), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if ! root

    ! gather data arrays on root, local size might be zero
    call goc%GatherV( self%data, data_all, status, nloc=self%npix )
    IF_NOTOK_RETURN(status=1)

    ! gather on root only ..
    if ( goc%root ) then

      ! number of selected pixels:
      nout = maxval(mapping)
      ! target array:
      allocate( data(self%nv,nout), source=self%fill_value, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over collected pixels:
      do iall = 1, npix_all
        ! target index:
        iout = mapping(iall)
        !! add contribution (IN FUTURE: APPLY WEIGHTS!)
        !data(:,iout) = data(:,iout) + data_all(:,iall)
        ! just copy, might overwrite previous copy:
        data(:,iout) = data_all(:,iall)
      end do ! iall

      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, data )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( data, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! root

    ! clear:
    deallocate( data_all, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Pixels_1D_NcPutGather


  ! ====================================================================
  ! ===
  ! === Pixels 2D
  ! ===
  ! ====================================================================


  subroutine Pixels_2D_Init( self, shp, status, &
                                  source, long_name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_2D), intent(out)   ::  self
    integer, intent(in)               ::  shp(3)  ! (/mv,nv,npix/)
    integer, intent(out)              ::  status
    
    real, intent(in), optional                ::  source
    character(len=*), intent(in), optional    ::  long_name
    character(len=*), intent(in), optional    ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_2D_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init base class:
    call self%Pixels_Init( shp(3), status, long_name=long_name, units=units )
    IF_NOTOK_RETURN(status=1)
    
    ! store:
    self%mv   = shp(1)
    self%nv   = shp(2)

    ! fill value:
    self%fill_value = - huge(1.0)
    ! storage:
    allocate( self%data(self%mv,self%nv,max(1,self%npix)), source=self%fill_value, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! replace if necessary:
    if ( present(source) ) self%data = source
    
    ! ok
    status = 0
    
  end subroutine Pixels_2D_Init


  ! ***


  subroutine Pixels_2D_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Pixels_2D), intent(inout)   ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Pixels_2D_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%data, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! empty:
    self%npix = 0
    self%nv   = 0
    self%mv   = 0

    ! ok
    status = 0
    
  end subroutine Pixels_2D_Done
  
  
  ! ***


  subroutine Pixels_2D_NcInit( self, ncid, varname, status, &
                                  nselect, select )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_2D), intent(out)   ::  self
    integer, intent(in)               ::  ncid
    character(len=*), intent(in)      ::  varname
    integer, intent(out)              ::  status
    
    integer, intent(in), optional     ::  nselect    ! npix or 0
    integer, intent(in), optional     ::  select(:)  ! (npix) indices in global arrays

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_2D_NcInit'
    
    integer, parameter     ::  ndim = 3
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    integer             ::  nd
    integer             ::  dimids(ndim)
    integer             ::  shp(ndim)
    integer             ::  idim
    integer             ::  npix
    integer             ::  ipix
    real, allocatable   ::  data(:,:,:)  ! (mv,nv,npix)
    
    ! --- begin ----------------------------------

    ! get variable id:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, ndims=nd, dimids=dimids )
    IF_NF90_NOTOK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, ndim
      ! scalar?
      if ( idim <= ndim-nd ) then
        shp(idim) = 1
      else
        ! get size:
        status = NF90_Inquire_Dimension( ncid, dimids(idim-(ndim-nd)), len=shp(idim) )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
    end do

    ! selection?
    if ( present(select) ) then
      ! number of selected pixels:
      if ( present(nselect) ) then
        ! check ...
        if ( (nselect < 0) .or. (nselect > size(select)) ) then
          write (gol,'("nselect=",i0," while size(select)=",i0)') nselect, size(select); call goErr
          TRACEBACK; status=1; return
        end if
        ! copy:
        npix = nselect
      else
        ! full size:
        npix = size(select)
      end if

      ! init storage, npix could be 0:
      call self%Init( (/shp(1),shp(2),npix/), status )
      IF_NOTOK_RETURN(status=1)

      ! read attributes:
      call self%NcGetAttrs( ncid, varid, status )
      IF_NOTOK_RETURN(status=1)

      ! any local pixels?
      if ( npix > 0 ) then
    
        ! number of selected pixels:
        npix = size(select)

        ! storage for all data:
        allocate( data(shp(1),shp(2),shp(3)), stat=status )
        IF_NOTOK_RETURN(status=1)

        ! read, seems needed to read exact number of dimensions ...
        if ( nd == 1 ) then
          status = NF90_Get_Var( ncid, varid, data(1,1,:) )
          IF_NF90_NOTOK_RETURN(status=1)
        else if ( nd == 2 ) then
          status = NF90_Get_Var( ncid, varid, data(1,:,:) )
          IF_NF90_NOTOK_RETURN(status=1)
        else
          status = NF90_Get_Var( ncid, varid, data(:,:,:) )
          IF_NF90_NOTOK_RETURN(status=1)
        end if

        !! testing ...
        !write (gol,*) rname//': long_name = ', trim(self%long_name), ' ; varname = ', trim(varname); call goPr
        !write (gol,*) rname//':   varid = ', varid; call goPr
        !write (gol,*) rname//':   shape = ', shape(data); call goPr
        !do ipix = 1, shp(3)
        !  write (gol,*) rname//':   pix ', ipix, ':', data(:,:,ipix); call goPr
        !end do
        write (gol,*) rname//':   range = ', minval(data), maxval(data); call goPr

        ! check ..
        if ( any(select < 1) .or. any(select > shp(ndim)) ) then
          write (gol,'("selection indices in range ",i0,",..,",i0," while number of pixels is ",i0)') &
                   minval(select), maxval(select), shp(ndim); call goErr
          TRACEBACK; status=1; return
        end if
        
        ! loop over selected pixels:
        do ipix = 1, npix
          ! copy value from global array:
          self%data(:,:,ipix) = data(:,:,select(ipix))
        end do

        ! clear:
        deallocate( data, stat=status )
        IF_NOTOK_RETURN(status=1)
        
      end if  ! npix > 0
    
    else

      ! init storage for all pixels:
      call self%Init( shp, status )
      IF_NOTOK_RETURN(status=1)

      ! read attributes:
      call self%NcGetAttrs( ncid, varid, status )
      IF_NOTOK_RETURN(status=1)

      ! read:
      status = NF90_Get_Var( ncid, varid, self%data )
      IF_NF90_NOTOK_RETURN(status=1)
      
    end if

    ! ok
    status = 0
    
  end subroutine Pixels_2D_NcInit


  ! ***


  subroutine Pixels_2D_NcDef( self, ncid, varname, dimids, status )
  
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_FLOAT
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_2D), intent(inout)     ::  self
    integer, intent(in)                   ::  ncid
    character(len=*), intent(in)          ::  varname
    integer, intent(in)                   ::  dimids(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_2D_NcDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! define variable:
    status = NF90_Def_Var( ncid, trim(varname), NF90_FLOAT, dimids, self%varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! add attributes:
    call self%NcPutAttrs( ncid, self%varid, status )
    IF_NOTOK_RETURN(status=1)
    ! define attribute:
    status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Pixels_2D_NcDef
  
  
  ! ***


  ! write global array from root, each pe has same copy
  
  subroutine Pixels_2D_NcPutGlb( self, ncid, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_2D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_2D_NcPutGlb'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, self%data )
      IF_NF90_NOTOK_RETURN(status=1)
    end if ! root
    
    ! ok
    status = 0
    
  end subroutine Pixels_2D_NcPutGlb
  
  
  ! ***


  ! write from root, selection from global array
  
  subroutine Pixels_2D_NcPutGlbSelect( self, ncid, mapping, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_2D), intent(in)      ::  self
    integer, intent(in)                 ::  ncid
    integer, intent(in)                 ::  mapping(:)  ! (nglb) target index, or <0 if not used
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_2D_NcPutGlbSelect'
    
    ! --- local ----------------------------------

    integer                 ::  nglb
    integer                 ::  n
    real, allocatable       ::  data(:,:,:)  ! (mv,nv,n)
    integer                 ::  k
    integer                 ::  iglb

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
    
      ! number of global pixels:
      nglb = size(mapping)
      ! number of selected pixels:
      n = maxval(mapping)
      ! target array:
      allocate( data(self%mv,self%nv,n), source=self%fill_value, stat=status )
      IF_NOTOK_RETURN(status=1)

      ! loop over global pixels:
      do iglb = 1, nglb
        ! target index:
        k = mapping(iglb)
        ! defined?
        if ( k > 0 ) then
          ! copy:
          data(:,:,k) = self%data(:,:,iglb)
        end if
      end do ! glb

      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, data )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( data, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! root

    ! ok
    status = 0
    
  end subroutine Pixels_2D_NcPutGlbSelect
  
  
  ! ***


  ! gather on root using mapping weights and write variable
  
  subroutine Pixels_2D_NcPutGather( self, ncid, mapping, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Pixels_2D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  mapping(:)   ! (npix_all) target index
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Pixels_2D_NcPutGather'
    
    ! --- local ----------------------------------

    integer                 ::  npix_all
    real, allocatable       ::  data_all(:,:,:)      ! (mv,nv,npix_all)
    integer                 ::  nout
    integer                 ::  iout
    real, allocatable       ::  data(:,:,:)          ! (mv,nv,n)
    integer                 ::  iall
    
    ! --- begin ----------------------------------
    
    ! gather on root only ..
    if ( goc%root ) then
      ! count:
      npix_all = size(mapping)
      ! storage for all data:
      allocate( data_all(self%mv,self%nv,npix_all), stat=status )
      IF_NOTOK_RETURN(status=1)
    else
      ! dummy ...
      allocate( data_all(1,1,1), stat=status )
      IF_NOTOK_RETURN(status=1)
    end if ! root

    ! gather data arrays on root, local size might be zero
    call goc%GatherV( self%data, data_all, status, nloc=self%npix )
    IF_NOTOK_RETURN(status=1)

    ! gather on root only ..
    if ( goc%root ) then

      ! number of selected pixels:
      nout = maxval(mapping)
      ! target array:
      allocate( data(self%mv,self%nv,nout), source=self%fill_value, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over collected pixels:
      do iall = 1, npix_all
        ! target index:
        iout = mapping(iall)
        !! add contribution (IN FUTURE: APPLY WEIGHTS!)
        !data(:,iout) = data(:,iout) + data_all(:,iall)
        ! just copy, might overwrite previous copy:
        data(:,:,iout) = data_all(:,:,iall)
      end do ! iall

      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, data )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( data, stat=status )
      IF_NOTOK_RETURN(status=1)

    end if ! root

    ! clear:
    deallocate( data_all, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Pixels_2D_NcPutGather


  ! ##########################################################################################
  ! ###
  ! ### Tracks
  ! ###
  ! ##########################################################################################


  subroutine Track_Init( self, ntx, nty, status, &
                             long_name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_Track), intent(out)       ::  self
    integer, intent(in)               ::  ntx, nty
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  long_name
    character(len=*), intent(in), optional    ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%ntx = ntx
    self%nty = nty
    
    ! defaults:
    self%long_name = ''
    self%units = ''
    
    ! attributes?
    if ( present(long_name) ) self%long_name = long_name
    if ( present(units) ) self%units = units
    
    ! ok
    status = 0
    
  end subroutine Track_Init


  ! ***


  subroutine Track_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Track), intent(inout)      ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Track_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! empty:
    self%ntx = 0
    self%nty = 0

    ! ok
    status = 0
    
  end subroutine Track_Done
  
  
  ! ***


  subroutine Track_NcGetAttrs( self, ncid, varid, status )
  
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_Track), intent(inout)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  varid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_NcGetAttrs'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! read attribute:
    status = NF90_Get_Att( ncid, varid, 'units', self%units )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Att( ncid, varid, 'long_name', self%long_name )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Track_NcGetAttrs
  
  
  ! ***


  subroutine Track_NcPutAttrs( self, ncid, varid, status )
  
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    class(T_Track), intent(in)       ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  varid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_NcPutAttrs'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! define attribute:
    status = NF90_Put_Att( ncid, varid, 'units', self%units )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid, 'long_name', self%long_name )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Track_NcPutAttrs


  ! ====================================================================
  ! ===
  ! === Track (1 value per pixel)
  ! ===
  ! ====================================================================


  subroutine Track_0D_Init( self, ntx, nty, status, &
                                  source, long_name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_Track_0D), intent(out)    ::  self
    integer, intent(in)               ::  ntx, nty
    integer, intent(out)              ::  status
    
    real, intent(in), optional                ::  source
    character(len=*), intent(in), optional    ::  long_name
    character(len=*), intent(in), optional    ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_0D_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init base class:
    call self%Track_Init( ntx, nty, status, long_name=long_name, units=units )
    IF_NOTOK_RETURN(status=1)

    ! fill value:
    self%fill_value = - huge(1.0)
    ! storage:
    allocate( self%data(self%ntx,self%nty), source=self%fill_value, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! replace if necessary:
    if ( present(source) ) self%data = source
    
    ! ok
    status = 0
    
  end subroutine Track_0D_Init


  ! ***


  subroutine Track_0D_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Track_0D), intent(inout)    ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Track_0D_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%data, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! empty:
    self%ntx = 0
    self%nty = 0

    ! ok
    status = 0
    
  end subroutine Track_0D_Done
  
  
  ! ***


  subroutine Track_0D_NcInit( self, ncid, varname, status )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
  
    ! --- in/out ---------------------------------
    
    class(T_Track_0D), intent(out)   ::  self
    integer, intent(in)               ::  ncid
    character(len=*), intent(in)      ::  varname
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_0D_NcInit'
    
    integer, parameter     ::  ndim = 2
    
    ! --- local ----------------------------------
    
    integer     ::  varid
    integer     ::  dimids(ndim)
    integer     ::  shp(ndim)
    integer     ::  idim
    
    ! --- begin ----------------------------------

    ! get variable id:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, dimids=dimids )
    IF_NF90_NOTOK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, ndim
      ! get size:
      status = NF90_Inquire_Dimension( ncid, dimids(idim), len=shp(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do

    ! init storage:
    call self%Init( shp(1), shp(2), status )
    IF_NOTOK_RETURN(status=1)

    ! read attributes:
    call self%NcGetAttrs( ncid, varid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read:
    status = NF90_Get_Var( ncid, varid, self%data )
    IF_NF90_NOTOK_RETURN(status=1)
    
      !! testing ...
      !write (gol,*) rname//': long_name = ', trim(self%long_name), ' ; varname = ', trim(varname); call goPr
      !write (gol,*) rname//':   varid = ', varid; call goPr
      !write (gol,*) rname//':   shape = ', shape(self%data); call goPr
      !write (gol,*) rname//':   range = ', minval(self%data), maxval(self%data); call goPr
      
    ! ok
    status = 0
    
  end subroutine Track_0D_NcInit


  ! ***


  subroutine Track_0D_NcDef( self, ncid, varname, dimids, status )
  
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_FLOAT
  
    ! --- in/out ---------------------------------
    
    class(T_Track_0D), intent(inout)      ::  self
    integer, intent(in)                   ::  ncid
    character(len=*), intent(in)          ::  varname
    integer, intent(in)                   ::  dimids(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_0D_NcDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! define variable:
    status = NF90_Def_Var( ncid, trim(varname), NF90_FLOAT, dimids, self%varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! add attributes:
    call self%NcPutAttrs( ncid, self%varid, status )
    IF_NOTOK_RETURN(status=1)
    ! define attribute:
    status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Track_0D_NcDef
  
  
  ! ***


  ! write global array from root, each pe has same copy
  
  subroutine Track_0D_NcPutGlb( self, ncid, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Track_0D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_0D_NcPutGlb'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, self%data )
      IF_NF90_NOTOK_RETURN(status=1)
    end if ! root

    ! ok
    status = 0
    
  end subroutine Track_0D_NcPutGlb


  ! ====================================================================
  ! ===
  ! === Track 1D
  ! ===
  ! ====================================================================


  subroutine Track_1D_Init( self, shp, status, &
                                  source, long_name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_Track_1D), intent(out)    ::  self
    integer, intent(in)               ::  shp(3)  ! (/nv,ntx,nty/)
    integer, intent(out)              ::  status
    
    real, intent(in), optional                ::  source
    character(len=*), intent(in), optional    ::  long_name
    character(len=*), intent(in), optional    ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_1D_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init base class:
    call self%Track_Init( shp(2), shp(3), status, long_name=long_name, units=units )
    IF_NOTOK_RETURN(status=1)
    
    ! store:
    self%nv   = shp(1)

    ! fill value:
    self%fill_value = - huge(1.0)
    ! storage:
    allocate( self%data(self%nv,self%ntx,self%nty), source=self%fill_value, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! replace if necessary:
    if ( present(source) ) self%data = source
    
    ! ok
    status = 0
    
  end subroutine Track_1D_Init


  ! ***


  subroutine Track_1D_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Track_1D), intent(inout)   ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Track_1D_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%data, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! empty:
    self%ntx = 0
    self%nty = 0
    self%nv  = 0

    ! ok
    status = 0
    
  end subroutine Track_1D_Done
  
  
  ! ***


  subroutine Track_1D_NcInit( self, ncid, varname, status )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
  
    ! --- in/out ---------------------------------
    
    class(T_Track_1D), intent(out)    ::  self
    integer, intent(in)               ::  ncid
    character(len=*), intent(in)      ::  varname
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_1D_NcInit'
    
    integer, parameter     ::  ndim = 3
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    integer             ::  nd
    integer             ::  dimids(ndim)
    integer             ::  shp(ndim)
    integer             ::  idim
    
    ! --- begin ----------------------------------

    ! get variable id:
    status = NF90_Inq_Varid( ncid, varname, varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, ndims=nd, dimids=dimids )
    IF_NF90_NOTOK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, ndim
      ! scalar?
      if ( idim <= ndim-nd ) then
        shp(idim) = 1
      else
        ! get size:
        status = NF90_Inquire_Dimension( ncid, dimids(idim-(ndim-nd)), len=shp(idim) )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
    end do

    ! init storage for all pixels:
    call self%Init( shp, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read attributes:
    call self%NcGetAttrs( ncid, varid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read:
    status = NF90_Get_Var( ncid, varid, self%data )
    IF_NF90_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Track_1D_NcInit


  ! ***


  subroutine Track_1D_NcDef( self, ncid, varname, dimids, status )
  
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    use NetCDF, only : NF90_FLOAT
  
    ! --- in/out ---------------------------------
    
    class(T_Track_1D), intent(inout)      ::  self
    integer, intent(in)                   ::  ncid
    character(len=*), intent(in)          ::  varname
    integer, intent(in)                   ::  dimids(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_1D_NcDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! define variable:
    status = NF90_Def_Var( ncid, trim(varname), NF90_FLOAT, dimids, self%varid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! add attributes:
    call self%NcPutAttrs( ncid, self%varid, status )
    IF_NOTOK_RETURN(status=1)
    ! define attribute:
    status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Track_1D_NcDef
  
  
  ! ***


  ! write global array from root, each pe has same copy
  
  subroutine Track_1D_NcPutGlb( self, ncid, status )
  
    use NetCDF, only : NF90_Put_Var
    use GO    , only : goc
  
    ! --- in/out ---------------------------------
    
    class(T_Track_1D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_1D_NcPutGlb'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( goc%root ) then
      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, self%data )
      IF_NF90_NOTOK_RETURN(status=1)
    end if ! root
    
    ! ok
    status = 0
    
  end subroutine Track_1D_NcPutGlb
  



end module Pixels
