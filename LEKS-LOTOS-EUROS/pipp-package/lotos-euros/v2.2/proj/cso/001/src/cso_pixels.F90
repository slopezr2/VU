!###############################################################################
!
! Pixels - storage for satellite data
!
!###############################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; csol=nf90_strerror(status); call csoErr; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!###############################################################################
  
module CSO_Pixels

  use CSO_Logging, only : csol, csoPr, csoErr
  use CSO_NcFile , only : T_NcAttrs
  use CSO_NcFile , only : T_NcDims
  use NetCDF     , only : NF90_StrError, NF90_NOERR

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private

  public  ::  T_PixelDatas
  
  public  ::  T_Track_0D
  public  ::  T_Track_1D
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Pixels'
  
  ! output real kind:
  integer, parameter            ::  wpr = 4


  ! --- types --------------------------------
  
  ! generic storage for all types
  type T_PixelData
    ! name:
    character(len=32)         ::  name
    ! type description: 0D_i,0D, 1D, 2D
    character(len=4)          ::  dtype
    ! global pixel list, same on all domains?
    ! this is used to store footprints for selection:
    logical                   ::  glb
    ! pixel dimension:
    integer                   ::  npix
    ! other dims:
    integer                   ::  ndim
    character(len=32)         ::  dimnames(7)
    integer                   ::  shp(7)
    ! storage:
    integer, pointer          ::  data0_i(:)    ! (npix)
    real, pointer             ::  data0(:)      ! (npix)
    real, pointer             ::  data1(:,:)    ! (nv,npix)
    real, pointer             ::  data2(:,:,:)  ! (mv,nv,npix)
    ! flag:
    logical                   ::  alloced
    ! source values:
    integer                   ::  source_i
    real(wpr)                 ::  source
    ! attributes:
    type(T_NcAttrs)           ::  attrs
    ! formula:
    character(:), allocatable ::  formula
    character(:), allocatable ::  formula_terms
    ! netcdf output:
    integer                   ::  varid
    integer                   ::  fill_value_i
    real(wpr)                 ::  fill_value
    !
  contains
    procedure :: Init            => PixelData_Init
    procedure :: Done            => PixelData_Done
    procedure :: Alloc           => PixelData_Alloc
    procedure :: NcInit          => PixelData_NcInit
    procedure :: Get             => PixelData_Get
    procedure :: Set             => PixelData_Set
    procedure :: SetPixelZero    => PixelData_SetPixelZero
    procedure ::                    PixelData_SetPixel_0D
    procedure ::                    PixelData_SetPixel_1D
    procedure ::                    PixelData_SetPixel_2D
    generic   :: SetPixel        => PixelData_SetPixel_0D, &
                                    PixelData_SetPixel_1D, &
                                    PixelData_SetPixel_2D
    procedure :: NcDef           => PixelData_NcDef
    procedure :: NcPutGlbSelect  => PixelData_NcPutGlbSelect
    procedure :: NcPutGather     => PixelData_NcPutGather
  end type T_PixelData
  
  ! *
  
  ! generic storage for all types
  type P_PixelData
    type(T_PixelData), pointer       ::  p
  end type P_PixelData
  
  ! *
  
  ! generic storage for multiple types:
  type T_PixelDatas
    ! storage:
    integer                         ::  n
    type(P_PixelData), pointer      ::  values(:)
!    ! list of unique dimensions:
!    integer                         ::  ndim
!    character(len=32), pointer      ::  dimnames(:)
!    integer, pointer                ::  dimlens(:)
    ! netcdf dimensions:
    type(T_NcDims)                  ::  ncdims
    ! flags:
    logical                         ::  initialized
    !
  contains
    procedure :: Init            => PixelDatas_Init
    procedure :: Done            => PixelDatas_Done
    procedure :: InqID           => PixelDatas_InqID
    procedure :: Get             => PixelDatas_Get
    procedure :: GetData         => PixelDatas_GetData
    procedure :: SetData         => PixelDatas_SetData
    procedure :: Def             => PixelDatas_Def
    procedure :: NcInit          => PixelDatas_NcInit
    procedure :: GetDim          => PixelDatas_GetDim
    procedure :: SetDim          => PixelDatas_SetDim
    procedure :: EndDef          => PixelDatas_EndDef
    procedure ::                    PixelDatas_SetAttr_i
    procedure ::                    PixelDatas_SetAttr_r
    procedure ::                    PixelDatas_SetAttr_c
    generic   :: SetAttr         => PixelDatas_SetAttr_i, &
                                    PixelDatas_SetAttr_r, &
                                    PixelDatas_SetAttr_c
    procedure ::                    PixelDatas_GetAttr_c
    generic   :: SetAttr         => PixelDatas_GetAttr_c
    procedure :: SetPixelZero    => PixelDatas_SetPixelZero
    procedure ::                    PixelDatas_SetPixel_0D
    procedure ::                    PixelDatas_SetPixel_1D
    procedure ::                    PixelDatas_SetPixel_2D
    generic   :: SetPixel        => PixelDatas_SetPixel_0D, &
                                    PixelDatas_SetPixel_1D, &
                                    PixelDatas_SetPixel_2D
    procedure :: SetFormula      => PixelDatas_SetFormula
    procedure :: GetFormulaData  => PixelDatas_GetFormulaData
    procedure :: ApplyFormulas   => PixelDatas_ApplyFormulas
    procedure :: NcDef           => PixelDatas_NcDef
    procedure :: NcPutGlbSelect  => PixelDatas_NcPutGlbSelect
    procedure :: NcPutGather     => PixelDatas_NcPutGather
  end type T_PixelDatas
  
  
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
    real, pointer         ::  data(:,:)  ! (ntx,nty)
    ! netcdf output:
    real(wpr)             ::  fill_value
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
    real, pointer         ::  data(:,:,:)  ! (nv,ntx,nty)
    ! netcdf output:
    real(wpr)             ::  fill_value
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
  ! === PixelData
  ! ===
  ! ====================================================================


  subroutine PixelData_Init( self, name, status, &
                                dtype, shp, glb, source, source_i )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(out)   ::  self
    character(len=*), intent(in)      ::  name
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  dtype
    integer, intent(in), optional             ::  shp(:)
    logical, intent(in), optional             ::  glb
    integer, intent(in), optional             ::  source_i
    real, intent(in), optional                ::  source

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/PixelData_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%name = trim(name)
    
    ! global arrays?
    self%glb = .false.
    if ( present(glb) ) self%glb = glb
    
    ! default data type is not defined:
    self%dtype = ''
    ! no dims:
    self%shp = -999
    
    ! fill values:
    self%fill_value_i = huge(1)
    self%fill_value   = huge(real(1.0,kind=wpr))
    
    ! no data allocated yet:
    self%alloced = .false.
    ! source values when allocated:    
    self%source_i = self%fill_value_i
    self%source   = self%fill_value
    ! reset?
    if ( present(source_i) ) self%source_i = source_i
    if ( present(source  ) ) self%source   = source
    
    ! init attribures:
    call self%attrs%Init( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine PixelData_Init


  ! ***


  subroutine PixelData_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelData_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(self%dtype) )
      !~
      case ( '0D_i' )
        deallocate( self%data0_i, stat=status )
        IF_NOT_OK_RETURN(status=1)
      !~
      case ( '0D' )
        deallocate( self%data0, stat=status )
        IF_NOT_OK_RETURN(status=1)
      !~
      case ( '1D' )
        deallocate( self%data1, stat=status )
        IF_NOT_OK_RETURN(status=1)
      !~
      case ( '2D' )
        deallocate( self%data2, stat=status )
        IF_NOT_OK_RETURN(status=1)
      !~
      case default
        write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
    end select
    ! empty:
    self%npix = 0
    
    ! done with attribures:
    call self%attrs%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine PixelData_Done
  
  
  ! ***


  subroutine PixelData_Alloc( self, shp, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(in)                 ::  shp(:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/PixelData_Alloc'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! number of dims (except pixel dimension):
    select case ( trim(self%dtype) )
      !~
      case ( '0D_i' )
        ! check ...
        if ( size(shp) /= 1 ) then
          write (csol,'("expected shape (npix), received size: ",i0)') size(shp); call csoErr
          TRACEBACK; status=1; return
        end if
        ! set dims::
        self%ndim = 0
        self%npix = shp(1)
        self%shp(1) = shp(1)
        ! storage:
        allocate( self%data0_i(max(1,self%npix)), source=self%source_i, stat=status )
        IF_NOT_OK_RETURN(status=1)
      !~
      case ( '0D' )
        ! check ...
        if ( size(shp) /= 1 ) then
          write (csol,'("expected shape (npix), received size: ",i0)') size(shp); call csoErr
          TRACEBACK; status=1; return
        end if
        ! set dims::
        self%ndim = 0
        self%npix = shp(1)
        self%shp(1) = shp(1)
        ! storage:
        allocate( self%data0(max(1,self%npix)), source=self%source, stat=status )
        IF_NOT_OK_RETURN(status=1)
      !~
      case ( '1D' )
        ! check ...
        if ( size(shp) /= 2 ) then
          write (csol,'("expected shape (npix), received size: ",i0)') size(shp); call csoErr
          TRACEBACK; status=1; return
        end if
        ! set dims::
        self%ndim = 1
        self%npix = shp(2)
        self%shp(1:2) = shp(1:2)
        ! storage:
        allocate( self%data1(shp(1),max(1,self%npix)), source=self%source, stat=status )
        IF_NOT_OK_RETURN(status=1)
      !~
      case ( '2D' )
        ! check ...
        if ( size(shp) /= 3 ) then
          write (csol,'("expected shape (npix), received size: ",i0)') size(shp); call csoErr
          TRACEBACK; status=1; return
        end if
        ! set dims:
        self%ndim = 2
        self%npix = shp(3)
        self%shp(1:3) = shp(1:3)
        ! storage:
        allocate( self%data2(shp(1),shp(2),max(1,self%npix)), source=self%source, stat=status )
        IF_NOT_OK_RETURN(status=1)
      !~
      case default
        write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
    end select
    
    ! reset flag:
    self%alloced = .true.

    ! ok
    status = 0
    
  end subroutine PixelData_Alloc
  
  
  ! ***


  subroutine PixelData_NcInit( self, ncid, varname, status, &
                                  nselect, select )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Inquire_Variable, NF90_Get_Var
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(in)                 ::  ncid
    character(len=*), intent(in)        ::  varname
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nselect    ! npix or 0
    integer, intent(in), optional       ::  select(:)  ! (npix) indices in global arrays

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/PixelData_NcInit'
    
    ! --- local ----------------------------------
    
    integer                 ::  varid
    integer                 ::  dimids(7)
    integer                 ::  xshp(7)
    integer                 ::  shp(7)
    integer                 ::  idim
    integer                 ::  ipix
    integer, allocatable    ::  data0_i(:)    ! (npix)
    real, allocatable       ::  data0(:)      ! (npix)
    real, allocatable       ::  data1(:,:)    ! (nv,npix)
    real, allocatable       ::  data2(:,:,:)  ! (mv,nv,npix)
    
    ! --- begin ----------------------------------
    
    ! number of dims (except pixel dimension):
    select case ( trim(self%dtype) )
      case ( '0D', '0D_i' ) ; self%ndim = 0
      case ( '1D'         ) ; self%ndim = 1
      case ( '2D'         ) ; self%ndim = 2
      case ( '3D'         ) ; self%ndim = 3
      case default
        write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
    end select

    ! get variable id:
    status = NF90_Inq_Varid( ncid, varname, varid )
    if (status/=NF90_NOERR) then
      csol=nf90_strerror(status); call csoErr
      write (csol,'("variable name: ",a)') trim(varname); call csoErr
      TRACEBACK; status=1; return
    end if

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, dimids=dimids(1:self%ndim+1) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, self%ndim+1
      ! get size:
      status = NF90_Inquire_Dimension( ncid, dimids(idim), len=xshp(idim) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! selection?
    if ( present(select) ) then

      ! number of selected pixels:
      if ( present(nselect) ) then
        ! check ...
        if ( (nselect < 0) .or. (nselect > size(select)) ) then
          write (csol,'("nselect=",i0," while size(select)=",i0)') nselect, size(select); call csoErr
          TRACEBACK; status=1; return
        end if
        ! copy:
        self%npix = nselect
      else
        ! full size:
        self%npix = size(select)
      end if
      
      ! target shape:
      shp = xshp
      shp(self%ndim+1) = self%npix
      ! storage:
      call self%Alloc( shp(1:self%ndim+1), status )
      IF_NOT_OK_RETURN(status=1)

      ! any local pixels?
      if ( self%npix > 0 ) then

        ! check ..
        if ( any(select < 1) .or. any(select > xshp(self%ndim+1)) ) then
          write (csol,'("selection indices in range ",i0,",..,",i0," while number of pixels is ",i0)') &
                   minval(select), maxval(select), xshp(self%ndim+1); call csoErr
          TRACEBACK; status=1; return
        end if

        ! switch:
        select case ( trim(self%dtype) )
          
          !~
          case ( '0D_i' )
            ! storage for all data:
            allocate( data0_i(xshp(1)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! read:
            status = NF90_Get_Var( ncid, varid, data0_i )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! loop over selected pixels:
            do ipix = 1, self%npix
              ! copy value from global array:
              self%data0_i(ipix) = data0_i(select(ipix))
            end do
            ! clear:
            deallocate( data0_i, stat=status )
            IF_NOT_OK_RETURN(status=1)
          
          !~
          case ( '0D' )
            ! storage for all data:
            allocate( data0(xshp(1)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! read:
            status = NF90_Get_Var( ncid, varid, data0 )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! loop over selected pixels:
            do ipix = 1, self%npix
              ! copy value from global array:
              self%data0(ipix) = data0(select(ipix))
            end do
            ! clear:
            deallocate( data0, stat=status )
            IF_NOT_OK_RETURN(status=1)
          
          !~
          case ( '1D' )
            ! storage for all data:
            allocate( data1(xshp(1),xshp(2)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! read:
            status = NF90_Get_Var( ncid, varid, data1 )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! loop over selected pixels:
            do ipix = 1, self%npix
              ! copy value from global array:
              self%data1(:,ipix) = data1(:,select(ipix))
            end do
            ! clear:
            deallocate( data1, stat=status )
            IF_NOT_OK_RETURN(status=1)
          
          !~
          case ( '2D' )
            ! storage for all data:
            allocate( data2(xshp(1),xshp(2),xshp(3)), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! read:
            status = NF90_Get_Var( ncid, varid, data2 )
            IF_NF90_NOT_OK_RETURN(status=1)
            ! loop over selected pixels:
            do ipix = 1, self%npix
              ! copy value from global array:
              self%data2(:,:,ipix) = data2(:,:,select(ipix))
            end do
            ! clear:
            deallocate( data2, stat=status )
            IF_NOT_OK_RETURN(status=1)

          !~
          case default
            write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
            TRACEBACK; status=1; return
        end select

      end if ! npix > 0
    
    else

      ! storage:
      call self%Alloc( xshp(1:self%ndim+1), status )
      IF_NOT_OK_RETURN(status=1)

      ! switch:
      select case ( trim(self%dtype) )
        !~
        case ( '0D_i' )
          ! read:
          status = NF90_Get_Var( ncid, varid, self%data0_i )
          IF_NF90_NOT_OK_RETURN(status=1)
        !~
        case ( '0D' )
          ! read:
          status = NF90_Get_Var( ncid, varid, self%data0 )
          IF_NF90_NOT_OK_RETURN(status=1)
        !~
        case ( '1D' )
          ! read:
          status = NF90_Get_Var( ncid, varid, self%data1 )
          IF_NF90_NOT_OK_RETURN(status=1)
        !~
        case ( '2D' )
          ! read:
          status = NF90_Get_Var( ncid, varid, self%data2 )
          IF_NF90_NOT_OK_RETURN(status=1)
        !~
        case default
          write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
          TRACEBACK; status=1; return
      end select
          
    end if ! select
    
    ! read attributes:
    call self%attrs%NcGet( ncid, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine PixelData_NcInit


  ! ***


  subroutine PixelData_Set( self, status, dtype )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(out)                ::  status
    
    character(len=*), intent(in), optional    ::  dtype
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelData_Set'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    if ( present(dtype) ) self%dtype = trim(dtype)

    ! ok
    status = 0
    
  end subroutine PixelData_Set


  ! ***


  subroutine PixelData_Get( self, status, &
                              data0, data1, data2, isfc_pressure )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(in)      ::  self
    integer, intent(out)                ::  status
    
    real, pointer, optional             ::  data0(:)       ! (npix)
    real, pointer, optional             ::  data1(:,:)     ! (nr,npix)
    real, pointer, optional             ::  data2(:,:,:)   ! (mr,nr,npix)
    integer, intent(out), optional      ::  isfc_pressure
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelData_Get'
    
    ! --- local ----------------------------------
    
    integer     ::  ipix
    
    ! --- begin ----------------------------------
    
    ! pointer to data?
    if ( present(data0) ) then  
      ! check ...
      if ( trim(self%dtype) /= '0D' ) then
        write (csol,'("argument `data0` requires dtype `0D`, not `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
      end if
      ! assign:
      data0 => self%data0
    end if ! data0
    
    ! pointer to data?
    if ( present(data1) ) then  
      ! check ...
      if ( trim(self%dtype) /= '1D' ) then
        write (csol,'("argument `data1` requires dtype `1D`, not `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
      end if
      ! assign:
      data1 => self%data1
    end if ! data1
    
    ! pointer to data?
    if ( present(data2) ) then  
      ! check ...
      if ( trim(self%dtype) /= '2D' ) then
        write (csol,'("argument `data2` requires dtype `2D`, not `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
      end if
      ! assign:
      data2 => self%data2
    end if ! data2
    
    ! level of surface pressure:
    if ( present(isfc_pressure) ) then 
      ! check ...
      if ( trim(self%dtype) /= '1D' ) then
        write (csol,'("argument `isfc_pressure` requires dtype `1D`, not `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
      end if
      ! define index of surface value;
      ! scan pixel until first filled pressure profile is found:
      do ipix = 1, self%npix
        ! no data? try next ...
        if ( self%data1(1,ipix) == self%fill_value ) cycle
        ! bottom-up or top-down ?
        if ( self%data1(1,ipix) > self%data1(self%shp(1),ipix) ) then
          isfc_pressure = 1
        else
          isfc_pressure = self%shp(1)
        end if
        ! found ...
        exit
      end do ! ipix
    end if ! isfc_pressure

    ! ok
    status = 0
    
  end subroutine PixelData_Get


  ! ***


  subroutine PixelData_SetPixelZero( self, ipix, status )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(in)                 ::  ipix
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelData_SetPixelZero'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( trim(self%dtype) )
      !~
      case ( '0D_i' )
        self%data0_i(ipix) = 0
      !~
      case ( '0D' )
        self%data0(ipix) = 0.0
      !~
      case ( '1D' )
        self%data1(:,ipix) = 0.0
      !~
      case ( '2D' )
        self%data2(:,:,ipix) = 0.0
      !~
      case default
        write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
  end subroutine PixelData_SetPixelZero


  ! ***


  subroutine PixelData_SetPixel_0D( self, ipix, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(in)                 ::  ipix
    real, intent(in)                    ::  value
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelData_SetPixel_0D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( trim(self%dtype) /= '0D' ) then
      write (csol,'("pixeldata type should be `0D`, found `",a,"`")') trim(self%dtype); call csoErr
      TRACEBACK; status=1; return
    end if
    
    !! testing ..
    !if ( trim(self%name) == 'mod_tcc' ) then
    !  if ( ipix == 27373 ) print *, 'sss1 ', ipix, value
    !end if

    ! fill:
    self%data0(ipix) = value

    ! ok
    status = 0
    
  end subroutine PixelData_SetPixel_0D


  ! ***


  subroutine PixelData_SetPixel_1D( self, ipix, values, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(in)                 ::  ipix
    real, intent(in)                    ::  values(:)
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelData_SetPixel_1D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( trim(self%dtype) /= '1D' ) then
      write (csol,'("pixeldata type should be `1D`, found `",a,"`")') trim(self%dtype); call csoErr
      TRACEBACK; status=1; return
    end if

    ! fill:
    self%data1(:,ipix) = values

    ! ok
    status = 0
    
  end subroutine PixelData_SetPixel_1D


  ! ***


  subroutine PixelData_SetPixel_2D( self, ipix, values, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(in)                 ::  ipix
    real, intent(in)                    ::  values(:,:)
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelData_SetPixel_2D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( trim(self%dtype) /= '2D' ) then
      write (csol,'("pixeldata type should be `2D`, found `",a,"`")') trim(self%dtype); call csoErr
      TRACEBACK; status=1; return
    end if

    ! fill:
    self%data2(:,:,ipix) = values

    ! ok
    status = 0
    
  end subroutine PixelData_SetPixel_2D
  

  ! ***
  
  
  subroutine PixelData_NcDef( self, ncid, ncdims, dimid_pixel, status )
  
    use NetCDF, only : NF90_INT
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att

    use CSO_NcFile, only : T_NcDims
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(inout)   ::  self
    integer, intent(in)                 ::  ncid
    type(T_NcDims), intent(in)          ::  ncdims
    integer, intent(in)                 ::  dimid_pixel
    integer, intent(out)                ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelData_NcDef'
    
    ! --- local ----------------------------------
    
    integer       ::  idim
    integer       ::  dimids(7)
    
    ! --- begin ----------------------------------
    
    ! loop over dimensions:
    do idim = 1, self%ndim
      ! get nc dimid:
      call ncdims%GetDim( trim(self%dimnames(idim)), status, dimid=dimids(idim) )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! pixel dimension:
    dimids(self%ndim+1) = dimid_pixel
    
    ! switch:
    select case ( trim(self%dtype) )

      !~
      case ( '0D_i' )
        ! define variable:
        status = NF90_Def_Var( ncid, trim(self%name), NF90_INT, dimids(1:self%ndim+1), self%varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! define sepecial attribute:
        status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value_i )
        IF_NF90_NOT_OK_RETURN(status=1)

      !~
      case ( '0D', '1D', '2D' )
        ! define variable:
        status = NF90_Def_Var( ncid, trim(self%name), NF90_FLOAT, dimids(1:self%ndim+1), self%varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! define sepecial attribute:
        status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
        IF_NF90_NOT_OK_RETURN(status=1)

      !~
      case default
        write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
    end select

    ! put attributes:
    call self%attrs%NcPut( ncid, self%varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine PixelData_NcDef
  
  
  ! ***


  ! write from root, this is a global array
  
  subroutine PixelData_NcPutGlbSelect( self, ncid, mapping, status )
  
    use NetCDF  , only : NF90_Put_Var
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(in)      ::  self
    integer, intent(in)                 ::  ncid
    integer, intent(in)                 ::  mapping(:)  ! (nglb) target index, or <0 if not used
    integer, intent(out)                ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/PixelData_NcPutGlbSelect'
    
    ! --- local ----------------------------------

    integer                 ::  nglb
    integer                 ::  n
    integer, allocatable    ::  data0_i(:)    ! (n)
    real, allocatable       ::  data0(:)      ! (n)
    real, allocatable       ::  data1(:,:)    ! (:,n)
    real, allocatable       ::  data2(:,:,:)  ! (:,:,n)
    integer                 ::  k
    integer                 ::  iglb

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( csoc%root ) then
    
      ! number of global pixels:
      nglb = size(mapping)
      ! number of selected pixels:
      n = maxval(mapping)

      ! switch:
      select case ( trim(self%dtype) )
      
        !~
        case ( '0D_i' )
        
          ! target array:
          allocate( data0_i(n), source=self%fill_value_i, stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! loop over global pixels:
          do iglb = 1, nglb
            ! target index:
            k = mapping(iglb)
            ! defined?
            if ( k > 0 ) then
              ! copy:
              data0_i(k) = self%data0_i(iglb)
            end if
          end do ! glb

          ! write variable:
          status = NF90_Put_Var( ncid, self%varid, data0_i )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data0_i, stat=status )
          IF_NOT_OK_RETURN(status=1)
      
        !~
        case ( '0D' )
        
          ! target array:
          allocate( data0(n), source=self%fill_value, stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! loop over global pixels:
          do iglb = 1, nglb
            ! target index:
            k = mapping(iglb)
            ! defined?
            if ( k > 0 ) then
              ! copy:
              data0(k) = self%data0(iglb)
            end if
          end do ! glb

          ! write variable:
          status = NF90_Put_Var( ncid, self%varid, data0 )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data0, stat=status )
          IF_NOT_OK_RETURN(status=1)
      
        !~
        case ( '1D' )
        
          ! target array:
          allocate( data1(self%shp(1),n), source=self%fill_value, stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! loop over global pixels:
          do iglb = 1, nglb
            ! target index:
            k = mapping(iglb)
            ! defined?
            if ( k > 0 ) then
              ! copy:
              data1(:,k) = self%data1(:,iglb)
            end if
          end do ! glb

          ! write variable:
          status = NF90_Put_Var( ncid, self%varid, data1 )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data1, stat=status )
          IF_NOT_OK_RETURN(status=1)
      
        !~
        case ( '2D' )
        
          ! target array:
          allocate( data2(self%shp(1),self%shp(2),n), source=self%fill_value, stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! loop over global pixels:
          do iglb = 1, nglb
            ! target index:
            k = mapping(iglb)
            ! defined?
            if ( k > 0 ) then
              ! copy:
              data2(:,:,k) = self%data2(:,:,iglb)
            end if
          end do ! glb

          ! write variable:
          status = NF90_Put_Var( ncid, self%varid, data2 )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data2, stat=status )
          IF_NOT_OK_RETURN(status=1)

        !~
        case default
          write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
          TRACEBACK; status=1; return
      end select

    end if ! root

    ! ok
    status = 0
    
  end subroutine  PixelData_NcPutGlbSelect
  
  
  ! ***


  !
  ! gather all pixels on root,
  ! add together into target array using "mapping" indices,
  ! and write variable
  !
  
  subroutine PixelData_NcPutGather( self, ncid, mapping, add, status )
  
    use NetCDF  , only : NF90_Put_Var
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_PixelData), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  mapping(:)   ! (npix_all) target index
    logical, intent(in)               ::  add          ! add contributions?
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/PixelData_NcPutGather'
    
    ! --- local ----------------------------------

    integer                 ::  npix_all
    integer, allocatable    ::  data0_i_all(:)    ! (npix_all)
    integer, allocatable    ::  data0_i(:)        ! (nout)
    real, allocatable       ::  data0_all(:)      ! (npix_all)
    real, allocatable       ::  data0(:)          ! (nout)
    real, allocatable       ::  data1_all(:,:)    ! (:,npix_all)
    real, allocatable       ::  data1(:,:)        ! (:,nout)
    real, allocatable       ::  data2_all(:,:,:)  ! (:,:,npix_all)
    real, allocatable       ::  data2(:,:,:)      ! (:,:,nout)
    integer                 ::  nout
    integer                 ::  iout
    integer                 ::  iall
    
    !logical                 ::  debug

    ! --- begin ----------------------------------
    
    !! testing ...
    !write (csol,'(a,": write `",a,"`")') rname, trim(self%name); call csoPr

    ! switch:
    select case ( trim(self%dtype) )

      !~
      case ( '0D_i' )

        ! gather on root only ..
        if ( csoc%root ) then
          ! count:
          npix_all = size(mapping)
          ! storage for all data:
          allocate( data0_i_all(npix_all), stat=status )
          IF_NOT_OK_RETURN(status=1)
        else
          ! dummy ...
          allocate( data0_i_all(1), stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if ! root

        ! gather data arrays on root, local size might be zero
        call csoc%GatherV( self%data0_i, data0_i_all, status, nloc=self%npix )
        IF_NOT_OK_RETURN(status=1)

        ! gather on root only ..
        if ( csoc%root ) then

          ! number of selected pixels:
          nout = maxval(mapping)
          ! target array:
          allocate( data0_i(nout), source=self%fill_value_i, stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop over collected pixels:
          do iall = 1, npix_all
            ! skip undefined source:
            if ( data0_i_all(iall) >= self%fill_value_i ) cycle
            ! target index:
            iout = mapping(iall)
            ! first contribution?
            if ( data0_i(iout) == self%fill_value_i ) then
              ! copy:
              data0_i(iout) = data0_i_all(iall)
            else if ( add ) then
              ! add contribution:
              data0_i(iout) = data0_i(iout) + data0_i_all(iall)
            end if
          end do ! iall

          ! write variable:
          status = NF90_Put_Var( ncid, self%varid, data0_i )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data0_i, stat=status )
          IF_NOT_OK_RETURN(status=1)

        end if ! root

        ! clear:
        deallocate( data0_i_all, stat=status )
        IF_NOT_OK_RETURN(status=1)

      !~
      case ( '0D' )

        ! gather on root only ..
        if ( csoc%root ) then
          ! count:
          npix_all = size(mapping)
          ! storage for all data:
          allocate( data0_all(npix_all), stat=status )
          IF_NOT_OK_RETURN(status=1)
        else
          ! dummy ...
          allocate( data0_all(1), stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if ! root

        ! gather data arrays on root, local size might be zero
        call csoc%GatherV( self%data0, data0_all, status, nloc=self%npix )
        IF_NOT_OK_RETURN(status=1)

        ! gather on root only ..
        if ( csoc%root ) then
        
          ! number of selected pixels:
          nout = maxval(mapping)
          ! target array:
          allocate( data0(nout), source=self%fill_value, stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop over collected pixels:
          do iall = 1, npix_all
            ! skip undefined source:
            if ( data0_all(iall) >= self%fill_value ) cycle
            ! target index:
            iout = mapping(iall)
            ! first contribution?
            if ( data0(iout) == self%fill_value ) then
              ! copy:
              data0(iout) = data0_all(iall)
            else if ( add ) then
              ! add contribution:
              data0(iout) = data0(iout) + data0_all(iall)
            end if
          end do ! iall

          ! write variable:
          status = NF90_Put_Var( ncid, self%varid, data0 )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data0, stat=status )
          IF_NOT_OK_RETURN(status=1)

        end if ! root

        ! clear:
        deallocate( data0_all, stat=status )
        IF_NOT_OK_RETURN(status=1)

      !~
      case ( '1D' )

        ! gather on root only ..
        if ( csoc%root ) then
          ! count:
          npix_all = size(mapping)
          ! storage for all data:
          allocate( data1_all(self%shp(1),npix_all), stat=status )
          IF_NOT_OK_RETURN(status=1)
        else
          ! dummy ...
          allocate( data1_all(self%shp(1),1), stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if ! root

        ! gather data arrays on root, local size might be zero
        call csoc%GatherV( self%data1, data1_all, status, nloc=self%npix )
        IF_NOT_OK_RETURN(status=1)

        ! gather on root only ..
        if ( csoc%root ) then

          ! number of selected pixels:
          nout = maxval(mapping)
          ! target array:
          allocate( data1(self%shp(1),nout), source=self%fill_value, stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop over collected pixels:
          do iall = 1, npix_all
            ! skip undefined source:
            if ( data1_all(1,iall) >= self%fill_value ) cycle
            ! target index:
            iout = mapping(iall)
            ! first contribution?
            if ( data1(1,iout) == self%fill_value ) then
              ! copy:
              data1(:,iout) = data1_all(:,iall)
            else if ( add ) then
              ! add contribution:
              data1(:,iout) = data1(:,iout) + data1_all(:,iall)
            end if
          end do ! iall

          ! write variable:
          status = NF90_Put_Var( ncid, self%varid, data1 )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data1, stat=status )
          IF_NOT_OK_RETURN(status=1)

        end if ! root

        ! clear:
        deallocate( data1_all, stat=status )
        IF_NOT_OK_RETURN(status=1)

      !~
      case ( '2D' )

        ! gather on root only ..
        if ( csoc%root ) then
          ! count:
          npix_all = size(mapping)
          ! storage for all data:
          allocate( data2_all(self%shp(1),self%shp(2),npix_all), stat=status )
          IF_NOT_OK_RETURN(status=1)
        else
          ! dummy ...
          allocate( data2_all(self%shp(1),self%shp(2),1), stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if ! root

        ! gather data arrays on root, local size might be zero
        call csoc%GatherV( self%data2, data2_all, status, nloc=self%npix )
        IF_NOT_OK_RETURN(status=1)

        ! gather on root only ..
        if ( csoc%root ) then

          ! number of selected pixels:
          nout = maxval(mapping)
          ! target array:
          allocate( data2(self%shp(1),self%shp(2),nout), source=self%fill_value, stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop over collected pixels:
          do iall = 1, npix_all
            ! skip undefined source:
            if ( data2_all(1,1,iall) >= self%fill_value ) cycle
            ! target index:
            iout = mapping(iall)
            ! first contribution?
            if ( data2(1,1,iout) == self%fill_value ) then
              ! copy:
              data2(:,:,iout) = data2_all(:,:,iall)
            else if ( add ) then
              ! add contribution:
              data2(:,:,iout) = data2(:,:,iout) + data2_all(:,:,iall)
            end if
          end do ! iall

          ! write variable:
          status = NF90_Put_Var( ncid, self%varid, data2 )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! clear:
          deallocate( data2, stat=status )
          IF_NOT_OK_RETURN(status=1)

        end if ! root

        ! clear:
        deallocate( data2_all, stat=status )
        IF_NOT_OK_RETURN(status=1)

      !~
      case default
        write (csol,'("unsupported data type `",a,"`")') trim(self%dtype); call csoErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
  end subroutine  PixelData_NcPutGather


  ! ====================================================================
  ! ===
  ! === PixelDatas
  ! ===
  ! ====================================================================


  subroutine PixelDatas_Init( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(out)   ::  self
    integer, intent(out)               ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/PixelDatas_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! no data yet:
    self%n = 0
    nullify( self%values )
    ! not initialized yet:
    self%initialized = .false.
    
!    ! no dimensions yet:
!    self%ndim = 0
!    nullify( self%dimnames )
!    nullify( self%dimlens )
    
    ! init netcdf dimension list:
    call self%ncdims%Init( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine PixelDatas_Init


  ! ***


  subroutine PixelDatas_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)   ::  self
    integer, intent(out)                 ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_Done'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! initialized?
    if ( self%initialized ) then
      ! defined?
      if ( self%n > 0 ) then
        ! loop:
        do i = 1, self%n
          ! done with attribute:
          call self%values(i)%p%Done( status )
          IF_NOT_OK_RETURN(status=1)
        end do ! i
      end if  ! n > 0
    end if ! initialized
    
!    ! dimensions defined?
!    if ( self%ndim > 0 ) then
!      ! clear:
!      deallocate( self%dimnames, stat=status )
!      IF_NOT_OK_RETURN(status=1)
!      deallocate( self%dimlens, stat=status )
!      IF_NOT_OK_RETURN(status=1)
!      ! no values anymore:
!      self%ndim = 0
!    end if

    ! done with netcdf dimension list:
    call self%ncdims%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_Done


  ! ***


  !
  ! Return index of pixel data set.
  ! If then name is not found, by default an error is raised and status>0 is returned,
  ! but with quiet=.true. only a warning status (<0) is returned.
  !

  subroutine PixelDatas_InqID( self, name, id, status, quiet )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)      ::  self
    character(len=*), intent(in)         ::  name
    integer, intent(out)                 ::  id
    integer, intent(out)                 ::  status
    
    logical, intent(in), optional        ::  quiet
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_InqID'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    logical     ::  verbose
    
    ! --- begin ----------------------------------
    
    ! shout?
    verbose = .true.
    if ( present(quiet) ) verbose = .not. quiet
    
    ! search:
    id = -999
    do i = 1, self%n
      if ( trim(self%values(i)%p%name) == trim(name) ) then
        id = i
        exit
      end if
    end do
    ! check ...
    if ( id < 0 ) then
      if ( verbose ) then
        write (csol,'("name `",a,"` not found in:")') trim(name); call csoErr
        do i = 1, self%n
          write (csol,'("  ",a)') trim(self%values(i)%p%name); call csoErr
        end do
        TRACEBACK; status=1; return
      else
        ! warning ...
        status = -1; return
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_InqID


  ! ***


  subroutine PixelDatas_Append_Empty( self, status )
    
    use CSO_String, only : CSO_SplitString, CSO_MatchValue
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_Append_Empty'
    
    ! --- local ----------------------------------
    
    type(P_PixelData), pointer      ::  new_values(:)
    integer                         ::  i
    
    ! --- begin ----------------------------------
    
    ! new storage:
    allocate( new_values(self%n+1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! allready values present?
    if ( self%n > 0 ) then
      ! loop over current values:
      do i = 1, self%n
        ! assign value:
        new_values(i)%p => self%values(i)%p
      end do
      ! clear current storage:
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! n > 0
    ! assign value list:
    self%values => new_values    
    
    ! increase counter:
    self%n = self%n + 1
    ! new value:
    allocate( self%values(self%n)%p, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_Append_Empty


  ! ***
  
  !
  ! Define new pixeldata:
  ! - vname: name for the dataset
  ! - dname: dimensions (apart from pixel dimension): 'layer' ;
  !   length read from netcdf file if ncid is present
  ! Output:
  ! - ivar: dataset index
  ! - status: 0, or >0 for error
  !
  ! Optional:
  ! - ncid : netcdf file id, used to read dimension lengths
  !

  subroutine PixelDatas_Def( self, vname, dnames, ivar, status, &
                               glb, xtype, source_i, source, ncid )
                                  
    use NetCDF, only : NF90_Inq_DimID
    use NetCDF, only : NF90_Inquire_Dimension
    
    use CSO_String, only : CSO_SplitString, CSO_MatchValue
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    character(len=*), intent(in)        ::  vname
    character(len=*), intent(in)        ::  dnames
    integer, intent(out)                ::  ivar
    integer, intent(out)                ::  status
    
    logical, intent(in), optional             ::  glb
    character(len=*), intent(in), optional    ::  xtype
    integer, intent(in), optional             ::  source_i
    real, intent(in), optional                ::  source
    integer, intent(in), optional             ::  ncid
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_Def'
    
    ! --- local ----------------------------------
    
    type(T_PixelData), pointer      ::  p
    character(len=8)                ::  xtype_curr
    character(len=8)                ::  dtype
    integer                         ::  idim
    integer                         ::  dimid
!    character(len=32), pointer      ::  new_dimnames(:)
!    integer, pointer                ::  new_dimlens(:)
!    integer                        ::  ndim
!    character(len=32)              ::  dimnames(7)
    integer                        ::  length
    
    ! --- begin ----------------------------------
    
    ! increase storage:
    call PixelDatas_Append_Empty( self, status )
    IF_NOT_OK_RETURN(status=1)
    ! variable id:
    ivar = self%n

    ! short:
    p => self%values(ivar)%p

    ! initialize with name, no type and allocation yet:
    call p%Init( vname, status, glb=glb, source_i=source_i, source=source )
    IF_NOT_OK_RETURN(status=1)
    
    ! split list with dimension names:
    call CSO_SplitString( dnames, p%ndim, p%dimnames, status, sep=' ' )
    IF_NOT_OK_RETURN(status=1)
!    call CSO_SplitString( dnames, ndim, dimnames, status, sep=' ' )
!    IF_NOT_OK_RETURN(status=1)
    !! info ...
    !write (csol,'(a,":  number of dimensions: ",i0)') rname, p%ndim; call csoPr
 
    ! int/real?
    xtype_curr = 'real'
    if ( present(xtype) ) xtype_curr = trim(xtype)

    ! type: 0D_i, 0D, 1D, ..
    select case ( trim(xtype_curr) )
      case ( 'int' )
        write (dtype,'(i1,"D_i")') p%ndim
!        write (dtype,'(i1,"D_i")') ndim
      case ( 'real' )
        write (dtype,'(i1,"D")') p%ndim
!        write (dtype,'(i1,"D")') ndim
      case default
        write (csol,'("unsupported xtype `",a,"`")') trim(xtype_curr); call csoErr
        TRACEBACK; status=1; return
    end select
    ! store:
    call p%Set( status, dtype=dtype )
    IF_NOT_OK_RETURN(status=1)
    
!    ! loop over dims:
!    do idim = 1, p%ndim
!      ! already some dims defined?
!      if ( self%ndim > 0 ) then
!        ! search, status<0 if not found, but do not shout error messages in that case:
!        call CSO_MatchValue( trim(p%dimnames(idim)), self%dimnames(1:self%ndim), dimid, status, quiet=.true. )
!        IF_ERROR_RETURN(status=1)
!        if ( status < 0 ) dimid = -1
!      else
!        ! not present yet:
!        dimid = -1
!      end if
!      ! not found?
!      if ( dimid < 0 ) then
!        !! info ..
!        !write (csol,'(a,":    new dim `",a,"`")') rname, trim(p%dimnames(idim)); call csoPr
!        ! new storage:
!        allocate( new_dimnames(self%ndim+1), stat=status )
!        IF_NOT_OK_RETURN(status=1)
!        allocate( new_dimlens(self%ndim+1), stat=status )
!        IF_NOT_OK_RETURN(status=1)
!        ! allready values present?
!        if ( self%ndim > 0 ) then
!          ! copy values:
!          new_dimnames(1:self%ndim) = self%dimnames(1:self%ndim)
!          new_dimlens (1:self%ndim) = self%dimlens (1:self%ndim)
!          ! clear current storage:
!          deallocate( self%dimnames, stat=status )
!          IF_NOT_OK_RETURN(status=1)
!          deallocate( self%dimlens , stat=status )
!          IF_NOT_OK_RETURN(status=1)
!        end if ! n > 0
!        ! assign values:
!        self%dimnames => new_dimnames
!        self%dimlens  => new_dimlens
!        ! increase counter:
!        self%ndim = self%ndim + 1
!        ! store name, no length yet:
!        self%dimnames(self%ndim) = p%dimnames(idim)
!        self%dimlens (self%ndim) = -999
!      end if ! idim < 0
!    end do ! dims
    ! loop over dims:
    do idim = 1, p%ndim
      ! read length?
      if ( present(ncid) ) then
        ! dim id:
        status = NF90_Inq_DimID( ncid, trim(p%dimnames(idim)), dimid )
        if ( status /= NF90_NOERR ) then
          write (csol,'("dimension `",a,"` not found in file")') trim(p%dimnames(idim)); call csoErr
          TRACEBACK; status=1; return
        end if
        ! length:
        status = NF90_Inquire_Dimension( ncid, dimid, len=length )
        IF_NF90_NOT_OK_RETURN(status=1)
      else
        ! dummy:
        length = -999
      end if
      ! append dimension, same name is stored once only,
      ! an error is raised if a name is used twice with different lengths:
      call self%ncdims%Append( p%dimnames(idim), length, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_Def


  ! ***
  
  
  subroutine PixelDatas_NcInit( self, name, dnames, ncid, varname, status, &
                                  nselect, select, glb, xtype )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  dnames
    integer, intent(in)                 ::  ncid
    character(len=*), intent(in)        ::  varname
    integer, intent(out)                ::  status
    
    integer, intent(in), optional             ::  nselect    ! npix or 0
    integer, intent(in), optional             ::  select(:)  ! (npix) indices in global arrays
    logical, intent(in), optional             ::  glb
    character(len=*), intent(in), optional    ::  xtype

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/PixelDatas_NcInit'
    
    ! --- local ----------------------------------
    
    integer     ::  ivar
!    integer     ::  idim
!    integer     ::  dimid

    ! --- begin ----------------------------------
    
    ! info ...
    write (csol,'(a,": init `",a,"` from nc variable ..")') rname, trim(varname); call csoPr
    
    ! define new pixeldata, read dimension lengths from netcdf file:
    call self%Def( name, dnames, ivar, status, glb=glb, xtype=xtype, ncid=ncid )
    IF_NOT_OK_RETURN(status=1)
    
!    ! loop:
!    do idim = 1, self%ndim
!      ! length not filed yet?
!      if ( self%dimlens(idim) < 0 ) then
!        ! dim id:
!        status = NF90_Inq_DimID( ncid, trim(self%dimnames(idim)), dimid )
!        if ( status /= NF90_NOERR ) then
!          write (csol,'("dimension `",a,"` not found in file")') trim(self%dimnames(idim)); call csoErr
!          TRACEBACK; status=1; return
!        end if
!        ! length:
!        status = NF90_Inquire_Dimension( ncid, dimid, len=self%dimlens(idim) )
!        IF_NF90_NOT_OK_RETURN(status=1)
!      end if ! inq dim
!    end do ! dims

    ! fill from ncfile:
    call self%values(self%n)%p%NcInit( ncid, varname, status, &
                                         nselect=nselect, select=select )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_NcInit


  ! ***


  subroutine PixelDatas_SetAttr_i( self, i, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  i
    character(len=*), intent(in)        ::  name
    integer, intent(in)                 ::  value
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetAttr_i'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check .. 
    if ( (i < 1) .or. (i > self%n) ) then
      write (csol,'("pixeldata index (",i0,") should be in range 1,..,",i0)') i, self%n; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! add attribute:
    call self%values(i)%p%attrs%Append( name, value, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine PixelDatas_SetAttr_i


  ! *


  subroutine PixelDatas_SetAttr_r( self, i, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  i
    character(len=*), intent(in)        ::  name
    real, intent(in)                    ::  value
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetAttr_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check .. 
    if ( (i < 1) .or. (i > self%n) ) then
      write (csol,'("pixeldata index (",i0,") should be in range 1,..,",i0)') i, self%n; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! add attribute:
    call self%values(i)%p%attrs%Append( name, value, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine PixelDatas_SetAttr_r
  
  ! *


  subroutine PixelDatas_SetAttr_c( self, id, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  id
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  value
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetAttr_c'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check .. 
    if ( (id < 1) .or. (id > self%n) ) then
      write (csol,'("pixeldata id (",i0,") should be in range 1,..,",i0)') id, self%n; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! add attribute:
    call self%values(id)%p%attrs%Append( name, value, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine PixelDatas_SetAttr_c


  ! ***


  subroutine PixelDatas_GetAttr_c( self, name, aname, value, status )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)      ::  self
    character(len=*), intent(in)         ::  name
    character(len=*), intent(in)         ::  aname
    character(len=*), intent(out)        ::  value
    integer, intent(out)                 ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_GetAttr_c'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! get id:
    call self%InqID( name, i, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! get from attributes:
    call self%values(i)%p%attrs%GetValue( aname, value, status )
    if ( status /= 0 ) then
      write (csol,'("variable `",a,"` has no `",a,"` attribute")') trim(name), trim(aname); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_GetAttr_c


  ! ***


  ! Return:
  !  - number of datasets
  !  - number of dims
  !  - number of user defined variables
  !  - user defined variable names
  !  - user defined variable units
    
  subroutine PixelDatas_Get( self, status, ndim, n, nuvar, uvarnames, uvarunits )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)      ::  self
    integer, intent(out)                 ::  status
    
    integer, intent(out), optional              ::  ndim
    integer, intent(out), optional              ::  n
    integer, intent(out), optional              ::  nuvar
    character(len=*), intent(out), optional     ::  uvarnames(:)
    character(len=*), intent(out), optional     ::  uvarunits(:)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_Get'
    
    ! --- local ----------------------------------
    
    integer                         ::  i
    type(T_PixelData), pointer      ::  p
    integer                         ::  iu
    
    ! --- begin ----------------------------------
    
    ! return number of user defined dimensions?
!    if ( present(ndim) ) ndim = self%ndim
    if ( present(ndim) ) then
      call self%ncdims%Get( status, n=ndim )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! number of values:
    if ( present(n) ) n = self%n
    
    ! number of user defined values:
    if ( present(nuvar) ) then
      ! init:
      nuvar = 0
      ! loop over values:
      do i = 1, self%n
        ! increase counter if not defined by formula:
        if ( len_trim(self%values(i)%p%formula) == 0 ) nuvar = nuvar + 1
      end do ! i
    end if ! nuvar
    
    ! names:
    if ( present(uvarnames) ) then
      ! init counter:
      iu = 0
      ! loop over values:
      do i = 1, self%n
        ! current:
        p => self%values(i)%p
        ! user defined?
        if ( len_trim(p%formula) == 0 ) then
          ! increase counter:
          iu = iu + 1
          ! check ...
          if ( len(uvarnames) < iu ) then
            write (csol,'("size of `uvarnames` (",i0,") should be at least ",i0)') size(uvarnames), iu; call csoErr
            TRACEBACK; status=1; return
          end if
          ! copy:
          uvarnames(iu) = trim(p%name)
        end if ! no formula
      end do ! i
    end if ! uvarnames
    
    ! units:
    if ( present(uvarunits) ) then
      ! init counter:
      iu = 0
      ! loop over values:
      do i = 1, self%n
        ! current:
        p => self%values(i)%p
        ! user defined?
        if ( len_trim(p%formula) == 0 ) then
          ! increase counter:
          iu = iu + 1
          ! check ...
          if ( len(uvarunits) < iu ) then
            write (csol,'("size of `uvarunits` (",i0,") should be at least ",i0)') size(uvarunits), iu; call csoErr
            TRACEBACK; status=1; return
          end if
          ! get from attributes:
          call p%attrs%GetValue( 'units', uvarunits(iu), status )
          if ( status /= 0 ) then
            write (csol,'("variable `",a,"` has no `units` attribute")') trim(p%name); call csoErr
            TRACEBACK; status=1; return
          end if
        end if ! no formula
      end do ! i
    end if ! uvarunits

    ! ok
    status = 0
    
  end subroutine PixelDatas_Get


  ! ***


  subroutine PixelDatas_GetData( self, status, &
                                   id, name, &
                                   data0, data1, data2, &
                                   units, isfc_pressure )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)      ::  self
    integer, intent(out)                 ::  status
    
    integer, intent(in), optional             ::  id
    character(len=*), intent(in), optional    ::  name
    real, pointer, optional                   ::  data0(:)
    real, pointer, optional                   ::  data1(:,:)
    real, pointer, optional                   ::  data2(:,:,:)
    character(len=*), intent(out), optional   ::  units
    integer, intent(out), optional            ::  isfc_pressure
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_GetData'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! get id:
    if ( present(id) ) then    
      ! check .. 
      if ( (id < 1) .or. (id > self%n) ) then
        write (csol,'("pixeldata id (",i0,") should be in range 1,..,",i0)') id, self%n; call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      i = id
    else if ( present(name) ) then
      ! from name:
      call self%InqID( name, i, status )
      IF_NOT_OK_RETURN(status=1)
    else
      write (csol,'("no `id` or `name` arguments provided")'); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! get:
    call self%values(i)%p%Get( status, data0=data0, data1=data1, data2=data2, isfc_pressure=isfc_pressure )
    IF_NOT_OK_RETURN(status=1)
    
    ! units:
    if ( present(units) ) then
      ! get from attributes:
      call self%values(i)%p%attrs%GetValue( 'units', units, status )
      if ( status /= 0 ) then
        write (csol,'("variable `",a,"` has no `units` attribute")') trim(self%values(i)%p%name); call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_GetData


  ! ***


  subroutine PixelDatas_SetData( self, status, &
                                   id, name, &
                                   data0, data1, data2 )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)      ::  self
    integer, intent(out)                 ::  status
    
    integer, intent(in), optional             ::  id
    character(len=*), intent(in), optional    ::  name
    real, intent(in), optional                ::  data0(:)
    real, intent(in), optional                ::  data1(:,:)
    real, intent(in), optional                ::  data2(:,:,:)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetData'
    
    ! --- local ----------------------------------
    
    real, pointer     ::  pdata0(:)
    real, pointer     ::  pdata1(:,:)
    real, pointer     ::  pdata2(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! data0:
    if ( present(data0) ) then
      ! get pointer:
      call self%GetData( status, id=id, name=name, data0=pdata0 )
      IF_NOT_OK_RETURN(status=1)
      ! check ..
      if ( any( shape(data0) /= shape(pdata0) ) ) then
        write (csol,'("input data has shape (",i0,") while storage has shape (",i0,")")') &
                            shape(data0), shape(pdata0); call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy values:
      pdata0 = data0
    end if ! data0
    
    ! data1:
    if ( present(data1) ) then
      ! get pointer:
      call self%GetData( status, id=id, name=name, data1=pdata1 )
      IF_NOT_OK_RETURN(status=1)
      ! check ..
      if ( any( shape(data1) /= shape(pdata1) ) ) then
        write (csol,'("input data has shape (",i0,",",i0,") while storage has shape (",i0,",",i0,")")') &
                            shape(data1), shape(pdata1); call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy values:
      pdata1 = data1
    end if ! data1
    
    ! data0:
    if ( present(data2) ) then
      ! get pointer:
      call self%GetData( status, id=id, name=name, data2=pdata2 )
      IF_NOT_OK_RETURN(status=1)
      ! check ..
      if ( any( shape(data2) /= shape(pdata2) ) ) then
        write (csol,'("input data has shape (",i0,",",i0,",",i0,") while storage has shape (",i0,",",i0,",",i0,")")') &
                            shape(data2), shape(pdata2); call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy values:
      pdata2 = data2
    end if ! data2
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_SetData


  ! ***


  subroutine PixelDatas_GetDim( self, id, status, name )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)               ::  self
    integer, intent(in)                           ::  id
    integer, intent(out)                          ::  status
    
    character(len=*), intent(out), optional       ::  name
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_GetDim'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
!    ! check ...
!    if ( (id < 1) .or. (id > self%ndim) ) then
!      write (csol,'("dimension ",i0," out of range 1:",i0)') id, self%ndim; call csoErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! dimension name:
!    if ( present(name) ) name = trim(self%dimnames(id))

    ! get dimension properties:
    call self%ncdims%GetDim( id, status, name=name )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_GetDim


  ! ***


  subroutine PixelDatas_SetDim( self, name, length, status )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)            ::  self
    character(len=*), intent(in)                  ::  name
    integer, intent(in)                           ::  length
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetDim'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! set dimension length:
    call self%ncdims%SetDim( name, length, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_SetDim


  ! ***


  subroutine PixelDatas_EndDef( self, npix, status )
  
    use CSO_String, only : CSO_MatchValue
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  npix
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_EndDef'
    
    ! --- local ----------------------------------
    
    type(T_PixelData), pointer      ::  p
    integer                         ::  id
    integer                         ::  idim
    integer                         ::  mv,nv
    
    ! --- begin ----------------------------------
    
    ! loop over data sets:
    do id = 1, self%n
      ! short:
      p => self%values(id)%p
      ! not needed if already allocated:
      if ( p%alloced ) cycle
      ! switch:
      select case ( trim(p%dtype) )
        !~
        case ( '0D_i', '0D' )
          ! storage:
          call p%Alloc( (/npix/), status )
          IF_NOT_OK_RETURN(status=1)
        !~
        case ( '1D' )
          ! dimension length:
          call self%ncdims%GetDim( p%dimnames(1), status, length=nv )
          IF_NOT_OK_RETURN(status=1)
          ! check ..
          if ( nv <= 0 ) then
            write (csol,'("dimension `",a,"` has length ",i0,"; not defined yet?")') trim(p%dimnames(1)), nv; call csoErr
            TRACEBACK; status=1; return
          end if
          ! init:
          call p%Alloc( (/nv,npix/), status )
          IF_NOT_OK_RETURN(status=1)
        !~
        case ( '2D' )
          ! dimension length:
          call self%ncdims%GetDim( p%dimnames(1), status, length=mv )
          IF_NOT_OK_RETURN(status=1)
          ! check ..
          if ( mv <= 0 ) then
            write (csol,'("dimension `",a,"` has length ",i0,"; not defined yet?")') trim(p%dimnames(1)), mv; call csoErr
            TRACEBACK; status=1; return
          end if
          ! dimension length:
          call self%ncdims%GetDim( p%dimnames(2), status, length=mv )
          IF_NOT_OK_RETURN(status=1)
          ! check ..
          if ( nv <= 0 ) then
            write (csol,'("dimension `",a,"` has length ",i0,"; not defined yet?")') trim(p%dimnames(2)), nv; call csoErr
            TRACEBACK; status=1; return
          end if
          ! init:
          call p%Alloc( (/mv,nv,npix/), status )
          IF_NOT_OK_RETURN(status=1)
        !~
        case default
          write (csol,'("unsupported `",a,"`")') trim(p%dtype); call csoErr
          TRACEBACK; status=1; return
      end select
    end do ! id

    ! now initialized:
    self%initialized = .true.

    ! ok
    status = 0
    
  end subroutine PixelDatas_EndDef


  ! ***


  subroutine PixelDatas_SetPixelZero( self, ipix, status )
    
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  ipix
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetPixelZero'
    
    ! --- local ----------------------------------
    
    ! --- local ----------------------------------
    
    integer                         ::  id
    
    ! --- begin ----------------------------------
    
    ! loop over data sets:
    do id = 1, self%n
      ! put out:
      call self%values(id)%p%SetPixelZero( ipix, status )
      IF_NOT_OK_RETURN(status=1)
    end do ! i

    ! ok
    status = 0
    
  end subroutine PixelDatas_SetPixelZero


  ! ***


  subroutine PixelDatas_SetPixel_0D( self, i, ipix, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  i
    integer, intent(in)                 ::  ipix
    real, intent(in)                    ::  value
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetPixel_0D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. self%initialized ) then
      write (csol,'("pixeldata not initialized yet ..")'); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( (i < 1) .or. (i > self%n) ) then
      write (csol,'("pixeldata ",i0," out of range 1:",i0)') i, self%n; call csoErr
      TRACEBACK; status=1; return
    end if

    ! fill:
    call self%values(i)%p%SetPixel( ipix, value, status )

    ! ok
    status = 0
    
  end subroutine PixelDatas_SetPixel_0D


  ! ***


  subroutine PixelDatas_SetPixel_1D( self, i, ipix, values, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  i
    integer, intent(in)                 ::  ipix
    real, intent(in)                    ::  values(:)
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetPixel_1D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. self%initialized ) then
      write (csol,'("pixeldata not initialized yet ..")'); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( (i < 1) .or. (i > self%n) ) then
      write (csol,'("pixeldata ",i0," out of range 1:",i0)') i, self%n; call csoErr
      TRACEBACK; status=1; return
    end if

    ! fill:
    call self%values(i)%p%SetPixel( ipix, values, status )

    ! ok
    status = 0
    
  end subroutine PixelDatas_SetPixel_1D


  ! ***


  subroutine PixelDatas_SetPixel_2D( self, i, ipix, values, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  i
    integer, intent(in)                 ::  ipix
    real, intent(in)                    ::  values(:,:)
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetPixel_2D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. self%initialized ) then
      write (csol,'("pixeldata not initialized yet ..")'); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( (i < 1) .or. (i > self%n) ) then
      write (csol,'("pixeldata ",i0," out of range 1:",i0)') i, self%n; call csoErr
      TRACEBACK; status=1; return
    end if

    ! fill:
    call self%values(i)%p%SetPixel( ipix, values, status )

    ! ok
    status = 0
    
  end subroutine PixelDatas_SetPixel_2D


  ! ***


  subroutine PixelDatas_SetFormula( self, i, formula, formula_terms, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    integer, intent(in)                 ::  i
    character(len=*), intent(in)        ::  formula
    character(len=*), intent(in)        ::  formula_terms
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_SetFormula'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (i < 1) .or. (i > self%n) ) then
      write (csol,'("pixeldata ",i0," out of range 1:",i0)') i, self%n; call csoErr
      TRACEBACK; status=1; return
    end if

    ! fill:
    self%values(i)%p%formula       = trim(formula)
    self%values(i)%p%formula_terms = trim(formula_terms)

    ! ok
    status = 0
    
  end subroutine PixelDatas_SetFormula


  ! ***


  subroutine PixelDatas_ApplyFormulas( self, pd, status )

    use CSO_Profile, only : T_ProfileMapping
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)  ::  self
    type(T_PixelDatas), intent(in)      ::  pd   ! data (kernels etc)
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_ApplyFormulas'

    ! gravity constant:
    real, parameter     ::  grav   = 9.80665 ! m/s2
    ! mole mass of air:
    real, parameter     ::  xm_air = 28.964e-3     ! kg/mol : ~80% N2, ~20% O2, ... 
     
    ! --- local ----------------------------------
    
    integer                         ::  id
    type(T_PixelData), pointer      ::  p
    integer                         ::  npix
    integer                         ::  ipix

    real, pointer                   ::  y_data(:,:)  ! (*,npix)
    character(len=64)               ::  y_units

    integer                         ::  nlayer
    integer                         ::  nz
    real, pointer                   ::  hp_data(:,:)  ! (nlayer,npix)
    character(len=64)               ::  hp_units
    integer                         ::  hp_isfc
    real, pointer                   ::  mod_hp_data(:,:)  ! (nz,npix)
    character(len=64)               ::  mod_hp_units
    integer                         ::  mod_hp_isfc
    real, pointer                   ::  mod_conc_data(:,:)  ! (nz,npix)
    character(len=64)               ::  mod_conc_units
    type(T_ProfileMapping)          ::  ProfileMapping
    real, allocatable               ::  mod_hpx(:)   ! (0:nz)
    real, allocatable               ::  hx(:)       ! (nlayer)

    real, pointer                   ::  A_data(:,:,:)  ! (nretr,nlayer,npix)
    character(len=64)               ::  A_units
    real, pointer                   ::  x_data(:,:)  ! (nlayer,npix)
    character(len=64)               ::  x_units

    ! --- begin ----------------------------------
    
    ! loop over datasets:
    do id = 1, self%n
      ! short:
      p => self%values(id)%p
      ! only if defined ...
      if ( len_trim(p%formula) == 0 ) cycle
      ! info ...
      write (csol,'(a,": fill `",a,"` using formula: ",a)') rname, trim(p%name), trim(p%formula); call csoPr
      ! switch ..
      select case ( trim(p%formula) )
        
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !~ vertical mapping
        case ( 'LayerAverage( hp, mod_hp, mod_conc )' )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! pointer to target array:
          call self%GetData( status, id=id, data1=y_data, units=y_units )
          IF_NOT_OK_RETURN(status=1)
          ! get pointers to source arrays:
          call self%GetFormulaData( p%formula_terms, 'hp', status, pd=pd, data1=hp_data, &
                                     units=hp_units, isfc_pressure=hp_isfc )
          IF_NOT_OK_RETURN(status=1)
          call self%GetFormulaData( p%formula_terms, 'mod_hp', status, pd=pd, data1=mod_hp_data, &
                                     units=mod_hp_units, isfc_pressure=mod_hp_isfc )
          IF_NOT_OK_RETURN(status=1)
          call self%GetFormulaData( p%formula_terms, 'mod_conc', status, pd=pd, data1=mod_conc_data, units=mod_conc_units )
          IF_NOT_OK_RETURN(status=1)
          
          ! check units:
          if ( trim(hp_units) /= trim(mod_hp_units) ) then
            write (csol,'("hp units `",a,"` should be equal to mod_hp units `",a,"`")') trim(hp_units), trim(mod_hp_units); call csoPr
            write (csol,'("  formula       : ",a)') trim(p%formula); call csoPr
            write (csol,'("  formula_terms : ",a)') trim(p%formula_terms); call csoPr
            write (csol,'("  variable      : ",a)') trim(p%name); call csoPr
            TRACEBACK; status=1; return
          end if

          ! number of apriori layers, half-level pressures on 1:nlayer+1
          nlayer = size(hp_data,1) - 1
          npix   = size(hp_data,2)
          ! number of model layers, half-level pressures on 1:nz+1
          nz = size(mod_hp_data,1) - 1
          ! init vertical mapping from "mod_hp" layers to "hp" layers:
          call ProfileMapping%Init( nz, nlayer, status )
          IF_NOT_OK_RETURN(status=1)
          
          ! storage:
          allocate( mod_hpx(0:nz), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( hx(nlayer), stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! select on unit conversion:
          select case ( trim(mod_conc_units)//' -> '//trim(y_units) )
              !~ mass mixing ratio to column density:
              case ( 'ppb -> mol m-2' )  
                ! loop over pixels:
                do ipix = 1, npix
                  ! skip if no-data:
                  if ( mod_hp_data(1,ipix) == p%fill_value ) cycle
                  ! scale model pressures to have same surface pressure as pixel:
                  mod_hpx = mod_hp_data(:,ipix) / mod_hp_data(mod_hp_isfc,ipix) * hp_data(hp_isfc,ipix)
                  ! compute weights for vertical mapping:
                  call ProfileMapping%Setup( mod_hpx, hp_data(:,ipix), status )
                  IF_NOT_OK_RETURN(status=1)
                  ! integral over pressure:
                  !   sum_i  c(i)  dp(i)
                  !          ppb    Pa
                  call ProfileMapping%Apply_WeightedSum( mod_conc_data(:,ipix), hx, status )
                  IF_NOT_OK_RETURN(status=1)
                  ! unit conversion:
                  y_data(:,ipix) = &         !  (mole tr)/m2 =
                                 hx &        ! ppb*Pa
                                 / grav &    ! / [g]   :  Pa/[g] = (kg air)/m2 
                                 / xm_air &  ! / ((kg air)/(mole air))
                                 * 1.0e-9    ! * (mole tr)/(mole air)/ppb
                end do ! ipix
            !~
            case default
              write (csol,'("unsupported conversion `",a,"`")') &
                         trim(mod_conc_units)//' -> '//trim(y_units); call csoErr
              TRACEBACK; status=1; return
          end select

          ! clear:
          deallocate( mod_hpx, stat=status )
          IF_NOT_OK_RETURN(status=1)
          deallocate( hx, stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! done with vertical mapping:
          call ProfileMapping%Done( status )
          IF_NOT_OK_RETURN(status=1)

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !~ kernel convolution, no apriori
        case ( 'A x' )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! pointer to target array:
          call self%GetData( status, id=id, data1=y_data, units=y_units )
          IF_NOT_OK_RETURN(status=1)
          ! get pointers to source arrays:
          call self%GetFormulaData( p%formula_terms, 'A', status, pd=pd, data2=A_data, units=A_units )
          IF_NOT_OK_RETURN(status=1)
          call self%GetFormulaData( p%formula_terms, 'x', status, pd=pd, data1=x_data, units=x_units )
          IF_NOT_OK_RETURN(status=1)
          ! check units:
          if ( trim(A_units) /= '1' ) then
            write (csol,'("A units `",a,"` should be `1`")') trim(A_units); call csoPr
            write (csol,'("  formula       : ",a)') trim(p%formula); call csoPr
            write (csol,'("  formula_terms : ",a)') trim(p%formula_terms); call csoPr
            write (csol,'("  variable      : ",a)') trim(p%name); call csoPr
            TRACEBACK; status=1; return
          end if
          if ( trim(x_units) /= trim(y_units) ) then
            write (csol,'("x units `",a,"` should be equal to y units `",a,"`")') trim(x_units), trim(y_units); call csoPr
            write (csol,'("  formula       : ",a)') trim(p%formula); call csoPr
            write (csol,'("  formula_terms : ",a)') trim(p%formula_terms); call csoPr
            write (csol,'("  variable      : ",a)') trim(p%name); call csoPr
            TRACEBACK; status=1; return
          end if
          ! size:
          npix = size(x_data,2)
          ! apply:
          do ipix = 1, npix
            ! fill:
            y_data(:,ipix) = matmul( A_data(:,:,ipix), x_data(:,ipix) )
          end do
        
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write (csol,'("unsupported formula `",a,"`")') trim(p%formula); call csoPr
          TRACEBACK; status=1; return
      end select
      
    end do  ! datasets

    ! ok
    status = 0
    
  end subroutine PixelDatas_ApplyFormulas
  
  
  ! ***
  
  
  !
  ! The  formula terms are used to obtain the actual data name for a standard name:
  !   A: A x:hx
  ! For term "x" the actual data name is "hx".
  ! Return properties for the actual variable:
  !   - data0, data1, data2 : pointers to data arrays
  !   - units : attribute
  !   - isfc  : index of surface level in case data represents pressures
  ! Eventuall supply extra "pd", this will be searched for the actual data if it is not found in self.
  !
  
  subroutine PixelDatas_GetFormulaData( self, formula_terms, term, status, &
                                          pd, data0, data1, data2, units, isfc_pressure )

    use CSO_String, only : CSO_ReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)   ::  self
    character(len=*), intent(in)      ::  formula_terms
    character(len=*), intent(in)      ::  term
    integer, intent(out)              ::  status

    type(T_PixelDatas), intent(in), optional   ::  pd
    real, pointer, optional                    ::  data0(:)      ! (npix)
    real, pointer, optional                    ::  data1(:,:)    ! (n,npix)
    real, pointer, optional                    ::  data2(:,:,:)  ! (m,n,npix)
    character(len=*), intent(out), optional    ::  units
    integer, intent(out), optional             ::  isfc_pressure

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_GetFormulaData'
    
    ! --- local ----------------------------------
    
    character(len=1024)     ::  line
    character(len=64)       ::  dterm
    character(len=64)       ::  aterm
    integer                 ::  n
    integer                 ::  id
    integer                 ::  i
    
    ! --- begin ----------------------------------

    ! copy:
    line = trim(formula_terms)
    
    ! init actual name to empty:
    aterm = ''
    ! loop over default/actual terms
    do
      ! done?
      if ( len_trim(line) == 0 ) exit

      ! extract default name:
      call CSO_ReadFromLine( line, dterm, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      ! check postfix:
      n = len_trim(dterm)
      if ( dterm(n:n) /= ':' ) then
        write (csol,'("default term `",a,"` should end with `:`")') trim(dterm); call csoErr
        write (csol,'("  formula terms: ",a)') trim(formula_terms); call csoErr
        TRACEBACK; status=1; return
      end if
      ! remove postfix:
      dterm = dterm(1:n-1)
      
      ! extract actual name:
      call CSO_ReadFromLine( line, aterm, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)

      ! match?
      if ( trim(dterm) == trim(term) ) then
        ! found:
        exit
      else
        ! reset:
        aterm = ''
      end if
    end do ! dterm/aterm pairs
    ! check ..
    if ( len_trim(aterm) == 0 ) then
      write (csol,'("could not find default term `",a,"` in formula terms: ",a)') trim(term), trim(formula_terms); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! search in self values, return status<0 if not found:
    call self%InqID( aterm, id, status, quiet=.true. )
    !~ found
    if ( status == 0 ) then
      ! get pointers and units from self:
      call self%GetData( status, id=id, data0=data0, data1=data1, data2=data2, &
                            units=units, isfc_pressure=isfc_pressure )
      IF_NOT_OK_RETURN(status=1)
    !~ not found
    else if ( status < 0 ) then
      ! alternative data?
      if ( present(pd) ) then
        ! search in additional data:
        call pd%InqID( aterm, id, status, quiet=.true. )
        !~ found 
        if ( status == 0 ) then
          ! get pointers and units from extra data:
          call pd%GetData( status, id=id, data0=data0, data1=data1, data2=data2, &
                              units=units, isfc_pressure=isfc_pressure )
          IF_NOT_OK_RETURN(status=1)
        !~ not found ... 
        else if ( status < 0 ) then
          write (csol,'("could not find term `",a,"` in pixel datas:")') trim(aterm); call csoErr
          do i = 1, self%n
            write (csol,'("    ",a)') trim(self%values(i)%p%name); call csoErr
          end do
          write (csol,'("  and neither in additional pixel datas:")'); call csoErr
          do i = 1, pd%n
            write (csol,'("    ",a)') trim(pd%values(i)%p%name); call csoErr
          end do
          write (csol,'("  formula terms: ",a)') trim(formula_terms); call csoErr
          TRACEBACK; status=1; return
        !~ error 
        else
          TRACEBACK; status=1; return
        end if
      else
        write (csol,'("could not find term `",a,"` in pixel datas:")') trim(aterm); call csoErr
        do i = 1, self%n
          write (csol,'("    ",a)') trim(self%values(i)%p%name); call csoErr
        end do
        write (csol,'("  and no additional pixel datas provide")'); call csoErr
        write (csol,'("  formula terms: ",a)') trim(formula_terms); call csoErr
        TRACEBACK; status=1; return
      end if ! pd
    !~ error 
    else
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine PixelDatas_GetFormulaData
  
  
  ! ***
  
  
  subroutine PixelDatas_NcDef( self, ncid, dimid_pixel, status )

    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(inout)    ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(in)                   ::  dimid_pixel
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_NcDef'
    
    ! --- local ----------------------------------
    
    integer                   ::  idim
    integer                   ::  i
    
    ! --- begin ----------------------------------
    
    ! define in nc file:
    call self%ncdims%Def( ncid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over data sets:
    do i = 1, self%n
      ! define variable in file:
      call self%values(i)%p%NcDef( ncid, self%ncdims, dimid_pixel, status )
      IF_NOT_OK_RETURN(status=1)
    end do ! i

    ! ok
    status = 0
    
  end subroutine PixelDatas_NcDef

  
  ! ***
  
  
  !
  ! "glb" data is available for all pixels and for all domains;
  ! from root, put subset to nc file
  !
  
  
  subroutine PixelDatas_NcPutGlbSelect( self, ncid, mapping, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)   ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  mapping(:)  ! (nglb) target index, or <0 if not used
    integer, intent(out)              ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_NcPutGlbSelect'
    
    ! --- local ----------------------------------
    
    integer                         ::  i
    
    ! --- begin ----------------------------------
    
    ! loop over data sets:
    do i = 1, self%n
      ! filter:
      if ( .not. self%values(i)%p%glb ) cycle
      ! put out:
      call self%values(i)%p%NcPutGlbSelect( ncid, mapping, status )
      IF_NOT_OK_RETURN(status=1)
    end do ! i
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_NcPutGlbSelect

  
  ! ***
  
  
  !
  ! Gather pixels from domains to root, put to nc file;
  ! skip 'glb' pixel data
  !
  
  
  subroutine PixelDatas_NcPutGather( self, ncid, mapping, add, status )
  
    ! --- in/out ---------------------------------
    
    class(T_PixelDatas), intent(in)   ::  self
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  mapping(:)   ! (npix_all) target index
    logical, intent(in)               ::  add          ! add contributions?
    integer, intent(out)              ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/PixelDatas_NcPutGather'
    
    ! --- local ----------------------------------
    
    integer                         ::  i
    
    ! --- begin ----------------------------------
    
    ! loop over data sets:
    do i = 1, self%n
      ! filter:
      if ( self%values(i)%p%glb ) cycle
      ! put out:
      call self%values(i)%p%NcPutGather( ncid, mapping, add, status )
      IF_NOT_OK_RETURN(status=1)
    end do ! i
    
    ! ok
    status = 0
    
  end subroutine PixelDatas_NcPutGather


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
    use NetCDF, only : NF90_ENOTATT
  
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
    if ( status == NF90_ENOTATT ) then
      write (csol,'("no `units` attribute found")'); call csoErr
      TRACEBACK; status=1; return
    else if ( status /= NF90_NOERR ) then
      csol=NF90_StrError(status); call csoErr
      TRACEBACK; status=1; return
    end if
    ! read attribute:
    status = NF90_Get_Att( ncid, varid, 'long_name', self%long_name )
    if ( status == NF90_ENOTATT ) then
      write (csol,'("no `long_name` attribute found")'); call csoErr
      TRACEBACK; status=1; return
    else if ( status /= NF90_NOERR ) then
      csol=NF90_StrError(status); call csoErr
      TRACEBACK; status=1; return
    end if
    
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
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( ncid, varid, 'long_name', self%long_name )
    IF_NF90_NOT_OK_RETURN(status=1)
    
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
    IF_NOT_OK_RETURN(status=1)

    ! fill value:
    self%fill_value = - huge(real(1.0,kind=wpr))
    ! storage:
    allocate( self%data(self%ntx,self%nty), source=real(self%fill_value), stat=status )
    IF_NOT_OK_RETURN(status=1)
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
    IF_NOT_OK_RETURN(status=1)
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
    if (status/=NF90_NOERR) then
      csol=nf90_strerror(status); call csoErr
      write (csol,'("variable name: ",a)') trim(varname); call csoErr
      TRACEBACK; status=1; return
    end if

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, dimids=dimids )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, ndim
      ! get size:
      status = NF90_Inquire_Dimension( ncid, dimids(idim), len=shp(idim) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do

    ! init storage:
    call self%Init( shp(1), shp(2), status )
    IF_NOT_OK_RETURN(status=1)

    ! read attributes:
    call self%NcGetAttrs( ncid, varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read:
    status = NF90_Get_Var( ncid, varid, self%data )
    IF_NF90_NOT_OK_RETURN(status=1)
    
      !! testing ...
      !write (csol,*) rname//': long_name = ', trim(self%long_name), ' ; varname = ', trim(varname); call csoPr
      !write (csol,*) rname//':   varid = ', varid; call csoPr
      !write (csol,*) rname//':   shape = ', shape(self%data); call csoPr
      !write (csol,*) rname//':   range = ', minval(self%data), maxval(self%data); call csoPr
      
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
    IF_NF90_NOT_OK_RETURN(status=1)

    ! add attributes:
    call self%NcPutAttrs( ncid, self%varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! define attribute:
    status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Track_0D_NcDef
  
  
  ! ***


  ! write global array from root, each pe has same copy
  
  subroutine Track_0D_NcPutGlb( self, ncid, status )
  
    use NetCDF  , only : NF90_Put_Var
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_Track_0D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_0D_NcPutGlb'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( csoc%root ) then
      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, self%data )
      IF_NF90_NOT_OK_RETURN(status=1)
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
    IF_NOT_OK_RETURN(status=1)
    
    ! store:
    self%nv   = shp(1)

    ! fill value:
    self%fill_value = - huge(real(1.0,kind=wpr))
    ! storage:
    allocate( self%data(self%nv,self%ntx,self%nty), source=real(self%fill_value), stat=status )
    IF_NOT_OK_RETURN(status=1)
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
    IF_NOT_OK_RETURN(status=1)
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
    if (status/=NF90_NOERR) then
      csol=nf90_strerror(status); call csoErr
      write (csol,'("variable name: ",a)') trim(varname); call csoErr
      TRACEBACK; status=1; return
    end if

    ! get dimension id's:
    status = NF90_Inquire_Variable( ncid, varid, ndims=nd, dimids=dimids )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! loop over dimensions:
    do idim = 1, ndim
      ! scalar?
      if ( idim <= ndim-nd ) then
        shp(idim) = 1
      else
        ! get size:
        status = NF90_Inquire_Dimension( ncid, dimids(idim-(ndim-nd)), len=shp(idim) )
        IF_NF90_NOT_OK_RETURN(status=1)
      end if
    end do

    ! init storage for all pixels:
    call self%Init( shp, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read attributes:
    call self%NcGetAttrs( ncid, varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read:
    status = NF90_Get_Var( ncid, varid, self%data )
    IF_NF90_NOT_OK_RETURN(status=1)

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
    IF_NF90_NOT_OK_RETURN(status=1)

    ! add attributes:
    call self%NcPutAttrs( ncid, self%varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! define attribute:
    status = NF90_Put_Att( ncid, self%varid, '_FillValue', self%fill_value )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Track_1D_NcDef
  
  
  ! ***


  ! write global array from root, each pe has same copy
  
  subroutine Track_1D_NcPutGlb( self, ncid, status )
  
    use NetCDF  , only : NF90_Put_Var
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_Track_1D), intent(in)    ::  self
    integer, intent(in)               ::  ncid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Track_1D_NcPutGlb'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! write from root only ..
    if ( csoc%root ) then
      ! write variable:
      status = NF90_Put_Var( ncid, self%varid, self%data )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if ! root
    
    ! ok
    status = 0
    
  end subroutine Track_1D_NcPutGlb
  



end module CSO_Pixels
