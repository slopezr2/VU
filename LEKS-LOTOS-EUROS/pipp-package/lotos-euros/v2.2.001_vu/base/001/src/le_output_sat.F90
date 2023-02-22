!###############################################################################
!
! LE_Output_Sat - simulations of satellite retrievals
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
  
module LE_Output_Sat

  use GO              , only : gol, goPr, goErr
  use NetCDF          , only : NF90_StrError, NF90_NOERR
  use LE_Output_Common, only : T_LE_Output_Common
  use Pixels          , only : T_Pixels_0D_i, T_Pixels_0D, T_Pixels_1D, T_Pixels_2D
  use Pixels          , only : T_Track_0D, T_Track_1D

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_LE_Output_Sat_Data
  public  ::  T_LE_Output_Sat_State
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'LE_Output_Sat'

  ! units used for (intermediate) vertical column densities:
  character(len=*), parameter   ::  vcd_units = '1e15 mlc/cm2'

  ! --- types --------------------------------
  
  !
  ! satelite data:
  !  footprints
  !  timestamps
  !  kernels
  !
  type T_LE_Output_Sat_Data
    ! key to identify this set:
    character(len=32)                 ::  name
    ! common stuff:
    type(T_LE_Output_Common)          ::  com
    ! tracer name:
    character(len=32)                 ::  tracer_name
    ! filename template:
    character(len=1024)               ::  filename
    ! variable names:
    character(len=32)                 ::  varname_hp
    character(len=32)                 ::  varname_ya
    character(len=32)                 ::  varname_yr
    character(len=32)                 ::  varname_sr
    character(len=32)                 ::  varname_vr
    ! flags:
    logical                           ::  with_corners
    logical                           ::  with_profile
    logical                           ::  with_apriori
    !
    !~ global pixel arrays
    !
    ! number of pixels in global doamin:
    integer                           ::  nglb
    ! number of corners in footprint:
    integer                           ::  ncorner
    ! number of layers in profiles:
    integer                           ::  nlayer
    ! number of layer interfaces:
    integer                           ::  nlayeri  ! nlayer+1
    ! base info (footprint)
    type(T_Pixels_0D)                 ::  glb_lon    ! (nglb)
    type(T_Pixels_0D)                 ::  glb_lat    ! (nglb)
    type(T_Pixels_1D)                 ::  glb_clons  ! (ncorner,nglb)
    type(T_Pixels_1D)                 ::  glb_clats  ! (ncorner,nglb)
    ! track info (2D footprint)
    type(T_Pixels_0D_i)               ::  glb_tpixel       ! (nglb)
    integer                           ::  ntx, nty         ! cross track, along track
    type(T_Track_0D)                  ::  glb_track_lon    ! (ntx,nty)
    type(T_Track_0D)                  ::  glb_track_lat    ! (ntx,nty)
    type(T_Track_1D)                  ::  glb_track_clons  ! (ncorner,ntx,nty)
    type(T_Track_1D)                  ::  glb_track_clats  ! (ncorner,ntx,nty)
    ! 
    !~ local pixels
    !  (overlap with local domain)
    !
    ! number of pixels loaded:
    integer                           ::  npix
    !
    ! index in global arrays:
    integer, allocatable              ::  iglb(:)   ! (npix)
    !
    ! for now: grid cell indices (i,j) with pixel center,
    ! to be changed with list with cell indices and weights 
    type(T_Pixels_0D_i)               ::  ilon    ! (npix)
    type(T_Pixels_0D_i)               ::  ilat    ! (npix)
    !
    ! analysis status flag:
    type(T_Pixels_0D_i)               ::  astatus    ! (npix)
    !
    ! number of retrieval values:
    ! (1 for column, >1 for profile)
    integer                           ::  nr
    ! half-level pressures:
    type(T_Pixels_1D)                 ::  hp   ! (nlayeri,npix)
    ! apriori profiles:
    type(T_Pixels_1D)                 ::  ya   ! (nlayer,npix)
    ! averaging kernels:
    type(T_Pixels_2D)                 ::  K    ! (nr,nlayer,npix)
    ! retrieval apriori (profiles):
    type(T_Pixels_1D)                 ::  yra   ! (nr,npix)
    ! retrieval (profiles):
    type(T_Pixels_1D)                 ::  yr    ! (nr,npix)
    ! retrieval error std.dev.
    type(T_Pixels_1D)                 ::  sr    ! (nr,npix)
    ! retrieval error covariance (profiles)
    type(T_Pixels_2D)                 ::  vr    ! (nr,nr,npix)
    !
    ! how to fill 'sr' ?
    !  'data'   : from file, eventually with extra scaling
    !  'frac'   : as fraction of observed value
    character(len=4)                  ::  r_type
    ! extra scaling factor for the 'data' values:
    real                              ::  r_data_scaling
    ! for 'frac':  r = frac * y, and bounded to range
    real                              ::  r_frac_factor
    real                              ::  r_frac_min
    real                              ::  r_frac_max
    !
  contains
    procedure :: Init            => LE_Output_Sat_Data_Init
    procedure :: Done            => LE_Output_Sat_Data_Done
    procedure :: Clear           => LE_Output_Sat_Data_Clear
    procedure :: GetPixel        => LE_Output_Sat_Data_GetPixel
    procedure :: Setup           => LE_Output_Sat_Data_Setup
    procedure :: PutOut          => LE_Output_Sat_Data_PutOut
  end type T_LE_Output_Sat_Data
  
  !
  ! state simulations
  !
  type T_LE_Output_Sat_State
    ! annote:
    character(len=16)                 ::  key
    character(len=256)                ::  description
    ! flags:
    logical                           ::  with_profile
    ! number of pixels stored:
    integer                           ::  npix
    ! number of levels:
    integer                           ::  nlevi
    integer                           ::  nlev
    ! retrievals:
    integer                           ::  nlayer
    integer                           ::  nr
    ! model profiles:
    type(T_Pixels_1D)                 ::  hp    ! (nlevi,npix)
    type(T_Pixels_1D)                 ::  vmr   ! (nlev,npix)
    type(T_Pixels_1D)                 ::  vcd   ! (nlev,npix)
    ! retrieval profiles:
    type(T_Pixels_1D)                 ::  hx   ! (nlayer,npix)
    ! retrieved profiles:
    type(T_Pixels_1D)                 ::  y    ! (nr,npix)
  contains
    procedure :: Init            => LE_Output_Sat_State_Init
    procedure :: Done            => LE_Output_Sat_State_Done
    procedure :: Clear           => LE_Output_Sat_State_Clear
    procedure :: GetPixel        => LE_Output_Sat_State_GetPixel
    procedure :: Setup           => LE_Output_Sat_State_Setup
    procedure :: PutOut          => LE_Output_Sat_State_PutOut
  end type T_LE_Output_Sat_State


  
contains


  ! ====================================================================
  ! ===
  ! === Data
  ! ===
  ! ====================================================================


  subroutine LE_Output_Sat_Data_Init( self, rcF, rcbase, typ, name, status )
  
    use GO              , only : TrcFile
    use LE_Output_Common, only : LE_Output_Common_Init

    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_Data), intent(out)    ::  self
    type(TrcFile), intent(in)                   ::  rcF
    character(len=*), intent(in)                ::  rcbase
    character(len=*), intent(in)                ::  typ
    character(len=*), intent(in)                ::  name
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Output_Sat_Data_Init'
    
    ! --- local ----------------------------------
    
    character(len=64)   ::  rckey
    
    ! --- begin ----------------------------------
    
    ! store:
    self%name = trim(name)
    
    ! prefix:
    rckey = trim(rcbase)//'.'//trim(typ)//'.'//trim(name)
    
    ! init common stuff:
    call LE_Output_Common_Init( self%com, rcF, rcbase, status )
    IF_NOTOK_RETURN(status=1)
    
    ! tracer name:
    call rcF%Get( trim(rckey)//'.tracer', self%tracer_name, status )
    IF_NOTOK_RETURN(status=1)
    
    ! flags:
    call rcF%Get( trim(rckey)//'.with_corners', self%with_corners, status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( trim(rckey)//'.with_profile', self%with_profile, status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( trim(rckey)//'.with_apriori', self%with_apriori, status )
    IF_NOTOK_RETURN(status=1)
    ! filename template:
    call rcF%Get( trim(rckey)//'.filename', self%filename, status )
    IF_NOTOK_RETURN(status=1)
    
    ! variable names:
    call rcF%Get( trim(rckey)//'.varname.hp', self%varname_hp, status )
    IF_NOTOK_RETURN(status=1)
    ! apri?
    if ( self%with_apriori ) then
      call rcF%Get( trim(rckey)//'.varname.ya', self%varname_ya, status )
    IF_NOTOK_RETURN(status=1)
    end if
    ! retrieval:
    call rcF%Get( trim(rckey)//'.varname.yr', self%varname_yr, status )
    IF_NOTOK_RETURN(status=1)
    ! profile or single value?
    if ( self%with_profile ) then
      call rcF%Get( trim(rckey)//'.varname.vr', self%varname_vr, status )
      IF_NOTOK_RETURN(status=1)
    else
      call rcF%Get( trim(rckey)//'.varname.sr', self%varname_sr, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! no pixels yet:
    self%nglb = 0
    self%npix = 0

    ! how to define representation error?
    call rcF%Get( trim(rckey)//'.r.type', self%r_type, status )
    IF_NOTOK_RETURN(status=1)
    ! switch:
    select case ( trim(self%r_type) )
      !~ from data?
      case ( 'data' )
        ! extra scaling factor:
        call rcF%Get( trim(rckey)//'.r.data.scaling', self%r_data_scaling, status, default=1.0 )
        IF_ERROR_RETURN(status=1)
      !~ fraction of data:
      case ( 'frac' )
        ! scaling and bounds:
        call rcF%Get( trim(rckey)//'.r.frac.factor', self%r_frac_factor, status )
        IF_NOTOK_RETURN(status=1)
        call rcF%Get( trim(rckey)//'.r.frac.min'   , self%r_frac_min   , status )
        IF_NOTOK_RETURN(status=1)
        call rcF%Get( trim(rckey)//'.r.frac.max'   , self%r_frac_max   , status )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unsupported r type `",a,"`")') trim(self%r_type); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_Data_Init


  ! ***


  subroutine LE_Output_Sat_Data_Done( self, status )

    use LE_Output_Common, only : LE_Output_Common_Done
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_Data), intent(inout)    ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_Data_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%Clear( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_Data_Done


  ! ***


  subroutine LE_Output_Sat_Data_Clear( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_Data), intent(inout)    ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_Data_Clear'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! any pixels in global domain?
    if ( self%nglb > 0 ) then
    
      ! clear:
      call self%glb_tpixel%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%glb_track_lon%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%glb_track_lat%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! corners defined?
      if ( self%with_corners ) then
        call self%glb_track_clons%Done( status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_track_clats%Done( status )
        IF_NOTOK_RETURN(status=1)
      end if
      
      ! clear:
      call self%glb_lon%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%glb_lat%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! corners defined?
      if ( self%with_corners ) then
        call self%glb_clons%Done( status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_clats%Done( status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! clear:
      deallocate( self%iglb, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! no pixels anymore:
      self%nglb = 0
    
      !~ arrays with local pixels, npix could be zero ...

      ! clear:
      call self%ilon%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%ilat%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      call self%astatus%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      call self%hp%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! apri?
      if ( self%with_apriori ) then
        call self%ya%Done( status )
        IF_NOTOK_RETURN(status=1)
        call self%yra%Done( status )
        IF_NOTOK_RETURN(status=1)
      end if
      call self%K%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%yr%Done( status )
      IF_NOTOK_RETURN(status=1)
      if ( self%with_profile ) then
        ! covariance:
        call self%vr%Done( status )
        IF_NOTOK_RETURN(status=1)
      else
        call self%sr%Done( status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! no pixels anymore:
      self%npix = 0
    end if

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_Data_Clear


  ! ***


  subroutine LE_Output_Sat_Data_GetPixel( self, ipix, status, &
                                            glbid, lon, lat, y, sigma, covar )
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_Data), intent(in)       ::  self
    integer, intent(in)                           ::  ipix
    integer, intent(out)                          ::  status
    
    integer, intent(out), optional                ::  glbid
    real, intent(out), optional                   ::  lon
    real, intent(out), optional                   ::  lat
    real, intent(out), optional                   ::  y
    real, intent(out), optional                   ::  sigma
    real, intent(out), optional                   ::  covar(:,:)  ! (nr,nr)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_Data_GetPixel'
    
    ! --- local ----------------------------------
    
    integer     ::  iglb
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ipix < 1) .or. (ipix > self%npix) ) then
      write (gol,'("pixel ",i0," out of range 1:",i0)') ipix, self%npix; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! global index:
    iglb = self%iglb(ipix)
    
    ! id?
    if ( present(glbid) ) glbid = iglb
    
    ! location:
    if ( present(lon) ) lon = self%glb_lon%data(iglb)
    if ( present(lat) ) lat = self%glb_lat%data(iglb)
    
    ! retrieval of total column: SHOULD BE REPLACED BY PROFILE!
    if ( present(y) ) then
      ! check ...
      if ( self%nr /= 1 ) then
        write (gol,'("profiles not supported yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! first value ...
      y = self%yr%data(1,ipix)
    end if
    
    ! std.dev.: SHOULD BE REPLACED BY COVARIANCE!
    if ( present(sigma) ) then
      ! check ...
      if ( self%nr /= 1 ) then
        write (gol,'("profiles not supported yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! first value ...
      sigma = self%sr%data(1,ipix)
    end if
    
    ! std.dev.: SHOULD BE REPLACED BY COVARIANCE!
    if ( present(covar) ) then
      ! check ..
      if ( any( shape(covar) /= (/self%nr,self%nr/) ) ) then
        write (gol,'("covar shape (",i0,",",i0,") while nr=",i0)') shape(covar), self%nr; call goPr
        TRACEBACK; status=1; return
      end if
      ! copy:
      covar = self%vr%data(:,:,ipix)
    end if

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_Data_GetPixel


  ! ***


  subroutine LE_Output_Sat_Data_Setup( self, tr, status )  
    
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    
    use Binas  , only : Avog, xm_air, grav

    use GO     , only : TDate, wrtgol
    use GO     , only : goUpCase, goLoCase
    use GO     , only : goReplace
    use LE_Grid, only : ugg
    use LE_Grid, only : dom
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_Data), intent(inout)    ::  self
    type(TDate), intent(in)                       ::  tr(2)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_Data_Setup'
    
    ! --- local ----------------------------------
    
    character(len=1024)     ::  fname
    logical                 ::  exist
    integer                 ::  ncid
    integer                 ::  dimid
    integer                 ::  varid
    integer, allocatable    ::  iglb_all(:)  ! (nglb)
    integer                 ::  ipix
    integer                 ::  iglb
    logical                 ::  indomain
    integer                 ::  off(2)
    
    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': setup for ', tr ); call goPr
    
    ! clear current storage:
    call self%Clear( status )
    IF_NOTOK_RETURN(status=1)

    ! start filename with template:
    fname = trim(self%filename)
    ! replace some values if necessary, use end time:
    call goReplace( fname, '%{yyyy}'  , '(i4.4)', tr(2)%year , status )
    IF_ERROR_RETURN(status=1)
    call goReplace( fname, '%{mm}'    , '(i2.2)', tr(2)%month, status )
    IF_ERROR_RETURN(status=1)
    call goReplace( fname, '%{dd}'    , '(i2.2)', tr(2)%day  , status ) 
    IF_ERROR_RETURN(status=1)
    call goReplace( fname, '%{hh}'    , '(i2.2)', tr(2)%hour , status ) 
    IF_ERROR_RETURN(status=1)
    call goReplace( fname, '%{mn}'    , '(i2.2)', tr(2)%min  , status ) 
    IF_ERROR_RETURN(status=1)
    call goReplace( fname, '%{TRACER}', goUpCase(trim(self%tracer_name)), status )
    IF_ERROR_RETURN(status=1)
    call goReplace( fname, '%{tracer}', goLoCase(trim(self%tracer_name)), status )
    IF_ERROR_RETURN(status=1)
   
    ! check ...
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      ! info ...
      write (gol,'(a,": WARNING - no file: ",a)') rname, trim(fname); call goPr
      
      ! no data:
      self%nglb = 0
      
    else

      ! info ...
      write (gol,'(a,": open file: ",a)') rname, trim(fname); call goPr

      ! open file:
      status = NF90_Open( trim(fname), NF90_NOWRITE, ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! pixel dimension:
      status = NF90_Inq_DimID( ncid, 'pixel', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! number of pixels
      status = NF90_Inquire_Dimension( ncid, dimid, len=self%nglb )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! track dimension:
      status = NF90_Inq_DimID( ncid, 'track_pixel', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! number of pixels
      status = NF90_Inquire_Dimension( ncid, dimid, len=self%ntx )
      IF_NF90_NOTOK_RETURN(status=1)
      ! track dimension:
      status = NF90_Inq_DimID( ncid, 'track_image', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! number of pixels
      status = NF90_Inquire_Dimension( ncid, dimid, len=self%nty )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! corners defined?
      if ( self%with_corners ) then
        ! corner dimension:
        status = NF90_Inq_DimID( ncid, 'corner', dimid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! number of pixels
        status = NF90_Inquire_Dimension( ncid, dimid, len=self%ncorner )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
      
      ! number of retrieval values:
      if ( self%with_profile ) then
        ! retrieval level dimension
        ! (SHOULD BE DEFINED WITH A DIFFERENT NAME ?)
        status = NF90_Inq_DimID( ncid, 'layer', dimid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! number of retrieval layers:
        status = NF90_Inquire_Dimension( ncid, dimid, len=self%nr )
        IF_NF90_NOTOK_RETURN(status=1)
      else
        ! column only:
        self%nr = 1
      end if

      ! layer dimension (used for kernel):
      status = NF90_Inq_DimID( ncid, 'layer', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! size:
      status = NF90_Inquire_Dimension( ncid, dimid, len=self%nlayer )
      IF_NF90_NOTOK_RETURN(status=1)

      ! number of half levels:
      status = NF90_Inq_DimID( ncid, 'layer_interface', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! size:
      status = NF90_Inquire_Dimension( ncid, dimid, len=self%nlayeri )
      IF_NF90_NOTOK_RETURN(status=1)

      ! info ...
      write (gol,'(a,":   shape of track: ",i0," x ",i0)') rname, self%ntx, self%nty; call goPr
      
      ! read global arrays from file:
      call self%glb_track_lon%NcInit( ncid, 'track_longitude', status )
      IF_NOTOK_RETURN(status=1)
      call self%glb_track_lat%NcInit( ncid, 'track_latitude', status )
      IF_NOTOK_RETURN(status=1)
      ! corners defined?
      if ( self%with_corners ) then
        call self%glb_track_clons%NcInit( ncid, 'track_corner_longitudes', status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_track_clats%NcInit( ncid, 'track_corner_latitudes', status )
        IF_NOTOK_RETURN(status=1)
      end if
      
      ! read compressed index in track:
      call self%glb_tpixel%NcInit( ncid, 'pixel', status )
      IF_NOTOK_RETURN(status=1)

      ! info ...
      write (gol,'(a,":   number of pixels in file: ",i0)') rname, self%nglb; call goPr
      
      ! read global arrays from file:
      call self%glb_lon%NcInit( ncid, 'longitude', status )
      IF_NOTOK_RETURN(status=1)
      call self%glb_lat%NcInit( ncid, 'latitude', status )
      IF_NOTOK_RETURN(status=1)
      ! corners defined?
      if ( self%with_corners ) then
        call self%glb_clons%NcInit( ncid, 'corner_longitudes', status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_clats%NcInit( ncid, 'corner_latitudes', status )
        IF_NOTOK_RETURN(status=1)
      end if
      
      ! temporary mapping at maximum size, init with no-data:
      allocate( iglb_all(self%nglb), source=-999, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over all pixels:
      do ipix = 1, self%nglb
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! for the moment, only select on on center location ...
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! check if center is in local domain:
        call ugg%InDomain( self%glb_lon%data(ipix), self%glb_lat%data(ipix), indomain, status )
        IF_NF90_NOTOK_RETURN(status=1)
        ! in local domain?
        if ( indomain ) then
          ! increase counter for local pixels:
          self%npix = self%npix + 1
          ! store index:
          iglb_all(self%npix) = ipix
        end if
      end do ! glb pixels
      
      ! storage for cell indices, npix could be zero ..
      call self%ilon%Init( self%npix, status, units='1', long_name='longitude cell index of pixel center' )
      IF_NOTOK_RETURN(status=1)
      call self%ilat%Init( self%npix, status, units='1', long_name='latitude cell index of pixel center' )
      IF_NOTOK_RETURN(status=1)

      ! any local pixels?
      if ( self%npix > 0 ) then

        ! storage for mapping:
        allocate( self%iglb(self%npix), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! copy from array with global length:
        self%iglb = iglb_all(1:self%npix)
        
        ! global offset for grid indices:
        call dom%Get( status, off=off )
        IF_NOTOK_RETURN(status=1)
        ! loop over local pixels:
        do ipix = 1, self%npix
          ! global index:
          iglb = self%iglb(ipix)
          ! fill ..
          call ugg%GetLocation( self%glb_lon%data(iglb), self%glb_lat%data(iglb), &
                                  self%ilon%data(ipix), self%ilat%data(ipix), status )
          IF_NOTOK_RETURN(status=1)
          ! add offset:
          self%ilon%data(ipix) = self%ilon%data(ipix) + off(1)
          self%ilat%data(ipix) = self%ilat%data(ipix) + off(2)
        end do

      else
    
        ! dummy ...
        allocate( self%iglb(1), stat=status )
        IF_NOTOK_RETURN(status=1)
        
      end if
        
      ! storage for status flags, default to 0 (analyse?)
      call self%astatus%Init( self%npix, status, source=0, units='1', long_name='assimilation status' )
      IF_NOTOK_RETURN(status=1)

      ! a-priori profiles in product?
      if ( self%with_apriori ) then
        write (gol,'("apriori profiles not supported yet")'); call goErr
        TRACEBACK; status=1; return
      else
        ! filled with fill_value ..
        call self%ya%Init( (/self%nlayer,self%npix/), status )
        IF_NOTOK_RETURN(status=1)
        ! filled with fill_value ..
        call self%yra%Init( (/self%nr,self%npix/), status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! read half-level pressures (local subset):
      call self%hp%NcInit( ncid, trim(self%varname_hp), status, &
                             nselect=self%npix, select=self%iglb )
      IF_NOTOK_RETURN(status=1)

      ! read kernels (local subset):
      call self%K%NcInit( ncid, 'kernel', status, &
                            nselect=self%npix, select=self%iglb )
      IF_NOTOK_RETURN(status=1)

      ! read retrieval profiles (local subset):
      call self%yr%NcInit( ncid, trim(self%varname_yr), status, &
                             nselect=self%npix, select=self%iglb )
      IF_NOTOK_RETURN(status=1)

      ! switch:
      select case ( trim(self%r_type) )
        !~ from data?
        case ( 'data' )
          ! profiles ?
          if ( self%with_profile ) then
            ! retrieval covariance (local subset):
            call self%vr%NcInit( ncid, trim(self%varname_vr), status, &
                                   nselect=self%npix, select=self%iglb )
            IF_NOTOK_RETURN(status=1)
            ! any local pixels?
            if ( self%npix > 0 ) then
              ! apply optional scaling:
              self%vr%data = self%vr%data * self%r_data_scaling**2
            end if
          else
            ! retrieval std.dev. (local subset):
            call self%sr%NcInit( ncid, trim(self%varname_sr), status, &
                                 nselect=self%npix, select=self%iglb )
            IF_NOTOK_RETURN(status=1)
            ! any local pixels?
            if ( self%npix > 0 ) then
              ! apply optional scaling:
              self%sr%data = self%sr%data * self%r_data_scaling
            end if
          end if
        !~ fraction of data:
        case ( 'frac' )
          ! profiles ?
          if ( self%with_profile ) then
            write (gol,'("profiles not supported yet")'); call goErr
            TRACEBACK; status=1; return
          else
            ! init storage:
            call self%sr%Init( (/self%nr,self%npix/), status )
            IF_NOTOK_RETURN(status=1)
            ! any local pixels?
            if ( self%npix > 0 ) then
              ! fraction of observed value, bounded:
              self%sr%data = min( max( self%r_frac_min, self%yr%data * self%r_frac_factor ), self%r_frac_max )
            end if
          end if
        !~
        case default
          write (gol,'("unsupported r type `",a,"`")') trim(self%r_type); call goErr
          TRACEBACK; status=1; return
      end select

      ! clear:
      deallocate( iglb_all, stat=status )
      IF_NOTOK_RETURN(status=1)

      ! close:
      status = NF90_Close( ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      
    end if  ! file found

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_Data_Setup


  ! ***


  subroutine LE_Output_Sat_Data_PutOut( self, t, status )

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim
    use NetCDF , only : NF90_EndDef
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER
  
    use GO              , only : TDate, wrtgol
    use GO              , only : goc
    use LE_Output_Common, only : PutOut_GlobalAttributes
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_Data), intent(inout)    ::  self
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_Data_PutOut'
    
    ! --- local ----------------------------------
    
    character(len=1024)       ::  fname
    integer                   ::  cmode
    integer                   ::  ncid
    integer                   ::  dimid_tx, dimid_ty
    integer                   ::  dimid_pixel
    integer                   ::  dimid_corner
    integer                   ::  dimid_layer
    integer                   ::  dimid_layeri
    integer                   ::  dimid_retr
    integer                   ::  ntot
    integer, allocatable      ::  iglb_all(:)
    
    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': put out for ', t ); call goPr
    
    ! any data?
    if ( self%nglb > 0 ) then
    
      ! total number of pixels handled by local domains;
      ! this will be >= nglb since some footprints cover multiple domains:
      call goc%ParInfo( self%npix, status, ntot=ntot )
      IF_NOTOK_RETURN(status=1)

      ! collect on root, this is much faster than parallel write ...
      if ( goc%root ) then

        ! target file:
        write (fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",2i2.2,".nc")') &
             trim(self%com%outdir), trim(self%com%model), trim(self%com%expid), &
             trim(self%name), t%year, t%month, t%day, t%hour, t%min

        ! info ..
        write (gol,'(a,":   create ",a," ...")') rname, trim(fname); call goPr

        ! set creation mode flag:
        if ( self%com%replace ) then
          cmode = NF90_CLOBBER       ! overwrite existing files
        else
          cmode = NF90_NOCLOBBER     ! do not overwrite existing files
        end if

        ! create file:
        status = NF90_Create( fname, cmode, ncid )
        if ( status /= 0 ) then
           write (gol,'("creating file :")'); call goErr
           write (gol,'("  ",a)') trim(fname); call goErr
           TRACEBACK; status=1; return
        end if

        ! write global attributes:
        call PutOut_GlobalAttributes( self%com, ncid, status )
        IF_NOTOK_RETURN(status=1)

        ! define dimensions:
        status = NF90_Def_Dim( ncid, 'pixel', self%nglb, dimid_pixel )
        IF_NF90_NOTOK_RETURN(status=1)
        if ( self%with_corners ) then
          status = NF90_Def_Dim( ncid, 'corner', self%ncorner, dimid_corner )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
        status = NF90_Def_Dim( ncid, 'track_pixel', self%ntx, dimid_tx )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'track_image', self%nty, dimid_ty )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'layer', self%nlayer, dimid_layer )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'layer_interface', self%nlayeri, dimid_layeri )
        IF_NF90_NOTOK_RETURN(status=1)
        ! retrieval layers:
        status = NF90_Def_Dim( ncid, 'retr', self%nr, dimid_retr )
        IF_NF90_NOTOK_RETURN(status=1)

        ! define variables:
        !~ global arrays:
        call self%glb_lon%NcDef( ncid, 'longitude', (/dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_lat%NcDef( ncid, 'latitude', (/dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        ! corners defined?
        if ( self%with_corners ) then
          call self%glb_clons%NcDef( ncid, 'corner_longitudes', (/dimid_corner,dimid_pixel/), status )
          IF_NOTOK_RETURN(status=1)
          call self%glb_clats%NcDef( ncid, 'corner_latitudes' , (/dimid_corner,dimid_pixel/), status )
          IF_NOTOK_RETURN(status=1)
        end if
        !~ global track arrays:
        call self%glb_tpixel%NcDef( ncid, 'pixel', (/dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_track_lon%NcDef( ncid, 'track_longitude', (/dimid_tx,dimid_ty/), status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_track_lat%NcDef( ncid, 'track_latitude', (/dimid_tx,dimid_ty/), status )
        IF_NOTOK_RETURN(status=1)
        ! corners defined?
        if ( self%with_corners ) then
          call self%glb_track_clons%NcDef( ncid, 'track_corner_longitudes', (/dimid_corner,dimid_tx,dimid_ty/), status )
          IF_NOTOK_RETURN(status=1)
          call self%glb_track_clats%NcDef( ncid, 'track_corner_latitudes' , (/dimid_corner,dimid_tx,dimid_ty/), status )
          IF_NOTOK_RETURN(status=1)
        end if
        !~ distributed arrays, will be collected on write:
        call self%hp%NcDef( ncid, 'hp' , (/dimid_layeri,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%yr%NcDef( ncid, 'yr' , (/dimid_retr,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        if ( self%with_profile ) then
          !call self%yr_vcd%NcDef( ncid, 'yr_vcd', (/dimid_retr,dimid_pixel/), status )
          !IF_NOTOK_RETURN(status=1)
          call self%vr%NcDef( ncid, 'covar', (/dimid_retr,dimid_retr,dimid_pixel/), status )
          IF_NOTOK_RETURN(status=1)
        else
          call self%sr%NcDef( ncid, 'sigma' , (/dimid_retr,dimid_pixel/), status )
          IF_NOTOK_RETURN(status=1)
        end if
        call self%K%NcDef( ncid, 'K' , (/dimid_retr,dimid_layer,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        !~ including a priori profiles?
        if ( self%with_apriori ) then
          call self%ya%NcDef( ncid, 'ya' , (/dimid_layer,dimid_pixel/), status )
          IF_NOTOK_RETURN(status=1)
          call self%yra%NcDef( ncid, 'yra', (/dimid_retr,dimid_pixel/), status )
          IF_NOTOK_RETURN(status=1)
        end if
        !~ testing ...
        call self%ilon%NcDef( ncid, 'ilon', (/dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%ilat%NcDef( ncid, 'ilat', (/dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        
        ! analysis status:
        call self%astatus%NcDef( ncid, 'astatus', (/dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)

        ! end defintion mode:
        status = NF90_EndDef( ncid )
        IF_NF90_NOTOK_RETURN(status=1)

        ! write arrays available as global arrays:
        call self%glb_lon%NcPutGlb( ncid, status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_lat%NcPutGlb( ncid, status )
        IF_NOTOK_RETURN(status=1)
        ! corners defined?
        if ( self%with_corners ) then
          call self%glb_clons%NcPutGlb( ncid, status )
          IF_NOTOK_RETURN(status=1)
          call self%glb_clats%NcPutGlb( ncid, status )
          IF_NOTOK_RETURN(status=1)
        end if
        
        ! write track arrays available as global arrays:
        call self%glb_tpixel%NcPutGlb( ncid, status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_track_lon%NcPutGlb( ncid, status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_track_lat%NcPutGlb( ncid, status )
        IF_NOTOK_RETURN(status=1)
        ! corners defined?
        if ( self%with_corners ) then
          call self%glb_track_clons%NcPutGlb( ncid, status )
          IF_NOTOK_RETURN(status=1)
          call self%glb_track_clats%NcPutGlb( ncid, status )
          IF_NOTOK_RETURN(status=1)
        end if
        
        ! storage for mapping:
        allocate( iglb_all(ntot), stat=status )
        IF_NOTOK_RETURN(status=1)
      
      else
      
        ! dummy ...
        allocate( iglb_all(1), stat=status )
        IF_NOTOK_RETURN(status=1)

      end if ! root

      ! gather on root:
      call goc%GatherV( self%iglb, iglb_all, status, nloc=self%npix )
      IF_NOTOK_RETURN(status=1)

      !! IN FUTURE: collect overlap fractions as weights ...
      !! storage:
      !allocate( w_all(ntot), stat=status )
      !IF_NOTOK_RETURN(status=1)
      !! gather:
      !call goc%Gather( self%w, w_all, status )
      !IF_NOTOK_RETURN(status=1)
      
      ! collect distributed arrays on root and write from there:
      call self%hp%NcPutGather( ncid, self%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%yr%NcPutGather( ncid, self%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%sr%NcPutGather( ncid, self%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%K%NcPutGather( ncid, self%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      !~ including a priori profiles?
      if ( self%with_apriori ) then
        call self%ya%NcPutGather( ncid, self%nglb, iglb_all, status )
        IF_NOTOK_RETURN(status=1)
        call self%yra%NcPutGather( ncid, self%nglb, iglb_all, status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! testing ..
      call self%ilon%NcPutGather( ncid, self%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%ilat%NcPutGather( ncid, self%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)

      ! analysis flags:
      call self%astatus%NcPutGather( ncid, self%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      
      ! written on root...
      if ( goc%root ) then
        ! close:
        status = NF90_Close( ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if  ! root

      ! clear:
      deallocate( iglb_all, stat=status )
      IF_NOTOK_RETURN(status=1)     
      
    else
    
      ! info ..
      write (gol,'(a,":   no data for this time ...")') rname; call goPr
      
    end if

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_Data_PutOut


  ! ====================================================================
  ! ===
  ! === State
  ! ===
  ! ====================================================================


  subroutine LE_Output_Sat_State_Init( self, status, key, description )
  
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_State), intent(out)   ::  self
    integer, intent(out)                        ::  status
    character(len=*), intent(in), optional      ::  key
    character(len=*), intent(in), optional      ::  description

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Output_Sat_State_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! no data yet:
    self%nlevi  = 0
    self%nlev   = 0
    self%nlayer = 0
    self%nr     = 0
    self%npix   = 0
    ! no flags ...
    self%with_profile = .false.
    
    ! defaults:
    self%key = 'state'
    self%description = 'simulated satellite retrievals'
    ! replace?
    if ( present(key) ) then
      if ( len_trim(key) > 0 ) self%key = trim(key)
    end if
    if ( present(description) ) self%description = trim(description)
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_State_Init


  ! ***


  subroutine LE_Output_Sat_State_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_State), intent(inout)   ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_State_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%Clear( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_State_Done


  ! ***


  subroutine LE_Output_Sat_State_Clear( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_State), intent(inout)   ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_State_Clear'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! any local pixels?
    if ( self%npix > 0 ) then
      ! clear:
      call self%hp%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%vmr%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%vcd%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      call self%hx%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      call self%y%Done( status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! no pixels anymore:
    self%npix = 0

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_State_Clear


  ! ***


  subroutine LE_Output_Sat_State_GetPixel( self, ipix, status, &
                                            y )
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_State), intent(in)      ::  self
    integer, intent(in)                           ::  ipix
    integer, intent(out)                          ::  status
    
    real, intent(out), optional                   ::  y
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_State_GetPixel'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ipix < 1) .or. (ipix > self%npix) ) then
      write (gol,'("pixel ",i0," out of range 1:",i0)') ipix, self%npix; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! retrieval of total column: SHOULD BE REPLACED BY PROFILE!
    if ( present(y) ) then
      ! check ...
      if ( self%nr /= 1 ) then
        write (gol,'("profiles not supported yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! first value ...
      y = self%y%data(1,ipix)
    end if
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_State_GetPixel


  ! ***
  
  
  ! allocate storage for simulations

  subroutine LE_Output_Sat_State_Setup( self, sdat, c, status )
  
    use GO             , only : goMatchValue
    use Binas          , only : Avog, xm_air, grav
    use Num            , only : IntervalSum
    use Num            , only : Interval
    use Indices        , only : specname
    use Indices        , only : specunit
    use LE_Grid        , only : dom
    use LE_Data_Common , only : nlev, nlev_top
    use LE_Data        , only : LE_Data_GetPointer
    use LE_Bound_Common, only : caloft
  
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_State), intent(inout)   ::  self
    class(T_LE_Output_Sat_Data), intent(in)       ::  sdat
    real, intent(in)                              ::  c(:,:,:,:)  ! (nx,ny,nlev,nspec)
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Output_Sat_State_Setup'
    
    ! --- local ----------------------------------
    
    integer                ::  ipix
    integer                ::  iglb
    integer                ::  i, j
    integer                ::  ilev
    integer                ::  ispec
    real, pointer          ::  hp(:,:,:)   ! (lon,lat,0:nlev)
    real                   ::  dm
    integer                ::  ilayer
    integer                ::  ilast
    integer                ::  ir
    integer                ::  off(2)
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a,": simulate ",i0," / ",i0," pixels ...")') rname, sdat%npix, sdat%nglb; call goPr

    ! clear current storage:
    call self%Clear( status )
    IF_NOTOK_RETURN(status=1)
    
    ! any pixels in global domain?
    if ( sdat%nglb > 0 ) then

      ! source tracer:
      call goMatchValue( trim(sdat%tracer_name), specname, ispec, status )
      IF_NOTOK_RETURN(status=1)

      ! copy flags:
      self%with_profile = sdat%with_profile
      ! local number of pixels:
      self%npix = sdat%npix
      ! model levels:
      self%nlevi  = nlev_top+1
      self%nlev   = nlev_top
      ! retrieval input layers:
      self%nlayer = sdat%nlayer
      self%nr     = sdat%nr

    ! init pixel arrays, if npix==0 no array is allocated
    ! allocate storage, filled with fill_value ..
    call self%hp%Init( (/nlev_top+1,self%npix/), status, &
                           long_name='model half level pressures', &
                           units=sdat%hp%units )
    IF_NOTOK_RETURN(status=1)
    ! allocate storage, filled with fill_value ..
    call self%vmr%Init( (/nlev_top,self%npix/), status, &
                           long_name='model volume mixing ratio', &
                           units=specunit(ispec) )
    IF_NOTOK_RETURN(status=1)
    ! allocate storage, filled with fill_value ..
    call self%vcd%Init( (/nlev_top,self%npix/), status, &
                           long_name='model vertical column density', &
                           units=sdat%yr%units )
    IF_NOTOK_RETURN(status=1)
    ! allocate storage, filled with fill_value ..
    call self%hx%Init( (/sdat%nlayer,self%npix/), status, &
                           long_name='model vertical column density at retrieval layers', &
                           units=sdat%yr%units )
    IF_NOTOK_RETURN(status=1)
    ! allocate storage, filled with fill_value ..
    call self%y%Init( (/sdat%nr,self%npix/), status, &
                           long_name='simulated retrieval', &
                           units=sdat%yr%units )
    IF_NOTOK_RETURN(status=1)

    ! any local pixels?
    if ( self%npix > 0 ) then

        ! pointer to meteo data:
        call LE_Data_GetPointer( 'hp', hp, status, check_units=sdat%hp%units )    
        IF_NOTOK_RETURN(status=1)
        
        ! global offset for grid indices:
        call dom%Get( status, off=off )
        IF_NOTOK_RETURN(status=1)

        ! loop over local pixels:
        do ipix = 1, self%npix
          ! global index:
          iglb = sdat%iglb(ipix)
          
          ! local grid cell:
          i = sdat%ilon%data(ipix) - off(1)
          j = sdat%ilat%data(ipix) - off(2)
          
          ! half level pressures:
          self%hp%data(:,ipix) = hp(i,j,0:nlev_top)
          
          ! fill full model profile (volume mixing ratios):
          self%vmr%data(     1:nlev    ,ipix) = c     (i,j,     1:nlev    ,ispec)
          self%vmr%data(nlev+1:nlev_top,ipix) = caloft(i,j,nlev+1:nlev_top,ispec)
          
          ! check units ...
          if ( (specunit(ispec) == 'ppb') .and. (sdat%yr%units == '1e15 mlc/cm2') ) then
            !
            ! Compute partial vertical column densities in (1e15 mlc trc)/cm2 :
            ! 
            !  1e15 mlc trc    mol trc    mol air  kg air  m2
            !  ---- -------  -----------  -------  ------  ---  / 1e15  = (1e15 mlc trc)/cm2
            !   1   mol trc  1e9 mol air  kg air     m2    cm2
            !
            ! loop over layers:
            do ilev = 1, nlev_top
              ! air mass density:
              dm = abs( self%hp%data(ilev,ipix) - self%hp%data(ilev+1,ipix) )/grav  ! (kg air)/m2
              ! convert:
              !                                 mol tracer ppb         1     mlc   mol air  kg air    m2
              !                                ---------------        ---    ---   -------  ------    ---
              !                                  mol air              ppb    mol    kg air    m2      cm2
              self%vcd%data(ilev,ipix) = self%vmr%data(ilev,ipix) * 1.0e-9 * Avog / xm_air *  dm   * 1.0e-4 / 1.0e15 ! (1e15 mlc trc)/cm2
            end do ! levels
          end if
          
          ! adhoc: lower bottom if pixel surface is below model ...
          ! note that total tracer mass is not changed since vcd has been computed already:
          self%hp%data(1,ipix) = max( self%hp%data(1,ipix), sdat%hp%data(1,ipix) )
          
          ! loop over retrieval layers:
          do ilayer = 1, sdat%nlayer
            ! fill each retrieval layer with sum of one or more fractions of model layers;
            ! use pressure axis to have mass-conservation; negate to have increasing axis:
            call IntervalSum( -1.0*self%hp%data(:,ipix), self%vcd%data(:,ipix), &
                               -1.0*sdat%hp%data(ilayer,ipix), -1.0*sdat%hp%data(ilayer+1,ipix), &
                               self%hx%data(ilayer,ipix), ilast, status )
            if ( status /= 0 ) then
              write (gol,'("mapping from model to retrieval layer:")'); call goErr
              write (gol,'("model layers (ph0,ph1,vmr,vcd):")'); call goErr
              do ilev = nlev_top, 1, -1
                write (gol,'(i4," ",2f12.4," ",2es12.4)') ilev, self%hp%data(ilev:ilev+1,ipix), &
                                self%vmr%data(ilev,ipix), self%vcd%data(ilev,ipix); call goErr
              end do
              write (gol,'("pixel: ",i6)') ipix; call goErr
              write (gol,'("retrieval layer:")'); call goErr
              write (gol,'(i4," ",2f12.4)') ilayer, sdat%hp%data(ilayer:ilayer+1,ipix); call goErr
              TRACEBACK; status=1; return
            end if
          end do ! retr layers
          
          ! apriori profile present?
          if ( sdat%with_apriori ) then
            write (gol,'("no apriori supported yet")'); call goErr
            TRACEBACK; status=1; return
          else
            ! loop over retrieval layers:
            do ir = 1, sdat%nr
              ! weighted sum:
              self%y%data(ir,ipix) = sum( sdat%K%data(ir,:,ipix) * self%hx%data(:,ipix) )
            end do ! retr layer
          end if
          
        end do  ! pixels
        
    end if  ! any local pixels
    
    end if  ! any global pixels

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_State_Setup


  ! ***


  subroutine LE_Output_Sat_State_PutOut( self, sdat, t, status )

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim
    use NetCDF , only : NF90_EndDef
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER
  
    use GO              , only : TDate, wrtgol
    use GO              , only : goc
    use LE_Output_Common, only : PutOut_GlobalAttributes
    use LE_Data_Common  , only : nlev_top
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_State), intent(inout)   ::  self
    class(T_LE_Output_Sat_Data), intent(in)       ::  sdat
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_State_PutOut'
    
    ! --- local ----------------------------------
    
    character(len=1024)       ::  fname
    integer                   ::  cmode
    integer                   ::  ncid
    integer                   ::  dimid_pixel
    integer                   ::  dimid_corner
    integer                   ::  dimid_mlay
    integer                   ::  dimid_mlayi
    integer                   ::  dimid_layer
    integer                   ::  dimid_retr
    integer                   ::  ntot
    integer, allocatable      ::  iglb_all(:)
    
    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': put out for ', t ); call goPr
    
    ! any data?
    if ( sdat%nglb > 0 ) then
    
      ! total number of pixels handled by local domains;
      ! this will be >= nglb since some footprints cover multiple domains:
      call goc%ParInfo( self%npix, status, ntot=ntot )
      IF_NOTOK_RETURN(status=1)

      ! collect on root, this is much faster than parallel write ...
      if ( goc%root ) then

        ! target file:
        write (fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",2i2.2,"_",a,".nc")') &
             trim(sdat%com%outdir), trim(sdat%com%model), trim(sdat%com%expid), &
             trim(sdat%name), t%year, t%month, t%day, t%hour, t%min, &
             trim(self%key)

        ! info ..
        write (gol,'(a,":   create ",a," ...")') rname, trim(fname); call goPr

        ! set creation mode flag:
        if ( sdat%com%replace ) then
          cmode = NF90_CLOBBER       ! overwrite existing files
        else
          cmode = NF90_NOCLOBBER     ! do not overwrite existing files
        end if

        ! create file:
        status = NF90_Create( fname, cmode, ncid )
        if ( status /= 0 ) then
           write (gol,'("creating file :")'); call goErr
           write (gol,'("  ",a)') trim(fname); call goErr
           TRACEBACK; status=1; return
        end if

        ! write global attributes:
        call PutOut_GlobalAttributes( sdat%com, ncid, status )
        IF_NOTOK_RETURN(status=1)

        ! define dimensions:
        status = NF90_Def_Dim( ncid, 'pixel', sdat%nglb, dimid_pixel )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'mlay', nlev_top, dimid_mlay )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'mlayi', nlev_top+1, dimid_mlayi )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'layer', sdat%nlayer, dimid_layer )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'retr', sdat%nr, dimid_retr )
        IF_NF90_NOTOK_RETURN(status=1)

        ! define variables:
        !~ distributed arrays, will be collected on write:
        call self%hp%NcDef( ncid, 'hp' , (/dimid_mlayi,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%vmr%NcDef( ncid, 'vmr' , (/dimid_mlay,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%vcd%NcDef( ncid, 'vcd' , (/dimid_mlay,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%hx%NcDef( ncid, 'hx' , (/dimid_layer,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%y%NcDef( ncid, 'y' , (/dimid_retr,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)

        ! end defintion mode:
        status = NF90_EndDef( ncid )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! storage for mapping:
        allocate( iglb_all(ntot), stat=status )
        IF_NOTOK_RETURN(status=1)
      
      else
      
        ! dummy ...
        allocate( iglb_all(1), stat=status )
        IF_NOTOK_RETURN(status=1)

      end if ! root

      ! gather on root:
      call goc%GatherV( sdat%iglb, iglb_all, status, nloc=self%npix )
      IF_NOTOK_RETURN(status=1)

      ! collect distributed arrays on root and write from there:
      call self%hp%NcPutGather( ncid, sdat%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%vmr%NcPutGather( ncid, sdat%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%vcd%NcPutGather( ncid, sdat%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%hx%NcPutGather( ncid, sdat%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%y%NcPutGather( ncid, sdat%nglb, iglb_all, status )
      IF_NOTOK_RETURN(status=1)
      
      ! written on root...
      if ( goc%root ) then
        ! close:
        status = NF90_Close( ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if  ! root

      ! clear:
      deallocate( iglb_all, stat=status )
      IF_NOTOK_RETURN(status=1)     
      
    else
    
      ! info ..
      write (gol,'(a,":   no data for this time ...")') rname; call goPr
      
    end if

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_State_PutOut





end module LE_Output_Sat

