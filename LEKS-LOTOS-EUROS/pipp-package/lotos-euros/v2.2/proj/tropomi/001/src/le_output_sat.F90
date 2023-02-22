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
    ! filename template(s):
    character(len=1024)               ::  filenames(10)
    integer                           ::  nfilename
    ! variable names:
    character(len=32)                 ::  varname_hp
    character(len=32)                 ::  varname_ya
    character(len=32)                 ::  varname_yra
    character(len=32)                 ::  varname_yr
    character(len=32)                 ::  varname_air_vcd
    character(len=32)                 ::  varname_sr
    character(len=32)                 ::  varname_vr
    ! flags:
    logical                           ::  with_corners
    logical                           ::  with_profile
    logical                           ::  with_apriori
    logical                           ::  with_track
    logical                           ::  with_air_vcd
    !
    !~ global pixel arrays,
    !  could be full track including pixels outside model domain
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
    ! on root only, collecting all "iglb" arrays:
    integer                           ::  npix_all
    integer, allocatable              ::  iglb_all(:)   ! (npix_all)
    ! number of pixels used at some domain:
    integer                           ::  nout
    ! mapping from 1:nglb to 1:nout:
    integer, allocatable              ::  iout_glb(:)   ! (nglb)
    ! mapping from 1:npix_all to 1:nout:
    integer, allocatable              ::  iout_all(:)   ! (npix_all)
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
    integer                           ::  nri  ! nr+1 for interfaces
    ! half-level pressures:
    type(T_Pixels_1D)                 ::  hp   ! (nlayeri,npix)
    type(T_Pixels_1D)                 ::  hpr  ! (nri,npix)
    ! apriori profiles:
    type(T_Pixels_1D)                 ::  ya   ! (nlayer,npix)
    ! averaging kernels:
    type(T_Pixels_2D)                 ::  K    ! (nr,nlayer,npix)
    ! retrieval apriori (profiles):
    type(T_Pixels_1D)                 ::  yra   ! (nr,npix)
    ! retrieval (profiles):
    type(T_Pixels_1D)                 ::  yr    ! (nr,npix)
    type(T_Pixels_1D)                 ::  yr_vcd    ! (nr,npix)
    ! air column (profiles):
    type(T_Pixels_1D)                 ::  air_vcd    ! (nr,npix)
    ! retrieval error std.dev. (no profiles)
    type(T_Pixels_1D)                 ::  sr    ! (nr,npix)
    ! retrieval error covariance (profiles)
    type(T_Pixels_2D)                 ::  vr    ! (nr,nr,npix)
    !
    ! index range 1,..,nlayer with valid values:
    type(T_Pixels_0D_i)               ::  ilayer1    ! (npix)
    type(T_Pixels_0D_i)               ::  ilayer2    ! (npix)
    type(T_Pixels_0D_i)               ::  ir1        ! (npix)
    type(T_Pixels_0D_i)               ::  ir2        ! (npix)
    !
    ! how to fill 'sr' or 'vr' ?
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
    procedure :: GetRandom       => LE_Output_Sat_Data_GetRandom
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
    type(T_Pixels_1D)                 ::  mod_hp    ! (nlevi,npix)
    type(T_Pixels_1D)                 ::  mod_hp0   ! (nlevi,npix)
    type(T_Pixels_1D)                 ::  mod_vmr   ! (nlev,npix)
    type(T_Pixels_1D)                 ::  mod_vcd   ! (nlev,npix)
    ! simulated input profiles:
    type(T_Pixels_1D)                 ::  hx   ! (nlayer,npix)
    ! simulated retrieved profiles:
    type(T_Pixels_1D)                 ::  y    ! (nr,npix)
    ! in case of profiles, also vcd for convenience:
    type(T_Pixels_1D)                 ::  y_vcd ! (nr,npix)
    !! TESTING ..
    !type(T_Pixels_1D)                 ::  test_dlogy           ! (nlayer,npix)
    !type(T_Pixels_1D)                 ::  test_Kdlogy          ! (nr,npix)
    !type(T_Pixels_1D)                 ::  test_logya_Kdlogy    ! (nr,npix)
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
    use GO              , only : goSplitString
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
    
    character(len=64)     ::  rckey
    character(len=2000)   ::  line
    
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
    call rcF%Get( trim(rckey)//'.with_track'  , self%with_track, status )
    IF_NOTOK_RETURN(status=1)

    ! filename template(s):
    call rcF%Get( trim(rckey)//'.filename', line, status )
    IF_NOTOK_RETURN(status=1)
    ! split:
    call goSplitString( line, self%nfilename, self%filenames, status )
    IF_NOTOK_RETURN(status=1)
    
    ! variable names:
    call rcF%Get( trim(rckey)//'.varname.hp', self%varname_hp, status )
    IF_NOTOK_RETURN(status=1)
    ! apri?
    if ( self%with_apriori ) then
      call rcF%Get( trim(rckey)//'.varname.yra', self%varname_yra, status )
      IF_NOTOK_RETURN(status=1)
      call rcF%Get( trim(rckey)//'.varname.ya' , self%varname_ya , status )
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

    ! read variable with air column?
    call rcF%Get( trim(rckey)//'.with_air_vcd'  , self%with_air_vcd, status, default=.false. )
    IF_ERROR_RETURN(status=1)
    ! enabled?
    if ( self%with_air_vcd ) then
      ! variable name:
      call rcF%Get( trim(rckey)//'.varname.air_vcd', self%varname_air_vcd, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! no pixels yet:
    self%nglb     = 0
    self%npix     = 0
    self%npix_all = 0
    self%nout     = 0

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
    
      ! track defined?
      if ( self%with_track ) then
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
        end if ! corners
      end if ! track

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
      deallocate( self%iglb_all, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%iout_glb, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%iout_all, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! no pixels anymore:
      self%nglb     = 0
      self%npix     = 0
      self%npix_all = 0
      self%nout     = 0
    
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
      call self%hpr%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! apri?
      if ( self%with_apriori ) then
        call self%ya%Done( status )
        IF_NOTOK_RETURN(status=1)
        call self%yra%Done( status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! air column?
      if ( self%with_air_vcd ) then
        call self%air_vcd%Done( status )
        IF_NOTOK_RETURN(status=1)
      end if
      call self%K%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%yr%Done( status )
      IF_NOTOK_RETURN(status=1)
      if ( self%with_profile ) then
        ! valid layer range (input profiles):
        call self%ilayer1%Done( status )
        IF_NOTOK_RETURN(status=1)
        call self%ilayer2%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! valid retrieved layer range:
        call self%ir1%Done( status )
        IF_NOTOK_RETURN(status=1)
        call self%ir2%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! column density:
        call self%yr_vcd%Done( status )
        IF_NOTOK_RETURN(status=1)
        ! covariance:
        call self%vr%Done( status )
        IF_NOTOK_RETURN(status=1)
      else
        call self%sr%Done( status )
        IF_NOTOK_RETURN(status=1)
      end if

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

    use GO     , only : goc
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
    
    integer                 ::  ifilename
    character(len=1024)     ::  fname
    logical                 ::  exist
    integer                 ::  ncid
    integer                 ::  dimid
    integer                 ::  varid
    integer, allocatable    ::  iglb_local(:)  ! (nglb)
    logical, allocatable    ::  used_glb(:)    ! (nglb)
    integer                 ::  ipix
    integer                 ::  iglb
    logical                 ::  indomain
    integer                 ::  off(2)
    integer                 ::  ilayer
    integer                 ::  ir
    character(len=64)       ::  conversion
    real                    ::  dm
    integer                 ::  i, j
    integer                 ::  k
    real, pointer           ::  hp(:,:,:)   ! (lon,lat,0:nlev)
    
    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': setup for ', tr ); call goPr
    
    ! clear current storage:
    call self%Clear( status )
    IF_NOTOK_RETURN(status=1)

    ! loop over filename templates:
    do ifilename = 1, self%nfilename
      ! start filename with template:
      fname = trim(self%filenames(ifilename))
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
      ! leave at first match:
      if ( exist ) exit
      ! testing ...
      write (gol,'(a,": WARNING - no file: ",a)') rname, trim(fname); call goPr
    end do  ! filenames
    
    ! found a file ?
    if ( exist ) then

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
      
      ! track defined?
      if ( self%with_track ) then
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
      end if
      
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
      ! interfaces:
      self%nri = self%nr + 1

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

      ! track defined?
      if ( self%with_track ) then

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

      end if ! track

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
      allocate( iglb_local(self%nglb), source=-999, stat=status )
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
          iglb_local(self%npix) = ipix
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
        self%iglb = iglb_local(1:self%npix)
        
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

      ! clear:
      deallocate( iglb_local, stat=status )
      IF_NOTOK_RETURN(status=1)

      ! total number of pixels handled by local domains, broadcasted to all;
      ! this might be >= nglb since some footprints cover multiple domains;
      call goc%ParInfo( self%npix, status, ntot=self%npix_all )
      IF_NOTOK_RETURN(status=1)
      ! storage for global indices for all local pixels, needed on root only:
      if ( goc%root ) then
        ! storage for mapping:
        allocate( self%iglb_all(self%npix_all), stat=status )
        IF_NOTOK_RETURN(status=1)
      else
        ! dummy ...
        allocate( self%iglb_all(1), stat=status )
        IF_NOTOK_RETURN(status=1)
      end if ! root
      ! gather on root:
      call goc%GatherV( self%iglb, self%iglb_all, status, nloc=self%npix )
      IF_NOTOK_RETURN(status=1)
      
      ! on root, info on selected pixels
      if ( goc%root ) then

        ! storage for flag per global pixel, disable by default:
        allocate( used_glb(self%nglb), source=.false., stat=status )
        IF_NOTOK_RETURN(status=1)
        ! loop over all locally used pixesl:
        do k = 1, self%npix_all
          ! global pixel index:
          iglb = self%iglb_all(k)
          ! enable:
          used_glb(iglb) = .true.
        end do

        ! count number of global pixels somewhere used:
        self%nout = count( used_glb )
        ! storage for mapping from 1:nglb to 1:nout
        allocate( self%iout_glb(self%nglb), source=-999, stat=status )
        IF_NOTOK_RETURN(status=1)
        ! init counter:
        ipix = 0
        ! loop over global pixels:
        do iglb = 1, self%nglb
          ! in use?
          if ( used_glb(iglb) ) then
            ! next:
            ipix = ipix + 1
            ! store mapping:
            self%iout_glb(iglb) = ipix
          end if
        end do ! iglb
        
        ! clear:
        deallocate( used_glb, stat=status )
        IF_NOTOK_RETURN(status=1)
        
        ! storage for mapping from "all" array to "out":
        allocate( self%iout_all(self%npix_all), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! loop over all locally used pixesl:
        do k = 1, self%npix_all
          ! global pixel index:
          iglb = self%iglb_all(k)
          ! mapping:
          self%iout_all(k) = self%iout_glb(iglb)
        end do ! k

      else
        ! dummy ..
        self%npix_all = -999
        ! dummy ...
        allocate( self%iout_glb(1), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! dummy ...
        allocate( self%iout_all(1), stat=status )
        IF_NOTOK_RETURN(status=1)
      end if ! root
      
      ! broadcast counter:
      call goc%BCast( goc%root_id, self%nout, status )
      IF_NOTOK_RETURN(status=1)
      
      ! ...
        
      ! storage for status flags, default to 0 (analyse?)
      call self%astatus%Init( self%npix, status, source=0, units='1', long_name='assimilation status' )
      IF_NOTOK_RETURN(status=1)

      ! read half-level pressures (local subset):
      call self%hp%NcInit( ncid, trim(self%varname_hp), status, &
                             nselect=self%npix, select=self%iglb )
      IF_NOTOK_RETURN(status=1)

      ! index range for input layers, here all valid
      ! IN FUTURE, USE FILL_VALUES FROM NETCDF FILES!
      !~ start index:
      call self%ilayer1%Init( self%npix, status )
      IF_NOTOK_RETURN(status=1)
      ! .. find first 
      do ipix = 1, self%npix
        do ilayer = 1, self%nlayer
          if ( self%hp%data(ilayer,ipix) >= 0.0 ) then
            self%ilayer1%data(ipix) = ilayer
            exit
          end if
        end do
      end do
      !~ end index:
      call self%ilayer2%Init( self%npix, status, source=self%nlayer )
      IF_NOTOK_RETURN(status=1)

      ! index range for retrieval layers:
      !~ start index:
      call self%ir1%Init( self%npix, status, source=1 )
      IF_NOTOK_RETURN(status=1)
      ! .. for profiles, reset to index of first valid input layer:
      if ( self%with_profile ) then
        ! ... reset:
        self%ir1%data = self%ilayer1%data
      end if
      !~ end index:
      call self%ir2%Init( self%npix, status, source=self%nr )
      IF_NOTOK_RETURN(status=1)

      ! a-priori profiles in product?
      if ( self%with_apriori ) then
        ! read 'true' apriori profiles (local subset):
        call self%ya%NcInit( ncid, trim(self%varname_ya), status, &
                               nselect=self%npix, select=self%iglb )
        IF_NOTOK_RETURN(status=1)
        ! read 'retrieval' apriori profiles (local subset):
        call self%yra%NcInit( ncid, trim(self%varname_yra), status, &
                               nselect=self%npix, select=self%iglb )
        IF_NOTOK_RETURN(status=1)
      end if

      ! air column?
      if ( self%with_air_vcd ) then
        ! read air column at retreival layers (local subset):
        call self%air_vcd%NcInit( ncid, trim(self%varname_air_vcd), status, &
                               nselect=self%npix, select=self%iglb )
        IF_NOTOK_RETURN(status=1)
      end if

      ! read kernels (local subset):
      call self%K%NcInit( ncid, 'kernel', status, &
                            nselect=self%npix, select=self%iglb )
      IF_NOTOK_RETURN(status=1)

      ! read retrieval profiles (local subset):
      call self%yr%NcInit( ncid, trim(self%varname_yr), status, &
                             nselect=self%npix, select=self%iglb )
      IF_NOTOK_RETURN(status=1)
      ! conversion to vcd if necessary:
      if ( self%with_profile ) then
        call self%yr_vcd%Init( (/self%nr,self%npix/), status, &
                                 units=vcd_units, long_name='vertical column density of retrieved profile' )
        IF_NOTOK_RETURN(status=1)
      end if

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

      ! no data ..
      do ipix = 1, self%npix
        ! loop over input layers:
        do ilayer = 1, self%nlayer
          ! outside valid range ?
          if ( (ilayer < self%ilayer1%data(ipix)) .or. (ilayer > self%ilayer2%data(ipix)) ) then
            ! reset pressure (half) levels:
            if ( ilayer < self%ilayer1%data(ipix) ) then
              self%hp%data(ilayer,ipix) = self%hp%fill_value
            else if ( ilayer > self%ilayer2%data(ipix) ) then
              self%hp%data(ilayer+1,ipix) = self%hp%fill_value
            end if
            ! a-priori profiles in product?
            if ( self%with_apriori ) then
              self%ya%data(ilayer,ipix) = self%ya%fill_value
            end if
            ! reset:
            self%K%data(:,ilayer,ipix) = self%K%fill_value
          end if
        end do ! ilayer
        ! loop over retrieval layers:
        do ir = 1, self%nr
          ! outside valid range?
          if ( (ir < self%ir1%data(ipix)) .or. (ir > self%ir2%data(ipix)) ) then
            ! a-priori profiles in product?
            if ( self%with_apriori ) then
              self%yra%data(ir,ipix) = self%yra%fill_value
            end if
            ! air column?
            if ( self%with_air_vcd ) then
              self%air_vcd%data(ir,ipix) = self%air_vcd%fill_value
            end if
            ! reset:
            self%K%data(ir,:,ipix) = self%K%fill_value
            self%yr%data(ir,ipix) = self%yr%fill_value
            ! single value or profile?
            if ( self%with_profile ) then
              self%vr%data(ir,:,ipix) = self%vr%fill_value
              self%vr%data(:,ir,ipix) = self%vr%fill_value
            else
              self%sr%data(ir,ipix) = self%sr%fill_value
            end if
          end if ! ir outside range
        end do ! ir
      end do ! ipix
      
      ! reset surface pressure to model surface;
      ! pointer to meteo data:
      call LE_Data_GetPointer( 'hp', hp, status, check_units=self%hp%units )    
      IF_NOTOK_RETURN(status=1)
      ! loop over pixels:
      do ipix = 1, self%npix
        ! local grid cell:
        i = self%ilon%data(ipix) - off(1)
        j = self%ilat%data(ipix) - off(2)
        ! lowset layer in input profile:
        ilayer = self%ilayer1%data(ipix)
        ! replace lowest value if below model ...
        self%hp%data(ilayer,ipix) = max( self%hp%data(ilayer,ipix), hp(i,j,0) )
      end do ! ipix

      ! half level pressures for retrieved profile:
      call self%hpr%Init( (/self%nri,self%npix/), status, &
                            units=self%hp%units, long_name='pressure interfaces for retrieved profile' )
      IF_NOTOK_RETURN(status=1)
      ! switch:
      if ( self%with_profile ) then
        ! copy:
        self%hpr%data = self%hp%data
      else
        ! loop over pixels:
        do ipix = 1, self%npix
          ! copy bottom from first input layer:
          self%hpr%data(1,ipix) = self%hp%data(self%ilayer1%data(ipix)  ,ipix)
          ! copy top from last input layer:
          self%hpr%data(2,ipix) = self%hp%data(self%ilayer2%data(ipix)+1,ipix)
        end do ! pixels
      end if ! profile or column
      
      ! compute vcd from retrieval?
      if ( self%with_profile ) then
        ! conversion from retrieval units to vcd units:
        conversion = trim(self%yr%units)//' -> '//trim(vcd_units)
        ! switch:
        select case ( trim(conversion) )
          !~ volume mixing ratio:
          case ( 'ppmv -> 1e15 mlc/cm2' )
            ! loop over pixels:
            do ipix = 1, self%npix
              ! loop over retrieval layers:
              do ir = self%ir1%data(ipix), self%ir2%data(ipix)
                ! air mass density:
                dm = abs( self%hpr%data(ir,ipix) - self%hpr%data(ir+1,ipix) )/grav  ! (kg air)/m2
                ! convert:
                !                              mol tracer ppm          1     mlc   mol air  kg air    m2
                !                              --------------         ---    ---   -------  ------    ---
                !                                 mol air             ppm    mol    kg air    m2      cm2
                self%yr_vcd%data(ir,ipix) = self%yr%data(ir,ipix) * 1.0e-6 * Avog / xm_air *  dm   * 1.0e-4 / 1.0e15 ! (1e15 mlc trc)/cm2
              end do ! levels
            end do  ! pixels
          !~ unknown ...
          case default
            write (gol,'("unsupported conversion `",a,"`")') trim(conversion); call goErr
            TRACEBACK; status=1; return
        end select
      end if ! with profile

      ! close:
      status = NF90_Close( ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      
    else

      ! no data:
      self%nglb = 0
      
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
    integer                   ::  dimid_retri
    
    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': put out for ', t ); call goPr
    
    ! any pixels to be put out?
    if ( self%nout > 0 ) then
    
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
        status = NF90_Def_Dim( ncid, 'pixel', self%nout, dimid_pixel )
        IF_NF90_NOTOK_RETURN(status=1)
        if ( self%with_corners ) then
          status = NF90_Def_Dim( ncid, 'corner', self%ncorner, dimid_corner )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
        ! track defined?
        if ( self%with_track ) then
          status = NF90_Def_Dim( ncid, 'track_pixel', self%ntx, dimid_tx )
          IF_NF90_NOTOK_RETURN(status=1)
          status = NF90_Def_Dim( ncid, 'track_image', self%nty, dimid_ty )
          IF_NF90_NOTOK_RETURN(status=1)
        end if ! track
        status = NF90_Def_Dim( ncid, 'layer', self%nlayer, dimid_layer )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'layer_interface', self%nlayeri, dimid_layeri )
        IF_NF90_NOTOK_RETURN(status=1)
        ! retrieval layers:
        status = NF90_Def_Dim( ncid, 'retr', self%nr, dimid_retr )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( ncid, 'retr_interface', self%nri, dimid_retri )
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
        ! track defined?
        if ( self%with_track ) then
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
        end if ! track
        !~ distributed arrays, will be collected on write:
        call self%hp%NcDef( ncid, 'hp' , (/dimid_layeri,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%hpr%NcDef( ncid, 'hpr' , (/dimid_retri,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%yr%NcDef( ncid, 'yr' , (/dimid_retr,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        if ( self%with_profile ) then
          call self%yr_vcd%NcDef( ncid, 'yr_vcd', (/dimid_retr,dimid_pixel/), status )
          IF_NOTOK_RETURN(status=1)
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
        !~ including air columns?
        if ( self%with_apriori ) then
          call self%air_vcd%NcDef( ncid, 'air_vcd', (/dimid_retr,dimid_pixel/), status )
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

      end if ! root

      ! write global arrays from root:
      if ( goc%root ) then

        ! write arrays available as global arrays:
        call self%glb_lon%NcPutGlbSelect( ncid, self%iout_glb, status )
        IF_NOTOK_RETURN(status=1)
        call self%glb_lat%NcPutGlbSelect( ncid, self%iout_glb, status )
        IF_NOTOK_RETURN(status=1)
        ! corners defined?
        if ( self%with_corners ) then
          call self%glb_clons%NcPutGlbSelect( ncid, self%iout_glb, status )
          IF_NOTOK_RETURN(status=1)
          call self%glb_clats%NcPutGlbSelect( ncid, self%iout_glb, status )
          IF_NOTOK_RETURN(status=1)
        end if
        
        ! track defined?
        if ( self%with_track ) then
          ! write track arrays available as global arrays:
          call self%glb_tpixel%NcPutGlbSelect( ncid, self%iout_glb, status )
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
          end if  ! corners
        end if  ! track

      end if ! root
      
      !! IN FUTURE: collect overlap fractions as weights ...
      !! storage:
      !allocate( w_all(ntot), stat=status )
      !IF_NOTOK_RETURN(status=1)
      !! gather:
      !call goc%Gather( self%w, w_all, status )
      !IF_NOTOK_RETURN(status=1)
      
      ! collect distributed arrays on root and write from there:
      call self%hp%NcPutGather( ncid, self%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%hpr%NcPutGather( ncid, self%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%yr%NcPutGather( ncid, self%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      if ( self%with_profile ) then
        call self%yr_vcd%NcPutGather( ncid, self%iout_all, status )
        IF_NOTOK_RETURN(status=1)
        call self%vr%NcPutGather( ncid, self%iout_all, status )
        IF_NOTOK_RETURN(status=1)
      else
        call self%sr%NcPutGather( ncid, self%iout_all, status )
        IF_NOTOK_RETURN(status=1)
      end if
      call self%K%NcPutGather( ncid, self%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      !~ including a priori profiles?
      if ( self%with_apriori ) then
        call self%ya%NcPutGather( ncid, self%iout_all, status )
        IF_NOTOK_RETURN(status=1)
        call self%yra%NcPutGather( ncid, self%iout_all, status )
        IF_NOTOK_RETURN(status=1)
      end if
      !~ including air columns?
      if ( self%with_air_vcd ) then
        call self%air_vcd%NcPutGather( ncid, self%iout_all, status )
        IF_NOTOK_RETURN(status=1)
      end if
      ! testing ..
      call self%ilon%NcPutGather( ncid, self%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%ilat%NcPutGather( ncid, self%iout_all, status )
      IF_NOTOK_RETURN(status=1)

      ! analysis flags:
      call self%astatus%NcPutGather( ncid, self%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      
      ! written on root...
      if ( goc%root ) then
        ! close:
        status = NF90_Close( ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if  ! root

    else
    
      ! info ..
      write (gol,'(a,":   no data for this time ...")') rname; call goPr
      
    end if

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_Data_PutOut


  ! ***


  !
  ! Fill array with local pixel size with random numbers out of N(0,1).
  ! An array with length 'nglb' will be filled with random numbers first;
  ! this array will be the same on all domains if the provided random generator
  ! was initialized with the same seed on all domains.
  ! The local pixel selection is then extracted from this array.
  !

  subroutine LE_Output_Sat_Data_GetRandom( self, rnd, v, status )
  
    use Num, only : T_Random
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_Sat_Data), intent(inout)    ::  self
    type(T_Random), intent(inout)                 ::  rnd
    real, intent(out)                             ::  v(:)  ! (npix)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_Sat_Data_GetRandom'
    
    ! --- local ----------------------------------
    
    real, allocatable     ::  vglb(:)  ! (nglb)
    integer               ::  iglb
    integer               ::  ipix
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(v) < self%npix ) then
      write (gol,'("size v is ",i0," while nix is ",i0)') size(v), ipix; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( vglb(self%nglb), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over global number of pixels:
    do iglb = 1, self%nglb
      ! generate random number out of N(0,1) distribution:
      call rnd%Get_Normal( vglb(iglb), status )
      IF_NOTOK_RETURN(status=1)
    end do ! i
    
    ! loop over local pixels:
    do ipix = 1, self%npix
      ! current global index:
      iglb = self%iglb(ipix)
      ! copy:
      v(ipix) = vglb(iglb)
    end do ! ipix
    
    ! clear:
    deallocate( vglb, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine LE_Output_Sat_Data_GetRandom


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
      call self%mod_hp%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%mod_hp0%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%mod_vmr%Done( status )
      IF_NOTOK_RETURN(status=1)
      call self%mod_vcd%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      call self%hx%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      call self%y%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear?
      if ( self%with_profile ) then
        call self%y_vcd%Done( status )
        IF_NOTOK_RETURN(status=1)
        !! testing ..
        !call self%test_dlogy%Done( status )
        !IF_NOTOK_RETURN(status=1)
        !call self%test_Kdlogy%Done( status )
        !IF_NOTOK_RETURN(status=1)
        !call self%test_logya_Kdlogy%Done( status )
        !IF_NOTOK_RETURN(status=1)
      end if
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
    integer                ::  i, j, k
    integer                ::  ilev
    integer                ::  ispec
    character(len=64)      ::  conversion
    character(len=64)      ::  p_units
    character(len=64)      ::  khx_units
    real, pointer          ::  hp(:,:,:)   ! (lon,lat,0:nlev)
    real                   ::  dm
    real                   ::  hx_min
    integer                ::  ilayer, k1, k2
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
      ! retrieved levels
      self%nr     = sdat%nr

      ! pressure units in input ...
      p_units = trim(sdat%hp%units)
      ! pointer to meteo data:
      call LE_Data_GetPointer( 'hp', hp, status, check_units=p_units )    
      IF_NOTOK_RETURN(status=1)

      ! testing
      call self%mod_hp0%Init( (/nlev_top+1,self%npix/), status, &
                           long_name='model half level pressures', &
                                 units=p_units )
      IF_NOTOK_RETURN(status=1)

      ! init pixel arrays, if npix==0 no array is allocated
      ! allocate storage, filled with fill_value ..
      call self%mod_hp%Init( (/nlev_top+1,self%npix/), status, &
                           long_name='model half level pressures', &
                                 units=p_units )
      IF_NOTOK_RETURN(status=1)
      ! allocate storage, filled with fill_value ..
      call self%mod_vmr%Init( (/nlev_top,self%npix/), status, &
                           long_name='model volume mixing ratio', &
                           units=specunit(ispec) )
      IF_NOTOK_RETURN(status=1)
      ! allocate storage, filled with fill_value ..
      call self%mod_vcd%Init( (/nlev_top,self%npix/), status, &
                                 long_name='model vertical column density', &
                                 units=vcd_units )
      IF_NOTOK_RETURN(status=1)
      ! allocate storage, filled with fill_value ..
      if ( sdat%with_apriori ) then
        call self%hx%Init( (/sdat%nlayer,self%npix/), status, &
                               long_name='model vertical column density at apriori layers', &
                               units=sdat%ya%units )
        IF_NOTOK_RETURN(status=1)
      else
        call self%hx%Init( (/sdat%nlayer,self%npix/), status, &
                               long_name='model vertical column density at retrieval layers', &
                               units=sdat%yr%units )
        IF_NOTOK_RETURN(status=1)
      end if
      ! allocate storage, filled with fill_value ..
      call self%y%Init( (/sdat%nr,self%npix/), status, &
                             long_name='simulated retrieval', &
                             units=sdat%yr%units )
      IF_NOTOK_RETURN(status=1)
      
      ! also convert to vertical column densities?
      if ( self%with_profile ) then
        ! allocate storage, filled with fill_value ..
        call self%y_vcd%Init( (/sdat%nr,self%npix/), status, &
                               long_name='simulated vertical column density', &
                               units=vcd_units )
        IF_NOTOK_RETURN(status=1)
        !! TESTING ...
        !call self%test_dlogy%Init( (/sdat%nlayer,self%npix/), status, &
        !                             long_name='log(hx) - log(ya)', &
        !                             units='log(ppmv)' )
        !IF_NOTOK_RETURN(status=1)
        !call self%test_Kdlogy%Init( (/sdat%nr,self%npix/), status, &
        !                             long_name='K[log(hx) - log(ya)]', &
        !                             units='log(ppmv)' )
        !IF_NOTOK_RETURN(status=1)
        !call self%test_logya_Kdlogy%Init( (/sdat%nr,self%npix/), status, &
        !                             long_name='log(ya) + K[log(hx) - log(ya)]', &
        !                             units='log(ppmv)' )
        !IF_NOTOK_RETURN(status=1)
      end if

      ! any local pixels?
      if ( self%npix > 0 ) then
        
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
          
          ! fill full model profile (volume mixing ratios):
          self%mod_vmr%data(     1:nlev    ,ipix) = c     (i,j,     1:nlev    ,ispec)
          self%mod_vmr%data(nlev+1:nlev_top,ipix) = caloft(i,j,nlev+1:nlev_top,ispec)
          
          ! half level pressures ;
          ! see also below for adhoc change of surface pressure:
          self%mod_hp%data(:,ipix) = hp(i,j,0:nlev_top)
          
          ! copy:
          self%mod_hp0%data(:,ipix) = hp(i,j,0:nlev_top)
          
          !! testing ... 
          !write (gol,*) 'xxx ilayer1  = ', sdat%ilayer1%data(ipix); call goPr
          !write (gol,*) 'xxx ilayer2  = ', sdat%ilayer2%data(ipix); call goPr
          !write (gol,*) 'xxx ir1      = ', sdat%ir1%data(ipix); call goPr
          !write (gol,*) 'xxx ir2      = ', sdat%ir2%data(ipix); call goPr
          !write (gol,*) 'xxx cris hp  = ', sdat%hp%data(:,ipix); call goPr

          ! compute vertical-column-densities per apriori layer in vcd units
          conversion = trim(specunit(ispec))//' -> '//trim(vcd_units)
          ! switch:
          select case ( trim(conversion) )
            !~ volume mixing ratio:
            case ( 'ppb -> 1e15 mlc/cm2' )
              ! loop over model layers:
              do ilev = 1, nlev_top
                ! air mass density:
                dm = abs( self%mod_hp%data(ilev,ipix) - self%mod_hp%data(ilev+1,ipix) )/grav  ! (kg air)/m2
                ! convert:
                !                                 mol tracer ppb         1     mlc   mol air  kg air    m2
                !                                ---------------        ---    ---   -------  ------    ---
                !                                  mol air              ppb    mol    kg air    m2      cm2
                self%mod_vcd%data(ilev,ipix) = self%mod_vmr%data(ilev,ipix) * 1.0e-9 * Avog / xm_air *  dm   * 1.0e-4 / 1.0e15 ! (1e15 mlc trc)/cm2
              end do ! levels
            !~ unknown ...
            case default
              write (gol,'("unsupported conversion `",a,"`")') trim(conversion); call goErr
              TRACEBACK; status=1; return
          end select
          
          ! replace surface pressure with retrieval input if that one is below the model surface;
          ! do this after computing mod_vcd to not change the total mass ...
          !~ index of retrieval input surface hp:
          ilayer = sdat%ilayer1%data(ipix)
          !~ index of model surface hp:
          ilev = 1
          ! scale pressure profile:
          self%mod_hp%data(:,ipix) = self%mod_hp%data(:,ipix) / self%mod_hp%data(ilev,ipix) * sdat%hp%data(ilayer,ipix)
          
          ! loop over retrieval input layers:
          do ilayer = sdat%ilayer1%data(ipix), sdat%ilayer2%data(ipix)
          
            ! fill each retrieval layer with sum of one or more fractions of model layers;
            ! use pressure axis to have mass-conservation; negate to have increasing axis;
            ! units of hx: (1e15 mlc trc)/cm2
            call IntervalSum( -1.0*self%mod_hp%data(:,ipix), self%mod_vcd%data(:,ipix), &
                               -1.0*sdat%hp%data(ilayer,ipix), -1.0*sdat%hp%data(ilayer+1,ipix), &
                               self%hx%data(ilayer,ipix), ilast, status )
            if ( status /= 0 ) then
              write (gol,'("mapping from model to retrieval layer:")'); call goErr
              write (gol,'("pixel           : ",i6)') ipix; call goErr
              write (gol,'("model layers (ph0,ph1,vmr,vcd):")'); call goErr
              do ilev = nlev_top, 1, -1
                write (gol,'(i4," ",2f12.4," ",2es12.4)') ilev, self%mod_hp%data(ilev:ilev+1,ipix), &
                                self%mod_vmr%data(ilev,ipix), self%mod_vcd%data(ilev,ipix); call goErr
              end do
              write (gol,'("retrieval layer : ",i4," ",2f12.4)') ilayer, sdat%hp%data(ilayer:ilayer+1,ipix); call goErr
              write (gol,'("local cell      : ",2i5)') i, j; call goErr
              TRACEBACK; status=1; return
            end if
          
            ! need to convert?
            if ( trim(sdat%yr%units) /= trim(vcd_units) ) then
              ! conversion needed:
              conversion = trim(vcd_units)//' -> '//trim(self%hx%units)
              ! switch:
              select case ( trim(conversion) )
                !~ vcd:
                case ( '1e15 mlc/cm2 -> molec/cm2' )
                  ! scale ..
                  self%hx%data(ilayer,ipix) = self%hx%data(ilayer,ipix) * 1.0e15
                !~ vcd:
                case ( '1e15 mlc/cm2 -> mol m-2' )
                  ! scale:                        1e15 mlc/cm2          *  mlc/(1e15 mlc)  * mole/mlc   cm2/m2 
                  self%hx%data(ilayer,ipix) = self%hx%data(ilayer,ipix) *      1.0e15        / Avog   * 1.0e4
                  !! testing ...
                  !write (gol,*) 'hx := ', self%hx%data(ilayer,ipix); call goErr
                !~ volume mixing ratio:
                case ( '1e15 mlc/cm2 -> ppmv' )
                  ! air mass density:
                  dm = abs( sdat%hp%data(ilayer,ipix) - sdat%hp%data(ilayer+1,ipix) )/grav  ! (kg air)/m2
                  ! convert:
                  !                                (1e15 mlc tr)            1    mol tr   cm2      m2     kg air
                  !                                -------------           ----  ------   ---    ------   -------   ppm
                  !                                    cm2                 1e15   mlc      m2    kg air   mol air
                  self%hx%data(ilayer,ipix) = self%hx%data(ilayer,ipix) * 1.0e15 /Avog * 1.0e4   / dm   * xm_air  * 1e6  ! ppmv
                !~ unknown ...
                case default
                  write (gol,'("unsupported conversion `",a,"`")') trim(conversion); call goErr
                  TRACEBACK; status=1; return
              end select
            end if

          end do ! retr input layers

          ! range of input layers:
          k1 = sdat%ilayer1%data(ipix)
          k2 = sdat%ilayer2%data(ipix)
          ! combined units of kernel and simulation:
          khx_units = trim(sdat%K%units)//' * '//trim(self%hx%units)
          ! loop over retrieval layers:
          do ir = sdat%ir1%data(ipix), sdat%ir2%data(ipix)
            ! apriori profile present?
            if ( sdat%with_apriori ) then
              ! check ...
              if ( trim(self%hx%units) /= trim(sdat%ya%units) ) then
                write (gol,'("simulation units `",a,"` do not match apriori units `",a,"`")') trim(self%hx%units), trim(sdat%ya%units); call goErr
                TRACEBACK; status=1; return
              end if
              ! weighted sum incl. apriori:
              select case ( trim(khx_units) )
                !~ mixing ratios:
                case ( 'ppmv/ppmv * ppmv', '1 * mol m-2' )
                  self%y%data(ir,ipix) = sdat%yra%data(ir,ipix) +      sum( sdat%K%data(ir,k1:k2,ipix) * (     self%hx%data(k1:k2,ipix) -      sdat%ya%data(k1:k2,ipix)  ) )
                !~ mixing ratios, log transform:
                case ( 'ln(ppmv)/ln(ppmv) * ppmv' )
                  ! use minium value to allow log tranforms ...
                  !hx_min = 1.0e-12  ! ppmv
                  hx_min = 1.0e-6   ! ppmv  = 1 ppt
                  ! apply convolution after log transforms:
                  !    ln(y) =      ln(yra) + K [ ln(hx) - ln(ya) ]
                  !       y  = exp( ln(yra) + K [ ln(hx) - ln(ya) ] )
                  self%y%data(ir,ipix) = exp( log(max(hx_min,sdat%yra%data(ir,ipix))) + sum( sdat%K%data(ir,k1:k2,ipix) * ( log(max(hx_min,self%hx%data(k1:k2,ipix))) - log(max(hx_min,sdat%ya%data(k1:k2,ipix))) ) ) )
                  !! testing ...
                  !self%test_dlogy%data    (k1:k2,ipix) =                                                                               log(max(hx_min,self%hx%data(k1:k2,ipix))) - log(max(hx_min,sdat%ya%data(k1:k2,ipix)))
                  !self%test_Kdlogy%data      (ir,ipix) =                                           sum( sdat%K%data(ir,k1:k2,ipix) * ( log(max(hx_min,self%hx%data(k1:k2,ipix))) - log(max(hx_min,sdat%ya%data(k1:k2,ipix))) ) )
                  !self%test_logya_Kdlogy%data(ir,ipix) = log(max(hx_min,sdat%yra%data(ir,ipix))) + sum( sdat%K%data(ir,k1:k2,ipix) * ( log(max(hx_min,self%hx%data(k1:k2,ipix))) - log(max(hx_min,sdat%ya%data(k1:k2,ipix))) ) )
                  !write (gol,*) 'xxx ir              = ', ir; call goPr
                  !write (gol,*) '  x k1,k2           = ', k1,k2; call goPr
                  !write (gol,*) '  x hx              = ', self%hx%data(k1:k2,ipix); call goPr
                  !write (gol,*) '  x ya              = ', sdat%ya%data(k1:k2,ipix); call goPr
                  !write (gol,*) '  x log(hx)         = ', log(max(hx_min,self%hx%data(k1:k2,ipix))); call goPr
                  !write (gol,*) '  x log(ya)         = ', log(max(hx_min,sdat%ya%data(k1:k2,ipix))); call goPr
                  !write (gol,*) '  x log(hx)-log(ya) = ', log(max(hx_min,self%hx%data(k1:k2,ipix))) - log(max(hx_min,sdat%ya%data(k1:k2,ipix))); call goPr
                  !write (gol,*) '  x K               = ', sdat%K%data(ir,k1:k2,ipix); call goPr
                  !write (gol,*) '  x K*(log-log)     = ', sum( sdat%K%data(ir,k1:k2,ipix) * ( log(max(hx_min,self%hx%data(k1:k2,ipix))) - log(max(hx_min,sdat%ya%data(k1:k2,ipix))) ) ); call goPr
                  !write (gol,*) '  x exp K*(log-log) = ', exp( sum( sdat%K%data(ir,k1:k2,ipix) * ( log(max(hx_min,self%hx%data(k1:k2,ipix))) - log(max(hx_min,sdat%ya%data(k1:k2,ipix))) ) ) ); call goPr
                  !write (gol,*) '  x yra             = ', sdat%yra%data(ir,ipix); call goPr
                  !write (gol,*) '  x y               = ', self%y%data(ir,ipix); call goPr
                !~ unknown ...
                case default
                  write (gol,'("unsupported kernel*simulation units `",a,"`")') trim(khx_units); call goErr
                  TRACEBACK; status=1; return
              end select

              ! postprocessing needed?
              if ( trim(sdat%yra%units) /= trim(self%y%units) ) then
                ! conversion:
                conversion = trim(sdat%yra%units)//' -> '//trim(self%y%units)
                ! switch:
                select case ( trim(conversion) )
                  !~ vertical column density to mixing ratio:
                  case ( 'mol m-2 -> ppb' )
                    ! devide by air column:  (mol tracer)/m2    /        (mol air)/m2        * ppb
                    self%y%data(ir,ipix) = self%y%data(ir,ipix) / sdat%air_vcd%data(ir,ipix) * 1e9
                  !~ 
                  case default
                    write (gol,'("unsupported conversion `",a,"`")') trim(conversion); call goErr
                    TRACEBACK; status=1; return
                end select
              end if

            !... no apriori ....
            else

              ! weighted sum:
              select case ( trim(khx_units) )
                !~ partical columns:
                case ( '1/1 * 1e15 mlc/cm2', '1 * 1e15 mlc/cm2', '1 * molec/cm2' )
                  self%y%data(ir,ipix) = sum( sdat%K%data(ir,k1:k2,ipix) * self%hx%data(k1:k2,ipix) )
                !~ unknown ...
                case default
                  write (gol,'("unsupported kernel*simulation units `",a,"`")') trim(khx_units); call goErr
                  TRACEBACK; status=1; return
              end select

              ! postprocessing needed?
              if ( trim(self%hx%units) /= trim(self%y%units) ) then
                ! conversion:
                conversion = trim(self%hx%units)//' -> '//trim(self%y%units)
                ! switch:
                select case ( trim(conversion) )
                  !~ 
                  case default
                    write (gol,'("unsupported conversion `",a,"`")') trim(conversion); call goErr
                    TRACEBACK; status=1; return
                end select
              end if

            end if ! apriori?

          end do ! retr layer
          
          ! for convenience, compute simulated vertical-column-densities for profiles:
          if ( self%with_profile ) then
            ! conversion from retrieval units to vcd units:
            conversion = trim(self%y%units)//' -> '//trim(vcd_units)
            ! switch:
            select case ( trim(conversion) )
              !~ volume mixing ratio:
              case ( 'ppmv -> 1e15 mlc/cm2' )
                ! loop over retrieval layers:
                do ir = sdat%ir1%data(ipix), sdat%ir2%data(ipix)
                  ! air mass density:
                  dm = abs( sdat%hpr%data(ir,ipix) - sdat%hpr%data(ir+1,ipix) )/grav  ! (kg air)/m2
                  ! convert:
                  !                            mol tracer ppm          1     mlc   mol air  kg air    m2
                  !                            --------------         ---    ---   -------  ------    ---
                  !                               mol air             ppm    mol    kg air    m2      cm2
                  self%y_vcd%data(ir,ipix) = self%y%data(ir,ipix) * 1.0e-6 * Avog / xm_air *  dm   * 1.0e-4 / 1.0e15 ! (1e15 mlc trc)/cm2
                end do ! levels
              !~ unknown ...
              case default
                write (gol,'("unsupported conversion `",a,"`")') trim(conversion); call goErr
                TRACEBACK; status=1; return
            end select
          end if ! with profile
          
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
    integer                   ::  dimid_mlay
    integer                   ::  dimid_mlayi
    integer                   ::  dimid_layer
    integer                   ::  dimid_retr
    integer                   ::  ntot
!    integer, allocatable      ::  iglb_all(:)
    
    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': put out for ', t ); call goPr
    
    ! any pixels to be put out?
    if ( sdat%nout > 0 ) then
    
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
        status = NF90_Def_Dim( ncid, 'pixel', sdat%nout, dimid_pixel )
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
        call self%mod_hp%NcDef( ncid, 'mod_hp' , (/dimid_mlayi,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%mod_hp0%NcDef( ncid, 'mod_hp0' , (/dimid_mlayi,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%mod_vmr%NcDef( ncid, 'mod_vmr' , (/dimid_mlay,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%mod_vcd%NcDef( ncid, 'mod_vcd' , (/dimid_mlay,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%hx%NcDef( ncid, 'hx' , (/dimid_layer,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        call self%y%NcDef( ncid, 'y' , (/dimid_retr,dimid_pixel/), status )
        IF_NOTOK_RETURN(status=1)
        if ( self%with_profile) then
          call self%y_vcd%NcDef( ncid, 'y_vcd' , (/dimid_retr,dimid_pixel/), status )
          IF_NOTOK_RETURN(status=1)
          !! TESTING ...
          !call self%test_dlogy%NcDef( ncid, 'test_dlogy' , (/dimid_layer,dimid_pixel/), status )
          !IF_NOTOK_RETURN(status=1)
          !call self%test_Kdlogy%NcDef( ncid, 'test_Kdlogy' , (/dimid_retr,dimid_pixel/), status )
          !IF_NOTOK_RETURN(status=1)
          !call self%test_logya_Kdlogy%NcDef( ncid, 'test_logya_Kdlogy' , (/dimid_retr,dimid_pixel/), status )
          !IF_NOTOK_RETURN(status=1)
        end if

        ! end defintion mode:
        status = NF90_EndDef( ncid )
        IF_NF90_NOTOK_RETURN(status=1)
        
      end if ! root

      ! collect distributed arrays on root and write from there:
      call self%mod_hp%NcPutGather( ncid, sdat%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%mod_hp0%NcPutGather( ncid, sdat%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%mod_vmr%NcPutGather( ncid, sdat%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%mod_vcd%NcPutGather( ncid, sdat%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%hx%NcPutGather( ncid, sdat%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      call self%y%NcPutGather( ncid, sdat%iout_all, status )
      IF_NOTOK_RETURN(status=1)
      if ( self%with_profile) then
        call self%y_vcd%NcPutGather( ncid, sdat%iout_all, status )
        IF_NOTOK_RETURN(status=1)
        !! TESTING ...
        !call self%test_dlogy%NcPutGather( ncid, sdat%nglb, iglb_all, status )
        !IF_NOTOK_RETURN(status=1)
        !call self%test_Kdlogy%NcPutGather( ncid, sdat%nglb, iglb_all, status )
        !IF_NOTOK_RETURN(status=1)
        !call self%test_logya_Kdlogy%NcPutGather( ncid, sdat%nglb, iglb_all, status )
        !IF_NOTOK_RETURN(status=1)
      end if
      
      ! written on root...
      if ( goc%root ) then
        ! close:
        status = NF90_Close( ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if  ! root

    else
    
      ! info ..
      write (gol,'(a,":   no data for this time ...")') rname; call goPr
      
    end if

    ! ok
    status = 0
    
  end subroutine LE_Output_Sat_State_PutOut


end module LE_Output_Sat

