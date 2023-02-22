!###############################################################################
!
! LE_Output_CSO - interface to CAMS Satellite Operator
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN (action) if (status> 0) then; TRACEBACK; action; return; end if
!
!#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################
  
module LE_Output_CSO

  use GO              , only : gol, goPr, goErr

  use CSO, only : T_CSO
  use CSO, only : T_CSO_RcFile
  use CSO, only : T_CSO_Listing
  use CSO, only : T_CSO_Sat_Data
  use CSO, only : T_CSO_Sat_State
  use CSO, only : T_GridMapping

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_LE_Output_CSO_Data
  public  ::  T_LE_Output_CSO_State
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'LE_Output_CSO'

!  ! units used for (intermediate) vertical column densities:
!  character(len=*), parameter   ::  vcd_units = '1e15 mlc/cm2'

  ! --- types --------------------------------
  
  !
  ! satelite data:
  !  footprints
  !  timestamps
  !  kernels
  !
  type T_LE_Output_CSO_Data
    ! key to identify this set:
    character(len=32)                 ::  name
    ! main object for CSO data:  
    type(T_CSO)                       ::  csod
    ! cso settings file:
    type(T_CSO_RcFile)                ::  cso_rcf
    ! basename for cso settings:
    character(len=64)                 ::  cso_rcbase
    ! orbit listing:
    type(T_CSO_Listing)               ::  listing
    ! current file, empty if none:
    character(len=1024)               ::  orbit_filename
    ! orbit data:
    type(T_CSO_Sat_Data)              ::  sdata
    ! mapping from regular grid to footprints:
    type(T_GridMapping)               ::  grid_mapping
    integer                           ::  grid_off(2)
    ! dims:
    integer                           ::  npix
    ! mapping arrays:
    real, pointer                     ::  areas(:)
    integer, pointer                  ::  ii(:), jj(:)
    real, pointer                     ::  ww(:)
    integer, pointer                  ::  iw0(:), nw(:)
    ! keys for output filenames:
    character(len=32)                 ::  id_model
    character(len=64)                 ::  id_expid
    ! tracer name and index:
    character(len=32)                 ::  tracer_name
    integer                           ::  ispec
    !
  contains
    procedure :: Init            => LE_Output_CSO_Data_Init
    procedure :: Done            => LE_Output_CSO_Data_Done
    procedure :: Clear           => LE_Output_CSO_Data_Clear
    !procedure :: GetPixel        => LE_Output_CSO_Data_GetPixel
    procedure :: Setup           => LE_Output_CSO_Data_Setup
    procedure :: PutOut          => LE_Output_CSO_Data_PutOut
  end type T_LE_Output_CSO_Data
  
  !
  ! state simulations
  !
  type T_LE_Output_CSO_State
    ! annote:
    character(len=16)                 ::  key
    character(len=256)                ::  description
    ! flag:
    logical                           ::  filled
    ! simulated state:
    type(T_CSO_Sat_State)             ::  sstate
    !
    ! user defined variables:
    integer                           ::  nuvar
    character(len=64), allocatable    ::  uvarnames(:)  ! (nuvar)
    character(len=64), allocatable    ::  uvarunits(:)  ! (nuvar)
    !
  contains
    procedure :: Init            => LE_Output_CSO_State_Init
    procedure :: Done            => LE_Output_CSO_State_Done
    procedure :: Clear           => LE_Output_CSO_State_Clear
    !procedure :: GetPixel        => LE_Output_CSO_State_GetPixel
    procedure :: Setup           => LE_Output_CSO_State_Setup
    procedure :: PutOut          => LE_Output_CSO_State_PutOut
  end type T_LE_Output_CSO_State


  
contains


  ! ====================================================================
  ! ===
  ! === Data
  ! ===
  ! ====================================================================


  subroutine LE_Output_CSO_Data_Init( self, rcF, rcbase, typ, name, status )
  
    use LE_Radiation    , only : tau  
    use GO     , only : TrcFile
    use GO     , only : goc
    use GO     , only : goMatchValue
    use LE_Grid, only : glb_ugg, ugg
    use LE_Grid, only : dom
    use Indices, only : specname
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_Data), intent(out)    ::  self
    type(TrcFile), intent(in)                   ::  rcF
    character(len=*), intent(in)                ::  rcbase
    character(len=*), intent(in)                ::  typ
    character(len=*), intent(in)                ::  name
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Output_CSO_Data_Init'
    
    ! --- local ----------------------------------
    
    character(len=2000)       ::  line
    character(len=1024)       ::  listing_filename
    real                      ::  glb_bbox(4)
    integer                   ::  mapping_levels
    
    ! --- begin ----------------------------------
    
    ! store:
    self%name = trim(name)
    
    ! read keys for output filenames:
    call rcF%Get( trim(rcbase)//'.model', self%id_model, status )
    IF_NOT_OK_RETURN(status=1)
    call rcF%Get( trim(rcbase)//'.expid', self%id_expid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! prefix:
    self%cso_rcbase = trim(rcbase)//'.'//trim(typ)//'.'//trim(name)
    
    ! tracer name:
    call rcF%Get( trim(self%cso_rcbase)//'.tracer', self%tracer_name, status )
    IF_NOT_OK_RETURN(status=1)
    ! source tracer:
    select case (trim(self%tracer_name))
       case('aod_440nm', 'aod_490nm', 'aod_550nm', 'aod_563nm', &
                   'aod_670nm', 'aod_675nm', 'aod_865nm', 'aod_870nm', &
                   'aod_1020nm' )
            !500 aod
            self%ispec=500
        case('aaod_440nm', 'aaod_490nm', 'aaod_550nm', 'aaod_563nm', &
                   'aaod_670nm', 'aaod_675nm', 'aaod_865nm', 'aaod_870nm', &
                   'aaod_1020nm' )
            !500 aod
            self%ispec=500
        case('angstrom_aeronet', 'angstrom_modis', 'angstrom_polder')       
             self%ispec=600
        case('ssa_440nm', 'ssa_490nm', 'ssa_550nm', 'ssa_563nm', &
                   'ssa_670nm', 'ssa_675nm', 'ssa_865nm', 'ssa_870nm', &
                   'ssa_1020nm' )      
             self%ispec=600
        case('extinction_440nm', 'extinction_490nm', 'extinction_550nm', 'extinction_563nm', &
                   'extinction_670nm', 'extinction_675nm', 'extinction_865nm', 'extinction_870nm', &
                   'extinction_1020nm' )      
             self%ispec=600
        
        case default
             call goMatchValue( trim(self%tracer_name), specname, self%ispec, status )
             IF_NOT_OK_RETURN(status=1)
   end select
    
    ! duplicate settings into cso structure ...
    call self%cso_rcf%Init( rcF%fname, status )
    IF_NOT_OK_RETURN(status=1)

    ! listing filename:
    call rcf%Get( trim(self%cso_rcbase)//'.listing', listing_filename, status )
    IF_NOT_OK_RETURN(status=1)
    ! info ...
    write (gol,'(a,": read listing file: ",a)') rname, trim(listing_filename); call goPr
    ! read listing file:
    call self%listing%Init( listing_filename, status )
    IF_NOT_OK_RETURN(status=1)

    ! init communication:
#if _MPI
    ! init CSO module including MPI communication:
    call self%csod%Init( status, comm=goc%comm )
    IF_NOT_OK_RETURN(status=1)
#else
    call self%csod%Init( status )
    IF_NOT_OK_RETURN(status=1)
#endif
    
    ! init grid mapping, currently for regular carthesian grid only ..
    select case ( trim(ugg%type) )

      ! regular lon/lat:
      case ( 'cartesian-regular' )

        ! bounding box (west,east,south,north) of global(!) domain:
        call glb_ugg%GetBoundingBox( glb_bbox, status )
        IF_NOT_OK_RETURN(status=1)

        ! offset in global grid:
        call dom%Get( status, off=self%grid_off )
        IF_NOT_OK_RETURN(status=1)

        ! mapping level:
        call rcf%Get( trim(self%cso_rcbase)//'.mapping.levels', mapping_levels, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init regular grid mapping:
        ! - west and south bound of global(!) grid;
        ! - grid spacing;
        ! - local grid cell offset in global grid;
        ! - local grid shape
        ! - recursion level
        call self%grid_mapping%Init( glb_bbox(1), glb_ugg%dlon, self%grid_off(1), ugg%nlon, &
                                     glb_bbox(3), glb_ugg%dlat, self%grid_off(2), ugg%nlat, &
                                      mapping_levels, status )
        IF_NOT_OK_RETURN(status=1)
        
      ! not yet ...
      case default
        write (gol,'("unsupported grid type `",a,"`")') trim(ugg%type); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_Data_Init


  ! ***


  subroutine LE_Output_CSO_Data_Done( self, status )

    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_Data), intent(inout)    ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_Data_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    call self%Clear( status )
    IF_NOT_OK_RETURN(status=1)
  
    ! done with grid mapping:
    call self%grid_mapping%Done( status )
    IF_NOT_OK_RETURN(status=1)
  
    ! done with listing:
    call self%listing%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! done with settings:
    call self%cso_rcf%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! done with CSO module:
    call self%csod%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_Data_Done


  ! ***


  subroutine LE_Output_CSO_Data_Clear( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_Data), intent(inout)    ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_Data_Clear'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! defined?
    if ( len_trim(self%orbit_filename) > 0 ) then
      
      ! done with data:
      call self%sdata%Done( status )
      IF_NOT_OK_RETURN(status=1)

    end if  ! orbit data present

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_Data_Clear


!  ! ***
!
!
!  subroutine LE_Output_CSO_Data_GetPixel( self, ipix, status, &
!                                            glbid, lon, lat, y, sigma, covar )
!    
!    ! --- in/out ---------------------------------
!    
!    class(T_LE_Output_CSO_Data), intent(in)       ::  self
!    integer, intent(in)                           ::  ipix
!    integer, intent(out)                          ::  status
!    
!    integer, intent(out), optional                ::  glbid
!    real, intent(out), optional                   ::  lon
!    real, intent(out), optional                   ::  lat
!    real, intent(out), optional                   ::  y
!    real, intent(out), optional                   ::  sigma
!    real, intent(out), optional                   ::  covar(:,:)  ! (nr,nr)
!  
!    ! --- const ----------------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_Data_GetPixel'
!    
!    ! --- local ----------------------------------
!    
!!    integer     ::  iglb
!    
!    ! --- begin ----------------------------------
!    
!    ! not yet ..
!    write (gol,'(a,": not implemented yet")'); call goPr
!    
!!    ! check ...
!!    if ( (ipix < 1) .or. (ipix > self%npix) ) then
!!      write (gol,'("pixel ",i0," out of range 1:",i0)') ipix, self%npix; call goErr
!!      TRACEBACK; status=1; return
!!    end if
!!    
!!    ! global index:
!!    iglb = self%iglb(ipix)
!!    
!!    ! id?
!!    if ( present(glbid) ) glbid = iglb
!!    
!!    ! location:
!!    if ( present(lon) ) lon = self%glb_lon%data(iglb)
!!    if ( present(lat) ) lat = self%glb_lat%data(iglb)
!!    
!!    ! retrieval of total column: SHOULD BE REPLACED BY PROFILE!
!!    if ( present(y) ) then
!!      ! check ...
!!      if ( self%nr /= 1 ) then
!!        write (gol,'("profiles not supported yet")'); call goErr
!!        TRACEBACK; status=1; return
!!      end if
!!      ! first value ...
!!      y = self%yr%data(1,ipix)
!!    end if
!!    
!!    ! std.dev.: SHOULD BE REPLACED BY COVARIANCE!
!!    if ( present(sigma) ) then
!!      ! check ...
!!      if ( self%nr /= 1 ) then
!!        write (gol,'("profiles not supported yet")'); call goErr
!!        TRACEBACK; status=1; return
!!      end if
!!      ! first value ...
!!      sigma = self%sr%data(1,ipix)
!!    end if
!!    
!!    ! std.dev.: SHOULD BE REPLACED BY COVARIANCE!
!!    if ( present(covar) ) then
!!      ! check ..
!!      if ( any( shape(covar) /= (/self%nr,self%nr/) ) ) then
!!        write (gol,'("covar shape (",i0,",",i0,") while nr=",i0)') shape(covar), self%nr; call goPr
!!        TRACEBACK; status=1; return
!!      end if
!!      ! copy:
!!      covar = self%vr%data(:,:,ipix)
!!    end if
!
!    ! ok
!    status = 0
!    
!  end subroutine LE_Output_CSO_Data_GetPixel


  ! ***


  subroutine LE_Output_CSO_Data_Setup( self, tr, status )  
    
    use GO     , only : TDate, wrtgol
    use LE_Grid, only : ugg, glb_ugg
    use CSO    , only : T_CSO_DateTime
    use CSO    , only : CSO_DateTime
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_Data), intent(inout)    ::  self
    type(TDate), intent(in)                       ::  tr(2)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_Data_Setup'
    
    ! --- local ----------------------------------
    
    type(T_CSO_DateTime)      ::  cso_t1, cso_t2
  
    integer                   ::  nglb
    real, pointer             ::  glb_lon(:)       ! (nglb)
    real, pointer             ::  glb_lat(:)       ! (nglb)
    real, pointer             ::  glb_clons(:,:)   ! (nglb,4)
    real, pointer             ::  glb_clats(:,:)   ! (nglb,4)
    logical, pointer          ::  glb_select(:)    ! (nglb)
    integer                   ::  iglb
    real                      ::  bbox(4)
    logical                   ::  indomain

    ! footprints:
    real, pointer             ::  clons(:,:), clats(:,:)   ! (ncorner,npix)
    
!    integer                   ::  npix
!    integer                   ::  ipix
!    real, pointer             ::  pix_clons(:), pix_clats(:)
!    real                      ::  pix_area
!    integer, pointer          ::  ii(:), jj(:)
!    real, pointer             ::  ww(:)
!    integer                   ::  nw
    
    !integer, pointer          ::  pix_ii(:), pix_jj(:)
    !real, pointer             ::  pix_ww(:)
    !integer                   ::  pix_nw
    !integer                   ::  iw
    !integer                   ::  i, j

    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': setup for ', tr ); call goPr
    
    ! clear current storage:
    call self%Clear( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! fill time structures:
    cso_t1 = CSO_DateTime( year=tr(1)%year, month=tr(1)%month, day=tr(1)%day, &
                             hour=tr(1)%hour, min=tr(1)%min, sec=tr(1)%sec )
    cso_t2 = CSO_DateTime( year=tr(2)%year, month=tr(2)%month, day=tr(2)%day, &
                             hour=tr(2)%hour, min=tr(2)%min, sec=tr(2)%sec )
    ! get name of CSO file for orbit with time average
    ! within requested interval:
    call self%listing%SearchFile( cso_t1, cso_t2, 'aver', self%orbit_filename, status )
    IF_NOT_OK_RETURN(status=1)

    ! defined?
    if ( len_trim(self%orbit_filename) > 0 ) then
      ! info ..
      write (gol,'(a,":  orbit file : ",a)') rname, trim(self%orbit_filename); call goPr
      
      ! ~ orbit data
      
      ! inititalize orbit data,
      ! read all footprints available in file:
      call self%sdata%Init( self%cso_rcf, self%cso_rcbase, trim(self%orbit_filename), status )
      IF_NOT_OK_RETURN(status=1)
      
      ! obtain from the orbit:
      ! - number of footprints (global, all pixels in file)
      ! - pointers to arrays with footprint centers
      ! - pointer to to selection flags, should be used to select pixels overlapping with local domain
      call self%sdata%Get( status, nglb=nglb, &
                               glb_lon=glb_lon, glb_lat=glb_lat, &
                               glb_clons=glb_clons, glb_clats=glb_clats, &
                               glb_select=glb_select )
      IF_NOT_OK_RETURN(status=1)
      
      ! select pixels that overlap with local domain;
      ! loop over global pixels:
      do iglb = 1, nglb
        !! check if center is in local domain:
        !call ugg%InDomain( glb_lon(iglb), glb_lat(iglb), glb_select(iglb), status )
        !IF_NOT_OK_RETURN(status=1)
        ! bounding box of pixel:
        bbox(1) = minval(glb_clons(:,iglb))
        bbox(2) = maxval(glb_clons(:,iglb))
        bbox(3) = minval(glb_clats(:,iglb))
        bbox(4) = maxval(glb_clats(:,iglb))
        !~ bounding box (west,east,south,north) should be entirely in global domain:
        call glb_ugg%InDomain( bbox(1), bbox(2), bbox(3), bbox(4), indomain, status )
        IF_NOT_OK_RETURN(status=1)
        ! outside? then skip:
        if ( .not. indomain ) cycle
        !~ check if bounding box (west,east,south,north) covers at least partly the local domain:
        call ugg%OverDomain( bbox(1), bbox(2), bbox(3), bbox(4), glb_select(iglb), status )
        IF_NOT_OK_RETURN(status=1)
      end do ! iglb
      
      ! read orbit, locally store only pixels that are flagged in 'glb_select'
      ! as having overplap with this domain:
      call self%sdata%Read( self%cso_rcf, self%cso_rcbase, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! obtain number of local pixels:
      call self%sdata%Get( status, npix=self%npix )
      IF_NOT_OK_RETURN(status=1)
      
      ! any local pixels?
      if ( self%npix > 0 ) then
      
        ! pointers to (*,pixel) arrays:
        ! - footprint corners
        call self%sdata%Get( status, clons=clons, clats=clats )
        IF_NOT_OK_RETURN(status=1)

        ! info ...
        write (gol,'(a,": compute mapping weights ...")') rname; call goPr
        ! get pointers to mapping arrays for all pixels:
        ! - areas(1:npix)            : pixel area [m2]
        ! - iw0(1:npix), nw(1:npix)  : offset and number of elements in ii/jj/ww
        ! - ii(:), jj(:), ww(:)      : cell and weight arrays for mapping to footprint,
        call self%grid_mapping%GetWeights( clons, clats, &
                                             self%areas, self%iw0, self%nw, self%ii, self%jj, self%ww, &
                                             status )
        IF_NOT_OK_RETURN(status=1)

        ! info ...
        write (gol,'(a,":   store weights ...")') rname; call goPr
        ! store mapping weights, might be saved and used to create gridded averages;
        ! cell indices ii/jj need to be the global index numbers!
        call self%sdata%SetMapping( self%areas, self%nw, &
                                      self%grid_off(1)+self%ii, self%grid_off(2)+self%jj, &
                                      self%ww, status )
        IF_NOT_OK_RETURN(status=1)

      end if  ! any local pixels
      
      ! clear:
      nullify( glb_lon )
      nullify( glb_lat )
      nullify( glb_select )

    end if  ! orbit found in listing

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_Data_Setup


  ! ***


  subroutine LE_Output_CSO_Data_PutOut( self, t, status )

    use GO              , only : TDate, wrtgol
    use LE_Config       , only : outputdir
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_Data), intent(inout)    ::  self
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_Data_PutOut'
    
    ! --- local ----------------------------------

    character(len=1024)       ::  fname
    
    ! --- begin ----------------------------------
    
    ! defined?
    if ( len_trim(self%orbit_filename) > 0 ) then

      ! target file:
      write (fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",2i2.2,"_data.nc")') &
           trim(outputdir), trim(self%id_model), trim(self%id_expid), &
           trim(self%name), t%year, t%month, t%day, t%hour, t%min

      ! info ..
      write (gol,'(a,":   create ",a," ...")') rname, trim(fname); call goPr
      
      ! write:
      call self%sdata%PutOut( trim(fname), status )
      IF_NOT_OK_RETURN(status=1)
    
    else
    
      ! info ..
      write (gol,'(a,":   no data for this time ...")') rname; call goPr
      
    end if
    
    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_Data_PutOut



  ! ====================================================================
  ! ===
  ! === State
  ! ===
  ! ====================================================================


  subroutine LE_Output_CSO_State_Init( self, sdat, status, key, description )
  
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_State), intent(out)   ::  self
    class(T_LE_Output_CSO_Data), intent(in)     ::  sdat
    integer, intent(out)                        ::  status
    character(len=*), intent(in), optional      ::  key
    character(len=*), intent(in), optional      ::  description

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Output_CSO_State_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! defaults:
    self%key = 'state'
    self%description = 'simulated satellite retrievals'
    ! replace?
    if ( present(key) ) then
      if ( len_trim(key) > 0 ) self%key = trim(key)
    end if
    if ( present(description) ) self%description = trim(description)
    
    ! no data yet:
    self%filled = .false.
    
    ! no user-defined variables yet:
    self%nuvar = 0

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_State_Init


  ! ***


  subroutine LE_Output_CSO_State_Done( self, sdat, status )
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_State), intent(inout)   ::  self
    class(T_LE_Output_CSO_Data), intent(in)       ::  sdat
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_State_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%Clear( sdat, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_State_Done


  ! ***


  subroutine LE_Output_CSO_State_Clear( self, sdat, status )
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_State), intent(inout)   ::  self
    class(T_LE_Output_CSO_Data), intent(in)       ::  sdat
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_State_Clear'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! any data?
    if ( self%filled ) then
      
      ! done with state:
      call self%sstate%Done( sdat%sdata, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! reset flag:
      self%filled = .false.
      
      ! user defined variables?
      if ( self%nuvar > 0 ) then
        ! clear:
        deallocate( self%uvarnames, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( self%uvarunits, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
    
    end if  ! orbit loaded

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_State_Clear


  ! ***


!  subroutine LE_Output_CSO_State_GetPixel( self, ipix, status, &
!                                            y )
!    
!    ! --- in/out ---------------------------------
!    
!    class(T_LE_Output_CSO_State), intent(in)      ::  self
!    integer, intent(in)                           ::  ipix
!    integer, intent(out)                          ::  status
!    
!    real, intent(out), optional                   ::  y
!  
!    ! --- const ----------------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_State_GetPixel'
!    
!    ! --- local ----------------------------------
!    
!    ! --- begin ----------------------------------
!    
!!    ! check ...
!!    if ( (ipix < 1) .or. (ipix > self%npix) ) then
!!      write (gol,'("pixel ",i0," out of range 1:",i0)') ipix, self%npix; call goErr
!!      TRACEBACK; status=1; return
!!    end if
!    
!    ! retrieval of total column: SHOULD BE REPLACED BY PROFILE!
!    if ( present(y) ) then
!      write (gol,'("TO BE IMPLEMENTED ...")'); call goErr
!      TRACEBACK; status=1; return
!!      ! check ...
!!      if ( self%nr /= 1 ) then
!!        write (gol,'("profiles not supported yet")'); call goErr
!!        TRACEBACK; status=1; return
!!      end if
!!      ! first value ...
!!      y = self%y%data(1,ipix)
!    end if
!    
!    ! ok
!    status = 0
!    
!  end subroutine LE_Output_CSO_State_GetPixel


  ! ***
  
  
  ! allocate storage for simulations

  subroutine LE_Output_CSO_State_Setup( self, sdat, c, status )
  
!    use Binas          , only : xm_air, grav
!    use Num            , only : IntervalSum
    use Indices        , only : specunit
!    use LE_Grid        , only : ugg
    use LE_Data_Common , only : nlev, nlev_top
    use LE_Data        , only : LE_Data_GetPointer
    use LE_Bound_Common, only : caloft
    use Dims   , only :  nz
    use LE_Radiation    , only : LE_Radiation_Calc
    use LE_Radiation    , only : tau, swbands,ssa,extinction
          ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_State), intent(inout)   ::  self
    class(T_LE_Output_CSO_Data), intent(in)       ::  sdat
    real, intent(in)                              ::  c(:,:,:,:)  ! (nx,ny,nlev,nspec)
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Output_CSO_State_Setup'
    
    ! --- local ----------------------------------
    
!    integer                   ::  npix
!    integer                   ::  ipix
!    integer                   ::  nlayer
!    real, pointer             ::  pix_clons(:), pix_clats(:)
!    real, pointer             ::  pix_hp(:)    ! (nlayer+1)
!    real, allocatable         ::  mod_conc(:)  ! (nlev_top)
!    real, allocatable         ::  mod_hp(:)    ! (0:nlev_top)
!    real, allocatable         ::  hx(:)        ! (nlayer)
!    real                      ::  pix_area
!    integer, pointer          ::  ii(:), jj(:)
!    real, pointer             ::  ww(:)
!    integer                   ::  nw
    integer                   ::  iw, inz
!    integer                   ::  i, j
!    character(len=32)         ::  p_units
!    character(len=32)         ::  y_units
!
!    integer                   ::  ilev
!    character(len=64)         ::  conversion
    real, pointer             ::   hp(:,:,:)   ! (lon,lat,0:nlev)
    real, pointer             ::   cc(:,:,:)   ! (lon,lat,1:nlev)
    real, pointer             ::  tcc(:,:,:)   ! (lon,lat,1)
!    real                      ::  dm
!    integer                   ::  ilayer
!    integer                   ::  ilast
  
    ! user defined dimensions:
    integer                           ::  nudim
    integer                           ::  iudim
    character(len=32)                 ::  udimname
    ! user defined variables:
    integer                           ::  iuvar
    ! pixel index:
    integer                           ::  npix
    integer                           ::  ipix 
    ! pixel arrays:
    real, pointer                     ::  data0(:)          ! (npix)
    real, pointer                     ::  data1(:,:)        ! (:,npix)
     ! optical properties:
    integer               ::  vlen
    integer               ::  lambda
    integer               ::  lambda1,lambda2
    integer               ::  iswband
    real, allocatable     ::  convfact(:,:)
    real, pointer          ::  dens(:,:,:)   ! (lon,lat,alt)    

    character(len=64)       ::  instrument
    ! --- begin ----------------------------------
    
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')    
    IF_NOT_OK_RETURN(status=1)   
    
    
    
    
    ! info ...
    write (gol,'(a,": simulate pixels ...")') rname; call goPr

    ! clear current storage:
    call self%Clear( sdat, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! orbit loaded?
    if ( len_trim(sdat%orbit_filename) > 0 ) then
      
      !
      ! initialize simulation state;
      ! optional arguments:
      !   description='long name'    : used for output attributes
      call self%sstate%Init( sdat%sdata, sdat%cso_rcf, sdat%cso_rcbase, status, &
                           description='simulated retrievals' )
      IF_NOT_OK_RETURN(status=1)
      
      ! number of user defined dimensions:
      call self%sstate%Get( status, nudim=nudim )
      IF_NOT_OK_RETURN(status=1)
      ! loop:
          
          
      do iudim = 1, nudim
        ! get dimension name to be defined:
        call self%sstate%GetDim( iudim, status, name=udimname )
        IF_NOT_OK_RETURN(status=1)
        ! switch:
        select case ( trim(udimname) )
          !~ model layers:
          case ( 'model_layer' )
            ! define:
            call self%sstate%SetDim( udimname, nlev_top, status )
            IF_NOT_OK_RETURN(status=1)
          !~ model layer interfaces:
           case ( 'model_total' )
            ! define:
            call self%sstate%SetDim( udimname,1, status )
            IF_NOT_OK_RETURN(status=1)
          !~Model two variables for aod wlenght and, and aod and aaod
          case ( 'model_sfc' )
            ! define:
            call self%sstate%SetDim( udimname,1, status )
            IF_NOT_OK_RETURN(status=1)
          !~Model two variables for aod wlenght and, and aod and aaod
          case ( 'model_total_ae' )
            ! define:
            call self%sstate%SetDim( udimname,2, status )
            IF_NOT_OK_RETURN(status=1)          
          !~ model layer interfaces:   
          case ( 'model_layeri' )
            ! define:
            call self%sstate%SetDim( udimname, nlev_top+1, status )
            IF_NOT_OK_RETURN(status=1)
          !~ unknown ...
          case default
            write (*,'("unsupported dimension `",a,"`")') trim(udimname)
            TRACEBACK; status=1; return
        end select
      end do ! idim
      
      ! dimension defined, allocate storage:
      call self%sstate%EndDef( status )
      IF_NOT_OK_RETURN(status=1)
      
      ! number of user defined variables:
      call self%sstate%Get( status, nuvar=self%nuvar )
      IF_NOT_OK_RETURN(status=1)
      ! any user defined?
      if ( self%nuvar > 0 ) then
        ! storage for variable names and units:
        allocate( self%uvarnames(self%nuvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%uvarunits(self%nuvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill:
        call self%sstate%Get( status, uvarnames=self%uvarnames, uvarunits=self%uvarunits )
        IF_NOT_OK_RETURN(status=1)
      end if ! nuvar > 0
      
 !This will calculate tau when concentrations changed     
!#ifdef with_radiation
    ! radiation properties:
    call LE_Radiation_Calc( c, status)
    IF_NOT_OK_RETURN(status=1)
!#endif
      ! loop over variables:
      do iuvar = 1, self%nuvar
        ! info ..
        write (*,'(a,": user defined variable: ",a)') rname, trim(self%uvarnames(iuvar))
        ! switch:
        select case ( trim(self%uvarnames(iuvar)) )

          !~ model concentrations
          case ( 'mod_conc' )
            ! check units:
            if ( trim(self%uvarunits(iuvar)) /= trim(specunit(sdat%ispec)) ) then
              write (*,'(a,": variable `",a,"` requires conversion to `",a,"` from `",a,"`")') &
                      rname, trim(self%uvarnames(iuvar)), trim(self%uvarunits(iuvar)), trim(specunit(sdat%ispec))
            end if
            ! get pointer to target array with shape (nz,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(:,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers and aloft:
                  data1(     1:nlev    ,ipix) = data1(     1:nlev    ,ipix) + c     (sdat%ii(iw),sdat%jj(iw),     1:nlev    ,sdat%ispec) * sdat%ww(iw)/sdat%areas(ipix)
                  data1(nlev+1:nlev_top,ipix) = data1(nlev+1:nlev_top,ipix) + caloft(sdat%ii(iw),sdat%jj(iw),nlev+1:nlev_top,sdat%ispec) * sdat%ww(iw)/sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix
            
            !~model column
            case ( 'mod_column' )
            ! get pointer to target array with shape (1,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:    
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(1,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers and aloft:
                  do inz=1,nlev    
                    data1(1,ipix) =data1(1,ipix)+ c     (sdat%ii(iw),sdat%jj(iw),     inz    ,sdat%ispec)
                  end do !inz
                   do inz=nlev+1,nlev_top    
                    data1(1,ipix) =data1(1,ipix)+  caloft(sdat%ii(iw),sdat%jj(iw),inz,sdat%ispec)
                  end do !inz
                  data1(1,ipix)= data1(1,ipix)* sdat%ww(iw) /sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix

            !~model surface
            case ( 'mod_surface' )
            ! get pointer to target array with shape (1,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:    
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(1,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model surface:
                    data1(1,ipix) = c     (sdat%ii(iw),sdat%jj(iw),     1    ,sdat%ispec)* sdat%ww(iw) /sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix
            
            
             !~ model aerosol extinction coefficient
          case ( 'mod_aec' )
            ! get pointer to target array with shape (nz,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            vlen = len_trim(sdat%tracer_name)
            read (sdat%tracer_name(12:vlen-2),*) lambda    
            ! band index:
            call swbands%FindBand( lambda*0.001, iswband, status )
            IF_NOT_OK_RETURN(status=1)
            
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(:,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers and aloft:
                  data1(     1:nlev    ,ipix) = data1(     1:nlev    ,ipix) + extinction   (sdat%ii(iw),sdat%jj(iw),     1:nlev    ,iswband) * sdat%ww(iw)/sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix
            
            !~ model aod
          case ( 'mod_aod' )
            ! get pointer to target array with shape (1,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            vlen = len_trim(sdat%tracer_name)
            read (sdat%tracer_name(5:vlen-2),*) lambda
            ! band index:
                
            call swbands%FindBand( lambda*0.001, iswband, status )
            IF_NOT_OK_RETURN(status=1)
    
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(1,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers and aloft:
                  do inz=1,nz    
                    data1(1,ipix) =data1(1,ipix)+ tau(sdat%ii(iw),sdat%jj(iw),inz,iswband)
                    
                  end do !inz
                  data1(1,ipix)= data1(1,ipix)* sdat%ww(iw) /sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix
            !print *,'Santiago mod_aod concentration mean', sum(c)/size(c)
            !print *,'Santiago mod_aod Tau', sum(tau)/size(tau)
         !~ model angstrom
          case ( 'mod_angstrom' )
            ! get pointer to target array with shape (2,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            vlen = len_trim(sdat%tracer_name)
            
            read (sdat%tracer_name(10:vlen),*) instrument
                  
            select case ( trim(instrument) )
                case('aeronet')
                    lambda1=440
                    lambda2=870
                case('modis')
                    lambda1=470
                    lambda2=650
                case('polder')
                    lambda1=443
                    lambda2=865
                case default
                    write (*,'("unsupported instrument name for angstrom exponent`",a,"`")') trim(instrument)
                                 TRACEBACK; status=1; return
               end select
            ! band index first AOD:
            call swbands%FindBand( lambda1*0.001, iswband, status )
            IF_NOT_OK_RETURN(status=1)
    
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(1,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers and aloft:
                  do inz=1,nz    
                    data1(1,ipix) =data1(1,ipix)+ tau(sdat%ii(iw),sdat%jj(iw),inz,iswband)
                    
                  end do !inz
                  data1(1,ipix)= data1(1,ipix)* sdat%ww(iw) /sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix

            call swbands%FindBand( lambda2*0.001, iswband, status )
            IF_NOT_OK_RETURN(status=1)
            !band index second AOD:
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(2,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers and aloft:
                  do inz=1,nz    
                    data1(2,ipix) =data1(2,ipix)+ tau(sdat%ii(iw),sdat%jj(iw),inz,iswband)
                    
                  end do !inz
                  data1(2,ipix)= data1(2,ipix)* sdat%ww(iw) /sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix

           !~ model aod
          case ( 'mod_aaod' )
            ! get pointer to target array with shape (1,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            vlen = len_trim(sdat%tracer_name)
            read (sdat%tracer_name(6:vlen-2),*) lambda
            ! band index:
                
            call swbands%FindBand( lambda*0.001, iswband, status )
            IF_NOT_OK_RETURN(status=1)
    
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(1,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers and aloft:
                  do inz=1,nz    
                    data1(1,ipix) =data1(1,ipix)+  tau(sdat%ii(iw),sdat%jj(iw),inz,iswband)*(1-ssa(sdat%ii(iw),sdat%jj(iw),inz,iswband))
                  end do !inz
                  data1(1,ipix)= data1(1,ipix)* sdat%ww(iw) /sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix

  
             !~ model Single Scatterig Albedo
          case ( 'mod_ssa' )
            ! get pointer to target array with shape (2,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            vlen = len_trim(sdat%tracer_name)
            read (sdat%tracer_name(5:vlen-2),*) lambda
            ! band index:
                
            call swbands%FindBand( lambda*0.001, iswband, status )
            IF_NOT_OK_RETURN(status=1)
    
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(1,ipix) = 0.0
                data1(2,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers and aloft:
                  do inz=1,nz    
                    data1(1,ipix) =data1(1,ipix)+ tau(sdat%ii(iw),sdat%jj(iw),inz,iswband)
                    data1(2,ipix) =data1(2,ipix)+  tau(sdat%ii(iw),sdat%jj(iw),inz,iswband)*(1-ssa(sdat%ii(iw),sdat%jj(iw),inz,iswband))
                  end do !inz
                  data1(1,ipix)= data1(1,ipix)* sdat%ww(iw) /sdat%areas(ipix)
                  data1(2,ipix)= data1(2,ipix)* sdat%ww(iw) /sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix

  
  
              
          !~ model half-level pressures:
          case ( 'mod_hp' )
            ! pointer to meteo data, check units:
            call LE_Data_GetPointer( 'hp', hp, status, check_units=trim(self%uvarunits(iuvar)) )    
            IF_NOT_OK_RETURN(status=1)
            ! get pointer to target array with shape (nz,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(:,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution:
                  data1(:,ipix) = data1(:,ipix) + hp(sdat%ii(iw),sdat%jj(iw),:) * sdat%ww(iw)/sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix

          !~ cloud cover: total:
          case ( 'mod_tcc' )
            ! pointer to meteo data, check units:
            call LE_Data_GetPointer( 'tcc', tcc, status, check_units=trim(self%uvarunits(iuvar)) )    
            IF_NOT_OK_RETURN(status=1)
            ! get pointer to target array with shape (npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data0=data0 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data0(ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution:
                  data0(ipix) = data0(ipix) + tcc(sdat%ii(iw),sdat%jj(iw),1) * sdat%ww(iw)/sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix

          !~ cloud cover: profile:
          case ( 'mod_cc' )
            ! pointer to meteo data, check units:
            call LE_Data_GetPointer( 'cc', cc, status, check_units=trim(self%uvarunits(iuvar)) )    
            IF_NOT_OK_RETURN(status=1)
            ! get pointer to target array with shape (nz,npix):
            call self%sstate%GetData( status, name=self%uvarnames(iuvar), data1=data1 )
            IF_NOT_OK_RETURN(status=1)
            ! loop over pixels:
            do ipix = 1, sdat%npix
              ! any source contributions?
              if ( sdat%nw(ipix) > 0 ) then
                ! init sum:
                data1(:,ipix) = 0.0
                ! loop over source contributions:
                do iw = sdat%iw0(ipix)+1, sdat%iw0(ipix)+sdat%nw(ipix)
                  ! add contribution from model layers, keep zero top:
                  data1(1:nlev,ipix) = data1(1:nlev,ipix) + cc(sdat%ii(iw),sdat%jj(iw),1:nlev) * sdat%ww(iw)/sdat%areas(ipix)
                end do ! iw
              end if ! sdat%nw > 0
            end do ! ipix

          !~
          case default
            write (*,'("unsupported variable name `",a,"`")') trim(self%uvarnames(iuvar))
            TRACEBACK; status=1; return
        end select

      end do ! user defined variables
      
      ! apply formula: kernel convolution etc:
      call self%sstate%ApplyFormulas( sdat%sdata, status )
      IF_NOT_OK_RETURN(status=1)
    
      ! reset flag:
      self%filled = .true.
    
    end if  ! any global pixels

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_State_Setup


  ! ***


  subroutine LE_Output_CSO_State_PutOut( self, sdat, t, status )
  
    use GO       , only : TDate, wrtgol
    use LE_Config, only : outputdir
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Output_CSO_State), intent(inout)   ::  self
    class(T_LE_Output_CSO_Data), intent(in)       ::  sdat
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Output_CSO_State_PutOut'
    
    ! --- local ----------------------------------
    
    character(len=1024)       ::  fname
    
    ! --- begin ----------------------------------
    
    ! info ...
    call wrtgol( rname//': put out for ', t ); call goPr
    
    ! any pixels to be put out?
    if ( len_trim(sdat%orbit_filename) > 0 ) then

      ! target file:
      write (fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",2i2.2,"_",a,".nc")') &
           trim(outputdir), trim(sdat%id_model), trim(sdat%id_expid), &
           trim(sdat%name), t%year, t%month, t%day, t%hour, t%min, &
           trim(self%key)

      ! info ..
      write (gol,'(a,":   create ",a," ...")') rname, trim(fname); call goPr

      ! write:
      call self%sstate%PutOut( sdat%sdata, trim(fname), status )
      IF_NOT_OK_RETURN(status=1)

    else
    
      ! info ..
      write (gol,'(a,":   no data for this time ...")') rname; call goPr
      
    end if

    ! ok
    status = 0
    
  end subroutine LE_Output_CSO_State_PutOut


end module LE_Output_CSO

