!###############################################################################
!
! NAME
!
!   LE_Output_Rad  -  LOTOS-EUROS output of 3D fields
!
! HISTORY
!
!   2007 may, Arjo Segers, TNO
!
!   2022 jun, Janot Tokaya, TNO, make the names of the optical params no 
!                                longer hard coded, but dependent on chosen swbs
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

module LE_Output_Rad

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf  
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_Rad

  public  ::  LE_Output_Rad_Init, LE_Output_Rad_Done
  public  ::  LE_Output_Rad_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Rad'


  ! maximum number of supported data sets:
  integer, parameter  ::  ndat = 200
  
  ! maximum number of levels in output:
  integer, parameter  ::  maxlev = 10

  ! list of keys for the optical output
  character(len=32),allocatable     ::  aod_keylist(:)    
  character(len=32),allocatable     ::  ssa_col_keylist(:) 
  character(len=32),allocatable     ::  extinction_keylist(:) 
  character(len=32),allocatable     ::  tau_keylist(:) 
  character(len=32),allocatable     ::  ssa_keylist(:) 
  character(len=32),allocatable     ::  asy_keylist(:)


  ! --- types ------------------------------

  type T_LE_Dat
    character(len=32)      ::  name
    character(len=32)      ::  unit
    integer                ::  rank
    logical                ::  const
  end type T_LE_Dat

  type T_LE_Output_Rad
    ! name for this file:
    character(len=16)           ::  typ
    character(len=16)           ::  name
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
    ! replace existing files ?
    logical                     ::  replace
    ! file opened ?
    logical                     ::  opened
    ! current time range:
    type(TDate)                 ::  tr(2)
    ! time resolution:
    real                        ::  dhour
    ! collect: daily, instant
    character(len=32)           ::  collect
    ! time record counter:
    integer                     ::  itrec
    ! state name:
    character(len=16)           ::  state
    ! file name:
    character(len=256)          ::  fname
    ! file handle:
    integer                     ::  ncid
    ! dimension handles:
    integer                     ::  dimid_lon
    integer                     ::  dimid_lat
    integer                     ::  dimid_lev
    integer                     ::  dimid_time
    ! dimension variables:
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    integer                     ::  varid_lev
    integer                     ::  varid_time
    !integer                     ::  varid_time_day
    integer                     ::  varid_time_dtg
    ! database with supported variables:
    type(T_LE_Dat)              ::  LE_Dat(ndat)
    ! tracer variables:
    integer                     ::  ndat
    integer, pointer            ::  idat(:)
    character(len=32), pointer  ::  name_dat(:)
    character(len=32), pointer  ::  unit_dat(:)
    real, pointer               ::  unitconv(:)
    integer, pointer            ::  varid_dat(:)
    ! level selection:
    character(len=16)           ::  levtype
    integer                     ::  nlev
    integer, pointer            ::  ilev(:)
    real                        ::  heights(maxlev)
    ! wavelength selection:
    integer                     ::  nband
    integer, allocatable        ::  iband(:)   !  (nband) band index in "swbands"
    ! grads ctl file ?
    logical                     ::  grads_ctl
    character(len=256)          ::  grads_ctl_file
    character(len=256)          ::  grads_ctl_base
    integer                     ::  grads_ctl_nt
    type(TDate)                 ::  grads_ctl_t1
    type(TIncrDate)             ::  grads_ctl_dt
    ! bounding box
    integer                     ::  i1, i2, ni
    integer                     ::  j1, j2, nj
    real                        ::  westb, southb
        
  end type T_LE_Output_Rad


contains


  ! ====================================================


  subroutine SetDat( d, name, unit, rank, const )

    ! --- in/out ----------------------------------

    type(T_LE_Dat), intent(out)     ::  d
    character(len=*), intent(in)    ::  name
    character(len=*), intent(in)    ::  unit
    integer, intent(in)             ::  rank
    logical, intent(in)             ::  const

    ! --- begin ----------------------------------

    d%name = name
    d%unit = unit
    d%rank = rank
    d%const = const

  end subroutine SetDat


  ! ====================================================


  subroutine LE_Output_Rad_Init( leo, rcF, rckey, typ, name, state, status )

    use GO     , only : TrcFile
    use GO     , only : goMatchValues, goSplitString, goReadFromLine
    use GO     , only : AnyDate
    use Grid   , only : GetDistribution
    use Dims   , only : nx, ny, nz
    use LE_Data_Common, only : nlev_top !change by JPT
    
    use LE_Grid, only : ugg
    use LE_Radiation_SWBands, only : T_SWBands
    use LE_Radiation, only : swbands

    ! --- in/out --------------------------------

    type(T_LE_Output_Rad), intent(out)    ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  state
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Rad_Init'

    ! rckey extensions:

    ! --- local ---------------------------------

    character(len=64)     ::  basekey,swb_key
    character(len=1024)   ::  field_names, lambda_names, field_names_rad_long
    integer               ::  idat
    character(len=32)     ::  level_names
    integer               ::  ilev
    integer               ::  i, i1, i2, iband

    character(len=512)    ::  key
    real                  ::  west, south, east, north
    real, pointer         ::  ff(:,:)
    
    ! wavelength bands:
    type(T_SWBands)                   ::  self
    real                              ::  lambda

    ! --- begin ---------------------------------

    ! store name:
    leo%typ  = typ
    leo%name = name

    ! init common stuff:
    call leo%com%Init( rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)

    ! replace existing files ?
    call rcF%Get( trim(rckey)//'.replace', leo%replace, status )
    IF_NOTOK_RETURN(status=1)

    ! write GrADS ctl file ?
    call rcF%Get( trim(rckey)//'.ctl', leo%grads_ctl, status )
    IF_NOTOK_RETURN(status=1)

    ! base key:
    write (basekey,'(a,".",a,".",a)') trim(rckey), trim(typ), trim(name)
    
    ! collect daily or instant
    call rcF%Get( trim(basekey)//'.collect', leo%collect, status )
    IF_NOTOK_RETURN(status=1)

    ! output time resolution:
    call rcF%Get( trim(basekey)//'.dhour', leo%dhour, status )
    IF_NOTOK_RETURN(status=1)

    ! tracer names:
    call rcF%Get( trim(basekey)//'.fields', field_names, status )
    IF_NOTOK_RETURN(status=1)

    ! level type::
    call rcF%Get( trim(basekey)//'.levtype', leo%levtype, status )
    IF_NOTOK_RETURN(status=1)

    ! level descriptions:
    call rcF%Get( trim(basekey)//'.levels', level_names, status )
    IF_NOTOK_RETURN(status=1)
    ! layers or heights ?
    select case ( trim(leo%levtype) )
      !~ model levels:
      case ( 'levels' )
        ! setup storage for level indices (surface + levels + upper boundary)
        allocate( leo%ilev(1:100) )
        ! set selected level indices:
                ! switch:
        if ( trim(level_names) == 'all' ) then
          ! all:
          leo%nlev = nz
          do i = 1, leo%nlev
            leo%ilev(i) = i
          end do
        else if ( trim(level_names) == 'top' ) then
          ! top change by JPT:
          leo%nlev = nlev_top
          do i = 1, leo%nlev
            leo%ilev(i) = i
          end do
        else if ( index(trim(level_names),':') > 0 ) then
          ! extract range:
          call goReadFromLine( level_names, i1, status, sep=':' )
          IF_NOTOK_RETURN(status=1)
          call goReadFromLine( level_names, i2, status )
          IF_NOTOK_RETURN(status=1)
          ! count:
          leo%nlev = i2 - i1 + 1
          ! store indices:
          do i = 1, leo%nlev
            leo%ilev(i) = i1 - 1 + i
          end do
        else
          ! set selected level indices:
          call goMatchValues( level_names, 0, nz+1, leo%nlev, leo%ilev, status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! info ...
        write (gol,'("selected levels for rad output:")'); call goPr
        do i = 1, leo%nlev
          ilev = leo%ilev(i)
          write (gol,'("  ",i3,"  ",i3)') i, ilev; call goPr
        end do
      
      !~ height levels:
      case ( 'heights', 'elevations' )
        ! extract heights:
        call goSplitString( trim(level_names), leo%nlev, leo%heights, status )
        IF_NOTOK_RETURN(status=1)
        ! dummy ...
        allocate( leo%ilev(leo%nlev) )
        leo%ilev = -999
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! line with wavelength description:
    call rcF%Get( trim(basekey)//'.wavelengths', lambda_names, status )
    IF_NOTOK_RETURN(status=1)
    
    ! maximum number of wavelengths:
    allocate( self%lambda(swbands%n), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! getting the swbands from the RCl
    
    ! switch:
    select case ( trim(lambda_names) )
      !~ all:
      case ( 'all' )
        ! copy all:
        self%n = swbands%n
        write (gol,'("selected wavelengths (nm) for rad output:")'); call goPr
        ! fill:
        do iband = 1, self%n
          self%lambda(iband) = swbands%lambda(iband)
          write (gol,'("  ",i3,"  ",i0)') iband, int(1.0e3 * swbands%lambda(iband)); call goPr
        end do
        
        
      !~ selection:
      case default
        ! init counter:
        self%n = 0
        write (gol,'("selected wavelengths (nm) for rad output:")'); call goPr
        ! loop over elements of line:
        do
          ! leave if empty:
          if ( len_trim(lambda_names) == 0 ) exit
          ! increase counter:
          self%n = self%n + 1
          ! check ..
          if ( self%n > swbands%n ) then
            write (gol,'("number of requested wavelengths exceeds maximum of ",i0)') swbands%n; call goErr
            TRACEBACK; status=1; return
          end if
          ! extract value
          call goReadFromLine( lambda_names, lambda, status, sep=' ')
          IF_NOTOK_RETURN(status=1)
          ! search:
          call swbands%FindBand( lambda*0.001, iband, status )
          IF_NOTOK_RETURN(status=1)
          
          self%lambda(self%n) = swbands%lambda(iband)
          
          write (gol,'("  ",i3,"  ",i0)') iband, int(1.0e3 * swbands%lambda(iband)); call goPr
          
        end do ! iband
      !
    end select
    

    ! init empty:
    do i = 1, size(leo%LE_Dat)
      call SetDat( leo%LE_Dat(i), '', '', 0, .false. )
    end do

    ! define supported data:
    !   name   :  used in rcfile to identify this field;
    !   unit   :  in which the corresponding field in LE is expressed;
    !             the model units are converted to 'cf'-conventions
    !             (SI-units) on output;
    !   rank   :  2 for 2D field, and 3 for 3D
    !   const  :  logical, T for constant fields
    !                            name               unit    rank  const
    i = 1
    ! AOD
    call SetDat( leo%LE_Dat(i), 'AOD2'            , '1'     , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'angstrom_modis'  , '1'     , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'angstrom_aeronet', '1'     , 2, .false. ); i=i+1
    call SetDat( leo%LE_Dat(i), 'angstrom_polder', '1'     , 2, .false. ); i=i+1
    
    !
    ! aerosol optical parameters for comparison with aeronet etc:
    !    

    allocate( aod_keylist(self%n))  
    allocate( ssa_col_keylist(self%n)) 
    allocate( extinction_keylist(self%n))     
    allocate( tau_keylist(self%n)) 
    allocate( ssa_keylist(self%n)) 
    allocate( asy_keylist(self%n))         
    
    field_names_rad_long = ''
    
    ! add variables that are not wavelength dependent to list of output.
    if (index(field_names, 'angstrom_modis ') /= 0) then    
        field_names_rad_long = trim(field_names_rad_long)//' '//'angstrom_modis'   
    end if
    
    if (index(field_names, 'angstrom_aeronet ') /= 0) then    
        field_names_rad_long = trim(field_names_rad_long)//' '//'angstrom_aeronet'   
    end if   
    
    do iband = 1, self%n
    !do iw = 1, nw    
    
      ! aod variables. We only include them if they are in field_names
      if (index(field_names, 'aod ') /= 0) then    
        write (swb_key,'("aod_",i0,"nm")') int(1.0e3 * self%lambda(iband))
        call SetDat( leo%LE_Dat(i), swb_key , '1', 2, .false. ); i=i+1
        aod_keylist(iband) = swb_key
        field_names_rad_long = trim(field_names_rad_long)//' '//trim(swb_key)
      end if
    
      ! ssa_col: single scatter albedo column integrated
      if (index(field_names, 'ssa_col ') /= 0) then          
        write (swb_key,'("ssa_col_",i0,"nm")') int(1.0e3 * self%lambda(iband))
        call SetDat( leo%LE_Dat(i), swb_key , '1', 2, .false. ); i=i+1
        ssa_col_keylist(iband) = swb_key
        field_names_rad_long = trim(field_names_rad_long)//' '//trim(swb_key)
      end if

      ! extinction coeff
      if (index(field_names, 'extinction ') /= 0) then    
        write (swb_key,'("extinction_",i0,"nm")') int(1.0e3 * self%lambda(iband))
        call SetDat( leo%LE_Dat(i), swb_key , 'm-1', 3, .false. ); i=i+1
        extinction_keylist(iband) = swb_key
        field_names_rad_long = trim(field_names_rad_long)//' '//trim(swb_key)
      end if

      ! tau
      if (index(field_names, 'tau ') /= 0) then    
        write (swb_key,'("tau_",i0,"nm")') int(1.0e3 * self%lambda(iband))
        call SetDat( leo%LE_Dat(i), swb_key , '1', 3, .false. ); i=i+1
        tau_keylist(iband) = swb_key
        field_names_rad_long = trim(field_names_rad_long)//' '//trim(swb_key)
      end if

      ! single scattering albedo
      if (index(field_names, 'ssa ') /= 0) then    
        write (swb_key,'("ssa_",i0,"nm")') int(1.0e3 * self%lambda(iband))
        call SetDat( leo%LE_Dat(i), swb_key , '1', 3, .false. ); i=i+1
        ssa_keylist(iband) = swb_key
        field_names_rad_long = trim(field_names_rad_long)//' '//trim(swb_key)
      end if

      ! assymetry
      if (index(field_names, 'asy ') /= 0) then    
        write (swb_key,'("asy_",i0,"nm")') int(1.0e3 * self%lambda(iband))
        call SetDat( leo%LE_Dat(i), swb_key , '1', 3, .false. ); i=i+1 
        asy_keylist(iband) = swb_key  
        field_names_rad_long = trim(field_names_rad_long)//' '//trim(swb_key)
      end if

    end do
    
    ! setup storage for tracer fields:
    allocate( leo%idat     (ndat) )
    allocate( leo%name_dat (ndat) )
    allocate( leo%unit_dat (ndat) )
    allocate( leo%unitconv (ndat) )
    allocate( leo%varid_dat(ndat) )

    ! search for requested names:
    call goMatchValues(  trim(adjustl(field_names_rad_long)), leo%LE_Dat(:)%name, &
                          leo%ndat, leo%name_dat, leo%idat, &
                          status )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected fields for data output:")'); call goPr
    do i = 1, leo%ndat
      idat = leo%idat(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3, a10," ",a10,")")') &
                  i, leo%name_dat(i), &
                  idat, trim(leo%LE_Dat(idat)%name), trim(leo%LE_Dat(idat)%unit); call goPr
    end do

    ! state name:
    leo%state = trim(state)

    ! files not open yet:
    leo%opened = .false.

    ! no time range set yet:
    leo%tr(1) = AnyDate()
    leo%tr(2) = AnyDate()

    ! init GrADS stuff:
    if ( leo%grads_ctl ) then
      ! no times written yet:
      leo%grads_ctl_nt = 0
    end if

    ! bounding box
    call rcF%Get( trim(basekey)//'.bounding_box', key, status )
    IF_NOTOK_RETURN(status=1)

    ! empty?
    if (len_trim(key) == 0) then
      ! full domain
      leo%i1 = 1
      leo%i2 = nx
      leo%ni = nx
      leo%j1 = 1
      leo%j2 = ny
      leo%nj = ny
      leo%westb = ugg%longitude_bnds(1,1,1)
      leo%southb = ugg%latitude_bnds(1,1,1)
    else
      
      ! not yet ...
      write (gol,'("no output subset supported for domain decomposition yet")'); call goErr
      TRACEBACK; status=1; return
      
      select case ( trim(ugg%type) ) 
        
        case ( 'cartesian-regular') 
          ! read domain from key
          read(key,*,iostat=status) west, east, south, north
          if(status/=0) then
            write (gol,'("could not read domain from key: ",a)') trim(key); call goErr
            TRACEBACK; status=1; return
          endif

          ! Check if bounding box is in run domain
          if ( west < ugg%longitude_bnds_1d(1,1) .or. east > ugg%longitude_bnds_1d(2,ugg%nlon) .or. &
               south < ugg%latitude_bnds_1d(1,1) .or. north > ugg%latitude_bnds_1d(2,ugg%nlat) ) then
            write( gol, '("Bounding box domain is (partly) outside run domain")' ) ; call goErr
            write( gol, '("Run domain: ", 4f8.2)' ) ugg%longitude_bnds_1d(1,1),ugg%longitude_bnds_1d(2,ugg%nlon),ugg%latitude_bnds_1d(1,1),ugg%latitude_bnds_1d(2,ugg%nlat); call goErr
            write( gol, '("Bounding Box domain: ", 4f8.2)' ) west, east, south, north ; call goErr
            TRACEBACK;status=1;return
          endif

          ! for safety
          nullify(ff)
          ! get cell range covered by box
          call ugg%GetDistribution(west,east,south,north,leo%i1,leo%i2,leo%j1,leo%j2,ff,status)
          IF_NOTOK_RETURN(status=1)
          !clear, fractions not used
          if ( associated(ff) ) deallocate(ff)
          ! set shape
          leo%ni = leo%i2-leo%i1+1
          leo%nj = leo%j2-leo%j1+1
          ! set west/south bounds
          leo%westb  = ugg%longitude_bnds_1d(1,leo%i1)
          leo%southb = ugg%latitude_bnds_1d(1,leo%j1)
        case default 
          write( gol, '("Definition of bounding box not clear for grid-type: ", a)' ) trim(ugg%type) ; call goErr
          TRACEBACK;status=1;return
      end select
    end if

    ! ok
    status = 0

  end subroutine LE_Output_Rad_Init


  ! ***


  subroutine LE_Output_Rad_Done( leo, status )

#ifdef with_netcdf  
    use NetCDF          , only : NF90_Close
#endif

    ! --- in/out --------------------------------

    type(T_LE_Output_Rad), intent(inout)   ::  leo
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Rad_Done'

    character(len=256) :: commandline

    ! --- begin ---------------------------------

    ! file opened ?
    if ( leo%opened ) then
      ! close:
#ifdef with_netcdf
      status = NF90_Close( leo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
#else
    stop 'not compiled with netcdf support'
#endif
      ! reset flag:
      leo%opened = .false.
      ! write GrADS ctl file if necessary:
      call Write_GrADS_Ctl( leo, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! clear storage for tracer fields:
    deallocate( leo%idat      )
    deallocate( leo%name_dat  )
    deallocate( leo%unit_dat  )
    deallocate( leo%unitconv  )
    deallocate( leo%varid_dat )

    ! layers or heights ?
    select case ( trim(leo%levtype) )
      !~ model levels:
      case ( 'levels' )
        ! clear storage for level indices:
        deallocate( leo%ilev )
      !~ height levels:
      case ( 'heights', 'elevations' )
        ! nothing to be done ..
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! done with common stuff ...
    call leo%com%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Rad_Done


  ! ***


  subroutine LE_Output_Rad_PutOut( leo, t, c, status )

    use Binas  , only : xm_air
    use GO     , only : goc
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>), rTotal, dTotal, iTotal
    use GO     , only : wrtgol, Precisely, MidNight
    use GO     , only : goMatchValue
    use Num    , only : Interp_Lin
    use LE_Grid, only : glb_ugg
    use C3PO   , only : T_Grid_NcDef

#ifdef with_netcdf
    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF , only : NF90_NETCDF4
    use NetCDF , only : NF90_GLOBAL
    use NetCDF , only : NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT
#endif

    use LE_CF_Conventions   , only : LE_CF_names

    use Dims   , only : nx, ny, nz
    use LE_Data_Common  , only : nlev_top !change by JPT
    use LE_Data      , only : LE_Data_GetPointer

    use LE_Output_Tools , only : LE_Output_Define_Dims_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Define_Dims_Lev
    use LE_Output_Tools , only : LE_Output_Define_Dims_Time

    use LE_Output_Tools , only : LE_Output_Define_Vars_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Define_Vars_Lev
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time

    use LE_Output_Tools , only : LE_Output_Put_Var_Lon_Lat
    use LE_Output_Tools , only : LE_Output_Put_Var_Lev
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains
    
    ! variables:
    use LE_Radiation    , only : LE_Radiation_Calc
    use LE_Radiation    , only : swbands
    use LE_Radiation    , only : tau, extinction, asy, ssa
    use LE_Radiation    , only : angstrom, i_ang_modis, i_ang_aeronet, i_ang_polder

    ! --- in/out --------------------------------

    type(T_LE_Output_Rad), intent(inout)  ::  leo
    type(TDate), intent(in)               ::  t
    real, intent(in)                      ::  c(:,:,:,:)  ! (nx,ny,nz,nspec) -> can also be (nx,ny,nlev_top,nspec) NO! top layers are not in c 
    integer, intent(out)                  ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Rad_PutOut'

    ! --- local ---------------------------------

    integer               ::  time6(6)
    integer               ::  time
#ifdef with_netcdf  
    integer               ::  cmode
#endif
    integer               ::  idat, ilev, ilev1
    integer               ::  varid
    type(T_Grid_NcDef)    ::  gncd
    type(TDate)           ::  t0
    integer               ::  i, j, k, l
    integer               ::  nz_top, with_aloft !change by JPT to deal with top layers in radiation output
    real                  ::  pat(nx,ny)
    real, allocatable     ::  field(:,:,:) !(nx, ny, 0:nz+1) or (nx, ny, 0:nlev_top+1) (if with top layers) 
    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units
    character(len=512)    ::  comment
    character(len=32)     ::  selection_key
    
    real, allocatable       ::  hh(:,:,:) !(nx, ny, 0:nz+1) or (nx, ny, 0:nlev_top+1) (if with top layers)     
    !real                  ::  hh(nx,ny,0:nz+1)
    real                  ::  hsamp(nx,ny)
    integer               ::  ilast

    integer               ::  icountry
    integer               ::  iemis
    integer               ::  itr
    
    character(len=256) :: commandline

    ! meteo data:
    real, pointer        ::   oro(:,:,:)      ! (lon,lat,1)
    real, pointer        ::   h_m(:,:,:)      ! (lon,lat,lev)

    ! optical properties:
    integer               ::  vlen
    integer               ::  lambda
    integer               ::  iswband
    real                  ::  sumext(nx, ny), sumscat(nx,ny)

    ! --- begin ---------------------------------

#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif

    
    ! setting nz_top to be either nz or nlev_top if we have top layers in our simulated concentrations
    if (leo%nlev > nz) then !also top layer output
      nz_top = nlev_top
      with_aloft = 0 ! aloft is part of the top layers so some loops should not go to nz+1 any more
      ! info ...
      write (gol,'("we have top layers:")'); call goPr
      write(*,*) shape(c) 
      write(*,*) nz_top 
      allocate(field(nx,ny,nz_top))
    else
      nz_top = nz    
      with_aloft = 1
      write (gol,'("we do not have top layers:")'); call goPr
      write(*,*) shape(c)
      write(*,*) nz_top 
      allocate(field(nx,ny,nz_top))        
    end if
    
    ! for multiples of dhour only ...
    if ( .not. Precisely(t,leo%dhour,'hour')  ) then
      status=0; return
    end if

    ! point to meteo data:
    call LE_Data_GetPointer( 'oro', oro, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)

    ! current time not in time range ?
    if ( (t < leo%tr(1)) .or. (leo%tr(2) < t) ) then
    
      ! extract time fields:
      call Get( t, time6=time6 )

      ! daily or less ?
      select case ( trim(leo%collect) )
        ! collect daily for [00,24)
        case ( 'daily' )
          ! set time range [00,24) for this day:
          leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
          leo%tr(2) = leo%tr(1) + IncrDate( day=1 ) - IncrDate(mili=1)
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), time6(1:3)
          if ( len_trim(leo%state) > 0 ) write (leo%fname,'(a,"_",a)') trim(leo%fname), trim(leo%state)
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! collect daily for (00,24]
        case ( 'daily24' )
          ! set time range (00,24] for this day:
          if ( Midnight(t) ) then
            leo%tr(1) = t - IncrDate(day=1)
            leo%tr(2) = t
          else
            leo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
            leo%tr(2) = leo%tr(1) + IncrDate( day=1 )
            leo%tr(1) = leo%tr(1) + IncrDate(mili=1)
          end if
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), leo%tr(1)%year, leo%tr(1)%month, leo%tr(1)%day
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! files with instant fields:
        case ( 'instant' )
          ! set time range for current instant time:
          leo%tr(1) = t
          leo%tr(2) = t
          ! new file name:
          write (leo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",2i2.2)') &
                    trim(leo%com%outdir), trim(leo%com%model), trim(leo%com%expid), &
                    trim(leo%name), time6(1:5)
          write (leo%fname,'(a,".nc")') trim(leo%fname)
        ! unknonw ...
        case default
          write (gol,'("unsupported collect key : ",a)') trim(leo%collect); call goErr
          TRACEBACK; status=1; return
      end select

      ! root only:
      if ( goc%root ) then

#ifdef with_netcdf
        ! set creation mode flag:
        if ( leo%replace ) then
          cmode = NF90_CLOBBER       ! overwrite existing files
        else
          cmode = NF90_NOCLOBBER     ! do not overwrite existing files
        end if

        ! enable large file support:
        cmode = or( cmode, NF90_NETCDF4 )

        ! create file:
        status = NF90_Create( trim(leo%fname), cmode, leo%ncid )
        if ( status /= NF90_NOERR ) then
          write (gol,'("creating file : ")'); call goErr
          write (gol,'("  file name  : ",a)') trim(leo%fname); call goErr
          write (gol,'("  nf90 error : ",a)') trim(nf90_strerror(status)); call goErr
          TRACEBACK; status=1; return
        end if
#else
        ! dummmy ...
        leo%ncid = -1
#endif

        ! reset flag:
        leo%opened = .true.

        ! write global attributes:
        call leo%com%PutOut_GlobalAttributes( leo%ncid, status )
        IF_NOTOK_RETURN(status=1)

#ifdef with_netcdf
        !
        ! define dimensions:
        !

        ! grid dimensions/variables
        call glb_ugg%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                    dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )
                                    !subset=(/leo%i1,leo%i2,leo%j1,leo%j2/) )
        IF_NOTOK_RETURN(status=1)

        ! level dimension
        call LE_Output_Define_Dims_Lev(leo%ncid, leo%dimid_lev, leo%nlev, trim(leo%com%CF_convention), status)
        IF_NOTOK_RETURN(status=1)

        ! time dimensions
        call LE_Output_Define_Dims_Time(leo%ncid, leo%dimid_time, status)
        IF_NOTOK_RETURN(status=1)

        !
        ! define variables:
        !
        ! level variables
        call LE_Output_Define_Vars_Lev(leo%ncid, leo%varid_lev,leo%dimid_lev, &
                                         trim(leo%levtype), trim(leo%com%CF_convention), status)
        IF_NOTOK_RETURN(status=1)

        ! time since t0
        t0 = leo%com%t0      
        ! time variables
        call LE_Output_Define_Vars_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                        leo%dimid_time, trim(leo%com%CF_convention), t0, status)
        IF_NOTOK_RETURN(status=1)
#endif

      end if  ! root

      ! loop over data fields to be written:
      do l = 1, leo%ndat

        ! global tracer index
        idat = leo%idat(l)

        ! CF standard name for concentration/mixing ratio/column:

        ! initial comment:
        select case ( trim(leo%LE_Dat(idat)%name) )
          case ( 'mpressure' )
            comment = 'pressure at mid of layer (full level)'
          case ( 'hpressure' )
            comment = 'pressure at top of layer (upper half level)'
          case default
            comment = ''
        end select

        ! get names following CF conventions;
        ! store conversion factor for later usage:
        call LE_CF_names( &
                     leo%LE_Dat(idat)%name, leo%LE_Dat(idat)%unit, &
                     cf_standard_name, cf_long_name, cf_units, &
                     leo%unitconv(l), comment, &
                     status )
        IF_NOTOK_RETURN(status=1)

        ! store units for later usage (GrADS ctl file):
        leo%unit_dat(l) = trim(cf_units)

        ! root only:
        if ( goc%root ) then

          ! define variable:
          select case ( leo%LE_Dat(idat)%rank )
            case ( 2 )
              if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf  
                status = NF90_Def_Var( leo%ncid, trim(leo%name_dat(l)), NF90_REAL, &
                                         (/leo%dimid_lon,leo%dimid_lat/), varid, &
                                         deflate_level=leo%com%deflate_level )
                IF_NF90_NOTOK_RETURN(status=1)
                status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'latitude longitude')
                IF_NF90_NOTOK_RETURN(status=1)
#endif
              else
#ifdef with_netcdf  
                status = NF90_Def_Var( leo%ncid, trim(leo%name_dat(l)), NF90_REAL, &
                                         (/leo%dimid_lon,leo%dimid_lat,leo%dimid_time/), varid, &
                                         deflate_level=leo%com%deflate_level )
                IF_NF90_NOTOK_RETURN(status=1)
                status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time latitude longitude')
                IF_NF90_NOTOK_RETURN(status=1)
#endif
              end if
            case ( 3 )
              if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf  
                status = NF90_Def_Var( leo%ncid, trim(leo%name_dat(l)), NF90_REAL, &
                                         (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev/), varid, &
                                         deflate_level=leo%com%deflate_level )
                IF_NF90_NOTOK_RETURN(status=1)
                status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'level latitude longitude')
                IF_NF90_NOTOK_RETURN(status=1)
#endif
              else
#ifdef with_netcdf  
                status = NF90_Def_Var( leo%ncid, trim(leo%name_dat(l)), NF90_REAL, &
                                         (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), varid, &
                                         deflate_level=leo%com%deflate_level )
                IF_NF90_NOTOK_RETURN(status=1)
                status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time level latitude longitude')
                IF_NF90_NOTOK_RETURN(status=1)
#endif
              end if
            case default
              write (gol,'("unsupported data rank : ",i4)') leo%LE_Dat(idat)%rank; call goErr
              TRACEBACK; status=1; return
          end select

          ! write attributes:
#ifdef with_netcdf  
          if ( len_trim(cf_standard_name) > 0 ) then
            status = nf90_put_att( leo%ncid, varid, 'standard_name', trim(cf_standard_name) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
          status = nf90_put_att( leo%ncid, varid, 'long_name', trim(cf_long_name) )
          IF_NF90_NOTOK_RETURN(status=1)
          ! write units:
          status = nf90_put_att( leo%ncid, varid, 'units', trim(cf_units) )
          IF_NF90_NOTOK_RETURN(status=1)
          call glb_ugg%DefCoor_NetCDF( gncd, varid, status )
          IF_NOTOK_RETURN(status=1)
          ! specials ...
          select case ( trim(leo%LE_Dat(idat)%name) )
            case ( 'dens' )
              ! moleweight of air:
              status = nf90_put_att( leo%ncid, varid, 'moleweight_air', xm_air )
              IF_NF90_NOTOK_RETURN(status=1)
              status = nf90_put_att( leo%ncid, varid, 'moleweight_air_unit', 'kg mole-1' )
              IF_NF90_NOTOK_RETURN(status=1)
          end select
          ! add comment:
          if ( len_trim(comment) > 0 ) then
            status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
#endif

          ! store variable id:
          leo%varid_dat(l) = varid
          
        end if  ! root

      end do  ! written tracers

      ! root only:
      if ( goc%root ) then

        ! end defintion mode:
#ifdef with_netcdf  
        status = NF90_EndDef( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif

      end if  ! root

      ! no records written yet:
      leo%itrec = 0

    end if

    ! next time record:
    leo%itrec = leo%itrec + 1

    ! root only:
    if ( goc%root ) then
    
      ! GrADS time counter:
      if ( leo%grads_ctl ) then
        ! increase counter:
        leo%grads_ctl_nt = leo%grads_ctl_nt + 1
        ! set times if necessary:
        if ( leo%grads_ctl_nt == 1 ) then
          leo%grads_ctl_t1 = t
          leo%grads_ctl_dt = IncrDate(day=1)   ! dummy ...
        end if
        if ( leo%grads_ctl_nt == 2 ) then
          leo%grads_ctl_dt = t - leo%grads_ctl_t1
        end if
      end if

      ! write dimension data only once ...
      if ( leo%itrec == 1 ) then

        ! write grid to netCDF file
        call glb_ugg%PutGrid_NetCDF( gncd, status )
        IF_NOTOK_RETURN(status=1)
#ifdef with_netcdf
        ! write level indices
        call LE_Output_Put_Var_Lev(leo%ncid, leo%varid_lev, leo%nlev, &
                                   trim(leo%levtype), leo%ilev, leo%heights, status)
#endif

      end if  ! first record

      ! date up to seconds:
      call Get( t, time6=time6 )

      ! time since t0
      t0 = leo%com%t0          
      time = iTotal( t - t0, 'sec' )

      ! write time record:
#ifdef with_netcdf
      call LE_Output_Put_Var_Time(leo%ncid, leo%varid_time, leo%varid_time_dtg, &
                                   time, time6, trim(leo%com%CF_convention), leo%itrec, status )
      IF_NOTOK_RETURN(status=1)
#endif

    end if  ! root

    ! sample ?    
    if ( (trim(leo%levtype) == 'heights'   ) .or. &
         (trim(leo%levtype) == 'elevations') ) then
      ! lowest is orography:
      hh(:,:,0) = oro(:,:,1)  ! m
      ! mid of layers:
      hh(:,:,1) = oro(:,:,1) + h_m(:,:,1) * 0.5  ! m
      do k = 2, nz_top
        hh(:,:,k) = oro(:,:,1) + ( h_m(:,:,k-1) + h_m(:,:,k) ) * 0.5  ! m
      end do
      ! dummy for aloft ... JPT: not a problem that this becomes too big with top layers
      hh(:,:,nz_top+1) = h_m(:,:,nz_top) + 1000.0   ! 1000 m above top
    end if

    ! loop over all tracer to be written:
    do l = 1, leo%ndat

      ! global tracer index:
      idat = leo%idat(l)

      ! constant fields written only once:
      if ( leo%LE_Dat(idat)%const .and. (leo%itrec > 1) ) cycle

      ! 2d or 3d ?
      select case ( leo%LE_Dat(idat)%rank )
      
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 2 )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! extract 2d field:
          selection_key =  trim(leo%LE_Dat(idat)%name) 
            
          if (selection_key == 'angstrom_aeronet') then
              pat = angstrom(:,:,i_ang_aeronet )
          else if ( selection_key == 'angstrom_modis' ) then
              pat = angstrom(:,:,i_ang_modis )
          else if (selection_key == 'angstrom_polder') then
               pat = angstrom(:,:,i_ang_polder ) 
            
          else if (any(aod_keylist == selection_key)) then
              ! extract wavelength:
              vlen = len_trim(leo%LE_Dat(idat)%name)
              read (leo%LE_dat(idat)%name(5:vlen-2),*) lambda
              ! band index:
              call swbands%FindBand( lambda*0.001, iswband, status )
              IF_NOTOK_RETURN(status=1)
              ! sum:
            pat = sum( tau(:,:,1:nz_top,iswband), dim=3 )

          else if (any(ssa_keylist == selection_key)) then
              ! extract wavelength:
              vlen = len_trim(leo%LE_Dat(idat)%name)
              read (leo%LE_dat(idat)%name(5:vlen-2),*) lambda
              ! band index:
              call swbands%FindBand( lambda*0.001, iswband, status )
              IF_NOTOK_RETURN(status=1)
              ! weighted sum:   ssa = sum ( ssa_k tau_k ) / sum( tau_k )
            sumext = sum( tau(:,:,1:nz_top,iswband), dim=3 )
              where ( sumext > 0.0 )
              pat = sum( ssa(:,:,1:nz_top,iswband) * tau(:,:,1:nz_top,iswband), dim=3 ) / sumext
              endwhere
              
          else
              write (gol,'("unsupported 2d data : ",a)') trim(leo%LE_Dat(idat)%name); call goErr
              TRACEBACK; status=1; return
          end if
          ! unit conversion:
          pat = pat * leo%unitconv(l)
          ! write data:
          if ( leo%LE_Dat(idat)%const ) then
            ! write 2D field, no level, no time index:
            call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), -1, -1, &
                                              pat, status )
            IF_NOTOK_RETURN(status=1)
          else
            ! write 2D field, no level, with time index:
            call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), -1, leo%itrec, &
                                              pat, status )
            IF_NOTOK_RETURN(status=1)
          end if

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 3 )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! loop over all model levels, including surface and aloft:
          do ilev = 0, nz_top+with_aloft
            ! corresponding model layer, extract surface (0) from layer 1:
            ilev1 = max(ilev,1)
            ! set to zero by default, might be used for ilev==nz+1
            pat = 0.0
            ! extract 2d field:
            !select case ( trim(leo%LE_Dat(idat)%name) )
            selection_key =  trim(leo%LE_Dat(idat)%name) 

			          ! extinction coefficient
            if (any(extinction_keylist == selection_key)) then

              ! extract wavelength:
              vlen = len_trim(leo%LE_Dat(idat)%name)
              read (leo%LE_dat(idat)%name(12:vlen-2),*) lambda
              ! band index:
              call swbands%FindBand( lambda*0.001, iswband, status )
              IF_NOTOK_RETURN(status=1)
              ! sum:
              if ( ilev <= nz_top ) then
                  pat = extinction(:,:,ilev1,iswband)
                else
                  pat = extinction(:,:,nz_top, iswband)     ! <-- extinction top layer
              end if
              
           
			        ! aerosol optical depth of grid cell):
            else if (any(tau_keylist == selection_key)) then

              ! extract wavelength:
              vlen = len_trim(leo%LE_Dat(idat)%name)
              read (leo%LE_dat(idat)%name(5:vlen-2),*) lambda
              ! band index:
              call swbands%FindBand( lambda*0.001, iswband, status )
              IF_NOTOK_RETURN(status=1)
              
                if ( ilev <= nz_top ) then
                  pat = tau(:,:,ilev1,iswband) ! ~~visible band, 550 nm
                else
                  pat = tau(:,:,nz_top, iswband)     ! <-- aod top layer
                end if
				                
                ! single scattering albedo of grid cell):
            else if (any(ssa_keylist == selection_key)) then

              ! extract wavelength:
              vlen = len_trim(leo%LE_Dat(idat)%name)
              read (leo%LE_dat(idat)%name(5:vlen-2),*) lambda
              ! band index:
              call swbands%FindBand( lambda*0.001, iswband, status )
              IF_NOTOK_RETURN(status=1)
              
                if ( ilev <= nz_top ) then
                  pat = ssa(:,:,ilev1,iswband)
                else
                  pat = ssa(:,:,nz_top, iswband)     ! <-- ssa top layer
                end if
                
                ! asymmetry parameter of grid cell):
            else if (any(asy_keylist == selection_key)) then
            
              vlen = len_trim(leo%LE_Dat(idat)%name)
              read (leo%LE_dat(idat)%name(5:vlen-2),*) lambda
              ! band index:
              call swbands%FindBand( lambda*0.001, iswband, status )
              IF_NOTOK_RETURN(status=1)
              
                !print *,'write asy'
                if ( ilev <= nz_top ) then
                  pat = asy(:,:,ilev1,iswband)
                else
                  pat = asy(:,:,nz_top,iswband)     ! <-- asy top layer
                end if

              else
                write (gol,'("unsupported 3d data : ",a)') trim(leo%LE_Dat(idat)%name); call goErr
                TRACEBACK; status=1; return
              end if
            ! unit conversion:
            pat = pat * leo%unitconv(l)
            
            write(*,*) shape(pat) 
            write(*,*) shape(field(:,:,ilev)) 
            ! store:
            field(:,:,ilev1) = pat
            
          end do  ! model layers

          ! * vertical selection or interpolation:
          
          ! which output levels ?
          select case ( trim(leo%levtype) )

            !~ selected model levels:
            case ( 'levels' )

              ! loop over selected layer:
              do k = 1, leo%nlev
                ! global level index:
                ilev = leo%ilev(k)
                ! extract 2D field:
                pat = field(:,:,ilev)
                ! write data:
                if ( leo%LE_Dat(idat)%const ) then
                  ! write 2D field for current level, no time index:
                  call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), k, -1, &
                                                    pat, status )
                  IF_NOTOK_RETURN(status=1)
                else
                  ! write 2D field for current level, with time index:
                  call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), k, leo%itrec, &
                                                    pat, status )
                  IF_NOTOK_RETURN(status=1)
                end if
              end do   ! selected model layers

            !~ height levels:
            case ( 'heights', 'elevations' )

              ! loop over target levels:
              do k = 1, leo%nlev
                ! sample height:
                if ( trim(leo%levtype) == 'heights' ) then
                  ! relative to orography:
                  hsamp = hh(:,:,0) + leo%heights(k)
                else if ( trim(leo%levtype) == 'elevations' ) then
                  ! absolute, interpolation will take surface value if below orography:
                  hsamp = leo%heights(k)
                else
                  write (gol,'("no sample height defined for level type : ",a)') trim(leo%levtype); call goErr
                  TRACEBACK; status=1; return
                end if
                ! loop over horizontal cells:
                do j = 1, ny
                  do i = 1, nx
                    ! vertical interpolation:
                    call Interp_Lin( hh(i,j,:), field(i,j,:), hsamp(i,j), pat(i,j), ilast, status )
                    IF_NOTOK_RETURN(status=1)
                  end do  ! i
                end do  ! j
                ! write data:
                if ( leo%LE_Dat(idat)%const ) then
                  ! write 2D field for current level, no time index:
                  call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), k, -1, &
                                                    pat, status )
                  IF_NOTOK_RETURN(status=1)
                else
                  ! write 2D field for current level, with time index:
                  call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), k, leo%itrec, &
                                                    pat, status )
                  IF_NOTOK_RETURN(status=1)
                end if
              end do   ! heights

            !~ unknown ...
            case default
              write (gol,'("unsupported level type : ",a)') trim(leo%levtype); call goErr
              TRACEBACK; status=1; return

          end select

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          write (gol,'("unsupported data rank : ",i4)') leo%LE_Dat(idat)%rank; call goErr
          TRACEBACK; status=1; return

      end select

    end do   ! variables

    ! root?
    if ( goc%root ) then

      ! next time exceeds interval ?
      if ( t+IncrDate(hour=int(leo%dhour)) > leo%tr(2) ) then
        ! close:
#ifdef with_netcdf
        status = NF90_Close( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif
        ! reset flag:
        leo%opened = .false.
        ! write GrADS ctl file if necessary:
        call Write_GrADS_Ctl( leo, status )
        IF_NOTOK_RETURN(status=1)
      end if
      
    end if ! root

    ! ok
    status = 0

  end subroutine LE_Output_Rad_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims, only : nx, ny
    use GrADS_Ctl
    use LE_Grid, only : glb_ugg

    ! --- in/out ---------------------------------

    type(T_LE_Output_Rad), intent(inout)    ::  leo
    integer, intent(out)                    ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/Write_GrADS_Ctl'

    ! --- local ----------------------------------

    type(T_GrADS_Ctl)     ::  ctl
    integer               ::  l, idat
    character(len=512)    ::  line

    ! --- begin ----------------------------------

    ! write ctl file ?
    if ( leo%grads_ctl ) then

      ! ctl file name:
      write (leo%grads_ctl_file,'(a,"_",a,"_",a)') &
                trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
      if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_file,'(a,"_",a)') trim(leo%grads_ctl_file), trim(leo%state)
      write (leo%grads_ctl_file,'(a,".ctl")') trim(leo%grads_ctl_file)

      ! daily or less ?
      select case ( trim(leo%collect) )
        ! collect daily for [00,24)
        case ( 'daily', 'daily24' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_base,'(a,"_",a)') trim(leo%grads_ctl_base), trim(leo%state)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! files with instant fields:
        case ( 'instant' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2_%h2%n2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          if ( len_trim(leo%state) > 0 ) write (leo%grads_ctl_base,'(a,"_",a)') trim(leo%grads_ctl_base), trim(leo%state)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! unknonw ...
        case default
          write (gol,'("unsupported collect key : ",a)') trim(leo%collect); call goErr
          TRACEBACK; status=1; return
      end select

      ! open ctl file:
      call GrADS_Ctl_Init( ctl, trim(leo%com%outdir), trim(leo%grads_ctl_file), status )
      IF_NOTOK_RETURN(status=1)
      ! comment ...
      call GrADS_Ctl_Comment( ctl, '', status )
      call GrADS_Ctl_Comment( ctl, 'GrADS Data Descriptor File', status )
      call GrADS_Ctl_Comment( ctl, '', status )
      ! data set:
      call GrADS_Ctl_DSet( ctl, trim(leo%grads_ctl_base), status )
      IF_NOTOK_RETURN(status=1)
      ! title:
      write (line,'("model: ",a,"; expid: ",a)') trim(leo%com%model), trim(leo%com%expid)
      call GrADS_Ctl_Title( ctl, trim(line), status )
      IF_NOTOK_RETURN(status=1)
      ! write xdef/ydef from grid/projection definition
      call glb_ugg%WriteCtlProjection( ctl, status)
      IF_NOTOK_RETURN(status=1)
      ! zdef:
      select case ( trim(leo%levtype) )
        case ( 'levels' )
          call GrADS_Ctl_ZDef( ctl, leo%ilev(1:leo%nlev), status )
          IF_NOTOK_RETURN(status=1)
        case ( 'heights', 'elevations' )
          call GrADS_Ctl_ZDef( ctl, nint(leo%heights(1:leo%nlev)), status )
          IF_NOTOK_RETURN(status=1)
        case default
          write (gol,'("do not know how to write zdef for levtyp : ",a)') trim(leo%levtype); call goPr
          TRACEBACK; status=1; return
      end select
      ! tdef:
      call GrADS_Ctl_TDef( ctl, leo%grads_ctl_nt, &
                               (/leo%grads_ctl_t1%year,leo%grads_ctl_t1%month,leo%grads_ctl_t1%day,leo%grads_ctl_t1%hour,leo%grads_ctl_t1%min/), &
                               (/                    0,                     0,leo%grads_ctl_dt%day,leo%grads_ctl_dt%hour,leo%grads_ctl_dt%min/), &
                               status )
      IF_NOTOK_RETURN(status=1)
      ! number of variables:
      call GrADS_Ctl_Vars( ctl, leo%ndat, status )
      IF_NOTOK_RETURN(status=1)
      ! loop over data fields to be written:
      do l = 1, leo%ndat
        ! global tracer index
        idat = leo%idat(l)
        ! set variable lineiption:
        write (line,'(a," [",a,"]")') trim(leo%name_dat(l)), trim(leo%unit_dat(l))
        ! define variable:
        select case ( leo%LE_Dat(idat)%rank )
          case ( 2 )
            if ( leo%LE_Dat(idat)%const ) then
              call GrADS_Ctl_Var( ctl, trim(leo%name_dat(l)), 1, 'y,x', trim(line), status )
              IF_NOTOK_RETURN(status=1)
            else
              call GrADS_Ctl_Var( ctl, trim(leo%name_dat(l)), 1, 't,y,x', trim(line), status )
              IF_NOTOK_RETURN(status=1)
            end if
          case ( 3 )
            if ( leo%LE_Dat(idat)%const ) then
              call GrADS_Ctl_Var( ctl, trim(leo%name_dat(l)), leo%nlev, 'z,y,x', trim(line), status )
              IF_NOTOK_RETURN(status=1)
            else
              call GrADS_Ctl_Var( ctl, trim(leo%name_dat(l)), leo%nlev, 't,z,y,x', trim(line), status )
              IF_NOTOK_RETURN(status=1)
            end if
          case default
            write (gol,'("unsupported data rank : ",i4)') leo%LE_Dat(idat)%rank; call goErr
            TRACEBACK; status=1; return
        end select
      end do
      ! end of variables section:
      call GrADS_Ctl_EndVars( ctl, status )
      IF_NOTOK_RETURN(status=1)
      ! close ctl file:
      call GrADS_Ctl_Done( ctl, status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! ok
    status = 0

  end subroutine Write_GrADS_Ctl



end module LE_Output_Rad
