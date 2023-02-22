!###############################################################################
!
! NAME
!
!   LE_Output_dat  -  LOTOS-EUROS output of 3D fields
!
! HISTORY
!
!   2007 may, Arjo Segers, TNO
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

module LE_Output_Dat2

  use GO, only : gol, goPr, goErr
  use GO, only : TDate, TIncrDate

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  use LE_Output_Common, only : T_LE_Output_Common

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output_Dat2

  public  ::  LE_Output_Dat2_Init, LE_Output_Dat2_Done
  public  ::  LE_Output_Dat2_PutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output_Dat2'

  ! maximum number of levels in output:
  integer, parameter  ::  maxlev = 10


  ! --- types ------------------------------

  type T_LE_Dat
    character(len=32)      ::  name
    character(len=256)     ::  long_name
    character(len=32)      ::  unit
    integer                ::  rank
    logical                ::  const
  end type T_LE_Dat

  type T_LE_Output_Dat2
    ! name for this file:
    character(len=16)                 ::  typ
    character(len=16)                 ::  name
    ! common stuff:
    type(T_LE_Output_Common)          ::  com
    ! replace existing files ?
    logical                           ::  replace
    ! file opened ?
    logical                           ::  opened
    ! current time range:
    type(TDate)                       ::  tr(2)
    ! time resolution:
    real                              ::  dhour
    ! collect: daily, instant
    character(len=32)                 ::  collect
    ! time record counter:
    integer                           ::  itrec
    ! file name:
    character(len=256)                ::  fname
    ! file handle:
    integer                           ::  ncid
    ! dimension handles:
    integer                           ::  dimid_lon
    integer                           ::  dimid_lat
    integer                           ::  dimid_lev
    integer                           ::  dimid_time
    integer                           ::  dimid_datelen
    ! dimension variables:
    integer                           ::  varid_lon
    integer                           ::  varid_lat
    integer                           ::  varid_lev
    integer                           ::  varid_time
    !integer                           ::  varid_time_day
    integer                           ::  varid_date
    integer                           ::  varid_time_dtg
    ! database with supported variables:
    type(T_LE_Dat), allocatable       ::  LE_Dat(:)
    ! tracer variables:
    integer                           ::  ndat
    integer, allocatable              ::  idat(:)       ! (ndat)
    character(len=1),  allocatable    ::  datcomp(:)
    character(len=32),  allocatable   ::  name_dat(:)
    character(len=32), allocatable    ::  unit_dat(:)
    !real, allocatable                 ::  unitconv(:)
    integer, allocatable              ::  varid_dat(:)
    ! level selection:
    character(len=16)                 ::  levtype
    integer                           ::  nlev
    integer, allocatable              ::  ilev(:)
    real                              ::  heights(maxlev)
    ! grads ctl file ?
    logical                           ::  grads_ctl
    character(len=256)                ::  grads_ctl_file
    character(len=256)                ::  grads_ctl_base
    integer                           ::  grads_ctl_nt
    type(TDate)                       ::  grads_ctl_t1
    type(TIncrDate)                   ::  grads_ctl_dt
    ! grid type:
    character(len=16)                 ::  gridtype
    ! bounding box
    integer                           ::  i1, i2, ni
    integer                           ::  j1, j2, nj
    real                              ::  westb, southb
  end type T_LE_Output_Dat2


  ! --- var --------------------------------------



contains


  ! ====================================================


  subroutine SetDat( d, name, long_name, unit, rank, const )

    ! --- in/out ----------------------------------

    type(T_LE_Dat), intent(out)     ::  d
    character(len=*), intent(in)    ::  name
    character(len=*), intent(in)    ::  long_name
    character(len=*), intent(in)    ::  unit
    integer, intent(in)             ::  rank
    logical, intent(in)             ::  const

    ! --- begin ----------------------------------

    d%name = name
    d%long_name = long_name
    d%unit = unit
    d%rank = rank
    d%const = const

  end subroutine SetDat


  ! ====================================================


  subroutine LE_Output_Dat2_Init( leo, rcF, rckey, typ, name, status )

    use GO     , only : TrcFile
    use GO     , only : goMatchValues, goSplitString, goReadFromLine
    use GO     , only : AnyDate
    use Dims   , only : nx, ny, nz
    
    use LE_Grid, only : ugg
    use LE_Output_Common, only : Init

    use LE_Data, only : LE_Data_Inquire
    use LE_Data, only : LE_Data_InqVar
    use LE_Data, only : LE_Data_Enable

    ! --- in/out --------------------------------

    type(T_LE_Output_Dat2), intent(out)   ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    character(len=*), intent(in)          ::  typ
    character(len=*), intent(in)          ::  name
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Dat2_Init'

    ! rckey extensions:

    ! --- local ---------------------------------

    character(len=64)     ::  basekey
    character(len=256)    ::  field_names
    integer               ::  idat
    character(len=32)     ::  level_names
    integer               ::  ilev
    integer               ::  i, i1, i2
    integer               ::  k
    
    integer               ::  nvar
    character(len=64)     ::  var_name, var_units
    character(len=256)    ::  var_long_name
    integer               ::  var_rank
    logical               ::  var_const
    character(len=1)      ::  var_comp
    character(len=64)     ::  var_levtype

    character(len=512)    ::  key
    real                  ::  west, south, east, north
    real, pointer         ::  ff(:,:)

    integer               ::  imf

    ! --- begin ---------------------------------

    ! store name:
    leo%typ  = typ
    leo%name = name

    ! init common stuff:
    call Init( leo%com, rcF, rckey, status )
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

    ! grid type:
    call rcF%Get( trim(basekey)//'.gridtype', leo%gridtype, status )
    IF_NOTOK_RETURN(status=1)

    ! level type::
    call rcF%Get( trim(basekey)//'.levtype', leo%levtype, status )
    IF_NOTOK_RETURN(status=1)

    ! define supported data:
    !   name   :  used in rcfile to identify this field;
    !   unit   :  in which the corresponding field in LE is expressed;
    !             the model units are converted to 'cf'-conventions
    !             (SI-units) on output;
    !   rank   :  2 for 2D field, and 3 for 3D
    !   const  :  logical, T for constant fields
    !
    ! number of variables:
    call LE_Data_Inquire( status, nvar=nvar )
    IF_NOTOK_RETURN(status=1)
    ! storage:
    allocate( leo%LE_Dat(nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! loop:
    do i = 1, nvar
      ! properties for this variable:
      call LE_Data_InqVar( i, status, name=var_name, long_name=var_long_name, &
                                units=var_units, rank=var_rank, const=var_const )
      IF_NOTOK_RETURN(status=1)
      ! define output properties:
      call SetDat( leo%LE_Dat(i), trim(var_name), trim(var_long_name), &
                                  trim(var_units), var_rank, var_const  )
    end do ! ivar
!    i = 0
!    ! switch:
!    select case ( trim(leo%gridtype) )
!      !
!      ! cells:
!      case ( 'cells' )
!        ! constant sfc
!        do imf = lbound(cmf_2d,1), ubound(cmf_2d,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(cmf_2d(imf)%name), trim(cmf_2d(imf)%units), 2, .true.  )
!        end do
!        ! temporal sfc
!        do imf = lbound(tmf_2d,1), ubound(tmf_2d,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(tmf_2d(imf)%name), trim(tmf_2d(imf)%units), 2, .false. )
!        end do
!        !! temporal 3d
!        !do imf = lbound(tmf_3d,1), ubound(tmf_3d,1)
!        !  i = i + 1
!        !  call SetDat( leo%LE_Dat(i), trim(tmf_3d(imf)%name), trim(tmf_3d(imf)%units), 3, .false. )
!        !end do
!        ! derived sfc
!        do imf = lbound(cf,1), ubound(cf,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(cf(imf)%name), trim(cf(imf)%units), 2, .false. )
!        end do
!        ! derived 3d
!        do imf = lbound(dmf_3d,1), ubound(dmf_3d,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(dmf_3d(imf)%name), trim(dmf_3d(imf)%units), 3, .false. )
!        end do
!      !
!      ! corners:
!      case ( 'corners' )
!        ! constant sfc
!        do imf = lbound(cmf_2dc,1), ubound(cmf_2dc,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(cmf_2dc(imf)%name), trim(cmf_2dc(imf)%units), 2, .true.  )
!        end do
!        ! temporal 3d at corners
!        do imf = lbound(tmf_3dc,1), ubound(tmf_3dc,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(tmf_3dc(imf)%name), trim(tmf_3dc(imf)%units), 3, .false. )
!        end do
!        ! derived 3d at corners
!        do imf = lbound(dmf_3dc,1), ubound(dmf_3dc,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(dmf_3dc(imf)%name), trim(dmf_3dc(imf)%units), 3, .false. )
!        end do
!      !
!      ! u-edges:
!      case ( 'u-edges' )
!        ! constant 2d at u-edges (some fields)
!        do imf = lbound(cf,1), ubound(cf,1)
!          select case ( trim(cf(imf)%name) )
!            case ( 'udy' )
!              i = i + 1
!              call SetDat( leo%LE_Dat(i), trim(cf(imf)%name), trim(cf(imf)%units), 2, .true.  )
!          end select
!        end do
!        ! derived 3d at u-edges
!        do imf = lbound(dmf_3du,1), ubound(dmf_3du,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(dmf_3du(imf)%name), trim(dmf_3du(imf)%units), 3, .false.  )
!        end do
!      !
!      ! v-edges:
!      case ( 'v-edges' )
!        ! constant 2d at v-edges (some fields)
!        do imf = lbound(cf,1), ubound(cf,1)
!          select case ( trim(cf(imf)%name) )
!            case ( 'vdx' )
!              i = i + 1
!              call SetDat( leo%LE_Dat(i), trim(cf(imf)%name), trim(cf(imf)%units), 2, .true.  )
!          end select
!        end do
!        ! derived 3d at v-edges
!        do imf = lbound(dmf_3dv,1), ubound(dmf_3dv,1)
!          i = i + 1
!          call SetDat( leo%LE_Dat(i), trim(dmf_3dv(imf)%name), trim(dmf_3dv(imf)%units), 3, .false.  )
!        end do
!      !
!      case default
!        write (gol,'("unsupported gridtype `",a,"`")') trim(leo%gridtype); call goErr
!        TRACEBACK; status=1; return
!    end select

    ! setup storage for output fields:
    allocate( leo%idat     (nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( leo%datcomp  (nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( leo%name_dat (nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( leo%unit_dat (nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    !allocate( leo%unitconv (nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( leo%varid_dat(nvar), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! match names, skip vector component extensions after dot:
    call goMatchValues( field_names, leo%LE_Dat(:)%name, &
                          leo%ndat, leo%name_dat, leo%idat, &
                          status, skipext='.' )
    IF_NOTOK_RETURN(status=1)
    ! info ...
    write (gol,'("selected fields for data output:")'); call goPr
    do i = 1, leo%ndat
      idat = leo%idat(i)
      write (gol,'("  ",i3,"  ",a10,"  (",i3, a20," ",a10,")")') &
                  i, leo%name_dat(i), &
                  idat, trim(leo%LE_Dat(idat)%name), trim(leo%LE_Dat(idat)%unit); call goPr
    end do

    ! enable variables requested as output:
    do i = 1, leo%ndat
      ! variable name:
      var_name = trim(leo%name_dat(i))
      var_comp = '-'
      ! name could be component of vector:  uv.u
      k = index(trim(var_name),'.')
      if ( k > 0 ) then
        var_comp = var_name(k+1:k+1)
        var_name = var_name(1:k-1)
      end if
      if ( len_trim(var_comp) == 0 ) stop
      ! enable:
      call LE_Data_Enable( trim(var_name), status )
      IF_NOTOK_RETURN(status=1)
      ! store component part:
      leo%datcomp(i) = var_comp
    end do    

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
      select case ( trim(leo%gridtype) )
        case ( 'cells' )
          leo%i1 = 1
          leo%i2 = nx
          leo%ni = nx
          leo%j1 = 1
          leo%j2 = ny
          leo%nj = ny
          leo%westb  = ugg%longitude_bnds(1,1,1)
          leo%southb = ugg%latitude_bnds (1,1,1)
        case ( 'corners' )
          leo%i1 = 0
          leo%i2 = nx
          leo%ni = nx+1
          leo%j1 = 0
          leo%j2 = ny
          leo%nj = ny+1
          leo%westb  = ugg%longitude_bnds(1,1,1)
          leo%southb = ugg%latitude_bnds (1,1,1)
        case ( 'u-edges' )
          leo%i1 = 0
          leo%i2 = nx
          leo%ni = nx+1
          leo%j1 = 1
          leo%j2 = ny
          leo%nj = ny
          leo%westb  = ugg%longitude_bnds(1,1,1)
          leo%southb = ugg%latitude_bnds (1,1,1)
        case ( 'v-edges' )
          leo%i1 = 1
          leo%i2 = nx
          leo%ni = nx
          leo%j1 = 0
          leo%j2 = ny
          leo%nj = ny+1
          leo%westb  = ugg%longitude_bnds(1,1,1)
          leo%southb = ugg%latitude_bnds (1,1,1)
        case default
          write (gol,'("unsupported gridtype `",a,"`")') trim(leo%gridtype); call goErr
          TRACEBACK; status=1; return
      end select
    else
      ! not yet ...
      write (gol,'("no output subset supported for domain decomposition yet")'); call goErr
      TRACEBACK; status=1; return
      ! subset
      select case ( trim(leo%gridtype) )
        case ( 'cells' )
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

    ! * levels

    ! level descriptions:
    call rcF%Get( trim(basekey)//'.levels', level_names, status )
    IF_NOTOK_RETURN(status=1)
    ! layers or heights ?
    select case ( trim(leo%levtype) )
      !
      !~ model levels:
      case ( 'levels' )
        ! setup storage for level indices (surface + levels + upper boundary)
        !allocate( leo%ilev(1+nz+1) )
        allocate( leo%ilev(1:100) )
        ! switch:
        if ( trim(level_names) == 'all' ) then
          ! all:
          leo%nlev = nz
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
        write (gol,'("selected levels for conc output:")'); call goPr
        do i = 1, leo%nlev
          ilev = leo%ilev(i)
          write (gol,'("  ",i3,"  ",i3)') i, ilev; call goPr
        end do

      case ( 'halflevels' )
        ! setup storage for level indices (surface + levels + upper boundary)
        !allocate( leo%ilev(1+nz+1) )
        allocate( leo%ilev(1:100) )
        ! switch:
        if ( trim(level_names) == 'all' ) then
          ! all:
          leo%nlev = nz+1
          do i = 1, leo%nlev
            leo%ilev(i) = i-1
          end do
        else if ( index(trim(level_names),':') > 0 ) then
          ! extract range:
          call goReadFromLine( level_names, i1, status, sep=':' )
          IF_NOTOK_RETURN(status=1)
          call goReadFromLine( level_names, i2, status )
          IF_NOTOK_RETURN(status=1)
          ! count:
          leo%nlev = i2 - i1 + 2
          ! store indices:
          do i = 0, leo%nlev-1
            leo%ilev(i) = i1 + i
          end do
        else
          ! set selected level indices:
          call goMatchValues( level_names, 0, nz+1, leo%nlev, leo%ilev, status )
          IF_NOTOK_RETURN(status=1)
        end if
        ! info ...
        write (gol,'("selected levels for conc output:")'); call goPr
        do i = 1, leo%nlev
          ilev = leo%ilev(i)
          write (gol,'("  ",i3,"  ",i3)') i, ilev; call goPr
        end do

      !~ model levels:
      case ( 'input_levels', 'input_halflevels' )

        ! loop over requested as output:
        do i = 1, leo%ndat
          ! variable name:
          var_name = trim(leo%name_dat(i))
          var_comp = '-'
          ! name could be component of vector:  uv.u
          k = index(trim(var_name),'.')
          if ( k > 0 ) then
            var_comp = var_name(k+1:k+1)
            var_name = var_name(1:k-1)
          end if
          if ( len_trim(var_comp) == 0 ) stop

          ! properties for this variable:
          call LE_Data_InqVar( trim(var_name), status, levtype=var_levtype )
          IF_NOTOK_RETURN(status=1)
          ! check ...
          if ( trim(var_levtype) /= trim(leo%levtype) ) then
            write (gol,'("requested output on level type `",a,"` while variable `",a,"` is defined on `",a,"`")') &
                              trim(leo%levtype), trim(var_name), trim(var_levtype); call goErr
            TRACEBACK; status=1; return
          end if
        end do ! output variables
        
        ! check ...
        if ( trim(level_names) /= 'all' ) then
          write (gol,'("output on input levels requires all, requested `",a,"`")') trim(level_names); call goErr
          TRACEBACK; status=1; return
        end if
        ! number of levels could only be defined on first output,
        ! set now to negative value to enable this:
        leo%nlev = -999
          
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! *

    ! ok
    status = 0

  end subroutine LE_Output_Dat2_Init


  ! ***


  subroutine LE_Output_Dat2_Done( leo, status )

#ifdef with_netcdf
    use NetCDF          , only : NF90_Close
#endif
    use LE_Output_Common, only : Done

    ! --- in/out --------------------------------

    type(T_LE_Output_Dat2), intent(inout)   ::  leo
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Dat2_Done'

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

    ! clear:
    deallocate( leo%LE_Dat, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! clear storage for tracer fields:
    deallocate( leo%idat     , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( leo%datcomp  , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( leo%name_dat , stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( leo%unit_dat , stat=status )
    IF_NOTOK_RETURN(status=1)
    !deallocate( leo%unitconv , stat=status )
    !IF_NOTOK_RETURN(status=1)
    deallocate( leo%varid_dat, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! layers or heights ?
    select case ( trim(leo%levtype) )
      !~ model levels:
      case ( 'levels', 'halflevels', 'input_levels', 'input_halflevels' )
        ! clear storage for level indices:
        deallocate( leo%ilev )
!      !~ height levels:
!      case ( 'heights', 'elevations' )
!        ! nothing to be done ..
      !~ unknown
      case default
        write (gol,'("unsupported level name : ",a)') trim(leo%levtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! done with common stuff ...
    call Done( leo%com, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Output_Dat2_Done


  ! ***


  subroutine LE_Output_Dat2_PutOut( leo, t, status )

    use Binas  , only : xm_air
    use GO     , only : goc
    use GO     , only : TDate, TIncrDate, NewDate, IncrDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>), rTotal, dTotal, iTotal
    use GO     , only : wrtgol, Precisely, MidNight
    use GO     , only : goMatchValue, goReplace
    use Num    , only : Interp_Lin
    use LE_Grid, only : glb_ugg, glb_ugg_ustag, glb_ugg_vstag
    use C3PO   , only : T_Grid_NcDef

#ifdef with_netcdf
    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER, NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT
#endif

    use Dims   , only : nx, ny, nz

    use LE_Output_Common, only : PutOut_GlobalAttributes

    use LE_Output_Tools , only : LE_Output_Define_Dims_Lev
    use LE_Output_Tools , only : LE_Output_Define_Dims_Time

    use LE_Output_Tools , only : LE_Output_Define_Vars_Lev
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time

    use LE_Output_Tools , only : LE_Output_Put_Var_Lev
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
    
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains

    use LE_Data         , only : LE_Data_GetPointer

    ! --- in/out --------------------------------

    type(T_LE_Output_Dat2), intent(inout)  ::  leo
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Output_Dat2_PutOut'

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
    real, allocatable     ::  pat(:,:)       ! (0:nx,0:ny)
!    real, allocatable     ::  field(:,:,:)   ! (0:nx,0:ny,nlev)
    integer               ::  nlev

    logical               ::  enabled
    character(len=32)     ::  varname
    character(len=256)    ::  cf_standard_name, cf_long_name, cf_units
    character(len=512)    ::  comment
    
    real                  ::  hh(nx,ny,0:nz+1)
    real                  ::  hsamp(nx,ny)
    integer               ::  ilast

    integer               ::  icountry
    integer               ::  iemis
    integer               ::  itr

    integer               ::  imf
    logical               ::  match

    real, pointer         ::  pdata(:,:,:)
    integer               ::  lbo(3), ubo(3)

    ! --- begin ---------------------------------
    
    !! info ...
    !write (gol,'(a,": start")') rname; call goPr

#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif

    ! for multiples of dhour only ...
    if ( .not. Precisely(t,leo%dhour,'hour')  ) then
      status=0; return
    end if

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

      ! number of levels not defined yet?
      if ( leo%nlev < 0 ) then
        ! check shape of first output field; select data id:
        l = 1
        idat = leo%idat(l)
        ! get pointer to data:
        call LE_Data_GetPointer( trim(leo%LE_Dat(idat)%name), pdata, status, &
                                   component=leo%datcomp(l) )
        IF_NOTOK_RETURN(status=1)
        ! reset shape:
        leo%nlev = size(pdata,3)
        ! storage:
        allocate( leo%ilev(leo%nlev), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! fill:
        do i = 1, leo%nlev
          leo%ilev(i) = i
        end do
        ! start at zero for half levels:
        if ( trim(leo%levtype) == 'input_halflevels' ) leo%ilev = leo%ilev - 1
      end if

      ! root only:
      if ( goc%root ) then

#ifdef with_netcdf
        ! set creation mode flag:
        if ( leo%replace ) then
          cmode = NF90_CLOBBER       ! overwrite existing files
        else
          cmode = NF90_NOCLOBBER     ! do not overwrite existing files
        end if

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
        call PutOut_GlobalAttributes( leo%com, leo%ncid, status )
        IF_NOTOK_RETURN(status=1)

        ! grid dimensions/variables
        select case ( trim(leo%gridtype) )
          !~ cell centers:
          case ( 'cells' )
            ! subset numbered as grid:
            call glb_ugg%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                        dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )!, &
                                        !subset=(/leo%i1,leo%i2,leo%j1,leo%j2/) )
            IF_NOTOK_RETURN(status=1)
          !~
          case ( 'u-edges' )
            ! grid numbered from 1, leo%i numbered from 0:
            call glb_ugg_ustag%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                        dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )!, &
                                        !subset=(/leo%i1+1,leo%i2+1,leo%j1,leo%j2/) )
            IF_NOTOK_RETURN(status=1)
          !~
          case ( 'v-edges' )
            ! grid numbered from 1, leo%j numbered from 0:
            call glb_ugg_vstag%DefGrid_NetCDF( gncd, leo%ncid, status, &
                                        dimid_lon=leo%dimid_lon, dimid_lat=leo%dimid_lat )!, &
                                        !subset=(/leo%i1,leo%i2,leo%j1+1,leo%j2+1/) )
            IF_NOTOK_RETURN(status=1)
          !~
          case default
            write (gol,'("could not define output grid for gridtype `",a,"`")') trim(leo%gridtype); call goErr
            TRACEBACK; status=1; return
        end select

        ! define level dimension:
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

          !! get names following CF conventions;
          !! store conversion factor for later usage:
          !call LE_CF_names( &
          !             leo%LE_Dat(idat)%name, leo%LE_Dat(idat)%unit, &
          !             cf_standard_name, cf_long_name, cf_units, &
          !             leo%unitconv(l), comment, &
          !             status )
          !IF_NOTOK_RETURN(status=1)
          ! asume variables follow cf already:
          cf_standard_name = ''
          cf_long_name = trim(leo%LE_Dat(idat)%long_name)
          if ( leo%datcomp(l) == 'u' ) cf_long_name = trim(cf_long_name)//' u-component'
          if ( leo%datcomp(l) == 'v' ) cf_long_name = trim(cf_long_name)//' v-component'
          cf_units     = trim(leo%LE_Dat(idat)%unit)
          comment      = ''

          ! variable name:
          varname = trim(leo%name_dat(l))
          ! vector component 'uv.u', replace dot:
          call goReplace( varname, '.', '_', status )
          IF_NOTOK_RETURN(status=1)

          ! define variable:
          select case ( leo%LE_Dat(idat)%rank )
            case ( 2 )
              if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf
                status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                         (/leo%dimid_lon,leo%dimid_lat/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
                status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'latitude longitude')
                IF_NF90_NOTOK_RETURN(status=1)
#endif
              else
#ifdef with_netcdf
                status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                         (/leo%dimid_lon,leo%dimid_lat,leo%dimid_time/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
                status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'time latitude longitude')
                IF_NF90_NOTOK_RETURN(status=1)
#endif
              end if
            case ( 3 )
              if ( leo%LE_Dat(idat)%const ) then
#ifdef with_netcdf
                status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                         (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev/), varid )
                IF_NF90_NOTOK_RETURN(status=1)
                status = nf90_put_att( leo%ncid, varid, '_CoordinateAxes', 'level latitude longitude')
                IF_NF90_NOTOK_RETURN(status=1)
#endif
              else
#ifdef with_netcdf
                status = NF90_Def_Var( leo%ncid, trim(varname), NF90_REAL, &
                                         (/leo%dimid_lon,leo%dimid_lat,leo%dimid_lev,leo%dimid_time/), varid )
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
          ! store units for later usage (GrADS ctl file):
          leo%unit_dat(l) = trim(cf_units)
          ! add comment:
          if ( len_trim(comment) > 0 ) then
            status = nf90_put_att( leo%ncid, varid, 'comment', trim(comment) )
            IF_NF90_NOTOK_RETURN(status=1)
          end if
#endif

          ! store variable id:
          leo%varid_dat(l) = varid

        end do  ! written tracers

        ! end defintion mode:
#ifdef with_netcdf
        status = NF90_EndDef( leo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
#endif

      end if ! root

      ! no records written yet:
      leo%itrec = 0

    end if  ! new file
    
    ! * write
    
    ! next time record:
    leo%itrec = leo%itrec + 1

    ! root ?
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

        ! grid dimensions/variables
        select case ( trim(leo%gridtype) )
          !~ cell centers:
          case ( 'cells' )
            ! write grid to netCDF file
            call glb_ugg%PutGrid_NetCDF( gncd, status )
            IF_NOTOK_RETURN(status=1)
          !~
          case ( 'u-edges' )
            ! write grid to netCDF file
            call glb_ugg_ustag%PutGrid_NetCDF( gncd, status )
            IF_NOTOK_RETURN(status=1)
          !~
          case ( 'v-edges' )
            ! write grid to netCDF file
            call glb_ugg_vstag%PutGrid_NetCDF( gncd, status )
            IF_NOTOK_RETURN(status=1)
          !~
          case default
            write (gol,'("could not define output grid for gridtype `",a,"`")') trim(leo%gridtype); call goErr
            TRACEBACK; status=1; return
        end select

        ! write level indices
        call LE_Output_Put_Var_Lev(leo%ncid, leo%varid_lev, leo%nlev, &
                                   trim(leo%levtype), leo%ilev, leo%heights, status)

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

    ! loop over all fields to be written:
    do l = 1, leo%ndat
    
      ! global tracer index:
      idat = leo%idat(l)

      !! info ..
      !write (gol,'(a,": write `",a,"` ...")') rname, trim(leo%LE_Dat(idat)%name); call goPr
      !write (gol,'(a,":   rank  ",i0)') rname, leo%LE_Dat(idat)%rank; call goPr
      !write (gol,'(a,":   const ",l1," itrec ",i0)') rname, leo%LE_Dat(idat)%const, leo%itrec; call goPr

      ! constant fields written only once:
      if ( leo%LE_Dat(idat)%const .and. (leo%itrec > 1) ) cycle
      
      ! get pointer to local data:
      call LE_Data_GetPointer( trim(leo%LE_Dat(idat)%name), pdata, status, &
                                 component=leo%datcomp(l) )
      IF_NOTOK_RETURN(status=1)
      
      ! 2d or 3d ?
      select case ( leo%LE_Dat(idat)%rank )

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 2 )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! write data:
          if ( leo%LE_Dat(idat)%const ) then
            ! write 2D field without level and time index:
            call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), -1, -1, &
                               pdata(leo%i1:leo%i2,leo%j1:leo%j2,1), status )
            IF_NOTOK_RETURN(status=1)
          else
            ! write 2D field without level but with time index:
            call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), -1, leo%itrec, &
                               pdata(leo%i1:leo%i2,leo%j1:leo%j2,1), status )
            IF_NOTOK_RETURN(status=1)
          end if

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 3 )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! * vertical selection or interpolation:
          
          ! index range:
          lbo = lbound(pdata)
          ubo = ubound(pdata)
          ! storage for 2D field:
          allocate( pat(lbo(1):ubo(1),lbo(2):ubo(2)), stat=status )
          IF_NOTOK_RETURN(status=1)

          !! testing ...
          !write (gol,'(a,":   levtype `",a,"` ...")') rname, trim(leo%levtype); call goPr
          !write (gol,'(a,":   nlev ",i0)') rname, leo%nlev; call goPr

          ! which output levels ?
          select case ( trim(leo%levtype) )

            !~ selected model levels:
            case ( 'levels', 'halflevels', 'input_levels', 'input_halflevels' )

              ! loop over selected layer:
              do k = 1, leo%nlev
                ! global level index:
                ilev = leo%ilev(k)
                ! extract 2D field:
                !pat = field(:,:,ilev)
                pat = pdata(:,:,ilev)
                ! write data:
                if ( leo%LE_Dat(idat)%const ) then
                  ! write 2D field with level but without time index:
                  call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), k, -1, &
                                     pat(leo%i1:leo%i2,leo%j1:leo%j2), status )
                  IF_NOTOK_RETURN(status=1)
                else
                  !! testing ..
                  !write (gol,'(a,": write (",i0,":",i0,",",i0,":",i0,") ",a," ",i0," record ",i0)') &
                  !               rname, leo%i1, leo%i2, leo%j1, leo%j2, &
                  !               trim(leo%levtype), k, leo%itrec; call goPr
                  ! write 2D field with level and time index:
                  call LE_Output_Put_Var_Domains( leo%ncid, leo%varid_dat(l), k, leo%itrec, &
                                     pat(leo%i1:leo%i2,leo%j1:leo%j2), status )
                 IF_NOTOK_RETURN(status=1)
                end if
              end do   ! selected model layers

!            !~ height levels:
!            case ( 'heights', 'elevations' )
!
!              ! loop over target levels:
!              do k = 1, leo%nlev
!                ! sample height:
!                if ( trim(leo%levtype) == 'heights' ) then
!                  ! relative to orography:
!                  hsamp = hh(:,:,0) + leo%heights(k)
!                else if ( trim(leo%levtype) == 'elevations' ) then
!                  ! absolute, interpolation will take surface value if below orography:
!                  hsamp = leo%heights(k)
!                else
!                  write (gol,'("no sample height defined for level type : ",a)') trim(leo%levtype); call goErr
!                  TRACEBACK; status=1; return
!                end if
!                ! loop over horizontal cells:
!                do j = 1, ny
!                  do i = 1, nx
!                    ! vertical interpolation:
!                    call Interp_Lin( hh(i,j,:), field(i,j,:), hsamp(i,j), pat(i,j), ilast, status )
!                    IF_NOTOK_RETURN(status=1)
!                  end do  ! i
!                end do  ! j
!                ! write data:
!                if ( leo%LE_Dat(idat)%const ) then
!                  status = NF90_Put_Var( leo%ncid, leo%varid_dat(l), pat(leo%i1:leo%i2,leo%j1:leo%j2), &
!                                 start=(/1,1,k/), count=(/leo%ni,leo%nj,1/) )
!                  IF_NF90_NOTOK_RETURN(status=1)
!                else
!                  status = NF90_Put_Var( leo%ncid, leo%varid_dat(l), pat(leo%i1:leo%i2,leo%j1:leo%j2), &
!                                 start=(/1,1,k,leo%itrec/), count=(/leo%ni,leo%nj,1,1/) )
!                 IF_NF90_NOTOK_RETURN(status=1)
!                end if
!              end do   ! heights

            !~ unknown ...
            case default
              write (gol,'("unsupported level type : ",a)') trim(leo%levtype); call goErr
              TRACEBACK; status=1; return

          end select
          
          ! clear:
          deallocate( pat, stat=status )
          IF_NOTOK_RETURN(status=1)

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

    !! clear:
    !deallocate( pat )
    !if ( allocated(field) ) deallocate( field )

    !! info ...
    !write (gol,'(a,": end")') rname; call goPr

    ! ok
    status = 0

  end subroutine LE_Output_Dat2_PutOut


  ! ***


  subroutine Write_GrADS_Ctl( leo, status )

    use dims   , only : nx, ny
    use GrADS_Ctl
    use LE_Grid, only : glb_ugg

    ! --- in/out ---------------------------------

    type(T_LE_Output_Dat2), intent(inout)    ::  leo
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
      write (leo%grads_ctl_file,'(a,".ctl")') trim(leo%grads_ctl_file)

      ! daily or less ?
      select case ( trim(leo%collect) )
        ! collect daily for [00,24)
        case ( 'daily', 'daily24' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
          write (leo%grads_ctl_base,'(a,".nc")') trim(leo%grads_ctl_base)
        ! files with instant fields:
        case ( 'instant' )
          ! template for data files:
          write (leo%grads_ctl_base,'("^",a,"_",a,"_",a,"_%y4%m2%d2_%h2%n2")') &
                    trim(leo%com%model), trim(leo%com%expid), trim(leo%name)
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
        case ( 'levels', 'halflevels', 'input_levels', 'input_halflevels' )
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



end module LE_Output_Dat2
