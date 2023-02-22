!###############################################################################
!
! NAME
!   LE_Bound_LE
!
! DESCRIPTION
!
!   Read boundary conditions from LOTOS-EUROS output in NetCDF form.
!
! HISTORY
!   Former 'boundary_le' module.
!
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
! deze versie leest elk uur een nieuwe cbound array  in
! ger, feb 2007
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!
! this module handles boundary data that is precalculated by the LE-model
! juli 2005
! adjusted april 2006
! adjusted Oct 2007, call scan_control
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

module LE_Bound_LE

  use GO, only : gol, goPr, goErr
  use Grid, only : TllGridInfo
#ifdef with_netcdf
  use NetCDF, only : NF90_NOERR, nf90_strerror
#endif

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  LE_Bound_LE_Init, LE_Bound_LE_Done
  public  ::  LE_Bound_LE_Get
  public  ::  LE_Bound_LE_Initial


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_LE'


  ! --- local ------------------------------------

  ! arachive description:
  character(len=512)    ::  arch_path
  character(len=32)     ::  arch_model
  character(len=64)     ::  arch_expid
  character(len=32)     ::  arch_name

  ! save for future calls:
  logical               ::  cache_setup
  ! dimensions of source fields:
  integer               ::  cache_nlon, cache_nlat, cache_nlev
  ! indices of covering grid cells in source fields:
  integer, pointer      ::  cache_ij_west(:,:,:), cache_ij_east(:,:,:)
  integer, pointer      ::  cache_ij_south(:,:,:), cache_ij_north(:,:,:)
  integer, pointer      ::  cache_ij(:,:,:)


contains


  ! ====================================================================


  subroutine LE_Bound_LE_Init( rcF, status )

    use GO, only : TrcFile, ReadRc
    use GO, only : goVarValue

    ! --- in/out ------------------------------

    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_LE_Init'

    ! --- local ------------------------------

    character(len=512)    ::  key

    ! --- begin -------------------------------

    !
    ! settings
    !

    ! archive description:
    call ReadRc( rcF, 'le.bound.le.path', arch_path, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'le.bound.le.key', key, status )
    IF_NOTOK_RETURN(status=1)

    ! extract parts from key:
    !
    !   model=LE;expid=base;name='conc'
    !
    arch_model = 'LE'
      call goVarValue( key, ';', 'model', '=', arch_model, status )
      IF_ERROR_RETURN(status=1)
    arch_expid = 'base'
      call goVarValue( key, ';', 'expid', '=', arch_expid, status )
      IF_ERROR_RETURN(status=1)
    arch_name = 'conc'
      call goVarValue( key, ';', 'name' , '=', arch_name , status )
      IF_ERROR_RETURN(status=1)

    !
    ! setup
    !

    ! empty cache:
    cache_setup = .false.
    nullify( cache_ij_west  )
    nullify( cache_ij_east  )
    nullify( cache_ij_south )
    nullify( cache_ij_north )
    nullify( cache_ij )

    !
    ! done
    !

    ! ok
    status = 0

  end subroutine LE_Bound_LE_Init


  ! ***


  subroutine LE_Bound_LE_Done( status )

    use Grid, only : Done

    ! --- in/out ------------------------------

    integer, intent(out)  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_LE_Done'

    ! --- begin -------------------------------

    ! clear:
    if ( cache_setup ) then
      ! clear:
      deallocate( cache_ij_west, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( cache_ij_east, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( cache_ij_south, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( cache_ij_north, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( cache_ij, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! reset flag:
      cache_setup = .false.
    end if

    ! ok
    status = 0

  end subroutine LE_Bound_LE_Done


  ! ***


  subroutine LE_Bound_LE_Get( t, bc_filled, status )

    use GO, only : TDate, TIncrDate, IncrDate, Get, ExpandTime
    use GO, only : Precisely, MidNight
    use GO, only : operator(-), operator(/=), operator(==), wrtgol
#ifdef with_netcdf
    use NetCDF
#endif
    use C3PO           , only : T_Grid_Ugg
    use dims           , only : nx, ny, nz, nspec
    use indices        , only : specname, specunit
    use LE_Bound_Common, only : bc_west, bc_east, bc_south, bc_north
    use LE_Grid        , only : ugg
    use LE_Grid        , only : ugg_west, ugg_east, ugg_south, ugg_north

    ! --- in/out ---------------------------------

    type(TDate), intent(in)     ::  t
    logical, intent(inout)      ::  bc_filled(nspec)
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_LE_Get'

    ! --- local ----------------------------------

    ! --- local ----------------------------------

    character(len=512)    ::  fname
    integer               ::  year, month, day
    integer               ::  nt    
    real(8), allocatable  ::  times(:)
    type(TDate)           ::  t_in
    integer               ::  irec
    integer               ::  it
    type(T_Grid_Ugg)      ::  ugg_bc
    integer               ::  ncid, dimid, varid
    real, allocatable     ::  lons(:), lats(:)
    real, allocatable     ::  lons_bnds(:,:), lats_bnds(:,:)
    real, allocatable     ::  cc(:,:,:)
    real                  ::  add_offset, scale_factor
    character(len=32)     ::  vname
    integer               ::  ispec
    integer               ::  ilon, ilat
    integer               ::  iz
    character(len=64)     ::  unit
    real                  ::  unit_fac
    integer               ::  slen

    ! --- begin -------------------------------------

    ! check implementation:
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read boundary conditions from LE NetCDF output.")'); call goErr
    write (gol,'("Define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! info ...
    call wrtgol( '        get_boundary_le_nc for ', t ); call goPr

    ! only at complete hour:
    if ( .not. Precisely( t, 1.0, 'hour' ) ) then
      write (gol,'("          only at whole hours ...")'); call goPr
      status=0; return
    end if

    ! extract time values for filename:
    call Get( t, year=year, month=month, day=day )

    ! filename :  <path>/<model>_<expid>_<name>_<yyyymmdd>.nc
    write (fname,'(a,"/",a,2("_",a),"_",i4.4,2i2.2,".nc")') &
            trim(arch_path), trim(arch_model), trim(arch_expid), trim(arch_name), &
            year, month, day
    
#ifdef with_netcdf
    ! open file (read-only):
    status = NF90_Open( trim(fname), NF90_NOWRITE, ncid )
    if (status/=NF90_NOERR) then
      write (gol,'("NF90 error: ",a)') trim(nf90_strerror(status)); call goErr
      write (gol,'("opening file : ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! setup grid etc ?
    if ( .not. cache_setup ) then

      ! get dimensions:
      status = NF90_Inq_DimID( ncid, 'longitude', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Inquire_Dimension( ncid, dimid, len=cache_nlon )
      IF_NF90_NOTOK_RETURN(status=1)
      ! get dimension:
      status = NF90_Inq_DimID( ncid, 'latitude', dimid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Inquire_Dimension( ncid, dimid, len=cache_nlat )
      IF_NF90_NOTOK_RETURN(status=1)

      ! storage:
      allocate( lons(cache_nlon), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( lats(cache_nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( lons_bnds(2,cache_nlon), stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( lats_bnds(2,cache_nlat), stat=status )
      IF_NOTOK_RETURN(status=1)

      ! get axis:
      status = NF90_Inq_VarID( ncid, 'longitude', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, lons )
      IF_NF90_NOTOK_RETURN(status=1)
      ! boundary axis
      status = NF90_Inq_VarID( ncid, 'longitude_bnds', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, lons_bnds )
      IF_NF90_NOTOK_RETURN(status=1)
      ! get axis:
      status = NF90_Inq_VarID( ncid, 'latitude', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, lats )
      IF_NF90_NOTOK_RETURN(status=1)
      ! boundary axis
      status = NF90_Inq_VarID( ncid, 'latitude_bnds', varid )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Get_Var( ncid, varid, lats_bnds )
      IF_NF90_NOTOK_RETURN(status=1)

      ! setup grid:
      call ugg_bc%Init( lons(:), lons_bnds(:,:), &
                        lats(:), lats_bnds(:,:), status )
      IF_NF90_NOTOK_RETURN(status=1)

      ! clear:
      deallocate( lons, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( lats, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( lons_bnds, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( lats_bnds, stat=status )
      IF_NOTOK_RETURN(status=1)

      ! simply assign concentrations from a source cell to a boundary condition cell;
      ! select and store source cell indices:
      allocate( cache_ij_west(2,ugg_west%nlon,ugg_west%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg_west%nlat
        do ilon = 1, ugg_west%nlon
          call ugg_bc%GetLocation( ugg_west%longitude_1d(ilon), ugg_west%latitude_1d(ilat), &
                              cache_ij_west(1,ilon,ilat), cache_ij_west(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do
      allocate( cache_ij_east(2,ugg_east%nlon,ugg_east%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg_east%nlat
        do ilon = 1, ugg_east%nlon
          call ugg_bc%GetLocation( ugg_east%longitude_1d(ilon), ugg_east%latitude_1d(ilat), &
                              cache_ij_east(1,ilon,ilat), cache_ij_east(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do
      allocate( cache_ij_south(2,ugg_south%nlon,ugg_south%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg_south%nlat
        do ilon = 1, ugg_south%nlon
          call ugg_bc%GetLocation( ugg_south%longitude_1d(ilon), ugg_south%latitude_1d(ilat), &
                              cache_ij_south(1,ilon,ilat), cache_ij_south(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do
      allocate( cache_ij_north(2,ugg_north%nlon,ugg_north%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg_north%nlat
        do ilon = 1, ugg_north%nlon
          call ugg_bc%GetLocation( ugg_north%longitude_1d(ilon), ugg_north%latitude_1d(ilat), &
                              cache_ij_north(1,ilon,ilat), cache_ij_north(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do
      ! core grid (used for initial condition)
      allocate( cache_ij(2,ugg%nlon,ugg%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      do ilat = 1, ugg%nlat
        do ilon = 1, ugg%nlon
          call ugg_bc%GetLocation( ugg%longitude(ilon,ilat), ugg%latitude(ilon,ilat), &
                              cache_ij(1,ilon,ilat), cache_ij(2,ilon,ilat), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end do

      ! get dimension:
      status = NF90_Inq_DimID( ncid, 'level', dimid )
      if (status /= NF90_NOERR ) then
        status = NF90_Inq_VarID( ncid, 'lev', varid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
      status = NF90_Inquire_Dimension( ncid, dimid, len=cache_nlev )
      IF_NF90_NOTOK_RETURN(status=1)

      ! check ...
      if ( cache_nlev /= nz ) then
        write (gol,'("number of levels in LE nc output does not match with current model:")'); call goErr
        write (gol,'("  nc file levels : ",i6)') cache_nlev; call goErr
        write (gol,'("  model levels   : ",i6)') nz; call goErr
        TRACEBACK; status=1; return
      end if

      ! reset flag:
      cache_setup = .true.

    end if   ! setup cache

    ! get time dimension:
    status = NF90_Inq_DimID( ncid, 'time', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=nt )
    IF_NF90_NOTOK_RETURN(status=1)
    ! storage for times:
    allocate( times(nt), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! time variable:
    status = NF90_Inq_VarID( ncid, 'time', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, times )
    IF_NF90_NOTOK_RETURN(status=1)
    ! time units:
    status = NF90_Get_Att( ncid, varid, 'units', unit )
    IF_NF90_NOTOK_RETURN(status=1)
    ! search:
    irec = -999
    do it = 1, nt
      ! expand:
      call ExpandTime( real(times(it)), unit, 'standard', t_in, status )
      IF_NOTOK_RETURN(status=1)
      ! compare:
      if ( t == t_in ) then
        irec = it
        exit
      end if
    end do
    if ( irec < 0 ) then
      write (gol,'("times in LE nc output do not match with requested time:")'); call goErr
      call wrtgol( '  model time      : ', t ); call goErr
      write (gol,'("  output times    ")'); call goErr
      do it = 1, nt
        call ExpandTime( real(times(it)), unit, 'standard', t_in, status )
        IF_NOTOK_RETURN(status=1)
        call wrtgol( '     ', t_in ); call goPr
      end do
      TRACEBACK; status=1; return
    end if
    ! clear:
    deallocate( times, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! storage for single field:
    allocate( cc(cache_nlon,cache_nlat,cache_nlev), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over all tracers:
    do ispec = 1, nspec

      ! variable name:
      vname = trim(specname(ispec))
 

      ! return variable id, or error status:
      status = NF90_Inq_VarID( ncid, trim(vname), varid )
      ! check status on error 'not a variable':

      ! outputname of PM in capitals
      if ((vname == 'ppm25') .and. ( status == NF90_ENOTVAR )) vname='PPM25'  
      if ((vname == 'ppm10') .and. ( status == NF90_ENOTVAR )) vname='PPM10'  

      ! get variable id for name:
      status = NF90_Inq_VarID( ncid, trim(vname), varid ) 
      ! not a variable?
      if ( status == NF90_ENOTVAR ) then
        ! info ...
        write (gol,'("        LE nc field not found : ",a)') trim(vname); call goPr
        ! not available; try next:
        cycle
      end if
      ! other errors ?
      IF_NF90_NOTOK_RETURN(status=1)

      ! info ..
      write (gol,'("        read LE nc field : ",a)') trim(vname); call goPr

      ! read record:
      status = NF90_Get_Var( ncid, varid, cc, start=(/1,1,1,irec/), &
                         count=(/cache_nlon,cache_nlat,cache_nlev,1/) )
      IF_NF90_NOTOK_RETURN(status=1)

      ! read packing numbers:
      status = NF90_Get_Att( ncid, varid, 'add_offset', add_offset )
      if ( status /= 0 ) then
        ! file not packed
        add_offset = 0.0
      end if
      status = NF90_Get_Att( ncid, varid, 'scale_factor', scale_factor )
      if ( status /= 0 ) then
        ! file not packed
        scale_factor = 1.0
      end if
    
      ! unpack
      cc = cc * scale_factor + add_offset      
      
      ! read unit:
      status = NF90_Get_Att( ncid, varid, 'units', unit )
      IF_NF90_NOTOK_RETURN(status=1)
      ! read length:
      status = NF90_Inquire_Attribute( ncid, varid, 'units', len=slen )
      IF_NF90_NOTOK_RETURN(status=1)
      ! truncate:
      unit = unit(1:slen)

      ! conversion factor:
      if ( (trim(unit) == 'unknown unit') .or. (trim(unit) == trim(specunit(ispec))) ) then
        unit_fac = 1.0
      else
        select case ( trim(unit) )
          case ( 'mole mole-1', 'mole mole**-1'   )
            select case ( trim(specunit(ispec)) )
              case ( 'ppm' ) ; unit_fac = 1.0e6
              case ( 'ppb' ) ; unit_fac = 1.0e9
              case default
                write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return
            end select
          case ( 'kg m-3', 'kg m**-3' )
            select case ( trim(specunit(ispec)) )
              case ( 'ug/m3', 'ug m-3' ) ; unit_fac = 1.0e9
              case default
                write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return
            end select
          case default
            write (gol,'("do not know how to convert `",a,"`")') trim(unit); call goErr
            TRACEBACK; status=1; return
        end select
      end if

      ! convert:
      cc = cc * unit_fac

      ! assume levels are the same; loop over levels:
      do iz = 1, nz
        ! copy concentrations from coarse grid to boundary conditions:
        do ilat = 1, ugg_west%nlat
          do ilon = 1, ugg_west%nlon
            bc_west(ilat,iz,ispec) = cc(cache_ij_west(1,ilon,ilat),cache_ij_west(2,ilon,ilat),iz)
          end do
        end do
        do ilat = 1, ugg_east%nlat
          do ilon = 1, ugg_east%nlon
            bc_east(ilat,iz,ispec) = cc(cache_ij_east(1,ilon,ilat),cache_ij_east(2,ilon,ilat),iz)
          end do
        end do
        do ilat = 1, ugg_south%nlat
          do ilon = 1, ugg_south%nlon
            bc_south(ilon,iz,ispec) = cc(cache_ij_south(1,ilon,ilat),cache_ij_south(2,ilon,ilat),iz)
          end do
        end do
        do ilat = 1, ugg_north%nlat
          do ilon = 1, ugg_north%nlon
            bc_north(ilon,iz,ispec) = cc(cache_ij_north(1,ilon,ilat),cache_ij_north(2,ilon,ilat),iz)
          end do
        end do
      end do  ! levels

      ! set flag:
      bc_filled(ispec) = .true.

    end do  ! specs

    ! clear:
    deallocate( cc, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine LE_Bound_LE_Get


  ! ***


  subroutine LE_Bound_LE_Initial( t, c, ic_filled, status )

    use GO, only : TDate, TIncrDate, Get, ExpandTime
    use GO, only : Precisely, MidNight, operator(-), operator(/=), operator(==)
    use GO, only : wrtgol
    use Grid, only : TllGridInfo, Init, GetLocation
#ifdef with_netcdf
    use NetCDF
#endif
    use dims           , only : nx, ny, nz, nspec
    use indices        , only : specname, specunit
    use LE_Bound_Common, only : bc_west, bc_east, bc_south, bc_north
    use LE_Grid        , only : ugg

    ! --- in/out ---------------------------------

    type(TDate), intent(in)     ::  t
    real, intent(inout)         ::  c(nx,ny,nz,nspec)
    logical, intent(inout)      ::  ic_filled(nspec)
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_LE_Initial'

    ! --- local ----------------------------------

    ! --- local ----------------------------------

    character(len=512)    ::  fname
    integer               ::  year, month, day
    integer               ::  nt
    real, allocatable     ::  times(:)
    type(TDate)           ::  t_in
    integer               ::  irec
    integer               ::  it
    integer               ::  ncid, dimid, varid
    real, allocatable     ::  cc(:,:,:)
    real                  ::  add_offset, scale_factor
    character(len=32)     ::  vname
    integer               ::  ispec
    integer               ::  ilon, ilat
    integer               ::  iz
    character(len=64)     ::  unit
    real                  ::  unit_fac
    integer               ::  slen

    ! --- begin -------------------------------------

    ! check implementation:
#ifndef with_netcdf
    write (gol,'("NetCDF support is required to read boundary conditions from LE NetCDF output.")'); call goErr
    write (gol,'("Define fpp macro `with_netcdf` to compile with NetCDF support ...")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! info ...
    call wrtgol( '        get_boundary_le_nc for ', t ); call goPr

    ! only at complete hour:
    if ( .not. Precisely( t, 1.0, 'hour' ) ) then
      write (gol,'("          only at whole hours ...")'); call goPr
      status=0; return
    end if

    ! extract time values for filename:
    call Get( t, year=year, month=month, day=day )

    ! filename :  <path>/<model>_<expid>_<name>_<yyyymmdd>.nc
    write (fname,'(a,"/",a,2("_",a),"_",i4.4,2i2.2,".nc")') &
            trim(arch_path), trim(arch_model), trim(arch_expid), trim(arch_name), &
            year, month, day

#ifdef with_netcdf
    ! open file (read-only):
    status = NF90_Open( trim(fname), NF90_NOWRITE, ncid )
    if (status/=NF90_NOERR) then
      write (gol,'("NF90 error: ",a)') trim(nf90_strerror(status)); call goErr
      write (gol,'("opening file : ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! get time dimension:
    status = NF90_Inq_DimID( ncid, 'time', dimid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Inquire_Dimension( ncid, dimid, len=nt )
    IF_NF90_NOTOK_RETURN(status=1)
    ! storage for dates:
    ! storage for times:
    allocate( times(nt), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! time variable:
    status = NF90_Inq_VarID( ncid, 'time', varid )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Get_Var( ncid, varid, times )
    IF_NF90_NOTOK_RETURN(status=1)
    ! time units:
    status = NF90_Get_Att( ncid, varid, 'units', unit )
    IF_NF90_NOTOK_RETURN(status=1)
    ! search:
    irec = -999
    do it = 1, nt
      ! expand:
      call ExpandTime( times(it), unit, 'standard', t_in, status )
      IF_NOTOK_RETURN(status=1)
      ! compare:
      if ( t == t_in ) then
        irec = it
        exit
      end if
    end do
    if ( irec < 0 ) then
      write (gol,'("times in LE nc output do not match with requested time:")'); call goErr
      call wrtgol( '  model time      : ', t ); call goErr
      write (gol,'("  output times    ")'); call goErr
      do it = 1, nt
        call ExpandTime( times(it), unit, 'standard', t_in, status )
        IF_NOTOK_RETURN(status=1)
        call wrtgol( '     ', t_in ); call goPr
      end do
      TRACEBACK; status=1; return
    end if
    ! clear:
    deallocate( times, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! storage for single field:
    allocate( cc(cache_nlon,cache_nlat,cache_nlev), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over all tracers:
    do ispec = 1, nspec

      ! variable name:
      vname = trim(specname(ispec))

      ! return variable id, or error status:
      status = NF90_Inq_VarID( ncid, trim(vname), varid )
      
      ! outputname of PM in capitals
      if ((vname == 'ppm25') .and. ( status == NF90_ENOTVAR )) vname='PPM25'  
      if ((vname == 'ppm10') .and. ( status == NF90_ENOTVAR )) vname='PPM10'  

      status = NF90_Inq_VarID( ncid, trim(vname), varid )
      ! check status on error 'not a variable':
      if ( status == NF90_ENOTVAR ) then
        ! info ...
        write (gol,'("        LE nc field not found : ",a)') trim(vname); call goPr
        ! not available; try next:
        cycle
      end if
      ! other errors ?
      IF_NF90_NOTOK_RETURN(status=1)

      ! info ..
      write (gol,'("        read LE nc field : ",a)') trim(vname); call goPr

      ! read record:
      status = NF90_Get_Var( ncid, varid, cc, start=(/1,1,1,irec/), &
                         count=(/cache_nlon,cache_nlat,cache_nlev,1/) )
      IF_NF90_NOTOK_RETURN(status=1)

      ! read packing numbers:
      status = NF90_Get_Att( ncid, varid, 'add_offset', add_offset )
      if ( status /= 0 ) then
        ! file not packed
        add_offset = 0.0
      end if
      status = NF90_Get_Att( ncid, varid, 'scale_factor', scale_factor )
      if ( status /= 0 ) then
        ! file not packed
        scale_factor = 1.0
      end if
    
      ! unpack
      cc = cc * scale_factor + add_offset      

      ! read unit:
      status = NF90_Get_Att( ncid, varid, 'units', unit )
      IF_NF90_NOTOK_RETURN(status=1)
      ! read length:
      status = NF90_Inquire_Attribute( ncid, varid, 'units', len=slen )
      IF_NF90_NOTOK_RETURN(status=1)
      ! truncate:
      unit = unit(1:slen)

      ! conversion factor:
      if ( (trim(unit) == 'unknown unit') .or. (trim(unit) == trim(specunit(ispec))) ) then
        unit_fac = 1.0
      else
        select case ( trim(unit) )
          case ( 'mole mole-1', 'mole mole**-1'   )
            select case ( trim(specunit(ispec)) )
              case ( 'ppm' ) ; unit_fac = 1.0e6
              case ( 'ppb' ) ; unit_fac = 1.0e9
              case default
                write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return
            end select
          case ( 'kg m-3', 'kg m**-3' )
            select case ( trim(specunit(ispec)) )
              case ( 'ug/m3', 'ug m-3' ) ; unit_fac = 1.0e9
              case default
                write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return
            end select
          case default
            write (gol,'("do not know how to convert `",a,"`")') trim(unit); call goErr
            TRACEBACK; status=1; return
        end select
      end if

      ! convert:
      cc = cc * unit_fac

      ! loop over levels:
      do iz = 1, nz
        ! loop over grid cells:
        do ilat = 1, ugg%nlat
          do ilon = 1, ugg%nlon
            ! copy concentrations:
            c(ilon,ilat,iz,ispec) = cc(cache_ij(1,ilon,ilat),cache_ij(2,ilon,ilat),iz)
          end do
        end do
        !
      end do  ! levels
      ! set flag:
      ic_filled(ispec) = .true.

    end do  ! specs

    ! clear:
    deallocate( cc, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine LE_Bound_LE_Initial


end module LE_Bound_LE
