!###############################################################################
!
! LE_Bound_Tools - LOTOS-EUROS boundary condition tools
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Bound_Tools

  use GO, only : gol, goErr, goPr
  
!  use CF, only : T_CF_AxisMapper

  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  LE_Bound_Fill
  public  ::  LE_Bound_Fill_Initial

  public  ::  PTQ_to_H

!  public  ::  T_Heights_to_LE
!  public  ::  Heights_to_LE_Init, Heights_to_LE_Done
!  public  ::  Heights_to_LE_Apply_Bounds
!  public  ::  Heights_to_LE_Apply_Initial
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_Tools'


  ! --- type --------------------------------------
  
!  ! mapping from 3D field defined at heights to model:
!  type T_Heights_to_LE
!    ! input dimensions:
!    integer                                 ::  nx_in, ny_in, nz_in
!    ! mappers for horizontal grid:
!    type(T_CF_AxisMapper)                   ::  Lon_Mapper
!    type(T_CF_AxisMapper)                   ::  Lat_Mapper
!    ! mappers for boundary condition grids:
!    type(T_CF_AxisMapper)                   ::  west_Lon_Mapper
!    type(T_CF_AxisMapper)                   ::  west_Lat_Mapper
!    type(T_CF_AxisMapper)                   ::  east_Lon_Mapper
!    type(T_CF_AxisMapper)                   ::  east_Lat_Mapper
!    type(T_CF_AxisMapper)                   ::  south_Lon_Mapper
!    type(T_CF_AxisMapper)                   ::  south_Lat_Mapper
!    type(T_CF_AxisMapper)                   ::  north_Lon_Mapper
!    type(T_CF_AxisMapper)                   ::  north_Lat_Mapper
!    ! mapping for model columns:
!    type(T_CF_AxisMapper), allocatable      ::  Lev_Mapper(:,:)  ! (nx,ny)
!    ! mapping for aloft columns:
!    type(T_CF_AxisMapper), allocatable      ::  aloft_Lev_Mapper(:,:)  ! (nx,ny)
!  end type
    
  
  
contains


  ! ====================================================================
  ! ===
  ! === boundary fields defined between pressure levels
  ! ===
  ! ====================================================================
  

  subroutine LE_Bound_Fill( ispec, mm, mm_air, &
                            bc_lons, bc_lats, bc_phlev, bc_ahlev, bc_vmr, bc_unit, &
                            status )
  
    use Binas, only : grav, twopi
    use JAQL , only : PotentialPressures
    use Num  , only : IntervalSum, Interval, Interval_Modulo, IntervalQuad_Const
    use Grid , only : TllGridInfo, Init, Done, FillGrid_AreaAverage, FillGrid_AA_Fast

    use dims   , only : runF
    use dims   , only : nx, ny, nz
    use indices, only : specname, specunit
    use dims   , only : bc_west, bc_east, bc_south, bc_north, caloft
    use LE_Data      , only : LE_Data_GetPointer    

    ! --- in/out ------------------------------
    
    integer, intent(in)             ::  ispec
    real, intent(in)                ::  mm, mm_air ! molemass
    real, intent(in)                ::  bc_lons(:)
    real, intent(in)                ::  bc_lats(:)
    real, intent(in)                ::  bc_phlev(:,:,:)  ! levels 1:bc_nlay+1
    real, intent(in)                ::  bc_ahlev(:,:,:)  ! levels 1:bc_nlay+1
    real, intent(in)                ::  bc_vmr(:,:,:)
    character(len=*), intent(in)    ::  bc_unit
    integer, intent(out)            ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Bound_Fill'
    
    ! number of horizontal hal cells at each boundary:
    integer, parameter    ::  nh = 1
    
    ! conversion factors:
    real, parameter   ::  vmr2ppb = 1.0e9   ! mole/mole -> mole/(1e9 mol)
    real, parameter   ::  kg2ug   = 1.0e9   ! kg -> ug
    
    ! --- local -------------------------------
    
    integer               ::  bc_nlon, bc_nlat, bc_nlay
    type(TllGridInfo)     ::  le_lli
    type(TllGridInfo)     ::  bc_lli
    real, allocatable     ::  bc_dm(:,:)
    real, allocatable     ::  bc_vcd__legrid(:,:,:)  ! vertical column density TM5 (kg tracer)/m2
    real, allocatable     ::  bc_phlev__legrid(:,:,:)  ! bc half level pressure at LE grid coordinates
    real, allocatable     ::  bc_ahlev__legrid(:,:,:)  ! bc half level altitude at LE grid coordinates
    integer               ::  l
    real                  ::  le_psurf, le_temp(nz+1), le_rh(nz+1)  ! pressure at surface, temperature and relative humidity
    real                  ::  le_hh(0:nz+1), le_ph(0:nz+1) 
    real                  ::  le_dm
    real                  ::  le_dh
    real                  ::  le_dp
    real                  ::  le_vcd                  ! vertical column density LOTOS-EUROS (kg tracer)/m2
    real, allocatable     ::  bc_ph(:), bc_prof(:)  ! TM5 pressure and concentration profile 
    integer               ::  ix, iy, iz
    integer               ::  lecore_ix, lecore_iy
    integer               ::  le_ix, le_iy
    integer               ::  ilast
    logical               ::  lfirst
    character(len=32)     ::  bc_vcd_units

    ! mappings and fracx added in orde to save computing time (VORtech)
    ! these mappings are w.r.t. a 'source' grid and a 'target' grid
    integer, allocatable :: target2source_i_l(:)
    integer, allocatable :: target2source_i_r(:)
    integer, allocatable :: target2source_j_l(:)
    integer, allocatable :: target2source_j_r(:)
    real,    allocatable :: fracxt(:,:,:,:)
    integer              :: ix1,ix2,jx1,jx2
    
    ! meteo data:
    real, pointer        ::  oro(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  h_m(:,:,:)     ! (lon,lat,1:nz)    
    real, pointer        ::  temp(:,:,:)     ! (lon,lat,1:nz)    

    ! --- begin -------------------------------
    
    ! point to meteo data:
    call LE_Data_GetPointer( 'oro', oro, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 't', temp, status, check_units='K' )
    IF_NOTOK_RETURN(status=1)   

    ! check unit ...
    select case ( trim(bc_unit) )
      case ( '1 = mole/mo', '1 = mole/mole', 'mole mole-1', 'mole mole**-1', 'kg kg**-1' )
        ! converted later on ..
      case ( 'ug/m3' )
        ! conversions are different
      case default
        write (gol,'("unsupported bc unit : ",a)') trim(bc_unit); call goErr
        TRACEBACK; status=1; return
    end select

    ! grid sizes:
    bc_nlon = size(bc_vmr,1)
    bc_nlat = size(bc_vmr,2)
    bc_nlay = size(bc_vmr,3)
    
    ! check ...
    if ( (size(bc_lons) /= bc_nlon) .or. &
         (size(bc_lats) /= bc_nlat) .or. &
         any(shape(bc_phlev) /= (/bc_nlon,bc_nlat,bc_nlay+1/)) ) then
      write (gol,'("one or more arguments have unexpected shape:")'); call goErr
      write (gol,'("  bc_lons   : ",i6)') shape(bc_lons); call goErr
      write (gol,'("  bc_lats   : ",i6)') shape(bc_lats); call goErr
      write (gol,'("  bc_phlev  : ",3i6)') shape(bc_phlev); call goErr
      write (gol,'("  bc_vmr    : ",3i6)') shape(bc_vmr); call goErr
      TRACEBACK; status=1; return
    end if

    ! define LE grid including nh halo cells at boundaries;
    !   core grid         :    1,2,...,nx-1,nx
    !   extended (nh=1)   :  1,2,3,...,nx  ,nx+1,nx+2
    call Init( le_lli, runF%westb +(0.5-nh)*runF%dlon, runF%dlon, nx+2*nh, &
                       runF%southb+(0.5-nh)*runF%dlat, runF%dlat, ny+2*nh,&
                       & status,'legrid')
    IF_NOTOK_RETURN(status=1)

    ! define bc grid:
    !source
    call Init( bc_lli, bc_lons(1), bc_lons(2)-bc_lons(1), bc_nlon, &
                       bc_lats(1), bc_lats(2)-bc_lats(1), bc_nlat, status&
                       &,'bcgrid')
    IF_NOTOK_RETURN(status=1)

    !   Vortech: allocate storage arrays
    ! grid is source
    ! fracxt is 4-dimensional. it should be allocated  when it is needed first;
    ! maybe the dimensions need not to be so large.
    ! start from west bound
    iX1 = 1
    call Interval_modulo( bc_lli%blon, le_lli%blon(0), twopi, iX1, status )
    IF_NOTOK_RETURN(status=1)
    ! start from cell found before
    iX2 = iX1
    call Interval_modulo( bc_lli%blon, le_lli%blon(le_lli%nlon), twopi, iX2, status )
    IF_NOTOK_RETURN(status=1)
    ! search lat index of model south bound:
    call Interval( bc_lli%blat, le_lli%blat(0), jX1, status )
    IF_NOTOK_RETURN(status=1)
    ! search lat index of model north bound:
    call Interval( bc_lli%blat, le_lli%blat(le_lli%nlat), jX2, status )
    IF_NOTOK_RETURN(status=1)
    ! storage:
    allocate(fracxt(le_lli%nlon,le_lli%nlat,ix1:ix2,jx1:jx2), stat=status)
    IF_NOTOK_RETURN(status=1)    

    ! grid is target: define mappings
    allocate(target2source_i_l(le_lli%nlon),stat=status)
    IF_NOTOK_RETURN(status=1)
    allocate(target2source_i_r(le_lli%nlon),stat=status)
    IF_NOTOK_RETURN(status=1)
    allocate(target2source_j_l(le_lli%nlat),stat=status)
    IF_NOTOK_RETURN(status=1)
    allocate(target2source_j_r(le_lli%nlat),stat=status)
    IF_NOTOK_RETURN(status=1)

    ! increasing pressure ax (top->down!) ?
    if ( bc_phlev(1,1,1) < bc_phlev(1,1,bc_nlay+1) ) then
      write (gol,'("boundary condition field defined on increasing pressure ax not supported yet ...")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( bc_dm(1:bc_nlon,bc_nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_vcd__legrid(1:le_lli%nlon,1:le_lli%nlat,1:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_phlev__legrid(1:le_lli%nlon,1:le_lli%nlat,0:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_ahlev__legrid(1:le_lli%nlon,1:le_lli%nlat,0:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_ph(0:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_prof(1:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)

    !
    !  >>> converting the whole grid is very slow ...
    !      should be replaced by conversions to lli_west etc
    !

    ! convert boundary model half level pressures and altitudes to LE horizontal grid;
    ! note difference in indices: 1:n+1 on input, 0:n on output:
    lfirst = .true.
    do l = 0, bc_nlay
      ! fill values "bc_phlev__legrid" defined on grid "le_lli" using
      !      values "bc_phlev" defined on grid "bc_lli"
      ! area weighted average:
      call FillGrid_AA_fast( le_lli, bc_phlev__legrid(:,:,l), &  ! Pa
                                bc_lli, bc_phlev(:,:,l+1), status,&
                                target2source_i_l,target2source_i_r,target2source_j_l,target2source_j_r,&
                                fracxt, le_lli%nlon,le_lli%nlat,ix1,ix2,jx1,jx2,lfirst )
      IF_NOTOK_RETURN(status=1)
      ! idem for altitudes:
      call FillGrid_AA_fast( le_lli, bc_ahlev__legrid(:,:,l), &  ! Pa
                                bc_lli, bc_ahlev(:,:,l+1), status,&
                                target2source_i_l,target2source_i_r,target2source_j_l,target2source_j_r,&
                                fracxt, le_lli%nlon,le_lli%nlat,ix1,ix2,jx1,jx2,lfirst )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! convert boundary model field to LE horizontal grid;
    do l = 1, bc_nlay
      ! switch over units:
      select case ( trim(bc_unit) )
        !
        ! mole mixing ratios (=volume mixing ratios):
        case ( '1 = mole/mo', '1 = mole/mole', 'mole mole-1', 'mole mole**-1' )
          ! mass per area; dimension bc_phlev(:,:,1:bc_nlay+1)
          bc_dm = abs( bc_phlev(:,:,l+1) - bc_phlev(:,:,l) )/grav  ! (kg air)/m2
          ! concentrations in volume mixing ratio;
          ! convert to masses; use that volume mixing ratio is mol mixing ratio:
          !  (mol trc) (kg trc)/(mol trc)
          !  --------- ------------------ (kg air)/m2 = (kg trc)/m2
          !  (mol air) (kg air)/(mol air)
          ! fill values "bc_vcd__legrid" defined on grid "le_lli" using
          !      values "bc_vmr*mm/mm_air*bc_dm" defined on grid "bc_lli"
          ! area weighted average:
          call FillGrid_AA_fast( le_lli, bc_vcd__legrid(:,:,l), &  ! (kg tracer)/m2
                                 bc_lli, bc_vmr(:,:,l)*mm/mm_air*bc_dm, status, &
                                 target2source_i_l, target2source_i_r, &
                                 target2source_j_l, target2source_j_r, &
                                 fracxt, &
                                 le_lli%nlon, le_lli%nlat, &
                                 ix1, ix2, jx1, jx2, lfirst )
          ! set units:
          bc_vcd_units = 'kg/m2'
        !
        ! mass mixing ratio's:
        case ( 'kg kg**-1' )
          ! mass per area; dimension bc_phlev(:,:,1:bc_nlay+1)
          bc_dm = abs( bc_phlev(:,:,l+1) - bc_phlev(:,:,l) )/grav  ! (kg air)/m2
          ! concentrations in mass mixing ratio; convert to masses:
          !  (kg trc)
          !  -------- (kg air)/m2 = (kg trc)/m2
          !  (kg air)
          ! fill values "bc_vcd__legrid" defined on grid "le_lli" using
          !      values "bc_mmr*bc_dm" defined on grid "bc_lli"
          ! (variable is named 'bc_vmr', but with this units it is actually 'bc_mmr')
          ! area weighted average:
          call FillGrid_AA_fast( le_lli, bc_vcd__legrid(:,:,l), &  ! (kg tracer)/m2
                                 bc_lli, bc_vmr(:,:,l)*bc_dm, status, &
                                 target2source_i_l, target2source_i_r, &
                                 target2source_j_l, target2source_j_r, &
                                 fracxt, &
                                 le_lli%nlon, le_lli%nlat, &
                                 ix1, ix2, jx1, jx2, lfirst )
          ! set units:
          bc_vcd_units = 'kg/m2'
        !
        ! mass concentrations:
        case ('ug/m3') 
          ! no extra conversions needed, units remain:
          call FillGrid_AA_fast( le_lli, bc_vcd__legrid(:,:,l), &     ! (ug tracer)/m3
                                 bc_lli, bc_vmr(:,:,l), status, &
                                 target2source_i_l, target2source_i_r, &
                                 target2source_j_l, target2source_j_r, &
                                 fracxt, &
                                 le_lli%nlon, le_lli%nlat, &
                                 ix1, ix2, jx1, jx2, lfirst )
          ! set units:
          bc_vcd_units = 'ug/m3'
        !
        case default
          write (gol,'("unsupported bc unit : ",a)') trim(bc_unit); call goErr
          TRACEBACK; status=1; return
      end select
    end do
        
    !
    !   <<<
    !
    
    ! western boundary conversion of pressure layering to km LE layering
    do iy = 1, ny
      ! most nearby LE cell index in 'normal' grid:
      lecore_ix = 1
      lecore_iy = iy
      ! LE cell index in extended le grid:
      le_ix =      1
      le_iy = nh + iy
      ! bc concentration profile:
      bc_prof = bc_vcd__legrid(le_ix,le_iy,:)
      ! bc pressure half levels;
      bc_ph = bc_phlev__legrid(le_ix,le_iy,:)
      ! LE surface pressure: copy from bc surface:
      le_psurf = bc_ph(0)
      ! LE temperature profile; nearest cell in core grid:
      le_temp(1:nz) = temp(lecore_ix,lecore_iy,1:nz)    ! K
      le_temp(nz+1) = le_temp(nz)      ! K  ; same as top level
      ! LE realitve humidity profile; use zero, no idea what it is in LE ...
      le_rh         = 0.0             ! (kg h2o)/(kg air)
      ! LE altitude profile; nearest cell in core grid:
      le_hh(0)    = oro(lecore_ix,lecore_iy,1)  ! m
      le_hh(1:nz) = le_hh(0) + h_m(lecore_ix,lecore_iy,:) ! m  
      le_hh(nz+1) = le_hh(nz) + 2.0e3    ! m ;  assume 2 km for aloft layer
      ! compute LE pressure half levels:
      call PotentialPressures( nz+1, le_psurf, le_temp, le_rh, le_hh, le_ph )
      if ( any(le_ph < 0.0) ) then
        write (gol,'("found negative pressures:")'); call goErr
        write (gol,*) ' le_psurf : ', le_psurf; call goErr
        write (gol,*) ' le_temp  : ', le_temp; call goErr
        write (gol,*) ' le_rh    : ', le_rh; call goErr
        write (gol,*) ' le_hh    : ', le_hh; call goErr
        write (gol,*) ' le_ph    : ', le_ph; call goErr
        TRACEBACK; status=1; return
      end if
      ! switch for original units:
      select case ( trim(bc_vcd_units) )
        !~ mass column density:
        case ( 'kg/m2' )
          ! map from TM5 profile to LE profile:
          ilast = 1
          do iz = 1, nz
            ! fill each LE layer with sum of one or more fractions of TM5 layers;
            ! use pressure axis to have mass-conservation; 
            ! use negative pressures to have increasing ax:
            call IntervalSum( -1.0*bc_ph, bc_prof, &
                              -1.0*le_ph(iz-1), -1.0*le_ph(iz), &
                             le_vcd, ilast, status )  ! (kg tracer)/m2
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target volume mixing ratio's:
              case ( 'ppb' )
                ! LE air mass per area:
                le_dm = ( le_ph(iz-1) - le_ph(iz) )/grav
                ! convert from (kg tracer)/m2 to volume mixing ratio:
                !              (kg air)/(mol air)      1
                !  (kg trc)/m2 ------------------ ----------- = (mol trc)/(mol air)
                !              (kg trc)/(mol trc) (kg air)/m2
                bc_west(iy,iz,ispec) = le_vcd * mm_air/mm / le_dm * vmr2ppb
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE layer thickness:
                le_dh = le_hh(iz) - le_hh(iz-1)
                ! convert from (kg tracer)/m2 to mass concentration (ug/m3):
                !  (kg trc)/m2 / m = (kg trc)/m3
                bc_west(iy,iz,ispec) = le_vcd / le_dh * kg2ug
              !~ not yet ...
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return                       
            end select
          end do  ! iz
        !~ input is mass concentrations, bc_prof is in (ug tr)/m3 :
        case ( 'ug/m3' )
          ! loop over levels:
          ilast = 1
          do iz = 1, nz
            ! fill model layer with integral over pressure range
            ! assuming that original values are constant within their layers;
            ! use negative pressures to have increasing ax:
            call IntervalQuad_Const( -1.0*bc_ph, bc_prof, &          ! Pa, (ug tracer)/m3
                                     -1.0*le_ph(iz-1), -1.0*le_ph(iz), &  ! Pa
                                    le_vcd, ilast, status)  ! (ug tracer)/m3*Pa
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE pressure gap:
                le_dp = le_ph(iz-1) - le_ph(iz)   ! Pa
                ! new mass concentration:
                bc_west(iy,iz,ispec) = le_vcd / le_dp  ! (ug tracer)
              !~ not yet ..
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return    
            end select 
          end do  ! levels
        !~ not yet ...
        case default
          write (gol,'("unsupported bc vcd units : ",a)') trim(bc_vcd_units); call goErr
          TRACEBACK; status=1; return
      end select
    end do ! iy

    ! eastern boundary conversion of pressure layering to LE layering
    do iy = 1, ny
      ! most nearby LE cell index in 'normal' grid:
      lecore_ix = nx
      lecore_iy = iy
      ! LE cell index in extended le grid:
      le_ix =      nx+2
      le_iy = nh + iy
      ! bc concentration profile:
      bc_prof = bc_vcd__legrid(le_ix,le_iy,:)
      ! bc pressure half levels:
      bc_ph = bc_phlev__legrid(le_ix,le_iy,:)
      ! LE surface pressure: copy from bc surface:
      le_psurf = bc_ph(0)
      ! LE temperature profile; nearest cell in core grid:
      le_temp(1:nz) = temp(lecore_ix,lecore_iy,1:nz)    ! K
      le_temp(nz+1) = le_temp(nz)      ! K  ; same as top level
      ! LE realitve humidity profile; use zero, no idea what it is in LE ...
      le_rh         = 0.0             ! (kg h2o)/(kg air)
      ! LE altitude profile; nearest cell in core grid:
      le_hh(0)    = oro(lecore_ix,lecore_iy,1)  ! m
      le_hh(1:nz) = le_hh(0) + h_m(lecore_ix,lecore_iy,:)  ! m
      le_hh(nz+1) = le_hh(nz) + 2.0e3    ! m ;  assume 2 km for aloft layer
      ! compute LE pressure half levels:
      call PotentialPressures( nz+1, le_psurf, le_temp, le_rh, le_hh, le_ph )
      if ( any(le_ph < 0.0) ) then
        write (gol,'("found negative pressures:")'); call goErr
        write (gol,*) ' le_psurf : ', le_psurf; call goErr
        write (gol,*) ' le_temp  : ', le_temp; call goErr
        write (gol,*) ' le_rh    : ', le_rh; call goErr
        write (gol,*) ' le_hh    : ', le_hh; call goErr
        write (gol,*) ' le_ph    : ', le_ph; call goErr
        TRACEBACK; status=1; return
      end if
      ! switch for original units:
      select case ( trim(bc_vcd_units) )
        !~ mass column density:
        case ( 'kg/m2' )
          ! map from TM5 profile to LE profile:
          ilast = 1
          do iz = 1, nz
            ! fill each LE layer with sum of one or more fractions of TM5 layers;
            ! use pressure axis to have mass-conservation; 
            ! use negative pressures to have increasing ax:
            call IntervalSum( -1.0*bc_ph, bc_prof, &
                              -1.0*le_ph(iz-1), -1.0*le_ph(iz), &
                             le_vcd, ilast, status )  ! (kg tracer)/m2
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target volume mixing ratio's:
              case ( 'ppb' )
                ! LE air mass per area:
                le_dm = ( le_ph(iz-1) - le_ph(iz) )/grav
                ! convert from (kg tracer)/m2 to volume mixing ratio:
                !              (kg air)/(mol air)      1
                !  (kg trc)/m2 ------------------ ----------- = (mol trc)/(mol air)
                !              (kg trc)/(mol trc) (kg air)/m2
                bc_east(iy,iz,ispec) = le_vcd * mm_air/mm / le_dm * vmr2ppb
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE layer thickness:
                le_dh = le_hh(iz) - le_hh(iz-1)
                ! convert from (kg tracer)/m2 to mass concentration (ug/m3):
                !  (kg trc)/m2 / m = (kg trc)/m3
                bc_east(iy,iz,ispec) = le_vcd / le_dh * kg2ug
              !~ not yet ...
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return                       
            end select
          end do  ! iz
        !~ input is mass concentrations, bc_prof is in (ug tr)/m3 :
        case ( 'ug/m3' )
          ! loop over levels:
          ilast = 1
          do iz = 1, nz
            ! fill model layer with integral over pressure range
            ! assuming that original values are constant within their layers;
            ! use negative pressures to have increasing ax:
            call IntervalQuad_Const( -1.0*bc_ph, bc_prof, &          ! Pa, (ug tracer)/m3
                                     -1.0*le_ph(iz-1), -1.0*le_ph(iz), &  ! Pa
                                    le_vcd, ilast, status)  ! (ug tracer)/m3*Pa
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE pressure gap:
                le_dp = le_ph(iz-1) - le_ph(iz)   ! Pa
                ! new mass concentration:
                bc_east(iy,iz,ispec) = le_vcd / le_dp  ! (ug tracer)
              !~ not yet ..
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return    
            end select 
          end do  ! levels
        !~ not yet ...
        case default
          write (gol,'("unsupported bc vcd units : ",a)') trim(bc_vcd_units); call goErr
          TRACEBACK; status=1; return
      end select
    end do  ! iy

    ! southern boundary conversion of pressure TM5 layering to km LE layering
    do ix = 1, nx
      ! most nearby LE cell index in 'normal' grid:
      lecore_ix = ix
      lecore_iy = 1
      ! LE cell index in extended le grid:
      le_ix = nh + ix
      le_iy =      1
      ! bc concentration profile:
      bc_prof = bc_vcd__legrid(le_ix,le_iy,:)
      ! bc pressure half levels:
      bc_ph = bc_phlev__legrid(le_ix,le_iy,:)
      ! LE surface pressure: copy from bc surface:
      le_psurf = bc_ph(0)
      ! LE temperature profile; nearest cell in core grid:
      le_temp(1:nz) = temp(lecore_ix,lecore_iy,1:nz)    ! K
      le_temp(nz+1) = le_temp(nz)      ! K  ; same as top level
      ! LE realitve humidity profile; use zero, no idea what it is in LE ...
      le_rh         = 0.0             ! (kg h2o)/(kg air)
      ! LE altitude profile; nearest cell in core grid:
      le_hh(0)    = oro(lecore_ix,lecore_iy,1)  ! m
      le_hh(1:nz) = le_hh(0) + h_m(lecore_ix,lecore_iy,:)  ! m
      le_hh(nz+1) = le_hh(nz) + 2.0e3    ! m ;  assume 2 km for aloft layer
      ! compute LE pressure half levels:
      call PotentialPressures( nz+1, le_psurf, le_temp, le_rh, le_hh, le_ph )
      if ( any(le_ph < 0.0) ) then
        write (gol,'("found negative pressures:")'); call goErr
        write (gol,*) ' le_psurf : ', le_psurf; call goErr
        write (gol,*) ' le_temp  : ', le_temp; call goErr
        write (gol,*) ' le_rh    : ', le_rh; call goErr
        write (gol,*) ' le_hh    : ', le_hh; call goErr
        write (gol,*) ' le_ph    : ', le_ph; call goErr
        TRACEBACK; status=1; return
      end if
      ! switch for original units:
      select case ( trim(bc_vcd_units) )
        !~ mass column density:
        case ( 'kg/m2' )
          ! map from TM5 profile to LE profile:
          ilast = 1
          do iz = 1, nz
            ! fill each LE layer with sum of one or more fractions of TM5 layers;
            ! use pressure axis to have mass-conservation; 
            ! use negative pressures to have increasing ax:
            call IntervalSum( -1.0*bc_ph, bc_prof, &
                              -1.0*le_ph(iz-1), -1.0*le_ph(iz), &
                             le_vcd, ilast, status )  ! (kg tracer)/m2
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target volume mixing ratio's:
              case ( 'ppb' )
                ! LE air mass per area:
                le_dm = ( le_ph(iz-1) - le_ph(iz) )/grav
                ! convert from (kg tracer)/m2 to volume mixing ratio:
                !              (kg air)/(mol air)      1
                !  (kg trc)/m2 ------------------ ----------- = (mol trc)/(mol air)
                !              (kg trc)/(mol trc) (kg air)/m2
                bc_south(ix,iz,ispec) = le_vcd * mm_air/mm / le_dm * vmr2ppb
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE layer thickness:
                le_dh = le_hh(iz) - le_hh(iz-1)
                ! convert from (kg tracer)/m2 to mass concentration (ug/m3):
                !  (kg trc)/m2 / m = (kg trc)/m3
                bc_south(ix,iz,ispec) = le_vcd / le_dh * kg2ug
              !~ not yet ...
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return                       
            end select
          end do  ! iz
        !~ input is mass concentrations, bc_prof is in (ug tr)/m3 :
        case ( 'ug/m3' )
          ! loop over levels:
          ilast = 1
          do iz = 1, nz
            ! fill model layer with integral over pressure range
            ! assuming that original values are constant within their layers;
            ! use negative pressures to have increasing ax:
            call IntervalQuad_Const( -1.0*bc_ph, bc_prof, &          ! Pa, (ug tracer)/m3
                                     -1.0*le_ph(iz-1), -1.0*le_ph(iz), &  ! Pa
                                    le_vcd, ilast, status)  ! (ug tracer)/m3*Pa
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE pressure gap:
                le_dp = le_ph(iz-1) - le_ph(iz)   ! Pa
                ! new mass concentration:
                bc_south(ix,iz,ispec) = le_vcd / le_dp  ! (ug tracer)
              !~ not yet ..
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return    
            end select 
          end do  ! levels
        !~ not yet ...
        case default
          write (gol,'("unsupported bc vcd units : ",a)') trim(bc_vcd_units); call goErr
          TRACEBACK; status=1; return
      end select
    enddo ! ix

    ! northern boundary conversion of pressure TM5 layering to km LE layering
    do ix = 1, nx
      ! most nearby LE cell index in 'normal' grid:
      lecore_ix = ix
      lecore_iy = ny
      ! LE cell index in extended le grid:
      le_ix = nh + ix
      le_iy =      ny+2
      ! TM5 concentration profile:
      bc_prof = bc_vcd__legrid(le_ix,le_iy,:)
      ! TM5 pressure half levels:
      bc_ph = bc_phlev__legrid(le_ix,le_iy,:)
      ! LE surface pressure: copy from bc surface:
      le_psurf = bc_ph(0)
      ! LE temperature profile; nearest cell in core grid:
      le_temp(1:nz) = temp(lecore_ix,lecore_iy,1:nz)    ! K
      le_temp(nz+1) = le_temp(nz)      ! K  ; same as top level
      ! LE realitve humidity profile; use zero, no idea what it is in LE ...
      le_rh         = 0.0             ! (kg h2o)/(kg air)
      ! LE altitude profile; nearest cell in core grid:
      le_hh(0)    = oro(lecore_ix,lecore_iy,1)  ! m
      le_hh(1:nz) = le_hh(0) + h_m(lecore_ix,lecore_iy,:)  ! m
      le_hh(nz+1) = le_hh(nz) + 2.0e3    ! m ;  assume 2 km for aloft layer
      ! compute LE pressure half levels:
      call PotentialPressures( nz+1, le_psurf, le_temp, le_rh, le_hh, le_ph )
      if ( any(le_ph < 0.0) ) then
        write (gol,'("found negative pressures:")'); call goErr
        write (gol,*) ' le_psurf : ', le_psurf; call goErr
        write (gol,*) ' le_temp  : ', le_temp; call goErr
        write (gol,*) ' le_rh    : ', le_rh; call goErr
        write (gol,*) ' le_hh    : ', le_hh; call goErr
        write (gol,*) ' le_ph    : ', le_ph; call goErr
        TRACEBACK; status=1; return
      end if
      ! switch for original units:
      select case ( trim(bc_vcd_units) )
        !~ mass column density:
        case ( 'kg/m2' )
          ! map from TM5 profile to LE profile:
          ilast = 1
          do iz = 1, nz
            ! fill each LE layer with sum of one or more fractions of TM5 layers;
            ! use pressure axis to have mass-conservation; 
            ! use negative pressures to have increasing ax:
            call IntervalSum( -1.0*bc_ph, bc_prof, &
                              -1.0*le_ph(iz-1), -1.0*le_ph(iz), &
                             le_vcd, ilast, status )  ! (kg tracer)/m2
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target volume mixing ratio's:
              case ( 'ppb' )
                ! LE air mass per area:
                le_dm = ( le_ph(iz-1) - le_ph(iz) )/grav
                ! convert from (kg tracer)/m2 to volume mixing ratio:
                !              (kg air)/(mol air)      1
                !  (kg trc)/m2 ------------------ ----------- = (mol trc)/(mol air)
                !              (kg trc)/(mol trc) (kg air)/m2
                bc_north(ix,iz,ispec) = le_vcd * mm_air/mm / le_dm * vmr2ppb
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE layer thickness:
                le_dh = le_hh(iz) - le_hh(iz-1)
                ! convert from (kg tracer)/m2 to mass concentration (ug/m3):
                !  (kg trc)/m2 / m = (kg trc)/m3
                bc_north(ix,iz,ispec) = le_vcd / le_dh * kg2ug
              !~ not yet ...
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return                       
            end select
          end do  ! iz
        !~ input is mass concentrations, bc_prof is in (ug tr)/m3 :
        case ( 'ug/m3' )
          ! loop over levels:
          ilast = 1
          do iz = 1, nz
            ! fill model layer with integral over pressure range
            ! assuming that original values are constant within their layers;
            ! use negative pressures to have increasing ax:
            call IntervalQuad_Const( -1.0*bc_ph, bc_prof, &          ! Pa, (ug tracer)/m3
                                     -1.0*le_ph(iz-1), -1.0*le_ph(iz), &  ! Pa
                                    le_vcd, ilast, status)  ! (ug tracer)/m3*Pa
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE pressure gap:
                le_dp = le_ph(iz-1) - le_ph(iz)   ! Pa
                ! new mass concentration:
                bc_north(ix,iz,ispec) = le_vcd / le_dp  ! (ug tracer)
              !~ not yet ..
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return    
            end select 
          end do  ! levels
        !~ not yet ...
        case default
          write (gol,'("unsupported bc vcd units : ",a)') trim(bc_vcd_units); call goErr
          TRACEBACK; status=1; return
      end select
    end do ! ix

    ! aloft:
    do iy = 1, ny
      do ix = 1, nx
        ! top layer:
        iz = nz + 1
        ! most nearby LE cell index in 'normal' grid:
        lecore_ix = ix
        lecore_iy = iy
        ! cell in extended le grid:
        le_ix = nh + ix
        le_iy = nh + iy
        ! bc concentration profile:
        bc_prof = bc_vcd__legrid(le_ix,le_iy,:)
        ! bc pressure half levels:
        bc_ph = bc_phlev__legrid(le_ix,le_iy,:)
        ! LE surface pressure: copy from bc surface:
        le_psurf = bc_ph(0)
        ! LE temperature profile; nearest cell in core grid:
        le_temp(1:nz) = temp(lecore_ix,lecore_iy,1:nz)    ! K
        le_temp(nz+1) = le_temp(nz)      ! K  ; same as top level
        ! LE realitve humidity profile; use zero, no idea what it is in LE ...
        le_rh         = 0.0             ! (kg h2o)/(kg air)
        ! LE altitude profile; nearest cell in core grid:
        le_hh(0)    = oro(lecore_ix,lecore_iy,1)  ! m
        le_hh(1:nz) = le_hh(0) + h_m(lecore_ix,lecore_iy,:)  ! m
        le_hh(nz+1) = le_hh(nz) + 2.0e3    ! m ;  assume 2 km for aloft layer
        ! compute LE pressure half levels:
        call PotentialPressures( nz+1, le_psurf, le_temp, le_rh, le_hh, le_ph )
        if ( any(le_ph < 0.0) ) then
          write (gol,'("found negative pressures:")'); call goErr
          write (gol,*) ' le_psurf : ', le_psurf; call goErr
          write (gol,*) ' le_temp  : ', le_temp; call goErr
          write (gol,*) ' le_rh    : ', le_rh; call goErr
          write (gol,*) ' le_hh    : ', le_hh; call goErr
          write (gol,*) ' le_ph    : ', le_ph; call goErr
          TRACEBACK; status=1; return
        end if
        ! switch for original units:
        select case ( trim(bc_vcd_units) )
          !~ mass column density:
          case ( 'kg/m2' )
            ! fill each LE layer with sum of one or more fractions of TM5 layers;
            ! use pressure axis to have mass-conservation; 
            ! use negative pressures to have increasing ax:
            call IntervalSum( -1.0*bc_ph, bc_prof, &
                              -1.0*le_ph(iz-1), -1.0*le_ph(iz), &
                             le_vcd, ilast, status )  ! (kg tracer)/m2
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target volume mixing ratio's:
              case ( 'ppb' )
                ! LE air mass per area:
                le_dm = ( le_ph(iz-1) - le_ph(iz) )/grav
                ! convert from (kg tracer)/m2 to volume mixing ratio:
                !              (kg air)/(mol air)      1
                !  (kg trc)/m2 ------------------ ----------- = (mol trc)/(mol air)
                !              (kg trc)/(mol trc) (kg air)/m2
                caloft(ix,iy,ispec) = le_vcd * mm_air/mm / le_dm * vmr2ppb
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE layer thickness:
                le_dh = le_hh(iz) - le_hh(iz-1)
                ! convert from (kg tracer)/m2 to mass concentration (ug/m3):
                !  (kg trc)/m2 / m = (kg trc)/m3
                caloft(ix,iy,ispec) = le_vcd / le_dh * kg2ug
              !~ not yet ...
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return                       
            end select
          !~ input is mass concentrations, bc_prof is in (ug tr)/m3 :
          case ( 'ug/m3' )
            ! fill model layer with integral over pressure range
            ! assuming that original values are constant within their layers;
            ! use negative pressures to have increasing ax:
            call IntervalQuad_Const( -1.0*bc_ph, bc_prof, &          ! Pa, (ug tracer)/m3
                                     -1.0*le_ph(iz-1), -1.0*le_ph(iz), &  ! Pa
                                    le_vcd, ilast, status)  ! (ug tracer)/m3*Pa
            IF_NOTOK_RETURN(status=1)
            ! switch for target units:
            select case ( specunit(ispec) )
              !~ target mass concentrations:
              case ( 'ug/m3' )
                ! LE pressure gap:
                le_dp = le_ph(iz-1) - le_ph(iz)   ! Pa
                ! new mass concentration:
                caloft(ix,iy,ispec) = le_vcd / le_dp  ! (ug tracer)
              !~ not yet ..
              case default
                write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                TRACEBACK; status=1; return    
            end select 
          !~ not yet ...
          case default
            write (gol,'("unsupported bc vcd units : ",a)') trim(bc_vcd_units); call goErr
            TRACEBACK; status=1; return
        end select
      end do  ! ix
    end do  ! iy
    
    ! clear:
    deallocate( bc_dm )
    deallocate( bc_vcd__legrid )
    deallocate( bc_phlev__legrid )
    deallocate( bc_ahlev__legrid )
    deallocate( bc_prof )
    deallocate( bc_ph )

    ! done with bc grid:
    call Done( bc_lli, status )
    IF_NOTOK_RETURN(status=1)
    
    ! done with le grid:
    call Done( le_lli, status )
    IF_NOTOK_RETURN(status=1)

    deallocate(target2source_i_l)
    deallocate(target2source_i_r)
    deallocate(target2source_j_l)
    deallocate(target2source_j_r)       
    deallocate(fracxt)
        
    ! ok
    status = 0

  end subroutine LE_Bound_Fill
  
  

  ! ===============================================================
  

  subroutine LE_Bound_Fill_Initial( c, ispec, mm, mm_air, &
                            bc_lons, bc_lats, bc_phlev, bc_vmr, bc_unit, &
                            status )
  
    use Binas, only : grav
    use JAQL , only : PotentialPressures
    use Num  , only : IntervalSum, IntervalQuad_Const
    use Grid , only : TllGridInfo, Init, Done, FillGrid_AreaAverage

    use dims   , only : runF
    use dims   , only : nx, ny, nz
    use indices, only : specname, specunit
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Grid, only : lli

    ! --- in/out ------------------------------
    
    real, intent(out)               ::  c(nx,ny,nz)
    integer, intent(in)             ::  ispec
    real, intent(in)                ::  mm, mm_air ! molemass
    real, intent(in)                ::  bc_lons(:)
    real, intent(in)                ::  bc_lats(:)
    real, intent(in)                ::  bc_phlev(:,:,:)  ! levels 1:bc_nlay+1
    real, intent(in)                ::  bc_vmr(:,:,:)
    character(len=*), intent(in)    ::  bc_unit
    integer, intent(out)            ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Bound_Fill_Initial'
    
    ! conversion factors:
    real, parameter   ::  vmr2ppb = 1.0e9   ! mole/mole -> mole/(1e9 mol)
    real, parameter   ::  kg2ug   = 1.0e9   ! kg -> ug
    
    ! --- local -------------------------------
    
    integer               ::  bc_nlon, bc_nlat, bc_nlay
    type(TllGridInfo)     ::  bc_lli
    real, allocatable     ::  bc_dm(:,:)
    real, allocatable     ::  bc_vcd__legrid(:,:,:)  ! vertical column density TM5 (kg tracer)/m2
    real, allocatable     ::  bc_phlev__legrid(:,:,:)  ! bc half level pressure at LE grid coordinates
    integer               ::  l
    real                  ::  le_psurf, le_temp(nz+1), le_rh(nz+1)  ! pressure at surface, temperature and relative humidity
    real                  ::  le_hh(0:nz+1), le_ph(0:nz+1) 
    real                  ::  le_dm
    real                  ::  le_dh
    real                  ::  le_dp
    real                  ::  le_vcd                  ! vertical column density LOTOS-EUROS (kg tracer)/m2
    real, allocatable     ::  bc_ph(:), bc_prof(:)  ! TM5 pressure and concentration profile 
    integer               ::  ix, iy, iz
    integer               ::  lecore_ix, lecore_iy
    integer               ::  le_ix, le_iy
    integer               ::  ilast
    character(len=32)     ::  bc_vcd_units
    
    ! meteo data:
    real, pointer        ::  oro(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  h_m(:,:,:)   ! (lon,lat,nz)
    real, pointer        ::  temp(:,:,:)   ! (lon,lat,nz)    
    
    ! --- begin -------------------------------
    
    ! point to meteo data:
    call LE_Data_GetPointer( 'oro', oro, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status,  check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 't', temp, status,  check_units ='K' )
    IF_NOTOK_RETURN(status=1)

    ! check unit ...
    select case ( trim(bc_unit) )
      case ( '1 = mole/mo', 'mole mole-1', 'mole mole**-1', 'kg kg**-1' )
        ! converted later on ..
      case ( 'ug/m3' )
        ! conversion is different
      case default
        write (gol,'("unsupported bc unit : ",a)') trim(bc_unit); call goErr
        TRACEBACK; status=1; return
    end select

    ! grid sizes:
    bc_nlon = size(bc_vmr,1)
    bc_nlat = size(bc_vmr,2)
    bc_nlay = size(bc_vmr,3)
    
    ! check ...
    if ( (size(bc_lons) /= bc_nlon) .or. &
         (size(bc_lats) /= bc_nlat) .or. &
         any(shape(bc_phlev) /= (/bc_nlon,bc_nlat,bc_nlay+1/)) ) then
      write (gol,'("one or more arguments have unexpected shape:")'); call goErr
      write (gol,'("  bc_lons   : ",i6)') shape(bc_lons); call goErr
      write (gol,'("  bc_lats   : ",i6)') shape(bc_lats); call goErr
      write (gol,'("  bc_phlev  : ",3i6)') shape(bc_phlev); call goErr
      write (gol,'("  bc_vmr    : ",3i6)') shape(bc_vmr); call goErr
      TRACEBACK; status=1; return
    end if

    ! define bc grid:
    call Init( bc_lli, bc_lons(1), bc_lons(2)-bc_lons(1), bc_nlon, &
                       bc_lats(1), bc_lats(2)-bc_lats(1), bc_nlat, status )
    IF_NOTOK_RETURN(status=1)

    ! increasing pressure ax (top->down!) ?
    if ( bc_phlev(1,1,1) < bc_phlev(1,1,bc_nlay+1) ) then
      write (gol,'("boundary condition field defined on increasing pressure ax not supported yet ...")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( bc_dm(1:bc_nlon,bc_nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_vcd__legrid(1:lli%nlon,1:lli%nlat,1:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_phlev__legrid(1:lli%nlon,1:lli%nlat,0:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_ph(0:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( bc_prof(1:bc_nlay), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! convert TM5 half level pressures to LE horizontal grid;
    ! note difference in indices: 1:n+1 on input, 0:n on output:
    do l = 0, bc_nlay
      ! area weighted average:
      call FillGrid_AreaAverage(    lli, bc_phlev__legrid(:,:,l), &  ! Pa
                                 bc_lli, bc_phlev(:,:,l+1), status )
       IF_NOTOK_RETURN(status=1)
    end do
    
    ! convert TM5 field to LE horizontal grid;
    do l = 1, bc_nlay
      select case ( trim(bc_unit) )
        !
        case ( '1 = mole/mo', 'mole mole-1', 'mole mole**-1' )
          ! mass per aera; dimension bc_phlev(:,:,1:bc_nlay+1)
          bc_dm = abs( bc_phlev(:,:,l+1) - bc_phlev(:,:,l) )/grav  ! (kg air)/m2
          ! concentrations in volume mixing ratio;
          ! convert to masses; use that volume mixing ratio is mol mixing ratio:
          !  (mol trc) (kg trc)/(mol trc)
          !  --------- ------------------ (kg air)/m2 = (kg trc)/m2
          !  (mol air) (kg air)/(mol air)
          ! area weighted average:
          call FillGrid_AreaAverage(    lli, bc_vcd__legrid(:,:,l), &  ! (kg tracer)/m2
                                   bc_lli, bc_vmr(:,:,l)*mm/mm_air*bc_dm, status )
          IF_NOTOK_RETURN(status=1) 
          ! set units:
          bc_vcd_units = 'kg/m2'
        !
        case ( 'kg kg**-1' )
          ! mass per aera; dimension bc_phlev(:,:,1:bc_nlay+1)
          bc_dm = abs( bc_phlev(:,:,l+1) - bc_phlev(:,:,l) )/grav  ! (kg air)/m2
          ! concentrations in mass mixing ratio; convert to masses:
          !  (kg trc)
          !  -------- (kg air)/m2 = (kg trc)/m2
          !  (kg air)
          ! (variable is named 'bc_vmr', but with this units it is actually 'bc_mmr')
          ! area weighted average:
          call FillGrid_AreaAverage(    lli, bc_vcd__legrid(:,:,l), &  ! (kg tracer)/m2
                                   bc_lli, bc_vmr(:,:,l)*bc_dm, status )
          IF_NOTOK_RETURN(status=1) 
          ! set units:
          bc_vcd_units = 'kg/m2'
        !
        case ( 'ug/m3' )
          ! area weighted average, units remain ug/m3 ;
          ! (name 'bc_vcd' is not correct therefore ...)
          call FillGrid_AreaAverage(    lli, bc_vcd__legrid(:,:,l), &  ! (ug tracer)/m3
                                     bc_lli, bc_vmr(:,:,l), status )
          IF_NOTOK_RETURN(status=1)  
          ! set units:
          bc_vcd_units = 'ug/m3'
        !
        case default
          write (gol,'("unsupported bc unit : ",a)') trim(bc_unit); call goErr
          TRACEBACK; status=1; return
      end select

      
    end do

    ! loop over grid cells:
    do iy = 1, ny
      do ix = 1, nx
        ! most nearby LE cell index in 'normal' grid:
        lecore_ix = ix
        lecore_iy = iy
        ! cell in le grid:
        le_ix = ix
        le_iy = iy
        ! bc concentration profile:
        bc_prof = bc_vcd__legrid(le_ix,le_iy,:)   ! (kg tr)/m2 or (ug tr)/m3
        ! bc pressure half levels:
        bc_ph = bc_phlev__legrid(le_ix,le_iy,:)
        ! LE surface pressure: copy from bc surface:
        le_psurf = bc_ph(0)
        ! LE temperature profile; nearest cell in core grid:
        le_temp(1:nz) = temp(lecore_ix,lecore_iy,1:nz)    ! K
        le_temp(nz+1) = le_temp(nz)      ! K  ; same as top level
        ! LE realitve humidity profile; use zero, no idea what it is in LE ...
        le_rh         = 0.0             ! (kg h2o)/(kg air)
        ! LE altitude profile; nearest cell in core grid:
        le_hh(0)    = oro(lecore_ix,lecore_iy,1)  ! m
        le_hh(1:nz) = le_hh(0) + h_m(lecore_ix,lecore_iy,:)  ! m
        le_hh(nz+1) = le_hh(nz) + 2.0e3    ! m ;  assume 2 km for aloft layer
        ! compute LE pressure half levels:
        call PotentialPressures( nz+1, le_psurf, le_temp, le_rh, le_hh, le_ph )
        if ( any(le_ph < 0.0) ) then
          write (gol,'("found negative pressures:")'); call goErr
          write (gol,*) ' le_psurf : ', le_psurf; call goErr
          write (gol,*) ' le_temp  : ', le_temp; call goErr
          write (gol,*) ' le_rh    : ', le_rh; call goErr
          write (gol,*) ' le_hh    : ', le_hh; call goErr
          write (gol,*) ' le_ph    : ', le_ph; call goErr
          TRACEBACK; status=1; return
        end if
        ! switch for original units:
        select case ( trim(bc_vcd_units) )
          !~ mass column density:
          case ( 'kg/m2' )
            ! loop over levels:
            ilast = 1
            do iz = 1, nz
              ! fill model layer with sum of one or more fractions of bc layers;
              ! use pressure axis to have mass-conservation; 
              ! use negative pressures to have increasing ax:
              call IntervalSum( -1.0*bc_ph, bc_prof, &          ! Pa, (kg tracer)/m2
                                -1.0*le_ph(iz-1), -1.0*le_ph(iz), &  ! Pa
                               le_vcd, ilast, status ) ! (kg tracer)/m2
              IF_NOTOK_RETURN(status=1)
              ! switch for target units:
              select case ( specunit(ispec) )
                !~ target volume mixing ratio's:
                case ( 'ppb' )
                  ! LE air mass per area:
                  le_dm = ( le_ph(iz-1) - le_ph(iz) )/grav
                  ! convert from (kg tracer)/m2 to volume mixing ratio:
                  !              (kg air)/(mol air)      1
                  !  (kg trc)/m2 ------------------ ----------- = (mol trc)/(mol air)
                  !              (kg trc)/(mol trc) (kg air)/m2
                  c(ix,iy,iz) = le_vcd * mm_air/mm / le_dm * vmr2ppb
                !~ target mass concentrations:
                case ( 'ug/m3' )
                  ! LE layer thickness:
                  le_dh = le_hh(iz) - le_hh(iz-1)
                  ! convert from (kg tracer)/m2 to mass concentration (ug/m3):
                  !  (kg trc)/m2 / m = (kg trc)/m3
                  c(ix,iy,iz) = le_vcd / le_dh * kg2ug
                !~ not yet ...
                case default
                  write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                  TRACEBACK; status=1; return
              end select
            end do  ! iz
          !~ input is mass concentrations too, bc_prof is in (ug tr)/m3 too:
          case ('ug/m3')
            ! loop over levels:
            ilast = 1
            do iz = 1, nz
              ! fill model layer with integral over pressure range
              ! assuming that original values are constant within their layers;
              ! use negative pressures to have increasing ax:
              call IntervalQuad_Const( -1.0*bc_ph, bc_prof, &          ! Pa, (ug tracer)/m3
                                       -1.0*le_ph(iz-1), -1.0*le_ph(iz), &  ! Pa
                                      le_vcd, ilast, status)  ! (ug tracer)/m3*Pa
              IF_NOTOK_RETURN(status=1)
              ! switch for target units:
              select case ( specunit(ispec) )
                !~ target mass concentrations:
                case ( 'ug/m3' )
                  ! LE pressure gap:
                  le_dp = le_ph(iz-1) - le_ph(iz)   ! Pa
                  ! new mass concentration:
                  c(ix,iy,iz) = le_vcd / le_dp  ! (ug tracer)/m3
                !~ not yet ...
                case default
                  write (gol,'("unsupported target unit : ",a)') trim(specunit(ispec)); call goErr
                  TRACEBACK; status=1; return   
              end select 
            end do  ! iz
          !~ not yet ..
          case default
            write (gol,'("unsupported bc vcd units : ",a)') trim(bc_vcd_units); call goErr
            TRACEBACK; status=1; return    
        end select 
      end do  ! ix
    end do  ! iy

    ! clear:
    deallocate( bc_dm )
    deallocate( bc_vcd__legrid )
    deallocate( bc_phlev__legrid )
    deallocate( bc_prof )
    deallocate( bc_ph )

    ! done with bc grid:
    call Done( bc_lli, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Bound_Fill_Initial
  
  

  ! ====================================================================
  ! ===
  ! === pressure, altitudes
  ! ===
  ! ====================================================================
  
  
  ! Given 3D half level pressures (1=surface,nz+1=top), 
  ! temperature, and specific humidity, estimate half level altitudes.
  ! Surface altitude is guessed from the gap between surface pressure
  ! and an assumed sea-level pressures, which is set to the maximu
  ! value in the surface pressures.
  
  subroutine PTQ_to_H( phlev, temper, shumid, ahlev, status )

    use Binas, only : p_sealevel
    use JAQL , only : PotentialHeight
  
    ! --- in/out ---------------------------------
    
    real, intent(in)        ::  phlev (:,:,:)   ! (nlon,nlat,nlev+1)  Pa
    real, intent(in)        ::  temper(:,:,:)   ! (nlon,nlat,nlev  )  K
    real, intent(in)        ::  shumid(:,:,:)   ! (nlon,nlat,nlev  )  kg/kg
    real, intent(out)       ::  ahlev(:,:,:)    ! (nlon,nlat,nlev+1)  m a.s.l.
    integer, intent(out)    ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/PTQ_to_H'
    
    ! --- local ----------------------------------
    
    real                      ::  psl
    integer                   ::  nlon, nlat, nlev
    integer                   ::  ilev
    real, allocatable         ::  dh(:,:)

    ! --- begin ----------------------------------
    
    ! dims:
    nlon = size(phlev,1)
    nlat = size(phlev,2)
    nlev = size(phlev,3) - 1
    
    ! storage:
    allocate( dh(nlon,nlat) )

    ! sea level pressure, ensure that it is at least the 
    ! highest pressure that occures in the input:
    psl = max( p_sealevel, maxval(phlev(:,:,1)) )

    ! guess surface altitude from pressure gap with ;
    call PotentialHeight( phlev(:,:,1), psl, &
                           temper(:,:,1), shumid(:,:,1), &
                           dh )
    ! fill:
    ahlev(:,:,1) = dh

    ! loop over layers:
    do ilev = 1, nlev
      ! height step in m for this layer;
      ! use minimum pressure at top since zero pressure 
      ! only reached at infinite distance:
      call PotentialHeight( max(0.1,phlev(:,:,ilev+1)), phlev(:,:,ilev), &
                             temper(:,:,ilev), shumid(:,:,ilev), &
                             dh )
      ! add to previous altitude:
      ahlev(:,:,ilev+1) = ahlev(:,:,ilev) + dh
    end do
    
    ! clear:
    deallocate( dh )
    
    ! ok
    status = 0
    
  end subroutine PTQ_to_H


!  ! ====================================================================
!  ! ===
!  ! === mapping
!  ! ===
!  ! ====================================================================
!
!
!  subroutine Heights_to_LE_Init( HtL, Lon_Axis_in, Lat_Axis_in, Lev_Axis_in, status )
!
!    use CF     , only : T_CF_Axis
!    use CF     , only : T_CF_CellAxis, CF_CellAxis_Init, CF_CellAxis_Done
!    use CF     , only : T_CF_AxisMapper, CF_AxisMapper_Init
!    use Dims   , only : nx, ny, nz
!    use LE_Data      , only : LE_Data_GetPointer
!    use LE_Grid, only : Lon_Axis, west_Lon_Axis, east_Lon_Axis, south_Lon_Axis, north_Lon_Axis
!    use LE_Grid, only : Lat_Axis, west_Lat_Axis, east_Lat_Axis, south_Lat_Axis, north_Lat_Axis
!
!    ! --- in/out -----------------------------------------
!    
!    type(T_Heights_to_LE), intent(out)    ::  HtL
!    type(T_CF_Axis), intent(in)           ::  Lon_Axis_in
!    type(T_CF_Axis), intent(in)           ::  Lat_Axis_in
!    type(T_CF_Axis), intent(in)           ::  Lev_Axis_in
!    integer, intent(out)                  ::  status
!    
!    ! --- const --------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/Heights_to_LE_Init'
!    
!    ! --- local -------------------------------
!    
!    type(T_CF_CellAxis)       ::  ML_Axis    ! model level axis
!    real, allocatable         ::  hhb(:,:,:)  ! (nx,ny,0:nz)  half level heights [m]
!    integer                   ::  i, j
!    ! meteo data:
!    real, pointer               ::  h_m(:,:,:)   ! (lon,lat,nz)
!      
!    ! --- begin -------------------------------
!
!    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m')    
!    IF_NOTOK_RETURN(status=1)
!  
!    ! fill dimensions:
!    call Lon_Axis_in%Get( status, length=HtL%nx_in )
!    IF_NOTOK_RETURN(status=1)
!    call Lat_Axis_in%Get( status, length=HtL%ny_in )
!    IF_NOTOK_RETURN(status=1)
!    call Lev_Axis_in%Get( status, length=HtL%nz_in )
!    IF_NOTOK_RETURN(status=1)
!  
!    ! setup conversion objects for horizontal axes:
!    call CF_AxisMapper_Init( HtL%Lon_Mapper, Lon_Axis_in, Lon_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Init( HtL%Lat_Mapper, Lat_Axis_in, Lat_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! setup conversion objects for horizontal axes:
!    call CF_AxisMapper_Init( HtL%west_Lon_Mapper, Lon_Axis_in, west_Lon_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Init( HtL%west_Lat_Mapper, Lat_Axis_in, west_Lat_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    ! setup conversion objects for horizontal axes:
!    call CF_AxisMapper_Init( HtL%east_Lon_Mapper, Lon_Axis_in, east_Lon_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Init( HtL%east_Lat_Mapper, Lat_Axis_in, east_Lat_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    ! setup conversion objects for horizontal axes:
!    call CF_AxisMapper_Init( HtL%south_Lon_Mapper, Lon_Axis_in, south_Lon_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Init( HtL%south_Lat_Mapper, Lat_Axis_in, south_Lat_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    ! setup conversion objects for horizontal axes:
!    call CF_AxisMapper_Init( HtL%north_Lon_Mapper, Lon_Axis_in, north_Lon_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Init( HtL%north_Lat_Mapper, Lat_Axis_in, north_Lat_Axis, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! storage for half level heights:
!    allocate( hhb(nx,ny,0:nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    ! fill:
!    hhb(:,:,0) = 0.0
!    hhb(:,:,1:nz) = h_m   ! m
!    
!    ! vertical mappers per cell:
!    allocate( HtL%Lev_Mapper(nx,ny), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    ! loop over cells:
!    do j = 1, ny
!      do i = 1, nx
!        ! define axis:
!        call CF_CellAxis_Init( ML_Axis, 'layer', 'height', nz, status, &
!                                 bvalues=hhb(i,j,:), axis='Z' )
!        IF_NOTOK_RETURN(status=1)
!        ! define mapping:
!        call CF_AxisMapper_Init( HtL%Lev_Mapper(i,j), Lev_Axis_in, ML_Axis, status )
!        IF_NOTOK_RETURN(status=1)
!        ! done:
!        call CF_CellAxis_Done( ML_Axis, status )
!        IF_NOTOK_RETURN(status=1)
!      end do
!    end do
!    
!    ! vertical mappers per aloft cell:
!    allocate( HtL%aloft_Lev_Mapper(nx,ny), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    ! loop over cells:
!    do j = 1, ny
!      do i = 1, nx
!        ! define axis, layer of 2km above top:
!        call CF_CellAxis_Init( ML_Axis, 'layer', 'height', 1, status, &
!                                 first_bound=hhb(i,j,nz), last_bound=hhb(i,j,nz)+2e3, &
!                                 axis='Z' )
!        IF_NOTOK_RETURN(status=1)
!        ! define mapping:
!        call CF_AxisMapper_Init( HtL%aloft_Lev_Mapper(i,j), Lev_Axis_in, ML_Axis, status )
!        IF_NOTOK_RETURN(status=1)
!        ! done:
!        call CF_CellAxis_Done( ML_Axis, status )
!        IF_NOTOK_RETURN(status=1)
!      end do
!    end do
!
!    ! clear:
!    deallocate( hhb, stat=status )
!    IF_NOTOK_RETURN(status=1)
!
!    ! ok
!    status = 0
!    
!  end subroutine Heights_to_LE_Init
!
!
!  ! ***
!
!
!  subroutine Heights_to_LE_Done( HtL, status )
!
!    use CF     , only : CF_AxisMapper_Done
!    use Dims   , only : nx, ny
!      
!    ! --- in/out -----------------------------------------
!    
!    type(T_Heights_to_LE), intent(inout)    ::  HtL
!    integer, intent(out)                    ::  status
!    
!    ! --- const --------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/Heights_to_LE_Done'
!    
!    ! --- local -------------------------------
!    
!    integer                   ::  i, j
!      
!    ! --- begin -------------------------------
!    
!    ! done with mapper:
!    call CF_AxisMapper_Done( HtL%Lon_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Done( HtL%Lat_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! done with mapper:
!    call CF_AxisMapper_Done( HtL%west_Lon_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Done( HtL%west_Lat_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    ! done with mapper:
!    call CF_AxisMapper_Done( HtL%east_Lon_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Done( HtL%east_Lat_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    ! done with mapper:
!    call CF_AxisMapper_Done( HtL%south_Lon_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Done( HtL%south_Lat_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    ! done with mapper:
!    call CF_AxisMapper_Done( HtL%north_Lon_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    call CF_AxisMapper_Done( HtL%north_Lat_Mapper, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! loop over cells:
!    do j = 1, ny
!      do i = 1, nx
!        ! done with aloft mapping:
!        call CF_AxisMapper_Done( HtL%Lev_Mapper(i,j), status )
!        IF_NOTOK_RETURN(status=1)
!      end do
!    end do
!    ! vertical mappers per cell:
!    deallocate( HtL%Lev_Mapper, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! loop over cells:
!    do j = 1, ny
!      do i = 1, nx
!        ! done with aloft mapping:
!        call CF_AxisMapper_Done( HtL%aloft_Lev_Mapper(i,j), status )
!        IF_NOTOK_RETURN(status=1)
!      end do
!    end do
!    ! vertical mappers per aloft cell:
!    deallocate( HtL%aloft_Lev_Mapper, stat=status )
!    IF_NOTOK_RETURN(status=1)
!
!    ! ok
!    status = 0
!    
!  end subroutine Heights_to_LE_Done
!
!
!  ! ***
!
!
!  ! 3D conversion factors from 'units_in' to spec units
!
!  subroutine LE_Bound_Conversion_Field( conv, units_in, ispec, status )
!
!    use Binas  , only : xm_air
!    use Indices, only : specname, specunit, specmolm
!
!    use LE_Data               , only : LE_Data_GetPointer
!
!    ! --- in/out -----------------------------------------
!    
!    real, intent(out)                       ::  conv(:,:,:)  ! (nx,ny,nz)
!    character(len=*), intent(in)            ::  units_in
!    integer, intent(in)                     ::  ispec
!    integer, intent(out)                    ::  status
!    
!    ! --- const --------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/Heights_to_LE_Apply_Bounds'
!    
!    ! --- local -------------------------------
!    
!    ! meteo data:
!    real, pointer       ::  dens(:,:,:)   ! (lon,lat,alt)            
!    
!    ! --- begin -------------------------------
!
!    ! point to meteo data:
!    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')            
!    IF_NOTOK_RETURN(status=1)
!
!    ! check input units:
!    select case ( trim(units_in) )
!      ! mass concentrations:
!      case ( 'kg m-3' )
!        ! target units:
!        select case ( trim(specunit(ispec)) )
!          ! mass concentrations:
!          case ( 'ug/m3' )
!            ! constant factor from kg -> ug :
!            conv = 1.0e9  ! ug/kg
!          ! volume mixing ratios:
!          case ( 'ppb' )
!            ! unit conversion from [kg m-3] to [ppb]:
!            ! (kg tracer)/(m3 air) / (kg tracer)/(mol tracer) / (kg air)/(m3 air) * (kg air)/(mol air) * 1e9
!            !       = 1e-9 (mol tracer)/(mol air) = ppb
!            conv = 1e9 / specmolm(ispec) * xm_air / dens  ! ppb/(kg/m3)
!          ! unknown ...
!          case default
!            write (gol,'("conversion of `",a,"` from `",a,"` to `",a,"` not implemented")') &
!                     trim(specname(ispec)), trim(units_in), trim(specunit(ispec)); call goErr
!            TRACEBACK; status=1; return
!        end select
!      ! unknown ...
!      case default
!        write (gol,'("conversion of `",a,"` from `",a,"` not implemented")') &
!                     trim(specname(ispec)), trim(units_in); call goErr
!        TRACEBACK; status=1; return
!    end select
!
!    ! ok
!    status = 0
!    
!  end subroutine LE_Bound_Conversion_Field
!  
!  
!  ! *
!    
!
!  subroutine Heights_to_LE_Apply_Bounds( HtL, bc_conc, bc_conc_units, ispec, spec_filled, status )
!
!    use Dims   , only : nx, ny, nz
!    use Dims   , only : bc_west, bc_east, bc_south, bc_north
!    use Dims   , only : caloft
!    use LE_Data               , only : LE_Data_GetPointer
!
!    ! --- in/out -----------------------------------------
!    
!    type(T_Heights_to_LE), intent(inout)    ::  HtL
!    real, intent(in)                        ::  bc_conc(:,:,:)
!    character(len=*), intent(in)            ::  bc_conc_units
!    integer, intent(in)                     ::  ispec
!    logical, intent(inout)                  ::  spec_filled(:)  ! (ispec) .true. if specie becomes filled
!    integer, intent(out)                    ::  status
!    
!    ! --- const --------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/Heights_to_LE_Apply_Bounds'
!    
!    ! --- local -------------------------------
!    
!    real, allocatable   ::  conv(:,:,:)
!    real, allocatable   ::  slab(:,:,:)
!    real, allocatable   ::  alfa(:,:,:)
!
!    ! meteo data:
!    real, pointer       ::  dens(:,:,:)   ! (lon,lat,alt)            
!    
!    ! --- begin -------------------------------
!
!    ! point to meteo data:
!    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')            
!    IF_NOTOK_RETURN(status=1)
!
!    ! Conversion array defined on core grid; values at edges are
!    ! used to convert boundary conditions.
!    ! storage:
!    allocate( conv(nx,ny,nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    ! fill:
!    call LE_Bound_Conversion_Field( conv, bc_conc_units, ispec, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! *
!    
!    ! east/west slab:
!    allocate( slab(1:1,1:ny,1:nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( alfa(1:1,1:ny,1:nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! fill boundary array:
!    call LE_Bound_Map_BC_to_Model( HtL%west_Lon_Mapper, HtL%west_Lat_Mapper, &
!                           HtL%Lev_Mapper, bc_conc, slab, alfa, status )
!    IF_NOTOK_RETURN(status=1)
!    ! add converted contribution:
!    call LE_Bound_Add_2D( ispec, bc_west, dens(1,:,:), alfa(1,:,:), slab(1,:,:)*conv(1,:,:), &
!                             spec_filled, status )
!    IF_NOTOK_RETURN(status=1)
!  
!    ! fill boundary array:
!    call LE_Bound_Map_BC_to_Model( HtL%east_Lon_Mapper, HtL%east_Lat_Mapper, &
!                           HtL%Lev_Mapper, bc_conc, slab, alfa, status )
!    IF_NOTOK_RETURN(status=1)
!    ! add converted contribution:
!    call LE_Bound_Add_2D( ispec, bc_east, dens(nx,:,:), alfa(1,:,:), slab(1,:,:)*conv(nx,:,:), spec_filled, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! clear:
!    deallocate( slab, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    deallocate( alfa, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! *
!    
!    ! north/south slab:
!    allocate( slab(1:nx,1:1,1:nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( alfa(1:nx,1:1,1:nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!  
!    ! fill boundary array:
!    call LE_Bound_Map_BC_to_Model( HtL%south_Lon_Mapper, HtL%south_Lat_Mapper, &
!                           HtL%Lev_Mapper, bc_conc, slab, alfa, status )
!    IF_NOTOK_RETURN(status=1)
!    ! add converted contribution:
!    call LE_Bound_Add_2D( ispec, bc_south, dens(:, 1,:), alfa(:,1,:), slab(:,1,:)*conv(:,1,:), spec_filled, status )
!    IF_NOTOK_RETURN(status=1)
!  
!    ! fill boundary array:
!    call LE_Bound_Map_BC_to_Model( HtL%north_Lon_Mapper, HtL%north_Lat_Mapper, &
!                           HtL%Lev_Mapper, bc_conc, slab, alfa, status )
!    IF_NOTOK_RETURN(status=1)
!    ! add converted contribution:
!    call LE_Bound_Add_2D( ispec, bc_north, dens(:,ny,:), alfa(:,1,:), slab(:,1,:)*conv(:,ny,:), spec_filled, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! clear:
!    deallocate( slab, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    deallocate( alfa, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! *
!    
!    ! aloft slab:
!    allocate( slab(1:nx,1:ny,1:1), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( alfa(1:nx,1:ny,1:1), stat=status )
!    IF_NOTOK_RETURN(status=1)
!  
!    ! fill aloft array:
!    call LE_Bound_Map_BC_to_Model( HtL%Lon_Mapper, HtL%Lat_Mapper, &
!                           HtL%aloft_Lev_Mapper, bc_conc, slab, alfa, status )
!    IF_NOTOK_RETURN(status=1)
!    ! add converted contribution:
!    call LE_Bound_Add_2D( ispec, caloft, dens(:,:,nz), alfa(:,:,1), slab(:,:,1)*conv(:,:,1), spec_filled, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! clear:
!    deallocate( slab, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    deallocate( alfa, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! *
!    
!    ! clear:
!    deallocate( conv, stat=status )
!    IF_NOTOK_RETURN(status=1)
!
!    ! ok
!    status = 0
!    
!  end subroutine Heights_to_LE_Apply_Bounds
!  
!  
!  ! *
!  
!
!  subroutine Heights_to_LE_Apply_Initial( HtL, bc_conc, bc_conc_units, ispec, c, spec_filled, status )
!
!    use Dims   , only : nx, ny, nz
!    use Dims   , only : bc_west, bc_east, bc_south, bc_north
!    use Dims   , only : caloft
!    use LE_Data, only : LE_Data_GetPointer
!
!    ! --- in/out -----------------------------------------
!    
!    type(T_Heights_to_LE), intent(inout)    ::  HtL
!    real, intent(in)                        ::  bc_conc(:,:,:)
!    character(len=*), intent(in)            ::  bc_conc_units
!    integer, intent(in)                     ::  ispec
!    real, intent(inout)                     ::  c(:,:,:,:)   ! (nx,ny,nz,nspec)
!    logical, intent(inout)                  ::  spec_filled(:)  ! (ispec) .true. if specie becomes filled
!    integer, intent(out)                    ::  status
!    
!    ! --- const --------------------------------
!    
!    character(len=*), parameter   ::  rname = mname//'/Heights_to_LE_Apply_Initial'
!    
!    ! --- local -------------------------------
!    
!    real, allocatable   ::  conv(:,:,:)
!    real, allocatable   ::  slab(:,:,:)
!    real, allocatable   ::  alfa(:,:,:)
!
!    ! meteo data:
!    real, pointer       ::  dens(:,:,:)   ! (lon,lat,alt)            
!    
!    ! --- begin -------------------------------
!
!    ! point to meteo data:
!    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')            
!    IF_NOTOK_RETURN(status=1)
!
!    ! Conversion array defined on core grid; values at edges are
!    ! used to convert boundary conditions.
!    ! 3D conversion array:
!    allocate( conv(nx,ny,nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    ! fill:
!    call LE_Bound_Conversion_Field( conv, bc_conc_units, ispec, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! *
!    
!    ! full slab:
!    allocate( slab(1:nx,1:ny,1:nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( alfa(1:nx,1:ny,1:nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!  
!    ! fill full array:
!    call LE_Bound_Map_BC_to_Model( HtL%Lon_Mapper, HtL%Lat_Mapper, &
!                           HtL%Lev_Mapper, bc_conc, slab, alfa, status )
!    IF_NOTOK_RETURN(status=1)
!    ! add converted contribution:
!    call LE_Bound_Add_3D( ispec, c, dens, alfa, slab*conv, spec_filled, status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! clear:
!    deallocate( slab, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    deallocate( alfa, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! *
!    
!    ! clear:
!    deallocate( conv, stat=status )
!    IF_NOTOK_RETURN(status=1)
!
!    ! ok
!    status = 0
!    
!  end subroutine Heights_to_LE_Apply_Initial
!  
!  
!  ! *
!  
!
!  subroutine LE_Bound_Map_BC_to_Model( Lon_Mapper, Lat_Mapper, Lev_Mapper, &
!                                 bc_conc, md_conc, md_alfa, status )
!
!    use CF, only : CF_AxisMapper_Apply
!
!    ! --- in/out ------------------------------
!
!    type(T_CF_AxisMapper), intent(in)     ::  Lon_Mapper
!    type(T_CF_AxisMapper), intent(in)     ::  Lat_Mapper
!    type(T_CF_AxisMapper), intent(in)     ::  Lev_Mapper(:,:)  ! (nx,ny)
!    real, intent(in)                      ::  bc_conc(:,:,:)   ! input concentrations
!    real, intent(out)                     ::  md_conc(:,:,:)   ! output concentrations
!    real, intent(out)                     ::  md_alfa(:,:,:)   ! fraction covered
!    integer, intent(out)                  ::  status
!
!    ! --- const -------------------------------
!    
!    character(len=*), parameter ::  rname = mname//'/LE_Bound_Map_BC_to_Model'
!    
!    ! --- local -------------------------------
!    
!    integer                   ::  bc_nx, bc_ny, bc_nz
!    integer                   ::  md_nx, md_ny, md_nz
!    real, allocatable         ::  md_conc_yz(:,:,:)
!    real, allocatable         ::  md_conc_z (:,:,:)
!    real, allocatable         ::  md_alfa_yz (:,:,:)
!    real, allocatable         ::  md_alfa2_yz(:,:,:)
!    real, allocatable         ::  md_alfa_z (:,:,:)
!    real, allocatable         ::  md_alfa2_z(:,:,:)
!    real, allocatable         ::  md_alfa2(:,:,:)
!    integer                   ::  i, j, k
!      
!    ! --- begin -------------------------------
!    
!    ! dimensions:
!    bc_nx = size(bc_conc,1) ; bc_ny = size(bc_conc,2) ; bc_nz = size(bc_conc,3)
!    md_nx = size(md_conc,1) ; md_ny = size(md_conc,2) ; md_nz = size(md_conc,3)
!    
!    ! partly interpolations:
!    allocate( md_conc_yz(md_nx,bc_ny,bc_nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( md_conc_z (md_nx,md_ny,bc_nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    ! partly coverages:
!    allocate( md_alfa_yz (md_nx,bc_ny,bc_nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( md_alfa2_yz(md_nx,bc_ny,bc_nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( md_alfa_z (md_nx,md_ny,bc_nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( md_alfa2_z(md_nx,md_ny,bc_nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    allocate( md_alfa2(md_nx,md_ny,md_nz), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! initial values:
!    md_conc_yz  = 0.0
!    md_conc_z   = 0.0
!    md_alfa_yz  = 0.0
!    md_alfa2_yz = 0.0
!    md_alfa_z   = 0.0
!    md_alfa2_z  = 0.0
!    md_alfa2    = 0.0
!  
!    ! loop over input layers:
!    do k = 1, bc_nz
!      ! loop over input latitudes:  
!      do j = 1, bc_ny
!        ! apply longitude conversion:
!        call CF_AxisMapper_Apply( Lon_Mapper, bc_conc(:,j,k), &
!                                    md_conc_yz(:,j,k), status, &
!                                    alfa=md_alfa_yz(:,j,k) )
!        IF_NOTOK_RETURN(status=1)
!      end do
!    end do
!    
!    ! loop over input layers:
!    do k = 1, bc_nz
!      ! loop over output longitudes:
!      do i = 1, md_nx
!        ! apply latitude conversion:
!        call CF_AxisMapper_Apply( Lat_Mapper, md_conc_yz(i,:,k), &
!                                      md_conc_z(i,:,k), status, &
!                                      alfa=md_alfa_z(i,:,k) )
!        IF_NOTOK_RETURN(status=1)
!        ! apply latitude conversion to fractions:
!        call CF_AxisMapper_Apply( Lat_Mapper, md_alfa_yz(i,:,k), &
!                                      md_alfa_z(i,:,k), status, &
!                                      alfa=md_alfa2_z(i,:,k) )
!        IF_NOTOK_RETURN(status=1)
!      end do
!    end do
!    ! combine fractions:
!    md_alfa_z = md_alfa_z * md_alfa2_z
!    
!    ! loop over output latitudes:  
!    do j = 1, md_ny
!      ! loop over output longitudes:
!      do i = 1, md_nx
!        ! apply level conversion:
!        call CF_AxisMapper_Apply( Lev_Mapper(i,j), md_conc_z(i,j,:), &
!                                       md_conc(i,j,:), status, &
!                                       alfa=md_alfa(i,j,:) )
!        IF_NOTOK_RETURN(status=1)
!        ! apply level conversion to fractions:
!        call CF_AxisMapper_Apply( Lev_Mapper(i,j), md_alfa_z(i,j,:), &
!                                       md_alfa(i,j,:), status, &
!                                       alfa=md_alfa2(i,j,:) )
!        IF_NOTOK_RETURN(status=1)
!      end do
!    end do
!    ! combine fractions:
!    md_alfa = md_alfa * md_alfa2
!    
!    ! clear:
!    deallocate( md_conc_yz )
!    deallocate( md_conc_z  )
!    deallocate( md_alfa_yz  )
!    deallocate( md_alfa2_yz )
!    deallocate( md_alfa_z   )
!    deallocate( md_alfa2_z  )
!    deallocate( md_alfa2    )
!        
!    ! ok
!    status = 0
!
!  end subroutine LE_Bound_Map_BC_to_Model
!  
!  
!  ! *
!  
!  
!  !
!  ! Fill new values into concentration array.
!  ! Which spec is determined by 'ispec' .
!  ! If this is an accumulated specie, distribute the new values over the
!  ! components maintaining the original ratio.
!  !
!  ! Example to fill boundary conditions:
!  !  call LE_Bound_Add_2D( ispec, bc_west, alfa(1,:,:), slab(1,:,:)*conv(1,:,:), status )
!  !
!
!  subroutine LE_Bound_Add_2D( ispec, bc, bdens, alfa, slab, spec_filled, status )
!
!    use Binas  , only : xm_air
!    use Dims   , only : nspec
!    use Indices, only : specname
!    use Indices, only : accum_n, accum_ii, accum_ww, accum_ppb_to_ugm3
!
!    ! --- in/out ------------------------------
!
!    integer, intent(in)                   ::  ispec
!    real, intent(inout)                   ::  bc  (:,:,:)   ! (xy,z,nspec) target concentrations
!    real, intent(in)                      ::  bdens(:,:)    ! (xy,z) air density slab [kg/m3]
!    real, intent(in)                      ::  alfa(:,:)     ! (xy,z) fraction covered
!    real, intent(in)                      ::  slab(:,:)     ! (xy,z) new values, converted!
!    logical, intent(inout)                ::  spec_filled(:)  ! (ispec) .true. if specie becomes filled
!    integer, intent(out)                  ::  status
!
!    ! --- const -------------------------------
!    
!    character(len=*), parameter ::  rname = mname//'/LE_Bound_Add_2D'
!    
!    ! --- local -------------------------------
!    
!    integer               ::  k
!    integer               ::  icomp
!    integer               ::  n1, n2
!    real, allocatable     ::  tot(:,:)
!    real, allocatable     ::  frac(:,:)
!    real, allocatable     ::  convfact(:,:)
!      
!    ! --- begin -------------------------------
!
!    ! support both single and accumulated species ...
!    if ( accum_n(ispec) == 1 ) then
!    
!      ! sinle species ...
!
!      ! add fraction of new values:
!      bc(:,:,ispec) = bc(:,:,ispec) * (1-alfa) + slab * alfa
!      ! set flag:
!      spec_filled(ispec) = .true.
!
!    else
!      
!      ! accumulated ...
!    
!      ! dimensions:
!      n1 = size(slab,1)
!      n2 = size(slab,2)
!
!      ! support both single and accumulated species ...
!
!      ! storage:
!      allocate( tot (n1,n2) )
!      allocate( frac(n1,n2) )
!      allocate( convfact(n1,n2) )
!
!      ! ~ compute total:
!      ! init to zero:
!      tot = 0.0
!      ! loop over components:
!      do k = 1, accum_n(ispec)
!        ! component index:
!        icomp = accum_ii(ispec,k)
!        ! conversion needed ?
!        convfact = 1.0
!        if ( accum_ppb_to_ugm3(ispec,k) ) convfact = bdens/xm_air
!        ! add weighted contribution to total:
!        tot = tot + bc(:,:,icomp) * accum_ww(ispec,k) * convfact
!      end do
!
!      ! ~ distribute new values according to original ratio:
!      ! loop over components:
!      do k = 1, accum_n(ispec)
!        ! component index:
!        icomp = accum_ii(ispec,k)
!        ! conversion needed ?
!        convfact = 1.0
!        if ( accum_ppb_to_ugm3(ispec,k) ) convfact = bdens/xm_air
!        ! could only fractionize if there was some tracer at all:
!        where ( tot > 0.0 )
!          ! fraction of this component in original total:
!          frac = bc(:,:,icomp) * accum_ww(ispec,k) * convfact / tot
!          ! add fraction of new values:
!          bc(:,:,icomp) = bc(:,:,icomp) * (1-alfa) + slab*frac * alfa
!        endwhere
!        ! set flag:
!        spec_filled(icomp) = .true.
!      end do
!
!      ! clear:
!      deallocate( tot  )
!      deallocate( frac )
!      deallocate( convfact )
!      
!    end if
!    
!    ! ok
!    status = 0
!
!  end subroutine LE_Bound_Add_2D
!  
!  ! *
!
!  subroutine LE_Bound_Add_3D( ispec, c, dens, alfa, slab, spec_filled, status )
!
!    use Binas  , only : xm_air
!    use Dims   , only : nspec
!    use Indices, only : specname
!    use Indices, only : accum_n, accum_ii, accum_ww, accum_ppb_to_ugm3
!
!    ! --- in/out ------------------------------
!
!    integer, intent(in)                   ::  ispec
!    real, intent(inout)                   ::  c   (:,:,:,:)   ! (x,y,z,nspec) target concentrations
!    real, intent(inout)                   ::  dens(:,:,:)     ! (x,y,z) air density (kg/m3)
!    real, intent(in)                      ::  alfa(:,:,:)     ! (x,y,z) fraction covered
!    real, intent(in)                      ::  slab(:,:,:)     ! (x,y,z) new values, converted!
!    logical, intent(inout)                ::  spec_filled(:)  ! (ispec) .true. if specie becomes filled
!    integer, intent(out)                  ::  status
!
!    ! --- const -------------------------------
!    
!    character(len=*), parameter ::  rname = mname//'/LE_Bound_Add_3D'
!    
!    ! --- local -------------------------------
!    
!    integer               ::  k
!    integer               ::  icomp
!    integer               ::  n1, n2, n3
!    real, allocatable     ::  tot(:,:,:)
!    real, allocatable     ::  frac(:,:,:)
!    real, allocatable     ::  convfact(:,:,:)
!      
!    ! --- begin -------------------------------
!    
!    ! support both single and accumulated species ...
!    if ( accum_n(ispec) == 1 ) then
!    
!      ! sinle species ...
!
!      ! add fraction of new values:
!      c(:,:,:,ispec) = c(:,:,:,ispec) * (1-alfa) + slab * alfa
!      ! set flag:
!      spec_filled(ispec) = .true.
!
!    else
!      
!      ! accumulated ...
!    
!      ! dimensions:
!      n1 = size(slab,1)
!      n2 = size(slab,2)
!      n3 = size(slab,3)
!
!      ! storage:
!      allocate( tot (n1,n2,n3) )
!      allocate( frac(n1,n2,n3) )
!      allocate( convfact(n1,n2,n3) )
!
!      ! ~ compute total:
!      ! init to zero:
!      tot = 0.0
!      ! loop over components:
!      do k = 1, accum_n(ispec)
!        ! component index:
!        icomp = accum_ii(ispec,k)
!        ! conversion needed ?
!        convfact = 1.0
!        if ( accum_ppb_to_ugm3(ispec,k) ) convfact = dens/xm_air
!        ! add weighted contribution to total:
!        tot = tot + c(:,:,:,icomp) * accum_ww(ispec,k) * convfact
!      end do
!
!      ! ~ distribute new values according to original ratio:
!      ! loop over components:
!      do k = 1, accum_n(ispec)
!        ! component index:
!        icomp = accum_ii(ispec,k)
!        ! conversion needed ?
!        convfact = 1.0
!        if ( accum_ppb_to_ugm3(ispec,k) ) convfact = dens/xm_air
!        ! could only fractionize if there was some tracer at all:
!        where ( tot > 0.0 )
!          ! fraction of this component in original total:
!          frac = c(:,:,:,icomp) * accum_ww(ispec,k) * convfact / tot
!          ! add fraction of new values:
!          c(:,:,:,icomp) = c(:,:,:,icomp) * (1-alfa) + slab*frac * alfa
!        endwhere
!        ! set flag:
!        spec_filled(icomp) = .true.
!      end do
!
!      ! clear:
!      deallocate( tot  )
!      deallocate( frac )
!      deallocate( convfact )
!      
!    end if   ! single or accum
!    
!    ! ok
!    status = 0
!
!  end subroutine LE_Bound_Add_3D
  

end module LE_Bound_Tools

