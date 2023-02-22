!#################################################################
!
! NAME
!
!   LE_MAORI  -  LOTOS-EUROS access to MAORI routines
!
! DESCRIPTION
!
!   Interface to LE output of various types.
!
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#################################################################

module LE_MAORI

  use GO, only : gol, goPr, goErr

  use MAORI, only : T_MAORI_State, T_MAORI_Data, T_MAORI_Output

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_MAORI_State, T_MAORI_Data, T_MAORI_Output
  public  ::  LE_MAORI_Data_Init, LE_MAORI_Data_Done, LE_MAORI_Data_Setup
  public  ::  LE_MAORI_State_Init, LE_MAORI_State_Done, LE_MAORI_State_Setup, LE_MAORI_State_Put
  public  ::  LE_MAORI_Output_Init, LE_MAORI_Output_Done, LE_MAORI_Output_Write


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_MAORI'




contains


  ! ====================================================================
  ! ===
  ! === maori data
  ! ===
  ! ====================================================================


  subroutine LE_MAORI_Data_Init( mdata, rcfile, t0, status )

    use GO, only : TDate

    use Dims           , only : nz
    use LE_Grid        , only : ugg
    use LE_BiasCorr    , only : biascorr__aod__factor

    use MAORI, only : T_MAORI_Data
    use MAORI, only : MAORI_Data_Init, MAORI_Data_Inquire, MAORI_Data_Start
    use MAORI, only : MAORI_Data_Set_Inquire
    use MAORI, only : MAORI_Data_Param_Inquire, MAORI_Data_Param_Put
    use MAORI, only : MAORI_TYPE_NAME, MAORI_INT, MAORI_REAL, MAORI_CHAR
    use MAORI, only : MAORI_LEN_NAME, MAORI_LEN_UNIT, MAORI_LEN_LINE

    ! --- in/out --------------------------------

    type(T_MAORI_Data), intent(out)     ::  mdata
    character(len=*), intent(in)        ::  rcfile
    type(TDate), intent(in)             ::  t0
    integer, intent(out)                ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_Data_Init'

    ! --- local ---------------------------------

    character(len=MAORI_LEN_NAME)   ::  sname
    integer                         ::  nset, iset
    integer                         ::  npar, ipar
    integer                         ::  ptype
    character(len=MAORI_LEN_NAME)   ::  pname
    character(len=MAORI_LEN_UNIT)   ::  punit
    character(len=MAORI_LEN_LINE)   ::  name_unit
    integer                         ::  ival
    real                            ::  rval
    character(len=MAORI_LEN_LINE)   ::  sval

    ! --- begin ---------------------------------

    ! initialise main data:
    call MAORI_Data_Init( mdata, trim(rcfile), t0, status )
    IF_NOTOK_RETURN(status=1)

    ! get number of sets to put out:
    call MAORI_Data_Inquire( mdata, status, nset=nset )
    IF_NOTOK_RETURN(status=1)

    ! loop over sets:
    do iset = 1, nset

      ! get number of requested parameters:
      call MAORI_Data_Set_Inquire( mdata, iset, status, name=sname, nparam=npar )
      IF_NOTOK_RETURN(status=1)

      ! info:
      write (gol,'("requested parameters for set ",i6," ",a," :")') iset, trim(sname); call goPr

      ! loop over requested parameters:
      do ipar = 1, npar

        ! get parameter type, name and required unit:
        call MAORI_Data_Param_Inquire( mdata, iset, ipar, &
                                        ptype, pname, punit, status )
       IF_NOTOK_RETURN(status=1)

        ! info:
        write (gol,'("  ",i6," `",a,"` (",a,")")') ipar, trim(pname), &
                                 trim(MAORI_TYPE_NAME(ptype)); call goPr

        ! fill parameter 'name [unit]'
        write (name_unit,'(a," [",a,"]")') trim(pname), trim(punit)
        ! extract per type:
        select case ( ptype )
          ! integer parameters
          case ( MAORI_INT )
            ! extract value from model data in requested unit:
            select case ( trim(name_unit) )
              case ( 'grid_nlon [1]' ) ; ival = ugg%nlon  ! grid size
              case ( 'grid_nlat [1]' ) ; ival = ugg%nlat  ! grid size
              case ( 'nlayer [1]'    ) ; ival = nz+1      ! all layers plus aloft
              case default
                write (gol,'("unsupported real parameter name [unit] : ",a)') trim(name_unit); call goErr
                TRACEBACK; status=1; return
            end select
            ! set parameter:
            call MAORI_Data_Param_Put( mdata, iset, ipar, ival, status )
            IF_NOTOK_RETURN(status=1)
          ! real parameters
          case ( MAORI_REAL )
            ! extract value from model data in requested unit:
            select case ( trim(name_unit) )
              case ( 'bound_west [degree_east]'   ) ; rval = ugg%longitude_bnds_1d(1,1)
              case ( 'bound_east [degree_east]'   ) ; rval = ugg%longitude_bnds_1d(2,ugg%nlon)
              case ( 'bound_south [degree_north]' ) ; rval = ugg%latitude_bnds_1d(1,1)
              case ( 'bound_north [degree_north]' ) ; rval = ugg%latitude_bnds_1d(2,ugg%nlat)
              case default
                write (gol,'("unsupported real parameter name [unit] : ",a)') trim(name_unit); call goErr
                TRACEBACK; status=1; return
            end select
            ! set parameter:
            call MAORI_Data_Param_Put( mdata, iset, ipar, rval, status )
            IF_NOTOK_RETURN(status=1)
          ! character parameters
          case ( MAORI_CHAR )
            ! extract value from model data in requested unit:
            select case ( trim(name_unit) )
              case ( 'bias_correction [free]'   )
                !write (sval,'("surface_ozone=",a,";total_pm_factor=",f4.1,";total_aod_factor=",f4.1)') &
                !         trim(biascorr__surface_ozone), biascorr__total_pm__factor, biascorr__aod__factor
                write (sval,'("total_aod_factor=",f4.1)') biascorr__aod__factor
              case default
                write (gol,'("unsupported char parameter name [unit] : ",a)') trim(name_unit); call goErr
                TRACEBACK; status=1; return
            end select
            ! set parameter:
            call MAORI_Data_Param_Put( mdata, iset, ipar, sval, status )
            IF_NOTOK_RETURN(status=1)
          case default
            write (*,'("ERROR - unsupported maori type ",i6," for parameter ",a)') &
                         ptype, trim(pname)
            TRACEBACK; status=1; return
        end select

      end do  ! parameters

    end do  ! sets

    ! start outut:
    call MAORI_Data_Start( mdata, status )
    IF_NOTOK_RETURN(status=1)

    !
    ! done
    !

    ! ok
    status = 0

  end subroutine LE_MAORI_Data_Init


  ! ***


  subroutine LE_MAORI_Data_Done( mdata, status )

    use MAORI, only : T_MAORI_Data
    use MAORI, only : MAORI_Data_Done

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)   ::  mdata
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_Data_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! done with maori stuff:
    call MAORI_Data_Done( mdata, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_MAORI_Data_Done


  ! ***


  subroutine LE_MAORI_Data_Setup( mdata, t1, t2, the_end, status )

    use GO   , only : TDate
    use MAORI, only : T_MAORI_Data
    use MAORI, only : MAORI_Data_Setup

    ! --- in/out --------------------------------

    type(T_MAORI_Data), intent(inout)   ::  mdata
    type(TDate), intent(in)             ::  t1, t2
    logical                             ::  the_end
    integer, intent(out)                ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_Data_Setup'

    ! --- local ---------------------------------

    ! --- begin ---------------------------------

    ! setup MAORI output for current interval:
    call MAORI_Data_Setup( mdata, t1, t2, the_end, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_MAORI_Data_Setup


  ! ====================================================================
  ! ===
  ! === state
  ! ===
  ! ====================================================================


  subroutine LE_MAORI_State_Init( mstate, mdata, name, status )

    use MAORI, only : T_MAORI_State, T_MAORI_Data
    use MAORI, only : MAORI_State_Init, MAORI_State_Start

    ! --- in/out --------------------------------

    type(T_MAORI_State), intent(out)      ::  mstate
    type(T_MAORI_Data), intent(in)        ::  mdata
    character(len=*), intent(in)          ::  name
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_State_Init'

    ! --- local ---------------------------------

    ! --- begin ---------------------------------

    ! init maori state:
    call MAORI_State_Init( mstate, mdata, name, status )
    IF_NOTOK_RETURN(status=1)

    ! could be started imediatately, sinde data has been started too:
    call MAORI_State_Start( mstate, mdata, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_MAORI_State_Init


  ! ***


  subroutine LE_MAORI_State_Done( mstate, mdata, status )

    use MAORI, only : T_MAORI_Data, T_MAORI_State
    use MAORI, only : MAORI_State_Done

    ! --- in/out --------------------------------

    type(T_MAORI_State), intent(inout)        ::  mstate
    type(T_MAORI_Data), intent(in)            ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_State_Done'

    ! --- local ------------------------------

    ! --- begin ---------------------------------

    ! done with maori stuff:
    call MAORI_State_Done( mstate, mdata, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_MAORI_State_Done


  ! ***


  subroutine LE_MAORI_State_Setup( mstate, mdata, status )

    use MAORI, only : T_MAORI_Data, T_MAORI_State
    use MAORI, only : MAORI_State_Setup

    ! --- in/out --------------------------------

    type(T_MAORI_State), intent(inout)        ::  mstate
    type(T_MAORI_Data), intent(in)            ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_State_Setup'

    ! --- local ------------------------------

    ! --- begin ---------------------------------

    ! done with maori stuff:
    call MAORI_State_Setup( mstate, mdata, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_MAORI_State_Setup


  ! ***


  subroutine LE_MAORI_State_Put( maori_state, maori_data, t, c, cg, aerh2o, status )

    use Binas            , only : xm_air
    use GO               , only : TDate
    use GO               , only : goMatchValue
    use Grid             , only : GetLocation
    use dims             , only : nx, ny, nz
    use Indices          , only : nspec, specname, specunit, specmolm
    use Indices          , only : accum_n, accum_ii, accum_ww
    use LE_Bound         , only : caloft
    use LE_Data          , only : LE_Data_GetPointer
    use LE_Grid          , only : ugg
    use LE_BiasCorr      , only : biascorr__aod__factor
    !use LE_AOD           , only : LE_AOD_Calc3D

    use MAORI, only : T_MAORI_Data, T_MAORI_State
    use MAORI, only : MAORI_TYPE_NAME, MAORI_SAMPLE, MAORI_SATELLITE
    use MAORI, only : MAORI_LEN_NAME, MAORI_LEN_UNIT
    use MAORI, only : MAORI_Data_Inquire
    use MAORI, only : MAORI_Data_Set_Inquire
    use MAORI, only : MAORI_Data_Loc_Inquire
    use MAORI, only : MAORI_Data_Var_Inquire
    use MAORI, only : MAORI_State_Setup
    use MAORI, only : MAORI_State_Put, MAORI_State_Post

    ! --- in/out --------------------------------

    type(T_MAORI_State), intent(inout)        ::  maori_state
    type(T_MAORI_Data), intent(in)            ::  maori_data
    type(TDate), intent(in)                   ::  t
    real, intent(in)                          ::  c(nx,ny,nz,nspec) ! tracer conc (ppb or ug/m3)
    real, intent(in)                          ::  cg(nx,ny,nspec)   ! ground conc (ppb or ug/m3)
    real, intent(in)                          ::  aerh2o(nx,ny,nz)  ! aerosol water conc (ug/m3)
    integer, intent(out)                      ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_State_Put'

    ! --- local ---------------------------------

    real                            ::  AOD(nx,ny,nz)
    real                            ::  extinction(nx,ny,nz)

    real, pointer                   ::  dens(:,:,:)
    real, pointer                   ::  Rh(:,:,:)
    real, pointer                   ::  temp(:,:,:)
    real, pointer                   ::  h(:,:,:)
    real, pointer                   ::  oro(:,:,:)

    integer                         ::  nset, iset
    integer                         ::  stype
    integer                         ::  nloc, iloc
    integer                         ::  nlev
    integer                         ::  nvar, ivar
    real                            ::  lon, lat
    integer                         ::  nlon, nlat
    integer                         ::  ilon, ilat
    integer                         ::  itr
    real                            ::  acw
    character(len=MAORI_LEN_NAME)   ::  vname
    character(len=MAORI_LEN_NAME)   ::  scode
    character(len=MAORI_LEN_UNIT)   ::  vunit
    real, allocatable               ::  profile(:)
    real, allocatable               ::  unitfacs(:)
    integer                         ::  icomp, ispec

    ! --- begin ---------------------------------
    
    ! meteo:
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'rh', Rh, status, check_units ='%' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 't', temp, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'oro', oro, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)

    ! AOD:
    !call LE_AOD_Calc3D( c, AOD, extinction, status )
    !IF_NOTOK_RETURN(status=1)
    AOD=0
    ! temporary storage;
    ! levels + aloft, and one extra for half levels:
    allocate( profile(1:nz+2) )
    allocate( unitfacs(1:nz+2) )

    !>>> this should be done only once, now  'LE_MAORI_State_Setup'
    !! setup MAORI output for current interval:
    !call MAORI_State_Setup( maori_state, maori_data, status )
    !IF_NOTOK_RETURN(status=1)
    !<<<

    ! get number of sets to put out:
    call MAORI_Data_Inquire( maori_data, status, nset=nset )
    IF_NOTOK_RETURN(status=1)

    ! loop over sets:
    do iset = 1, nset

      ! get type of output:
      call MAORI_Data_Set_Inquire( maori_data, iset, status, type=stype )
      IF_NOTOK_RETURN(status=1)

      ! do something given the type:
      select case ( stype )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( MAORI_SAMPLE )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! get number of locations and number of variables to be written:
          call MAORI_Data_Set_Inquire( maori_data, iset, status, nloc=nloc, nvar=nvar )
          IF_NOTOK_RETURN(status=1)

          ! loop over stations:
          do iloc = 1, nloc

            ! get (lon,lat) of location:
            call MAORI_Data_Loc_Inquire( maori_data, iset, iloc, status, &
                                            lon=lon, lat=lat, code=scode )
            IF_NOTOK_RETURN(status=1)

            ! get cell indices:
            call ugg%GetLocation( lon, lat, ilon, ilat, status )
            if (status/=0) then
              write (gol,'("sample location not in grid : ",2f8.2)') lon, lat; call goErr
              TRACEBACK; status=1; return
            end if

            ! set default:
            profile = 0.0
            unitfacs = 1.0

            ! loop over variables:
            do ivar = 1, nvar

              ! get variable name and unit:
              call MAORI_Data_Var_Inquire( maori_data, iset, ivar, status, &
                                         name=vname, unit=vunit, nlev=nlev )
              IF_NOTOK_RETURN(status=1)

              ! check ...
              if ( (nlev < 1) .or. (nlev > nz+2) ) then
                write (gol,'("unsupported nlev ",i4,"; should be in 1..nz+2 with nz=",i4)') nlev, nz; call goErr
                TRACEBACK; status=1; return
              end if

              ! search variable name in tracer names:
              call goMatchValue( trim(vname), specname, ispec, status, quiet=.true. )
              IF_ERROR_RETURN(status=1)
              ! found matching specname ?
              if ( status == 0 ) then
                ! might be accumulated tracer; init sum:
                profile = 0.0
                ! loop over contributions:
                do icomp = 1, accum_n(ispec)
                  ! tracer index:
                  itr = accum_ii(ispec,icomp)
                  ! weight of tracer:
                  acw = accum_ww(ispec,icomp)
                  ! add contributions:
                  if ( nlev == 1 ) then
                    profile(1)    = profile(1) + cg(ilon,ilat,itr) * acw
                  else
                    profile(1:nz) = profile(1:nz) +      c(ilon,ilat,1:nz,itr) * acw
                    profile(nz+1) = profile(nz+1) + caloft(ilon,ilat,nz+1,itr) * acw
                  end if
                end do
                ! convert ?
                if ( trim(specunit(ispec)) /= trim(vunit) ) then
                  ! unit conversion factor:
                  select case ( trim(specunit(ispec))//' -> '//trim(vunit) )
                    ! scale:
                    case ( 'ppb -> mole mole-1' ) ; unitfacs = 1.0e-9   ! ppb -> mole/mole
                    case ( 'ug/m3 -> kg m-3'    ) ; unitfacs = 1.0e-9   ! ug/m3 -> kg/m3
                    ! volume mixing ratio to mass concentration:
                    case ( 'ppb -> kg m-3' )
                      ! (kg/m3)/ppb  = ((mole tr/mole air)/ppb) * (kg tr/mole tr) *   (kg air/m3 air)    / (kg air/mole air)
                      unitfacs(1:nz) =       1e-9               * specmolm(ispec) * dens(ilon,ilat,1:nz) /     xm_air
                      unitfacs(nz+1) = unitfacs(nz)
                    ! unkown ...
                    case default
                      write (gol,'("do not know how to convert from `",a,"` to `",a,"`")') &
                              trim(specunit(ispec)), trim(vunit); call goErr
                      TRACEBACK; status=1; return
                  end select
                  ! convert:
                  profile = profile * unitfacs
                end if  ! different units
              else
                ! not a tracer value ...
                ! extract value from model data:
                select case ( trim(vname)//' ['//trim(vunit)//']' )
                  case ( 'aod [1]', 'aod_biascorr [1]' )
                    ! total column or profile ?
                    if ( nlev == 1 ) then
                      profile(1) = sum(AOD(ilon,ilat,1:nz))  ! 1
                    else
                      profile(1:nz) = AOD(ilon,ilat,1:nz)  ! 1
                      profile(nz+1) = 0.0                  ! 1
                    end if
                    ! bias corrected ?
                    if ( trim(vname) == 'aod_biascorr' ) profile = profile * biascorr__aod__factor
                  case ( 'aerh2o [kg m-3]' )
                    profile(1:nz) = aerh2o(ilon,ilat,1:nz) * 1e-9  ! ug/m3 -> kg/m3
                    profile(nz+1) = 0.0
                  case ( 'rhumid [%]' )
                    profile(1:nz) = Rh(ilon,ilat,1:nz)   ! %
                    profile(nz+1) = 0.0                  ! %
                  case ( 'temper [K]' )
                    profile(1:nz) = temp(ilon,ilat,1:nz)  ! K
                    profile(nz+1) = temp(ilon,ilat,  nz)  ! K
                  case ( 'halflevel_altitude [m]' )
                    profile(1     ) = oro(ilon,ilat,1)     ! m
                    profile(2:nz+1) = h(ilon,ilat,1:nz)    ! m
                    profile(  nz+2) = 200.0e3   ! 200 km -> m                  ! kg/kg
                  case ( 'dens [kg m-3]' )
                    profile(1:nz) = dens(ilon,ilat,1:nz)  ! kg/m3
                    profile(nz+1) = dens(ilon,ilat,  nz)  ! kg/m3
                  case default
                    write (gol,'("unsupported variable name or unit : ",a," [",a,"]")') trim(vname), trim(vunit); call goErr
                    TRACEBACK; status=1; return
                end select
              end if

              ! put profile:
              call MAORI_State_Put( maori_state, maori_data, iset, ivar, status, &
                                      iloc=iloc, values=profile(1:nlev) )
              IF_NOTOK_RETURN(status=1)

            end do   ! variables

          end do  ! stations

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( MAORI_SATELLITE )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! get number of locations and number of variables to be written:
          call MAORI_Data_Set_Inquire( maori_data, iset, status, &
                                        nlon=nlon, nlat=nlat, nloc=nloc, nvar=nvar )
          IF_NOTOK_RETURN(status=1)

          ! loop over points:
          do iloc = 1, nloc

            ! get location:
            call MAORI_Data_Loc_Inquire( maori_data, iset, iloc, status, &
                                                 lon=lon, lat=lat )
            IF_NOTOK_RETURN(status=1)
            ! get cell location:
            call ugg%GetLocation( lon, lat, ilon, ilat, status )
            if (status/=0) then
              write (gol,'("sample location not in grid : ",2f8.2)') lon, lat; call goErr
              TRACEBACK; status=1; return
            end if

            ! set default:
            profile = 0.0

            ! loop over variables:
            do ivar = 1, nvar

              ! get variable name and unit:
              call MAORI_Data_Var_Inquire( maori_data, iset, ivar, status, &
                                         name=vname, unit=vunit, nlev=nlev )
              IF_NOTOK_RETURN(status=1)

              ! check ...
              if ( nlev /= 1 ) then
                write (gol,'("unsupported nlev ",i4,"; should be 1 ...")') nlev; call goErr
                TRACEBACK; status=1; return
              end if

              ! extract value from model data:
              select case ( trim(vname)//' ['//trim(vunit)//']' )
                case ( 'aod [1]' )
                  ! total column:
                  profile(1) = sum(AOD(ilon,ilat,1:nz))  ! m
                case ( 'aod_biascorr [1]' )
                  ! total column :
		          profile(1) = sum(AOD(ilon,ilat,1:nz))   * biascorr__aod__factor
                case default
                  write (gol,'("unsupported variable name or unit : ",a," [",a,"]")') trim(vname), trim(vunit); call goErr
                  TRACEBACK; status=1; return
              end select

              ! put profile:
              call MAORI_State_Put( maori_state, maori_data, iset, ivar, status, &
                                      iloc=iloc, values=profile(1:nlev) )
              IF_NOTOK_RETURN(status=1)

            end do   ! variables

          end do  ! iloc

          ! loop over grid cells:
          do ilat = 1, nlat
            do ilon = 1, nlon

              ! set default:
              profile = 0.0

              ! loop over variables:
              do ivar = 1, nvar

                ! get variable name and unit:
                call MAORI_Data_Var_Inquire( maori_data, iset, ivar, status, &
                                           name=vname, unit=vunit, nlev=nlev )
                IF_NOTOK_RETURN(status=1)

                ! check ...
                if ( nlev /= 1 ) then
                  write (gol,'("unsupported nlev ",i4,"; should be 1 ...")') nlev; call goErr
                  TRACEBACK; status=1; return
                end if

                ! extract value from model data:
                select case ( trim(vname)//' ['//trim(vunit)//']' )
                  case ( 'aod [1]' )
                    ! total column:
                    profile(1) = sum(AOD(ilon,ilat,1:nz))  ! 1
                  case ( 'aod_biascorr [1]' )
                    ! total column :
		            profile(1) = sum(AOD(ilon,ilat,1:nz))   * biascorr__aod__factor  ! 1
                  case default
                    write (gol,'("unsupported variable name or unit : ",a," [",a,"]")') trim(vname), trim(vunit); call goErr
                    TRACEBACK; status=1; return
                end select

                ! put profile:
                call MAORI_State_Put( maori_state, maori_data, iset, ivar, status, &
                                        ilon=ilon, ilat=ilat, values=profile(1:nlev) )
                IF_NOTOK_RETURN(status=1)

              end do   ! variables

            end do   ! ipx
          end do  ! ipy

        ! ~ unknown ...
        case default

          write (gol,'("ERROR - unknown maori type ",i6," (`",a,"`) for set ",i6)') &
                                   stype, trim(MAORI_TYPE_NAME(stype)), iset; call goPr
          TRACEBACK; status=1; return

      end select

    end do   ! sets

    ! update simulations, eventually write to output files:
    call MAORI_State_Post( maori_state, maori_data, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( profile )

    ! ok
    status = 0

  end subroutine LE_MAORI_State_Put




  ! ====================================================================
  ! ===
  ! === output
  ! ===
  ! ====================================================================


  subroutine LE_MAORI_Output_Init( moutput, mdata, name, status )

    use MAORI, only : T_MAORI_Output, T_MAORI_Data
    use MAORI, only : MAORI_Output_Init, MAORI_Output_Start

    ! --- in/out --------------------------------

    type(T_MAORI_Output), intent(out)      ::  moutput
    type(T_MAORI_Data), intent(in)        ::  mdata
    character(len=*), intent(in)          ::  name
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_Output_Init'

    ! --- local ---------------------------------

    ! --- begin ---------------------------------

    ! init maori state:
    call MAORI_Output_Init( moutput, mdata, name, status )
    IF_NOTOK_RETURN(status=1)

    ! could be started imediatately, sinde data has been started too:
    call MAORI_Output_Start( moutput, mdata, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_MAORI_Output_Init


  ! ***


  subroutine LE_MAORI_Output_Done( moutput, mdata, status )

    use MAORI, only : T_MAORI_Data, T_MAORI_Output
    use MAORI, only : MAORI_Output_Done

    ! --- in/out --------------------------------

    type(T_MAORI_Output), intent(inout)       ::  moutput
    type(T_MAORI_Data), intent(in)            ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_Output_Done'

    ! --- local ------------------------------

    ! --- begin ---------------------------------

    ! done with maori stuff:
    call MAORI_Output_Done( moutput, mdata, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_MAORI_Output_Done


  ! ***


  subroutine LE_MAORI_Output_Write( moutput, mstate, mdata, status )

    use MAORI, only : T_MAORI_Data, T_MAORI_Output
    use MAORI, only : MAORI_Output_Write

    ! --- in/out ---------------------------------

    type(T_MAORI_Output), intent(inout)       ::  moutput
    type(T_MAORI_State), intent(inout)        ::  mstate
    type(T_MAORI_Data), intent(in)            ::  mdata
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_MAORI_Output_Write'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! done with maori stuff:
    call MAORI_Output_Write( moutput, mstate, mdata, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_MAORI_Output_Write


end module LE_MAORI
