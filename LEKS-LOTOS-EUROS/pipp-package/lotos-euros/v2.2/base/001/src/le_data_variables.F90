!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Data_Variables

  use GO  , only : gol, goPr, goErr, goLabel
  use LE_Data_Variable, only : T_Variable

  ! --- in/out -----------------------------------
  
  private

  public  ::  T_Variables


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Data_Variables'
  
  ! maximum number of dependencies:
  integer, parameter          ::  maxdep = 10


  ! --- types ------------------------------------

  type T_Variables
    ! number of variables:
    integer                                   ::  nvar
    ! list:
    type(T_Variable), pointer                 ::  var(:)
  contains
    procedure :: Init            => Variables_Init
    procedure :: Done            => Variables_Done
    procedure :: GetID           => Variables_GetID
    procedure :: GetVariable     => Variables_GetVariable
    procedure :: SetVariable     => Variables_SetVariable
    procedure :: Setup           => Variables_Setup
    !
    procedure ::                    Variables_GetValidPtr
    procedure ::                    Variables_GetValidPtr_uv
    generic   :: GetValidPtr     => Variables_GetValidPtr, &
                                    Variables_GetValidPtr_uv
  end type T_Variables



contains


  ! ====================================================================
  ! ===
  ! === Variable
  ! ===
  ! ====================================================================


  subroutine Variables_Init( self, rcF, status )
  
    use GO              , only : TrcFile, ReadRc
    use GO              , only : goSplitString
    use GO              , only : T_Field
    use LE_Data_Variable, only : LEN_VARNAME
  
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(out)         ::  self
    type(TrcFile), pointer                  ::  rcF
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Variables_Init'
    
    integer, parameter    ::  maxvar = 500
    
    ! --- local ------------------------------------------
    
    character(len=4000)           ::  line
    character(len=LEN_VARNAME)    ::  varnames(maxvar)
    integer                       ::  ivar
    class(T_Field), pointer       ::  F
    
    ! --- begin ------------------------------------------
    
    ! read line with supported variables:
    call ReadRc( rcF, 'data.vars', line, status )
    IF_NOTOK_RETURN(status=1)
    ! split:
    call goSplitString( line, self%nvar, varnames, status )
    IF_NOTOK_RETURN(status=1)

    ! storage:
    allocate( self%var(self%nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! loop:
    do ivar = 1, self%nvar
      ! init variable:
      call self%var(ivar)%Init( rcF, trim(varnames(ivar)), status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! ok
    status = 0
    
  end subroutine Variables_Init


  ! ***


  subroutine Variables_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Done'
    
    ! --- local ------------------------------------------
    
    integer                       ::  ivar
      
    ! --- begin ----------------------------------

    ! loop:
    do ivar = 1, self%nvar
      ! init variable:
      call self%var(ivar)%Done( status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! clear:
    deallocate( self%var, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Variables_Done

  ! ***


  subroutine Variables_GetID( self, name, id, status, silent )
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(in)                ::  self
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  id
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  silent
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_GetID'
    
    ! --- local ------------------------------------------
    
    integer                       ::  ivar
    logical                       ::  shout
      
    ! --- begin ----------------------------------
    
    ! shout?
    shout = .true.
    if ( present(silent) ) shout = .not. silent
    
    ! init result:
    id = -1
    ! loop:
    do ivar = 1, self%nvar
      ! compare:
      if ( trim(name) == trim(self%var(ivar)%name) ) then
        ! store:
        id = ivar
        ! leave:
        exit
      end if  ! match
    end do  ! variables
    ! not found?
    if ( id < 0 ) then
      ! message?
      if ( shout ) then
        write (gol,'("no variable with name `",a,"`")') trim(name); call goErr
        TRACEBACK; status=1; return
      end if
      ! warning status:
      status = -1; return
    end if  ! found
    
    ! ok
    status = 0
    
  end subroutine Variables_GetID


  ! ***


  subroutine Variables_GetVariable( self, name, status, V, &
                                      units, check_units, &
                                      Field, VectorField, Constant_Field, Accumulated_Field_Series, &
                                        Instant_Field, Instant_Field_Series, &
                                        Instant_VectorField, &
                                      pdata, levtype, gridtype, lbo, ubo, check_lbo, check_ubo, enabled )

    use GO            , only : T_Field
    use GO            , only : T_VectorField
    use GO            , only : T_Constant_Field
    use GO            , only : T_Accumulated_Field_Series
    use GO            , only : T_Instant_Field
    use GO            , only : T_Instant_Field_Series
    use GO            , only : T_Instant_VectorField
    use LE_Data_Common, only : LE_Data_CompareUnits
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(in)                ::  self
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status
    
    character(len=*), intent(out), optional               ::  units
    character(len=*), intent(in), optional                ::  check_units
    class(T_Variable), pointer, intent(out), optional     ::  V
    class(T_Field), pointer, optional                     ::  Field
    class(T_VectorField), pointer, optional               ::  VectorField
    class(T_Constant_Field), pointer, optional            ::  Constant_Field
    class(T_Accumulated_Field_Series), pointer, optional  ::  Accumulated_Field_Series
    class(T_Instant_Field), pointer, optional             ::  Instant_Field
    class(T_Instant_Field_Series), pointer, optional      ::  Instant_Field_Series
    class(T_Instant_VectorField), pointer, optional       ::  Instant_VectorField
    real, pointer, optional                               ::  pdata(:,:,:)
    character(len=*), intent(out), optional               ::  levtype
    logical, intent(out), optional                        ::  enabled
    character(len=*), intent(out), optional               ::  gridtype
    integer, intent(out), optional                        ::  lbo(3), ubo(3)
    integer, intent(in), optional                         ::  check_lbo(3), check_ubo(3)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_GetVariable'
    
    ! --- local ----------------------------------
    
    integer                       ::  id
      
    ! --- begin ----------------------------------
    
    ! get id by name:
    call self%GetID( trim(name), id, status )
    IF_NOTOK_RETURN(status=1)
    
    ! attributes:
    call self%var(id)%Get( status, enabled=enabled, units=units, &
                              levtype=levtype, gridtype=gridtype, &
                              lbo=lbo, ubo=ubo, &
                              check_lbo=check_lbo, check_ubo=check_ubo )
    IF_NOTOK_RETURN(status=1)

    ! check units?
    if ( present(check_units) ) then
      ! compare, status<0 if not equivalent:
      call LE_Data_CompareUnits( trim(self%var(id)%units), trim(check_units), status )
      IF_ERROR_RETURN(status=1)
      if ( status < 0 ) then
        write (gol,'("variable `",a,"` has units `",a,"` while expected `",a,"`")') &
                        trim(name), trim(self%var(id)%units), trim(check_units); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! return variable pointer?
    if ( present(V) ) then
      ! assign pointer:
      V => self%var(id)
    end if ! V
    
    ! pointers etc:
    call self%var(id)%Get( status, pdata=pdata, &
                                 Field=Field, &
                                 VectorField=VectorField, &
                                 Constant_Field=Constant_Field, &
                                 Accumulated_Field_Series=Accumulated_Field_Series, &
                                 Instant_Field=Instant_Field, &
                                 Instant_Field_Series=Instant_Field_Series, &
                                 Instant_VectorField=Instant_VectorField )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Variables_GetVariable


  ! ***


  recursive subroutine Variables_SetVariable( self, name, status, enabled )
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  enabled
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_SetVariable'
    
    ! --- local ----------------------------------
    
    class(T_Variable), pointer    ::  V
    integer                       ::  idep
     
    ! --- begin ----------------------------------
    
    ! pointer to variable:
    call self%GetVariable( trim(name), status, V=V )
    IF_NOTOK_RETURN(status=1)
    
    ! enable variable?
    if ( present(enabled) ) then
      !! testing ...
      !write (gol,'("enable `",a,"`")') trim(name); call goPr
      ! set:
      call V%Set( status, enabled=enabled )
      IF_NOTOK_RETURN(status=1)
      ! loop over dependent variables:
      do idep = 1, V%ndep
        ! skip if not in use:
        if ( trim(V%deps(idep)) == 'None' ) cycle
        if ( index(trim(V%deps(idep)),'=') > 0 ) cycle        
        ! recursive call:
        call self%SetVariable( V%deps(idep), status, enabled=enabled )
        IF_NOTOK_RETURN(status=1)
      end do ! dependencies
    end if ! enable flag?
    
    ! ok
    status = 0
    
  end subroutine Variables_SetVariable


  ! ***


  recursive subroutine Variables_Setup( self, varname, tref, t1, t2, status, nodeps )
  
    use GO, only : TDate, operator(+), operator(-), operator(*)

    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    type(TDate), intent(in)                       ::  tref
    type(TDate), intent(in)                       ::  t1, t2
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  nodeps
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Setup'
    
    ! --- local ----------------------------------
    
    type(TDate)                 ::  tmid
    integer                     ::  varid
    class(T_Variable), pointer  ::  V
    logical                     ::  enabled
    logical                     ::  valid
    logical                     ::  with_deps
    integer                     ::  idep
      
    ! --- begin ----------------------------------
    
    ! label:
    call goLabel(rname)
    
    ! mid time:
    tmid = t1 + (t2-t1)*0.5

    ! pointer to variable given the name:
    call self%GetVariable( trim(varname), status, V=V )
    IF_NOTOK_RETURN(status=1)
    
    ! get flag:
    call V%Get( status, enabled=enabled )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    if ( .not. enabled ) then
      ! ok, just leave:
      status=0; return
    end if

    !! info ...
    !write (gol,'(a,": setup variable `",a,"` ...")') rname, trim(varname); call goPr
    
    ! already valid? return status <0 if not:
    call V%Valid( t1, t2, status, silent=.true. )
    IF_ERROR_RETURN(status=1)
    ! already valid?
    if ( status == 0 ) then
      !! info ...
      !write (gol,'(a,":   field already valid for current time step ...")') rname; call goPr
    else

      ! input defined ?
      if ( len_trim(V%input) > 0 ) then
        ! setup variable:
        call V%Setup( tref, t1, t2, status )
        IF_NOTOK_RETURN(status=1)
      end if

      ! apply routine ?
      if ( len_trim(V%routine) > 0 ) then

        ! update depended variables first, unless explicitly requested
        ! not to do that (for example because field depends on previous
        ! fields that should ot be updated to current time) ;
        ! check optional argument:
        with_deps = .true.  ! default with deps
        if ( present(nodeps) ) with_deps = .not. nodeps
        ! update dependencies?
        if ( with_deps ) then
          ! loop over dependencies:
          do idep = 1, V%ndep
            ! skip if dummy:
            if ( trim(V%deps(idep)) == 'None' ) cycle
            if ( index(trim(V%deps(idep)),'=') > 0 ) cycle        
            ! recursive call:
            call Variables_Setup( self, trim(V%deps(idep)), tref, t1, t2, status )
            IF_NOTOK_RETURN(status=1)
          end do ! deps
        end if ! with deps

        !! info ...
        !write (gol,'(a,":   call routine `",a,"` ...")') rname, trim(V%routine); call goPr
        ! switch ...
        select case ( trim(V%routine) )

          !~ dummy ..
          case ( 'Dummy' )
            ! this is only used to define dependencies,
            ! no action needs to be taken ...

          !~ fill area:
          case ( 'Area' )
            ! arguments: ()
            call Variables_Calc_Area( self, trim(varname), status )
            IF_NOTOK_RETURN(status=1)
          
          !~ edge lengths:
          case ( 'EdgeLengthU' )
            ! arguments: ()
            call Variables_Calc_EdgeLengthU( self, trim(varname), status )
            IF_NOTOK_RETURN(status=1)
          case ( 'EdgeLengthV' )
            ! arguments: ()
            call Variables_Calc_EdgeLengthV( self, trim(varname), status )
            IF_NOTOK_RETURN(status=1)
          
          !~ fill normal vectors at u-edges:
          case ( 'UNormal' )
            ! arguments: ()
            call Variables_Calc_UNormal( self, trim(varname), status )
            IF_NOTOK_RETURN(status=1)
          
          !~ fill normal vectors at v-edges:
          case ( 'VNormal' )
            ! arguments: ()
            call Variables_Calc_VNormal( self, trim(varname), status )
            IF_NOTOK_RETURN(status=1)
          
          !~ array like aother
          case ( 'SetValue' )
            ! arguments: ( value=3.14 )
            call Variables_SetValue( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          
          !~ like other variable
          case ( 'Like' )
            ! arguments: ( t, value=3.14 )
            call Variables_Like( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          
          !~ fill lon
          case ('FillLonFromRc' )
            ! arguments: ()
            call Variables_FillLonFromRC( self, status )
            IF_NOTOK_RETURN(status=1)

          !~ fill lat
          case ('FillLatFromRc' )
            ! arguments: ()
            call Variables_FillLatFromRC( self, status )
            IF_NOTOK_RETURN(status=1)

          !~ geopotential (m2/s2) to orography (m)
          case ( 'GpToOro' )
            ! arguments: ( gp )
            call Variables_Calc_GpToOro( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ compute vertex half-level heights:
          case ( 'MixlayerVertexHeights' )
            ! arguments: ( oro_crnr, blh_crnr )
            call Variables_Calc_MixlayerVertexHeights( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !!~ add top layer:
          !case ( 'MixlayerHAltBounds' )
          !  ! arguments: ( halt )
          !  call Variables_Calc_MixlayerHAltBounds( self, trim(varname), V%deps(1:V%ndep), tmid, status )
          !  IF_NOTOK_RETURN(status=1)

          !~ compute hybride half-level pressures:
          case ( 'HybrideHalflevelPressure' )
            ! arguments: ( sp )
            call Variables_Calc_HybrideHalflevelPressure( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ compute layer thickness at corners:
          case ( 'Thickness' )
            ! arguments: ( halt_crnr )
            call Variables_Calc_Thickness( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ compute level top height (above surface)
          case ( 'LevelTopHeight' )
            ! arguments: ( halt )
            call Variables_Calc_LevelTopHeight( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ compute grid cell volumes:
          case ( 'Thickness_to_Volume' )
            ! arguments: ( area, dh )
            call Variables_Calc_Thickness_to_Volume( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ compute average of corner values:
          case ( 'CornerAverage' )
            ! arguments: ( x_crnr )
            call Variables_Calc_CornerAverage( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          case ( 'CornerAverageU' )
            ! arguments: ( x_crnr )
            call Variables_Calc_CornerAverageU( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          case ( 'CornerAverageV' )
            ! arguments: ( x_crnr )
            call Variables_Calc_CornerAverageV( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          ! select subset of input levels:
          case ( 'SelectLevels' )
            ! arguments: ( *_met )
            call Variables_Calc_SelectLevels( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1) 

          !~ compute average of half levels:
          case ( 'HalfLevelAverage' )
            ! arguments: ( hx )
            call Variables_Calc_HalfLevelAverage( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ vertex heights: 
          case ( 'HalflevelAltitudes' )
            ! arguments: ( oro, t_met, q_met, hp_met )
            call Variables_Calc_HalflevelAltitudes( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ vertex heights: 
          case ( 'HalflevelAltitudes_GPH' )
            ! arguments: ( gphbase,gphpert )
            call Variables_Calc_HalflevelAltitudes_GPH( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          !~ potential temperature from temperature:
          case ('T_from_potentialtemperature')
            call T_from_potentialtemperature( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          !~ half level pressures from surface pressure, layer thickness, T and Q:
          case ( 'HalflevelPressures') 
            call Variables_Calc_HalflevelPressures( self,trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ air mass:
          case ( 'AirMass' )
            ! arguments: ( area, hp )
            call Variables_Calc_AirMass( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ air mass:
          case ( 'AirDensity' )
            ! arguments: ( airm, vol )
            call Variables_Calc_AirDensity( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ vector length:
          case ( 'VectorLength' )
            ! arguments: ( u, v )
            call Variables_Calc_VectorLength( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ meteorological wind direction:
          case ( 'MeteoWindDir' )
            ! arguments: ( u, v )
            call Variables_Calc_MeteoWindDir( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !!~ volume flux through u-edges:
          !case ( 'VolumeFlux_u' )
          !  ! arguments: ( u_crnr_met, halt_crnr_met, halt_crnr )
          !  call Variables_Calc_VolumeFlux_u( self, trim(varname), V%deps(1:V%ndep), tmid, status )
          !  IF_NOTOK_RETURN(status=1)
          !!~ volume flux through v-edges:
          !case ( 'VolumeFlux_v' )
          !  ! arguments: ( v_crnr_met, halt_crnr_met, halt_crnr )
          !  call Variables_Calc_VolumeFlux_v( self, trim(varname), V%deps(1:V%ndep), tmid, status )
          !  IF_NOTOK_RETURN(status=1)

          !~ volume flux through u-edges from vector field:
          case ( 'VolumeFlux_uv2u' )
            ! arguments: ( uv_crnr_met, halt_crnr_met, halt_crnr )
            call Variables_Calc_VolumeFlux_uv2u( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          !~ volume flux through v-edges from vector field:
          case ( 'VolumeFlux_uv2v' )
            ! arguments: ( uv_crnr_met, halt_crnr_met, halt_crnr )
            call Variables_Calc_VolumeFlux_uv2v( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ volume flux through w-edges:
          case ( 'VolumeFlux_w' )
            ! arguments: ( uflux, vflux, pvol, nvol )
            call Variables_Calc_VolumeFlux_w( self, trim(varname), V%deps(1:V%ndep), tmid, t1, t2, status )
            IF_NOTOK_RETURN(status=1)
          !~ volume flux through w-edges for constant volume:
          case ( 'VolumeFlux_w0' )
            ! arguments: ( uflux, vflux )
            call Variables_Calc_VolumeFlux_w0( self, trim(varname), V%deps(1:V%ndep), tmid, t1, t2, status )
            IF_NOTOK_RETURN(status=1)

          !~ compute mass average over layers:
          case ( 'LayerAverage' )
            ! arguments: ( hp_met, x_met, hp )
            call Variables_Calc_LayerAverage( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          !~ idem for vectorfield:
          case ( 'LayerAverageUV' )
            ! arguments: ( hp_met, uv_met, hp )
            call Variables_Calc_LayerAverageUV( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ vertical interpolation
          case ( 'LayerInterpol' )
            ! arguments: ( x_in, y_in, x )
            call Variables_Calc_LayerInterpol( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ vertical interpolation, special case for Kz
          case ( 'LayerInterpol_MixLayer_Kz' )
            ! arguments: ( halt_met, kz_met, halt, blh )
            call Variables_Calc_LayerInterpol_MixLayer_Kz( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          case ('SH_from_wvmixingratio')
            call Variables_Calc_SH_from_wvmixingratio( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ relative humidity from temperature and dewpoint:
          case ( 'RelativeHumidityTD' )
            ! arguments: ( T, D )
            call Variables_Calc_RelativeHumidityTD( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ relative humidity from pressure, temperature, and specific humidity:
          case ( 'RelativeHumidityPTQ' )
            ! arguments: ( p, T, Q )
            call Variables_Calc_RelativeHumidityPTQ( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ Total cloud coverage at surface based on clouds on layers:
          case ( 'TotalCloudCoverage' )
            ! arguments: ( cc)
            call Variables_Calc_TotalCloudCoverage( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ windspeed
          case ( 'WindSpeed' )
            ! arguments: ( u, v )
            call Variables_Calc_WindSpeed( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ tendency, used to convert from accumulated to time average fields;
          !  this is available as internal variable for temporal interpolated fields:
          case ( 'Tendency' )
            ! arguments: ( varname )
            call Variables_Calc_Tendency( self, trim(varname), V%deps(1:V%ndep), t1, t2, status )
            IF_NOTOK_RETURN(status=1)

          !~ total rain
          case ( 'TotalRain' )
            ! arguments: ( lsp, cp, ... )
            call Variables_Calc_TotalRain( self, trim(varname), V%deps(1:V%ndep), t1, t2, status )
            IF_NOTOK_RETURN(status=1)

          !~ cloud profiles
          case ( 'CloudProfiles' )
            ! target variables: ( icc, bcc, occ, clwc, praini )
            ! arguments       : ( hp_met, tcc_met, cc_met, clwc_met, hp )
            call Variables_Calc_CloudProfiles( self, V%outs(1:V%nout), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ 3D rain intensity
          case ( 'RainIntensity' )
            ! arguments: ( rain, praini )
            call Variables_Calc_RainIntensity( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ in-cloud (liquid or ice) water content
          case ( 'InCloudWaterContent' )
            ! arguments: ( icc, cwc )
            call Variables_Calc_InCloudWaterContent( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          !~ snow water equivalent
          case ( 'Snowdepth_lw' )
            ! arguments: (physical snow depth in m )
            call Variables_Calc_Snowdepth_lw( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          ! rain rate [m/s] from rain amount [kg/m2/s]
          case ( 'LWE' )
            ! arguments: ( rain_amount )
            call Variables_Calc_LWE( self, varname, V%deps(1:V%ndep), t1, t2, status )
            IF_NOTOK_RETURN(status=1)

          !soil moisture
          case ( 'Extract_swvl' )
            ! arguments: ( swvl, soil layer )
            ! where some might be 'None':
            call Variables_Calc_Extract_swvl( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ volumetric soil water average over 2 or 3 layers:
          case ( 'VolumetricSoilWater' )
            ! arguments: ( swvl1, swvl2, swvl3 )
            ! where some might be 'None':
            call Variables_Calc_VolumetricSoilWater( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ gravimetric soil water:
          case ( 'GravimetricSoilWater' )
            ! arguments: ( swvl1 )
            call Variables_Calc_GravimetricSoilWater( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ leaf area index:
          case ( 'LeafAreaIndex' )
            ! arguments: ()
            call Variables_Calc_LeafAreaIndex( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ surface area index:
          case ( 'SurfaceAreaIndex' )
            ! arguments: ()
            call Variables_Calc_SurfaceAreaIndex( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ surface roughness for mass per landuse class:
          case ( 'SurfaceRoughness_Mass_LU' )
            ! arguments: ( lai, wsurf )
            call Variables_Calc_SurfaceRoughness_Mass_LU( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ surface roughness for heat per landuse class:
          case ( 'SurfaceRoughness_Heat_LU' )
            ! arguments: ( lai, wsurf )
            call Variables_Calc_SurfaceRoughness_Heat_LU( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ canopy height:
          case ( 'CanopyHeight_LU' )
            ! arguments: ( z0h, lai )
            call Variables_Calc_CanopyHeight_LU( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ surface roughness for dust emission landuse class:
          case ( 'SurfaceRoughness_Dust_LU' )
            ! arguments: ( lai )
            call Variables_Calc_SurfaceRoughness_Dust_LU( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          !~ cell averaged surface roughness:
          case ( 'SurfaceRoughness_CellAverage' )
            ! arguments: ( z0_lu )
            call Variables_Calc_SurfaceRoughness_CellAverage( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          !~ soil moisture index:
          case ( 'SoilMoistureIndex' )
            ! arguments: ( swvl, lsm, slt )
            call Variables_Calc_SoilMoistureIndex( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          
          !~ ustar:
          case ( 'UstarJacobsen' )
            ! arguments: ( monin_inv, wsurf )
            call Variables_Calc_UstarJacobsen( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          
          !~ ustar:
          case ( 'UstarGrs2Aver' )
            ! arguments: ( zust, monin_inv )
            call Variables_Calc_UstarGrs2Aver( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          
          !~ monin Obukhov length
          case ( 'MoninObukhovLengthIFS' )
            ! arguments: ustar, t, Q0v 
            call Variables_Calc_MoninObukhovLengthIFS( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
         
          !~ monin Obukhov length
          case ( 'MoninObukhovLengthExposure' )
            ! arguments: radd, wsurf
            call Variables_Calc_MoninObukhovLengthExposure( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
         
          !~ stability class (as also computed in MoninObukhovLengthExposure)
          case ( 'StabilityClass' )
            ! arguments: radd, wsurf, etc
            call Variables_Calc_StabilityClass( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          !~ Inverse of Monin Obukhov length
          case ( 'InverseMonin' )
            ! arguments: monin
            call Variables_Calc_InverseMonin( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          !~ Monin from Inverse of Monin Obukhov length
          case ( 'InverseMoninInverse' )
            ! arguments: monin_inv
            call Variables_Calc_InverseMoninInverse( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          ! ~ Scale Parameter SpecificHumidity (used for calculation of kz)
          case ( 'ScaleParameterSpecificHumidity' )
            ! arguments: sshf, dens, ustar 
            call Variables_Calc_ScaleParameterSpecificHumidity( self, trim(varname), V%deps(1:V%ndep), tmid, status )
           !~ soil moisture index:
          case ( 'SoilMoistureIndex_WRF' )
            ! arguments: ( swvl, lsm, slt )
            call Variables_Calc_SoilMoistureIndex_WRF( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          ! ~ Scale Parameter Dry Static Energy (used for calculation of kz)
          case ( 'ScaleParameterDryStaticEnergy' )
            ! arguments: slhf, t_met, dens, ustar 
            call Variables_Calc_ScaleParameterDryStaticEnergy( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          
          ! ~ Half level gradient
          case ( 'HalfLevelGradient' )
            ! arguments: uu_ifs, dh_met
            call Variables_Calc_HalfLevelGradient( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          
          ! ~ Bottom windspeed
          case ( 'BottomWindSpeed' )
            ! arguments: t, u, Q0v
            call Variables_Calc_BottomWindSpeed( self, trim(varname), V%deps(1:V%ndep), tmid, status )
          !~ Exponential moving average: (temp/ssrd/...)
          case ( 'CalcEMA' )
            ! arguments: ( varname, alpha=0.1 )
            call Variables_CalcEMA( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          case ( 'ChangeSign' )
            ! arguments: ( varname)
            call Variables_ChangeSign( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)

          ! ~ Virtual temperature flux
          case ( 'VirtualTemperatureFlux' )
            ! arguments: ustar, sstar, t_met, qstar
            call Variables_Calc_VirtualTemperatureFlux( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          ! ~ Dry static energy
          case ( 'DryStaticEnergy' )
            ! arguments: h, t, tsurf
            call Variables_Calc_DryStaticEnergy( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
          
          !~ Richardson number
          case ( 'RichardsonNumber' )
            ! arguments: h,uv,q,t,dse
            call Variables_Calc_RichardsonNumber( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
 
          !~ Half level kz-ifs
          case ( 'HalfLevelKzIFS' )
            ! arguments: temp, dse, uu_ifs, duudz, ri, monin, halt, blh, ustar
            call Variables_Calc_HalfLevelKzIFS( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          !~ Half level kz-msp
          case ( 'HalfLevelKzMSP' )
            ! arguments: h, monin_inv, blh, ustar
            call Variables_Calc_HalfLevelKzMSP( self, trim(varname), V%deps(1:V%ndep), tmid, status )
            IF_NOTOK_RETURN(status=1)
            
          !~ not yet ..
          case default
            write (gol,'("unsupported routine `",a,"` for variable `",a,"`")') &
                    trim(V%routine), trim(varname); call goErr
            TRACEBACK; status=1; return
        end select

      end if ! apply routine
      
      ! truncate if filled:
      if ( (len_trim(V%input) > 0) .or. (len_trim(V%routine) > 0) ) then
        ! truncate if necessary:
        call V%Truncate( status )
        IF_NOTOK_RETURN(status=1)
      end if
      
    end if ! valid?
    
    ! unlabel:
    call goLabel()
    
    ! ok
    status = 0
    
  end subroutine Variables_Setup


  ! ***


  subroutine Variables_GetValidPtr( self, varname, t1, t2, pdata, status, &
                                     units, lbo, ubo, &
                                     check_units, check_lbo, check_ubo )
  
    use GO            , only : TDate
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    type(TDate), intent(in)                       ::  t1, t2
    real, pointer, intent(out)                    ::  pdata(:,:,:)
    integer, intent(out)                          ::  status
    character(len=*), intent(out), optional       ::  units
    integer, intent(out), optional                ::  lbo(3), ubo(3)
    character(len=*), intent(in), optional        ::  check_units
    integer, intent(in), optional                 ::  check_lbo(3), check_ubo(3)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_GetValidPtr'
    
    ! --- local ----------------------------------
    
    class(T_Variable), pointer           ::  V
      
    ! --- begin ----------------------------------
    
    ! get pointer to variable given name:
    call self%GetVariable( trim(varname), status, V=V )
    IF_NOTOK_RETURN(status=1)
    
    ! should be valid:
    call V%Valid( t1, t2, status )
    IF_NOTOK_RETURN(status=1)
    
    ! get pointer to data etc:
    call V%Get( status, pdata=pdata, &
                  units=units, lbo=lbo, ubo=ubo, &
                  check_units=check_units, check_lbo=check_lbo, check_ubo=check_ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_GetValidPtr


  ! ***


  subroutine Variables_GetValidPtr_uv( self, varname, t1, t2, pudata, pvdata, status, &
                                     units, lbo, ubo, &
                                     check_units, check_lbo, check_ubo )
  
    use GO            , only : TDate
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    type(TDate), intent(in)                       ::  t1, t2
    real, pointer, intent(out)                    ::  pudata(:,:,:)
    real, pointer, intent(out)                    ::  pvdata(:,:,:)
    integer, intent(out)                          ::  status
    character(len=*), intent(out), optional       ::  units
    integer, intent(out), optional                ::  lbo(3), ubo(3)
    character(len=*), intent(in), optional        ::  check_units
    integer, intent(in), optional                 ::  check_lbo(3), check_ubo(3)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_GetValidPtr_uv'
    
    ! --- local ----------------------------------
    
    class(T_Variable), pointer           ::  V
      
    ! --- begin ----------------------------------
    
    ! get pointer to variable given name:
    call self%GetVariable( trim(varname), status, V=V )
    IF_NOTOK_RETURN(status=1)
    
    ! should be valid:
    call V%Valid( t1, t2, status )
    IF_NOTOK_RETURN(status=1)
    
    ! get pointer to data etc:
    call V%Get( status, pudata=pudata, pvdata=pvdata, &
                  units=units, lbo=lbo, ubo=ubo, &
                  check_units=check_units, check_lbo=check_lbo, check_ubo=check_ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_GetValidPtr_uv

  
  ! ***
  
  subroutine Variables_FillLonFromRC( self, status )
    
    use GO     , only : T_Field
    use C3PO   , only : T_Grid_Ugg
    use LE_Grid, only : ugg, ugg_bnds
    
    ! ---in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_FillLonFromRC'
    
    ! --- local ----------------------------------
    
    class(T_Field), pointer     ::  F
    character(len=32)           ::  gridtype
    type(T_Grid_Ugg), pointer   ::  gri
    real, allocatable           ::  values(:,:,:)
    integer                     ::  i,j
      
    ! --- begin ----------------------------------
    
    ! pointer to variable given the name:
    call self%GetVariable( 'lon', status, Field=F, check_units='degrees_east', &
                             gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! storage:
    allocate( values(gri%nlon,gri%nlat,1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    do j = 1, gri%nlat
      do i = 1, gri%nlon
        values(i,j,1) = gri%longitude(i,j)
      end do
    end do

    ! store:
    call F%Put( values, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
  end subroutine Variables_FillLonFromRC  


  ! ***
  
  
  subroutine Variables_FillLatFromRC( self, status )
    
    use GO     , only : T_Field
    use C3PO   , only : T_Grid_Ugg
    use LE_Grid, only : ugg, ugg_bnds
    
    ! ---in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_FillLatFromRC('
    
    ! --- local ----------------------------------
    
    class(T_Field), pointer     ::  F
    character(len=32)           ::  gridtype
    type(T_Grid_Ugg), pointer   ::  gri
    real, allocatable           ::  values(:,:,:)
    integer                     ::  i,j
      
    ! --- begin ----------------------------------
    
    ! pointer to variable given the name:
    call self%GetVariable( 'lat', status, Field=F, check_units='degrees_north', &
                             gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! storage:
    allocate( values(gri%nlon,gri%nlat,1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    do i = 1, gri%nlon
      do j = 1, gri%nlat
        values(i,j,1) = gri%latitude(i,j)
      end do
    end do

    ! store:
    call F%Put( values, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
  end subroutine Variables_FillLatFromRC  


  ! ***

  subroutine Variables_Calc_Area( self, varname, status )
  
    use GO     , only : TDate
    use GO     , only : T_Field
    use C3PO   , only : T_Grid_Ugg
    use LE_Grid, only : ugg, ugg_bnds
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_Area'
    
    ! --- local ----------------------------------
    
    class(T_Field), pointer     ::  F
    character(len=32)           ::  gridtype
    type(T_Grid_Ugg), pointer   ::  gri
    real, allocatable           ::  values(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! pointer to variable given the name:
    call self%GetVariable( trim(varname), status, Field=F, check_units='m2', &
                             gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! switch:
    select case ( trim(gridtype) )
      case ( 'cell' )
        gri => ugg
      case ( 'cell_bnds' )
        gri => ugg_bnds
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select

    ! storage:
    allocate( values(gri%nlon,gri%nlat,1), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill with area:
    call gri%AreaOper( values, '=', 'm2', status )
    IF_NOTOK_RETURN(status=1)
    ! store:
    call F%Put( values, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_Area


  ! ***


  !
  ! Fill field with edge lengths
  !

  subroutine Variables_Calc_EdgeLengthU( self, varname, status )
  
    use GO     , only : T_Field
    use C3PO   , only : T_Grid_Ugg, EDGE_LEFT, EDGE_RIGHT
    use LE_Grid, only : ugg
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_EdgeLengthU'
    
    ! --- local ----------------------------------
    
    class(T_Field), pointer             ::  F
    character(len=64)                   ::  gridtype
    type(T_Grid_Ugg), pointer           ::  gri
    real, allocatable                   ::  dy(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Field=F, &
                              check_units='m', gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! grid definition of cells:
    select case ( trim(gridtype) )
      case ( 'u-edge' )
        gri => ugg
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! storage for edge lengths:
    allocate( dy(0:gri%nlon,1:gri%nlat,1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over lat bands:
    do j = 1, gri%nlat
      ! first:
      i = 1
      call gri%GetEdgeLength( i, j, EDGE_LEFT, dy(i-1,j,1), status )
      IF_NOTOK_RETURN(status=1)
      ! right edges:
      do i = 1, gri%nlon
        call gri%GetEdgeLength( i, j, EDGE_RIGHT, dy(i,j,1), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j
    
    ! store:
    call F%Put( dy, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( dy, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_EdgeLengthU
  
  ! *

  subroutine Variables_Calc_EdgeLengthV( self, varname, status )
  
    use GO     , only : T_Field
    use C3PO   , only : T_Grid_Ugg, EDGE_LOWER, EDGE_UPPER
    use LE_Grid, only : ugg
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_EdgeLengthV'
    
    ! --- local ----------------------------------
    
    class(T_Field), pointer             ::  F
    character(len=64)                   ::  gridtype
    type(T_Grid_Ugg), pointer           ::  gri
    real, allocatable                   ::  dx(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Field=F, &
                              check_units='m', gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! grid definition of cells:
    select case ( trim(gridtype) )
      case ( 'v-edge' )
        gri => ugg
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! storage for vector components:
    allocate( dx(1:gri%nlon,0:gri%nlat,1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over lon bands:
    do i = 1, gri%nlon
      ! first:
      j = 1
      call gri%GetEdgeLength( i, j, EDGE_LOWER, dx(i,j-1,1), status )
      IF_NOTOK_RETURN(status=1)
      ! upper edges:
      do j = 1, gri%nlat
        call gri%GetEdgeLength( i, j, EDGE_UPPER, dx(i,j,1), status )
        IF_NOTOK_RETURN(status=1)
      end do ! j
    end do ! i
    
    ! store:
    call F%Put( dx, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( dx, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_EdgeLengthV


  ! ***


  !
  ! Fill vector field with normal vectors at mid of u-edges.
  !

  subroutine Variables_Calc_UNormal( self, varname, status )
  
    use GO     , only : T_VectorField
    use C3PO   , only : T_Grid_Ugg
    use LE_Grid, only : ugg
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_UNormal'
    
    ! --- local ----------------------------------
    
    class(T_VectorField), pointer       ::  VF
    character(len=64)                   ::  gridtype
    type(T_Grid_Ugg), pointer           ::  gri
    real, allocatable                   ::  n_lon(:,:,:)
    real, allocatable                   ::  n_lat(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, VectorField=VF, &
                              check_units='1', lbo=lbo, ubo=ubo, &
                              gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! grid definition of cells:
    select case ( trim(gridtype) )
      case ( 'u-edge' )
        gri => ugg
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! storage for vector components:
    allocate( n_lon(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( n_lat(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! copy:
    n_lon(:,:,1) = gri%unormal(1,:,:)
    n_lat(:,:,1) = gri%unormal(2,:,:)
    
    ! store:
    call VF%Put( n_lon, n_lat, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( n_lon, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( n_lat, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_UNormal


  ! ***


  !
  ! Fill vector field with normal vectors at mid of v-edges.
  !

  subroutine Variables_Calc_VNormal( self, varname, status )
  
    use GO     , only : T_VectorField
    use C3PO   , only : T_Grid_Ugg
    use LE_Grid, only : ugg
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_VNormal'
    
    ! --- local ----------------------------------
    
    class(T_VectorField), pointer       ::  VF
    character(len=64)                   ::  gridtype
    type(T_Grid_Ugg), pointer           ::  gri
    real, allocatable                   ::  n_lon(:,:,:)
    real, allocatable                   ::  n_lat(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, VectorField=VF, &
                              check_units='1', lbo=lbo, ubo=ubo, &
                              gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! grid definition of cells:
    select case ( trim(gridtype) )
      case ( 'v-edge' )
        gri => ugg
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! storage for vector components:
    allocate( n_lon(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( n_lat(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! copy:
    n_lon(:,:,1) = gri%vnormal(1,:,:)
    n_lat(:,:,1) = gri%vnormal(2,:,:)
    
    ! store:
    call VF%Put( n_lon, n_lat, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( n_lon, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( n_lat, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_VNormal

  ! ***

  !
  ! Set constant value.
  ! Arguments:
  !   input :  'value=0.0'
  !

  subroutine Variables_SetValue( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use GO          , only : GoVarValue
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_SetValue'
    
    ! argument names:
    integer, parameter        ::  i_value = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    real                                ::  value
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! read value
    call GoVarValue( args(i_value), ';', 'value','=', value, status )
    if ( status /= 0 ) then
      write (gol,'("Format should be `value=3.14`, found `",a,"`")') trim(args(i_value)) ; call GoErr
      TRACEBACK;status=1;return
    end if
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    values = value

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_SetValue
  
  ! ***

  !
  ! Fille variable like other variable.
  ! Arguments:
  !   input :  'vname', 'value=0.0'
  !

  subroutine Variables_Like( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use GO          , only : GoVarValue
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Like'
    
    ! argument names:
    integer, parameter        ::  i_var   = 1
    integer, parameter        ::  i_value = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  source(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    real                                ::  value
      
    ! --- begin ----------------------------------
    
    ! dependencies:
    call self%GetValidPtr( args(i_var), t,t, source, status, &
                               lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! read value
    call GoVarValue( args(i_value), ';', 'value','=', value, status )
    if ( status /= 0 ) then
      write (gol,'("Format should be `value=3.14`, found `",a,"`")') trim(args(i_value)) ; call GoErr
      TRACEBACK;status=1;return
    end if
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    values = value

    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF )
    IF_NOTOK_RETURN(status=1)    
    ! store:
    call iF%PutSample( values, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Like

  ! ***


  !
  ! oro = gp/grav
  !

  subroutine Variables_Calc_GpToOro( self, varname, args, t, status )
  
    use Binas, only : grav
    use GO   , only : TDate
    use GO   , only : T_Constant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_GpToOro'
    
    ! argument names:
    integer, parameter        ::  i_gp = 1
    
    ! --- local ----------------------------------
    
    class(T_Constant_Field), pointer    ::  cF
    real, pointer                       ::  gp(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Constant_Field=cF, &
                              check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_gp), t,t, gp, status, check_units='m2/s2', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! convert:
    values = gp / grav

    ! store:
    call cF%PutSample( values, (/t,t/), status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Variables_Calc_GpToOro


  ! ***


  !
  ! Compute vertex heights:
  !   target variable : halt (m)
  !   input variables :  oro (m), blh (m)
  !

  subroutine Variables_Calc_MixlayerVertexHeights( self, varname, args, t, status )
  
    use GO     , only : TDate
    use GO     , only : T_Instant_Field
    use Num    , only : StdDev
    
    use LE_Data_Common, only : nlev
    use LE_Data_Common, only : mixlayer_surf_top
    use LE_Data_Common, only : mixlayer_mix_topmin
    use LE_Data_Common, only : mixlayer_dmin, mixlayer_sdofac
    use LE_Data_Common, only : mixlayer_top
    use LE_Data_Common, only : mixlayer_daloft
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_MixlayerVertexHeights'
    
    ! argument names:
    integer, parameter        ::  i_oro = 1
    integer, parameter        ::  i_blh = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  oro(:,:,:)
    real, pointer                       ::  blh(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real                                ::  mu
    real, allocatable                   ::  sdo (:,:)
    real, allocatable                   ::  dmin(:,:)
    real, allocatable                   ::  halt(:,:,:)
    integer                             ::  i, j
    integer                             ::  ilev
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! check ...
    if ( (lbo(3) /= 0) .or. ((ubo(3) /= nlev) .and. (ubo(3) /= nlev+1)) ) then
      write (gol,'("expected vertical index space 0:",i0," or 0:",i0," for `",a,"`, found ",i0,":",i0)') &
                nlev, nlev+1, trim(varname), lbo(3), ubo(3); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! index space for surface fields:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)
    
    ! dependencies:
    call self%GetValidPtr( args(i_oro), t,t, oro, status, check_units='m', &
                              check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_blh), t,t, blh, status, check_units='m', &
                              check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)

    ! storage:
    allocate( sdo(lbo(1):ubo(1),lbo(2):ubo(2)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( dmin(lbo(1):ubo(1),lbo(2):ubo(2)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( halt(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! reset boundary layer height in regions with strong orography;
    ! init std.dev. of oro:
    sdo = 0.0
    ! loop over cells inside boundary:
    do j = lbo(2)+1, ubo(2)-1
      do i = lbo(1)+1, ubo(1)-1
        ! compute std.dev. over 3x3 stencil:
        call StdDev( oro(i-1:i+1,j-1:j+1,1), sdo(i,j), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j
    
    ! loop over cells in target space:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
      
        ! higher minium above rough orography:
        dmin = mixlayer_dmin + mixlayer_sdofac(4) * sdo
        ! init with orography:
        halt(i,j,0) = oro(i,j,1)
        ! surface layer:
        halt(i,j,1) = oro(i,j,1) + mixlayer_surf_top
        ! mix layer with minium top:
        halt(i,j,2) = max( oro(i,j,1) + mixlayer_mix_topmin, oro(i,j,1) + blh(i,j,1) )
        ! top of second residual layer, minium thickness for each:
        halt(i,j,4) = max( halt(i,j,2) + 2*dmin(i,j), mixlayer_top(4) )
        ! top of first residual layer is in between:
        halt(i,j,3) = 0.5*( halt(i,j,2) + halt(i,j,4) )
        ! extra layers:
        if ( nlev > 4 ) then
          ! loop:
          do ilev = 5, nlev
            ! higher minium above rough orography:
            dmin = mixlayer_dmin + mixlayer_sdofac(ilev) * sdo
            ! add layer:
            halt(i,j,ilev) = max( halt(i,j,ilev-1) + dmin(i,j), mixlayer_top(ilev) )
          end do ! levels
        end if ! nlev
        ! aloft?
        if ( ubo(3) > nlev ) then
          ! extra layer:
          halt(i,j,nlev+1) = halt(i,j,nlev) + mixlayer_daloft
        end if
        
      end do ! i
    end do ! j
    
    ! store:
    call iF%PutSample( halt, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( sdo, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( dmin, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( halt, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_MixlayerVertexHeights


!  ! ***
!
!
!  !
!  ! Add boundary values to half-level altitudes;
!  ! extra layer added on top,
!  ! horizontal boundaries are copied from core.
!  !   target variable : halt_bnds (m)
!  !   input variables : halt (m)
!  !
!
!  subroutine Variables_Calc_MixlayerHAltBounds( self, varname, args, t, status )
!  
!    use GO     , only : TDate
!    use GO     , only : T_Instant_Field
!    
!    use LE_Data_Common, only : nlev
!    use LE_Data_Common, only : mixlayer_dtop
!    
!    ! --- in/out ---------------------------------
!    
!    class(T_Variables), intent(inout)             ::  self
!    character(len=*), intent(in)                  ::  varname
!    character(len=*), intent(in)                  ::  args(2)
!    type(TDate), intent(in)                       ::  t
!    integer, intent(out)                          ::  status
!  
!    ! --- const ----------------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_MixlayerHAltBounds'
!    
!    ! argument names:
!    integer, parameter        ::  i_halt = 1
!    
!    ! --- local ----------------------------------
!    
!    class(T_Instant_Field), pointer     ::  iF
!    real, pointer                       ::  halt(:,:,:)
!    integer                             ::  lbo(3), ubo(3)
!    real, allocatable                   ::  halt_bnds(:,:,:)
!      
!    ! --- begin ----------------------------------
!    
!    ! pointer to target field given the name:
!    call self%GetVariable( trim(varname), status, Instant_Field=iF, check_units='m' )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! dependencies:
!    call self%GetValidPtr( args(i_halt), t,t, halt, status, check_units='m' )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! dimensions:
!    lbo = lbound( halt )
!    ubo = ubound( halt )
!
!    ! storage with halo of 1 cell:
!    allocate( halt_bnds(lbo(1)-1:ubo(1)+1,lbo(2)-1:ubo(2)+1,0:nlev+1), stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! copy levels form core:
!    halt_bnds(lbo(1):ubo(1),lbo(2):ubo(2),0:nlev) = halt
!    ! add top:
!    halt_bnds(:,:,nlev+1) = halt_bnds(:,:,nlev) + mixlayer_dtop
!    ! copy boundaries:
!    halt_bnds(lbo(1)-1,:,:) = halt_bnds(lbo(1),:,:)  ! west
!    halt_bnds(ubo(1)+1,:,:) = halt_bnds(ubo(1),:,:)  ! east
!    halt_bnds(:,lbo(2)-1,:) = halt_bnds(:,lbo(2),:)  ! south
!    halt_bnds(:,ubo(2)+1,:) = halt_bnds(:,ubo(2),:)  ! north
!    
!    ! store:
!    call iF%PutSample( halt_bnds, t, status )
!    IF_NOTOK_RETURN(status=1)
!
!    ! clear:
!    deallocate( halt_bnds, stat=status )
!    IF_NOTOK_RETURN(status=1)
!    
!    ! ok
!    status = 0
!    
!  end subroutine Variables_Calc_MixlayerHAltBounds


  ! ***


  !
  ! Compute halflevel pressure:
  !   target variable : hp (Pa)
  !   input variables : sp (Pa)
  !

  subroutine Variables_Calc_HybrideHalflevelPressure( self, varname, args, t, status )
  
    use GO              , only : TDate
    use GO              , only : T_Instant_Field
    use LE_Data_Common  , only : nlev, nlev_top
    use LE_Data_Common  , only : hyb, hyb_top
    use LE_Data_Variable, only : LEN_LEVTYPE
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_HybrideHalflevelPressure'
    
    ! argument names:
    integer, parameter        ::  i_sp = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    character(len=LEN_LEVTYPE)          ::  levtype
    real, pointer                       ::  sp(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  hp(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             levtype=levtype, check_units='Pa' )
    IF_NOTOK_RETURN(status=1)

    ! dependencies:
    call self%GetValidPtr( args(i_sp), t,t, sp, status, check_units='Pa' )
    IF_NOTOK_RETURN(status=1)
    
    ! dimensions, used for horizontal grid:
    lbo = lbound( sp )
    ubo = ubound( sp )

    ! should be half levels, eventually incl top:
    select case ( trim(levtype) )
      !~ bounds of model levels:
      case ( 'halflevels' )
        ! storage:
        allocate( hp(lbo(1):ubo(1),lbo(2):ubo(2),0:nlev), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! half-level pressure from surface pressure and hybride coeff:
        call hyb%Ps_to_Ph( sp(:,:,1), hp, status )
        IF_NOTOK_RETURN(status=1)
      !~ bounds of model levels incl top layers:
      case ( 'halflevels_top' )
        ! storage:
        allocate( hp(lbo(1):ubo(1),lbo(2):ubo(2),0:nlev_top), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! half-level pressure from surface pressure and hybride coeff:
        call hyb_top%Ps_to_Ph( sp(:,:,1), hp, status )
        IF_NOTOK_RETURN(status=1)
      !~
      case default
        write (gol,'("unsupported levtype `",a,"`")') trim(levtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! store:
    call iF%PutSample( hp, t, status, lbo=lbound(hp), ubo=ubound(hp) )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( hp, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_HybrideHalflevelPressure


  ! ***


  !
  ! Thickness from half-level altitudes ('halt'):
  !
  !   call Variables_Calc_Thickness( self, 'dh', (/'halt'/), t, status )
  !

  subroutine Variables_Calc_Thickness( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_Thickness'
    
    ! argument names:
    integer, parameter        ::  i_halt = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  halt(:,:,:)  ! (:,:,0:nz)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_hlv(3), ubo_hlv(3)
    real, allocatable                   ::  dh(:,:,:)
    integer                             ::  ilev
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! dependencies:
    call self%GetValidPtr( args(i_halt), t,t, halt, status, check_units='m', &
                             lbo=lbo_hlv, ubo=ubo_hlv )
    IF_NOTOK_RETURN(status=1)
    
    ! target shape
    lbo = lbo_hlv+(/0,0,1/)
    ubo = ubo_hlv
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='m' )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( dh(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over levels:
    do ilev = lbo(3), ubo(3)
      ! loop over cells:
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! fill:
          dh(i,j,ilev) = abs( halt(i,j,ilev) - halt(i,j,ilev-1) )
        end do ! i
      end do ! j
    end do ! levels
    
    ! store:
    call iF%PutSample( dh, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( dh, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_Thickness


  ! ***


  !
  ! Level top height above surface ('h'):
  !
  !   call Variables_Calc_LevelTopHeight( self, 'h', (/'halt'/), t, status )
  !

  subroutine Variables_Calc_LevelTopHeight( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_LevelTopHeight'
    
    ! argument names:
    integer, parameter        ::  i_halt = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  halt(:,:,:)  ! (:,:,0:nz)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  h(:,:,:)
    integer                             ::  ilev
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_halt), t,t, halt, status, check_units='m', &
                              check_lbo=lbo+(/0,0,-1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for level field:
    allocate( h(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over levels:
    do ilev = lbo(3), ubo(3)
      ! loop over cells in target space:
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! fill:
          h(i,j,ilev) = halt(i,j,ilev) - halt(i,j,lbo(3)-1)
        end do ! i
      end do ! j
    end do ! levels
    
    ! store:
    call iF%PutSample( h, t, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( h, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_LevelTopHeight


  ! ***


  !
  ! Cell average ('x') from corners ('x_crnr') :
  !
  !   call Variables_Calc_CornerAverage( self, 'x', (/'x_crnr'/), t, status )
  !

  subroutine Variables_Calc_CornerAverage( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_CornerAverage'
    
    ! argument names:
    integer, parameter        ::  i_crnr = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  values_crnr(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field, get size (number of levels might be undefined yet):
    call self%GetVariable( varname, status, Instant_Field=iF, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! set level index space if not defined yet:
    if ( ubo(3) < lbo(3) ) then
      ! pointer to source array:
      call self%GetValidPtr( args(i_crnr), t,t, values_crnr, status )
      IF_NOTOK_RETURN(status=1)
      ! copy level bounds:
      lbo(3) = lbound(values_crnr,3)
      ubo(3) = ubound(values_crnr,3)
    end if
    
    ! pointer to source array, now check bounds:
    call self%GetValidPtr( args(i_crnr), t,t, values_crnr, status, &
                             check_lbo=lbo+(/-1,-1,0/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for corner averages:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! average:
        values(i,j,:) = 0.25 * ( values_crnr(i-1,j-1,:) + &
                                 values_crnr(i  ,j-1,:) + &
                                 values_crnr(i-1,j  ,:) + &
                                 values_crnr(i  ,j  ,:)     )
      end do ! i
    end do ! j
    
    ! store, provide bounds:
    call iF%PutSample( values, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_CornerAverage

  !
  ! U-edge average ('x') from corners ('x_crnr') :
  !
  !   call Variables_Calc_CornerAverageU( self, 'x', (/'x_crnr'/), t, status )
  !

  subroutine Variables_Calc_CornerAverageU( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_CornerAverageU'
    
    ! argument names:
    integer, parameter        ::  i_crnr = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  values_crnr(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field, get size (number of levels might be undefined yet):
    call self%GetVariable( varname, status, Instant_Field=iF, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! set level index space if not defined yet:
    if ( ubo(3) < lbo(3) ) then
      ! pointer to source array:
      call self%GetValidPtr( args(i_crnr), t,t, values_crnr, status )
      IF_NOTOK_RETURN(status=1)
      ! copy level bounds:
      lbo(3) = lbound(values_crnr,3)
      ubo(3) = ubound(values_crnr,3)
    end if
    
    ! pointer to source array, check bounds for corners of u-edge:
    call self%GetValidPtr( args(i_crnr), t,t, values_crnr, status, &
                              check_lbo=lbo+(/-1,0,0/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for corner averages:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! average of two corners:
        values(i,j,:) = 0.5 * ( values_crnr(i-1,j,:) + values_crnr(i,j,:) )
      end do ! i
    end do ! j
    
    ! store, provide bounds:
    call iF%PutSample( values, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_CornerAverageU

  !
  ! U-edge average ('x') from corners ('x_crnr') :
  !
  !   call Variables_Calc_CornerAverageV( self, 'x', (/'x_crnr'/), t, status )
  !

  subroutine Variables_Calc_CornerAverageV( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_CornerAverageV'
    
    ! argument names:
    integer, parameter        ::  i_crnr = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  values_crnr(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field, get size (number of levels might be undefined yet):
    call self%GetVariable( varname, status, Instant_Field=iF, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! set level index space if not defined yet:
    if ( ubo(3) < lbo(3) ) then
      ! pointer to source array:
      call self%GetValidPtr( args(i_crnr), t,t, values_crnr, status )
      IF_NOTOK_RETURN(status=1)
      ! copy level bounds:
      lbo(3) = lbound(values_crnr,3)
      ubo(3) = ubound(values_crnr,3)
    end if
    
    ! pointer to source array, check bounds for corners of u-edge:
    call self%GetValidPtr( args(i_crnr), t,t, values_crnr, status, &
                              check_lbo=lbo+(/0,-1,0/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for corner averages:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! average of two corners:
        values(i,j,:) = 0.5 * ( values_crnr(i,j-1,:) + values_crnr(i,j,:) )
      end do ! i
    end do ! j
    
    ! store, provide bounds:
    call iF%PutSample( values, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_CornerAverageV


  ! ***


  !
  ! Select levels from input fields and copy into target field:
  !
  !   call Variables_SelectLevels( self, 'x_met', (/'x'/), t, status )
  !

  subroutine Variables_Calc_SelectLevels( self, varname, args, t, status )
  
    use GO            , only : TDate
    use GO            , only : T_Instant_Field
    use LE_Data_Common, only : metlevel_select, metlevel_select_top
    use LE_Data_Variable, only : LEN_LEVTYPE
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_SelectLevels'
    
    ! argument names:
    integer, parameter        ::  i_arg = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    character(len=LEN_LEVTYPE)          ::  levtype
    integer, pointer                    ::  mselect(:)
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  values(:,:,:)
    character(len=32)                   ::  units
    real, pointer                       ::  values_in(:,:,:)
    integer                             ::  i, j
    integer                             ::  k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( varname, status, Instant_Field=iF, &
                             levtype=levtype, units=units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( values, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! dependencies:
    call self%GetValidPtr( args(i_arg), t,t, values_in, status, &
                               check_units=units )
    IF_NOTOK_RETURN(status=1)

    ! should be half levels, eventually incl top:
    select case ( trim(levtype) )
      !~ bounds of model levels:
      case ( 'halflevels' )
        ! check ...
        if ( .not. allocated(metlevel_select) ) then
          write (gol,'("array `metlevel_select` not allocated")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! set pointer:
        mselect => metlevel_select
      !~ bounds of model levels incl top layers:
      case ( 'halflevels_top' )
        ! check ...
        if ( .not. allocated(metlevel_select_top) ) then
          write (gol,'("array `metlevel_select_top` not allocated")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! set pointer:
        mselect => metlevel_select_top
      !~
      case default
        write (gol,'("unsupported levtype `",a,"`")') trim(levtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! check ...
    if ( (lbo(3) < lbound(mselect,1)) .or. &
         (ubo(3) > ubound(mselect,1)) ) then
      write (gol,'("target level range (",i0,":",i0,") should be within metlevel indices (",i0,":",i0,")")') &
              lbo(3), ubo(3), lbound(mselect,1), ubound(mselect,1); call goErr
      TRACEBACK; status=1; return
    end if

    ! check ...
    if ( (mselect(lbo(3)) < lbound(values_in,3)) .or. &
         (mselect(ubo(3)) > ubound(values_in,3)) ) then
      write (gol,'("metlevel selection (",i0,":",i0,") should be within source range (",i0,":",i0,")")') &
              mselect(lbo(3)), mselect(ubo(3)), &
              lbound(values_in,3), ubound(values_in,3); call goErr
      TRACEBACK; status=1; return
    end if
    
    !! info ...
    !write (gol,'("    xxx copy levels ...")'); call goPr
    ! loop over target levels:
    do k = lbo(3), ubo(3)
      !! testing ...
      !write (gol,'("      x fill level ",i2," from ",i2)') k, mselect(k); call goPr
      ! copy:
      values(lbo(1):ubo(1),lbo(2):ubo(2),k) = values_in(lbo(1):ubo(1),lbo(2):ubo(2),mselect(k))
    end do
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_SelectLevels


  ! ***
  !
  ! Cell average ('x') from half-levels ('hx') :
  !
  !   call Variables_Calc_HalfLevelAverage( self, 'x', (/'hx'/), t, status )
  !

  subroutine Variables_Calc_HalfLevelAverage( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_HalfLevelAverage'
    
    ! argument names:
    integer, parameter        ::  i_hlv = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  values_hlv(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
    integer                             ::  k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_Field=iF, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_hlv), t,t, values_hlv, status, &
                              check_lbo=lbo+(/0,0,-1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for full-levels:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over levels:
    do k = lbo(3), ubo(3)
      ! loop over cells:
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! average:
          values(i,j,k) = 0.5 * ( values_hlv(i,j,k-1) + values_hlv(i,j,k) )
        end do ! i
      end do ! j
    end do ! k
    
    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_HalfLevelAverage


  ! ***


  !
  ! Cell volume ('vol') from area ('area') and cell heights ('dh') :
  !
  !   call Variables_Calc_Thickness_to_Volume( self, 'vol', (/'area','dh'/), t, status )
  !

  subroutine Variables_Calc_Thickness_to_Volume( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_Thickness_to_Volume'
    
    ! argument names:
    integer, parameter        ::  i_area = 1
    integer, parameter        ::  i_dh   = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  area(:,:,:)
    real, pointer                       ::  dh(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  volume(:,:,:)
    integer                             ::  i, j
    integer                             ::  ilev
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              lbo=lbo, ubo=ubo, check_units='m3' )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for surface fields:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)

    ! dependencies:
    call self%GetValidPtr( args(i_area), t,t, area, status, check_units='m2', &
                              check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_dh  ), t,t, dh  , status, check_units='m', &
                              check_lbo=lbo, check_ubo=ubo  )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( volume(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over levels:
    do ilev = lbo(3), ubo(3)
      ! loop over target cells:
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! fill:
          volume(i,j,ilev) = area(i,j,1) * dh(i,j,ilev)
        end do ! i
      end do ! j
    end do ! levels
    
    ! store:
    call iF%PutSample( volume, t, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( volume, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_Thickness_to_Volume


  ! ***


  !
  ! Vertex half-level heights ('halt') from orography ('oro'), temperature ('t'),
  ! specific humidity ('q'), and half-level pressures ('hp') :
  !
  !   call Variables_Calc_Thickness( self, 'vol', (/'oro','t','q','hp'/), t, status )
  !
  ! Default values if name is 'None':
  !    temperature 293 K
  !    specific humidity 0 kg/kg
  !

  subroutine Variables_Calc_HalflevelAltitudes( self, varname, args, t, status )
  
    use GO  , only : TDate
    use GO  , only : T_Instant_Field
    use JAQL, only : GeoPotentialHeightB
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(4)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_HalflevelAltitudes'
    
    ! argument names:
    integer, parameter        ::  i_oro    = 1
    integer, parameter        ::  i_temper = 2
    integer, parameter        ::  i_shumid = 3
    integer, parameter        ::  i_hpress = 4
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  oro(:,:,:)
    real, pointer                       ::  temper(:,:,:)
    real, pointer                       ::  shumid(:,:,:)
    real, pointer                       ::  hpress(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  halt(:,:,:)
    integer                             ::  i, j
    logical                             ::  change_levs
      
    ! --- begin ----------------------------------
    
    ! get pointer to half-level pressure to obtain number of levels if necessary:
    call self%GetValidPtr( args(i_hpress), t,t, hpress, status )
    IF_NOTOK_RETURN(status=1)

    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! set level index space if neceesary:
    if ( lbo(3) /= lbound(hpress,3) .or. ubo(3) /= ubound(hpress,3) ) then
      lbo(3) = lbound(hpress,3)
      ubo(3) = ubound(hpress,3)
      change_levs = .true.
    else 
      change_levs = .false.
    end if
    
    ! index space for surface fields:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)
    
    ! arguments:
    call self%GetValidPtr( args(i_oro), t,t, oro, status, check_units='m', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_hpress), t,t, hpress, status, check_units='Pa', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! default values or actual data?
    if ( args(i_temper) == 'None' ) then
      ! storage:
      allocate( temper(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3)+1:ubo(3)), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! default value:
      temper = 293.0 ! K
    else
      ! pointer to data:
      call self%GetValidPtr( args(i_temper), t,t, temper, status, check_units='K', &
                               check_lbo=lbo+(/0,0,1/), check_ubo=ubo )
      IF_NOTOK_RETURN(status=1)
    end if

    ! default values or actual data?
    if ( args(i_shumid) == 'None' ) then
      ! storage:
      allocate( shumid(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3)+1:ubo(3)), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! default value:
      shumid = 0.0 ! kg/kg
    else
      ! pointer to data:
      call self%GetValidPtr( args(i_shumid), t,t, shumid, status, check_units='kg/kg', &
                               check_lbo=lbo+(/0,0,1/), check_ubo=ubo )
      IF_NOTOK_RETURN(status=1)
    end if

    ! storage for half-level variable:
    allocate( halt(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! compute height profile ;
        ! zero pressure is at infinite distance, so use minimum value:
        call GeoPotentialHeightB( ubo(3)-lbo(3), & ! number of full levels
                  max( 0.01, hpress(i,j,:) ), & ! half lev pres
                  temper(i,j,:), & ! temperature
                  shumid(i,j,:), & ! rel.humid.
                  oro(i,j,1), & ! orography
                  halt(i,j,:) )
      end do ! i
    end do ! j
    
    ! store, use index space if necessary:
    call iF%PutSample( halt, t, status, lbo=lbo, ubo=ubo, change_levs=change_levs )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( halt, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! default values?
    if ( args(i_temper) == 'None' ) then
      ! clear:
      deallocate( temper, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! default values?
    if ( args(i_shumid) == 'None' ) then
      ! clear:
      deallocate( shumid, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine Variables_Calc_HalflevelAltitudes


  ! ***

  subroutine Variables_Calc_HalflevelAltitudes_GPH( self, varname, args, t, status )
  
    use GO  , only : TDate
    use GO  , only : T_Instant_Field
    
    use Binas, only : grav  
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_HalflevelAltitudes_GPH'
    
    ! argument names:
    
    integer, parameter        ::  i_gphbase= 1
    integer, parameter        ::  i_gphpert =2
    !integer, parameter        ::  i_hpress = 4
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  gphbase(:,:,:)
    real, pointer                       ::  gphpert(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  halt(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! get pointer to half-level pressure to obtain number of levels if necessary:
    call self%GetValidPtr( args(i_gphbase), t,t, gphbase, status, &
                              check_units='m2/s2', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! get pointer to half-level pressure to obtain number of levels if necessary:
    call self%GetValidPtr( args(i_gphpert), t,t, gphpert, status, &
                              check_units='m2/s2', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, check_units='m' )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target variable:
    allocate( halt(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! compute height profile ;
    ! zero pressure is at infinite distance, so use minimum value:
    halt = ( gphbase + gphpert )/grav
    
    !! testing ...
    !print *, 'halt_GPH, gphbase, gphper, halt', gphbase(10,10,:), gphpert(10,10,:), halt(10,10,:)
    
    ! store, use index space if necessary:
    call iF%PutSample( halt, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( halt, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Variables_Calc_HalflevelAltitudes_GPH
  

  ! ***


  !
  ! calculate temperature from potential temperature using reference temperature, reference pressure and kappa
  ! T=(potentialtemperature+tref)*(p/p0)**kappa
  !
  
  subroutine T_from_potentialtemperature( self, varname, args, t, status )
  
    use GO  , only : TDate
    use GO  , only : T_Instant_Field
    
    use Binas, only : p0, kappa
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'T_from_potentialtemperature/'
    
    ! argument names:    
    integer, parameter        ::  i_tpot= 1
    integer, parameter        ::  i_pbase =2
    integer, parameter        ::  i_ppert = 3
    
    ! reference temperature for WRF simulation;
    ! not to be confused with T00=290 in WRF output ..
    real, parameter           ::  Tbase=300.0

    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  tpot(:,:,:)
    real, pointer                       ::  ppert(:,:,:)
    real, pointer                       ::  pbase(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  temperature(:,:,:)
     
    ! --- begin ----------------------------------
    
    ! get pointer to  obtain number of levels if necessary:
    call self%GetValidPtr( args(i_tpot), t,t, tpot, status, &
                              check_units='K', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! get pointer to half-level pressure to obtain number of levels if necessary:
    call self%GetValidPtr( args(i_pbase), t,t, pbase, status, &
                              check_units='Pa', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
   ! get pointer to half-level pressure to obtain number of levels if necessary:
    call self%GetValidPtr( args(i_ppert), t,t, ppert, status, &
                              check_units='Pa', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, check_units='K' )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target variable:
    allocate( Temperature(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! compute height profile ;
    ! zero pressure is at infinite distance, so use minimum value:
    temperature= (tpot+Tbase)*((pbase+ppert)/p0)**(kappa)
    
    !! testing ...
    !print *, 'tfrompotentialtemperature, temperature tpot', temperature(10,10,:), tpot(10,10,:)
    !print *, 'tfrompotentialtemperature, pbase, ppert', pbase(10,10,:), ppert(10,10,:)
    
    ! store, use index space if necessary:
    call iF%PutSample( Temperature, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( Temperature, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine T_from_potentialtemperature
  
  
  ! ***
  

  subroutine Variables_Calc_HalflevelPressures( self, varname, args, t, status )
  
    use GO  , only : TDate
    use GO  , only : T_Instant_Field
    use JAQL, only : PotentialPressures
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(4)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_HalflevelPressures'
    
    ! argument names:
    integer, parameter        ::  i_psurf  = 1
    integer, parameter        ::  i_temper = 2
    integer, parameter        ::  i_shumid = 3
    integer, parameter        ::  i_halt = 4
   
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  psurf(:,:,:)
    real, pointer                       ::  temper(:,:,:)
    real, pointer                       ::  shumid(:,:,:)
    real, pointer                       ::  halt(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  hpress(:,:,:)
    integer                             ::  i, j
    integer                             ::  nlayer
    real                                ::  dh  
    ! --- begin ----------------------------------
    
    ! get pointer to half-level altitude to obtain number of levels if necessary:
    call self%GetValidPtr( args(i_halt), t,t, halt, status, &
                           lbo=lbo, ubo=ubo, check_units='m')
    IF_NOTOK_RETURN(status=1)

    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='Pa' )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for surface fields:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)
    
    ! argument:
    call self%GetValidPtr( args(i_psurf), t,t, psurf, status, check_units='Pa', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    ! argument:
    call self%GetValidPtr( args(i_temper), t,t, temper, status, check_units='K', &
                             check_lbo=lbo+(/0,0,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! argument:
    call self%GetValidPtr( args(i_shumid), t,t, shumid, status, check_units='kg/kg', &
                             check_lbo=lbo+(/0,0,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for half-level variable:
    allocate( hpress(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! number of layers (number of half levels minus 1 ...)
    nlayer = ubo(3) - lbo(3)
    ! loop over grid cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)       
        ! compute pressure profile:                     
        call PotentialPressures( nlayer, psurf(i,j,1), &
                                  temper(i,j,:), shumid(i,j,:), halt(i,j,:), &
                                  hpress(i,j,:) ) 
      end do ! i
    end do ! j
    
    ! store, use index space if necessary:
    call iF%PutSample( hpress, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( hpress, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Variables_Calc_Halflevelpressures


  ! ***


  !
  ! Airmass from area and half-level pressures:
  !
  !   call Variables_Calc_Area( self, 'airm', (/'area','hp'/), t, status )
  !

  subroutine Variables_Calc_AirMass( self, varname, args, t, status )
  
    use Binas, only : grav
    use GO   , only : TDate
    use GO   , only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_AirMass'
    
    ! argument names:
    integer, parameter        ::  i_area = 1
    integer, parameter        ::  i_hp   = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  area(:,:,:)
    real, pointer                       ::  hp(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  airm(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='kg', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for surface fields:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)
    
    ! arguments:
    call self%GetValidPtr( args(i_area), t,t, area, status, check_units='m2', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_hp  ), t,t, hp  , status, check_units='Pa', &
                               check_lbo=lbo+(/0,0,-1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for full-level variable:
    allocate( airm(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over target levels:
    do k = lbo(3), ubo(3)
      ! loop over grid cells:
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! airmass:         m2              Pa (kg m/s2 / m2)       / m/s2
          airm(i,j,k) = area(i,j,1) * abs( hp(i,j,k) - hp(i,j,k-1) ) / grav  ! kg
        end do ! i
      end do ! j
    end do ! k
      
    ! store:
    call iF%PutSample( airm, t, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( airm, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_AirMass


  ! ***


  !
  ! Air density from air mass and cell volume:
  !
  !   call Variables_Calc_Area( self, 'dens', (/'airm','vol'/), t, status )
  !

  subroutine Variables_Calc_AirDensity( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_AirMass'
    
    ! argument names:
    integer, parameter        ::  i_airm = 1
    integer, parameter        ::  i_vol  = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  airm(:,:,:)
    real, pointer                       ::  vol(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  dens(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='kg/m3', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! arguments:
    call self%GetValidPtr( args(i_airm), t,t, airm, status, check_units='kg', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_vol ), t,t, vol , status, check_units='m3', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for full-level variable:
    allocate( dens(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do k = lbo(3), ubo(3)
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! air density:
          !               kg           m3
          dens(i,j,k) = airm(i,j,k) / vol(i,j,k)  ! kg/m3
        end do ! i
      end do ! j
    end do ! k

    ! store:
    call iF%PutSample( dens, t, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( dens, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_AirDensity


  ! ***


  !
  ! Vector length (wind speed):
  !
  !   call Variables_Calc_VectorLength( self, 'wsp', (/'uv'/), t, status )
  !

  subroutine Variables_Calc_VectorLength( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_VectorLength'
    
    ! argument names:
    integer, parameter        ::  i_uv = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer           ::  iF
    real, pointer                             ::  u(:,:,:)
    real, pointer                             ::  v(:,:,:)
    character(len=64)                         ::  units
    integer                                   ::  lbo(3), ubo(3)
    real, allocatable                         ::  wsp(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! input variable, extract units and shape:
    call self%GetValidPtr( args(i_uv), t,t, u, v, status, &
                               units=units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to target field, should have same units:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units=trim(units) )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for full-level variable:
    allocate( wsp(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! compute:
    wsp = sqrt( u**2 + v**2 )
    
    ! store:
    call iF%PutSample( wsp, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( wsp, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_VectorLength


  ! ***


  !
  ! Meteorological wind direction:
  !    wind from north =   0 degrees
  !              east  =  90
  !              south = 180
  !              west  = 270
  !
  !    From:
  !      http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
  !    Meteorological wind direction:
  !      wdir = atan2(-u,-v) * 180/pi    # values in [-180,180]
  !    and
  !      where wdir < 0 :  wdir += 360   # values in [0,360]
  !
  !   call Variables_Calc_MeteoWindDir( self, 'wdir', (/'uv'/), t, status )
  !

  subroutine Variables_Calc_MeteoWindDir( self, varname, args, t, status )
  
    use Binas, only : deg2rad
    use GO   , only : TDate
    use GO   , only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_MeteoWindDir'
    
    ! argument names:
    integer, parameter        ::  i_uv = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer           ::  iF
    real, pointer                             ::  u(:,:,:)
    real, pointer                             ::  v(:,:,:)
    integer                                   ::  lbo(3), ubo(3)
    real, pointer                             ::  wdir(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! input variable, extract units and shape:
    call self%GetValidPtr( args(i_uv), t,t, u, v, status, &
                              lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to target field, check units:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='degrees' )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( wdir, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! compute:
    wdir = atan2( -u, -v ) / deg2rad   ! [-180,180]
    ! normalize:
    where ( wdir < 0.0 )
      wdir = wdir + 360.0   ! [0,360]
    end where
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_MeteoWindDir


  ! ***


  !
  ! Compute u-flux from (u,v) windvector and normal vector at mid of u-edge,
  ! and half-level altitudes at corners that define the area of the edge:
  !
  !   call Variables_Calc_VolumeFlux_uv2u( self, 'uflux', &
  !                              (/'uv_u_met','udy','unormal','halt_crnr_met','halt_crnr'/), t, status )
  !

  subroutine Variables_Calc_VolumeFlux_uv2u( self, varname, args, t, status )
  
    use GO     , only : TDate
    use GO     , only : T_Instant_Field
    use C3PO   , only : T_Grid_Ugg
    use Num    , only : IntervalQuad_Const
    use LE_Grid, only : ugg, EDGE_LEFT, EDGE_RIGHT
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(5)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_VolumeFlux_uv2u'
    
    ! argument names:
    integer, parameter        ::  i_uv_u_met      = 1
    integer, parameter        ::  i_udy           = 2
    integer, parameter        ::  i_unormal       = 3
    integer, parameter        ::  i_halt_crnr_met = 4
    integer, parameter        ::  i_halt_crnr     = 5
    
    ! --- local ----------------------------------
    
    class(T_Variable), pointer          ::  V
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  uv_u_met_lon(:,:,:)
    real, pointer                       ::  uv_u_met_lat(:,:,:)
    real, pointer                       ::  halt_crnr_met(:,:,:)
    real, pointer                       ::  halt_crnr(:,:,:)
    real, pointer                       ::  unormal_lon(:,:,:)
    real, pointer                       ::  unormal_lat(:,:,:)
    real, pointer                       ::  udy(:,:,:)
    character(len=64)                   ::  gridtype
    type(T_Grid_Ugg), pointer           ::  gri
    real                                ::  dy
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_met(3), ubo_met(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  halt_met(:)
    real, allocatable                   ::  u_met(:)
    real, allocatable                   ::  uflux(:,:,:)
    integer                             ::  i, j, k
    integer                             ::  ilast
    real                                ::  a, b
    
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='m3/s', lbo=lbo, ubo=ubo, &
                              gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! grid definition of cells:
    select case ( trim(gridtype) )
      case ( 'u-edge' )
        gri => ugg
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! get pointers to wind vector components,
    ! only used to obtain vertical index space:
    call self%GetValidPtr( args(i_uv_u_met), t,t, uv_u_met_lon, uv_u_met_lat, status )
    IF_NOTOK_RETURN(status=1)

    ! index space for meteo on u-edge and original levels:
    lbo_met = (/ lbo(1), lbo(2), lbound(uv_u_met_lon,3) /)
    ubo_met = (/ ubo(1), lbo(2), ubound(uv_u_met_lon,3) /)
    
    ! index space for surface field on u-edge:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), lbo(2), 1 /)
    
    ! get pointers to wind vector components, check bounds:
    call self%GetValidPtr( args(i_uv_u_met   ), t,t, uv_u_met_lon, uv_u_met_lat, status, &
                            check_units='m/s', check_lbo=lbo_met, check_ubo=ubo_met )
    IF_NOTOK_RETURN(status=1)
    ! get pointers to edge lengths:
    call self%GetValidPtr( args(i_udy), t,t, udy, status, &
                             check_units='m', check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    ! get pointers to normal vector components:
    call self%GetValidPtr( args(i_unormal), t,t, unormal_lon, unormal_lat, status, &
                             check_units='1', check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to meteo half-level altitude at corners (extra point in 2nd and 3rd dimension):
    call self%GetValidPtr( args(i_halt_crnr_met), t,t, halt_crnr_met, status, &
                             check_units='m', check_lbo=lbo_met+(/0,-1,-1/), check_ubo=ubo_met )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to target half-level altitude at corners (extra point in 2nd and 3rd dimension):
    call self%GetValidPtr( args(i_halt_crnr    ), t,t, halt_crnr    , status, &
                             check_units='m', check_lbo=lbo+(/0,-1,-1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for half level atitudes at mid of u-edge:
    allocate( halt_met(lbo_met(3)-1:ubo_met(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! storage for wind vector project on normal:
    allocate( u_met(lbo_met(3):ubo_met(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! storage for target array:
    allocate( uflux(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over u-edges:
    do j = lbo(2), ubo(2)      ! 1:ny
      do i = lbo(1), ubo(1)    ! 0:nx

        ! half-level altitudes at mid of u-edge, compute as average of corners:
        halt_met = 0.5*( halt_crnr_met(i,j-1,:) + halt_crnr_met(i,j,:) )
        
        !
        ! Wind vector v=(vx,vy) projected on normal n=(nx,ny) :
        ! Angle between the two vectors, use that length of n is 1:
        !   cos(alfa) = dot(v,n)/len(v)/len(n) = dot(v,n)/len(v)
        ! Length of projection of v on n is:
        !   len(v on n) = len(v) cos(alfa) = len(v) dot(v,n)/len(v) = dot(v,n)
        !
        u_met = uv_u_met_lon(i,j,:) * unormal_lon(i,j,1) + &
                uv_u_met_lat(i,j,:) * unormal_lat(i,j,1)  ! m/s
                
        ! loop over model levels:
        ilast = 1
        do k = lbo(3), ubo(3)
          ! integration interval is between target half-level altiudes at mid of u-edge;
          ! compute as average of corners:
          a = 0.5*( halt_crnr(i,j-1,k-1) + halt_crnr(i,j,k-1) )    ! m
          b = 0.5*( halt_crnr(i,j-1,k  ) + halt_crnr(i,j,k  ) )    ! m
          ! sometimes model is slightly lower due to interpolations;
          ! reset value to bottom of meteo half-level altitudes:
          a = max( halt_met(lbo_met(3)-1), a )
          ! vertcal integral of projected wind vecor over model levels,
          ! keep u constant in meteo layer:
          !                          m        m/s   m  m      m2/s
          call IntervalQuad_Const( halt_met, u_met, a, b, uflux(i,j,k), ilast, status )
          IF_NOTOK_RETURN(status=1)
        end do ! k

        ! edge length:
        dy = udy(i,j,1)
        ! multiply with edge length:
        !                   m2/s      m
        uflux(i,j,:) = uflux(i,j,:) * dy   ! m3/s

      end do ! i
    end do ! j
    
    ! store:
    call iF%PutSample( uflux, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( halt_met, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( u_met, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( uflux, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_VolumeFlux_uv2u


  ! ***


  !
  ! Compute v-flux from (u,v) windvector and normal vector at mid of v-edge,
  ! and half-level altitudes at corners that define the area of the edge:
  !
  !   call Variables_Calc_VolumeFlux_uv2v( self, 'vflux', &
  !                              (/'uv_v_met','vdx','vnormal','halt_crnr_met','halt_crnr'/), t, status )
  !

  subroutine Variables_Calc_VolumeFlux_uv2v( self, varname, args, t, status )
  
    use GO     , only : TDate
    use GO     , only : T_Instant_Field
    use C3PO   , only : T_Grid_Ugg
    use Num    , only : IntervalQuad_Const
    use LE_Grid, only : ugg, EDGE_LOWER, EDGE_UPPER
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(5)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_VolumeFlux_uv2v'
    
    ! argument names:
    integer, parameter        ::  i_uv_v_met      = 1
    integer, parameter        ::  i_vdx           = 2
    integer, parameter        ::  i_vnormal       = 3
    integer, parameter        ::  i_halt_crnr_met = 4
    integer, parameter        ::  i_halt_crnr     = 5
    
    ! --- local ----------------------------------
    
    class(T_Variable), pointer          ::  V
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  uv_v_met_lon(:,:,:)
    real, pointer                       ::  uv_v_met_lat(:,:,:)
    real, pointer                       ::  halt_crnr_met(:,:,:)
    real, pointer                       ::  halt_crnr(:,:,:)
    real, pointer                       ::  vnormal_lon(:,:,:)
    real, pointer                       ::  vnormal_lat(:,:,:)
    real, pointer                       ::  vdx(:,:,:)
    character(len=64)                   ::  gridtype
    type(T_Grid_Ugg), pointer           ::  gri
    real                                ::  dx
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_met(3), ubo_met(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  halt_met(:)
    real, allocatable                   ::  v_met(:)
    real, allocatable                   ::  vflux(:,:,:)
    integer                             ::  i, j, k
    integer                             ::  ilast
    real                                ::  a, b
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='m3/s', lbo=lbo, ubo=ubo, &
                              gridtype=gridtype )
    IF_NOTOK_RETURN(status=1)
    
    ! grid definition of cells:
    select case ( trim(gridtype) )
      case ( 'v-edge' )
        gri => ugg
      case default
        write (gol,'("unsupported grid type: ",a)') trim(gridtype); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! get pointers to wind vector components,
    ! only used to obtain vertical index space:
    call self%GetValidPtr( args(i_uv_v_met), t,t, uv_v_met_lon, uv_v_met_lat, status )
    IF_NOTOK_RETURN(status=1)

    ! index space for meteo on u-edge and original levels:
    lbo_met = (/ lbo(1), lbo(2), lbound(uv_v_met_lon,3) /)
    ubo_met = (/ ubo(1), lbo(2), ubound(uv_v_met_lon,3) /)
    
    ! index space for surface field on u-edge:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), lbo(2), 1 /)
    
    ! get pointers to wind vector components, check bounds:
    call self%GetValidPtr( args(i_uv_v_met), t,t, uv_v_met_lon, uv_v_met_lat, status, &
                            check_units='m/s', check_lbo=lbo_met, check_ubo=ubo_met )
    IF_NOTOK_RETURN(status=1)
    ! get pointers to edge lengths:
    call self%GetValidPtr( args(i_vdx), t,t, vdx, status, &
                             check_units='m', check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    ! get pointers to normal vector components:
    call self%GetValidPtr( args(i_vnormal), t,t, vnormal_lon, vnormal_lat, status, &
                             check_units='1', check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to meteo half-level altitude at corners (extra point in 1st and 3rd dimension):
    call self%GetValidPtr( args(i_halt_crnr_met), t,t, halt_crnr_met, status, &
                             check_units='m', check_lbo=lbo_met+(/-1,0,-1/), check_ubo=ubo_met )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to target half-level altitude at corners (extra point in 1st and 3rd dimension):
    call self%GetValidPtr( args(i_halt_crnr), t,t, halt_crnr, status, &
                             check_units='m', check_lbo=lbo+(/-1,0,-1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for half level atitudes at mid of v-edge:
    allocate( halt_met(lbo_met(3)-1:ubo_met(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! storage for wind vector project on normal:
    allocate( v_met(lbo_met(3):ubo_met(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! storage for target array:
    allocate( vflux(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over v-edges:
    do j = lbo(2), ubo(2)      ! 0:ny
      do i = lbo(1), ubo(1)    ! 1:nx

        ! half-level altitudes at mid of v-edge, compute as average of corners:
        halt_met = 0.5*( halt_crnr_met(i-1,j,:) + halt_crnr_met(i,j,:) )
        
        !
        ! Wind vector v=(vx,vy) projected on normal n=(nx,ny) :
        ! Angle between the two vectors, use that length of n is 1:
        !   cos(alfa) = dot(v,n)/len(v)/len(n) = dot(v,n)/len(v)
        ! Length of projection of v on n is:
        !   len(v on n) = len(v) cos(alfa) = len(v) dot(v,n)/len(v) = dot(v,n)
        !
        v_met = uv_v_met_lon(i,j,:) * vnormal_lon(i,j,1) + &
                uv_v_met_lat(i,j,:) * vnormal_lat(i,j,1)  ! m/s
        
        ! loop over model levels:
        ilast = 1
        do k = lbo(3), ubo(3)
          ! integration interval is between target half-level altiudes at mid of u-edge;
          ! compute as average of corners:
          a = 0.5*( halt_crnr(i-1,j,k-1) + halt_crnr(i,j,k-1) )    ! m
          b = 0.5*( halt_crnr(i-1,j,k  ) + halt_crnr(i,j,k  ) )    ! m
          ! sometimes model is slightly lower due to interpolations;
          ! reset value to bottom of meteo half-level altitudes:
          a = max( halt_met(lbo_met(3)-1), a )
          ! vertcal integral of projected wind vecor over model levels,
          ! keep v constant in meteo layer:
          !                          m        m/s   m  m      m2/s
          call IntervalQuad_Const( halt_met, v_met, a, b, vflux(i,j,k), ilast, status )
          IF_NOTOK_RETURN(status=1)
        end do ! k

        ! edge length:
        dx = vdx(i,j,1)
        ! multiply with edge length:
        !                   m2/s      m
        vflux(i,j,:) = vflux(i,j,:) * dx   ! m3/s

      end do ! i
    end do ! j
    
    ! store:
    call iF%PutSample( vflux, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( halt_met, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( v_met, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( vflux, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_VolumeFlux_uv2v


  ! ***


  !
  ! Volume flux through half-levels from uv-flux and volume change:
  !
  !   call Variables_Calc_VolumeFlux_w( self, 'wflux', (/'uflux','vflux','pvol','nvol'/), t1, t2, status )
  !

  subroutine Variables_Calc_VolumeFlux_w( self, varname, args, t, t1, t2, status )
  
    use GO, only : TDate, operator(-), rTotal
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(4)
    type(TDate), intent(in)                       ::  t, t1, t2
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_VolumeFlux_w'
    
    ! argument names:
    integer, parameter        ::  i_uflux = 1
    integer, parameter        ::  i_vflux = 2
    integer, parameter        ::  i_pvol  = 3
    integer, parameter        ::  i_nvol  = 4
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  uflux(:,:,:)
    real, pointer                       ::  vflux(:,:,:)
    real, pointer                       ::  pvol(:,:,:)
    real, pointer                       ::  nvol(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  dwflux(:,:,:)
    real, allocatable                   ::  wflux(:,:,:)
    integer                             ::  i, j, k
    real                                ::  dt_sec
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                               check_units='m3/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! arguments:
    call self%GetValidPtr( args(i_uflux), t ,t , uflux, status, check_units='m3/s', &
                               check_lbo=lbo+(/-1,0,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_vflux), t ,t , vflux, status, check_units='m3/s', &
                               check_lbo=lbo+(/0,-1,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_pvol ), t1,t1, pvol , status, check_units='m3', &
                               check_lbo=lbo+(/0,0,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_nvol ), t2,t2, nvol , status, check_units='m3', &
                               check_lbo=lbo+(/0,0,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! step in seconds:
    dt_sec = rTotal( t2 - t1, 'sec' )

    ! storage for net wflux:
    allocate( dwflux(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3)+1:ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do k = lbo(3)+1, ubo(3)
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! fill with first contribution: negative volume change:
          !          (      m3     )/  s
          dwflux(i,j,k) = -( nvol(i,j,k) - pvol(i,j,k) ) / dt_sec
          ! add net horizontal flux:   leaving through    enter through
          !                              'upper' edge  +  'lower' edge
          !                                   m3/s             m3/s
          dwflux(i,j,k) = dwflux(i,j,k) - uflux(i,j,k) + uflux(i-1,j  ,k)
          dwflux(i,j,k) = dwflux(i,j,k) - vflux(i,j,k) + vflux(i  ,j-1,k)
        end do
      end do
    end do
    
    ! storage on half levels:
    allocate( wflux(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! no flux through surface:
    wflux(:,:,lbo(3)) = 0.0
    ! loop over upper half levels:
    do k = lbo(3)+1, ubo(3)
      ! from:
      !  delta_flux = upper_flux - lower_flux
      ! use:
      !  upper_flux = lower_flux + delta_flux
      wflux(:,:,k) = wflux(:,:,k-1) + dwflux(:,:,k)
    end do
    
    ! store:
    call iF%PutSample( wflux, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( dwflux, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( wflux, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_VolumeFlux_w


  !
  ! Volume flux through half-levels from uv-flux (asuming constant volume):
  !
  !   call Variables_Calc_VolumeFlux_w0( self, 'wflux', (/'uflux','vflux'/), t1, t2, status )
  !

  subroutine Variables_Calc_VolumeFlux_w0( self, varname, args, t, t1, t2, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t, t1, t2
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_VolumeFlux_w0'
    
    ! argument names:
    integer, parameter        ::  i_uflux = 1
    integer, parameter        ::  i_vflux = 2
    
    ! conversion factors:
    real, parameter  ::  km3_per_m3 = 1.0e-9   ! km3/m3
    real, parameter  ::  s_per_min  = 60.0     ! s/min
    real, parameter  ::  m3s_to_km3min = km3_per_m3 * s_per_min ! (km3/min)/(m3/s)
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  uflux(:,:,:)
    real, pointer                       ::  vflux(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  dwflux(:,:,:)
    real, allocatable                   ::  wflux(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m3/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! arguments:
    call self%GetValidPtr( args(i_uflux), t ,t , uflux, status, check_units='m3/s', &
                               check_lbo=lbo+(/-1,0,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_vflux), t ,t , vflux, status, check_units='m3/s', &
                               check_lbo=lbo+(/0,-1,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for net wflux:
    allocate( dwflux(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3)+1:ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! initialize for no volume change:
    dwflux = 0.0  ! m3/s
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! add net horizontal flux:   leaving through    enter through
        !                              'upper' edge  +  'lower' edge
        !                                   m3/s             m3/s
        dwflux(i,j,:) = dwflux(i,j,:) - uflux(i,j,:) + uflux(i-1,j  ,:)
        dwflux(i,j,:) = dwflux(i,j,:) - vflux(i,j,:) + vflux(i  ,j-1,:)
      end do
    end do
    
    ! storage on half levels:
    allocate( wflux(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! no flux through surface:
    wflux(:,:,lbo(3)) = 0.0
    ! loop over upper half levels:
    do k = lbo(3)+1, ubo(3)
      ! from:
      !  delta_flux = upper_flux - lower_flux
      ! use:
      !  upper_flux = lower_flux + delta_flux
      wflux(:,:,k) = wflux(:,:,k-1) + dwflux(:,:,k)
    end do
    
    ! store:
    call iF%PutSample( wflux, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( dwflux, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( wflux, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_VolumeFlux_w0


  ! ***
  
  !
  ! Vertical interpolation.
  ! Input:   x_in, y_in, x
  ! Output:  y
  !

  subroutine Variables_Calc_LayerInterpol( self, varname, args, t, status )
  
    use GO , only : TDate
    use GO , only : T_Instant_Field
    use Num, only : Interp_Lin, EXTRAPOL_LINEAR
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_LayerInterpol'
    
    ! argument names:
    integer, parameter        ::  i_x_in = 1
    integer, parameter        ::  i_y_in = 2
    integer, parameter        ::  i_x    = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  x_in(:,:,:)
    real, pointer                       ::  y_in(:,:,:)
    real, pointer                       ::  x(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_in(3), ubo_in(3)
    real, allocatable                   ::  y(:,:,:)
    character(len=64)                   ::  x_units
    character(len=64)                   ::  y_units
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              units=y_units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! set pointer to input to obtain vertical index space:
    call self%GetValidPtr( args(i_x_in), t,t, x_in, status )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for input:
    lbo_in = (/ lbo(1), lbo(2), lbound(x_in,3) /)
    ubo_in = (/ ubo(1), ubo(2), ubound(x_in,3) /)
    
    ! input:
    call self%GetValidPtr( args(i_x_in), t,t, x_in, status, units=x_units, &
                               check_lbo=lbo_in, check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_y_in), t,t, y_in, status, check_units=y_units, &
                               check_lbo=lbo_in, check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_x   ), t,t, x   , status, check_units=x_units, &
                               check_lbo=lbo, check_ubo=ubo  )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( y(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over grid cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! linear interpolation:
        call Interp_Lin( x_in(i,j,:), y_in(i,j,:), x(i,j,lbo(3):ubo(3)), y(i,j,:), status, &
                           extrapol=EXTRAPOL_LINEAR )
        IF_NOTOK_RETURN(status=1)
      end do  ! i
    end do  ! j
    
    ! store:
    call iF%PutSample( y, t, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( y, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_LayerInterpol


  ! ***
  
  !
  ! Vertical interpolation, special case for interpolation
  ! At the interface between mixing-layer and reservoir layer
  !  pick the top value, not interpolation
  ! Input:   halt_met, kz_ifs, halt, blh
  ! Output:  kz
  !

  subroutine Variables_Calc_LayerInterpol_MixLayer_Kz( self, varname, args, t, status )
  
    use GO , only : TDate
    use GO , only : T_Instant_Field
    use Num, only : Interp_Lin, EXTRAPOL_LINEAR
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(4)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_LayerInterpol_MixLayer_Kz'
    
    ! argument names:
    integer, parameter        ::  i_halt_met = 1
    integer, parameter        ::  i_kz_ifs   = 2
    integer, parameter        ::  i_halt     = 3
    integer, parameter        ::  i_blh      = 4
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  halt_met(:,:,:)
    real, pointer                       ::  kz_ifs(:,:,:)
    real, pointer                       ::  halt(:,:,:)
    real, pointer                       ::  blh(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_in(3), ubo_in(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  kz(:,:,:)
    character(len=64)                   ::  x_units
    character(len=64)                   ::  y_units
    integer                             ::  i, j, l
    integer                             ::  ilev_blh
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              units=y_units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! set pointer to input to obtain vertical index space:
    call self%GetValidPtr( args(i_halt_met), t,t, halt_met, status )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for input:
    lbo_in = (/ lbo(1), lbo(2), lbound(halt_met,3) /)
    ubo_in = (/ ubo(1), ubo(2), ubound(halt_met,3) /)
    
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)
    ! input:
    call self%GetValidPtr( args(i_halt_met), t,t, halt_met, status, units=x_units, &
                               check_lbo=lbo_in, check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_kz_ifs), t,t, kz_ifs, status, check_units=y_units, &
                               check_lbo=lbo_in, check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_halt   ), t,t, halt   , status, check_units=x_units, &
                               check_lbo=lbo, check_ubo=ubo  )
    IF_NOTOK_RETURN(status=1)
    
    call self%GetValidPtr( args(i_blh    ), t,t, blh    , status, check_units=x_units, &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc  )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( kz(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over grid cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)

        ! linear interpolation:
        call Interp_Lin( halt_met(i,j,:), kz_ifs(i,j,:), halt(i,j,lbo(3):ubo(3)), kz(i,j,:), status, &
                           extrapol=EXTRAPOL_LINEAR )
        IF_NOTOK_RETURN(status=1)

        ! change value at interface mixing/reservoir layer
        ! second layer is defined as mixing layer, so at top of this layer almost no diffusion 
        ilev_blh = 2
        ! Take value from layer above
        kz(i,j,ilev_blh) = kz(i,j,ilev_blh+1)
      end do  ! i
    end do  ! j
    
    ! store:
    call iF%PutSample( kz, t, status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( kz, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_LayerInterpol_MixLayer_Kz


  ! ***


  !
  ! Layer average of 'values' :
  ! - mass weighted using ph (Pa) afterp mass/m2 = p/g
  ! - height weighted using hh (m)
  ! Example:
  !   call Variables_Calc_LayerAverage( self, 'y', (/'hx_in','values_in','hx'/), t, status )
  !

  subroutine Variables_Calc_LayerAverage( self, varname, args, t, status )
  
    use GO              , only : TDate
    use GO              , only : T_Instant_Field
    use Num             , only : IntervalQuad_Const
    use LE_Data_Variable, only : LEN_UNITS
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_LayerAverage'
    
    ! argument names:
    integer, parameter        ::  i_hx_in     = 1
    integer, parameter        ::  i_values_in = 2
    integer, parameter        ::  i_hx        = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  hx_in(:,:,:)
    real, pointer                       ::  values_in(:,:,:)
    real, pointer                       ::  hx(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_in(3), ubo_in(3)
    real, allocatable                   ::  hx_inx(:)
    real, allocatable                   ::  values(:,:,:)
    character(len=LEN_UNITS)            ::  x_units
    character(len=LEN_UNITS)            ::  y_units
    integer                             ::  i, j, k
    integer                             ::  isfc_in, isfc
    real                                ::  sgn
    integer                             ::  ilast
    integer                             ::  ilev
      
    ! --- begin ----------------------------------

    ! set pointer to input to obtain vertical index space:
    call self%GetValidPtr( args(i_values_in), t,t, values_in, status )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_Field=iF, &
                              units=y_units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for input:
    lbo_in = (/ lbo(1), lbo(2), lbound(values_in,3) /)
    ubo_in = (/ ubo(1), ubo(2), ubound(values_in,3) /)

    ! dependencies:
    call self%GetValidPtr( args(i_hx_in    ), t,t, hx_in    , status, units=x_units, &
                               check_lbo=lbo_in+(/0,0,-1/), check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_values_in), t,t, values_in, status, check_units=y_units, &
                               check_lbo=lbo_in, check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_hx       ), t,t, hx       , status, check_units=x_units, &
                               check_lbo=lbo+(/0,0,-1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! surface level of model fields:
    if ( hx(lbo(1),lbo(2),lbo(3)-1) > hx(lbo(1),lbo(2),ubo(3)) ) then
      if ( trim(x_units) == 'm' ) then
        isfc = ubo(3)
      else if ( trim(x_units) == 'Pa' ) then
        isfc = lbo(3)-1
      else
        write( gol, '("Unsupported units to define sfc level: ", a)' ) trim(x_units) ; call GoErr
        TRACEBACK;status=1;return
      end if
    else
      if (trim(x_units) == 'm') then
        isfc = lbo(3)-1
      else if ( trim(x_units) == 'Pa') then
        isfc = ubo(3)
      else
        write( gol, '("Unsupported units to define sfc level: ", a)' ) trim(x_units) ; call GoErr
        TRACEBACK;status=1;return
      end if
    end if
    
    ! surface level of meteo fields:
    if ( hx_in(lbo_in(1),lbo_in(2),lbo_in(3)-1) > hx_in(lbo_in(1),lbo_in(2),ubo_in(3)) ) then
      if ( trim( x_units) == 'm' ) then
        isfc_in = ubo_in(3)
      else if ( trim(x_units) == 'Pa') then
        isfc_in = lbo_in(3)-1
      else
        write( gol, '("Unsupported units to define input sfc level: ", a)' ) trim(x_units) ; call GoErr
        TRACEBACK;status=1;return
      end if
      sgn = -1.0  ! negate levels to have increasing order
    else
      if ( trim( x_units) == 'm' ) then
        isfc_in = lbo_in(3)-1
      else if ( trim(x_units) == 'Pa') then
        isfc_in = ubo_in(3)
      else
        write( gol, '("Unsupported units to define input sfc level: ", a)' ) trim(x_units) ; call GoErr
        TRACEBACK;status=1;return
      end if
      sgn = 1.0  ! keep input levels
    end if
    ! storage for copy with synchronized surface pressure:
    allocate( hx_inx(lbo_in(3)-1:ubo_in(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! storage for field on levels:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! copy:
        hx_inx = hx_in(i,j,:)
        ! reset input surface to target surface if that is lower:
        if ( sgn > 0.0 ) then
          hx_inx(isfc_in) = min( hx_inx(isfc_in), hx(i,j,isfc) )
        else
          hx_inx(isfc_in) = max( hx_inx(isfc_in), hx(i,j,isfc) )
        end if
        ! init index:
        ilast = 1
        ! loop over levels:
        do k = lbo(3), ubo(3)
          ! integral over ax:
          call IntervalQuad_Const( sgn*hx_inx, values_in(i,j,:), &
                             sgn*hx(i,j,k-1), sgn*hx(i,j,k), &
                             values(i,j,k), ilast, status )
          if ( status /= 0 ) then
            write (gol,'("cell (",i0,",",i0,")")') i, j; call goErr
            write (gol,'("swap sign       : ",f7.1)') sgn; call goErr
            write (gol,'("input axes `",a,"` (surface at ",i0,")")') trim(args(i_hx_in)), isfc_in; call goErr
            do ilev = lbound(hx_inx,1), ubound(hx_inx,1)
              write (gol,'(i6,2e15.6)') ilev, sgn*hx_inx(ilev), hx_in(i,j,ilev); call goErr
            end do
            write (gol,'("input values `",a,"`")') trim(args(i_values_in)); call goErr
            write (gol,'("target ax `",a,"` (surface at ",i0,")")') trim(args(i_hx)), isfc; call goErr
            do ilev = lbo(3), ubo(3)
              write (gol,'(i6,2e15.6)') ilev, sgn*hx(i,j,ilev-1), sgn*hx(i,j,ilev); call goErr
            end do
            write (gol,'("target values `",a,"`")') trim(varname); call goErr
            TRACEBACK; status=1; return
          end if
          ! weight:
          values(i,j,k) = values(i,j,k) / ( sgn*hx(i,j,k) - sgn*hx(i,j,k-1) )
        end do ! levels
      end do ! i
    end do ! j
    
    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( hx_inx, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_LayerAverage

  ! *

  subroutine Variables_Calc_LayerAverageUV( self, varname, args, t, status )
  
    use GO , only : TDate
    use GO , only : T_Instant_VectorField
    use Num, only : IntervalQuad_Const
    use LE_Data_Variable, only : LEN_UNITS
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_LayerAverageUV'
    
    ! argument names:
    integer, parameter        ::  i_hx_in = 1
    integer, parameter        ::  i_uv_in = 2
    integer, parameter        ::  i_hx    = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_VectorField), pointer     ::  iVF
    real, pointer                             ::  hx_in(:,:,:)
    real, pointer                             ::  u_in(:,:,:)
    real, pointer                             ::  v_in(:,:,:)
    real, pointer                             ::  hx(:,:,:)
    integer                                   ::  lbo(3), ubo(3)
    integer                                   ::  lbo_in(3), ubo_in(3)
    real, allocatable                         ::  hx_inx(:)
    real, allocatable                         ::  u(:,:,:)
    real, allocatable                         ::  v(:,:,:)
    character(len=LEN_UNITS)                  ::  x_units
    character(len=LEN_UNITS)                  ::  uv_units
    integer                                   ::  i, j, k
    integer                                   ::  isfc_in, isfc
    real                                      ::  sgn
    integer                                   ::  ilast
    integer                                   ::  ilev
      
    ! --- begin ----------------------------------

    ! set pointer to input to obtain vertical index space:
    call self%GetValidPtr( args(i_uv_in), t,t, u_in, v_in, status )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_VectorField=iVF, &
                              units=uv_units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for input:
    lbo_in = (/ lbo(1), lbo(2), lbound(u_in,3) /)
    ubo_in = (/ ubo(1), ubo(2), ubound(u_in,3) /)

    ! dependencies:
    call self%GetValidPtr( args(i_hx_in), t,t, hx_in, status, units=x_units, &
                               check_lbo=lbo_in+(/0,0,-1/), check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_uv_in), t,t, u_in, v_in, status, check_units=uv_units, &
                               check_lbo=lbo_in, check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_hx), t,t, hx, status, check_units=x_units, &
                               check_lbo=lbo+(/0,0,-1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! surface level of model fields:
    if ( hx(lbo(1),lbo(2),lbo(3)-1) > hx(lbo(1),lbo(2),ubo(3)) ) then
      isfc = lbo(3)-1
    else
      isfc = ubo(3)
    end if
    
    ! surface level of meteo fields:
    if ( hx_in(lbo_in(1),lbo_in(2),lbo_in(3)-1) > hx_in(lbo_in(1),lbo_in(2),ubo_in(3)) ) then
      isfc_in = lbo_in(3)-1
      sgn = -1.0  ! negate levels to have increasing order
    else
      isfc_in = ubo_in(3)
      sgn = 1.0  ! keep input levels
    end if
    ! storage for copy with synchronized surface pressure:
    allocate( hx_inx(lbo_in(3)-1:ubo_in(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! storage for field on levels:
    allocate( u(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( v(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! copy:
        hx_inx = hx_in(i,j,:)
        ! reset input surface to target surface if that is lower:
        if ( sgn > 0.0 ) then
          hx_inx(isfc_in) = min( hx_inx(isfc_in), hx(i,j,isfc) )
        else
          hx_inx(isfc_in) = max( hx_inx(isfc_in), hx(i,j,isfc) )
        end if
        ! init index:
        ilast = 1
        ! loop over levels:
        do k = lbo(3), ubo(3)
          ! integral over ax:
          call IntervalQuad_Const( sgn*hx_inx, u_in(i,j,:), &
                             sgn*hx(i,j,k-1), sgn*hx(i,j,k), &
                             u(i,j,k), ilast, status )
          if ( status /= 0 ) then
            write (gol,'("cell (",i0,",",i0,")")') i, j; call goErr
            write (gol,'("swap sign       : ",f7.1)') sgn; call goErr
            write (gol,'("input axes `",a,"` (surface at ",i0,")")') trim(args(i_hx_in)), isfc_in; call goErr
            do ilev = lbound(hx_inx,1), ubound(hx_inx,1)
              write (gol,'(i6,2e15.6)') ilev, sgn*hx_inx(ilev), hx_in(i,j,ilev); call goErr
            end do
            write (gol,'("input values `",a,"`")') trim(args(i_uv_in)); call goErr
            write (gol,'("target ax `",a,"` (surface at ",i0,")")') trim(args(i_hx)), isfc; call goErr
            do ilev = lbo(3), ubo(3)
              write (gol,'(i6,2e15.6)') ilev, sgn*hx(i,j,ilev-1), sgn*hx(i,j,ilev); call goErr
            end do
            write (gol,'("target values `",a,"`")') trim(varname); call goErr
            TRACEBACK; status=1; return
          end if
          ! weight:
          u(i,j,k) = u(i,j,k) / ( sgn*hx(i,j,k) - sgn*hx(i,j,k-1) )
          !
          !~ idem for second field:
          call IntervalQuad_Const( sgn*hx_inx, v_in(i,j,:), &
                             sgn*hx(i,j,k-1), sgn*hx(i,j,k), &
                             v(i,j,k), ilast, status )
          ! weight:
          v(i,j,k) = v(i,j,k) / ( sgn*hx(i,j,k) - sgn*hx(i,j,k-1) )
        end do ! levels
      end do ! i
    end do ! j
    
    ! store:
    call iVF%PutSample( u, v, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( u, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( v, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( hx_inx, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_LayerAverageUV

  ! ***


  !
  ! windspeed
  !

  subroutine Variables_Calc_WindSpeed( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_WindSpeed'
    
    ! argument names:
    integer, parameter        ::  i_u = 1
    integer, parameter        ::  i_v = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  u(:,:,:)
    real, pointer                       ::  v(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_Field=iF, &
                              check_units='m/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_u), t,t, u, status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_v), t,t, v, status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells in target space:
    do k = lbo(3), ubo(3)
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! compute:
          values(i,j,k) = sqrt( u(i,j,k)**2 + v(i,j,k)**2 )
        end do ! i
      end do ! j
    end do ! k
    
    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_WindSpeed


  ! ***


  !
  ! tendency to convert from accumulated fields to temporal average;
  ! for example:
  !     precip  =  d(accum_precip)/dt
  !      m/s             m        /s
  !

  subroutine Variables_Calc_Tendency( self, varname, args, t1,t2, status )
  
    use GO  , only : TDate
    use GO  , only : T_Accumulated_Field_Series, T_Constant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t1,t2
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_Tendency'
    
    ! argument names:
    integer, parameter        ::  i_var = 1
    
    ! --- local ----------------------------------
    
    class(T_Constant_Field), pointer              ::  cF
    class(T_Accumulated_Field_Series), pointer    ::  aFS
    character(len=32)                             ::  units
    integer                                       ::  lbo(3), ubo(3)
    real, allocatable                             ::  ddata_dt(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Constant_Field=cF, &
                              units=units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to accumulated series as specified in arguments ;
    ! check if the units match (sum over time):
    call self%GetVariable( args(i_var), status, Accumulated_Field_Series=aFS, &
                              check_units=trim(units)//' s', &
                              check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( ddata_dt(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! get tendency, check if valid for requested interval:
    call aFS%Get_Tendency( (/t1,t2/), ddata_dt, status )
    IF_NOTOK_RETURN(status=1)
    ! store:
    call cF%PutSample( ddata_dt, (/t1,t2/), status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( ddata_dt, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_Tendency


  ! ***


  !
  ! total rain:
  !  - sum of large-scale precip, convective precip, ..
  !  - only substantial amount (> 0.1 mm/hr ?)
  !

  subroutine Variables_Calc_TotalRain( self, varname, args, t1,t2, status )
  
    use GO  , only : TDate
    use GO  , only : T_Constant_Field
    use Dims, only : substantial_rain__mps
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(:)
    type(TDate), intent(in)                       ::  t1,t2
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_TotalRain'
       
    ! --- local ----------------------------------
    
    class(T_Constant_Field), pointer    ::  cF
    real, pointer                       ::  p(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  precip(:,:,:)
    integer                             ::  ivar
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Constant_Field=cF, &
                              check_units='m/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to target array (allocate if necessary), set target time:
    call cF%GetTargetPtr( precip, (/t1,t2/), status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! init sum:
    precip = 0.0
    ! loop over input variables:
    do ivar = 1, size(args)

      ! pointer to source array:
      call self%GetValidPtr( args(ivar), t1,t2, p, status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
      IF_NOTOK_RETURN(status=1)

      ! loop over cells:
      do k = lbo(3), ubo(3)
        do j = lbo(2), ubo(2)
          do i = lbo(1), ubo(1)
            ! add contribution:
            precip(i,j,k) = precip(i,j,k) + p(i,j,k)
          end do ! i
        end do ! j
      end do ! k

    end do ! input variables

    ! only substantial amount:
    where ( precip <= substantial_rain__mps ) precip = 0.0   ! m/s
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_TotalRain


  ! ***


  !
  ! 3D rain from total (2D) and potential intensity (3D)
  !

  subroutine Variables_Calc_RainIntensity( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_RainIntensity'
    
    ! argument names:
    integer, parameter        ::  i_rain   = 1
    integer, parameter        ::  i_praini = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  rain  (:,:,:)   ! (nlon,nlat,1)
    real, pointer                       ::  praini(:,:,:)   ! (nlon,nlat,nlev)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  values(:,:,:)   ! (nlon,nlat,nlev)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_Field=iF, &
                              check_units='m/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for surface field:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)
    
    ! dependencies:
    call self%GetValidPtr( args(i_rain  ), t,t, rain  , status, check_units='m/s', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_praini), t,t, praini, status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do k = lbo(3), ubo(3)
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! compute:
          values(i,j,k) = rain(i,j,1) * praini(i,j,k)
        end do ! i
      end do ! j
    end do ! k
    
    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_RainIntensity


  ! ***


  !
  ! Relative humidity ('RH'):
  !
  !   call Variables_Calc_RelativeHumidityTD( self, 'RH', (/'T','D'/), t, status )
  !

  subroutine Variables_Calc_RelativeHumidityTD( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_RelativeHumidityTD'
    
    ! argument names:
    integer, parameter        ::  i_temper   = 1
    integer, parameter        ::  i_dewpoint = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  temper(:,:,:)
    real, pointer                       ::  dewpoint(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_Field=iF, &
                              check_units='%', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_temper  ), t,t, temper  , status, check_units='K', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_dewpoint), t,t, dewpoint, status, check_units='K', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do k = lbo(3), ubo(3)
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! compute relative humidity from temperature and dewpoint:
          !                       K              K                K              K
          values(i,j,k) = ES( dewpoint(i,j,k) - 273.15 ) / ES( temper(i,j,k) - 273.15 ) * 100.0  ! %
        end do ! i
      end do ! j
    end do ! k
    
    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_RelativeHumidityTD

  
  ! ***


  subroutine  Variables_Calc_SH_from_wvmixingratio( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field_series
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'Variables_calc_SH_from_wvmixingratio'
    
    ! argument names:
    integer, parameter        ::  i_wvmix   = 1
    
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field_series), pointer     ::  iFs
    real, pointer                       ::  wvmix(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_Field_series=iFs, &
                              check_units='kg/kg' )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_wvmix ), t,t, wvmix , status, check_units='kg/kg', &
                               lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

   
    !from water vapour mixing ratio in kg water/kg air to specific humidity in kg/kg
    values=wvmix/(1+wvmix)
  
    ! store:
    call iFs%PutSample( values, t, status , lbo=lbo, ubo=ubo)
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine  Variables_calc_SH_from_wvmixingratio 
  ! *
  ! ***


  !
  ! *

  real elemental function ES( t )

    ! --- in/out ---------------------------

    real, intent(in)   ::  t  ! degrees Celcius

    ! --- begin ----------------------------

    ES = 6.1365 * exp( 17.502 * t / (240.97 + t) )

  end function ES


  ! ***


  !
  ! Relative humidity ('RH'):
  !
  !   call Variables_Calc_RelativeHumidityPTQ( self, 'RH', (/'p','T','Q'/), t, status )
  !

  subroutine Variables_Calc_RelativeHumidityPTQ( self, varname, args, t, status )
  
    use GO  , only : TDate
    use GO  , only : T_Instant_Field
    use JAQL, only : RelativeHumidity_from_SpecificHumidity
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_RelativeHumidityPTQ'
    
    ! argument names:
    integer, parameter        ::  i_p = 1
    integer, parameter        ::  i_T = 2
    integer, parameter        ::  i_Q = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  pressure(:,:,:)
    real, pointer                       ::  temper(:,:,:)
    real, pointer                       ::  shumid(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_Field=iF, &
                              check_units='%', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_p), t,t, pressure, status, check_units='Pa', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_T), t,t, temper  , status, check_units='K', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_Q), t,t, shumid  , status, check_units='kg/kg', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do k = lbo(3), ubo(3)
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! compute relative humidity from pressure, temperature, and specific humidity:
          values(i,j,k) = RelativeHumidity_from_SpecificHumidity( shumid(i,j,k), temper(i,j,k), pressure(i,j,k) ) ! %
        end do ! i
      end do ! j
    end do ! k
    
    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_RelativeHumidityPTQ


  ! ***


  subroutine Variables_Calc_TotalCloudCoverage( self, varname, args, t, status )
  
    use GO  , only : TDate
    use GO  , only : T_Instant_Field
    use jaql_cloud, only: CloudCoverOverhead
    
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_TotalCloudCoverage'
    
    ! argument names:
    integer, parameter        ::  i_cc    = 1
    
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  cc(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  tcc(:,:,:)
    real, allocatable                   ::  occ(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
 

    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='1' )
    IF_NOTOK_RETURN(status=1)

    ! get pointer to cc to obtain number of levels if necessary:
    call self%GetValidPtr( args(i_cc), t,t, cc, status, &
                lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for surface fields:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)

    ! storage for tcc
    allocate( tcc(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( occ(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

     ! loop:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! compute height profile ;
        ! zero pressure is at infinite distance, so use minimum value:
        call CloudCoverOverhead( ubo(3)-lbo(3)+1, & ! number of full levels
                  cc(i,j,:), & !cloud fractions on layers
                  occ(i,j,:), & ! overhead cloud cover
                  status)
        ! copy total cloud cover from overhead at surface:
        tcc(i,j,1) = occ(i,j,lbo(3))
      end do ! i
    end do ! j
    
    !! testing ...
    !print *, 'totalcloudcover, tcc, occ', tcc(10,10,1), occ(10,10,:)
    
    ! store, use index space if necessary:
    call iF%PutSample( tcc, t, status, lbo=lbo_sfc, ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( tcc, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( occ, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Variables_Calc_TotalCloudCoverage


  ! ***


  !
  ! Cloud profiles
  !   output: icc, bcc, occ, clwc, pri
  !   input : hp_met, tcc_met, cc_met, clwc_met, hp
  !

  subroutine Variables_Calc_CloudProfiles( self, varnames, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use LE_Data_Calc, only : Get_Cloud_Profiles
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varnames(6)
    character(len=*), intent(in)                  ::  args(6)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_CloudProfiles'
    
    ! output names:
    integer, parameter        ::  o_icc  = 1
    integer, parameter        ::  o_bcc  = 2
    integer, parameter        ::  o_occ  = 3
    integer, parameter        ::  o_clwc = 4
    integer, parameter        ::  o_ciwc = 5
    integer, parameter        ::  o_pri  = 6
    ! argument names:
    integer, parameter        ::  i_hp_met   = 1
    integer, parameter        ::  i_tcc_met  = 2
    integer, parameter        ::  i_cc_met   = 3
    integer, parameter        ::  i_clwc_met = 4
    integer, parameter        ::  i_ciwc_met = 5
    integer, parameter        ::  i_hp       = 6
    
    ! --- local ----------------------------------
    
    logical                             ::  enabled
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::    hp_met(:,:,:)
    real, pointer                       ::   tcc_met(:,:,:)
    real, pointer                       ::    cc_met(:,:,:)
    real, pointer                       ::  clwc_met(:,:,:)
    real, pointer                       ::  ciwc_met(:,:,:)
    real, pointer                       ::        hp(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    integer                             ::  lbo_met(3), ubo_met(3)
    real, allocatable                   ::   icc(:,:,:)
    real, allocatable                   ::   bcc(:,:,:)
    real, allocatable                   ::   occ(:,:,:)
    real, allocatable                   ::  clwc(:,:,:)
    real, allocatable                   ::  ciwc(:,:,:)
    real, allocatable                   ::   pri(:,:,:)
    integer                             ::  i, j
    integer                             ::  nlev_met
    integer                             ::  nlev
    real, allocatable                   ::   cc_met1(:)
    real, allocatable                   :: clwc_met1(:)
    real, allocatable                   :: ciwc_met1(:)
      
    ! --- begin ----------------------------------

    ! pointer to target field given the name;
    ! not sure what happends if this is not enabled,
    ! but we need the index space for checking ..
    call self%GetVariable( varnames(o_icc ), status, Instant_Field=iF, &
                              check_units='1', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! set pointer to input field to obtain vertical index space:
    call self%GetValidPtr( args(i_cc_met  ), t,t, cc_met, status )
    IF_NOTOK_RETURN(status=1)
    
    ! index space for input:
    lbo_met = (/ lbo(1), lbo(2), lbound(cc_met,3) /)
    ubo_met = (/ ubo(1), ubo(2), ubound(cc_met,3) /)
    ! index space for surface input:
    lbo_sfc = (/ lbo(1), lbo(2), 1 /)
    ubo_sfc = (/ ubo(1), ubo(2), 1 /)
    
    ! dependencies:
    call self%GetValidPtr( args(i_hp_met  ), t,t, hp_met  , status, check_units='Pa', &
                               check_lbo=lbo_met+(/0,0,-1/), check_ubo=ubo_met    )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_tcc_met ), t,t, tcc_met , status, check_units='1', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc     )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_cc_met  ), t,t, cc_met  , status, check_units='1', &
                               check_lbo=lbo_met, check_ubo=ubo_met     )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_clwc_met), t,t, clwc_met, status, check_units='kg/kg', &
                               check_lbo=lbo_met, check_ubo=ubo_met )
    IF_NOTOK_RETURN(status=1)
    if ( args(i_ciwc_met) == 'None' ) then
      ! not available on input, allocate dummy array with zeros:
      allocate( ciwc_met(size(clwc_met,1),size(clwc_met,2),size(clwc_met,3)), stat=status, source=0.0 )
      IF_NOTOK_RETURN(status=1)
    else
      call self%GetValidPtr( args(i_ciwc_met), t,t, ciwc_met, status, check_units='kg/kg', &
                                 check_lbo=lbo_met, check_ubo=ubo_met )
      IF_NOTOK_RETURN(status=1)
    end if
    call self%GetValidPtr( args(i_hp      ), t,t, hp      , status, check_units='Pa', &
                               check_lbo=lbo+(/0,0,-1/), check_ubo=ubo    )
    IF_NOTOK_RETURN(status=1)
    
    ! number of input layers:
    nlev_met = ubo_met(3) - lbo_met(3) + 1
    ! number of model layers:
    nlev = ubo(3) - lbo(3) + 1
    
    ! storage for input columns, will be truncated to avoid rounding errors:
    allocate(   cc_met1(nlev_met), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( clwc_met1(nlev_met), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( ciwc_met1(nlev_met), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! storage for target fields (full levels):
    allocate(  icc(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate(  bcc(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate(  occ(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( clwc(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( ciwc(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate(  pri(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! extract:
          cc_met1 =   cc_met(i,j,:)
        clwc_met1 = clwc_met(i,j,:)
        ciwc_met1 = ciwc_met(i,j,:)
        ! truncate very small cloud fractions (1% or less):
        where ( cc_met1 < 0.01 )
            cc_met1 = 0.0
          clwc_met1 = 0.0
          ciwc_met1 = 0.0
        end where
        ! convert cloud profile to LE layers:
        call Get_Cloud_Profiles( nlev_met, hp_met(i,j,:), &
                                  tcc_met(i,j,1), cc_met1, clwc_met1, ciwc_met1, &
                                 nlev, hp(i,j,:), &
                                  icc(i,j,:), bcc(i,j,:), occ(i,j,:), clwc(i,j,:), ciwc(i,j,:), pri(i,j,:), &
                                  status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j
    
    ! get flag:
    call self%GetVariable( varnames(o_icc), status, enabled=enabled )
    IF_NOTOK_RETURN(status=1)
    ! enabled?
    if ( enabled ) then
      ! pointer to target field given the name:
      call self%GetVariable( varnames(o_icc ), status, Instant_Field=iF, check_units='1' )
      IF_NOTOK_RETURN(status=1)
      ! store:
      call iF%PutSample( icc, t, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! get flag:
    call self%GetVariable( varnames(o_bcc), status, enabled=enabled )
    IF_NOTOK_RETURN(status=1)
    ! enabled?
    if ( enabled ) then
      ! pointer to target field given the name:
      call self%GetVariable( varnames(o_bcc ), status, Instant_Field=iF, check_units='1' )
      IF_NOTOK_RETURN(status=1)
      ! store:
      call iF%PutSample( bcc, t, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! get flag:
    call self%GetVariable( varnames(o_occ), status, enabled=enabled )
    IF_NOTOK_RETURN(status=1)
    ! enabled?
    if ( enabled ) then
      ! pointer to target field given the name:
      call self%GetVariable( varnames(o_occ), status, Instant_Field=iF, check_units='1' )
      IF_NOTOK_RETURN(status=1)
      ! store:
      call iF%PutSample( occ, t, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! get flag:
    call self%GetVariable( varnames(o_clwc), status, enabled=enabled )
    IF_NOTOK_RETURN(status=1)
    ! enabled?
    if ( enabled ) then
      ! pointer to target field given the name:
      call self%GetVariable( varnames(o_clwc), status, Instant_Field=iF, check_units='kg/kg' )
      IF_NOTOK_RETURN(status=1)
      ! store:
      call iF%PutSample( clwc, t, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! put out ice water?
    if ( args(o_ciwc) /= 'None' ) then
      ! get flag:
      call self%GetVariable( varnames(o_ciwc), status, enabled=enabled )
      IF_NOTOK_RETURN(status=1)
      ! enabled?
      if ( enabled ) then
        ! pointer to target field given the name:
        call self%GetVariable( varnames(o_ciwc), status, Instant_Field=iF, check_units='kg/kg' )
        IF_NOTOK_RETURN(status=1)
        ! store:
        call iF%PutSample( ciwc, t, status )
        IF_NOTOK_RETURN(status=1)
      end if
    end if

    ! get flag:
    call self%GetVariable( varnames(o_pri), status, enabled=enabled )
    IF_NOTOK_RETURN(status=1)
    ! enabled?
    if ( enabled ) then
      ! pointer to target field given the name:
      call self%GetVariable( varnames(o_pri ), status, Instant_Field=iF, check_units='1' )
      IF_NOTOK_RETURN(status=1)
      ! store:
      call iF%PutSample( pri, t, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! clear:
    deallocate(  icc, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(  bcc, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(  occ, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( clwc, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( ciwc, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate(  pri, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate(   cc_met1, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( clwc_met1, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( ciwc_met1, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    if ( args(i_ciwc_met) == 'None' ) then
      ! clear:
      deallocate( ciwc_met, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine Variables_Calc_CloudProfiles


  ! ***


  !
  ! Cloud profiles
  !   input : icc, clwc (or ciwc)
  !   output: iclwc  (or iciwc)
  !

  subroutine Variables_Calc_InCloudWaterContent( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_InCloudWaterContent'
    
    ! argument names:
    integer, parameter        ::  i_icc = 1
    integer, parameter        ::  i_cwc = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  icc(:,:,:)
    real, pointer                       ::  cwc(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='kg/kg', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_icc ), t,t, icc , status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_cwc), t,t, cwc, status, check_units='kg/kg', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target fields (full levels):
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do k = lbo(3), ubo(3)
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          ! in-cloud liquid water content:
          if ( icc(i,j,k) > 0.0 ) then
            values(i,j,k) = cwc(i,j,k) / icc(i,j,k)
          else
            values(i,j,k) = 0.0
          end if
        end do ! i
      end do ! j
    end do ! k

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_InCloudWaterContent


  ! ***


  !
  ! Average volumetric soil water:
  !   output: sd [m water eqv.]
  !   input : sd [m snow]
  !

  subroutine Variables_Calc_Snowdepth_lw( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_Snowdepth_lw'
    
    ! argument names:
    integer, parameter        ::  i_sdphys  = 1

    ! densities:
    real, parameter  :: rho_snow  =  109.0  ! kg/m3 first order approximation from IFS
    real, parameter  :: rho_water = 1000.0  ! kg/m3

    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  sdphys(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! init to zero:
    values = 0.0

    ! input data:
    call self%GetValidPtr( args(i_sdphys), t,t, sdphys, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! loop over grid cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        !                 (m snow)      kg/(m3 snow) / kg/(m3 water)
        values(i,j,1) = sdphys(i,j,1) *   rho_snow   /   rho_water   ! m water
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_Snowdepth_lw


  ! ***


  !
  ! liquid water equivalent (devide by rho_water)
  !   input :  kg/m2/s
  !   output:  m/s
  !

  subroutine Variables_Calc_LWE( self, varname, args, t1, t2, status )
  
    use GO, only : TDate
    use GO, only : T_Constant_Field
    use GO, only : GoVarValue
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t1, t2
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_LWE'
    
    ! argument names:
    integer, parameter        ::  i_rain_amount  = 1

    ! densities:
    real, parameter  :: rho_water = 1000.0  ! kg/m3

    ! --- local ----------------------------------
    
    class(T_Constant_Field), pointer    ::  cF
    real, pointer                       ::  water(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Constant_Field=cF, &
                              check_units='m/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! input data, should be valid for requested time interval:
    call self%GetValidPtr( args(i_rain_amount), t1,t2, water, status, &
                             check_units='kg/m2/s', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! loop over target grid cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        !                  kg/m2     /   kg/m3  
        values(i,j,1) = water(i,j,1) / rho_water   ! m
      end do ! i
    end do ! j

    ! store:
    call cF%PutSample( values, (/t1,t2/), status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_LWE


  ! ***


  !
  ! Average volumetric soil water:
  !   output: swvl1t2 or swvl1t3
  !   input : swvl1, swvl2, swvl3 or None
  !

  subroutine Variables_Calc_VolumetricSoilWater( self, varname, args, t, status )
  
    use GO, only : TDate
    use GO, only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_VolumetricSoilWater'
    
    ! hard-coded layer depths, should be passed with data in future:
    real, parameter  ::  depths(0:3) = (/ 0.0, 7.0, 28.0, 100.0 /) ! cm
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  swvl(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, k
    real                                ::  thickness
    real                                ::  tdepth
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='m3/m3', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! init to zero:
    values = 0.0
    ! init total depth:
    tdepth = 0.0
    ! loop:
    do k = 1, size(args)
      ! skip if not defined:
      if ( args(k) == 'None' ) cycle
      ! input data:
      call self%GetValidPtr( args(k), t,t, swvl, status, check_units='m3/m3', &
                               check_lbo=lbo, check_ubo=ubo )
      IF_NOTOK_RETURN(status=1)
      ! current thickness:
      thickness = depths(k) - depths(0)
      ! add contribution:
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)
          values(i,j,1) = values(i,j,1) + thickness * swvl(i,j,1)
        end do ! i
      end do ! j
      tdepth = tdepth + thickness
    end do
    ! average, could not be zero since at least first value is checked to be defined:
    values = values / tdepth

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_VolumetricSoilWater


  ! ***


  !
  ! Extract volumetric soil water in layers
  !   
  !   input : 'swvl', 'layer=1'
  !   output : swvl1 or swvl2 or swvl3 or swvl4

  subroutine Variables_Calc_Extract_swvl( self, varname, args, t, status )
 
    use GO, only : TDate
    use GO, only : T_Instant_Field
    use GO          , only : GoVarValue
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_Extract_swvl'
    
    ! arguments:
    integer, parameter        ::  i_swvl= 1
    integer, parameter        ::  i_Nlayer =2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  swvl(:,:,:)
    integer                             ::  nlayer
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
   
    ! --- begin ----------------------------------

    ! get pointer to source array:
    call self%GetValidPtr( args(i_swvl), t,t, swvl, status, &
                              check_units='m3/m3', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                              check_units='m3/m3', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! extract layer index:
    call GoVarValue( args(i_Nlayer), ';', 'nlayer','=', nlayer, status )
    if ( status /= 0 ) then
      write( gol, '("Format should be: nlayer=N, found: ", a )' ) trim(args(i_Nlayer)) ; call GoErr
      TRACEBACK;status=1;return
    end if
   
    ! copy:    
    values(:,:,1) = swvl(:,:,nlayer)

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
     
  end subroutine Variables_Calc_Extract_swvl


  ! ***

  !
  ! Gravimetric soil water:
  !   output: swg
  !   input : swv, lsm
  !

  subroutine Variables_Calc_GravimetricSoilWater( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use Binas          , only : rho_water   ! kg/m3
    use LE_LandUse_Soil, only : rho_soil    ! kg/m3
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_GravimetricSoilWater'
    
    ! --- begin ----------------------------------
    
    ! argument names:
    integer, parameter        ::  i_swv = 1
    integer, parameter        ::  i_lsm = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  swv(:,:,:)
    real, pointer                       ::  lsm (:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, k
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( varname, status, Instant_Field=iF, &
                              check_units='kg/kg', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! dependencies:
    call self%GetValidPtr( args(i_swv), t,t, swv, status, check_units='m3/m3', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_lsm ), t,t, lsm , status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do k = lbo(3), ubo(3)
      do j = lbo(2), ubo(2)
        do i = lbo(1), ubo(1)

          ! Use ECMWF soil volumetric water content (m3/m3) and the assumed soil density.
          ! Note that here we use the ECMWF soil water, which is valid for a
          ! a layer of about 25 cm (need to check this).
          ! For dust we need the real top only however, which is expected
          ! to contain much less water in case of sun, and much more in case of rain.
          !  kg/kg           m3/m3         kg/m3          (kg water)/m3                 (kg soil)/m3
          values(i,j,k) =  swv(i,j,k)  * rho_water / ( swv(i,j,k) * rho_water  + (1.0-swv(i,j,k))*rho_soil )

          ! ADHOC: add water part, by default no soil over sea so swv is zero ...
          ! lsm is 1 over land and 0 over sea, at coast lines in between:
          !                 kg/kg       land-fraction  kg/kg    sea-fraction
          values(i,j,k) = values(i,j,k) * lsm(i,j,1) +  1.0  * (1.0-lsm(i,j,1))

        end do ! i
      end do ! j
    end do ! k

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_GravimetricSoilWater


  ! ***


  !
  ! Average volumetric soil water:
  !   output: smi
  !   input : swvl, lsm, slt
  !

  subroutine Variables_Calc_SoilMoistureIndex( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use LE_Data_Calc, only : SoilMoistureIndex
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_SoilMoistureIndex'
    
    ! argument names:
    integer, parameter        ::  i_swvl = 1
    integer, parameter        ::  i_lsm  = 2
    integer, parameter        ::  i_slt  = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  swvl(:,:,:)
    real, pointer                       ::  slt (:,:,:)
    real, pointer                       ::  lsm (:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='1', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_swvl), t,t, swvl, status, check_units='m3/m3', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_lsm ), t,t, lsm , status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_slt ), t,t, slt , status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! soil moisture index for this cell ;
        ! pass integer soil type and time value 
        ! to decide on interpretation of soil type:
        call SoilMoistureIndex( swvl(i,j,1), lsm(i,j,1), nint(slt(i,j,1)), t, &
                                 values(i,j,1), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_SoilMoistureIndex
  

  ! ***


  !
  ! Average volumetric soil water:
  !   output: smi
  !   input : swvl, lsm, slt
  !

  subroutine Variables_Calc_SoilMoistureIndex_WRF( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use LE_Data_Calc, only : SoilMoistureIndex_wrf
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_SoilMoistureIndex_WRF'
    
    ! argument names:
    integer, parameter        ::  i_swvl = 1
    integer, parameter        ::  i_lsm  = 2
    integer, parameter        ::  i_slt  = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  swvl(:,:,:)
    real, pointer                       ::  slt (:,:,:)
    real, pointer                       ::  lsm (:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='1', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_swvl), t,t, swvl, status, check_units='m3/m3', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_lsm ), t,t, lsm , status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_slt ), t,t, slt , status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! soil moisture index for this cell ;
        ! pass integer soil type and time value 
        ! to decide on interpretation of soil type:
        call SoilMoistureIndex_WRF( swvl(i,j,1), lsm(i,j,1), nint(slt(i,j,1)), t, &
                                 values(i,j,1), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_SoilMoistureIndex_WRF
  
  ! ***

  !
  ! Exponential moving average ;
  ! new average is combination of previous average and newe value:
  !    var_ema[t]  =  (1-alpha) * var_ema[t-1]  +  alpha * var[t-1]
  ! Arguments:
  !   input :  var, alpha
  !   output:  var_ema
  !

  subroutine Variables_CalcEMA( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Constant_Field
    use GO          , only : GoVarValue
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_CalcEMA'
    
    ! argument names:
    integer, parameter        ::  i_param = 1
    integer, parameter        ::  i_alpha = 2
    
    ! --- local ----------------------------------
    
    class(T_Constant_Field), pointer    ::  cF
    real, pointer                       ::  field(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  values_o(:,:,:)
    real, allocatable                   ::  values(:,:,:)
    logical                             ::  filled
    integer                             ::  i, j
    character(len=32)                   ::  units
    real                                ::  alpha
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Constant_Field=cF, &
                            units=units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_param), t,t, field, status, check_units=trim(units), &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! read alpha:
    call GoVarValue( args(i_alpha), ';', 'alpha','=', alpha, status )
    if ( status /= 0 ) then
      write( gol, '("format of argument should be `alpha=0.123456`, found `",a,"`")' ) trim(args(i_alpha)) ; call GoErr
      TRACEBACK;status=1;return
    end if
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! get current data
    call cF%Get( status, pdata=values_o, filled=filled )
    IF_NOTOK_RETURN(status=1)

    ! If not filled yet (first hour) point to current temperature
    if ( .not. filled ) then
      call self%GetValidPtr( args(i_param), t,t, values_o, status, check_units=trim(units), &
                               check_lbo=lbo, check_ubo=ubo )
      IF_NOTOK_RETURN(status=1)
    end if                               
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! update Exponential moving average 
        values(i,j,1) = (1-alpha) * values_o(i,j,1) + alpha* field(i,j,1)
      end do ! i
    end do ! j

    ! store:
    call cF%PutSample( values, (/t,t/), status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_CalcEMA


  !
  ! Change sign ;
  ! Arguments:
  !   input :  var
  !   output:  var * -1.0
  !

  subroutine Variables_ChangeSign( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use GO          , only : GoVarValue
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_ChangeSign'
    
    ! argument names:
    integer, parameter        ::  i_param = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  field(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  values_o(:,:,:)
    real, allocatable                   ::  values(:,:,:)
    logical                             ::  filled
    integer                             ::  i, j
    character(len=32)                   ::  units
    real                                ::  alpha
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            units=units, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_param), t,t, field, status, check_units=trim(units), &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
        
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! get current data
    call iF%Get( status, pdata=values_o, filled=filled )
    IF_NOTOK_RETURN(status=1)

    ! If not filled yet (first hour) point to current temperature
    if ( .not. filled ) then
      call self%GetValidPtr( args(i_param), t,t, values_o, status, check_units=trim(units), &
                               check_lbo=lbo, check_ubo=ubo )
      IF_NOTOK_RETURN(status=1)
    end if                               
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! change sign:
        values(i,j,1) = -1.0 * values_o(i,j,1)
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_ChangeSign
  
 
  !
  ! LeafAreaIndex
  !   output: lai_lu
  !   input : lat
  !

  subroutine Variables_Calc_LeafAreaIndex( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_Data, only : start_gs, end_gs
    use LE_Landuse     , only : Setup_GS
    use LE_Landuse     , only : Setup_LAI
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_LeafAreaIndex'
 
    ! argument names:
    integer, parameter        ::  i_lat    = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  lai_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  lat(:,:,:)      ! (nx,ny,1)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            check_units='m2/m2', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! (re)set level dimension:
    lbo(3) = 1
    ubo(3) = nlu
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( lai_lu, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_lat), t,t, lat, status, &
                             check_units='degrees_north', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)

    ! update information on growing seasons:
    call Setup_GS( t, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! fill:
        call Setup_LAI( lai_lu(i,j,:), lat(i,j,1), t, &
                        lu_fracs(i,j,:), start_gs(i,j,:), end_gs(i,j,:), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Variables_Calc_LeafAreaIndex
  
 
  !
  ! SurfaceAreaIndex
  !   output: sai_lu
  !   input : lai_lu, lat
  !

  subroutine Variables_Calc_SurfaceAreaIndex( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_Data, only : start_gs, end_gs
    use LE_Landuse     , only : Setup_GS
    use LE_Landuse     , only : Setup_SAI
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_SurfaceAreaIndex'
 
    ! argument names:
    integer, parameter        ::  i_lai_lu = 1
    integer, parameter        ::  i_lat    = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  sai_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  lai_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  lat(:,:,:)      ! (nx,ny,1)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            check_units='m2/m2', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! (re)set level dimension:
    lbo(3) = 1
    ubo(3) = nlu
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( sai_lu, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_lai_lu), t,t, lai_lu, status, &
                             check_units='m2/m2', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_lat), t,t, lat, status, &
                             check_units='degrees_north', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)

    ! update information on growing seasons:
    call Setup_GS( t, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! fill:
        call Setup_SAI( sai_lu(i,j,:), lai_lu(i,j,:), lat(i,j,1), t, &
                        lu_fracs(i,j,:), start_gs(i,j,:), end_gs(i,j,:), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Variables_Calc_SurfaceAreaIndex
  
 
  !
  ! SurfaceRoughness_Mass_LU
  !   output: z0m_lu
  !   input : lai_lu, wsurf, sd
  !

  subroutine Variables_Calc_SurfaceRoughness_Mass_LU( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_z0  , only : Setup_z0m
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_SurfaceRoughness_Mass_LU'
 
    ! argument names:
    integer, parameter        ::  i_lai_lu = 1
    integer, parameter        ::  i_wsurf  = 2
    integer, parameter        ::  i_sd     = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  z0m_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  lai_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  wsurf(:,:,:)    ! (nx,ny,1)
    real, pointer                       ::  sd   (:,:,:)    ! (nx,ny,1)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! (re)set level dimension:
    lbo(3) = 1
    ubo(3) = nlu
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( z0m_lu, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_lai_lu), t,t, lai_lu, status, &
                             check_units='m2/m2', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_wsurf), t,t, wsurf, status, &
                             check_units='m/s', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_sd), t,t, sd, status, &
                             check_units='m', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! fill:
        call Setup_z0m( z0m_lu(i,j,:), lai_lu(i,j,:), wsurf(i,j,1), sd(i,j,1), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Variables_Calc_SurfaceRoughness_Mass_LU
  
 
  !
  ! SurfaceRoughness_Heat_LU
  !   output: z0h_lu
  !   input : lai_lu, wsurf, sd
  !

  subroutine Variables_Calc_SurfaceRoughness_Heat_LU( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_z0  , only : Setup_z0h
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_SurfaceRoughness_Heat_LU'
 
    ! argument names:
    integer, parameter        ::  i_lai_lu = 1
    integer, parameter        ::  i_wsurf  = 2
    integer, parameter        ::  i_sd     = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  z0h_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  lai_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  wsurf(:,:,:)    ! (nx,ny,1)
    real, pointer                       ::  sd   (:,:,:)    ! (nx,ny,1)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! (re)set level dimension:
    lbo(3) = 1
    ubo(3) = nlu
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( z0h_lu, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_lai_lu), t,t, lai_lu, status, &
                             check_units='m2/m2', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_wsurf), t,t, wsurf, status, &
                             check_units='m/s', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_sd), t,t, sd, status, &
                             check_units='m', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! fill:
        call Setup_z0h( z0h_lu(i,j,:), lai_lu(i,j,:), wsurf(i,j,1), sd(i,j,1), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Variables_Calc_SurfaceRoughness_Heat_LU
  
 
  !
  ! canopy height
  !   output: zcanopy_lu
  !   input : z0h_lu, lai_lu
  !

  subroutine Variables_Calc_CanopyHeight_LU( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_z0  , only : Setup_zcanopy
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_CanopyHeight_LU'
 
    ! argument names:
    integer, parameter        ::  i_z0h_lu = 1
    integer, parameter        ::  i_lai_lu = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  zcanopy_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  z0h_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  lai_lu(:,:,:)   ! (nx,ny,nlu)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! (re)set level dimension:
    lbo(3) = 1
    ubo(3) = nlu
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( zcanopy_lu, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_z0h_lu), t,t, z0h_lu, status, &
                             check_units='m', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_lai_lu), t,t, lai_lu, status, &
                             check_units='m2/m2', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! fill:
        call Setup_zcanopy( zcanopy_lu(i,j,:), z0h_lu(i,j,:), lai_lu(i,j,:), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Variables_Calc_CanopyHeight_LU
  
 
  !
  ! SurfaceRoughness_Dust_LU
  !   output: z0dust_lu
  !   input : lai_lu
  !

  subroutine Variables_Calc_SurfaceRoughness_Dust_LU( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_z0  , only : Setup_z0dust
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_SurfaceRoughness_Dust_LU'
 
    ! argument names:
    integer, parameter        ::  i_lai_lu = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  z0dust_lu(:,:,:)   ! (nx,ny,nlu)
    real, pointer                       ::  lai_lu(:,:,:)   ! (nx,ny,nlu)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! (re)set level dimension:
    lbo(3) = 1
    ubo(3) = nlu
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( z0dust_lu, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_lai_lu), t,t, lai_lu, status, &
                             check_units='m2/m2', check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! fill:
        call Setup_z0dust( z0dust_lu(i,j,:), lai_lu(i,j,:), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Variables_Calc_SurfaceRoughness_Dust_LU
  
 
  !
  ! cell averaged SurfaceRoughness
  !   output: z0
  !   input : z0_lu
  !

  subroutine Variables_Calc_SurfaceRoughness_CellAverage( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_z0  , only : Average_z0
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_SurfaceRoughness_CellAverage'
 
    ! argument names:
    integer, parameter        ::  i_z0_lu = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  z0(:,:,:)      ! (nx,ny,1)
    real, pointer                       ::  z0_lu(:,:,:)   ! (nx,ny,nlu)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                            check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( z0, t, status )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_z0_lu), t,t, z0_lu, status, &
                             check_units='m', &
                             check_lbo=(/lbo(1),lbo(2),nlu/), &
                             check_ubo=(/ubo(1),ubo(2),nlu/) )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! fill:
        call Average_z0( z0_lu(i,j,:), lu_fracs(i,j,:), z0(i,j,1), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Variables_Calc_SurfaceRoughness_CellAverage
  
 
  !
  !
  ! Stability class (1-6)
  !   output: index
  !   input : lon, wsurf, tcc, radd, sd
  !

  subroutine Variables_Calc_StabilityClass( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use binas       , only : grav, vkarman
    
    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_Data, only : ilu_water_sea
  
    use dims, only   : coszen
  
    use JAQL_stability, only : exposure_new
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(5)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_StabilityClass'
 
    ! argument names:
    integer, parameter        ::  i_lon   = 1
    integer, parameter        ::  i_wsurf = 2
    integer, parameter        ::  i_tcc   = 3
    integer, parameter        ::  i_radd  = 4
    integer, parameter        ::  i_sd    = 5
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  lon(:,:,:)
    real, pointer                       ::  wsurf(:,:,:)
    real, pointer                       ::  tcc(:,:,:)
    real, pointer                       ::  radd(:,:,:)
    real, pointer                       ::  sd(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  stabcls(:,:,:)
    integer                             ::  i, j
    real                                ::  lon_local
    integer                             ::  dt_local, h_local
    integer                             ::  stab_class
    character(len=1)                    ::  stab_class_letter
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='1', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( stabcls, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_lon), t,t, lon, status, check_units='degrees_east', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_wsurf), t,t, wsurf, status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_tcc ), t,t, tcc , status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)    
    call self%GetValidPtr( args(i_radd ), t,t, radd , status, check_units='J/m2/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_sd ), t,t, sd , status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)

        ! find hour of local time (needed to check sunrise/sunset )
        lon_local = lon(i,j,1)
        dt_local  = floor(lon_local/15.)
        h_local   = modulo(t%hour + dt_local, 24)

        ! Exposure class via radiation/wspd and clouds (table 16.3 in S&P )
        call exposure_new( radd(i,j,1), wsurf(i,j,1), tcc(i,j,1), coszen(i,j), sd(i,j,1), h_local, &
                             stab_class_letter, stab_class, status )
        IF_NOTOK_RETURN(status=1)

        ! Restrict stability class above water to slightly stable, neutral or slighty unstable
        if ( lu_fracs(i,j,ilu_water_sea) > 0.5 .and. ( stab_class_letter == 'A' .or. stab_class_letter == 'B' ) ) then
          stab_class_letter = 'C'
          stab_class        =  3
        else if ( lu_fracs(i,j,ilu_water_sea) > 0.5 .and. stab_class_letter == 'F' ) then
          stab_class_letter = 'E'
          stab_class        =  5
        end if
        
        ! store:
        stabcls(i,j,1) = stab_class
              
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( stabcls, t, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Variables_Calc_StabilityClass
  
 
  !
  !
  ! Monin-Obukhov-Length:
  !   output: monin (grass)
  !   input : z0m, stabcls
  !

  subroutine Variables_Calc_MoninObukhovLengthExposure( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use binas       , only : grav, vkarman
    
    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_Data, only : ilu_water_sea
  
    use dims, only   : coszen
  
    use JAQL_stability, only : exposure_new, CalcMonin_stabclass
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(6)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_MoninObukhovLengthExposure'
 
    ! argument names:
    integer, parameter            ::  i_z0m     = 1
    integer, parameter            ::  i_stabcls = 2 

    ! Stability calculations with z0 above 0.75 are cut-off to 0.75
    ! Curves from S&P only valid for low z0-values, above keep stability regime
    real, parameter               ::  z0_thr = 0.75
    
    !! minimum L to avoid extreem stable conditions:
    !real, parameter               ::  Lmin = 100.0
    
    ! letter codes:
    character(len=1), parameter   ::  stab_class_letter(6) = (/'A','B','C','D','E','F'/)
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  z0m(:,:,:)
    real, pointer                       ::  stabcls(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  MOL(:,:,:)
    integer                             ::  i, j
    integer                             ::  stab_class
    real                                ::  L
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( MOL, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_z0m), t,t, z0m, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_stabcls), t,t, stabcls, status, check_units='1', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)

        ! integer stability class:
        stab_class = nint( stabcls(i,j,1) )

        ! Calculate monin obuhkov length via table 16.4 in S&P
        call CalcMonin_StabClass( stab_class_letter(stab_class), min(z0m(i,j,1),z0_thr), L, status )
        IF_NOTOK_RETURN(status=1)

        !! stable? then use minimum value:
        !if ( L > 0.0 ) then
        !  ! not too stable ...
        !  L = max( Lmin, L )
        !end if
        
        ! store:
        MOL(i,j,1) = L
              
      end do ! i
    end do ! j

    ! ok
    status = 0
    
  end subroutine Variables_Calc_MoninObukhovLengthExposure
  
  !
  ! Ustar-average:
  !   output: ustar
  !   input : z0m_lu, z0m, monin_inv (calculated by exposure classes), wsurf
  !

  subroutine Variables_Calc_UstarJacobsen( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use binas          , only : grav, vkarman, kappa_stab
    use LE_Landuse_Data, only : ilu_grass
    use JAQL_Stability , only : f_m_stability
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(4)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_UstarJacobsen'
      
    ! argument names:
    integer, parameter        ::  i_z0m_lu    = 1
    integer, parameter        ::  i_z0m       = 2
    integer, parameter        ::  i_monin_inv = 3
    integer, parameter        ::  i_wsurf     = 4
  
    real, parameter :: zref  = 10.0    ! zref = 10m, the wind speed is given at 10m height
                                       
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  z0m_lu(:,:,:)
    real, pointer                       ::  z0m(:,:,:)
    real, pointer                       ::  monin_inv(:,:,:)
    real, pointer                       ::  wsurf(:,:,:)    
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
    real                                ::  fm  
    
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input:
    call self%GetValidPtr( args(i_z0m_lu), t,t, z0m_lu, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_z0m), t,t, z0m, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_monin_inv), t,t, monin_inv , status, check_units='1/m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_wsurf), t,t, wsurf , status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)

        fm = f_m_stability( zref, z0m(i,j,1), monin_inv(i,j,1) )      
        values(i,j,1) = kappa_stab * wsurf(i,j,1) / fm
                       
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_UstarJacobsen
  
  !
  !
  ! Monin-Obukhov-Length:
  !   output: monin (grass)
  !   input : ustar_grass,t,Q0v
  !

  subroutine Variables_Calc_MoninObukhovLengthIFS( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use binas       , only : grav, vkarman
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_MoninObukhovLengthIFS'
    
    ! argument names:
    integer, parameter        ::  i_zust = 1
    integer, parameter        ::  i_temp = 2
    integer, parameter        ::  i_Q0v  = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  zust(:,:,:)
    real, pointer                       ::  tsurf(:,:,:)
    real, pointer                       ::  Q0v(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_zust), t,t, zust, status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_Q0v ), t,t, Q0v , status, check_units='K m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! surface bound checks
    lbo_sfc = (/lbo(1),lbo(2), 1/)
    ubo_sfc = (/ubo(1),ubo(2), 1/)
    
    call self%GetValidPtr( args(i_temp ), t,t, tsurf , status, check_units='K', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! monin obukhov length for this cell ; valid for grass, because zust: u* for grass
        values(i,j,1) = -zust(i,j,1)**3 / ( ( grav*vkarman)/ tsurf(i,j,1) * Q0v(i,j,1) )
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_MoninObukhovLengthIFS
  
  !
  ! Inverse of Monin-Obukhov-Length 
  !   output: monin_inv 
  !   input : monin  
  !

  subroutine Variables_Calc_InverseMonin( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_InverseMonin'
    
    ! argument names:
    integer, parameter        ::  i_monin = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  monin(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='1/m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_monin), t,t, monin, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! Inverse of monin obukhov length for this cell ;
        if ( monin(i,j,1) == 0 ) then
          values(i,j,1) = 1.0e6
        else 
          values(i,j,1) = 1.0/monin(i,j,1)
        end if
        
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_InverseMonin

  !
  ! Calculate Monin from Inverse Monin-Obukhov-Length 
  !   output: monin
  !   input : monin_inv  
  !

  subroutine Variables_Calc_InverseMoninInverse( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(1)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_InverseMoninInverse'
    
    ! argument names:
    integer, parameter        ::  i_monin_inv = 1
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  monin_inv(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_monin_inv), t,t, monin_inv, status, check_units='1/m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! Inverse of monin obukhov length for this cell ;
        if ( monin_inv(i,j,1) > 1e6 .or. monin_inv(i,j,1) < -1e6 ) then
          values(i,j,1) = 0.0
        else 
          values(i,j,1) = 1.0/monin_inv(i,j,1)
        end if
        
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_InverseMoninInverse

  ! Ustar-average:
  !   output: ustar
  !   input : z0m_lu, z0m, ustar_grass, monin_inv for grass 
  !

  subroutine Variables_Calc_UstarGrs2Aver( self, varname, args, t, status )
  
    use GO             , only : TDate
    use GO             , only : T_Instant_Field
    use binas          , only : grav, vkarman, kappa_stab
    use LE_Landuse_Data, only : ilu_grass
    use JAQL_Stability , only : f_m_stability

    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(4)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_UstarGrs2Aver'
      
    ! argument names:
    integer, parameter        ::  i_z0m_lu    = 1
    integer, parameter        ::  i_z0m       = 2
    integer, parameter        ::  i_zust      = 3
    integer, parameter        ::  i_monin_inv = 4
  
    real, parameter :: zref50 = 50.0   ! zref50 = 50 m is introduced as
                                       !   the reference height at which
                                       !   the logaritmic velocity profile
                                       !   is assumed to be independent
                                       !   upon land use class.
                                       
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  z0m_lu(:,:,:)
    real, pointer                       ::  z0m(:,:,:)
    real, pointer                       ::  zust(:,:,:)    ! ustar for grass
    real, pointer                       ::  monin_inv(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
    real                                ::  fm  
    real                                ::  u_zref50
    
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_z0m_lu), t,t, z0m_lu, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_z0m), t,t, z0m, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_zust), t,t, zust, status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_monin_inv), t,t, monin_inv , status, check_units='1/m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
      
        ! calculate free wind from ustar/monin for grass
        fm = f_m_stability( zref50, z0m_lu(i,j,ilu_grass), monin_inv(i,j,1) )
        u_zref50 = zust(i,j,1) * fm / kappa_stab
        
        ! calculate ustar on other z0 (with monin from grass)
        fm = f_m_stability( zref50, z0m(i,j,1), monin_inv(i,j,1) )
        values(i,j,1) = kappa_stab * u_zref50 / fm

      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_UstarGrs2Aver

  !
  ! Scale parameter Specific humidity
  !   output: qstar
  !   input : slhf,dens,ustar
  !

  subroutine Variables_Calc_ScaleParameterSpecificHumidity( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use Binas       , only : lvap
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_ScaleParameterSpecificHumidity'
    
    ! argument names:
    integer, parameter        ::  i_slhf  = 1
    integer, parameter        ::  i_dens  = 2
    integer, parameter        ::  i_zust  = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  slhf(:,:,:)
    real, pointer                       ::  dens(:,:,:)
    real, pointer                       ::  zust(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='kg/kg', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_slhf), t,t, slhf, status, check_units='J/m2/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_dens ), t,t, dens , status, check_units='kg/m3', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_zust ), t,t, zust , status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! Scale parameter specific humidity;
        values(i,j,1) = -slhf(i,j,1) / lvap / dens(i,j,1) / zust(i,j,1)
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_ScaleParameterSpecificHumidity

  !
  ! Scale parameter Dry static energy
  !   output: sstar
  !   input : sshf, dens,ustar
  !

  subroutine Variables_Calc_ScaleParameterDryStaticEnergy( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_ScaleParameterDryStaticEnergy'
    
    ! argument names:
    integer, parameter        ::  i_sshf  = 1
    integer, parameter        ::  i_dens  = 2
    integer, parameter        ::  i_zust  = 3
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  sshf(:,:,:)
    real, pointer                       ::  dens(:,:,:)
    real, pointer                       ::  zust(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='J/kg', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_sshf), t,t, sshf, status, check_units='J/m2/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    call self%GetValidPtr( args(i_dens), t,t, dens , status, check_units='kg/m3', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_zust), t,t, zust , status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),1), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        ! Scale parameter dry static energy:
        values(i,j,1) = -sshf(i,j,1) /dens(i,j,1)/ zust(i,j,1)
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_ScaleParameterDryStaticEnergy
  
  !
  ! Half level gradient
  !   output: duudz
  !   input : uu_ifs, dh
  !

  subroutine Variables_Calc_HalfLevelGradient( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field

    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_HalfLevelGradient'
    
    ! argument names:
    integer, parameter        ::  i_uu_ifs = 1
    integer, parameter        ::  i_dh_met = 2
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  uu_ifs(:,:,:)
    real, pointer                       ::  dh_met(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_out(3), ubo_out(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, l
      
    ! --- begin ----------------------------------
    
    ! pointers to first input dat, get bounds:
    call self%GetValidPtr( args(i_uu_ifs), t,t, uu_ifs, status, check_units='m/s', &
                                 lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_dh_met), t,t, dh_met, status, check_units='m', &
                                 check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! Add level 0 (surface to output
    lbo_out = (/ lbo(1), lbo(2), 0 /)
    ubo_out = (/ ubo(1), ubo(2), ubo(3)/)
    
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='1/s')
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo_out(1):ubo_out(1),lbo_out(2):ubo_out(2),lbo_out(3):ubo_out(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo_out(2), ubo_out(2)
      do i = lbo_out(1), ubo_out(1)
        do l = lbo_out(3), ubo_out(3)
          if ( l == lbo_out(3) ) then
            ! No flux through surface
            values(i,j,l) = 0.0
          else if ( l == ubo_out(3) ) then
            ! No flux through top of input
            values(i,j,l) = 0.0
          else
            ! Half level gradient
            values(i,j,l) = ( uu_ifs(i,j,l+1) - uu_ifs(i,j,l) ) / dh_met(i,j,l+1)
          end if
        end do ! l
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status, lbo=lbo_out, ubo=ubo_out )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_HalfLevelGradient
  
  !
  ! Bottom windspeed
  !   output: uu_ifs
  !   input : t, uu, Q0v
  !
  ! Absolute windspeed (squared) at the bottom
  ! model layer following eq. (3.17):
  !
  !    |U_n|^2 = u_n^2 + v_n^2 + w_*^2
  !
  ! with w_* the free convection velocity scale defined by:
  !
  !                g       1/3
  !   w_* = ( z_i --- Q0v )
  !               T_n
  !
  ! Note that Q0v could be negative (which gives L>0, stable conditions).
  ! Therefore it is better to formulate the convection velocity scale as:
  !
  !        3        g
  !   (w_*)  = z_i --- Q0v
  !                T_n
  ! and
  !           |      g      | 2/3
  !   w_*^2 = | z_i --- Q0v |
  !           |     T_n     |
  !
  ! Parameters:
  !    z_i   : scale height of the boundary layer depth (1000 m)
  !    g     : acceleration of gravity
  !
  ! Arguments:
  !    u_n, v_n  : horizontal wind components at lowest model layer
  !    T_n       : temperature                at lowest model layer
  !    Q0v       : virtual temperature flux in the surface layer [K m/s ?]
  !

  subroutine Variables_Calc_BottomWindSpeed( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use binas       , only : grav

    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(3)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_BottomWindSpeed'
    
    ! argument names:
    integer, parameter        ::  i_temp  = 1
    integer, parameter        ::  i_uu    = 2
    integer, parameter        ::  i_Q0v   = 3
    
    real, parameter           ::  z_i = 1.e3  ! scale height of the boundary layer depth and is set to constant 
                                              ! value of 1000 m, since only the order of magnitude matters
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  temper(:,:,:)
    real, pointer                       ::  uu(:,:,:)
    real, pointer                       ::  Q0v(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, l
    real                                ::  wstar ! free convective scale (m/s)
    
      
    ! --- begin ----------------------------------
    
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_uu), t,t, uu, status, check_units='m/s', &
                               lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_temp), t,t, temper, status, check_units='K', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! surface bound checks
    lbo_sfc = (/lbo(1),lbo(2), 1/)
    ubo_sfc = (/ubo(1),ubo(2), 1/)
    
    call self%GetValidPtr( args(i_Q0v), t,t, Q0v, status, check_units='K m/s', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m/s' )
    IF_NOTOK_RETURN(status=1)
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
        do l = lbo(3), ubo(3)

          ! which layer ?
          if ( l == lbo(3) ) then
            ! free convective velocity:
            wstar = abs( z_i * grav / temper(i,j,l) * Q0v(i,j,l) )**(1.0/3.0) 
            ! absolute windspeed at bottom layer incl. free convective velocity:
            values(i,j,l) = sqrt( uu(i,j,l)**2 + wstar**2 )
          else
            ! absolute windspeed:
            values(i,j,l) = uu(i,j,l)
          end if

        end do ! l
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_BottomWindSpeed
  
  !
  ! Virtual Temperature flux
  !   output: Q0v
  !   input : ustar,sstar,t_met,qstar
  !

  subroutine Variables_Calc_VirtualTemperatureFlux( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use binas       , only : cp_air, eps1

    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(4)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_VirtualTemperatureFlux'
    
    ! argument names:
    integer, parameter        ::  i_zust   = 1
    integer, parameter        ::  i_sstar  = 2
    integer, parameter        ::  i_temp   = 3
    integer, parameter        ::  i_qstar  = 4
    
    ! maximum negative Q0v value, used for rounding values close to zero:
    real, parameter           ::  Q0v_maxneg = -1.0e-6
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  zust(:,:,:)
    real, pointer                       ::  sstar(:,:,:)
    real, pointer                       ::  tsurf(:,:,:)
    real, pointer                       ::  qstar(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j
    
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='K m/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_zust), t,t, zust, status, check_units='m/s', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_sstar), t,t, sstar, status, check_units='J/kg', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_temp), t,t, tsurf, status, check_units='K', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_qstar), t,t, qstar, status, check_units='kg/kg', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)    

        ! compute following eq. (3.8) :
        !                   m/s         J/kg   / J/kg/K    1         K          m/s      kg/kg            
        values(i,j,1) = zust(i,j,1)*sstar(i,j,1)/cp_air + eps1*tsurf(i,j,1)*zust(i,j,1)*qstar(i,j,1)  ! K m/s

        ! avoid Q0v==0 which leads to division by zero; 
        ! round small numbers to negatve threshold (which gives L>0, stable conditions):
        if ( abs(values(i,j,1)) < abs(Q0v_maxneg) ) then
          values(i,j,1) = Q0v_maxneg
        end if

      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_VirtualTemperatureFlux
  
  !
  ! Dry static energy
  !   output: dse
  !   input : h,t
  !

  subroutine Variables_Calc_DryStaticEnergy( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use binas       , only : cp_air, grav

    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(2)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_DryStaticEnergy'
    
    ! argument names:
    integer, parameter        ::  i_temp  = 1
    integer, parameter        ::  i_halt  = 2

    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  halt(:,:,:)
    real, pointer                       ::  temper(:,:,:)
    integer                             ::  lbo_in(3), ubo_in(3)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, l
    real                                ::  oro
    real                                ::  z
      
    ! --- begin ----------------------------------
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_temp), t,t, temper, status, check_units='K', &
                               lbo=lbo_in, ubo=ubo_in  )
    IF_NOTOK_RETURN(status=1)
    ! half-level altitudes on different grid:
    call self%GetValidPtr( args(i_halt), t,t, halt, status, check_units='m', &
                               check_lbo=lbo_in-(/0,0,1/), check_ubo=ubo_in )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to target field given the name,
    ! get target shape (horizontal is used):
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='J/kg', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! replace target levels:
    lbo(3) = lbo_in(3)
    ubo(3) = ubo_in(3)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)
      
        ! orography is minimum (needed in case of swapped layers ...)
        oro = min( halt(i,j,0), halt(i,j,ubo(3)) )
        ! loop over target layers:
        do l = lbo(3), ubo(3)
          ! height above surface in mid of layer:
          z = ( halt(i,j,l-1) + halt(i,j,l) )/2.0 - oro
          ! compute:
          !               m/s2   m   J/kg/K    K
          values(i,j,l) = grav * z + cp_air * temper(i,j,l)   ! [J/kg]
        end do ! l

      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status, lbo=lbo,ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_DryStaticEnergy
  
  !
  ! Richardson Number
  !   output: ri
  !   input : h,u,v,q,t,dse
  !

  subroutine Variables_Calc_RichardsonNumber( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use binas       , only : eps, grav, cp_air

    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(6)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_RichardsonNumber'
    
    ! argument names:
    integer, parameter        ::  i_dse   = 1
    integer, parameter        ::  i_halt  = 2
    integer, parameter        ::  i_uv    = 3
    integer, parameter        ::  i_q     = 4
    integer, parameter        ::  i_temp  = 5
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  dse(:,:,:)
    real, pointer                       ::  halt(:,:,:)
    real, pointer                       ::  u(:,:,:)
    real, pointer                       ::  v(:,:,:)
    real, pointer                       ::  q(:,:,:)
    real, pointer                       ::  temper(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, l
    real                                ::  delta_u2  ! [m2/s2]
    real                                ::  ds_cpT    ! [-]
      
    ! --- begin ----------------------------------
    
    ! pointers to input data:
    
    ! input variable, extract units and shape:    
    call self%GetValidPtr( args(i_uv), t,t, u, v, status, check_units='m/s' , &
                               lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_dse), t,t, dse, status, check_units='J/kg', &
                               check_lbo=lbo, check_ubo=ubo)
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_halt), t,t, halt, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_q), t,t, q, status, check_units='kg/kg', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_temp), t,t, temper, status, check_units='K', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='1' )
    IF_NOTOK_RETURN(status=1)
    

    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)       
        do l = lbo(3), ubo(3)
          if ( l == lbo(3) ) then
            values(i,j,l) = 0.0
          else if ( l == ubo(3) ) then     
            values(i,j,l) = 0.0
          else       
            delta_u2 =  ( u(i,j,l+1) - u(i,j,l) )**2 + ( v(i,j,l+1)-v(i,j,l) )**2  
            ds_cpT   =  ( dse(i,j,l+1)-dse(i,j,l) ) / ( cp_air* (temper(i,j,l+1)+temper(i,j,l))/2.) + &
                          eps * ( q(i,j,l+1)-q(i,j,l) )   
                
            values(i,j,l) = grav * ( halt(i,j,l+1) - halt(i,j,l) ) * ds_cpT / delta_u2  ! [-]
          end if
        end do
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_RichardsonNumber
  
  !
  ! HalfLevel Kz IFS
  !   output: kz_ifs
  !   input : temp, dse, uu_ifs, duudz, ri, monin, halt, blh, ustar
  !

  subroutine Variables_Calc_HalfLevelKzIFS( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use JAQL        , only : HalfLevelKzIFS
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(7)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_HalfLevelKzIFS'
    
    ! argument names:
    integer, parameter        ::  i_temp   = 1
    integer, parameter        ::  i_duudz  = 2
    integer, parameter        ::  i_ri     = 3
    integer, parameter        ::  i_monin  = 4
    integer, parameter        ::  i_halt   = 5
    integer, parameter        ::  i_blh    = 6
    integer, parameter        ::  i_ustar  = 7
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    !real, pointer                       ::  dse(:,:,:)
    real, pointer                       ::  temper(:,:,:)
    !real, pointer                       ::  uu_ifs(:,:,:)
    real, pointer                       ::  duudz(:,:,:)
    real, pointer                       ::  ri(:,:,:)
    real, pointer                       ::  monin(:,:,:)
    real, pointer                       ::  halt(:,:,:)
    real, pointer                       ::  blh(:,:,:)
    real, pointer                       ::  ustar(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    integer                             ::  lbo_sfc(3), ubo_sfc(3)
    real, allocatable                   ::  values(:,:,:)
    integer                             ::  i, j, l
    real                                ::  delta_u2  ! [m2/s2]
    real                                ::  ds_cpT    ! [-]
      
    ! --- begin ----------------------------------
    
    
    ! pointers to input data:
    call self%GetValidPtr( args(i_duudz), t,t, duudz, status, check_units='1/s', &
                               lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
        
    !call self%GetValidPtr( args(i_dse), t,t, dse, status, check_units='J/kg', &
    !                           check_lbo=lbo, check_ubo=ubo )
    !IF_NOTOK_RETURN(status=1)
    
    call self%GetValidPtr( args(i_temp), t,t, temper, status, check_units='K', &
                               check_lbo=lbo+(/0,0,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    !call self%GetValidPtr( args(i_uu_ifs), t,t, uu_ifs, status, check_units='m/s', &
    !                           check_lbo=lbo+(/0,0,1/), check_ubo=ubo  )
    !IF_NOTOK_RETURN(status=1)
    
    call self%GetValidPtr( args(i_ri), t,t, ri, status, check_units='1', &
                               check_lbo=lbo+(/0,0,1/), check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_halt), t,t, halt, status, check_units='m', &
                               check_lbo=lbo, check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! surface bound checks
    lbo_sfc = (/lbo(1),lbo(2), 1/)
    ubo_sfc = (/ubo(1),ubo(2), 1/)
    
    call self%GetValidPtr( args(i_monin), t,t, monin, status, check_units='m', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_blh), t,t, blh, status, check_units='m', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)
    call self%GetValidPtr( args(i_ustar), t,t, ustar, status, check_units='m/s', &
                               check_lbo=lbo_sfc, check_ubo=ubo_sfc )
    IF_NOTOK_RETURN(status=1)

    ! pointer to target field given the name:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m2/s')
    IF_NOTOK_RETURN(status=1)
    
    ! storage for target field:
    allocate( values(lbo(1):ubo(1),lbo(2):ubo(2),lbo(3):ubo(3)), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)       
        do l = lbo(3), ubo(3)
          if ( l == lbo(3) ) then
            values(i,j,l) = 0.0
          else 
            !call HalfLevelKzIFS( temper(i,j,l), dse(i,j,l), uu_ifs(i,j,l), duudz(i,j,l), ri(i,j,l), &
            !                     monin(i,j,1), halt(i,j,l), blh(i,j,1), ustar(i,j,1), values(i,j,l), &
            !                     status )
            !IF_NOTOK_RETURN(status=1)
            call HalfLevelKzIFS( temper(i,j,l), duudz(i,j,l), ri(i,j,l), &
                                 monin(i,j,1), halt(i,j,l), blh(i,j,1), ustar(i,j,1), values(i,j,l), &
                                 status )
            IF_NOTOK_RETURN(status=1)
          end if
        end do
      end do ! i
    end do ! j

    ! store:
    call iF%PutSample( values, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_HalfLevelKzIFS
  
  !
  ! HalfLevel Kz MSP (corrected kz value for large difference between level heigths
  !   output: kz
  !   input : h, monin, blh, ustar
  ! 

  subroutine Variables_Calc_HalfLevelKzMSP( self, varname, args, t, status )
  
    use GO          , only : TDate
    use GO          , only : T_Instant_Field
    use JAQL        , only : HalfLevelKzMSP
    
    ! --- in/out ---------------------------------
    
    class(T_Variables), intent(inout)             ::  self
    character(len=*), intent(in)                  ::  varname
    character(len=*), intent(in)                  ::  args(4)
    type(TDate), intent(in)                       ::  t
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Variables_Calc_HalfLevelKzMSP'
    
    ! argument names:   
    integer, parameter        ::  i_h         = 1
    integer, parameter        ::  i_monin_inv = 2
    integer, parameter        ::  i_blh       = 3
    integer, parameter        ::  i_ustar     = 4
    
    ! --- local ----------------------------------
    
    class(T_Instant_Field), pointer     ::  iF
    real, pointer                       ::  monin_inv(:,:,:)
    real, pointer                       ::  h(:,:,:)
    real, pointer                       ::  blh(:,:,:)
    real, pointer                       ::  ustar(:,:,:)
    integer                             ::  lbo(3), ubo(3)
    real, pointer                       ::  Kz(:,:,:)
    integer                             ::  i, j, l
      
    ! --- begin ----------------------------------
    
    ! pointer to target field given the name,
    ! check units, and return index space:
    call self%GetVariable( trim(varname), status, Instant_Field=iF, &
                             check_units='m2/s', lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    
    ! expected 3d dimension (0:nz)
    if ( lbo(3) /= 0 ) then
      write (gol,'("expected Kz defined on half levels with 3rd lower bound 0, found ",i0)') lbo(3); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! get pointer to target array (allocate if necessary), set target time:
    call iF%GetTargetPtr( Kz, t, status, lbo=lbo, ubo=ubo )
    IF_NOTOK_RETURN(status=1)

    ! pointer to input data:
    call self%GetValidPtr( args(i_h), t,t, h, status, &
                             check_units='m', &
                             check_lbo=lbo+(/0,0,1/), &
                             check_ubo=ubo )
    IF_NOTOK_RETURN(status=1)
    ! pointer to input data:
    call self%GetValidPtr( args(i_monin_inv), t,t, monin_inv, status, &
                             check_units='1/m', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)
    ! pointer to input data:
    call self%GetValidPtr( args(i_blh), t,t, blh, status, &
                             check_units='m', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)
    ! pointer to input data:
    call self%GetValidPtr( args(i_ustar), t,t, ustar, status, &
                             check_units='m/s', &
                             check_lbo=(/lbo(1),lbo(2),1/), &
                             check_ubo=(/ubo(1),ubo(2),1/) )
    IF_NOTOK_RETURN(status=1)

    ! loop over cells:
    do j = lbo(2), ubo(2)
      do i = lbo(1), ubo(1)       
        ! kz value equal to zero at surface
        Kz(i,j,0) = 0.0
        ! calculate dh corrected kz for other layers
        call HalfLevelKzMSP(  h(i,j,1:ubo(3)), monin_inv(i,j,1), blh(i,j,1), ustar(i,j,1), &
                              Kz(i,j,1:ubo(3)), status )
        IF_NOTOK_RETURN(status=1)
      end do ! i
    end do ! j
    
    ! ok
    status = 0
    
  end subroutine Variables_Calc_HalfLevelKzMSP
  
  
end module LE_Data_Variables
