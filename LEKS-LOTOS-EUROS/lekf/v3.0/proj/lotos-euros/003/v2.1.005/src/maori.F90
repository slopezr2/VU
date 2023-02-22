!#######################################################################
!
! NAME
!
!   MAORI - Model And Output Routine Interface
!
!
! DESCRIPTION
!
!   Mode independent output routines.
!
!
! USAGE
!
!   ! ##################################################################
!
!   program everything_you_never_wanted
!
!   use GO, only : TDate, NewDate
!   use MAORI
!
!   integer             ::  status
!   type(T_MAORI_Data)  ::  mdata
!   type(T_MAORI_State) ::  mstate
!   type(TDate)         ::  t
!   integer             ::  nset, iset
!   integer             ::  itype
!   real, allocatable   ::  values(:)
!
!   ! assimilation: observed variables:
!   real                ::  y, r, alfa
!   real                ::  Ht_xi(nmodes)
!
!   !
!   ! * init
!   !
!
!   ! read output settings:
!   call MAORI_Data_Init( mdata, 'settings.rc', status )
!   if (status/=0) stop
!
!   ! init state:
!   call MAORI_State_Init( mstate, mdata, status )
!   if (status/=0) stop
!
!   !
!   ! * provide information about the model
!   !
!
!   ! get number of sets to put out:
!   call MAORI_Data_Inquire( mdata, status, nset=nset )
!   if (status/=0) stop
!
!   ! loop over sets:
!   do iset = 1, nset
!
!     ! get number of requested parameters:
!     call MAORI_Data_Set_Inquire( mdata, iset, status, &
!                               name=sname, nparam=nparam )
!     if (status/=0) stop
!
!     ! info:
!     write (*,'("requested parameters for set ",i6," ",a," :")') iset, trim(sname)
!
!     ! loop over requested parameters:
!     do ipar = 1, npar
!
!       ! get parameter type, name and required unit:
!       call MAORI_Data_Param_Inquire( mdata, iset, ipar, &
!                                       itype, pname, punit, status )
!       if (status/=0) stop
!
!       ! info:
!       write (*,'("  ",i6," `",a,"` (",a,")") ipar, trim(pname), &
!                                trim(MAORI_TYPE_NAME(itype))
!
!       ! extract per type:
!       select case ( itype )
!         case ( MAORI_INT )
!           ! extract value from model data in requested unit:
!           ival = ...  pname, punit, ...
!           ! set parameter:
!           call MAORI_Data_Param_Put( mdata, iset, ipar, ival, status )
!           if (status/=0) stop
!         case ( MAORI_REAL )
!           ! extract value from model data in requested unit:
!           rval = ...  pname, punit, ...
!           ! set parameter:
!           call MAORI_Data_Param_Put( mdata, iset, ipar, rval, status )
!           if (status/=0) stop
!         case default
!           write (*,'("ERROR - unsupported maori type ",i6," for parameter ",a)')
!                        itype, trim(pname)
!           stop
!       end select
!
!     end do  ! parameters
!
!   end do  ! sets
!
!   ! start output:
!   call MAORI_Data_Start( mdata, status )
!   if (status/=0) stop
!
!   ! start output for state:
!   call MAORI_State_Start( mstate, mdata, status )
!   if (status/=0) stop
!
!   !
!   ! * loop
!   !
!
!   ! time loop:
!   do
!
!     ! ...
!
!     ! current time interval:
!     t1 = NewDate( 2001, 02, 03, 12, 00 )
!     t2 = NewDate( 2001, 02, 03, 13, 00 )
!
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ! setup for current time
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     ! setup for current time interval:
!     call MAORI_Data_Setup( mdata, t1, t2, status )
!     if (status/=0) stop
!
!     ! setup state for current time interval:
!     call MAORI_State_Setup( mstate, mdata, status )
!     if (status/=0) stop
!
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ! put simulations
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     ! get number of sets to put out:
!     call MAORI_Data_Inquire( mdata, status, nset=nset )
!     if (status/=0) stop
!
!     ! loop over sets:
!     do iset = 1, nset
!
!       ! get type of output:
!       call MAORI_Data_Set_Inquire( mdata, iset, status, type=stype )
!       if (status/=0) stop
!
!       ! do something given the type:
!       select case ( stype )
!
!         case ( MAORI_SAMPLE )
!
!           ! get number of locations and number of variables to be written:
!           call MAORI_Data_Inquire( mdata, iset, status, nloc=nloc, nvar=nvar )
!           if (status/=0) stop
!
!           ! loop over stations:
!           do iloc = 1, nloc
!
!             ! get location and altitude:
!             call MAORI_Data_Loc_Inquire( mdata, iset, iloc, status, &
!                                       lon=lon, lat=lat, alt=alt )
!             if (status/=0) stop
!
!             ! loop over variables:
!             do ivar = 1, nvar
!
!               ! get variable name and unit:
!               call MAORI_Data_Var_Inquire( mdata, iset, ivar, status, &
!                                          name=vname, unit=vunit, nlev=nlev )
!               if (status/=0) stop
!
!               ! storage:
!               allocate( values(1:nlev) )
!
!               ! extract values in requested unit on this location :
!               values(1:nlev) =   ! ... lon, lat, alt, vname, vunit, ..
!
!               ! put single layer or profile:
!               call MAORI_State_Put( mstate, mdata, iset, ivar, status, &
!                                                iloc=iloc, values=values )
!               if (status/=0) stop
!
!               ! clear:
!               deallocate( values )
!
!             end do   ! variables
!
!           end do  ! stations
!
!         case default
!           write (*,'("ERROR - unknown maori type ",i6," (`",a,"`) for set ",i6)') &
!                                    stype, trim(MAORI_TYPE_NAME(stype)), iset
!           stop
!
!       end select
!
!     end do   ! sets
!
!     ! all values are put; postprocess:
!     call MAORI_State_Post( mstate, mdata, status )
!     if (status/=0) stop
!
!     ! write if necessary:
!     call MAORI_State_Write( mstate, mdata, status )
!     if (status/=0) stop
!
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ! assimilation
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     ! get number of sets to put out:
!     call MAORI_Data_Inquire( mdata, status, nset=nset )
!     if (status/=0) stop
!
!     ! loop over sets:
!     do iset = 1, nset
!
!       ! extract flag:
!       call MAORI_Data_Set_Inquire( mdata, iset, status, assim_analyse_now=assim_analyse_now )
!       if (status/=0) stop
!       ! not necessary to continue ?
!       if ( .not. assim_analyse_now ) cycle
!
!       ! extract other assimilation parameters:
!       call MAORI_Data_Set_Inquire( mdata, iset, status, &
!                                      assim_rho=rho )
!       if (status/=0) stop
!
!       ! get type of output:
!       call MAORI_Data_Set_Inquire( mdata, iset, status, type=stype )
!       if (status/=0) stop
!
!       ! do something given the type:
!       select case ( stype )
!
!         case ( MAORI_SAMPLE )
!
!           ! get number of locations and number of variables to be written:
!           call MAORI_Data_Set_Inquire( mdata, iset, status, &
!                                              nloc=nloc, obs_nvar=obs_nvar )
!           if (status/=0) stop
!
!           ! loop over stations:
!           do iloc = 1, nloc
!
!             ! loop over variables:
!             do obs_ivar = 1, obs_nvar
!
!               ! get observed value, error std.dev., screening factor, and localisation parameter:
!               call MAORI_Data_Obs_Get( mdata, iset, iloc, obs_ivar, status, 
!                                               y=y, r=r, alfa=alfa )
!               if (status/=0) stop
!
!               ! given some tests values, set status bits;
!               ! see the section 'PARAMETERS' for possible bits:
!               if ( y < 0.0 ) then
!                 call MAORI_Data_Obs_Put( mdata, iset, iloc, obs_ivar, status, 
!                                                  astat_ibset=MAORI_ASTAT_NODATA )
!                 if (status/=0) stop
!               end if
!
!               ! extract simulations of the observed values from ensemble members:
!               do j = 1, nmodes
!                 call MAORI_State_Obs_Get( ens(j)%mstate, mdata, iset, obs_ivar, status, 
!                                                   iloc=iloc, value=Ht_xi(j) )
!                 if (status/=0) stop
!               end do
!
!               ! analysis with localisation:
!               call Analyse( y, r, Ht_xi, alfa, rho, status )
!               if (status/=0) stop
!
!             end do   ! variables
!
!           end do  ! stations
!
!         case default
!           write (*,'("ERROR - unknown maori type ",i6," (`",a,"`) for set ",i6)') &
!                                    stype, trim(MAORI_TYPE_NAME(stype)), iset
!           stop
!
!       end select
!
!     end do
!
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ! operations on state values
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     ! total number of elements in a state:
!     call MAORI_Data_Inquire( mdata, status, nvalue=nvalue )
!     if (status/=0) stop
!
!     ! storage:
!     allocate( values(nvalue) )
!
!     ! extract rank-1 array with state values:
!     call MAORI_State_Values_Get( mstate, mdata, values, status )
!     if (status/=0) stop
!
!     ! do something useful ...
!     values = values + 0.0
!
!     ! restore:
!     call MAORI_State_Values_Put( mstate, mdata, values, status )
!     if (status/=0) stop
!
!     ! clear:
!     deallocate( values )
!     
!
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ! other
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     ! ...
!
!   end do   ! time loop
!
!   ! done with state:
!   call MAORI_State_Done( mstate, mdata, status )
!   if (status/=0) stop
!
!   ! done with output:
!   call MAORI_Data_Done( mdata, status )
!   if (status/=0) stop
!
!   end program everything_you_never_wanted
!
!
!   ! ##################################################################
!
!
!   subroutine Analyse( lon0, lat0, y, r, Htx, alfa, rho, status )
!
!     real, intent(in)      ::  y, r, Htx, alfa, ro
!     integer, intent(out)  ::  status
!
!     ! get number of sets to put out:
!     call MAORI_Data_Inquire( mdata, status, nset=nset )
!     if (status/=0) stop
!
!     ! loop over sets:
!     do iset = 1, nset
!
!       ! get number of locations and number of state values per location:
!       call MAORI_Data_Set_Inquire( mdata, iset, status, nloc=nloc, nlocvalue=nlocvalue )
!       if (status/=0) stop
!
!       ! storage:
!       allocate( xi_maori(nlocvalue,nmodes) )
!       allocate( Sf_maori(nlocvalue,nmodes) )
!       allocate( xf_maori(nlocvalue)        )
!       allocate(  K_maori(nlocvalue)        )
!
!       ! loop over locations:
!       do iloc = 1, nloc
!
!         ! get location:
!         call MAORI_Data_Loc_Inquire( mdata, iset, iloc, status, lon=lon, lat=lat )
!         if (status/=0) stop
!
!         ! distance in m:
!         call Distance( lon0, lat0, lon, lat, dist, status )
!         if (status/=0) stop
!
!         ! outside localisation range ? then skip:
!         if ( dist > 3.5 * rho ) cycle
!         ! ... else compute correlation factor:
!         fcc = exp( -(dist/rho)**2 )
!
!         ! loop over ensemble members::
!         do j = 1, nmodes
!           ! extract state on current location:
!           call MAORI_State_LocValues_Get( ens(j)%mstate, mdata, iset, iloc, xi_maori(:,j), status )
!           if (status/=0) stop
!         end do
!
!         ! ensemble mean:
!         xf_maori = sum(Sf_maori,2) / real(nmodes)
!         ! covariance square root:
!         do j = 1, nmodes
!           Sf_maori(:,j) = ( xi_maori(:,j) - xf_maori ) / sqrt(nmodes-1.0)
!         end do
!
!         ! gain for this location: K = S(:,1:nmodes) * Theta(1:nmodes,ip)
!         K_maori = 0.0
!         do j = 1, nmodes
!           K_maoir = K_maori + Sf_maori(:,j) * Theta(j)
!         end do
!         ! bound gain:
!         K_maori = K_maori * fcc
!       
!         ! loop over ensemble members:
!         do j = 1, nmodes
!           ! analyse ensemble member:
!           xi_maori(:,j) = xi_maori(:,j) + K_maori * dd(j)
!         end do
!       
!         ! loop over ensemble members:
!         do j = 1, nmodes
!           ! restore ensemble members:
!           call MAORI_State_LocValues_Put( ens(j)%mstate, mdata, iset, iloc, xi_maori(:,j), status )
!           if (status/=0) stop
!         end do
!
!       end do  ! locations
!
!       ! clear:
!       deallocate( xi_maori )
!       deallocate( Sf_maori )
!       deallocate( xf_maori )
!       deallocate(  K_maori )
!
!     end do  ! sets
!
!     ! ok
!     status = 0
!
!   end subroutine Analyse
!     
!
!   ! ##################################################################
!
!
! PARAMETERS
!
!   Assimilation status parameters:
!    MAORI_ASTAT_OUTSIDE       ! location outside domain
!    MAORI_ASTAT_NODATA     
!    MAORI_ASTAT_VALIDATION 
!    MAORI_ASTAT_SCREENED   
!    MAORI_ASTAT_ANALYSED   
!
!
! RCFILE
!
!
!    ! list with moari sets to put out (space seperated, empty for none):
!    !maori.sets                        :  
!    maori.sets                        :  emep emep24 earlinet aeronet
!
!    ! general stuff:
!    maori.path                        :  ${OUTPUTDIR}/${PROJECT}/${RUNID}
!    maori.author.name                 :  The LOTOS-EUROS team.
!    maori.author.inst                 :  TNO, RIVM, MNP, KNMI
!    maori.model.name                  :  LE
!    maori.model.version               :  ${VERSION_NR}
!    maori.experiment.id               :  ${RUNID}
!    maori.data.version                :  -
!
!    ! settings for ground stations:
!    !
!    !   type      :  sample
!    !
!    !   tres      :
!    !                step=1.0;unit=hour                ! hourly average
!    !                step=1.0;unit=hour;hist=T         ! hourly average, store history during day
!    !                step=24.0;unit=hour;start=0.0     ! daily average
!    !
!    !   loc.file  :
!    !                type=lml;file=/data/input/station.list
!    !                type=emep;file=/modas/shared/observations/EMEP/stations.txt
!    !                type=earlinet;file=${PROJECTS}/GEOmon/data/EARLINET_stations.csv
!    !                type=aeronet;file=${PROJECTS}/SATELLIETEN/AERONET/aeronet_locations_europe_2006.txt
!    !
!    !   observations:
!    !     obs.type  : emep-daily-csv        ! EMEP dail averages, csv file by HvdB
!    !       obs.var       :  PM10      
!    !       obs.file      :  /data/emep_%{var}_${year4}.csv
!    !       obs.err.frac  :  0.20           ! assumed error std.dev. is fraction of value
!    !
!    !   assimilation:
!    !     assim.analyse   :  T     ! (T|F)
!    !     assim.alfa      :  5.0   ! screening factor
!    !     assim.rho       :  50.0  ! localisation scale (km)
!    !
!    ! aerosol variables:
!    #eval var_aero="tpm10 tpm25 bc pm25 pm10 so4a no3a nh4a dust_f dust_c na_f na_c aod"
!    #eval var_profile="temper rhumid halflevel_altitude"
!    !
!    ! EMEP : hourly samples
!    maori.emep.type                 :  sample
!    maori.emep.loc.file             :  type=emep;file=/modas/shared/observations/EMEP/stations.txt
!    maori.emep.tres                 :  step=1.0;unit=hour
!    maori.emep.profile              :  F
!    maori.emep.var                  :  ${var_aero}
!    maori.emep.obs.var              :  
!    maori.emep.obs.type             :  
!    maori.emep.obs.file             :   
!    maori.emep.obs.err.frac         :  
!    maori.emep.assim.analyse        :  F
!    maori.emep.assim.alfa           :  
!    maori.emep.assim.rho            : 
!    !
!    ! EMEP24 : daily total pm10
!    maori.emep24.type               :  sample
!    maori.emep24.loc.file           :  type=emep;file=/modas/shared/observations/EMEP/stations.txt
!    maori.emep24.tres               :  step=24.0;unit=hour
!    maori.emep24.profile            :  F
!    maori.emep24.var                :  ${var_aero}
!    maori.emep24.obs.type           :  emep-daily-csv
!    maori.emep24.obs.var            :  PM10
!    maori.emep24.obs.file           :  ${PROJECTS}/GEOmon/data/Stations_DailyTimeSeries_%{var}_%{year4}.csv
!    maori.emep24.obs.err.frac       :  0.20
!    maori.emep24.assim.analyse      :  T
!    maori.emep24.assim.alfa         :  5.0
!    maori.emep24.assim.rho          :  50.0
!    !
!    ! EARLINET : aerosol lidars
!    maori.earlinet.type             :  sample
!    maori.earlinet.loc.file         :  type=earlinet;file=${PROJECTS}/GEOmon/data/EARLINET_stations.csv
!    maori.earlinet.tres             :  step=1.0;unit=hour
!    maori.earlinet.profile          :  T
!    maori.earlinet.var              :  ${var_aero} ${var_profile}
!    maori.earlinet.obs.var          :  
!    maori.earlinet.obs.type         :  
!    maori.earlinet.obs.file         :  
!    !
!    ! AERONET : AOD profiles
!    maori.aeronet.type              :  sample
!    maori.aeronet.loc.file          :  type=aeronet;file=${PROJECTS}/SATELLIETEN/AERONET/aeronet_locations_europe_2006.txt
!    maori.aeronet.tres              :  step=1.0;unit=hour
!    maori.aeronet.profile           :  T
!    maori.aeronet.var               :  ${var_aero} ${var_profile}
!    maori.aeronet.obs.var           :  
!    maori.aeronet.obs.type          :  
!    maori.aeronet.obs.file          :  
!
!
! IMPLEMENTATION
!
!   Module tree:
!
!     MAORI
!       MAORI_Set
!         MAORI_Set_Sample
!         MAORI_Set_Satellite
!
!   The main data types have names based on the module name:
!
!     T_MAORI[_Set[_Sample]]_Data
!     T_MAORI[_Set[_Sample]]_State
!     T_MAORI[_Set[_Sample]]_Output
!
!   Idem for the routine names:
!
!     MAORI[_Set[_Sample]]_Data_Init
!     MAORI[_Set[_Sample]]_Data_Start
!     MAORI[_Set[_Sample]]_Data_Setup
!     MAORI[_Set[_Sample]]_Data_Done
!     MAORI[_Set[_Sample]]_Data_Inquire
!     MAORI               _Data_Set_Inquire
!     MAORI[_Set[_Sample]]_Data_Param_Inquire
!     MAORI[_Set[_Sample]]_Data_Param_Put
!     MAORI[_Set[_Sample]]_Data_Var_Inquire
!     MAORI[_Set[_Sample]]_Data_Loc_Inquire
!     MAORI[_Set[_Sample]]_Data_Obs_Get
!     MAORI[_Set[_Sample]]_Data_Obs_Put
!
!     T_MAORI_State
!     MAORI_State_Init
!     MAORI_State_Start
!     MAORI_State_Setup
!     MAORI_State_Done
!     MAORI_State_Put
!     MAORI_State_Post
!     MAORI_State_Obs_Get
!     MAORI_State_Values_Get
!     MAORI_State_Values_Put
!     MAORI_State_LocValues_Get
!     MAORI_State_LocValues_Put
!
!     T_MAORI_Output
!     MAORI_Output_Init
!     MAORI_Output_Start
!     MAORI_Output_Done
!     MAORI_Output_Write
!
!
! HISTORY
!
!   2009-2011, Arjo Segers, TNO
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  use MAORI_Param
  use MAORI_Set, only : T_MAORI_Set_Data, T_MAORI_Set_State, T_MAORI_Set_Output

  implicit none


  ! --- in/out ---------------------------

  private

  public  ::  MAORI_TYPE_NAME, MAORI_INT, MAORI_REAL, MAORI_CHAR
  public  ::  MAORI_SAMPLE, MAORI_SATELLITE
  public  ::  MAORI_LEN_NAME, MAORI_LEN_UNIT, MAORI_LEN_LINE
  
  public  ::  MAORI_ASTAT_OUTSIDE, MAORI_ASTAT_NODATA, MAORI_ASTAT_VALIDATION, &
                MAORI_ASTAT_SCREENED, MAORI_ASTAT_ANALYSED
  
  public  ::  T_MAORI_Data
  public  ::  MAORI_Data_Init, MAORI_Data_Start,  MAORI_Data_Setup, MAORI_Data_Done
  public  ::  MAORI_Data_Inquire
  public  ::  MAORI_Data_Set_Inquire
  public  ::  MAORI_Data_Param_Inquire, MAORI_Data_Param_Put
  public  ::  MAORI_Data_Var_Inquire
  public  ::  MAORI_Data_Loc_Inquire
  public  ::  MAORI_Data_Obs_Get, MAORI_Data_Obs_Put, MAORI_Data_ObsValid_Get

  public  ::  T_MAORI_State
  public  ::  MAORI_State_Init, MAORI_State_Start,  MAORI_State_Setup, MAORI_State_Done
  public  ::  MAORI_State_Put, MAORI_State_Post
  public  ::  MAORI_State_Obs_Get
  public  ::  MAORI_State_Values_Get, MAORI_State_Values_Put
  public  ::  MAORI_State_LocValues_Get, MAORI_State_LocValues_Put
  public  ::  MAORI_State_ObsValid_Get
  
  public  ::  T_MAORI_Output
  public  ::  MAORI_Output_Init, MAORI_Output_Start, MAORI_Output_Done
  public  ::  MAORI_Output_Write



  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI'


  ! --- types ------------------------------

  type T_MAORI_Data
    ! number of data sets:
    integer                             ::  nset
    ! names:
    character(len=MAORI_LEN_NAME)       ::  setname(MAORI_MAX_SET)
    ! output sets:
    type(T_MAORI_Set_Data), pointer     ::  set(:)
    ! flags:
    logical                             ::  started
    ! current time interval:
    type(TDate)                         ::  t1, t2
  end type T_MAORI_Data

  type T_MAORI_State
    ! name:
    character(len=MAORI_LEN_NAME)       ::  name
    ! number of state sets:
    integer                             ::  nset
    ! output states:
    type(T_MAORI_Set_State), pointer    ::  set(:)
  end type T_MAORI_State

  type T_MAORI_Output
    ! name:
    character(len=MAORI_LEN_NAME)       ::  name
    ! number of output sets:
    integer                             ::  nset
    ! output states:
    type(T_MAORI_Set_Output), pointer   ::  set(:)
  end type T_MAORI_Output


  ! --- interfaces -------------------------------
  
  interface MAORI_Data_Param_Put
    module procedure MAORI_Data_Param_Put_i
    module procedure MAORI_Data_Param_Put_r
    module procedure MAORI_Data_Param_Put_s
  end interface


contains


  ! ********************************************************************
  ! ***
  ! *** main routines
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Init( mdata, rcfile, t, status )

    use GO, only : TRcFile, Init, Done, ReadRc
    use GO, only : goSplitString
    use GO, only : TDate, NewDate
    use MAORI_Param , only : MAORI_LEN_LINE
    use MAORI_Common, only : MAORI_Common_Init
    use MAORI_Set   , only : MAORI_Set_Data_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(out)       ::  mdata
    character(len=*), intent(in)          ::  rcfile
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Init'

    ! --- local -----------------------------------

    type(TRcFile)                         ::  rcF
    character(len=MAORI_LEN_LINE)         ::  line
    integer                               ::  iset

    ! --- begin -----------------------------------

    ! open rcfile:
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read set names:
    call ReadRc( rcF, 'maori.sets', line, status )
    IF_NOTOK_RETURN(status=1)

    ! extract list of set names:
    call goSplitString( line, mdata%nset, mdata%setname, status )
    IF_NOTOK_RETURN(status=1)

    ! any request ?
    if ( mdata%nset > 0 ) then

      ! init common variables:
      call MAORI_Common_Init( rcF, status )
      IF_NOTOK_RETURN(status=1)

      ! storage:
      allocate( mdata%set(1:mdata%nset) )

      ! loop over sets:
      do iset = 1, mdata%nset
        ! init set:
        call MAORI_Set_Data_Init( mdata%set(iset), rcF, mdata%setname(iset), t, status )
        IF_NOTOK_RETURN(status=1)
      end do  ! sets

    end if

    ! close:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)

    ! not started yet:
    mdata%started = .false.
    
    ! no setup yet, set dummy times:
    mdata%t1 = NewDate(0001,1,1)
    mdata%t2 = NewDate(0001,1,1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Init


  ! ***


  subroutine MAORI_Data_Start( mdata, status )

    use MAORI_Set, only : MAORI_Set_Data_Start

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)      ::  mdata
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Start'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    write (gol,'("MAORI:  start data ...")'); call goPr

    ! loop over sets:
    do iset = 1, mdata%nset
      ! init set:
      call MAORI_Set_Data_Start( mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! reset flag:
    mdata%started = .true.

    ! ok
    status = 0

  end subroutine MAORI_Data_Start


  ! ***


  subroutine MAORI_Data_Done( mdata, status )

    use MAORI_Common, only : MAORI_Common_Done
    use MAORI_Set   , only : MAORI_Set_Data_Done
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)      ::  mdata
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Done'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! init set:
      call MAORI_Set_Data_Done( mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! done with common variables:
    call MAORI_Common_Done( status )
    IF_NOTOK_RETURN(status=1)

    ! reset flag:
    mdata%started = .false.

    ! ok
    status = 0

  end subroutine MAORI_Data_Done


  ! ***


  !
  !  call MAORI_Data_Check( mdata, status, iset=12 )
  !
  !    ! check if iset is in valid range
  !    ! check if type of the set corresponds with allocated data
  !
  !  call MAORI_Data_Check( mdata, status, iset=12, type=MAORI_SAMPLE )
  !
  !    ! check iset as described above
  !    ! check if type is in valid range
  !    ! check if the set has the specified type
  !

  subroutine MAORI_Data_Check( mdata, status, iset, type )

    use MAORI_Param, only : MAORI_TYPE_MIN, MAORI_TYPE_MAX, MAORI_TYPE_NAME
    use MAORI_Param, only : MAORI_SAMPLE
    use MAORI_Set  , only : MAORI_Set_Data_Check

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(in)         ::  mdata
    integer, intent(out)              ::  status

    integer, intent(in), optional     ::  type
    integer, intent(in), optional     ::  iset

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Check'

    ! --- begin -----------------------------------

    ! check set index:
    if ( present(iset) ) then

      ! in range ?
      if ( (iset < 1) .or. (iset > mdata%nset) ) then
        write (gol,'("iset = ",i6,", should be in 1 .. ",i6)') mdata%nset; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! set specific tests:
      call MAORI_Set_Data_Check( mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)

    end if

    ! check if set has specified type:
    if ( present(type) ) then

      ! requires iset ..
      if ( .not. present(iset) ) then
        write (gol,'("argument type requires iset")'); call goErr
        TRACEBACK; status=1; return
      end if

      ! check range:
      if ( (type < MAORI_TYPE_MIN) .or. (type > MAORI_TYPE_MAX) ) then
        write (gol,'("type id ",i6," out of range ",i6," .. ",i6)') &
                        type, MAORI_TYPE_MIN, MAORI_TYPE_MAX; call goErr
        TRACEBACK; status=1; return
      end if

      ! check type:
      if ( mdata%set(iset)%type /= type ) then
        write (gol,'("MAORI Type mismatch :")'); call goErr
        write (gol,'("  requested  : ",i6," (",a,")")') &
                           type, trim(MAORI_TYPE_NAME(type))
        write (gol,'("  found      : ",i6," (",a,")")') &
                           mdata%set(iset)%type, trim(MAORI_TYPE_NAME(mdata%set(iset)%type))
        TRACEBACK; status=1; return
      end if

    end if  ! itype

    ! ok
    status = 0

  end subroutine MAORI_Data_Check


  ! ***


  subroutine MAORI_Data_Inquire( mdata, status, nset, nobsvalue, nobsvalid, nvalue )

    use MAORI_Set, only : MAORI_Set_Data_Inq
    use MAORI_Set, only : MAORI_Set_Data_Obs_Get
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(in)          ::  mdata
    integer, intent(out)                    ::  status
    
    integer, intent(out), optional          ::  nset
    integer, intent(out), optional          ::  nobsvalue    ! nr. of observations
    integer, intent(out), optional          ::  nobsvalid    ! nr. of observations >= 0.0
    integer, intent(out), optional          ::  nvalue       ! nr. state values (observations plus history etc)

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Inquire'

    ! --- local -----------------------------------
    
    integer     ::  iset
    integer     ::  stype
    integer     ::  nval
    integer     ::  nloc, iloc
    integer     ::  obs_nvar, obs_ivar
    real        ::  y

    ! --- begin -----------------------------------

    ! return number of output sets:
    if ( present(nset) ) nset = mdata%nset
    
    ! number of observation values:
    if ( present(nobsvalue) ) then
      ! nothing yet:
      nobsvalue = 0
      ! loop over sets:
      do iset = 1, mdata%nset
        ! inquire number of values in set:
        call MAORI_Set_Data_Inq( mdata%set(iset), status, nloc=nloc, obs_nvar=obs_nvar )
        IF_NOTOK_RETURN(status=1)
        ! add contribution:
        nobsvalue = nobsvalue + nloc * obs_nvar
      end do   ! sets
    end if
    
    ! number of observation values:
    if ( present(nobsvalid) ) then
      ! nothing yet:
      nobsvalid = 0
      ! loop over sets:
      do iset = 1, mdata%nset
        ! get type of output:
        call MAORI_Set_Data_Inq( mdata%set(iset), status, type=stype )
        IF_NOTOK_RETURN(status=1)
        ! check ...
        if ( stype /= MAORI_SAMPLE ) stop 'MAORI_Data_Inquire only implemented for samples yet'
        ! inquire number of values in set:
        call MAORI_Set_Data_Inq( mdata%set(iset), status, nloc=nloc, obs_nvar=obs_nvar )
        IF_NOTOK_RETURN(status=1)
        ! loop over stations:
        do iloc = 1, nloc
          ! loop over variables:
          do obs_ivar = 1, obs_nvar
            ! extract measurement:
            call MAORI_Set_Data_Obs_Get( mdata%set(iset), obs_ivar, status, y=y, iloc=iloc )
            IF_NOTOK_RETURN(status=1)
            ! valid data ?
            if ( y >= 0.0 ) nobsvalid = nobsvalid + 1
          end do  ! observed variables
        end do  ! station locations
      end do   ! sets
    end if
    
    ! number of state values used in assimilation:
    if ( present(nvalue) ) then
      ! nothing yet:
      nvalue = 0
      ! loop over sets:
      do iset = 1, mdata%nset
        ! inquire number of values in set:
        call MAORI_Set_Data_Inq( mdata%set(iset), status, nvalue=nval )
        IF_NOTOK_RETURN(status=1)
        ! add contribution:
        nvalue = nvalue + nval
      end do   ! sets
    end if

    ! ok
    status = 0

  end subroutine MAORI_Data_Inquire


  ! ***


  subroutine MAORI_Data_Setup( mdata, t1, t2, the_end, status )

    use GO       , only : TDate, IsAnyDate, wrtgol, operator(/=)
    use MAORI_Set, only : MAORI_Set_Data_Setup

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)       ::  mdata
    type(TDate), intent(in)                 ::  t1, t2
    logical, intent(in)                     ::  the_end
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Data_Setup'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    ! different ?
    if ( (t1 /= mdata%t1) .or. (t2 /= mdata%t2) ) then
      ! info ...
      call wrtgol( 'MAORI Setup for ', t1, ' - ', t2 ); call goPr
      ! loop over sets:
      do iset = 1, mdata%nset
        ! init set:
        call MAORI_Set_Data_Setup( mdata%set(iset), t1, t2, the_end, status )
        IF_NOTOK_RETURN(status=1)
      end do
    else
      ! info ..
      call wrtgol( 'MAORI Setup already done for ', t1, ' - ', t2 ); call goPr
    end if

    ! store current time interval:
    mdata%t1 = t1
    mdata%t2 = t2

    ! ok
    status = 0

  end subroutine MAORI_Data_Setup


  ! ***


  subroutine MAORI_Data_Set_Inquire( mdata, iset, status, &
                                  name, type, &
                                  nparam, nloc, loc_id_range, nlon, nlat, &
                                  nvar, obs_nvar, &
                                  assim_analyse, assim_analyse_now, assim_rho, &
                                  nlocvalue, varname )

    use MAORI_Set, only : MAORI_Set_Data_Inq

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(in)            ::  mdata
    integer, intent(in)                       ::  iset
    integer, intent(out)                      ::  status

    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  type
    integer, intent(out), optional            ::  nparam
    integer, intent(out), optional            ::  nloc
    integer, intent(out), optional            ::  loc_id_range(2)
    integer, intent(out), optional            ::  nlon, nlat
    integer, intent(out), optional            ::  nvar
    integer, intent(out), optional            ::  obs_nvar
    logical, intent(out), optional            ::  assim_analyse
    logical, intent(out), optional            ::  assim_analyse_now
    real, intent(out), optional               ::  assim_rho
    integer, intent(out), optional            ::  nlocvalue
    character(len=*), intent(out), optional   ::  varname

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Set_Inquire'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! inquire output set:
    call MAORI_Set_Data_Inq( mdata%set(iset), status, &
                          name=name, type=type, &
                          nparam=nparam, nloc=nloc, loc_id_range=loc_id_range, nlon=nlon, nlat=nlat, &
                          nvar=nvar, obs_nvar=obs_nvar, &
                          assim_analyse=assim_analyse, assim_analyse_now=assim_analyse_now, &
                          assim_rho=assim_rho, &
                          nlocvalue=nlocvalue, varname=varname )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Set_Inquire


  ! ********************************************************************
  ! ***
  ! *** parameters
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Param_Inquire( mdata, iset, ipar, type, name, unit, status )

    use MAORI_Param, only : MAORI_SAMPLE
    use MAORI_Set  , only : MAORI_Set_Data_Param_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)         ::  mdata
    integer, intent(in)                       ::  iset
    integer, intent(in)                       ::  ipar
    integer, intent(out)                      ::  type
    character(len=*), intent(out)             ::  name
    character(len=*), intent(out)             ::  unit
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Param_Inquire'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! inquire output set:
    call MAORI_Set_Data_Param_Inq( mdata%set(iset), ipar, type, name, unit, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Param_Inquire


  ! ***
  
  
  subroutine MAORI_Data_Param_Put_i( mdata, iset, ipar, ival, status )

    use MAORI_Param, only : MAORI_SAMPLE
    use MAORI_Set  , only : MAORI_Set_Data_Param_Put
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)         ::  mdata
    integer, intent(in)                       ::  iset
    integer, intent(in)                       ::  ipar
    integer, intent(in)                       ::  ival
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Param_Put_r'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! put parameter to output set:
    call MAORI_Set_Data_Param_Put( mdata%set(iset), ipar, ival, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Param_Put_i


  ! ***
  
  
  subroutine MAORI_Data_Param_Put_r( mdata, iset, ipar, rval, status )

    use MAORI_Param, only : MAORI_SAMPLE
    use MAORI_Set  , only : MAORI_Set_Data_Param_Put
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)         ::  mdata
    integer, intent(in)                       ::  iset
    integer, intent(in)                       ::  ipar
    real, intent(in)                          ::  rval
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Param_Put_r'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! put parameter to output set:
    call MAORI_Set_Data_Param_Put( mdata%set(iset), ipar, rval, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Param_Put_r


  ! ***
  
  
  subroutine MAORI_Data_Param_Put_s( mdata, iset, ipar, sval, status )

    use MAORI_Param, only : MAORI_SAMPLE
    use MAORI_Set  , only : MAORI_Set_Data_Param_Put
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)         ::  mdata
    integer, intent(in)                       ::  iset
    integer, intent(in)                       ::  ipar
    character(len=*), intent(in)              ::  sval
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Param_Put_s'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! put parameter to output set:
    call MAORI_Set_Data_Param_Put( mdata%set(iset), ipar, sval, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Param_Put_s


  ! ********************************************************************
  ! ***
  ! *** locations
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Loc_Inquire( mdata, iset, iloc, status, &
                                      loc_id, name, code, lon, lat, alt, &
                                      varname )

    use MAORI_Set  , only : MAORI_Set_Data_Loc_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(in)              ::  mdata
    integer, intent(in)                         ::  iset
    integer, intent(in)                         ::  iloc
    integer, intent(out)                        ::  status

    integer, intent(out), optional              ::  loc_id
    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  code
    real, intent(out), optional                 ::  lon
    real, intent(out), optional                 ::  lat
    real, intent(out), optional                 ::  alt
    character(len=*), intent(out), optional     ::  varname

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Loc_Inquire'

    ! --- begin ----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! inquire per type:
    call MAORI_Set_Data_Loc_Inq( mdata%set(iset), iloc, status, &
                                   loc_id=loc_id, name=name, code=code, lon=lon, lat=lat, alt=alt, &
                                   varname=varname )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Loc_Inquire


  ! ********************************************************************
  ! ***
  ! *** variables
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Var_Inquire( mdata, iset, ivar, status, &
                                          name, unit, nlev )

    use MAORI_Set  , only : MAORI_Set_Data_Var_Inq
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(in)                   ::  mdata
    integer, intent(in)                         ::  iset
    integer, intent(in)                         ::  ivar
    integer, intent(out)                        ::  status
    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  unit
    integer, intent(out), optional              ::  nlev

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Var_Inquire'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! inquire per type:
    call MAORI_Set_Data_Var_Inq( mdata%set(iset), ivar, status, &
                                 name=name, unit=unit, nlev=nlev )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Var_Inquire


  ! ********************************************************************
  ! ***
  ! *** observations
  ! ***
  ! ********************************************************************


  subroutine MAORI_Data_Obs_Get( mdata, iset, obs_ivar, status, iloc, ilon, ilat, y, r, alfa )

    use MAORI_Set, only : MAORI_Set_Data_Obs_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(in)        ::  mdata
    integer, intent(in)                   ::  iset
    integer, intent(in)                   ::  obs_ivar
    integer, intent(out)                  ::  status
    
    integer, intent(in), optional         ::  iloc
    integer, intent(in), optional         ::  ilon, ilat
    real, intent(out), optional           ::  y       ! measured value
    real, intent(out), optional           ::  r       ! error std.dev.
    real, intent(out), optional           ::  alfa    ! screening factor

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Obs_Get'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! get observation values:
    call MAORI_Set_Data_Obs_Get( mdata%set(iset), obs_ivar, status, &
                                       iloc=iloc, ilon=ilon, ilat=ilat, y=y, r=r, alfa=alfa )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Obs_Get


  ! ***


  subroutine MAORI_Data_ObsValid_Get( mdata, status, y, r )

    use MAORI_Set, only : MAORI_Set_Data_Inq
    use MAORI_Set, only : MAORI_Set_Data_Obs_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(in)        ::  mdata
    integer, intent(out)                  ::  status
    real, intent(out), optional           ::  y(:)
    real, intent(out), optional           ::  r(:)

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_ObsValid_Get'

    ! --- local -----------------------------------
    
    integer     ::  iset
    integer     ::  stype
    integer     ::  i
    integer     ::  nloc, iloc
    integer     ::  obs_nvar, obs_ivar
    real        ::  y_curr, r_curr

    ! --- begin -----------------------------------
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! get type of output:
      call MAORI_Set_Data_Inq( mdata%set(iset), status, type=stype )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( stype /= MAORI_SAMPLE ) stop 'MAORI_Data_ObsValid_Get only implemented for samples yet'
      ! inquire number of locations and observed values in set:
      call MAORI_Set_Data_Inq( mdata%set(iset), status, nloc=nloc, obs_nvar=obs_nvar )
      IF_NOTOK_RETURN(status=1)
      ! loop over stations:
      do iloc = 1, nloc
        ! loop over variables:
        do obs_ivar = 1, obs_nvar
          ! extract measurement:
          call MAORI_Set_Data_Obs_Get( mdata%set(iset), obs_ivar, status, iloc=iloc, y=y_curr, r=r_curr )
          IF_NOTOK_RETURN(status=1)
          ! valid data ?
          if ( y_curr >= 0.0 ) then
            ! increase counter:
            i = i + 1
            ! store observation:
            if ( present(y) ) then
              ! storage should be large enough ...
              if ( i > size(y) ) then
                write (gol,'("size of output array y (",i6,") too small for number of values (",i6,") up to set ",i4)') &
                                     size(y), i, iset; call goErr
                TRACEBACK; status=1; return
              end if
              ! fill:
              y(i) = y_curr
            end if
            ! store observation error std.dev.:
            if ( present(r) ) then
              ! storage should be large enough ...
              if ( i > size(r) ) then
                write (gol,'("size of output array r (",i6,") too small for number of values (",i6,") up to set ",i4)') &
                                     size(r), i, iset; call goErr
                TRACEBACK; status=1; return
              end if
              ! fill:
              r(i) = r_curr
            end if
          end if  ! valid data
        end do  ! observed variables
      end do  ! station locations
    end do   ! sets

    ! ok
    status = 0

  end subroutine MAORI_Data_ObsValid_Get


  ! ***
  
  
  subroutine MAORI_Data_Obs_Put( mdata, iset, obs_ivar, status, iloc, ilon, ilat, astat_ibset )

    use MAORI_Set, only : MAORI_Set_Data_Obs_Put

    ! --- in/out ---------------------------------

    type(T_MAORI_Data), intent(inout)     ::  mdata
    integer, intent(in)                   ::  iset
    integer, intent(in)                   ::  obs_ivar
    integer, intent(out)                  ::  status
    
    integer, intent(in), optional         ::  iloc
    integer, intent(in), optional         ::  ilon, ilat
    integer, intent(in), optional         ::  astat_ibset

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Data_Obs_Put'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! put observation entities:
    call MAORI_Set_Data_Obs_Put( mdata%set(iset), obs_ivar, status, &
                                       iloc=iloc, ilon=ilon, ilat=ilat, astat_ibset=astat_ibset )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_Data_Obs_Put


  ! ********************************************************************
  ! ***
  ! *** state routines
  ! ***
  ! ********************************************************************


  subroutine MAORI_State_Init( mstate, mdata, name, status )

    use MAORI_Set, only : MAORI_Set_State_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(out)  ::  mstate
    type(T_MAORI_Data), intent(in)         ::  mdata
    character(len=*), intent(in)      ::  name
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_State_Init'

    ! --- local -----------------------------------

    integer       ::  iset

    ! --- begin -----------------------------------
    
    write (gol,'("MAORI: init state `",a,"` ...")') trim(name); call goPr
    
    ! store:
    mstate%name = trim(name)

    ! number of requested sets:
    mstate%nset = mdata%nset

    ! any request ?
    if ( mstate%nset > 0 ) then

      ! storage:
      allocate( mstate%set(1:mstate%nset) )

      ! loop over sets:
      do iset = 1, mstate%nset
        ! init set:
        call MAORI_Set_State_Init( mstate%set(iset), mdata%set(iset), trim(name), status )
        IF_NOTOK_RETURN(status=1)
      end do  ! sets

    end if

    ! ok
    status = 0

  end subroutine MAORI_State_Init


  ! ***


  subroutine MAORI_State_Start( mstate, mdata, status )

    use MAORI_Set, only : MAORI_Set_State_Start

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(inout)  ::  mstate
    type(T_MAORI_Data), intent(in)           ::  mdata
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Startup'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! init set:
      call MAORI_Set_State_Start( mstate%set(iset), mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! ok
    status = 0

  end subroutine MAORI_State_Start


  ! ***


  subroutine MAORI_State_Done( mstate, mdata, status )

    use MAORI_Set, only : MAORI_Set_State_Done
    
    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(inout)      ::  mstate
    type(T_MAORI_Data), intent(in)               ::  mdata
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_State_Done'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! init set:
      call MAORI_Set_State_Done( mstate%set(iset), mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! ok
    status = 0

  end subroutine MAORI_State_Done


  ! ***


  subroutine MAORI_State_Setup( mstate, mdata, status )

    use MAORI_Set, only : MAORI_Set_State_Setup

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(inout)      ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Setup'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! init set:
      call MAORI_Set_State_Setup( mstate%set(iset), mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! ok
    status = 0

  end subroutine MAORI_State_Setup


  ! ***


  subroutine MAORI_State_Put( mstate, mdata, iset, ivar, status, ilon, ilat, iloc, values )

    use MAORI_Param, only : MAORI_SAMPLE
    use MAORI_Set  , only : MAORI_Set_State_Put

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(inout)    ::  mstate
    type(T_MAORI_Data), intent(in)        ::  mdata
    integer, intent(in)                   ::  iset
    integer, intent(in)                   ::  ivar
    integer, intent(out)                  ::  status
    
    integer, intent(in), optional         ::  ilon, ilat, iloc
    real, intent(in), optional            ::  values(:)

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_State_Put'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! put simulated value:
    call MAORI_Set_State_Put( mstate%set(iset), mdata%set(iset), ivar, status, &
                                    ilon=ilon, ilat=ilat, iloc=iloc, values=values )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_State_Put


  ! ***

  
  subroutine MAORI_State_Post( mstate, mdata, status )

    use MAORI_Param, only : MAORI_SAMPLE
    use MAORI_Set  , only : MAORI_Set_State_Post

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(inout)      ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_State_Post'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! post processing after put:
      call MAORI_Set_State_Post( mstate%set(iset), mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! ok
    status = 0

  end subroutine MAORI_State_Post

  
  ! ***
  
  
  ! simulation of the observed value in the state:
  
  subroutine MAORI_State_Obs_Get( mstate, mdata, iset, obs_ivar, status, iloc, ilon, ilat, value )

    use MAORI_Param, only : MAORI_SAMPLE
    use MAORI_Set  , only : MAORI_Set_State_Obs_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(in)       ::  mstate
    type(T_MAORI_Data), intent(in)        ::  mdata
    integer, intent(in)                   ::  iset
    integer, intent(in)                   ::  obs_ivar
    integer, intent(out)                  ::  status
    
    integer, intent(in), optional         ::  iloc
    integer, intent(in), optional         ::  ilon, ilat
    real, intent(out), optional           ::  value

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_State_Obs_Get'

    ! --- begin -----------------------------------

    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! get simulated value:
    call MAORI_Set_State_Obs_Get( mstate%set(iset), mdata%set(iset), obs_ivar, status, &
                                       iloc=iloc, ilon=ilon, ilat=ilat, value=value )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_State_Obs_Get

  
  ! ***


  subroutine MAORI_State_Values_Get( mstate, mdata, values, status )

    use MAORI_Set, only : MAORI_Set_Data_Inq
    use MAORI_Set, only : MAORI_Set_State_Values_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(in)         ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    real, intent(out)                       ::  values(:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Values_Get'

    ! --- local -----------------------------------
    
    integer     ::  iset
    integer     ::  i
    integer     ::  nval

    ! --- begin -----------------------------------
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! inquire number of values in set:
      call MAORI_Set_Data_Inq( mdata%set(iset), status, nvalue=nval )
      IF_NOTOK_RETURN(status=1)
      ! some values available:
      if ( nval > 0 ) then
        ! check ...
        if ( size(values) < i+nval ) then
          write (gol,'("size of output array (",i6,") too small for number of values (",i6,") up to set ",i4)') &
                               size(values), i+nval, iset; call goErr
          TRACEBACK; status=1; return
        end if
        ! get values from set and insert in array:
        call MAORI_Set_State_Values_Get( mstate%set(iset), mdata%set(iset), values(i+1:i+nval), status )
        IF_NOTOK_RETURN(status=1)
        ! increase counter:
        i = i + nval
      end if  ! some values
    end do   ! sets

    ! ok
    status = 0

  end subroutine MAORI_State_Values_Get


  ! ***


  subroutine MAORI_State_Values_Put( mstate, mdata, values, status )

    use MAORI_Set, only : MAORI_Set_Data_Inq
    use MAORI_Set, only : MAORI_Set_State_Values_Put

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(inout)      ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    real, intent(in)                        ::  values(:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Values_Put'

    ! --- local -----------------------------------
    
    integer     ::  iset
    integer     ::  i
    integer     ::  nval

    ! --- begin -----------------------------------
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! inquire number of values in set:
      call MAORI_Set_Data_Inq( mdata%set(iset), status, nvalue=nval )
      IF_NOTOK_RETURN(status=1)
      ! some values available:
      if ( nval > 0 ) then
        ! check ...
        if ( size(values) < i+nval ) then
          write (gol,'("size of input array (",i6,") too small for number of values (",i6,") up to set ",i4)') &
                               size(values), i+nval, iset; call goErr
          TRACEBACK; status=1; return
        end if
        ! get values from set and insert in array:
        call MAORI_Set_State_Values_Put( mstate%set(iset), mdata%set(iset), values(i+1:i+nval), status )
        IF_NOTOK_RETURN(status=1)
        ! increase counter:
        i = i + nval
      end if  ! some values
    end do   ! sets

    ! ok
    status = 0

  end subroutine MAORI_State_Values_Put


  ! ***


  subroutine MAORI_State_LocValues_Get( mstate, mdata, iset, values, status, iloc, ilon, ilat )

    use MAORI_Set, only : MAORI_Set_Data_Inq
    use MAORI_Set, only : MAORI_Set_State_LocValues_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(in)         ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    integer, intent(in)                     ::  iset
    real, intent(out)                       ::  values(:)
    integer, intent(out)                    ::  status
    integer, intent(in), optional           ::  iloc
    integer, intent(in), optional           ::  ilon, ilat

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_Values_Get'

    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! get values for requested location:
    call MAORI_Set_State_LocValues_Get( mstate%set(iset), mdata%set(iset), values, status, &
                                           iloc=iloc, ilon=ilon, ilat=ilat )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_State_LocValues_Get


  ! ***


  subroutine MAORI_State_LocValues_Put( mstate, mdata, iset, values, status, iloc, ilon, ilat )

    use MAORI_Set, only : MAORI_Set_Data_Inq
    use MAORI_Set, only : MAORI_Set_State_LocValues_Put

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(inout)      ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    integer, intent(in)                     ::  iset
    real, intent(in)                        ::  values(:)
    integer, intent(out)                    ::  status
    integer, intent(in), optional           ::  iloc
    integer, intent(in), optional           ::  ilon, ilat

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_LocValues_Put'

    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! check ...
    call MAORI_Data_Check( mdata, status, iset=iset )
    IF_NOTOK_RETURN(status=1)

    ! put values for requested location:
    call MAORI_Set_State_LocValues_Put( mstate%set(iset), mdata%set(iset), values, status, &
                                           iloc=iloc, ilon=ilon, ilat=ilat )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine MAORI_State_LocValues_Put  

  
  ! ***


  subroutine MAORI_State_ObsValid_Get( mstate, mdata, values, status )

    use MAORI_Set, only : MAORI_Set_Data_Inq
    use MAORI_Set, only : MAORI_Set_Data_Obs_Get
    use MAORI_Set, only : MAORI_Set_State_Obs_Get

    ! --- in/out ---------------------------------

    type(T_MAORI_State), intent(in)         ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    real, intent(out)                       ::  values(:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_State_ObsValid_Get'

    ! --- local -----------------------------------
    
    integer     ::  iset
    integer     ::  stype
    integer     ::  i
    integer     ::  nval
    integer     ::  nloc, iloc
    integer     ::  obs_nvar, obs_ivar
    real        ::  y

    ! --- begin -----------------------------------
    
    ! latest filled value; nothing filled yet:
    i = 0
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! get type of output:
      call MAORI_Set_Data_Inq( mdata%set(iset), status, type=stype )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( stype /= MAORI_SAMPLE ) then
        write (gol,'(a," only implemented for samples yet")') rname; call goErr
        TRACEBACK; status=1; return
      end if
      ! inquire number of locations and observed values in set:
      call MAORI_Set_Data_Inq( mdata%set(iset), status, nloc=nloc, obs_nvar=obs_nvar )
      IF_NOTOK_RETURN(status=1)
      ! loop over stations:
      do iloc = 1, nloc
        ! loop over variables:
        do obs_ivar = 1, obs_nvar
          ! extract measurement:
          call MAORI_Set_Data_Obs_Get( mdata%set(iset), obs_ivar, status, y=y, iloc=iloc )
          IF_NOTOK_RETURN(status=1)
          ! valid data ?
          if ( y >= 0.0 ) then
            ! increase counter:
            i = i + 1
            ! storage should be large enough ...
            if ( i > size(values) ) then
              write (gol,'("size of output array (",i6,") too small for number of values (",i6,") up to set ",i4)') &
                                   size(values), i, iset; call goErr
              TRACEBACK; status=1; return
            end if
            ! get corresponding simulation:
            call MAORI_Set_State_Obs_Get( mstate%set(iset), mdata%set(iset), obs_ivar, status, iloc=iloc, &
                                             value=values(i) )
            IF_NOTOK_RETURN(status=1)
          end if  ! valid data
        end do  ! observed variables
      end do  ! station locations
    end do   ! sets

    ! ok
    status = 0

  end subroutine MAORI_State_ObsValid_Get


  ! ********************************************************************
  ! ***
  ! *** output routines
  ! ***
  ! ********************************************************************


  subroutine MAORI_Output_Init( mstate, mdata, name, status )

    use MAORI_Set, only : MAORI_Output_Set_Init

    ! --- in/out ---------------------------------

    type(T_MAORI_Output), intent(out)  ::  mstate
    type(T_MAORI_Data), intent(in)         ::  mdata
    character(len=*), intent(in)      ::  name
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Output_Init'

    ! --- local -----------------------------------

    integer       ::  iset

    ! --- begin -----------------------------------
    
    write (gol,'("MAORI: init output `",a,"` ...")') trim(name); call goPr
    
    ! store:
    mstate%name = trim(name)

    ! number of requested sets:
    mstate%nset = mdata%nset

    ! any request ?
    if ( mstate%nset > 0 ) then

      ! storage:
      allocate( mstate%set(1:mstate%nset) )

      ! loop over sets:
      do iset = 1, mdata%nset
        ! init set:
        call MAORI_Output_Set_Init( mstate%set(iset), mdata%set(iset), trim(name), status )
        IF_NOTOK_RETURN(status=1)
      end do  ! sets

    end if

    ! ok
    status = 0

  end subroutine MAORI_Output_Init


  ! ***


  subroutine MAORI_Output_Start( mstate, mdata, status )

    use MAORI_Set, only : MAORI_Output_Set_Start

    ! --- in/out ---------------------------------

    type(T_MAORI_Output), intent(inout)  ::  mstate
    type(T_MAORI_Data), intent(in)      ::  mdata
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Startup'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! start set:
      call MAORI_Output_Set_Start( mstate%set(iset), mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! ok
    status = 0

  end subroutine MAORI_Output_Start


  ! ***


  subroutine MAORI_Output_Done( mstate, mdata, status )

    use MAORI_Set, only : MAORI_Output_Set_Done
    
    ! --- in/out ---------------------------------

    type(T_MAORI_Output), intent(inout)      ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/MAORI_Output_Done'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------

    ! loop over sets:
    do iset = 1, mdata%nset
      ! done with set:
      call MAORI_Output_Set_Done( mstate%set(iset), mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! ok
    status = 0

  end subroutine MAORI_Output_Done


  ! ***


  subroutine MAORI_Output_Write( moutput, mstate, mdata, status )

    use MAORI_Set, only : MAORI_Output_Set_Write

    ! --- in/out ---------------------------------

    type(T_MAORI_Output), intent(inout)     ::  moutput
    type(T_MAORI_State), intent(in)         ::  mstate
    type(T_MAORI_Data), intent(in)          ::  mdata
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Output_Write'

    ! --- local -----------------------------------
    
    integer     ::  iset

    ! --- begin -----------------------------------
    
    ! loop over sets:
    do iset = 1, mdata%nset
      ! write set:
      call MAORI_Output_Set_Write( moutput%set(iset), mstate%set(iset), mdata%set(iset), status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! ok
    status = 0

  end subroutine MAORI_Output_Write
  

end module MAORI
