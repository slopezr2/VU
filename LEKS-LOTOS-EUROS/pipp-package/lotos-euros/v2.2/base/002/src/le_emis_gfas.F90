!###############################################################################
!
! LE_Emis_GFAS - GFAS fire emissions
!
! See description of GFAS service:
!   Copernicus Knowledge Base / Global Fire Assimilation System (GFAS) data documentation
!      https://confluence.ecmwf.int/display/CKB/CAMS%3A+Global+Fire+Assimilation+System+%28GFAS%29+data+documentation
!
!
! REFERENCES
! ----------
!
! Remy, S., Veira, A., Paugam, R., Sofiev, M., Kaiser, J. W., Marenco, F., Burton, S. P., Benedetti, A., 
!   Engelen, R. J., Ferrare, R., and Hair, J. W.: 
!   Two global data sets of daily fire emission injection heights since 2003, 
!   Atmos. Chem. Phys., 17, 2921-2942, https://doi.org/10.5194/acp-17-2921-2017, 2017.
!

! HISTORY
! -------
!
! 2011-11-03, Richard Kranenburg (TNO)
!   Changed routines such that fire emissions can be used for grids
!   other than MACC Emission grid. Emissions in grid cell from fire-files
!   distributed over underlying grid cell from run-domain
!
! 2012-04-11, Henk Eskes (KNMI)
!   Support files with list of fire sources as produced by IDL script.
!
! 2012-07-02, Arjo Segers (TNO)
!   Support files with compressed coordinates following CF-conventions
!   as produced by LEIP scripts.
!
! 2015-10, Arjo Segers (TNO)
!   Support CO2 emissions.
!   Bugs fixed:
!   - Only latest MACC/FIRE species that contributes to a
!     model CBM4 lumped species was taken into account, now changed.
!   - In "voc_to_cbm4_ole" an input molemass for "C3H6" of "56.0e-3"
!     was used, now changed to correct value.
!
! 2018-04, John Douros (KNMI)
!   Read and use hourly emissions files. Implementation is still a bit slopy, as 
!   it is based on a number of if-checks for the tres key.
!
!   NOTE: current implementation with hardcoded indices is tricky ..
!     In future better consider:
!      - read macc-to-model species mapping from csv file
!      - use loops over all species instead of subsets
!      - use specunits to convert to target units
!     See template settings in:
!       proj/macc3-co2/005/rc/lotos-euros-emissions-gfas.rc
!
! 2019-12, AMM
!   Fixed bug: emitted total PM is not just ppm, but also includes EC and POM;
!   corrected assignment to LE tracers.
!
! 2020-01, Marine Desmons, TNO
!   Distribute emisisons over height following (Remy, 2017)
!   using the latest height variables:
!
!      float mami(time, point) ;
!        mami:units = "m" ;
!        mami:long_name = "Mean altitude of maximum injection" ;
!      float apt(time, point) ;
!        apt:units = "m" ;
!        apt:long_name = "Altitude of plume top" ;
!      float apb(time, point) ;
!        apb:units = "m" ;
!        apb:long_name = "Altitude of plume bottom" ;
!
!   Alternative is available, but not used yet:
!
!      float injh(time, point) ;
!        injh:units = "m" ;
!        injh:long_name = "Injection height (from IS4FIRES)" ;
!    
! 2021-05, Arjo Segers, TNO
!   Re-implemented based on "le_emis_fire_macc.F90".
!   Uses composition table to form and assign emissions.
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_GFAS

  use NetCDF , only : NF90_StrError
  
  use GO     , only : gol, goPr, goErr
  use LE_Data, only : T_VarSum
  
  use GO, only : TDate

  implicit none


  ! --- in/out -----------------------------------

  private
  
  public  ::  T_Emis_GFAS


  ! --- const ------------------------------------

  ! module name:
  character(len=*), parameter ::  mname = 'LE_Emis_GFAS'


  ! --- types ------------------------------------
  
  ! linear combination of tracers:
  type, extends(T_VarSum) :: T_IndxVarSum
    ! tracer indices:
    integer, allocatable              ::  indx(:)     ! (n)
  contains
    procedure :: Init          => IndxVarSum_Init
    procedure :: Done          => IndxVarSum_Done
    procedure :: Show          => IndxVarSum_Show
    procedure :: AssignIndx    => IndxVarSum_AssignIndx
  end type T_IndxVarSum
  
  ! *
  
  type T_Composition
    ! filename:
    character(len=1024)                ::  filename
    ! number of substancances:
    integer                            ::  n
    ! emission properties:
    character(len=32), allocatable     ::  name(:)      ! (n)
    real, allocatable                  ::  molemass(:)  ! (n)
    type(T_IndxVarSum), allocatable    ::  emform(:)    ! (n)
    type(T_IndxVarSum), allocatable    ::  trform(:)    ! (n)
    ! list of unique emission variables:
    integer                            ::  nemvar
    character(len=16), allocatable     ::  emvars(:)    ! (n)
    ! limitted list of elements with data:
    integer                            ::  nd
    integer, allocatable               ::  irec(:)        ! (nd)
    !
  contains
    procedure   ::  Init          => Composition_Init
    procedure   ::  Done          => Composition_Done
!    procedure   ::  Inq_SubstID   => Composition_Inq_SubstID
  end type T_Composition
  
  ! *
  
  type T_Emis_GFAS
    ! label assigned to this emission:
    character(len=32)                 ::  label
    ! list of filename templates:
    character(len=1024)               ::  filename_templates
    ! allow missing files ?
    logical                           ::  allow_missing
    ! compmosition table:
    type(T_Composition)               ::  compo
    ! number of persistency days in forecast mode:
    integer                           ::  fc_nday_persist
    ! info?
    logical                           ::  verbose
    ! current
    character(len=1024)               ::  filename_curr
    ! dims:
    integer                           ::  npoint
    integer                           ::  nlon, nlat
    integer                           ::  ntime
    ! coordinates:
    real, allocatable                 ::  lons(:)  ! (nlon)
    real, allocatable                 ::  lon_bounds(:,:)  ! (2,nlon)
    real, allocatable                 ::  lats(:)  ! (nlat)
    real, allocatable                 ::  lat_bounds(:,:)  ! (2,nlat)
    integer, allocatable              ::  points(:)  ! (npoint)
    type(TDate), allocatable          ::  times(:)              ! (ntime)
    type(TDate), allocatable          ::  time_bounds(:,:)      ! (2,ntime)
    character(len=32)                 ::  time_bounds_from
    ! emission values:
    real, allocatable                 ::  emis(:,:,:)  ! (npoint,ntime,nd)
    character(len=32), allocatable    ::  units(:)  ! (nd)
    ! hourly profile:
    real                              ::  hour_fracs(0:23)
    ! injection heights:
    real                              ::  default_mami
    real                              ::  default_apb
    real                              ::  default_apt
    real, allocatable                 ::  mami(:,:)  ! (npoint,ntime)
    real, allocatable                 ::  apb(:,:)  ! (npoint,ntime)
    real, allocatable                 ::  apt(:,:)  ! (npoint,ntime)
    !
  contains
    procedure   ::  Init                => LE_Emis_GFAS_Init
    procedure   ::  Done                => LE_Emis_GFAS_Done
    procedure   ::  Clear               => LE_Emis_GFAS_Clear
    procedure   ::  GetFilename         => LE_Emis_GFAS_GetFilename
    procedure   ::  Setup               => LE_Emis_GFAS_Setup
    procedure   ::  AddEmis             => LE_Emis_GFAS_AddEmis
  end type T_Emis_GFAS
    


contains


  ! ===============================================================
  ! ===
  ! === GFAS
  ! ===
  ! ===============================================================


  subroutine LE_Emis_GFAS_Init( self, label, rcF, rckey, status )

    use GO     , only : TrcFile, ReadRc

    ! --- in/out ------------------------------

    class(T_Emis_GFAS), intent(out)       ::  self
    character(len=*), intent(in)          ::  label
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_GFAS_Init'

    ! --- local ------------------------------
    
    character(len=1024)       ::  compofile
    character(len=32)         ::  htype
    
    ! --- begin -------------------------------

    ! info?
    call rcf%Get( trim(rckey)//'.verbose', self%verbose, status )
    IF_NOTOK_RETURN(status=1)
    
    ! info ...
    if ( self%verbose ) then
      write (gol,'(a,": initialize GFAS fire emissions ...")'); call goPr
    end if
    
    ! store label:
    self%label = trim(label)

    ! filename templates:
    call rcf%Get( trim(rckey)//'.files', self%filename_templates, status )
    IF_NOTOK_RETURN(status=1)

    ! allow missing files ?
    call rcf%Get( trim(rckey)//'.allow_missing', self%allow_missing, status )
    IF_NOTOK_RETURN(status=1)
    
    ! optional: how to form time bounds if not in file
    call rcf%Get( trim(rckey)//'.time_bounds.from', self%time_bounds_from, status, default='' )
    IF_NOTOK_RETURN(status=1)
    
    ! defaults for injection heights:
    call rcf%Get( trim(rckey)//'.default_mami', self%default_mami, status )
    IF_NOTOK_RETURN(status=1)
    call rcf%Get( trim(rckey)//'.default_apt', self%default_apt, status )
    IF_NOTOK_RETURN(status=1)
    call rcf%Get( trim(rckey)//'.default_apb', self%default_apb, status )
    IF_NOTOK_RETURN(status=1)

    ! composition file:
    call rcf%Get( trim(rckey)//'.composition', compofile, status )
    IF_NOTOK_RETURN(status=1)
    ! read:
    call self%compo%Init( compofile, status )
    IF_NOTOK_RETURN(status=1)
    
    ! number of persistency days in forecast mode:
    call rcf%Get( trim(rckey)//'.fc_nday_persist', self%fc_nday_persist, status )
    IF_NOTOK_RETURN(status=1)

    ! hourly profile:
    call rcf%Get( trim(rckey)//'.hour_fracs', htype, status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    select case ( trim(htype) )
      !~ flat
      case ( 'constant' )
        ! fill:
        self%hour_fracs = 1.0
      !~
      ! Johannes Kaiser, ECMWF Tech. Memo. 596, August 2009, p4
      ! (May be suitable for the tropics, but probably too little night fires for Europe )
      case ( 'Kaiser' )
        ! fill:
        self%hour_fracs = (/ 0.2000, 0.2000, 0.2000, 0.2000, 0.2005, 0.2034, &
                             0.2195, 0.2873, 0.5047, 1.0283, 1.9534, 3.0909, &
                             3.9120, 3.9120, 3.0909, 1.9534, 1.0283, 0.5047, &
                             0.2873, 0.2195, 0.2034, 0.2005, 0.2000, 0.2000 /)
      !~
      ! The following form corresponds roughly to the results found for the US by
      ! Zhang et al, RSE 2008, doi: 10.1016/j.rse.2008.02.006, Fig.10
      case ( 'Zhang' )
        ! fill:
        self%hour_fracs = (/ 0.4000, 0.4000, 0.4000, 0.4000, 0.4003, 0.4025, &
                             0.4146, 0.4655, 0.6285, 1.0212, 1.7151, 2.5682, &
                             3.1840, 3.1840, 2.5682, 1.7151, 1.0212, 0.6285, &
                             0.4655, 0.4146, 0.4025, 0.4003, 0.4000, 0.4000 /)
      !~
      case default
        write (gol,'("unsupportetd hour profile `",a,"`")') trim(htype); call goErr
        TRACEBACK; status=1; return
    end select

    ! no file read yet:
    self%filename_curr = ''
             
    ! ok
    status = 0

  end subroutine LE_Emis_GFAS_Init


  ! ***


  subroutine LE_Emis_GFAS_Done( self, status )
    
    ! --- in/out ------------------------------

    class(T_Emis_GFAS), intent(inout)       ::  self
    integer, intent(out)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_GFAS_Done'

    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! done:
    call self%compo%Done( status )
    IF_NOTOK_RETURN(status=1)
      
    ! clear:
    call self%Clear( status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Emis_GFAS_Done


  ! ***


  subroutine LE_Emis_GFAS_Clear( self, status )
    
    ! --- in/out ------------------------------

    class(T_Emis_GFAS), intent(inout)       ::  self
    integer, intent(out)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_GFAS_Clear'

    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    if ( allocated(self%emis) ) then
      deallocate( self%lons, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%lon_bounds, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%lats, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%lat_bounds, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%points, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%times, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%time_bounds, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( self%emis, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%units, stat=status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( self%mami, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%apt, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( self%apb, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0

  end subroutine LE_Emis_GFAS_Clear


  ! ***


  subroutine LE_Emis_GFAS_GetFilename( self, t1, t2, filename, status )

    use GO     , only : TDate, Get
    use GO     , only : operator(-), operator(+), operator(<), operator(<=)
    use GO     , only : IsAnyDate, IncrDate, wrtgol
    use GO     , only : goReplace
    use GO     , only : goSplitString
    use Dims   , only : runF
 
    ! --- in/out ---------------------------

    class(T_Emis_GFAS), intent(inout)       ::  self
    type(TDate), intent(in)                 ::  t1, t2
    character(len=*), intent(out)           ::  filename
    integer, intent(out)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_GFAS_GetFilename'

    ! max number of file templates:
    integer, parameter    ::  maxfile = 10
    
    ! --- local ----------------------------
    
    integer                          ::  year, month, day, hour
    character(len=1024)              ::  filenames(maxfile)
    integer                          ::  nfile, ifile
    logical                          ::  exist
    
    ! --- begin ----------------------------

    ! forecast base specified?
    if ( IsAnyDate(runF%t_base) ) then
      ! target time:
      call Get( t1, year=year, month=month, day=day, hour=hour )
    else
      ! before forecast base ?
      if ( t1 < runF%t_base ) then
        ! as usual:
        call Get( t1, year=year, month=month, day=day, hour=hour )
        !
      ! within 2 days after forecast base?
      else if ( t1 <= runF%t_base+IncrDate(day=self%fc_nday_persist) ) then
        ! info ...
        write (gol,'("WARNING - no fire emission forecast, use data from day before forecast base!")'); call goPr
        ! use analysis day:
        call Get( runF%t_base-IncrDate(day=1), year=year, month=month, day=day, hour=hour )
        !
      ! no data ..
      else
        ! info ...
        write (gol,'("WARNING - no fire emission forecast, keep zero at this day!")'); call goPr
        ! ok
        status=0; return
      end if
    end if

    ! path might consist of several values,
    ! should be tested until the first match:
    call goSplitString( trim(self%filename_templates), nfile, filenames, status )
    IF_NOTOK_RETURN(status=1)
    ! dummy:
    filename = '/no/file/found/yet'
    ! loop over directories:
    do ifile = 1, nfile
      ! replace some values:
      call goReplace( filenames(ifile), '%{year}', '(i4.4)', year, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filenames(ifile), '%{month}', '(i2.2)', month, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filenames(ifile), '%{day}', '(i2.2)', day, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( filenames(ifile), '%{hour}', '(i2.2)', hour, status )
      IF_NOTOK_RETURN(status=1)
      ! check, if present then store match and leave:
      inquire( file=trim(filenames(ifile)), exist=exist )
      if ( exist ) then
        filename = trim(filenames(ifile))
        exit
      end if
    end do
    ! check ..
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      if ( self%allow_missing) then
        write (gol,'("WARNING - missing MACC fire emissions file; continue ...")'); call goPr
        status=-1; return
      else
        write (gol,'("None of the specified MACC fire emission files found:")'); call goErr
        do ifile = 1, nfile
          write (gol,'("  ",a)') trim(filenames(ifile)); call goErr
        end do
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine LE_Emis_GFAS_GetFilename


  ! ***


  subroutine LE_Emis_GFAS_Setup( self, emis_a, t1, t2, status )

    use GO     , only : TDate, IncrDate, wrtgol
    use GO     , only : operator(-), operator(+), operator(<), operator(<=)
    use GO     , only : goReplace
    use GO     , only : goMatchValue
    use GO     , only : goSplitString
    use C3PO   , only : T_File_Nc
    use Dims   , only : nx, ny, nz
    use Indices, only : nspec
 
    ! --- in/out ---------------------------

    class(T_Emis_GFAS), intent(inout)       ::  self
    real, intent(inout)                     ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)                 ::  t1, t2
    integer, intent(out)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_GFAS_Setup'

    ! --- local ----------------------------
    
    character(len=1024)              ::  filename
    type(T_File_Nc)                  ::  ncfile
    character(len=256)               ::  description
    character(len=64)                ::  vname
    integer                          ::  nvar, ivar
    integer                          ::  varid
    real, allocatable                ::  emis_in(:,:,:)   ! (npoint,ntime,nemvar)
    character(len=32), allocatable   ::  units_in(:)      ! (nemvar)
    integer                          ::  iemvar
    character(len=32)                ::  units
    integer                          ::  itime
    integer                          ::  id
    integer                          ::  irec
    integer                          ::  j
    real, allocatable                ::  emis_fire(:,:,:,:)   ! (nx,ny,nz,nspec)
    ! --- begin ----------------------------

    ! info:
    if ( self%verbose ) then
      call wrtgol( rname//': setup GFAS fire emissions for ', (/t1,t2/) ); call goPr
    end if

    ! check ...
    if ( self%compo%nemvar == 0 ) then
      ! info ...
      if ( self%verbose ) then
        write (gol,'("warning - no gfas fire emission variables selected")'); call goPr
      end if
      ! leave, status ok
      status=0; return
    end if
    
    ! get filename for current time,
    ! status<0 is warning that no file is found but that is allowed ...
    call self%GetFilename( t1, t2, filename, status )
    if ( status < 0 ) return
    IF_NOTOK_RETURN(status=1)
    
    ! new file?
    if ( trim(filename) /= trim(self%filename_curr) ) then
      ! store:
      self%filename_curr = trim(filename)

      ! info:
      if ( self%verbose ) then
        write (gol,'(a,":   read from file : ",a)') trim(filename); call goPr
      end if

      ! clear:
      call self%Clear( status )
      IF_NOTOK_RETURN(status=1)

      ! open:
      call ncfile%Open( filename, status )
      IF_NOTOK_RETURN(status=1)

      ! dims:
      call ncfile%Inquire_Dimension( 'point', status, length=self%npoint )
      IF_NOTOK_RETURN(status=1)
      call ncfile%Inquire_Dimension( 'lon', status, length=self%nlon )
      IF_NOTOK_RETURN(status=1)
      call ncfile%Inquire_Dimension( 'lat', status, length=self%nlat )
      IF_NOTOK_RETURN(status=1)
      call ncfile%Inquire_Dimension( 'time', status, length=self%ntime )
      IF_NOTOK_RETURN(status=1)
      
      ! storage for coordinate:
      allocate( self%lons(self%nlon), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call ncfile%Get_Var( 'var_name=lon', self%lons, units, status )
      IF_NOTOK_RETURN(status=1)
      
      ! storage for coordinate:
      allocate( self%lon_bounds(2,self%nlon), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call ncfile%Get_Var( 'var_name=lon_bounds', self%lon_bounds, units, status )
      IF_NOTOK_RETURN(status=1)
      
      ! storage for coordinate:
      allocate( self%lats(self%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call ncfile%Get_Var( 'var_name=lat', self%lats, units, status )
      IF_NOTOK_RETURN(status=1)
      
      ! storage for coordinate:
      allocate( self%lat_bounds(2,self%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call ncfile%Get_Var( 'var_name=lat_bounds', self%lat_bounds, units, status )
      IF_NOTOK_RETURN(status=1)
      
      ! storage for coordinate:
      allocate( self%points(self%npoint), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call ncfile%Get_Var( 'var_name=point;units=1', self%points, units, status )
      IF_NOTOK_RETURN(status=1)
      
      ! storage for coordinate:
      allocate( self%times(self%ntime), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call ncfile%Get_Var( 'var_name=time', self%times, status )
      IF_NOTOK_RETURN(status=1)
      
      ! time bounds:
      allocate( self%time_bounds(2,self%ntime), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! variable description:
      description = 'var_name=time_bounds'
      ! in file?
      call ncfile%Inq_VarID( description, varid, status, quiet=.true. )
      if ( status > 0 ) then
        TRACEBACK; status=1; return
      else if ( status == 0 ) then
        ! read:
        call ncfile%Get_Var( description, self%time_bounds, status )
        IF_NOTOK_RETURN(status=1)
      else
        ! no time bounds in file, how to form?
        select case ( trim(self%time_bounds_from) )
          !~ time in file is end of day:
          case ( 'end-of-day' )
            ! loop:
            do itime = 1, self%ntime
              ! start time:
              self%time_bounds(1,itime) = self%times(itime) - IncrDate(day=1)
              ! end time:
              self%time_bounds(2,itime) = self%times(itime)
            end do ! itime
          !~ time in file is end of hour:
          case ( 'end-of-hour' )
            ! loop:
            do itime = 1, self%ntime
              ! start time:
              self%time_bounds(1,itime) = self%times(itime) - IncrDate(hour=1)
              ! end time:
              self%time_bounds(2,itime) = self%times(itime)
            end do ! itime
          !~
          case default
            write (gol,'("unsupported time_bounds keyword `",a,"`")') trim(self%time_bounds_from); call goErr
            TRACEBACK; status=1; return
        end select
      end if
      
      ! height info:
      allocate( self%mami(self%npoint,self%ntime), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! variable description:
      description = 'long_name=Mean altitude of maximum injection;var_name=mami'
      ! in file?
      call ncfile%Inq_VarID( description, varid, status, quiet=.true. )
      if ( status > 0 ) then
        TRACEBACK; status=1; return
      else if ( status == 0 ) then
        ! read:
        call ncfile%Get_Var( description, self%mami, units, status )
        IF_NOTOK_RETURN(status=1)
      else
        ! info ...
        if ( self%verbose ) then
          write (gol,'("    no `mami`, use default value ",f12.4)') self%default_mami; call goPr
        end if
        ! default:
        self%mami = self%default_mami
      end if
      
      ! height info:
      allocate( self%apt(self%npoint,self%ntime), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! variable description:
      description = 'long_name=Altitude of plume top;var_name=apt'
      ! in file?
      call ncfile%Inq_VarID( description, varid, status, quiet=.true. )
      if ( status > 0 ) then
        TRACEBACK; status=1; return
      else if ( status == 0 ) then
        ! read:
        call ncfile%Get_Var( description, self%apt, units, status )
        IF_NOTOK_RETURN(status=1)
      else
        ! info ...
        if ( self%verbose ) then
          write (gol,'("    no `apt` , use default value ",f12.4)') self%default_apt; call goPr
        end if
        ! default:
        self%apt = self%default_apt
      end if
      
      ! height info:
      allocate( self%apb(self%npoint,self%ntime), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! variable description:
      description = 'long_name=Altitude of plume bottom;var_name=apb'
      ! in file?
      call ncfile%Inq_VarID( description, varid, status, quiet=.true. )
      if ( status > 0 ) then
        TRACEBACK; status=1; return
      else if ( status == 0 ) then
        ! read:
        call ncfile%Get_Var( description, self%apb, units, status )
        IF_NOTOK_RETURN(status=1)
      else
        ! info ...
        if ( self%verbose ) then
          write (gol,'("    no `apb` , use default value ",f12.4)') self%default_apb; call goPr
        end if
        ! default:
        self%apb = self%default_apb
      end if
      
      ! storage for emissions from file:
      allocate( emis_in(self%npoint,self%ntime,self%compo%nemvar), source=0.0, stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( units_in(self%compo%nemvar), stat=status )
      IF_NOTOK_RETURN(status=1)

      ! number of variables:
      call ncfile%Inquire( status, nvar=nvar )
      IF_NOTOK_RETURN(status=1)
      ! loop over variables:
      do ivar = 1, nvar
        ! current name:
        call ncfile%Inq_Var( ivar, status, name=vname )
        IF_NOTOK_RETURN(status=1)
        ! check ..
        select case ( trim(vname) )
          !~ coordinates:
          case ( 'point', 'time', 'lon', 'lat', 'lon_bounds', 'lat_bounds' )
            ! ok
          !~ plume height variables:
          case ( 'mami', 'apt', 'injh', 'apb' )
            ! ok
          !~ emission:
          case default
            ! probably an emission variable, search in list of variables to be used:
            call goMatchValue( vname, self%compo%emvars(1:self%compo%nemvar), iemvar, status, quiet=.true. )
            if ( status < 0 ) then
              ! testing ..
              if ( self%verbose ) then
                write (gol,'(a,": emission `",a,"` not used ...")') rname, trim(vname); call goPr
              end if
            else
              ! testing ..
              if ( self%verbose ) then
                write (gol,'(a,": emission `",a,"` read ...")') rname, trim(vname); call goPr
              end if
              ! read:
              call ncfile%Get_Var( 'var_name='//trim(vname), emis_in(:,:,iemvar), units_in(iemvar), status )
              IF_NOTOK_RETURN(status=1)
            end if
        end select ! vname
      end do ! ivar
      
      ! target emissions, combination of input:
      allocate( self%emis(self%npoint,self%ntime,self%compo%nd), source=0.0, stat=status )
      IF_NOTOK_RETURN(status=1)
      allocate( self%units(self%compo%nd), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over target emissions:
      do id = 1, self%compo%nd
        ! compo record:
        irec = self%compo%irec(id)
        ! testing ...
        write (gol,'(a,": fire emis ",i2," from composition ",i0," `",a,"`")') &
                    rname, id, irec, trim(self%compo%name(irec)); call goPr
        ! no units yet:
        self%units(id) = ''
        ! loop over contributions:
        do j = 1, self%compo%emform(irec)%n
          ! corresponding emisison record:
          iemvar = self%compo%emform(irec)%indx(j)
          ! check units:
          if ( len_trim(self%units(id)) == 0 ) then
            self%units(id) = trim(units_in(iemvar))
          else if ( trim(units_in(iemvar)) /= trim(self%units(id)) ) then
            write (gol,'("current contribution has units `",a,"` while previous had `",a,"`")') &
                     trim(units_in(iemvar)), trim(self%units(id)); call goErr
            TRACEBACK; status=1; return
          end if
          ! testing ...
          write (gol,'(a,":  add ",f5.2," emvar ",i2," `",a,"`")') &
                  rname, self%compo%emform(irec)%w(j), iemvar, trim(self%compo%emvars(iemvar)); call goPr
          ! add contribution:
          self%emis(:,:,id) = self%emis(:,:,id) + self%compo%emform(irec)%w(j) * emis_in(:,:,iemvar)
        end do ! irec
      end do ! id
      
      ! clear:
      deallocate( emis_in, stat=status )
      IF_NOTOK_RETURN(status=1)
      deallocate( units_in, stat=status )
      IF_NOTOK_RETURN(status=1)

      ! close netcdf file
      call ncfile%Close( status )
      IF_NOTOK_RETURN(status=1)

    end if  ! new file
    
    ! storage for fire emissions:
    allocate( emis_fire(nx,ny,nz,nspec), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! fill current emissions:
    call self%AddEmis( emis_fire, t1, t2, status )
    IF_NOTOK_RETURN(status=1)
    
    ! add:
    emis_a = emis_a + emis_fire

    
    ! clear:
    deallocate( emis_fire, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_GFAS_Setup


  ! ***


  subroutine LE_Emis_GFAS_AddEmis( self, emis_fire, t1, t2, status )

    use GO     , only : TDate, wrtgol
    use GO     , only : operator(<=)
    use Grid   , only : ll_area_deg_m2
    use LE_Grid, only : ugg
    use Dims   , only : nx, ny, nz
    use Indices, only : nspec, specunit, specname
    use LE_Data, only : LE_Data_GetPointer
 
    ! --- in/out ---------------------------

    class(T_Emis_GFAS), intent(inout)       ::  self
    real, intent(out)                       ::  emis_fire(nx,ny,nz,nspec)
    type(TDate), intent(in)                 ::  t1, t2
    integer, intent(out)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_GFAS_AddEmis'

    ! conversion factors:
    real, parameter    ::  ug_per_kg   = 1.0e9   ! ug/kg
    real, parameter    ::  sec_per_min = 60.0    ! sec/min

    ! --- local ----------------------------
    
    logical                          ::  found
    integer                          ::  itime
    real                             ::  hfrac
    integer                          ::  n_in_domain
    integer                          ::  npoint
    integer                          ::  ipoint
    integer                          ::  ilon, ilat
    real, pointer                    ::  fracs(:,:)
    real                             ::  west, east, south, north
    real                             ::  area
    integer                          ::  i1, i2, j1, j2
    integer                          ::  i, j, k
    real, pointer                    ::  halt(:,:,:)   ! (lon,lat,0:lev)
    real, allocatable                ::  fraction(:)
    integer                          ::  id
    integer                          ::  irec
    integer                          ::  iw
    integer                          ::  ispec
    real                             ::  hour_frac
    character(len=64)                ::  conversion
    
    ! --- begin ----------------------------

    ! Meteo pointers:
    call LE_Data_GetPointer( 'halt', halt, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    
    ! search time record ...
    found = .false.
    ! loop over records:
    do itime = 1, self%ntime
      ! interval within bounds?
      found = (self%time_bounds(1,itime) <= t1) .and. (t2 <= self%time_bounds(2,itime))
      ! leave?
      if ( found ) exit
    end do
    ! not found?
    if ( .not. found ) then
      call wrtgol( 'no time record for interval ', (/t1,t2/) ); call goErr
      do itime = 1, self%ntime
        call wrtgol( '  record: ', self%time_bounds(:,itime) ); call goErr
      end do
      TRACEBACK; status=1; return
    end if
    ! testing ...
    call wrtgol( rname//': add time record valid for ', self%time_bounds(:,itime) ); call goPr
    
    ! hour fraction, might be 1.0 in case of hourly records:
    hour_frac = self%hour_fracs(t1%hour)
    
    ! init result:
    emis_fire = 0.0

    ! profile:    
    allocate( fraction(nz), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! counter:
    n_in_domain = 0
    ! init fractions:
    nullify( fracs )
    ! loop over points:
    do ipoint = 1, self%npoint
      ! point gives the zero-based 1D index of the compressed coordinates,
      ! expand to 1-based:
      ilat = self%points(ipoint) / self%nlon + 1
      ilon = modulo( self%points(ipoint), self%nlon ) + 1
      ! boundingbox of emission:
      west  = self%lon_bounds(1,ilon)
      east  = self%lon_bounds(2,ilon)
      south = min( self%lat_bounds(1,ilat), self%lat_bounds(2,ilat) )
      north = max( self%lat_bounds(1,ilat), self%lat_bounds(2,ilat) )
      ! emission area:
      area = ll_area_deg_m2( west, east, south, north ) ! m2
      ! Given a box with sides [west,east]x[south,north] in degrees;
      ! determine how to distribute it's area over a regular lat/lon grid.
      ! Returns indices ranges i1..i2 and j1..j2 and fractions
      !   fracs(i1:i2,j1:j2)
      ! Return status:
      !   -1          : location outside domain
      !   0           : ok
      !   >0          : error
      call ugg%GetDistribution( west, east, south, north, &
                                 i1, i2, j1, j2, fracs, status )
      if ( status == -1 ) cycle  ! cell outside of target domain
      IF_NOTOK_RETURN(status=1)
      ! increase counter:
      n_in_domain = n_in_domain + 1
      ! loop over covered cells:
      do i = i1, i2
        do j = j1, j2

          ! injection profile:
          call GetInjectionProfile( self%mami(ipoint,itime), self%apb(ipoint,itime), self%apt(ipoint,itime), &
                                      halt(i,j,0:nz), fraction, status )
          IF_NOTOK_RETURN(status=1)
          ! testing ...
          !write (gol,*) 'xxx1 ', ipoint, i, j, ';', self%mami(ipoint,itime), self%apb(ipoint,itime), self%apt(ipoint,itime), ';', fraction; call goPr

          ! loop over emission datasets:
          do id = 1, self%compo%nd
            ! compo record:
            irec = self%compo%irec(id)
            ! loop over tracers:
            do iw = 1, self%compo%trform(irec)%n

              ! corresponding tracer index:
              ispec = self%compo%trform(irec)%indx(iw)

              !! testing ...
              !write (gol,'(a,":  add ",f5.2," `",a,"` to `",a,"`")') &
              !        rname, self%compo%trform(irec)%w(iw), trim(self%compo%name(irec)), trim(specname(ispec)); call goPr

              ! conversion:
              if ( specunit(ispec) == 'ppb' ) then
                conversion = trim(self%units(id))//' -> mol/min'
              else if ( specunit(ispec) == 'ug/m3' ) then
                conversion = trim(self%units(id))//' -> ug/min'
              else
                write (gol,'("unsupported units `",a,"` for target tracer `",a,"`")') trim(specunit(ispec)), trim(specname(ispec)); call goErr
                TRACEBACK; status=1; return
              end if
              ! check units:
              select case ( trim(conversion) )
                !~
                case ( 'kg m**-2 s**-1 -> mol/min' )
                  ! add contribution:
                  emis_fire(i,j,:,ispec) = emis_fire(i,j,:,ispec) + &   ! mol/min
                         fraction                        &  ! 1        (height profile)
                         * self%compo%trform(irec)%w(iw) &  ! 1        (fraction assigned to this tracer)
                         * self%emis(ipoint,itime,id)    &  ! kg/m2/s
                         / self%compo%molemass(irec)     &  ! mol/kg
                         * area * fracs(i,j)             &  ! m2       (fraction of cell overlapping with emission area)
                         * sec_per_min                   &  ! s/min
                         * hour_frac                        ! 1
                !~
                case ( 'kg m**-2 s**-1 -> ug/min' )
                  ! add contribution:
                  emis_fire(i,j,:,ispec) = emis_fire(i,j,:,ispec) + &   ! ug/min
                         fraction                        &  ! 1        (height profile)
                         * self%compo%trform(irec)%w(iw) &  ! 1        (fraction assigned to this tracer)
                         * self%emis(ipoint,itime,id)    &  ! kg/m2/s
                         * ug_per_kg                     &  ! ug/kg
                         * area * fracs(i,j)             &  ! m2       (fraction of cell overlapping with emission area)
                         * sec_per_min                   &  ! s/min
                         * hour_frac                        ! 1
                !~
                case default
                  write (gol,'("unsupported emission conversion `",a,"`")') trim(conversion); call goErr
                  TRACEBACK; status=1; return
              end select

            end do ! irec
          end do ! id

        end do  ! i
      end do   ! j

    end do  ! points

    ! info ...
    if ( self%verbose ) then
       write (gol,'(a,":   nr of fire locations = ",i5,",  inside domain = ",i5)' ) &
                                      rname, self%npoint, n_in_domain ; call goPr
    end if

    ! clear:
    deallocate( fraction, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    if ( associated(fracs) ) then
      deallocate( fracs, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_GFAS_AddEmis
  
  
  ! ***
  
  !
  ! Emission injection profile:
  !
  !     !
  !     |\\         apt  = altitude of plume top
  !     |\\
  !     |\\
  !     |////\\      mami = mean altitude of maximum injection
  !     |////       apb  = altiude of plume bottom
  !     |
  !  ---*----------------
  !
  !  The profile assigns:
  !     \\   50% to the layers >= mami
  !     //   50% to the layers <= mami
  !
  ! When mami=0 then a smouldering fire is assumed,
  ! and all emissions are injected in the first layer.
  !
  
  subroutine GetInjectionProfile( mami, apb, apt, hhb, fraction, status )
    
    use Dims, only : nz

    ! --- in/out ---------------------------------
    
    real, intent(inout)                           ::  mami ! m
    real, intent(in)                              ::  apb, apt ! m
    real, intent(in)                              ::  hhb(0:nz)   ! half-level altitudes [m]
    real, intent(out)                             ::  fraction(nz) ! 1
    integer, intent(out)                          ::  status
      
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GetInjectionProfile'
    
    ! --- local ----------------------------------
    
    integer     ::  k
    integer     ::  maxinjlev, minlev, maxlev
    real        ::  laysum
    real        ::  fsum

    ! --- begin ----------------------------------
    
    ! init profile:
    fraction = 0.0
  
    ! check if injectoin is at surface leve.
    if ( mami <= 0.0 ) then

      ! smouldering, inject in surface layer:
      fraction(1) = 1.0

    else
    
      ! dummy to avoid warnings about possibly unitialized variables ...
      minlev = -999
      ! layer with bottom:
      if ( apb < hhb(0) ) then
        ! below surface, take first layer:
        minlev = 1
      else if ( apb < hhb(nz) ) then
        ! loop:
        do k = 1, nz
          ! "altitude of plume bottom" in this layer?
          if ((hhb(k-1) <= apb) .and. (apb < hhb(k)))  then 
            minlev = k
            exit
          end if
        end do ! k
      else
        ! plume bottom above model top ...
        minlev = nz+1
      end if
      
      ! dummy to avoid warnings about possibly unitialized variables ...
      maxlev = -999
      ! layer with top:
      if ( apt < hhb(0) ) then
        ! below surface, take first layer:
        maxlev = 1
      else if ( apt < hhb(nz) ) then
        ! loop:
        do k = 1, nz
          ! "altitude of plume bottom" in this layer?
          if ((hhb(k-1) <= apt) .and. (apt < hhb(k)))  then 
            maxlev = k
            exit
          end if
        end do ! k
      else
        ! plume top above model top, assign to top anyway ...
        maxlev = nz
      end if

      !! check ...
      if ( (mami < apb) .or. (mami > apt) ) then      
      !  write (gol,'("mean altitude of injection should be between bottom and top variables:")'); call goErr
      !  write (gol,'("  altitude of plume top              : ",f16.8," m")') apt; call goErr
      !  write (gol,'("  mean altitude of maximum injection : ",f16.8," m")') mami; call goErr
      !  write (gol,'("  altitude of plume bottom           : ",f16.8," m")') apb; call goErr
      !  TRACEBACK; status=1; return
        mami = apb + (apt-apb)/2.0
      end if
      ! dummy to avoid warnings about possibly unitialized variables ...
      maxinjlev = -999
      ! layer with maximum:
      if ( mami < hhb(0) ) then
        ! below surface, take first layer:
        maxinjlev = 1
      else if ( mami < hhb(nz) ) then
        ! loop:
        do k = 1, nz
          ! "mean altitude of maximum injection" in this layer?
          if ((hhb(k-1) <= mami) .and. (mami < hhb(k)))  then 
            maxinjlev = k
            exit
          end if
        end do ! k
      else
        ! plume bottom above model top ..
        maxinjlev = nz+1
      end if
      
      ! plume below model top?
      if ( minlev <= nz ) then

        ! asssign 50% to layers <= mami
        if ( minlev == maxinjlev ) then
          fraction(maxinjlev) = 0.5
        else

          laysum=0
          do k = minlev, maxinjlev
            laysum = laysum+(hhb(k)-hhb(k-1))
          end do 

          do k = minlev, maxinjlev
            fraction(k) = 0.5*(hhb(k)-hhb(k-1))/laysum
          end do 
        endif

        ! assign 50% to layers >= mami
        if ( maxinjlev == maxlev ) then
          fraction(maxinjlev) = fraction(maxinjlev) + 0.5
        else
          laysum=0
          do k = maxinjlev, maxlev
            laysum = laysum + (hhb(k)-hhb(k-1))
          end do  
          do k = maxinjlev, maxlev
            fraction(k)=fraction(k)+0.5*(hhb(k)-hhb(k-1))/laysum 
          end do 
        end if
        
      end if  ! plume below model top

    end if  ! surface or profile
    
    ! check ...
    fsum = sum(fraction)
    if (  (fsum > 0.0) .and. (abs(fsum-1.0) > 1.0e-4) ) then
      write (gol,'("sum over fraction not 1.0 but: ",f16.8)') fsum; call goPr
      write (gol,'("injection variables:")'); call goErr
      write (gol,'("  altitude of plume top              : ",f16.8," m")') apt; call goErr
      write (gol,'("  mean altitude of maximum injection : ",f16.8," m")') mami; call goErr
      write (gol,'("  altitude of plume bottom           : ",f16.8," m")') apb; call goErr
      write (gol,'("estimated profile:")'); call goErr
      do k = 1, nz
        write (gol,'("  layer ",i3," [",f9.2,",",f9.2,"]  : ",f8.2)') k, hhb(k-1), hhb(k), fraction(k); call goErr
      end do
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine GetInjectionProfile


  ! ===============================================================
  ! ===
  ! === IndxVarSum
  ! ===
  ! ===============================================================


  subroutine IndxVarSum_Init( self, input, status, verbose )

    use LE_Data_VarSum, only : VarSum_Init
    
    ! --- in/out ---------------------------------
    
    class(T_IndxVarSum), intent(out)              ::  self
    character(len=*), intent(in)                  ::  input
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  verbose
      
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/IndxVarSum_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init ancestor:
    call VarSum_Init( self, input, status, verbose=verbose )
    IF_NOTOK_RETURN(status=1)
    
    ! any elements?
    if ( self%n > 0 ) then
      ! extra storage:
      allocate( self%indx(self%n), source=-999, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine IndxVarSum_Init
  
  
  ! ***


  subroutine IndxVarSum_Done( self, status )
  
    use LE_Data_VarSum, only : VarSum_Done
    
    ! --- in/out ---------------------------------
    
    class(T_IndxVarSum), intent(inout)            ::  self
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/IndexVarSum_Done'
      
    ! --- begin ----------------------------------
    
    ! any storage?
    if ( self%n > 0 ) then
      ! clear:
      deallocate( self%indx, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! init ancestor:
    call VarSum_Done( self, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine IndxVarSum_Done


  ! ***


  subroutine IndxVarSum_Show( self, prompt, status )
    
    ! --- in/out ---------------------------------
    
    class(T_IndxVarSum), intent(in)               ::  self
    character(len=*), intent(in)                  ::  prompt
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/IndxVarSum_Done'
      
    ! --- local ----------------------------------
      
    integer     ::  i
      
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a,"number of elements: ",i0)') prompt, self%n; call goPr
    ! any?
    if ( self%n > 0 ) then
      ! loop:
      do i = 1, self%n
        write (gol,'(a,"..",i2," ",f16.4," ",i6," ",a)') prompt, i, self%w(i), self%indx(i), trim(self%name(i)); call goPr
      end do ! i
    end if
    
    ! ok
    status = 0
    
  end subroutine IndxVarSum_Show
  
  
  ! ***
  
  ! for each element,
  ! search index of "name" in "values(:)" list;
  ! return status:
  !   -2 : all undefined
  !   -1 : some undefined
  !    0 : all found (or none defined)
  !   >0 : error

  subroutine IndxVarSum_AssignIndx( self, values, status, quiet )
  
    use GO, only : goMatchValue
    
    ! --- in/out ---------------------------------
    
    class(T_IndxVarSum), intent(inout)            ::  self
    character(len=*), intent(in)                  ::  values(:)
    integer, intent(out)                          ::  status
    logical, intent(in), optional                 ::  quiet

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/IndexVarSum_AssignIndx'
      
    ! --- local ----------------------------------
    
    integer     ::  i
    integer     ::  nmissing
      
    ! --- begin ----------------------------------
    
    ! any?
    if ( self%n > 0 ) then

      ! init counter:
      nmissing = 0
      ! loop over elements:
      do i = 1, self%n
        ! assign index:
        call goMatchValue( self%name(i), values, self%indx(i), status, quiet=quiet )
        ! not found?
        if ( status < 0 ) then
          ! set warning value:
          self%indx(i) = -999
          ! increse counter:
          nmissing = nmissing + 1
        else if ( status > 0 ) then
          TRACEBACK; status=1; return
        end if
      end do
      
      ! set status:
      if ( nmissing == self%n ) then
        status = -2
      else if ( nmissing > 0 ) then
        status = -1
      else
        status = 0
      end if
      
    else
      ! always ok
      status = 0
    end if
    
  end subroutine IndxVarSum_AssignIndx


  ! ===============================================================
  ! ===
  ! === composition table
  ! ===
  ! ===============================================================
  

  subroutine Composition_Init( self, filename, status )

    use GO     , only :  goGetFU
    use GO     , only :  goVarValue, goSplitString, goMatchValue
    use LE_Data, only :  T_VarSum
    use Indices, only :  specname

    ! --- in/out ------------------------------
    
    class(T_Composition), intent(out)           ::  self
    character(len=*), intent(in)                ::  filename
    integer, intent(out)                        ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Composition_Init'
    
    ! maximum number of rows and columns in csv file:
    integer, parameter    ::  maxrow = 100
    integer, parameter    ::  maxcol = 40
    
    ! --- local -------------------------------
    
    logical                 ::  exist
    character(len=1)        ::  comment
    character(len=1)        ::  sep

    integer                 ::  ifield_name
    integer                 ::  ifield_vars
    integer                 ::  ifield_molm
    integer                 ::  ifield_tracers

    integer                 ::  fu
    character(len=1024)     ::  line
    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    integer                 ::  nfield
    integer                 ::  ifield
    character(len=64)       ::  fields(maxcol)
    integer                 ::  iline
    type(T_VarSum)          ::  varsum
    integer                 ::  i
    integer                 ::  j
    character(len=16)       ::  name
    integer                 ::  indx
    
    ! --- begin -------------------------------
    
    ! testing ...
    write (gol,'(a,": read composition table: ",a)') rname, trim(filename); call goPr
    
    ! store:
    self%filename = trim(filename)
    
    ! check ..
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("composition file not found: ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! seperation character:
    sep = ';'
    ! comment character:
    comment = '#'

    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( fu, file=trim(filename), status='old', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! line counter:          
    iline = 0

    ! read header line after first comment:
    do
      ! read line:
      read (fu,'(a)',iostat=status) line
      if (status/=0) then
        write (gol,'("reading header line from file : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle
      ! found non-comment line, leave loop:
      exit
    end do
    ! split:
    call goSplitString( line, nheader, headers, status, sep=sep )
    IF_NOTOK_RETURN(status=1)
    ! search columns:
    call goMatchValue( 'name', headers, ifield_name, status )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( 'variables', headers, ifield_vars, status )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( 'molemass', headers, ifield_molm, status )
    IF_NOTOK_RETURN(status=1)
    call goMatchValue( 'tracers', headers, ifield_tracers, status )
    IF_NOTOK_RETURN(status=1)

    ! storage for tracers, fill with warning value:
    allocate( self%name(maxrow), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%molemass(maxrow), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%emform(maxrow), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%trform(maxrow), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! no records yet:
    self%n = 0

    ! loop over records:
    do

      ! increase record counter:
      iline = iline + 1
      ! try to read line:
      read (fu,'(a)',iostat=status) line
      ! eof ?
      if (status<0) exit
      ! error ?
      if (status>0) then
        write (gol,'("reading line ",i6," from file : ",a)') iline, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle

      ! split into records:
      call goSplitString( line, nfield, fields, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( nfield /= nheader ) then
        write (gol,'("number of fields (",i0,") in line ",i0," :")') nfield, iline; call goPr
        write (gol,'("  ",a)') trim(line); call goErr
        write (gol,'("fields:")'); call goErr
        do ifield = 1, nfield
          write (gol,'(i6," : ",a)') ifield, trim(fields(ifield)); call goErr
        end do
        write (gol,'("differs from number of headers (",i0,") in file ",a)') nheader, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! increase counter:
      self%n = self%n + 1

      ! store:
      self%name(self%n) = trim(fields(ifield_name))
      
      ! formula for source variables:
      call self%emform(self%n)%Init( fields(ifield_vars), status )
      IF_NOTOK_RETURN(status=1)

      ! evaluate molemass description:
      call varsum%Init( fields(ifield_molm), status )
      IF_NOTOK_RETURN(status=1)
      ! interpreted as molemass ...
      call varsum%Get( status, molemass=self%molemass(self%n) )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      call varsum%Done( status )
      IF_NOTOK_RETURN(status=1)
      
      ! formula for source variables:
      call self%trform(self%n)%Init( fields(ifield_tracers), status )
      IF_NOTOK_RETURN(status=1)

    end do  ! lines

    ! close
    close( fu, iostat=status )
    if (status/=0) then
      write (gol,'("closing file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! info ...
    write (gol,'(a,": create list of emisison variables to be read ..")') rname; call goPr
    ! storage for list of unique emisison variables that should be read:
    allocate( self%emvars(maxrow), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! no values yet:
    self%nemvar = 0
    ! loop over records:
    do i = 1, self%n
      ! any tracers defined?
      if ( self%trform(i)%n > 0 ) then
        ! assign tracer indices; return status:
        !   -2 : all undefined
        !   -1 : some undefined
        !    0 : all found (or none defined)
        !   >0 : error
        call self%trform(i)%AssignIndx( specname, status, quiet=.true. )
        if ( status > 0 ) then
          ! error ..
          TRACEBACK; status=1; return
        else if ( status == -2 ) then
          ! all undefined, skip
        else if ( status >= -1 ) then

          ! testing ...
          write (gol,'(a,":   emission `",a,"` ..")') rname, trim(self%name(i)); call goPr
          do j = 1, self%trform(i)%n
            if ( self%trform(i)%indx(j) < 0 ) then
              write (gol,'(a,":     (undefined: ",a,")")') rname, trim(self%trform(i)%name(j)); call goPr
            else
              write (gol,'(a,":     ",f5.2," ",a)') rname, self%trform(i)%w(j), trim(self%trform(i)%name(j)); call goPr
            end if
          end do
          
          ! loop over source variables:
          do j = 1, self%emform(i)%n
            ! current name:
            name = trim(self%emform(i)%name(j))
            ! could be present already?
            if ( self%nemvar > 0 ) then
              ! search, do not complain if not found:
              call goMatchValue( name, self%emvars(1:self%nemvar), indx, status, quiet=.true. )
              if ( status < 0 ) then
                ! not found
              else if ( status == 0 ) then
                ! found, continue with next:
                cycle
              else
                TRACEBACK; status=1; return
              end if
            end if
            ! increase counter:
            self%nemvar = self%nemvar + 1
            ! check ...
            if ( self%nemvar > size(self%emvars) ) then
              write (gol,'("number of unique emission variables (",i0,") exceeds storage (",i0,")")') &
                        self%nemvar, size(self%emvars); call goErr
              TRACEBACK; status=1; return
            end if
            ! store:
            self%emvars(self%nemvar) = trim(name)
          end do ! j

        end if ! some or all tracers enabled in model
        
      end if ! tracer assignment defined
    end do ! records
    
    ! storage for mapping to records with data:
    allocate( self%irec(maxrow), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! no records selected yet:
    self%nd = 0    
    ! any variables to be read?
    if ( self%nemvar == 0 ) then
      ! testing ...
      write (gol,'(a,": no emission variables selected ...")') rname; call goPr
    else
      ! testing ...
      write (gol,'(a,": selected emission variables:")') rname; call goPr
      do i = 1, self%nemvar
        write (gol,'(a,":   ",i3," ",a)') rname, i, trim(self%emvars(i)); call goPr
      end do
      ! testing ...
      write (gol,'(a,": emission data sets:")') rname; call goPr
      ! loop over records:
      do i = 1, self%n
        ! assign emission variable indices to elements of formula:
        call self%emform(i)%AssignIndx( self%emvars(1:self%nemvar), status, quiet=.true. )
        if ( status > 0 ) then
          ! error
          TRACEBACK; status=1; return
        else if ( status == 0 ) then
          ! this record is composed form emissions that are read:
          ! increase counter:
          self%nd = self%nd + 1
          ! fill index:
          self%irec(self%nd) = i
          ! testing ...
          write (gol,'(a,":   data set ",i4," defined by record ",i4)') rname, self%nd, i; call goPr
        else
          ! no emisison variables assign, so no need to use this record ...
        end if
      end do ! i
    end if ! nemvar > 0
   
    ! ok
    status = 0

  end subroutine Composition_Init
  
  
  ! ***
  

  subroutine Composition_Done( self, status )
  
    ! --- in/out ------------------------------
    
    class(T_Composition), intent(inout)    ::  self
    integer, intent(out)                   ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Composition_Done'
    
    ! --- local -------------------------------
    
    integer     ::  i
    
    ! --- begin -------------------------------
    
    ! loop over variables:
    do i = 1, self%n
      ! done:
      call self%emform(i)%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! done:
      call self%trform(i)%Done( status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! clear:
    deallocate( self%name, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%molemass, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%emform, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%trform, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( self%emvars, stat=status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    deallocate( self%irec, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Composition_Done
  
  
end module LE_Emis_GFAS
