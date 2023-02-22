!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_Data

  use GO, only : gol, goPr, goErr

  use LE_Data_Variable, only : LEN_VARNAME

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Emis_Data


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Data'
  

  ! --- types --------------------------------
  
  ! hourly profile:
  type T_Profile_Hourly
    ! factors for local times [00,01], [01,02], ...
    real                                ::  factor(24)
  contains
    procedure :: Init          => Profile_Hourly_Init
    procedure :: Done          => Profile_Hourly_Done
    procedure :: GetFactor     => Profile_Hourly_GetFactor
  end type T_Profile_Hourly
  
  ! *

  ! emission data base
  type T_Emis_Data
    ! list with data variables:
    integer                               ::  nvar
    character(len=LEN_VARNAME)            ::  varnames(10)   ! just a maximum ...
    ! target spec indices:
    integer, allocatable                  ::  ispecs(:)   ! (nvar)
    ! hourly profiles:
    type(T_Profile_Hourly), allocatable   ::  profile_hourly(:)  ! (nvar)
  contains
    procedure :: Init          => Emis_Data_Init
    procedure :: Done          => Emis_Data_Done
    procedure :: Setup         => Emis_Data_Setup
  end type T_Emis_Data



contains


  ! ====================================================================
  ! ===
  ! === profile hourly
  ! ===
  ! ====================================================================
  
  !
  ! Example content:
  !
  !   LT_start;LT_end;factor
  !   1;2;0,008317338451695
  !   2;3;0,003710812539987
  !   3;4;0,009596928982726
  !   :
  !

  subroutine Profile_Hourly_Init( self, profile, status )
  
    use GO     , only : TTextFile, Init, Done, ReadLine
    use GO     , only : goReadFromLine
    
    ! --- in/out ---------------------------------
    
    class(T_Profile_Hourly), intent(out)      ::  self
    character(len=*), intent(in)              ::  profile
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Profile_Hourly_Init'
      
    ! --- local ----------------------------------
    
    type(TTextFile)               ::  F
    character(len=1)              ::  sep
    character(len=1024)           ::  line
    character(len=32)             ::  header
    integer                       ::  LT_start, LT_end
    integer                       ::  ihour
    real                          ::  factor
    integer                       ::  i
    real                          ::  fsum, fave
    
    ! --- begin ----------------------------------
    
    ! default?
    if ( len_trim(profile) == 0 ) then

      ! no factors ..
      self%factor = 1.0

    else

      ! init commented textfile:
      call Init( F, profile, status, status='old', comment='#' )
      IF_NOTOK_RETURN(status=1)

      ! formatting:
      sep = ';'
      ! header line:
      call ReadLine( F, line, status )
      IF_NOTOK_RETURN(status=1)
      ! check on expected header ...
      call goReadFromLine( line, header, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      if ( trim(header) /= 'LT_start' ) then
        write (gol,'("unexpected header `",a,"`")') trim(header); call goErr
        TRACEBACK; status=1; return
      end if
      ! check on expected header ...
      call goReadFromLine( line, header, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      if ( trim(header) /= 'LT_end' ) then
        write (gol,'("unexpected header `",a,"`")') trim(header); call goErr
        TRACEBACK; status=1; return
      end if
      ! check on expected header ...
      call goReadFromLine( line, header, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      if ( trim(header) /= 'factor' ) then
        write (gol,'("unexpected header `",a,"`")') trim(header); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! init factors:
      self%factor = -999
      ! loop over records:
      do
        ! read next line:
        call ReadLine( F, line, status )
        if ( status < 0 ) exit  ! eof
        IF_NOTOK_RETURN(status=1)
        
        ! read:
        call goReadFromLine( line, lt_start, status, sep=sep )
        IF_NOTOK_RETURN(status=1)
        call goReadFromLine( line, lt_end, status, sep=sep )
        IF_NOTOK_RETURN(status=1)
        ! check ...
        if ( lt_end /= lt_start+1 ) then
          write (gol,'("invalid LT interval [",i0,",",i0,"]")') lt_start, lt_end; call goErr
          TRACEBACK; status=1; return
        end if
        ! index is end hour:
        ihour = lt_end
        ! check ...
        if ( self%factor(ihour) >= 0 ) then
          write (gol,'("record ",i0," already filled")') ihour; call goErr
          TRACEBACK; status=1; return
        end if
        
        ! read:
        call goReadFromLine( line, factor, status, sep=sep )
        IF_NOTOK_RETURN(status=1)
        ! check ...
        if ( factor < 0.0 ) then
          write (gol,'("found negative factor in file: ",a)') trim(profile); call goErr
          TRACEBACK; status=1; return
        end if
        print *, 'xxx', lt_start, lt_end, factor
        ! fill:
        self%factor(ihour) = factor
      
      end do ! records
      
      ! close:
      call Done( F, status )
      IF_NOTOK_RETURN(status=1)
      
      ! check ...
      if ( any(self%factor < 0) ) then
        write (gol,'("not all factors defined")'); call goErr
        do i = 1, size(self%factor)
          write (gol,'("  ",i2," ",f12.6)') self%factor(i); call goErr
        end do
        TRACEBACK; status=1; return
      end if
      
      ! sum and average:
      fsum = sum(self%factor)
      fave = fsum / size(self%factor)
      ! if sum is close to 1, the factors are actually fractions ...
      if ( abs(fsum-1.0) < 1.0e-4 ) then
        ! convert to deviations around nominal value:
        self%factor = self%factor / fave
        ! new sum and average:
        fsum = sum(self%factor)
        fave = fsum / size(self%factor)
      end if
      ! check ...
      if ( abs(fave-1.0) > 1.0e-4 ) then
        write (gol,'("average of factors should be 1.0, found: ",es16.10)') fave; call goErr
        TRACEBACK; status=1; return
      end if
      ! scale to have average of exactly 1.0:
      self%factor = self%factor / fave

    end if  ! default or csv
    
    ! ok
    status = 0
    
  end subroutine Profile_Hourly_Init


  ! ***


  subroutine Profile_Hourly_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Profile_Hourly), intent(inout)    ::  self
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Profile_Hourly_Done'
      
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
    
  end subroutine Profile_Hourly_Done


  ! ***


  subroutine Profile_Hourly_GetFactor( self, lhour, factor, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Profile_Hourly), intent(inout)    ::  self
    integer, intent(in)                       ::  lhour(:,:)  ! local hour 0,..,23
    real, intent(out)                         ::  factor(:,:)
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Profile_Hourly_GetFactor'
      
    ! --- local ----------------------------------
    
    integer           ::  i, j
    integer           ::  irec
    
    ! --- begin ----------------------------------
    
    ! loop over cells:
    do i = 1, size(lhour,1)
      do j = 1, size(lhour,2)
        ! record:
        irec = lhour(i,j) + 1  ! 1,..,24
        ! copy:
        factor(i,j) = self%factor(irec)
      end do ! j
    end do ! i
    
    ! ok
    status = 0
    
  end subroutine Profile_Hourly_GetFactor


  ! ====================================================================
  ! ===
  ! === Emis Data
  ! ===
  ! ====================================================================


  subroutine Emis_Data_Init( self, rcF, rckey, status )
  
    use GO     , only : TrcFile, ReadRc
    use GO     , only : goSplitString, goMatchValue
    use LE_Data, only : LE_Data_Enable
    use Indices, only : specname
    
    ! --- in/out ---------------------------------
    
    class(T_Emis_Data), intent(out)           ::  self
    type(TrcFile), intent(in)                 ::  rcF
    character(len=*), intent(in)              ::  rckey
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Emis_Data_Init'
      
    ! --- local ----------------------------------
    
    character(len=1024)           ::  line
    integer                       ::  ivar
    character(len=32)             ::  spec
    character(len=1024)           ::  profile
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("setup emissions from data variables ...")'); call goPr

    ! enable data:
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lon', status )
    IF_NOTOK_RETURN(status=1)
    
    ! read line with supported variables:
    call ReadRc( rcF, trim(rckey)//'.vars', line, status )
    IF_NOTOK_RETURN(status=1)
    ! split:
    call goSplitString( line, self%nvar, self%varnames, status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( self%ispecs(self%nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%profile_hourly(self%nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over variables:
    do ivar = 1, self%nvar

      ! enable data:
      call LE_Data_Enable( trim(self%varnames(ivar)), status )
      IF_NOTOK_RETURN(status=1)

      ! target spec:
      call ReadRc( rcF, trim(rckey)//'.'//trim(self%varnames(ivar))//'.spec', spec, status )
      IF_NOTOK_RETURN(status=1)
      ! global index:
      call goMatchValue( spec , specname, self%ispecs(ivar), status )
      IF_NOTOK_RETURN(status=1)
      
      ! csv file with hourly profile, empty for none:
      call ReadRc( rcF, trim(rckey)//'.'//trim(self%varnames(ivar))//'.profile_hourly', profile, status )
      IF_NOTOK_RETURN(status=1)
      ! init:
      call self%profile_hourly(ivar)%Init( profile, status )
      IF_NOTOK_RETURN(status=1)

    end do
    
    ! ok
    status = 0
    
  end subroutine Emis_Data_Init


  ! ***


  subroutine Emis_Data_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Emis_Data), intent(inout)         ::  self
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Emis_Data_Done'
      
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%ispecs, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%profile_hourly, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Emis_Data_Done


  ! ***


  subroutine Emis_Data_Setup( self, t1, t2, status )
  
    use GO     , only : TDate
    use LE_Data, only : LE_Data_GetPointer
    use Dims   , only : nx, ny
    use dims   , only : emis_a
    use Indices, only : specname, specunit, specmolm
  
    ! --- in/out ---------------------------------
    
    class(T_Emis_Data), intent(inout)         ::  self
    type(TDate), intent(in)                   ::  t1, t2
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Emis_Data_Setup'
    
    ! conversions:
    real, parameter   ::  sec_per_min  = 60.0
    real, parameter   ::  min_per_year = 365.0 * 24.0 * 60.0
      
    ! --- local ----------------------------------
    
    integer                 ::  ivar
    character(len=64)       ::  units
    integer                 ::  ispec
    real, pointer           ::  area(:,:,:)
    real, pointer           ::  lon(:,:,:)
    real, pointer           ::  emis(:,:,:)
    character(len=256)      ::  conversion
    integer, allocatable    ::  zone(:,:)  ! (nx,ny)
    integer, allocatable    ::  lhour(:,:)   ! (nx,ny)
    real, allocatable       ::  factor(:,:)  ! (nx,ny)
    
    ! --- begin ----------------------------------

    ! grid cell area:
    call LE_Data_GetPointer( 'area', area, status, check_units='m2' )
    IF_NOTOK_RETURN(status=1)
    ! longitude, used for local time:
    call LE_Data_GetPointer( 'lon', lon, status, check_units='degrees_east' )
    IF_NOTOK_RETURN(status=1)
    
    ! hour 0,..,23 in local time based on 15 degrees lon band:
    !       -11                0   1                12
    !     +-----+--- ... ---+----+----+--- ... ---+----+
    !   -180  -165        -15    0   15         165   180
    ! storage:
    allocate( zone(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( lhour(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( factor(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! round 15 degree bands:
    zone = int(ceiling(lon(1:nx,1:ny,1)/15.0))
    ! add current hour in UTC, convert to 0,..,23:
    lhour = modulo( t1%hour + zone, 24 )
        
    ! loop over variables:
    do ivar = 1, self%nvar
    
      ! info ...
      write (gol,'(" add `",a,"` emissions ...")') trim(self%varnames(ivar)); call goPr

      ! pointer to emission array:
      call LE_Data_GetPointer( trim(self%varnames(ivar)), emis, status, &
                                 units=units )
      IF_NOTOK_RETURN(status=1)

      ! field with hour factors:
      call self%profile_hourly(ivar)%GetFactor( lhour, factor, status )
      IF_NOTOK_RETURN(status=1)

      ! index of target spec:
      ispec = self%ispecs(ivar)
      
      ! required conversion:
      select case ( specunit(ispec) )
        !~ emission for tracers expressed as volume mixing ratio's:
        case ( 'ppb' )
          ! unit emis_a for gases   : mol/min
          conversion = trim(units)//' -> mol/min'
        !~ emission for aerosol tracers:
        case ( 'ug/m3' )
           conversion = trim(units)//' -> ug/min'
        !~ unknown ...
        case default
          write (gol,'("unsupported unit `",a,"` for tracer ",i4," (",a,")")') &
                           trim(specunit(ispec)), ispec, trim(specname(ispec)); call goErr
          TRACEBACK; status=1; return
      end select

      ! fill emission array:
      select case ( trim(conversion) )

        !~ emission of trace gasses:
        case ( 'kg/m2/s -> mol/min' )
          !       mol/min                     mol/min          
          emis_a(1:nx,1:ny,1,ispec) = emis_a(1:nx,1:ny,1,ispec) &
                   !   kg/m2/s   *         m2        /     (kg/mol)    *  (s/min)    *   1
                   + emis(:,:,1) * area(1:nx,1:ny,1) / specmolm(ispec) * sec_per_min * factor
	
        !~ emission of aerosols:
        case ( 'Mg/year -> ug/min' )
          !        ug/min                      ug/min
          emis_a(1:nx,1:ny,1,ispec) = emis_a(1:nx,1:ny,1,ispec) &
                   !  Mg/year    *  ug/Mg  /  (min/year)  *    1
                   + emis(:,:,1) * 1.0e-12 / min_per_year * factor

        !~ unknown ...
        case default
          write (gol,'("unsupported conversion `",a,"` for tracer ",i4," (",a,") emitted from `",a,"`")') &
                           trim(conversion), ispec, trim(specname(ispec)), &
                           trim(self%varnames(ivar)); call goErr
          TRACEBACK; status=1; return
      end select
      
    end do  ! variables
    
    ! clear:
    deallocate( zone, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( lhour, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( factor, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Emis_Data_Setup


end module LE_Emis_Data

