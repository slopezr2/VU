!###############################################################################
!
! LE_Bound_EMEP_Monthly - LOTOS-EUROS boundary conditions from CIFS
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

module LE_Bound_EMEP_Monthly

  use GO, only : gol, goErr, goPr
  !use GO,  only: gomem
  use GO, only : TDate
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_LE_Bound_EMEP_Monthly
  public  ::  LE_Bound_EMEP_Monthly_Init, LE_Bound_EMEP_Monthly_Done
  public  ::  LE_Bound_EMEP_Monthly_Get
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_EMEP_Monthly'
  
  integer, parameter     :: max_archives = 10


  ! --- local --------------------------------
  
  type T_EMEP_Monthly_Arch
    type(TDate)               ::  tstart
    character(len=512)        ::  path
    integer                   ::  fixyear
  end type T_EMEP_Monthly_Arch
  
  type T_EMEP_Monthly_Field
    ! varable name:
    character(len=32)   ::  vname
    ! filled?
    logical             ::  filled
    ! time range:
    type(TDate)         ::  t1, t2
    ! axis:
    integer             ::  nlon, nlat, nlay
    real, allocatable   ::  lons(:)   ! (nlon)
    real, allocatable   ::  lats(:)   ! (nlat)
    ! data:
    real, allocatable   ::  data(:,:,:)  ! (lon,lat,lev)
    real, allocatable   ::  phlev(:,:,:)  ! (lon,lat,0:lev)
    real, allocatable   ::  ahlev(:,:,:)  ! (lon,lat,0:lev)
    character(len=32)   ::  units
  end type T_EMEP_Monthly_Field
  
  ! *
  
  type T_LE_Bound_EMEP_Monthly
    integer                                   ::  n
    integer                                   ::  current
    type(T_EMEP_Monthly_Arch), allocatable    ::  arch(:)
    ! number of variables:
    integer                                   ::  nvar
    ! maping from model ispec to variable:
    integer, allocatable                      ::  ivar(:)  ! (nspec)
    ! fields read:
    type(T_EMEP_Monthly_Field), allocatable   ::  f(:)  ! (nvar)
  end type T_LE_Bound_EMEP_Monthly
  

contains


  ! ===============================================================
  

  subroutine LE_Bound_EMEP_Monthly_Init( self, rcfile, rckey, status )
  
    use GO, only : TrcFile, Init, Done, ReadRc
    !use GO, only: gomem
    use GO, only : goVarValue, goSplitString, goTranslate
    use GO, only : operator(>=), NewDate

    use Dims   , only : nspec
    use Indices, only : i_dust_f, i_dust_c

    ! --- in/out ------------------------------
    
    type(T_LE_Bound_EMEP_Monthly), intent(out) ::  self
    character(len=*), intent(in)        ::  rcfile
    character(len=*), intent(in)        ::  rckey
    integer, intent(out)                ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Bound_EMEP_Monthly_Init'
    
    ! --- local ------------------------------
    
    type(TRcFile)         ::  rcF
    character(len=512)    ::  keys
    character(len=512)    ::  archives(max_archives)

    integer               ::  n_arch, iarch
    character(len=32)     ::  stime
    integer               ::  yy, mm, dd
    integer               ::  fixyear
    integer               ::  ivar
    
    ! --- begin -------------------------------
    
    ! open rcfile:
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)

    ! read from rc file: path, tstart, expid, dhour, allow_missing
    ! e.g. le.bound.macc-cifs.keys   :  path=${MODELS}/MACC/CIFS;tstart=2014-09-18;expid=D1_raq;dhour=3;allow_missing=F 
    !    note: there is one such line for each archive, with it's own starting date

    call ReadRc( rcF, trim(rckey)//'.keys', keys, status )
    IF_NOTOK_RETURN(status=1)
    
    call goSplitString( trim(keys), n_arch, archives, status )
    IF_NOTOK_RETURN(status=1)
    
    self%n = n_arch
    self%current = -1
    allocate( self%arch(n_arch) )
    
    ! Fill properties for each archive
    do iarch = 1, n_arch

      ! start time, should be present:
      call goVarValue( archives(iarch), ';', 'tstart', '=', stime, status )
      IF_NOTOK_RETURN(status=1)
      ! replace seperation characters by white space:
      call goTranslate( stime, '/-:', ' ', status )
      IF_NOTOK_RETURN(status=1)
      ! read time elements:
      read (stime,*,iostat=status) yy, mm, dd
      if ( status /= 0 ) then
        write (gol,'("could not read start time from : ",a, " in EMEP_Monthly boundary archives")') trim(stime); call goErr
        TRACEBACK; status=1; return
      end if
      ! fill:
      self%arch(iarch)%tstart = NewDate( year=yy, month=mm, day=dd )
      
      ! input path, should be present:
      call goVarValue( archives(iarch), ';', 'path', '=', self%arch(iarch)%path, status )
      IF_NOTOK_RETURN(status=1)

      ! default value:
      self%arch(iarch)%fixyear = -999
      ! reset if defined; warning status returned if not defined, so only check on error:
      call goVarValue( archives(iarch), ';', 'fixyear', '=', self%arch(iarch)%fixyear, status )
      IF_ERROR_RETURN(status=1)

    end do

    ! storage:
    allocate( self%ivar(nspec) )
    ! by default no variable for any spec:
    self%ivar = -999
    
    ! mapping:
    ! now hardcoded, in future from rcfile settings ...
    self%nvar = 2
    ! storage:
    allocate( self%f(self%nvar) )
    ! loop:
    do ivar = 1, self%nvar
      ! switch:
      select case ( ivar )
        case ( 1 ) 
          self%ivar(i_dust_f) = ivar
          self%f(ivar)%vname = 'DUST_WB_F'
        case ( 2 ) 
          self%ivar(i_dust_c) = ivar
          self%f(ivar)%vname = 'DUST_WB_C'
        case default
          write (gol,'("unsupported ivar ",i0)') ivar; call goErr
          TRACEBACK; status=1; return
      end select
      ! no data yet:
      self%f(ivar)%filled = .false.
    end do

    ! close rcfile:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Bound_EMEP_Monthly_Init
  
  
  ! ***
  

  subroutine LE_Bound_EMEP_Monthly_Done( self, status )
  
    ! --- in/out ------------------------------
    
    type(T_LE_Bound_EMEP_Monthly), intent(inout) ::  self
    integer, intent(out)                         ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Bound_EMEP_Monthly_Done'
    
    ! --- local -------------------------------
    
    integer    ::  ivar
    
    ! --- begin -------------------------------
    
    ! clear:
    deallocate( self%arch )
    
    ! loop:
    do ivar = 1, self%nvar
      ! data filled ?
      if ( self%f(ivar)%filled ) then
        ! clear:
        deallocate( self%f(ivar)%lons )
        deallocate( self%f(ivar)%lats )
        deallocate( self%f(ivar)%data )
        deallocate( self%f(ivar)%phlev )
        deallocate( self%f(ivar)%ahlev )
      end if
    end do
    ! clear
    deallocate( self%f )
    deallocate( self%ivar )
    
    
    ! ok
    status = 0

  end subroutine LE_Bound_EMEP_Monthly_Done
  

  ! ***
  

  !
  ! fields are read for time t
  ! pointers are allocated on output
  !
  ! return status:
  !   -2   : file not found
  !   -1   : variable not found
  !    0   : ok
  !   else : error
  !
  
  subroutine LE_Bound_EMEP_Monthly_Get( self, t, nhour, ispec, &
                  lons, lats, phlev, ahlev, vmr, vmr_unit, &
                  moleweight_tracer, moleweight_air, status, &
                  verbose )
  
    use GO, only : TDate, wrtgol, Get_Begin_Of, Get_End_Of
    use GO, only : operator(<), operator(==)
    
    use EMEP_Monthly, only : T_EMEP_Monthly
    use EMEP_Monthly, only : EMEP_Monthly_Init, EMEP_Monthly_Done
    use EMEP_Monthly, only : EMEP_Monthly_Get, EMEP_Monthly_GetData
    
    use indices, only : specname
    use dims   , only : nspec
    use LE_Bound_Tools, only : PTQ_to_H

    ! --- in/out ------------------------------
    
    type(T_LE_Bound_EMEP_Monthly), intent(inout) ::  self
    type(TDate), intent(in)         ::  t
    integer, intent(in)             ::  nhour    ! hours since start of forecst (or midnight)
    integer, intent(in)             ::  ispec
    real, pointer                   ::  lons(:), lats(:)
    real, pointer                   ::  phlev(:,:,:)
    real, pointer                   ::  ahlev(:,:,:)
    real, pointer                   ::  vmr(:,:,:)
    character(len=*), intent(out)   ::  vmr_unit
    real, intent(out)               ::  moleweight_tracer
    real, intent(out)               ::  moleweight_air
    integer, intent(out)            ::  status
    logical, intent(in), optional   ::  verbose
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Bound_EMEP_Monthly_Get'
    
    ! --- local -------------------------------
    
    logical                   ::  do_verbose
    integer                   ::  ivar
    logical                   ::  do_read
    type(T_EMEP_Monthly)      ::  emp
    real, allocatable         ::  temper(:,:,:)
    real, allocatable         ::  humid(:,:,:)

    ! --- begin -------------------------------

    ! messages ?
    do_verbose = .false.
    if ( present(verbose) ) do_verbose = verbose
    
    !! testing ...
    !do_verbose = .true.
    
    ! check ...
    if ( (ispec < 1) .or. (ispec > nspec) ) then
      write (gol,'("tracer index not in range:")'); call goErr
      write (gol,'("  ispec   : ",i6)') ispec; call goErr
      write (gol,'("  nspec   : ",i6)') nspec; call goErr
      TRACEBACK; status=1; return
    end if

    ! data provided by this set ?
    if ( self%ivar(ispec) < 0 ) then
      ! not available ...
      status=-1; return
    end if
    
    ! info ...
    if ( do_verbose ) then
      call wrtgol( rname//': get '//trim(specname(ispec))//' for time: ', t ); call goPr
    end if

    ! current:
    ivar = self%ivar(ispec)
    ! already data read?
    if ( self%f(ivar)%filled ) then
      ! target time outside interval for which data is valid?
      if ( (t < self%f(ivar)%t1) .or. (self%f(ivar)%t2 < t) ) then
        ! need to read:
        do_read = .true.
      else
        ! data already valid:
        do_read = .false.
      end if
    else
      ! not data yet, so read:
      do_read = .true.
    end if
    
    ! read data ?
    if ( do_read ) then

      ! info ...
      if ( do_verbose ) then
        write (gol,'(a,":   read new field for ",a)') rname, trim(specname(ispec)); call goPr
      end if

      ! new archive ?
      if ( self%current < 0 ) self%current = 1
      ! switch archive ?
      do
        ! last archive ?
        if ( self%current == self%n ) exit
        ! check if next should be used ...
        if ( t < self%arch(self%current+1)%tstart ) then
          ! keep current number
          exit
        end if
        ! try next:
        self%current = self%current + 1
      end do

      ! info ...
      if ( do_verbose ) then
        call wrtgol( rname//': selected archive for time: ', t ); call goPr
        write (gol,'(a,":   archive number : ",i0)') rname, self%current; call goPr
        write (gol,'(a,":   path           : ",a)') rname, trim(self%arch(self%current)%path); call goPr
        write (gol,'(a,":   fixyear        : ",i0)') rname, self%arch(self%current)%fixyear; call goPr
        write (gol,'(a,":   variable       : ",a)') rname, trim(self%f(ivar)%vname); call goPr
      end if
      !call goMem( 'MEMORY - monthly_init:', status )
      !IF_NOTOK_RETURN(status=1)
      ! open file:
      call EMEP_Monthly_Init( emp, self%arch(self%current)%path, 'DUST', t, status, &
                                fixyear=self%arch(self%current)%fixyear, &
                                verbose=do_verbose )
      IF_NOTOK_RETURN(status=1)
      !call goMem( 'MEMORY - monthly get:', status )
      !IF_NOTOK_RETURN(status=1)
      ! dims:
      call EMEP_MOnthly_Get( emp, status, &
                               nlon=self%f(ivar)%nlon, &
                               nlat=self%f(ivar)%nlat, &
                               nlay=self%f(ivar)%nlay )
      IF_NOTOK_RETURN(status=1)
      
      ! clear existing:
      if ( self%f(ivar)%filled ) then
        deallocate( self%f(ivar)%lons )
        deallocate( self%f(ivar)%lats )
        deallocate( self%f(ivar)%data )
        deallocate( self%f(ivar)%phlev )
        deallocate( self%f(ivar)%ahlev )
      end if
      !call goMem( 'MEMORY - allocate:', status )
      !IF_NOTOK_RETURN(status=1)
      ! storage:
      allocate( self%f(ivar)%lons(self%f(ivar)%nlon) )
      allocate( self%f(ivar)%lats(self%f(ivar)%nlat) )
      allocate( self%f(ivar)%data(self%f(ivar)%nlon,self%f(ivar)%nlat,self%f(ivar)%nlay) )
      allocate( self%f(ivar)%phlev(self%f(ivar)%nlon,self%f(ivar)%nlat,0:self%f(ivar)%nlay) )
      allocate( self%f(ivar)%ahlev(self%f(ivar)%nlon,self%f(ivar)%nlat,0:self%f(ivar)%nlay) )
      !call goMem( 'MEMORY - montlhyl get axis:', status )
      !IF_NOTOK_RETURN(status=1)
      ! axis:
      call EMEP_Monthly_Get( emp, status, lons=self%f(ivar)%lons, lats=self%f(ivar)%lats )
      IF_NOTOK_RETURN(status=1)
      !call goMem( 'MEMORY - montlhyl get data', status )
      !IF_NOTOK_RETURN(status=1)
      ! data:
      call EMEP_Monthly_GetData( emp, trim(self%f(ivar)%vname), &
                               self%f(ivar)%data, self%f(ivar)%units, &
                               self%f(ivar)%phlev, status )
      IF_NOTOK_RETURN(status=1)
      !call goMem( 'MEMORY - montlhyl done', status )
      !IF_NOTOK_RETURN(status=1)
      ! close:
      call EMEP_Monthly_Done( emp, status )
      IF_NOTOK_RETURN(status=1)
      !call goMem( 'MEMORY - montlhyl done2', status )
      !IF_NOTOK_RETURN(status=1)
      ! storage:
      allocate( temper(self%f(ivar)%nlon,self%f(ivar)%nlat,self%f(ivar)%nlay) )
      allocate( humid (self%f(ivar)%nlon,self%f(ivar)%nlat,self%f(ivar)%nlay) )
      ! no meteo, fill arrays with default values:
      temper = 293.0  ! K
      humid  = 0.0
      ! estimate half level altitudes:
      call PTQ_to_H( self%f(ivar)%phlev, temper, humid, self%f(ivar)%ahlev, status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( temper )
      deallocate( humid )
      
      ! field is valid for month:
      self%f(ivar)%t1 = Get_Begin_Of( t, 'month' )
      self%f(ivar)%t2 = Get_End_Of( t, 'month' )
      ! set flag:
      self%f(ivar)%filled = .true.
      
    else

      ! info ...
      if ( do_verbose ) then
        write (gol,'(a,":   keep current ",a)') rname, trim(specname(ispec)); call goPr
      end if
      
    end if  ! do_read
    
    ! clear if necessary:
    if ( associated(lons ) ) deallocate( lons  )
    if ( associated(lats ) ) deallocate( lats  )
    if ( associated(vmr  ) ) deallocate( vmr   )
    if ( associated(phlev) ) deallocate( phlev )
    if ( associated(ahlev) ) deallocate( ahlev )
    ! storage:
    allocate( lons(self%f(ivar)%nlon) )
    allocate( lats(self%f(ivar)%nlat) )
    allocate( vmr(self%f(ivar)%nlon,self%f(ivar)%nlat,self%f(ivar)%nlay) )
    allocate( phlev(self%f(ivar)%nlon,self%f(ivar)%nlat,0:self%f(ivar)%nlay) )
    allocate( ahlev(self%f(ivar)%nlon,self%f(ivar)%nlat,0:self%f(ivar)%nlay) )
    ! copy:
    lons  = self%f(ivar)%lons
    lats  = self%f(ivar)%lats
    vmr   = self%f(ivar)%data
    phlev = self%f(ivar)%phlev
    ahlev = self%f(ivar)%ahlev
    vmr_unit = trim(self%f(ivar)%units)
    ! not relevant for aerosols ...but EMEP gives mixing ratio, 200 g/mole
    !moleweight of dry air =29 g/mole
    moleweight_tracer = 200
    moleweight_air    = 29 
    
    !call goMem( 'MEMORY - after LE-bound_storage:', status )
    !IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Bound_EMEP_Monthly_Get
  

end module LE_Bound_EMEP_Monthly

