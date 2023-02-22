!###############################################################################
!
! NAME
!
!   LE_Output  -  LOTOS-EUROS output
!
! DESCRIPTION
!
!   Interface to LE output of various types.
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

module LE_Output

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  use LE_Output_Grid         , only : T_LE_Output_Grid
  use LE_Output_Dat          , only : T_LE_Output_Dat
  use LE_Output_Dat2         , only : T_LE_Output_Dat2
  use LE_Output_conc         , only : T_LE_Output_conc
#ifdef with_labeling
  use LE_Output_Label        , only : T_LE_Output_Label
  use LE_Output_Labelled_conc, only : T_LE_Output_Labelled_conc
  use LE_Output_Emis_Label   , only : T_LE_Output_Emis_Label
  use LE_Output_drydep_Label , only : T_LE_Output_drydep_Label
  use LE_Output_wetdep_Label , only : T_LE_Output_wetdep_Label
  use LE_Output_budget_Label , only : T_LE_Output_budget_Label
#endif    
#ifdef with_radiation
  use LE_Output_rad          , only : T_LE_Output_rad
#endif
  use LE_Output_emis         , only : T_LE_Output_emis
  use LE_Output_vd           , only : T_LE_Output_vd
  use LE_Output_vd_diag      , only : T_LE_Output_vd_diag
  use LE_Output_drydep       , only : T_LE_Output_drydep
  use LE_Output_wetdep       , only : T_LE_Output_wetdep
  use LE_Output_budget       , only : T_LE_Output_budget
  use LE_Output_column       , only : T_LE_Output_column
  use LE_Output_Common       , only : T_LE_Output_Common
  use LE_Output_CSO          , only : T_LE_Output_CSO_Data, T_LE_Output_CSO_State
  
  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  T_LE_Output, T_LE_OutputState
  public  ::  T_LE_Output_Dat
  public  ::  T_LE_Output_Dat2
  public  ::  T_LE_Output_conc
#ifdef with_labeling
  public  ::  T_LE_Output_Label
  public  ::  T_LE_Output_Labelled_conc
  public  ::  T_LE_Output_Emis_Label
  public  ::  T_LE_Output_Drydep_Label
  public  ::  T_LE_Output_Wetdep_Label
  public  ::  T_LE_Output_Budget_Label
#endif    
#ifdef with_radiation
  public  ::  T_LE_Output_rad
#endif
  public  ::  T_LE_Output_emis
  public  ::  T_LE_Output_vd
  public  ::  T_LE_Output_vd_diag
  public  ::  T_LE_Output_drydep
  public  ::  T_LE_Output_wetdep
  public  ::  T_LE_Output_budget
  public  ::  T_LE_Output_column
  public  ::  T_LE_Output_Common
  public  ::  T_LE_Output_CSO_Data, T_LE_Output_CSO_State

  public  ::  Init, Done
  public  ::  Setup
  public  ::  PutOut
  public  ::  EndPutOut


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Output'

  ! maximum number of output files defined in rcfile:
  integer, parameter    ::  maxfile = 10


  ! --- types ------------------------------

  type T_LE_Output
    ! settings
    character(len=64)                     ::  rckey
    ! setup interval:
    type(TDate)                           ::  t1, t2
    ! lockfile:
    logical                               ::  lockfile_apply
    real                                  ::  lockfile_dhour
    character(len=1024)                   ::  lockfile_template
    ! grid output:
    type(T_LE_Output_Grid)                ::  file_grid
    ! dat output
    integer                               ::  ndat
    character(len=16)                     ::  file_dat(maxfile)
    type(T_LE_Output_Dat), pointer        ::  output_dat(:)
    ! dat2 output
    integer                               ::  ndat2
    character(len=16)                     ::  file_dat2(maxfile)
    type(T_LE_Output_Dat2), pointer       ::  output_dat2(:)
    ! conc output
    integer                               ::  nconc
    character(len=16)                     ::  file_conc(maxfile)
#ifdef with_labeling
    integer                               ::  nlabel
    character(len=16)                     ::  file_label(maxfile)   

    integer                               ::  nlabelled_conc
    character(len=32)                     ::  file_labelled_conc(maxfile)   

    integer                               ::  nemislabel
    character(len=16)                     ::  file_emis_label(maxfile)   
    
    integer                               ::  ndrydeplabel
    character(len=16)                     ::  file_drydep_label(maxfile)   
    
    integer                               ::  nwetdeplabel
    character(len=16)                     ::  file_wetdep_label(maxfile) 
    
    integer                               ::  nbudgetlabel
    character(len=16)                     ::  file_budget_label(maxfile) 
#endif
#ifdef with_radiation
    integer                               ::  nrad
    character(len=16)                     ::  file_rad(maxfile)
#endif    
    ! emis output
    integer                               ::  nemis
    character(len=16)                     ::  file_emis(maxfile)
    ! vd output
    integer                               ::  nvd
    character(len=16)                     ::  file_vd(maxfile)
    ! diagnostic vd output
    integer                               ::  nvd_diag
    character(len=16)                     ::  file_vd_diag(maxfile)
    ! drydep output
    integer                               ::  ndrydep
    character(len=16)                     ::  file_drydep(maxfile)
    ! wetdep output
    integer                               ::  nwetdep
    character(len=16)                     ::  file_wetdep(maxfile)
    ! budget output
    integer                               ::  nbudget
    character(len=16)                     ::  file_budget(maxfile)
    ! column output
    integer                               ::  ncolumn
    character(len=16)                     ::  file_column(maxfile)
    ! cso output
    integer                               ::  ncso
    character(len=16)                     ::  file_cso(maxfile)
    type(T_LE_Output_CSO_Data), pointer   ::  cso_data(:)
  end type T_LE_Output
  
  ! simulation per state:
  type T_LE_OutputState
    character(len=32)                     ::  name
    ! conc output
    type(T_LE_Output_conc), pointer       ::  output_conc(:)
#ifdef with_labeling
    ! label output
    type(T_LE_Output_Label), pointer      ::  output_Label(:)
    ! labelled conc output
    type(T_LE_Output_Labelled_Conc), pointer  ::  output_Labelled_conc(:)
    ! emis label output
    type(T_LE_Output_Emis_Label), pointer     ::  output_emis_Label(:)
    ! drydep label output
    type(T_LE_Output_drydep_Label), pointer   ::  output_drydep_Label(:)
    ! wetdep label output
    type(T_LE_Output_wetdep_Label), pointer   ::  output_wetdep_Label(:)
    ! budget label output
    type(T_LE_Output_budget_Label), pointer   ::  output_budget_Label(:)
#endif    
#ifdef with_radiation
    ! radiation output
    type(T_LE_Output_rad), pointer        ::  output_rad(:)
#endif
    ! emis output
    type(T_LE_Output_emis), pointer       ::  output_emis(:)
    ! vd output:
    type(T_LE_Output_vd), pointer         ::  output_vd(:)
    ! diagnostic vd output:
    type(T_LE_Output_vd_diag), pointer    ::  output_vd_diag(:)
    ! drydep output
    type(T_LE_Output_drydep), pointer     ::  output_drydep(:)
    ! wetdep output
    type(T_LE_Output_wetdep), pointer     ::  output_wetdep(:)
    ! budget output
    type(T_LE_Output_budget), pointer     ::  output_budget(:)
    ! column output
    type(T_LE_Output_column), pointer     ::  output_column(:)
    ! cso output
    type(T_LE_Output_CSO_State), pointer  ::  cso_state(:)
  end type T_LE_OutputState


  ! --- interfaces -------------------------

  interface Init
    module procedure leo_Init
    module procedure leos_Init
  end interface

  interface Done
    module procedure leo_Done
    module procedure leos_Done
  end interface

  interface Setup
    module procedure leo_Setup
  end interface

  interface PutOut
    module procedure leo_PutOut
    module procedure leos_PutOut
  end interface

  interface EndPutOut
    module procedure leo_EndPutOut
  end interface


contains


  ! ====================================================================
  ! ===    
  ! === output data
  ! ===    
  ! ====================================================================


  subroutine leo_Init( leo, rcF, rckey, status )

    use GO     , only : TrcFile, ReadRc
    use GO     , only : goSplitString
    use GO     , only : TDate

    use LE_Output_Dat  , only : LE_Output_Dat_Init
    use LE_Output_Dat2 , only : LE_Output_Dat2_Init

    ! --- in/out --------------------------------

    type(T_LE_Output), intent(out)      ::  leo
    type(TrcFile), intent(in)           ::  rcF
    character(len=*), intent(in)        ::  rckey
    integer, intent(out)                ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_Init'

    ! --- local ---------------------------------

    character(len=512)              ::  list
    integer                         ::  i

    ! --- begin ---------------------------------

    ! store:
    leo%rckey  = trim(rckey)

    ! lockfile to be written ?
    call rcF%Get( trim(leo%rckey)//'.lockfile.apply', leo%lockfile_apply, status )
    IF_NOTOK_RETURN(status=1)
    if ( leo%lockfile_apply ) then
       call rcF%Get( trim(leo%rckey)//'.lockfile.dhour', leo%lockfile_dhour, status )
       IF_NOTOK_RETURN(status=1)
       call rcF%Get( trim(leo%rckey)//'.lockfile.template', leo%lockfile_template, status )
       IF_NOTOK_RETURN(status=1)
    end if
    
    ! create grid definition file:
    call leo%file_grid%Init( rcF, trim(leo%rckey), status )
    IF_NOTOK_RETURN(status=1)

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.dat.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%ndat, leo%file_dat, status )
      IF_NOTOK_RETURN(status=1)

      ! any dat output ?
      if ( leo%ndat > 0 ) then
        ! storage for output stuff:
        allocate( leo%output_dat(leo%ndat) )
        ! loop over files to be written:
        do i = 1, leo%ndat
          call LE_Output_Dat_Init( leo%output_dat(i), rcF, &
                                     rckey, 'dat', leo%file_dat(i), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end if

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.dat2.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%ndat2, leo%file_dat2, status )
      IF_NOTOK_RETURN(status=1)

      ! any dat output ?
      if ( leo%ndat2 > 0 ) then
        ! storage for output stuff:
        allocate( leo%output_dat2(leo%ndat2) )
        ! loop over files to be written:
        do i = 1, leo%ndat2
          call LE_Output_Dat2_Init( leo%output_dat2(i), rcF, &
                                      rckey, 'dat2', leo%file_dat2(i), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end if

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.conc.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%nconc, leo%file_conc, status )
      IF_NOTOK_RETURN(status=1)

#ifdef with_labeling
    ! read files:
    call rcF%Get( trim(leo%rckey)//'.label.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%nlabel, leo%file_label, status )
      IF_NOTOK_RETURN(status=1)   
  
    ! read files:
    call rcF%Get( trim(leo%rckey)//'.labelled_conc.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%nlabelled_conc, leo%file_labelled_conc, status )
      IF_NOTOK_RETURN(status=1)   
  
    ! read files
    call rcF%Get( trim(leo%rckey)//'.emis_label.files', list, status )
    IF_NOTOK_RETURN(status=1)
    
      ! extract parts:
      call goSplitString( list, leo%nemislabel, leo%file_emis_label, status )
      IF_NOTOK_RETURN(status=1)      
      
    ! read files
    call rcF%Get( trim(leo%rckey)//'.drydep_label.files', list, status )
    IF_NOTOK_RETURN(status=1)
    
      ! extract parts:
      call goSplitString( list, leo%ndrydeplabel, leo%file_drydep_label, status )
      IF_NOTOK_RETURN(status=1)      
      
    ! read files
    call rcF%Get( trim(leo%rckey)//'.wetdep_label.files', list, status )
    IF_NOTOK_RETURN(status=1)
    
      ! extract parts:
      call goSplitString( list, leo%nwetdeplabel, leo%file_wetdep_label, status )
      IF_NOTOK_RETURN(status=1)
    
    ! read files
    call rcF%Get( trim(leo%rckey)//'.budget_label.files', list, status )
    IF_NOTOK_RETURN(status=1)
      
      ! extract parts:
      call goSplitString( list, leo%nbudgetlabel, leo%file_budget_label, status )
      IF_NOTOK_RETURN(status=1)
#endif
!    ! read files:
!    call rcF%Get( trim(leo%rckey)//'.dep.files', list, status )
!    IF_NOTOK_RETURN(status=1)
!
!      ! extract parts:
!      call goSplitString( list, leo%ndep, leo%file_dep, status )
!      IF_NOTOK_RETURN(status=1)
!      
!      ! any dep output ?
!      if ( leo%ndep > 0 ) then
!        ! storage for output stuff:
!        allocate( leo%output_dep(leo%ndep) )
!        ! loop over files to be written:
!        do i = 1, leo%ndep
!          call LE_Output_dep_Init( leo%output_dep(i), rcfile, &
!                       rckey, 'dep',  leo%file_dep(i), status )
!          IF_NOTOK_RETURN(status=1)
!        end do
!      end if

#ifdef with_radiation
    ! read files:
    call rcF%Get( trim(leo%rckey)//'.rad.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%nrad, leo%file_rad, status )
      IF_NOTOK_RETURN(status=1)
#endif      

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.emis.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%nemis, leo%file_emis, status )
      IF_NOTOK_RETURN(status=1)

    ! * vd

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.vd.files', list, status )
    IF_NOTOK_RETURN(status=1)

    ! extract parts:
    call goSplitString( list, leo%nvd, leo%file_vd, status )
    IF_NOTOK_RETURN(status=1)

    ! * diagnostic vd

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.vd-diag.files', list, status )
    IF_NOTOK_RETURN(status=1)

    ! extract parts:
    call goSplitString( list, leo%nvd_diag, leo%file_vd_diag, status )
    IF_NOTOK_RETURN(status=1)

    ! *

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.drydep.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%ndrydep, leo%file_drydep, status )
      IF_NOTOK_RETURN(status=1)

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.wetdep.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%nwetdep, leo%file_wetdep, status )
      IF_NOTOK_RETURN(status=1)

    ! read files:
    call rcF%Get( trim(leo%rckey)//'.budget.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%nbudget, leo%file_budget, status )
      IF_NOTOK_RETURN(status=1)
      
    ! read files:
    call rcF%Get( trim(leo%rckey)//'.column.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%ncolumn, leo%file_column, status )
      IF_NOTOK_RETURN(status=1)

    ! read files:
    call ReadRc( rcF, trim(leo%rckey)//'.cso.files', list, status )
    IF_NOTOK_RETURN(status=1)

      ! extract parts:
      call goSplitString( list, leo%ncso, leo%file_cso, status )
      IF_NOTOK_RETURN(status=1)
    
      ! info ...
      write (gol,'(a,":   number of `cso` files: ",i0)') rname, leo%ncso; call goPr

      ! any cso output ?
      if ( leo%ncso > 0 ) then
        ! storage for output stuff:
        allocate( leo%cso_data(leo%ncso), stat=status )
        IF_NOTOK_RETURN(status=1)
        ! loop over files to be written:
        do i = 1, leo%ncso
          ! info ...
          write (gol,'(a,":     `cso` file: ",a)') rname, trim(leo%file_cso(i)); call goPr
          ! setup data:
          call leo%cso_data(i)%Init( rcF, trim(rckey), 'cso', trim(leo%file_cso(i)), status )
          IF_NOTOK_RETURN(status=1)
        end do  ! output files
      end if  ! any output?
    
    ! info ...
    write (gol,'(a,": end")') rname; call goPr

    ! ok
    status = 0

  end subroutine leo_Init


  ! ***


  subroutine leo_Done( leo, status )

    use LE_Output_dat  , only : LE_Output_Dat_Done
    use LE_Output_Dat2 , only : LE_Output_Dat2_Done

    ! --- in/out ---------------------------------

    type(T_LE_Output), intent(inout)    ::  leo
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_Done'

    ! --- local ----------------------------------

    integer                         ::  i

    ! --- begin ----------------------------------

    ! info ...
    write (gol,'(a,": clear ...")') rname; call goPr
    
    ! done grid definition file:
    call leo%file_grid%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! any dat output ?
    if ( leo%ndat > 0 ) then
      ! loop over files:
      do i = 1, leo%ndat
        ! done with output file:
        call LE_Output_Dat_Done( leo%output_dat(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leo%output_dat, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! any dat2 output ?
    if ( leo%ndat2 > 0 ) then
      ! loop over files:
      do i = 1, leo%ndat2
        ! done with output file:
        call LE_Output_Dat2_Done( leo%output_dat2(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leo%output_dat2, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! any cso output ?
    if ( leo%ncso > 0 ) then
      ! loop over files to be written:
      do i = 1, leo%ncso
        ! done with data:
        call leo%cso_data(i)%Done( status )
        IF_NOTOK_RETURN(status=1)
      end do  ! output files
      ! clear:
      deallocate( leo%cso_data, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if  ! any output?

    ! ok
    status = 0

  end subroutine leo_Done


  ! ***


  subroutine leo_Setup( self, t1, t2, status )

    use GO, only : TDate, Midnight, IncrDate, operator(+), wrtgol
    
    ! --- in/out --------------------------------

    type(T_LE_Output), intent(inout)    ::  self
    type(TDate), intent(in)             ::  t1, t2
    integer, intent(out)                ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_Setup'

    ! --- local ---------------------------------

    integer       ::  i
    
    ! --- begin ---------------------------------

    ! info ...
    call wrtgol( rname//': setup for ', (/t1,t2/) ); call goPr

    ! store interval:
    self%t1 = t1
    self%t2 = t2
    
    ! put out any cso simulations ?
    if ( self%ncso > 0 ) then
      ! loop:
      do i = 1, self%ncso
        ! info ...
        write (gol,'(a,":   setup cso `",a,"` ...")') rname, trim(self%file_cso(i)); call goPr
        ! load observation info for current interval:
        call self%cso_data(i)%Setup( (/t1,t2/), status )
        IF_NOTOK_RETURN(status=1)
      end do ! output files
    end if ! any output
    
    ! ok
    status = 0

  end subroutine leo_Setup


  ! ***


  subroutine leo_PutOut( self, t, status )

    use GO            , only : TDate, MidNight, wrtgol
    use dims          , only : nx, ny, nz, nspec
    use LE_Output_Dat , only : LE_Output_Dat_PutOut
    use LE_Output_Dat2, only : LE_Output_Dat2_PutOut

    ! --- in/out --------------------------------

    type(T_LE_Output), intent(inout)    ::  self
    type(TDate), intent(in)             ::  t
    integer, intent(out)                ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_PutOut'

    ! --- local ---------------------------------

    integer        ::  i

    ! --- begin ---------------------------------

    ! info ...
    call wrtgol( rname//': put out at ', t ); call goPr

    ! put out grid definition file:
    call self%file_grid%PutOut( status )
    IF_NOTOK_RETURN(status=1)

    ! put out dat files ?
    if ( self%ndat > 0 ) then
      do i = 1, self%ndat
        call LE_Output_Dat_PutOut( self%output_dat(i), t, status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! put out dat2 files ?
    if ( self%ndat2 > 0 ) then
      do i = 1, self%ndat2
        call LE_Output_Dat2_PutOut( self%output_dat2(i), t, status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! put out any cso simulations ?
    if ( self%ncso > 0 ) then
      ! loop:
      do i = 1, self%ncso
        ! info ...
        write (gol,'(a,":   put out cso `",a,"` ...")') rname, trim(self%file_cso(i)); call goPr
        ! load observation info for current interval:
        call self%cso_data(i)%PutOut( t, status )
        IF_NOTOK_RETURN(status=1)
        ! info ...
        write (gol,'(a,":     done")') rname; call goPr
      end do ! output files
    end if ! any output

    ! ok
    status = 0

  end subroutine leo_PutOut


  ! ***


  subroutine leo_EndPutOut( leo, t, status )

    use GO, only : TDate, Get, Precisely
    use GO, only : goGetFU
    use GO, only : goReplace

    ! --- in/out --------------------------------

    type(T_LE_Output), intent(inout)    ::  leo
    type(TDate), intent(in)             ::  t
    integer, intent(out)                ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/leo_EndPutOut'

    ! --- local ---------------------------------

    character(len=1024) ::  lockfile
    integer             ::  time6(6)
    character(len=8)    ::  ccyymmdd
    character(len=4)    ::  hhmm
    integer             ::  fu

    ! --- begin ---------------------------------

    ! write lockfile ?
    if ( leo%lockfile_apply ) then
       ! for multiples of dhour only ...
       if ( Precisely(t,leo%lockfile_dhour,'hour')  ) then
          ! copy from template:
          lockfile = trim(leo%lockfile_template)
          ! get time:
          call Get( t, time6=time6 )
          ! replace some keys:
          write (ccyymmdd,'(i4.4,2i2.2)') time6(1:3)
          write (hhmm,'(2i2.2)') time6(4:5)
          call goReplace( lockfile, '%{ccyymmdd}', ccyymmdd, status )
          IF_NOTOK_RETURN(status=1)
          call goReplace( lockfile, '%{hhmm}', hhmm, status )
          IF_NOTOK_RETURN(status=1)
          ! obtain free file unit:
          call goGetFU( fu, status )
          IF_NOTOK_RETURN(status=1)
          ! open:
          open( unit=fu, file=trim(lockfile), status='new', form='formatted', iostat=status )
          if ( status /= 0 ) then
             write (gol,'("writing lockfile `",a,"`; already present ?")') trim(lockfile); call goErr
             TRACEBACK; status=1; return
          end if
          ! write dummy text into file:
          write (fu,'("LOTOS-EUROS output written for ",i4,2("-",i2.2)," ",i2.2,2(":",i2.2))',iostat=status) time6
          IF_NOTOK_RETURN(status=1)
          ! close:
          close( fu, iostat=status )
          IF_NOTOK_RETURN(status=1)
       end if
    end if

    ! ok
    status = 0

  end subroutine leo_EndPutOut


  ! ====================================================================
  ! ===    
  ! === output state
  ! ===    
  ! ====================================================================


  subroutine leos_Init( leos, leo, rcF, name, status )

    use GO                     , only : TrcFile
    use LE_Output_Conc         , only : LE_Output_Conc_Init
#ifdef with_labeling
    use LE_Output_Label        , only : LE_Output_Label_Init
    use LE_Output_Labelled_conc, only : LE_Output_Labelled_conc_Init
    use LE_Output_Emis_Label   , only : LE_Output_Emis_Label_Init
    use LE_Output_drydep_Label , only : LE_Output_drydep_Label_Init
    use LE_Output_wetdep_Label , only : LE_Output_wetdep_Label_Init
    use LE_Output_budget_Label , only : LE_Output_budget_Label_Init
#endif        
#ifdef with_radiation
    use LE_Output_Rad    , only : LE_Output_Rad_Init
#endif
    use LE_Output_Emis   , only : LE_Output_Emis_Init
    use LE_Output_vd     , only : LE_Output_vd_Init
    use LE_Output_vd_diag, only : LE_Output_vd_diag_Init
    use LE_Output_Drydep , only : LE_Output_Drydep_Init
    use LE_Output_Wetdep , only : LE_Output_Wetdep_Init
    use LE_Output_budget , only : LE_Output_budget_Init
    use LE_Output_Column , only : LE_Output_Column_Init
    
    ! --- in/out --------------------------------

    type(T_LE_OutputState), intent(out)   ::  leos
    type(T_LE_Output), intent(in)         ::  leo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  name
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/leos_Init'

    ! --- local ---------------------------------

    integer                         ::  i

    ! --- begin ---------------------------------
    
    ! store:
    leos%name = trim(name)

    ! any conc output ?
    if ( leo%nconc > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_conc(leo%nconc) )
      ! loop over files to be written:
      do i = 1, leo%nconc
        call LE_Output_Conc_Init( leos%output_conc(i), rcF, &
                     trim(leo%rckey), 'conc', leo%file_conc(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

#ifdef with_labeling
    ! any label output ?
    if ( leo%nlabel > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_label(leo%nlabel) )
      ! loop over files to be written:
      do i = 1, leo%nlabel
        call LE_Output_label_Init( leos%output_label(i), rcF, &
                     trim(leo%rckey), 'label', leo%file_label(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if
    ! any labelled conc output ?
    if ( leo%nlabelled_conc > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_labelled_conc(leo%nlabelled_conc) )
      ! loop over files to be written:
      do i = 1, leo%nlabelled_conc
        call LE_Output_labelled_Conc_Init( leos%output_labelled_conc(i), rcF, &
                     trim(leo%rckey), 'labelled_conc', leo%file_labelled_conc(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if
    ! any labelled emissions output ?
    if ( leo%nemislabel > 0  ) then
      ! storage for output stuff:
      allocate( leos%output_emis_label(leo%nemislabel) )
      ! loop over files to be written:
      do i = 1, leo%nemislabel
        call LE_Output_Emis_label_Init( leos%output_emis_label(i), rcF, &
                     trim(leo%rckey), 'emis_label', leo%file_emis_label(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! any labelled dry deposition output ?
    if ( leo%ndrydeplabel > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_drydep_label(leo%ndrydeplabel) )
      ! loop over files to be written:
      do i = 1, leo%ndrydeplabel
        call LE_Output_Drydep_label_Init( leos%output_drydep_label(i), rcF, &
                     trim(leo%rckey), 'drydep_label', leo%file_drydep_label(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! any labelled wet deposition output ?
    if ( leo%nwetdeplabel > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_wetdep_label(leo%nwetdeplabel) )
      ! loop over files to be written:
      do i = 1, leo%nwetdeplabel
        call LE_Output_Wetdep_label_Init( leos%output_wetdep_label(i), rcF, &
                     trim(leo%rckey), 'wetdep_label', leo%file_wetdep_label(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! any labelled wet deposition output ?
    if ( leo%nbudgetlabel > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_budget_label(leo%nbudgetlabel) )
      ! loop over files to be written:
      do i = 1, leo%nbudgetlabel
        call LE_Output_Budget_label_Init( leos%output_budget_label(i), rcF, &
                     trim(leo%rckey), 'budget_label', leo%file_budget_label(i), &
                     status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if
    
    !
#endif
    
#ifdef with_radiation
    ! any radiation output ?
    if ( leo%nrad > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_rad(leo%nrad) )
      ! loop over files to be written:
      do i = 1, leo%nrad
        call LE_Output_Rad_Init( leos%output_rad(i), rcF, &
                     trim(leo%rckey), 'rad', leo%file_rad(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if
#endif

    ! any emis output ?
    if ( leo%nemis > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_emis(leo%nemis) )
      ! loop over files to be written:
      do i = 1, leo%nemis
        call LE_Output_Emis_Init( leos%output_emis(i), rcF, &
                     trim(leo%rckey), 'emis', leo%file_emis(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! any vd output ?
    if ( leo%nvd > 0 ) then
       ! storage for output stuff:
       allocate( leos%output_vd(leo%nvd) )
       ! loop over files to be written:
       do i = 1, leo%nvd
          call LE_Output_vd_Init( leos%output_vd(i), rcF, &
                       trim(leo%rckey), 'vd', leo%file_vd(i), &
                       trim(leos%name), status )
          IF_NOTOK_RETURN(status=1)
       end do
    end if
    
    ! any diagnostic vd output ?
    if ( leo%nvd_diag > 0 ) then
       ! storage for output stuff:
       allocate( leos%output_vd_diag(leo%nvd_diag) )
       ! loop over files to be written:
       do i = 1, leo%nvd_diag
          call LE_Output_vd_diag_Init( leos%output_vd_diag(i), rcF, &
                       trim(leo%rckey), 'vd-diag', leo%file_vd_diag(i), &
                       trim(leos%name), status )
          IF_NOTOK_RETURN(status=1)
       end do
    end if
    
    ! any drydep output ?
    if ( leo%ndrydep > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_drydep(leo%ndrydep) )
      ! loop over files to be written:
      do i = 1, leo%ndrydep
        call LE_Output_DryDep_Init( leos%output_drydep(i), rcF, &
                     trim(leo%rckey), 'drydep', leo%file_drydep(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! any wetdep output ?
    if ( leo%nwetdep > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_wetdep(leo%nwetdep) )
      ! loop over files to be written:
      do i = 1, leo%nwetdep
        call LE_Output_WetDep_Init( leos%output_wetdep(i), rcF, &
                     trim(leo%rckey), 'wetdep', leo%file_wetdep(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! any budget output ?
    if ( leo%nbudget > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_budget(leo%nbudget) )
      ! loop over files to be written:
      do i = 1, leo%nbudget
        call LE_Output_budget_Init( leos%output_budget(i), rcF, &
                     trim(leo%rckey), 'budget',  leo%file_budget(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! any column output ?
    if ( leo%ncolumn > 0 ) then
      ! storage for output stuff:
      allocate( leos%output_column(leo%ncolumn) )
      ! loop over files to be written:
      do i = 1, leo%ncolumn
        call LE_Output_Column_Init( leos%output_column(i), rcF, &
                     trim(leo%rckey), 'column', leo%file_column(i), &
                     trim(leos%name), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end if

    ! any cso output ?
    if ( leo%ncso > 0 ) then
      ! storage for output stuff:
      allocate( leos%cso_state(leo%ncso), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! loop over files to be written:
      do i = 1, leo%ncso
        ! setup state:
        call leos%cso_state(i)%Init( leo%cso_data(i), status, key=trim(leos%name), &
                                       description='simulated '//trim(leo%file_cso(i))//' retrievals' )
        IF_NOTOK_RETURN(status=1)
      end do  ! output files
    end if  ! any output?

    ! ok
    status = 0

  end subroutine leos_Init


  ! ***


  subroutine leos_Done( leos, leo, status )

    use LE_Output_Conc   , only : LE_Output_Conc_Done
#ifdef with_labeling
    use LE_Output_Label   , only : LE_Output_Label_Done
    use LE_Output_Labelled_Conc   , only : LE_Output_Labelled_conc_Done
    use LE_Output_Emis_Label   , only : LE_Output_Emis_Label_Done
    use LE_Output_Drydep_Label  , only : LE_Output_Drydep_Label_Done
    use LE_Output_Wetdep_Label  , only : LE_Output_Wetdep_Label_Done
    use LE_Output_Budget_Label  , only : LE_Output_Budget_Label_Done
#endif
#ifdef with_radiation
    use LE_Output_Rad    , only : LE_Output_Rad_Done
#endif
    use LE_Output_Emis   , only : LE_Output_Emis_Done
    use LE_Output_vd     , only : LE_Output_vd_Done
    use LE_Output_vd_diag, only : LE_Output_vd_diag_Done
    use LE_Output_Drydep , only : LE_Output_Drydep_Done
    use LE_Output_Wetdep , only : LE_Output_Wetdep_Done
    use LE_Output_budget , only : LE_Output_budget_Done
    use LE_Output_Column , only : LE_Output_Column_Done

    ! --- in/out --------------------------------

    type(T_LE_OutputState), intent(inout)     ::  leos
    type(T_LE_Output), intent(in)             ::  leo
    integer, intent(out)                      ::  status

    ! --- const ----------------------------

    character(len=*), parameter   ::  rname = mname//'/leos_Done'

    ! --- local ------------------------------

    integer        ::  i

    ! --- begin ---------------------------------

    ! any conc output ?
    if ( leo%nconc > 0 ) then
      ! loop over files:
      do i = 1, leo%nconc
        ! done with output file:
        call LE_Output_Conc_Done( leos%output_conc(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_conc )
    end if

#ifdef with_labeling
    ! any label output ?
    if ( leo%nlabel > 0 ) then
      ! loop over files:
      do i = 1, leo%nlabel
        ! done with output file:
        call LE_Output_label_Done( leos%output_label(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_label )
    end if
  
    ! any label output ?
    if ( leo%nlabelled_conc > 0 ) then
      ! loop over files:
      do i = 1, leo%nlabelled_conc
        ! done with output file:
        call LE_Output_labelled_conc_Done( leos%output_labelled_conc(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_labelled_conc )
    end if
  
    ! any label output ?
    if ( leo%nemislabel > 0 ) then
      ! loop over files:
      do i = 1, leo%nemislabel
        ! done with output file:
        call LE_Output_Emis_label_Done( leos%output_emis_label(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_emis_label )
    end if
    
    ! any dry deposition label output ?
    if ( leo%ndrydeplabel > 0 ) then
      ! loop over files:
      do i = 1, leo%ndrydeplabel
        ! done with output file:
        call LE_Output_Drydep_label_Done( leos%output_drydep_label(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_drydep_label )
    end if

    ! any wet deposition label output ?
    if ( leo%nwetdeplabel > 0 ) then
      ! loop over files:
      do i = 1, leo%nwetdeplabel
        ! done with output file:
        call LE_Output_Wetdep_label_Done( leos%output_wetdep_label(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_wetdep_label )
    end if

    ! any budget label output ?
    if ( leo%nbudgetlabel > 0 ) then
      ! loop over files:
      do i = 1, leo%nbudgetlabel
        ! done with output file:
        call LE_Output_Budget_label_Done( leos%output_budget_label(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_budget_label )
    end if

#endif        

#ifdef with_radiation
    ! any radiation output ?
    if ( leo%nrad > 0 ) then
      ! loop over files:
      do i = 1, leo%nrad
        ! done with output file:
        call LE_Output_Rad_Done( leos%output_rad(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_rad )
    end if
#endif

    ! any emis output ?
    if ( leo%nemis > 0 ) then
      ! loop over files:
      do i = 1, leo%nemis
        ! done with output file:
        call LE_Output_Emis_Done( leos%output_emis(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_emis )
    end if

    ! any vd output ?
    if ( leo%nvd > 0 ) then
       ! loop over files:
       do i = 1, leo%nvd
          ! done with output file:
          call LE_Output_vd_Done( leos%output_vd(i), status )
          IF_NOTOK_RETURN(status=1)
       end do
       ! clear:
       deallocate( leos%output_vd )
    end if
    
    ! any diagnostic vd output ?
    if ( leo%nvd_diag > 0 ) then
       ! loop over files:
       do i = 1, leo%nvd_diag
          ! done with output file:
          call LE_Output_vd_diag_Done( leos%output_vd_diag(i), status )
          IF_NOTOK_RETURN(status=1)
       end do
       ! clear:
       deallocate( leos%output_vd_diag )
    end if
    
    ! any drydep output ?
    if ( leo%ndrydep > 0 ) then
      ! loop over files:
      do i = 1, leo%ndrydep
        ! done with output file:
        call LE_Output_DryDep_Done( leos%output_drydep(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_drydep )
    end if

    ! any wetdep output ?
    if ( leo%nwetdep > 0 ) then
      ! loop over files:
      do i = 1, leo%nwetdep
        ! done with output file:
        call LE_Output_WetDep_Done( leos%output_wetdep(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_wetdep )
    end if

    ! any budget output ?
    if ( leo%nbudget > 0 ) then
      ! loop over files:
      do i = 1, leo%nbudget
        ! done with output file:
        call LE_Output_budget_Done( leos%output_budget(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_budget )
    end if

    ! any column output ?
    if ( leo%ncolumn > 0 ) then
      ! loop over files:
      do i = 1, leo%ncolumn
        ! done with output file:
        call LE_Output_Column_Done( leos%output_column(i), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! clear:
      deallocate( leos%output_column )
    end if

    ! any cso output ?
    if ( leo%ncso > 0 ) then
      ! loop over files to be written:
      do i = 1, leo%ncso
        ! done with state:
        call leos%cso_state(i)%Done( leo%cso_data(i), status )
        IF_NOTOK_RETURN(status=1)
      end do  ! output files
      ! clear:
      deallocate( leos%cso_state, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if  ! any output?

    ! ok
    status = 0

  end subroutine leos_Done


  ! ***


  subroutine leos_PutOut( leos, leo, t, c, cg, bud, status, last, without_data )

    use GO               , only : TDate, MidNight, operator(>)
    use GO               , only : goMatchValue
    use dims             , only : runF
    use dims             , only : nx, ny, nz, nspec
    use indices          , only : i_o3, ispec_o3
    use LE_Output_Conc   , only : LE_Output_Conc_PutOut
#ifdef with_labeling
    use LE_Output_Label   , only : LE_Output_Label_PutOut
    use LE_Output_Labelled_Conc   , only : LE_Output_Labelled_Conc_PutOut
    use LE_Output_Emis_Label   , only : LE_Output_Emis_Label_PutOut
    use LE_Output_Drydep_Label  , only : LE_Output_Drydep_Label_PutOut
    use LE_Output_Wetdep_Label  , only : LE_Output_Wetdep_Label_PutOut
    use LE_Output_Budget_Label  , only : LE_Output_Budget_Label_PutOut
    ! reset labelled daily dry deposition budgets
    use SA_Labeling             , only : SA_Depos_Reset_day
#endif        
#ifdef with_radiation
    use LE_Output_Rad    , only : LE_Output_Rad_PutOut
#endif

    use LE_Output_Emis   , only : LE_Output_Emis_PutOut
    use LE_Output_vd     , only : LE_Output_vd_PutOut
    use LE_Output_vd_diag, only : LE_Output_vd_diag_PutOut
    use LE_Output_Drydep , only : LE_Output_Drydep_PutOut
    use LE_Output_Wetdep , only : LE_Output_Wetdep_PutOut
    use LE_Output_budget , only : LE_Output_budget_PutOut
    use LE_Output_Column , only : LE_Output_Column_PutOut
    use LE_Budget         , only : T_Budget, Budget_Update_O3max, Budget_Reset_Day
    use LE_Budget_WetDepos, only : ex_wetdepo

    ! --- in/out --------------------------------

    type(T_LE_OutputState), intent(inout)     ::  leos
    type(T_LE_Output), intent(in)             ::  leo
    type(TDate), intent(in)                   ::  t
    real, intent(in)                          ::  c(nx,ny,nz,nspec)
    real, intent(in)                          ::  cg(nx,ny,nspec)
    type(T_Budget), intent(inout)             ::  bud
    integer, intent(out)                      ::  status
    
    logical, intent(in), optional             ::  last
    logical, intent(in), optional             ::  without_data

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/leos_PutOut'

    ! --- local ---------------------------------

    integer            ::  i
    logical            ::  with_data

    ! --- begin ---------------------------------
    
    ! flag to skip output of model data, used in LEKF:
    with_data = .true.
    if ( present(without_data) ) with_data = .not. without_data

    ! ozone enabled ?
    if ( ispec_o3 > 0 ) then
      ! update ozone maximum:
      call Budget_Update_O3max( bud, cg(:,:,i_o3), status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! put out conc files ?
    if ( leo%nconc > 0 ) then

      ! loop over concentration files to be written:
      do i = 1, leo%nconc
        ! write instantenous fields:
        call LE_Output_Conc_PutOut( leos%output_conc(i), t, c, cg, status )
        IF_NOTOK_RETURN(status=1)
      end do

    end if
    
#ifdef with_labeling
    ! put out label files ?
    if ( leo%nlabel > 0 ) then

      ! loop over labelentration files to be written:
      do i = 1, leo%nlabel
        ! write instantenous fields:
        call LE_Output_label_PutOut( leos%output_label(i), t, c, cg, status )
        IF_NOTOK_RETURN(status=1)
      end do

    end if
    
    ! put out labelled conc files ?
    if ( leo%nlabelled_conc > 0 ) then

      ! loop over labelentration files to be written:
      do i = 1, leo%nlabelled_conc
        ! write instantenous fields:
        call LE_Output_labelled_conc_PutOut( leos%output_labelled_conc(i), t, c, cg, status )
        IF_NOTOK_RETURN(status=1)
      end do

    end if
    
    ! put out emission label files ?
    if ( leo%nemislabel > 0 ) then

      ! emis files only after first timestep:
      if ( t > runF%t_start ) then
        ! loop over labelentration files to be written:
        do i = 1, leo%nemislabel
          ! write instantenous fields:
          call LE_Output_Emis_label_PutOut( leos%output_emis_label(i), t, c, cg, status )
          IF_NOTOK_RETURN(status=1)
        end do
      end if
    end if
    
    ! put out dry deposition label files ?
    if ( leo%ndrydeplabel > 0 ) then

      ! budget files only after first new concentrations are computed:
      if ( t > runF%t_start ) then
        ! loop over labelentration files to be written:
        do i = 1, leo%ndrydeplabel
          ! write instantenous fields:
          call LE_Output_Drydep_label_PutOut( leos%output_drydep_label(i), t, status )
          IF_NOTOK_RETURN(status=1)
        end do
      end if
    end if
    
    ! put out wet deposition label files ?
    if ( leo%nwetdeplabel > 0 ) then

      ! budget files only after first new concentrations are computed:
      if ( t > runF%t_start ) then
        ! loop over labelentration files to be written:
        do i = 1, leo%nwetdeplabel
          ! write instantenous fields:
          call LE_Output_Wetdep_label_PutOut( leos%output_wetdep_label(i), t, status )
          IF_NOTOK_RETURN(status=1)
        end do
      end if
    end if

    ! put out budget label files ?
    if ( leo%nbudgetlabel > 0 ) then

      ! budget files are only written at midnight,
      ! but skip the start time (might be midnight too):
      if ( (t > runF%t_start) .and. MidNight(t) ) then

        ! loop over labelentration files to be written:
        do i = 1, leo%nbudgetlabel
          ! write instantenous fields:
          call LE_Output_Budget_label_PutOut( leos%output_budget_label(i), t, status )
          IF_NOTOK_RETURN(status=1)
        end do

        ! reset daily budgets:
        call SA_Depos_Reset_day( status )
        IF_NOTOK_RETURN(status=1)

      end if  ! midnight

    end if
#endif        
    
#ifdef with_radiation
    ! put out radiation files ?
    if ( leo%nrad > 0 ) then

      ! loop over radiation files to be written:
      do i = 1, leo%nrad
        ! write instantenous fields:
        call LE_Output_Rad_PutOut( leos%output_rad(i), t, c, status )
        IF_NOTOK_RETURN(status=1)
      end do  ! rad files
    end if ! something to put out
#endif  

    ! put out emis files ?
    if ( with_data .and. (leo%nemis > 0) ) then

      ! emis files only after first timestep:
      if ( t > runF%t_start ) then
        ! loop over emisentration files to be written:
        do i = 1, leo%nemis
          ! write instantenous fields:
          call LE_Output_Emis_PutOut( leos%output_emis(i), t, status )
          IF_NOTOK_RETURN(status=1)
        end do  ! emis files
      end if  ! after init
    end if ! something to put out
    
    ! put out vd files ?
    if ( leo%nvd > 0 ) then
       do i = 1, leo%nvd
          call LE_Output_vd_PutOut( leos%output_vd(i), t, c, bud%drydepos, status )
          IF_NOTOK_RETURN(status=1)
       end do
    end if
    
    ! put out diagnostic vd files ?
    if ( leo%nvd_diag > 0 ) then
       do i = 1, leo%nvd_diag
          call LE_Output_vd_diag_PutOut( leos%output_vd_diag(i), t, c, bud%drydepos, status )
          IF_NOTOK_RETURN(status=1)
       end do
    end if
    
    ! put out drydep files ?
    if ( leo%ndrydep > 0 ) then
      ! budget files only after first new concentrations are computed:
      if ( t > runF%t_start ) then
        ! loop over drydep files to be written:
        do i = 1, leo%ndrydep
          ! write current budgets:
          call LE_Output_DryDep_PutOut( leos%output_drydep(i), t, &
                           bud%drydepos%ex_hour, status )
          IF_NOTOK_RETURN(status=1)
        end do
      end if  ! after start ?
    end if
    
    ! put out wetdep files ?
    if ( leo%nwetdep > 0 ) then
      ! budget files only after first new concentrations are computed:
      if ( t > runF%t_start ) then
        ! loop over wetdep files to be written:
        do i = 1, leo%nwetdep
          ! write current budgets:
          call LE_Output_WetDep_PutOut( leos%output_wetdep(i), t, &
                            bud%wetdepos%ex_hour(:,:,:,:,ex_wetdepo), status )
          IF_NOTOK_RETURN(status=1)
        end do
      end if  ! after start ?
    end if

    ! put out budget files ?
    if ( leo%nbudget > 0 ) then
      ! budget files are only written at midnight,
      ! but skip the start time (might be midnight too):
      if ( (t > runF%t_start) .and. MidNight(t) ) then

        ! loop over budgetentration files to be written:
        do i = 1, leo%nbudget
          ! write instantenous fields:
          call LE_Output_budget_PutOut( leos%output_budget(i), t, bud, status )
          IF_NOTOK_RETURN(status=1)
        end do

        ! reset daily budgets:
        call Budget_Reset_Day( bud, status )
        IF_NOTOK_RETURN(status=1)

      end if   ! midnight
    end if
    
    ! put out column files ?
    if ( leo%ncolumn > 0 ) then

      ! loop over columnentration files to be written:
      do i = 1, leo%ncolumn
        ! write instantenous fields:
        call LE_Output_Column_PutOut( leos%output_column(i), t, c, status )
        IF_NOTOK_RETURN(status=1)
      end do

    end if

    ! any cso output ?
    if ( leo%ncso > 0 ) then
      ! loop over files to be written:
      do i = 1, leo%ncso
        ! simulate retrievals:
        call leos%cso_state(i)%Setup( leo%cso_data(i), c, status )
        IF_NOTOK_RETURN(status=1)
        ! put out:
        call leos%cso_state(i)%PutOut( leo%cso_data(i), t, status )
        IF_NOTOK_RETURN(status=1)
      end do  ! output files
    end if  ! any output?

    ! ok
    status = 0

  end subroutine leos_PutOut



end module LE_Output
