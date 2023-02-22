!
! go print : tools for standard output
!
! Example:
!
!   ! messages printed by root only:
!   call GO_Print_Init( status, apply=myid==root, &
!                           prompt_pe=npes>1, pe=myid, &
!                           trace=.false. )
!   if (status/=0) stop
!
!   ! set routine label:
!   call goLabel( 'mymod/myroutine' )
!
!   ! write single message (including processor prompt?) :
!   !   [00] This is number  3
!   write (gol,'("This is number ",i2)') 3; call goPr
!
!   ! write error message and traceback using the
!   ! previous defined routine label:
!   !   [00] ERROR - Something wrong.
!   !   [00] ERROR in mymod/myroutine
!   write (gol,'("Something wrong.")'); call goErr
!   call goErr
!
!   ! close label
!   call goLabel()
!
!   ! done
!   call GO_Print_Done( status )
!   if (status/=0) stop
!  

! Nedit macro's:
!
!  o change error traceback:
!       write (*,'("ERROR in ",a)') rname
!       call goErr
!
!  o change error traceback:
!       write (*,'("ERROR in ",a)') rname
!       write (gol,'("in ",a)') rname; call goErr
!
!  o change error message:
!       write \(\*,'\("ERROR - (.*$)
!       write (gol,'("\1; call goErr
!
!  o change other message:
!       write \(\*,(.*$)
!       write (gol,\1; call goPr
!
!  o change error time messages:
!       (goprdt.*ERROR.*$)
!       \1; call goErr
!
!  o change time messages:
!       printdate2
!       wrtgol
!       call printdate2\( 'ERROR - (.*$)
!       call wrtgol( '\1; call goErr
!       printdate
!       wrtgol
!       call printdate\( 'ERROR - (.*$)
!       call wrtgol( '\1; call goErr
!       
!       
!
!  o change error messages:
!       (ERROR.*; call )goPr
!       \1goErr
!
!  o change time messages:
!       (call goprdt.*$)
!       \1; call goPr
!


module GO_Print

  implicit none
  
  ! --- in/out ---------------------------------
  
  private
  
  public  ::  gol
  public  ::  GO_Print_Init, GO_Print_Done
  public  ::  goPr, goErr, goBug
  public  ::  goLabel
  
  
  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'GO_Print'

  
  ! --- var ------------------------------------
  
  ! buffer for standard output
  character(len=1024)  ::  gol
  
  ! stack with labels:
  integer, parameter   ::  mstack = 400
  character(len=64)    ::  labels(0:mstack)
  integer              ::  istack = 0
  
  ! initialized ? 
  ! some errors might be printed before initialization ...
  logical              ::  pr_initialized = .false.
  
  ! destination file unit:
  integer              ::  pr_fu
  
  ! flags etc
  logical              ::  pr_apply
  logical              ::  pr_trace
  
  ! processor prompt
  logical              ::  pr_prompt_pe
  integer              ::  pr_pe
  
  ! white space for indents:
  integer, parameter          ::  dindent = 2
  integer                     ::  indent = 0
  
  ! writ to file ?
  logical              ::  pr_file
  character(len=256)   ::  pr_file_name
  
  
contains


  ! ***************************************************************************
  ! ***
  ! *** module init/done
  ! ***
  ! ***************************************************************************
  
  
  subroutine GO_Print_Init( status, apply, prompt_pe, pe, trace, file, file_name )
  
    use go_fu, only : goStdOut
  
    ! --- in/out ----------------------------
    
    integer, intent(out)                     ::  status
    logical, intent(in), optional            ::  apply
    logical, intent(in), optional            ::  prompt_pe
    integer, intent(in), optional            ::  pe
    logical, intent(in), optional            ::  trace
    logical, intent(in), optional            ::  file
    character(len=*), intent(in), optional   ::  file_name
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GO_Print_Init'
    
    ! --- local -----------------------------
    
    logical   ::  opened
    
    ! --- begin -----------------------------
    
    ! print or not ?
    pr_apply = .true.
    if ( present(apply) ) pr_apply = apply
    
    ! processor number
    pr_pe = 0
    if ( present(pe) ) pr_pe = pe
    
    ! prompt processor number ?
    pr_prompt_pe = .false.
    if ( present(prompt_pe) ) pr_prompt_pe = prompt_pe
    
    ! trace labels ?
    pr_trace = .false.
    if ( present(trace) ) pr_trace = trace
    
    ! write to file ?
    pr_file = .false.
    if ( present(file) ) pr_file = file
    pr_file_name = 'go.out'
    if ( present(file_name) ) pr_file_name = file_name
    
    ! init label stack:
    labels(0) = '<no-label>'
    istack = 0
    
    ! no indent yet
    indent = 0
    if ( .not. pr_trace ) indent = -2
    
    ! write messages to file ?
    if ( pr_file ) then
      ! select free file unit:
      pr_fu = 10
      do
        inquire( pr_fu, opened=opened )
        if ( .not. opened ) exit
        pr_fu = pr_fu + 1
      end do
      ! open requested output file:
      open( unit=pr_fu, file=pr_file_name, status='replace', iostat=status )
      if ( status/=0 ) then
        write (*,'("ERROR - opening file for output:")')
        write (*,'("ERROR -   unit : ",i6)') pr_fu
        write (*,'("ERROR -   file : ",a)') trim(pr_file_name)
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
    else
      ! write to standard output:
      pr_fu = goStdOut
    end if

    ! now the module is initialized ...
    pr_initialized = .true.

    ! ok
    status = 0
    
  end subroutine GO_Print_Init
  
  
  ! ***
  
  
  subroutine GO_Print_Done( status )
  
    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GO_Print_Done'
    
    ! --- begin -----------------------------
    
    ! output to file ?
    if ( pr_file ) then
      ! close file:
      close( pr_fu, iostat=status )
      if ( status/=0 ) then
        write (*,'("ERROR - closing output file:")')
        write (*,'("ERROR -   unit : ",i6)') pr_fu
        write (*,'("ERROR -   file : ",a)') trim(pr_file_name)
        write (*,'("ERROR in ",a)') rname; status=1; return
      end if
    end if

    
    ! ok
    status = 0
    
  end subroutine GO_Print_Done
  
  
  ! ***************************************************************************
  ! ***
  ! *** printing
  ! ***
  ! ***************************************************************************
  
  
  subroutine goPr

    ! --- local --------------------------------
    
    character(len=16)   ::  prompt
    integer             ::  nind
  
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goPr'
    
    ! --- begin --------------------------------
    
    ! not initialized yet ? then print to standard output:
    if ( .not. pr_initialized ) then
      write (*,'(a)') trim(gol)
      gol = ''
      return
    end if
    
    ! print go line ?
    if ( pr_apply ) then
    
      ! number of spaces to indent:
      nind = max( 0, indent )

      ! processor prompt ?
      if ( pr_prompt_pe ) then
        write (prompt,'("[",i2.2,"]")') pr_pe
        nind = nind + 1
      else
        prompt = ''
      end if

      ! write prompt, indention, go line:
      if ( nind > 0 ) then  
        write (pr_fu,'(a,a,a)') trim(prompt), repeat(' ',nind), trim(gol)
      else
        write (pr_fu,'(a,a)') trim(prompt), trim(gol)
      end if

    end if
    
    ! clear output line:
    gol = ''

  end subroutine goPr


  ! ***
  
  ! Print error message.
  ! Now printed to standard output, in future to standard error ?
  ! Make gol empty before leaving.
  ! If still empty in next call, this is a trace back
  !   (print error label, one label back)
  

  subroutine goErr
  
    ! --- local -------------------------------
    
    integer     ::  ilab
  
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goErr'
    
    ! --- local ----------------------------
    
    logical                   ::  save_pr_apply
    character(len=len(gol))   ::  gol2
    
    ! --- begin --------------------------------
    
    ! store original apply flag:
    save_pr_apply = pr_apply
    ! always print error messages:
    pr_apply = .true.
    
    ! message in buffer ?
    if ( len_trim(gol) > 0 ) then
    
      ! error message;
      ! make a copy of the message to avoid problems with re-writing on some machiens:
      gol2 = trim(gol)
      write (gol,'("ERROR - ",a)') trim(gol2); call goPr
      
    else
    
      ! label index:
      ilab = min( istack, mstack )

      ! write error message:
      write (gol,'("ERROR in ",a)') trim(labels(ilab)); call goPr

      ! one level back:
      call goLabel()

    end if 
    
    ! restore apply flag:
    pr_apply = save_pr_apply 

  end subroutine goErr


  ! ***
  
  
  subroutine goBug

    ! --- local ----------------------------

    logical   ::  save_pr_apply

    ! --- begin --------------------------------

    ! store original apply flag:
    save_pr_apply = pr_apply
    ! always print bug messages:
    pr_apply = .true.

    ! write message
    write (gol,'("BUG - ",a)') trim(gol); call goPr

    ! restore apply flag:
    pr_apply = save_pr_apply

  end subroutine goBug

  
  ! ***************************************************************************
  ! ***
  ! *** routine labels
  ! ***
  ! ***************************************************************************
  
  
  subroutine goLabel( label )
  
    ! --- in/out -------------------------------
    
    character(len=*), intent(in), optional  ::  label

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goLabel'
    
    ! --- begin --------------------------------
    
    ! add new label to stack ?
    if ( present(label) ) then
      istack = istack + 1
      if ( istack > mstack ) then
        write (gol,'("BUG - stack too small; please increase mstack in go_print")'); call goPr
      else
        labels(istack) = label
      end if
      if (pr_trace) then
        write (gol,'("<",a,">")') trim(labels(istack)); call goPr
      end if
      indent = indent + dindent
    else
      indent = indent - dindent
      if (pr_trace) then
        write (gol,'("(",a,")")') trim(labels(istack)); call goPr
      end if
      istack = max( 0, istack - 1 )
    end if
    
  end subroutine goLabel
    

end module go_print



! #############################################################################
! ###
! ### test program
! ###
! #############################################################################
!
!
!module testmod
!
!  implicit none
!  
!  public
!  
!contains
!
!  subroutine subr( i, status )
!  
!    use go_print, only : goLabel, gol, goPr, goErr
!    
!    ! --- in/out ----------------------------------------
!    
!    integer, intent(in)           ::  i
!    integer, intent(out)          ::  status
!    
!    ! --- begin -----------------------------------------
!    
!    call goLabel( 'subr' )
!    
!    write (gol,'("welcome to subr !")'); call goPr
!    
!    select case ( i )
!    
!      case ( 0 )
!        write (gol,'("testing i : ",i2)') i; call goPr
!      
!      case ( 1 )
!        call subr2( 0, status )
!        if (status/=0) then; call goErr; status=1; return; end if
!      
!      case ( 2 )
!        call subr2( 1, status )
!        if (status/=0) then; call goErr; status=1; return; end if
!      
!      case default
!        write (gol,'("unsupported i : ",i2)') i; call goErr
!        call goErr; status=1; return
!        
!    end select
!        
!    call goLabel(); status=0
!    
!  end subroutine subr
!    
!
!  ! ***
!
!    
!  subroutine subr2( i, status )
!  
!    use go_print, only : goLabel, gol, goPr, goErr
!    
!    ! --- in/out ----------------------------------------
!    
!    integer, intent(in)           ::  i
!    integer, intent(out)          ::  status
!    
!    ! --- begin -----------------------------------------
!    
!    call goLabel('subr2')
!
!    write (gol,'("testing subr2")'); call goPr
!    
!    select case ( i )
!      case ( 0 )
!      case default
!        write (gol,'("wrong i : ",i2)') i; call goErr
!        call goErr; status=1; return
!    end select
!    
!    call goLabel; status=0
!    
!  end subroutine subr2
!    
!
!
!end module testmod
!
!
! ################################################################
!
!
!program test
!
!  use go_print
!  use testmod
!  
!  ! --- local -----------------------------------------
!  
!  integer           ::  status
!  
!  ! --- begin ------------------------------------------
!
!  call GO_Print_Init( status, trace=.false. )
!  call goLabel('test prog')
!  
!  write (gol,'("begin of program")'); call goPr
!  
!  call Subr( 2, status )
!  if (status/=0) then; call goErr; call exit(1); end if
!  
!  write (gol,'("end of program")'); call goPr
!
!  call goLabel()
!  call GO_Print_Done( status )
!
!end program test
!
