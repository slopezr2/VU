!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!######################################################################

module C3PO_Spectral_Gradient_Main

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  ParseArguments
  public  ::  Spectral_Gradient
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Spectral_Gradient_Main'
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** arguments
  ! ***
  ! ********************************************************************


  !
  ! Return status:
  !     -1 : usage displayed
  !      0 : ok
  !  other : error
  !
  
  subroutine ParseArguments( input_filename, output_filename, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(out)    ::  input_filename
    character(len=*), intent(out)    ::  output_filename
    integer, intent(out)             ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ParseArguments'
    
    ! --- local ----------------------------------
    
    integer                         ::  narg, iarg
    integer                         ::  arglength
    character(len=1024)             ::  argvalue
    
    ! --- begin ----------------------------------
    
    ! dummy values at start:
    input_filename = ''
    output_filename = ''
    
    ! count:
    narg = Command_Argument_Count()

    ! loop:
    do iarg = 1, narg
      ! get length:
      call Get_Command_Argument( iarg, length=arglength, status=status )
      if ( status /= 0 ) then
        write (gol,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength; call goErr
        TRACEBACK; status=1; return
      end if
      ! get value: 
      call Get_Command_Argument( iarg, value=argvalue, status=status )
      if ( status /= 0 ) then
        write (gol,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength; call goErr
        TRACEBACK; status=1; return
      end if
      ! switch:
      select case ( trim(argvalue) )
        case ( '-h', '--help' )
          write (gol,'("Usage:")'); call goPr
          write (gol,'("  c3po_spectral_gradient.x inputfile outputfile")'); call goPr
          write (gol,'("  c3po_select.x -h|--help")'); call goPr
          status = -1; return
        case default
          if ( len_trim(input_filename) == 0 ) then
            input_filename = trim(argvalue)
          else if ( len_trim(output_filename) == 0 ) then
            output_filename = trim(argvalue)
          else
            write (*,'("unsupported argument `",a,"`")') trim(argvalue)
            TRACEBACK; status=1; return
          end if
      end select
    end do   ! arguments

    ! check ...
    if ( len_trim(input_filename) == 0 ) then
      write (*,'("no input filename argument specified")')
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( len_trim(output_filename) == 0 ) then
      write (*,'("no output filename argument specified")')
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
  
  end subroutine ParseArguments
  
  
  ! ***
  
  subroutine Spectral_Gradient( input_filename, output_filename, status )
  
    use C3PO_File_Spg, only : T_File_Spg
    use C3PO_Grid_Spg, only : T_Grid_Spg
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)    ::  input_filename
    character(len=*), intent(in)    ::  output_filename
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Spectral_Gradient'
    
    ! --- local ----------------------------------
    
    logical                 ::  exist
    type(T_File_Spg)        ::  spgfile_in
    type(T_File_Spg)        ::  spgfile_out
    type(T_Grid_Spg)        ::  spg
    character(len=1024)     ::  message
    integer                 ::  ndim, idim
    character(len=32)       ::  dimname
    integer                 ::  nvar, ivar
    character(len=32)       ::  varname
    character(len=32)       ::  dvarnames(2)
    logical                 ::  spectral
    integer, allocatable    ::  varids(:)
    integer, allocatable    ::  varids_nabla(:,:)
    
    ! --- begin ----------------------------------
    
    ! check ..
    inquire( file=trim(input_filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("input file not found : ",a)') trim(input_filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! info ...
    write (gol,'("open ",a," ...")') trim(input_filename); call goPr
    
    ! open file for reading:
    call spgfile_in%Open( trim(input_filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("read grid definition ...")'); call goPr
    
    ! fill grid definition from file variables:
    call spgfile_in%Get_Grid( spg, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("create target file ",a," ...")') trim(output_filename); call goPr
    
    ! create target file:
    call spgfile_out%Create( trim(output_filename), status, formatNum=spgfile_in%formatNum )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("copy attributes ...")'); call goPr
    
    ! copy global attributes:
    call spgfile_in%Copy_Atts( '', spgfile_out, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("extend history ...")'); call goPr
    
    ! write informative message:
    write (message,'("c3po: computed gradient of spectral variables from ",a)') &
                               trim(input_filename)
    ! add to history attribute:
    call spgfile_out%Extend_History( trim(message), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("write grid definition ...")'); call goPr
    
    ! define grid definition from file variables:
    call spgfile_out%Def_Grid( spg, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of dimensions:
    call spgfile_in%Inquire( status, ndim=ndim )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do idim = 1, ndim
      ! get name:
      call spgfile_in%Inquire_Dimension( idim, status, name=dimname )
      IF_NOT_OK_RETURN(status=1)
      ! skip grid dimensions:
      if ( trim(dimname) == 'nc2' ) cycle
      if ( trim(dimname) == 'nsp' ) cycle
      ! info ...
      write (gol,'("copy dimension ",a," ...")') trim(dimname); call goPr
      ! other:
      call spgfile_in%Copy_Dimension( trim(dimname), spgfile_out, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! number of variables:
    call spgfile_in%Inquire( status, nvar=nvar )
    IF_NOT_OK_RETURN(status=1)
    ! storage for variable id's in output file:
    allocate( varids(nvar) )
    allocate( varids_nabla(nvar,2) )
    ! init with dummy values, used later on to check on need to copy:
    varids       = -999
    varids_nabla = -999
    ! loop:
    do ivar = 1, nvar
      ! get name:
      call spgfile_in%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! skip grid variables:
      !if ( trim(varname) == 'rgrid'     ) cycle
      ! get flag:
      call spgfile_in%Inquire_Variable_Spectral( ivar, status, spectral=spectral )
      IF_NOT_OK_RETURN(status=1)
      ! spectral field ?
      if ( spectral ) then
        ! info ...
        write (gol,'("define gradients for variable ",a," ...")') trim(varname); call goPr
        !
        ! define as copy:
        call spgfile_in%Define_Variable_Copy( varname, spgfile_out, varids_nabla(ivar,1), status )
        IF_NOT_OK_RETURN(status=1)
        ! rename:
        write (dvarnames(1),'("d",a,"_dlambda")') trim(varname)
        call spgfile_out%Rename_Variable( varname, trim(dvarnames(1)), status )
        IF_NOT_OK_RETURN(status=1)
        ! extend attribute:
        call spgfile_out%Extend_Variable_Att( dvarnames(1), 'derivative of ', 'long_name', ' to lambda (longitude in radians)', status )
        IF_NOT_OK_RETURN(status=1)
        !
        ! define as copy:
        call spgfile_in%Define_Variable_Copy( varname, spgfile_out, varids_nabla(ivar,2), status )
        IF_NOT_OK_RETURN(status=1)
        ! rename:
        write (dvarnames(2),'("d",a,"_dmu")') trim(varname)
        call spgfile_out%Rename_Variable( varname, trim(dvarnames(2)), status )
        IF_NOT_OK_RETURN(status=1)
        ! extend attribute:
        call spgfile_out%Extend_Variable_Att( dvarnames(2), 'derivative of ', 'long_name', ' to mu (sin(latitude) in radians)', status )
        IF_NOT_OK_RETURN(status=1)
        !
      else
        ! info ...
        write (gol,'("define variable ",a," ...")') trim(varname); call goPr
        ! define copy:
        call spgfile_in%Define_Variable_Copy( varname, spgfile_out, varids(ivar), status )
        IF_NOT_OK_RETURN(status=1)
      end if
    end do
    
    ! end of definition phase:
    call spgfile_out%EndDef( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write grid definition:
    call spgfile_out%Put_Grid( spg, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy data of non-spectral variables:
    do ivar = 1, nvar
      ! not defined ? then skip:
      if ( varids(ivar) < 0 ) cycle
      ! get name:
      call spgfile_in%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'("copy data of ",a," ...")') trim(varname); call goPr
      ! copy selection:
      call spgfile_in%Copy_Variable( varname, spgfile_out, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! write gradients of spectral fields:
    do ivar = 1, nvar
      ! not defined ? then skip:
      if ( any(varids_nabla(ivar,:) < 0) ) cycle
      ! get name:
      call spgfile_in%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'("write gradients of ",a," ...")') trim(varname); call goPr
      ! copy selection:
      call spgfile_in%Put_Variable_Gradients( spg, varname, spgfile_out, dvarnames, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! info ...
    write (gol,'("done ...")'); call goPr
    
    ! clear:
    deallocate( varids )
    deallocate( varids_nabla )
    
    ! close file:
    call spgfile_out%Close( status )
    IF_NOT_OK_RETURN(status=1) 

    ! clear:
    call spg%Final()
    
    ! close file:
    call spgfile_in%Close( status )
    IF_NOT_OK_RETURN(status=1)    
    
    ! ok
    status = 0
  
  end subroutine Spectral_Gradient


end module C3PO_Spectral_Gradient_Main
