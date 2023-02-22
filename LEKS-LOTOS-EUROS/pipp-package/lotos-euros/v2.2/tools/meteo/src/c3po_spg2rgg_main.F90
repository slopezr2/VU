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

module C3PO_SPG2RGG_Main

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  ParseArguments
  public  ::  SPG2RGG
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_SPG2RGG_Main'
  
  

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
  
  subroutine ParseArguments( input_filename, template_filename, output_filename, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(out)    ::  input_filename
    character(len=*), intent(out)    ::  template_filename
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
    template_filename = ''
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
          write (gol,'("  c3po_SPG2RGG.x inputfile templatefile outputfile")'); call goPr
          write (gol,'("  c3po_select.x -h|--help")'); call goPr
          status = -1; return
        case default
          if ( len_trim(input_filename) == 0 ) then
            input_filename = trim(argvalue)
          else if ( len_trim(template_filename) == 0 ) then
            template_filename = trim(argvalue)
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
    if ( len_trim(template_filename) == 0 ) then
      write (*,'("no template filename argument specified")')
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
  
  subroutine SPG2RGG( input_filename, template_filename, output_filename, status )
  
    use C3PO_Grid_Spg    , only : T_Grid_Spg
    use C3PO_Grid_Rgg    , only : T_Grid_Rgg
    use C3PO_File_Spg    , only : T_File_Spg
    use C3PO_File_Rgg    , only : T_File_Rgg
    use C3PO_File_Convert, only : Spg_Define_Variable_Copy_Rgg, Spg_Define_Variable_RegriddedCopy_Rgg
    use C3PO_File_Convert, only : Spg_Regrid_Variable_Rgg, Spg_Copy_Variable_Rgg
    use C3PO_File_Convert, only : Spg_Copy_Atts_Rgg
    use C3PO_File_Convert, only : Spg_Copy_Dimension_Rgg
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)    ::  input_filename
    character(len=*), intent(in)    ::  template_filename
    character(len=*), intent(in)    ::  output_filename
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SPG2RGG'
    
    ! --- local ----------------------------------
    
    logical                 ::  exist
    type(T_File_Spg)        ::  spgfile
    type(T_File_Rgg)        ::  rggfile_in
    type(T_File_Rgg)        ::  rggfile
    type(T_Grid_Spg)        ::  spg
    type(T_Grid_Rgg)        ::  rgg
    character(len=1024)     ::  message
    integer                 ::  ndim, idim
    character(len=32)       ::  dimname
    integer                 ::  nvar, ivar
    character(len=32)       ::  varname
    logical                 ::  spectral
    integer, allocatable    ::  varids(:)
    integer, allocatable    ::  varids_spectral(:)
    
    ! --- begin ----------------------------------
    
    ! check ..
    inquire( file=trim(input_filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("input file not found : ",a)') trim(input_filename); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    inquire( file=trim(template_filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("template file not found : ",a)') trim(template_filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! info ...
    write (gol,'("read template file ",a," ...")') trim(template_filename); call goPr
    
    ! open file for reading:
    call rggfile_in%Open( trim(template_filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! fill grid definition from file variables:
    call rggfile_in%Get_Grid( rgg, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! done:
    call rggfile_in%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("open input file ",a," ...")') trim(input_filename); call goPr
    
    ! open file for reading:
    call spgfile%Open( trim(input_filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("read grid definition ...")'); call goPr
    
    ! fill grid definition from file variables:
    call spgfile%Get_Grid( spg, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("create target file ",a," ...")') trim(output_filename); call goPr
    
    ! create target file:
    call rggfile%Create( trim(output_filename), status, formatNum=spgfile%formatNum )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("copy attributes ...")'); call goPr
    
    ! copy global attributes:
    !call spgfile%Copy_Atts( '', rggfile, status )
    call Spg_Copy_Atts_Rgg( spgfile, '', rggfile, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("extend history ...")'); call goPr
    
    ! write informative message:
    write (message,'("c3po: evaluated on reduced gaussian grid from spectral fields in ",a)') &
                               trim(input_filename)
    ! add to history attribute:
    call rggfile%Extend_History( trim(message), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("write grid definition ...")'); call goPr
    
    ! number of dimensions:
    call spgfile%Inquire( status, ndim=ndim )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do idim = 1, ndim
      ! get name:
      call spgfile%Inquire_Dimension( idim, status, name=dimname )
      IF_NOT_OK_RETURN(status=1)
      ! skip grid dimensions:
      if ( trim(dimname) == 'nc2' ) cycle
      if ( trim(dimname) == 'nsp' ) cycle
      ! info ...
      write (gol,'("copy dimension ",a," ...")') trim(dimname); call goPr
      ! other:
      !call spgfile%Copy_Dimension( trim(dimname), rggfile, status )
      call Spg_Copy_Dimension_Rgg( spgfile, trim(dimname), rggfile, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! define grid definition from file variables:
    call rggfile%Def_Grid( rgg, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of variables:
    call spgfile%Inquire( status, nvar=nvar )
    IF_NOT_OK_RETURN(status=1)
    ! storage for variable id's in output file:
    allocate( varids(nvar) )
    allocate( varids_spectral(nvar) )
    ! init with dummy values, used later on to check on need to copy:
    varids          = -999
    varids_spectral = -999
    ! loop:
    do ivar = 1, nvar
      ! get name:
      call spgfile%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! skip grid variables:
      if ( trim(varname) == 'rgrid'     ) cycle
      ! get flag:
      call spgfile%Inquire_Variable_Spectral( ivar, status, spectral=spectral )
      IF_NOT_OK_RETURN(status=1)
      ! spectral field ?
      if ( spectral ) then
        ! info ...
        write (gol,'("define variable ",a," from sepectral to rgrid ...")') trim(varname); call goPr
        ! define as copy:
        call Spg_Define_Variable_RegriddedCopy_Rgg( spgfile, varname, 'rgrid', rggfile, varids_spectral(ivar), status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! info ...
        write (gol,'("define variable ",a," ...")') trim(varname); call goPr
        ! define copy:
        call Spg_Define_Variable_Copy_Rgg( spgfile, varname, rggfile, varids(ivar), status )
        IF_NOT_OK_RETURN(status=1)
      end if
    end do ! variables:
    
    ! end of definition phase::
    call rggfile%EndDef( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write grid definition:
    call rggfile%Put_Grid( rgg, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy data of non-spectral variables:
    do ivar = 1, nvar
      ! not defined ? then skip:
      if ( varids(ivar) < 0 ) cycle
      ! get name:
      call spgfile%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'("copy data of ",a," ...")') trim(varname); call goPr
      ! copy selection:
      call Spg_Copy_Variable_Rgg( spgfile, varname, rggfile, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! evaluate spectral fields:
    do ivar = 1, nvar
      ! not defined ? then skip:
      if ( varids_spectral(ivar) < 0 ) cycle
      ! get name:
      call spgfile%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'("evaluate spectral field ",a," ...")') trim(varname); call goPr
      ! copy selection:
      call Spg_Regrid_Variable_Rgg( spgfile, spg, varname, rggfile, rgg, status )
      IF_NOT_OK_RETURN(status=1)
    end do  ! variables
    
    ! info ...
    write (gol,'("done ...")'); call goPr
    
    ! clear:
    deallocate( varids )
    deallocate( varids_spectral )
    
    ! close file:
    call rggfile%Close( status )
    IF_NOT_OK_RETURN(status=1) 

    ! clear:
    call spg%Final()
    
    ! close file:
    call spgfile%Close( status )
    IF_NOT_OK_RETURN(status=1)    
    
    ! ok
    status = 0
  
  end subroutine SPG2RGG


end module C3PO_SPG2RGG_Main
