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

module C3PO_CoarseLevs_Main

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  ParseArguments
  public  ::  CoarsenLevels
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_CoarseLevs_Main'
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** arguments
  ! ***
  ! ********************************************************************


  !
  ! Expected command line arguments:
  !     level  hlevel  3/3/3/3/3/1  input.nc output.nc
  !
  ! Return status:
  !     -1 : usage displayed
  !      0 : ok
  !  other : error
  !
  
  subroutine ParseArguments( levelname, hlevelname, n, input_filename, output_filename, status )
  
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    integer, allocatable             ::  n(:)
    character(len=*), intent(out)    ::  levelname
    character(len=*), intent(out)    ::  hlevelname
    character(len=*), intent(out)    ::  input_filename
    character(len=*), intent(out)    ::  output_filename
    integer, intent(out)             ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ParseArguments'
    
    ! --- local ----------------------------------
    
    integer                         ::  narg, iarg
    integer                         ::  arglength
    character(len=1024)             ::  argvalue
    integer                         ::  nval, k
    
    ! --- begin ----------------------------------
    
    ! dummy values at start:
    levelname = ''
    hlevelname = ''
    if ( allocated(n) ) deallocate( n )
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
        TRACEBACK; stop 1
      end if
      ! get value: 
      call Get_Command_Argument( iarg, value=argvalue, status=status )
      if ( status /= 0 ) then
        write (gol,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength; call goErr
        TRACEBACK; stop 1
      end if
      ! switch:
      select case ( trim(argvalue) )
        case ( '-h', '--help' )
          write (gol,'("Usage:")'); call goPr
          write (gol,'("  c3po_coarselevs.x level hlevel 4/4/4/1 inputfile outputfile")'); call goPr
          write (gol,'("  c3po_coarselevs.x -h|--help")'); call goPr
          status = -1; return
        case default
          if ( len_trim(levelname) == 0 ) then
            levelname = trim(argvalue)
          else if ( len_trim(hlevelname) == 0 ) then
            hlevelname = trim(argvalue)
          else if ( .not. allocated(n) ) then
            ! count:
            nval = 1
            do k = 1, len_trim(argvalue)
              if ( argvalue(k:k) == '/' ) nval = nval + 1
            end do
            ! storage:
            allocate( n(nval), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! read:
            do k = 1, nval
              call goReadFromLine( argvalue, n(k), status, sep='/' )
              IF_NOT_OK_RETURN(status=1)
            end do
          else if ( len_trim(input_filename) == 0 ) then
            input_filename = trim(argvalue)
          else if ( len_trim(output_filename) == 0 ) then
            output_filename = trim(argvalue)
          else
            write (*,'("unsupported argument `",a,"`")') trim(argvalue)
            TRACEBACK; stop 1
          end if
      end select
    end do   ! arguments

    ! check ...
    if ( len_trim(input_filename) == 0 ) then
      write (*,'("no input filename argument specified")')
      TRACEBACK; stop 1
    end if
    ! check ...
    if ( len_trim(output_filename) == 0 ) then
      write (*,'("no output filename argument specified")')
      TRACEBACK; stop 1
    end if
    
    ! ok
    status = 0
  
  end subroutine ParseArguments
  
  
  ! ***
  
  subroutine CoarsenLevels( levelname, hlevelname, nn, input_filename, output_filename, status )
  
    use C3PO_File_Nc , only : T_File_Nc  
    use C3PO_Levs_Hyb, only : T_Levs_Hyb
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)    ::  levelname
    character(len=*), intent(in)    ::  hlevelname
    integer, intent(in)             ::  nn(:)
    character(len=*), intent(in)    ::  input_filename
    character(len=*), intent(in)    ::  output_filename
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/CoarsenLevels'
    
    ! --- local ----------------------------------
    
    logical                 ::  exist
    type(T_File_Nc)         ::  ncfile_in
    type(T_File_Nc)         ::  ncfile_out
    type(T_Levs_Hyb)        ::  levs_in
    type(T_Levs_Hyb)        ::  levs_out
    integer, allocatable    ::  mapping(:,:)  ! (nlev,2)
    character(len=256)      ::  message
    integer                 ::  ndim, idim
    character(len=32)       ::  dimname
    integer                 ::  nvar, ivar
    character(len=32)       ::  varname
    integer, allocatable    ::  varids(:)
    
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
    call ncfile_in%Open( trim(input_filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("read level definition ...")'); call goPr
    
    ! fill grid definition from file variables:
    call ncfile_in%Get_Levs( levelname, hlevelname, levs_in, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("coarse level definition ...")'); call goPr
    ! init output level definition as coarsened version of input:
    call levs_out%Init( levs_in, nn, status, mapping=mapping )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("create target file ",a," ...")') trim(output_filename); call goPr
    
    ! create target file:
    call ncfile_out%Create( trim(output_filename), status, formatNum=ncfile_in%formatNum )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("copy attributes ...")'); call goPr
    
    ! copy global attributes:
    call ncfile_in%Copy_Atts( '', ncfile_out, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("extend history ...")'); call goPr
    
    ! write informative message:
    write (message,'("c3po: coarsened to ",i0," levels from ",a)') &
                            size(nn), trim(input_filename)
    ! add to history attribute:
    call ncfile_out%Extend_History( trim(message), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("write level definition ...")'); call goPr
    
    ! define levels definition in ouput file:
    call ncfile_out%Def_Levs( levs_out, levelname, hlevelname, status )
   IF_NOT_OK_RETURN(status=1)
    
    ! number of dimensions:
    call ncfile_in%Inquire( status, ndim=ndim )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do idim = 1, ndim
      ! get name:
      call ncfile_in%Inquire_Dimension( idim, status, name=dimname )
      IF_NOT_OK_RETURN(status=1)
      ! skip level dimensions:
      if ( trim(dimname) == trim(levelname ) ) cycle
      if ( trim(dimname) == trim(hlevelname) ) cycle
      ! info ...
      write (gol,'("copy dimension ",a," ...")') trim(dimname); call goPr
      ! other:
      call ncfile_in%Copy_Dimension( trim(dimname), ncfile_out, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! number of variables:
    call ncfile_in%Inquire( status, nvar=nvar )
    IF_NOT_OK_RETURN(status=1)
    ! storage for variable id's in output file:
    allocate( varids(nvar) )
    ! init with dummy values, used later on to check on need to copy:
    varids = -999
    ! loop:
    do ivar = 1, nvar
      ! get name:
      call ncfile_in%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! skip level variables:
      if ( varname      == trim(levelname ) ) cycle
      if ( varname      == trim(hlevelname) ) cycle
      if ( varname(1:3) == 'hya'            ) cycle
      if ( varname(1:3) == 'hyb'            ) cycle
      ! info ...
      write (gol,'("define variable copy ",a," ...")') trim(varname); call goPr
      ! define copy, the lenght of the level dimension is different now:
      call ncfile_in%Define_Variable_Copy( varname, ncfile_out, varids(ivar), status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! end of definition phase:
    call ncfile_out%EndDef( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write level definition:
    call ncfile_out%Put_Levs( levs_out, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy data of variables:
    do ivar = 1, nvar
      ! not defined ? then skip:
      if ( varids(ivar) < 0 ) cycle
      ! get name:
      call ncfile_in%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'("coarsen ",a," ...")') trim(varname); call goPr
      ! copy selection:
      call ncfile_in%Copy_Variable_CoarsenLevs( varname, ncfile_out, &
                        levelname, hlevelname, levs_in, levs_out, mapping, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! info ...
    write (gol,'("done ...")'); call goPr
    
    ! clear:
    deallocate( varids )
    
    ! close file:
    call ncfile_out%Close( status )
    IF_NOT_OK_RETURN(status=1)  
    
    ! clear level defintion:
    call levs_out%Done( status )
    IF_NOT_OK_RETURN(status=1)
     
    ! clear level defintion:
    call levs_in%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close file:
    call ncfile_in%Close( status )
    IF_NOT_OK_RETURN(status=1) 
    
    ! ok
    status = 0
  
  end subroutine CoarsenLevels


end module C3PO_CoarseLevs_Main
