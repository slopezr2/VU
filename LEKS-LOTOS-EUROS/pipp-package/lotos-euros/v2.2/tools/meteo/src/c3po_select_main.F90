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

module C3PO_Select_Main

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  ParseArguments
  public  ::  ExtractDomain
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Select_Main'
  
  

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
  
  subroutine ParseArguments( domain, input_filename, output_filename, status )
  
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    real, intent(out)                ::  domain(4) ! (w,e,s,n) [deg]
    character(len=*), intent(out)    ::  input_filename
    character(len=*), intent(out)    ::  output_filename
    integer, intent(out)             ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ParseArguments'
    
    ! --- local ----------------------------------
    
    integer                         ::  narg, iarg
    integer                         ::  arglength
    character(len=1024)             ::  argvalue
    integer                         ::  k
    
    ! --- begin ----------------------------------
    
    ! dummy values at start:
    domain = -999.9
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
          write (gol,'("  c3po_select.x west,easth,south,north inputfile outputfile")'); call goPr
          write (gol,'("  c3po_select.x -h|--help")'); call goPr
          status = -1; return
        case default
          if ( any( domain < -900.0 ) ) then
            ! extract domain:
            do k = 1, size(domain)
              call goReadFromLine( argvalue, domain(k), status, sep=',' )
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
  
  subroutine ExtractDomain( domain, input_filename, output_filename, status )
  
    use C3PO_File_Rgg, only : T_File_Rgg    
    use C3PO_Grid_RGG, only : T_Grid_RGG
  
    ! --- in/out ---------------------------------
    
    real, intent(in)                ::  domain(4) ! (w,e,s,n) [deg]
    character(len=*), intent(in)    ::  input_filename
    character(len=*), intent(in)    ::  output_filename
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ExtractDomain'
    
    ! --- local ----------------------------------
    
    logical                 ::  exist
    type(T_File_Rgg)         ::  rggfile_in
    type(T_File_Rgg)         ::  rggfile_out
    type(T_Grid_RGG)        ::  rgg
    type(T_Grid_RGG)        ::  rgg_out
    integer                 ::  nselected
    integer, allocatable    ::  selection(:)
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
    call rggfile_in%Open( trim(input_filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("read grid definition ...")'); call goPr
    
    ! fill grid definition from file variables:
    call rggfile_in%Get_Grid( rgg, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("select domain ...")'); call goPr
    
    ! storage for selection indices:
    allocate( selection(rgg%npoint) )
    ! extract definition:
    call rgg%GetRegion( domain, rgg_out, nselected, selection, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("create target file ",a," ...")') trim(output_filename); call goPr
    
    ! create target file:
    call rggfile_out%Create( trim(output_filename), status, formatNum=rggfile_in%formatNum )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("copy attributes ...")'); call goPr
    
    ! copy global attributes:
    call rggfile_in%Copy_Atts( '', rggfile_out, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("extend history ...")'); call goPr
    
    ! write informative message:
    write (message,'("c3po: extracted region [",f7.2,",",f7.2,"] x ["f6.2,",",f6.2,"] from ",a)') &
                            domain, trim(input_filename)
    ! add to history attribute:
    call rggfile_out%Extend_History( trim(message), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("write grid definition ...")'); call goPr
    
    ! define grid definition from file variables:
    call rggfile_out%Def_Grid( rgg_out, status )
   IF_NOT_OK_RETURN(status=1)
    
    ! number of dimensions:
    call rggfile_in%Inquire( status, ndim=ndim )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do idim = 1, ndim
      ! get name:
      call rggfile_in%Inquire_Dimension( idim, status, name=dimname )
      IF_NOT_OK_RETURN(status=1)
      ! skip grid dimensions:
      if ( trim(dimname) == 'ulat'  ) cycle
      if ( trim(dimname) == 'ulon'  ) cycle
      if ( trim(dimname) == 'rgrid' ) cycle
      if ( trim(dimname) == 'nv'    ) cycle
      ! info ...
      write (gol,'("copy dimension ",a," ...")') trim(dimname); call goPr
      ! other:
      call rggfile_in%Copy_Dimension( trim(dimname), rggfile_out, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! number of variables:
    call rggfile_in%Inquire( status, nvar=nvar )
    IF_NOT_OK_RETURN(status=1)
    ! storage for variable id's in output file:
    allocate( varids(nvar) )
    ! init with dummy values, used later on to check on need to copy:
    varids = -999
    ! loop:
    do ivar = 1, nvar
      ! get name:
      call rggfile_in%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! skip grid variables:
      if ( trim(varname) == 'rgrid'     ) cycle
      if ( varname(1:9)  == 'longitude' ) cycle
      if ( varname(1:8)  == 'latitude'  ) cycle
      if ( varname(1:4)  == 'ulat'      ) cycle
      ! info ...
      write (gol,'("define variable copy ",a," ...")') trim(varname); call goPr
      ! define copy:
      call rggfile_in%Define_Variable_Copy( varname, rggfile_out, varids(ivar), status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! write grid definition:
    call rggfile_out%EndDef( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write grid definition:
    call rggfile_out%Put_Grid( rgg_out, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy data of variables:
    do ivar = 1, nvar
      ! not defined ? then skip:
      if ( varids(ivar) < 0 ) cycle
      ! get name:
      call rggfile_in%Inquire_Variable( ivar, status, name=varname )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'("copy selection of ",a," ...")') trim(varname); call goPr
      ! copy selection:
      call rggfile_in%Copy_Variable_Selection( varname, rggfile_out, &
                                                selection(1:nselected), status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! info ...
    write (gol,'("done ...")'); call goPr
    
    ! clear:
    deallocate( varids )
    
    ! close file:
    call rggfile_out%Close( status )
    IF_NOT_OK_RETURN(status=1)  
    
    ! clear grid defintion:
    call rgg_out%Done( status )
    IF_NOT_OK_RETURN(status=1)
     
    ! clear grid defintion:
    call rgg%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( selection, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close file:
    call rggfile_in%Close( status )
    IF_NOT_OK_RETURN(status=1)    
    
    ! ok
    status = 0
  
  end subroutine ExtractDomain


end module C3PO_Select_Main
