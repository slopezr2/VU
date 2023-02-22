!###############################################################################
!
! NAME
!   LE_Logging  -  LOTOS-EUROS logging stuff
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

module LE_Logging

  use GO, only : gol, goPr, goErr
  use GO, only : goStdOut

  implicit none


  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Logging_Init, LE_Logging_Done
  public  ::  LE_Logging_SetFU
  public  ::  PrintDebug

  public  ::  ident1, ident2
  public  ::  u_log
  !public  ::  u_log_c
  public  ::  u_err


  ! --- const ------------------------------------
    
  character(len=*), parameter   ::  mname = 'LE_Logging'
  
  ! identations
  character(len=3), parameter :: ident1 = '   '
  character(len=6), parameter :: ident2 = '      '
  

  ! --- var --------------------------------------

  ! file units, set both to standard output:
  integer     ::  u_log = goStdOut
  integer     ::  u_err = goStdOut


contains


  ! ====================================================================
  
  
  subroutine LE_Logging_Init( rcF, status )
  
    use GO       , only : lgr
    use GO       , only : goc
    use GO       , only : TrcFile, ReadRc
    !use GO       , only : goGetFU
    use GO       , only : pathsep
    use C3PO     , only : show_filename_on_open
    use Dims     , only : outF
    use LE_Data  , only : LE_Data_Enable
    use LE_Config, only : outputdir
      
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)         ::  rcF
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Logging_Init'
    
    ! --- local ----------------------------------

    logical              ::  log2files
    logical              ::  root2stdout
    character(len=3)     ::  yesno
    character(len=64)    ::  values
    character(len=1024)  ::  filename

    ! --- begin ----------------------------------
    
    ! enable on all processes:
    call lgr%Set( status, apply=.true. )
    IF_NOTOK_RETURN(status=1)
    
    ! write to log files?   
    call ReadRc( rcF, 'log.files', log2files, status )
    IF_NOTOK_RETURN(status=1)
    ! apply?
    if ( log2files ) then
      ! echo root ?
      call ReadRc( rcF, 'log.files.root_to_stdout', root2stdout, status )
      IF_NOTOK_RETURN(status=1)
      ! log file name, including domain number if decomposed:
      if ( goc%npes > 1 ) then
        write (filename,'(a,a,"le-log-d",i2.2,".out")') trim(outputdir), pathsep, goc%id
      else
        write (filename,'(a,a,"le-log.out")') trim(outputdir), pathsep
      end if
      ! switch to log file ;
      ! eventually let root also write to standard output:
      call lgr%LogFile( status, file=trim(filename), echo=(goc%root .and. root2stdout) )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! *
    
    call ReadRc( rcF, 'log.suppress', yesno, status )
    IF_NOTOK_RETURN(status=1)
    call process_yesno(yesno,outF%suppress, &
                       'suppressing most of the LOTOS output', &
                       'printing all LOTOS output')

    call ReadRc( rcF, 'log.debug', yesno, status )
    IF_NOTOK_RETURN(status=1)
    call process_yesno(yesno,outF%debug_print, &
                       'printing debug infomation (per process)', &
                       'do not debug infomation')
    ! enabled?
    if (outF%debug_print) then
      ! selected cell:
      call ReadRc( rcF, 'log.debug.cell', values, status )
      IF_NOTOK_RETURN(status=1)
      read (values,*) outF%debug_ix, outF%debug_iy
      ! enable meteo:
      call LE_Data_Enable( 'dh', status )
      IF_NOTOK_RETURN(status=1)
      !call LE_Data_Enable( 'vol', status )
      !IF_NOTOK_RETURN(status=1)
    endif

    !! output directory:
    !call ReadRc( rcF, 'le.output.outdir', outpath, status )
    !IF_NOTOK_RETURN(status=1)

    !! full path:
    !write (filename,'(a,a,"log.txt")') trim(outpath), pathsep
    !! new file unit:
    !call goGetFU( u_log, status )
    !IF_NOTOK_RETURN(status=1)
    !! open:
    !open( u_log, file=trim(filename), iostat=status )
    !if (status/=0) then
    !  write (gol,'("openening file `",a,"`")') trim(filename); call goErr
    !  TRACEBACK; status=1; return
    !end if

    !! full path:
    !write (filename,'(a,a,"log_chem.txt")') trim(outpath), pathsep
    !! new file unit:
    !call goGetFU( u_log_c, status )
    !IF_NOTOK_RETURN(status=1)
    !! open:
    !open( u_log, file=trim(filename), iostat=status )
    !if (status/=0) then
    !  write (gol,'("openening file `",a,"`")') trim(filename); call goErr
    !  TRACEBACK; status=1; return
    !end if
    
    ! write message when file is opened with c3po_file_nc ?
    call ReadRc( rcF, 'log.show_filename_on_open', show_filename_on_open, status )
    IF_NOTOK_RETURN(status=1)

    ! set file units:
    call LE_Logging_SetFU( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Logging_Init
  
  ! ***
  
  
  subroutine LE_Logging_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Logging_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    !! info ..
    !write (u_log,*) 'model simulation finished'

    !! close logfile:
    !close( u_log, iostat=status )
    !IF_NOTOK_RETURN(status=1)

    !! close logfile:
    !close( u_log_c, iostat=status )
    !IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Logging_Done


  ! ***
  
  
  subroutine LE_Logging_SetFU( status )
  
    use GO, only : lgr

    ! --- in/out ---------------------------------
    
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Logging_SetFU'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! copy file units:
    u_log = lgr%PrFU()
    u_err = lgr%ErrFU()
    
    ! ok
    status = 0
  
  end subroutine LE_Logging_SetFU


  ! ***


  subroutine process_yesno( yesno, flag, stringyes, stringno )
  
    ! --- in/out ---------------------------------

    character(len=3), intent(in)  ::  yesno
    logical, intent(out)          ::  flag
    character(len=*), intent(in)  ::  stringyes
    character(len=*), intent(in)  ::  stringno
    
    ! --- begin ----------------------------------

    ! set flag:
    select case ( yesno )
      case ( 'yes', 'YES', 'y', 'Y' )
        flag = .true.
      case default
        flag = .false.
    end select

    ! info ..
    if (flag) then
      write (gol,'(a)') trim(stringyes); call goPr
    else
      write (gol,'(a)') trim(stringno); call goPr
    endif

  end subroutine process_yesno


  ! ***


  ! print concentration to screen for debugging purposes

  subroutine PrintDebug( c, ca, name, status )

    use dims         , only : outF
    use Dims         , only : nx, ny, nz, nspec
    use Indices      , only : specname
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Meteo_Data, only : ovolume, volume
    use LE_Meteo_Data, only : oairmass, airmass
    use LE_Meteo_Data, only : ohpres
    
    !use LE_Grid, only : debug_i, debug_j
    !use LE_Grid, only : debug_vars
    !use indices, only : ispec_o3
    !use GO, only : goReadFromLine

    ! --- in/out ---------------------------------

    real, intent(in)               ::  c(nx,ny,nz,nspec)  ! concentration
    real, intent(in)               ::  ca(nx,ny,nspec)  ! concentration aloft
    character(len=*), intent(in)   ::  name      ! name of related subroutine (PrintDebug is called after this routine)
    integer, intent(out)           ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/PrintDebug'    

    ! --- local ----------------------------------
    
    integer            :: ix,iy,k,ispec
    character(len=128) :: fmt
    character(len=10)  :: nam1
    real               :: hh(nz)
    
    ! meteo data:
    real, pointer                         ::  dh(:,:,:)
    !real, pointer                         ::  pdata(:,:,:)
    !character(len=16)                     ::  varname
    !character(len=len_trim(debug_vars))   ::  line
    
    ! --- begin ----------------------------------

    !! useful for OpenMP debugging ...
    !write (gol,*) '        xxx ', name, sum(c), minval(c), maxval(c); call goPr
    
    !! idem per component:
    !do ispec = 1, nspec
    !  write (gol,*) '        xxx ', name, ' ', trim(specname(ispec)), sum(c(:,:,:,ispec)), minval(c(:,:,:,ispec)), maxval(c(:,:,:,ispec)); call goPr
    !end do
    
    !! testing decompo ...
    !do ix = 1, nx
    !  if ( ix == debug_i ) then
    !    write (gol,*) 'xxx '; call goPr
    !    write (gol,*) 'xxx ', trim(name); call goPr
    !    do iy = ny, 1, -1
    !      write (gol,*) 'xxx ', iy, c(ix,iy,1,ispec_o3); call goPr
    !    end do
    !  end if
    !end do
    
!    ! testing adjust ...
!    do iy = 1, ny
!      do ix = 1, nx
!        if ( (ix == debug_i) .and. (iy == debug_j) ) then
!
!          ! info ...
!          write (gol,*) 'xxx '; call goPr
!          write (gol,*) 'xxx ', trim(name); call goPr
!
!          ! concentrations:
!          do ispec = 1, nspec
!            write (gol,*) 'xxx c ', trim(specname(ispec)), c(ix,iy,:,ispec), ca(ix,iy,ispec); call goPr
!          end do
!
!          ! module variables:
!          write (gol,*) 'xxx oairmass      ', oairmass(ix,iy,:); call goPr
!          write (gol,*) 'xxx  airmass      ',  airmass(ix,iy,:); call goPr
!          write (gol,*) 'xxx ovolume       ', ovolume(ix,iy,:); call goPr
!          write (gol,*) 'xxx  volume       ',  volume(ix,iy,:); call goPr
!          write (gol,*) 'xxx ohpres        ', ohpres(ix,iy,:); call goPr 
!          
!          ! data:
!          line = trim(debug_vars)
!          do
!            ! empty?
!            if ( len_trim(line) == 0 ) exit
!            ! next:
!            call goReadFromLine( line, varname, status, sep=' ' )
!            IF_NOTOK_RETURN(status=1)
!            ! set pointer to data:
!            call LE_Data_GetPointer( trim(varname), pdata, status )    
!            IF_NOTOK_RETURN(status=1)
!            ! show:
!            write (gol,*) 'xxx ', varname, pdata(ix,iy,:); call goPr
!          end do
!          write (gol,*) 'xxx '; call goPr
!
!        end if
!      end do
!    end do
    
    ! .....................................................

    nam1 = name
    if (outF%debug_print) then

      ! point to meteo data:
      call LE_Data_GetPointer( 'dh', dh, status, check_units ='m' )
      IF_NOTOK_RETURN(status=1)
      !call LE_Data_GetPointer( 'vol', volume, status, check_units ='m3' )
      !IF_NOTOK_RETURN(status=1)
    
      ! selected cell:
      ix = outF%debug_ix 
      iy = outF%debug_iy 

      ! cell output or grid statistics?
      if ( (ix == 0) .or. (iy == 0) ) then

        ! print min,max,average over each layer
        fmt = '(a,a8,1x,a10,a,xx(1x,e10.3))'
        write(fmt(16:17),'(i2.2)') nz
        do ispec = 1, nspec
          !write (gol,fmt) 'DB ', specname(ispec), nam1, ' min ', (minval(c(:,:,k,ispec)), k = 1,nz); call goPr
          !write (gol,fmt) 'DB ', specname(ispec), nam1, ' max ', (maxval(c(:,:,k,ispec)), k = 1,nz); call goPr
          write (gol,fmt) 'DB ', specname(ispec), nam1, ' ave ', (sum(c(:,:,k,ispec))/(nx*ny), k = 1,nz); call goPr
        end do

      else

        ! print concentration at cell 
        if ( (ix < 1) .or. (ix > nx) .or. (iy < 1) .or. (iy > ny) ) then
          write (gol,'("wrong value of debug cell ix,iy  = ",i0," ",i0)') ix,iy; call goErr
          TRACEBACK; status=1; return
        end if
        fmt = '(a,a8,1x,a10,a,2i5,xx(1x,e10.3))'
        write(fmt(20:21),'(i2.2)') nz
        ! compute the depths of the cells
        hh = dh(ix,iy,:)  ! m
        ! loop over specs:
        do ispec = 1,nspec
          write (gol,fmt) 'DB ', specname(ispec), nam1, ' hh___cell ', ix,iy, (hh(k)                           , k = 1,nz); call goPr
          write (gol,fmt) 'DB ', specname(ispec), nam1, ' conc_cell ', ix,iy, (c(ix,iy,k,ispec)                , k = 1,nz); call goPr
          write (gol,fmt) 'DB ', specname(ispec), nam1, ' mass_cell ', ix,iy, (volume(ix,iy,k)*c(ix,iy,k,ispec), k = 1,nz); call goPr
        end do
          
       end if ! cell output?

    end if !debug_print
    
    ! ok
    status = 0

  end subroutine PrintDebug

end module LE_Logging


