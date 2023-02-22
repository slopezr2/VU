!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################

module LE_Emis_TNO_Timeprof_Selection

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  T_Emis_Time_Prof_Selection
  public  ::  LE_Emis_TNO_Timeprof_Selection_Init
  public  ::  LE_Emis_TNO_Timeprof_Selection_Done
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_TNO_Timeprof_Selection'
  
  
  ! --- types --------------------------------
  
  ! Time profile selection
  type T_Emis_Time_Prof_Selection
    integer                           ::  ncat
    integer                           ::  nfuel
    integer                           ::  ntracer
    character(len=1024), allocatable  ::  tprof_fname(:,:)  ! (cat,tracer)
    character(len=32), allocatable    ::  tprof_type(:,:)   ! (cat,tracer)
    logical, allocatable              ::  filled(:,:)       ! (cat,tracer)

    character(len=512), allocatable   ::  fname_hourly(:)
    character(len=512), allocatable   ::  fname_hourly_gridded(:)
    integer                           ::  ntype_hourly
    integer                           ::  ntype_hourly_gridded
  end type


contains



  ! ===============================================================
  ! ===
  ! === module init/done
  ! ===
  ! ===============================================================
  

  subroutine LE_Emis_TNO_Timeprof_Selection_Init( em_tps, cat_codes, tracer_names, query, year, status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue
    use GO, only :  goReplace
    use Go, only :  calc_DayNumber, days_in_year
    use MDF, only : MDF_Open, MDF_Close
    use MDF, only : MDF_Inquire
    use MDF, only : MDF_Inq_DimID, MDF_Inquire_Dimension
    use MDF, only : MDF_Inq_VarID, MDF_Inquire_Variable, MDF_Get_Var
    use MDF, only : MDF_Get_Att
    use MDF, only : MDF_NETCDF, MDF_READ
    
    use LE_Emis_Tools, only : ShortSNAP_to_Code
    use LE_Emis_Tools, only : MDF_Get_StrArr

    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof_Selection), intent(out) ::  em_tps
    character(len=*), intent(in)                  ::  cat_codes(:)
    character(len=*), intent(in)                  ::  tracer_names(:)
    character(len=512), intent(in)                ::  query
    integer, intent(out)                          ::  year
    integer, intent(out)                          ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Timeprof_Selection_Init'
    integer, parameter          ::  maxcol = 10
    
    integer, parameter          ::  max_hourly_types           = 50
    integer, parameter          ::  max_hourly_gridded_types   = 30
    
    ! --- local -------------------------------
        
    character(len=512)      ::  filename
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    logical                 ::  exist

    integer                 ::  fu
    character(len=1024)     ::  line
    integer                 ::  iline
    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    integer                 ::  nfield
    integer                 ::  ifield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    
    integer                 ::  ifield_cat, ifield_tracer
    integer                 ::  ifield_file, ifield_type
    
    integer                 ::  icat, itrac
    integer                 ::  idummy
    
    ! --- begin -------------------------------
    
    ! store currently known dimensions:
    em_tps%ncat    = size(cat_codes)
    em_tps%ntracer = size(tracer_names)
    
    ! init:
    allocate( em_tps%tprof_fname( em_tps%ncat, em_tps%ntracer), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( em_tps%tprof_type( em_tps%ncat, em_tps%ntracer), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( em_tps%filled( em_tps%ncat, em_tps%ntracer), stat=status )
    IF_NOTOK_RETURN(status=1)
    em_tps%filled = .false.
    
    ! Filenames
    allocate( em_tps%fname_hourly(max_hourly_types), stat=status)
    IF_NOTOK_RETURN(status=1)      
    allocate( em_tps%fname_hourly_gridded(max_hourly_gridded_types), stat=status)
    IF_NOTOK_RETURN(status=1)  

    ! set counters
    em_tps%ntype_hourly = 0
    em_tps%ntype_hourly_gridded = 0
  
    ! extract filename, by default empty:
    filename = ''
    call goVarValue( trim(query), ';', 'file', '=', filename, status )
    IF_ERROR_RETURN(status=1)
    
    ! replace year string if necessary
    call goReplace( filename, '%{year}', '(i4.4)', year, status )
    IF_NOTOK_RETURN(status=1)
    
    ! seperation character:
    sep = ';'
    call goVarValue( trim(query), ';', 'sep', '=', sep, status )
    IF_ERROR_RETURN(status=1)
    ! comment character:
    comment = '#'
    call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    IF_ERROR_RETURN(status=1)

    ! file should be present:
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( fu, file=trim(filename), status='old', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    iline = 0
    ! read header line after first comment
    do 
      ! read line:
      read (fu, '(a)', iostat=status) line
      if (status/=0 ) then
        write( gol, '("Reading header line from file: ",a)') trim(filename); call goErr
        TRACEBACK;status=1;return
      end if
      ! skip empty line:
      if (len_trim(line) == 0 ) cycle
      ! skip comment line
      if (line(1:1) == comment ) cycle
      ! found non comment line (must be header)
      exit
    end do
    ! split:
    call goSplitString( line, nheader, headers, status, sep=sep )
    IF_NOTOK_RETURN(status=1)
    
    ! set header indices
    ifield_cat    = -1
    ifield_tracer = -1
    ifield_file   = -1
    ifield_type   = -1
    ! loop over fields
    do ifield = 1, nheader
      ! current
      header = headers(ifield)
      ! which column?
      select case ( trim(header) )
        ! emission category
        case ( 'Sector_Fuel' )
          ifield_cat = ifield
        ! emission tracer
        case( 'Tracer', 'tracer' )
          ifield_tracer = ifield
        ! profile file
        case( 'ProfileFile', 'profilefile' )
          ifield_file = ifield
        ! profile type
        case ( 'ProfileType', 'profiletype' )
          ifield_type = ifield
        ! Unknown
        case default
          write( gol, '("Unsupported header found: `", a,"`")') trim(header); call GoErr
          TRACEBACK;status=1;return
      end select
    end do
    
    ! check all headers found
    if ( ifield_cat < 0 ) then
      write( gol, '("Emission category index not defined yet")' ) ; call GoErr
      TRACEBACK;status=1;return
    end if
    if ( ifield_tracer < 0 ) then
      write( gol, '("Emission tracer index not defined yet")' ) ; call GoErr
      TRACEBACK;status=1;return
    end if
    if ( ifield_file < 0 ) then
      write( gol, '("Emission time profile file index not defined yet")' ) ; call GoErr
      TRACEBACK;status=1;return
    end if
    if ( ifield_type < 0 ) then
      write( gol, '("Emission time profile type index not defined yet")' ) ; call GoErr
      TRACEBACK;status=1;return
    end if
    
    ! Loop over records
    do
      ! increase record counter
      iline = iline + 1
      ! try to read line
      read (fu,'(a)', iostat=status) line
      ! eof?
      if ( status < 0 ) exit
      ! error
      if ( status > 0 ) then
        write(gol, '("reading line ",i6," from file: ",a)') iline, trim(filename); call GoErr
        TRACEBACK;status=1;return
      end if
      
      ! skip empty line:
      if (len_trim(line) == 0 ) cycle
      ! skip comment line
      if (line(1:1) == comment ) cycle
            
      ! split into records
      call GoSplitString( line, nfield, fields, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      ! check number of fields
      if ( nfield /= nheader ) then
        write (gol,'("number of fields (",i2,") in line ",i2," :")') nfield, iline; call goPr
        write (gol,'("  ",a)') trim(line); call goErr
        write (gol,'("fields:")'); call goErr
        do ifield = 1, nfield
          write (gol,'(i6," : ",a)') ifield, trim(fields(ifield)); call goErr
        end do
        write (gol,'("differs from number of headers ",i2," in file ",a)') nheader, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
    
      ! get category index
      call GoMatchValue( trim(fields(ifield_cat)), cat_codes, icat, status )
      IF_NOTOK_RETURN(status=1)
      ! get tracer index
      ! Specific case for PM25 and PM10
      if ( trim(fields(ifield_tracer)) == 'PM25' ) then
        call GoMatchValue( 'pm2_5', tracer_names, itrac, status )
        IF_NOTOK_RETURN(status=1)
      else if ( trim(fields(ifield_tracer)) == 'PM10' ) then
        call GoMatchValue( 'pm25_pm10', tracer_names, itrac, status )
        IF_NOTOK_RETURN(status=1)
      else          
        call GoMatchValue( trim(fields(ifield_tracer)), tracer_names, itrac, status )
        IF_NOTOK_RETURN(status=1)
      end if
      
      ! fill profile file/type
      em_tps%tprof_fname( icat,itrac) = trim(fields(ifield_file))
      em_tps%tprof_type( icat,itrac)  = trim(fields(ifield_type))
      
      ! COunt number of profiles per type
      select case ( trim(fields(ifield_type)) )
        
        case ( 'hourly' )
        
          ! ~ hourly timeprofile
          ! Check if already available
          call goMatchValue( trim(fields(ifield_file)), em_tps%fname_hourly(1:em_tps%ntype_hourly), idummy, status, quiet=.true. )
          IF_ERROR_RETURN(status=1)
          if ( status < 0 ) then
            em_tps%ntype_hourly = em_tps%ntype_hourly + 1
            em_tps%fname_hourly(em_tps%ntype_hourly) = trim(fields(ifield_file))
          end if
          
        case ( 'hourly-gridded' )
          
          ! ~ hourly gridded timeprodile
          ! Check if already available
          call goMatchValue( trim(fields(ifield_file)), em_tps%fname_hourly_gridded(1:em_tps%ntype_hourly_gridded), idummy, status, quiet=.true. )
          IF_ERROR_RETURN(status=1)
          if ( status < 0 ) then
            em_tps%ntype_hourly_gridded = em_tps%ntype_hourly_gridded + 1
            em_tps%fname_hourly_gridded(em_tps%ntype_hourly_gridded) = trim(fields(ifield_file))
          end if
          
        case default
          ! ~ Unknown profiles
          write( gol, '("Unknown timeprofile type found: ",a)' ) trim(fields(ifield_type)) ; call GoErr
          write( gol, '("In selection table: ",a)' ) trim(filename) ; call GoErr
          write( gol, '("Line number: ",i6)' ) iline ; call GoErr
          TRACEBACK;status=1;return
      end select
      
      ! Put flag 
      em_tps%filled(icat,itrac) = .true.
      
    end do  ! lines
      
    ! Check if all filles
    do icat = 1, em_tps%ncat
    do itrac = 1, em_tps%ntracer
      if ( .not. em_tps%filled(icat,itrac) ) then
        write(gol, '("Not all sectors/tracers fille for time prof selection: ", a," ",a)' ) trim(cat_codes(icat)), trim(tracer_names(itrac)) ; call GoErr
        TRACEBACK;status=1;return
      end if
    end do
    end do
    
    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Timeprof_Selection_Init
  
  
  ! ***
  

  subroutine LE_Emis_TNO_Timeprof_Selection_Done( em_tps, status )
  
    ! --- in/out ------------------------------
    
    type(T_Emis_Time_Prof_Selection), intent(out) ::  em_tps
    integer, intent(out)                  ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_TNO_Timeprof_Selection_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    ! clear:
    if ( allocated(em_tps%tprof_fname) ) deallocate( em_tps%tprof_fname  )
    if ( allocated(em_tps%tprof_type) ) deallocate( em_tps%tprof_type  )
    
    if ( allocated(em_tps%fname_hourly) ) deallocate( em_tps%fname_hourly )
    if ( allocated(em_tps%fname_hourly_gridded) ) deallocate( em_tps%fname_hourly_gridded )
    
    ! ok
    status = 0

  end subroutine LE_Emis_TNO_Timeprof_Selection_Done


end module LE_Emis_TNO_Timeprof_Selection

