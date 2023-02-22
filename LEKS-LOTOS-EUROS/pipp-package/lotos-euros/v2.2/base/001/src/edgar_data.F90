!###############################################################################
!
! EDGAR emission inventories.
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!###############################################################################

module EDGAR_Data

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_EDGAR_Data

  public  ::  EDGAR_Data_Init, EDGAR_Data_Done
  public  ::  EDGAR_Data_Get


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'EDGAR_Data'
  

  ! --- local --------------------------------

  ! emission data base
  type T_EDGAR_Data
    ! label assigned to this emission:
    character(len=64)                 ::  label
    ! settings:
    character(len=64)                 ::  rckey
    ! input file template:
    character(len=2024)               ::  fname_template
    ! current year:
    integer                           ::  year
    integer                           ::  month
    ! emisison field:
    real, pointer                     ::  emis(:,:)  ! (nx,ny)
    character(len=32)                 ::  emis_units
    character(len=64)                 ::  tracer
    character(len=64)                 ::  source
  end type T_EDGAR_Data



contains


  ! ===============================================================


  subroutine EDGAR_Data_Init( edg, fname_template, tracer, source, status )

    ! --- in/out ------------------------------

    type(T_EDGAR_Data), intent(out)     ::  edg
    character(len=*), intent(in)        ::  fname_template
    character(len=*), intent(in)        ::  tracer
    character(len=*), intent(in)        ::  source
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/EDGAR_Data_Init'

    ! --- local ------------------------------

    ! --- begin -------------------------------
    
    ! store:
    edg%fname_template = trim(fname_template)
    edg%tracer    = trim(tracer)
    edg%source    = trim(source)

    ! no year/month read yet:
    edg%year  = -999
    edg%month = -999
    
    ! no data stored yet:
    nullify( edg%emis )

    ! ok
    status = 0

  end subroutine EDGAR_Data_Init


  ! ***


  subroutine EDGAR_Data_Done( edg, status )

    ! --- in/out ------------------------------

    type(T_EDGAR_Data), intent(inout)   ::  edg
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/EDGAR_Data_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------
    
    ! clear ?
    if ( associated(edg%emis) ) deallocate( edg%emis )

    ! no stored:
    edg%year  = -999
    edg%month = -999

    !  ok
    status = 0

  end subroutine EDGAR_Data_Done


  ! ***


  subroutine EDGAR_Data_Get( edg, year, month, status )

    use GO              , only : TDate, Get, operator(+), operator(-), operator(/)
    use GO              , only : goReplace, goLoCase
         
    use LE_Grid, only : ugg
    use C3PO   , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors

    ! --- in/out ---------------------------

    type(T_EDGAR_Data), intent(inout)     ::  edg
    integer, intent(in)                   ::  year, month
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/EDGAR_Data_Get'
    
    ! --- local ----------------------------

    
    character(len=1024)                 ::  fname
    character(len=1024)                 ::  vname
    character(len=1024)                 ::  description
    integer                             ::  varid 
    
    character(len=64)                   ::  gridtype_in
    type(T_File_Ugg)                    ::  file_in
    type(T_Grid_Ugg)                    ::  grid_in
    integer                             ::  start_ind(2)
    integer                             ::  count_ind(2)
    real, allocatable                   ::  values_in(:,:)

    
    ! --- begin ----------------------------
    

    ! info ...
    write (gol,'("    read ",a," emission from EDGAR source ", a, " for year ",i4)') trim(edg%tracer), trim(edg%source), year; call goPr

    ! target array:
    if ( .not. associated(edg%emis) ) allocate( edg%emis(ugg%nlon,ugg%nlat) )

    ! copy:
    fname = trim(edg%fname_template)
    ! replace some values:
    call goReplace( fname, '%{year}', '(i4.4)', year, status )
    IF_NOTOK_RETURN(status=1)
    ! replace some values:
    call goReplace( fname, '%{month}', '(i2.2)', month, status )
    IF_NOTOK_RETURN(status=1)
    ! replace some values:
    call goReplace( fname, '%{tracer}', trim(edg%tracer), status )
    IF_NOTOK_RETURN(status=1)
    ! replace some values:
    call goReplace( fname, '%{source}', trim(edg%source), status )
    IF_NOTOK_RETURN(status=1)


    ! variable name:
    write (vname,'("emi_",a)') goLoCase(trim(edg%tracer))
    description = 'var_name='//trim(vname)
    if ( trim(edg%tracer) == 'PM2.5_bio') then
      description = 'var_name=emi_pm2.5'
    else if ( trim(edg%tracer) == 'PM2.5_fossil' ) then
      description = 'var_name=emi_pm2.5'
    end if

    ! open file:
    call file_in%Open( trim(fname), status )
    IF_NOTOK_RETURN(status=1)

    ! variable id:
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! init grid definition
    call file_in%Get_Grid( varid, grid_in, status, &
                          ugg_to=ugg, start_ind=start_ind, count_ind=count_ind )
    IF_NOTOK_RETURN(status=1)

    ! storage
    allocate( values_in(grid_in%nlon, grid_in%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! read:
    call file_in%Get_Var( trim(description), values_in, edg%emis_units, status, &
                           start=start_ind, count=count_ind )
    IF_NOTOK_RETURN(status=1)

    !! testing ...
    !print *, ' emission units = ', trim(edg%emis_units)

    select case ( trim(edg%emis_units) )  
      ! EDGAR units ...
      case ( 'kg m-2 s-1',  'kg(N) m-2 s-1', 'kg(S) m-2 s-1')   
        ! do nothing

      ! MEIC units
      case ( 'Mg/cell/year' )
        ! divide by area to get emission per m2
        call grid_in%AreaOper( values_in, '/', 'm2', status )
        IF_NOTOK_RETURN(status=1)
        ! reset units:
        edg%emis_units= 'Mg/m2/year'

      case ( 'Mg/cell/month' )
        ! divide by area to get emission per m2
        call grid_in%AreaOper( values_in, '/', 'm2', status )
        IF_NOTOK_RETURN(status=1)
        ! reset units:
        edg%emis_units= 'Mg/m2/month'

      ! MEIC units            
      case ( 'Mmol/cell/year' ) 
        ! divide by area to get emission per m2
        call grid_in%AreaOper( values_in, '/', 'm2', status )
        IF_NOTOK_RETURN(status=1) 
        ! reset units:
        edg%emis_units= 'Mmol/m2/year'
      case ( 'Mmol/cell/month' ) 
        ! divide by area to get emission per m2
        call grid_in%AreaOper( values_in, '/', 'm2', status )
        IF_NOTOK_RETURN(status=1) 
        ! reset units:
        edg%emis_units= 'Mmol/m2/month'

      case default
        write (gol,'("Unknown emission unit ")'); call goErr
        TRACEBACK; status=1; return
    end select

    ! ~ regrid
    call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, edg%emis, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( values_in)
    call grid_in%Done(status)
    IF_NOTOK_RETURN(status=1)
    call file_in%Close(status)
    IF_NOTOK_RETURN(status=1)

    ! store current year, month:
    edg%year = year
    edg%month = month         
    
    ! ok
    status = 0

  end subroutine EDGAR_Data_Get


end module EDGAR_Data

