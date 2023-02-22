!###############################################################################
!
! Access to MEGAN input files.
!
!### macro's ##################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_MEGAN_Input

  use GO, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_LE_Emis_MEGAN_Input

  public  ::  LE_Emis_MEGAN_Input_Init, LE_Emis_MEGAN_Input_Done
  public  ::  LE_Emis_MEGAN_Input_Setup


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_MEGAN_Input'
  

  ! --- local --------------------------------

  ! emission data base
  type T_LE_Emis_MEGAN_Input
    ! input file template:
    character(len=1024)               ::  fname_template
    character(len=1024)               ::  vname_template
    ! latest file and variable actually read:
    character(len=1024)               ::  fname
    character(len=1024)               ::  vname
    ! emisison field:
    real, allocatable                 ::  data(:,:)  ! (nx,ny)
    character(len=32)                 ::  units
  end type T_LE_Emis_MEGAN_Input



contains


  ! ===============================================================


  subroutine LE_Emis_MEGAN_Input_Init( meg, fname_template, vname_template, status )

    ! --- in/out ------------------------------

    type(T_LE_Emis_MEGAN_Input), intent(out) ::  meg
    character(len=*), intent(in)             ::  fname_template
    character(len=*), intent(in)             ::  vname_template
    integer, intent(out)                     ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_MEGAN_Input_Init'

    ! --- local ------------------------------

    ! --- begin -------------------------------
    
    ! store:
    meg%fname_template = trim(fname_template)
    meg%vname_template = trim(vname_template)

    ! no data stored yet:
    meg%fname = '-'
    meg%vname = '-'

    ! ok
    status = 0

  end subroutine LE_Emis_MEGAN_Input_Init


  ! ***


  subroutine LE_Emis_MEGAN_Input_Done( meg, status )

    ! --- in/out ------------------------------

    type(T_LE_Emis_MEGAN_Input), intent(inout)  ::  meg
    integer, intent(out)                        ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_MEGAN_Input_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------
    
    ! clear ?
    if ( allocated(meg%data) ) then
      deallocate( meg%data, stat=status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! no stored:
    meg%fname = '-'
    meg%vname = '-'

    !  ok
    status = 0

  end subroutine LE_Emis_MEGAN_Input_Done


  ! ***


  subroutine LE_Emis_MEGAN_Input_Setup( meg, status, &
                                   year, month, component, comp, mt, &
                                   units )

    use GO            , only : goReplace
    use LE_Grid       , only : ugg
    use C3PO          , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors
    
    ! --- in/out ---------------------------

    type(T_LE_Emis_MEGAN_Input), intent(inout)  ::  meg
    integer, intent(out)                        ::  status

    integer, intent(in), optional             ::  year
    integer, intent(in), optional             ::  month
    character(len=*), intent(in), optional    ::  component   ! component long name "limone" etc
    character(len=*), intent(in), optional    ::  comp        ! component short name "lim" etc
    character(len=*), intent(in), optional    ::  mt   ! '' or 'MT' for directory part

    character(len=*), intent(in), optional    ::  units   ! target units

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_MEGAN_Input_Setup'
    
    ! month names as used in variable names:
    character(len=3), parameter  ::  monthname(12) = &
         (/ 'Jan' , 'Feb' , 'Mar' , 'Apr' , 'May' , 'Jun' , &
            'Jul' , 'Aug' , 'Sep' , 'Oct' , 'Nov' , 'Dec'   /)
    
    ! --- local ----------------------------

    character(len=1024)   ::  fname
    character(len=1024)   ::  vname
    character(len=1024)   ::  description
    integer               ::  varid 

    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    real, allocatable     ::  values_in(:,:)
    character(len=64)     ::  units_in
    real                  ::  missing_value
    
    integer               ::  start_ind(2)
    integer               ::  count_ind(2)

    ! --- begin ----------------------------
      
    ! start with file and variable names from templates:
    fname = trim(meg%fname_template)
    vname = trim(meg%vname_template)
    ! replace year part ?
    if ( present(year) ) then
      call goReplace( fname, '%{year}', '(i4.4)', year, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( vname, '%{year}', '(i4.4)', year, status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! replace month part ?
    if ( present(month) ) then
      call goReplace( fname, '%{month}', '(i2.2)', month, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( vname, '%{month}', '(i2.2)', month, status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( fname, '%{monthname}', monthname(month), status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( vname, '%{monthname}', monthname(month), status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! component part in filename ?
    if ( present(component) ) then
      call goReplace( fname, '%{component}', trim(component), status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( vname, '%{component}', trim(component), status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! component short name part in filename ?
    if ( present(component) ) then
      call goReplace( fname, '%{comp}', trim(comp), status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( vname, '%{comp}', trim(comp), status )
      IF_NOTOK_RETURN(status=1)
    end if
    ! monoterp dir?
    if ( present(mt) ) then
      call goReplace( fname, '%{mt}', trim(mt), status )
      IF_NOTOK_RETURN(status=1)
      call goReplace( vname, '%{mt}', trim(mt), status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! new file and/or variable ?
    if ( (trim(fname) /= trim(meg%fname)) .or. &
         (trim(vname) /= trim(meg%vname))         ) then
    
      ! store current file and variable name:
      meg%fname = trim(fname)
      meg%vname = trim(vname)
    
      ! info ...
      write (gol,'("      read MEGAN field `",a,"` ...")') trim(vname); call goPr

      description = 'var_name='//trim(vname)
      ! open file:
      call file_in%Open( trim(fname), status )
      IF_NOTOK_RETURN(status=1)
      
      ! variable id:
      call file_in%Inq_VarID( trim(description), varid, status )
      IF_NOTOK_RETURN(status=1)
      ! init grid definition, limit to current domain:
      call file_in%Get_Grid( varid, grid_in, status, &
                                  ugg_to=ugg, start_ind=start_ind, count_ind=count_ind )
      IF_NOTOK_RETURN(status=1)
      
      ! storage
      allocate( values_in(grid_in%nlon, grid_in%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! read:
      call file_in%Get_Var( trim(description), values_in, units_in, status, &
                             start=start_ind, count=count_ind, &
                             missing_value=missing_value )
      IF_NOTOK_RETURN(status=1)
      ! replace missing values by zero:
      where( values_in == missing_value )
        values_in = 0.0
      endwhere      
      
      ! ~ unit conversion
      
      !~ According to the "MEGAN Community Data Portal (CDP) User Guide",
      !   "II. INPUT FOLDER", "LAI" :
      !      `... The unist are m2 per 1000 m2. Deivde by 1000 to get units of m2 m-2
      !      which is required for the [EFMAP_data.csv input file] for the MEGAN code.'
      !  The "units" attribute says "Degree" however ...
      if ( (index(vname,'LAI') > 0) .and. (index(vname,'(m2_per_m2)') > 0) &
             .and. (units_in == 'Degree') ) then
        ! convert from "m2 per 1000 m2" to "m2 m-2" :
        values_in = values_in / 1000.0
        ! reset units:
        meg%units = 'm2 m-2'
      !~ PFT fractions seem to be in '%' while units say 'Degree' :
      else if ( (index(vname,'cover_fraction') > 0) .and. (index(vname,'(m2_per_m2)') > 0) &
                   .and. (units_in == 'Degree') ) then
        ! reset units:
        meg%units = '%'
      !~ strange units for EF fields ..
      else if ( (index(vname,'(microgram_per_m2_per_hr)') > 0) &
                 .and. (units_in == 'Degree') ) then
        ! reset units:
        meg%units = 'ug m-2 h-1'
      !~ ensure that all variable are checked for need of unit conversion ..
      else
        write (gol,'("no unit conversion for variable `",a,"` supported yet")') trim(vname); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! target array:
      if ( .not. allocated(meg%data ) ) then
        allocate( meg%data(ugg%nlon,ugg%nlat), stat=status )
        IF_NOTOK_RETURN(status=1)
      end if
           
      ! ~ regrid
      call Grid_Convertors%Ugg_AreaAver( grid_in, values_in, ugg, meg%data, status )
      IF_NOTOK_RETURN(status=1)
      
      ! clear:
      deallocate( values_in )
      IF_NOTOK_RETURN(status=1)
      ! done:
      call grid_in%Done(status)
      IF_NOTOK_RETURN(status=1)
      ! close:
      call file_in%Close( status )
      IF_NOTOK_RETURN(status=1)

    end if  ! new year

    ! ok
    status = 0

  end subroutine LE_Emis_MEGAN_Input_Setup



end module LE_Emis_MEGAN_Input
