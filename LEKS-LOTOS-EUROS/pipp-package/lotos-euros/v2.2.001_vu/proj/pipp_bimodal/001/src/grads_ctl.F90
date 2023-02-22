!#######################################################################
!
! NAME
!   GrADS_Ctl  -  interface to direct access data with GrADS control file
!
!
! GrADS FILES
!
!   GrADS is a visualisation package for geographical data (grads.iges.org).
!   A data set is read using an ascii control file with suffix '.ctl',
!   which defines the filename(s) of the actual data. 
!   Current implementation writes the data to unformatted direct access files (.dat).
!   Future versions should be able to write control files for other data formats
!   supported by GrADS such as NetCDF.
!
!
! PROCECURES
!
!    type(T_GrADS_Ctl)   ::  ctl
!    integer             ::  status   ! non zero in case of error
!
!    ! open ctl file:
!    call GrADS_Ctl_Init( ctl, '/data/output/', 'LE_ctl_conc-sfc.ctl', status )
!
!    ! comment ...
!    call GrADS_Ctl_Comment( ctl, '', status )
!    call GrADS_Ctl_Comment( ctl, 'GrADS Data Descriptor File', status )
!    call GrADS_Ctl_Comment( ctl, '', status )
!
!    ! data set:
!    call GrADS_Ctl_DSet( ctl, '^LE_ctl_conc-sfc_%y4%m2%d2.nc', status )
!
!    ! title:
!    call GrADS_Ctl_Title( ctl, 'This is output., status )
!
!    ! xdef and ydef:
!    call GrADS_Ctl_XYDef( ctl, 'XDEF', 100, -9.750, 0.50, status )
!    call GrADS_Ctl_XYDef( ctl, 'YDEF', 140, 35.125, 0.25, status )
!
!    ! zdef:
!    call GrADS_Ctl_ZDef( ctl, (/1,2,3,4/), status )
!
!    ! number of time steps, start time, and time step:
!    call GrADS_Ctl_TDef( ctl, 25, (/2006,07,01,00,00/), (/0000,00,00,01,00/), status )
!
!    ! number of variables:
!    call GrADS_Ctl_Vars( ctl, 3, status )
!
!    ! add variable description:
!    call GrADS_Ctl_Var( ctl, 'oro', 0,     'y,x', 'oro [m]', status )
!    call GrADS_Ctl_Var( ctl, 'blh', 0,   't,y,x', 'blh [m]', status )
!    call GrADS_Ctl_Var( ctl, 'o3' , 4, 't,z,y,x', 'O3 [ppm]', status )
!
!    ! end of variables section:
!    call GrADS_Ctl_EndVars( ctl, status )
!
!    ! close ctl file:
!    call GrADS_Ctl_Done( ctl, status )
!
!
! EXAMPLE CTL FILE
!
!    DSET ^LE_ctl_conc-sfc_%y4%m2%d2.nc
!    DTYPE netcdf
!    OPTIONS template
!    TITLE model: LE; expid: ctl
!    UNDEF  -0.3402823E+39
!    XDEF     64 LINEAR    -1.875000    0.250000
!    YDEF     80 LINEAR    47.062500    0.125000
!    ZDEF 1 LEVELS 0
!    TDEF 2 LINEAR 00Z01jul2006  1hr
!    VARS     17
!    o3=>o3                   1 t,z,y,x o3 [mole mole-1]
!    no2=>no2                 1 t,z,y,x no2 [mole mole-1]
!    no=>no                   1 t,z,y,x no [mole mole-1]
!    oh=>oh                   1 t,z,y,x oh [mole mole-1]
!    hno3=>hno3               1 t,z,y,x hno3 [mole mole-1]
!    nh3=>nh3                 1 t,z,y,x nh3 [mole mole-1]
!    so2=>so2                 1 t,z,y,x so2 [mole mole-1]
!    so4a=>so4a               1 t,z,y,x so4a [kg m-3]
!    no3a=>no3a               1 t,z,y,x no3a [kg m-3]
!    nh4a=>nh4a               1 t,z,y,x nh4a [kg m-3]
!    bc=>bc                   1 t,z,y,x bc [kg m-3]
!    pm25=>pm25               1 t,z,y,x pm25 [kg m-3]
!    pm10=>pm10               1 t,z,y,x pm10 [kg m-3]
!    na_f2=>naf               1 t,z,y,x na_f2 [kg m-3]
!    na_c3=>nac               1 t,z,y,x na_c3 [kg m-3]
!    dust_f=>dustf            1 t,z,y,x dust_f [kg m-3]
!    dust_c=>dustc            1 t,z,y,x dust_c [kg m-3]
!    ENDVARS
!
!
!### GO routines #######################################################
!
!
!  ! substitutes for message routines from GO modules
!  
!  ! message line:
!  character(len=256)            ::  gol
!
!  ! display message:
!  subroutine goPr
!    write (*,'(a)') trim(gol)
!  end subroutine goPr
!
!  ! display error message:
!  subroutine goErr
!    write (*,'("ERROR - ",a)') trim(gol)
!  end subroutine goErr
!  
!  ! free file unit:
!  subroutine goGetFU( fu, status )
!    integer, intent(out)    ::  fu
!    integer, intent(out)    ::  status
!    logical                 ::  opened
!    fu = 456
!    do
!      inquire( unit=fu, opened=opened )
!      if ( .not. opened ) exit
!      fu = fu + 1
!    end do
!    status = 0
!  end subroutine goGetFU
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if

#include "le.inc"
!
!#######################################################################


module GrADS_Ctl

  use GO, only : gol, goPr, goErr

  implicit none
  

  ! --- in/out ---------------------------------------

  private

  public  :: T_GrADS_Ctl 
  public  :: GrADS_Ctl_Init, GrADS_Ctl_Done
  public  :: GrADS_Ctl_DSet, GrADS_Ctl_Title, GrADS_Ctl_Comment
  public  :: GrADS_Ctl_XYDef, GrADS_Ctl_ZDef, GrADS_Ctl_TDef
  public  :: GrADS_Ctl_PDef
  public  :: GrADS_Ctl_Vars, GrADS_Ctl_Var, GrADS_Ctl_EndVars
  

  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'GrADS_Ctl'
  
  ! value for undefined data:
  real, parameter   ::  GRADS_UNDEF = -1.0 * huge(1.0)

    
  ! --- types ----------------------------------------

  type T_GrADS_Ctl
    ! full file name:
    character(len=512)          ::  fname
    ! file unit
    integer                     ::  fu
  end type T_GrADS_Ctl


contains


  ! ********************************************************************
  ! ***
  ! *** Control file
  ! ***
  ! ********************************************************************

  
  subroutine GrADS_Ctl_Init( ctl, path, file, status )
  
    use GO, only : goGetFU
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(out)          ::  ctl
    character(len=*), intent(in)            ::  path
    character(len=*), intent(in)            ::  file
    integer, intent(out)                    ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_Init'
    
    ! --- local --------------------------------------

    ! --- begin -------------------------------------
    
    ! file unit not in use yet:
    call goGetFU( ctl%fu, status )
    IF_NOTOK_RETURN(status=1)

    ! filename of ctl file:
    write (ctl%fname,'(a,"/",a)') trim(path), trim(file)
    
    ! open file:
    open( ctl%fu, file=trim(ctl%fname), status='replace', &
                            form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening grads control file :")'); call goErr
      write (gol,'("  file    : ",a)') trim(ctl%fname); call goErr
      write (gol,'("  iostat  : ",i6)') status; call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine GrADS_Ctl_Init
  
  
  ! ***
  
  
  subroutine GrADS_Ctl_Done( ctl, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)          ::  ctl
    integer, intent(out)                      ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_Done'
    
    ! --- local --------------------------------------

    ! --- begin -------------------------------------
    
    ! close:
    close( ctl%fu, iostat=status )
    if (status/=0) then
      write (gol,'("closing grads control file :")'); call goErr
      write (gol,'("  ",a)') trim(ctl%fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine GrADS_Ctl_Done
  
  
  ! ***
  
  
  subroutine GrADS_Ctl_Comment( ctl, line, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    character(len=*), intent(in)          ::  line
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_Comment'

    ! --- local --------------------------------------
    
    ! --- begin -------------------------------------
    
    ! add comment line:
    write (ctl%fu,'("* ",a)') trim(line)

    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_Comment
  
  
  ! ***
  
  
  subroutine GrADS_Ctl_DSet( ctl, dset_template, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    character(len=*), intent(in)          ::  dset_template
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_DSet'

    ! --- local --------------------------------------
    
    ! --- begin -------------------------------------
    
    ! add name of data file template:
    write (ctl%fu,'("DSET ",a)') trim(dset_template)
    
    ! set option to note that DSET is in template format:
    write (ctl%fu,'("OPTIONS template")')
    
    ! data files are in NetCDF:
    write (ctl%fu,'("DTYPE netcdf")')

    ! add value for no-data:
    write (ctl%fu,'("UNDEF ",e15.7)') GRADS_UNDEF
    
    ! add unpack description:
    write (ctl%fu,'("UNPACK scale_factor add_offset")')
    
    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_DSet
  
  
  ! ***
  
  
  subroutine GrADS_Ctl_Title( ctl, title, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    character(len=*), intent(in)          ::  title
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_Title'

    ! --- local --------------------------------------
    
    ! --- begin -------------------------------------
    
    ! add title:
    write (ctl%fu,'("TITLE ",a)') trim(title)

    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_Title
  
  
  ! ***
  
  
  subroutine GrADS_Ctl_PDef( ctl, name, type, nx, ny, lonref, latref, iref, jref, STrueLat, NTrueLat, slon,dx,dy, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  type
    integer, intent(in)                   ::  nx
    integer, intent(in)                   ::  ny
    
    real, intent(in)                      ::  lonref
    real, intent(in)                      ::  latref
    
    integer, intent(in)                   ::  iref
    integer, intent(in)                   ::  jref
    
    real, intent(in)                      ::  STrueLat
    real, intent(in)                      ::  NTrueLat
    real, intent(in)                      ::  slon
    
    real, intent(in)                      ::  dx
    real, intent(in)                      ::  dy

    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_PDef'

    ! --- local --------------------------------------
    
    ! --- begin -------------------------------------
    
    ! x/y definition:
    if (trim(type) == 'LAMBERT') then
      write (ctl%fu,'(a,2i6," LCCR ",2f12.6,2i6,5f12.6)',iostat=status) trim(name), nx,ny,lonref, latref, iref, jref, STrueLat, NTrueLat, slon,dx,dy
      IF_NOTOK_RETURN(status=1)
    else
      write( gol, '("Unknown projection type to define PDEF for projection in GrADS control file: ", a)' ) trim(type) ; call GoErr
      TRACEBACK;status=1;return
    end if

    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_PDef
  
  
  ! ***
  
  
  subroutine GrADS_Ctl_XYDef( ctl, name, n, first, step, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    character(len=*), intent(in)          ::  name
    integer, intent(in)                   ::  n
    real, intent(in)                      ::  first, step
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_XYDef'

    ! --- local --------------------------------------
    
    ! --- begin -------------------------------------
    
    ! x/y definition:
    write (ctl%fu,'(a,i6," LINEAR ",2f12.6)',iostat=status) trim(name), n, first, step
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_XYDef
  
  
  ! ***
  
  
  subroutine GrADS_Ctl_ZDef( ctl, levels, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    integer, intent(in)                   ::  levels(:)
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_ZDef'

    ! --- local --------------------------------------
    
    ! --- begin -------------------------------------
    
    ! z definition:
    write (ctl%fu,'("ZDEF ",i6," LEVELS ",100i6)',iostat=status) size(levels), levels
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_ZDef


  ! ***
  
  
  subroutine GrADS_Ctl_TDEF( ctl, nt, t1, dt, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    integer, intent(in)                   ::  nt
    integer, intent(in)                   ::  t1(5)
    integer, intent(in)                   ::  dt(5)
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_TDEF'
    
    ! 3-character month names:
    character(len=3), parameter  ::  month(12) = &
        (/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/)

    ! --- local --------------------------------------
    
    integer                           ::  nv
    character(len=4)                  ::  tincr

    ! --- begin -------------------------------------
    
    ! set time increment:
    if ( all(dt(2:5) == 0) ) then
      nv = 1
      write (tincr,'(i2,"yr")') dt(1)
    else if ( (dt(1) == 0) .and. all(dt(3:5) == 0) ) then
      nv = 2
      write (tincr,'(i2,"mo")') dt(2)
    else if ( all(dt(1:2) == 0) .and. all(dt(4:5) == 0) ) then
      nv = 3
      write (tincr,'(i2,"dy")') dt(3)
    else if ( all(dt(1:3) == 0) .and. (dt(5) == 0) ) then
      nv = 4
      write (tincr,'(i2,"hr")') dt(4)
    else if ( all(dt(1:4) == 0) ) then
      nv = 5
      write (tincr,'(i2,"mn")') dt(5)
    else
      !write (gol,'("WARNING - unsupported time step : ",5i6)') dt; call goPr
      !tincr = ' 1mn'
      write (gol,'("unsupported time step : ",5i6)') dt; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! write ..
    if ( nv == 2 ) then
      write (ctl%fu,'("TDEF ",i6," LINEAR ",a3,i4.4," ",a)') &
                nt, month(t1(2)), t1(1), tincr
    else if ( nv == 3 ) then
      write (ctl%fu,'("TDEF ",i6," LINEAR ",i2.2,a3,i4.4," ",a)') &
                nt, t1(3), month(t1(2)), t1(1), tincr
    else if ( nv == 4 ) then
      write (ctl%fu,'("TDEF ",i6," LINEAR ",i2.2,"Z",i2.2,a3,i4.4," ",a)') &
                nt, t1(4), t1(3), month(t1(2)), t1(1), tincr
    else if ( nv == 5 ) then
      write (ctl%fu,'("TDEF ",i6," LINEAR ",i2.2,":",i2.2,"Z",i2.2,a3,i4.4," ",a)') &
                nt, t1(4), t1(5), t1(3), month(t1(2)), t1(1), tincr
    else
      !write (gol,'("WARNING - number of time values outside 2-5 : ",i4)') nv; call goPr
      write (gol,'("number of time values outside 2-5 : ",i4)') nv; call goErr
      TRACEBACK; status=1; return
    end if
   
    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_TDEF
  
  
  ! ***
  
  
  subroutine GrADS_Ctl_Vars( ctl, nvar, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    integer, intent(in)                   ::  nvar
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_Vars'

    ! --- local --------------------------------------
    
    ! --- begin -------------------------------------
    
    ! z definition:
    write (ctl%fu,'("VARS ",i6)',iostat=status) nvar
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_Vars


  ! ***
  
  
  subroutine GrADS_Ctl_Var( ctl, name, nz, & 
                            dims, descr, &
                            status &
#ifdef with_labeling  
                            ,label_name &
#endif                            
                           )   
                           
    use GO, only : goReplace
    
    ! --- in/out ---------------------------------
    
    type(T_GrADS_Ctl), intent(inout)          ::  ctl
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  nz
    character(len=*), intent(in)              ::  dims
    character(len=*), intent(in)              ::  descr
    integer, intent(out)                      ::  status
#ifdef with_labeling  
    character(len=*), intent(in), optional    ::  label_name
#endif                            
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_Var'
    
    ! --- local ----------------------------------
    
    character(len=30)   ::  gname
    character(len=50)   ::  s

    ! --- begin ----------------------------------
    
    ! GrADS variable could not contain all characters ...
    gname = trim(name)
    call goReplace( gname, '_', '', status )
    IF_NOTOK_RETURN(status=1)
    
    ! special rename fo cz-canopytop
    if ( trim(gname) == 'czcanopytop' ) then
      call goReplace( gname, 'czcanopytop', 'czctop', status)
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! name part:
    s = repeat(' ',50)
    write (s,'(a,"=>",a)') trim(name), trim(gname)
#ifdef with_labeling
    ! overwrite string s if label is available 
    s = repeat(' ',50)
    write (s,'(a,"=>",a,a)') trim(name), trim(gname), trim(label_name)
#endif
    ! add variable line:
    write (ctl%fu,'(a50,i6," ",a," ",a)',iostat=status) &
              s, nz, trim(dims), trim(descr)
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine GrADS_Ctl_Var


  ! ***
  
  
  subroutine GrADS_Ctl_EndVars( ctl, status )
  
    ! --- in/out -------------------------------------
    
    type(T_GrADS_Ctl), intent(inout)      ::  ctl
    integer, intent(out)                  ::  status
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/GrADS_Ctl_EndVars'

    ! --- local --------------------------------------
    
    ! --- begin -------------------------------------
    
    ! z definition:
    write (ctl%fu,'("ENDVARS")',iostat=status)
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
     
  end subroutine GrADS_Ctl_EndVars



end module GrADS_Ctl



