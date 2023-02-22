!###############################################################################
!
! NAME
!
!   KF_Output_DC  -  write noise stuff
!
! HISTORY
!
!   2007 may, Arjo Segers, TNO
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
!###############################################################################


module LEKF_Output_DC

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  use NetCDF, only : NF90_StrError, NF90_NOERR
  
  use LE_Output_Common, only : T_LE_Output_Common
  
  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  T_LEKF_Output_DC_State
  
  public  ::  ReadDC


  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LEKF_Output_DC'
  
  ! number of nhist to put out:
  integer, parameter            ::  nhist_out = 1
  
    
  ! --- types ------------------------------
  
  type T_LEKF_Output_DC_State
    ! short name: "xb", "xf", "xa", ...
    character(len=32)           ::  name
    ! common stuff:
    type(T_LE_Output_Common)    ::  com
    ! file opened ?
    logical                     ::  opened
    ! current time range:
    type(TDate)                 ::  tr(2)
    ! reference time:
    type(TDate)                 ::  t0
    ! time resolution:
    real                        ::  dhour
    ! latest time:
    type(TDate)                 ::  tprev
    ! time record counter:
    integer                     ::  itrec
    ! file name:
    character(len=512)          ::  fname
    ! file handle:
    integer                     ::  ncid
    ! dimension handles:
    integer                     ::  dimid_namelen
    integer                     ::  dimid_lon
    integer                     ::  dimid_lat
    integer                     ::  dimid_noise
    integer                     ::  dimid_hist
    integer                     ::  dimid_time
    integer                     ::  dimid_datelen
    ! dimension variables:
    integer                     ::  varid_name
    integer                     ::  varid_time
    integer                     ::  varid_time_dtg
    integer                     ::  varid_lon
    integer                     ::  varid_lat
    ! data set:
    integer                     ::  varid_dc
  contains
    procedure ::  Init        => LEKF_Output_DC_State_Init
    procedure ::  Done        => LEKF_Output_DC_State_Done
    procedure ::  PutOut      => LEKF_Output_DC_State_PutOut
  end type T_LEKF_Output_DC_State
  
  
  
contains


  ! ====================================================================
  ! ===
  ! === dc output file
  ! ===
  ! ====================================================================
  
  
  subroutine LEKF_Output_DC_State_Init( self, rcF, rckey, name, status )

    use GO              , only : TrcFile
    use GO              , only : AnyDate
    use LE_Output_Common, only : Init
    
    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_DC_State), intent(out)    ::  self
    type(TrcFile), intent(in)                     ::  rcF
    character(len=*), intent(in)                  ::  rckey
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_DC_State_Init'
    
    ! --- local ---------------------------------
    
    ! --- begin ---------------------------------
    
    ! store name:
    self%name = trim(name)
    
    ! init common stuff:
    call Init( self%com, rcF, rckey, status )
    IF_NOTOK_RETURN(status=1)
    
    ! output time resolution:
    call rcF%Get( 'kf.output.dhour.dc', self%dhour, status )
    IF_NOTOK_RETURN(status=1)
    
    ! files not open yet:
    self%opened = .false.
    
    ! no time range set yet:
    self%tr(1) = AnyDate()
    self%tr(2) = AnyDate()
    
    ! ok
    status = 0
    
  end subroutine LEKF_Output_DC_State_Init
  
  
  ! ***
  

  subroutine LEKF_Output_DC_State_Done( self, status )
  
    use NetCDF          , only : NF90_Close
    
    use GO              , only : goc
    use LE_Output_Common, only : Done
  
    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_DC_State), intent(inout)    ::  self
    integer, intent(out)                            ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_DC_State_Done'
    
    ! --- begin ---------------------------------
    
    ! file opened ?
    if ( self%opened ) then
      ! root only:
      if ( goc%root ) then
        ! close:
        status = NF90_Close( self%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
      end if
      ! reset flag:
      self%opened = .false.
    end if
    
    ! done with common stuff ...
    call Done( self%com, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LEKF_Output_DC_State_Done
  
  
  ! ***
  
  
  subroutine LEKF_Output_DC_State_PutOut( self, x, t, status )
  
    use GO     , only : TDate, IncrDate, NewDate, AnyDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>)
    use GO     , only : iTotal, MidNight
    use GO     , only : wrtgol, Precisely
    use GO     , only : goc

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim, NF90_Def_Var, NF90_EndDef
    use NetCDF , only : NF90_Put_Var, NF90_Put_Att
    use NetCDF , only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF , only : NF90_GLOBAL, NF90_UNLIMITED
    use NetCDF , only : NF90_REAL, NF90_INT, NF90_CHAR, NF90_BYTE
    
    use LE_Grid, only : glb_ugg
    use C3PO   , only : T_Grid_NcDef

    use LE_Output_Common, only : PutOut_GlobalAttributes
    use LE_Output_Tools , only : LE_Output_Define_Dims_Time
    use LE_Output_Tools , only : LE_Output_Define_Vars_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_Time
    use LE_Output_Tools , only : LE_Output_Put_Var_4D
    
    use LEKF_Dims , only : nx, ny
    use LEKF_Data , only : lekfo_replace
    use LEKF_Noise, only : nnoise, nhist
    use LEKF_State, only : TState
    use LEKF_noise, only : noise_namelen, noise_name
 
    ! --- in/out --------------------------------
    
    class(T_LEKF_Output_DC_State), intent(inout)    ::  self
    type(TState), intent(in)                        ::  x
    type(TDate), intent(in)                         ::  t
    integer, intent(out)                            ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LEKF_Output_DC_State_PutOut'
    
    ! --- local ---------------------------------
    
    type(TDate)           ::  tday
    integer               ::  time6(6)
    integer               ::  time

    integer               ::  cmode

    type(T_Grid_NcDef)    ::  gncd

    integer               ::  i, j

    integer               ::  inoise

    ! --- begin ---------------------------------

    ! at certain hours only ..
    if ( .not. Precisely(t,self%dhour,'hour') ) then
      status=0; return
    end if
    
    ! info ...
    call wrtgol('KF: put out dc for : ',t); call goPr
    
    ! current time not in time range ?
    if ( (t < self%tr(1)) .or. (self%tr(2) < t) ) then

      ! file opened ?
      if ( self%opened ) then
        ! root only:
        if ( goc%root ) then
          ! close:
          status = NF90_Close( self%ncid )
          IF_NF90_NOTOK_RETURN(status=1)
        end if
        ! reset flag:
        self%opened = .true.
      end if ! file opened

      ! day is defined for (00,24]
      tday = t
      if ( MidNight(t) ) tday = tday - IncrDate(day=1)

      ! extract time fields:
      call Get( tday, time6=time6 )

      ! time range for this file is (00,24]
      self%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
      self%tr(2) = self%tr(1) + IncrDate( day=1 )

      ! reference time for "seconds since ..."
      self%t0 = NewDate( time6=(/self%tr(1)%year,01,01,00,00,00/) )

      ! root only:
      if ( goc%root ) then

        ! new file name:
        write (self%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,"_",a,".nc")') &
                  trim(self%com%outdir), &
                  trim(self%com%model), trim(self%com%expid), 'dc', time6(1:3), &
                  trim(self%name)

        ! set creation mode flag:
        if ( lekfo_replace ) then
          cmode = NF90_CLOBBER       ! overwrite existing files
        else
          cmode = NF90_NOCLOBBER     ! do not overwrite existing files
        end if

        ! create file:
        status = NF90_Create( self%fname, cmode, self%ncid )
        IF_NF90_NOTOK_RETURN(status=1)

        ! write global attributes:
        call PutOut_GlobalAttributes( self%com, self%ncid, status )
        IF_NOTOK_RETURN(status=1)

        ! define dimensions:

        status = NF90_Def_Dim( self%ncid, 'namelen', noise_namelen, self%dimid_namelen )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( self%ncid, 'noise', nnoise, self%dimid_noise )
        IF_NF90_NOTOK_RETURN(status=1)
        status = NF90_Def_Dim( self%ncid, 'hist', nhist_out, self%dimid_hist )
        IF_NF90_NOTOK_RETURN(status=1)

        ! grid dimensions/variables
        call glb_ugg%DefGrid_NetCDF( gncd, self%ncid, status, &
                                    dimid_lon=self%dimid_lon, dimid_lat=self%dimid_lat )
        IF_NOTOK_RETURN(status=1)

        ! time dimensions:
        call LE_Output_Define_Dims_Time(self%ncid, self%dimid_time, status)
        IF_NOTOK_RETURN(status=1)    
        ! time variables
        call LE_Output_Define_Vars_Time(self%ncid, self%varid_time, self%varid_time_dtg, &
                                        self%dimid_time, trim(self%com%CF_convention), self%t0, status)
        IF_NOTOK_RETURN(status=1)

        ! define noise description variable:
        status = NF90_Def_Var( self%ncid, 'noise_name', NF90_CHAR, (/self%dimid_namelen,self%dimid_noise/), self%varid_name )
        IF_NF90_NOTOK_RETURN(status=1)

        ! define variable:
        status = NF90_Def_Var( self%ncid, 'dc', NF90_REAL, &
            (/self%dimid_lon,self%dimid_lat,self%dimid_noise,self%dimid_hist,self%dimid_time/), self%varid_dc )
        IF_NF90_NOTOK_RETURN(status=1)

        ! end defintion mode:
        status = NF90_EndDef( self%ncid )
        IF_NF90_NOTOK_RETURN(status=1)

        ! write noise description:
        status = NF90_Put_Var( self%ncid, self%varid_name, noise_name(1:nnoise) )

        ! write grid to netCDF file
        call glb_ugg%PutGrid_NetCDF( gncd, status )
        IF_NOTOK_RETURN(status=1)

      end if ! root
      
      ! reset flag:
      self%opened = .true.

      ! no records written yet:
      self%itrec = 0
      self%tprev = AnyDate()
    
    end if  ! new file
    
    ! next time record ?
    if ( t > self%tprev ) self%itrec = self%itrec + 1

    ! store this time:
    self%tprev = t
    
    ! root only:
    if ( goc%root ) then

      ! time since yyyy-01-01 00:00
      time = iTotal( t - self%t0, 'sec' )
      ! date up to seconds:
      call Get( t, time6=time6 )

      ! write time record:
      call LE_Output_Put_Var_Time(self%ncid, self%varid_time, self%varid_time_dtg, &
                                   time, time6, trim(self%com%CF_convention), self%itrec, status )
      IF_NOTOK_RETURN(status=1)

    end if ! root
    
    !! testing ..
    !do j = 1, ny
    !  write (gol,*) 'ccc1 ', j, x%dc(1,j,1,1); call goPr
    !end do

    ! write 4D field with time index:
    call LE_Output_Put_Var_4D( self%ncid, self%varid_dc, self%itrec, x%dc(:,:,:,1:nhist_out), status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LEKF_Output_DC_State_PutOut



  ! ====================================================================
  ! ===
  ! === read dc field
  ! ===
  ! ====================================================================
  
  subroutine ReadDC( fdir, model, expid, year, month, day, lag, nt, dc, status )

    use GO        , only : TDate, NewDate, IncrDate, operator(+)
    use LE_Grid   , only : dom
    use LEKF_Dims , only : nx, ny
    use LEKF_Noise, only : nnoise, nhist
    use NetCDF

    ! --- in/out -------------------------------

    character(len=*), intent(in)          ::  fdir
    character(len=32), intent(in)         ::  model
    character(len=32), intent(in)         ::  expid
    integer, intent(in)                   ::  year, month, day
    integer, intent(in)                   ::  lag
    integer, intent(in)                   ::  nt
    real, intent(out)                     ::  dc(nx,ny,nnoise,nt)
    integer, intent(out)                  ::  status

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/ReadDC'

    ! --- local ---------------------------------

    character(len=512)      ::  fname
    type(TDate)             ::  tfile
    logical                 ::  exist
    integer                 ::  ncid, varid
    integer                 ::  off(2), shp(2)

    ! --- begin ---------------------------------

    write (gol,'("read dc from day ",i4,2("-",i2.2)," ; lag ",i2)') year, month, day, lag; call goPr

    ! file 'LE_ner-ax_dc_20030101.nc' contains dc
    !   20030101   ! ihist=1
    !   20021231   ! ihist=2
    !   20021230   ! ihist=3

    ! set day in file name:
    tfile = NewDate( year=year, month=month, day=day ) + IncrDate(day=lag-1)

    ! set file name:
    write (fname,'(a,"/",a,"_",a,"_dc_",i4,2i2.2,"_xa.nc")') &
             trim(fdir), trim(model), trim(expid), &
             tfile%year, tfile%month, tfile%day

    ! check if file exist ...
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("WARNING - dc file not found (keep current dc) :")'); call goPr
      write (gol,'("WARNING -   ",a)') trim(fname); call goPr
      status=-1; return
    else
      write (gol,'("  dc file : ",a)') trim(fname); call goPr
      write (gol,'("  ihist   : ",i2)') lag; call goPr
    end if

    ! open file:
    status = NF90_Open( trim(fname), NF90_NOWRITE, ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! get variable id's:
    status = NF90_Inq_Varid( ncid, 'dc', varid   )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! offset in global index space, shape:
    call dom%Get( status, off=off, shp=shp )
    IF_NOTOK_RETURN(status=1)

    ! read dc, fill in lag-1
    status = NF90_Get_Var( ncid, varid, dc, &
                            start=(/off(1)+1,off(2)+1,     1,lag, 1/), &
                            count=(/shp(1)  ,shp(2)  ,nnoise,  1,nt/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ReadDC



end module LEKF_Output_DC
