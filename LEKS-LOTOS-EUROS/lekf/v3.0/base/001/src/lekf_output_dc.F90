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
  
  public  ::  T_LEKF_Output_DC
  
  public  ::  Init, Done
  public  ::  PutOut
  public  ::  ReadDC


  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LEKF_Output_DC'
  
  ! filter states to be written:
  integer, parameter            ::  nfs = 3
  character(len=2), parameter   ::  fsname(nfs) = (/'xb','x ','s '/)
  
  ! filter moments:
  integer, parameter            ::  nfm = 2
  
  ! number of nhist to put out:
  integer, parameter            ::  nhist_out = 1
  
    
  ! --- types ------------------------------
  
  type T_LEKF_Output_DC
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
    ! data sets:
    integer                     ::  varid(nfs,nfm)
  end type T_LEKF_Output_DC
  
  
  ! --- interfaces -------------------------
  
  interface Init
    module procedure kfo_dc_Init
  end interface
  
  interface Done
    module procedure kfo_dc_Done
  end interface

  interface PutOut
    module procedure kfo_dc_PutOut
  end interface

  
contains


  ! ====================================================
  
  
  subroutine kfo_dc_Init( kfo, rcF, rckey, status )

    use GO              , only : TrcFile
    use GO     , only : AnyDate
    use LE_Output_Common, only : Init
    
    ! --- in/out --------------------------------
    
    type(T_LEKF_Output_DC), intent(out)     ::  kfo
    type(TrcFile), intent(in)             ::  rcF
    character(len=*), intent(in)          ::  rckey
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/kfo_dc_Init'
    
    
    ! --- local ---------------------------------
    

    ! --- begin ---------------------------------
    
    ! init common stuff:
    call Init( kfo%com, rcF%fname, rckey, status )
    IF_NOTOK_RETURN(status=1)
    
    ! output time resolution:
    call rcF%Get( 'kf.output.dhour.dc', kfo%dhour, status )
    IF_NOTOK_RETURN(status=1)
    
    ! files not open yet:
    kfo%opened = .false.
    
    ! no time range set yet:
    kfo%tr(1) = AnyDate()
    kfo%tr(2) = AnyDate()
    
    ! ok
    status = 0
    
  end subroutine kfo_dc_Init
  
  
  ! ***
  

  subroutine kfo_dc_Done( kfo, status )
  
    use NetCDF          , only : NF90_Close
    use LE_Output_Common, only : Done
  
    ! --- in/out --------------------------------
    
    type(T_LEKF_Output_DC), intent(inout)   ::  kfo
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/kfo_dc_Done'
    
    ! --- begin ---------------------------------
    
    ! file opened ?
    if ( kfo%opened ) then
      ! close:
      status = NF90_Close( kfo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! reset flag:
      kfo%opened = .true.
    end if
    
    ! done with common stuff ...
    call Done( kfo%com, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine kfo_dc_Done
  
  
  ! ***
  
  
  subroutine kfo_dc_PutOut( kfo, key, t, status )
  
    use GO     , only : TDate, IncrDate, NewDate, AnyDate, Get
    use GO     , only : operator(+), operator(-), operator(<), operator(>)
    use GO     , only : iTotal, MidNight
    use GO     , only : wrtgol, Precisely

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
    use LE_Output_Tools , only : LE_Output_Put_Var_Domains
    
    use LEKF_Dims , only : nx, ny
    use LEKF_Data , only : lekfo_replace
    use LEKF_Noise, only : nnoise, nhist
    use LEKF_State, only : kf_with_xb, kf_with_xm
    use LEKF_State, only : xb, x, sigma
!    use LEKF_State, only : dc_s, dc_e
    use LEKF_State, only : TState
    use LEKF_noise, only : noise_namelen, noise_name
 
    ! --- in/out --------------------------------
    
    type(T_LEKF_Output_DC), intent(inout)   ::  kfo
    character(len=*), intent(in)          ::  key
    type(TDate), intent(in)               ::  t
    integer, intent(out)                  ::  status
    
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/kfo_dc_PutOut'
    
    ! --- local ---------------------------------
    
    type(TDate)           ::  tday
    integer               ::  time6(6)
    integer               ::  time

    integer               ::  cmode

    integer               ::  varid
    type(T_Grid_NcDef)    ::  gncd

    integer               ::  ifs
    character(len=16)     ::  varname
    character(len=512)    ::  descr
    character(len=512)    ::  tunits

    integer               ::  i, j

    ! dimensions lons(nx), lats(ny)
    real, dimension(:), allocatable  ::  lons
    real, dimension(:), allocatable  ::  lats

    real, allocatable     ::  dc(:,:,:,:)

    integer               ::  ifm

    ! --- begin ---------------------------------

    ! storage
    allocate( lons(nx) )
    allocate( lats(ny) )
    
    ! at certain hours only ..
    if ( .not. Precisely(t,kfo%dhour,'hour') ) then
      status=0; return
    end if
    
    ! info ...
    call wrtgol('KF: put out dc for : ',t); call goPr
    
    ! current time not in time range ?
    if ( (t < kfo%tr(1)) .or. (kfo%tr(2) < t) ) then

      ! file opened ?
      if ( kfo%opened ) then
        ! close:
        status = NF90_Close( kfo%ncid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! reset flag:
        kfo%opened = .true.
      end if

      ! day is defined for (00,24]
      tday = t
      if ( MidNight(t) ) tday = tday - IncrDate(day=1)

      ! extract time fields:
      call Get( tday, time6=time6 )

      ! time range for this file is (00,24]
      kfo%tr(1) = NewDate( year=time6(1), month=time6(2), day=time6(3), hour=00 )
      kfo%tr(2) = kfo%tr(1) + IncrDate( day=1 )

      ! reference time for "seconds since ..."
      kfo%t0 = NewDate( time6=(/kfo%tr(1)%year,01,01,00,00,00/) )

      ! new file name:
      write (kfo%fname,'(a,"/",a,"_",a,"_",a,"_",i4.4,2i2.2,".nc")') &
                trim(kfo%com%outdir), &
                trim(kfo%com%model), trim(kfo%com%expid), 'dc', time6(1:3)

      ! set creation mode flag:
      if ( lekfo_replace ) then
        cmode = NF90_CLOBBER       ! overwrite existing files
      else
        cmode = NF90_NOCLOBBER     ! do not overwrite existing files
      end if

      ! create file:
      status = NF90_Create( kfo%fname, cmode, kfo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
      
      ! reset flag:
      kfo%opened = .true.

      ! write global attributes:
      call PutOut_GlobalAttributes( kfo%com, kfo%ncid, status )
      IF_NOTOK_RETURN(status=1)
      
      ! define dimensions:

      status = NF90_Def_Dim( kfo%ncid, 'namelen', noise_namelen, kfo%dimid_namelen )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfo%ncid, 'noise', nnoise, kfo%dimid_noise )
      IF_NF90_NOTOK_RETURN(status=1)
      status = NF90_Def_Dim( kfo%ncid, 'hist', nhist_out, kfo%dimid_hist )
      IF_NF90_NOTOK_RETURN(status=1)

      ! grid dimensions/variables
      call glb_ugg%DefGrid_NetCDF( gncd, kfo%ncid, status, &
                                  dimid_lon=kfo%dimid_lon, dimid_lat=kfo%dimid_lat )
      IF_NOTOK_RETURN(status=1)

      ! time dimensions:
      call LE_Output_Define_Dims_Time(kfo%ncid, kfo%dimid_time, status)
      IF_NOTOK_RETURN(status=1)    
      ! time variables
      call LE_Output_Define_Vars_Time(kfo%ncid, kfo%varid_time, kfo%varid_time_dtg, &
                                      kfo%dimid_time, trim(kfo%com%CF_convention), kfo%t0, status)
      IF_NOTOK_RETURN(status=1)
      
      ! define noise description variable:
      status = NF90_Def_Var( kfo%ncid, 'noise_name', NF90_CHAR, (/kfo%dimid_namelen,kfo%dimid_noise/), varid )
      IF_NF90_NOTOK_RETURN(status=1)
      ! store:
      kfo%varid_name = varid

      ! loop over filter states:
      do ifs = 1, nfs
      
        ! skip ?
        select case ( trim(fsname(ifs)) )
          case ( 'xb'     ) ; if ( .not. kf_with_xb ) cycle
          case ( 'x', 's' ) ; if ( .not. kf_with_xm ) cycle
        end select
        
        ! filter state description:
        select case ( trim(fsname(ifs)) )
          case ( 'xb' ) ; descr = 'model background run'
          case ( 'x'  ) ; descr = 'ensemble mean'
          case ( 's'  ) ; descr = 'ensemble standard deviation'
          case default
            write (*,'("unsupported var name : ",a)') trim(fsname(ifs)); call goErr
            TRACEBACK; status=1; return
        end select
        
        ! loop over filter moments (1=analysis,2=forecast)
        do ifm = 1, nfm
        
          ! skip ?
          select case ( trim(fsname(ifs)) )
            case ( 'xb' ) ; if ( ifm/=1 ) cycle   ! write xb 'analysis' only
          end select
        
          ! variable name:
          write (varname,'(a,"_",a)') 'dc', trim(fsname(ifs))
          
          ! extend for forecast moment:
          if ( ifm == 2 ) varname = trim(varname)//'_f'

          ! define variable:
          status = NF90_Def_Var( kfo%ncid, trim(varname), NF90_REAL, &
              (/kfo%dimid_lon,kfo%dimid_lat,kfo%dimid_noise,kfo%dimid_hist,kfo%dimid_time/), varid )
          IF_NF90_NOTOK_RETURN(status=1)

          ! description:
          status = NF90_Put_Att( kfo%ncid, varid, 'description', descr )
          IF_NF90_NOTOK_RETURN(status=1)

          ! store variable id:
          kfo%varid(ifs,ifm) = varid

        end do  ! filter moments
        
      end do  ! filter states

      ! end defintion mode:

      status = NF90_EndDef( kfo%ncid )
      IF_NF90_NOTOK_RETURN(status=1)
    
      ! write noise description:
      status = NF90_Put_Var( kfo%ncid, kfo%varid_name, noise_name(1:nnoise) )

      ! write grid to netCDF file
      call glb_ugg%PutGrid_NetCDF( gncd, status )
      IF_NOTOK_RETURN(status=1)
      
      ! no records written yet:
      kfo%itrec = 0
      kfo%tprev = AnyDate()
    
    end if
    
    ! next time record ?
    if ( t > kfo%tprev ) kfo%itrec = kfo%itrec + 1

    ! store this time:
    kfo%tprev = t

    ! time since yyyy-01-01 00:00
    time = iTotal( t - kfo%t0, 'sec' )
    ! date up to seconds:
    call Get( t, time6=time6 )
    
    ! write time record:
    call LE_Output_Put_Var_Time(kfo%ncid, kfo%varid_time, kfo%varid_time_dtg, &
                                   time, time6, trim(kfo%com%CF_convention), kfo%itrec, status )
    IF_NOTOK_RETURN(status=1)  

    ! temporary storage:    
    allocate( dc(nx,ny,nnoise,nhist_out) )
    
    ! loop over filter states:
    do ifs = 1, nfs
    
      ! loop over filter moments:
      do ifm = 1, nfm
    
        ! skip ?
        select case ( key )
          ! forecast written to ifm=2
          case ( 'forecast' )
            if ( ifm /= 2 ) cycle
          ! analysis written to ifm=1
          case ( 'analysis' )
            if ( ifm /= 1 ) cycle
          ! error ...
          case default
            write (gol,'("unsupported key : ",a)') trim(key); call goErr
            TRACEBACK; status=1; return
        end select
    
        ! extract simulated concentrations:
        select case ( trim(fsname(ifs)) )
          case ( 'xb' )
            if ( .not. kf_with_xb ) cycle
            if ( ifm /= 1 ) cycle   ! write xb 'analysis' only
            dc(:,:,:,1:nhist_out) = xb%dc(:,:,:,1:nhist_out)
          case ( 'x' )
            if ( .not. kf_with_xm ) cycle
            dc(:,:,:,1:nhist_out) = x%dc(:,:,:,1:nhist_out)
          case ( 's' )
            if ( .not. kf_with_xm ) cycle 
            dc (:,:,:,1:nhist_out)= sigma%dc(:,:,:,1:nhist_out)
          case default
            write (*,'("unsupported filter state name : ",a)') trim(fsname(ifs)); call goErr
            TRACEBACK; status=1; return
        end select

        ! write concentrations:
        status = NF90_Put_Var( kfo%ncid, kfo%varid(ifs,ifm), dc, &
                       start=(/1,1,1,1,kfo%itrec/), count=(/nx,ny,nnoise,nhist_out,1/) )
        IF_NF90_NOTOK_RETURN(status=1)
        
      end do  ! filter momeents
      
    end do  ! filter states
    
    ! clear:
    deallocate( dc )

    ! clear
    deallocate( lons )
    deallocate( lats )

    ! ok
    status = 0
    
  end subroutine kfo_dc_PutOut


  ! ***


  subroutine ReadDC( fdir, model, expid, year, month, day, lag, nt, dc, status )

    use GO     , only : TDate, NewDate, IncrDate, operator(+)
    use LEKF_Dims, only : nx, ny
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

    ! --- begin ---------------------------------

    write (gol,'("read dc from day ",i4,2("-",i2.2)," ; lag ",i2)') year, month, day, lag; call goPr

    ! file 'LE_ner-ax_dc_20030101.nc' contains dc
    !   20030101   ! ihist=1
    !   20021231   ! ihist=2
    !   20021230   ! ihist=3

    ! set day in file name:
    tfile = NewDate( year=year, month=month, day=day ) + IncrDate(day=lag-1)

    ! set file name:
    write (fname,'(a,"/",a,"_",a,"_dc_",i4,2i2.2,".nc")') &
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
    status = NF90_Inq_Varid( ncid, 'dc_x', varid   )
    IF_NF90_NOTOK_RETURN(status=1)

    ! read dc, fill in lag-1
    status = NF90_Get_Var( ncid, varid, dc, &
                            start=(/1,1,1,lag,1/), count=(/nx,ny,nnoise,1,nt/) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! close file:
    status = nf90_close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ReadDC



end module LEKF_Output_DC
