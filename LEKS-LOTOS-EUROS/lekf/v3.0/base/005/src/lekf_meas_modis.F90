!###############################################################################
!
! NAME
!
!   KF_Meas_MODIS  -  interface to MODIS data
!
!
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


module LEKF_Meas_MODIS

  use GO, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out ----------------------------
  
  private

  public  ::  InitMeas_MODIS, DoneMeas_MODIS
  public  ::  GetMeas_MODIS, CalcMeas_MODIS, MeasUpdate_MODIS
  
  
  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LEKF_Meas_MODIS'
  
  

contains


  ! ======================================================================
  ! ===
  ! === input, update
  ! ===
  ! ======================================================================
  
  
  ! opens the measurement file (if present)
  ! end sets some parameters

  subroutine InitMeas_MODIS( status )

    use GO, only : TDate
    
    ! --- in/out --------------------------------
    
    integer, intent(out)            ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/initmeas_MODIS'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    write (gol,'("KF:     setup MODIS measurements ...")'); call goPr
    
    ! ok
    status = 0

  end subroutine InitMeas_MODIS


  ! ***
  

  subroutine DoneMeas_MODIS( status )

    ! --- in/out --------------------------------
    
    integer, intent(out)            ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/DoneMeas_MODIS'
    
    ! --- local ------------------------------
    
    ! --- begin ------------------------------
    
    ! ok
    status = 0

  end subroutine DoneMeas_MODIS
  
  
  ! ***
  
 
  subroutine GetMeas_MODIS( status )

    ! --- in/out ----------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Getmeas_MODIS'
    
    ! --- local -----------------------------

    ! --- begin -----------------------------
    
    ! ok
    status = 0

  end subroutine GetMeas_MODIS


  ! ***
  
  
  ! simulate MODIS AOD given concentrations
  
  subroutine CalcMeas_MODIS( c, modis_aod, n, status )
  
    use Dims           , only : nx, ny, nz, nspec
    use LE_Grid        , only : lli
    use LE_AOD         , only : LE_AOD_calc
    use LE_Output_MODIS, only : T_Simulation_MODIS, Simulate
    use LEKF_Dims      , only : modis_maxpix
    use LEKF_Data      , only : leo
  
    ! --- in/out -----------------------------

    real, intent(in)          ::  c(nx,ny,nz,nspec)   ! volume ppb
    real, intent(out)         ::  modis_aod(modis_maxpix)  ! aod
    integer, intent(out)      ::  n
    integer, intent(out)      ::  status
    
    ! --- const --------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/CalcMeas_MODIS'
    
    ! --- local --------------------------------
    
    real                        ::  AOD(nx,ny)
    type(T_Simulation_MODIS)    ::  sim_modis
    
    ! --- local --------------------------------
    
    write (gol,'("KF:     update modis simulation ...")'); call goPr
    
    ! simulate AOD field:
    call LE_AOD_calc( c, AOD, status )
    IF_NOTOK_RETURN(status=1)
    
    ! simulate modis fields:
    call Simulate( leo%data_modis, lli, AOD, sim_modis, status )
    IF_NOTOK_RETURN(status=1)

    ! extract:
    n = sim_modis%npix
    modis_aod = 0.0
    modis_aod(1:n) = sim_modis%aod(1:n)

    ! ok
    status = 0
    
  end subroutine CalcMeas_MODIS


  ! ***
  
  
  !
  ! analyse all observations in (t1,t2]
  !

  subroutine MeasUpdate_MODIS( nmodes, status )

    use LE_Output_MODIS, only : bpos_nodata, bpos_validation, bpos_screened, bpos_analysed
    use LEKF_Dims      , only : maxmodes
    use LEKF_State     , only : substate_modis
    use LEKF_Data      , only : leo
    use LEKF_Meas      , only : Meas_Update_Scalar_Cell

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  nmodes     ! actual number of modes
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MeasUpdate_MODIS'

    ! --- local ----------------------------------
    
    integer             ::  npix, ipix
    real                ::  y, r
    integer             ::  inds(0:7)
    integer             ::  ix, iy
    real                ::  alfa
    
    ! --- begin ----------------------------------
    
    ! validation only ?
    if ( .not. leo%data_modis%analyse ) then
      ! info ...
      write (gol,'("LEKF:   do not analyse MODIS, validation only ...")'); call goPr
      ! loop over pixels:
      do ipix = 1, leo%data_modis%npix
        ! set validation flag:
        leo%data_modis%pix_assim_status(ipix) = &
              IBSet( leo%data_modis%pix_assim_status(ipix), bpos_validation )
      end do
      ! ok
      status = 0; return
    end if
    
    ! info ...
    write (gol,'("LEKF:   measurement update MODIS ...")'); call goPr

    ! check ...
    if ( nmodes > maxmodes ) then
      write (gol,'("LEKF:     nmodes > maxmodes : ",2i6)') nmodes, maxmodes; call goErr
      TRACEBACK; status=1; return
    end if

    ! get current number of observed values:
    npix = leo%data_modis%npix
    
    ! no measurements ? then leave
    if ( npix < 1 ) then
      write (gol,'("LEKF:     no measurements available ...")'); call goPr
      status=0; return
    end if
    
    write (gol,'("LEKF:     analyse ",i6," measurements")') npix; call goPr
    
    ! loop over measurements:
    do ipix = 1, npix
    
      ! already flagged ? then skip:
      if ( leo%data_modis%pix_assim_status(ipix) > 0 ) cycle
    
      ! extract measurment:
      y = leo%data_modis%pix_Optical_Depth_Land_And_Ocean(ipix)  ! m
      ! assumed error std.dev.:
      r = leo%data_modis%pix_Optical_Depth_Land_And_Ocean_sigma(ipix)  ! m
      
      ! screening factor:
      alfa = leo%data_modis%screening_factor
    
      ! nan ? then skip
      if ( y < 0 ) then
        ! set no-data bit:
        leo%data_modis%pix_assim_status(ipix) = &
              IBSet( leo%data_modis%pix_assim_status(ipix), bpos_nodata )
        ! next:
        cycle
      end if

      !! info ...
      !iana = iana + 1
      !if ( modulo(iana,1000) == 0 ) then
      !  write (gol,'("LEKF:     analyse ",2i6)') iana, nana; call goPr
      !end if
      
      ! define how to extract value from state:
      inds(0:1) = (/ substate_modis, ipix /)
      ! corresponding grid cell:
      ix = leo%data_modis%pix_ix(ipix)
      iy = leo%data_modis%pix_iy(ipix)

      ! analyse measurement (only modes are adjusted, 
      ! so don't forget to update mean later on):
      call meas_update_scalar_cell( y, r**2, inds, ix, iy, alfa, nmodes, status )
      ! rejected or error ?
      if ( status == -1 ) then
        ! set screened flag:
        leo%data_modis%pix_assim_status(ipix) = &
              IBSet( leo%data_modis%pix_assim_status(ipix), bpos_screened )
      else if ( status == 0 ) then
        ! analysed ok; set flag:
        leo%data_modis%pix_assim_status(ipix) = &
              IBSet( leo%data_modis%pix_assim_status(ipix), bpos_analysed )
      else
        ! some error ...
        TRACEBACK; status=1; return
      end if
      
            !! testing ...
            !print *,  'break after first pixel'
            !if (ipix==5) exit
      
    end do  ! measurements

    ! ok
    status = 0

  end subroutine MeasUpdate_MODIS


end module LEKF_Meas_MODIS



