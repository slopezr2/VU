!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Data_Calc

  use GO  , only : gol, goPr, goErr
  use GO  , only : TDate

  implicit none


  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Data_Calc_Init, LE_Data_Calc_Done
  public  ::  Get_Cloud_Profiles
  public  ::  SoilMoistureIndex
  public  ::  SoilMoistureIndex_wrf


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Data_Calc'


  ! --- local ------------------------------------
  
  ! time that soiltype files swith definition
  type(TDate)     ::  t_soiltype_switch


contains


  ! ====================================================================
  ! ===
  ! === module
  ! ===
  ! ====================================================================
  

  subroutine LE_Data_Calc_Init( status )
  
    use GO, only : NewDate
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Data_Calc_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------
    
    ! switch of soil types in ECMWF archive (2007-11-06)  
    t_soiltype_switch = NewDate( year=2007, month = 11, day = 6)
    
    ! ok
    status = 0
    
  end subroutine LE_Data_Calc_Init
  
  
  ! ***
  
  
  subroutine LE_Data_Calc_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Calc_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
    
  end subroutine LE_Data_Calc_Done
  
  
  ! ***
  
  !
  ! Compute cloud properties at LOTOS-EUROS layers given
  ! cloud profiles from meteo model.
  !
  ! Levels in profiles are ordered surface to top.
  !
  ! Copy from original version in "le_meteo_ecmwf.F90".
  !
  
  subroutine Get_Cloud_Profiles( mm_n, mm_ph_in, mm_tcc, mm_icc, mm_clwc, mm_ciwc, &
                             le_n, le_ph, &
                             le_icc, le_bcc, le_occ, le_clwc, le_ciwc, le_pri, status )
                             !debug )
  
    use JAQL, only : CloudCoverOverhead
    use Num , only : IntervalQuad_Const, IntervalQuad_Lin
    use Binas,only : Grav
  
    ! --- in/out ---------------------------------
    
    integer, intent(in)   ::  mm_n             ! number of meteo model layers
    real, intent(in)      ::  mm_ph_in(0:mm_n) ! meteo model half-level pressure [Pa]
    real, intent(in)      ::  mm_tcc           ! meteo model total cloud cover
    real, intent(in)      ::  mm_icc(1:mm_n)   ! meteo model (in)cloud cover fractions [0-1]
    real, intent(in)      ::  mm_clwc(1:mm_n)  ! meteo model cloud liquid water content [(kg water)/(kg air)]
    real, intent(in)      ::  mm_ciwc(1:mm_n)  ! meteo model cloud ice water content [(kg water)/(kg air)]
    integer, intent(in)   ::  le_n             ! number of LE layers
    real, intent(in)      ::  le_ph(0:le_n)    ! LE half-level pressure [Pa]
    real, intent(out)     ::  le_icc(1:le_n)   ! LE in-cloud cover [0-1]
    real, intent(out)     ::  le_bcc(1:le_n)   ! LE below cloud cover [0-1]
    real, intent(out)     ::  le_occ(1:le_n)   ! LE overhead cloud cover [0-1]
    real, intent(out)     ::  le_clwc(1:le_n)  ! LE cloud liquid water content [(kg water)/(kg air)]
    real, intent(out)     ::  le_ciwc(1:le_n)  ! LE cloud ice water content [(kg water)/(kg air)]
    real, intent(out)     ::  le_pri(1:le_n)   ! LE Potential Rain intensity profile [1]
    integer, intent(out)  ::  status
    
    !logical, intent(in), optional   ::  debug
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/Get_Cloud_Profiles'
    
    ! --- local ----------------------------------

    real      ::  mm_ph(0:mm_n)     ! meteo model half-level pressure [Pa]
    real      ::  mm_cco(mm_n)      ! cloud cover overhead [0-1]
    real      ::  mm_bcc(mm_n)      ! below cloud cover [0-1]
    integer   ::  k
    integer   ::  ilast
    
    real      ::  mm_pri(0:mm_n)
    real      ::  mm_clwc_overhead(0:mm_n)
    real      ::  mm_ciwc_overhead(0:mm_n)
    real      ::  dm
    real      ::  tclwc             ! total liquid water column [kg/m2]
    
    !logical   ::  dbg
    
    ! --- begin ----------------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! Compute cloud-cover-overhead fractions, e.g. fraction of cell with some clouds above;
    ! for the distirbution of the clouds ECMWF maximum/random overlap approach is used:
    call CloudCoverOverhead( mm_n, mm_icc, mm_cco, status )
    IF_NOTOK_RETURN(status=1)
    
    ! correct for the case that part of the model top is missing;
    ! the 2D field of "total cloud cover" is then higher 
    ! than cco at the bottom:
    if ( mm_tcc > mm_cco(1) ) then
      ! add the difference:
      mm_cco = mm_cco + (mm_tcc - mm_cco(1))
    end if
    ! clip values outside [0,1] that might be present due to rounding errors:
    mm_cco = min( max( 0.0, mm_cco ), 1.0 )
    
    ! Compute below-cloud-coverage: part of the cell that has no cloud,
    ! but is below a cloud in a higher layer.
    ! Here set to the difference between the cloud-cover-overhead
    ! and the in-cloud-coverage; clip to [0,1] to remove rounding errors:
    mm_bcc = min( max( 0.0, mm_cco - mm_icc ), 1.0 )
    
    ! adhoc fix: surface of meteo profile should be equal or lower than model surface ...
    mm_ph = mm_ph_in
    mm_ph(0) = max( mm_ph(0), le_ph(0) )
    
    ! total water column;
    ! use that air mass column density 'dm' is : dp/g [kg/m2] ,
    ! thus water column is integral over clwc:
    !     tclwc = int clwc dm    [(kg water)/m2]
    ! accumulated water column..  
    mm_clwc_overhead(mm_n) = 0.0
    mm_ciwc_overhead(mm_n) = 0.0
    do k = mm_n-1, 0, -1
      ! layer air mass column density [kg air/m2]
      dm = ( mm_ph(k) - mm_ph(k+1) ) / Grav
      ! total accumulated water column [kg water/m2]
      !                         (kg water)/m2       (kg water)/(kg air) * (kg air)/m2
      mm_clwc_overhead(k) = mm_clwc_overhead(k+1) +      mm_clwc(k+1)   *      dm
      mm_ciwc_overhead(k) = mm_ciwc_overhead(k+1) +      mm_ciwc(k+1)   *      dm
    end do
    ! store total column
    tclwc = mm_clwc_overhead(0)+mm_ciwc_overhead(0)
    if ( tclwc > 0.0 ) then
      ! compute accumulated potential rain intensity
      !          (kg water)/m2    (kg water)/m2
      mm_pri = (mm_clwc_overhead + mm_ciwc_overhead(0) )/      tclwc     ! 1
    else
      ! no water in total column
      mm_pri = 0.0  ! 1
    end if
    
    !! info ...
    !if ( dbg ) then
    !  print *, 'ppp0 mm_n ', mm_n
    !  print *, ''
    !  print *, 'ppp1 *** meteo model layers:'
    !  print *, 'ppp1 lv      ptop      pbot       icc       cco       bcc        clwc        ciwc    mm_pri'
    !  do k = mm_n, 1, -1
    !    write (*,'(" ppp1 ",i2,2f10.2,3f10.6,3E12.4)') &
    !              k, mm_ph(k), mm_ph(k-1), mm_icc(k), mm_cco(k), mm_bcc(k), mm_clwc(k), mm_ciwc(k),mm_pri(k)
    !  end do
    !  k=0
    !  write (*,'(" ppp1 ",i2,f10.2,"                                                    ",E12.4)') &
    !            k, mm_ph(k),mm_pri(k)
    !  print *, ''
    !  print *, 'ppp1 mm tclwc = ', tclwc, ' kg/m2'
    !  print *, ''
    !end if
    
    ! average values over LE layers:
    ilast = 1
    do k = 1, le_n

      !~~ in-cloud cover
      ! integral over pressure ax; negative pressure to have increasing ax:
      call IntervalQuad_Const( -1.0*mm_ph, mm_icc, &
                               -1.0*le_ph(k-1), -1.0*le_ph(k), &
                               le_icc(k), ilast, status )
      IF_NOTOK_RETURN(status=1)
      ! devide for average:
      le_icc(k) = le_icc(k) / abs(le_ph(k-1)-le_ph(k))

      !~~ below-cloud cover
      ! integral over pressure ax; negative pressure to have increasing ax:
      call IntervalQuad_Const( -1.0*mm_ph, mm_bcc, &
                               -1.0*le_ph(k-1), -1.0*le_ph(k), &
                               le_bcc(k), ilast, status )
      IF_NOTOK_RETURN(status=1)
      ! devide for average:
      le_bcc(k) = le_bcc(k) / abs(le_ph(k-1)-le_ph(k))
      
      !~~ cloud liquid water content (grid cell average);
      ! integral over pressure ax is equivalent to integral over air mass,
      ! thus result is average clwc weighted with air mass ;
      ! negative pressure to have increasing ax:
      call IntervalQuad_Const( -1.0*mm_ph, mm_clwc, &
                               -1.0*le_ph(k-1), -1.0*le_ph(k), &
                               le_clwc(k), ilast, status )
      IF_NOTOK_RETURN(status=1)
      ! devide for average:
      le_clwc(k) = le_clwc(k) / abs(le_ph(k-1)-le_ph(k))
      
      !~~ cloud ice water content (grid cell average);
      ! integral over pressure ax is equivalent to integral over air mass,
      ! thus result is average clwc weighted with air mass ;
      ! negative pressure to have increasing ax:
      call IntervalQuad_Const( -1.0*mm_ph, mm_ciwc, &
                               -1.0*le_ph(k-1), -1.0*le_ph(k), &
                               le_ciwc(k), ilast, status )
      IF_NOTOK_RETURN(status=1)
      ! devide for average:
      le_ciwc(k) = le_ciwc(k) / abs(le_ph(k-1)-le_ph(k))
      
      ! Average potential rain intensity
      call IntervalQuad_Lin( -1.0*mm_ph, mm_pri, &
                             -1.0*le_ph(k-1), -1.0*le_ph(k), &
                               le_pri(k), ilast, status )  ! 1*Pa
      IF_NOTOK_RETURN(status=1)
      ! scale:       1*Pa   /           Pa
      le_pri(k) = le_pri(k) / abs(le_ph(k-1)-le_ph(k))   ! 1

    end do ! levels
    
    ! overhead cloud coverage is in-cloud + below-cloud:
    le_occ = le_icc + le_bcc

    !! info ...
    !if ( dbg ) then
    !  print *, 'ppp2 *** LE layers:'
    !  print *, 'ppp2 lv      ptop      pbot       icc       bcc        clwc      ciwc       ri'
    !  do k = le_n, 1, -1
    !    write (*,'(" ppp2 ",i2,2f10.2,2f10.6,2f12.8)') &
    !                  k, le_ph(k), le_ph(k-1), le_icc(k), le_bcc(k), le_clwc(k), le_ciwc(k), le_pri(k)
    !  end do
    !  print *, ''
    !  print *, 'Break after debugging Get_Cloud_Profiles'
    !  TRACEBACK; status=1; return
    !end if
    
    ! ok
    status = 0
    
  end subroutine Get_Cloud_Profiles

  ! ***


  subroutine SoilMoistureIndex_WRF( swvl, lsm, islt, t, smi, status )
 
 !AM: 31-3-2017 WRF implementation 
 !First order approach: translate soil types  (NOAH) to soil size classes (H-tessel, dominant land use used by ECMWF next to 
 ! other soil type classification, see Balsano et al).
    use GO, only : operator(<)
       
    ! --- in/out --------------
    
    real, intent(in)                  ::  swvl  ! volumetric soil water (m3 water)/(m3 soil)
    real, intent(in)                  ::  lsm   ! land-sea mask [0,1]
    integer, intent(in)               ::  islt  ! soil type index {0,1,..}
    type(TDate), intent(in)           ::  t
    real, intent(out)                 ::  smi   ! soil moisture index [0,1]
    integer, intent(out)              ::  status
    
    ! --- const ----------------

    character(len=*), parameter ::  rname = mname//'/SoilMoistureIndex'
    
    ! parameter database for calculate smi:parameters for permanent wilting point
    real, parameter   :: SW_min(6) = (/ 0.059, 0.151, 0.133, 0.279, 0.335, 0.267 /)
    real, parameter   :: SW_max(6) = (/ 0.242, 0.346, 0.382, 0.448, 0.541, 0.662 /)
    
    ! parameter database for constant soil type
    ! (in the ECMWF archive until 2007-11-06)
    real, parameter   :: SW_min_const = 0.171
    real, parameter   :: SW_max_const = 0.323
    
    ! why no soil type parameters ... ?
    
    ! --- local ----------------
    
    ! --- begin ----------------
    
    ! if cell is sea, smi = 1
    if ( lsm == 0.0 ) then
      ! all water:
      smi = 1.0
    else

      ! select on soil type:
      !~ water ?
      if ( islt == 14 ) then
        ! all water?
        smi = 1.0
      ! 
      else if ( islt == 15.or.islt==16 ) then
        ! bedrock or ice
        smi = 0.0
      !~ types 1:6 :
      else if ( islt == 1) then
      !sand->coarse
        ! switch in soil type:          
          smi = min( max( 0.0, (swvl - SW_min(1)) / (SW_max(1) - SW_min(1)) ), 1.0 )  
      else if ( islt==2.or.islt==3.or.islt==6.or.islt==7 ) then  
      !loamy->medium
          smi = min( max( 0.0, (swvl - SW_min(2)) / (SW_max(2) - SW_min(2)) ), 1.0 )
      else if ( islt==4.or.islt==5 ) then  
      !silty->medium fine
          smi = min( max( 0.0, (swvl - SW_min(3)) / (SW_max(3) - SW_min(3)) ), 1.0 )
      else if ( islt==8.or.islt==9.or.islt==10.or.islt==11 ) then  
      ! clayey->fine
          smi = min( max( 0.0, (swvl - SW_min(4)) / (SW_max(4) - SW_min(4)) ), 1.0 )
      else if ( islt==12) then 
      ! clay-> very fine 
          smi = min( max( 0.0, (swvl - SW_min(5)) / (SW_max(5) - SW_min(5)) ), 1.0 )
      else if (islt==13) then
      ! organic
          smi = min( max( 0.0, (swvl - SW_min(6)) / (SW_max(6) - SW_min(6)) ), 1.0 )
      ! 
      
      else
        write (gol,'("unsupported soil type ",i0)') islt; call goErr
        TRACEBACK; status=1; return
      end if  ! select soil type
        
    end if ! 100% sea
      
    ! ok
    status = 0
  
  end subroutine SoilMoistureIndex_WRF


  ! ***


  subroutine SoilMoistureIndex( swvl, lsm, islt, t, smi, status )
  
    use GO, only : operator(<)
       
    ! --- in/out --------------
    
    real, intent(in)                  ::  swvl  ! volumetric soil water (m3 water)/(m3 soil)
    real, intent(in)                  ::  lsm   ! land-sea mask [0,1]
    integer, intent(in)               ::  islt  ! soil type index {0,1,..}
    type(TDate), intent(in)           ::  t
    real, intent(out)                 ::  smi   ! soil moisture index [0,1]
    integer, intent(out)              ::  status
    
    ! --- const ----------------

    character(len=*), parameter ::  rname = mname//'/SoilMoistureIndex'
    
    ! parameter database for calculate smi
    real, parameter   :: SW_min(6) = (/ 0.059, 0.151, 0.133, 0.279, 0.335, 0.267 /)
    real, parameter   :: SW_max(6) = (/ 0.242, 0.346, 0.382, 0.448, 0.541, 0.662 /)
    
    ! parameter database for constant soil type
    ! (in the ECMWF archive until 2007-11-06)
    real, parameter   :: SW_min_const = 0.171
    real, parameter   :: SW_max_const = 0.323
    
    ! why no soil type parameters ... ?
    
    ! --- local ----------------
    
    ! --- begin ----------------
    
    ! if cell is sea, smi = 1
    if ( lsm == 0.0 ) then
      ! all water:
      smi = 1.0
    else

      ! select on soil type:
      !~ water ?
      if ( islt == 0 ) then
        ! all water?
        smi = 1.0
      ! 
      !~ types 1:6 :
      else if ( islt <= 6 ) then
        ! switch in soil type:
        if ( t < t_soiltype_switch ) then
          smi = min( max( 0.0, (swvl - SW_min_const) / (SW_max_const - SW_min_const) ), 1.0 )
        else   
          smi = min( max( 0.0, (swvl - SW_min(islt)) / (SW_max(islt) - SW_min(islt)) ), 1.0 )
        end if
      !
      !~ (which type is this?) treat as soil type 1 :
      else if ( islt == 7 ) then
        ! switch in soil type:
        if ( t < t_soiltype_switch ) then
          smi = min( max( 0.0, (swvl - SW_min_const) / (SW_max_const - SW_min_const) ), 1.0 )
        else 
          smi = min( max( 0.0, (swvl - SW_min(1)   ) / (SW_max(1)    - SW_min(1)   ) ), 1.0 )
        end if
      !
      else
        write (gol,'("unsupported soil type ",i0)') islt; call goErr
        TRACEBACK; status=1; return
      end if  ! select soil type
        
    end if ! 100% sea
      
    ! ok
    status = 0
  
  end subroutine SoilMoistureIndex


end module LE_Data_Calc
