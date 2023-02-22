!###############################################################################
!
! NAME
!
!   le_biascorr  -  apply bias corrections
!
!       updated for LuKwa v4 (LE v1.8.001)
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

module LE_BiasCorr

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------

  private

  !public  ::  biascorr__surface_ozone
  !public  ::  biascorr__total_pm__factor
  public  ::  biascorr__aod__factor
  
  public  ::  LE_BiasCorr_Init, LE_BiasCorr_Done
  public  ::  LE_BiasCorr_Fill
  

  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_BiasCorr'


  ! --- var --------------------------------

  ! description of current bias correction:
  character(len=64)   ::  biascorr__surface_ozone
  real                ::  biascorr__total_pm__factor
  real                ::  biascorr__aod__factor


contains


  ! ============================================================


  subroutine LE_BiasCorr_Init( rcf, status )

    use GO              , only : TrcFile, ReadRc
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_BiasCorr_Init'

    ! --- local -------------------------------

    ! --- begin -------------------------------
    
    ! read settings:
    call ReadRc( rcF, 'le.biascorr.surface_ozone', biascorr__surface_ozone, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'le.biascorr.total_pm.factor', biascorr__total_pm__factor, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'le.biascorr.aod.factor', biascorr__aod__factor, status )
    IF_NOTOK_RETURN(status=1)

    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_BiasCorr_Init


  ! ***


  subroutine LE_BiasCorr_Done( status )
  
    ! --- in/out ------------------------------

    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_BiasCorr_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------
    
    ! ok
    status = 0

  end subroutine LE_BiasCorr_Done


  ! ***


  subroutine LE_BiasCorr_Fill( c, cg, status )
  
    use Dims   , only : nx, ny, nz, nspec
    use Indices, only : specname
    use Indices, only : n_biascorr, ispecs_biascorr
    use Indices, only : ispec_o3_biascorr
    use Indices, only : i_o3
    use Indices, only : ispec_tpm25_biascorr, ispec_tpm10_biascorr
    use Indices, only : ispec_tpm25, ispec_tpm10
    use Indices, only : accum_n, accum_ii, accum_ww
    use Indices, only : tracer_is_dust, tracer_is_seasalt
  
    ! --- in/out ------------------------------

    real, intent(inout)                   ::  c(:,:,:,:)  ! (nx,ny,nz,nspec)
    real, intent(inout)                   ::  cg(:,:,:)   ! (nx,ny,nspec)
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_BiasCorr_Fill'

    ! --- local -------------------------------
    
    integer            ::  itr
    integer            ::  ispec
    integer            ::  ispec_tmp
    integer            ::  it, is
    real               ::  w
    integer            ::  iz

    ! --- begin -------------------------------
    
    ! loop over biascorrected tracers:
    do itr = 1, n_biascorr
      ! global index:
      ispec = ispecs_biascorr(itr)
      ! switch:
      select case ( ispec )
      
        ! biascorrected ozone for surface:
        case ( ispec_o3_biascorr )
        
          ! check ...
          if ( i_o3 < 0 ) then
            write (gol,'("no ozone tracer present to form bias corrected field")'); call goErr
            TRACEBACK; status=1; return
          end if
          ! check ...
          if ( ispec > size(c,4) ) then
            write (gol,'("concentration array has only ",i0," tracers and does not have entry for bias corrected tracer ",i0)') &
                             size(c,4), ispec; call goErr
            TRACEBACK; status=1; return
          end if
          ! copy:
          c (:,:,:,ispec) = c (:,:,:,i_o3)
          cg(:,:  ,ispec) = cg(:,:  ,i_o3)
          ! switch:
          select case ( biascorr__surface_ozone )
            !
            case ( 'none' )
              ! nothing ...
            !
            case ( 'o3_smogprog_v19' )
              !! info ..
              !write (gol,'("apply ozone biascorrection, smoprog v1.9, tracer o3_biascorr")'); call goPr
              ! apply ozone correction to surface field:
              call biascorr__o3_smogprog_v19( cg(:,:,ispec), status )
              IF_NOTOK_RETURN(status=1)
              ! also apply to levels to have at least some values there:
              do iz = 1, nz
                call biascorr__o3_smogprog_v19( c(:,:,iz,ispec), status )
                IF_NOTOK_RETURN(status=1)
              end do
            !
            case ( 'o3_smogprog_v110' )
              !! info ..
              !write (gol,'("apply ozone biascorrection, smoprog v1.8, tracer o3_biascorr")'); call goPr
              ! apply ozone correction:
              call biascorr__o3_smogprog_v110( cg(:,:,ispec), status )
              IF_NOTOK_RETURN(status=1)
              ! also apply to levels to have at least some values there:
              do iz = 1, nz
                call biascorr__o3_smogprog_v110( c(:,:,iz,ispec), status )
                IF_NOTOK_RETURN(status=1)
              end do
            !
            case default
              write (gol,'("unsupported surface_ozone bias correction : ",a)') trim(biascorr__surface_ozone); call goErr
              TRACEBACK; status=1; return
          end select
          
        ! biascorrected tpm10: multiply with factor, except for dust and seasalt parts:
        case ( ispec_tpm25_biascorr, ispec_tpm10_biascorr )

          ! check ...
          if ( ispec > size(c,4) ) then
            write (gol,'("concentration array has only ",i0," tracers and does not have entry for bias corrected tracer ",i0)') &
                             size(c,4), ispec; call goErr
            TRACEBACK; status=1; return
          end if
          ! init sum:
          c (:,:,:,ispec) = 0.0
          cg(:,:  ,ispec) = 0.0
          ! template accumulated tracer:
          if ( ispec == ispec_tpm25_biascorr ) then
            ispec_tmp = ispec_tpm25
          else if ( ispec == ispec_tpm10_biascorr ) then
            ispec_tmp = ispec_tpm10
          else
            write (gol,'("no template accumulated tracer defined for ispec ",i0," `",a,"`")') ispec, trim(specname(ispec)); call goErr
            TRACEBACK; status=1; return
          end if
          ! check ...
          if ( ispec_tmp < 0 ) then
            write (gol,'("originating tpmNN tracer not present to form bias corrected field")'); call goErr
            TRACEBACK; status=1; return
          end if
          ! loop over components of accumulated tracer:
          do it = 1, accum_n(ispec_tmp)
            ! global index:
            is = accum_ii(ispec_tmp,it)
            ! extra factor:
            if ( tracer_is_dust(is) .or. tracer_is_seasalt(is) ) then
              ! no extra factor:
              w = 1.0
            else
              ! factor to compensate for missing aerosols:
              w = biascorr__total_pm__factor
            end if
            ! add contributions:
            c (:,:,:,ispec) = c (:,:,:,ispec) + c (:,:,:,is) * accum_ww(ispec_tmp,it) * w
            cg(:,:  ,ispec) = cg(:,:  ,ispec) + cg(:,:  ,is) * accum_ww(ispec_tmp,it) * w
          end do ! accumulation
          
        ! unkown ...
        case default
          write (gol,'("unsupported biascorr ispec ",i0," `",a,"`")') ispec, trim(specname(ispec)); call goErr
          TRACEBACK; status=1; return
      end select
      
    end do  ! biascorrected tracers
    
    ! ok
    status = 0

  end subroutine LE_BiasCorr_Fill


  ! ====================================================


  !
  ! Conversion from ug/m3 -> ppb :
  !
  !   ug O3    m3    (kg air)  (kg O3) (mol O3)   (mol O3)
  !   ----- -------- --------- ------- -------- = ---------
  !     m3  (kg air) (mol air) (ug O3) (kg O3)    (mol air)
  !
  !   ugm3 /  dens  *  xm_air  * 1e-9 /  xm_O3 * 1e9   ! ppb
  !
  !   ugm3 /  dens  *  xm_air  /  xm_O3   ! ppb
  !
  !   ugm3 /  dens  *  xm_air  /  xm_O3   ! ppb
  !
  ! Density is computed in meteo module using density at 1000 hPa :
  !
  !    dens(i,j,l) =  p_hPa(l) / 1000.0 * 1.293
  !
  ! Example for ozone:
  !
  !   ugm3 /  1.293 * 28.964 / 48.0 = ugm3 * 0.47
  !
  ! Formula suggested by SCJ :
  !
  !   ppb2ug =  (48.00*0.001*101300)/(1.3807*283*6.022)
  !   ug2ppb = (1.3807*283*6.022) / (48.00*0.001*101300)
  !
  ! Usage in LE:
  !
  !   ! single layer, will with 10.00 ug/m3 O3 :
  !   csurf = ugm3_to_ppb( 10.0, xm_O3, dens(:,:,1) )
  !
  
  elemental real function ugm3_to_ppb( ugm3, xm_tracer, airdens )   ! concentration (ppb)
  
    use Binas, only : xm_air      ! 28.964e-3   kg/mol
  
    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  ugm3         ! concentration (ug/m3)
    real, intent(in)      ::  xm_tracer    ! tracer mol mass (kg/mol)
    real, intent(in)      ::  airdens      ! air density (kg/m3)
    
    ! --- begin ----------------------------------
    
    ugm3_to_ppb = ugm3 /  airdens  *  xm_air  /  xm_tracer   ! ppb
    
  end function ugm3_to_ppb
  
  
  ! ====================================================    
  ! 
  ! Version for LE v1.9.000
  !
  !  Tsurf <= 20.0 C       -> c = 4.20 ug/m3
  !  20 C < Tsurf < 40 C   -> c = (0.051081 Tsurf^2 -1.2418 Tsurf + 8.6004)  ug/m3
  !  Tsurf >= 40 C         -> c = 40.66 ug/m3
  !
  ! N.B. de biascorrectie is hier in ug/m3, in LE/mix2ground.F90 moet 
  ! dat nog eerst omgerekend worden naar ppb.
  !
  ! Ref: "Biascorrection LOTOS-EUROS v1.9"
  !        Ferd Sauter, 2013-04-16
  !
  !   "LOTOS-EUROS bias in daily maxima ozone for rural stations in the Netherlands"
  !   2011-02-23, Ferd Sauter, RIVM
  !
  ! 2012-03-23, Henk Eskes: 
  !   An extra multiplication factor 1.02 is added to account for 
  !   the new RIVM ozone sensors which have higher ozone readings than the old ones
  !

  subroutine biascorr__o3_smogprog_v19( cg_o3, status )
  
    use Binas, only : xm_O3
    use dims,  only : nx, ny
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------------
    
    real, intent(inout)       ::  cg_o3(nx,ny)
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/biascorr__o3_smogprog'
    
    ! offset for degree Kelvin to degree Celcius
    real, parameter       ::  kelvin_to_celcius = -273.15
    
    ! --- local ----------------------------------
    
    real          ::  c0
    integer       ::  i, j
    real, pointer               ::  tsurf(:,:,:)   ! (lon,lat,1)
    real, pointer               ::  dens(:,:,:)   ! (lon,lat,alt)    
    
    ! --- begin ----------------------------------

    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units ='K')
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')   
    IF_NOTOK_RETURN(status=1)

    ! Apply ozone bias correction, Ferd Sauter    

    do i = 1, nx
       do j = 1, ny

          if ( tsurf(i,j,1)+kelvin_to_celcius  < 20.0 ) then
             cg_o3(i,j) = cg_o3(i,j) + ugm3_to_ppb( 4.20, xm_O3, dens(i,j,1) )
          end if

          if ( (20.0 <= tsurf(i,j,1)+kelvin_to_celcius) .and. (tsurf(i,j,1)+kelvin_to_celcius < 40.0) ) then
             cg_o3(i,j) = cg_o3(i,j) + &
                  ugm3_to_ppb(   0.051081 * (tsurf(i,j,1)+kelvin_to_celcius)**2 &
                               - 1.2418  * (tsurf(i,j,1)+kelvin_to_celcius) &
                               + 8.6004 , xm_O3, dens(i,j,1)                )
          end if

          if ( tsurf(i,j,1)+kelvin_to_celcius  >= 40.0 ) then
             cg_o3(i,j) = cg_o3(i,j) + ugm3_to_ppb( 40.66, xm_O3, dens(i,j,1) )
          end if

       end do
    end do

    ! Apply extra correction for new RIVM ozone sensors

    cg_o3 = cg_o3 * 1.02

    ! ok
    status = 0
    
  end subroutine biascorr__o3_smogprog_v19


  ! ====================================================    
  ! 
  ! Version for LE v1.10.000
  !
  !  Tsurf <= 15.0 C       -> c = 10.23 ug/m3
  !  15 C < Tsurf < 35 C   -> c = (0.132443 Tsurf^2 -4.23369 Tsurf + 43.9356)  ug/m3
  !  Tsurf >= 35 C         -> c = 58.00 ug/m3
  !
  ! N.B. de biascorrectie is hier in ug/m3, in LE/mix2ground.F90 moet 
  ! dat nog eerst omgerekend worden naar ppb.
  !
  ! Ref: "Biascorrection LOTOS-EUROS v1.10"
  !        Ferd Sauter, 2014-03-20
  !
  !   "LOTOS-EUROS bias in daily maxima ozone for rural stations in the Netherlands"
  !   2011-02-23, Ferd Sauter, RIVM
  !
  ! 2014-04-01, Ujjwal Kumar: 
  !Up till now, an extra correction for the new RIVM monitors has been supplied
  !inside the LOTOS-EUROS bias correction module. From v1.10, this should be kept
  !outside the code, since the corrections have already been used in the bias
  !correction method described in this report. (Ferd Sauter, 2014-03-20)
  !

  subroutine biascorr__o3_smogprog_v110( cg_o3, status )
  
    use Binas, only : xm_O3
    use dims,  only : nx, ny
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------------
    
    real, intent(inout)       ::  cg_o3(nx,ny)
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/biascorr__o3_smogprog'
    
    ! offset for degree Kelvin to degree Celcius
    real, parameter       ::  kelvin_to_celcius = -273.15
    
    ! --- local ----------------------------------
    
    real          ::  c0
    integer       ::  i, j

    real, pointer               ::  tsurf(:,:,:)   ! (lon,lat,1)
    real, pointer               ::  dens(:,:,:)   ! (lon,lat,alt)
    
    ! --- begin ----------------------------------

    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units ='K')
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')             
    IF_NOTOK_RETURN(status=1)

    ! Apply ozone bias correction, Ferd Sauter    

    do i = 1, nx
       do j = 1, ny

          if ( tsurf(i,j,1)+kelvin_to_celcius  < 15.0 ) then
             cg_o3(i,j) = cg_o3(i,j) + ugm3_to_ppb( 10.23, xm_O3, dens(i,j,1) )
          end if

          if ( (15.0 <= tsurf(i,j,1)+kelvin_to_celcius) .and. (tsurf(i,j,1)+kelvin_to_celcius < 35.0) ) then
             cg_o3(i,j) = cg_o3(i,j) + &
                  ugm3_to_ppb(   0.132443 * (tsurf(i,j,1)+kelvin_to_celcius)**2 &
                               - 4.23369  * (tsurf(i,j,1)+kelvin_to_celcius) &
                               + 43.9356 , xm_O3, dens(i,j,1)                )
          end if

          if ( tsurf(i,j,1)+kelvin_to_celcius  >= 35.0 ) then
             cg_o3(i,j) = cg_o3(i,j) + ugm3_to_ppb( 58.00, xm_O3, dens(i,j,1) )
          end if

       end do
    end do


    ! ok
    status = 0
    
  end subroutine biascorr__o3_smogprog_v110


end module LE_BiasCorr
