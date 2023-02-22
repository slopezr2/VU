!#######################################################################
!
! LE_Emis_Fires_Data
!
! Emissions from 'MODIS' .
!
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!#######################################################################


module LE_Emis_Fires_Data

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  public
  

  ! --- const ------------------------------------

  ! >>> these factors are copied from the SILAM forest fires module;
  !     since absolute emission is unknown anyway, these numbers are of for ratio
  ! update following AQMEII-2 protocol, which is based on Andrea and Merlet (but differs from 
  ! Transphorm settings)
  
  ! fraction NO from NOx   from  Andrea and Merlet 2001 
  real, parameter   ::  NO_rfire = 0.9

  ! Unit of emissions is (kg PM2.5)/sec 
  ! Ratios between PM2.5 emission factors and other constituants (ref ?)
  ! factor results
  real, parameter   ::  rEF_CO   = 6.08 !8.15
  real, parameter   ::  rEF_CH4  = 0.37
  real, parameter   ::  rEF_NOx  = 0.17 !0.23
  real, parameter   ::  rEF_NO   = rEF_NOx *  NO_rfire
  real, parameter   ::  rEF_NO2  = rEF_NOx * (1.0-NO_rfire)

  real, parameter   ::  rEF_VOC  = 0.32!0.44
  real, parameter   ::  rEF_SO2  = 0.057!0.08
  real, parameter   ::  rEF_NH3  = 0.08 !0.11
  ! 74% of mass emitted as pm2.5
  !part of the pm25 is emitted as black-carbon or oc:
  real, parameter   ::  rEF_EC_F   = 0.043           !  4.3% of pm2.5
  real, parameter   ::  rEF_POM_F =  0.52            !70.4% of PM2.5
  real, parameter   ::  rEF_PPM_f   = 0.187           !25.3% of pm2.5
  real, parameter   ::  rEF_PPM_C = 0.26             !26% of tpm
  ! VOC partionning UNIT [mol/kg VOC] 
  !   REF ?
  real, parameter   ::  xfire_OLE  =  0.45
  real, parameter   ::  xfire_PAR  =  7.079
  real, parameter   ::  xfire_TOL  =  0.218
  real, parameter   ::  xfire_XYL  =  0.094
  real, parameter   ::  xfire_FORM =  1.036
  real, parameter   ::  xfire_ALD  =  1.698
 !real, parameter   ::  xfire_KET  =  0.
 !real, parameter   ::  xfire_ACET =  0.
  real, parameter   ::  xfire_ETH  =  5.355
  real, parameter   ::  xfire_UNR  = 38.0
  
  !<<<

  
  ! --- var --------------------------------------
  
  ! by default not enabled ...
  logical           ::  with_ffire = .false.

  ! 2d field with released bio-fuel emission:
  real(4), allocatable          ::  ffire_release(:,:)  ! kg/cell/s
  !real(4), allocatable          ::  ffire_proba  (:,:)  ! confidence
  
  !! collect emissions:
  !real, allocatable             ::  emis_ffire(:,:,:,:) ! nx, ny, nz, nspec


end module LE_Emis_Fires_Data
