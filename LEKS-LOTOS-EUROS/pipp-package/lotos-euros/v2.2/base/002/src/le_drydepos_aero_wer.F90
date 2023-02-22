!###############################################################################
!
! Aerosol dry deposition following 
!   Wesely (1985), Erisman (1994), and Ruijgrok (1994)
!
! Calculation of dry deposition velocity for fine particles according
! to Wesely et al. (1985) for low vegetation and other areas with
! a roughness length (znul) < 0.5m,
! and to Erisman et al (1994) and Ruijgrok et al. (1994) for forest 
! and other areas with a roughness length above 0.5m.
!
! Calculation of the collecting efficiency E for different components
! for different relative humidity and surface wetness.
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


module LE_DryDepos_Aero_WER

  implicit none
  
  private
  
  public  ::  vdpart_fine
  public  ::  vdpart_coarse


contains


  subroutine vdpart_fine (ust, ol , rh, nwet, znul, rstot)
 
    ! vegetation height: h=10*z0:
    integer :: nwet
    real :: ust, ol, rh, E, Uh, znul, rstot, vegh, d

    vegh = 10*znul
    d = 0.7*vegh
    if (ol.ne.0.)then
      Uh=(ust/0.41)*(log((vegh-d)/znul)- & 
             fpsih((vegh-d)/ol)+fpsih(znul/ol))                         
    else
      Uh=(ust/0.41)*log((vegh-d)/znul)
    endif

    if (rh.le.80.) then
            if (nwet.ge.1) then
                        E = 0.08*ust**0.45
            else
                        E = 0.05*ust**0.28
            endif
    else
            if (nwet.ge.1) then
                        E = 0.08*ust**0.45*(1.+0.37*EXP((rh-80.)/20.))                 
            else
                        E = 0.05*ust**0.28*(1.+0.18*EXP((rh-80.)/20.))
            endif
    endif

    ! Calculation of the deposition velocity
    ! for low vegetation and other areas with a roughness length
    ! below 0.5m by Wesely et al (1985) (for stable/neutral conditions
    ! (ol<0) and stable conditions (ol>0)), and for forest and other
    ! areas with a roughness length above 0.5m by Erisman et al (1994),

    if (znul.lt.0.5) then
          if (ol.lt.0.) then
              rstot = 1/ ((ust/500.)*(1.+((300./(-1.*ol))**2./3.)) )
          else
              rstot = 500/ust
          endif
    else
          rstot = 1./(ust*ust/Uh*E) 
    endif

  end subroutine


  ! ***       


  ! Calculation of dry deposition velocity for coarse particles: 
  ! Base cations BC (Na,Mg,Ca,K)
  !
  ! Calculation of the collecting efficiency E for different components
  ! for different relative humidity and surface wetness

  subroutine vdpart_coarse (ust, ol, rh, nwet, znul, rstot)

    !c vegetation height: h=10*z0:
    integer :: nwet
    real :: ust, ol, rh, E, Uh, znul, rstot, vegh, d

     vegh = 10*znul
     d = 0.7*vegh
     if (ol.ne.0.)then
     Uh=(ust/0.41)*(alog((vegh-d)/znul)- &
           fpsih((vegh-d)/ol)+fpsih(znul/ol))
     else
     Uh=(ust/0.41)*alog((vegh-d)/znul)
     endif

     if (rh.le.80.) then
                 if (nwet.ge.1) then
                         E = 0.679*ust**0.56
                 else
                         E = 0.140*ust**0.12
                 endif
     else
                 if (nwet.ge.1) then
                         E = 0.679*ust**0.56*(1.+0.37*EXP((rh-80.)/20.))
                 else
                         E = 0.140*ust**0.12*(1.-0.09*EXP((rh-80.)/20.))
                 endif
     endif
     rstot = 1./(ust*ust/Uh*E) 

  end subroutine


  !
  ! STABILITY CORRECTION FUNCTION IN THE SURFACE LAYER TEMPERATURE
  ! PROFILE
  !
  !  INPUT:  ETA   : STABILITY PARAMETER Z/L
  !  OUTPUT: FPSIH : CORRECTION IN LOGARITHMIC TEMPERATURE PROFILE
  !
  ! THE PRESENT MODEL IS AN EMPIRICAL FIT BY HOLTSLAG AND
  ! DE BRUIN(1987)OF DATA BY HICKS(1976, QUART. J. R. METEOR. SOC., 102,
  ! 535-551); SEE ALSO HOLTSLAG(1984, BLM, 29, 225-250)
  !

  real function fpsih(eta)

    real :: y, eta

    if (eta .lt. 0.) then
       y=sqrt(1.-16.*eta)
       fpsih=2*log((1+y)/2)
    else
         fpsih=-0.7*eta -(0.75*eta-14.29)*exp(-0.35*eta) -10.72
    endif

  end function fpsih


end module LE_DryDepos_Aero_WER
