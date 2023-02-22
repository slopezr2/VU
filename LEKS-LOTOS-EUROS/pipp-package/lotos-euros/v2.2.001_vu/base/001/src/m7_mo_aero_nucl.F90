MODULE mo_aero_nucl

  ! Astrid Manders, July 2013,  extension with Activation parameterization
  ! based on  Kulmala et al (2006, ACP), since the Vehkamaki scheme gives 
  ! too little nucleation in boundary layer, intended for free troposphere.
  ! Approach of Vehkamaki scheme used, no analytical integration.

  ! *mo_aero_nucl* holds routines calculating the nucleation rate
  !
  ! Author:
  ! -------
  ! Philip Stier, MPI-Met, Hamburg,                 01/2003
  !
  ! Purpose:                                                           
  ! --------
  ! This module holds the routines for the calculation of the nucleation
  ! rate for the binary nucleation in the sulfate water system. 
  ! These routines are called and applied within HAM from the routine
  ! m7_nuck.
  ! The old parameterization of Kulmala (1998), that apparently contains
  ! some errors is kept for consistency. It is recommended to use the 
  ! new scheme of Vehkamaeki et al. (2002). However, this parameterization
  ! is no longer analytically integrable, the effect of this is to be 
  ! scrutinized.
  ! The modularized version of the Kulmala parameterization has been tested
  ! to give bit-identical results with the old version.
  ! 
  ! References:
  ! -----------
  ! Vehkamaeki et al. (2002), An improved parameterization for sulfuric
  !    acid/water nucleation rates for tropospheric and stratospheric
  !    conditions, J. Geophys. Res, 107, D22, 4622
  ! Kulmala et al. (1998), Parameterizations for sulfuric acid/water
  !    nucleation rates. J. Geophys. Res., 103, No D7, 8301-8307
  ! Vignati, E. (1999), Modelling Interactions between Aerosols and 
  !    Gaseous Compounds in the Polluted Marine Atmosphere. PhD-Thesis,
  !    RISO National Laborartory Copenhagen, Riso-R-1163(EN)


  IMPLICIT NONE


CONTAINS

  subroutine nucl_kulmala2006(kproma,  kbdim,  klev,    &  ! ECHAM5 dims
                           pso4g,   ptp1,   &
                           pxtrnucr, pntot                )  !nucleation rate, number of molecules in the
                                                            ! critical cluster
  !J_act=A[H2SO4], formation rate in cm-3 s-1, Kulmala et al ACP 2006
  !critical cluster size=1nm       Following Fountoukis et al ACPD 2012                 
  
  implicit none
  INTEGER :: kproma, kbdim, klev
  
  INTEGER :: jk, jl,jrow
  
  REAL    ::   ptp1(kbdim,klev), prhd(kbdim,klev), &             
               pso4g(kbdim,klev), &            
               pxtrnucr(kbdim,klev), pntot(kbdim,klev)
               
  real, parameter ::  Act=1e-6  !nucleation rate coefficient, Spracklen et al 2009, Kulmala 2006, Fountoukis 2012
                                !should be in range 1e-6 to 1e-7 per second (Kulmala 2006), 2e-6 maybe too large
                                ! (integration step too large-> too much production) 
  real, parameter ::  critsize= 1e-9 !critical cluster size=1nm=1e-9m Fountoukis et al 2012
	
  real       :: zjnuc  !   
    DO jk=1, klev
     DO jl=1,kproma
       
         
         !convert total number of molecules in the critical cluster
        ! to number of sulfate molecules:
        zjnuc=Act*pso4g(jl,jk)
        !pntot(jl,jk)=zntot*zxmole  

        ! limitation to 1E+10 [1/cm3s]
       
        zjnuc=MIN(zjnuc,1.e10)
        IF(zjnuc<1.e-7)  zjnuc=0
        pxtrnucr(jl,jk) = zjnuc
        pntot(jl,jk)=1 !critical clusters contain 1 molecule of H2SO4, Kulmala et al 2006
        ! convert total number of molecules in the critical cluster
        ! to number of sulfate molecules:
        ! critical size =1nm, volume=(4/3)*pi*critsize**3 
        !pntot(jl,jk)=zntot*zxmole
        !assume 
        

     end do
    end do
  
  end subroutine nucl_kulmala2006

  SUBROUTINE nucl_kulmala(kproma,  kbdim,  klev,    &  ! ECHAM5 dims
                          pso4g,   ptp1,   prelhum, & 
                          pbnrate, palpha, pbeta    )

    
    !  Authors:
    !  --------
    !  P. Stier, MPI-Met, Hamburg,    from the original f77 code
    !                                 in the routine m7_nuck        2001-2003
    !  J. Wilson, E. Vignati, JRC/EI, original source                 09/2000
    !
    !  Purpose:                                                           
    !  --------                                                           
    !  This routine calculates the instananeous nucleation rate            
    !  znucrate [molec. cm-3 s-1] from a given gas phase H2SO4 concentration      
    !  pso4g [molec. cm-3]. It also calculates the integrated change of 
    !  H2SO4 gas phase mass over one timestep due to nucleation 
    !  pa4delt(:,:,1) [molec. cm-3] as well as the number of nucleated 
    !  particles panew [1] during the timestep.
    !
    !  Interface:
    !  ----------
    !  *nucl_kulmala* is called from *m7_nuck*
    !
    !  Method:
    !  -------
    !  Kulmala et al. (1998) 's formula for binary nucleation is 
    !  rewritten to take the form znucrate = exp[zalpha+ln(pso4g)*beta]. 
    !  Equation numbers are taken from Kulmala et al. (1998).
    !  After the calculation of the nucleation rate znucrate, it is 
    !  integrated in 2) analytically over one timestep, i.e.:
    !
    !  Integration of:
    ! 
    !  znucrate=d(critn*znav)/dt=exp[zalpha + zbeta*ln(znav)]
    !
    !  gives znav(t0+dt, znav(t0)) where znav(t0) is pso4g.
    !  znav is temporarily stored in zso4g_new and confined between
    !  0 and pso4g. 
    !  The number of nucleated particles is then calculated by 
    !  assuming a fixed critical mass critn of newly nucleated 
    !  particles and dividing the total mass of nucleated sulfate 
    !  by this critical mass of one nucleated particle. 
    !
    !  Externals:
    !  ----------
    !  None

    USE mo_aero_m7, ONLY: bk

    IMPLICIT NONE

    INTEGER :: kproma,      kbdim,        klev

    REAL    :: pso4g(kbdim,klev),         ptp1(kbdim,klev),      &
               prelhum(kbdim,klev),       pbnrate(kbdim,klev),   &
               palpha(kbdim,klev),        pbeta(kbdim,klev)

    INTEGER :: jl,          jk

    REAL    :: znwv,        zln_nac,      ztk,         zsupsat,  &
               zpeh2o,      zpeh2so4,     zra,         zxal,     &
               ztkn,        zssn,         zdelta


    !---1) Calculation of the nucleation rate: ----------------------------
    
    DO jk=1,klev
       DO jl=1,kproma

          IF (pso4g(jl,jk) .GT. 1e-5) THEN
             ztk=ptp1(jl,jk)
             zsupsat=prelhum(jl,jk)
             !
             !--- 1.1) Restrict t, and rh to limits where the parameterization ok:
             !
             ztkn = MAX(ztk, 220.0)
             zssn = MIN(zsupsat, 0.90)

             !
             !--- 1.2) Equlibrium vapour pressures (Jaeker-Mirabel (1995), JGR):
             !
             !--- H2O equlibrium vapour pressure (Tabata):
             !
             zpeh2o=0.750064*(10.**(8.42926609-1827.17843/          &
                  ztkn-71208.271/ztkn/ztkn))*1333./bk/ztkn
             !
             !--- H2SO4 equlibrium vapour pressure at 360 
             !@@@ Check source: which Ayers?
             !
             zpeh2so4=EXP(-10156./ztkn+16.259)*7.6e2*1333./bk/ztkn
             !
             !--- H2SO4 equlibrium vapour pressure - correction of ayers
             !    by kulmala - currently not used
             !
             !     payers=exp(-10156/360+16.259)*7.6e2
             !     zpeh2so4=exp(log(payers)+10156*(-1./ztkn+1./360.+0.38/(905-360) * &
             !              (1+log(360./ztkn)-360./ztkn)))*1333/bk/ztkn
             !
             !--- 1.3) Relative acidity (0.0 -1.0):
             ! 
             zra=pso4g(jl,jk)/zpeh2so4
             !
             !--- 1.4) Water vapour molecule concentration [cm-3]:
             !
             znwv=zsupsat*zpeh2o
             !
             !--- 1.5) Factor delta in Eq. 22:
             ! 
             zdelta=1.0+(ztkn-273.15)/273.15
             !
             !--- 1.6) Molefraction of H2SO4 in the critical cluster 
             !         minus the H2SO4(g) term in Eq. 17:
             !
             zxal = 1.2233-0.0154*zra/(zra+zssn)-0.0415*LOG(znwv)+ 0.0016*ztkn 
             !
             !--- 1.7) Exponent of the critical cluster (Eq. 18):
             !
             zln_nac = -14.5125+0.1335*ztkn-10.5462*zssn+1958.4*zssn/ztkn
             !
             !--- 1.8) Sum of all terms in Eq. 20 containing H2SO4(g):
             !
             pbeta(jl,jk) = 25.1289 - 4890.8/ztkn + 7643.4*0.0102/ztkn - &
                            2.2479*zdelta*zssn - 1.9712*0.0102*zdelta/zssn
             ! 
             !--- 1.9) Sum all terms in Eq. 20 not containing H2SO4(g):
             !
             palpha(jl,jk) = zln_nac*(-25.1289 + 4890.8/ztkn + 2.2479*zdelta*zssn) - &
                             1743.3/ztkn + zxal*(7643.4/ztkn - 1.9712*zdelta/zssn)
             !
             !--- 1.10) Nucleation rate [cm-3 s-1] (Kulmala et al., 1998):
             !
             pbnrate(jl,jk) = EXP(palpha(jl,jk)+LOG(pso4g(jl,jk))*pbeta(jl,jk))

          ELSE

             palpha(jl,jk) =0. 
             pbeta(jl,jk)  =0.
             pbnrate(jl,jk)=0.

          END IF ! pso4g(jl,jk) .GT. 1e-5
       END DO ! kproma
    END DO !klev


  END SUBROUTINE nucl_kulmala



  SUBROUTINE nucl_vehkamaeki(kproma,   kbdim, klev,        &  ! ECHAM5 dimensions
                             ptp1,     prhd,  pmolecH2SO4, &  ! ECHAM5 temperature, relative humidity
                             pxtrnucr, pntot               )  ! nucleation rate, number of molecules in the
                                                              ! critical cluster
			 
    !
    !   Authors:
    !   ---------
    !   C. TIMMRECK, MPI HAMBURG                                             2002
    !
    !   Purpose
    !   ---------
    !   Calculation of classical nucleation rate
    !               
    !   calculation of the nucleation rate after Vehkamaeki et al. (2002)
    !   The calculation of the nucrate ZKNH2SO4 is in cm^-3 s^-1
    !   and a coarse approxmation for the first class
    !
    !   Modifications:
    !   --------------
    !   R. Hommel; rewrite in f90, adopted to ECHAM5;        MPI HAMBURG;      Dec. 2002
    !   P. Stier; bugfixes, modularisation and optimization; MPI HAMBURG;      2003-2004
    !
    !   H2SO4 still fixed to xxx molc/cm3, no sulfur cycle coupling yet
    !
    !   References:
    !   -----------
    !   Vehkamaeki et al. (2002), An improved parameterization for sulfuric
    !      acid/water nucleation rates for tropospheric and stratospheric
    !      conditions, J. Geophys. Res, 107, D22, 4622
    !
    !   Parameters
    !   ----------
    !   prho = prhop_neu in *sam*
    !
    !   prhd = relative humidity in sam_aeroprop & sam_nucl
    !
    !   pxtrnucr = nucleation rate in [1/m3s]
    !   xrhoc    = density of the critical nucleus in kg/m^3
    !   zrxc = ?
    !
			     
  USE mo_control,  ONLY: nrow


  !----------------------------------------------------

  IMPLICIT NONE

  !----------------------------------------------------
  
  INTEGER :: kproma, kbdim, klev
  
  INTEGER :: jk, jl,jrow
  
  !----------------------------------------------------
  !
  
  REAL    ::   ptp1(kbdim,klev), prhd(kbdim,klev), &
               pxtrnucr(kbdim,klev),  &
               pmolecH2SO4(kbdim,klev), &            ! revisited, ok
               pntot(kbdim,klev)
	       
  !----------------------------------------------------  
  ! Local Arrays
  
  REAL    ::   zrxc(kbdim) 

  REAL    ::   zrhoa, zrh, zt, x, zjnuc, zrc, &
               zxmole, zntot


  !--- 0) Initializations:

  jrow=nrow(2)

  DO jk=1, klev
     DO jl=1,kproma

  !----1.) Parameterization of  nucleation rate after Vehkamaeki et al. (2002)

        ! t: temperature in K (190.15-300.15K)                                  
        ! zrh: saturatio ratio of water (0.0001-1)                               
        ! zrhoa: sulfuric acid concentration in 1/cm3 (10^4-10^11 1/cm3)         
        ! jnuc: nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)                  
        ! ntot: total number of molecules in the critical cluster (ntot>4)      
        ! x: molefraction of H2SO4 in the critical cluster                      
        ! rc: radius of the critical cluster in nm                              

        ! Calculate nucleation only for valid thermodynamic conditions:

        IF( (pmolecH2SO4(jl,jk)>=1.E+4)                      .AND. &
            (prhd(jl,jk) >=1.E-4)                            .AND. &
            (ptp1(jl,jk)>=190.15 .AND. ptp1(jl,jk)<=300.15)      ) THEN

        zrhoa=MIN(pmolecH2SO4(jl,jk),1.e11)
        zrh=MIN(prhd(jl,jk),1.0)
        zt=ptp1(jl,jk)

        ! Equation (11) - molefraction of H2SO4 in the critical cluster

        x=0.7409967177282139 - 0.002663785665140117*zt   &
          + 0.002010478847383187*LOG(zrh)    &
          - 0.0001832894131464668*zt*LOG(zrh)    &
          + 0.001574072538464286*LOG(zrh)**2        &  
          - 0.00001790589121766952*zt*LOG(zrh)**2    &
          + 0.0001844027436573778*LOG(zrh)**3     &
          -  1.503452308794887e-6*zt*LOG(zrh)**3    &
          - 0.003499978417957668*LOG(zrhoa)   &
          + 0.0000504021689382576*zt*LOG(zrhoa)   

        zxmole=x

        ! Equation (12) - nucleation rate in 1/cm3s
	  
        zjnuc=0.1430901615568665 + 2.219563673425199*zt -   &
              0.02739106114964264*zt**2 +     &
              0.00007228107239317088*zt**3 + 5.91822263375044/x +     &
              0.1174886643003278*LOG(zrh) + 0.4625315047693772*zt*LOG(zrh) -     &
              0.01180591129059253*zt**2*LOG(zrh) +     &
              0.0000404196487152575*zt**3*LOG(zrh) +    &
              (15.79628615047088*LOG(zrh))/x -     &
              0.215553951893509*LOG(zrh)**2 -    &
              0.0810269192332194*zt*LOG(zrh)**2 +     &
              0.001435808434184642*zt**2*LOG(zrh)**2 -    & 
              4.775796947178588e-6*zt**3*LOG(zrh)**2 -     &
              (2.912974063702185*LOG(zrh)**2)/x -   &
              3.588557942822751*LOG(zrh)**3 +     &
              0.04950795302831703*zt*LOG(zrh)**3 -     &
              0.0002138195118737068*zt**2*LOG(zrh)**3 +    & 
              3.108005107949533e-7*zt**3*LOG(zrh)**3 -     &
              (0.02933332747098296*LOG(zrh)**3)/x +     &
              1.145983818561277*LOG(zrhoa) -    &
              0.6007956227856778*zt*LOG(zrhoa) +    &
              0.00864244733283759*zt**2*LOG(zrhoa) -    &
              0.00002289467254710888*zt**3*LOG(zrhoa) -    &
              (8.44984513869014*LOG(zrhoa))/x +    &
              2.158548369286559*LOG(zrh)*LOG(zrhoa) +   & 
              0.0808121412840917*zt*LOG(zrh)*LOG(zrhoa) -    &
              0.0004073815255395214*zt**2*LOG(zrh)*LOG(zrhoa) -   & 
              4.019572560156515e-7*zt**3*LOG(zrh)*LOG(zrhoa) +    &
              (0.7213255852557236*LOG(zrh)*LOG(zrhoa))/x +    &
              1.62409850488771*LOG(zrh)**2*LOG(zrhoa) -    &
              0.01601062035325362*zt*LOG(zrh)**2*LOG(zrhoa) +   & 
              0.00003771238979714162*zt**2*LOG(zrh)**2*LOG(zrhoa) +    &
              3.217942606371182e-8*zt**3*LOG(zrh)**2*LOG(zrhoa) -    &
              (0.01132550810022116*LOG(zrh)**2*LOG(zrhoa))/x +    &
              9.71681713056504*LOG(zrhoa)**2 -    &
              0.1150478558347306*zt*LOG(zrhoa)**2 +    &
              0.0001570982486038294*zt**2*LOG(zrhoa)**2 +    &
              4.009144680125015e-7*zt**3*LOG(zrhoa)**2 +    &
              (0.7118597859976135*LOG(zrhoa)**2)/x -    &
              1.056105824379897*LOG(zrh)*LOG(zrhoa)**2 +    &
              0.00903377584628419*zt*LOG(zrh)*LOG(zrhoa)**2 -    &
              0.00001984167387090606*zt**2*LOG(zrh)*LOG(zrhoa)**2 +    &
              2.460478196482179e-8*zt**3*LOG(zrh)*LOG(zrhoa)**2 -    &
              (0.05790872906645181*LOG(zrh)*LOG(zrhoa)**2)/x -    &
              0.1487119673397459*LOG(zrhoa)**3 +    &
              0.002835082097822667*zt*LOG(zrhoa)**3 -    &
              9.24618825471694e-6*zt**2*LOG(zrhoa)**3 +    &
              5.004267665960894e-9*zt**3*LOG(zrhoa)**3 -    &
              (0.01270805101481648*LOG(zrhoa)**3)/x
        
        zjnuc=EXP(zjnuc)      !   add. Eq. (12) [1/(cm^3s)]      


        ! Equation (13) - total number of molecules in the critical cluster

        zntot=-0.002954125078716302 - 0.0976834264241286*zt +   &
               0.001024847927067835*zt**2 - 2.186459697726116e-6*zt**3 -    &
               0.1017165718716887/x - 0.002050640345231486*LOG(zrh) -   &
               0.007585041382707174*zt*LOG(zrh) +    &
               0.0001926539658089536*zt**2*LOG(zrh) -   &
               6.70429719683894e-7*zt**3*LOG(zrh) -    &
               (0.2557744774673163*LOG(zrh))/x +   &
               0.003223076552477191*LOG(zrh)**2 +   &
               0.000852636632240633*zt*LOG(zrh)**2 -    &
               0.00001547571354871789*zt**2*LOG(zrh)**2 +   &
               5.666608424980593e-8*zt**3*LOG(zrh)**2 +    &
               (0.03384437400744206*LOG(zrh)**2)/x +   &
               0.04743226764572505*LOG(zrh)**3 -    &
               0.0006251042204583412*zt*LOG(zrh)**3 +   &
               2.650663328519478e-6*zt**2*LOG(zrh)**3 -    &
               3.674710848763778e-9*zt**3*LOG(zrh)**3 -   &
               (0.0002672510825259393*LOG(zrh)**3)/x -    &
               0.01252108546759328*LOG(zrhoa) +   &
               0.005806550506277202*zt*LOG(zrhoa) -    &
               0.0001016735312443444*zt**2*LOG(zrhoa) +   &
               2.881946187214505e-7*zt**3*LOG(zrhoa) +    &
               (0.0942243379396279*LOG(zrhoa))/x -   &
               0.0385459592773097*LOG(zrh)*LOG(zrhoa) -   & 
               0.0006723156277391984*zt*LOG(zrh)*LOG(zrhoa) +   &
               2.602884877659698e-6*zt**2*LOG(zrh)*LOG(zrhoa) +    &
               1.194163699688297e-8*zt**3*LOG(zrh)*LOG(zrhoa) -   &
               (0.00851515345806281*LOG(zrh)*LOG(zrhoa))/x -    &
               0.01837488495738111*LOG(zrh)**2*LOG(zrhoa) +   &
               0.0001720723574407498*zt*LOG(zrh)**2*LOG(zrhoa) -   & 
               3.717657974086814e-7*zt**2*LOG(zrh)**2*LOG(zrhoa) -    &
               5.148746022615196e-10*zt**3*LOG(zrh)**2*LOG(zrhoa) +    &
               (0.0002686602132926594*LOG(zrh)**2*LOG(zrhoa))/x -   &
               0.06199739728812199*LOG(zrhoa)**2 +    &
               0.000906958053583576*zt*LOG(zrhoa)**2 -   &
               9.11727926129757e-7*zt**2*LOG(zrhoa)**2 -    &
               5.367963396508457e-9*zt**3*LOG(zrhoa)**2 -   &
               (0.007742343393937707*LOG(zrhoa)**2)/x +    &
               0.0121827103101659*LOG(zrh)*LOG(zrhoa)**2 -   &
               0.0001066499571188091*zt*LOG(zrh)*LOG(zrhoa)**2 +    &
               2.534598655067518e-7*zt**2*LOG(zrh)*LOG(zrhoa)**2 -    &
               3.635186504599571e-10*zt**3*LOG(zrh)*LOG(zrhoa)**2 +    &
               (0.0006100650851863252*LOG(zrh)*LOG(zrhoa)**2)/x +   &
               0.0003201836700403512*LOG(zrhoa)**3 -    &
               0.0000174761713262546*zt*LOG(zrhoa)**3 +   &
               6.065037668052182e-8*zt**2*LOG(zrhoa)**3 -    &
               1.421771723004557e-11*zt**3*LOG(zrhoa)**3 +   &
               (0.0001357509859501723*LOG(zrhoa)**3)/x
	       
        zntot=EXP(zntot)  !  add. Eq. (13)

          
        ! Equation (14) - radius of the critical cluster in nm

        zrc=EXP(-1.6524245+0.42316402*x+0.33466487*LOG(zntot))    ! [nm]

        ! Conversion [nm -> m]

        zrxc(jl)=zrc*1e-9       
	  
        !----1.2) Limiter

        IF(zjnuc<1.e-7 .OR. zntot<4.0) zjnuc=0.0

        ! limitation to 1E+10 [1/cm3s]
      
        zjnuc=MIN(zjnuc,1.e10)

        pxtrnucr(jl,jk) = zjnuc

        ! convert total number of molecules in the critical cluster
        ! to number of sulfate molecules:

        pntot(jl,jk)=zntot*zxmole

        ELSE ! pmolecH2SO4, ptp1 , prhd out of range

        pntot(jl,jk)   =0.0
        pxtrnucr(jl,jk)=0.0

        END IF

     END DO ! kproma
  END DO ! klev

END SUBROUTINE nucl_vehkamaeki


END MODULE mo_aero_nucl
