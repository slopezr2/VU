SUBROUTINE m7_equil (kproma, kbdim,  klev,   prelhum, paerml, paernl, &
                     pm6rp,  pm6dry, phplus, pww,     prhop,  ptp1    )
  !Modified by AM, relative humidity not larger than 1.00000 because of round-off errors in zmo(3)
  !
  !**** *m7_equil* calculates the ambient radii of accumulation 
  !                and coarse mode particles in presence of sea salt
  !     
  !
  !  Authors:
  !  --------
  !  E. Vignati, JRC/EI (original source)                    01/2000
  !  P. Stier, MPI      (f90-version, changes, comments)        2001
  !
  !  Purpose:
  !  --------
  !  This routine calculates the equilibrium radius of sea salt
  !  particles, and of sea salt particles mixed with sulphate for 
  !  accumulation and coarse modes. 
  !
  !  Method:
  !  -------
  !  from the mass of sea salt and sulphate (in ug/m3),
  !   
  !
  !**Interface:
  !  ----------
  !  *m7_equil* is called from *m7*
  !
  !  Externals:
  !  ----------
  !  none
  !
  !  References:
  !  -----------
  !  Jacobson, M.Z., Tabazadeh, A., Turco, R.P., (1996). Simulating
  !     equilibrium within aerosols and nonequilibrium between gases
  !     and aerosols. 
  !  Tang, I.N., (1997). Thermodynamic and optical properties of 
  !     mixed-salt aerosols of atmospheric importance. 
  !
  !@@@  ToDo:  Rewrite with identifiers of ions and ion pairs!!!!!
  !
  !--- Local variables: to be completed

  USE mo_aero_m7, ONLY: naermod, nsol,    nss,     nmod,      &
                        wnacl,   wna2so4, wnahso4, wh2so4,    &
                        ddust,   dbc,     doc,     dnacl,     &
                        dna2so4, dnahso4, dh2so4,  dh2o,      &
                        ram2cmr, pi,      r_kcal,  crh,       &
                        avo

  IMPLICIT NONE
  !
  !--- Parameter list:
  !
  ! phplus = molality of H+ [mol/kg(water)]
  !
  !
  !--- Local variables:
  ! 
  !  zna(:,:,i) = sodium concentration [ug m-3]
  !               i=1: accumulation mode; i=2: coarse mode
  !  zcl(:,:,i) = chlorine concentration [ug m-3]
  !               i=1: accumulation mode; i=2: coarse mode
  !
  ! !@@@ To be completed!
  !
  ! Indices for the arrays zmm and zmmr:
  !  -----------------------------------
  ! |     ions               ion pair   |
  ! |  (mole m-3)            (mole m-3) |
  ! | i   zmm(i)             zmmr(i)    |
  ! | 1   Na+                NaCl       |
  ! | 2   Cl-                NaHSO4     |
  ! | 3   SO4--              Na2SO4     |
  ! | 4   HSO4-              H2-SO4     |
  ! | 5   H+                            |
  !  -----------------------------------

  !--- Parameters:

  INTEGER :: kproma, kbdim, klev

  REAL    :: prelhum(kbdim,klev),        ptp1(kbdim,klev)
  
  REAL    :: paerml(kbdim,klev,naermod), paernl(kbdim,klev,nmod), &
             pm6rp(kbdim,klev,nmod),     pm6dry(kbdim,klev,nsol), &
             prhop(kbdim,klev,nmod),     pww(kbdim,klev,nmod)

  REAL    :: phplus(kbdim,klev,nss)
             
  !--- Local variables:

  INTEGER :: i,           jl,            jk,           jmod

  REAL    :: zaw,         zvolw,         zdryvol,      zdrymass,  &
             zkw,         zdryvol_mean,  zambvol_mean

  REAL    :: zmm(5),      zmmr(5),       zmol(5),      zmo(4),    &
             zmmt(4)
  !   
  !
  !print*, 'inequil',kproma, kbdim,  klev,   prelhum, paerml, paernl, &
  !                   pm6rp,  pm6dry, phplus, pww,     prhop,  ptp1 
  DO 50 jmod=1,nss
     DO 60 jk=1,klev
        DO 70 jl=1,kproma

           IF ((paerml(jl,jk,jmod+12) > 1.E-15) .AND. (paernl(jl,jk,jmod+2)>1.E-10))  THEN
              !
              !--- 1) Dry Calculations: --------------------------------------------------
              !
              !--- 1.1) Calculate initial concentrations of the compounds in mole/m+3:
              !
              !--- Na, Cl:
              !
              zmm(1)=paerml(jl,jk,jmod+12)*1.E-6 / wnacl       ! n(Na)    [mole m-3]
              zmm(2)=paerml(jl,jk,jmod+12)*1.E-6 / wnacl       ! n(Cl)    [mole m-3]
              !      [         g m-3           ] / [g mole-1]           = [mole m-3]
              !
              !--- SO4--:
              !
              zmm(3)=paerml(jl,jk,jmod+2) *1.E+6 / avo         ! n(H2SO4) [mole m-3]
              !      [      m-3                ] / [mole-1]             = [mole m-3]
              !
              !--- HSO4-:
              !
              zmm(4)=0.                                        ! n(HSO4)  [mole m-3]
              !
              !--- 1.2) Calculation of the concentration of the different species:
              !         The ions are supposed to be arranged such that
              !         sodium is associated to sulphate first in
              !         Na2SO4, the remaining sodium is associated to Cl
              !         in form of NaCl.
              !
              zmmt(1)=zmm(1)                        ! n(Na)  
              zmmt(3)=zmm(3)                        ! n(SO4)    
              !
              zmmr(3)=MIN(zmmt(1)/2. , zmmt(3))     ! n(Na2SO4) 
              zmmt(1)=zmmt(1)-2*zmmr(3)             ! Remaining n(Na) after association
                                                    ! with Na2SO4: n(Na)=n(Na)-2*n(Na2SO4)
              !
              zmmr(1)=MIN(zmm(2),zmmt(1))           ! n(NaCl) 
              !
              zmm(2)=zmmr(1)                        ! n(Cl) bound in NaCl, the rest is
                                                    ! assumed to evaporate in presence of SO4
              !
              zmmr(2)=0.                            ! n(NaHSO4)
              !
              zmmr(4)=zmm(3)-zmmr(2)-zmmr(3)        ! n(H2-SO4)(t)=n(H2SO4)(t0)-n(NaHSO4)-n(Na2SO4)
              !                                     ! as n(H2SO4)(t0)=n(SO4--)(t0)
              !
              !--- 1.3) Total aerosol dry volume [cm3/m3]:
              !
              zdryvol= zmmr(1)*wnacl/dnacl               + &
                       zmmr(2)*wnahso4/dnahso4           + &
                       zmmr(3)*wna2so4/dna2so4           + &
                       zmmr(4)*wh2so4/dh2so4             + &
                       paerml(jl,jk,jmod+14)*1.e-6/ddust + &
                       paerml(jl,jk,jmod+5) *1.e-6/dbc   + &
                       paerml(jl,jk,jmod+9) *1.e-6/doc
              !
              !--- 1.4) Mean aerosol dry volume [cm+3]:
              zdryvol_mean = zdryvol    / (paernl(jl,jk,jmod+2)*1.E6)
              ! [cm+3]     = [cm+3/m+3] /  [         m-3            ]
              !
              !--- 1.5) Dry radius [cm]:
              !
              pm6dry(jl,jk,jmod+2)=((3./(4.*pi))*zdryvol_mean)**(1./3.) * ram2cmr(jmod+2)
              !
              !--- 1.6) Total aerosol dry mass [gr/m3]:
              !
              zdrymass= zmmr(1)*wnacl                    + &
                        zmmr(2)*wnahso4                  + &
                        zmmr(3)*wna2so4                  + &
                        zmmr(4)*wh2so4                   + &
                        paerml(jl,jk,jmod+14)*1.e-6      + & ! Dust
                        paerml(jl,jk,jmod+5)*1.e-6       + & ! Black Carbon
                        paerml(jl,jk,jmod+9)*1.e-6           ! Organic Carbon
           !    print *,'zmmr', zmmr(:)
            !    print *,'zmmt', zmmt(:)
            !    print *, 'prelhum', prelhum(jl,jk)
                        
              !
              !
              !--- 2) Wet calculations: --------------------------------------------------
              !
              !--- Set threshold for relative humidity:
              !    If RH is smaller than the Critical Relative Humidity 
              !    (currently crh=0.45) the equilibrium radius is set to the dry radius:

              IF (prelhum(jl,jk) < crh) THEN
              
                
                 !
                 phplus(jl,jk,jmod)  = 0. 
                 !
                 pww(jl,jk,jmod+2)     = 0.
                 !
                 pm6rp(jl,jk,jmod+2) = pm6dry(jl,jk,jmod+2)
                 !
                 prhop(jl,jk,jmod+2) = zdrymass/zdryvol
                 
                 
                 !
              ELSE
              
              
              
              
                
                 !
                 !--- 2.1) Calculate thermodynamic properties under ambient conditions
                 !
                 !--- 2.1.1) Water activity:
                 !
                 zaw=min(prelhum(jl,jk),1.00000)
                 
                 !
                 !--- 2.1.2) Molality as function of the water activity:
                 !           Currently sulfate is assumed to be fully dissociated,
                 !           i.e. zmmr(2)=0. and zmo(2) is not calculated.
                 !
                 !           Changed reference to Jacobson et al. (1996):
                 !
                 !--- NaCl:
                 
                 zmo(1)=(-1.918004E2+2.001540E3*zaw-8.557205E3*zaw**2           &
                         +1.987670E4*zaw**3-2.717192E4*zaw**4+2.187103E4*zaw**5 &
                         -9.591577E3*zaw**6+1.763672E3*zaw**7                    )**2 

                 !--- NaHSO4:

                 zmo(2)=(+4.662777E0-1.128472E1*zaw+7.049464E1*zaw**2           &
                         -2.788050E2*zaw**3+6.103105E2*zaw**4-7.409417E2*zaw**5 & 
                         +4.614577E2*zaw**6-1.150735E2*zaw**7                    )**2 

                 !--- Na2SO4:

                 zmo(3)=(-3.295311E3+3.188349E4*zaw-1.305168E5*zaw**2           &
                         +2.935608E5*zaw**3-3.920423E5*zaw**4+3.109519E5*zaw**5 &
                         -1.356439E5*zaw**6+2.510249E4*zaw**7                    )**2

                 !--- H2-SO4:
                 
                 zmo(4)=(+5.611895-1.387446E1*zaw+1.750682E1*zaw**2             &
                         +7.138146E1*zaw**3-3.109173E2*zaw**4+4.662288E2*zaw**5 &
                         -3.128612E2*zaw**6+7.76097E1*zaw**7                     )**2
                 !
                 !
                 !--- 2.2) Calculation of the water content in kg water/m3 air:
                 !         (zmmr[mole/m3(air)]/zmo[mole/kg(water)] =zww[kg(water)/m3(air)]
                 !
                 pww(jl,jk,jmod+2)=zmmr(1)/zmo(1)+zmmr(2)/zmo(2)+zmmr(3)/zmo(3)+zmmr(4)/zmo(4)
                 !print *, 'zaw, zmo, pww',zaw, zmo(:), pww
                 !
                 !--- 2.3) Calculate the molality of the ions
                 !
                 !--- 2.3.1) For Na+, Cl-, SO4--, HSO4- :
                 !
                 DO i=1,4
                    zmol(i)=zmm(i)/pww(jl,jk,jmod+2)
                 END DO
                 !
                 !--- 2.3.2) For h+ :
                 !
                 !    [H+] = -[Na+] + [Cl-] + [OH-] + 2[SO4--] +[HSO4-]
                 !
                 !    with [OH-] = kw/[H+]
                 !
                 !--- Calculate autodissociation constant (kw) for water:
                 !    (Seinfeld &Pandis (1998): Eq. (6.5) + Table 6.5)

                 zkw=1.0E-14*exp( (13.35/r_kcal) * (1./298. - 1./ptp1(jl,jk)) )

                 !--- Calculate molality of H+:

                 zmol(5)=( (-zmol(1)+zmol(2)+2.*zmol(3)+zmol(4)) +                          &
                            SQRT( (zmol(1)-zmol(2)-2.*zmol(3)-zmol(4))**2 + 4.*zkw ) ) / 2.
                 !
                 !
                 !print *,'zmol, pww,zkw', zmol(:),pww(:,:,:),zkw
                 zmm(5)=zmol(5)*pww(jl,jk,jmod+2)
                 !
                 phplus(jl,jk,jmod)=zmol(5)
                 !
                 !
                 !--- 2.4) Calculation of the wet radius
                 !
                 !--- 2.4.1) Total water volume  [cm3/m3]:
                 !
                 zvolw=pww(jl,jk,jmod+2)/(dh2o*1.E-3)  ![cm3/m3]=[kg/m3]/([g/cm3]*[1.E-3 kg/g])
                 !
                 !
                 !--- 2.4.2) Mean aerosol ambient volume:
                 zambvol_mean = (zdryvol+zvolw) / (paernl(jl,jk,jmod+2)*1.E6)
                 ! [cm+3]     = [  cm+3/m+3   ] / [         m-3             ]                 
                 
                 !--- 2.4.3) Equilibrium wet count mean radius [cm]:
                 !
                 pm6rp(jl,jk,jmod+2)=((3./(4.*pi))*zambvol_mean)**(1./3.) * ram2cmr(jmod+2)
                 !
                 !--- 2.4.4) Calculation of the particle density (g cm-3):
                 !
                 prhop(jl,jk,jmod+2)=(zdrymass+zvolw*dh2o)/(zvolw+zdryvol)
                 !

              END IF !(prelhum(jl,jk) < crh)
             !print *, 'zmm,zmo',zmm(:),zmo(:)
           
           END IF !((paerml(jl,jk,jmod+12) > 1.E-15) .AND. (paernl(jl,jk,jmod+2)>1.E-10)) 

70      END DO
60   END DO
50 END DO
  !
 ! print*, 'outequil',kproma, kbdim,  klev,   prelhum, paerml, paernl, &
!                     pm6rp,  pm6dry, phplus, pww,     prhop,  ptp1 
END SUBROUTINE m7_equil


