SUBROUTINE m7_equiz(kproma,  kbdim, klev,   &
                    papp1,   pttn,  ptp1,   &
                    prelhum, pm6rp, pm6dry, &
                    prhop,   pww,   paernl  )
  !
  !   *m7_equiz*   calculates the ambient radii of the sulphate particles
  !
  !    Authors:
  !    --------
  !    J. Wilson, E. Vignati, JRC/EI (original source)                05/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)    2001
  !
  !    Purpose:
  !    --------
  !    This routine calculates the ambient radii for sulfate particles
  !    with mass of ttn molecules, converts them to count mean radii and
  !    stores them in the array with address pm6rp.
  !    It additionally calculates the ambient particle density.
  !
  !    Method:
  !    -------
  !    The calculations of the ambient particle properties are based on 
  !    parameterisations of the mass of sulfate and density derived 
  !    by Julian Wilson from a regression analysis of results of solving 
  !    the generalised Kelvin equation using (F. J. Zeleznik, J. Phys. Chem.
  !    Ref. Data 20, 1157, 1991), for an H2SO4-H2O mixture, in the 
  !    following parameter ranges: 

  !       1e2   < pttn    < 1E11  [molecules]
  !       0.2   < prelhum < 0.9   [1]
  !       240   < ptp1    < 330   [K]
  !       10000 < papp1   < 100000  [Pa]
  !
  !    Due to the limitations of the parametrisation, the ambient temperature
  !    is restricted to a minimum of 240 K within this subroutine. 
  !      
  !    Interface:
  !    ----------
  !    *m7_equiz* is called from *m7*
  !
  !    Externals:
  !    ----------
  !    none
  !
  !
  USE mo_aero_m7, ONLY: naermod, nmod,    nsol,           &
                        ibcks,   iocks,   issas,   ibcas, &
                        iocas,   iduas,   isscs,   ibccs, &
                        ioccs,   iducs,                   &
                        wvb,     gmb,     avo,     wh2so4,&
                        pi,      ram2cmr, dh2so4,  dh2o
  !
  ! pttn      = average mass for single compound in each mode 
  !             [in molec. for sulphate and in ug for bc, oc, ss, and dust]
  ! pm6rp     = count mean radius under ambient conditions [cm]
  ! pm6dry    = count mean radius under dry conditions [cm]
  ! paernl    = aerosol number for each mode [cm-3]
  ! pww       = aerosol water content for each mode [kg(water) m-3(air)]
  ! zwso4     = percentage by mass of sulfate in a H2O-H2SO4 particle 
  !             containing pttn molecules of sulfate under ambient conditions
  ! zvso4     = volume of pttn molecules of sulfate [cm3]
  ! zmso4     = mass of pttn molecules of sulfate [g]
  ! zdso4h2o  = density of sulfate-h2o fraction of a particle with average 
  !             mass [g.cm-3]
  ! zmso4h2o  = mass of sulfate-h2o fraction of a particle with average mass [g]
  ! zvso4h2o  = volume of sulfate-h2o fraction of a particle with average 
  !             mass [cm3]

  IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL    :: papp1(kbdim,klev),        ptp1(kbdim,klev),       &
             prelhum(kbdim,klev)

  REAL    :: pttn(kbdim,klev,naermod), prhop(kbdim,klev,nmod), &
             pm6dry(kbdim,klev,nsol),  pm6rp(kbdim,klev,nmod), &
             pww(kbdim,klev,nmod),     paernl(kbdim,klev,nmod)

  !--- Local variables:

  INTEGER :: jk, jl, jmod

  REAL    :: zaerelse,                                         &
             zwso4,       zvso4,       zmso4,                  &
             zvso4h2o,    zmso4h2o,    zdso4h2o,               &
             zapp1,       ztk,         zrh

  REAL    :: ztk2,        zln3,        zln32,                  &
             zlnm,        zss2,        zlnm2

!CDIR unroll=5
  DO 100 jmod=1,nsol
     DO 90 jk=1,klev
        DO 80 jl=1,kproma

           !--- 1) Determine mass of non sulfate compounds in a particle: ---------

           SELECT CASE (jmod)

              CASE (1)
                 zaerelse=0.
              CASE (2)
                 zaerelse = pttn(jl,jk,ibcks)+pttn(jl,jk,iocks)
              CASE (3)
                 zaerelse = pttn(jl,jk,issas)+pttn(jl,jk,ibcas)+  &
                            pttn(jl,jk,iocas)+pttn(jl,jk,iduas)
              CASE (4)
                 zaerelse = pttn(jl,jk,isscs)+pttn(jl,jk,ibccs)+  &
                            pttn(jl,jk,ioccs)+pttn(jl,jk,iducs)
           END SELECT

           !--- 2) Calculation of the particle properties in the absense of -------
           !       other compounds than sulfate:

           IF (pttn(jl,jk,jmod) > 0.0 .AND. zaerelse < 1.E-15) THEN

	      !
              !--- 2.1) Calculation of the ambient particle properties: -----------
              !      
              !--- Constrain ambient temperature to conditions for which the 
              !    parametrisation of the liquid water content works:
              
              ! Temperature:
              ztk = ptp1(jl,jk)
              ztk = MAX(ztk , 240.)

              ! Relative Humidity:
              zrh = prelhum(jl,jk)
              zrh = MAX(zrh , 0.05)
              zrh = MIN(zrh , 0.90)

              !--- Assign auxiliary variables:

              zapp1=papp1(jl,jk)
              zlnm = LOG(pttn(jl,jk,jmod))
              zlnm2 = zlnm*zlnm
              zss2 = zrh**2
              ztk2 = ztk*ztk
              zln3 = zlnm/3.0
              zln32 = zln3*zln3
              !
              !--- Percentage by weight of sulfate in the particle [%]:
              !    (Here we ignore any insoluble mass.)
              !
              zwso4 = wvb(1) + wvb(2)*zlnm + wvb(3)*zrh*zlnm + wvb(4)*ztk*zlnm + &
                      wvb(5)*zrh/ztk + wvb(6)*zlnm2*zrh + wvb(7)*zlnm2*ztk +     &
                      wvb(8)*zlnm*zss2 + wvb(9)*zlnm*ztk2 + wvb(10)*zlnm2*zss2 + &
                      wvb(11)*zlnm2*ztk2 + wvb(12)*zss2/ztk2 + wvb(13)*zlnm2 +   &
                      wvb(14)*zlnm2*zlnm + wvb(15)*zlnm2*zlnm2 +                 &
                      wvb(16)*zss2*zrh/(ztk2*ztk) + wvb(17)*LOG(zrh*ztk/zapp1)

              !--- Dry mass of sulfate in an average particle [g]:

              zmso4 = pttn(jl,jk,jmod)*wh2so4/avo

              !--- Dry volume of sulfate in an average particle [cm3]:
              !    Any temperature or pressure dependency of the 
              !    sulfate density is ingored.

              zvso4 = zmso4/dh2so4

              !--- Mass of sulfate + water in an average particle [g]:

              zmso4h2o = zmso4/(zwso4/100.0)

              !--- Density of the sulfate-water fraction of an average particle [g cm-3]:
              !@@@ Check: changed zwvso4 into zwso4 (now the mass!)

              zdso4h2o = gmb(1) + gmb(2)*zwso4 + gmb(3)*zln3 + gmb(4)*zrh +        &
                         gmb(5)*ztk + gmb(6)*zln32 + gmb(7)*zln3/zrh +             &
                         gmb(8)*zln3/ztk + gmb(9)*ztk2
         
              !--- Limits for zdso4h2o: H2O(0.99) and pure H2SO4 (1.841):
              !
              zdso4h2o=MAX(zdso4h2o,dh2o)
              zdso4h2o=MIN(zdso4h2o,dh2so4)
              
              !--- Volume of sulfate-water fraction of an average particle [cm3]:

              zvso4h2o = zmso4h2o/zdso4h2o

              !--- 2.2) Calculatiion of the particle radii: ----------------------------

              !--- 2.2.1) Dry count mean radius [cm]:
 
              pm6dry(jl,jk,jmod)=((zvso4)*0.75/pi)**(1./3.)*ram2cmr(jmod)

              !--- 2.2.2) Equilibrium wet count mean radius [cm]:
             
              pm6rp(jl,jk,jmod) =((zvso4h2o)*0.75/pi)**(1./3.)*ram2cmr(jmod)
              
              !--- 2.3) Assignment of the particle density [g cm-3]: -------------------
 
              prhop(jl,jk,jmod)=zdso4h2o

              !--- 2.4) Store aerosol water for each mode [kg(water) m-3(air)]:

              pww(jl,jk,jmod)=(zmso4h2o-zmso4)*paernl(jl,jk,jmod)*1.E3

           END IF
80      END DO
90   END DO
100 END DO

END SUBROUTINE m7_equiz
