SUBROUTINE m7_delcoa(kproma,   kbdim,  klev,     paerml,   &
                     paernl,   pm6rp,  pa4delt,  panew,    &
                     pa,       pb,     pbfract1, pbfract2, &
                     pbfract5, pso4_5, pso4_6,   pso4_7, ptime    )
  ! 
  !    Authors: 
  !    --------- 
  !    E. Vignati and J. Wilson, JRC/EI (original source)                09/2000
  !    P. Stier, MPI                    (f90-version, changes, comments)    2001 
  ! 
  !    Version/History: 
  !    ----------------
  !    equivalent to the version delco_n2 of the m7 boxmodel 
  !    + use of the analytical solution 
  ! 
  !    Purpose 
  !    --------- 
  !    This routine calculates changes in number concentration of 
  !    each aerosol mode over the time step, due to coagulation with 
  !    the current mode and all higher ones. 
  ! 
  !    Method: 
  !    ----------- 
  !    *delcoa*  integrates for each mode dn/dt=c -a*n^2 -b*n  over ztmst  
  ! 
  !    The resulting particles are assumed to reside in the 
  !    mode of highest mode of the pair of particles colliding. 
  !    1+1=>1 1+2=>2 1+3=>3, 2+2=>2 2+3=>3, 3+3=>3. 
  !    zc is now non zero for mode 1 only (nucleation). 
  !    All formation of higher mode particles is handled in dconc. 
  ! 
  !    For climatological studies, 5 day accumulation mode concs are  
  !    within a factor of 2 of the full model.  
  ! 
  !    Interface 
  !    ----------- 
  !    *m7_delcoa* is called from *m7_dnum* 
  ! 
  !    Externals 
  !    ----------- 
  !    none 
  ! 
 
!  USE mo_time_control,  ONLY: time_step_len
  USE mo_aero_m7,       ONLY: nmod,    nsol,    naermod,             & 
                              iaiti,   ibcki,   iocki,   ibcks,      & 
                              iocks,   ibcas,   iocas,   ibccs,      &
                              ioccs,   iduci,   iduai,   iacci,      &
                              icoai,   iaits

  USE mo_aero_m7,        ONLY: m7_coat

                                  
  IMPLICIT NONE 

  !--- Parameter list:
  !
  ! paerml          = total aerosol mass for each compound 
  !                   [molec. cm-3 for sulphate; ug m-3 for others]
  ! paernl          = aerosol number for each mode [cm-3]
  ! pm6rp           = mean mode actual radius (wet radius for soluble modes 
  !                   and dry radius for insoluble modes) [cm]
  ! pa4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one timstep 
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in m7_nuck)
  !                      - coagulation (calculated in m7_concoag)
  ! panew           = number of nucleated particles (during 1 timestep) [1] 
  ! pa              = unimodal coagulation coefficient (zcom(mod)) []
  ! pb              = inter-modal coagulation with higher modes
  !                   (zcom(mod) * n(jmod+1)) []
  ! pbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x that is moved to mode y+1 [1]
  ! pso4_x          = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  ! 
  !--- Local variables:
  !
  ! zansum          = aerosol number in the respective mode [cm-3]
  ! zxxsum          = aerosol mass for compound xx in the respective
  !                   mode, e.g. xx = bc, oc, a4 (sulfate) 
  !                   [g cm-3 for bc,oc and molec. cm-3 for sulfate]
  ! zxxav           = average mass of a sulfate particle in the respective
  !                   mode [molecules]
  ! zxxavy          = average mass of species xx in mode y []
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  !                   [molecules for sulfate and ug for others]
  ! zanli(:,:,:)    = Number of particles moved by the inter-modal 
  !                   coagulation []
  ! zansq(:,:,:)    = Number of particles moved by the intra-modal 
  !                   coagulation []
  ! zaernt(:,:,:)   = New particle number n(t+dt) after the integration
  !                   of the aerosol dynamics equation [cm-3]

  !--- Parameters:

  INTEGER :: kproma, kbdim, klev

  REAL    :: ptime, time_step_len
 
  REAL    :: paerml(kbdim,klev,naermod),   paernl(kbdim,klev,nmod),         & 
             pa4delt(kbdim,klev,naermod),  panew(kbdim,klev),               &
             pm6rp(kbdim,klev,nmod)
 
  REAL    :: pbfract1(kbdim,klev,nmod-1),  pbfract2(kbdim,klev,nmod-1),     & 
             pbfract5(kbdim,klev,3)
 
  REAL    :: pa(kbdim,klev,nmod),          pb(kbdim,klev,nmod)
 
  REAL    :: pso4_5(kbdim,klev),           pso4_6(kbdim,klev),              & 
             pso4_7(kbdim,klev)
 
  ! Local variables: 
 
  INTEGER :: jmod, jl, jk, kmod
 
  REAL    :: za4av1(kbdim,klev),           za4av2(kbdim,klev),              & 
             zbcav2(kbdim,klev),           zocav2(kbdim,klev),              &
             zbcav3(kbdim,klev),           zocav3(kbdim,klev),              & 
             zbcav4(kbdim,klev),           zocav4(kbdim,klev),              & 
             zbcav5(kbdim,klev),           zocav5(kbdim,klev),              & 
             zduav6(kbdim,klev),           zduav7(kbdim,klev) 
 
  REAL    :: zanli(kbdim,klev,nmod),       zansq(kbdim,klev,nmod),          &
             zaernt(kbdim,klev,nmod) 

  REAL    :: zm6rp(nmod),                  zcrtcst(nmod)

  ! Auxiliary variables: 
 
  REAL    :: zansum, zbcsum, zocsum, za4sum, za4av,                         & 
             ztop,   zbot,   zatot,  zanse,  zanle,                         & 
             ze1,    zf1,    zf2,    zf3,    zf4,    zr1,                   &
             ztmst,  zc

  REAL    :: zamloss5,      zamloss6,      zamloss7,      zanloss5,         &
             zanloss6,      zanloss7,      ztloss,        zbfofac,   zms,   &
             zbfnfac,       zbftot,        zaerni,        zanytnt,          &
             zanytni,       zanytns,       zanytnm,       ztotn,            &
             zaerns

  !--- m7_box:
  REAL    :: zaernl_diff(kbdim,klev,nmod), zaerml_diff(kbdim,klev,naermod), &
             zaernl_0(kbdim,klev,nmod),    zaerml_0(kbdim,klev,naermod)

  REAL    :: zcrit_5, zcrit_6, zcrit_7


  
  zcrit_5=0.
  zcrit_6=0.
  zcrit_7=0.
  !--- END m7_box

  !--- 0) Initialisations: ------------------------------------------------ 
 
  time_step_len = ptime

  ztmst  = time_step_len 
 
  za4av       = 0. 
  za4av1(:,:) = 0.
  za4av2(:,:) = 0.
  zbcav2(:,:) = 0. 
  zocav2(:,:) = 0. 
  zbcav3(:,:) = 0.
  zocav3(:,:) = 0.
  zbcav4(:,:) = 0.
  zocav4(:,:) = 0.
  zbcav5(:,:) = 0.
  zocav5(:,:) = 0.
 
  !--- m7_box:
  zaernl_0=paernl
  zaerml_0=paerml
  !--- END m7_box

  !--- 1) Insoluble modes 
  ! 
 
  DO jk=1,klev 
     DO jl=1,kproma 

        !--- Dust modes insoluble: 
        !    Calculate average mass (used in concoag)

        IF(paernl(jl,jk,iacci) > 0.0) THEN
           zduav6(jl,jk) = (paerml(jl,jk,iduai)*1.E-6) / paernl(jl,jk,iacci)
        ELSE
           zduav6(jl,jk)=0.
        END IF
        IF(paernl(jl,jk,icoai) > 0.0) THEN
           zduav7(jl,jk) = (paerml(jl,jk,iduci)*1.E-6) / paernl(jl,jk,icoai)
        ELSE
           zduav7(jl,jk)=0. 
        END IF

        !--- Aitken mode insoluble:
        !    Only considered process: 
        !    Coagulation and transfer from the insoluble aitken mode 
        !    to the soluble aitken and accumulation modes.

        zansum=paernl(jl,jk,iaiti) 
        zbcsum=paerml(jl,jk,ibcki)*1.e-6 
        zocsum=paerml(jl,jk,iocki)*1.e-6 
        zaernt(jl,jk,iaiti)=zansum 
        zanli(jl,jk,iaiti)=0.0 
        zansq(jl,jk,iaiti)=0.0 

        !--- Calculations only in presence of sufficient particles:

        IF (zansum > 1.E-10) THEN 
 
           ! --- Average mass for a bc/oc particle in the aitken mode [ug]:

           zbcav5(jl,jk)=zbcsum/zansum 
           zocav5(jl,jk)=zocsum/zansum 

           !--- 1.1) Case of no coagulation:
           !
           !         (pa(jl,jk,iaiti) < 1.e-15 .AND. pb(jl,jk,iaiti) < 1.e-15) 
           !
           !         => Nothing to be done
              
           !--- 1.2) Case with coagulation:

           IF (pa(jl,jk,iaiti) >= 1.e-15 .OR. pb(jl,jk,iaiti) >= 1.e-15) THEN
                 
              !--- 1.2.1) Case of no inter-modal coagulation:
              !           dn/dt = -a*n**2 => 
              !           n(t)  = n0/(1 + n0*a*(t-t0))

              IF (pb(jl,jk,iaiti) < 1.e-15) THEN 
                 zaernt(jl,jk,iaiti)=zansum/(1.0+zansum*pa(jl,jk,iaiti)*ztmst) 
                 zanli(jl,jk,iaiti)=0.0 
                 zansq(jl,jk,iaiti)=zansum-zaernt(jl,jk,iaiti) 

              !--- 1.2.2) Case with inter- and intra-modal coagulation:
              !           dn/dt = -a*n**2 - b*n => 
              !           n(t)  = (b*n0*exp(-b(t-t0)))/((n0*a)(1-exp(-b(t-t0)))+b)
              ELSE   
                 !--- Calculate n(t+dt):
                 ze1=EXP(-pb(jl,jk,iaiti)*ztmst) 
                 ztop=pb(jl,jk,iaiti)*zansum*ze1 
                 zbot=zansum*pa(jl,jk,iaiti)*(1.0-ze1)+pb(jl,jk,iaiti) 
                 zaernt(jl,jk,iaiti)=ztop/zbot 
                 !--- Limit n(t+dt) to available particle in the mode:
                 zaernt(jl,jk,iaiti)=MIN(zaernt(jl,jk,iaiti), zansum) 
                 !--- Total change in particle numbers of the mode due to coagulation:
                 zatot=zansum-zaernt(jl,jk,iaiti) 
                 !--- Contribution of the intra-modal coagulation:
                 zanse=zansum*zansum*pa(jl,jk,iaiti)
                 !--- Contribution of the inter-modal coagulation:
                 zanle=zansum*pb(jl,jk,iaiti) 
                 !--- Number of particles moved by the inter-modal coagulation:
                 zanli(jl,jk,iaiti)=zatot*zanle/(zanse+zanle) 
                 !--- Number of particles moved by the intra-modal coagulation:
                 zansq(jl,jk,iaiti)=zatot*zanse/(zanse+zanle) 
              END IF 

              !--- 1.2.3) Change masses of the insoluble aitken mode due to 
              !           intra-modal coagulation and the coagulation with the
              !           nucleation mode (transfers to the soluble modes
              !           of the particles coagulating with the nucleation mode 
              !           are done in m7_concoag):
 
              paerml(jl,jk,ibcki)=( zaernt(jl,jk,iaiti) +                    &
                                    zansq(jl,jk,iaiti)  +                    &
                                    pbfract5(jl,jk,1)*zanli(jl,jk,iaiti) ) * & 
                                    zbcav5(jl,jk)*1.e6 

              paerml(jl,jk,iocki)=( zaernt(jl,jk,iaiti) +                    &
                                    zansq(jl,jk,iaiti)  +                    &
                                    pbfract5(jl,jk,1)*zanli(jl,jk,iaiti) ) * & 
                                    zocav5(jl,jk)*1.e6 
              
              !--- 1.2.4) Change the numbers of the insoluble aitken mode due to 
              !           intra-modal coagulation:
 
              paernl(jl,jk,iaiti)=zaernt(jl,jk,iaiti)   +                    &
                                  pbfract5(jl,jk,1)*zanli(jl,jk,iaiti)
 
              !--- 1.2.5) Store changes in masses of compounds in the insoluble 
              !           aitken mode due to inter-modal coagulation:
              !           (zanli(:,:,x)   = total number of particles moved from mode x
              !            pbfract5(:,:,x)= fraction of the total number of particles
              !                             moved from mode 5 that is moved to mode x   )
 
              pa4delt(jl,jk,ibcks)=pbfract5(jl,jk,2)*zanli(jl,jk,iaiti)*zbcav5(jl,jk)*1.e6 
              pa4delt(jl,jk,iocks)=pbfract5(jl,jk,2)*zanli(jl,jk,iaiti)*zocav5(jl,jk)*1.e6 
              pa4delt(jl,jk,ibcas)=pbfract5(jl,jk,3)*zanli(jl,jk,iaiti)*zbcav5(jl,jk)*1.e6 
              pa4delt(jl,jk,iocas)=pbfract5(jl,jk,3)*zanli(jl,jk,iaiti)*zocav5(jl,jk)*1.e6 
 
           END IF
        END IF
     END DO
  END DO
  ! 
  !--- 2) Soluble modes: --------------------------------------------------
  ! 
!CDIR unroll=5
   DO jmod=1,nsol !mode :
     DO jk=1,klev !level
         DO jl=1,kproma !longitude

           !--- Nucleation mode:
           IF (jmod .EQ. 1) THEN 
              zansum=paernl(jl,jk,jmod)+panew(jl,jk) 
              za4sum=paerml(jl,jk,jmod)+pa4delt(jl,jk,1)
           !--- Others:
           ELSE 
              zansum=paernl(jl,jk,jmod) 
              za4sum=paerml(jl,jk,jmod)
           END IF 

           IF (jmod.EQ.2) THEN 
	      zbcsum=paerml(jl,jk,ibcks)*1.e-6 
	      zocsum=paerml(jl,jk,iocks)*1.e-6 
           END IF 
           IF (jmod.EQ.3) THEN 
	      zbcsum=paerml(jl,jk,ibcas)*1.e-6 
	      zocsum=paerml(jl,jk,iocas)*1.e-6 
           END IF 
           IF (jmod.EQ.4) THEN 
	      zbcsum=paerml(jl,jk,ibccs)*1.e-6 
	      zocsum=paerml(jl,jk,ioccs)*1.e-6 
           END IF 

           zaernt(jl,jk,jmod)=zansum 
           zanli(jl,jk,jmod)=0.0 
           zansq(jl,jk,jmod)=0.0 

           !--- Calculations only in presence of sufficient particles:

           IF (zansum > 1.E-10) THEN 

              za4av=za4sum/zansum 

              IF (jmod.EQ.1) THEN 
                 za4av1(jl,jk)=za4av 
              ELSE IF (jmod.EQ.2) THEN 
                 zbcav2(jl,jk)=zbcsum/zansum 
                 zocav2(jl,jk)=zocsum/zansum 
                 za4av2(jl,jk)=za4av 
              ELSE IF (jmod.EQ.3) THEN 
                 zbcav3(jl,jk)=zbcsum/zansum 
                 zocav3(jl,jk)=zocsum/zansum                  
              ELSE IF (jmod.EQ.4) THEN 
                 zbcav4(jl,jk)=zbcsum/zansum 
                 zocav4(jl,jk)=zocsum/zansum 
              END IF

              !--- 2.1) Case of no coagulation:
              !         
              IF (pa(jl,jk,jmod) < 1.e-15 .AND. pb(jl,jk,jmod) < 1e-15) THEN 
 
                 !--- Nucleation in mode 1 only. 
                 !    Nothing to be done for other modes.

                 IF(jmod.EQ.1) THEN 
                    paerml(jl,jk,jmod)=za4sum
                    paernl(jl,jk,jmod)=zansum
                 END IF 
 
              !--- 2.2) Case with coagulation:

              ELSE 

                 !--- 2.2.1) Case of no nucleation:

                 !--- Not Mode 1 or Nucleation rate below 1/s:

                 IF ( (jmod .NE. 1) .OR. (panew(jl,jk)/ztmst < 1.0) ) THEN

                    paernl(jl,jk,jmod)=zansum 

                    !--- 2.2.1a) Case of no inter-modal coagulation:
                    !            dn/dt = -a*n**2 => 
                    !            n(t)  = n0/(1 + n0*a*(t-t0))

                    IF (pb(jl,jk,jmod) < 1.e-15) THEN 
                       zaernt(jl,jk,jmod)=zansum/(1.0+zansum*pa(jl,jk,jmod)*ztmst) 
                       zanli(jl,jk,jmod)=0.0 
                       zansq(jl,jk,jmod)=zansum-zaernt(jl,jk,jmod) 

                    !--- 2.2.1b) Case with inter- and intra-modal coagulation:
                    !            dn/dt = -a*n**2 - b*n => 
                    !            n(t)  = (b*n0*exp(-b(t-t0)))/((n0*a)(1-exp(-b(t-t0)))+b)

                    ELSE            
                       !--- Calculate n(t+dt):
                       ze1=EXP(-pb(jl,jk,jmod)*ztmst) 
                       ztop=pb(jl,jk,jmod)*zansum*ze1 
                       zbot=zansum*pa(jl,jk,jmod)*(1.0-ze1)+pb(jl,jk,jmod) 
                       zaernt(jl,jk,jmod)=ztop/zbot 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum) 
                       !--- Total change in particle numbers of the mode due to coagulation:
                       zatot=zansum-zaernt(jl,jk,jmod) 
                       !--- Contribution of the intra-modal coagulation:
                       zanse=zansum*zansum*pa(jl,jk,jmod)
                       !--- Contribution of the inter-modal coagulation:
                       zanle=zansum*pb(jl,jk,jmod) 
                       !--- Number of particles moved by the inter-modal coagulation:
                       zanli(jl,jk,jmod)=zatot*zanle/(zanse+zanle) 
                       !--- Number of particles moved by the intra-modal coagulation:
                       zansq(jl,jk,jmod)=zatot*zanse/(zanse+zanle) 
                    END IF 
                 
                 !--- 2.2.2) Case with nucleation:

                 ELSE IF ( (jmod .EQ. 1) .AND. (panew(jl,jk)/ztmst >= 1.0) ) THEN

                    !--- 2.2.2a) Nucleation, inter- and intra-modal coagulation:
                    !            dn/dt = -a*n**2 - b*n + c => 
                    !            n(t)  = -(b/(2a)) + 
                    !                    R/2a * [ ((1 - (-2ax0-b+R)/(+2ax0+b+R))exp(-Rt)) /
                    !                             ((1 + (-2ax0-b+R)/(+2ax0+b+R))exp(-Rt))  ]
                    !            where:  R=SQRT(b**2+4ac)
                    !
                    !            If b/=0 then always a/=0. The only case where a would be 0
                    !            and b unequal zero is the case of no pre-existing particles 
                    !            in the nucleation mode but pre-existing particles in other
                    !            modes. For this case a is calculated for an assumed radius
                    !            of a critical cluster in m7_coaset. 

                    IF (pb(jl,jk,jmod) >= 1.E-15) THEN
                       !--- Calculate n(t):
                       !--- c:
                       zc=panew(jl,jk)/ztmst
                       !--- R:
                       zf1=pb(jl,jk,jmod)*pb(jl,jk,jmod)+4.0*pa(jl,jk,jmod)*zc 
                       zr1=SQRT(zf1) 
                       !--- exp(-Rt):
                       ze1=EXP(-zr1*ztmst) 
                       !--- 2ax0+b:
                       zf2=2.0*pa(jl,jk,jmod)*paernl(jl,jk,jmod)+pb(jl,jk,jmod) 
                       !--- Term in squared bracket:
                       zf3=ze1*(zr1-zf2)/(zr1+zf2) 
                       zf4=(1.0-zf3)/(1.0+zf3) 
                       !--- n(t):
                       zaernt(jl,jk,jmod)=(zr1*zf4-pb(jl,jk,jmod))/2.0/pa(jl,jk,jmod) 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum) 
                       !--- Total change in particle numbers of the mode due to coagulation:
                       zatot=zansum-zaernt(jl,jk,jmod) 
                       !--- Contribution of the intra-modal coagulation:
                       zanse=zansum*zansum*pa(jl,jk,jmod)
                       !--- Contribution of the inter-modal coagulation:
                       zanle=zansum*pb(jl,jk,jmod) 
                       !--- Number of particles moved by the inter-modal coagulation:
                       zanli(jl,jk,jmod)=zatot*zanle/(zanse+zanle) 
                       !--- Number of particles moved by the intra-modal coagulation:
                       zansq(jl,jk,jmod)=zatot*zanse/(zanse+zanle) 

                    !--- 2.2.2b) Nucleation and intra-modal coagulation:
                    !            dn/dt = -a*n**2 - b*n + c with b=0 =>
                    !            dn/dt = -a*n**2 + c => 
                    !            n(t)  = R/2a * [ ((1 - (-2ax0+R)/(+2ax0+R))exp(-Rt)) /
                    !                             ((1 + (-2ax0+R)/(+2ax0+R))exp(-Rt))  ]
                    !            where:  R=SQRT(4ac)
                    !
                    !            Can be shown to be equivalent to:
                    !
                    !            n(t)  = R1*((x0+R1)/(x0-R1)+exp(-SQRT(-4ac)t)) / 
                    !                       ((x0+R1)/(x0-R1)-exp(-SQRT(-4ac)t))
                    !            where R1=SQRT(c/a)

                    ELSE IF (pb(jl,jk,jmod) < 1.E-15) THEN 
                       !--- c:
                       zc=panew(jl,jk)/ztmst
                       !--- R1:
                       zr1=SQRT(zc/pa(jl,jk,jmod)) 
                       !--- exp(-Rt):
                       ze1=EXP(-zr1*2.0*pa(jl,jk,jmod)*ztmst)
                       !--- n(t):
                       zf1=(paernl(jl,jk,jmod)+zr1)/(paernl(jl,jk,jmod)-zr1) 
                       ztop=zr1*(zf1+ze1) 
                       zbot=zf1-ze1 
                       IF (zbot < 1.E-15) THEN 
                          zaernt(jl,jk,jmod)=zansum 
                       ELSE 
                          zaernt(jl,jk,jmod)=ztop/zbot 
                       END IF 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt(jl,jk,jmod)=MIN(zaernt(jl,jk,jmod), zansum) 
                       !--- Number of particles moved by the inter-modal coagulation:
                       zanli(jl,jk,jmod)=0.0 
                       !--- Number of particles moved by the intra-modal coagulation:
                       zansq(jl,jk,jmod)=zansum-zaernt(jl,jk,jmod) 
                    END IF 
                 END IF 
                 !---2.2.3 New bit for insoluble/souble coagulation
                 !--- sum total insoluble+soluble paticles in mode jmod JJNW
                 IF (jmod .EQ. 1 .AND. zanli(jl,jk,jmod)>0.0) THEN
                 zaerni=paernl(jl,jk,iaiti)+paernl(jl,jk,iacci)+paernl(jl,jk,icoai)
                 zaerns=zansum+paernl(jl,jk,iaits)
!                 zaerns=zansum
                 ztotn=zaerns+zaerni
                 IF (zaerns .gt. zaerni .and. zaerni .gt. 0.0) THEN
                    ! calculate analytical solution no of mixed particles for coagulation
                    ! between paernl(jl,jk,jmod) soluble particles and zaerni insouble of
                    ! the same dimensions
                    IF (zaerni .gt. 1.0) then
                       zanytni=4.0*zaerni/((2.0+pa(jl,jk,jmod)*ztmst*ztotn)*                     &
                              (2.0+pa(jl,jk,jmod)*ztmst*(ztotn-zaerni)))
                    ELSE
                      zanytni = 0.0
                    ENDIF
                    zanytnt=2.0*ztotn/(2.0+pa(jl,jk,jmod)*ztmst*ztotn)
                    zanytns=4.0*zaerns/((2.0+pa(jl,jk,jmod)*ztmst*ztotn)*            &
                            (2.0+pa(jl,jk,jmod)*ztmst*(ztotn-zaerns)))
                    zanytnm=zanytnt-(zanytni+zanytns)
                    zanytnm=min(zanytnm,zaerni)
                    zanytni=zaerni-zanytnm
                    zanytns=zaernt(jl,jk,jmod)
!CDIR UNROLL=7
                    DO kmod=1,nmod
                       zm6rp(kmod)=pm6rp(jl,jk,kmod)
                    END DO
                    CALL m7_coat(zm6rp,zcrtcst)
                    zamloss5=paernl(jl,jk,5)/zaerni*zanytnm*zcrtcst(5)
                    zanloss5=zamloss5/za4av
                    zamloss6=paernl(jl,jk,6)/zaerni*zanytnm*zcrtcst(6)
                    zanloss6=zamloss6/za4av
                    zamloss7=paernl(jl,jk,7)/zaerni*zanytnm*zcrtcst(7)
                    zanloss7=zamloss7/za4av
                    ztloss=zanloss5+zanloss6+zanloss7
                    zms=zansq(jl,jk,jmod)*0.95

                    ztloss=min(ztloss,zansq(jl,jk,jmod)*0.95)                               
                    zbfofac=zanli(jl,jk,jmod)/(zanli(jl,jk,jmod)+ztloss)
                    zbfnfac=ztloss/(zanli(jl,jk,jmod)+ztloss)
                    zanli(jl,jk,jmod)=zanli(jl,jk,jmod)+ztloss
                    zansq(jl,jk,jmod)=zansq(jl,jk,jmod)-ztloss
                    zbftot=0.0
!CDIR UNROLL=7
                    DO kmod=1,nmod
                       IF(kmod>jmod) THEN
                          pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)*zbfofac
                          IF (kmod.GE.5) THEN
                             pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)+                  &
                                                       zbfnfac*paernl(jl,jk,kmod)/zaerni
                          END IF
                          zbftot=zbftot+pbfract1(jl,jk,kmod-jmod)
                       END IF
                    END DO       
!CDIR UNROLL=7
                    DO kmod=1,nmod
                       IF (kmod>jmod) THEN
                          pbfract1(jl,jk,kmod-jmod)=pbfract1(jl,jk,kmod-jmod)/zbftot
                       END IF
                    END DO
                 END IF
                 END IF
                 !---- End of new inslouble/soluble caogulation routine JJNW
                 !--- 2.3) Change masses and numbers of the respective modes to account-----------
                 !         for intra-modal coagulation (zansq) and coagulation with
                 !         higher modes (zaernt):
                 !
                 !--- 2.3.1) Change mass of the sulfur compounds:

                 paerml(jl,jk,jmod)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*za4av 
 
                 !--- 2.3.2) Change mass of the carbon compounds:

                 IF (jmod.EQ.2) THEN 
                    paerml(jl,jk,ibcks)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*zbcav2(jl,jk)*1.e6 
                    paerml(jl,jk,iocks)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*zocav2(jl,jk)*1.e6 
                 ELSE IF (jmod.EQ.3) THEN 
                    paerml(jl,jk,ibcas)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*zbcav3(jl,jk)*1.e6 
                    paerml(jl,jk,iocas)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*zocav3(jl,jk)*1.e6 
                 ELSE IF (jmod.EQ.4) THEN 
                    paerml(jl,jk,ibccs)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*zbcav4(jl,jk)*1.e6 
                    paerml(jl,jk,ioccs)=(zaernt(jl,jk,jmod)+zansq(jl,jk,jmod))*zocav4(jl,jk)*1.e6 
                 END IF 

                 !--- 2.3.3) Particle numbers:
 
                 paernl(jl,jk,jmod)=zaernt(jl,jk,jmod)!@@@+zansq(jl,jk,jmod)/2.0 
 
                 !--- 2.4) Calculate changes in particle masses due to inter-modal --------------
                 !         coagulation:

                 !--- 2.4.1) Transfer of mass from mode 1 to other modes:

                 IF (jmod .EQ. 1) THEN 
                     
                    ! Mass from 1 to 2: 
 
                    pa4delt(jl,jk,2)=pbfract1(jl,jk,1)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 2 due to coag. with 5:
 
                    pa4delt(jl,jk,2)=pa4delt(jl,jk,2)+pbfract1(jl,jk,4)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 3:
 
                    pa4delt(jl,jk,3)=pbfract1(jl,jk,2)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 3 due to coag. with 6:
 
                    pa4delt(jl,jk,3)=pa4delt(jl,jk,3)+pbfract1(jl,jk,5)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 4: 
 
                    pa4delt(jl,jk,4)=pbfract1(jl,jk,3)*zanli(jl,jk,1)*za4av 
 
                    ! Mass from 1 to 4 due to coag. with 7:
 
                    pa4delt(jl,jk,4)=pa4delt(jl,jk,4)+pbfract1(jl,jk,6)*zanli(jl,jk,1)*za4av 

                 !---  2.4.2) Transfer of mass from mode 2 to other modes:
 
                 ELSE IF (jmod .EQ. 2) THEN 
 
                    ! Mass from 2 to 3: 
                     
                    pa4delt(jl,jk,3)=pa4delt(jl,jk,3)+pbfract2(jl,jk,2)*zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibcas)=pa4delt(jl,jk,ibcas)+                                & 
                                         pbfract2(jl,jk,2)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.e6 
                    pa4delt(jl,jk,iocas)=pa4delt(jl,jk,iocas)+                                & 
                                         pbfract2(jl,jk,2)*zanli(jl,jk,2)*zocav2(jl,jk)*1.e6 
 
                    ! Mass from 2 to 3 due to coag. with 6: 
 
                    pa4delt(jl,jk,3)=pa4delt(jl,jk,3)+pbfract2(jl,jk,5)*zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibcas)=pa4delt(jl,jk,ibcas)+                                & 
                                         pbfract2(jl,jk,5)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.e6 
                    pa4delt(jl,jk,iocas)=pa4delt(jl,jk,iocas)+                                & 
                                         pbfract2(jl,jk,5)*zanli(jl,jk,2)*zocav2(jl,jk)*1.e6 

                    ! Mass from 2 to 4: 
 
                    pa4delt(jl,jk,4)=pa4delt(jl,jk,4)+pbfract2(jl,jk,3)*zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibccs)=pa4delt(jl,jk,ibccs)+                                & 
                                         pbfract2(jl,jk,3)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.E6 
                    pa4delt(jl,jk,ioccs)=pa4delt(jl,jk,ioccs)+                                & 
                                         pbfract2(jl,jk,3)*zanli(jl,jk,2)*zocav2(jl,jk)*1.E6 

                    ! Mass from 2 to 4 due to coag. with 7: 
 
                    pa4delt(jl,jk,4)=pa4delt(jl,jk,4)+pbfract2(jl,jk,6)*zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibccs)=pa4delt(jl,jk,ibccs)+                                & 
                                         pbfract2(jl,jk,6)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.E6 
                    pa4delt(jl,jk,ioccs)=pa4delt(jl,jk,ioccs)+                                & 
                                         pbfract2(jl,jk,6)*zanli(jl,jk,2)*zocav2(jl,jk)*1.E6 

                    ! Mass from 2 due to coagulation of 2 with 5 remains in 2:
                    !
                    !@@@ No effect as pbfract2(:,:,4)=0.!
                    !@@@ (Not needed as it is assumed that 5 coagulates with 2 
                    !@@@ and therefor the masses in 2 remain unchanged!)

                    pa4delt(jl,jk,2)=pa4delt(jl,jk,2)+pbfract2(jl,jk,4)*zanli(jl,jk,2)*za4av
                    pa4delt(jl,jk,ibcks)=pa4delt(jl,jk,ibcks)+                                &
                                         pbfract2(jl,jk,4)*zanli(jl,jk,2)*zbcav2(jl,jk)*1.E6
                    pa4delt(jl,jk,iocks)=pa4delt(jl,jk,iocks)+                                &
                                         pbfract2(jl,jk,4)*zanli(jl,jk,2)*zocav2(jl,jk)*1.E6
 
                 END IF
              END IF
           END IF
        END DO !longitude 
     END DO !level 
  END DO !mode 

  !--- 3) Calculate transfer from the insoluble to the soluble modes: -------------------------
  CALL m7_concoag (kproma,   kbdim,   klev,                       &
                   paerml,   paernl,  pm6rp,  pa4delt, zanli,     & 
                   za4av1,   za4av2,  zbcav5, zocav5,  zduav6,    & 
                   zduav7,   pso4_5,  pso4_6, pso4_7,             &
                   pbfract1, pbfract2,                            &
                   zcrit_5,  zcrit_6,  zcrit_7                    ) !--- m7_box: Added for diagnostics
 
 
  !--- 4) Final change of the aerosol masses due to nucleation, -------------------------------
  !       inter-modal coagulation and condensation on the insoluble modes:
  !       (Nucleation mode already done above.)

  DO jmod=2,naermod 
     paerml(1:kproma,:,jmod)=paerml(1:kproma,:,jmod)+pa4delt(1:kproma,:,jmod) 
  END DO 

  !--- m7_box: Include diagnostics for conservation of numbers:

  DO jmod=1,nmod 
     DO jk=1,klev 
        DO jl=1,kproma

           SELECT CASE(jmod)
              CASE(1)
                 zaernl_diff(jl,jk,jmod) = zaernl_0(jl,jk,jmod) -                   &
                                          (paernl(jl,jk,jmod) - panew(jl,jk) +      &
                                           zansq(jl,jk,jmod)  + zanli(jl,jk,jmod))
              CASE(2)
                 zaernl_diff(jl,jk,jmod) = zaernl_0(jl,jk,jmod) -                   &
                                          (paernl(jl,jk,jmod) +                     &
                                           zansq(jl,jk,jmod)  + zanli(jl,jk,jmod) - &
                                           zcrit_5)
              CASE(3)
                 zaernl_diff(jl,jk,jmod) = zaernl_0(jl,jk,jmod) -                   &
                                          (paernl(jl,jk,jmod) +                     &
                                           zansq(jl,jk,jmod)  + zanli(jl,jk,jmod) - &
                                           zcrit_6)
              CASE(4)
                 zaernl_diff(jl,jk,jmod) = zaernl_0(jl,jk,jmod) -                   &
                                          (paernl(jl,jk,jmod) +                     &
                                           zansq(jl,jk,jmod)  + zanli(jl,jk,jmod) - &
                                           zcrit_7)
              CASE (5)
                 zaernl_diff(jl,jk,jmod) = zaernl_0(jl,jk,jmod) -                   &
                                          (paernl(jl,jk,jmod) +                     &
                                           zansq(jl,jk,jmod)  + zanli(jl,jk,jmod) + &
                                           zcrit_5)
              CASE (6)
                 zansq(jl,jk,jmod)=0.
                 zanli(jl,jk,jmod)=0.
                 zaernl_diff(jl,jk,jmod) = zaernl_0(jl,jk,jmod) -                   &
                                          (paernl(jl,jk,jmod) +                     &
                                           zansq(jl,jk,jmod)  + zanli(jl,jk,jmod) + &
                                           zcrit_6)
              CASE (7)
                 zansq(jl,jk,jmod)=0.
                 zanli(jl,jk,jmod)=0.
                 zaernl_diff(jl,jk,jmod) = zaernl_0(jl,jk,jmod) -                   &
                                          (paernl(jl,jk,jmod) +                     &
                                           zansq(jl,jk,jmod)  + zanli(jl,jk,jmod) + &
                                           zcrit_7)
              END SELECT
                 
              zaerml_diff(jl,jk,jmod) = zaerml_0(jl,jk,jmod) - &
                                       (paerml(jl,jk,jmod)-panew(jl,jk)+zansq(jl,jk,jmod)+zanli(jl,jk,jmod))
        END DO
     END DO
  END DO

  
   
END SUBROUTINE m7_delcoa 
