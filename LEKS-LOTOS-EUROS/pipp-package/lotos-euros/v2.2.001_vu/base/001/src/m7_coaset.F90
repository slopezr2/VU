SUBROUTINE m7_coaset(kproma, kbdim, klev,  paernl, ptp1, &
                     papp1,  pm6rp, prhop, pcom          )
  !
  ! *m7_coaset*  calculates the coagulation kernels between the modes
  !
  ! Authors:
  ! ---------
  ! J. Wilson  and E. Vignati, JRC/EI (original source)                09/2000
  ! P. Stier, MPI                     (f90-version, changes, comments)    2001 
  !
  ! Modifications:
  ! --------------
  ! Philip Stier, MPI                             2001
  !
  ! Purpose
  ! ---------
  ! This routine calculates the coaglation kernels between particles
  ! with the count median radii of the three modes.
  ! Coagulation allowed between the following modes:
  ! soluble modes:   1+1=1, 2+2=2, 1+2=2, 1+3=3, 1+4=4, 2+3=3, 2+4=4
  ! insoluble modes: 2i+2i=2i
  ! mixed modes:     1+2i=2, 1+4i=4, 2+2i=2, 2+4i=4, 3+2i=3, 4+2i=4
  !
  ! Interface:
  ! -----------
  ! *m7_coaset* is called from *m7_dnum*
  !
  ! Externals:
  ! -----------
  ! none
  !
  ! Reference:
  ! -----------
  ! The calculations are based on:
  ! Fuchs, N.A. (1964). The Mechanics of Aerosols. Pergamon Press. Oxford. 
  ! (Chapter VII, Section 49)
  !
  !  Warning:
  !  --------
  !  For optimization purposes currently only "physically reasonable" elements of the
  !  coagulation kernel pcom are calculated in m7_concoag. These elements are specified
  !  in the matrix locoagmask in mo_aero_m7. Check carefully and adapt locoagmask 
  !  accordingly  before changes in the code below.
  !
  USE mo_aero_m7, ONLY: pi, bk, sqrt2, nmod, locoagmask
  !
  IMPLICIT NONE
  !
  !--- Parameter list:
  !
  !  paernl            = aerosol number for each mode [cm-3]
  !  ptp1              = atmospheric temperature at time t+1 [K]
  !  papp1             = atmospheric pressure at time t+1 [Pa]
  !  pm6rp             = mean mode actual radius (wet radius for soluble modes 
  !                      and dry radius for insoluble modes) [cm]
  !  prhop             = mean mode particle density [g cm-3]
  !  pcom(:,:,jm1,jm2) = Coagulation coefficient for modes jm1 and jm2 []
  !
  !--- List of local variables:
  !
  ! zwlc              = mean free pathlength []
  ! zairvisc          = air viscosity []
  ! zrpav             = average radius of the two interacting modes
  ! zpvx              = volume of the xth interacting mode
  ! zpmx              = mass of the xth interacting mode
  !
  ! zrknudx           = knudsen number of the xth interacting mode
  ! zpd2x             = particle diffusion of the xth interacting mode
  !
  !--- Parameters:
  !
  INTEGER :: kproma, kbdim, klev

  REAL    :: ptp1(kbdim,klev),          papp1(kbdim,klev)

  REAL    :: pm6rp(kbdim,klev,nmod),    prhop(kbdim,klev,nmod),                 &
             paernl(kbdim,klev,nmod)

  REAL    :: pcom(kbdim,klev,nmod,nmod)

  !--- Local variables:
  
  INTEGER :: jm2, jm1, jl, jk
  !
  REAL    :: zpbyone,     zwlc,        zairvisc,    zbtketc,   &
             zrpvm1,      zrpvm2,      zrpav,       zpv1,      &
             zpm1,        zcv21,       zpv2,        zpm2,      &
             zcv22,       zcv2av,      zrknud1,     zrknud2,   &
             ze1,         ze2,         zpd21,	    zpd22,     &
             zpdav,       zxd1,        zh21,        zxd2,      &
             zh22,        zh2av,       zcoc,        zhu1,      &
             zhu2,        zhu

  !--- 1) Calculation of the coagulation coefficient: ---------------------------
  !
  DO jm2=1,nmod
     DO jm1=jm2,nmod

        IF (locoagmask(jm1,jm2)) THEN

        DO jk=1,klev
           DO jl=1,kproma
	      IF (paernl(jl,jk,jm1) > 1.E-10 .AND. paernl(jl,jk,jm2) > 1.E-10 .AND. &
                  pm6rp(jl,jk,jm1)  > 1.E-10 .AND. pm6rp(jl,jk,jm2)  > 1.E-10 .AND. &
                  prhop(jl,jk,jm1)  > 1.E-10 .AND. prhop(jl,jk,jm2)  > 1.E-10        ) THEN        
           
                 !--- 1.1) Calculate ambient properties:

                 !--- Mean free pathlength ? (from Knudsen Number below):
                 !    Parametrisation?
                 zpbyone=1000.0 / (papp1(jl,jk)/100.0)
                 zwlc=6.6e-6 * ptp1(jl,jk) / 293.15 * zpbyone

                 !--- Viscosity:
                 zairvisc=1.827e-4 * (ptp1(jl,jk) / 291.15)**0.74

                 !--- 
                 zbtketc=bk * ptp1(jl,jk) / 6.0 / pi / zairvisc

                 !--- Count median radii of the respective modes:
                 zrpvm1=pm6rp(jl,jk,jm1)
                 zrpvm2=pm6rp(jl,jk,jm2)

                 !--- Average radius of the modes:
                 zrpav=(zrpvm1 + zrpvm2) / 2.0

                 !--- Volume and mass of mode jm1:
                 zpv1=4.0 / 3.0 * pi * zrpvm1**3.0
                 zpm1=zpv1 * prhop(jl,jk,jm1)

                 !--- Squared mean particle velocity of mode jm1:
                 zcv21=8.0 * bk * ptp1(jl,jk) / (pi * zpm1)

                 !--- Volume and mass of particles in mode jm2:
                 zpv2=4.0 / 3.0 * pi * zrpvm2**3.0
                 zpm2=zpv2 * prhop(jl,jk,jm2)

                 !--- Squared mean particle velocity of mode jm2:
                 zcv22=8.0 * bk * ptp1(jl,jk) / (pi * zpm2)

                 !--- Fuchs: G_r (below Eq. 49.27):
                 zcv2av=SQRT(zcv21 + zcv22)

                 !--- Knudsen numbers of the modes:
                 !@@@ Check: is zwlc mean free path then zrknud=zwlc/zrpvm!
                 zrknud1=zwlc/zrpvm1/2.0
                 zrknud2=zwlc/zrpvm2/2.0

                 !--- Diffusivities of the respective modes:
                 !    !@@@ Parameterisation?
                 ze1=EXP(-0.43/zrknud1)
                 ze2=EXP(-0.43/zrknud2)
                 zpd21=zbtketc * (1.0 + zrknud1*2.492 + zrknud1*0.84*ze1) / zrpvm1
                 zpd22=zbtketc * (1.0 + zrknud2*2.492 + zrknud2*0.84*ze2) / zrpvm2

                 !--- Average diffusivity of the modes:
                 zpdav=(zpd21 + zpd22) / 2.0

                 !--- Average mean free path of particles in jm1:
                 zxd1=8.0 * zpd21 / (pi*SQRT(zcv21))

                 !--- Mean distance from surface after mean free path (Eq. 49.13):
                 zh21=(((2.0*zrpvm1 + zxd1)**3.0 -                        &
                       (4.0*zrpvm1*zrpvm1 + zxd1*zxd1)**1.5) /            &
                       (6.0*zrpvm1*zxd1) - 2*zrpvm1           ) * sqrt2

                 !--- Average mean free path of particles in jm2:
                 zxd2=8.0 * zpd22 / (pi*SQRT(zcv22))

                 !--- Mean distance from surface after mean free path (Eq. 49.13):

                 zh22=(((2.0*zrpvm2 + zxd2)**3.0 -                        &
                       (4.0*zrpvm2*zrpvm2 + zxd2*zxd2)**1.5) /            &
                       (6.0*zrpvm2*zxd2) - 2*zrpvm2           ) * sqrt2

                 !--- Fuchs: delta_r !@@@ (why division by sqrt2?)
                 zh2av=SQRT(zh21*zh21 + zh22*zh22) / sqrt2
                 
                 !--- 1.2) Calculation of the coagulation coefficient pcom (Eq. 49.26):
                 !         Factor 16 instead factor 8 as in Fuchs as his formulation
                 !         applies for the inter-modal coagulation. This is taken into
                 !         account in the assignment of the inter-modal coagulation
                 !         coefficient.

                 zcoc=16.0 * pi * zpdav * zrpav

                 !--- Calculation of beta=1/zhu (Eq. 49.27):
                 zhu1=4.0 * zpdav / (zcv2av * zrpav)
                 zhu2=zrpav / (zrpav + zh2av /2.0)
                 zhu=zhu1 +  zhu2

                 !--- Coagulation coefficient following (Eq.49.26):
                 pcom(jl,jk,jm1,jm2)=zcoc / zhu

	      ELSE
                 pcom(jl,jk,jm1,jm2)=0.
	      END IF

              !--- 2) Mirror the coagulation matrix (symmetric): -----------------

              pcom(jl,jk,jm2,jm1)=pcom(jl,jk,jm1,jm2)

           END DO
        END DO

        ELSE

           pcom(1:kproma,:,jm1,jm2)=0.

        END IF ! locoagmask

     END DO
  END DO

END SUBROUTINE m7_coaset
