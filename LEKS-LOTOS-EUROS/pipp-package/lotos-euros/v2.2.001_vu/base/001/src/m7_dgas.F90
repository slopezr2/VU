SUBROUTINE m7_dgas(kproma, kbdim, klev,  pso4g,  paerml, paernl, &
                   ptp1,   papp1, pm6rp, pso4_5, pso4_6, pso4_7, &
                   ptime                                           )
  !
  !**** *m7_dgas*  calculates the transfer of mass due to 
  !                sulfate condensation on soluble and 
  !                insoluble modes.
  !
  !    Authors:
  !    -----------
  !    J. Wilson, E. Vignati, JRC/EI (original source)                05/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)    2001   
  !
  !    Purpose:
  !    -----------
  !    This routine calculates the changes in aerosol mass and
  !    gas phase sulfate due to sulfate condensation.
  !
  !**  Interface:
  !    -----------
  !    *m7_dgas* is called from *m7*
  !
  !    Method:
  !    -----------------
  !    The transfer of sulfate to the particles is based on
  !    Fuchs (1959). Soluble and insoluble particles are distinguished
  !    by a different accomodation coefficient "caccso4" defined
  !    in mo_aero_m7.
  !    (Currently 1.0 for soluble and 0.3 for insoluble modes). 
  !
  !    Externals
  !    -----------
  !    none
  !
  !    References:
  !    -----------
  !    Fuchs, N.A. (1959). Evaporation and droplet growth in gaseous media;
  !       Pergamon, New York, pp72.

 ! USE mo_time_control,  ONLY: delta_time
  USE mo_control,       ONLY: nrow
  USE mo_aero_m7,       ONLY: pi,      wh2so4,  rerg,    avo,            &
                              fmax,    caccso4,                          &
                              nmod,    naermod,                          &
                              iso4ns,  iso4ks,  iso4as,  iso4cs
!  USE mo_aero_mem,      ONLY: d_cond_so4

  IMPLICIT NONE
  !
  !--- Parameter list:
  !
  ! pso4g     = mass of gas phase sulfate [molec. cm-3]
  ! pm6rp     = mean mode actual radius (wet radius for soluble modes 
  !             and dry radius for insoluble modes) [cm]
  ! pso4_x    = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  ! 
  !--- Local Variables:
  !
  ! zde       = molecular diffusion []
  ! zvelb     = velocity []
  ! zcondo    = condensation coefficient []
  ! zc2(nmod) = flux of sulfate condensing on the respective mode 
  !             per sulfate gas phase concentration []
  ! zcondo    = total flux of condensing sulfate 
  !             per sulfate gas phase concentration []
  ! zfcond    = total mass of condensing sulfate for one timestep []

  INTEGER :: kproma, kbdim, klev

  REAL    :: ptime, time_step_len

  REAL    :: ptp1(kbdim,klev),        papp1(kbdim,klev),                 &
             pso4g(kbdim,klev),       pso4_5(kbdim,klev),                &
             pso4_6(kbdim,klev),      pso4_7(kbdim,klev) 
  
  REAL    :: paernl(kbdim,klev,nmod), paerml(kbdim,klev,naermod),        &
             pm6rp(kbdim,klev,nmod)
  !
  ! Local variables:

  INTEGER :: jl,         jk,         jmod,       jrow
  
  REAL    :: zfcond,     zftot,      zpbyone,    zde2,                   & 
             zvelb,      zxibc,      zm6rp,      zf1,                    &
             zqtmst

  REAL    :: zcondo(kbdim,klev)

  REAL    :: zc2(kbdim,klev,nmod)

  !--- 0) Initialisations: -------------------------------------------------
  !
  jrow=nrow(2)

  zcondo(:,:)=0.0
      
  zc2(:,:,:) = 0.0  

  time_step_len = ptime

  zqtmst=1/time_step_len
  
  !--- 1) Calculate condensation rate for cm diameter sulphate aerosols: ---
  !
  DO jmod=1,nmod
     DO jk=1,klev
        DO jl=1,kproma
           IF (pm6rp(jl,jk,jmod).GT.0.) THEN

              !--- Diffusion coefficient (Reference???):

              zpbyone=1000.0 / (papp1(jl,jk)/100.0)

              zde2=0.073 * zpbyone * (ptp1(jl,jk) / 298.15)**1.5  

              !--- Mean molecule velocity (Moore, 1962 (S+P equ. 8.2)):

              zvelb=SQRT(8.0 * rerg * ptp1(jl,jk) / pi / wh2so4)  

              !--- ???Fuchs???

              zxibc=8.0 * zde2 / pi / zvelb
              !
              ! Use count median radius:

              zm6rp=pm6rp(jl,jk,jmod)

              !--- Distance from particle up to which the kinetic regime applies:

              zf1=( (zm6rp + zxibc)**3.0 - (zm6rp**2.0 + zxibc**2.0)**1.5 ) / &
                  (3.0 * zm6rp * zxibc) - zm6rp

              !--- Diffusive flux to single particle surface:
              !    (Elisabetta's thesis: fraction in equ. 2.26)

              zc2(jl,jk,jmod)=(4.0 * pi * zde2 * zm6rp ) /                      &
                              ((4.0 * zde2) / (zvelb * zm6rp * caccso4(jmod)) +   &
                               (zm6rp/(zm6rp+zf1))                              )

              !--- Total diffusive flux to all particles in the respective mode:
              !    (per concentration of gas phase sulfate)

              zc2(jl,jk,jmod)=zc2(jl,jk,jmod) * paernl(jl,jk,jmod)

              !--- Total diffusive flux to all particles in all modes:
              !    (per concentration of gas phase sulfate)

              zcondo(jl,jk)=zcondo(jl,jk)+ zc2(jl,jk,jmod)  

           END IF
        END DO
     END DO
  END DO
  !
  !--- 2) Calculation of the new sulfate aerosol masses and of the ---------
  !       mass of sulfate condensing on the respective modes:
  !
  DO jk=1,klev
     DO jl=1,kproma
        IF(zcondo(jl,jk).GT.0..AND.pso4g(jl,jk).GT.1.e-10) THEN

           !--- Total diffusive flux to all particles in all modes:

           zfcond=zcondo(jl,jk)*pso4g(jl,jk)

           !--- Total mass of sulfate condensing during 1 timestep:

           zftot=zfcond*time_step_len
           
           !--- Limit condensing sulfate to 
           !    fmax[%] x (available gas-phase sulfate) :

           zfcond=MIN(zftot,(pso4g(jl,jk)*fmax)) 

           !--- Remaining gas phase sulfate:
           
           pso4g(jl,jk)=pso4g(jl,jk)-zfcond
           
           !--- Add mass to sulfate compounds:
           !    zc2(:,:,jmod)/zcondo = fraction of total sulfate flux
           !                           that condenses on mode jmod
           !    => (   "    )*zfcond = mass of sulfate condensing on
           !                           the respective mode

           paerml(jl,jk,iso4ns)=paerml(jl,jk,iso4ns)+ &
                                zc2(jl,jk,iso4ns)/zcondo(jl,jk)*zfcond
           paerml(jl,jk,iso4ks)=paerml(jl,jk,iso4ks)+ &
                                zc2(jl,jk,iso4ks)/zcondo(jl,jk)*zfcond
           paerml(jl,jk,iso4as)=paerml(jl,jk,iso4as)+ &
                                zc2(jl,jk,iso4as)/zcondo(jl,jk)*zfcond
           paerml(jl,jk,iso4cs)=paerml(jl,jk,iso4cs)+ &
                                zc2(jl,jk,iso4cs)/zcondo(jl,jk)*zfcond

           !--- Mass of sulphate condensing on the insoluble modes:
           !    (Transfer from insoluble to soluble modes 
           !     calculated in m7_concoag.)
           
           pso4_5(jl,jk)=zc2(jl,jk,5)/zcondo(jl,jk)*zfcond
           pso4_6(jl,jk)=zc2(jl,jk,6)/zcondo(jl,jk)*zfcond
           pso4_7(jl,jk)=zc2(jl,jk,7)/zcondo(jl,jk)*zfcond


        ELSE 

           pso4_5(jl,jk)=0.
           pso4_6(jl,jk)=0.
           pso4_7(jl,jk)=0.

        END IF
     END DO
  END DO

!  write(*,*) 'dgas', 'so4gnew', pso4g(2100,1), 'condso4', zcondo(2100,1)*pso4g(2100,1)*time_step_len, 'limit=', pso4g(2100,1)*fmax
!  write(*,*) 'cond5= ', pso4_5(2100,1), 'rad5= ', pm6rp(2100,1,5)

END SUBROUTINE m7_dgas
