SUBROUTINE m7_averageproperties(kproma, kbdim, klev, paernl, paerml, pttn, pm6rp, prhop)
  !
  !  Author:
  !  --------
  !  E. Vignati, JRC/EI (original source)                10/2000
  !  P. Stier, MPI      (f90-version, changes, comments)    2001
  !
  !  Purpose:                                                           
  !  ---------                                                           
  !  Calculation of the mean particle mass (pttn).
  !     [molecules cm-3] for the sulphate mass
  !     [ug m-3] for the other compounds
  !
  !  Calculation of the (dry) radius and the density 
  !  of the particles of the insoluble modes.
  !
  !  Interface:
  !  ----------
  !  m7_averageproperties is called from m7
  !
  !  Externals:
  !  ----------
  !  none
  !
  USE mo_aero_m7, ONLY: dbc,    doc,     ddust,  pi,      &  
                        critn,  ram2cmr, nmod,   naermod, &
                        ibcks,  ibcas,   ibccs,  ibcki,   &
                        iocks,  iocas,   ioccs,  iocki,   &
                        issas,  isscs,                    &
                        iduas,  iducs,   iduai,  iduci,   &
                        iaiti,  iacci,   icoai,           &
                        nsol

  IMPLICIT NONE
  ! 
  !--- Parameter list:
  !
  ! paerml(kbdim,klev,naermod) = total aerosol mass for each compound 
  !                             [molec. cm-3 for sulfate and ug m-3 for others]
  ! paernl(kbdim,klev,nmod)    = aerosol number for each mode [cm-3]
  ! pttn(kbdim,klev,naermod)   = average mass for single compound in each mode 
  !                             [in molec. for sulphate and in ug for others]
  ! pm6rp(kbdim,klev,nmod)     = mean mode actual radius (wet radius for soluble
  !                             modes and dry radius for insoluble modes) [cm]
  ! prhop(kbdim,klev,nmod)     = mean mode particle density [g cm-3]
  !
  !--- Local variables:
  !
  ! zinsvol                   = average volume for single particle in the 
  !                             insolulbe mode [cm3]
  ! zinsmas                   = average mass for single particle in the 
  !                             insolulbe mode [g]

  !--- Parameters:

  INTEGER :: kproma, kbdim, klev

  REAL    :: paerml(kbdim,klev,naermod), paernl(kbdim,klev,nmod), & 
             pttn(kbdim,klev,naermod),   pm6rp(kbdim,klev,nmod),  & 
             prhop(kbdim,klev,nmod)
  
  !--- Local variables:

  INTEGER :: jmod,         jk,          jl
  
  REAL    :: zinsvol,      zinsmas,     zeps

  !--- 0) Initialization:

!  zeps=EPSILON(1.0)
   zeps=1.e-20
   
  !--- 1) Calculate mean particle masses at start of timestep: ---------------------------
  !
  !       To be able to compute a intra-modal coagulation coefficient for the nucleation
  !       mode for the case of no pre-existing particles but coagulation of freshly formed
  !       particles during the timestep, pttn is set to the mass of the critical cluster
  !       for this case. This allows to calculate an ambient radius of the 
  !       freshly formed particles and subsequently the calculation of the coagulation 
  !       coefficient. This mass is "virtual" as it is not added to the mode but used 
  !       only for the described computation of the coagulation coefficient. 
  !       !@@@ Check whether this is always fulfilled. 
  
  DO jmod=1,nsol
     DO jk=1,klev
        DO jl=1,kproma
           IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,jmod) .GT. 1e-15) THEN

              pttn(jl,jk,jmod)=paerml(jl,jk,jmod)/paernl(jl,jk,jmod)

           ELSE IF (jmod == 1 .AND. paernl(jl,jk,jmod) <= 1e-10 .AND. paerml(jl,jk,jmod) <= 1e-15) THEN

              pttn(jl,jk,jmod)=critn

           END IF
        END DO
     END DO
  END DO
  !
  !--- 3) Calculation of the mean mass pttn [ug] for each compound in the modes: ---------
  !       [Factor 1.E-6 to convert(ug m-3)/cm-3 into ug]
  !
  DO jmod=2,nmod
     DO jk=1,klev
        DO jl=1,kproma
           IF (jmod.EQ.2) THEN
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,ibcks) .GT. 1e-15) THEN
                 pttn(jl,jk,ibcks)=paerml(jl,jk,ibcks)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ibcks)=0.
              END IF
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,iocks) .GT. 1e-15) THEN
                 pttn(jl,jk,iocks)=paerml(jl,jk,iocks)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iocks)=0.
              END IF
           END IF
           !	    
           IF (jmod.EQ.3) THEN
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,ibcas) .GT. 1e-15) THEN
                 pttn(jl,jk,ibcas)=paerml(jl,jk,ibcas)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ibcas)=0.
              END IF
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,iocas) .GT. 1e-15) THEN
                 pttn(jl,jk,iocas)=paerml(jl,jk,iocas)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iocas)=0.
              END IF
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,issas) .GT. 1e-15) THEN
                 pttn(jl,jk,issas)=paerml(jl,jk,issas)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,issas)=0.
              END IF
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,iduas) .GT. 1e-15) THEN
                 pttn(jl,jk,iduas)=paerml(jl,jk,iduas)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iduas)=0.
              END IF
           END IF
           !	    
           IF (jmod.EQ.4) THEN
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,ibccs) .GT. 1e-15) THEN
                 pttn(jl,jk,ibccs)=paerml(jl,jk,ibccs)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ibccs)=0.
              END IF
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,ioccs) .GT. 1e-15) THEN
                 pttn(jl,jk,ioccs)=paerml(jl,jk,ioccs)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ioccs)=0.
              END IF
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,isscs) .GT. 1e-15) THEN
                 pttn(jl,jk,isscs)=paerml(jl,jk,isscs)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,isscs)=0.
              END IF
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,iducs) .GT. 1e-15) THEN
                 pttn(jl,jk,iducs)=paerml(jl,jk,iducs)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iducs)=0.
              END IF
           END IF
           !	    
           IF (jmod.EQ.5) THEN
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,ibcki) .GT. 1e-15) THEN
                 pttn(jl,jk,ibcki)=paerml(jl,jk,ibcki)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,ibcki)=0.
              END IF
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,iocki) .GT. 1e-15) THEN
                 pttn(jl,jk,iocki)=paerml(jl,jk,iocki)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iocki)=0.
              END IF
           END IF
           !	    
           IF (jmod.EQ.6) THEN
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,iduai) .GT. 1e-15) THEN
                 pttn(jl,jk,iduai)=paerml(jl,jk,iduai)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iduai)=0.
              END IF
           END IF
           !
           IF (jmod.EQ.7) THEN
              IF (paernl(jl,jk,jmod) .GT. 1e-10 .AND. paerml(jl,jk,iduci) .GT. 1e-15) THEN
                 pttn(jl,jk,iduci)=paerml(jl,jk,iduci)/paernl(jl,jk,jmod)*1.E-6
              ELSE
                 pttn(jl,jk,iduci)=0.
              END IF
           END IF
        END DO
     END DO
  END DO

  !
  !--- 4) Calculate count median radii for lognormal distribution from -------------------
  !       mass for insoluble modes:

  DO jk=1,klev
     DO jl=1,kproma

        !--- 4.1) Aitken mode insoluble:

        zinsmas=1.e-6*(pttn(jl,jk,ibcki)+pttn(jl,jk,iocki))
        zinsvol=1.e-6*(pttn(jl,jk,ibcki)/dbc+pttn(jl,jk,iocki)/doc)
        IF (zinsvol > zeps) THEN
           prhop(jl,jk,iaiti)=zinsmas/zinsvol
           pm6rp(jl,jk,iaiti)=(0.75/pi*1.e-6*                                              &
                              (pttn(jl,jk,ibcki)/dbc+pttn(jl,jk,iocki)/doc))**(1./3.)*ram2cmr(iaiti)
        ELSE
           prhop(jl,jk,iaiti)=0.
           pm6rp(jl,jk,iaiti)=0.
        END IF
      
        !--- 4.2) Accumulation mode insoluble:

        IF (pttn(jl,jk,iduai) > zeps) THEN 
           prhop(jl,jk,iacci)=ddust
           pm6rp(jl,jk,iacci)=(0.75/pi*1.e-6*pttn(jl,jk,iduai)/ddust)**(1./3.)*ram2cmr(iacci)
        ELSE
           prhop(jl,jk,iacci)=0.
           pm6rp(jl,jk,iacci)=0.  
        END IF

        !--- 4.3) Coarse mode insoluble:

        IF (pttn(jl,jk,iduci) > zeps) THEN
           prhop(jl,jk,icoai)=ddust
           pm6rp(jl,jk,icoai)=(0.75/pi*1.e-6*pttn(jl,jk,iduci)/ddust)**(1./3.)*ram2cmr(icoai)
        ELSE
           prhop(jl,jk,icoai)=0.
           pm6rp(jl,jk,icoai)=0.
        END IF

     END DO
  END DO

!  write(*,*) 'averprop', 'zinsvol= ', zinsvol, 'zeps= ', zeps, 'rad=', pm6rp(2100,1,5), 'massbc= ', pttn(2100,1,8), 'ram2cmr= ', ram2cmr(5)      

END SUBROUTINE m7_averageproperties
