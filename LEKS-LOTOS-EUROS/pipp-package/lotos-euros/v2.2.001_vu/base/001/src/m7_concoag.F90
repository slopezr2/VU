SUBROUTINE m7_concoag (kproma,   kbdim,    klev,                      &
                       paerml,   paernl,   pm6rp,   pa4delt, panli,   &
                       pa4av1,   pa4av2,   pbcav5,  pocav5,  pduav6,  &
                       pduav7,   pso4_5,   pso4_6,  pso4_7,           &
                       pbfract1, pbfract2,                            &
                       zcrit_5,  zcrit_6, zcrit_7                     ) !--- m7_box: Added for diagnostics 
                                                                        !            in m7_delcoa
  !
  !   *m7_concoag*
  !
  !   Author:
  !   ----------
  !   E. Vignati, JRC/EI     (original source)                09/2000
  !   P. Stier, MPI          (f90-version, changes, comments)    2001 

  !   Version:
  !   ----------
  !   This version is equivalent to the version concoa_n of the boxmodel. 
  !
  !   Purpose
  !   ----------
  !   m7_concoag transfers aerosol mass and numbers from the insoluble
  !   to the soluble modes.
  !
  !   Interface:
  !   ----------
  !   *m7_concoag* is called from *m7_delcoa*
  !
  !   Externals
  !   ----------
  !   none

  USE mo_aero_m7, ONLY: m7_coat, nmod,    naermod,         &
                        ibcks,   ibcki,   iocks,  iocki,   &
                        iduas,   iducs,   iduai,  iduci,   &
                        iaiti,   iacci,   icoai,           &
                        iaits,   iaccs,   icoas

  IMPLICIT NONE 

  !--- Parameters:
  !
  ! paerml          = total aerosol mass for each compound 
  !                   [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and dust]
  ! paernl          = aerosol number for each mode [cm-3]
  ! pm6rp           = mean mode actual radius (wet radius for soluble modes 
  !                   and dry radius for insoluble modes) [cm]
  ! pa4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one timstep 
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in m7_nuck)
  !                      - coagulation (calculated here in m7_concoag)
  ! pxxavy          = average mass of species xx in mode y []!@@@
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  ! panli(:,:,x)    = total number of particles moved by inter-modal 
  !                   coagulation from mode x [cm-3]
  ! pbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x that is moved to mode y+1 [1] 
  !                   !@@@ Clumsy notation! Should be moved to mode y !!!
  ! pso4_x          = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  !
  !--- Local variables / Constants:
  !
  ! zso4x    = available mass of sulfate from mode 1 and 2 
  !            condensing and coagulating on mode x (x = insoluble modes 5,6,7).  
  !
  ! zcrtcst  = Critical constant, i.e. number of sulfate molecules to cover 
  !            an average particle of the mode with a layer of the thickness
  !            determined by cLayerThickness in mo_aero_m7. Calculated by
  !            m7_coat.
  !
  ! =>         zso4x/zcrtcst is the total number of particles that could be moved
  !            from insoluble mode x to soluble modes.
  !
  ! zcrit_x  = total available number of particles in mode x that are moved from 
  !            insoluble mode x to the corresponding soluble mode.

  INTEGER :: kproma, kbdim, klev

  REAL    :: pso4_5(kbdim,klev),          pso4_6(kbdim,klev),                   &
             pso4_7(kbdim,klev),                                                &
             pa4av1(kbdim,klev),          pa4av2(kbdim,klev),                   &
             pbcav5(kbdim,klev),          pocav5(kbdim,klev),                   &
             pduav6(kbdim,klev),          pduav7(kbdim,klev)
              
  REAL    :: paerml(kbdim,klev,naermod),  paernl(kbdim,klev,nmod),              &
             pbfract1(kbdim,klev,nmod-1), pbfract2(kbdim,klev,nmod-1),          &
             panli(kbdim,klev,nmod),      pa4delt(kbdim,klev,naermod),          &
             pm6rp(kbdim,klev,nmod)


  ! Local variables:

  INTEGER :: jl,          jk,            jmod   

  REAL    :: zcrit_5,     zcrit_6,       zcrit_7,                               &
             zso45,       zso46,         zso47,                                 &
             zeps

  REAL    :: zm6rp(nmod), zcrtcst(nmod)

  !--- 0) Initializations:

  zeps=EPSILON(1.)

  !--- 1) Redistribution of mass and numbers after nucleation, coagulation ----
  !       and coagulation calculated in the preceeding subroutines:

  DO jk=1,klev
     DO jl=1,kproma

        !--- 1.1) Sum mass of sulphate added to modes 5, 6, and 7 due to 
        !         coagulation with modes 1 and 2 (1st term) and the mass
        !         of sulfate condensed on the insoluble mode x (pso4_x):
        
        zso45=panli(jl,jk,1)*pbfract1(jl,jk,4)*pa4av1(jl,jk)+pso4_5(jl,jk)

        zso46=panli(jl,jk,1)*pbfract1(jl,jk,5)*pa4av1(jl,jk)+                &
              panli(jl,jk,2)*pbfract2(jl,jk,5)*pa4av2(jl,jk)+pso4_6(jl,jk)

        zso47=panli(jl,jk,1)*pbfract1(jl,jk,6)*pa4av1(jl,jk)+                &
              panli(jl,jk,2)*pbfract2(jl,jk,6)*pa4av2(jl,jk)+pso4_7(jl,jk)

        !--- 1.2) Determine number of particles that can be sufficiently coated
        !         by the available sulfate to be transfered to the soluble modes:

        !    Optimization of the call of m7_coat to allow for unroll and 
        !    subsequent vectorization.

!CDIR UNROLL=7
        DO jmod = 1, nmod
          zm6rp(jmod) = pm6rp(jl,jk,jmod)
        END DO

        CALL m7_coat(zm6rp,zcrtcst)

        !@@@ Changed security check to allow for inconsistent radii:

        IF(paernl(jl,jk,iaiti) >= 1.E-5 .AND. zcrtcst(5)>zeps) THEN
           zcrit_5=MIN(paernl(jl,jk,iaiti), zso45/zcrtcst(5))
        ELSE
           zcrit_5=0.
        END IF
        IF(paernl(jl,jk,iacci) >= 1.E-5 .AND. zcrtcst(6)>zeps) THEN
           zcrit_6=MIN(paernl(jl,jk,iacci), zso46/zcrtcst(6))
        ELSE
           zcrit_6=0.
        END IF
        IF(paernl(jl,jk,icoai) >= 1.E-5 .AND. zcrtcst(7)>zeps) THEN
           zcrit_7=MIN(paernl(jl,jk,icoai), zso47/zcrtcst(7))
        ELSE
           zcrit_7=0.
        END IF

        !--- 1.3) Number of particles moved from the mode 5 to 2 due to
        !         interaction with 1 and due to condensation:
        
        paernl(jl,jk,iaits)=paernl(jl,jk,iaits)+zcrit_5
        paernl(jl,jk,iaiti)=paernl(jl,jk,iaiti)-zcrit_5
        
        !--- 1.4) Mass moved from mode 5 to 2:
        
        pa4delt(jl,jk,2)=pa4delt(jl,jk,2)+pso4_5(jl,jk)

        ! JadB: I use an 'own' zero-concentration cap for the masses
        ! instead of using the zero-concentration cap for the numbers.
        ! Those gave rounding errors and negative concentrations,
        ! especially with 8-byte floating points.
        ! The same will be done for modi 6 and 7.

        ! pa4delt(jl,jk,ibcks)=pa4delt(jl,jk,ibcks)+zcrit_5*pbcav5(jl,jk)*1.e6
        ! pa4delt(jl,jk,iocks)=pa4delt(jl,jk,iocks)+zcrit_5*pocav5(jl,jk)*1.e6
        IF(zcrit_5 .NE. 0.) THEN
          ! Only transport mass if the same conditions
          ! as in the case of number transport are met.
          pa4delt(jl,jk,ibcks)=pa4delt(jl,jk,ibcks)+                         &
            MIN((zso45/zcrtcst(5))*pbcav5(jl,jk)*1.e6,paerml(jl,jk,ibcki))
          pa4delt(jl,jk,iocks)=pa4delt(jl,jk,iocks)+                         &
            MIN((zso45/zcrtcst(5))*pocav5(jl,jk)*1.e6,paerml(jl,jk,iocki))

        !--- 1.5) Mass remaining in mode 5:
         
          paerml(jl,jk,ibcki)=paerml(jl,jk,ibcki)-                           &
            MIN((zso45/zcrtcst(5))*pbcav5(jl,jk)*1.e6,paerml(jl,jk,ibcki))
          paerml(jl,jk,iocki)=paerml(jl,jk,iocki)-                           &
            MIN((zso45/zcrtcst(5))*pocav5(jl,jk)*1.e6,paerml(jl,jk,iocki))
        END IF
        ! paerml(jl,jk,ibcki)=paerml(jl,jk,ibcki)-zcrit_5*pbcav5(jl,jk)*1.e6
        ! paerml(jl,jk,iocki)=paerml(jl,jk,iocki)-zcrit_5*pocav5(jl,jk)*1.e6

        !--- 1.6) Number of particles moved from the mode 6 to 3:
        
        paernl(jl,jk,iaccs)=paernl(jl,jk,iaccs)+zcrit_6
        paernl(jl,jk,iacci)=paernl(jl,jk,iacci)-zcrit_6
        
        !--- 1.7) Mass moved from mode 6 to 3:

        pa4delt(jl,jk,3)=pa4delt(jl,jk,3)+pso4_6(jl,jk)
        ! pa4delt(jl,jk,iduas)=pa4delt(jl,jk,iduas)+zcrit_6*pduav6(jl,jk)*1.e6
        IF(zcrit_6 .NE. 0.) THEN
          pa4delt(jl,jk,iduas)=pa4delt(jl,jk,iduas)+                         &
            MIN((zso46/zcrtcst(6))*pduav6(jl,jk)*1.e6,paerml(jl,jk,iduai))

        !--- 1.8) Mass remaining in mode 6:

          paerml(jl,jk,iduai)=paerml(jl,jk,iduai)-                           &
            MIN((zso46/zcrtcst(6))*pduav6(jl,jk)*1.e6,paerml(jl,jk,iduai))
        END IF
        ! paerml(jl,jk,iduai)=paerml(jl,jk,iduai)-zcrit_6*pduav6(jl,jk)*1.e6
        
        !--- 1.9) Number of particles moved from the mode 7 to 4:

        paernl(jl,jk,icoas)=paernl(jl,jk,icoas)+zcrit_7
        paernl(jl,jk,icoai)=paernl(jl,jk,icoai)-zcrit_7

        !--- 1.10) Mass moved from mode 7 to 4:
        
        pa4delt(jl,jk,4)=pa4delt(jl,jk,4)+pso4_7(jl,jk)
        ! pa4delt(jl,jk,iducs)=pa4delt(jl,jk,iducs)+zcrit_7*pduav7(jl,jk)*1.e6
        IF(zcrit_7 .NE. 0.) THEN
          pa4delt(jl,jk,iducs)=pa4delt(jl,jk,iducs)+                         &
            MIN((zso47/zcrtcst(7))*pduav7(jl,jk)*1.e6,paerml(jl,jk,iduci))

        !--- 1.11) Mass remaining in mode 7:

          paerml(jl,jk,iduci)=paerml(jl,jk,iduci)-                           &
            MIN((zso47/zcrtcst(7))*pduav7(jl,jk)*1.e6,paerml(jl,jk,iduci))
        END IF
        ! paerml(jl,jk,iduci)=paerml(jl,jk,iduci)-zcrit_7*pduav7(jl,jk)*1.e6

     END DO
  END DO

END SUBROUTINE m7_concoag
