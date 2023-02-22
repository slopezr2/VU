SUBROUTINE m7_dnum(kproma, kbdim,   klev,          &
                   pso4g,  paerml,  paernl, ptp1,  &
                   papp1,  prelhum, pm6rp,  prhop, &
                   pso4_5, pso4_6,  pso4_7, ptime    )
  !
  !  *m7_dnum*  changes gas phase sulfate, aerosol numbers and masses
  !             due to nucleation and coagulation
  !
  !  Authors:
  !  ---------
  !  J. Wilson  and E. Vignati, JRC/EI (original source)                09/2000
  !  P. Stier, MPI                     (f90-version, changes, comments)    2001 
  !
  !  Version: 
  !  --------- 
  !  This version is equivalent to the version dnum2 of the m7 boxmodel. 
  !
  !  Purpose
  !  ---------
  !  This routine calculates new gas phase sulfate and aerosol
  !  numbers and masses after timestep ztmst.
  !
  !  Interface
  !  -----------
  !  *m7_dnum* is called from *m7*
  !
  !  Externals
  !  -----------
  !  *m7_coaset*   calculates the coagulation kernels
  !  *m7_nuck*     calculates the integral mass nucleated sulfate and
  !                the number of nucleated particles over one timestep
  !  *m7_delcoa*   integrates equations for the changes in aerosol numbers
  !                dn/dt=c -an2 -bn over one timestep and calculates the 
  !                corresponding changes in aerosol masses
  !  *m7_concoag*  calculates particle numbers and mass moved from the 
  !                insoluble to the mixed modes due to coagulation with
  !                smaller mixed particles and condensation of H2SO4.
  !
  !  Warning:
  !  --------
  !  For optimization purposes currently only "physically reasonable" elements of the
  !  coagulation kernel zcom are calculated in m7_concoag. These elements are specified
  !  in the matrix locoagmask in mo_aero_m7. Check carefully and adapt locoagmask 
  !  accordingly  before changes in the code below.


  USE mo_aero_m7, ONLY: naermod, nmod,  lsnucl, lscoag

  !use m7_delcoa_mod, only : m7_delcoa

  IMPLICIT NONE

  !--- Parameters:
  !
  !  pso4g      = mass of gas phase sulfate [molec. cm-3]
  !  paerml     = total aerosol mass for each compound 
  !               [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and dust]
  !  paernl     = aerosol number for each mode [cm-3]
  !  ptp1       = atmospheric temperature at time t+1 [K]
  !  papp1      = atmospheric pressure at time t+1 [Pa]
  !  prelhum    = atmospheric relative humidity [%]
  !  pm6rp      = mean mode actual radius (wet radius for soluble modes 
  !               and dry radius for insoluble modes) [cm]
  !  prhop      = mean mode particle density [g cm-3]
  !  pso4_x     = mass of sulphate condensed on insoluble mode x []
  !
  !--- Local variables:
  !
  ! zcom            = general coagulation coefficient []
  !                   (calculated in m7_coaset)
  ! za              = effectively used coagulation coefficient for 
  !                   unimodal coagulation []
  ! zb              = effectively used coagulation coefficient for 
  !                   inter-modal coagulation []
  ! zbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x (finally calculated in m7_delcoa)
  !                   that is moved to mode y / y+1 (modes 5,6,7 and mode 1,2 resp.) [1]
  !                   !@@@ Careful! Inconsistent usage!!!
  ! za4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one timstep 
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in m7_nuck)
  !                      - coagulation (calculated in m7_concoag)
  ! zanew           = number of nucleated particles
  !                   zanew=za4delt/critn i.e. mass of formed sulfate 
  !                   divided by an assumed mass of a nucleus. []
  !                   (calculated in m7_nuck)
  ! zxxavy          = average mass of species xx in mode y []
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  !                   [molecules for sulfate and ug for others]
  !
  !@@@ to be continued!

  !--- Parameters:

  INTEGER :: kproma, kbdim, klev

  REAL    :: ptime

  REAL    :: pso4g(kbdim,klev),          ptp1(kbdim,klev),                 &
             papp1(kbdim,klev),          prelhum(kbdim,klev),              &
             pso4_5(kbdim,klev),         pso4_6(kbdim,klev),               &
             pso4_7(kbdim,klev)

  REAL    :: paerml(kbdim,klev,naermod), paernl(kbdim,klev,nmod),          &
             pm6rp(kbdim,klev,nmod),     prhop(kbdim,klev,nmod)

             
  ! Local Variables:

  INTEGER :: jl, jk

  REAL    :: zanew(kbdim,klev)

  REAL    :: za(kbdim,klev,nmod),        zb(kbdim,klev,nmod),              &
             za4delt(kbdim,klev,naermod)

  REAL    :: zbfract1(kbdim,klev,nmod-1),zbfract2(kbdim,klev,nmod-1),      &
             zbfract5(kbdim,klev,3),     zbfract6(kbdim,klev,2),           &
             zbfract7(kbdim,klev,2)
             

  REAL    :: zcom(kbdim,klev,nmod,nmod)


  !--- 0) Initialisations: ----------------------------------------------

  za4delt(:,:,:) = 0.  ! Mode 1 changed by m7_nuck if lsnucl=TRUE .
                       ! Has to be initialized for the other modes.
  zanew(:,:)     = 0.  ! Changed by m7_nuck if lsnucl=TRUE .

  !--- 1) Calculate  coagulation coefficients: --------------------------
  !
 
  IF (lscoag) CALL m7_coaset(kproma, kbdim, klev,  paernl, ptp1, &
                             papp1,  pm6rp, prhop, zcom          )
  !
  !--- 2) Calculate nucleation rate, number of nucleated particles ------
  !       and changes in gas phase sulfate mass.
  !
 
  IF (lsnucl) CALL m7_nuck(kproma, kbdim,   klev,  pso4g,   &
                           ptp1,   prelhum, zanew, za4delt, &
                           ptime                               )

  !
  !--- 3) Assign coagulation coefficients (unimodal coagulation)---------
  !       and the normalised fraction of particle numbers moved 
  !       between the modes (inter-modal coagulation):
  !
  !       The general equation for dn/dt for each mode is:
  !
  !       dn/dt=-za*n^2 - zb*n + zc
  !
  !       where za=unimodal coagulation coefficient (zcom(mod))
  !             zb=inter-modal coagulation with higher modes
  !                (zcom(mod) * n(jmod+1))
  !             zc=particle formation rate 
  !                (=zanew/ztmst if jmod=1, zero for higher modes) 
  !
  !             zb is zero when n(jmod+1)=zero, or jmod=naermod
  !
  IF (lscoag) THEN

     DO jk=1,klev
        DO jl=1,kproma 

           !---  3.1) Unimodal coagulation coefficients:
           !@@@Coag:
           za(jl,jk,1)=zcom(jl,jk,1,1)/2.    ! Unimodal coagulation
           za(jl,jk,2)=zcom(jl,jk,2,2)/2.    ! only allowed for modes
           za(jl,jk,3)=zcom(jl,jk,3,3)/2.    ! 1,2,3 and 5.
           za(jl,jk,4)=0.
           za(jl,jk,5)=zcom(jl,jk,5,5)/2.
           za(jl,jk,6)=0.
           za(jl,jk,7)=0.
           !
           !---  3.2) Inter-modal coagulation - soluble modes
           ! 
           !--- Sum all higher mode coagulation terms for 
           !    soluble modes 1,2,3,4:
           !
           !--- 3.2.1) Mode 1:
           
           !--- Number of particles (zbfract1(:,:,x)) that are moved 
           !    from mode 1 to the mode x+1 :
           !    !@@@ Clumsy! Change to x - also in concoag!!!
           zbfract1(jl,jk,1)=zcom(jl,jk,2,1)*paernl(jl,jk,2)
           zbfract1(jl,jk,2)=zcom(jl,jk,3,1)*paernl(jl,jk,3)
           zbfract1(jl,jk,3)=zcom(jl,jk,4,1)*paernl(jl,jk,4)
           zbfract1(jl,jk,4)=zcom(jl,jk,5,1)*paernl(jl,jk,5)
           zbfract1(jl,jk,5)=zcom(jl,jk,6,1)*paernl(jl,jk,6)
           zbfract1(jl,jk,6)=zcom(jl,jk,7,1)*paernl(jl,jk,7)
           !
           !--- Sum of all particles that are moved from mode 1:
           !
           zb(jl,jk,1)=zbfract1(jl,jk,1)+zbfract1(jl,jk,2)+            &
                       zbfract1(jl,jk,3)+zbfract1(jl,jk,4)+            &
                       zbfract1(jl,jk,5)+zbfract1(jl,jk,6)
           !
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 1:
           !
           IF (zb(jl,jk,1).GT.0.0) THEN
              zbfract1(jl,jk,1)=zbfract1(jl,jk,1)/zb(jl,jk,1)
              zbfract1(jl,jk,2)=zbfract1(jl,jk,2)/zb(jl,jk,1)
              zbfract1(jl,jk,3)=zbfract1(jl,jk,3)/zb(jl,jk,1)
              zbfract1(jl,jk,4)=zbfract1(jl,jk,4)/zb(jl,jk,1)
              zbfract1(jl,jk,5)=zbfract1(jl,jk,5)/zb(jl,jk,1)
              zbfract1(jl,jk,6)=zbfract1(jl,jk,6)/zb(jl,jk,1)
           END IF
           !
           !--- 3.2.2) Mode 2:
           !
           !--- Number of particles (zbfract1(:,:,x)) that are moved 
           !    from mode 2 to the mode x+1 :
           !
           zbfract2(jl,jk,2)=zcom(jl,jk,3,2)*paernl(jl,jk,3)
           zbfract2(jl,jk,3)=zcom(jl,jk,4,2)*paernl(jl,jk,4)
           !zbfract2(jl,jk,4)=zcom(jl,jk,5,2)*paernl(jl,jk,5)
           zbfract2(jl,jk,4)=0.
           zbfract2(jl,jk,5)=zcom(jl,jk,6,2)*paernl(jl,jk,6)
           zbfract2(jl,jk,6)=zcom(jl,jk,7,2)*paernl(jl,jk,7)
           
           !--- Sum of all particles that are moved from mode 2:
           
           zb(jl,jk,2)=zbfract2(jl,jk,2)+zbfract2(jl,jk,3)+            &
                       zbfract2(jl,jk,4)+zbfract2(jl,jk,5)+            &
                       zbfract2(jl,jk,6)

           !--- Normalize particle numbers by the total number of 
           !    particles moved from mode 2:

           IF (zb(jl,jk,2).GT.0.0) THEN
              zbfract2(jl,jk,2)=zbfract2(jl,jk,2)/zb(jl,jk,2)
              zbfract2(jl,jk,3)=zbfract2(jl,jk,3)/zb(jl,jk,2)
              zbfract2(jl,jk,4)=zbfract2(jl,jk,4)/zb(jl,jk,2)
              zbfract2(jl,jk,5)=zbfract2(jl,jk,5)/zb(jl,jk,2)
              zbfract2(jl,jk,6)=zbfract2(jl,jk,6)/zb(jl,jk,2)
           END IF
        
           !--- 3.2.3) Mode 3 and Mode 4 - considered ineffective:

           zb(jl,jk,3)=0.0

           zb(jl,jk,4)=0.0

           !
           !--- 3.3) Inter-modal coagulation - insoluble modes
           !
           !         For the insoluble modes coagulation with soluble modes
           !         is a sink as they are transfered to the corresponding
           !         mixed/solublemode. Therefore, terms with a lower mode 
           !         number or the same mode number are included. 
           !         (!@@@ There are still some switches for testing.)
           !
           !--- 3.3.1) Mode 5:
           !
           !--- Number of particles (zbfract5(:,:,x)) that are moved 
           !    from mode 5 to the mode x:
           
           !@@@ zbfract5(jl,jk,1)= zcom(jl,jk,1,5)*paernl(jl,jk,1)    
           zbfract5(jl,jk,1)=0.
           zbfract5(jl,jk,2)=zcom(jl,jk,2,5)*paernl(jl,jk,2)
           zbfract5(jl,jk,3)=zcom(jl,jk,3,5)*paernl(jl,jk,3)
           
           !--- Sum of all particles that are moved from mode 5:
           
           zb(jl,jk,5)=zbfract5(jl,jk,1)+zbfract5(jl,jk,2)+zbfract5(jl,jk,3)
           
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 5:
           
           IF (zb(jl,jk,5).GT.0.0) THEN
              zbfract5(jl,jk,1)=zbfract5(jl,jk,1)/zb(jl,jk,5)
              zbfract5(jl,jk,2)=zbfract5(jl,jk,2)/zb(jl,jk,5)
              zbfract5(jl,jk,3)=zbfract5(jl,jk,3)/zb(jl,jk,5)
           END IF
           
           !--- 3.3.2) Mode 6:
           
           !--- Number of particles (zbfract6(:,:,x)) that are moved 
           !    from mode 6 to the mode x:
           
           !@@@zbfract6(jl,jk,1)=zcom(jl,jk,1,6)*paernl(jl,jk,1)
           !@@@zbfract6(jl,jk,2)=zcom(jl,jk,2,6)*paernl(jl,jk,2)
           zbfract6(jl,jk,1)=0.
           zbfract6(jl,jk,2)=0.
           
           !--- Sum of all particles that are moved from mode 6:
           
           zb(jl,jk,6)=zbfract6(jl,jk,1)+zbfract6(jl,jk,2)
           
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 6:
           
           IF (zb(jl,jk,6).GT.0.0) THEN
              zbfract6(jl,jk,1)=zbfract6(jl,jk,1)/zb(jl,jk,6)
              zbfract6(jl,jk,2)=zbfract6(jl,jk,2)/zb(jl,jk,6)
           END IF
           
           !--- 3.3.3) Mode 7:
           
           !--- Number of particles (zbfract7(:,:,x)) that are moved 
           !    from mode 7 to the mode x:
           
           !@@@ zbfract7(jl,jk,1)=zcom(jl,jk,1,7)*paernl(jl,jk,1)
           !@@@ zbfract7(jl,jk,2)=zcom(jl,jk,2,7)*paernl(jl,jk,2)
           zbfract7(jl,jk,1)=0.
           zbfract7(jl,jk,2)=0.  
           
           !--- Sum of all particles that are moved from mode 7:
           
           zb(jl,jk,7)=zbfract7(jl,jk,1)+zbfract7(jl,jk,2)
           
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 7:
           
           IF (zb(jl,jk,7).GT.0.0) THEN
              zbfract7(jl,jk,1)=zbfract7(jl,jk,1)/zb(jl,jk,7)
              zbfract7(jl,jk,2)=zbfract7(jl,jk,2)/zb(jl,jk,7)
           END IF

        END DO
     END DO
     !
  ELSE

     za(:,:,:)       = 0.
     zb(:,:,:)       = 0.
     zbfract1(:,:,:) = 0.
     zbfract2(:,:,:) = 0.
     zbfract5(:,:,:) = 0.
     zbfract6(:,:,:) = 0.
     zbfract7(:,:,:) = 0.

  END IF !(lscoag)
  !
  !
  
  CALL m7_delcoa(kproma,   kbdim,  klev,     paerml,   &
                 paernl,   pm6rp,  za4delt,  zanew,    &
                 za,       zb,     zbfract1, zbfract2, &
                 zbfract5, pso4_5, pso4_6,   pso4_7, ptime    )
 ! stop 'break after delcoa'

END SUBROUTINE m7_dnum
