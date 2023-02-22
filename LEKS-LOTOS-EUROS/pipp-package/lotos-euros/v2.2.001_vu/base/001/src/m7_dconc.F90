SUBROUTINE m7_dconc (kproma, kbdim, klev, paerml, paernl, pm6dry)
  !
  !    *m7_dconc*  changes aerosol numbers and masses to account for
  !                condensational growth of the mode mean radii
  !  
  !    Authors:
  !    --------
  !    J. Wilson and E. Vignati, JRC (original source)            May 2000
  !    P. Stier, MPI-MET (f90 version, changes, comments)             2001
  !
  !    Purpose:
  !    --------
  !    This routine repartitions aerosol number and mass between the
  !    the modes to account for condensational growth and the formation
  !    of an accumulation mode from the upper tail of the aitken mode.
  !
  !    Interface:
  !    ----------
  !    *m7_dconc* is called from *m7*
  !
  !    Method:
  !    -------
  !    The routine calculates the cumulativ number and mass distribution of the
  !    modes up to the respective mode boundary:
  !
  !                        / x                              _
  !                 N      |       1           1   ln(R)-ln(R)  2
  !    N(0,x) = ---------  |   --------  exp(- - ( ----------- )   ) d ln(R) 
  !             ln(sigma)  |   sqrt(2PI)       2    ln(sigma)
  !                        / 0 
  !                         
  !                         /tx                   2
  !                        |        1            t
  !           =     N      |     --------  exp(- - ) d t 
  !                        |     sqrt(2PI)       2 
  !                        /-inf
  ! 
  !    where:                   
  !
  !                        _
  !               ln(R)-ln(R)
  !    t      =   -----------
  !                ln(sigma)
  !
  !    and:
  !                        _
  !               ln(x)-ln(R)
  !    tx     =   -----------
  !                ln(sigma)
  !    _
  !    R is the Count Mean Radius or the Mass Mean Radius.
  !
  !    Practically, the routine m7_cumnor calculates the fraction of the number and
  !    mass distribution for each mode lying below the respective upper mode boundary (1).
  !    In a next step the net fraction of each mode lying between the upper and lower
  !    mode boundaries are summed up (2) and the numbers and masses exceeding the mode
  !    boundaries are transfered to the neighboring larger mode (3). 
  !    Finally, these quantities are stored in the respective arrays
  !    paernl and paerml (4).
  !    The repartititioning is currently only done for the soluble modes as it is
  !    assumed that insoluble modes are rather transfered to the soluble modes 
  !    and grow as soluble particles.
  !
  !    Externals:
  !    ----------
  !    None
  !
  !--- Parameter list:
  !
  !    paerml(kbdim,klev,naermod)= total aerosol mass for each compound 
  !                                [molec. cm-3 for sulfate and ug m-3 for others]
  !    paernl(kbdim,klev,nmod)   = aerosol number for each mode [cm-3]  !
  !    sigma(jmod)               = standard deviation of mode jmod [1]
  !    crdiv                     = threshold radii between the different modes [cm]
  !                                crdiv(jmod) is the lower bound and crdiv(jmod+1) is 
  !                                the upper bound of the respective mode
  !
  !--- Local Variables:
  !
  !    zfconn(:,:,jnum,jmod)     = absolute fraction of the number of particles in mode jmod, 
  !                                with CMD=2*pm6dry(jmod) and a geometric standard 
  !                                deviation zrcsig, that are smaller than crdiv(jnum+1).
  !                                I.e. with 0 < CMD < crdiv(jnum+1) [1]
  !    zfconm(:,:,jnum,jmod)     = absolute fraction of the mass in mode jmod, 
  !                                with CMD=2*pm6dry(jmod) and a geometric standard 
  !                                deviation zrcsig, that are smaller than crdiv(jnum+1).
  !                                I.e. with 0 < CMD < crdiv(jnum+1) [1]

  USE mo_kind,    ONLY: wp
  USE mo_aero_m7, ONLY: sigma, crdiv,          &
                        nmod,  naermod, nsol,  &
                        iso4ns,iso4ks,  iso4as,&
                        ibcks, ibcas,   ibccs, &
                        iocks, iocas,   ioccs, &
                        issas, isscs,          &
                        iduas, iducs,          &
                        pi,    avo,     wh2so4,&
                        dh2so4,dbc,     doc,   &
                        dnacl, ddust,          &
                        cmr2ram,cmedr2mmedr

  IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev

  REAL    :: paerml(kbdim,klev,naermod),   paernl(kbdim,klev,nmod),              &
             pm6dry(kbdim,klev,nsol)

  ! Local variables:
  !
  INTEGER :: jmod,          jnum,          jl,           jk

  REAL    :: zrcsig,        zarg1,         zarg2,        zdpmm,         zdpcm,   &
             zarg3,         zdpam,         zcongn,       zcongm,        zdummy,  &
             zr1,           zr2,           zttnj,        zavnj,         zmrj,    &
             zmt,           znt,           zavmt,        zmcr,          zfconmj, &
             zntnew,        zmtnew,        zdm,          zeps

  REAL    :: zambc2(4),     zambc3(4),     zambc4(4),                            &
             zamoc2(4),     zamoc3(4),     zamoc4(4),                            &
                            zamss3(4),     zamss4(4),                            &
                            zamdu3(4),     zamdu4(4)

  REAL    :: ztotmass(kbdim,klev)

  REAL    :: zsumn(kbdim,klev,nmod),       zsumm(kbdim,klev,nmod),               &
             zsumbc(kbdim,klev,3),         zsumoc(kbdim,klev,3),                 &
             zsumss(kbdim,klev,2),         zsumdu(kbdim,klev,2)
 
  REAL    :: zfconn(kbdim,klev,nsol,nsol), zfconm(kbdim,klev,nsol,nsol)


  !--- 0) Initialisations: ----------------------------------------------------------

  zeps=EPSILON(1.0_wp)
 

  zsumn(:,:,:)    = 0.
  zsumm(:,:,:)    = 0.
  zsumbc(:,:,:)   = 0.
  zsumoc(:,:,:)   = 0.
  zsumss(:,:,:)   = 0.
  zsumdu(:,:,:)   = 0.

  zfconm(:,:,:,:) = 0.
  zfconn(:,:,:,:) = 0.

  !
  !--- 1) Identify how much the mode jmod has grown into the next higher mode -------
  !
  DO jmod=1,nsol-1

     !--- Total mass of the mode in equivalent molecules of sulfate:

     SELECT CASE(jmod)
     CASE(1)
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4ns)
     CASE(2)
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4ks) + &
                               (paerml(1:kproma,:,ibcks)/dbc+paerml(1:kproma,:,iocks)/doc)  &
                               *dh2so4/wh2so4*avo*1.E-12
     CASE(3)
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4as) + &
                               (paerml(1:kproma,:,ibcas)/dbc+paerml(1:kproma,:,iocas)/doc+  &
                                paerml(1:kproma,:,issas)/dnacl+paerml(1:kproma,:,iduas)/ddust) &
                               *dh2so4/wh2so4*avo*1.E-12
     END SELECT

     DO jnum=jmod,nsol-1
        DO jk=1,klev
           DO jl=1,kproma
              IF (paernl(jl,jk,jmod) .GT. zeps .AND. pm6dry(jl,jk,jmod) .GT. 0.0) THEN

                 !--- 1.1) Calculate necessary parameters:                 

                 !--- Geometric Standard Deviation:
                 zrcsig=LOG(sigma(jmod))

                 !--- Mass Median Radius:
                 zarg1=pm6dry(jl,jk,jmod)*cmedr2mmedr(jmod)

                 !--- Count Median Radius:
                 zarg2=pm6dry(jl,jk,jmod)

                 !--- Threshold radius between the modes:
                 zarg3=crdiv(jnum+1)
            
                 !--- Transfer to logarithmic scale:
                 zdpmm=LOG(zarg1)
                 zdpcm=LOG(zarg2)
                 zdpam=LOG(zarg3)

                 !--- Distance of the CMD of the mode from the threshold mode  
                 !    diameter in terms of geometric standard deviations:

                 zcongn=(zdpam-zdpcm)/zrcsig
                 

                 !--- Distance of the MMD of the mode from the threshold mode  
                 !    diameter in terms of geometric standard deviations (t):

                 zcongm=(zdpam-zdpmm)/zrcsig
                 
                 !--- Calculate the cumulative of the log-normal number distribution:

                 CALL m7_cumnor(zcongn,zfconn(jl,jk,jnum,jmod),zdummy)
                

                 !--- Limit transfer only to adjacent modes:

                 IF (jnum .GT. jmod) THEN
                    zfconn(jl,jk,jnum,jmod)= 1.0
                 END IF

                 !--- Set minimum radius and maximum radius:

                 zr1 = crdiv(jmod)
                 zr2 = crdiv(jmod+1)

                 !--- Radius of average mass for a lognormal distribution

                 zdm = EXP((LOG(zr1)+LOG(zr2))/2.0)*cmr2ram(jmod)

                 !--- Average mass contained in the mode

                 zttnj = ztotmass(jl,jk)/paernl(jl,jk,jmod)

                 !--- Average number of sulfate molecules or equivalent for mixed modes,
                 !    for a particle with radius zdm

                 zavnj=zdm**3.0*pi*avo*dh2so4/wh2so4/0.75

                 !--- If the average mass contained in the mode is larger than the average mass
                 !    for a particle with radius zdm, the transfer of number and mass is done,
                 !    else there is no transfer

                 IF (zttnj .GT. zavnj .AND. jnum .EQ. jmod) THEN

                    !--- Mass remaining in the mode

                    zmrj=zfconn(jl,jk,jnum,jmod)*paernl(jl,jk,jmod)*zavnj

                    !--- Mass transferred

                    zmt=ztotmass(jl,jk)-zmrj

                    !--- Numbers transferred 

                    znt=(1.0-zfconn(jl,jk,jnum,jmod))*paernl(jl,jk,jmod)
                      
                    !--- Average mass of particles transferred                                      

                    IF(znt>zeps) THEN
                       zavmt=zmt/znt
                    ELSE
                       zavmt=0.
                    END IF

                    !--- Average mass of particles of radius zr2

                    zmcr=(zr2*cmr2ram(jmod))**3.0*pi*avo*dh2so4/wh2so4/0.75

                    !--- If the average mass of particle transferred is smaller than the average mass
                    !    mass of particles with radius zr2 then reduce the particles transferred
                    !    so that zavmt=zmcr, else calculate the mass fraction transferred zfconmj

                    IF (zavmt .GE. zmcr) THEN
                       zfconmj=zmrj/ztotmass(jl,jk)
                    ELSE
                       zntnew = znt/(1.0 + (zmcr-zavmt)/(zavmt-zavnj))
                       zmtnew = zntnew*zmcr
                       zfconmj = 1.0 - zmtnew/ztotmass(jl,jk)
                       zfconn(jl,jk,jnum,jmod) = 1.0 - zntnew/paernl(jl,jk,jmod)
                    END IF
                    zfconm(jl,jk,jnum,jmod)=zfconmj
                 ELSE
                    zfconn(jl,jk,jnum,jmod)=1.
                    zfconm(jl,jk,jnum,jmod)=1.
                 END IF
              ELSE
                 zfconn(jl,jk,jnum,jmod)=1.
                 zfconm(jl,jk,jnum,jmod)=1.
              END IF
           END DO
        END DO
     END DO
  END DO

  DO jmod=1,nsol
     DO jk=1,klev
        DO jl=1,kproma           

           !--- 2) Calculate the net fraction of mode jmod that is transfered -------
           !       to the mode jnew zfconn(:,:,jnew,jmod) :
           
           !--- Numbers:
           zfconn(jl,jk,4,jmod)=1.0                 -zfconn(jl,jk,3,jmod)
           zfconn(jl,jk,3,jmod)=zfconn(jl,jk,3,jmod)-zfconn(jl,jk,2,jmod)
           zfconn(jl,jk,2,jmod)=zfconn(jl,jk,2,jmod)-zfconn(jl,jk,1,jmod)

           !--- Mass:
           zfconm(jl,jk,4,jmod)=1.0                 -zfconm(jl,jk,3,jmod)
           zfconm(jl,jk,3,jmod)=zfconm(jl,jk,3,jmod)-zfconm(jl,jk,2,jmod)
           zfconm(jl,jk,2,jmod)=zfconm(jl,jk,2,jmod)-zfconm(jl,jk,1,jmod)

           !--- 3) Sum the net masses and numbers transfered between the modes: -----

           !--- 3.1) Soluble mode numbers and sulfate mass:

           zsumn(jl,jk,1)=zsumn(jl,jk,1)+paernl(jl,jk,jmod)*zfconn(jl,jk,1,jmod)
           zsumm(jl,jk,1)=zsumm(jl,jk,1)+paerml(jl,jk,jmod)*zfconm(jl,jk,1,jmod)
           zsumn(jl,jk,2)=zsumn(jl,jk,2)+paernl(jl,jk,jmod)*zfconn(jl,jk,2,jmod)
           zsumm(jl,jk,2)=zsumm(jl,jk,2)+paerml(jl,jk,jmod)*zfconm(jl,jk,2,jmod)
           zsumn(jl,jk,3)=zsumn(jl,jk,3)+paernl(jl,jk,jmod)*zfconn(jl,jk,3,jmod)
           zsumm(jl,jk,3)=zsumm(jl,jk,3)+paerml(jl,jk,jmod)*zfconm(jl,jk,3,jmod)
           zsumn(jl,jk,4)=zsumn(jl,jk,4)+paernl(jl,jk,jmod)*zfconn(jl,jk,4,jmod)
           zsumm(jl,jk,4)=zsumm(jl,jk,4)+paerml(jl,jk,jmod)*zfconm(jl,jk,4,jmod)

        END DO
     END DO
  END DO

  DO jk=1,klev
     DO jl=1,kproma

        !--- 3.2) Non-sulfate masses:

        zambc2(2)=paerml(jl,jk,ibcks)*zfconm(jl,jk,2,2)
        zambc2(3)=paerml(jl,jk,ibcks)*zfconm(jl,jk,3,2)
        zambc2(4)=paerml(jl,jk,ibcks)*zfconm(jl,jk,4,2)
        zamoc2(2)=paerml(jl,jk,iocks)*zfconm(jl,jk,2,2)
        zamoc2(3)=paerml(jl,jk,iocks)*zfconm(jl,jk,3,2)
        zamoc2(4)=paerml(jl,jk,iocks)*zfconm(jl,jk,4,2)
        zambc3(3)=paerml(jl,jk,ibcas)*zfconm(jl,jk,3,3)
        zambc3(4)=paerml(jl,jk,ibcas)*zfconm(jl,jk,4,3)
        zamoc3(3)=paerml(jl,jk,iocas)*zfconm(jl,jk,3,3)
        zamoc3(4)=paerml(jl,jk,iocas)*zfconm(jl,jk,4,3)
        zambc4(4)=paerml(jl,jk,ibccs)
        zamoc4(4)=paerml(jl,jk,ioccs)
        zamss3(3)=paerml(jl,jk,issas)*zfconm(jl,jk,3,3)
        zamss3(4)=paerml(jl,jk,issas)*zfconm(jl,jk,4,3)
        zamdu3(3)=paerml(jl,jk,iduas)*zfconm(jl,jk,3,3)
        zamdu3(4)=paerml(jl,jk,iduas)*zfconm(jl,jk,4,3)
        zamss4(4)=paerml(jl,jk,isscs)
        zamdu4(4)=paerml(jl,jk,iducs)
        
        zsumbc(jl,jk,1)=zambc2(2)
        zsumbc(jl,jk,2)=zambc2(3)+zambc3(3)
        zsumbc(jl,jk,3)=zambc2(4)+zambc3(4)+zambc4(4)
        zsumoc(jl,jk,1)=zamoc2(2)
        zsumoc(jl,jk,2)=zamoc2(3)+zamoc3(3)
        zsumoc(jl,jk,3)=zamoc2(4)+zamoc3(4)+zamoc4(4)
        zsumss(jl,jk,1)=zamss3(3)
        zsumss(jl,jk,2)=zamss3(4)+zamss4(4)
        zsumdu(jl,jk,1)=zamdu3(3)
        zsumdu(jl,jk,2)=zamdu3(4)+zamdu4(4)

     END DO
  END DO

  !--- 4) Store final masses and numbers of the modes: ------------------------

  DO jmod=1,nsol
     DO jk=1,klev
        DO jl=1,kproma

           !--- Particle numbers:

           paernl(jl,jk,jmod)=zsumn(jl,jk,jmod)

           !--- Sulfate mass:

           paerml(jl,jk,jmod)=zsumm(jl,jk,jmod)

        END DO
     END DO
  END DO

  !--- Non sulfate masses:

  DO jk=1,klev
     DO jl=1,kproma

        paerml(jl,jk,ibcks)=zsumbc(jl,jk,1)
        paerml(jl,jk,ibcas)=zsumbc(jl,jk,2)
        paerml(jl,jk,ibccs)=zsumbc(jl,jk,3)
        paerml(jl,jk,iocks)=zsumoc(jl,jk,1)
        paerml(jl,jk,iocas)=zsumoc(jl,jk,2)
        paerml(jl,jk,ioccs)=zsumoc(jl,jk,3)
        paerml(jl,jk,issas)=zsumss(jl,jk,1)
        paerml(jl,jk,isscs)=zsumss(jl,jk,2)
        paerml(jl,jk,iduas)=zsumdu(jl,jk,1)
        paerml(jl,jk,iducs)=zsumdu(jl,jk,2)

     END DO
  END DO

END SUBROUTINE m7_dconc
