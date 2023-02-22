SUBROUTINE m7_nuck(kproma, kbdim,   klev,  pso4g,   &
                   ptp1,   prelhum, panew, pa4delt, &
                   ptime                              )
  !      
  !  Authors:
  !  --------
  !  J. Wilson, E. Vignati, JRC/EI, original source                 09/2000
  !  P. Stier, MPI                  f90-version, 
  !                                 changes, 
  !                                 comments,
  !                                 modularisation and 
  !                                 implementation of 
  !                                 Vehkamaeki (2002)             2001-2003
  !
  !  Purpose:                                                           
  !  --------
  !  This routine calls the routines for the computation of the 
  !  nucleation rate znucrate [molec. cm-3 s-1] and, for 
  !  Vehkamaeki (2002), of the number of molecules in the critical
  !  cluster from a given gas phase H2SO4 concentration 
  !  pso4g [molec. cm-3]. It also calculates the integrated change of 
  !  H2SO4 gas phase mass over one timestep due to nucleation 
  !  pa4delt(:,:,1) [molec. cm-3] as well as the number of nucleated 
  !  particles panew [1] during the timestep. Whilst this is done 
  !  analytically for the old Kulmala (1998) parameterization, it has
  !  to be done numerically for the new Vehkamaeki (2002) scheme.
  !
  !  Interface:
  !  ----------
  !  *m7_nuck* is called from *m7_dnum*
  !
  !  Method:
  !  -------
  !
  !  1) Kulmala et al. (1998):
  !
  !  For the Kulmala et al. (1998) scheme the formula for binary
  !  nucleation is rewritten to take the form 
  !  znucrate = exp[zalpha+ln(pso4g)*beta]. 
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
  !
  !  2) Vehkamaeki et al. (2002):
  !  
  !  An analytical integration of the nucleation equation is not
  !  possible, therefor the nucleation rate is simply multiplied
  !  by dt, implying a fixed gas-phase H2SO4 concentration over 
  !  the timestep.
  !@@@ The dependency of the results on the used timestep 
  !@@@ should be checked carefully in a sensitivity study! 
  !  The number of nucleated particles is then calculated by 
  !  taking the calculated critical mass critn of newly nucleated 
  !  particles and dividing the total mass of nucleated sulfate 
  !  by this critical mass of one nucleated particle. 
  !  
  !  Externals:
  !  ----------
  !  nucl_kulmala
  !  nucl_vehkamaeki
  !
  !  References:
  !  -----------
  !  Vehkamaeki et al. (2002), An improved parameterization for sulfuric
  !     acid/water nucleation rates for tropospheric and stratospheric
  !     conditions, J. Geophys. Res, 107, D22, 4622
  !  Kulmala et al. (1998), Parameterizations for sulfuric acid/water
  !     nucleation rates. J. Geophys. Res., 103, No D7, 8301-8307
  !  Vignatti, E. (1999), Modelling Interactions between Aerosols and 
  !     Gaseous Compounds in the Polluted Marine Atmosphere. PhD-Thesis,
  !     RISO National Laborartory Copenhagen, Riso-R-1163(EN)


  !USE mo_time_control,  ONLY: delta_time
  USE mo_control,       ONLY: nrow
  USE mo_kind,          ONLY: wp
  USE mo_aero_m7,       ONLY: critn, naermod, nnucl, avo
!  USE mo_aero_mem,      ONLY: d_nuc_so4
  USE mo_aero_nucl,     ONLY: nucl_kulmala, nucl_vehkamaeki, nucl_kulmala2006

  IMPLICIT NONE

  !--- Parameters:
  !
  ! pso4g          = mass of gas phase sulfate [molec. cm-3]
  ! ptp1           = atmospheric temperature at time t+1 [K]
  ! prelhum        = atmospheric relative humidity [%]
  ! pa4delt(:,:,1) = mass of H2SO4 added to the nucleation mode due 
  !                  to nucleation of H2SO4 over ztmst. 
  !                  Equilvalent to the integral of H2SO4 gas loss 
  !                  due to nucleation over timestep ztmst. [molec. cm-3]
  ! panew          = number of nucleated particles over timestep ztmst
  !                  panew=pa4delt/critn i.e. mass of formed sulfate 
  !                  divided by an assumed mass of a nucleus. [cm-3]
  !        
  !--- Local variables:
  !
  ! zso4g_new      = temporay storage of gas phase sulfate [molec. cm-3]
  ! zncrit         = number of molecules in the critical cluster [1]
  !
  ! See comments!

  INTEGER :: kproma, kbdim, klev

  REAL    :: ptime, time_step_len

  REAL    :: pso4g(kbdim,klev),          ptp1(kbdim,klev),     &
             prelhum(kbdim,klev),        panew(kbdim,klev)

  REAL    :: pa4delt(kbdim,klev,naermod)

  ! Local variables:

  INTEGER :: jk,          jl,           jrow

  REAL    :: ztmst,       zqtmst,       zf1,          zeps

  REAL    :: znucrate(kbdim,klev),  & ! nucleation rate [m-3 s-1] !@@@ Make consistent with Kulmala
             zso4g_new(kbdim,klev), & ! new gas phase sulfate concentration [molec. cm-3]
             zalpha(kbdim,klev),    & ! auxiliary coefficient for the analytical integration of Kulmala
             zbeta(kbdim,klev),     & ! auxiliary coefficient for the analytical integration of Kulmala
             zncrit(kbdim,klev)       ! number of molecules in the critical cluster [1]

  !--- 0) Initialisations: ------------------------------------------------------

  jrow=nrow(2)

  time_step_len = ptime
  
  ztmst=time_step_len

  zqtmst=1/ztmst

  zeps=EPSILON(1.0_wp)

  !--- 1) Calculate nucleation rate:

  IF(nnucl==1) THEN

     CALL nucl_vehkamaeki(kproma,   kbdim,   klev,  & ! ECHAM5 dimensions
                          ptp1,     prelhum, pso4g, & ! ECHAM5 temperature, relative humidity
                          znucrate, zncrit          )  


     !--- Calculate updated gas phase concentration:
     !
     !    N(t)   = N(0)   - znucrate   * zncrit * dt
     !    [cm-3] = [cm-3] - [cm-3 s-1] * [1]    * [s] 

     DO jk=1, klev
        DO jl=1, kproma

           zso4g_new(jl,jk)=pso4g(jl,jk)-(znucrate(jl,jk)*zncrit(jl,jk)*ztmst)

        END DO
     END DO

  ELSE IF (nnucl==2) THEN

     zncrit(:,:)=critn

     CALL nucl_kulmala(kproma,  kbdim,   klev,    &
                       pso4g,   ptp1,    prelhum, &
                       znucrate, zalpha, zbeta    )

     !--- 2) Analytical integration of the nucleation rate (Eq. 19) ----------------
     !       over timestep ztmst assuming no new H2SO4(g) production:
     !
     !       d(N_av/critn)/dt=exp(alpha + beta*ln(N_av)) => ... =>
     !           
     !       N_av(t0+dt)=
     !       [N_av(t0)**(1-beta) + critn*(1-beta)exp(alpha)*dt]**(1/(1-beta)
     !

     DO jk=1, klev
        DO jl=1, kproma

           IF (znucrate(jl,jk) .GT. 1e-10) THEN
              zf1 = pso4g(jl,jk)**(1.0-zbeta(jl,jk))-critn*EXP(zalpha(jl,jk))*(1.0-zbeta(jl,jk))*ztmst
              zso4g_new(jl,jk) = EXP(LOG(zf1)/(1.0 - zbeta(jl,jk)))
           ELSE
              zso4g_new(jl,jk) = pso4g(jl,jk)
           END IF

        END DO
     END DO
     
    ELSE IF(nnucl==3) THEN

     CALL nucl_kulmala2006(kproma,   kbdim,   klev,  & ! ECHAM5 dimensions
                           pso4g, ptp1,  & ! ECHAM5 temperature, relative humidity
                          znucrate, zncrit          )  


     !--- Calculate updated gas phase concentration:
     !
     !    N(t)   = N(0)   - znucrate   * zncrit * dt
     !    [cm-3] = [cm-3] - [cm-3 s-1] * [1]    * [s] 

     DO jk=1, klev
        DO jl=1, kproma

           zso4g_new(jl,jk)=pso4g(jl,jk)-(znucrate(jl,jk)*zncrit(jl,jk)*ztmst)

        END DO
     END DO

  END IF

  !--- 3) Calculate changes in gas-phase and aerosol mass and aerosol numbers: --

  DO jk=1, klev
     DO jl=1, kproma

        IF( znucrate(jl,jk) > zeps  ) THEN

  !--- 3.1) Security check:

           zso4g_new(jl,jk) = MAX(zso4g_new(jl,jk), 0.0)
           zso4g_new(jl,jk) = MIN(zso4g_new(jl,jk), pso4g(jl,jk))
   
  !--- 3.2) Calculate mass of nucleated H2SO4 (equals the net 
  !         gas phase H2SO4 loss):
  !
           ! pa4delt(jl,jk,1) = (pso4g(jl,jk)-zso4g_new(jl,jk))
           ! JadB: A problem occurs. We calculate gas_new, wich is gas_old - nucleation, 
           ! and then calculate the nucleation again with gas_old - gas_new. 
           ! Huge errors when the time step is small.
           ! That's why we take the nucleation as fresh as possible (from znucrate) and apply the limiters.
           pa4delt(jl,jk,1) = (znucrate(jl,jk)*zncrit(jl,jk)*ztmst)
           pa4delt(jl,jk,1) = Max(pa4delt(jl,jk,1),0.0)
           pa4delt(jl,jk,1) = Min(pa4delt(jl,jk,1),pso4g(jl,jk))
  !
  !--- 3.3) Calculate the number of nucleated particles (nucleated mass 
  !         divided by the assumed mass of a critical cluster critn):

           panew(jl,jk)=pa4delt(jl,jk,1)/zncrit(jl,jk)      
  !
  !--- 3.4) Calculate changes in gas phase H2SO4 due to nucleation:
  !
           pso4g(jl,jk)=pso4g(jl,jk)-pa4delt(jl,jk,1)

  !--- 4) Vertically integrate mass of nucleated sulfate for diagnostics
  !       Convert [molec. cm-3] to [kg(S) m-2]:
  

        END IF

     END DO
  END DO
!  write(*,*) 'nuck', 'nucl part= ',panew(2100,1), 'so4gn', pso4g(2100,1), 'deltag', pa4delt(2100,1,1)
END SUBROUTINE m7_nuck
