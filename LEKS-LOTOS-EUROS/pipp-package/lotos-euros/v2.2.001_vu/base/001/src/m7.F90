SUBROUTINE m7(kproma, kbdim,   klev,           &  ! TM5  indices
              papp1,  prelhum, ptp1,           &  !   "   thermodynamics
              pso4g,  paerml, paernl,     &       !  M7   tracers
              prhop,  pww,    pm6rp,  pm6dry,  &  !   "   aerosol properties
              ptime )                             ! TM5  time step
  !
  !   ****m7* Aerosol model for the system so4,bc,oc,ss,dust in 7 modes.
  !
  !   Authors:
  !   ---------
  !   E. Vignati, JRC/EI    (original source)					    2000
  !   P. Stier, MPI         (f90-version, changes, comments)		2001 
  !   E. Vignati, JRC/IES   (so2 is not required in this version)   2005
  !
  !   Purpose
  !   ---------
  !   Aerosol model for the system so4,bc,oc,ss,dust in 7 modes.
  !
  !   Externals
  !   ---------
  !
  !   *m7_averageproperties* 
  !       calculates the average mass for all modes and the particle 
  !       dry radius and density for the insoluble modes.
  !    
  !   *m7_equiz*   
  !       calculates the ambient radius of sulphate particles
  !
  !   *m7_equimix* 
  !       calculates the ambient radius of so4,bc,oc (dust) particles
  !
  !   *m7_equil* 
  !       calculates the ambient radius of so4,ss particles 
  !
  !   *m7_dgas*    
  !       calculates the sulfate condensation on existing particles
  !
  !   *m7_dnum*    
  !       calculates new gas phase sulfate and aerosol numbers and masses
  !       after condensation, nucleation and coagulation over one timestep
  !
  !   *m7_dconc*   
  !       repartitions aerosol number and mass between the
  !       the modes to account for condensational growth and the formation
  !       of an accumulation mode from the upper tail of the aitken mode and 
  !       of a coarse mode from the upper tail of the accumulation mode
  !

  USE mo_aero_m7, ONLY: lsnucl, lscoag, lscond,         &
                        nmod,   nss,    nsol,   naermod

  IMPLICIT NONE 

  !--- Parameter list:
  !
  !  papp1      = atmospheric pressure at time t+1 [Pa]
  !  prelhum    = atmospheric relative humidity [% (0-1)]
  !  ptp1       = atmospheric temperature at time t+1 [K]
  !  pso4g      = mass of gas phase sulfate [molec. cm-3]
  !  paerml     = total aerosol mass for each compound 
  !               [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and dust]
  !  paernl     = aerosol number for each mode [cm-3]
  !  prhop      = mean mode particle density [g cm-3]
  !  pm6rp      = mean mode actual radius (wet radius for soluble modes 
  !               and dry radius for insoluble modes) [cm]
  !  pm6dry     = dry radius for soluble modes [cm]
  !  pww        = aerosol water content for each mode [kg(water) m-3(air)]
  !
  !--- Local variables:
  !
  !  zttn       = average mass for single compound in each mode 
  !               [in molec. for sulphate and in ug for bc, oc, ss, and dust]
  !  zhplus     = number of h+ in mole [???] (kg water)-1
  !  zso4_x     = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  !               (calculated in dgas used in concoag)


  ! Parameters:

  INTEGER :: kproma, kbdim, klev

  REAL    :: ptime

  REAL    :: prelhum(kbdim,klev),        papp1(kbdim,klev),           &
             ptp1(kbdim,klev),                                        &
             pso4g(kbdim,klev)
 
  REAL    :: paerml(kbdim,klev,naermod), paernl(kbdim,klev,nmod),     &
             pm6rp(kbdim,klev,nmod),     pm6dry(kbdim,klev,nsol),     &
             prhop(kbdim,klev,nmod),     pww(kbdim,klev,nmod)

  ! Local variables:

  !--- m7_box: 
  INTEGER :: i
  LOGICAL :: lofirst=.TRUE.
  !--- END m7_box

  REAL    :: zso4_5(kbdim,klev),         zso4_6(kbdim,klev),          &
             zso4_7(kbdim,klev)

  REAL    :: zhplus(kbdim,klev,nss)

  REAL    :: zttn(kbdim,klev,naermod)

  !
  !--- 0) Initialisations: -------------------------------------------------
  !
  zhplus(:,:,:) = 0.
  pm6dry(:,:,:) = 0. 
  pm6rp(:,:,:)  = 0.
  zttn(:,:,:)   = 0.
  prhop(:,:,:)  = 0.
  pww(:,:,:)    = 0. 
  zso4_5(:,:)   = 0.
  zso4_6(:,:)   = 0.
  zso4_7(:,:)   = 0.
  
 ! print *, 'entering m7'
 ! print *, 'paernl: ', paernl(6950,1,:)
 ! print *, 'paerml: ', paerml(6950,1,:)
  !
  !--- 1) Calculation of particle properties under ambient conditions: -----
  !
  !--- 1.1) Calculate mean particle mass for all modes 
  !         and dry radius and density for the insoluble modes.
  !
 ! print *, 'm7_averageproperties' 
 ! print *, 'paernl: ', paernl(6950,1,:)
 ! print *, 'paerml: ', paerml(6950,1,:)
  CALL m7_averageproperties(kproma, kbdim, klev, paernl, paerml, zttn, pm6rp, prhop)
  !
  !--- 1.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  !
  !         Sulfate particles:
  !
  ! print *, 'm7_equiz'
  ! print *, 'paernl: ', paernl(6950,1,:)
  !print *, 'paerml: ', paerml(6950,1,:)
  
  CALL m7_equiz(kproma,  kbdim, klev,   &
                papp1,   zttn,  ptp1,   &
                prelhum, pm6rp, pm6dry, &
                prhop,   pww,   paernl  )
  !         
  !         Mixed particles with sulfate, b/o carbon and dust: 
  !
   !print *, 'm7_equimix'
   !print *, 'paernl: ', paernl(6950,1,:)
  !print *, 'paerml: ', paerml(6950,1,:)
  CALL m7_equimix(kproma,  kbdim, klev,   &
                  papp1,   zttn,  ptp1,   &
                  prelhum, pm6rp, pm6dry, &
                  prhop,   pww,   paernl  )
  !
  !         Accumulation and coarse mode particles in presence of
  !         sea salt particles:
  !
  !print *, 'm7_equil'
  !print *, 'paernl: ', paernl(6950,1,:)
  !print *, 'paerml: ', paerml(6950,1,:)
  CALL m7_equil(kproma, kbdim,  klev,   prelhum, paerml, paernl, &
                pm6rp,  pm6dry, zhplus, pww,     prhop,  ptp1    )
  !
  !--- m7_box: Write initial particle distribution -------------------------
  IF (lofirst) THEN

     lofirst=.FALSE.

     !AJS: for debugging only ?
     !WRITE(16,FMT='(8(1x,e12.4))') 0.0, (pm6rp(844,1,i), i=1,7)
     !WRITE(20,FMT='(1x,f5.1,7x,7(1x,e12.4))') 0.0, (pm6dry(844,1,i), i=1,4), (pm6rp(1,1,i), i=5,7)

  END IF
   
  !--- END m7_box ----------------------------------------------------------
  !
  !--- 2) Calculate changes in aerosol mass and gas phase sulfate ----------
  !       due to sulfate condensation:
  !       No change in particle mass/number relationships.
  !
  !print *, 'm7_dgas'
  !print *, 'paernl: ', paernl(6950,1,:)
  !print *, 'paerml: ', paerml(6950,1,:)
  IF (lscond) CALL m7_dgas(kproma, kbdim,  klev,  pso4g,  paerml, paernl, &
                           ptp1,   papp1,  pm6rp, zso4_5, zso4_6, zso4_7, &
                           ptime) 

  
  !
  !
  !--- 3) Calculate change in particle number concentrations ---------------
  !       due to nucleation and coagulation:
  !       Change particle mass/number relationships.
  !
  ! JadB: Removed "If (lsnucl .OR. lscoag)".
  ! If only lscond is set, the m7_dnum is required for storing sulfuric acid on insoluble aerosols (making the soluble).
  ! Without m7_dnum, the sulfuric acid condensed on insoluble particles is turned into void.
  ! IF (lsnucl.OR.lscoag) CALL m7_dnum(kproma, kbdim,   klev,          &
  !                                    pso4g,  paerml,  paernl, ptp1,  &
  !                                    papp1,  prelhum, pm6rp,  prhop, &
  !                                    zso4_5, zso4_6,  zso4_7, ptime   )
  
  !print *, 'm7_dnum'
  !print *, 'paernl: ', paernl(6950,1,:)
  !print *, 'paerml: ', paerml(6950,1,:)
  CALL m7_dnum(kproma, kbdim,   klev,          &
               pso4g,  paerml,  paernl, ptp1,  &
               papp1,  prelhum, pm6rp,  prhop, &
               zso4_5, zso4_6,  zso4_7, ptime   )
               
   
  !
  !
  !--- 4) Recalculation of particle properties under ambient conditions: ---
  !
  !--- 4.1) Recalculate mean masses for all modes 
  !         and dry radius and density for the insoluble modes.
  !
  !print *, 'after m7_dnum'
  !print *, 'paernl: ', paernl(6950,1,:)
  !print *, 'paerml: ', paerml(6950,1,:)
  CALL m7_averageproperties(kproma, kbdim, klev, paernl, paerml, zttn, pm6rp, prhop)
  !
  !--- 4.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  !
  !         Sulfate particles:
  !
  
 
  CALL m7_equiz(kproma,  kbdim, klev,   &
                papp1,   zttn,  ptp1,   &
                prelhum, pm6rp, pm6dry, &
                prhop,   pww,   paernl  )
  !
  !         Mixed particles with sulfate, b/o carbon and dust:
  !
  
  CALL m7_equimix(kproma,  kbdim, klev,   &
                  papp1,   zttn,  ptp1,   &
                  prelhum, pm6rp, pm6dry, &
                  prhop,   pww,   paernl  )
  !
  !         Accumulation and coarse mode particles in presence of
  !         sea salt particles:  
  ! 
 
  CALL m7_equil(kproma, kbdim,  klev,   prelhum, paerml, paernl,  &
                pm6rp,  pm6dry, zhplus, pww,     prhop,  ptp1     )
  !
  !--- 5) Repartitition particles among the modes: -------------------------
  !
  IF (lscond.OR.lscoag) THEN
  
  
     CALL m7_dconc(kproma, kbdim, klev, paerml, paernl, pm6dry)
 
  END IF
  !
  !--- 6) Recalculation of particle properties under ambient conditions: ---
  !
  !--- 6.1) Calculate mean particle mass for all modes 
  !         and dry radius and density for the insoluble modes:
  !
  
  CALL m7_averageproperties(kproma, kbdim, klev, paernl, paerml, zttn, pm6rp, prhop)
  !
  !--- 6.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  !
  !         Sulfate particles:
  !
  !print *, 'm7_equiz'
  !print *, 'paernl: ', paernl(6950,1,:)
  !print *, 'paerml: ', paerml(6950,1,:)
  CALL m7_equiz(kproma,  kbdim, klev,   &
                papp1,   zttn,  ptp1,   &
                prelhum, pm6rp, pm6dry, &
                prhop,   pww,   paernl  )
  !
  !         Mixed particles with sulfate, b/o carbon and dust: 
  !
  
  CALL m7_equimix(kproma,  kbdim, klev,   &
                  papp1,   zttn,  ptp1,   &
                  prelhum, pm6rp, pm6dry, &
                  prhop,   pww,   paernl  )
  !
  
  !         Accumulation and coarse mode particles in presence of
  !         sea salt particles:
  !
  CALL m7_equil(kproma, kbdim,  klev, prelhum, paerml, paernl, &
                pm6rp,  pm6dry, zhplus, pww,   prhop,  ptp1    )

  
  !print *, 'return to interface'
  
END SUBROUTINE m7
