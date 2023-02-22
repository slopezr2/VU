!#######################################################################
!
!
!  Heights
!  -------
!
!     z50m ----- to compute landuse-independed windspeed
!
!                   -----------------  ztop : no impact of surface on 
!                    ^    ^     ^         concentration profile here
!                    |    |     |  
!                    |  +--+  +--+ 
!                    |  |  |  |  | 
!                    |  |  |  |Ra|      atm. resistance [zobs,ztop]
!                    |  |Ra|  |  |      atm. resistance [z0  ,ztop]
!                    |  |  |  +--+
!                    |  |  |    |   - -  zcanopy top (used for O3 deposition)
!      z10m -----    |  +--+    v
!                    |    |   - - - - -   zobs = z0_lu + dhobs
!                    |    v             ) dhobs : observation height offset to z0_lu
!          - - - -   |  - - - - - - - -   z0_lu : roughness length (landuse dependend)
!      z0  (         |                               wind speed profile is zero here
!    -------------  --------------------  z=0 surface (all heights are relative to orography)
!    /////////////  ////////////////////
!
!       overall    ;    landuse dependend
!
!         
!
!  Deposition velocity
!  -------------------
!
!  Deposition is modelled as a constant flux through a resistance,
!  with zero concentration in the bottom (inside vegetation),
!  and concentration c(ztop,t) at the top. Leads to a flux Phi0
!  that is constant in time and height:
!
!     d(ch)            1
!     ----- = Phi0 = ---- c(ztop,t)       (kg/m2/s)
!      dt            Rtot
!
!  with solution:
!
!     c(z,t) = exp( - (1/Rtot)/h dt )
!
!  Consider three resistances in line:
!   o atmospheric resistance Ra over [z0,htop]
!   o viscousious resistance Rb
!   o vegetation resistance Rc
!
!  Total resistance is then:
!
!       1                 1
!     ---- = ----------------------------            [m/s]
!     Rtot   Ra([z0,ztop]) + Rb(?) + Rc(?)
!
!  Reformulate this in terms of a deposition velocity:
!
!                          1
!     vd_tot = ----------------------------            [m/s]
!              Ra([z0,ztop]) + Rb(?) + Rc(?)
!
!  Concentration is assumed zero inside the vegetation,
!  and c(ztop,t) at z=ztop. Resistance model then leads to:
!  Flux is zero if concentration is zero, assume linear dependency:
!
!  Reformulate in terms of resistances
!
!                          1
!     vd_tot = ----------------------------            [m/s]
!              Ra([z0,ztop]) + Rb(?) + Rc(?)
!
!     c(z,t) = exp( - vd_tot/h dt )
!
!
!  Concentration profile
!  ---------------------
!
!  Describe flux using diffusion coefficient:
!
!    Phi(z) = - Kz(z) dc/dz
!
!  Assume same downward flux -Phi0 through every surface:
!
!     - Kz(z) dc/dz = -Phi0
!     [m2/s] [kg/m4]  [kg/m2/s]
!
!             dc/dz = Phi0 1/Kz(z)
!
!              c(z) = A + Phi0 int(s=-inf,z) 1/Kz(s) ds
!
!  Assume known concentration ctop at z=ztop :
!
!                                    c(ztop) =  ctop
!      A + Phi0 int(s=-inf,ztop) 1/Kz(s) ds  =  ctop
!
!  gives solution:
!
!     c(z) = A + Phi0 int(s=-inf,z) 1/Kz(s) ds
!
!          = ctop - Phi0 int(s=-inf,ztop) 1/Kz(s) ds + Phi0 int(s=-inf,z) 1/Kz(s) ds
!
!          = ctop - Phi0 int(s=z,ztop) 1/Kz(s) ds
!
!  or:
!
!    c(z) =  ctop - Phi0 Ra([z,ztop])
!
!  where Ra is the 'atmospheric resistance' between z and ztop :
!
!                    ztop
!                     (   1
!    Ra([z,ztop])  =  | ----- ds        [s/m]
!                     ) Kz(s)
!                    s=z
!
!  Using the formula for Phi0 this becomes:
!
!    c(z) = c(ztop,t) - vd_tot c(ztop,t) Ra([z,ztop])
!
!         = c(ztop,t) { 1 - vd_tot Ra([z,ztop]) }
!
!
!  Diffusion coefficient
!  ---------------------
!
!  Common formulation of diffusion coefficient:
!
!            kappa u* z
!    Kz(z) = ----------    [m2/s]
!            phi_h(z/L)
!
!  where:
!     kappa   :  Von Karman constant (0.35)
!     L       :  Monin-Obukhov length                   [m]
!     u*      :  friction velocity; also depends on L   [m/s]
!     phi_h   :  shear function for 'heat'
!
!
!  Monin-Obukhov length L
!  ----------------------
!
!  For the Monin-Obukhov length different implementations are used:
!
!   o Parameterisation of 1/L using Pasquill classes and z0
!     (Seinfeld and Pandis, 2006, section 16.4.4) .
!     First determine the Pasquill class given solar radiation, cloud cover, 
!     and surface wind.
!
!     A parameterisation for the Monin-Obukov length L
!     has been defined by Golder (1972),
!     see (Seinfeld and Pandis, 2006, section 16.4.4) :
!
!        1/L = a + b log(z0)
!
!     where z0 is the roughness length (m) and a,b are coefficents 
!     given the Pasquill Stabillity Classes:
!
!                                              coefficients
!        Pasquill Stability Class            a            b
!        A (extremely unstable)           - 0.096       0.029
!        B (moderately unstable)          - 0.037       0.029
!        C (slightly unstable)            - 0.002       0.018
!        D (neutral)                        0.000       0.000
!        E (slightly stable)                0.004     - 0.018
!        F (moderately stable)              0.035     - 0.036
!    
!     The stability class is defined in:
!
!     --------------------------------------------------------------------------------------
!                                            daytime                     nighttime
!                                        --------------------------  -------------------------
!     cloud fraction                           1-7/8         8/8      8/8     4-7/8    1-3/8
!                               #)       -----------------  -------  ------  -------   -------
!     "incoming solar radiation"  (W/m2) >700 350-700 <350
!                                        ---- ------- ----
!     windspeed (m/s)
!       0-2                               A     A-B     B      D       D        F*)       F*)
!       2-3                              A-B     B      C      D       D        E         F
!       3-5                               B     B-C     C      D       D        D         E
!       5-6                               C     C-D     D      D       D        D         D
!       6-                                C      D      D      D       D        D         D
!     --------------------------------------------------------------------------------------
!       *) In the original table, these fields are filled with '-' ;
!          no idea how this should be interpreted ...
!       #) We assume that we can use ECMWF field 'ssrd' (surface solar radiation downward, W/m2)
!          for this. There is also an ECMWF field 'ssr'; probably 'ssrd' is what reaches
!          the surface, and ssr is the fraction that is scatterd back directly plus the thermal
!          radiance.
!     --------------------------------------------------------------------------------------
!
!   o Definition using surface heat flux:
!
!       ...
!
!
!  Friction velocity u*
!  --------------------
!
!  The friction velocity u* follows from an integration of the Monin-Obukhov similarity:
! 
!                  z
!            u*   / phi_m(s/L)         u*
!    u(z) = ----  | ---------- ds  = ----- f_m(z0,z)
!           kappa /     s            kappa
!                s=z0
!
!  where phi_m is the shear function for momentum transport 
!  (phi_m is also called a 'Businger function')
!
!                    {   1 + beta_m  (z/L)            ,  z/L >= 0   (stable)
!       phi_m(z/L) = {
!                    { ( 1 - gamma_m (z/L) )^(-1/4)   ,  z/L < 0    (unstable)
!
!  where
!
!       beta_m  =  4.7
!       gamma_m = 15.0
!
!  and f_m the corresponding 'stability function' defined by:
!
!                   z
!                   / phi_m(s/L)   
!     f_m(z0,z) =   | ---------- ds  =  F_m(z) - F_m(z0)
!                   /     s      
!                  s=z0
!
!  For the above defined shear function phi_m the integrant F evaluates to (approximates to)
!  (see Jacobson 2005, 8.4, eq. 8.31)
!
!     F_m(z)   =  ln(z) + beta_m z/L                                             ,  z/L >= 0   (stable)
!
!     F_m(z)   ~  ln( [1/phi_m(z/L)-1][1/phi_m(z/L)+1] ) + 2 atan(1/phi_m(z/L))  ,  z/L < 0   (unstable)
!
!  Evaluate these using the windspeed at z=10m to get u* :
!
!      u(z10m) = u* / kappa f_m(z0,z10m) 
!      u10m    = u* / kappa [ F_m(z10m) - F_m(z0) ]
!
!      u*  =  kappa u10m / [ F_m(z10m) - F_m(z0) ]
!
!  NOTE: In current implementation, first an overall (landuse independend) u*
!  is computed based on u10m using the previous formula.
!  Reasoning: u10m comes from the meteo model and is an average too.
!  With this u*, an u50m at z=z50m is computed from the wind speed profile.
!  Then u(ztop) is used to compute a land-use dependend friction velocity:
!      u*_lu = kappa u50m / [ F_m(z50m) - F_m(z0_lu) ]
!
!
!  Atmospheric resistance
!  ----------------------
!
!  The shear function for 'heat' requires Prandtl number Pr :
!
!                     { Pr + beta_h (z/L)                 ,  z/L >= 0   (stable)
!       phi_h(z/L) =  {
!                     { Pr ( 1 - gamma_h (z/L) )^(-1/2)   ,  z/L < 0    (unstable)
!
!  where
!
!      Pr   : Prandtl number, seems to be 0.74 for kappa=0.35
!
!  and with the 'stability function' defined by:
!
!                   z
!                   ( phi_h(s/L)   
!     f_h(z0,z) =   | ---------- ds  =  F_h(z) - F_h(z0)
!                   )     s      
!                  s=z0
!
!  For the above defined shear function phi_h the integrant F evaluates to (approximates to):
!  (see Jacobson 2005, $8.4, eq. 8.38)
!
!     F_h(z)   =  ln(z) + gamma_h z/L    ,  z/L >= 0   (stable)
!
!     F_h(z)   ~ ....                    ,  z/L < 0   (unstable)
!
!  Now evaluate Ra :
!
!                    ztop                      ztop
!                     (   1              1      ( phi_h(s)    
!     Ra([z,ztop]) =  | ----- ds   =  --------  | -------- ds 
!                     ) Kz(s)         kappa u*  )    s        
!                    s=z                       s=z
!
!                 1      
!           =  --------  [ F_h(ztop) - F_h(z) ]
!              kappa u* 
!
!                   1   
!           =  ------------ [ F_m(z10m) - F_m(z0) ][ F_h(ztop) - F_h(z) ]
!              kappa^2 u10m
!
!  We seem to use z=z0, thus resistance over [z0,ztop]
!
!
!  References
!  ----------
!
!    Jacobson, 2005
!    
!    TvN, note on the resistance modelling in LE, see LE techdoc.
!
!    J.H. Seinfeld, S.N. Pandis
!    Atmospheric Chemistry and Physics; From Air Pollution to Climate Change; Second Edition
!    2006
!    Sections 16.4.3-4
!
!
!  History
!  -------
!
!  2010-12-10 Astrid Manders
!    Added land use dependent deposition (based on a version by ECJH).
!    Canopy height (depoparameters) should be set to zero, 
!    otherwise results are inconsistent with iKz.
!    Reference height is the top of the 25 m layer
!
!  2011-02-21 Arjo Segers
!    Set minimum canopy height ('ctop' in module depoparameters) to 'z0_lu',
!    otherwise problems with negative range [z0,zcanopytop] .
!    Added comment in top.
!    Added templates for shear and stability functions.
!
!  2012-01 Arjo Segers, TNO
!    Compute Ra([zobs,zref50]) per landuse, since used in mix2ground in combination
!    with landuse depended vd.
!    Re-structured the different heights used in the varios formula;
!    now assume a stack   z0+dhobs < htop ,
!    with a canopy top  z0 <= zcanopytop < htop
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Stability

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Stability_Init, LE_Stability_Done
  public  ::  monin_ustar
  public  ::  Atmospheric_Resistance
  public  ::  f_h_stability

  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'Stability'

  ! maximum heigth at which concentration profile starts
  ! so this profile is valid for z0 < z < max(20.0, min(L, htopmx) )
  real, parameter :: htopmx = 50.0 ! layer height for nz=3, for nz=4 see code

  real, parameter :: zref   = 10.0        ! zref = 10 since wsurf = U10
  real, parameter :: zref50 = 50.0        ! zref50 = 50 m is introduced as
                                          !   the reference height at which
                                          !   the logaritmic velocity profile
                                          !   is assumed to be independent
                                          !   upon land use class.
  real, parameter :: a1=0.004349, a2=0.00324, &     ! a1 to b3 are fit parameter
                     b1=0.5034, b2=0.2310, b3=0.0325
                     
  real, parameter :: dhobs = 2.5          ! old name: "hrfub"

! --- PARAMETERS beta, gamma, kappa, pr_t shifted to Binas Routine

!  ! turbulent Prandtl number, value corresponds to kappa=0.35
!  ! to be consistent throughout the code, we should use kappa=0.4
!  ! constant describing the dimensionless temperature profile in the surface layer under unstable conditions
!  real, parameter :: pr_t = 0.74 
!
!  real, parameter :: kappa  = 0.35
!  ! stability constants:
!  real, parameter ::  beta_m  =  4.7   ! mass, stable
!  real, parameter ::  gamma_m = 15.0   ! mass, unstable
!  real, parameter ::  beta_h  =  4.7   ! heat, stable
!  real, parameter ::  gamma_h =  9.0   ! heat, unstable

  ! --- var --------------------------------------
                                                                                
  real, allocatable   :: u_zref50(:,:)


contains


  ! ====================================================================
                                                                                

  subroutine LE_Stability_Init( status )
  
    use dims   , only : nx, ny
    use LE_Data, only : LE_Data_Enable
    
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Stability_Init'
    
    ! --- begin ----------------------------------

    ! storage:
    allocate( u_zref50(nx,ny) )

    ! enable data:
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'blh', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'occ', status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Stability_Init
  
  
  ! ***
  
  
  subroutine LE_Stability_Done( status )

    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Stability_Done'
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate(u_zref50)
   
    ! ok
    status = 0
    
  end subroutine LE_Stability_Done
  
  
  ! ====================================================================


  ! routine computes ustar and the monin obukov length based on
  ! meteo variables

  subroutine monin_ustar( lfirst, status )

    use LE_Logging, only : ident2
    use Binas, only  : kappa_stab, pr_t
    use Binas, only  : gamma_h, gamma_m,beta_h, beta_m
    use dims, only   : outF
    use dims, only   : nx, ny, nz
    use dims, only   : kz
    use dims, only   : coszen
    use dims, only   : ustar
    use dims, only   : monin, kstabl,expcls
    
    use LE_Landuse_Data, only : z0m, z0h, z0m_lu, z0h_lu, zcanopy_lu, z0dust_emis_lu
    
    use LE_Landuse_Data, only : Ra_h_htop_zobs_lu
    use LE_Landuse_Data, only : Ra_h_htop_z0_lu
    use LE_Landuse_Data, only : Ra_h_htop_zcanopytop_lu
    use LE_Landuse_Data, only : Ra_h_zcanopytop_z0_lu
    
    use LE_Landuse_Data, only : ustar_lu, u_zcanopytop_lu
    use LE_Landuse_Data, only : ustar_lu_dust_emis
    use LE_Landuse_Data, only : nlu, zcanopytop, Ld, lu_name

    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out ---------------------------------

    logical, intent(in)   ::  lfirst
    integer, intent(out)  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Stability_Init'
    
    integer ::  i, j, ilu
    real    ::  elinv, windmg 
    real    ::  ce, cw, s, ssqd, se
    
    real    ::  fm, fh
    
    real    ::  zobs, htop

    ! meteo data:
    real, pointer               ::  h_m  (:,:,:)   ! (lon,lat,nz)
    real, pointer               ::  occ  (:,:,:)   ! (lon,lat,nz)
    real, pointer               ::  wsurf(:,:,:)   ! (lon,lat,1)  

    ! --- begin ----------------------------
    
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'occ', occ, status, check_units ='1' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'wsurf', wsurf, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)

    !if (.NOT. outF%suppress) print *, ident2,'<constructing ustar and monin-obuhkov length>'

    ustar_lu              = 0.
    u_zcanopytop_lu       = 0.
    
    Ra_h_htop_z0_lu         = 0.
    Ra_h_htop_zobs_lu       = 0.
    Ra_h_htop_zcanopytop_lu = 0.
    Ra_h_zcanopytop_z0_lu   = 0.
    
    do i=1,nx
    do j=1,ny

      ! instantaneaous z0-values that also account for sea surface rougness
      ! are stored in 'z0_lu(i,j,ilu)'

      ! first determine the exposure class for the present cell (surface)
      call exposure( expcls(i,j), occ(i,j,1), wsurf(i,j,1), kstabl(i,j), coszen(i,j), lfirst )
      ! calculate the stability function s, from the wind speed and
      ! the exposure class

      ce = expcls(i,j)
      windmg = max( 0.01, wsurf(i,j,1) )
      CW=min(0.5*windmg, 4.0)    ! moet er wel u10 in voor stab class?
      ! set elinv small if s=0
      elinv=1.0e-6
      ! skip calculation if ce is zero
      if (abs(ce) > 0.0) then
          S=-0.5*(3.0-CW+abs(ce))*sign(1.0,ce)
          if(abs(S) > 0.0) then
             ! CALCULATE 1/(MONIN-OBUKHV LENGTH) --elinv
             SSQD=S**2
             SE=B2*abs(S)-B3*SSQD-B1
             elinv=S*(A1+A2*SSQD)*z0m(i,j)**SE
          endif
      endif

      if ( (h_m(i,j,1) > 25.01) .or. (h_m(i,j,1) < 24.99)  ) then !trap numerical precision
        write (gol,'("height of lowest model layer is not 25 meter")'); call goErr
        write (gol,'("define what to do with htop to get concentration on observation height")'); call goErr
        write (gol,'("height is: ",f12.4," in grid cell: (",i0,",",i0,")")') h_m(i,j,1), i, j; call goErr
        TRACEBACK; status=1; return
      else 
        htop = 25.0
      end if

      ! CALCULATE THE FRICTION VELOCITY, ustar(i,j), and THE SPATIALLY VARYING
      ! COMPONENT OF THE SURFACE DEPOSITION RATE, ftop(i,j), ACCORDING TO THE
      ! STABILITY FUNCTION (Jacobsen-2005, eq 8.31)
      
      fm = f_m_stability( zref, z0m(i,j), elinv )      
      ustar(i,j) = calc_ustar( kappa_stab, windmg, fm )
      
      fm = f_m_stability( zref50, z0m(i,j), elinv )     
      u_zref50(i,j) = ustar(i,j) * fm / kappa_stab

      ! loop over landuse classes:
      do ilu = 1, nlu
        
        ! landuse specific ustar, via stability function for momentum
        fm = f_m_stability( zref50, z0m_lu(i,j,ilu), elinv )
        ustar_lu(i,j,ilu) = calc_ustar(kappa_stab, u_zref50(i,j), fm)
        
        if ( z0dust_emis_lu(i,j,ilu) < 0.0 ) then
          ustar_lu_dust_emis(i,j,ilu) = -999.0 
        else
          ! landuse specific ustar for dust emissions (roughness parameters for local scale used)
          fm = f_m_stability( zref50, z0dust_emis_lu(i,j,ilu), elinv )
          ustar_lu_dust_emis(i,j,ilu) = calc_ustar(kappa_stab, u_zref50(i,j), fm)
        end if
        
        ! windspeed at canopy top
        ! (for calculation of Rb for ozone after McNaughton and Van der Hurk (1995)) ;
        ! only applied for leaf-dimension Ld(ilu) > 0.0 ,
        ! this is used to pre-select landuses, since for water the z0m is computed online
        ! given wind stress and therefore could exceed the fixed canopy top:
        if ( Ld(ilu) > 0.0 ) then
          ! check ...
          if ( zcanopytop(ilu) < z0m_lu(i,j,ilu) ) then
            write (gol,'("found canopy top lower than z0:")'); call goErr
            write (gol,'("  cell        : ",2i4)') i, j; call goErr
            write (gol,'("  depac class : ",i6," (",a,")")') ilu, trim(lu_name(ilu)); call goErr
            write (gol,'("  canopy top  : ",f12.6)') zcanopytop(ilu); call goErr
            write (gol,'("  z0m         : ",f12.6)') z0m_lu(i,j,ilu); call goErr
            TRACEBACK; status=1; return
          end if
          ! fill:             
          fm = f_m_stability( zcanopytop(ilu), z0m_lu(i,j,ilu), elinv )
          u_zcanopytop_lu(i,j,ilu) = ustar_lu(i,j,ilu) * fm / kappa_stab
        else
          ! no leaves, no canopy ...
          u_zcanopytop_lu(i,j,ilu) = 0.0
        end if

        ! check ...
        if ( (htop                < zcanopy_lu(i,j,ilu)      ) .or. &
             (htop                <     z0h_lu(i,j,ilu)      ) .or. &
             (htop                <     z0h_lu(i,j,ilu)+dhobs) .or. &
             (zcanopy_lu(i,j,ilu) <     z0h_lu(i,j,ilu)      )      ) then
          write (gol,'("found unrealistic heights for stability:")'); call goErr
          write (gol,'("  cell        : ",2i4)') i, j; call goErr
          write (gol,'("  depac class : ",i6," (",a,")")') ilu, trim(lu_name(ilu)); call goErr
          write (gol,'("  htop        : ",f12.6," (should exceed zcanopy, z0h, and z0h+dhobs)")') htop; call goErr
          write (gol,'("  zcanopy     : ",f12.6," (should exceed z0h)")') zcanopy_lu(i,j,ilu); call goErr
          write (gol,'("  z0h+dhobs   : ",f12.6)') z0h_lu(i,j,ilu)+dhobs; call goErr
          write (gol,'("  z0h         : ",f12.6)') z0h_lu(i,j,ilu); call goErr
          TRACEBACK; status=1; return
        end if

        ! Calculate aerodynamic resistance (Ra) from evaluation height (htop) to canopy top (zcanopytop)  
        ! note that the model is assumed to calculate relative to displacement height + z0
        fh = f_h_stability( htop, zcanopy_lu(i,j,ilu), elinv )
        Ra_h_htop_zcanopytop_lu(i,j,ilu)  = atmospheric_resistance( kappa_stab, ustar_lu(i,j,ilu), fh )
        
        ! Ra from evaluation height to z0 (used for deposition ), based on stability for heat
        fh = f_h_stability( htop, z0h_lu(i,j,ilu), elinv )        
        Ra_h_htop_z0_lu(i,j,ilu) = atmospheric_resistance(kappa_stab, ustar_lu(i,j,ilu), fh )
        
        ! Ra from evaluation height to zobs (used for mix2ground (conc-sfc files) ), based on stability for heat
        fh = f_h_stability( htop, z0h_lu(i,j,ilu)+dhobs, elinv )
        Ra_h_htop_zobs_lu(i,j,ilu) = atmospheric_resistance( kappa_stab, ustar_lu(i,j,ilu), fh )
        
        ! Ra from canopytop to z0, based on stability for heat ;
        ! the (over water) dynamically computed z0h might exceed the  
        fh =  f_h_stability( zcanopy_lu(i,j,ilu), z0h_lu(i,j,ilu), elinv )
        Ra_h_zcanopytop_z0_lu(i,j,ilu)  = atmospheric_resistance( kappa_stab, ustar_lu(i,j,ilu), fh )
        
      end do
        
      ! fill:
      monin(i,j) = elinv       

     ! ???? shouldn't this be monin(i,j)= 1/elinv ?? dd jan 2009 ECJH
     ! in subroutine comp_kz the monin(i,j) definition is used as elinv. it would 
     ! be logical if the name is changed, especially in the case that monin is used 
     ! outside stability.F90.... dd 4 feb 2009 ECJH 
     ! outside stability.F90: it is used correctly in depos.F90.... dd 4 feb 2009 ECJH 


    enddo  ! i
    enddo  ! j
    
    if (nz > 3) then
      call comp_kz( status )
      IF_NOTOK_RETURN(status=1)
    else
      ! we only have the mixing layer as one layer
      do i=1,nx
      do j=1,ny
         !if (monin(i,j) > 0) then
         !  kz(j,j,1)    = 0.1
         !  kz(i,j,2:nz) = 0.01
         !else
          !! lower Kz for low mixing heights (h in km, such that kz=1 if h=1 and 0.1 if 50.0)
          !kz(i,j,1)    = min(1.0, max(0.0, (h(i,j,1)-0.05)/0.95 ))
           kz(i,j,1)    = 1.0
           kz(i,j,2:nz) = 0.1
         !endif
      enddo
      enddo
    endif
    
    ! ok
    status = 0

  end subroutine monin_ustar
  
  
  ! ***
 
  
  !
  !  'Stability function' for mass defined by:
  !
  !                   z
  !                   / phi_m(s/L)   
  !     f_m(z0,z) =   | ---------- ds  =  F_m(z) - F_m(z0)
  !                   /     s      
  !                  s=z0
  !
  ! where 'phi_m' is the shear function for momentum transport ('Businger function'):
  !
  !                    {   1 + beta_m  (z/L)            ,  z/L >= 0   (stable)
  !       phi_m(z/L) = {
  !                    { ( 1 - gamma_m (z/L) )^(-1/4)   ,  z/L < 0    (unstable)
  !
  !  For the above defined shear function phi_m the integrant F evaluates to (approximates to)
  !  (see Jacobson 2005, 8.4, eq. 8.31)
  !
  !     F_m(z)   =  ln(z) + beta_m z/L                                             ,  z/L >= 0   (stable)
  !
  !     F_m(z)   ~  ln( [1/phi_m(z/L)-1][1/phi_m(z/L)+1] ) + 2 atan(1/phi_m(z/L))  ,  z/L < 0   (unstable)
  !
  !  Here define 'F_m(z)', which is the same as 'f_m(0,z)'
  !
          
  elemental function f_m_stability( z, znul, monin_inv ) result (f_m)
  
    use Binas, only : beta_m, gamma_m
    
    ! --- in/out ---------------------------------

    real, intent(in) :: z
    real, intent(in) :: znul
    real, intent(in) :: monin_inv
    real             :: f_m  ! output variable
     
    ! --- const ----------------------------------
    
    ! arbitrary small number to trap numerical errors
    real, parameter  ::  eps = 1.0e-6
    ! threshold for interval [delta,0.0]
    real, parameter  ::  delta = 1.0 - (1.0+eps)**4

    ! --- local ----------------------------------
    
    real  ::  phi_mz
    real  ::  phi_mznul
    
    ! --- begin ----------------------------------
    
    if ( monin_inv >= 0.0 ) then     ! stable case
      
      ! take stable case from formula (8.31, Jacobsen-2005)
      f_m = log(z/znul) + beta_m * (z-znul) * monin_inv

    else if ( monin_inv < 0.0 ) then   ! unstable case
      
      ! trap numerical errors for monin_inv too close to zero;
      ! use that z >= z0 :
      !
      !          phi_mz                    - 1 <        eps
      !     (1-gamma_m*z*monin_inv)**0.25  - 1 <        eps
      !      1-gamma_m*z*monin_inv             <     (1+eps)**4
      !        gamma_m*z*monin_inv             > 1 - (1+eps)**4
      !
      if ( gamma_m*znul*monin_inv > delta ) then
        
        ! take neutral case from formula (8.31, Jacobsen-2005)
        f_m = log(z/znul)
        
      else
        
        ! take unstable case from formula (8.31, Jacobsen-2005)
        phi_mz    = ( 1.0 - gamma_m*z   *monin_inv )**0.25
        phi_mznul = ( 1.0 - gamma_m*znul*monin_inv )**0.25
        f_m = ( log( ( phi_mz   -1 ) / (phi_mz    + 1 ) ) + 2*atan( phi_mz   ) ) - &
              ( log( ( phi_mznul-1 ) / (phi_mznul + 1 ) ) + 2*atan( phi_mznul) )

      end if
    
    end if
  
  end function f_m_stability  
        
 
  ! ***
  
  
  ! Stability function for heat 
  
  elemental function f_h_stability( z, znul, monin_inv ) result (f_h)
  
    use Binas, only : beta_h, gamma_h, pr_t
    
    ! --- in/out ---------------------------------
    
    real, intent(in) :: znul
    real, intent(in) :: z
    real, intent(in) :: monin_inv
    real             :: f_h
     
    ! --- const ----------------------------------
    
    ! arbitrary small number to trap numerical errors
    real, parameter  ::  eps = 1.0e-6
    ! threshold:
    real, parameter  ::  delta = 1.0 - (1.0+eps)**2
  
    ! --- local ----------------------------------
    
    real  ::  phi_hz
    real  ::  phi_hznul
    
    ! --- begin ----------------------------------
    
    if ( monin_inv >= 0.0 ) then     ! stable case
      
      ! take stable case from formula (8.38, Jacobsen-2005)
      f_h = pr_t * log(z/znul) + beta_h * (z-znul) * monin_inv
      
    else if ( monin_inv < 0.0 ) then   ! unstable case
      
      ! trap numerical errors for monin_inv too close to zero;
      ! use that z >= z0 :
      !
      !        phi_hz0                  -  1  <      eps
      !  (1-gamma_h*z0*monin_inv)**0.5  -  1  <      eps
      !   1-gamma_h*z0*monin_inv              <   (1+eps)**2
      !     gamma_h*z0*monin_inv              > 1-(1+eps)**2
      !
      if ( gamma_h*znul*monin_inv > delta ) then
        
        ! take neutral case from formula (8.38, Jacobsen-2005)
        f_h = pr_t * log(z/znul)
        
      else

        ! take unstable case from formula (8.38, Jacobsen-2005)
        phi_hz    = ( 1.0 - gamma_h*z   *monin_inv )**0.5
        phi_hznul = ( 1.0 - gamma_h*znul*monin_inv )**0.5
        f_h = pr_t * ( log( ( phi_hz   -1 ) / (phi_hz    + 1 ) ) - &
                       log( ( phi_hznul-1 ) / (phi_hznul + 1 ) )   )

      end if

    end if
  
  end function f_h_stability  
  
  
  ! ***
  
  
  ! function to calculate u* from stability for mass
  elemental function calc_ustar( kappa_stab, windmg, fm ) result ( ust )
           
    ! --- arguments -----------

    real, intent(in)    ::  kappa_stab
    real, intent(in)    ::  windmg
    real, intent(in)    ::  fm
    real                ::  ust
    
    !------------------------------
    
    ust = kappa_stab * windmg / fm
    
  end function calc_ustar  
  

  ! ***

  
  elemental function Atmospheric_Resistance( kappa_stab, ustar, fh ) result ( Ra )
  
    ! --- arguments -----------

    real, intent(in)  ::  kappa_stab
    real, intent(in)  ::  ustar
    real, intent(in)  ::  fh
    real              ::  Ra

    ! --- begin ---------------------------------
           
    Ra = fh / kappa_stab / ustar
    
  end function Atmospheric_Resistance
  
  
  ! ***
  

  subroutine exposure( expcls, cloud, ws, kstabl, zcoef, lfirst )

    real, intent(inout)     ::  expcls
    real, intent(in)        ::  cloud
    real, intent(in)        ::  ws
    integer, intent(out)    ::  kstabl
    real, intent(in)        ::  zcoef
    logical, intent(in)     ::  lfirst

    integer             ::  ice
    integer             ::  iws
    real                ::  exvary
    real                ::  exclas

    integer, dimension(6,6) :: stab_class

    ! lookup table to get kstabl from exposure class
    data stab_class /6,6,5,5,4,4, &   ! kstabl= stab_class(iws,ice) 
               6,5,4,4,4,4, &
               4,4,4,4,4,4, &
               2,3,3,3,4,4, &
               2,2,3,3,4,4, &
               1,2,2,2,3,3 /

    ! calculate exposure class

    if (cloud.gt.0.95) then     
      ! total overcast expcls=0.0 regardless
      exclas=0.0

    
    elseif (zcoef.ge.0.866) then      ! zcoef=cos(sza)
      ! solar zenith angle > 60 degrees  !RT 11-3-10 should be solar elevation angle > 60 deg and sza<30
      !exclas=3.0  
    !  make the exclas fluent as function of sza instead of stepwise
      exclas=2.5 + 0.5 * ((zcoef-0.866)/ (1-0.866))  ! sza=0 then exclas=3, sza =30 then exclas =2.5
      !if (cloud.GT.0.50) exclas=2.0
      if (cloud.GT.0.50) exclas=1.5 + 0.5 * ((zcoef-0.866)/ (1-0.866))

    elseif (zcoef.ge.0.574) then     
      ! solar zenith angle 15 < theta < 60   !RT 11-3-10 should be solar elevation angle between 35 and 60 deg, sza between 30 and 55
!      exclas=2.0
    !  make the exclas fluent as function of sza instead of stepwise
      exclas=1 + 1.5 * ((zcoef-0.574)/ (0.866-0.574))  ! sza=30 then exclas=2.5, sza =55 then exclas =1      
    !  if (cloud.GT.0.50) exclas=1.0
      if (cloud.GT.0.50) exclas=0.5 + 1 * ((zcoef-0.574)/ (0.866-0.574)) 

   
   elseif (zcoef.ge.0.259) then     
      ! 15 < theta < 35 degrees   !RT 11-3-10 should be solar elevation angle between 15 and 35 deg, sza between 55 and 75
      !exclas=0.0
    !  make the exclas fluent as function of sza instead of stepwise
      exclas=-1 + 2 * ((zcoef-0.259)/ (0.574-0.259))  ! sza=55 then exclas=1, sza =75 then exclas =-1      
      !if (cloud.GT.0.50) exclas=0.0
      if (cloud.GT.0.50) exclas=-0.5 + 1 * ((zcoef-0.259)/ (0.574-0.259))

    else     
       ! solar zenith angle < 15 degrees   !RT 11-3-10 should be solar elevation angle < 15 deg, sza > 75
      ! exclas=-2.0
      !  make the exclas fluent as function of sza instead of stepwise
      exclas=-2 + 1 * ((zcoef-0.0)/ (0.259-0.0))  ! sza=75 then exclas=1, sza =-1 then exclas =-2          
      !if (cloud.gt.0.60) exclas=-1.0
      if (cloud.gt.0.60) exclas=-1 + 0.5 * ((zcoef-0.0)/ (0.259-0.0))
    endif

    if (lfirst) expcls=exclas

    ! do not allow exposure class to vary more than
    ! 1/2 a class per time step
    exvary=exclas-expcls
    if (abs(exvary).le.0.5) then    
       expcls=exclas
    else 
       expcls=expcls+sign(0.5,exvary)
    endif

    ice=INT(expcls+3.1)
    iws=INT(ws)
    if(iws.EQ.0) iws=1
    if(iws.GT.6) iws=6

    ! now look in lookuptable for kstabl 
    kstabl= stab_class(iws,ice) 

    ! check for unstable at night
    if (zcoef.le.0.0 .and. kstabl.lt.4) kstabl=4

  end subroutine exposure
  
  
  ! ***


  ! routine computes kz based on u* and L

  subroutine comp_kz( status )

    use dims, only : nx, ny, nz
    use dims, only : kz
    use dims, only : ustar, monin

    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out ---------------------------

    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/comp_kz'
    
    ! --- local ----------------------------------

    !real :: helinv, hx, busfac
    !integer :: k, iz, iz_surf, ilayer
    integer :: i,j
    !real    ::  busfc, stab, wstar
    real    ::  z1, z2
    real    ::  hm
    real    ::  ufac
    real    ::  z1fac, z2fac
    real    ::  p
    real    ::  r12
    real    ::  alpha
    real    ::  qfac, q1fac, q2fac
    
    ! meteo data:
    real, pointer               ::  h_m(:,:,:)   ! (lon,lat,lev)
    real, pointer               ::  blh(:,:,:)   ! (lon,lat,lev)    

    ! --- begin ----------------------------

    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m')    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'blh', blh, status, check_units ='m')    
    IF_NOTOK_RETURN(status=1)
    
    do i = 1, nx
      do j = 1, ny

        ! at the mixing layer
        kz(i,j,2) = 1.0

        ! above mixing layer
        kz(i,j,3:nz) = 0.1
!    monin(i,j) as it is here, seems to be the inverse of the Monin-Obukhov length ?? Check!! dd jan 2009 ECJH

        ! both 'h' and 'blh' are in m:
        z1    = 0.5 * h_m(i,j,1)
        z2    = 0.5 * ( h_m(i,j,2) + h_m(i,j,1) )
        hm    = max( h_m(i,j,2), blh(i,j,1) )

        ufac  = ustar(i,j)*0.4/0.74
        if (monin(i,j) > 0.0) then
          z1fac = (1.0-z1/hm)**0.5
          z2fac = (1.0-z2/hm)**0.5
          P     = (1.0-z2fac)*(1.0+z1fac)/( (1.0-z1fac)*(1.0+z2fac) )
          r12   = (log(P) + 2*(1.0+6.3*hm*monin(i,j))*(1.0/z2fac-1.0/z1fac))/ufac
        else
          alpha=9.0*abs(monin(i,j))
          qfac =(1.0+alpha*hm)**0.5
          q1fac=(1.0+alpha*z1)**0.5
          q2fac=(1.0+alpha*z2)**0.5
          r12  = (log((qfac+q2fac)*(qfac-q1fac)/((qfac+q1fac)*(qfac-q2fac)))/qfac - &
                 log((1.0+q2fac)*(1.0-q1fac)/((1.0+q1fac)*(1.0-q2fac))) )/ufac
        endif
        ! kz changed by correction factor for dc/dz by MvL and MSP
        !  kz(i,j,1) = kz(i,j,1) * 0.5*h(i,j,2)/h(i,j,1)
        kz(i,j,1) = (z2-z1)/R12
        ! new expression for Kz to be checked, what is the reference here?
        ! TvN
        
        ! truncate Kz between surface and mixing layer to [1,100] :
        kz(i,j,1) = min( max( 1.0, kz(i,j,1) ), 100.0 )
        
      end do
    end do

    ! ok
    status = 0

  end subroutine comp_kz

end module LE_Stability

