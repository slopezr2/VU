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

  use GO,   only : gol, goPr, goErr
  use dims, only : nx,ny,nz

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Stability_Init, LE_Stability_Done
  public  ::  monin_ustar
  public  ::  dhobs

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

    ! enable roughness lengths:
    call LE_Data_Enable( 'z0m_lu', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'z0m', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'z0h_lu', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'zcanopy_lu', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'z0dust_lu', status )
    IF_NOTOK_RETURN(status=1)

    ! enable diffusion data:
    call LE_Data_Enable( 'monin_inv', status )  
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_Enable( 'ustar', status )      
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

  subroutine monin_ustar(lfirst, tmid, status )

    use LE_Logging, only : ident2
    use Binas, only  : kappa_stab, pr_t
    use Binas, only  : gamma_h, gamma_m,beta_h, beta_m
    use dims, only   : outF
    use dims, only   : nx, ny, nz
    use dims, only   : substantial_snowdepth
    
    use LE_Landuse_Data, only : Ra_h_htop_zobs_lu
    use LE_Landuse_Data, only : Ra_h_htop_z0_lu
    use LE_Landuse_Data, only : Ra_h_htop_zcanopytop_lu
    use LE_Landuse_Data, only : Ra_h_zcanopytop_z0_lu
    
    use LE_Landuse_Data, only : ustar_lu, u_zcanopytop_lu
    use LE_Landuse_Data, only : ustar_lu_dust_emis
    use LE_Landuse_Data, only : nlu, zcanopytop, Ld, lu_name

    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer
    use go, only : TDate
    
    use JAQL_Stability, only : f_m_stability, f_h_stability
    use JAQL_Stability, only : calc_ustar
    use JAQL_Stability, only : Atmospheric_Resistance

    ! --- in/out ---------------------------------

    logical, intent(in)       ::  lfirst
    type(TDate), intent(in)   ::  tmid
    integer, intent(out)      ::  status

    ! --- const -------------------------------    

    character(len=*), parameter ::  rname = mname//'/monin_ustar'

    !-------local----------------------------------

    integer ::  i, j, ilu
    real    ::  elinv, windmg 
    real    ::  ce, cw, s, ssqd, se
    
    real    ::  fm, fh
    
    real    ::  zobs, htop
    
    ! meteo data:
    real, pointer               ::  z0m_lu(:,:,:)   ! (lon,lat,:)
    real, pointer               ::  z0m   (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  z0h_lu(:,:,:)   ! (lon,lat,:)
    real, pointer               ::  zcanopy_lu(:,:,:)   ! (lon,lat,:)
    real, pointer               ::  z0dust_lu(:,:,:)    ! (lon,lat,:)
    real, pointer               ::  h_m   (:,:,:)   ! (lon,lat,nz)
    real, pointer               ::  occ   (:,:,:)   ! (lon,lat,nz)
    real, pointer               ::  wsurf (:,:,:)   ! (lon,lat,1)
    real, pointer               ::  ustar (:,:,:)   ! (lon,lat,1)    
    real, pointer               ::  monin_inv(:,:,:)   ! (lon,lat,1)        


    ! --- begin ----------------------------    

    ! point to data:
    call LE_Data_GetPointer( 'z0m_lu', z0m_lu, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'z0m', z0m, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'z0h_lu', z0h_lu, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'zcanopy_lu', zcanopy_lu, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'z0dust_lu', z0dust_lu, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'ustar', ustar, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'monin_inv', monin_inv, status, check_units ='1/m' )
    IF_NOTOK_RETURN(status=1)
    

    ustar_lu              = 0.
    u_zcanopytop_lu       = 0.
    
    Ra_h_htop_z0_lu         = 0.
    Ra_h_htop_zobs_lu       = 0.
    Ra_h_htop_zcanopytop_lu = 0.
    Ra_h_zcanopytop_z0_lu   = 0.
    
    do i=1,nx
     do j=1,ny

       elinv = monin_inv(i,j,1) ! [1/m]

      ! top of lowest layer:
      htop = h_m(i,j,1)
      ! check 
      if ( (htop < 10.0) .or. (htop > 60.0) ) then
        write (gol,'("expected height of lowest model layer ~25 (mixlayer approach) or ~20 (hybride levels),")'); call goErr
        write (gol,'("found ",f12.4," in grid cell (",i0,",",i0,");")') htop, i, j; call goErr
        write (gol,'("define appropriate htop to get concentration on observation height")'); call goErr
        TRACEBACK; status=1; return
      end if

       ! CALCULATE THE FRICTION VELOCITY, ustar(i,j), and THE SPATIALLY VARYING
       ! COMPONENT OF THE SURFACE DEPOSITION RATE, ftop(i,j), ACCORDING TO THE
       ! STABILITY FUNCTION (Jacobsen-2005, eq 8.31)

       fm = f_m_stability( zref50, z0m(i,j,1), elinv )     
       u_zref50(i,j) = ustar(i,j,1) * fm / kappa_stab

       ! loop over landuse classes:
       do ilu = 1, nlu

         ! landuse specific ustar, via stability function for momentum
         fm = f_m_stability( zref50, z0m_lu(i,j,ilu), elinv )
         ustar_lu(i,j,ilu) = calc_ustar(kappa_stab, u_zref50(i,j), fm)

         if ( z0dust_lu(i,j,ilu) < 0.0 ) then
           ustar_lu_dust_emis(i,j,ilu) = -999.0 
         else
           ! landuse specific ustar for dust emissions (roughness parameters for local scale used)
           fm = f_m_stability( zref50, z0dust_lu(i,j,ilu), elinv )
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

     enddo  ! i
    enddo  ! j
    
    ! ok
    status = 0

  end subroutine monin_ustar  
  
end module LE_Stability

