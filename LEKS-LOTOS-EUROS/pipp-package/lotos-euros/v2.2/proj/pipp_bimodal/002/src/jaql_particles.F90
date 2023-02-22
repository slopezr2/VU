!#######################################################################
!
! JAQL_SFC_Dust
!   Joint Air Quality Library
!   Particle modelling
!
!
!#######################################################################


module JAQL_Particles

  implicit none
  

  ! --- in/out --------------------------------
  
  private
  
  public  ::  Particle_Volume_Rp
  public  ::  Particle_Volume_D
  public  ::  Particle_CrossSection_D
  public  ::  Normal_Number_Density
  public  ::  Particle_dNp_dDp, Particle_dNp_dlnDp

  public  ::  ParticleMode_AMeanDiameter_Moment

contains


  ! ====================================================================
  
  ! 
  ! Particle volume as function of radius:
  !
  !    4/3 pi R**3 
  !
  
  elemental function Particle_Volume_Rp( Rp ) result( Vp )
  
    use Binas, only : pi

    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  Rp  ! radius (m)
    real                  ::  Vp  ! volume (m3)
    
    ! --- begin ----------------------------------
    
    Vp = 4.0/3.0 * pi * Rp**3
    
  end function Particle_Volume_Rp
  
  ! 
  ! Particle volume as function of diameter:
  !
  !    4/3 pi R**3 = 4/3 pi (D/2)**3 = 1/6 pi D**3
  !
  
  elemental function Particle_Volume_D( Dp ) result( Vp )
  
    use Binas, only : pi

    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  Dp  ! diameter (m)
    real                  ::  Vp  ! volume (m3)
    
    ! --- begin ----------------------------------
    
    Vp = pi * Dp**3 / 6.0
    
  end function Particle_Volume_D
  
  ! 
  ! Particle cross-section as function of diameter:
  !
  !    pi R**2 = pi (D/2)**3 = 1/4 pi D**2
  !
  
  elemental function Particle_CrossSection_D( Dp ) result( Cp )
  
    use Binas, only : pi

    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  Dp  ! diameter (m)
    real                  ::  Cp  ! cross-section (m2)
    
    ! --- begin ----------------------------------
    
    Cp = pi * Dp**2 / 4.0
    
  end function Particle_CrossSection_D


  ! ====================================================================


  ! Densitiy of normal distribution:
  !
  !    dN/dx  ~  G(x,mu,sigma)
  !
  
  elemental function Normal_Number_Density( x, mu, sigma ) result( dN_dx )
  
    use Binas, only : pi
    
    ! --- in/out ---------------------------------
    
    real, intent(in)        ::  x
    real, intent(in)        ::  mu
    real, intent(in)        ::  sigma
    real                    ::  dN_dx
    
    ! --- begin ----------------------------------
    
    dN_dx = exp( - 0.5 * (x-mu)**2 / sigma**2 ) / (sqrt(2*pi)*sigma)
    
  end function Normal_Number_Density

  !
  ! Number of particles with diameter Dp follows log-normal distribution:
  !
  !      dNp
  !    -------  ~  G(ln(Dp),ln(gmean),ln(gsigma))
  !    dln(Dp)
  !
  
  elemental function Particle_dNp_dlnDp( Dp, mu, sigma ) result( dNp_dlnDp )
  
    ! --- in/out ---------------------------------
    
    real, intent(in)        ::  Dp
    real, intent(in)        ::  mu
    real, intent(in)        ::  sigma
    real                    ::  dNp_dlnDp
    
    ! --- begin ----------------------------------
    
    dNp_dlnDp = Normal_Number_Density( log(Dp), log(mu), log(sigma) )
    
  end function Particle_dNp_dlnDp
  
  
  ! For integral over all diameters it is sometimes usefull to have:
  !
  !     dNp       dNp    dln(Dp)  
  !     ---  =  ------- -------  =  G(ln(Dp),ln(gmean),ln(gsigma)) / Dp
  !     dDp     dln(Dp)   dDp     
  !  
  
  elemental function Particle_dNp_dDp( Dp, mu, sigma ) result( dNp_dDp )
  
    ! --- in/out ---------------------------------
    
    real, intent(in)        ::  Dp
    real, intent(in)        ::  mu
    real, intent(in)        ::  sigma
    real                    ::  dNp_dDp
    
    ! --- begin ----------------------------------
    
    dNp_dDp = Particle_dNp_dlnDp( Dp, mu, sigma ) / Dp
    
  end function Particle_dNp_dDp


  ! ====================================================================


!  !
!  ! Mean diameter:
!  !
!  !   Dp = Dpg exp( 0.5 * ln(sigma_g)**2 )
!  !
!
!  elemental function ParticleDiameter_Median2Mean( Dpg, sigma ) result( Dp )
!
!    ! --- in/out ---------------------------------
!    
!    real, intent(in)        ::  Dpg
!    real, intent(in)        ::  sigma
!    real                    ::  Dp
!    
!    ! --- begin ----------------------------------
!    
!    Dp = Dpg * exp( 0.5 * log(sigma) )
!    
!  end function ParticleDiameter_Median2Mean
!
!  ! *
!  
!  elemental function ParticleDiameter_Mean2Median( Dp, sigma ) result( Dpg )
!
!    ! --- in/out ---------------------------------
!    
!    real, intent(in)        ::  Dp
!    real, intent(in)        ::  sigma
!    real                    ::  Dpg
!    
!    ! --- begin ----------------------------------
!    
!    Dpg = Dp * exp( -0.5 * log(sigma) )
!    
!  end function ParticleDiameter_Mean2Median
!  
!  
!  !
!  ! Number median to surface or volume median:
!  !
!  !   DpgS = Dpg exp( 2 (ln(sigma))**2 )
!  !   DpgV = Dpg exp( 3 (ln(sigma))**2 )
!  !
!
!  elemental function ParticleDiameter_Median2SurfaceMedian( Dpg, sigma ) result( DpgS )
!
!    ! --- in/out ---------------------------------
!    
!    real, intent(in)        ::  Dpg
!    real, intent(in)        ::  sigma
!    real                    ::  DpgS
!    
!    ! --- begin ----------------------------------
!    
!    DpgV = Dpg * exp( 2 * log(sigma)**2 )
!    
!  end function ParticleDiameter_Median2SurfaceMedian
!
!  ! *
!  
!  elemental function ParticleDiameter_Median2VolumeMedian( Dpg, sigma ) result( DpgV )
!
!    ! --- in/out ---------------------------------
!    
!    real, intent(in)        ::  Dpg
!    real, intent(in)        ::  sigma
!    real                    ::  DpgV
!    
!    ! --- begin ----------------------------------
!    
!    DpgV = Dpg * exp( 3 * log(sigma)**2 )
!    
!  end function ParticleDiameter_Median2VolumeMedian


  !
  ! Particle mode with log-normal distribution:
  !   ln(D)  ~  N( lnDg, ln2sigmag )
  ! Arithmetic mean (average) of diameter moment:
  !
  !   E[D**m] =  Dam**m
  !
  ! where the moment arithmetic mean diameter is:
  !
  !   Dam = Dg exp( m/2 ln(sigmag)**2 )
  !

  elemental function ParticleMode_AMeanDiameter_Moment( Dg, sigmag, m ) result( Dam )

    ! --- in/out ---------------------------------
    
    real, intent(in)        ::  Dg
    real, intent(in)        ::  sigmag
    integer, intent(in)     ::  m
    real                    ::  Dam
    
    ! --- begin ----------------------------------
    
    Dam = Dg * exp( 0.5*m * log(sigmag)**2 )
    
  end function ParticleMode_AMeanDiameter_Moment
    


end module JAQL_Particles
