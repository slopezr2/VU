!#######################################################################
!
! JAQL_Humidity
!   Joint Air Quality Library
!     Humidity tools
!
!
!#######################################################################


module JAQL_Humidity

  implicit none
  

  ! --- in/out --------------------------------
  
  private
  
  public  ::  WaterVaporPartialPressure
  public  ::  SaturationVaporPressure
  public  ::  RelativeHumidity_from_SpecificHumidity


contains


  ! ====================================================================
  
  !
  !  Compute water vapor partial pressure (e_w)
  !  given specific humidity Q [(kg water)/(kg air)].
  !
  !  Use that gas law for volume V with temperature T
  !  holds for the total mixture as well as the water part:
  !
  !    R T / V = p_air / n_air = p_water / n_water
  !
  !  thus:
  !
  !    p_water = p_air n_water / n_air
  !
  !  Use:
  !    n_air =   m_air   /        xm_air
  !            [kg air]  /  [(kg air)/(mole air)]
  !  and:
  !    n_water =  m_air * Q  /     xm_water
  !              [kg water]  /  [(kg water)/(mole water)]
  !  thus:
  !    p_water = p_air Q / (xm_water/xm_air)
  !
  
  elemental function WaterVaporPartialPressure( Q, p ) result ( p_w )
  
    use Binas, only : xm_air, xm_H2O
    
    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  Q   ! specific humidity [(kg water)/(kg air)]
    real, intent(in)      ::  p   ! air pressure [Pa]
    real                  ::  p_w  ! water vapor partial pressure [Pa]
    
    ! --- const ----------------------------------
    
    ! mole mass ratio:
    real, parameter  ::  eps = xm_H2O / xm_air  !  ~ 0.622
    
    ! --- begin ----------------------------------
    
    ! partial pressure of water vapor:
    p_w = p * Q / eps
    
  end function WaterVaporPartialPressure
  

  !
  !  Saturation vapor pressure.
  !  From (Stull 1988, eq. 7.5.2d):
  !
  !      e_sat = p0 exp( 17.67 * (T-273.16) / (T-29.66) )     [Pa]
  !
  !  where:
  !      p0 = 611.2 [Pa]   : reference pressure
  !
  !  Arguments:
  !      T  [K]  : air temperature
  !  Result:
  !      e_sat_w  [Pa]  : saturation vapor pressure
  !
  !  References:
  !      Roland B. Stull, 1988
  !      An introduction to boundary layer meteorology.
  !
  
  elemental function SaturationVaporPressure( T ) result( e_sat_w )

    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  T        ! temperature [K]
    real                  ::  e_sat_w  ! saturation vapor pressure  [Pa]
    
    ! --- const ----------------------------------
    
    ! base pressure:
    real, parameter  ::  p0 = 611.2   ! [Pa]
    
    ! --- begin ----------------------------------
    
    ! saturation vapor pressure:
    e_sat_w = p0 * exp( 17.67 * (T-273.16) / (T-29.66) )   ! [Pa]
    
  end function SaturationVaporPressure


  !
  !  Relative humidity RH [%] is by definition:
  !
  !           e_w            # water vapor partial pressure
  !    Rh = -------- * 100
  !         e_sat_w          # saturation vapor pressure
  !
  !  Compute here from:
  !    Q specific humidity [(kg water)/(kg air)]
  !    T temperature [K]
  !    P pressure [Pa]
  !
  
  elemental function RelativeHumidity_from_SpecificHumidity( Q, T, p ) result( Rh )
  
    ! --- in/out ---------------------------------
    
    real, intent(in)      ::  Q   ! specific humidity [(kg water)/(kg air)]
    real, intent(in)      ::  T   ! temperature [K]
    real, intent(in)      ::  p   ! air pressure [Pa]
    real                  ::  Rh  ! relative humidity [%]
    
    ! --- begin ----------------------------------
    
    ! relative humidity:
    Rh = WaterVaporPartialPressure( Q, p ) / SaturationVaporPressure( T ) * 100.0
    
  end function RelativeHumidity_from_SpecificHumidity


end module JAQL_Humidity
