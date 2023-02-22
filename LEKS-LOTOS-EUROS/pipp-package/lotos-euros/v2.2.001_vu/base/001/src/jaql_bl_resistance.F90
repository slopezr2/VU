!#######################################################################
!
! JAQL_BL_Resistance 
!   Joint Air Quality Library
!   Boundary Layer tools
!   Resistance tools
!
!
!#######################################################################


module JAQL_BL_Resistance

  implicit none
  

  ! --- in/out --------------------------------
  
  private

  public  ::  JAQL_BL_Rb
  public  ::  JAQL_BL_Rb_McNaughton_VdHurk
  

contains


  ! ===============================================================
  
  !
  ! Leaf-level quasi-laminar boundary layer resistance.
  ! Input:
  !   V(h)  : windspeed at development height
  !   diffr : ratio between diffusivity of heat and tracer,
  !           for example: 
  !              diffc_heat = 0.20e-4
  !              diffc_o3   = 0.13e-4
  !           then 
  !              diffr = diffc_heat / diffc_o3 = 1.54
  ! Output:
  !    Rb = 50.0 / V * diffr**0.67
  !
  elemental function JAQL_BL_Rb( V, diffr ) result( Rb )
  
    ! --- in/out ---------------------------------
    
    real, intent(in)        ::  V     ! windspeed at development height
    real, intent(in)        ::  diffr ! diffusivity ratio
    real                    ::  Rb
    
    ! --- const ----------------------------------

    ! --- begin -----------------------------------
    
    ! trap division by zero:
    if ( V > 0.01 ) then
      ! apply parameterization:
      Rb = 5.0 / V * diffr**0.67
    else
      ! set to extreme large value:
      Rb = 100.0
    end if
     
  end function JAQL_BL_Rb
  
  !
  ! Leaf-level quasi-laminar boundary layer resistance.
  ! Formulation by McNaughton and van der Hurk (1995)
  ! Input:
  !   Ld    : cross-wind leaf dimension
  !   V(h)  : windspeed at development height
  !   diffr : ratio between diffusivity of heat and tracer,
  !           for example: diffr_o3 = 1.3
  ! Output:
  !    Rb = diffr * 150 * sqrt{ Ld / V(h) }
  !
  elemental function JAQL_BL_Rb_McNaughton_VdHurk( Ld, V, diffr ) result( Rb )
  
    ! --- in/out ---------------------------------
    
    real, intent(in)        ::  Ld    ! cross-wind leaf dimension
    real, intent(in)        ::  V     ! windspeed at development height
    real, intent(in)        ::  diffr ! ratio between  in diffusivity between heat 
                                      ! and tracer; example: diffr_o3 = 1.3
    real                    ::  Rb
    
    ! --- const ----------------------------------

    real, parameter  ::  mcn = 150.0  ! unspecified constant
    
    ! --- begin -----------------------------------
    
    ! trap division by zero:
    if ( (V > 0.01) .and. (Ld > 0.0) ) then
      ! apply parameterization:
      Rb = diffr * mcn * sqrt( Ld / V )
    else
      ! set to extreme large value:
      Rb = 100.0
    end if
     
  end function JAQL_BL_Rb_McNaughton_VdHurk


end module JAQL_BL_Resistance

