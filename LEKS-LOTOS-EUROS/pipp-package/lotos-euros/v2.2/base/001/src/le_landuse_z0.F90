!###############################################################################
!
! Landuse : surface roughness
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


module LE_LandUse_z0

  use GO, only : gol, goPr, goErr
 
  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  Setup_z0m
  public  ::  Setup_z0h
  public  ::  Setup_zcanopy
  public  ::  Setup_z0dust
  public  ::  Average_z0
  

  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_LandUse_z0'


  ! --- var -------------------------------


contains


  ! ========================================================================
  

  !
  ! Fill surface roughness for mass.
  !

  subroutine Setup_z0m( z0m_lu, lai_lu, wsurf, sd, status )
  
    use LE_Landuse_Data, only : z0m_depac
    use LE_Landuse_Data, only : z0m_snow_surface
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_Data, only : ilu_water_sea, ilu_arable, ilu_desert
    use LE_Landuse_Data, only : ilu_grass, ilu_semi_natural_veg
    use Dims           , only : substantial_snowdepth

    ! --- in/out ---------------------------------
    
    real, intent(out)           ::  z0m_lu(nlu)
    real, intent(in)            ::  lai_lu(nlu)
    real, intent(in)            ::  wsurf   ! windspeed (at 10m?) [m/s]
    real, intent(in)            ::  sd      ! snow depth [m]
    integer, intent(out)        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Setup_z0m'
    
    ! --- local ------------------------------------------
    
    integer     ::  iclass

    ! --- begin ------------------------------------------

    ! fill land-use depended z0 ;
    ! loop over landuse classes:
    do iclass = 1, nlu

      ! default:
      z0m_lu(iclass) = z0m_depac(iclass)

      ! check whether sea or inland water 
      ! for wind dependent rougnhess length z0
      if ( iclass == ilu_water_sea ) then

        z0m_lu(iclass) = (0.0617/30.) * wsurf**1.52
        ! set minimum:
        if ( z0m_lu(iclass) < z0m_depac(iclass) ) then
          z0m_lu(iclass) = z0m_depac(iclass)
        end if

      else if ( iclass == ilu_arable ) then

        if ( lai_lu(iclass) < 0.01 ) then
          !z0m_lu(iclass)  = z0m_lu(ilu_desert)
          z0m_lu(iclass)  = z0m_depac(ilu_desert)
        end if

      else if ( iclass == ilu_grass .or. iclass == ilu_semi_natural_veg ) then

        if ( sd > substantial_snowdepth ) then
          z0m_lu(iclass) = z0m_snow_surface
        end if

      end if  ! class

    end do  ! landuse classes

    ! ok
    status = 0
    
  end subroutine Setup_z0m
   
  ! *
  
  !
  ! Fill surface roughness for heat.
  !

  subroutine Setup_z0h( z0h_lu, lai_lu, wsurf, sd, status )
  
    use LE_Landuse_Data, only : z0h_depac
    use LE_Landuse_Data, only : z0h_snow_surface
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_Data, only : ilu_water_sea, ilu_arable, ilu_desert
    use LE_Landuse_Data, only : ilu_grass, ilu_semi_natural_veg
    use Dims           , only : substantial_snowdepth

    ! --- in/out ---------------------------------
    
    real, intent(out)           ::  z0h_lu(nlu)
    real, intent(in)            ::  lai_lu(nlu)
    real, intent(in)            ::  wsurf   ! windspeed (at 10m?) [m/s]
    real, intent(in)            ::  sd      ! snow depth [m]
    integer, intent(out)        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Setup_z0h'
    
    ! --- local ------------------------------------------
    
    integer     ::  iclass

    ! --- begin ------------------------------------------

    ! fill land-use depended z0 ;
    ! loop over landuse classes:
    do iclass = 1, nlu

      ! default:
      z0h_lu(iclass) = z0h_depac(iclass)

      ! check whether sea or inland water 
      ! for wind dependent rougnhess length z0
      if ( iclass == ilu_water_sea ) then

        z0h_lu(iclass) = (0.0617/30.) * wsurf**1.52
        ! set minimum:
        if ( z0h_lu(iclass) < z0h_depac(iclass) ) then
          z0h_lu(iclass) = z0h_depac(iclass)
        end if

      else if ( iclass == ilu_arable ) then

        if ( lai_lu(iclass) < 0.01 ) then
          !z0h_lu(iclass)  = z0h_lu(ilu_desert)
          z0h_lu(iclass)  = z0h_depac(ilu_desert)
        end if

      else if ( iclass == ilu_grass .or. iclass == ilu_semi_natural_veg ) then

        if ( sd > substantial_snowdepth ) then
          z0h_lu(iclass) = z0h_snow_surface
        end if

      end if  ! class

    end do  ! landuse classes

    ! ok
    status = 0
    
  end subroutine Setup_z0h
   
  ! *
  
  !
  ! Fill canopy height
  !

  subroutine Setup_zcanopy( zcanopy_lu, z0h_lu, lai_lu, status )
  
    use LE_Landuse_Data, only : z0h_depac
    use LE_Landuse_Data, only : zcanopytop
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_Data, only : ilu_water_sea, ilu_arable, ilu_desert
  
    ! --- in/out ---------------------------------
    
    real, intent(out)           ::  zcanopy_lu(nlu)
    real, intent(in)            ::  z0h_lu(nlu)
    real, intent(in)            ::  lai_lu(nlu)
    integer, intent(out)        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Setup_zcanopy'
    
    ! --- local ------------------------------------------
    
    integer     ::  iclass

    ! --- begin ------------------------------------------

    ! fill land-use depended z0 ;
    ! loop over landuse classes:
    do iclass = 1, nlu

      ! default:
      zcanopy_lu(iclass) = zcanopytop(iclass)

      ! check whether sea or inland water 
      ! for wind dependent rougnhess length z0
      if ( iclass == ilu_water_sea ) then

        ! dummy:
        zcanopy_lu(iclass) = 3.33 * z0h_lu(iclass)
        ! set minimum:
        if ( z0h_lu(iclass) < z0h_depac(iclass) ) then
          zcanopy_lu(iclass) = 3.33 * z0h_lu(iclass)
        end if

      else if ( iclass == ilu_arable ) then

        if ( lai_lu(iclass) < 0.01 ) then
          !zcanopy_lu(iclass) = zcanopy_lu(ilu_desert)
          zcanopy_lu(iclass) = zcanopytop(ilu_desert)
        end if

      end if  ! class

    end do  ! landuse classes

    ! ok
    status = 0
    
  end subroutine Setup_zcanopy
   
  ! *
  
  !
  ! Fill z0 for wind blown dust emission
  !

  subroutine Setup_z0dust( z0dust_lu, lai_lu, status )
  
    use LE_Landuse_Data, only : z0dust_emis
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_Data, only : ilu_arable
  
    ! --- in/out ---------------------------------
    
    real, intent(out)           ::  z0dust_lu(nlu)
    real, intent(in)            ::  lai_lu(nlu)
    integer, intent(out)        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Setup_z0dust'
    
    ! --- local ------------------------------------------
    
    integer     ::  iclass

    ! --- begin ------------------------------------------

    ! fill land-use depended z0 ;
    ! loop over landuse classes:
    do iclass = 1, nlu

      ! default:
      z0dust_lu(iclass) = z0dust_emis(iclass)

      if ( iclass == ilu_arable ) then

        if ( lai_lu(iclass) < 0.01 ) then
          !! Skip production of dust from arable land for now
          !!z0dust_emis_lu(iclass) = z0dust_emis_lu(ilu_desert)
          !z0dust_emis_lu(iclass) = z0dust_emis(ilu_desert)
        end if

      end if  ! class

    end do  ! landuse classes

    ! ok
    status = 0
    
  end subroutine Setup_z0dust
  
  
  ! ***
  

  subroutine Average_z0( z0_lu, lu_fracs, z0, status )
  
    use LE_Landuse_Data, only : nlu
  
    ! --- in/out ---------------------------------
    
    real, intent(in)            ::  z0_lu(nlu)
    real, intent(in)            ::  lu_fracs(nlu)
    real, intent(out)           ::  z0
    integer, intent(out)        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Average_z0'
    
    real, parameter ::  kappa = 0.41
    real, parameter ::  zr = 60
    
    ! --- local ------------------------------------------
    
    integer     ::  iclass
    
    real        ::  cd
    real        ::  cd_ave
    

    ! --- begin ------------------------------------------

    !! overall z0 ; log weighted over landuse classes:
    !! for momentum
    !z0m = 0.0
    !do iclass = 1, nlu
    !  z0m = z0m + lu_fracs(:,:,iclass) * log( z0m_lu(:,:,iclass) )
    !end do
    !z0m = exp( z0m )

    !
    ! Calculate average z0 via drag coefficients following:
    !   Agterberg, R & Wieringa, J. (1989)
    !   Mesoscale terrain roughness mapping of The Netherlands.
    !
    ! init average:
    cd_ave = 0.0
    ! loop over landuse clases:
    do iclass = 1, nlu
      cd = ( kappa / log(zr/z0_lu(iclass)) )**2
      cd_ave = cd_ave + lu_fracs(iclass) * cd
    end do
    ! average roughness length:
    z0 = zr * exp( - kappa / sqrt(cd_ave) )

    ! ok
    status = 0
    
  end subroutine Average_z0
  
  
end module LE_LandUse_z0
