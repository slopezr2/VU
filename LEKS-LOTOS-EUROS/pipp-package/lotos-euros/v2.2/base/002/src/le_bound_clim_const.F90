!###############################################################################
!
! NAME
!   LE_Bound_Clim_Const
!
! DESCRIPTION
!
!   Set some constant boundary conditions.
!
! HISTORY
!   Taken from form 'boundary' module in 'bound.F90' .
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

module LE_Bound_Clim_Const

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  
  ! --- in/out -----------------------------
  
  private
  
  public    ::  LE_Bound_Clim_Const_Setup
  
  
  ! --- const ------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_Clim_Const'



contains


  ! ==================================================
  
  ! fills the boundary conditions
  !
  !   spec_filled(1:nspec)  : set to .true. if bc is filled, unchanged otherwise
  !

  subroutine LE_Bound_Clim_Const_Setup( spec_filled )

    use dims   , only : runF
    use dims   , only : nspec
    use indices
    use LE_Data_Common , only : nlev
    use LE_Bound_Common, only : bc_west, bc_east, bc_south, bc_north, caloft

    ! --- in/out ----------------------------
    
    logical, intent(inout)      ::  spec_filled(nspec)
    
    ! --- local -----------------------------
    
    integer     ::  ispec, itr

    ! --- begin -----------------------------

    ! loop over all species:
    do ispec = 1, nspec
      select case ( ispec )
        case ( ispec_so4a_f )
          bc_east (:,:,ispec) = 0.6
          bc_west (:,:,ispec) = 0.6
          bc_north(:,:,ispec) = 0.6
          bc_south(:,:,ispec) = 0.6
          caloft  (:,:,nlev+1,ispec) = 0.6
          spec_filled(ispec) = .true.
        case ( ispec_nh4a_f )
          bc_east (:,:,ispec) = 0.6 / 48.0 * 18.0
          bc_west (:,:,ispec) = 0.6 / 48.0 * 18.0
          bc_north(:,:,ispec) = 0.6 / 48.0 * 18.0
          bc_south(:,:,ispec) = 0.6 / 48.0 * 18.0
          caloft  (:,:,nlev+1,ispec) = 0.6 / 48.0 * 18.0
          spec_filled(ispec) = .true.
        case ( ispec_nh3 )
          bc_east (:,:,ispec) = 0.2
          bc_west (:,:,ispec) = 0.2
          bc_north(:,:,ispec) = 0.2
          bc_south(:,:,ispec) = 0.2
          caloft  (:,:,nlev+1,ispec) = 0.2
          spec_filled(ispec) = .true.
        case ( ispec_pb_c )
          bc_east (:,:,ispec) = 0.9*2e-3
          bc_west (:,:,ispec) = 0.9*1e-3
          bc_north(:,:,ispec) = 0.9*0.6e-3
          bc_south(:,:,ispec) = 0.9*1.5e-3
          caloft  (:,:,nlev+1,ispec) = 0.9*0.1e-3
          spec_filled(ispec) = .true.
        case ( ispec_pb_f )
          bc_east (:,:,ispec) = 0.1*2e-3
          bc_west (:,:,ispec) = 0.1*1e-3
          bc_north(:,:,ispec) = 0.1*0.6e-3
          bc_south(:,:,ispec) = 0.1*1.5e-3
          caloft  (:,:,nlev+1,ispec) = 0.1*0.1e-3
          spec_filled(ispec) = .true.
        case ( ispec_cd_c )
          bc_east (:,:,ispec) = 0.9*0.05e-3
          bc_west (:,:,ispec) = 0.9*0.03e-3
          bc_north(:,:,ispec) = 0.9*0.02e-3
          bc_south(:,:,ispec) = 0.9*0.04e-3
          caloft  (:,:,nlev+1,ispec) = 0.9*0.003e-3
          spec_filled(ispec) = .true.
        case ( ispec_cd_f )
          bc_east (:,:,ispec) = 0.1*0.05e-3
          bc_west (:,:,ispec) = 0.1*0.03e-3
          bc_north(:,:,ispec) = 0.1*0.02e-3
          bc_south(:,:,ispec) = 0.1*0.04e-3
          caloft  (:,:,nlev+1,ispec) = 0.1*0.003e-3
          spec_filled(ispec) = .true.
      end select
    end do

  end subroutine LE_Bound_Clim_Const_Setup
  

end module LE_Bound_Clim_Const
