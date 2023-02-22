!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
!###############################################################################

module LEKF_Dims

  use dims, only : nx, ny, nz, nspec, nspec_all
  
  ! maximum number of modis pixels:
  use LE_Output_MODIS, only : modis_maxpix => maxp
  
  public
  
  !! number of noise fields
  !integer, parameter :: nnoise  = 4
  
  !! size of the noise history
  !!integer, parameter :: nhist   = 8
  !integer, parameter :: nhist   = 2
  
  ! maximum number of (scalar) measurements applied for allocation of hourly measurements
  ! and 24h measurements; the latter have no history
  integer, parameter :: maxmeas = 150
  
  !! maximum number of components for 24h averages
  !integer, parameter :: maxcomp24 = 1

!>>> use dynamic allocation  
!#ifdef with_kf_meas_omi_trc
!  !! maximum number of omi_no2 measurements:
!  !integer, parameter :: maxomi = 30e3
!  ! synthetic tropomi data
!  integer, parameter :: maxomi = 400e3
!#endif
!<<<
  
  ! the number of AOD fields (different wavelengths)
  !integer, parameter :: nz_aod = 1           ! total column
  !integer, parameter :: nz_aod = RAL_nlay    ! profiles
  !integer            :: nz_aod = -1   ! init with dummy value, set in kf_meas_aod
 
  ! the number of modes
!#ifdef with_many_modes
!  integer, parameter :: maxmodes =  120
!#else
  integer, parameter :: maxmodes =  15
!#endif
  !! the number of modes
  !integer, parameter :: maxmodesnoise = maxmodes + nhist*nnoise
 

end module LEKF_Dims
