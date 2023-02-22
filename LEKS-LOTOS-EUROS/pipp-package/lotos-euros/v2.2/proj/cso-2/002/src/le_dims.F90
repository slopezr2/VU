!###############################################################################
!
! dimensions, global data, ...
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


module dims

  use GO, only : TDate

  ! for backwards compatibility ..
  use Indices, only : nspec, nspec_all

  implicit none
  
  
  ! --- const ------------------------------------


  ! --- types ------------------------------------

  ! NOW DEFINE TYPES FOR INPUT AND OUTPUT:
  type runspec
    !
    ! time fields (should move to LE_Time module)
    !
    ! the start time for the simulation (in GMT!)
    integer :: yy_s, mm_s, dd_s, hh_s
    ! the end time for the simulation (in GMT!)
    integer     ::  yy_e, mm_e, dd_e, hh_e
    ! base time for forecasts:
    integer     ::  yy_b, mm_b, dd_b, hh_b
    ! time structures
    type(TDate) ::  t_base, t_start, t_end
    ! forecast run ?
    logical ::  forecast_mode
    !
    ! flags
    !
    ! first time step
    logical :: first
    !
    ! restart stuff:
    logical             ::  restart
    character(len=256)  ::  restart_path
    character(len=256)  ::  restart_key
    logical             ::  restart_save
    real                ::  restart_save_dhour
    character(len=256)  ::  restart_save_path
    character(len=256)  ::  restart_save_key                                                                                     
    !
  end type runspec

  type outspec
    ! SCREEN OUTPUT:
    logical :: suppress
    ! DEBUG OUTPUT:
    logical :: debug_print           ! print debug info (after each process)
    integer :: debug_ix,debug_iy     ! cell where debug info is written 
                                     ! = (0,0) ==> min,max,average of grid is written
  end type outspec


  ! --- var --------------------------------------

  ! the number of grid cells:
  integer :: nx, ny, nz
  
  real :: soldec
  real, allocatable, dimension(:,:,:,:) :: emis_a

  real, allocatable, dimension(:,:,:) :: kz
  real, allocatable, dimension(:,:,:) :: tmass_old
  real, allocatable, dimension(:,:)   :: expcls, rk1rk3
  integer, allocatable, dimension(:,:) ::  kstabl
  integer, allocatable, dimension(:,:) :: lst

  ! sun: 
  real, allocatable, dimension(:,:) :: zenith ! deduced meteo field ?
  real, allocatable, dimension(:,:) :: coszen ! deduced meteo field

  ! pH (for in-cloud chemistry, in/below cloud scavenging)
  real, allocatable     ::  pH(:,:,:) ! deduced meteo field

  ! set threshold for substantial rain:
  ! units of rain is m/s, limit at 0.1 mm/hr:
  !                                                mm/hr / (mm/m) / (s/hr)
  real, parameter       ::  substantial_rain__mps = 0.1  / 1000.0 / 3600.0  ! m/s

  ! set threshold for substantial snowdepth to 1 mm water-equivalent
  real, parameter       ::  substantial_snowdepth = 0.001  ! m of water equivalent

  !! from FUB meteo:
  !real, allocatable, dimension(:,:,:) :: cvf
  
  !<<<

#ifdef with_hdiff
  real, allocatable, dimension(:,:,:) :: khx, khy
#endif

  ! OH concentration for methane/sulphur only runs:
  real, allocatable, dimension(:,:,:) :: OH

  real,allocatable, dimension(:) :: wash
  !real, allocatable, dimension(:,:,:) :: rb
  !real, allocatable, dimension(:,:) :: ra

  real, allocatable, dimension(:,:,:) :: OPS_mass

  ! the data type with the run flags
  type (runspec) :: runF
  type (outspec) :: outF


contains


  ! ====================================================================
  

  subroutine LE_Dims_Alloc

    ! this routine allocates the dimensions of the global arrays in lotos:
    allocate(lst(nx,ny)  )
    allocate(wash(nspec) )
!    allocate(ra(nx,ny) )
!    allocate(rb(nx,ny,nspec) )

    !EMISSIONS:
    allocate( emis_a(nx,ny,nz,nspec) )
    !allocate(flux_dust(nx,ny,3) )
    
    allocate(pH(nx,ny,nz) )              

    !<<<
     
    allocate(zenith(nx,ny) )
    allocate(coszen(nx,ny) ) 
#ifdef with_hdiff
    allocate(khx(0:nx,ny,nz) ) 
    allocate(khy(nx,0:ny,nz) )
#endif
    allocate(kz(nx,ny,nz) ) 
    allocate(kstabl(nx,ny) ) 
    allocate( expcls(nx,ny) )
    ! save the old concentration array
    allocate(tmass_old(nx,ny,nspec) ) 
    !allocate(aerh2o(nx,ny,nz) ) 
    allocate(rk1rk3(nx,ny) ) 
    allocate(oh(nx,ny,nz) ) 

    allocate(OPS_mass(nx,ny,nz) )

  end subroutine LE_Dims_Alloc

end module dims

