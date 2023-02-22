!###############################################################################
!
! ISORROPIA2 - secondary inorganic aerosol physics
!
!
! OPENMP IMPLEMENTATION
!
!   Not trivial due to use of common blocks in the original f77 code.
!
!   Example code:
!     http://map.nasa.gov/GEOS_CHEM_f90toHTML/html_code/
!   In particular:
!     o for the main module:
!         http://map.nasa.gov/GEOS_CHEM_f90toHTML/html_code/src/isoropia_mod.f.html
!     o see the 'THREADPRIVATE' statements at the end of:
!         http://map.nasa.gov/GEOS_CHEM_f90toHTML/html_code/inc/isoropia.h.html
!
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

module LE_SIA_Isorropia2

  use GO, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out ----------------------------------------------
  
  private
  
  public    ::  LE_SIA_Isorropia2_Init, LE_SIA_Isorropia2_Done
  public    ::  LE_SIA_Isorropia2_Step
  
  
  ! --- const ------------------------------------------------
  
  ! module name:
  character(len=*), parameter :: mname = 'LE_SIA_Isorropia2_Step'

  ! some parameters
  real, parameter :: mw_so4 = 96.0
  real, parameter :: mw_nh4 = 18.0
  real, parameter :: mw_no3 = 62.0
  real, parameter :: mw_h2o = 18.0  ! g/mol
  
  ! working precission:
  integer, parameter  ::  wp = 8
  
  ! ISORROPIA2 indices for 'wi' and 'wt' arrays
  ! with gas+aerosol concentrations:
  !integer, parameter    ::  igpa_Na  = 1   ! sodium
  integer, parameter    ::  igpa_SO4 = 2   ! sulfate
  integer, parameter    ::  igpa_NH4 = 3   ! ammonium
  integer, parameter    ::  igpa_NO3 = 4   ! nitrate
  !integer, parameter    ::  igpa_Cl  = 5   ! chloride
  !integer, parameter    ::  igpa_Ca  = 6   ! calcium
  !integer, parameter    ::  igpa_K   = 7   ! potassium
  !integer, parameter    ::  igpa_Mg  = 8   ! magnesium
  
  ! ISORROPIA2 indices for 'gas' array
  ! with gas concentrations:
  integer, parameter    ::  igas_NH3  = 1
  integer, parameter    ::  igas_HNO3 = 2
  !integer, parameter    ::  igas_HCl  = 3


contains


  ! ====================================================================
  

  subroutine LE_SIA_Isorropia2_Init( status )

    ! --- in/out ------------------------------

    integer, intent(out)            ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_SIA_Isorropia2_Init'

    ! --- external -------------------------------
    
    external  ::  LE_SIA_Isorropia2_IsoCom_Init

    ! --- local ------------------------------

    ! --- begin -------------------------------
    
    ! init data:
    call LE_SIA_Isorropia2_IsoCom_Init()
    
    ! ok
    status = 0

  end subroutine LE_SIA_Isorropia2_Init


  ! ***


  subroutine LE_SIA_Isorropia2_Done( status )

    ! --- in/out ------------------------------

    integer, intent(out)  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_SIA_Isorropia2_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! ok
    status = 0

  end subroutine LE_SIA_Isorropia2_Done


  ! ***


  subroutine LE_SIA_Isorropia2_Step( c, ah2o, status )

    use dims   , only : nx, ny, nz
    use Indices, only : nspec
    use Indices, only : i_HNO3, i_NH3
    use Indices, only : i_SO4a_f, i_NO3a_f, i_NH4a_f
    !use LE_Budget_Chem, only : budget_chem
    !use dims   , only :  h

    use LE_Data      , only : LE_Data_GetPointer

    ! --- in/out ---------------------------------

    real, intent(inout)             ::  c(:,:,:,:)    ! ppb or ug/m3
    real, intent(out)               ::  ah2o(:,:,:)   ! ug/m3
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_SIA_Isorropia2_Step'

    ! --- external -------------------------------
    
    external  ::  isorropia
    
    ! --- local ----------------------------------
    
    ! the stuff needed for isorropia
    real(wp)            ::  wi(8)
    real(wp)            ::  rhi
    real(wp)            ::  tk
    real(wp)            ::  cntrl(2)
    real(wp)            ::  wt(8)
    real(wp)            ::  gas(3)
    real(wp)            ::  aerliq(15)
    real(wp)            ::  aersld(19)
    character(len=15)   ::  scase
    real(wp)            ::  other(9)
    ! loop counters:
    integer             ::  i, j, k
    !! layer thinkness:
    !real                ::  dh

    ! meteo data:
    real, pointer        ::   temp(:,:,:)   ! (lon,lat,nz)
    real, pointer        ::   rh  (:,:,:)   ! (lon,lat,nz)

    ! --- begin ---------------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 't', temp, status )    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'rh', rh, status )    
    IF_NOTOK_RETURN(status=1)

    ! control parameters for isorropia
    cntrl(1) = 0.0d0   ! forward problem; wi should contain gas+aerosol
    cntrl(2) = 1.0d0   ! aerosol is in only liquid state (metastable aerosol)

    ! loop over grid cells:
    do k = 1, nz
      !$OMP parallel &
      !$OMP   default ( none ) &
      !$OMP   private ( i, j ) &
      !$OMP   private ( wi, rhi, tk, wt, gas, aerliq, aersld, scase, other ) &
      !$OMP   shared ( nx, ny, nz ) &
      !$OMP   shared ( k ) &
      !$OMP   shared ( temp, rh ) &
      !$OMP   shared ( c, ah2o ) &
      !$OMP   shared ( i_nh3, i_hno3, i_so4a_f, i_no3a_f, i_nh4a_f ) &
      !$OMP   shared ( cntrl )
      
      ! Set accuracy parameters.
      ! Run time was found to be most sensitive to the 'eps' and 'epsact' values,
      ! especially the later one strongly determines the speed.
      ! Will be copied to threadprivate common blocks, therefore need
      ! to set this within the parallel region.
      ! For negative values the default remains.
      ! Note that the eps numbers should be in double precission.
      call SetParm ( -1    , &    ! wftyp  : weighting algorithm
                     -1    , &    ! iacalc : (activity coeff. calc) 0=compute online, 1=lookup-tables (default)
                     1.0d-2, &    ! eps    : convergence criterion (default 1.0d-6; DOUBLE!)
                     -1    , &    ! maxit  : max. number of iterations (default 100)
                     -1    , &    ! nsweep : (activity coeff. calc) max. number of iterations (default 4)
                     5.0d-1, &    ! epsact : (activity coeff. calc) convergence criterion (default 5.0d-2; DOUBLE!)
                     -1    , &    ! ndiv   : 
                     -1      )    ! nadj   :
                     
      ! init input concentrations to zero, non-zero will be reset:
      wi = 0.0

      !$OMP   do
      do j = 1, ny
        do i = 1, nx

          tk  = temp(i,j,k)
          rhi = min(1.0, rh(i,j,k)/100.0)

          ! compute total concentrations in moles/m3
          wi(igpa_SO4) =                             c(i,j,k,i_so4a_f) / mw_so4 * 1.0e-6
          wi(igpa_NH4) = c(i,j,k,i_nh3 )*0.0409e-6 + c(i,j,k,i_nh4a_f) / mw_nh4 * 1.0e-6
          wi(igpa_NO3) = c(i,j,k,i_hno3)*0.0409e-6 + c(i,j,k,i_no3a_f) / mw_no3 * 1.0e-6
          
          ! input ready, compute aerosol concentrations
          call isorropia( wi, rhi, tk, cntrl, &                    ! in
                          wt, gas, aerliq, aersld, scase, other )  ! out

          ! Adjust gasfaseconcentrations(ppb) and aerosolconcentrations(microgram/m**3)
          c(i,j,k,i_so4a_f) = max( 0.0, real( wt(igpa_SO4)                  )*mw_so4*1.0e6 )
          c(i,j,k,i_nh4a_f) = max( 0.0, real( wt(igpa_NH4) - gas(igas_NH3 ) )*mw_nh4*1.0e6 )
          c(i,j,k,i_no3a_f) = max( 0.0, real( wt(igpa_NO3) - gas(igas_HNO3) )*mw_no3*1.0e6 )

          !! set dh (for budget calculations)
          !if (k==1) then
          !  dh = h(i,j,1)
          !else
          !  dh = h(i,j,k) - h(i,j,k-1)
          !endif
          !dh = dh/h(i,j,nz)

          !! budget for HNO3 loss/production from NO3
          !budget_chem(i,j,6) = budget_chem(i,j,6) + dh* &
          !            ( gas(2)/0.0409 * 1.0e6 - c(i,j,k,i_hno3) )

          c(i,j,k,i_nh3 )  = max( 0.0, real( gas(igas_NH3 ) )/0.0409 * 1.0e6 )
          c(i,j,k,i_hno3)  = max( 0.0, real( gas(igas_HNO3) )/0.0409 * 1.0e6 )

          ! set the aerosol water
          !                    mole/m3      g/mole   ug/g
          ah2o(i,j,k) = real( aerliq(8) ) * mw_h2o * 1.0e6   ! ug/m3

        end do  ! i
      end do  ! j
      !$OMP   end do
      !$OMP end parallel
    end do  ! k

    ! ok
    status = 0

  end subroutine LE_SIA_Isorropia2_Step


end module LE_SIA_Isorropia2
