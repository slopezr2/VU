!###############################################################################
!
! LE_Emis - LOTOS-EUROS emission routines
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

module LE_Emis_Bio_NO

  use GO, only : gol, goErr, goPr
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Emis_Bio_NO_Init, LE_Emis_Bio_NO_Done
  public  ::  LE_Emis_Bio_NO_Setup
  

  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Bio_NO'

contains


  ! ==============================================
  
  
  subroutine LE_Emis_Bio_NO_Init( status )

    use LE_Data, only : LE_Data_Enable
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Bio_NO_Init'
    
    ! --- begin ----------------------------------

    ! enable data:
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_Bio_NO_Init
  
  
  ! ***
  
  
  subroutine LE_Emis_Bio_NO_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Bio_NO_Done'
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_Bio_NO_Done
  
  
  ! ***


  ! this routine adds biogenic VOC emissions to the emissions array
  ! input: the emission array

  subroutine LE_Emis_Bio_NO_Setup( emis, status )

    use dims, only : runF
    use dims, only : nx, ny, nz, nspec
    use Binas, only : xm_N, T0
    use LE_Landuse_data, only : lu_fracs, nlu
    use LE_Landuse_data, only : A_bio_no
    use indices

    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out ---------------------------------
    
    real, intent(inout)       ::  emis(nx,ny,nz,nspec)
    integer, intent(out)      ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Bio_NO_Setup'
    
    
    ! --- local ---------------------------------
    
    real    :: Ta, Ts
    real    :: no_emis
    integer :: ix, iy, ilu
    real    :: delta_emis

    ! meteo data:
    real, pointer               ::  temp(:,:,:)   ! (lon,lat,1)
    real, pointer               ::  area(:,:,:)   ! (lon,lat,1)

    
    ! --- begin ----------------------------------
    
    ! point to meteo data:
    call LE_Data_GetPointer( 't', temp, status, check_units='K' )    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'area', area, status, check_units ='m2' )
    IF_NOTOK_RETURN(status=1)
    

    ! biogene nox emissions
    do ix=1,nx
      do iy=1,ny

        ! temperature (K):
        Ta = temp(ix,iy,1) - T0

        ! initialize 
        no_emis = 0.0

        ! loop over smiatek landuse classes
        do ilu = 1, nlu
            
            call get_temp_factor_no(ilu, Ta, Ts, status)
            IF_NOTOK_RETURN(status=1)

            ! formula given in ng N /m2/s, convert to  ( (mol N)/m2/s ) 
            ! factor 0.071 from Yienger and Levy 1995
            ! (mol N)/m2/s  =         +    kg/ng / (kg/(mol N)) * (ng N)/m2/s   *       1          * fraction of grid cell
            no_emis         = no_emis + ( (1e-12 / xm_N )       * A_bio_no(ilu) * exp(0.071*Ts) )  * lu_fracs(ix,iy,ilu)
        end do    

        ! set the emissions (ug/m2/hr) in the emission array in mol/min!
        ! mol / min = sec/min *  m2  * ( (mol N)/m2/s ) 
        delta_emis  = 60.0    * area(ix,iy,1) *     no_emis
        emis(ix,iy,1,i_no)  = emis(ix,iy,1,i_no) + delta_emis
                

      enddo
    enddo
  
    ! ok
    status = 0

  end subroutine LE_Emis_Bio_NO_Setup


  ! ***


  subroutine get_temp_factor_no(ilu, Ta, Ts, status)
  
    use LE_Landuse_Data, only : nlu
    use LE_Landuse_Data, only : temp_coeff_no, base_coeff_no

    ! --- in/out -----------------------
    integer, intent(in)         ::  ilu
    real, intent(in)            ::  Ta
    real, intent(out)           ::  Ts
    integer, intent(out)        ::  status

    ! --- const ------------------------

    character(len=*), parameter ::  rname = mname//'/get_temp_factor_no'

    ! --- begin ------------------------

    Ts = temp_coeff_no(ilu) * Ta + base_coeff_no(ilu)

    ! ok
    status = 0

  end subroutine get_temp_factor_no
  
end module
