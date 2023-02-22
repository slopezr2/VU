!###############################################################################
!
! LE_Emis_HM_Natural - heavy-metal emissions from natural sources
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') mname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_HM_Natural

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -------------------------------
  
  private
  
  public  ::  LE_Emis_HM_Natural_Init
  public  ::  get_natural_hm

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_emis_hm_natural'

contains


  ! ====================================================================
  
  subroutine LE_Emis_HM_Natural_Init( status )
    
    use LE_Data, only : LE_Data_Enable
    
    ! --- in/out ---
    
    integer, intent(out)    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_HM_Natural_SILAM_Init'
    
    ! --- begin ---

    ! enable data:
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'sd', status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine
  
  ! ***
  
  subroutine get_natural_hm( emis_a, status )

    use dims, only : nspec

    use LE_Data, only : LE_Data_GetPointer
    use Dims   , only : nx,ny,nz
    use Dims   , only : substantial_snowdepth
    
    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_Data, only : ilu_water_sea
    
    use indices
    
    !  --- in/out ---

    real, intent(inout)   ::  emis_a(nx,ny,nz,nspec)
    integer, intent(out)  ::  status
    
    ! --- const -----

    real, parameter :: fraction_c=.9, fraction_f=.1, conv_ug=1e6 
    real, parameter :: conv_min=365.*24.*60.
    
    ! --- local -----

    integer   ::  ix,iy
    real      ::  c_snow
    integer   ::  ilu

    real, pointer          ::  area(:,:,:)   ! (lon,lat,1)
    real, pointer          ::  sd(:,:,:)   ! (lon,lat,1)
    
    ! --- begin -----
    
    ! access meteo:
    call LE_Data_GetPointer( 'area', area, status, check_units ='m2' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'sd', sd, status, check_units ='1' )
    IF_NOTOK_RETURN(status=1)
   
    ! determine landuse type sea (less natural hm-emissions over sea)
    ilu = ilu_water_sea

    do ix=1,nx
      do iy=1,ny
         ! if snowcover no natural emissions
         c_snow=1.

         if ( sd(ix,iy,1) > substantial_snowdepth ) c_snow = 0.0 

         ! in g/km2/y * fraction fine of coarse/conv naar minuut*conv km2 naar gridcel opp/1e6 naar ug  
         emis_a(ix,iy,1,i_pb_f) = emis_a(ix,iy,1,i_pb_f) + fraction_f*area(ix,iy,1)*conv_ug/conv_min*&
                               (160*lu_fracs(ix,iy,ilu)+220*c_snow*(1-lu_fracs(ix,iy,ilu)) )
         emis_a(ix,iy,1,i_pb_c) = emis_a(ix,iy,1,i_pb_c) + fraction_c*area(ix,iy,1)*conv_ug/conv_min*&
                               (160*lu_fracs(ix,iy,ilu)+220*c_snow*(1-lu_fracs(ix,iy,ilu)) )
         emis_a(ix,iy,1,i_cd_f) = emis_a(ix,iy,1,i_cd_f) + fraction_f*area(ix,iy,1)*conv_ug/conv_min*&
                               (8  *lu_fracs(ix,iy,ilu)+12 *c_snow*(1-lu_fracs(ix,iy,ilu)) )
         emis_a(ix,iy,1,i_cd_c) = emis_a(ix,iy,1,i_cd_c) + fraction_c*area(ix,iy,1)*conv_ug/conv_min*&
                               (8  *lu_fracs(ix,iy,ilu)+12 *c_snow*(1-lu_fracs(ix,iy,ilu)) )
      enddo
    enddo
  
    ! ok
    status = 0

  end subroutine

end module LE_Emis_HM_Natural
