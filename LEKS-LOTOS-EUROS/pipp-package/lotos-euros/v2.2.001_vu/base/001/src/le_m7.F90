!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_M7

  use GO, only : gol, goPr, goErr

  use MO_Aero_M7   , only : sigmaln
  use LE_M7_Data   , only : rdrym7modes ,rwetm7modes
  use LE_M7_Data   , only : densm7modes
  use LE_M7_DryDepo, only : LE_M7_DryDepo_Setup, LE_M7_DryDepo_vd
  use LE_M7_DryDepo, only : m7_vd_table

  implicit none


  ! --- in/out -----------------------------------

  private

  public    ::  LE_M7_Init, LE_M7_Done
  public    ::  LE_M7_Apply
  public    ::  LE_M7_UpdateModes
  public    ::  LE_M7_DryDepo_Setup, LE_M7_DryDepo_vd

  public    ::  sigmaln
  public    ::  rdrym7modes
  public    ::  rwetm7modes
  public    ::  densm7modes


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_M7'



contains


  ! ====================================================================
  ! ===
  ! === module init/done
  ! ===
  ! ====================================================================


  subroutine LE_M7_Init( status )

    use indices      , only : n_m7
    use mo_aero_m7   , only : m7_initialize
    use LE_M7_Data   , only : LE_M7_Data_Init
    use LE_M7_DryDepo, only : LE_M7_DryDepo_Init
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
     
    ! enabled ?
    if ( n_m7 > 0 ) then
      ! init data:
      call LE_M7_Data_Init( status )
      IF_NOTOK_RETURN(status=1)
      ! init drydepo:
      call LE_M7_DryDepo_Init( status )
      IF_NOTOK_RETURN(status=1)
      ! initialize parameters for M7
      call m7_initialize()
      ! enable data:
      call LE_Data_Enable( 't', status )
      IF_NOTOK_RETURN(status=1)
      call LE_Data_Enable( 'rh', status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine LE_M7_Init


  ! ***


  subroutine LE_M7_Done( status )

    use indices      , only : n_m7
    use LE_M7_Data   , only : LE_M7_Data_Done
    use LE_M7_DryDepo, only : LE_M7_DryDepo_Done

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_Done'

    ! --- begin ----------------------------------
     
    ! enabled ?
    if ( n_m7 > 0 ) then
      ! done with drydepo:
      call LE_M7_DryDepo_Done( status )
      IF_NOTOK_RETURN(status=1)
      ! done with data:
      call LE_M7_Data_Done( status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine LE_M7_Done


  ! ====================================================================
  ! ===
  ! === step
  ! ===
  ! ====================================================================


  ! Call to m7 with corresponding conversions of input and output
  ! input: concentrations in
  ! mass concentrations: ug/m3 for BC, OC, sea salt, dust
  ! number concentrations: particles/cm3
  ! h2so4 gas and SO4 aerosol concentrations in molecules/cm3
  ! Astrid Manders 2009

  subroutine LE_M7_Apply( c, dtime, status )

    use binas, only: Rgas 
    use dims,  only : nx, ny, nz, nspec
    use LE_M7_Data, only : rdrym7modes, rwetm7modes, waterm7modes, densm7modes
    use mo_aero_m7, only: naermod, nmod, nsol
    use indices
    
    use LE_Data               , only : LE_Data_GetPointer
    
    ! --- in/out ----------------------------------

    real, intent(inout)       ::  c(nx,ny,nz,nspec)
    real, intent(in)          ::  dtime !in minutes, to be converted to seconds in call to m7
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_Apply'

    ! --- local ----------------------------------

    real, dimension(1,1):: presdum, tempdum,rhdum
    real,  dimension(1,1,naermod)::densm7, waterm7, radiusm7, dryradiusm7
    real,  dimension(1,1,naermod):: concentm
    real, dimension(1,1,nmod):: concentn
    real,  dimension(1,1)::concentso4g
    integer :: lev, sp1,sp2,j,i
    ! meteo data:
    real, pointer        ::  temp(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  rh  (:,:,:)   ! (lon,lat,nz)

    ! --- begin ----------------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 't', temp, status, check_units='K' )    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'rh', rh, status, check_units='%' )    
    IF_NOTOK_RETURN(status=1)
    
    ! enabled ?
    if ( n_m7 > 0 ) then

      presdum=101300. !Should be improved, variation with height
      radiusm7=0.0    !radius and density are determined in M7
      dryradiusm7=0.0
      densm7=0.0
      densm7modes=0.0
      waterm7modes=0.0
      rdrym7modes=0.0
      rwetm7modes=0.0
      !temp=288.
      !rh=80.

      !restructure arrays from LOTOS-EUROS (nx,ny,nz)  to M7 (nx*ny,nz) conventions
      !call m7, then restructure m7 output to LOTOS-EUROS again.
      do lev=1,nz

         do j=1,ny
            do i=1, nx

              do sp1=1, naermod !species mass
                 concentm(1,1,sp1)=c(i,j,lev,i_so4ns+sp1-1)
              end do
              do sp2=1, nmod !number modes
                 concentn(1,1,sp2)=c(i,j,lev,i_nnus+sp2-1)
              end do
              concentso4g(1,1)=  c(i,j,lev,i_h2so4)  
              tempdum(1,1)=  temp(i,j,lev)
              rhdum(1,1)=  rh(i,j,lev)/100.0


              call m7(1, 1,1,presdum(1,1), rhdum(1,1), tempdum(1,1),&
                 concentso4g(1,1), concentm(1,1,:), &
                 concentn(1,1,:), &
                 densm7, waterm7,radiusm7, dryradiusm7,&
                 dtime*60)


              do sp1=1, naermod
                 c(i,j,lev,i_so4ns+sp1-1)=concentm(1,1,sp1)
              end do
              do sp2=1, nmod
                  c(i,j,lev,i_nnus+sp2-1)=concentn(1,1,sp2)
                  densm7modes(i,j,lev,sp2)=densm7(1,1,sp2)
                  waterm7modes(i,j,lev,sp2)=waterm7(1,1,sp2)
                  rwetm7modes(i,j,lev,sp2)=radiusm7(1,1,sp2)
                  if (sp2.le.nsol) then
                      rdrym7modes(i,j,lev,sp2)=dryradiusm7(1,1,sp2)      
                  endif
              end do
              c(i ,j,lev,i_h2so4)  = concentso4g(1,1)
            end do !nx
         end do !ny
      end do !levels
      
    end if  ! m7 tracers enabled ?
    
    ! ok
    status = 0

  end subroutine LE_M7_Apply

  !Call to m7 routines which calculate count median diameters and densities of the mode
  !depending on mass and number concentrations.
  !No additional aerosol processes, only change in distribution (diameter, density), 
  !not in concentration.
  !July 2009 AM 

  subroutine LE_M7_UpdateModes( c, status, warning )

    use binas, only: Rgas 
    use dims,  only : nx, ny, nz, nspec
    use LE_M7_Data,  only: rdrym7modes, rwetm7modes, waterm7modes, densm7modes
    use mo_aero_m7, only: naermod, nmod, nsol, nss
    use indices
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out ---------------------------------

    real, intent(inout)       ::  c(nx,ny,nz,nspec)
    integer, intent(out)      ::  status

    logical, intent(in), optional  ::  warning

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_UpdateModes'

    ! --- local ----------------------------------

    real,  dimension(1,1):: presdum, tempdum,rhdum
    real,  dimension(1,1,naermod)::densm7, waterm7, radiusm7, dryradiusm7
    real,  dimension(1,1,naermod):: concentm
    real, dimension(1,1,nmod)::concentn
    real, dimension(1,1)::concentso4g
    real, dimension(1,1,naermod):: avmasscomp
    real, dimension(1,1,nss)::hplus
    integer ::  lev, sp1,sp2,i,j

    ! meteo data:
    real, pointer               ::  temp(:,:,:)   ! (lon,lat,nz)
    real, pointer               ::  rh  (:,:,:)   ! (lon,lat,nz)

    ! --- begin ----------------------------------
    
    ! enabled ?
    if ( n_m7 > 0 ) then

      ! point to meteo data:
      call LE_Data_GetPointer( 't', temp, status, check_units='K' )    
      IF_NOTOK_RETURN(status=1)
      call LE_Data_GetPointer( 'rh', rh, status, check_units='%' )    
      IF_NOTOK_RETURN(status=1)
    
      !  avmasscomp      = average mass for single compound in each mode 
      !               [in molec. for sulphate and in ug for bc, oc, ss, and dust]
      !  hplus     = number of h+ in mole [???] (kg water)-1

      presdum=101300. !should be improved, height dependent pressure
      !print *,'begin update, densm7mode, rdry',densm7modes(8,10,2,3),rdrym7modes(8,10,2,3) 
      do lev=1,nz
         do j=1,ny
           do i=1,nx
              do sp1=1, naermod !species mass
               concentm(1,1,sp1)=c(i,j,lev,i_so4ns+sp1-1)
              end do
              do sp2=1, nmod !number modes
               concentn(1,1,sp2)=c(i,j,lev,i_nnus+sp2-1)
              end do
              concentso4g(1,1)=  c(i,j,lev,i_h2so4)  
              tempdum(1,1)=  temp(i,j,lev)
              rhdum(1,1)=  rh(i,j,lev)/100.0



         avmasscomp(1,1,:)=0.
        !--- 1) Calculation of particle properties under ambient conditions: -----
        !
        !--- 1.1) Calculate mean particle mass for all modes 
        !         and dry radius and density for the insoluble modes.
         call m7_averageproperties(1,1,1,&
                 concentn(1,1,:), concentm(1,1,:),avmasscomp(1,1,:),radiusm7,densm7)

        !--- 1.2) Calculate ambient count median radii and density 
        !         for lognormal distribution of particles.
        !
        !         Sulfate particles:
        !
        CALL m7_equiz(1,1,1,presdum(1,1),   &
                      avmasscomp(1,1,:),  tempdum(1,1),  &
                      rhdum(:,1), radiusm7, dryradiusm7, &
                      densm7,  waterm7,  concentn(1,1,:)   )

        !         
        !         Mixed particles with sulfate, b/o carbon and dust: 
        !
        CALL m7_equimix(1, 1,1,presdum(1,1),   &
                      avmasscomp(:,1,:),  tempdum(1,1),  &
                      rhdum(1,1), radiusm7, dryradiusm7, &
                      densm7,  waterm7,  concentn(1,1,:)   )

        !         Accumulation and coarse mode particles in presence of
        !         sea salt particles:
        !
        CALL m7_equil(1, 1,1, rhdum(1,1),concentm(1,1,:), concentn(1,1,:),  &
                      radiusm7, dryradiusm7, hplus, waterm7,    densm7,  tempdum(1,1)    )
         !

            do sp1=1, naermod
              c(i,j,lev,i_so4ns+sp1-1)=concentm(1,1,sp1)
            end do
            do sp2=1, nmod
              c(i,j,lev,i_nnus+sp2-1)=concentn(1,1,sp2)
              densm7modes(i,j,lev,sp2)=densm7(1,1,sp2)
              waterm7modes(i,j,lev,sp2)=waterm7(1,1,sp2)
              rwetm7modes(i,j,lev,sp2)=radiusm7(1,1,sp2)
              if (sp2.le.nsol) then
                 rdrym7modes(i,j,lev,sp2)=dryradiusm7(1,1,sp2)      
              endif
            end do
            c(i,j,lev,i_h2so4)  = concentso4g(1,1)
         end do
       end do
      end do 
      
    end if  ! m7 tracers enabled ?

    ! ok
    status = 0

  end subroutine LE_M7_UpdateModes


end module LE_M7
