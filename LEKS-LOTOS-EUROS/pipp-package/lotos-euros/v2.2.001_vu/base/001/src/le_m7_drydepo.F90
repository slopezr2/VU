!###############################################################################
!
! calculates deposition velocities for aerosol mass and number modes
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

module LE_M7_DryDepo

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private

  public  ::  LE_M7_DryDepo_Init, LE_M7_DryDepo_Done
  public  ::  LE_M7_DryDepo_vd
  public  ::  LE_M7_DryDepo_Setup

  public  ::  m7_vd_table

  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_M7_DryDepo'


  ! lookup table for deposition velocities for fixed set of radii (in um)
  integer, parameter :: nvdep=23
  real, parameter :: vdepradius(nvdep) = &
                   (/  0.001,  0.01,   0.05,  0.1,  0.3,  &
                       0.5  ,  0.7 ,   0.8 ,  0.9,  1.0,  &
                       1.2  ,  1.5 ,   2.0 ,  3.0,  4.0,  &
                       5.0  ,  6.0 ,   8.0 , 10.0, 15.0,  &
                      20.0  , 50.0 , 100.0                   /)

  real, parameter :: nmin=0.01 !minimum number of particles/cm for deposition

  ! should be function of grid cell ...
  real, parameter :: psurf=100000 !Pa
  real, parameter :: rhoref=1800. !kg/m3

  ! --- var ------------------------------------

  ! lookup table:
  real, allocatable   ::  m7_vd_table(:,:,:,:)   ! (nx,ny,nvdep,nlu)


contains


  ! ====================================================================
  ! ===
  ! === module init/done
  ! ===
  ! ====================================================================


  subroutine LE_M7_DryDepo_Init( status )

    use dims            , only : nx, ny, nz
    use indices         , only : n_m7
    use LE_Landuse_Data , only : nlu
    use LE_Data         , only : LE_Data_Enable

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_Init'

    ! --- begin ----------------------------------

    ! enabled ?
    if ( n_m7 > 0 ) then
      ! lookup table:
      allocate( m7_vd_table(nx,ny,nvdep,nlu) )
      ! enable meteo:
      call LE_Data_Enable( 'tsurf', status )
      IF_NOTOK_RETURN(status=1)
      call LE_Data_Enable( 't', status )
      IF_NOTOK_RETURN(status=1)
      call LE_Data_Enable( 'p', status )
      IF_NOTOK_RETURN(status=1)
      call LE_Data_Enable( 'dens', status )
      IF_NOTOK_RETURN(status=1)
      call LE_Data_Enable( 'lai_lu', status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine LE_M7_DryDepo_Init


  ! ***


  subroutine LE_M7_DryDepo_Done( status )

    use indices, only : n_m7

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_DryDepo_Done'

    ! --- begin ----------------------------------

    ! enabled ?
    if ( n_m7 > 0 ) then

      ! clear:
      deallocate( m7_vd_table )

    end if

    ! ok
    status = 0

  end subroutine LE_M7_DryDepo_Done


  ! ====================================================================
  ! ===
  ! === depos velocity
  ! ===
  ! ====================================================================


  ! ***
  
  
  ! table with surface deposition velocities for reference density.
  
  subroutine LE_M7_DryDepo_Setup( status )
  
    use Binas                 , only : grav, pi
    use JAQL_DryDeposition    , only : Slip_Correction_Factor, Sedimentation_Velocity
    use dims                  , only : nx, ny
    use LE_Data               , only : LE_Data_GetPointer

    use indices               , only : n_m7
    use Indices               , only : rhopart_default ! 1.14e3 kg/m3, to be consequent with le_sedim
    use LE_meteo_data         , only:  viscos, freepathlen
    use LE_Landuse_Data       , only : lu_fracs, nlu
    use LE_Landuse_Data       , only : ustar_lu, Ra_h_htop_z0_lu, gamma_lu, A_lu, alfa_lu
    use LE_DryDepos_Aero_Zhang, only : LE_DryDepos_Aero_Zhang_vd

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_DryDepo_Setup'

    

    ! --- local ----------------------------------

    integer ::  nwet
    integer ::  luc, nluc, i,j,k,n, mode
    real    ::  vd_table, vs_table, partsize_table, slipcor_table
    real    ::  Rs ! sedimentation resistance (not needed here)
    ! meteo data:
    real, pointer        ::  tsurf (:,:,:)   ! (lon,lat,1)
    real, pointer        ::  temp  (:,:,:)   ! (lon,lat,1)    
    real, pointer        ::  pres  (:,:,:)   ! (lon,lat,alt)        
    real, pointer        ::  dens  (:,:,:)   ! (lon,lat,alt)            
    real, pointer        ::  lai_lu(:,:,:)   ! (lon,lat,nlu)          

    ! --- begin ----------------------------------
    
    ! point to meteo data:
    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units='K' )          
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 't', temp, status, check_units='K' )            
    IF_NOTOK_RETURN(status=1)  
    call LE_Data_GetPointer( 'p', pres, status, check_units='Pa' )            
    IF_NOTOK_RETURN(status=1)      
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')            
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lai_lu', lai_lu, status, check_units='m2/m2' )
    IF_NOTOK_RETURN(status=1)

    ! only necessary if M7 is enabled ...
    if ( n_m7 > 0 ) then
    
      ! init table to zero:
      m7_vd_table = 0.0

      ! loop over aerosol radius' in lookup table:
      do k = 1, nvdep
        
        ! calculate vs, slip correction for sizes in table
        do i=1,nx
          do j=1,ny

            ! fom radius in um to diameter in m for vd_particles
            partsize_table = 2.0e-6 * vdepradius(k)
            slipcor_table = Slip_Correction_Factor( freepathlen(i,j,1), partsize_table )

            ! calculate vs, slip correction for sizes in table
            vs_table =  Sedimentation_Velocity( rhopart_default, partsize_table, &
                                                            slipcor_table, viscos(i,j,1))
            
            ! loop over landuse classes:
            do luc=1,nlu
              ! skip ?
              if ( lu_fracs(i,j,luc) <= 0.0 ) cycle
              ! for this landuse class only ;
              ! sedimentation velocity in lowest layer is input:
              call LE_DryDepos_Aero_Zhang_vd( vd_table, Rs, &
                              vs_table, partsize_table, slipcor_table, &
                              nwet, tsurf(i,j,1), dens(i,j,1), viscos(i,j,1), &
                              luc, lai_lu(i,j,luc), Ra_h_htop_z0_lu(i,j,luc), ustar_lu(i,j,luc), &
                              status )
              IF_NOTOK_RETURN(status=1)
              ! store: lookup table per landuse class:
              m7_vd_table(i,j,k,luc) = vd_table
            end do !luc

          end do !ny
        end do !nx

      end do !nvdep
      
    end if  ! with M7
    
    ! ok
    status = 0

  end subroutine LE_M7_DryDepo_Setup


  ! ***
  
  
  ! Deposition velocity of particles, in m/s
  ! Sedimentation is accounted for elsewhere
  ! Method from TM5 (Maarten Krol, personal communication), 
  ! adapted for LOTOS-EUROS by Astrid Manders

  subroutine LE_M7_DryDepo_vd( vdmass, vdnumber, &
                                mode, rwetm7, densm7, ntot, vdtab, &
                                status )

    use Binas     , only : pi
    !use Dims      , only : nx, ny, nz, nspec
    !use indices   , only : n_m7!, i_nnus
    use mo_aero_m7, only : sigmaln
    !use mo_aero_m7, only : nmod
    !use LE_M7_Data, only : rwetm7modes, rdrym7modes, densm7modes  

    ! --- in/out ---------------------------------

    real, intent(out)         ::  vdmass        ! m/s
    real, intent(out)         ::  vdnumber      ! m/s
    integer, intent(in)       ::  mode          ! 1..7
    real, intent(in)          ::  rwetm7        ! cm
    real, intent(in)          ::  densm7        ! g/cm3
    real, intent(in)          ::  ntot          ! total number density
    real, intent(in)          ::  vdtab(nvdep)  ! row from lookup table
    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_M7_DryDepo_vd'
    
    ! --- local ----------------------------------
    
    real  ::  vdshift(nvdep),vdeprad2(nvdep), d_aer(nvdep), logdp(nvdep), logdp2(nvdep)!nog allocatable maken?!
    real  :: dlogdp(nvdep+1), ddp(nvdep+1), ddpi(nvdep)
    real  :: nnumb(nvdep), nvolume(nvdep)
    real  :: lnsigma,density,vt, diampg, radius

    integer  :: i,j,ir, ir1,ir2, n,nrd, nshift
    
    ! --- begin ----------------------------------
    
    ! Calculate lookup table with precalculated deposition velocities for reference 
    ! density but actual meteo parameters and seasonal dependencies.
    ! For different densities the curve may shift, this will be taken into account
    ! further on in the code

    !All radii and diameters in um

    ! calculate the binsizes (um) around the radii of the pre-calculated vd's
    logdp(:)    = log10(2*vdepradius(:)) ! log(diameter)

    ! mean radius of the mode in um, density in kg/m3
    radius  = 1.0e4  * rwetm7 !cm->um
    density = 1000.0 * densm7 !g/cm3->kg/m3

    ! initialize array with deposition velocities for increasing radia:      
    vdshift(:) = vdtab
    lnsigma = sigmaln(mode)

    !Calculate deposition velocity for nonzero radius and density, nonzero number
    ! concentrations
    if ( (radius > 1e-4) .and. (density > 1e-2) .and. (ntot > nmin) ) then
       vdeprad2(:) = vdepradius(:)
       logdp2(:) = logdp(:)

       ! account for density different than density_ref of the look-up table (lur --> vdi):
       do ir = 2, nvdep
          if(vdshift(ir) > vdshift(ir-1)) exit  ! for bigger r's :  impaction dominates (density effects)
          if ( ir == nvdep ) exit   ! trap upper boundary
       end do      
       do ir1 = ir, nvdep
          vdeprad2(ir1) = vdepradius(ir1)*sqrt(rhoref/density)
          logdp2(ir1) = log10(2*vdeprad2(ir1))
       end do   
       ! compress look-up table such that radii are increasing monotonically:
       nshift = 0
       ir1 = ir
       do
          if ( logdp2(ir1) > logdp2(ir-1) ) exit
          nshift = nshift + 1
          ir = ir -1
          if(ir == 1) exit
       end do

       nrd = nvdep - nshift
       if (nshift > 0) then 
          do ir1 = ir, nrd
             dlogdp(ir1) = dlogdp(ir1+nshift)
             vdeprad2(ir1) = vdeprad2(ir1+nshift)
             vdshift(ir1) = vdshift(ir1+nshift)
          end do   
       endif
       ! Do the integration of the shifted lookup table
       ! prepare integration intervals

       dlogdp(1) = -3.0 
       do ir=2,nrd 
          dlogdp(ir) = 0.5*(logdp2(ir-1)+logdp2(ir))  !take middle of the log scale
       end do
       dlogdp(nrd+1) = 3.0   ! 1000 um
       ddp(1:nrd+1)  = 10**dlogdp(1:nrd+1)
       ddpi(1:nrd) = ddp(2:nrd+1)-ddp(1:nrd)   ! integration intervals (um)
       d_aer(1:nrd) = 2.0*vdeprad2(1:nrd)

       ! perform convolution with log-normal distribution:
       diampg = 2*radius  ! median diameter of mode in um
       ! calculate the distribution (number and mass) over the deposition bins:
       if(ntot > nmin .and. radius > tiny(radius) ) then !integrate
          do n=1,nrd
             nnumb(n) = &
             ntot/(sqrt(2.*pi)*d_aer(n)*lnsigma)*exp(-(log(d_aer(n))-log(diampg))**2/(2*lnsigma**2))
             nvolume(n) = nnumb(n)*(pi/6.0)*d_aer(n)**3
          end do   
          vt = sum(nnumb(1:nrd)*ddpi(1:nrd)*vdshift(1:nrd))/sum(nnumb(1:nrd)*ddpi(1:nrd))
       else   
          vt = 0.0
       endif   
       ! deposition velocity for number
       vdnumber = vt

       ! for mass:
       if(ntot > nmin .and. radius > tiny(radius) ) then !integrate
          vt = sum(nvolume(1:nrd)*ddpi(1:nrd)*vdshift(1:nrd))/sum(nvolume(1:nrd)*ddpi(1:nrd))
       else
          vt = 0.0
       endif   
       vdmass = vt

    else
       vdmass   = 0.0
       vdnumber = 0.0
    endif

    ! ok
    status = 0

  end subroutine LE_M7_DryDepo_vd


end module LE_M7_DryDepo
