!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_M7_Emis

  implicit none
  
  private
  public  ::  emism7
  
contains


subroutine emism7(c, dt_min)

!Converts emissions from LOTOS-EUROS sectional conventions to M7 lognormal modes
!AM July 2009

!emissions of BC, SO4a from emission database
use dims, only: nx,ny,nz, nspec, emis_a, volume
use indices
use binas, only: pi
use mo_aero_m7, only: naermod, nmod, sigmaln, avo, wso4

implicit none
real    :: c(nx,ny,nz,nspec) ! for m7 aerosol concentration in ug/cm3
real    :: dt_min

!densities and radii of emitted species
real, parameter:: rad_so4ait=0.03e-6 !(m)
real, parameter:: rad_so4acc=0.075e-6 !(m) 
real, parameter:: rad_ssa=0.079e-6 !(m)
real, parameter:: rad_ssc=0.63e-6 !(m)
real, parameter:: rad_bcai=0.034e-6 !(m)
real, parameter:: rad_ocai=0.025e-6  !(m)
real, parameter:: rad_ssac=0.079e-6 !(m)
real, parameter:: rad_ssco=0.63e-6 !(m)
real, parameter:: rad_duai=0.079e-6 !(m)
real, parameter:: rad_duci=0.63e-6 !(m)
real, parameter:: rad_duas=0.079e-6 !(m)
real, parameter:: rad_ducs=0.63e-6 !(m)

real, parameter:: dens_so4=1841. !kg/m3
real, parameter:: dens_ss=2165. !kg/m3
real, parameter:: dens_bc=2000.  !kg/m3
real, parameter:: dens_oc=2000 !kg/m3
real, parameter:: dens_du=2650 !kg/m3
real, parameter:: fracsol_oc=0.65
!real, allocatable :: radavmass(:,:,:), mass(:,:,:), number(:,:,:)
real, allocatable :: mass(:,:,:), number(:,:,:)

!conversion from ug/m3 to ug/cm3 and ug/cm3 to molec/cm3
real, parameter:: convcm=1.e-6
real, parameter:: kgtoug=1.e9
real convso4
!simple conversion to relate oc to bc, very inaccurate, depends on region and season
!should be improved
real, parameter:: bc2oc=3.0

real radavmass
!allocate(radavmass(nx,ny,nz))
allocate(mass(nx,ny,nz))
allocate(number(nx,ny,nz))

print *, 'constructing emism7 emissions'

!SO4 half of mass in aitken, half in accumulation mode 
!conversion from ug/m3 to molec/cm3 for so4a
convso4=avo/(wso4*1.e6)


!SO4ait
  mass=0.5*convcm*dt_min*1.0e-9*emis_a(:,:,:,i_so4as)/volume(:,:,:) !ug/m3
  radavmass=rad_so4ait*exp(1.5*(sigmaln(2))**2) 
  number=mass*3./(4*pi*(radavmass**3)*dens_so4*kgtoug) !#/cm3
  mass=mass*convso4 !molec/cm3
  !mass=number of moleculs! special treatment SO4
  c(:,:,:,i_so4ks) = c(:,:,:,i_so4ks) + mass 
  c(:,:,:,i_nais)=  c(:,:,:,i_nais)+number  
!SO4acc
  mass=0.5*convcm*dt_min*1.0e-9*emis_a(:,:,:,i_so4as)/volume(:,:,:) !ug/m3
  radavmass=rad_so4acc*exp(1.5*(sigmaln(3))**2) 
  number=mass*3./(4*pi*(radavmass**3)*dens_so4*kgtoug)  !#/cm3  
  mass=mass*convso4 !molec/cm3 
  !mass=number of molecules!special treatment SO4  
  c(:,:,:,i_so4as) = c(:,:,:,i_so4as) + mass
  c(:,:,:,i_nacs) = c(:,:,:,i_nacs) + number
!Black Carbon: all in aitken insoluble
  mass=dt_min*1.0e-9*emis_a(:,:,:,i_bcki)/volume(:,:,:) !ug/m3
  radavmass=rad_bcai*exp(1.5*(sigmaln(5))**2) 
  number=convcm*mass*3./(4.*pi*(radavmass**3)*dens_bc*kgtoug )  !#/cm3 
  c(:,:,:,i_bcki) = c(:,:,:,i_bcki) + mass
  c(:,:,:,i_naii) = c(:,:,:,i_naii) + number
 
!Organic Carbon: 
!all in aitken mode, fraction soluble, fraction insoluble
!organic carbon not in emission database.
!For now: organic carbon crudely related to black carbon
!soluble
  mass=bc2oc*fracsol_oc*dt_min*1.0e-9*emis_a(:,:,:,i_bcki)/volume(:,:,:) !ug/m3
  radavmass=rad_ocai*exp(1.5*(sigmaln(2))**2) 
  number=convcm*mass*3./(4.*pi*(radavmass**3)*dens_oc *kgtoug) !#/cm3  
  c(:,:,:,i_ocks) = c(:,:,:,i_ocks) + mass
  c(:,:,:,i_nais) = c(:,:,:,i_nais) + number
!insoluble
  mass=bc2oc*(1.-fracsol_oc)*dt_min*1.0e-9*emis_a(:,:,:,i_bcki)/volume(:,:,:) !ug/cm3
  radavmass=rad_ocai*exp(1.5*(sigmaln(5))**2) 
  number=convcm*mass*3./(4*pi*radavmass**3*dens_oc*kgtoug )   
  c(:,:,:,i_ocki) = c(:,:,:,i_ocki) + mass
  c(:,:,:,i_naii) = c(:,:,:,i_naii) + number
     
!Sea salt
!fine fraction in accumulation, coarse fraction in coarse mode, soluble
   emis_a(:,:,1,i_ssas)=emis_a(:,:,1,i_na_f)
   emis_a(:,:,1,i_sscs)=emis_a(:,:,1,i_na_c)
!Fine
  mass=dt_min*1.0e-9*emis_a(:,:,:,i_ssas)/volume(:,:,:) !ug/m3
  radavmass=rad_ssac*exp(1.5*(sigmaln(3))**2) 
  number=convcm*mass*3./(4*pi*radavmass**3*dens_ss*kgtoug ) !#/cm3  
  c(:,:,:,i_ssas) = c(:,:,:,i_ssas) + mass
  c(:,:,:,i_nacs) = c(:,:,:,i_nacs) + number
!coarse
  mass=dt_min*1.0e-9*emis_a(:,:,:,i_sscs)/volume(:,:,:) !g/cm3
  radavmass=rad_ssco*exp(1.5*(sigmaln(4))**2) 
  number=convcm*mass*3./(4*pi*radavmass**3*dens_ss*kgtoug ) !#/cm3   
  c(:,:,:,i_sscs) = c(:,:,:,i_sscs) + mass
  c(:,:,:,i_ncos) = c(:,:,:,i_ncos) + number

!dust not yet implemented in this version
  c(:,:,:,i_ducs) = 0.
  c(:,:,:,i_duas) = 0.
  c(:,:,:,i_duci) = 0.
  c(:,:,:,i_duai) = 0.
  
print *, 'after emissions', c(18,24,1,i_so4ns:i_duci)
print *, 'number conc',c(18,24,1,i_nnus:i_ncoi)
!deallocate(radavmass)
deallocate(mass)
deallocate(number)

 end subroutine emism7


end module LE_M7_Emis

