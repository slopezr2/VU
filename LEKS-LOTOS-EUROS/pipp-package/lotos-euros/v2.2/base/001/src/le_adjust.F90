!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Adjust

  use GO, only : gol, goPr, goErr
  
  implicit none
  

  ! --- in/out -----------------------------------
  
  private

  public  ::  LE_Adjust_Apply


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Adjust'


contains


  ! ====================================================================


  ! Redistributes the concentrations according to the
  ! volumes in the array "volume".
  ! On inpot, the concentrations are still diveded over
  ! the layers according to "ovolume" 

  subroutine LE_Adjust_Apply( c, status)

    use Dims         , only : outF
    use LE_Logging   , only : ident2
    use Dims         , only : nx, ny, nz
    use LE_Meteo_Data, only : ovolume, volume
    use LE_Meteo_Data, only : oairmass, airmass
    use LE_Meteo_Data, only : ohpres
    use LE_Bound     , only : caloft  ! (nx,ny,nz+1:nz_top,nspec)
    use Indices      , only : nspec

    ! --- in/out ---------------------------------

    real, intent(inout)       ::  c(nx,ny,nz,nspec)
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter :: rname = mname//'/LE_Adjust_Apply'
    
    ! --- local ----------------------------------
    
    integer     ::  i,j
    integer     ::  status_par
    
    real, pointer   ::  hp
    
    ! --- begin ----------------------------------
    
    ! info:
    if ( .not. outF%suppress ) then
      write (gol,'(a,"adjusting layer depths ...")') ident2; call goPr
    end if

    ! init parallel status:
    status_par = 0

    ! loop over the grid cells
    !$OMP parallel &
    !$OMP   default( none ) &
    !$OMP   shared( nx,ny ) &
    !$OMP   shared( volume, ovolume ) &
    !$OMP   shared( airmass, oairmass ) &
    !$OMP   shared( ohpres ) &
    !$OMP   shared( c ) &
    !$OMP   shared( caloft ) &
    !$OMP   private( i,j ) &
    !$OMP   private( status ) &
    !$OMP   reduction( + : status_par )
    !$OMP   do SCHEDULE( DYNAMIC )
    do i = 1, nx
      do j = 1, ny
    
      ! apply on column and all tracers:
      call LE_Adjust_Apply_1d_all( i, j, ovolume(i,j,:), volume(i,j,:), &
                                     oairmass(i,j,:), airmass(i,j,:), &
                                     ohpres(i,j,:), &
                                     caloft(i,j,nz+1,:), c(i,j,:,:), status )
      ! check on non-zero status:
      status_par = status_par + abs(status)

      end do ! j
    end do ! i
    !$OMP   end do
    !$OMP end parallel
    
    ! check ...
    if ( status_par /= 0 ) then
      write (gol,'("received non-zero return status from omp loop: ",i0)') status_par; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine LE_Adjust_Apply
  
  ! *
  
  subroutine LE_Adjust_Apply_1d_all( iloc, jloc, &
                                      ovolume, volume, oairm, airm, ohpr, &
                                       caloft, c, status )

    use Binas     , only : xm_air
    use Dims      , only : nz
    use Indices   , only : nspec, specunit, specmolm, i_o3

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  iloc, jloc  ! needed for labeling
    real, intent(in)          ::  ovolume(nz)    ! km3 (zucht)
    real, intent(in)          ::   volume(nz)    ! km3 (zucht)
    real, intent(in)          ::  oairm(nz)      ! kg
    real, intent(in)          ::   airm(nz)      ! kg
    real, intent(in)          ::  ohpr(0:nz)     ! Pa
    real, intent(in)          ::  caloft(nspec)  ! ppb or ug/m3
    real, intent(inout)       ::  c(nz,nspec)    ! ppb or ug/m3
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter :: rname = mname//'/LE_Adjust_Apply_1d_all'
    
    ! --- local ----------------------------------
    
    real        ::  v_int(nz), v_int_o(0:nz+1)
    integer     ::  l
    integer     ::  ispec
    real        ::  c_conc(nz), ca_conc

    ! --- begin ----------------------------------
  
    ! determine cumulative volumes by vertical integration
    v_int_o(0) = 0.0
    v_int_o(1) = ovolume(1)
    v_int(1)   =  volume(1)
    do l=2,nz
       v_int(l)   =  volume(l) + v_int  (l-1)
       v_int_o(l) = ovolume(l) + v_int_o(l-1)
    end do
    ! create extra data value above top of model:
    v_int_o(nz+1) = max( v_int_o(nz), v_int(nz) ) + 10.0

    ! deterime the cumulative mass for each species seperatly
    do ispec = 1, nspec
    
      ! check units, convert to (mass) concentration (per volume)
      ! since adjust is based on volume change
      select case ( trim(specunit(ispec)) )
        ! mole (volume) mixing ratio?
        case ( 'ppb' )
          ! conversion to mass concentation in ug/m3:
          !         1.0e-9             ! ((mole c)/(mole air)) / ppb
          !         1.0e9              ! ug/kg
          c_conc = c(:,ispec)                        &  ! ppb (mole c)/(mole air)/ppb ug/kg
                    * specmolm(ispec)/xm_air         &  ! (kg c)/(mole c) (mole air)/(kg air)
                    * oairm / (ovolume*1e9)            ! (kg air) / m3
          ! use top layer properties for aloft:
          ca_conc = caloft(ispec)                    &  ! (mole c)/(mole air) ug/kg
                    * specmolm(ispec)/xm_air         &  ! (kg c)/(mole c) (mole air)/(kg air)
                    * oairm(nz) / (ovolume(nz)*1e9)    ! (kg air) / m3

        !~ mass concentration
        case ( 'ug/m3', 'grns/m3' )
          c_conc  = c(:,ispec)     ! ug/m3
          ca_conc = caloft(ispec)  ! ug/m3
        case ( 'mlc/cm3', '1/cm3' )
          c_conc  = c(:,ispec)     ! 
          ca_conc = caloft(ispec)  ! 
        !~
        case default
          write (gol,'("volume concentrations not implemented for units `",a,"`")') &
                            trim(specunit(ispec)); call goErr
          TRACEBACK; status=1; return
      end select
    
      ! apply for this spec:
      call LE_Adjust_Apply_1d( iloc, jloc, ispec, &
                                ovolume, volume, &
                                ohpr, &
                                v_int, v_int_o, &
                                ca_conc, c_conc, &
                                status )
      IF_NOTOK_RETURN(status=1)
    
      ! convert back using new volume and airmass,
      ! no need to reset the aloft:
      select case ( trim(specunit(ispec)) )
        ! mole (volume) mixing ratio?
        case ( 'ppb' )
          ! conversion from mass concentation in ug/m3:
          !         1.0e9              ! ppb / ((mole c)/(mole air))
          !         1.0e-9             ! kg/ug
          c(:,ispec) = c_conc                        &  ! ug/m3 kg/ug ppb/((mole air)/(mole c))
                    / specmolm(ispec)*xm_air         &  ! (mole c)/(kg c) (kg air)/(mole air)
                    / airm * (volume*1e9)              ! m3 / (kg air)
        !~ mass concentration
        case ( 'ug/m3', 'grns/m3' )
          c(:,ispec) = c_conc    ! ug/m3
        case ( 'mlc/cm3', '1/cm3' )
          c(:,ispec) = c_conc    ! ug/m3
        !~
        case default
          write (gol,'("volume concentrations not implemented for units `",a,"`")') &
                            trim(specunit(ispec)); call goErr
          TRACEBACK; status=1; return
      end select

    end do  ! spec
       
    ! ok
    status = 0
       
  end subroutine LE_Adjust_Apply_1d_all
  
  ! *
  
  !
  ! apply adjust on 1D column (single spec)
  ! in future, return changes to pass to labeling routines
  ! instead of passing (iloc,jloc) and ispec
  !
  
  subroutine LE_Adjust_Apply_1d( iloc, jloc, ispec, &
                                   ovolume, volume, &
                                   hp_o, &
                                   v_int, v_int_o, &
                                   caloft, c, status )

    use Dims      , only : nz
    use Indices   , only : nspec, i_o3

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  iloc, jloc, ispec  ! needed for labeling
    real, intent(in)          ::  ovolume(nz)   ! km3 ?
    real, intent(in)          ::   volume(nz)   ! km3 ?
    real, intent(in)          ::  hp_o(0:nz)
    real, intent(in)          ::  v_int(nz)
    real, intent(in)          ::  v_int_o(0:nz+1)
    real, intent(in)          ::  caloft   ! ug/m3 ?
    real, intent(inout)       ::  c(nz)    ! ug/m3 ?
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter :: rname = mname//'/LE_Adjust_Apply_1d'
    
    ! --- local ----------------------------------
    
    integer     ::  l
    integer     ::  ix
    real        ::  int_m(0:nz+1)
    real        ::  dvratio, dpvratio, mass
    real        ::  hp_o_col(0:nz+1)
    
    ! --- begin ----------------------------------
      
    ! integrated mass=0 at bottom
    int_m(0) = 0.0
    do l=1,nz
      int_m(l) = ovolume(l)*c(l) + int_m(l-1)
    enddo
    
    ! aloft:
    int_m(nz+1) = int_m(nz) + (v_int_o(nz+1)-v_int_o(nz))*caloft
    
    ! pressure including aloft 
    hp_o_col(0:nz) = hp_o
    hp_o_col(nz+1) = hp_o(nz) - (hp_o(nz-1)-hp_o(nz))


    ! data points have now been created for the integrated concentrations
    ! interpolared linearly for positive redistribution
    ! NB. data points are int_m(i), v_int_o(i), for i=0, nz+1

    ! construct new conc. for each layer
    ! first compute the masses
    mass = 0.0
    do l=1,nz
      ix = 1
      ! determine between which layers of the old distribution
      ! interpolation must be done
      do
        if (v_int(l) > v_int_o(ix)) then
          ix=ix+1
        else
          exit
       endif
      enddo
      ! interpolation takes place between the old levels ix-1 and ix

      !! by definition, ix > 0, only check if ix has not become too large
      !! (should be impossible, only possible due to implementation error)
      !if (ix > nz+1) then
      !   write (gol,'("ix ",i0," above nz+1 ",i0)') ix, nz+1; call goErr
      !   TRACEBACK; status=1; return
      !end if
        
       ! interpolate linearly ratio between volumes
       dvratio = (v_int(l)-v_int_o(ix-1))/(v_int_o(ix)-v_int_o(ix-1))
       
       ! interpolate linearly ratio between volume*pressure ( with ideal gas law this estimates the mass flux)
       dpvratio = (0.5*(hp_o_col(ix) - hp_o_col(ix-1))*dvratio + hp_o_col(ix-1)) * dvratio*(v_int_o(ix)-v_int_o(ix-1)) / &
                  (0.5*(hp_o_col(ix) + hp_o_col(ix-1))*(v_int_o(ix)-v_int_o(ix-1)) )
       
       ! assign integrated mass to the concentration array of this
       ! species and layer
       ! which is equal to the interpolated value, minus the already
       ! redistributed mass in the lower layers
       c(l) = int_m(ix-1) + (int_m(ix)-int_m(ix-1))*dpvratio - mass

       ! trap round-off errors:
       c(l) = max( 0.0, c(l) )

       ! update the cumulative mass assigned so far
       mass = mass + c(l)

       ! divide by volume in order to get concentration 
       c(l) = c(l) / volume(l)

     ! end loop over layer
     enddo


    ! ok
    status = 0
       
  end subroutine LE_Adjust_Apply_1d
  

end module LE_Adjust

