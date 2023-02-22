!###############################################################################
!
! Wet deposition.
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

module LE_WetDepos_Emep

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private
  
  public  ::  LE_WetDepos_Emep_Init, LE_WetDepos_Emep_Done
  public  ::  LE_WetDepos_Emep_Apply


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'LE_WetDepos_Emep'



contains


  ! ====================================================================


  subroutine LE_WetDepos_Emep_Init( status )
  
    use indices
    use dims, only : wash
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    integer, intent(out)                  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Emep_Init'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------

    ! set washout ratios
    wash        = 0.0
    if (i_hno3>0) wash(i_hno3) = 5.0e5
    if (i_so2 >0) wash(i_so2 ) = 1.5e5
    if (i_nh3 >0) wash(i_nh3 ) = 5.0e5
    if (i_h2o2>0) wash(i_h2o2) = 5.0e5
    if (i_form>0) wash(i_form) = 0.5e5

    ! enable data:
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_WetDepos_Emep_Init


  ! ***


  subroutine LE_WetDepos_Emep_Done( status )

    ! --- in/out ------------------------------

    integer, intent(out)                  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Emep_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------

    ! ok
    status = 0
    
  end subroutine LE_WetDepos_Emep_Done


  ! ***


  !***********************************************************************
  !                                 copyright by
  !       National Institute of Public Health and Environmental Protection
  !                          Laboratory for Air Research
  !                                  (RIVM/LLO)
  !                               The Netherlands
  ! No part of this software may be used, copied or distributed without permission
  ! of RIVM/LLO (1994)
  !
  !  NAME                 : %M%
  !  SCCS(SOURCE)         : %P%
  !  RELEASE - LEVEL      : %R% - %L%
  !  BRANCH -SEQUENCE     : %B% - %S%
  !  DATE - TIME          : %E% - %U%
  !  MAKE                 : mk%M% (without extension)
  !  WHAT                 : %W%:%E%
  !  AUTHOR               : Liesbeth de Waal/Anton van Weel
  !  FIRM/INSTITUTE       : RIVM-LLO-APS
  !  LANGUAGE             : FORTRAN(HP-UX, HP-F77)
  !  USAGE                : %M%
  !  DESCRIPTION  : calculates new concentration after wet deposition.
  !  PRECONDITION(S)      :
  !  POSTCONDITION(S)     :
  !  EXIT CODES           :
  !  REFERENCE            :
  !  FILES AND OTHER
  !       I/O DEVICES     :
  !  SYSTEM DEPENDENCIES  :
  !  UPDATE HISTORY       :
  !  2004-7-20 Ferd Sauter: test versie voor LOTOS-EUROS
  !  2004-9-21 Martijn Schaap: Combined routine LOTOS-EUROS
  !                            new aerosol parameterisation (LOTOS)
  !                            New washout ratio's, plus allocatable
  ! 2009-8 Astrid Manders : Washout for modal aerosol distribution (M7), 
  !                         collection efficiencies according to
  !                         UNI-AERO (EMEP)
  !*********************************************************************

  subroutine LE_WetDepos_Emep_Apply( c, cwet, update_cwet, deltat, status )

    use dims        , only : nx, ny, nz, nspec
    !use dims        , only : h
    use dims        , only : outF, runF
    use dims        , only : wash
    use indices     , only : specmode
    use indices     , only : NO_AEROSOL_MODE, AEROSOL_FINE_MODES, AEROSOL_COARSE_MODE
    use indices     , only : AEROSOL_FF_MODES, AEROSOL_CCC_MODE, AEROSOL_CC_MODE
#ifdef with_m7
    use indices     , only : AEROSOL_NUCL_MODE, AEROSOL_AITKEN_MODE, AEROSOL_ACCUM_MODE
#endif
    use indices     , only : AEROSOL_ULTRA_FINE_MODE, AEROSOL_ULTRA_FINE_FINE_MODE, AEROSOL_FINE_MEDIUM_MODE, AEROSOL_MEDIUM_COARSE_MODE
    use indices     , only : specname
    use LE_Logging  , only : ident2
    use LE_Data      , only : LE_Data_GetPointer

    ! --- in/out ------------------------------

    !  c         : concentrations for ncomp components and nz layers in ug/m3
    !  cwet      : budget term
    !  deltat       : timestep (min)
    real, intent(inout)   ::  c(nx,ny,nz,nspec)
    real, intent(inout)   ::  cwet(nx,ny,nz,nspec)
    logical, intent(in)   ::  update_cwet
    real, intent(in)      ::  deltat
    integer, intent(out)  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Emep_Apply'

    real, parameter   ::  A = 5.2
    real, parameter   ::  Vrd = 5.0
    real, parameter   ::  rho = 1e3

    ! --- local -------------------------------

    !  wetd       : change in concentration due to wet deposition (ug/m3)
    !  nldt         : number of times to run the loop / dt (1/h)
    !  ct           : concentration in air per layer
    !  wrt          : wash-out ratio times precipitation amount per loop
    !  r            : rain concentration per layer
    !  d            : change of concentration in air per layer
    !  ispec        : component routine is called for
    !  dh           : depth of the nz layers (m)
    !  wash         : wash-out ratio (=-99 if no wet deposition)
    !  preco        : precipation (m/h)
    !  atmos        : 0; impossible for comp. to return from rain to atmosphere
    !               : 1; possible for comp. to return from rain to atmosphere
    real        ::  preco, dt
    real        ::  ri, rate, dh(nz), rate_a, rate_g
    integer     ::  atmos
    integer     ::  ix,iy,ispec

    real        ::  scav_efficiency
    real        ::  r_sub
    integer     ::  l
    
    real        ::  delta_cwet(nx,ny,nz,nspec)
    real        ::  dcwet_lay

    ! meteo data:
    real, pointer       ::  temp(:,:,:)   ! (lon,lat,nz)
    real, pointer       ::  delta_h(:,:,:)   ! (lon,lat,nz)
    real, pointer       ::  rain(:,:,:)   ! (lon,lat,1)
    
    ! --- begin ----------------------------------

    call LE_Data_GetPointer( 't', temp, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'rain', rain, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dh', delta_h, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)

    ! info ...
    if (.NOT.outF%suppress) then
      write (gol,'(a,"<wet deposition>")') ident2; call goPr
    end if
    ! initialize
    delta_cwet = 0.0

    !set atmos such that evap is not possible
    dt = deltat*60 ! convert minutes to sec
    atmos = 0

    ! loop over grid cells:
    DO ix = 1, nx
      DO iy = 1,ny

        ! for now: no wetdep if snow or ice rain.....
        ! some rain: > 0.3 mm/3hr  =  0.1e-3/3600 m/s
        if (rain(ix,iy,1) > 0.1e-3/3600.0 ) then

          ! rain intensity in mm/s
          !      m/s            mm/m
          ri = rain(ix,iy,1) * 1000.0 ! mm/s
          ! preco in m/s
          preco = rain(ix,iy,1) ! m/s
          ! set dh (the layer depth) for all layers in meters
          dh = delta_h(ix,iy,1:nz)

          ! Loop over all tracers:
          do ispec = 1, nspec

            ! aerosol or gas ?
            if ( specmode(ispec) == NO_AEROSOL_MODE ) then
              ! gas tracer; not washed out ?
              if ( wash(ispec) <= 0.0 ) cycle
              ! some resistance parameterization, needs rain intens in mm/s
              r_sub =  (wash(ispec)*ri/1000./rho)
            else
              ! aerosols; set scav. efficiency:
              select case ( specmode(ispec) )
#ifdef with_m7
                case ( AEROSOL_NUCL_MODE   )          ; scav_efficiency = 0.4
                case ( AEROSOL_AITKEN_MODE )          ; scav_efficiency = 0.3
                case ( AEROSOL_ACCUM_MODE  )          ; scav_efficiency = 0.1
#endif
                case ( AEROSOL_FINE_MODES  )          ; scav_efficiency = 0.1              
                case ( AEROSOL_FF_MODES  )            ; scav_efficiency = 0.1
                case ( AEROSOL_CCC_MODE )             ; scav_efficiency = 0.4
                case ( AEROSOL_CC_MODE )              ; scav_efficiency = 0.4
                case ( AEROSOL_COARSE_MODE )          ; scav_efficiency = 0.4
                case ( AEROSOL_ULTRA_FINE_MODE )      ; scav_efficiency = 0.1
                case ( AEROSOL_ULTRA_FINE_FINE_MODE ) ; scav_efficiency = 0.1
                case ( AEROSOL_FINE_MEDIUM_MODE )     ; scav_efficiency = 0.4
                case ( AEROSOL_MEDIUM_COARSE_MODE )   ; scav_efficiency = 0.4
                case default
                  write (gol,'("unexpected aerosol mode ",i6," for tracer ",i6," (",a,")")') &
                                 specmode(ispec), ispec, trim(specname(ispec)); call goErr
                  TRACEBACK; status=1; return
              end select
              ! resistance parameterization:
              r_sub = (A*ri/Vrd) * scav_efficiency
            end if

            ! loop over layers:
            do l = 1, nz
              ! removal rate, seems the same for each layer:
              rate = exp(-(r_sub)*dt)
              ! change:
              dcwet_lay = (1-rate)*c(ix,iy,l,ispec)*dh(l)
              ! acumulated change:
              delta_cwet(ix,iy,l,ispec) = delta_cwet(ix,iy,l,ispec) + dcwet_lay
              ! update budgets ?
              if ( update_cwet ) then
                ! update budget:
                cwet(ix,iy,l,ispec) = cwet(ix,iy,l,ispec) + dcwet_lay
              end if
     
             ! update concentration:
              c   (ix,iy,l,ispec) =                        rate *c(ix,iy,l,ispec)
            end do  ! layers

          end do  ! tracers

        end if   ! rain

      end do   ! loops over x and y
    end do

    ! ok
    status = 0

  end subroutine LE_WetDepos_Emep_Apply


end module LE_WetDepos_Emep
