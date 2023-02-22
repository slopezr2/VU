!###############################################################################
!
! LE_Bound - LOTOS-EUROS boundary condition routines
!            Storage to save top layers of parent model, 
!            useful for simulation of satellite retrievals.
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

module LE_Bound_Top

  use GO, only : gol, goErr, goPr

  implicit none

  ! --- in/out --------------------------------

  private

  public  ::  LE_Bound_Top_MergeProfiles
  public  ::  LE_Bound_Top_GetProfile


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_Top'


contains

  ! ***


  !
  ! Merge LE profile with top profile.
  ! The original layer that includes the top of the LE layers 
  ! is split into two parts, the upper half and the higher layers
  ! are stacked on the LE arrays.
  ! Concentrations should be vetrical column densities (kg/m2),
  ! since this is the common unit used for interpolation of boundary conditions
  ! and simulations of satellite profiles.
  !

  subroutine LE_Bound_Top_MergeProfiles( ix, iy, ispec,&
                                         top_ph_entry, top_trc_entry, &
                                         le_phlev, le_vcd, le_units, &
                                         phlev, vcd, &
                                         status & 
#ifdef with_labeling
                                         ,label_aloft &
#endif                                         
                                        )
                                        
                                         
    use Num  , only : Interval        
    use Binas, only : grav, xm_air
    use Indices, only : specmolm
    use LE_Grid, only : lli

    ! --- in/out ------------------------------

    integer, intent(in)                   ::  ix, iy
    integer, intent(in)                   ::  ispec
    character(len=*), intent(in)          ::  top_ph_entry
    character(len=*), intent(in)          ::  top_trc_entry
    real, intent(in)                      ::  le_phlev(:)   ! (1:nz+1) [Pa]
    real, intent(in)                      ::  le_vcd(:)     ! (nz) [units]
    character(len=*), intent(in)          ::  le_units
    real, pointer                         ::  phlev(:)   ! (0:nlev) [Pa]
    real, pointer                         ::  vcd(:)     ! (nlev) [units]
    integer, intent(out)                  ::  status
#ifdef with_labeling
    logical, intent(in),optional          ::  label_aloft  ! needed for check aloft label
#endif

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_Top_MergeProfiles'

    ! --- local -------------------------------
    
    integer             ::  le_nlev
    real                ::  le_ptop
    real, pointer       ::  bt_phlev(:) ! 0:bt_nlev
    real, pointer       ::  bt_ahlev(:) ! 0:bt_nlev
    real, pointer       ::  bt_vcd(:)   ! 1:bt_nlev
    character(len=32)   ::  bt_units
    integer             ::  bt_nlev
    integer             ::  ibt
    integer             ::  nlev
    integer             ::  k
    character(len=32)   ::  conversion
    integer             ::  ilev
    
    ! --- begin -------------------------------
    
    ! init results:
    nullify( bt_phlev ) 
    nullify( bt_ahlev ) 
    nullify( bt_vcd )
    
    ! get profile for requested location:
    call LE_Bound_Top_GetProfile( ix, iy, &
                                  top_ph_entry, top_trc_entry, &
                                  bt_phlev, bt_vcd, bt_units, &
                                  status )
    IF_NOTOK_RETURN(status=1)
    
    
    ! level
    bt_nlev = size(bt_vcd)
    
    ! unit conversion ..
    conversion = trim(bt_units)//' ; '//trim(le_units)
    
    select case ( trim(conversion) )
        
      case ( 'kg/m2 ; kg/m2' )
        ! no conversion needed
        bt_vcd = bt_vcd
      
      case ( 'mole/mole ; kg/m2')
        
        do ilev = 1, bt_nlev
          ! kg tr /m2    mole tr/mole air   kg air/mole air          kg air m-2  m s-2               m s-2     kg tr / mole tr
          bt_vcd(ilev) = bt_vcd(ilev)    /  xm_air           * ( bt_phlev(ilev-1) - bt_phlev(ilev) )/grav  * specmolm(ispec)  
        end do

      case ( 'kg/kg ; kg/m2')
        do ilev = 1, bt_nlev 
          ! kg tr /m2    kg tr /kg air           kg air m-2 m s-2               m s-2
          bt_vcd(ilev) = bt_vcd(ilev)   * ( bt_phlev(ilev-1) - bt_phlev(ilev) )/grav 
        end do
      case default 
        
        write( gol, '("Define unit conversion for ",a, " to ", a)' ) trim(bt_units), trim(le_units) ; call goErr
        TRACEBACK;status=1;return
    end select
    
                
    ! dimensions:
    le_nlev = size(le_vcd)
    
    ! check ... 
    if ( size(le_phlev) /= le_nlev+1 ) then
      write (gol,'("LE concentration array has ",i4," layers while phlev has ",i4)') le_nlev, size(le_phlev); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! search layer with top of LE pressures;
    ! negate pressures to ensure increasing ax:
    ibt = 1  ! init search
    le_ptop = le_phlev(le_nlev+1)
    call Interval( -1.0*bt_phlev, -1.0*le_ptop, ibt, status )
    if ( status /= 0 ) then
      write (gol,'("could not find interval with value ",f12.2," in:")') le_ptop; call goErr
      do k = 0, bt_nlev
        write (gol,'("  ",i6,f12.2)') k, bt_phlev(k); call goErr
      end do
      TRACEBACK; status=1; return
    end if
    
    ! combined number of layers:
    nlev = le_nlev + (bt_nlev-ibt+1)
    
    ! clear current storage ?
    if ( associated(phlev) .and. (size(phlev) /= nlev+1) ) deallocate( phlev )
    if ( associated(vcd  ) .and. (size(vcd  ) /= nlev  ) ) deallocate( vcd   )
    ! storage:
    if ( .not. associated(phlev) ) allocate( phlev(0:nlev) )
    if ( .not. associated(vcd  ) ) allocate( vcd  (1:nlev) )
    
    ! copy pressures:
    phlev(0:le_nlev) = le_phlev
    phlev(le_nlev+1:nlev) = bt_phlev(ibt:bt_nlev)
    
    ! copy LE concentrations:
    vcd(1:le_nlev) = le_vcd
    
#ifdef with_labeling
    ! fill top boundary if necessary
    if ( present(label_aloft) ) then
      if ( label_aloft ) then
#endif    
        ! fill fraction of interface layer:
        vcd(le_nlev+1) =  bt_vcd(ibt) * (bt_phlev(ibt)-le_ptop)/(bt_phlev(ibt)-bt_phlev(ibt-1))
        ! copy top concentrations:
        vcd(le_nlev+2:nlev) = bt_vcd(ibt+1:bt_nlev) 
#ifdef with_labeling
      else
        ! Fill top concentrations with zero for not label aloft
        vcd(le_nlev+1) = 0.0
        vcd(le_nlev+2:nlev) = 0.0
      end if
    else ! no label aloft (probably normal OMI-output wanted in a label run)
      ! fill fraction of interface layer:
      vcd(le_nlev+1) =  bt_vcd(ibt) * (bt_phlev(ibt)-le_ptop)/(bt_phlev(ibt)-bt_phlev(ibt-1))
      ! copy top concentrations:
      vcd(le_nlev+2:nlev) = bt_vcd(ibt+1:bt_nlev) 
    endif
#endif    

    ! clear:
    deallocate( bt_phlev )
    deallocate( bt_vcd )

    ! ok:
    status = 0

  end subroutine LE_Bound_Top_MergeProfiles

  ! ***

  subroutine LE_Bound_Top_GetProfile( ix, iy, &
                                      input_ph_entry, input_trc_entry, &
                                      phlev, conc, units, &
                                      status )

    use LE_Data      , only : LE_Data_GetPointer

    ! --- in/out ------------------------------

    integer, intent(in)                   ::  ix, iy
    character(len=*), intent(in)          ::  input_ph_entry
    character(len=*), intent(in)          ::  input_trc_entry
    real, pointer                         ::  phlev(:)   ! (0:nlev) [Pa]
    real, pointer                         ::  conc(:)    ! (nlev) [units]
    character(len=*), intent(out)         ::  units
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_Top_GetProfile'

    ! --- local -------------------------------
    
    integer         ::  nlev
        
    real, pointer   ::  ph(:,:,:)
    real, pointer   ::  ctrc(:,:,:)
    
    ! --- begin -------------------------------
        
    ! get variable pointer
    call LE_Data_GetPointer( input_ph_entry, ph, status, check_units='Pa')
    IF_NOTOK_RETURN(status=1)
        
    call LE_Data_GetPointer( input_trc_entry, ctrc, status, units=units )
    IF_NOTOK_RETURN(status=1)
    
    ! levels in boundary file
    nlev = size(ctrc,3)
    if ( associated(phlev) .and. (size(phlev) /= nlev+1) ) deallocate( phlev )
    if ( associated(conc ) .and. (size(conc ) /= nlev) ) deallocate( conc  )
    ! new storage ?
    if ( .not. associated(phlev) ) allocate( phlev(0:nlev) )
    if ( .not. associated(conc ) ) allocate( conc (1:nlev) )

    ! put variables
    phlev(:) = ph(ix,iy,:)
    conc(:)  = ctrc(ix,iy,:)
    
    ! ok:
    status = 0

  end subroutine LE_Bound_Top_GetProfile


end module LE_Bound_Top

