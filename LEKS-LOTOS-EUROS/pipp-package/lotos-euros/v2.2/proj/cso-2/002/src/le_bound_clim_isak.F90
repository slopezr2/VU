!###############################################################################
!
! NAME
!   LE_Bound_Clim_Isak - LOTOS-EUROS boundary conditions from Isaksen model
!
! DESCRIPTION
!
!   Read concentrations from 2D Isaksen model (very ancient ...).
!
!    spec    Isak 
!    ------- ---- 
!    so4a         
!    nh4a         
!    no3a         
!    bc           
!    pm25     
!    pm10     
!    na_f         
!    na_c         
!    dust_f       
!    dust_c       
!    no2      i   
!    no       i   
!    o3       i   
!    eth      i   
!    ole      i   
!    par      i   
!    ald      i   
!    form     i   
!    xyl      i
!    tol      i
!    co       i   
!    ch4      i   
!    so2      i   
!    pan      i   
!    mgly     i
!    cres     i
!    hno2     i
!    hno3         
!    nh3          
!    h2o2         
!    open     
!    to2      
!    iso      
!    ispd     
!    terp     
!    no3      
!    oh           
!    ho2      
!    n2o5     
!    c2o3     
!    xo2      
!    xo2n     
!    cro      
!    -------------
!
! HISTORY
!   Former module 'boundary' in 'bound.F90' .
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') mname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Bound_Clim_Isak

  use GO, only : gol, goErr, goPr
  
  use dims, only : nx, ny, nz, nspec
  use indices

  implicit none
  
  
  ! --- const ------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_Clim_Isak'

  ! the dimension of the boundary conditions in the files
  ! the number of species
  integer, parameter :: bc_spec = 18

  ! the number of layers and zones
  integer, parameter :: bc_lay = 7, bc_zone = 5
  
  ! the heights of the Isaksen levels
  real, parameter :: bc_height(0:bc_lay) = &
                     (/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5/) ! km
  
  ! --- var --------------------------------

  ! the array with the isaksen boundary conditions
  real :: bc_isak(bc_spec, 0:bc_lay, bc_zone)

  ! array pointing to the actual species for each species in
  ! the boundary condition file
  integer, dimension(bc_spec) :: l_bc_spec

  ! pointers to the zones for each y-coordinate
  integer, allocatable, dimension(:) :: l_bc_zone


contains


  ! ==================================================
  

  subroutine LE_Bound_Clim_Isak_Init( status )

    use LE_Data, only : LE_Data_Enable
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_Clim_Isak_Init'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)

    ! pointers to the zones for each y-coordinate
    allocate( l_bc_zone(ny) )
    
    ! ok
    status = 0

  end subroutine LE_Bound_Clim_Isak_Init
  
  
  ! ***
  

  subroutine read_bc( im, path, status )

    ! subroutine reads the Isaksen boundary conditions

    use GO,   only : goGetFU
    use LE_Logging, only : ident1
    use indices
    use LE_Grid, only : ugg

    ! --- in/out --------------------
    
    integer, intent(in)             ::  im
    character(len=*), intent(in)    ::  path
    integer, intent(out)            ::  status
    
    ! --- const ---------------------

    character(len=*), parameter ::  rname = mname//'/read_bc'
    
    ! abbreviations for each month
    character(len=3), parameter :: ext(12) = &
           (/'jan','feb','mar','apr','may','jun', &
             'jul','aug','sep','oct','nov','dec' /)

    ! --- local ---------------------
    
    integer             ::  u_bound
    character(len=1024) ::  fname
    integer             ::  id1, id2, id3, id4, izone, ispec, iy
    real                ::  y
    
    ! --- begin ---------------------

    ! array pointing to the actual species for each species in
    ! the boundary condition file
    ! Note that i_so4a is replaced by 0 since the SO4a boundary condition is not set here
    ! but in bound_logan
    l_bc_spec = &
          (/i_no2, i_o3, i_mgly, i_eth, i_ole, i_par, i_ald, i_form, i_xyl, &
            i_tol, i_pan, i_co, i_cres, i_hno2, i_so2, 0, i_no, i_ch4 /)
    ! in cb99, ald is called ald2 ...

    ! init the array
    bc_isak = 0.0

    ! print message
    !print *, '<reading b.c. for the month ', im,'>'

    ! determine file name
    fname = trim(path)//'bc_isak_'//ext(im)//'.txt'

    !! print message
    !print *,ident1, '>> reading file ',trim(fname)
    
    ! new file unit:
    call goGetFU( u_bound, status )
    IF_NOTOK_STOP

    open (u_bound,file=fname,status='old')

    ! read the file
    do ispec=1,bc_spec
      do izone=1,bc_lay
        !!! NB: zones in the file are ordered from north the south, here
        !!! we reverse the order !!!!!!
        read (u_bound,*) id1,id2,id3,(bc_isak(id1,id2,id4),id4=bc_zone,1,-1)
    !   read (u_bound,18) id1,id2,id3,(bc_isak(id1,id2,id4),id4=bc_zone,1,-1)
    !   18  format (3i3,5f15.7)
      enddo
      ! set boundary conditions for SO4 to zero, which is
      ! much more realistic value than those from Isaksen
    ! if (l_bc_spec(ispec)==i_so4) bc_isak(ispec, : ,: ) = 0.0
      if (l_bc_spec(ispec)==i_so2) bc_isak(ispec, : ,: ) = 0.0
    enddo

    ! copy the first layer to layer zero
    bc_isak(:,0,:)=bc_isak(:,1,:)

    ! close the file
    close (u_bound)

    select case( trim(ugg%type) )
      case ( 'cartesian-regular', 'cartesian' ) 
        ! construct the pointer array to the zones
        ! zones are
        do iy = 1, ny
          y = ugg%latitude_bnds_1d(1,1) + (iy-0.5)*ugg%dlat
          if (y < 40) then
            l_bc_zone(iy) = 1
          else if (y < 50) then
            l_bc_zone(iy) = 2
          else if (y < 60) then
            l_bc_zone(iy) = 3
          else if (y < 70) then
            l_bc_zone(iy) = 4
          else if (y < 80) then
            l_bc_zone(iy) = 5
          else
            print *,'no Isaksen data for latitude = ', y
            stop
          endif
        enddo

        ! copy the first layer to layer zero
        bc_isak(:,0,:)=bc_isak(:,1,:)


        ! construct the pointer array to the zones
        ! zones are
        !do iy = 1, ny
        !   y = runF%southb + (iy-0.5)*runF%dlat
        !   if (y < 40) then
        !       l_bc_zone(iy) = 1
        !   else if (y < 50) then
        !       l_bc_zone(iy) = 2
        !   else if (y < 60) then
        !       l_bc_zone(iy) = 3
        !   else if (y < 70) then
        !       l_bc_zone(iy) = 4
        !   else if (y < 80) then
        !       l_bc_zone(iy) = 5
        !   else
        !       print *, 'error: outside domain of boundary conditions, lat = ', y
        !   endif
        !enddo
      case default
        write( gol, '("Conversion not yet made for boundary ISAK to universal grid")' ) ; call goErr
        write( gol, '("Make input file netCDF compliant to be used in le-bound-data format")' ) ; call goErr
        TRACEBACK;status=1;return

    end select
    
    ! ok
    status = 0

  end subroutine read_bc


  ! ***
  
  
  ! fills the boundary conditions
  !
  !   spec_filled(1:nspec)  : set to .true. if bc is filled, unchanged otherwise
  !

  subroutine makebc( spec_filled, status )

    use dims           , only : nspec
    use LE_Logging     , only : ident2
    use LE_Data        , only : LE_Data_GetPointer
    use LE_Data_Common , only : nlev
    use LE_Bound_Common, only : bc_west, bc_east, bc_north, bc_south, caloft

    ! --- in/out ----------------------------
    
    logical, intent(inout)      ::  spec_filled(nspec)
    integer, intent(out)        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/makebc'
    
    ! --- local -----------------------------
    
    integer             :: ix, iy
    integer             :: ind_lat, i_lat
    integer             :: ind_spec
    integer             :: k
    integer             :: ind
    integer             :: iz 
    integer             :: i_lon
    integer             :: i, j
    real                :: frac, hh
    real, pointer       ::  h_m(:,:,:)   ! (lon,lat,nz)
    
    ! --- begin -----------------------------

    !print message
    !print *,ident2,'<constructing boundary conditions>'

    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m' )    
    IF_NOTOK_RETURN(status=1)

    !
    ! ~~ Isaksen specs
    !

    ! the boundary conditions in longitue direction (i.e. north and south)
    do i=1,nx
    do j=1,nz

    ! i_lon =1 refers to south, i_lon=1 to north
    do i_lon = 1, 2
       ! compute level half-way model layer
       if ( i_lon == 1 ) then
         iy = 1
       else
         iy = ny
       end if

       if (j > 1) then
          hh = 0.5*(h_m(i,iy,j)+h_m(i,iy,j-1))*1e-3  ! km
       else
          hh = 0.5*h_m(i,iy,1)*1e-3  ! km
       endif

       ! determine the layer
       do iz = 1, bc_lay
         ind = iz
         if (hh < bc_height(iz) ) exit
       enddo
       ! determine interpolation factors
       frac = (bc_height(ind)-hh)/(bc_height(ind)-bc_height(ind-1))
       frac = min( 1.0, max(0.0, frac ))

       !loop over the species
       do k=1,bc_spec
          ind_spec = l_bc_spec(k)
          ind_lat  = l_bc_zone(iy)
          if ( ind_spec > 0 ) then
            if (i_lon == 1) then
              bc_south(i,j,ind_spec) =  frac*bc_isak(k,ind-1,ind_lat) + &
                                  (1.0-frac)*bc_isak(k,ind  ,ind_lat) 
            else
              bc_north(i,j,ind_spec) =  frac*bc_isak(k,ind-1,ind_lat) + &
                                  (1.0-frac)*bc_isak(k,ind  ,ind_lat) 
            endif
          endif
       enddo

    ! end north-south loop
    enddo

    ! end loop over cells
    enddo
    enddo

    ! the boundary conditions in latitude direction, i.e. west and east
    do i=1,ny
    do j=1,nz

       do i_lat = 1, 2
       ! computations at west boundary
       if ( i_lat == 1 ) then
         ix = 1
       else
         ix = nx
       end if

       ! compute level half-way model layer
       if (j > 1) then
          hh = 0.5*(h_m(ix,i,j)+h_m(ix,i,j-1))*1e-3  ! km
       else
          hh = 0.5*h_m(ix,i,1)*1e-3  ! km
       endif

       ! determine the layer
       do iz = 1, bc_lay
         ind = iz
         if (hh < bc_height(iz) ) exit
       enddo
       ! determine interpolation factors
       frac = (bc_height(ind)-hh)/(bc_height(ind)-bc_height(ind-1))
       frac = min( 1.0, max(0.0, frac ))

       !loop over the species
       do k=1,bc_spec
          ind_spec = l_bc_spec(k)
          ind_lat  = l_bc_zone(i)
          if ( ind_spec > 0 ) then
            if (i_lat == 1) then
              bc_west(i,j,ind_spec) =  frac*bc_isak(k,ind-1,ind_lat) + &
                                 (1.0-frac)*bc_isak(k,ind  ,ind_lat) 
            else
              bc_east(i,j,ind_spec) =  frac*bc_isak(k,ind-1,ind_lat) + &
                                 (1.0-frac)*bc_isak(k,ind  ,ind_lat) 
            endif
          endif
       enddo

    !end east-west loop
    enddo

    ! end loop over cells
    enddo
    enddo

    ! loop over the species
    do k = 1, bc_spec
      ! local index:
      ind_spec = l_bc_spec(k)
      ! skip ?
      if ( ind_spec <= 0 ) cycle
      ! set flag:
      spec_filled(ind_spec) = .true.
    end do
    
    ! the aloft concentration
    ! loop over the species
    do k = 1, bc_spec
      ! local index:
      ind_spec = l_bc_spec(k)
      ! skip ?
      if ( ind_spec <= 0 ) cycle
      ! loop over lons:
      do i=1,nx
        ! fill lat row:
        caloft(i,:,nlev+1,ind_spec) = bc_west(:,nz,ind_spec)
      end do
    end do
    
    !
    ! ~~ done
    !
    
    ! ok
    status = 0

  end subroutine


end module LE_Bound_Clim_Isak
