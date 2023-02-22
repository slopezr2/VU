!###############################################################################
!
! NAME
!   LE_Bound_Clim_Logan
!
! DESCRIPTION
!
!   For handling the Logan boundary conditions.
!   Also parameters from the emep report for O3
!
! HISTORY
!   Former module 'logan' in 'bound_logan.F90' .
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

module LE_Bound_Clim_Logan

  use GO, only : gol, goErr, goPr
  
  use dims, only : nx, ny, nz, nspec
  use LE_Logging, only : ident1

  implicit none


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_Bound_Clim_Logan'

  real, allocatable   ::  data_logan(:,:,:,:)
  real, allocatable   ::  levels_logan(:)
  real                ::  dx_logan, dy_logan, west_logan, south_logan
  integer             ::  nx_logan, ny_logan, nz_logan

  

contains

  
  ! ==========================================================
  
  
  subroutine LE_Bound_Clim_Logan_Init( path, status )

    use GO     , only : goGetFU
    use indices, only : i_o3
    use LE_Data, only : LE_Data_Enable
    use LE_Grid, only : ugg
    use dims   , only : nx, ny
    
    ! --- in/out ------------------------
    
    character(len=*), intent(in)    ::  path
    integer, intent(out)            ::  status
    
    ! --- const -------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Bound_Clim_Logan_Init'
    
    ! --- local -------------------------

    integer             ::  u_tmp
    character(len=1024) ::  fname
    integer             ::  imonth, ilayer, j
    
    ! --- begin -------------------------

    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    
    select case( trim(ugg%type) )
      case ( 'cartesian-regular', 'cartesian' ) 
        if ( i_o3 > 0 ) then

          !! info ...
          !print *, ident1,'reading logan BC values for O3'

          ! new file unit:
          call goGetFU( u_tmp, status )
          IF_NOTOK_RETURN(status=1)

          fname = trim(path)//'o3_month_logan.txt'
          open (u_tmp,file=trim(fname),status='old')
          !read the dimensions of the logan data
          read (u_tmp,*)
          read (u_tmp,*) nx_logan
          read (u_tmp,*) ny_logan
          read (u_tmp,*) nz_logan
          ! allocate arrays
          allocate(data_logan(nx_logan,ny_logan,nz_logan,12))
          allocate(levels_logan(nz_logan))
          ! read levels
          read (u_tmp,*) levels_logan(1:nz_logan)
          ! read grid spec
          read (u_tmp,*) dx_logan
          read (u_tmp,*) dy_logan
          read (u_tmp,*) west_logan
          read (u_tmp,*) south_logan

          ! check ...
          if ( (ugg%longitude_bnds_1d(1,       1) < west_logan-           0.5 *dx_logan) .or. &
               (ugg%longitude_bnds_1d(2,ugg%nlon) > west_logan +(nx_logan-0.5)*dx_logan) .or. &
               (ugg%latitude_bnds_1d (1,       1) < south_logan-          0.5 *dy_logan) .or. &
               (ugg%latitude_bnds_1d (2,ugg%nlat) > south_logan+(ny_logan-0.5)*dy_logan) ) then
            write (gol,'("LOGAN grid too small for LE grid : ")'); call goErr
            write (gol,'("  logan west , dx, nx : ",2f8.2,i6)') west_logan, dx_logan, nx_logan; call goErr
            write (gol,'("  logan south, dy, ny : ",2f8.2,i6)') south_logan, dy_logan, ny_logan; call goErr
            write (gol,'("  logan w,e,s,n       : ",4f8.2)') &
                west_logan -          0.5 *dx_logan, &
                west_logan +(nx_logan-0.5)*dx_logan, &
                south_logan-          0.5 *dy_logan, &
                south_logan+(ny_logan-0.5)*dy_logan; call goErr
            write (gol,'("  model w,e,s,n       : ",4f8.2)') ugg%longitude_bnds_1d(1,1), ugg%longitude_bnds_1d(2,ugg%nlon), &
                                                             ugg%latitude_bnds_1d(1,1), ugg%latitude_bnds_1d(2,ugg%nlat); call goErr
            TRACEBACK; status=1; return
          end if

          ! do the reading
          do imonth=1,12
             do ilayer=1,nz_logan
                read (u_tmp,*)
                ! the file is ordered from North to South
                do j=ny_logan,1,-1
                   read(u_tmp,*) data_logan(1:nx_logan,j,ilayer,imonth)
                enddo
             enddo
          enddo

          close (u_tmp)

          !write(*,*) 'Logan data level 3, july:'
          !do j=1,ny_logan
          !write(*,*) data_logan(:,j,3,7)
          !enddo
          !print *, 'nx,ny,nzlogan', nx_logan, ny_logan, nz_logan
        endif ! tracer o3?

      case default
        write( gol, '("Conversion not yet made for boundary Logan to universal grid")' ) ; call goErr
        write( gol, '("Make input file netCDF compliant to be used in le-bound-data format")' ) ; call goErr
        TRACEBACK;status=1;return

    end select

    
    ! ok
    status = 0

  end subroutine LE_Bound_Clim_Logan_Init


  ! ***
  
  
  ! purpose: get O3 boundary date from the Logan data set
  ! put it in photostationary equilibrium with NO2 and NO
  ! and insert the results into the BC array.
  !
  !   spec_filled(1:nspec)  : set to .true. if bc is filled, unchanged otherwise
  !

  subroutine get_logan( iyear, imonth, spec_filled, status )

    use dims, only : nspec
    use indices
    use LE_Data_Common , only : nlev
    use LE_Bound_Common, only : bc_west, bc_east, bc_south, bc_north, caloft

    use LE_Data, only : LE_Data_GetPointer
    use LE_Grid, only : ugg
    
    ! --- in/out ----------------------------
    
    integer, intent(in)         ::  iyear, imonth
    logical, intent(inout)      ::  spec_filled(nspec)
    integer, intent(out)        ::  status
    
    ! --- const -------------------------
    
    character(len=*), parameter ::  rname = mname//'/get_logan'
    
    ! ----MACE HEAD correction tables from EMEP model
  ! Mace Head ozone concentrations for backgroudn sectors
  ! from Fig 5.,  Derwent et al., 1998, AE Vol. 32, No. 2, pp 145-157
  real, parameter      ::maceheadlat=53.3
  real, parameter      ::maceheadlon=-9.9
  integer, parameter :: MH_YEAR1 = 1990, MH_YEAR2 = 2012
  real, dimension(12,MH_YEAR1:MH_YEAR2), parameter :: macehead_year=reshape(&
   [35.3,36.3,38.4,43.0,41.2,33.4,35.1,27.8,33.7,36.2,28.4,37.7,& !1990
    36.1,38.7,37.7,45.8,38.8,36.3,29.6,33.1,33.4,35.7,37.3,36.7,& !1991
    36.1,37.3,41.8,39.6,41.2,31.5,28.3,30.3,31.3,34.2,36.1,34.9,& !1992
    37.6,40.4,44.4,42.6,43.4,29.2,28.5,29.6,32.2,37.3,37.3,38.3,& !1993
    38.6,37.3,45.7,43.8,42.9,35.1,30.8,30.5,33.8,36.5,34.0,37.3,& !1994
    37.5,37.1,41.6,42.4,41.1,33.1,29.1,28.7,33.7,34.8,35.0,36.0,& !1995
    37.0,40.1,42.9,44.6,41.3,38.3,29.3,29.4,35.6,38.4,37.8,38.4,& !1996
    36.2,41.9,41.8,40.4,40.6,34.4,26.2,29.3,31.3,35.2,25.7,39.5,& !1997
    38.6,42.0,44.6,45.1,44.2,33.0,29.7,32.9,35.7,38.8,39.7,40.4,& !1998
    39.9,44.5,49.4,45.0,42.8,34.3,29.0,30.0,31.8,36.9,39.6,39.2,& !1999
    39.5,42.1,41.8,43.8,43.4,34.5,28.0,27.3,33.6,37.4,35.6,35.8,& !2000
    37.3,38.0,42.2,44.8,42.6,34.9,28.9,29.4,29.9,35.3,37.3,37.5,& !2001
  ! Preliminary BCs generated using Mace Head CFC and other greenhouse gases
  ! data to define clean air masses. Data cover all of 2002 and 9 months
  ! of 2003. What to do for Oct-Dec 2003?
  ! Could use (1) 2002 data or (2) 10-year average?
  ! Simmonds paper would support (1), simplicity (2).
  ! After seeing earlier 2003 plots, chose (2).
  !!Values in ug/m3
    42.4,44.4,45.5,45.0,45.9,39.8,32.5,28.7,37.7,39.3,40.5,42.3,& !2002
    39.8,40.1,44.7,45.4,45.7,41.7,33.3,31.0,35.7,37.9,40.9,38.1,& !2003
    40.8,42.0,48.3,46.6,39.9,31.9,32.4,32.1,33.9,36.7,40.2,39.8,& !2004
    40.9,41.4,44.1,45.6,42.7,32.9,26.7,30.0,33.2,37.7,39.5,38.0,& !2005
  ! 2006 and 2007 are calculated with using IE31 O3 data and
  ! trajectory sectors (based on PARLAM-PS and HIRLAM20 met) for resp. year
    39.8,42.4,44.2,48.3,41.3,39.0,31.9,29.5,34.8,37.4,41.9,39.9,& !2006
    40.7,38.2,46.1,46.4,40.9,34.5,31.2,28.8,33.3,36.1,40.6,41.7,& !2007
  ! 2008 Mace Head correction calculated using IE31 O3 data and
  ! trajectory sectors (based on HIRLAM20 met) for 2008
    41.0,45.1,48.0,46.3,44.2,37.1,30.8,31.3,34.3,37.5,37.9,40.0,& !2008
  ! 2009 to 2011 Mace Head correction calculated using IE31 O3 data and
  ! trajectory sectors (based on ECMWF met) for respective year
    37.7,43.3,46.5,46.2,41.6,39.1,31.0,29.0,34.5,34.4,40.5,38.4,& !2009
    36.8,38.9,43.9,46.4,41.7,35.5,31.0,31.3,35.6,36.7,33.4,33.8,& !2010
    36.5,42.4,43.3,44.5,40.2,34.6,30.1,30.8,32.0,34.7,37.7,38.1,& !2011
    35.0,40.2,41.0,46.8,43.1,34.0,29.6,33.8,34.9,33.3,37.9,38.7]& !2012
    ,[12,MH_YEAR2-MH_YEAR1+1])
  real, dimension(12), parameter :: macehead_default=&
  ! Defaults from 1998-2010 average
    (/39.8,41.9,45.4,46.5,43.2,36.2,30.5,30.1,34.1,37.0,39.0,38.5/)
  real, dimension(12):: macehead_O3 
     
    
    ! --- local -----------------------------


    integer :: i,j,k,ix,iy,ind,kk
    real    :: press
    real    :: macehead_cor
    !real, parameter::  ppb=1.92 !conversion ug/m3->ppb
    integer  :: maceheadindexlon, maceheadindexlat

    ! meteo data:
    real, pointer    ::  temp(:,:,:)   ! (lon,lat,lev)
    real, pointer    ::  dens(:,:,:)   ! (lon,lat,lev)                                

    ! --- begin -----------------------------

    call LE_Data_GetPointer( 't', temp, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3' )
    IF_NOTOK_RETURN(status=1)

    select case( trim(ugg%type) )
      case ( 'cartesian-regular', 'cartesian' ) 
        if ( i_o3 > 0 ) then

          if((iyear>=MH_YEAR1).and.(iyear<=MH_YEAR2))then
            macehead_O3=macehead_year(:,iyear)
                 !apply mace head correction
             !average over quadrant south west of Mace Head

             maceheadindexlon = int((maceheadlon- west_logan)/dx_logan) +1
             maceheadindexlat = int((maceheadlat-south_logan)/dy_logan) +1
             !print *, 'maeheadindices', maceheadindexlon, maceheadindexlat, maceheadlat, maceheadlon, south_logan, west_logan
             !print *, data_logan(maceheadindexlon, maceheadindexlat, :, imonth)
             macehead_cor=macehead_O3(imonth)-sum(data_logan(1:maceheadindexlon, 1:maceheadindexlat, 1, imonth))/ &
                    (maceheadindexlon*maceheadindexlat)
             !print *, 'sum, elements',  sum(data_logan(1:maceheadindexlon, 1:maceheadindexlat, 1, imonth)) ,   (maceheadindexlon*maceheadindexlat)
             !print *, 'macehead_cor', macehead_cor, macehead_O3(imonth), imonth, iyear
             !print *, 'sectorave', sum(data_logan(1:maceheadindexlon, 1:maceheadindexlat, 1, imonth))/ &
             !       (maceheadindexlon*maceheadindexlat)

          else !
                macehead_cor=0.!macehead_default(imonth)/1.92 !default logan used
          end if
          !print *, ident1,'getting logan b.c. values for O3 for this timestep'

          ! west and east boundaries
          do j=1,ny
          do k=1,nz
             ! west
             !search for appropriate height using dens to get pressure at level
             press = dens(1,j,k)*temp(1,j,k)*8.31451*10.0/29.0
             ! logan levels are descending, so if another index
             ! is found in the test below, that's OK, otherwise the
             ! pressure is higher than the first level and then
             ! the first logan layer is the best choice
             ind = 1
             do kk = 1, nz_logan
                if (press < levels_logan(kk)) ind = kk
             enddo
             ix = int( (ugg%longitude_bnds_1d(1,1)-west_logan)/dx_logan) + 1
             iy = int( (ugg%latitude_bnds_1d(1,1)+(j-0.5)*ugg%dlat-south_logan)/dy_logan) + 1
             bc_west(j,k,i_o3)  = data_logan(ix,iy,ind,imonth)+macehead_cor

             ! east
             press = dens(nx,j,k)*temp(nx,j,k)*8.31451*10.0/29.0
             ind = 1
             do kk = 1, nz_logan
                if (press < levels_logan(kk)) ind = kk
             enddo
             ix = int( (ugg%longitude_bnds_1d(1,1)+nx*ugg%dlon-west_logan)/dx_logan) + 1
             iy = int( (ugg%latitude_bnds_1d(1,1)+(j-0.5)*ugg%dlat-south_logan)/dy_logan) + 1
             bc_east(j,k,i_o3) = data_logan(ix,iy,ind,imonth)+macehead_cor
          enddo
          enddo

          ! south and north boundary
          do i=1,nx
          do k=1,nz
             ! south
             press = dens(i,1,k)*temp(i,1,k)*8.31451*10.0/29.0
             ind = 1
             do kk = 1, nz_logan
                if (press < levels_logan(kk)) ind = kk
             enddo
             ix = int( (ugg%longitude_bnds_1d(1,1)+(i-0.5)*ugg%dlon-west_logan)/dx_logan) + 1
             iy = int( (ugg%latitude_bnds_1d(1,1)-south_logan)/dy_logan) + 1
             bc_south(i,k,i_o3) = data_logan(ix,iy,ind,imonth)+macehead_cor

             ! north
             press = dens(i,ny,k)*temp(i,ny,k)*8.31451*10.0/29.0
             ind = 1
             do kk = 1, nz_logan
                if (press < levels_logan(kk)) ind = kk
             enddo
             ix = int( (ugg%longitude_bnds_1d(1,1)+(i-0.5)*ugg%dlon-west_logan)/dx_logan) + 1
             iy = int( (ugg%latitude_bnds_1d(1,1)+ny*ugg%dlat-south_logan)/dy_logan) + 1
             bc_north(i,k,i_o3) = data_logan(ix,iy,ind,imonth)+macehead_cor
          enddo
          enddo

          ! upper boundary
          do i=1,nx
          do j=1,ny
             press = dens(i,j,nz)*temp(i,j,nz)*8.31451*10.0/29.0
             ind = 1
             do kk = 1, nz_logan
                if (press < levels_logan(kk)) ind = kk
             enddo
             ix = int( (ugg%longitude_bnds_1d(1,1)+(i-0.5)*ugg%dlon-west_logan)/dx_logan) + 1
             iy = int( (ugg%latitude_bnds_1d(1,1)+(j-0.5)*ugg%dlat-south_logan)/dy_logan) + 1
             caloft(i,j,nlev+1,i_o3) = data_logan(ix,iy,ind,imonth)+macehead_cor

          enddo
          enddo


          ! set flag:
          spec_filled(i_o3) = .true.

        end if ! o3 ttracer?
      case default
        write( gol, '("Conversion not yet made for boundary Logan to universal grid")' ) ; call goErr
        write( gol, '("Make input file netCDF compliant to be used in le-bound-data format")' ) ; call goErr
        TRACEBACK;status=1;return

    end select
    
    ! ok
    status = 0

  end subroutine get_logan
  
  
end module LE_Bound_Clim_Logan

