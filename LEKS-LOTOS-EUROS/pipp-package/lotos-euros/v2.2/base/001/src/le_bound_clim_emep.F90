!###############################################################################
!
! NAME
!   LE_Bound_Clim_EMEP
!
! DESCRIPTION
!
!   "Parameters from the emep report for O3."
!
! HISTORY
!   Taken from form 'logan' module in 'bound_logan.F90' .
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

module LE_Bound_Clim_EMEP

  use GO, only : gol, goErr, goPr
  
  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  LE_Bound_Clim_EMEP_Init, get_emep_bc
  

  ! --- const ------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_Clim_EMEP'
  

  ! --- local ------------------------------

  ! emep table
  real, allocatable :: cmean(:), daymax(:), dc(:), hz(:), cmin(:)
  real, allocatable :: latdep(:,:)
  real, allocatable :: latbc(:)

  ! latitute dependance
  integer, parameter :: nlat_bc = 9


contains

  
  ! ==========================================================
  
  
  subroutine LE_Bound_Clim_EMEP_Init( status )

    use dims, only : nspec
    use indices
    use LE_Data, only : LE_Data_Enable
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_Clim_EMEP_Init'

    ! --- local ------------------------------

    ! --- begin -------------------------------------

    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)

    ! first allocate
    allocate(cmean(nspec) )
    allocate(daymax(nspec) )
    allocate(dc(nspec) )
    allocate(hz(nspec) )
    allocate(cmin(nspec) )
    allocate( latbc(nlat_bc) )
    allocate( latdep(nlat_bc,nspec) )

    ! set the parameters for the EMEP BC
    cmean = 0.0
    daymax= 0.0
    dc    = 0.0
    hz    = 0.0
    cmin  = 0.0

    if ( i_so2 > 0 ) then
      cmean(i_so2)  = 0.15
      daymax(i_so2) = 15.0
      dc(i_so2)     = 0.05
      hz(i_so2)     = 1.0e10
      cmin(i_so2)   = 0.03
    endif

    if ( i_no > 0 ) then
      cmean(i_no)  = 0.1
      daymax(i_no) = 15.0
      dc(i_no)     = 0.03
      hz(i_no)     = 4.0
      cmin(i_no)   = 0.02
    end if
    
    if ( i_no2 > 0 ) then
      cmean(i_no2)  = 0.1
      daymax(i_no2) = 15.0
      dc(i_no2)     = 0.03
      hz(i_no2)     = 4.0
      cmin(i_no2)   = 0.04
    end if
    
    if ( i_pan > 0 ) then
      cmean(i_pan)  = 0.2
      daymax(i_pan) = 120.0
      dc(i_pan)     = 0.15
      hz(i_pan)     = 1.0e10
      cmin(i_pan)   = 0.1
    end if
    
    if ( i_hno3 > 0 ) then
      cmean(i_hno3)  = 0.1
      daymax(i_hno3) = 15.0
      dc(i_hno3)     = 0.03
      hz(i_hno3)     = 1.0e10
      cmin(i_hno3)   = 0.05
    end if
    
    if ( i_co > 0 ) then
      cmean(i_co)  = 125.0
      daymax(i_co) = 75.0
      dc(i_co)     = 35.0
      hz(i_co)     = 25.0
      cmin(i_co)   = 30.0
    end if
    
    if ( i_eth > 0 ) then
      cmean(i_eth)  = 2.0
      daymax(i_eth) = 75.0
      dc(i_eth)     = 1.0
      hz(i_eth)     = 10.0
      cmin(i_eth)   = 0.05
    end if
    
    if ( i_form > 0 ) then
      cmean(i_form)  = 0.7
      daymax(i_form) = 180.0
      dc(i_form)     = 0.3
      hz(i_form)     = 6.0
      cmin(i_form)   = 0.05
    end if
    
    if ( i_ald > 0 ) then
      cmean(i_ald)  = 0.7
      daymax(i_ald) = 180.0
      dc(i_ald)     = 0.3
      hz(i_ald)     = 6.0
      cmin(i_ald)   = 0.05
    end if
    
    if ( i_ald2 > 0 ) then
      cmean(i_ald2)  = 0.7
      daymax(i_ald2) = 180.0
      dc(i_ald2)     = 0.3
      hz(i_ald2)     = 6.0
      cmin(i_ald2)   = 0.05
    end if
    
    if ( i_so4a_f > 0 ) then
      ! in our model SO4 is in ug/m3!
      cmean(i_so4a_f)  = 0.15*0.0409*96
      daymax(i_so4a_f) = 180.0
      dc(i_so4a_f)     = 0.0
      hz(i_so4a_f)     = 1.6
      cmin(i_so4a_f)   = 0.03*0.0409*96
    end if
    
    if ( i_nh4a_f > 0 ) then
      ! put NH4 in balance with SO4
      cmean(i_nh4a_f)  = 0.15*0.0409*96*0.375
      daymax(i_nh4a_f) = 180.0
      dc(i_nh4a_f)     = 0.0
      hz(i_nh4a_f)     = 1.6
      cmin(i_nh4a_f)   = 0.03*0.0409*96*0.375
    endif

    ! the array for latitude dependance
    latbc(1:nlat_bc) = (/30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0/)
    ! default:
    latdep = 1.0
    ! overwrite if necessary:
    if ( i_so2    > 0 ) latdep(:,i_so2   ) = (/0.05,0.15,0.30,0.80,1.00,0.60,0.20,0.12,0.05/)
    if ( i_hno3   > 0 ) latdep(:,i_hno3  ) = (/1.00,1.00,1.00,0.85,0.70,0.55,0.40,0.30,0.20/)
    if ( i_pan    > 0 ) latdep(:,i_pan   ) = (/0.15,0.33,0.50,0.80,1.00,0.75,0.50,0.30,0.10/)
    if ( i_co     > 0 ) latdep(:,i_co    ) = (/0.60,0.70,0.80,0.90,1.00,1.00,0.95,0.85,0.80/)
    if ( i_no2    > 0 ) latdep(:,i_no2   ) = (/0.05,0.15,0.30,0.80,1.00,0.60,0.20,0.12,0.05/)
    if ( i_no     > 0 ) latdep(:,i_no    ) = (/0.05,0.15,0.30,0.80,1.00,0.60,0.20,0.12,0.05/)
    if ( i_form   > 0 ) latdep(:,i_form  ) = (/0.15,0.33,0.50,0.80,1.00,0.75,0.50,0.30,0.10/)
    if ( i_ald    > 0 ) latdep(:,i_ald   ) = (/0.15,0.33,0.50,0.80,1.00,0.75,0.50,0.30,0.10/)
    if ( i_ald2   > 0 ) latdep(:,i_ald2  ) = (/0.15,0.33,0.50,0.80,1.00,0.75,0.50,0.30,0.10/)
    if ( i_so4a_f > 0 ) latdep(:,i_so4a_f) = (/0.05,0.15,0.30,0.80,1.00,0.60,0.20,0.12,0.05/)
    if ( i_nh4a_f > 0 ) latdep(:,i_nh4a_f) = (/0.05,0.15,0.30,0.80,1.00,0.60,0.20,0.12,0.05/)
    if ( i_so2    > 0 ) latdep(:,i_so2   ) = (/0.05,0.15,0.30,0.80,1.00,0.60,0.20,0.12,0.05/)
    
    ! ok
    status = 0

  end subroutine LE_Bound_Clim_EMEP_Init


  ! ***
  
  
  ! Put NO and NO2 boundary condition data from EMEP report
  !
  !   spec_filled(1:nspec)  : set to .true. if bc is filled, unchanged otherwise
  !

  subroutine get_emep_bc(  iyear,iday, spec_filled, status )

    use GO             , only : goc
    use LE_Logging     , only : ident1
    use dims           , only : nx, ny, nz, nspec
    use indices
    use LE_Data        , only : LE_Data_GetPointer
    use LE_Data_Common , only : nlev
    use LE_Bound_Common, only : bc_west, bc_east, bc_south, bc_north, caloft
    use LE_Grid        , only : ugg
    
    ! --- in/out ----------------------------
    
    integer, intent(in)         ::  iday, iyear
    logical, intent(inout)      ::  spec_filled(nspec)
    integer, intent(out)        ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/get_emep_bc'
    
    ! --- local -----------------------------

    integer :: j, ilat, k, ispec, i, i0,i1
    real :: hh, arg, cosarg, topos, topon, topow, topoe, f0, f1

    real    :: hhe, hhw, hhn, hhs, delta, delta_s, delta_n, top
    integer :: ind, ind_n, ind_s

    real, pointer        ::  h_m(:,:,:)   ! (lon,lat,nz)
    
    !----------------------------------------------------------

  !Trends 1980-2003 derived from EPA emissions of so2,nox.
  ! nh4 derived from 2/3so3+1/3nox
  !Support for SO2 can be found in Hicks, Artz, Meyer and Hosker, 2002
  ! Figure 7 (Eastern US) which show 'close' correspondance between national
  ! emissions and concentration trend

  !1920-1970 BCs derived from:
  !  NH4: nh3 emissions
  !  SOx: winter ice cores, Col du dome
  !  NOx: winter ice cores
  !1890-1920: trends from emissions for SOx,NOx,NH3, Aardenne USA
  ! Updated: April 2013
  ! - use data above to 1980, then EPA download of April 2013
  ! - then IIASA/ECLAIRE/ECLIPSE
  type :: SIAfac ! trends in boundary conditions
    integer :: year
    real:: so2,nox,nh4
  end type SIAfac
   type(SIAfac), dimension(37), save :: SIAtrends = (/ &
      SIAfac(1890,   0.12 ,   0.15 ,   0.44 )&
     ,SIAfac(1900,   0.18 ,   0.20 ,   0.48 )&
     ,SIAfac(1910,   0.27 ,   0.27 ,   0.52 )&
     ,SIAfac(1920,   0.32 ,   0.33 ,   0.59 )&
     ,SIAfac(1930,   0.35 ,   0.33 ,   0.55 )&
     ,SIAfac(1940,   0.46 ,   0.25 ,   0.59 )&
     ,SIAfac(1950,   0.59 ,   0.33 ,   0.69 )&
     ,SIAfac(1960,   0.76 ,   0.50 ,   0.76 )&
     ,SIAfac(1970,   0.95 ,   0.75 ,   0.90 )&
     ,SIAfac(1980,   1.000,   1.000,   1.000)&     
     ,SIAfac(1985,   0.899,   0.951,   0.989)&
     ,SIAfac(1990,   0.890,   0.943,   0.920)&
     ,SIAfac(1991,   0.863,   0.930,   0.934)&
     ,SIAfac(1992,   0.852,   0.933,   0.947)&
     ,SIAfac(1993,   0.840,   0.936,   0.963)&
     ,SIAfac(1994,   0.823,   0.936,   0.978)&
     ,SIAfac(1995,   0.718,   0.922,   0.993)&
     ,SIAfac(1996,   0.709,   0.915,   1.007)&
     ,SIAfac(1997,   0.727,   0.912,   1.027)&
     ,SIAfac(1998,   0.731,   0.899,   1.052)&
     ,SIAfac(1999,   0.677,   0.844,   1.035)&
     ,SIAfac(2000,   0.631,   0.835,   1.046)&
     ,SIAfac(2001,   0.615,   0.796,   0.786)&
     ,SIAfac(2002,   0.570,   0.781,   0.880)&
     ,SIAfac(2003,   0.568,   0.753,   0.877)&
     ,SIAfac(2004,   0.565,   0.726,   0.874)&
     ,SIAfac(2005,   0.572,   0.703,   0.870)&
     ,SIAfac(2006,   0.514,   0.681,   0.885)&
     ,SIAfac(2007,   0.456,   0.658,   0.900)&
     ,SIAfac(2008,   0.399,   0.635,   0.930)&
     ,SIAfac(2009,   0.320,   0.579,   0.928)&
     ,SIAfac(2010,   0.292,   0.543,   0.925)&
     ,SIAfac(2011,   0.265,   0.488,   0.921)&
     ,SIAfac(2012,   0.213,   0.421,   0.917)&
     ! Default here from IIASA ECLAIRE/ECLIPSE
     ! related to 2005 emissions as base
     ! (Created by mk.UStrends, April 2013)
     ,SIAfac(2030,   0.155,   0.252,   0.953)&
     ,SIAfac(2050,   0.225,   0.276,   0.977)&! Last year which works
     ,SIAfac(2200,   0.225,   0.276,   0.977)&! FAKE for interp
   /)
    type(SIAfac), save :: SIAtrend

    ! --- begin -----------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 'h', h_m, status, check_units='m' )
    IF_NOTOK_RETURN(status=1)

    if (iyear < SIAtrends(1)%year .or. iyear >= SIAtrends(37)%year ) then
      print  *, 'Unspecified trend BCs for this year:', iyear
      stop
    end if

  ! Interpolate between boundary condition years if needed
    i0 = 1
    do i = 1, size(SIAtrends(:)%year) 
      if ( iyear >= SIAtrends(i)%year ) i0 = i
      !if(MasterProc) print "(a,5i6)", "USAsrch: ", i, i0, BCtrend(i)%year, BCtrend(i0)%year
    end do

    i1= i0 + 1
    f0 =     (SIAtrends(i1)%year - iyear)/&
         real(SIAtrends(i1)%year - SIAtrends(i0)%year )
    f1 =     (iyear        - SIAtrends(i0)%year )/&
         real(SIAtrends(i1)%year - SIAtrends(i0)%year )

    SIAtrend%so2 =f0*SIAtrends(i0)%so2 + f1*SIAtrends(i1)%so2
    SIAtrend%nox =f0*SIAtrends(i0)%nox + f1*SIAtrends(i1)%nox
    SIAtrend%nh4 =f0*SIAtrends(i0)%nh4 + f1*SIAtrends(i1)%nh4

    ! !default for 1990-now
    !trend_co=1.0
    !trend_voc=1.0
    !if (iyear.lt.1990) then
    !   trend_co = exp(-0.01*0.85*(1990-iyear)) ! Zander:CO
     !  trend_voc= exp(-0.01*0.85*(1990-iyear)) ! Zander,1975-1990
    !end if

    !print *, ident1,'getting EMEP b.c. values for several species'

    ! Check type of input grid
    select case ( trim(ugg%type) )

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'cartesian-regular', 'cartesian' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! topomax
        top = maxval( h_m(:,:,nz) * 1e-3 )  ! km
        ! maximum over all domains:
        call goc%AllReduce( 'max', top, status )
        IF_NOTOK_RETURN(status=1)

        !
        ! west and east
        !
        
        ! loop from north to south:
        do j = 1, ny

          ! get the latitude index
          ind = -1
          do ilat = 1,nlat_bc
            delta  = abs(ugg%latitude_bnds_1d(1,1) + (j-0.5)*ugg%dlat - latbc(ilat))
            if (delta <=  2.5) then
              ind  = ilat
              exit
            end if
          end do
          if ( ind == -1 ) then
            write (gol,*) 'wrong index west/east'; call goErr
            write (gol,*) j, ugg%latitude_bnds_1d(1,1) + (j-0.5)*ugg%dlat; call goErr
            write (gol,*) latbc(1:nlat_bc); call goErr
            TRACEBACK; status=1; return
          end if

          ! loop over layers:
          do k = 1, nz

            ! get (by a trick) the topography, i..
            topow = top - h_m( 1,j,nz) * 1e-3  ! km
            topoe = top - h_m(nx,j,nz) * 1e-3  ! km
            ! get the height
            if ( k == 1 ) then
              hhw = 0.5*h_m( 1,j,1)*1e-3 + topow  ! km
              hhe = 0.5*h_m(nx,j,1)*1e-3 + topoe  ! km
            else
              hhw = 0.5*(h_m( 1,j,k)+h_m( 1,j,k-1))*1e-3 + topow  ! km
              hhe = 0.5*(h_m(nx,j,k)+h_m(nx,j,k-1))*1e-3 + topoe  ! km
            end if

            ! apply formula for all relevant species
            do ispec = 1, nspec
              ! only if a mean larger than zero is specified
              if ( cmean(ispec) > 0.0 ) then
                ! argument for the cosine - used in BC computations for NO and NO2
                arg = 4.0*atan(1.0)*(iday-daymax(ispec))/365.0
                cosarg = cos(arg)
                bc_west(j,k,ispec) = (cmean(ispec) + dc(ispec)*cosarg)*exp(-hhw/hz(ispec))*latdep(ind,ispec)
                bc_east(j,k,ispec) = (cmean(ispec) + dc(ispec)*cosarg)*exp(-hhe/hz(ispec))*latdep(ind,ispec)
              end if
            end do ! ispec

          end do  ! k

        end do ! j

        !
        ! south and north
        !

        ! get latitude indices for north and south boundary
        ind_s = -1
        ind_n = -1
        do ilat = 1,nlat_bc
           delta_s  = abs(ugg%latitude_bnds_1d(1,1) + (j-0.5)*ugg%dlat - latbc(ilat))
           delta_n  = abs(ugg%latitude_bnds_1d(1,1) + ny*ugg%dlat - latbc(ilat))
           if (delta_s <=  2.5) ind_s  = ilat
           if (delta_n <=  2.5) ind_n  = ilat
        enddo
        if ( ind_s == -1 .or. ind_n == -1 ) then
          write (gol,*) 'wrong delta_s/n'; call goErr
          TRACEBACK; status=1; return
        end if

        ! loop from west to east:
        do i = 1, nx
          ! loop over layers:
          do k = 1, nz

            ! get (by a trick) the topography, i..
            topos = top - h_m(i, 1,nz)*1e-3  ! km
            topon = top - h_m(i,ny,nz)*1e-3  ! km
            if (k==1) then
              hhs = 0.5*h_m(i, 1,1)*1e-3 + topos  ! km
              hhn = 0.5*h_m(i,ny,1)*1e-3 + topon  ! km
            else
              hhs = 0.5*(h_m(i, 1,k)+h_m(i, 1,k-1))*1e-3 + topos  ! km
              hhn = 0.5*(h_m(i,ny,k)+h_m(i,ny,k-1))*1e-3 + topon  ! km
            end if

            ! apply formula for all relevant species
            do ispec=1,nspec
              ! only if a mean larger than zero is specified
              if (cmean(ispec) > 0.0) then
                !argument for the cosine - used in BC computations for NO and NO2
                arg = 4.0*atan(1.0)*(iday-daymax(ispec))/365.0
                cosarg = cos(arg)
                bc_south(i,k,ispec) = (cmean(ispec) + dc(ispec)*cosarg)*exp(-hhs/hz(ispec))*latdep(ind_s,ispec)
                bc_north(i,k,ispec) = (cmean(ispec) + dc(ispec)*cosarg)*exp(-hhn/hz(ispec))*latdep(ind_n,ispec)
              end if

            end do  ! ispec

          end do  ! k
        end do  ! i

        !
        ! the upper boundary arrays
        !

        ! height is the maximum height (this corrects for topography)
        hh = maxval( h_m(:,:,nz)*1e-3 ) ! km
        ! maximum over all domains:
        call goc%AllReduce( 'max', hh, status )
        IF_NOTOK_RETURN(status=1)
        
        ! loop from south to north:
        do j = 1 ,ny

          ! get latitude index
          ind = -1
          do ilat = 1, nlat_bc
             delta  = abs(ugg%latitude_bnds_1d(1,1) + (j-0.5)*ugg%dlat - latbc(ilat))
             if (delta <=  2.5) then
                ind  = ilat
                exit
             end if
          end do
          if (ind == -1) then
             print *, 'wrong index caloft'
             print *, j, ugg%latitude_bnds_1d(1,1) + (j-0.5)*ugg%dlat
             stop 'wrong ind'
          endif

          ! loop from west to east:
          do i = 1, nx
            ! loop over all specs:
            do ispec = 1, nspec
              ! boundary conditions defined?
              if ( cmean(ispec) > 0.0 ) then
                ! time within year:
                arg = 4.0*atan(1.0)*(iday-daymax(ispec))/365.0
                cosarg = cos(arg)
                ! fill:
                caloft(i,j,nlev+1,ispec) = ( cmean(ispec) + dc(ispec)*cosarg ) * exp(-hh/hz(ispec)) * latdep(ind,ispec)
              end if
            end do  ! ispec
          end do ! i

        end do  ! j
        
        
        !
        ! apply trends
        !

        ! account for trends in species 
        do ispec=1,nspec
          ! switch:
          if ( (ispec == i_so2) .or. (ispec == i_so4a_f) ) then
            bc_west (:,:, ispec) = bc_west (:,:,ispec)*SIAtrend%so2
            bc_east (:,:, ispec) = bc_east (:,:,ispec)*SIAtrend%so2
            bc_south(:,:, ispec) = bc_south(:,:,ispec)*SIAtrend%so2
            bc_north(:,:, ispec) = bc_north(:,:,ispec)*SIAtrend%so2
            caloft  (:,:,nlev+1,ispec) = caloft  (:,:,nlev+1,ispec)*SIAtrend%so2
          else if (ispec == i_nh4a_f) then
            bc_west (:,:, ispec) = bc_west (:,:,ispec)*SIAtrend%nh4
            bc_east (:,:, ispec) = bc_east (:,:,ispec)*SIAtrend%nh4
            bc_south(:,:, ispec) = bc_south(:,:,ispec)*SIAtrend%nh4
            bc_north(:,:, ispec) = bc_north(:,:,ispec)*SIAtrend%nh4
            caloft  (:,:,nlev+1,ispec) = caloft  (:,:,nlev+1,ispec)*SIAtrend%nh4

          else if ( ispec == i_no3a_f .or. ispec == i_hno3 .or. &
                    ispec == i_no2    .or. ispec == i_no   .or. &
                    ispec == i_pan ) then
            bc_west (:,:, ispec) = bc_west (:,:,ispec)*SIAtrend%nox
            bc_east (:,:, ispec) = bc_east (:,:,ispec)*SIAtrend%nox
            bc_south(:,:, ispec) = bc_south(:,:,ispec)*SIAtrend%nox
            bc_north(:,:, ispec) = bc_north(:,:,ispec)*SIAtrend%nox
            caloft  (:,:,nlev+1,ispec) = caloft  (:,:,nlev+1,ispec)*SIAtrend%nox
           !trends in co, voc only pre 1990     
           !else if (ispec ==  i_co) then
           !     bc_west(:,:, ispec) = bc_west(:,:,ispec)*trend_co
           !     bc_east(:,:, ispec) = bc_east(:,:,ispec)*trend_co
           !     bc_south(:,:, ispec) = bc_south(:,:,ispec)*trend_co
           !     bc_north(:,:, ispec) = bc_north(:,:,ispec)*trend_co
           !     caloft(:,:,nlev+1,ispec) = caloft(:,:,nlev+1,ispec)*trend_co
           !else if (ispec == i_form.or.ispec == i_ald.or.ispec ==  i_ald2  ) then
           !     bc_west(:,:, ispec) = bc_west(:,:,ispec)*trend_voc
           !     bc_east(:,:, ispec) = bc_east(:,:,ispec)*trend_voc
           !     bc_south(:,:, ispec) = bc_south(:,:,ispec)*trend_voc
           !     bc_north(:,:, ispec) = bc_north(:,:,ispec)*trend_voc
           !     caloft(:,:,nlev+1,ispec) = caloft(:,:,nlev+1,ispec)*trend_voc
          end if
        end do
        
        
        !
        ! constraint on the boundary conditions
        !
        
        ! loop over specs:
        do ispec = 1, nspec
          ! only if species taken into account in this approach
          if ( cmean(ispec) > 0.0 ) then
            ! apply minimum:
            where ( bc_west (:,:,ispec) < cmin(ispec) ) bc_west (:,:,ispec) = cmin(ispec)
            where ( bc_east (:,:,ispec) < cmin(ispec) ) bc_east (:,:,ispec) = cmin(ispec)
            where ( bc_south(:,:,ispec) < cmin(ispec) ) bc_south(:,:,ispec) = cmin(ispec)
            where ( bc_north(:,:,ispec) < cmin(ispec) ) bc_north(:,:,ispec) = cmin(ispec)
            where ( caloft  (:,:,nlev+1,ispec) < cmin(ispec) ) caloft  (:,:,nlev+1,ispec) = cmin(ispec)
            ! set flag:
            spec_filled(ispec) = .true.
            !! info ...
            !print *, '    -->  filled emep bc for ', ispec, ' ', trim(specname(ispec))
            
          end if  ! defined?
        end do ! ispec


      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
        write( gol, '("Conversion not yet made for boundary Emep-climatology to universal grid")' ) ; call goErr
        write( gol, '("Make input file netCDF compliant to be used in le-bound-data format")' ) ; call goErr
        TRACEBACK;status=1;return
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    end select
    
    ! ok
    status = 0

  end subroutine get_emep_bc


end module LE_Bound_Clim_EMEP

