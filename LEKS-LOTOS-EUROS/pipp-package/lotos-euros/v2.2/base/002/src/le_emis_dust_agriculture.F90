!###############################################################################
!
! LE_Emis - LOTOS-EUROS dust emission routines
!
! For a full description of dust emisisons, see the report of the BOP project:
!  "Mineral Dust as a component of particulate matter"
! Download from:
!  http://www.pbl.nl/en/publications/2010/Mineral-Dust-component-particulate-matter
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

module LE_Emis_Dust_Agriculture

    use dims, only : nx,ny, nz,nspec, runF, outF, &
                           nspec
    use LE_LandUse_BCatIon, only : bcat_frac
    use LE_Time, only : day_of_the_week, local_time2
    use indices

    implicit none


    ! --- in/out -----------------------------------

    private

    public  ::  init_agriacti_dust, agriacti_dust_em
    

    ! --- const ------------------------------------

    character(len=*), parameter ::  mname = 'LE_Emis_Dust_Agriculture'

    !real, parameter    :: g = 9.8067       ! gravity constant
    integer, parameter          :: nsize=2

    integer, parameter :: ii_ca =1
    integer, parameter :: ii_mg =2
    integer, parameter :: ii_k =3
    integer, parameter :: ii_na =4

    integer, parameter :: n_aaperiod=5


    ! --- local ------------------------------------

    !integer, parameter :: DP = kind(0.0d0)
    
    integer :: aaperiod(12), aaperiod_len(n_aaperiod)
    integer :: days(12)
    
contains


  subroutine init_agriacti_dust( status )
  
    use LE_Data, only : LE_Data_Enable
    use LE_Country, only : read_countries
    
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/init_agriacti_dust'

    ! --- local ----------------------------------

    integer ::  im

    ! --- begin ----------------------------------

    ! enable data:
    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'rain', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'ustar', status )
    IF_NOTOK_RETURN(status=1)

    ! fill the array with days
    days      = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    ! check on leap year
    if ( mod(runF%yy_s ,4)*4 == runF%yy_s ) then
      days(2)        = days(2) + 1
    endif

    ! set the aaperiod factor
    aaperiod=(/5,5,1,1,2,2,2,3,3,3,4,4/)
    ! where 
    !  1 :  cultivation and seeding period
    !  2 :  growing period 
    !  3 :  harvesting period
    !  4 :  before winter cultivation 
    !  5 :  winter period

    ! define length of aaperiods
    aaperiod_len = 0
    do im=1,12
      aaperiod_len(aaperiod(im)) = aaperiod_len(aaperiod(im)) + days(im)
    enddo
    
    ! Put flag to read country fractions (needed for timezone information )
    read_countries = .true.
    
    ! ok
    status = 0

  end subroutine init_agriacti_dust


  ! ***


  subroutine agriacti_dust_em( emis, t, status )

    use GO             , only : TDate
    use GO             , only : GoMatchValue
    use JAQL_SFC_Soil  , only : rho_mineral_soil
    use JAQL_Particles , only : Particle_Volume_Rp
    use LE_LandUse_Data, only : ilu_arable, lu_fracs      
    use LE_Country     , only : country_map_code, country_map_frac
    use LE_Data        , only : LE_Data_GetPointer

    ! --- in/out ---------------------------------

    real, intent(inout)       ::  emis(nx,ny,nz,nspec)
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/agriacti_dust_em'

    ! --- local ----------------------------------

    real              ::  flux_dust(nsize)
    real              ::  temperatureAir2m, u_1, density, conv
    integer           ::  ilu, i, j, k, yyh, mmh, ddh, hhh, weekday
    logical           ::  SnowCoverage, Precipitation, PrecipitationLast48hours
    character(len=3)  ::  country
    integer           ::  itr, ispec, icountry

    real              ::  delta_emis
    real              ::  max_val

    ! meteo data:
    real, pointer        ::  dens (:,:,:)   ! (lon,lat,1)
    real, pointer        ::  tsurf(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  rain (:,:,:)   ! (lon,lat,1)
    real, pointer        ::  area(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  ustar(:,:,:)   ! (lon,lat,1)


    ! --- begin ----------------------------------

    ! info ..
    write (gol,'("Constructing dust emissions from agricultural activity")'); call goPr

    ! point to meteo data:
    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'rain', rain, status, check_units ='m/s')
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'ustar', ustar, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'area', area, status, check_units ='m2' )
    IF_NOTOK_RETURN(status=1)


    ! info ...
    write (gol,'("<constructing dust emissions>")'); call goPr

    do i=1,nx
      do j=1,ny
        flux_dust(:) = 0.0
        temperatureAir2m = tsurf(i,j,1)
        density = dens(i,j,1)
        SnowCoverage = .False.
        Precipitation = .False.
        !if(snowcover(i,j).gt.0.0) SnowCoverage = .true.
        ! threshold at 0.1 mm/hr = 0.1 * 1.0e-3/3600.0 m/s
        if( rain(i,j,1) > 0.1 * 1.0e-3/3600.0 ) Precipitation= .true.
        PrecipitationLast48hours = .False.
        u_1 = ustar(i,j,1) 

        ! Find maximum country contribution in each cell
        max_val = maxval(country_map_frac(i,j,:))

        ! Match country code
        call goMatchValue( max_val, country_map_frac(i,j,:), icountry, status )
        IF_NOTOK_RETURN(status=1)
        country = country_map_code(icountry)

        ! Find local time for this country      
        call local_time2(country, t%year, t%month, t%day, t%hour, &
                           yyh, mmh, ddh, hhh, weekday, status)
        IF_NOTOK_RETURN(status=1)                           

        ! We only need to know the agricultural emissions from arable land
        ! old landuse data:
        ! for areas  basedefined as arable land the dust contribution from agricultural activity is calculated:
        ! new landuse data:
        ! for areas where landuse is defined as agricultural land the dust contribution from agricultural activity is calculated

        !if(ilu.eq.2.) then - arable/agricultural land
        ilu = ilu_arable
        
        if(lu_fracs(i,j,ilu).gt.0.0) Call agriacti_dust_calc(ilu, temperatureAir2m, SnowCoverage, Precipitation,  &
                                  PrecipitationLast48hours, u_1, density, & 
                                  flux_dust, mmh, weekday, hhh)

        ! fill the emission array:
        conv = 60*1e6*area(i,j,1) !/s -> /min, and g -> ug

        ! add to dust tracers:
        do itr = 1, n_dust
          ispec = ispecs_dust(itr)
          select case ( ispec )
            case (ispec_dust_ff)
              delta_emis=0 !nothing specified yet
              emis(i,j,1,ispec) =  emis(i,j,1,ispec) + delta_emis
            case (ispec_dust_f)
              delta_emis = flux_dust(1) * conv*lu_fracs(i,j,ilu)
              emis(i,j,1,ispec) =  emis(i,j,1,ispec) + delta_emis
            case (ispec_dust_ccc)
              delta_emis = 0
              emis(i,j,1,ispec) =  emis(i,j,1,ispec) + delta_emis
            case (ispec_dust_cc)
              delta_emis = 0
              emis(i,j,1,ispec) =  emis(i,j,1,ispec) + delta_emis
            case (ispec_dust_c)
              delta_emis = flux_dust(2) * conv*lu_fracs(i,j,ilu)
              emis(i,j,1,ispec) =  emis(i,j,1,ispec) + delta_emis
              
            case default
              write (gol,'("unsupported tracer index ",i6)') ispec; call goErr
              TRACEBACK; status=1; return
          end select
        end do
        
        ! base-cat-ions enabled ? just Na is not enough ...
        if ( n_basecation > 5 ) then
          ! add to basecation tracers:
          do itr = 1, n_basecation
            ! tracer index:
            ispec = ispecs_basecation(itr)
            ! switch:
            select case ( ispec )
              case(ispec_Na_ff) ; emis(i,j,1,ispec) = emis(i,j,1,ispec) 
              case(ispec_Na_f)  ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(1) * conv * bcat_frac(i,j,ii_Na)*lu_fracs(i,j,ilu)
              case(ispec_Na_ccc); emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(2) * conv * bcat_frac(i,j,ii_Na)*lu_fracs(i,j,ilu)
              case(ispec_Na_cc) ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(2) * conv * bcat_frac(i,j,ii_Na)*lu_fracs(i,j,ilu)
              case(ispec_Na_c)  ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(2) * conv * bcat_frac(i,j,ii_Na)*lu_fracs(i,j,ilu)
              case(ispec_Ca_f)  ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(1) * conv * bcat_frac(i,j,ii_Ca)*lu_fracs(i,j,ilu)
              case(ispec_Ca_c)  ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(2) * conv * bcat_frac(i,j,ii_Ca)*lu_fracs(i,j,ilu)
              case(ispec_Mg_f)  ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(1) * conv * bcat_frac(i,j,ii_Mg)*lu_fracs(i,j,ilu)
              case(ispec_Mg_c)  ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(2) * conv * bcat_frac(i,j,ii_Mg)*lu_fracs(i,j,ilu)
              case(ispec_K_f )  ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(1) * conv * bcat_frac(i,j,ii_K )*lu_fracs(i,j,ilu)
              case(ispec_K_c )  ; emis(i,j,1,ispec) = emis(i,j,1,ispec) + flux_dust(2) * conv * bcat_frac(i,j,ii_K )*lu_fracs(i,j,ilu)
              case default
                write (gol,'("unsupported tracer index ",i6)') ispec; call goErr
                TRACEBACK; status=1; return
            end select
          end do
        end if

      enddo  ! ny
    enddo ! nx  

    !OK
    status = 0

  end subroutine agriacti_dust_em


  ! ***
  
  
  subroutine agriacti_dust_calc( landUse_Class, temperatureAir2m, SnowCoverage, &
          Precipitation, PrecipitationLast48hours, u, density, &
          flux_dust, im, weekday, hhh)  

      implicit none

      integer, intent(IN) :: landUse_Class

      real :: flux_dust(nsize)
      real,   intent(IN) :: temperatureAir2m
      real,   intent(IN) :: u               ! friction velocity
      real,   intent(IN) :: density         ! density air [kg/m3]       
      logical , intent(IN) :: SnowCoverage, Precipitation, PrecipitationLast48hours

      logical :: DustCondition

      integer :: im, hhh, weekday, iaclass
      integer, parameter :: naclass=5          

     ! agricultural activity categories (iaclass) with (total) emission factors after Ottl et al (2005) 
     ! and Funk et al (2007):
     !     iaclass = 1     :   harrowing
     !     iaclass = 2     :   disking
     !     iaclass = 3     :   cultivating
     !     iaclass = 4     :   ploughing
     ! and with an emission factor for corn harvesting after Gaffney et al (2006)
     !     iaclass = 5     :   harvesting
     ! class(es) for harvesting of other types of crops and/or other activities need to be added ... 
     !     iaclass = 6

     ! total emission factor for all agricultural activities in g/m2 
      ! Original dry emission factor ploughing real, dimension(naclass),parameter ::  emfac=(/0.082,0.137,0.186,1.045,0.191/)
      real, dimension(naclass),parameter ::  emfac=(/0.082,0.137,0.186,0.120,0.191/)

     ! ---------- not used here now -----------
     ! total emission factor of agricultural activity in g/m2.h (as used for the calculation of emission fluxes)
     ! is converted from emfac [g/m2] as follows:
     !    With assumed machine track width and driving speed of 4m, 10km/h respectively 1m2 of land is 
     !    cultivated/harvested within time period delta_t=area/(speed*width)=1m2/(10000m/h*4m)=2.5e-05 h
     !    Therefore, during 1 hour of cultivation 40000m2 of land is cultivated with a corresponding emission
     !    of 190.7*40000 mg =7.628e+03 mg   
     !    The conversion factor for emfac [g/m2]-> emfac [g/m2.h] is thus given by (emfacconv):
     ! real, dimension(naclass), parameter ::  efconv=(/4.0e04,4.0e04,4.0e04,4.0e04,4.0e04/)
     !--------- end of not used here now ------

      !hour factor
      real, dimension(24), parameter :: h_fac=(/0.0,0.0,0.0,0.0,0.0,0.0,&
           0.0,1.6,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,1.6,0.0,&
           0.0,0.0,0.0,0.0,0.0,0.0/)
      !day factor
      real, dimension(7), parameter :: d_fac=(/1.1,1.1,1.1,1.1,1.1,0.75,0.75/)

     ! Additional timing conversion
     ! For each period of agricultural activity  emfac[g/m2] is spread over the entire period of 
     ! activity instead of applying a 'peak emission factor' during a very short period of time 
     ! within the period of activity. Thus, emfac [g/m2] is therefore divided by the activity period 
     ! length (aaperiod_len, the number of days during period, set in the init_agriact routine) in hours: 
     ! aaperiod_len * 24h

     ! Fine fraction corresponding to agricultural activity (total) emission factors   
      real, dimension(naclass), parameter ::  fine_fr=(/0.35,0.09,0.1,0.12,0.1/)

      integer :: isize
      real, dimension(nsize,naclass) ::  mass_fr              

      ! Assign fine and coarse fractions to a mass fraction array
      do isize=1,nsize
       do iaclass=1,naclass
         if (isize.eq. 1) mass_fr(isize,iaclass)=fine_fr(iaclass)
         if (isize.eq. 2) mass_fr(isize,iaclass)=1-fine_fr(iaclass)
       enddo
      enddo

      flux_dust(:)=0.

     ! NoDust condition (cases for which dust from agricultural activity is set to zero)
      DustCondition = .not. (temperatureAir2m < 0. .OR. SnowCoverage .OR. Precipitation .OR. PrecipitationLast48hours)


      If (DustCondition ) Then

        !flux of agriculturally activated dust in g/m2.s
        !product of total agriculaturally activated dust emission factor [g/m2.h] and timing conversion factors
        ! (h_fac, d_fac and conversion from /h -> /s, and aaperiod_len*24)   

        ! Selection of agricultural activity, time and aaperiod dependent

          do isize=1,nsize
          ! spring case
          if (aaperiod(im).eq.1) then                
            do iaclass=1,2   
              flux_dust(isize)=flux_dust(isize) + mass_fr(isize,iaclass) &
                           * emfac(iaclass)*h_fac(hhh+1)*d_fac(weekday)/(3600*aaperiod_len(aaperiod(im))*24)  
            enddo
            iaclass=4
            flux_dust(isize)=flux_dust(isize) + mass_fr(isize,iaclass) &
                           * emfac(iaclass)*h_fac(hhh+1)*d_fac(weekday)/(3600*aaperiod_len(aaperiod(im))*24)
          endif
          ! harvesting case
          if (aaperiod(im).eq. 3) then
            iaclass=5
            flux_dust(isize)=flux_dust(isize) +mass_fr(isize,iaclass) &
                           * emfac(iaclass)*h_fac(hhh+1)*d_fac(weekday)/(3600*aaperiod_len(aaperiod(im))*24)  
          endif
          ! before winter case - not applied right now
          if (aaperiod(im).eq.4) then                
            iaclass=3
            flux_dust(isize)=flux_dust(isize) +mass_fr(isize,iaclass) &
                           * emfac(iaclass)*h_fac(hhh+1)*d_fac(weekday)/(3600*aaperiod_len(aaperiod(im))*24)  
          endif
       ! if (aaperiod(im).eq.1 .or. aaperiod(im).eq.3) then 
       !  print *, 'hhh, h_fac, weekday, d_fac: ',hhh, h_fac(hhh), weekday, d_fac(weekday)
       !  print *, 'size mode (1=fine, 2=coarse): ', isize
       !  print *, 'agricultural activity period: ', aaperiod(im)
       !  print *, 'length of this period (in days): ', aaperiod_len(im)
       !  print *, 'emfactor of period', iaclass, '= ', emfac(iaclass)
       !  print *, 'mass_fr: ', mass_fr(isize,iaclass)
       !  print *, 'agricultural dust flux: ', flux_dust(isize)
       ! endif
        enddo ! isize           

       ! pause
      endif ! DustCondition

  end subroutine agriacti_dust_calc
    

end module LE_Emis_Dust_Agriculture
