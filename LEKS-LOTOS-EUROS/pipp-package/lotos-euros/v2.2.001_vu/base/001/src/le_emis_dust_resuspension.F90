!###############################################################################
!
! LE_Emis - LOTOS-EUROS dust dust resuspension routines
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
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_Dust_Resuspension

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private

  public  ::  LE_Emis_Dust_Resuspension_Init
  public  ::  roadresp_dust_em
  public  ::  Get_SWG_Data
  public  ::  correction_swg


  ! --- const --------------------------

  character(len=*), parameter   ::  mname = 'LE_Emis_Dust_Resuspension'
  

  ! --- var -----------------------------
  
  real, allocatable     ::  swg_avg(:,:)  ! average gravimetric soil moisture map



contains


  ! ====================================================================
  
  
  subroutine LE_Emis_Dust_Resuspension_Init( rcF, rckey, status )

    use GO                , only : TrcFile, ReadRc
    use Dims              , only : nx, ny, nspec
    use LE_Landuse_Traffic, only : n_traffic_map
    use LE_Country        , only : read_countries

    use LE_Data           , only : LE_Data_Enable
    
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)       ::  rcF
    character(len=*), intent(in)    ::  rckey
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_Emis_Dust_Resuspension_Init'

    ! --- local ----------------------------------
    
    character(len=512)    ::  fname

    ! --- begin ----------------------------------
    
    ! check maps ...
    if ( n_traffic_map < 1 ) then
      write (gol,'("no traffic maps loaded; check if:")'); call goErr
      write (gol,'("  - module Landuse_Traffic is initialized;")'); call goErr
      write (gol,'("  - file(s) with traffic maps are specified in rcfile")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! storage for soil water volume:    
    allocate( swg_avg(nx,ny) )

    ! input file:
    call ReadRc( rcF, trim(rckey)//'.swg.file', fname, status )
    IF_NOTOK_RETURN(status=1)
    ! read, check expected units:
    call Get_SWG_Data( trim(fname), 'kg kg-1', swg_avg, status )
    IF_NOTOK_RETURN(status=1)

    ! to correct for cells with low average swv, i.e. coastal cells
    call correction_swg( swg_avg, status )
    IF_NOTOK_RETURN(status=1)
    
    ! Enable data
    call LE_Data_Enable( 'rain', status )
    IF_NOTOK_RETURN(status=1)

    ! Put flag to read country fractions (needed for timezone information )
    read_countries = .true.

    ! ok
    status = 0

  end subroutine LE_Emis_Dust_Resuspension_Init
  
  
  ! ***
  

  subroutine Get_SWG_Data( fname, units, swg, status )

    use UDUnits       , only : UDUNITS_NOERR, UDUnits_StrError
    use UDUnits       , only : UDUnits_ConversionFactor
    use LE_Grid, only : ugg
    use C3PO   , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors
    
    ! --- in/out ------------------------
    
    character(len=*), intent(in)    ::  fname
    character(len=*), intent(in)    ::  units
    real, intent(out)               ::  swg(:,:)
    integer, intent(out)            ::  status

    ! --- const --------------------------

    character(len=*), parameter   ::  rname = mname//'/Get_SWG_Data'

    character(len=32), parameter  ::    description='var_name=swg'
    
    ! --- local --------------------------

    integer               ::  varid
    
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    real, allocatable     ::  swg_in(:,:)
    character(len=64)     ::  units_in
    real                  ::  factor
    
    ! --- begin -------------------------
    
    ! open file:
    call file_in%Open( trim(fname), status )
    IF_NOTOK_RETURN(status=1)

    ! variable id:
    call file_in%Inq_VarID( trim(description), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! init grid definition
    call file_in%Get_Grid( varid, grid_in, status )
    IF_NOTOK_RETURN(status=1)

    ! storage:
    allocate( swg_in(grid_in%nlon,grid_in%nlat), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! read:
    call file_in%Get_Var( trim(description), swg_in, units_in, status, &
                           start=(/1,1/), count=(/grid_in%nlon,grid_in%nlat/) )
    IF_NOTOK_RETURN(status=1)                     

    ! check ...
    if ( trim(units_in) /= trim(units) ) then
      ! get conversion factor:
      call UDUnits_ConversionFactor( trim(units_in), trim(units), factor, status )
      if ( status /= UDUNITS_NOERR ) then
        gol=trim(UDUnits_StrError(status)); call goErr
        write (gol,'("could not obtain conversion factor from `",a,"` to `",a,"`")') &
                         trim(units_in), trim(units); call goErr
        TRACEBACK; status=1; return
      end if
      ! apply:
      swg_in = swg_in * factor
    end if

    ! regrid:
    call Grid_Convertors%Ugg_AreaAver( grid_in, swg_in, ugg, swg, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( swg_in )
    call grid_in%Done( status )
    IF_NOTOK_RETURN(status=1)
    ! clear:
    call file_in%Close( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Get_SWG_Data


  ! ***


  ! this routine adds resuspension dust emissions to the emissions array

  subroutine roadresp_dust_em( emis, t, status )

    use GO             , only : TDate
    use GO             , only : GoMatchValue
    use LE_Time        , only : local_time2
    use Dims           , only : nx, ny, nz, nspec
    use LE_Data        , only : LE_Data_GetPointer
    use Dims           , only : substantial_rain__mps
    use indices
    use JAQL_SFC_Soil     , only : rho_mineral_soil
    use JAQL_Particles    , only : Particle_Volume_Rp
    use LE_Landuse_Traffic, only : n_traffic_map, traffic_map
    use LE_Landuse_Traffic, only : ROADTYPE_HIGHWAYPLUS
    use LE_Landuse_Traffic, only : ROADTYPE_NONHIGHWAY
    use LE_Landuse_Traffic, only : ROADTYPE_MAX
    use LE_Landuse_Traffic, only : VEHICLETYPE_TRUCK
    use LE_Landuse_Traffic, only : VEHICLETYPE_CAR
    use LE_Landuse_Traffic, only : VEHICLETYPE_MOTORCYCLE
    use LE_Landuse_Traffic, only : VEHICLETYPE_MAX
    use LE_LandUse_BCatIon, only : bcat_frac
    use LE_Country        , only : country_map_code, country_map_frac

    ! --- in/out ---------------------------------

    real, intent(inout)       ::  emis(nx,ny,nz,nspec)
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/windblwn_dust_em'

    real, dimension(24),parameter :: h_fac = (/0.19,0.09,0.06,0.05, &
          0.09,0.22,0.86,1.84,1.86,1.41,1.24,1.2,1.32,1.44,1.45, &
          1.59,2.03,2.08,1.51,1.06,0.74,0.62,0.61,0.44/)
    real, dimension(7), parameter :: d_fac=(/1.02,1.06,1.08,1.1,1.14,0.81,0.79/)

    !Definitions for regional dependence

    !rc and constant applied in multiplication factor (clim_fac) applied  
    !to correct for regional difference in climatic condtions based upon soil moisture 
    !Based upon hmix ratio of 750/550, pm10conc_obs ratio 5.0/2.1 and soil 
    !moisture ratio 0.225/0.3 for southern European region (SEU) versus central European 
    !region (CEU), the best guess for rc and const are:
    real, parameter :: rc=31.00          ! rc for soil moisture correlation with 
                                         ! emission
    real, parameter :: const=10.3        ! constant for soil moisture correlation with 
                                         ! emission

    real, parameter :: swv_ref1=0.3
    real, parameter :: swv_ref2=0.225

    ! month dependent factor for emission intensity due to road sanding
    real, parameter :: rsfactor(12) = (/2,2,2,2,1,1,1,1,1,1,1,1/)

    integer, parameter :: ii_Ca = 1
    integer, parameter :: ii_Mg = 2
    integer, parameter :: ii_K  = 3
    integer, parameter :: ii_Na = 4

    ! --- local ----------------------------------

    real                    ::  emfac_road(ROADTYPE_MAX)
    real                    ::  finefr_road(ROADTYPE_MAX)
    real                    ::  emfac_vehicle(VEHICLETYPE_MAX)

    real                    ::  emis_f, emis_c

    !regional climatic factor and correction factor for sanding and studded tyre usage
    real, dimension(nx,ny) :: clim_fac
    real, dimension(nx,ny) :: cnordic

    integer :: i,j,k
    integer :: hhh,yyh,mmh,ddh,weekday      !local time definitions


    integer                 ::  roadtype
    integer                 ::  vehicletype
    real                    ::  emis_pm
    character(len=3)        ::  country
    integer                 ::  itr, ispec, icountry

    real                    ::  max_val
    real                    ::  delta_emis

    ! meteo data:
    real, pointer        ::  rain(:,:,:)   ! (lon,lat,1)

    ! --- begin ----------------------------------

    !! info ...
    !call wrtgol( rname//': road dust resuspension emission for ', t ); call goPr

    ! point to meteo data:
    call LE_Data_GetPointer( 'rain', rain, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)
      
    ! emission factors (mg/km)
    !
    !          higway         rural        urban
    !         heavy light  heavy light  heavy  light
    ! emfac   198.0, 22.0, 432.0, 48.0, 432.0, 48.0 
    ! fine_fr   0.1,  0.1,   0.1,  0.1,   0.1,  0.1 
    !
    ! emission factor per road type:
    emfac_road(ROADTYPE_HIGHWAYPLUS) = 220.0  ! mg/km
    emfac_road(ROADTYPE_NONHIGHWAY ) = 480.0  ! mg/km
    !
    ! emission factor per vehicle type:
    emfac_vehicle(VEHICLETYPE_TRUCK     ) = 0.90   ! scale factor
    emfac_vehicle(VEHICLETYPE_CAR       ) = 0.10   ! scale factor
    emfac_vehicle(VEHICLETYPE_MOTORCYCLE) = 0.10   ! scale factor
    !
    ! fine fraction part per road type:
    finefr_road(ROADTYPE_HIGHWAYPLUS) = 0.1   ! fraction
    finefr_road(ROADTYPE_NONHIGHWAY ) = 0.1   ! fraction

    ! traffic resuspension  emissions
    do i=1,nx
      do j=1,ny

        ! skip if rains:
        ! units of rain and substantial_rain are both m/s:
        !       m/s     
        if ( rain(i,j,1) >= substantial_rain__mps ) cycle

        !initialize regional factor (clim_fac) and road sanding factor (cnordic)
        clim_fac(i,j)=1.0
        cnordic(i,j)=1.0

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
        
        ! -----------------------------------------------------------
        ! temporary !! to test whether clim_fac application works..
        !    if (lat .lt. 41) then
        !      swg_avg(i,j)=0.225
        !    elseif (lat .lt. 55 .and. lat .ge. 41) then 
        !      swg_avg(i,j)=0.3
        !    elseif (lat .ge. 55) then
        !      swg_avg(i,j)=0.325
        !    endif
        ! -----------------------------------------------------------

        ! set climatic factor clim_fac(i,j)
        !print*, 'i=',i, 'j= ',j,'swg_avg',swg_avg(i,j)
        if (swg_avg(i,j) .ge. swv_ref1) clim_fac(i,j)=const - rc*swv_ref1  ! constant value = 1
        if (swg_avg(i,j) .le. swv_ref2)  clim_fac(i,j)=const - rc*swv_ref2 ! constant value at upper level

        if (swg_avg(i,j) .gt. swv_ref2 .and. swg_avg(i,j) .lt. swv_ref1) clim_fac(i,j)=const - rc*swg_avg(i,j)

        if (clim_fac(i,j) .lt. 1.)  clim_fac(i,j)=1 ! as we do not extrapolate above the threshold SWC


        ! set nordic factor cnordic(i,j)    
        if(country.eq.'SWE'.OR.country.eq.'NOR'.OR.country.eq.'FIN') then
           cnordic(i,j)=rsfactor(mmh)
        endif
        !if (lat .le. 55.0 .and. elevation.ge 1200) then
        !   cnordic(i,j)=rsfactor(mmh)
        !endif

        ! loop over traffic maps:
        do k = 1, n_traffic_map
       
          ! extract:
          roadtype    = traffic_map(k)%roadtype
          vehicletype = traffic_map(k)%vehicletype

          ! tabel-traffic_km map coupling
          ! total pm emission:
          emis_pm =                             &  ! ug/min =
              traffic_map(k)%km(i,j) * 1.e6        &  ! (mlj km)/year  km/(mlj km)
                * emfac_road(roadtype) * 1.e3      &  ! mg/km    ug/mg
                * emfac_vehicle(vehicletype)       &  ! 1
                / (8760*60)                        &  ! year/min
                * h_fac(hhh+1) * d_fac(weekday)    &  ! 1
                * cnordic(i,j) * clim_fac(i,j)        ! 1
          ! distribute over tracers:
          emis_f = emis_pm *     finefr_road(roadtype) 
          emis_c = emis_pm *(1.0-finefr_road(roadtype))

          ! add to dust tracers:
          do itr = 1, n_dust
            ispec = ispecs_dust(itr)
            select case ( ispec )
              case (ispec_dust_ff)
                 !nothing specified yet
              case (ispec_dust_f)
                emis(i,j,1,ispec) =  emis(i,j,1,ispec) + emis_f
              case (ispec_dust_ccc) 
                 !nothing specified yet
              case (ispec_dust_cc) 
                 !nothing specified yet
              case (ispec_dust_c)
                emis(i,j,1,ispec) =  emis(i,j,1,ispec) + emis_c

#ifdef with_m7
              ! M7 species:for mass take usual convention, 
              !   for number convert using density in kg/m3, radius in m
              !   and take into account that emissions are in ug (=1e-9 kg)           
              case (ispec_duai)
                delta_emis = emis_f 
                emis(i,j,1,ispec ) = emis(i,j,1,ispec ) + delta_emis !mass 
                emis(i,j,1,i_naci) = emis(i,j,1,i_naci) + delta_emis*1e-9/(Particle_Volume_Rp(0.17e-6)*rho_mineral_soil)
              case (ispec_duci)
                delta_emis = emis_c
                emis(i,j,1,ispec ) = emis(i,j,1,ispec ) + delta_emis !mass
                emis(i,j,1,i_ncoi) = emis(i,j,1,i_ncoi) + delta_emis*1e-9/(Particle_Volume_Rp(2.0e-6)*rho_mineral_soil)
              case (ispec_duas) 
                ! skip, assume all wind-blown dust insoluble
                delta_emis = 0.0
              case (ispec_ducs) 
                ! skip, assume all wind-blown dust insoluble
                delta_emis = 0.0
#endif                
              case default
                write (gol,'("unsupported tracer index ",i6)') ispec; call goErr
                TRACEBACK; status=1; return
            end select
          end do

          ! base-cat-ions enabled ? just Na is not enought ...
          if ( n_basecation > 5 ) then
            ! add to basecation tracers:
            do itr = 1, n_basecation
              ispec = ispecs_basecation(itr)
              select case ( ispec )
                case(ispec_Ca_f); emis(i,j,1,ispec) = emis(i,j,1,ispec) + emis_f * bcat_frac(i,j,ii_Ca)
                case(ispec_Ca_c); emis(i,j,1,ispec) = emis(i,j,1,ispec) + emis_c * bcat_frac(i,j,ii_Ca)
                case(ispec_Na_f); emis(i,j,1,ispec) = emis(i,j,1,ispec) + emis_f * bcat_frac(i,j,ii_Na)
                case(ispec_Na_ff); emis(i,j,1,ispec) =  emis(i,j,1,ispec) !nothing yet
                case(ispec_Na_c); emis(i,j,1,ispec) = emis(i,j,1,ispec) + emis_c * bcat_frac(i,j,ii_Na)
                case(ispec_Na_cc); emis(i,j,1,ispec) = emis(i,j,1,ispec) !nothing yet
                case(ispec_Na_ccc); emis(i,j,1,ispec) = emis(i,j,1,ispec) !nothing yet
                case(ispec_Mg_f); emis(i,j,1,ispec) = emis(i,j,1,ispec) + emis_f * bcat_frac(i,j,ii_Mg)
                case(ispec_Mg_c); emis(i,j,1,ispec) = emis(i,j,1,ispec) + emis_c * bcat_frac(i,j,ii_Mg)
                case(ispec_K_f ); emis(i,j,1,ispec) = emis(i,j,1,ispec) + emis_f * bcat_frac(i,j,ii_K )
                case(ispec_K_c ); emis(i,j,1,ispec) = emis(i,j,1,ispec) + emis_c * bcat_frac(i,j,ii_K )
                case default
                  write (gol,'("unsupported tracer index ",i6)') ispec; call goErr
                  TRACEBACK; status=1; return
              end select
            end do
          end if

          ! if(runF%do_hm) etc

        end do  ! traffic maps

      end do   ! i
    end do   ! j


    !!open dat-file for testing:
    !open(u_tmp, file=trim(locF%output)//'test_swv.dat', status='replace', &
    !              access='direct', recl=nx*ny*4, iostat=status)
    !    if (status/=0) then
    !      write (*,'("ERROR(1) - opening file:")')
    !      write (*,'("ERROR(1) -   file    : ",a)') trim(locF%output)//'test_swv.dat'
    !      write (*,'("ERROR(1) -   iostat  : ",i6)') status
    !      stop 'ERROR(1) in print_swv'
    !    end if
    !    
    !    write (u_tmp,rec=1) ((swg_avg(i,j),i=1,nx),j=1,ny)
    !    write (u_tmp,rec=2) ((clim_fac(i,j),i=1,nx),j=1,ny)
    !    write (u_tmp,rec=3) ((cnordic(i,j),i=1,nx),j=1,ny)
    !close(u_tmp)
    !
    !stop

    !OK
    status = 0

  end subroutine roadresp_dust_em


  ! ***
  

  subroutine correction_swg( swg, status )
  
    use Dims      , only : nx, ny
    use LE_Landuse_Data, only : lu_fracs
    use LE_Landuse_Data, only : ilu_water_sea

    ! --- in/out ---------------------------------
    
    real, intent(inout)   ::  swg(:,:)
    integer, intent(out)  ::  status
    
    ! --- local ----------------------------------
     
    ! --- begin ----------------------------------
     
    ! correction for coastal low swv values due to mapping on LE grid via landuse class 
    ! (as a result of the use of fraction of land..), 
    ! necessary for the application for proper use of soil water content information for 
    ! regional factor (climate correction), and soil moisture inhibition as applied in Fecan 
    ! formulation.
    !   soil water corrected = soil water / (1 - landuse_class_fraction_sea)  

    where ( lu_fracs(:,:,ilu_water_sea) < 0.98 )
      swg = swg / ( 1.0 - lu_fracs(:,:,ilu_water_sea) )
    elsewhere
      swg = 1.0 ! set to very high soilwater content in case 
                ! some KM are given in these "wet" cells.
    endwhere

    ! ok
    status = 0

  end subroutine correction_swg


end module LE_Emis_Dust_Resuspension
