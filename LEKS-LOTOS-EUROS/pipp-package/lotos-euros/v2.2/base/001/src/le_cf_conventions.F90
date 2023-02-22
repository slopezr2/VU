!#################################################################
!
! NAME
!
!   LE_CF_Conventions
!
! DISCRIPTION
!
!   Tools to form standard names for use in output.
!
! USAGE
!
!   ! storage for cf names; chosen length is an indication ...
!   character(len=64)     ::  cf_standard_name
!   character(len=64)     ::  cf_long_name
!   character(len=16)     ::  cf_units
!   real                  ::  cf_units_conversion_factor
!   character(len=1024)   ::  comment
!   integer               ::  status
!
!   ! get names following cf conventions given tracer name and unit:
!   call LE_CF_names( &
!                     'O3', 'ppb', &
!                     cf_standard_name, cf_long_name, cf_units, &
!                     cf_units_conversion_factor, comments, &
!                     status )
!
!   ! NOTES:
!   !  o many different tracer names and units are accepted ...
!   !  o use the unit convesion factor to obtain data in
!   !    in CF units:
!   !      data_in_cf_units = data * cf_units_conversion_factor
!   !  o comments (in/out) is extended if necessary,
!   !    new comments are seperated using the ';' character
!
!
! SEE ALSO
!
!   http://wiki.esipfed.org/index.php/Air_Quality/Chemistry_Naming_Conventions

! HISTORY
!
!   2007 may, Arjo Segers, TNO
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

module LE_CF_Conventions


  use GO, only : gol, goPr, goErr

  !use NetCDF, only : NF90_StrError, NF90_NOERR

  implicit none


  ! --- in/out -----------------------------

  private

  public  ::  LE_CF_names


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_CF_Conventions'


contains


  ! ====================================================


  subroutine LE_CF_names( name, unit, &
                     cf_standard_name, cf_long_name, cf_units, &
                     cf_units_conversion_factor, comments, &
                     status, &
                     halflevel, &
                     cf_tracer_name )

    ! --- in/out --------------------------------

    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  unit
    character(len=*), intent(out)         ::  cf_standard_name
    character(len=*), intent(out)         ::  cf_long_name
    character(len=*), intent(out)         ::  cf_units
    real, intent(out)                     ::  cf_units_conversion_factor
    character(len=*), intent(inout)       ::  comments
    integer, intent(out)                  ::  status
    logical, intent(out), optional        ::  halflevel
    character(len=*), intent(out), optional    ::  cf_tracer_name

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/LE_CF_names'

    ! --- local ---------------------------------

    character(len=32)     ::  cf_type
    character(len=64)     ::  cf_medium_stnd, cf_medium_long
    character(len=64)     ::  cf_enti_stnd, cf_enti_long
    character(len=64)     ::  cf_spec_stnd, cf_spec_long
    !character(len=256)    ::  cf_name_stnd, cf_name_long, cf_name_unit
    character(len=512)    ::  comment

    ! --- begin ---------------------------------

    ! defaults:
    cf_standard_name           = ''
    cf_long_name               = 'unknown cf long name'
    cf_units                   = 'unknown cf units'
    cf_units_conversion_factor = 1.0
    if ( present(halflevel) ) halflevel = .false.

    ! no comments by default ...
    comment = ''

    ! standard and long names for a tracer according to CF conventions:
    select case ( name )

      case ( 'latitude', 'lat' )
        cf_type          = 'default'
        cf_standard_name = 'latitude'
        cf_long_name     = 'latitude'
        cf_units         = 'degrees_north'
      case ( 'longitude', 'lon' )
        cf_type          = 'default'
        cf_standard_name = 'longitude'
        cf_long_name     = 'longitude'
        cf_units         = 'degrees_south'
      case ( 'h' )
        cf_type          = 'default'
        cf_standard_name = 'height'
        cf_long_name     = 'layer height relative to orography'
        cf_units         = 'm'
      case ( 'oro' )
        cf_type          = 'default'
        cf_standard_name = 'surface_altitude'
        cf_long_name     = 'orography above sea level'
        cf_units         = 'm'
      case ( 'blh' )
        cf_type          = 'default'
        cf_standard_name = 'atmosphere_boundary_layer_thickness'
        cf_long_name     = 'boundary layer height'
        cf_units         = 'm'
      case ( 'tcc' )
        cf_type          = 'default'
        cf_standard_name = 'cloud_area_fraction'
        cf_long_name     = 'total cloud cover'
        cf_units         = '0-1'
      case ( 'rad' )
        cf_type          = 'default'
        cf_standard_name = ''
        cf_long_name     = 'surface solar radiation'
        cf_units         = 'W m-2'
      case ( 'radd' )
        cf_type          = 'default'
        cf_standard_name = ''
        cf_long_name     = 'surface solar radiation downward'
        cf_units         = 'W m-2'
      case ( 'tsurf' )
        cf_type          = 'default'
        cf_standard_name = 'surface_temperature'
        cf_long_name     = 'surface temperature'
        cf_units         = 'K'
      case ( 'slpress' )
        cf_type          = 'default'
        cf_standard_name = 'sealevel_pressure'
        cf_long_name     = 'air pressure at sea level'
        cf_units         = 'Pa'
      case ( 'rain' )
        cf_type          = 'default'
        cf_standard_name = 'lwe_precipitation_rate'
        cf_long_name     = 'liquid water equivalent precipitation rate'
        cf_units         = 'm s-1'
      case ( 'dens' )
        cf_type          = 'default'
        cf_standard_name = 'air_density'
        cf_long_name     = 'air density'
        cf_units         = 'kg m-3'
      case ( 'rhumid', 'srh' )
        cf_type          = 'default'
        cf_standard_name = 'relative_humidity'
        cf_long_name     = 'relative humidity'
        cf_units         = '%'
      case ( 'shumid' )
        cf_type          = 'default'
        cf_standard_name = 'specific_humidity'
        cf_long_name     = 'specific humidity'
        cf_units         = 'kg kg-1'
      case ( 'temper' )
        cf_type          = 'default'
        cf_standard_name = 'air_temperature'
        cf_long_name     = 'air temperature '
        cf_units         = 'K'
      case ( 'u' )
        cf_type          = 'default'
        cf_standard_name = 'eastward_wind'
        cf_long_name     = 'eastward wind '
        cf_units         = 'm s-1'
      case ( 'v' )
        cf_type          = 'default'
        cf_standard_name = 'northward_wind'
        cf_long_name     = 'northward wind '
        cf_units         = 'm s-1'
      case ( 'ustar' )
        cf_type          = 'default'
        cf_standard_name = 'friction_velocity'
        cf_long_name     = 'friction velocity'
        cf_units         = 'm s-1'
      case ( 'kz' )
        cf_type          = 'default'
        cf_standard_name = ''
        cf_long_name     = 'vertical diffusion velocity'
        cf_units         = 'm2 s-1'
      case ( 'wsurf' )
        cf_type          = 'default'
        cf_standard_name = '10m_wind_speed'
        cf_long_name     = '10m wind speed'
        cf_units         = 'm s-1'
      case ( 'psurf' )
        cf_type          = 'default'
        cf_standard_name = 'surface_air_pressure'
        cf_long_name     = 'surface pressure'
        cf_units         = 'Pa'
      case ( 'pressure' )
        cf_type          = 'default'
        cf_standard_name = 'air_pressure'
        cf_long_name     = 'pressure'
        cf_units         = 'Pa'
      case ( 'swvl1' )
        cf_type          = 'default'
        cf_standard_name = 'volumetric_soil_water_content_layer_1'
        cf_long_name     = 'volumetric soil water content layer 1'
        cf_units         = 'm3 m-3'
      case ( 'swvl2' )
        cf_type          = 'default'
        cf_standard_name = 'volumetric_soil_water_content_layer_2'
        cf_long_name     = 'volumetric soil water content layer 2'
        cf_units         = 'm3 m-3'
      case ( 'swvl3' )
        cf_type          = 'default'
        cf_standard_name = 'volumetric_soil_water_content_layer_3'
        cf_long_name     = 'volumetric soil water content layer 3'
        cf_units         = 'm3 m-3'
      case ( 'swvl4' )
        cf_type          = 'default'
        cf_standard_name = 'volumetric_soil_water_content_layer_4'
        cf_long_name     = 'volumetric soil water content layer 4'
        cf_units         = 'm3 m-3'
      case ( 'smi2' )
        cf_type          = 'default'
        cf_standard_name = 'soil_moisture_index_upper_2_layers'
        cf_long_name     = 'soil moisture index upper 2 layers'
        cf_units         = '-'
      case ( 'smi3' )
        cf_type          = 'default'
        cf_standard_name = 'soil_moisture_index_upper_3_layers'
        cf_long_name     = 'soil moisture index upper 3 layers'
        cf_units         = '-'
      case ( 'area' )
        cf_type          = 'default'
        cf_standard_name = 'cell_area'
        cf_long_name     = 'area '
        cf_units         = 'm2'
      case ( 'halflevel_altitude' )
        cf_type          = 'default'
        cf_standard_name = ''      ! not available yet ...
        cf_long_name     = 'half level altitude above sea level'      ! <-- guess ...
        cf_units         = 'm'
        if (present(halflevel)) halflevel = .true.
      case ( 'tau' )
        cf_type          = 'default'
        cf_standard_name = 'atmosphere_optical_thickness_per_layer'
        cf_long_name     = 'atmosphere_optical_thickness_due_to_ambient_aerosol_per_layer'
        cf_units         = '-'
      case ( 'extinction' )
        cf_type          = 'default'
        cf_standard_name = 'atmosphere_extinction_coefficient'
        cf_long_name     = 'atmosphere_extinction_coefficient_due_to_ambient_aerosol'
        cf_units         = '-'
      case ( 'ssa' )
        cf_type          = 'default'
        cf_standard_name = 'single_scattering_albedo'
        cf_long_name     = 'atmosphere_single_scattering_albedo_due_to_ambient_aerosol'
        cf_units         = '-'
      case ( 'asy' )
        cf_type          = 'default'
        cf_standard_name = 'atmosphere_asymmetry_parameter'
        cf_long_name     = 'atmosphere_asymmetry_parameter_due_to_ambient_aerosol'
        cf_units         = '-'
      case ( 'AOD2' )
        cf_type          = 'default'
        cf_standard_name = 'atmosphere_optical_thickness'
        cf_long_name     = 'atmosphere_optical_thickness_due_to_ambient_aerosol'
        cf_units         = '-'
      case ( 'CO', 'co'   )
        cf_type      = 'tracer'
        cf_spec_stnd = 'carbon_monoxide'
        cf_spec_long = 'CO'
        cf_units     = 'mole mole-1'
      case ( 'CO2', 'co2'   )
        cf_type      = 'tracer'
        cf_spec_stnd = 'carbon_dioxide'
        cf_spec_long = 'CO2'
        cf_units     = 'mole mole-1'
      case ( 'O3', 'o3', 'o3_biascorr' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'ozone'
        cf_spec_long = 'O3'
        cf_units     = 'mole mole-1'
      case ( 'NO', 'no'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'nitrogen_monoxide'
        cf_spec_long = 'NO'
        cf_units     = 'mole mole-1'
      case ( 'NO2', 'no2', 'no2_obs'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'nitrogen_dioxide'
        cf_spec_long = 'NO2'
        cf_units     = 'mole mole-1'
      case ( 'N2O5', 'n2o5'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'dinitrogen_pentoxide'
        cf_spec_long = 'N2O5'
        cf_units     = 'mole mole-1'
      case ( 'NOy', 'noy'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''             ! not available yet ...
        cf_spec_long = 'NOy'
        comment      = 'NOy = NOx + HNO3 + PAN + org.ntr., '// &
                         'with NOx = NO + NO2 + NO3 + HNO4 + N2O5 ;'
        cf_units     = 'mole mole-1'
      case ( 'CH2O', 'ch2o', 'CHOH', 'choh', 'FORM', 'form' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'formaldehyde'
        cf_spec_long = 'CH2O'
        cf_units     = 'mole mole-1'
      case ( 'SO2', 'so2'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'sulfur_dioxide'
        cf_spec_long = 'SO2'
        cf_units     = 'mole mole-1'
      case ( 'CH4', 'ch4'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'methane'
        cf_spec_long = 'CH4'
        cf_units     = 'mole mole-1'
      case ( 'OH', 'oh'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'hydroxyl_radical'
        cf_spec_long = 'OH'
        cf_units     = 'mole mole-1'
      case ( 'H2O2', 'h2o2'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'hydrogen_peroxide'
        cf_spec_long = 'H2O2'
        cf_units     = 'mole mole-1'
      case ( 'NH3', 'nh3' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'ammonia'
        cf_spec_long = 'NH3'
        cf_units     = 'mole mole-1'
      case ( 'HNO3', 'hno3'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'nitric_acid'
        cf_spec_long = 'HNO3'
        cf_units     = 'mole mole-1'
      case ( 'PAN', 'pan'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'peroxyacetyl_nitrate'
        cf_spec_long = 'PAN'
        cf_units     = 'mole mole-1'
      case ( 'OLE', 'ole'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''             ! not available yet ...
        cf_spec_long = 'olefin'
        cf_units     = 'mole mole-1'
      case ( 'ALD', 'ald'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'aldehyde'
        cf_units     = 'mole mole-1'
      case ( 'MGLY', 'mgly'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'methylglyoxal'
        cf_units     = 'mole mole-1'
      case ( 'XYL', 'xyl'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'xylene'
        cf_units     = 'mole mole-1'
      case ( 'TOL', 'tol'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'toluene'
        cf_spec_long = 'toluene'
        cf_units     = 'mole mole-1'
      case ( 'PAR', 'par'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'parrafin bond'
        cf_units     = 'mole mole-1'
      case ( 'ISO', 'iso'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'isoprene'
        cf_units     = 'mole mole-1'
      case ( 'TERP', 'terp'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'terpenes'
        cf_units     = 'mole mole-1'
      case ( 'ETH', 'eth'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'ethane'
        cf_units     = 'mole mole-1'
      case ( 'NO3', 'no3'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'nitrate_radical'
        cf_spec_long = 'nitrate radical'
        cf_units     = 'mole mole-1'
      case ( 'Rn', 'rn', 'Radon', 'radon'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'radon'
        cf_spec_long = 'Rn'
        cf_units     = 'mole mole-1'
      case ( 'SF6', 'sf6'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'sulfur hexafluoride'
        cf_spec_long = 'SF6'
        cf_units     = 'mole mole-1'
      case ( 'Pb', 'pb', 'Lead', 'lead'  )
        cf_type      = 'tracer'
        cf_spec_stnd = ''    ! not available yet ...
        cf_spec_long = 'Pb'
        cf_units     = 'mole mole-1'
      case ( 'SO4', 'so4', 'SO4a', 'so4a', 'so4a_f', 'so4a_c' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'sulfate_dry_aerosol'
        cf_spec_long = 'SO4'
        cf_units     = 'kg m-3'
      case ( 'NO3_A', 'no3_a', 'NO3A', 'no3a', 'no3a_f', 'no3a_c'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'nitrate_dry_aerosol'
        cf_spec_long = 'NO3'
        cf_units     = 'kg m-3'
      case ( 'ec_f'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'black_carbon_dry_aerosol'
        cf_spec_long = 'ec_f'
        cf_units     = 'kg m-3'
      case ( 'ec_c'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'black_carbon_dry_aerosol'
        cf_spec_long = 'ec_c'
        cf_units     = 'kg m-3'
      case ( 'pom_f'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'primary_organic_matter_dry_aerosol'
        cf_spec_long = 'primary organic matter dry aerosol'
        cf_units     = 'kg m-3'
      case ( 'pom_c'  )
        cf_type      = 'tracer'
        cf_spec_stnd = 'primary_organic_matter_dry_aerosol'
        cf_spec_long = 'primary organic matter dry aerosol'
        cf_units     = 'kg m-3'
      case ( 'nh4a_f', 'nh4a_c' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'ammonium_dry_aerosol'
        cf_spec_long = 'ammonium dry aerosol'
        cf_units     = 'kg m-3'
      case ( 'ppm_f', 'ppm_c' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'primary_particulate_matter_dry_aerosol'
        cf_spec_long = 'PPM2.5'
        cf_units     = 'kg m-3'
      case ( 'na_f', 'na_c', 'na_ff','na_ccc','na_cc' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'sodium'
        cf_units     = 'kg m-3'
      case ( 'dust_f', 'dust_c', 'dust_ff', 'dust_ccc','dust_cc' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'dust_dry_aerosol'
        cf_spec_long = 'dust dry aerosol'
        cf_units     = 'kg m-3'
      case ( 'vbs_soa1', 'vbs_soa2', 'vbs_soa3', &
             'vbs_soa4', 'vbs_soa5', 'vbs_soa6', &
             'vbs_soa7', 'vbs_soa8', 'vbs_soa9' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'soa_aerosol'
        cf_spec_long = 'soa aerosol'
        cf_units     = 'kg m-3'
      case ( 'vbs_cg1', 'vbs_cg2', 'vbs_cg3', &
             'vbs_cg4', 'vbs_cg5', 'vbs_cg6', &
             'vbs_cg7', 'vbs_cg8', 'vbs_cg9' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'condensable gas'
        cf_spec_long = 'condensable gas'
        cf_units     = 'mole mole-1'
      case ( 'tpm10', 'tpm10_biascorr', 'tpm10_comb' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TPM10'
        cf_units     = 'kg m-3'
      case ( 'tpm25', 'tpm25_biascorr', 'tpm25_comb' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TPM25'
        cf_units     = 'kg m-3'
      case ( 'TSS', 'tss' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TSS'
        cf_units     = 'kg m-3'
      case ( 'TNMVOC', 'tnmvoc' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TNMVOC'
        cf_units     = 'kg m-3'
      case ( 'TSOA', 'tsoa' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TSOA'
        cf_units     = 'kg m-3'
      case ( 'TPOA', 'tpoa' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TPOA'
        cf_units     = 'kg m-3'
      case ( 'TSISOA', 'tsisoa' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TSISOA'
        cf_units     = 'kg m-3'
      case ( 'TASOA', 'tasoa' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TASOA'
        cf_units     = 'kg m-3'
      case ( 'TBSOA', 'tbsoa' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TBSOA'
        cf_units     = 'kg m-3'
      case ( 'TCG', 'tcg' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TCG'
        cf_units     = 'kg m-3'
      case ( 'TPOG', 'tpog' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TPOG'
        cf_units     = 'kg m-3'
      case ( 'TSISOG', 'tsisog' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TSISOG'
        cf_units     = 'kg m-3'
      case ( 'TASOG', 'tasog' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TASOG'
        cf_units     = 'kg m-3'
      case ( 'TBSOG', 'tbsog' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''         ! not available yet ...
        cf_spec_long = 'TBSOG'
        cf_units     = 'kg m-3'
      case ( 'TDUST', 'tdust' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'dust_dry_aerosol'
        cf_spec_long = 'TDUST'
        cf_units     = 'kg m-3'
      case ( 'SIA25', 'sia25' )
        cf_type      = 'tracer'
        cf_spec_stnd = 'secondary_inorganic_aerosol <2.5um'
        cf_spec_long = 'SIA25'
        cf_units     = 'kg m-3'

      case ( 'drysox' )
        cf_type      = 'integrated-tracer'
        cf_standard_name = ''
        cf_long_name = 'drysox'
        cf_units     = 'ug m-2'
      case ('wetsox' )
        cf_type      = 'integrated-tracer'
        cf_standard_name = 'wetsox'
        cf_long_name = 'wetsox'
        cf_units     = 'ug m-2'
      case ( 'drynoy' )
        cf_type      = 'integrated-tracer'
        cf_standard_name = ''
        cf_long_name = 'drynoy'
        cf_units     = 'ug m-2'
      case ( 'wetnoy' )
        cf_type      = 'integrated-tracer'
        cf_standard_name = ''
        cf_long_name = 'wetnoy'
        cf_units     = 'ug m-2'
      case ( 'wetnhx' )
        cf_type      = 'integrated-tracer'
        cf_standard_name = ''
        cf_long_name = 'wetnhx'
        cf_units     = 'ug m-2'
      case ( 'drynhx' )
        cf_type      = 'integrated-tracer'
        cf_standard_name = ''
        cf_long_name = 'drynhx'
        cf_units     = 'ug m-2'
      case ( 'AOD', 'aod' )
        cf_type          = 'integrated-tracer'
        cf_standard_name = 'atmosphere_optical_thickness_due_to_aerosol'
        cf_long_name     = 'atmosphere optical thickness due to aerosol'
        cf_units         = '1'
      case ( 'aod_biascorr' )
        cf_type          = 'integrated-tracer'
        cf_standard_name = ''
        cf_long_name     = 'atmosphere optical thickness due to aerosol (bias corrected)'
        cf_units         = '1'
      case ( 'aerh2o' )
        cf_type          = 'tracer'
        cf_standard_name = ''
        cf_long_name     = 'aerosol water concentration'
        cf_units         = 'kg m-3'
      case( 'nh3ave' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''
        cf_spec_long = 'nh3ave'
        cf_units     = 'mole mole-1'
      case( 'so2ave' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''
        cf_spec_long = 'so2ave'
        cf_units     = 'mole mole-1'
      case( 'o3max' )
        cf_type      = 'tracer'
        cf_spec_stnd = ''
        cf_spec_long = 'o3max'
        cf_units     = 'mole mole-1'
      case default

        !write (gol,'("unsupported tracer name : ",a)') trim(name); call goErr
        !TRACEBACK; status=1; return

        ! set dummy values:
        cf_type          = 'default'

        !cf_units         = trim(unit)
        !cf_units_conversion_factor = 1.0
        !write (gol,'("LE_CF_names - WARNING - unsupported tracer name for CF conventions : ",a)') trim(name); call goPr
        ! continue with units and conversion factor ...
        !status=0; return

    end select

    ! further process given type ...
    select case ( cf_type )

      case ( 'default' )

        ! unit conversion factors:
        if ( trim(unit) == trim(cf_units) ) then
          cf_units_conversion_factor = 1.0   ! no conversion
        else if ( unit == '<cf>' ) then
          cf_units_conversion_factor = 1.0   ! assume unit is the prefered cf unit ..
        else
          select case ( unit )
            case ( '%' )
              select case ( cf_units )
                case ( '0-1', 'kg kg-1', 'kg kg**-1', 'kg/kg' )
                  cf_units_conversion_factor = 0.01  ! [0,100] -> [0,1]
                case default
                  write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(cf_units); call goErr
                  TRACEBACK; status=1; return
              end select
            case ( 'km' )
              select case ( cf_units )
                case ( 'm' ) ; cf_units_conversion_factor = 1.0e3   ! km -> m
                case default
                  write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(cf_units); call goErr
                  TRACEBACK; status=1; return
              end select
            case ( 'km2', 'km**2' )
              select case ( cf_units )
                case ( 'm2' ) ; cf_units_conversion_factor = 1.0e6   ! km2 -> m2
                case default
                  write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(cf_units); call goErr
                  TRACEBACK; status=1; return
              end select
            case ( 'km min-1', 'km min**-1', 'km/min' )
              select case ( cf_units )
                case ( 'm s-1' ) ; cf_units_conversion_factor = 1.0e3/60.0   ! km/min -> m/s
                case default
                  write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(cf_units); call goErr
                  TRACEBACK; status=1; return
              end select
            case ( 'm 3hr-1' )
              select case ( cf_units )
                case ( 'mm hr-1' ) ; cf_units_conversion_factor = 1.0e3/3.0   ! m/3hr -> mm/hr
                case default
                  write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(cf_units); call goErr
                  TRACEBACK; status=1; return
              end select
            case default
              !write (gol,'("unsupported default unit : ",a)') trim(unit); call goErr
              !TRACEBACK; status=1; return
              ! just continue ...
              cf_units = trim(unit)
              cf_units_conversion_factor = 1.0   ! no conversion
          end select
        end if

      case ( 'tracer' )

        if ( present(cf_tracer_name) ) cf_tracer_name = trim(cf_spec_stnd)

        select case ( cf_units )
          case ( 'mole mole-1' )
            cf_enti_stnd = 'mole_fraction'
            cf_enti_long = 'volume mixing ratio'
            select case ( unit )
              case ( '<cf>'          ) ; cf_units_conversion_factor = 1.0      ! assume unit is the prefered cf unit ..
              case ( 'mole mole-1', 'mole mole**-1', 'mole/mole' )
                cf_units_conversion_factor = 1.0      ! already cf ...
              case ( 'ppm'           ) ; cf_units_conversion_factor = 1.0e-6   ! ppm -> mole/mole
              case ( 'ppb'           ) ; cf_units_conversion_factor = 1.0e-9   ! ppb -> mole/mole
              case ( 'ppt'           ) ; cf_units_conversion_factor = 1.0e-12  ! ppt -> mole/mole
              case default
                write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(cf_units); call goErr
                TRACEBACK; status=1; return
            end select
          case ( 'kg m-3' )
            cf_enti_stnd = 'mass_concentration'
            cf_enti_long = 'mass concentration'
            select case ( unit )
              case ( '<cf>'              ) ; cf_units_conversion_factor = 1.0      ! assume unit is the prefered cf unit ..
              case ( 'kg m**-3', 'kg/m3' ) ; cf_units_conversion_factor = 1.0      ! already cf ...
              case ( 'ug/m3'             ) ; cf_units_conversion_factor = 1.0e-9   ! ug/m3 -> kg/m3
              case default
                write (gol,'("do not know how to convert `",a,"` from `",a,"` to `",a,"`")') &
                             trim(name), trim(unit), trim(cf_units); call goErr
                TRACEBACK; status=1; return
            end select
          case ('#/m3' )
            cf_enti_stnd = 'number_concentration'
            cf_enti_long = 'number concentration'
            select case( unit )
              case ( '<cf>'              ) ; cf_units_conversion_factor = 1.0      ! assume unit is the prefered cf unit ..
              case ( 'grns/m3'           ) ; cf_units_conversion_factor = 1.0      ! grains conc. is the number conc.
              case default
                write (gol,'("do not know how to convert `",a,"` from `",a,"` to `",a,"`")') &
                             trim(name), trim(unit), trim(cf_units); call goErr
                TRACEBACK; status=1; return
            end select
          case default
            write (gol,'("unsupported cf unit : ",a)') trim(cf_units); call goErr
            write (gol,'("for name,unit : ",a,", ",a)') trim(name), trim(unit); call goErr
            TRACEBACK; status=1; return
        end select

        ! CF standard name for medium:
        cf_medium_stnd = 'in_air'
        cf_medium_long = 'in humid air'

        ! total names:
        if ( len_trim(cf_spec_stnd) > 0 ) then
          write (cf_standard_name,'(a,"_of_",a,"_",a)') trim(cf_enti_stnd), trim(cf_spec_stnd), trim(cf_medium_stnd)
        else
          cf_standard_name = ''
        end if
        write (cf_long_name    ,'(a," of ",a," ",a)') trim(cf_enti_long), trim(cf_spec_long), trim(cf_medium_long)

      case ( 'integrated-tracer' )

        ! unit conversion factors:
        if ( (trim(unit) == trim(cf_units)) .or. (trim(unit) == '<cf>') ) then
          cf_units_conversion_factor = 1.0   ! no conversion
        else
          ! unit and conversion factors:
          !select case ( unit )
          !  case ( 'm' )
          !    select case ( cf_units )
          !      case ( 'm' ) ; cf_units_conversion_factor = 1.0   ! m -> m
          !      case default
          !        write (gol,'("do not know how to convert `",a,"` to `",a,"`")') trim(unit), trim(cf_units); call goErr
          !        TRACEBACK; status=1; return
          !    end select
          !  case default
              write (gol,'("unsupported ",a," unit : ",a)') trim(cf_type), trim(unit); call goErr
              TRACEBACK; status=1; return
          !end select
        end if

      case default

        !write (gol,'("unsupported cf type : ",a)') trim(cf_type); call goErr
        !TRACEBACK; status=1; return

    end select

    ! add comments if necessary:
    if ( len_trim(comment) > 0 ) then
      if ( len_trim(comments) > 0 ) then
        comments = trim(comments)//' ; '//trim(comment)
      else
        comments = trim(comment)
      end if
    end if

    ! ok
    status = 0

  end subroutine LE_CF_names


end module LE_CF_Conventions

