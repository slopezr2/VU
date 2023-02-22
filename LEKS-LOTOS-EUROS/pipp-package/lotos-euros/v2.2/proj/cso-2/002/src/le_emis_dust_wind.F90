!###############################################################################
!
! LE_Emis - LOTOS-EUROS wind blown dust emission routines
!
! Reports:
!
!   For a full description of dust emisisons as used in LOTOS-EUROS,
!   see the report of the BOP project:
!
!   [Schaap 2009]
!     Schaap M; Manders AMM; Hendriks ECJ; Cnossen JM; Segers AJS; Denier van der Gon HAC;
!     Jozwicka M; Sauter F; Velders G; Matthijsen J; Builtjes PJH
!     "Regional modelling of particulate matter for the Netherlands"
!     PBL Report no.    500099008
!     http://www.pbl.nl/bibliotheek/rapporten/500099008.pdf
!
!   or the summary report:
!
!   [Schaap et al., 2010]
!    "Mineral Dust as a component of particulate matter"
!    http://www.pbl.nl/en/publications/2010/Mineral-Dust-component-particulate-matter
!
! References:
!
!   [Alfaro&Gomes 2001]
!     St\'ephae C. Alfaro, Laurent Gomes
!     Modeling mineral aerosol production by wind erosion:
!       Emission intensities and aerosol size distribution in source areas
!     JGR 106 D16, p 18.075-18.084
!
!   [Marticorena and Bergametti, 1995]
!     Modeling the atmoshperic dust cycle: 1. Design of a soil-derived dust emission scheme
!     JGR, vol. 100, no. D8, p. 16,415-16,430 .
!
!
! Notes by Elise Hendriks, 2011-08-30
! ===================================
!
! WINDBLOWN DUST contributions to PM10 calculations for BOP project
! (periods 2005, 2007-2008)
!
! Basic approach
! --------------
!
! Windblown dust emission after Marticorena and Bergametti (1995) and
! applying soil texture dependence  of Afaro et al.(1997, 1998) with:
!
! * Soiltexture settings on basis of JRC soil texture and Chatenet etal
! (1996) classification (see BOP  technical report).
!
! * Soil moisture (and other meteo, except u* which is inferred from
! surface wind velocity and  meteorological stability parameters in
! stability.F90 subroutine) from ECMWF.
!
! * Roughness length and smooth roughness length definitions (see below).
!
! * Threshold friction velocity definitions (lower boundary condition of
! local friction velocity for uplift  of sand particles).
!
!
! Settings and switches
! ----------------------
!
! Roughness length:
! * Standard applied aerodynamic roughness length are the DEPAC defined
! and landuse dependent values  (z0_depac(ilu) in the code).  For this
! purpose though, a different approach is taken to calculate the
! windblowndust production:  It is assumed there is only uplift over
! arable land (ilu=2, depac) and the uplift is vegetation cover
! dependent (monthly variation is assumed, bare land in winter, covered
! land in summer).
! * Vegetation cover dependent roughness length (z0_1 in the code)
!   ~ (OLD) defined as 1e-5 x z0_factor(mm),
!     where z0_factor(mm) is the vegetation factor,
!       z0_factor(12) = (/1.0,1.0,1.0,1.0e4,2.5e4,2.5e4,2.5e4,2.5e4,1.0e4, 1.0, 1.0, 1.0/)  and the value 1e-5m is
!     chosen such that z0=z0s, hence f_eff=1 in case of no vegetation. In
!     other  words, largest possible effective uplift is assumed when there
!     is no vegetation cover!
!   ~ New implementation in 2001-08:
!     monthly z0 values with more or less realistic values.
! * z0s is the smooth roughness length, being set to 1e-5m (values in
!   literature are 4e-6 in Afaro et al,  (1997, 1998) 1e-5m in Marticorena
!   and Bergametti (1995)).
!
! Rain and snow conditions:
! * no uplift when snow, nor with rain>0.3mm.
! * a setup for no uplift of particles within 48 hours after rain is made
! (PrecipitationLast48hours flag),  but not used (value .false.).
!
! Threshold friction velocity:
! Reference and uncorrected threshold friction velocity is set to
! 0.25m/s. This value is hence, corrected for soilmoisture (after Fecan
! etal (1999), f_w), surface roughness and vegetation (f_eff,
! Marticorena  and Bergametti (1995) and surface roughness and
! vegetation assumption above).
!
! ==============================================================
!
! Reformulation by Vincent Kamphuis and Arjo Segers.
! September 2011 .
!
! Classification of aerosols released from surface following
! [Alfaro&Gomes,2001] into 3 modes:
!
!   Values and labels:
!     imode                          :   1     2       3
!     meaning                        : fine, coarse, extreme
!   Settings following Table 1 in [Alfaro&Gomes,2001]
!   or similar Table 5.4 in BOP report:
!     size distribution parameters   : Dp_a(imode), sigma_a(imode)
!     binding energy                 : e_a(imode)
!
! Dust is released from erodible soils only.
! The texture of an erodible soil is one of the parameters that determine
! the actual erosion. Here the following texture classes are distinguished:
!
!  itext = texture type
!     possible values   :  course, medium, medium fine, fine, very fine
!
! Emission in a cell is the sum of vertical fluxes launched from
! the different erodible soils within a cell:
!
!                 ntext
!  emw(imode) =    sum   texture_area(itext)  Fv(imode,itext)
!                itext=1
!
! Texture implies soil size distribution (particle sizes):
!
!  issd  = soil size mode
!     possible values                :  FFS, FS, MS, CS
!     size distribution parameters   :  Dp_s(issd), sigma_s(issd)
!
! Vertical flux from a certian soil is sum of contributions
! from the present soil size modes:
!
!                      nssd
!   Fv(imode,itext) =  sum   texture_frac(issd)  Fv_ssd(imode,itext,issd)
!                     issd=1
!
! Assume log-normal number distribution for the soil particles within
! a soil size mode. Distribution is defined by the median Dp_s(issd)
! and the width sigma_s(issd), see Table 5.2 in report.
!
!      dNp
!    ------- ~ Gaussian( Dp_s(issd), sigma_s(issd) )
!    dln(Dp)
!
! Vertical flux of aerosol tracers is the result of sand blasting,
! i.e. a horizontal flux of particles smashes aerosols from the soil
! into the air. Per soil particle diameter one can compute the horizontal
! flux Fh and the multiplication factor alpha that gives the resulting
! vertical flux. Integrate over all soil particles:
!
!                                                                                      dNp
!   Fv_ssd(imode,itext,issd,u*)  = int Fh(u*,u*t(Dp,z0,z0s,itext)) alpha(Dp,u*,imode) ----- dlnDp
!                                 lnDp                                                 dlnDp
!
! Horizontal flux:
!
!                          { 0               ,  u* <= u*t(Dp,..)
!   Fh(u*,u*t(Dp,z0,z0s) = {
!                          { rho_air/g (u*)**3 (1-u*t(Dp,..)/u*) (1+u*t(Dp,..)/u*)**2
!
! Threshold for friction velocity:
!
!   u*t(Dp,z0,z0s,itext) = u*tu(Dp) f_w(itext)  /  f_eff(z0,z0s)
!
! with the uncorrected friction velocity:
!
!   u*tu(Dp) = ... Dp, rho_soil, rho_air, ...   (formulas p41 of report)
!
! Threshold factor depending on water and clay content:
!
!   f_w(itext) = ..., clay_frac(itext), soilwater_frac, ...
!
! Threshold scaling depending on roughness:
!
!                            ln(z0/z0s)
!   f_eff(z0,z0s) = 1 - ( --------------- )
!                         ln(a(X/z0s)**p)
!
! Amount of vertical flux relative to horizontal flux:
!
!                                          p(imode,e_c(Dp,u*)) Dp_a(imode)**3
!    alpha(Dp,u*,imode) =  pi/6 rho_soil B ----------------------------------
!                                                       e_a(imode)
!
!                                                                   p(imode,e_c(Dp,u*))
!                       =  [4/3 pi (Dp_a(imode)/2)**3]  rho_soil B -------------------
!                                                                      e_a(imode)
!
! with energy :
!
!    e_c(Dp,u*) = 0.5 (pi/6 Dp**3) rho_soil (20 u*)**2
!
! and parameters p from table 5.5 :
!
!    p(imode,ec(Dp,u*)) = ...
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! all above is revised by S.L. Janson (July, 2014) > internship
! report found in folder Users/Sjoerd/Verslag inc. figuren/
! in SV folder (\\tsn.tno.nl\Data\SV)
! methods mostly based on article by Mokhtari et al. 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

Module LE_Emis_Dust_Wind

  use GO, only : gol, goPr, goErr
  use GO, only : T_LUT

  implicit none


  ! --- in/out -----------------------------------

  private

  public    ::  LE_Emis_Dust_Wind_Init, LE_Emis_Dust_Wind_Done
  public    ::  LE_Emis_Dust_Wind_Add
  
  public    :: ustarthreshold


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Dust_Wind'


  ! --- local ------------------------------------

  ! Lookup table for dust flux:
  type(T_LUT)          ::  lut
  ! Function of (u*,u*_thr,issd,imode) ;
  ! - friction velocities within [0,6] m/s
  !   5 m/s prooved to be too small somewhere in July 2008
  integer, parameter   ::  n_ustar = 120  
  real, parameter      ::  d_ustar = 0.05
  !! - friction velocity threshold factor in [0,30]
  !integer, parameter   ::  n_ustar_thrf = 120
  !real, parameter      ::  d_ustar_thrf = 0.25
  ! - friction velocity threshold factor in [0,30]
  integer, parameter   ::  n_ustar_thrf = 600
  real, parameter      ::  d_ustar_thrf = 0.05
  real, parameter      ::  max_ustar_thrf = n_ustar_thrf * d_ustar_thrf

  real, allocatable     ::  pot_source(:,:)

  integer, parameter    ::  nbins = 5       ! number of aerosol bins
  ! aerosol bins maximum  and minimum particle sizes
  real, parameter       ::  binsize_max(nbins) = (/ 1.0e-6,     2.5e-6,     4.0e-6,     7.0e-6,    10.0e-6 /)
  real, parameter       ::  binsize_min(nbins) = (/ 0.01e-6,    1.0e-6,     2.5e-6,     4.0e-6,     7.0e-6 /)

  
  real, allocatable     ::  soilwater_coast_correction(:,:)
  real, allocatable     ::  ustarthreshold(:,:)


contains


  ! ========================================================================

  subroutine dust_wind_potential_sources( rcF, status )
  
    use GO            , only : TrcFile, ReadRc
    use LE_Grid       , only : ugg

    use C3PO          , only : T_File_Ugg, T_Grid_Ugg
    use LE_Data_Common, only : Grid_Convertors
      
    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/Dust_Wind_Potential_Sources' 
 

    character(len=32), parameter  ::    description='var_name=level'
    
    ! --- local ---------------------------------

    character(len=512)    ::  fname
    integer               ::  varid
    
    type(T_File_Ugg)      ::  file_in
    type(T_Grid_Ugg)      ::  grid_in
    real, allocatable     ::  pot_sources_in(:,:)
    character(len=64)     ::  units_in
    
    
    ! --- begin ---------------------------------
  
    call ReadRc( rcF, 'landuse.pot_sources.file', fname, status )
    IF_NOTOK_RETURN(status=1)
    
    
    ! check if potential sources is required
    if ( trim(fname) /= '' ) then 
        
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
      allocate( pot_sources_in(grid_in%nlon,grid_in%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
    
      ! read:
      call file_in%Get_Var( trim(description), pot_sources_in, units_in, status, &
                             start=(/1,1/), count=(/grid_in%nlon,grid_in%nlat/) )
      IF_NOTOK_RETURN(status=1)                     
      
      ! allocate on target grid        
      allocate( pot_source(ugg%nlon,ugg%nlat), stat=status )
      IF_NOTOK_RETURN(status=1)
      
      ! regrid:
      call Grid_Convertors%Ugg_AreaAver( grid_in, pot_sources_in, ugg, pot_source, status )
      IF_NOTOK_RETURN(status=1)
    
      ! clear
      deallocate( pot_sources_in )
      call grid_in%Done( status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      call file_in%Close( status )
      IF_NOTOK_RETURN(status=1)
    
    else
        ! no special potential sources file defined
        allocate( pot_source(ugg%nlon,ugg%nlat), stat = status )
        IF_NOTOK_RETURN(status=1)
        
        pot_source = 1.0
    end if
    
    ! ok
    status = 0 
  
  end subroutine dust_wind_potential_sources
  
  
  ! ***
  

  subroutine LE_Emis_Dust_Wind_Init( rcF, status )

    use GO           , only : goGetFU, pathsep
    use GO           , only : LUT_Init, LUT_Def_Axis, LUT_Coordinate, LUT_Set_Value
    use LE_Comm      , only : goc_dom
    use Binas        , only : rho_dry_air_20C_1013hPa
    use JAQL_SFC_Soil, only : nssd, ssd_gmean, ssd_gsigma, ssd_mfrac
    use JAQL_SFC_Soil, only : rho_mineral_soil
    use JAQL_SFC_Dust, only : nmode
    use JAQL_SFC_Dust, only : Friction_Velocity_Threshold_Factor_Soilwater_Demo
    use GO           , only : TrcFile, ReadRc
    use Dims         , only : nx, ny
    use Binas        , only : pi
    use LE_Data      , only : LE_Data_Enable
    use LE_Config    , only : outputdir
    
    ! --- in/out ---------------------------------

    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_Emis_Dust_Wind_Init'

    ! --- local ------------------------------------------

    integer               ::  imode
    integer               ::  issd
    integer               ::  i_thrf
    integer               ::  i_ustar
    real                  ::  ustar_thrf
    real                  ::  ustar
    !real                  ::  em(nmode)
    real                  ::  em
    character(len=1024)   ::  filename
    integer               ::  fu
    logical               ::  exist
    
    ! --- begin ------------------------------------------

    ! enable data:
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'swg', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'sd', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'ustar', status )
    IF_NOTOK_RETURN(status=1)
    
    ! output file with sampled threshold factor 
    ! related to clay fraction and soil water;
    ! same on all domains, so root only:
    if ( goc_dom%root ) then
      ! output file:
      write (filename,'(3a)') trim(outputdir), pathsep, 'dust-threshold-factor-clay.csv'
      ! file unit:
      call goGetFU( fu, status )
      IF_NOTOK_RETURN(status=1)
      ! create and write table:
      call Friction_Velocity_Threshold_Factor_Soilwater_Demo( fu, filename, status )
      IF_NOTOK_RETURN(status=1)
    end if ! root

    ! read potential dust sources file and regrid to LE
    write (gol, '("  read potential dust sources file ...")'); call goPr
    call dust_wind_potential_sources( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! init 4D look-up-table for evaluated integrals:
    call LUT_Init( lut, 4, status, name='Fv', units='kg/m2/s' )
    IF_NOTOK_RETURN(status=1)

    ! define axis:
    call LUT_Def_Axis( lut, 1, 'ustar'     , n_ustar     , 0.0, d_ustar     , status )
    IF_NOTOK_RETURN(status=1)
    call LUT_Def_Axis( lut, 2, 'ustar_thrf', n_ustar_thrf, 0.0, d_ustar_thrf, status )
    IF_NOTOK_RETURN(status=1)
    call LUT_Def_Axis( lut, 3, 'ssd' , nssd , 0.5, 1.0, status )
    IF_NOTOK_RETURN(status=1)
    call LUT_Def_Axis( lut, 4, 'mode', 1, 0.5, 1.0, status )
    IF_NOTOK_RETURN(status=1)
    
    ! data file:
    call ReadRc( rcF, 'le.emis.dust-wind.lut', filename, status )
    IF_NOTOK_RETURN(status=1)

    ! present ?
    inquire( file=trim(filename), exist=exist )
    !! testing, always re-compute:
    !exist = .false.
    ! read from file ?
    if ( exist ) then

      ! info ...
      write (gol,'("  read wind blown dust look-up-table ...")'); call goPr
      write (gol,'("    table file : ",a)') trim(filename); call goPr
    
      ! read:
      call LUT_Read( lut, trim(filename), status )
      IF_NOTOK_RETURN(status=1)
      
    else

      ! info ...
      write (gol,'("  fill wind blown dust look-up-table ...")'); call goPr

      ! loop over soil size modes:
      do issd = 1, nssd

        ! loop over u* threshold factors:
        do i_thrf = 1, n_ustar_thrf

          ! current value:
          call LUT_Coordinate( lut, 2, i_thrf, ustar_thrf, status )
          IF_NOTOK_RETURN(status=1)

          ! loop over u* values:
          do i_ustar = 1, n_ustar

            ! current value:
            call LUT_Coordinate( lut, 1, i_ustar, ustar, status )
            IF_NOTOK_RETURN(status=1)

            ! evaluate integral over soil particle diameters ;
            ! assume constant soil and air densities:
            call Emis_Dust_Fv_ssd( issd, rho_mineral_soil, &
                                    rho_dry_air_20C_1013hPa, ustar, ustar_thrf, &
                                    em, status )
            IF_NOTOK_RETURN(status=1)

            ! safety check ...           
            if ( em < 0.0 ) then
              write (gol,'("found negative dust emission")'); call goErr
              TRACEBACK; status=1; return
            end if

            ! emission is not aerosol size dependent (at this step)
            call LUT_Set_Value( lut, (/i_ustar,i_thrf,issd,1/), em, status )
            IF_NOTOK_RETURN(status=1)

          end do  ! u*

        end do  ! u* threshold factors

      end do  ! soil size modes
      
      ! info ...
      write (gol,'("    write to file : ",a)') trim(filename); call goPr
      
      ! write:
      call LUT_Write( lut, trim(filename), status )
      IF_NOTOK_RETURN(status=1)
      
    end if  ! read or compute
    
    !! output dir:
    !call ReadRc( rcF, 'le.output.outdir', outdir, status )  
    !IF_NOTOK_RETURN(status=1)
    !! dump:
    !call nc_dump( trim(outdir)//'/erodible_soil_texture.nc', erodible_soil_texture, &
    !                'texture', (/'lon','lat','lbg','texture'/), status ) 
    !IF_NOTOK_RETURN(status=1)
    !! testing ...
    !write (gol,'("break after dumping texture")'); call goPr
    !TRACEBACK; status=1; return
    
    ! storage:
    allocate( ustarthreshold(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Dust_Wind_Init


  ! ***


  subroutine LE_Emis_Dust_Wind_Done( status )

    use GO, only : LUT_Done

    ! --- in/out ---------------------------------

    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Dust_Wind_Done'

    ! --- begin ----------------------------------

    ! clear:
    call LUT_Done( lut, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( ustarthreshold, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Dust_Wind_Done


  ! ***
  
  
  subroutine LUT_Write( lut, filename, status )

    use GO, only : T_LUT  
    use MDF
    
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(in)         ::  lut
    character(len=*), intent(in)    ::  filename
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Write'
    
    ! --- local ----------------------------------
    
    integer                 ::  ncid
    integer                 ::  idim
    integer, allocatable    ::  dimids(:)
    integer, allocatable    ::  varids(:)
    integer                 ::  varid
    real, allocatable       ::  values(:)
    integer                 ::  i, n
    character(len=256)      ::  name
  
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. lut%defined ) then
      write (gol,'("table not completely defined yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( dimids(lut%ndim), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( varids(lut%ndim), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! create file:
    call MDF_Create( filename, MDF_NETCDF, MDF_REPLACE, ncid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, lut%ndim
      ! define dimension:
      call MDF_Def_Dim( ncid, trim(lut%dimname(idim)), lut%shp(idim), dimids(idim), status )
      IF_NOTOK_RETURN(status=1)
      ! define coordinate variable:
      call MDF_Def_Var( ncid, trim(lut%dimname(idim)), MDF_FLOAT, &
                          dimids(idim:idim), varids(idim), status )
      IF_NOTOK_RETURN(status=1)
    end do ! dims
    
    ! table name:
    if ( len_trim(lut%name) > 0 ) then
      name = trim(lut%name)
    else
      name = 'data'
    end if
    ! define table:
    call MDF_Def_Var( ncid, trim(name), MDF_FLOAT, dimids(1:lut%ndim), varid, status )
    IF_NOTOK_RETURN(status=1)
    if ( len_trim(lut%long_name) > 0 ) then
      call MDF_Put_Att( ncid, varid, 'long_name', trim(lut%long_name), status )
      IF_NOTOK_RETURN(status=1)
    end if
    if ( len_trim(lut%units) > 0 ) then
      call MDF_Put_Att( ncid, varid, 'units', trim(lut%units), status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! end of definition phase:
    call MDF_EndDef( ncid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, lut%ndim
      ! size:
      n = lut%shp(idim)
      ! storage:
      allocate( values(n), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! fill:
      do i = 1, n
        values(i) = lut%offset(idim) + lut%scale(idim) * (i-1)
      end do ! values
      ! write:
      call MDF_Put_Var( ncid, varids(idim), values, status )
      IF_NOTOK_RETURN(status=1)
      ! clear:
      deallocate( values, stat=status )
      IF_NOTOK_RETURN(status=1)
    end do ! dims
    
    ! write data:
    select case ( lut%ndim )
      case ( 1 )
        call MDF_Put_Var( ncid, varid, lut%values(:,1,1,1,1,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 2 )
        call MDF_Put_Var( ncid, varid, lut%values(:,:,1,1,1,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 3 )
        call MDF_Put_Var( ncid, varid, lut%values(:,:,:,1,1,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 4 )
        call MDF_Put_Var( ncid, varid, lut%values(:,:,:,:,1,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 5 )
        call MDF_Put_Var( ncid, varid, lut%values(:,:,:,:,:,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 6 )
        call MDF_Put_Var( ncid, varid, lut%values(:,:,:,:,:,:,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 7 )
        call MDF_Put_Var( ncid, varid, lut%values(:,:,:,:,:,:,:), status )
        IF_NOTOK_RETURN(status=1)
      case default
        write (gol,'("unsupported ndim ",i0)'); call goErr
        TRACEBACK; status=1; return
    end select

    ! close file:
    call MDF_Close( ncid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( dimids, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( varids, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LUT_Write


  ! ***
  
  
  !
  ! Read LUT from file.
  ! On input, lut should have the dimensions defined.
  !
  
  subroutine LUT_Read( lut, filename, status )

    use GO, only : T_LUT  
    use MDF
    
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(inout)      ::  lut
    character(len=*), intent(in)    ::  filename
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Read'
    
    ! --- local ----------------------------------
    
    integer                 ::  ncid
    integer                 ::  idim
    integer                 ::  dimid
    integer                 ::  varid
    real, allocatable       ::  values(:)
    integer                 ::  n
    character(len=256)      ::  name
  
    ! --- begin ----------------------------------
    
    ! open file:
    call MDF_Open( trim(filename), MDF_NETCDF, MDF_READ, ncid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, lut%ndim
      ! check ...
      if ( .not. lut%axis_defined(idim) ) then
        write (gol,'("axis ",i0," not defined yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! inquire dimension id:
      call MDF_Inq_DimID( ncid, trim(lut%dimname(idim)), dimid, status )
      IF_NOTOK_RETURN(status=1)
      ! properties:
      call MDF_Inquire_Dimension( ncid, dimid, status, length=n )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( n /= lut%shp(idim) ) then
        write (gol,'("lut axis ",i0," (",a,") has length ",i0," while found ",i0," in file ",a)') &
                 idim, trim(lut%dimname(idim)), lut%shp(idim), n, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      allocate( values(n), stat=status )
      IF_NOTOK_RETURN(status=1)
      ! get variable id:
      call MDF_Inq_VarID( ncid, trim(lut%dimname(idim)), varid, status )
      IF_NOTOK_RETURN(status=1)
      ! read:
      call MDF_Get_Var( ncid, varid, values, status )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( values(1) /= lut%offset(idim) ) then
        write (gol,'("axis definition does not match with file:")'); call goErr
        write (gol,'("  dimension     : ",i0," (",a,")")') idim, trim(lut%dimname(idim)); call goErr
        write (gol,'("  offset        : ",f12.6)') lut%offset(idim); call goErr
        write (gol,'("  coordinate    : ",f12.6)') values(1); call goErr
        write (gol,'("  filename      : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( n > 1 ) then
        if ( abs( (values(2)-values(1)) - lut%scale(idim) ) > 1.0e-4 ) then
          write (gol,'("axis definition does not match with file:")'); call goErr
          write (gol,'("  dimension     : ",i0," (",a,")")') idim, trim(lut%dimname(idim)); call goErr
          write (gol,'("  scale         : ",f12.6)') lut%scale(idim); call goErr
          write (gol,'("  coordinate    : ",f12.6)') values(2)-values(1); call goErr
          write (gol,'("  filename      : ",a)') trim(filename); call goErr
          TRACEBACK; status=1; return
        end if
      end if
      ! clear:
      deallocate( values, stat=status )
      IF_NOTOK_RETURN(status=1)
    end do ! dims
    
    ! table name:
    if ( len_trim(lut%name) > 0 ) then
      name = trim(lut%name)
    else
      name = 'data'
    end if
    ! data:
    call MDF_Inq_VarID( ncid, trim(name), varid, status )
    IF_NOTOK_RETURN(status=1)
    ! storage:
    allocate( lut%values(lut%shp(1),lut%shp(2),lut%shp(3),lut%shp(4),lut%shp(5),lut%shp(6),lut%shp(7)), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( lut%filled(lut%shp(1),lut%shp(2),lut%shp(3),lut%shp(4),lut%shp(5),lut%shp(6),lut%shp(7)), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! read:
    select case ( lut%ndim )
      case ( 1 )
        call MDF_Get_Var( ncid, varid, lut%values(:,1,1,1,1,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 2 )
        call MDF_Get_Var( ncid, varid, lut%values(:,:,1,1,1,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 3 )
        call MDF_Get_Var( ncid, varid, lut%values(:,:,:,1,1,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 4 )
        call MDF_Get_Var( ncid, varid, lut%values(:,:,:,:,1,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 5 )
        call MDF_Get_Var( ncid, varid, lut%values(:,:,:,:,:,1,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 6 )
        call MDF_Get_Var( ncid, varid, lut%values(:,:,:,:,:,:,1), status )
        IF_NOTOK_RETURN(status=1)
      case ( 7 )
        call MDF_Get_Var( ncid, varid, lut%values(:,:,:,:,:,:,:), status )
        IF_NOTOK_RETURN(status=1)
      case default
        write (gol,'("unsupported ndim ",i0)'); call goErr
        TRACEBACK; status=1; return
    end select
    ! set flags:
    lut%filled  = .true.
    lut%defined = .true.

    ! close file:
    call MDF_Close( ncid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LUT_Read
  

  ! ***


  subroutine LE_Emis_Dust_Wind_Add( emis, t, status )

    use GO, only : TDate, ReadRc

    use JAQL_SFC_Soil  , only : rho_mineral_soil, ssd_gmean
    use JAQL_Particles , only : Particle_Volume_Rp
    use JAQL_SFC_Dust  , only : nmode
    use JAQL_SFC_Dust  , only : Friction_Velocity_Threshold_Scale_Roughness
    use JAQL_SFC_Dust  , only : Friction_Velocity_Threshold_Factor_Soilwater
    use JAQL_SFC_Dust  , only : mode_gmean, mode_gsigma, mode_massfrac
 
    use Indices        , only : n_dust, ispecs_dust, ispec_dust_f2, ispec_dust_f1, ispec_dust_c3, ispec_dust_c1, ispec_dust_c2
    use Indices        , only : ispec_duas,ispec_ducs,ispec_duai,ispec_duci
    use Indices        , only : i_naci, i_ncoi
    use Indices        , only : n_seasalt, ispecs_seasalt, ispec_Na_f2, ispec_Na_c3, ispec_Na_f1, ispec_Na_c1, ispec_Na_c2
    use Indices        , only : n_basecation, ispecs_basecation, ispec_K_f, ispec_K_c, ispec_Ca_f, ispec_Ca_c, ispec_Mg_f, ispec_Mg_c
#ifdef with_m7
    use Indices        , only : n_m7
#endif
    use Binas          , only : grav  ! m/s2

    use Dims           , only : runF
    use Dims           , only : nx, ny, nz, nspec

    use LE_LandUse_Data, only : ustar_lu
    use LE_LandUse_Data, only : ustar_lu_dust_emis
    use LE_LandUse_Data, only : ilu_desert, ilu_arable, ilu_water_sea
    use LE_LandUse_Data, only : lu_fracs

    use LE_LandUse_BCatIon, only : bcat_frac, ibcat_Ca, ibcat_Mg, ibcat_K, ibcat_Na
    use LE_LandUse_Soil   , only : n_soiltexture
    use LE_LandUse_Soil   , only : n_erodible
    use LE_LandUse_Soil   , only : erodible_soil, erodible_soil_texture
    use LE_LandUse_Soil   , only : clay_frac
    use LE_LandUse_Soil   , only : i_erodible_arable, i_erodible_bare

    use file_nc
#ifdef with_labeling
    use SA_Labeling       , only : SA_Emis_Setup_Dust
#endif    
    use LE_Data, only : LE_Data_GetPointer
    use Dims   , only : substantial_snowdepth

    ! --- in/out ---------------------------------

    real, intent(inout)       ::  emis(nx,ny,nz,nspec)
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Dust_Wind_Add'

    real, parameter  ::  kg_s_to_ug_min = 1e9  *  60

    ! --- local ----------------------------------

    real        ::  z0s
    integer     ::  imonth
    integer     ::  i, j
    real        ::  frac
    real        ::  rho_air
    real        ::  u_star
    real        ::  u_star_thrf
    real        ::  f_w
    real        ::  f_eff
    real        ::  soil_water
    integer     ::  i_soiltexture
    integer     ::  i_erodible
    real        ::  z0_eff
    real        ::  em
    real        ::  em_mode(nbins)
    real        ::  em_frac(nbins)
    real        ::  frac_transport(nbins)
    integer     ::  itr, ispec
    real        ::  delta_emis
    real        ::  pots, alpha, constant_c
    
    integer     ::  ibin, imode
    integer, parameter :: imode_f1  = 1
    integer, parameter :: imode_f2  = 2
    integer, parameter :: imode_c1  = 3
    integer, parameter :: imode_c2  = 4
    integer, parameter :: imode_c3  = 5

    !real, parameter    ::  arable_factor = 0.02
    real, parameter     ::  arable_factor = 1.0
    
    ! meteo data:
    real, pointer          ::  dens(:,:,:)   ! (lon,lat,alt)
    real, pointer          ::  swg(:,:,:)   ! (lon,lat,1)  
    real, pointer          ::  sd(:,:,:)   ! (lon,lat,1)  
    real, pointer          ::  area(:,:,:)   ! (lon,lat,1)  
    real, pointer          ::  z0dust_lu(:,:,:)   ! (lon,lat,:)  

    ! --- begin ----------------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3')      
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'swg', swg, status, check_units ='kg/kg' )    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'sd', sd, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'area', area, status, check_units ='m2' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'z0dust_lu', z0dust_lu, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    
    ! info ...
    write (gol,'("<LE_Emis_Dust_Wind_Add>")'); call goPr
      
    ! time indices:
    imonth = t%month
    
    ! tuning parameter
    constant_c = 0.5
    
    frac_transport = 0.0
    do ibin = 1, nbins  
        do imode = 1, nmode
            frac_transport(ibin) = frac_transport(ibin) + 0.5 * &
                (   erf( log(binsize_max(ibin)/mode_gmean(imode)) / ( 2**0.5 * log(mode_gsigma(imode)) ) )      -   &
                    erf( log(binsize_min(ibin)/mode_gmean(imode)) / ( 2**0.5 * log(mode_gsigma(imode)) ) )  )   *   &
                    mode_massfrac(imode)
        end do
    end do
   
    !! testing ..
    !write (gol,'("WARNING - dummy u*")'); call goPr
    !ustar_lu = ( (t%day-1)*24 + t%hour ) * 0.01

    ! check:
    if ( any(clay_frac < 0.0) .or. any(clay_frac > 1.0) ) then
      write (gol,'("clay_fracs should bin within [0,1], found range : ",2e16.6)') &
          minval(clay_frac), maxval(clay_frac); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop over grid cells:
    do i = 1, nx
      do j = 1, ny

        ! not erodible at all ? then skip:
        if ( .not. erodible_soil(i,j) ) cycle


        ! soil covered by snow ? then skip;
        if ( sd(i,j,1) > substantial_snowdepth ) cycle

        ! potential sources
        pots = pot_source(i,j)
        
        ! extract air density:
        rho_air = dens(i,j,1)   !  kg/m3

        ! soil water at coastal areas set to 1.0 if grid is at least 10% water
        soil_water = swg(i,j,1)
        if ( lu_fracs(i,j,ilu_water_sea) >= 0.1 ) soil_water = 1.0
        
        ! round:
        if ( (soil_water < 0.0) .and. (abs(soil_water    ) < 1.0e-4) ) soil_water = 0.0
        if ( (soil_water > 1.0) .and. (abs(soil_water-1.0) < 1.0e-4) ) soil_water = 1.0
        ! check ...
        if ( (soil_water < 0.0) .or. (soil_water > 1.0) ) then
          write (gol,'("found strange soil water fraction ",f12.6," in cell (",i0,",",i0,")")') soil_water, i, j; call goPr
          TRACEBACK; status=1; return
        end if

        ! init to zero:
        em_mode = 0.0

        if ( pots == 0.0 ) cycle
        
        ! loop over erodible soils:
        do i_erodible = 1, n_erodible

          ! arable or bare? roughness length and u*
          if ( i_erodible == i_erodible_arable ) then 
            z0_eff = z0dust_lu(i,j,ilu_arable)
            u_star = ustar_lu_dust_emis(i,j,ilu_arable)
          end if
          if ( i_erodible == i_erodible_bare ) then 
            z0_eff = z0dust_lu(i,j,ilu_desert)
            u_star = ustar_lu_dust_emis(i,j,ilu_desert)
          end if
        
          ! Skip if no dust emission (arable land during growing season)
          if ( z0_eff < 0 ) cycle 
          
          ! loop over all textures for this soil:
          do i_soiltexture = 1, n_soiltexture

            ! fraction of this texture in current cell:
            frac = erodible_soil_texture(i,j,i_erodible,i_soiltexture)
            
            ! restirct emissions over arable land
            if ( i_erodible == i_erodible_arable ) then
              frac = frac * arable_factor
            end if
            
            ! skip if not present anyway:
            if ( frac <= 0.0 ) cycle
          
            !
            ! threshold for u* :
            !             u*_s                f_w
            !    u*_thr = ----- f_w  =  u*_s -----
            !             f_eff              f_eff
            ! where:
            !    u*_s    :  small-scale friction velocity u*_s
            !    f_eff   :  scaling number from small-scale friction velocity, <= 1
            !    f_w     : correction for soil water: [1,inf)
            !
            ! See usage in:
            !    JAQL_SFC_Dust_Emis
            !
            
            ! smooth roughness equal to desert roughness
            z0s = 30e-6
            
            ! scaling number for u* due to local roughness elements:
            !   (1 for small z0, towards 0 for rough terain)
            f_eff = Friction_Velocity_Threshold_Scale_Roughness( z0_eff, z0s )
            
            ! u* threshold correction for soil water: [1,inf)
            f_w = Friction_Velocity_Threshold_Factor_Soilwater( soil_water, clay_frac(i_soiltexture) )

            ! total factor:
            u_star_thrf = f_w / f_eff
            
            ! store for output:
            ustarthreshold(i,j) = u_star_thrf

            ! exceeds maximum ? then skip:
            if ( u_star_thrf > max_ustar_thrf ) cycle

            ! get emission for current aerosol mode
            ! for a soil with the given texture:
            call Emis_Dust_Fv( i_soiltexture, rho_mineral_soil, &
                                 rho_air, u_star, u_star_thrf, &
                                 em, &  ! kg/m2/s
                                 status )
            IF_NOTOK_RETURN(status=1)
            
            ! put emissions into aerosol bins
            do ibin = 1, nbins
                em_frac(ibin) = frac_transport(ibin) * em
            end do
            
            ! contribution       kg/m2/s   0-1    m2
            em_mode = em_mode + em_frac * frac * area(i,j,1) ! kg/s
            
          end do  ! textures
        
        end do  ! erodible soils

        ! emission = EM * potential sources * tuning parameter
        em_mode = em_mode * pots * constant_c

    
        ! add to dust tracers:
        do itr = 1, n_dust
          ispec = ispecs_dust(itr)
          select case ( ispec )
            case (ispec_dust_f1)
               delta_emis = em_mode(imode_f1  ) * kg_s_to_ug_min
              emis(i,j,1,ispec) = emis(i,j,1,ispec) + delta_emis
#ifdef with_labeling
              call SA_Emis_Setup_Dust(i,j,ispec,delta_emis,'wind',status)
              IF_NOTOK_RETURN(status=1)
#endif
            case (ispec_dust_f2)
              delta_emis = em_mode(imode_f2  ) * kg_s_to_ug_min
              emis(i,j,1,ispec) = emis(i,j,1,ispec) + delta_emis
#ifdef with_labeling
              call SA_Emis_Setup_Dust(i,j,ispec,delta_emis,'wind',status)
              IF_NOTOK_RETURN(status=1)
#endif
            case (ispec_dust_c1)
                delta_emis = em_mode(imode_c1  ) * kg_s_to_ug_min
                emis(i,j,1,ispec) = emis(i,j,1,ispec) + delta_emis
#ifdef with_labeling
              call SA_Emis_Setup_Dust(i,j,ispec,delta_emis,'wind',status)
              IF_NOTOK_RETURN(status=1)
#endif
            case (ispec_dust_c2)
               delta_emis = em_mode(imode_c2  ) * kg_s_to_ug_min
              emis(i,j,1,ispec) = emis(i,j,1,ispec) + delta_emis
#ifdef with_labeling
              call SA_Emis_Setup_Dust(i,j,ispec,delta_emis,'wind',status)
              IF_NOTOK_RETURN(status=1)
#endif
            case (ispec_dust_c3)
              delta_emis = em_mode(imode_c3) * kg_s_to_ug_min
              emis(i,j,1,ispec) = emis(i,j,1,ispec) + delta_emis
#ifdef with_labeling
              call SA_Emis_Setup_Dust(i,j,ispec,delta_emis,'wind',status)
              IF_NOTOK_RETURN(status=1)
#endif

#ifdef with_m7
            ! M7 species: for mass take usual convention, 
            !   for number convert using density in kg/m3, radius in m
            !   and take into account that emissions are in ug (=1e-9 kg)
            case (ispec_duai)
              delta_emis = (em_mode(imode_f2  ) + em_mode(imode_f1  ) ) * kg_s_to_ug_min
              emis(i,j,1,ispec ) = emis(i,j,1,ispec ) + delta_emis !mass 
              emis(i,j,1,i_naci) = emis(i,j,1,i_naci) + delta_emis*1.0e-9/(Particle_Volume_Rp(0.17e-6)*rho_mineral_soil)
            case (ispec_duci)
              delta_emis = (em_mode(imode_c3  )  + em_mode(imode_c2  ) + &
                            em_mode(imode_c1  ) )* kg_s_to_ug_min
              emis(i,j,1,ispec ) = emis(i,j,1,ispec ) + delta_emis !mass
              emis(i,j,1,i_ncoi) = emis(i,j,1,i_ncoi) + delta_emis*1.0e-9/(Particle_Volume_Rp(2.0e-6)*rho_mineral_soil)
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

        ! base-cat-ions enabled ? only Na is not enough ...
        if ( n_basecation > 5 ) then
          ! add to basecation tracers:
          do itr = 1, n_basecation
            ! tracer index:
            ispec = ispecs_basecation(itr)
            ! switch:
            select case ( ispec )
              case(ispec_Na_f1); emis(i,j,1,ispec) = emis(i,j,1,ispec) + em_mode(imode_f1) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Na)
              case(ispec_Na_f2); emis(i,j,1,ispec) = emis(i,j,1,ispec) + em_mode(imode_f2) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Na)
              case(ispec_Na_c1); emis(i,j,1,ispec) = emis(i,j,1,ispec) + em_mode(imode_c1) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Na)
              case(ispec_Na_c2); emis(i,j,1,ispec) = emis(i,j,1,ispec) + em_mode(imode_c2) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Na)
              case(ispec_Na_c3); emis(i,j,1,ispec) = emis(i,j,1,ispec) + em_mode(imode_c3) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Na)
              
              case(ispec_Ca_f); emis(i,j,1,ispec) = emis(i,j,1,ispec) + (em_mode(imode_f2  ) + em_mode(imode_f1  ) ) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Ca)
              case(ispec_Ca_c); emis(i,j,1,ispec) = emis(i,j,1,ispec) + (em_mode(imode_c3  )  + em_mode(imode_c2  ) + &
                            em_mode(imode_c1  ) ) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Ca)
              case(ispec_Mg_f); emis(i,j,1,ispec) = emis(i,j,1,ispec) + (em_mode(imode_f2  ) + em_mode(imode_f1  ) ) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Mg)
              case(ispec_Mg_c); emis(i,j,1,ispec) = emis(i,j,1,ispec) + (em_mode(imode_c3  )  + em_mode(imode_c2  ) + &
                            em_mode(imode_c1  ) ) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_Mg)
              case(ispec_K_f ); emis(i,j,1,ispec) = emis(i,j,1,ispec) + (em_mode(imode_f2  ) + em_mode(imode_f1  ) ) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_K )
              case(ispec_K_c ); emis(i,j,1,ispec) = emis(i,j,1,ispec) + (em_mode(imode_c3  )  + em_mode(imode_c2  ) + &
                            em_mode(imode_c1  ) ) * kg_s_to_ug_min * bcat_frac(i,j,ibcat_K )
              case default
                write (gol,'("unsupported tracer index ",i6)') ispec; call goErr
                TRACEBACK; status=1; return
            end select
          end do
        end if

      end do    ! i
    end do  ! j
    
    ! info ...
    write (gol,'("<end LE_Emis_Dust_Wind_Add>")'); call goPr

    ! ok
    status = 0

  end subroutine LE_Emis_Dust_Wind_Add


  ! ***


  !
  ! Emission of dust following Alfaro&Gomes in single model cell.
  ! Input:
  !   imode         : aerosol mode (fine,coarse,extreme)
  !   i_soiltexture  : soil texture (coarse,medium,medium fine,fine,very fine)
  ! Ouptut:
  !   em(nmode)     :  emission (kg/m2/s) per aerosol mode fine/course/extreme
  !   status        :  non-zero in case of error
  !

  subroutine Emis_Dust_Fv( i_soiltexture, rho_s, &
                           rho_a, ustar, ustar_thrf, &
                           em, status )

    use GO           , only : LUT_Index, LUT_Get_Value
    use JAQL_SFC_Dust, only : nmode
    use JAQL_SFC_Soil, only : nssd, ssd_gmean, ssd_gsigma
    use LE_LandUse_Soil , only : ssd_frac, clay_frac

    ! --- in/out ---------------------------------

    integer, intent(in )      ::  i_soiltexture   ! 1:n_soiltexture
    real, intent(in)          ::  rho_s           ! soil density (kg/m3)
    real, intent(in)          ::  rho_a           ! air density (kg/m3)
    real, intent(in)          ::  ustar           ! friction velocity (m/s)
    real, intent(in)          ::  ustar_thrf   ! friction velocity threshold factor
    !real, intent(out)         ::  em(nmode)       ! kg/m2/s
    real, intent(out)         ::  em
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Emis_Dust_Fv'

    ! --- local ----------------------------------

    integer     ::  issd
    !real        ::  em_ssd(nmode)
    real        ::  em_ssd
    integer     ::  imode
    integer     ::  i_ustar
    integer     ::  i_thrf
    real        ::  alpha

    ! --- begin ----------------------------------

    ! nearest indices:
    call LUT_Index( lut, 1, ustar, i_ustar, status )
    IF_NOTOK_RETURN(status=1)
    call LUT_Index( lut, 2, ustar_thrf, i_thrf, status )
    IF_NOTOK_RETURN(status=1)

    ! init result:
    em = 0.0
    call LUT_Get_Value( lut, (/i_ustar,i_thrf,i_soiltexture,1/), em_ssd, status )
    em = em_ssd

    ! ok
    status = 0

  end subroutine Emis_Dust_Fv


  ! *


  !
  ! Emission is integral over all possible soil particle diameters:
  !
  !                             dNp
  !  Fv = int Fh(..) alpha(..) ----- dlnDp
  !       lnDp                 dlnDp
  !

  subroutine Emis_Dust_Fv_ssd( issd, rho_s, &
                                rho_a, ustar, ustar_thrf, &
                                em, status )

    use JAQL_Particles, only : Particle_dNp_dlnDp
    use JAQL_SFC_Dust , only : JAQL_SFC_Dust_Emis
    use JAQL_SFC_Soil, only  : nmodes
    use JAQL_SFC_Soil, only  : ssd_gmean, ssd_gsigma, ssd_mfrac
    use Binas,        only   : pi

    ! --- in/out ---------------------------------

    integer, intent(in)       ::  issd           ! geometric std.dev.  (-)
    real, intent(in)          ::  rho_s           ! soil density (kg/m3)
    real, intent(in)          ::  rho_a           ! air density (kg/m3)
    real, intent(in)          ::  ustar           ! friction velocity (m/s)
    real, intent(in)          ::  ustar_thrf      ! friction velocity threshold factor
    real, intent(out)         ::  em
    integer, intent(out)      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/Emis_Dust_Fv_ssd'

    ! number of bins:
    integer, parameter :: maxclass = 196

    ! --- local ----------------------------------

    integer       ::  ibin, idp, imodes
    integer       ::  iclass, nclass
    real          ::  lnDp_offset
    real          ::  dlnDp
    real          ::  Dp
    real          ::  dNp_dlnDp
    !real          ::  em_Dp(nmode)
    real          ::  em_Dp
    !real          ::  dln_Dp
    real          ::  S_total, dM, Dstep
    real          ::  dS
    !real          ::  B
    real          ::  gmean, gsigma, mfrac
    real          ::  S_class(maxclass), S_rel(maxclass)
    

    ! --- begin ----------------------------------


    ! 
    ! Following Ina Tegen (2002)
    !
    
    ! init storage:
    S_class = 0.0
    ! init total:
    S_total = 0.0
    !B = 0.0
    ! step size to next bin:
    !dln_Dp = ( log(2e-3) - log(1e-8) ) / 1000.0
    Dstep = 0.06527 ! [m]  ! this seems a dlnDp !
    ! diameter of first bin:
    Dp = 1e-8       ! [m]
    ! init bin counter:
    nclass = 0
    do
        ! increase counter:
        nclass = nclass + 1
        ! check ...
        if ( nclass > maxclass ) then
          write (gol,'("array size `nclass` too small; use dynamic allocation instead")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! init sum:
        dM = 0.0
        ! loop over size modes:
        do imodes = 1, nmodes
            ! parameters of paricle mass distribution:
            gmean = ssd_gmean(issd, imodes)
            gsigma = ssd_gsigma(issd, imodes)
            ! fraction of mass in this mode:
            mfrac = ssd_mfrac(issd, imodes)
            ! Contribution to total mass?
            ! Mass distribution as function of diameter:
            !   dm(Dp)/dlnDp = rho * dV(Dp)/dlnDp
            !                = rho * f( lnDp ; log(gmean), log(gsigma)**2 )
            ! with
            !   f(x;mu,sigma2) =  1/sqrt(2 pi sigma2) exp(-(x-mu)**2/(2 sigma2))
            ! Here something is computed that seems equal to:
            !   dm(Dp)  3/rho/Dp
            ! This and the name "dS" below suggest some conversion to surface
            ! (factor 2 missing, but finally only fracion of total is used)
            dM = dM + mfrac / ( sqrt(2.0*pi) * log(gsigma) ) *          &
                   exp( ((log( Dp ) - log( gmean))**2) / (-2.0 * (log(gsigma))**2) ) &
                    * Dstep / ( 2.0/3.0 * rho_s * Dp/2.0)
        end do
        ! copy (why?)
        dS = dM 
        ! store:
        S_class(nclass) = dS
        ! update total:
        S_total = S_total + dS !* (1e-8 * exp( idp*dln_Dp ) - Dp)
        ! next bin:
        Dp = Dp * exp( Dstep )   !  lnDp = lnDp + Dstep
        ! leave?
        if ( Dp > 2e-3 ) exit
    end do 

    ! init sum:
    em = 0.0
    ! init array to zero:
    S_rel = 0.0
    ! first diameter: DANGEROUS, should be same as above ...
    Dp = 1e-8
    ! loop over diameter bins:
    do iclass = 1, nclass
        ! fraction of total "surface" (?) that is in this bin ;
        ! shouldn't this be the number fraction instead?
        S_rel(iclass) = S_class(iclass) / S_total
        ! dust emission for this diamater:
        call JAQL_SFC_Dust_Emis(issd, Dp, rho_s, rho_a, ustar, ustar_thrf, &
                                   em_Dp )!, Fh_Dp, alpha_Dp )
        ! add tot total emission:
        em = em + em_Dp * S_rel(iclass)
        ! next:
        Dp = Dp * exp( Dstep )
    end do
    
    ! ok
    status = 0

  end subroutine Emis_Dust_Fv_ssd
  

end module LE_Emis_Dust_Wind
