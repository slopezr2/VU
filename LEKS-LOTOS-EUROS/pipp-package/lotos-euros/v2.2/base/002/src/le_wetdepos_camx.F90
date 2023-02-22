!###############################################################################
!
! Wet scavenging, module from CAMX
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


module LE_WetDepos_Camx

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private
  
  public  ::  LE_WetDepos_Camx_Init, LE_WetDepos_Camx_Done
  public  ::  LE_WetDepos_Camx_Apply


  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'LE_WetDepos_Camx'

  real, parameter :: xl_threshold    = 1.0e-5  ! new value after problems with RACMO meteo;
                                               ! alternative would be to change the cloud conversion
  
  real, allocatable :: henry0(:)
  real, allocatable :: Tfact(:)  
  real, allocatable :: difrat(:)                                             
  real, allocatable :: rscale(:)
  real, allocatable :: bdnl(:)


contains
  
  
  ! ====================================================================


  subroutine LE_WetDepos_Camx_Init( rcF, status )
    
    use GO     , only : TrcFile, ReadRc
    use GO     , only : GoSplitString, GoGetFu, GoMatchValue
    use dims, only : nspec
    use indices, only : specname
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------
    
    type(TrcFile), intent(in)             ::  rcF
    integer, intent(out)                  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Camx_Init'
    
    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 10
    
    ! --- local -------------------------------
    
    character(len=512)      ::  wetscav_datafile
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    integer                 ::  fu
    logical                 ::  exist
    
    real                    ::  value
    integer                 ::  iline, iheader, ispec
    character(len=1024)     ::  line

    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    
    integer                 ::  nfield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    integer                 ::  ifield
    
    integer                 ::  ifield_tracer
    integer                 ::  ifield_lbnd
    integer                 ::  ifield_hlaw
    integer                 ::  ifield_Tfac
    integer                 ::  ifield_diffrat
    integer                 ::  ifield_rscale
    integer                 ::  ifield_reactivity

    ! --- begin -------------------------------
    
    ! allocate
    allocate ( henry0(nspec) )
    allocate ( Tfact(nspec) )
    allocate ( difrat(nspec) )
    allocate ( rscale(nspec) )
    allocate ( bdnl(nspec) )
    
    henry0 = 0.0
    Tfact  = 0.0
    difrat = 0.0
    rscale = 1.0
    bdnl   = 0.0
    
    ! enable data:
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'rh', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'p', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'dens', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'icc', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'bcc', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'iclwc', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'raini', status )
    IF_NOTOK_RETURN(status=1)
    
    ! data file for wet scavenging parameters
    call ReadRc( rcF, 'le.wet_depos.CAMx.datafile', wetscav_datafile, status )
    IF_NOTOK_RETURN(status=1)
    
    ! file should be present:
    inquire( file=trim(wetscav_datafile), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found : ",a)') trim(wetscav_datafile); call goErr
      TRACEBACK; status=1; return
    end if

    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( fu, file=trim(wetscav_datafile), status='old', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening file : ",a)') trim(wetscav_datafile); call goErr
      TRACEBACK; status=1; return
    end if

    ! line counter:          
    iline = 0
    
    ! comment character:
    comment = '#'
    
    ! seperation character:
    sep = ';'

    ! read header line after first comment:
    do
      ! read line:
      read (fu,'(a)',iostat=status) line
      if (status/=0) then
        write (gol,'("reading header line from file : ",a)') trim(wetscav_datafile); call goErr
        TRACEBACK; status=1; return
      end if
      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle
      ! found non-comment line, leave loop:
      exit
    end do

    ! split:
    call goSplitString( line, nheader, headers, status, sep=sep )
    IF_NOTOK_RETURN(status=1)
    
    ! dummy values to avoid warnings ...
    ifield_tracer     = -999
    ifield_lbnd       = -999
    ifield_hlaw       = -999
    ifield_Tfac       = -999
    ifield_diffrat    = -999
    ifield_rscale     = -999
    ifield_reactivity = -999

    ! loop over headers
    do iheader = 1, nheader
      ! current:
      header = headers(iheader)
      ! which column ?
      select case ( trim(header) )
        ! tracer
        case ( 'Tracer' )
          ifield_tracer = iheader
        ! lower bnd
        case ( 'lower_bnd' )
          ifield_lbnd = iheader
        ! Henry laws constant
        case ( 'H_Law' )
          ifield_hlaw = iheader
        ! Temperature Fact
        case ( 'T_Fact' )
          ifield_Tfac = iheader
        ! Diffusion ratio
        case ( 'Diffrat' )
          ifield_diffrat = iheader
        ! reactivity
        case ( 'Reactivity' )
          ifield_reactivity = iheader
        ! Rscale
        case ( 'Rscale' )
          ifield_Rscale = iheader 
        ! unknown 
        case default
          write (gol,'("unsupported header `",a,"`")') trim(header); call goErr
          TRACEBACK; status=1; return
      end select        
    end do 

    ! loop over records:
    do

      ! increase record counter:
      iline = iline + 1
      ! try to read line:
      read (fu,'(a)',iostat=status) line
      ! eof ?
      if (status<0) exit
      ! error ?
      if (status>0) then
        write (gol,'("reading line ",i6," from file : ",a)') iline, trim(wetscav_datafile); call goErr
        TRACEBACK; status=1; return
      end if

      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle

      ! split into records:
      call goSplitString( line, nfield, fields, status, sep=sep )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( nfield /= nheader ) then
        write (gol,'("number of fields (",i2,") in line ",i2," :")') nfield, iline; call goPr
        write (gol,'("  ",a)') trim(line); call goErr
        write (gol,'("fields:")'); call goErr
        do ifield = 1, nfield
          write (gol,'(i6," : ",a)') ifield, trim(fields(ifield)); call goErr
        end do
        write (gol,'("differs from number of headers ",i2," in file ",a)') nheader, trim(wetscav_datafile); call goErr
        TRACEBACK; status=1; return
      end if
    
      ! get tracer code
      call goMatchValue( trim(fields(ifield_tracer)), specname, ispec, status, quiet=.true. )
      if ( status /= 0 ) then
        write ( gol, '("Tracer in data file:, ",a, ", not in model run")' ) fields(ifield_tracer) ; call goPr
        cycle
      end if
      
      do ifield = 1, nfield
        
        ! skip tracer field
        if (ifield == ifield_tracer ) cycle

        field = fields(ifield)
        ! read field:
        read (field,*,iostat=status) value
        if (status/=0) then
          write (gol,'("reading fraction from `",a,"`")') field; call goErr
          TRACEBACK; status=1; return
        end if
        
        !~ henrys law
        if ( ifield == ifield_hlaw ) then
          henry0(ispec) = value
        !~ Temperature factor
        else if ( ifield == ifield_Tfac ) then
          Tfact(ispec)  = value
        ! diffusivity rate
        else if ( ifield == ifield_diffrat ) then
          difrat(ispec) = value
        ! rscale value  
        else if ( ifield == ifield_rscale ) then  
          rscale(ispec) = value
        ! lower bound
        else if ( ifield == ifield_lbnd ) then
          bdnl(ispec)   = value

        end if ! which header
      end do ! fields
    end do ! lines
    
    ! close file
    close( fu, iostat=status)
    if (status/=0) then
      write (gol,'("closing file : ",a)') trim(wetscav_datafile); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! which tracers are not filled with henry constant?
    do ispec = 1, nspec
      if ( henry0(ispec) == 0.0 ) then
        write( gol, '(" No henry diffusion coeff for tracer: ",a)') specname(ispec) ; call goPr
      endif
    enddo  

    ! ok
    status = 0

  end subroutine LE_WetDepos_Camx_Init

  
  ! ====================================================================


  subroutine LE_WetDepos_Camx_Done( status )
  

    ! --- in/out ------------------------------

    integer, intent(out)                  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Camx_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
    deallocate( henry0 )
    deallocate( Tfact )
    deallocate( difrat )
    deallocate( rscale )
    deallocate( bdnl )
    
    ! ok
    status = 0

  end subroutine LE_WetDepos_Camx_Done

  
  ! =====================================================================
  
  subroutine LE_WetDepos_Camx_Apply( c, cwet, update_cwet, dt, n, status )
  
    use Binas,    only  : xm_air, Rgas
    use Indices,  only  : rhopart, tracer_is_aerosol, ispec2aerosol
    use Indices,  only  : specname, specunit
    use dims,     only  : nx, ny, nz, nspec
    use dims,     only  : pH
    use JAQL,     only  : get_hplus
    use JAQL,     only  : dissociation_henry
    use JAQL,     only  : henry_func
    use JAQL,     only  : Rh_Growth
    
    use LE_Particle_Data, only : partsize
    
        
    use indices , only : i_nh3, i_so2, i_hno3
    use indices , only : i_so4a_f, i_nh4a_f, i_no3a_f
    use indices , only : i_h2o2, i_o3

    use LE_Data      , only : LE_Data_GetPointer
    
    ! --- in/out ------------------------------
    
    real, intent(inout)   ::  c(nx,ny,nz,nspec)
    real, intent(inout)   ::  cwet(nx,ny,nz,nspec)
    logical, intent(in)   ::  update_cwet
    real, intent(in)      ::  dt   
    integer, intent(in)   ::  n  ! step number of wet deposition 
    integer, intent(out)  ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_WetDepos_Camx_Apply'
    real, parameter             ::  thr_rain__mm_hr = 0.1  ! thres hold for rain (mm/hr)

    ! assumed CO2 concentration (320 ppm)
    real, parameter             ::     co2 = 320.0e-6       ! mole/mole
    ! TO BE TESTED: SHOULDN'T THIS BE USED ?
    !real, parameter             ::  ch_co2 = 320.0e3        ! ppb

    real, parameter             ::  Tmin  = 243.       ! minimum temperature for clud water
    real, parameter             ::  densfac = 44.9
    real, parameter             ::  rhoh2o = 1.e6      ! water density [g/m3]
    
    ! default pH of rain water, somewhat acid due to presence of CO2:
    real, parameter             ::  pH_rain = 5.0
    ! minimum value for pH, this is already extreme acid ...
    real, parameter             ::  pH_min = 1.0
    
    ! --- local -------------------------------
    
    ! 
    integer        ::  ix, iy, iz, ispec, iaerosol
    logical        ::  ltop
    real           ::  dh
    
    real           ::  f_in, f_below
    real           ::  loss_prev, loss_curr
    real           ::  rate
    real           ::  cov_nc
    real           ::  volrat
    real           ::  ch, ch_old
    real           ::  delc, delr_in, delr_below
    real           ::  delc_in, delc_below
    real           ::  cwat
    real           ::  pwat
    real           ::  c0(nspec)
    real           ::  convfac
    real           ::  cmin, ceq
    real           ::  totc, totw
    real           ::  delc0
    real           ::  gdiff
    real           ::  psize
    
    ! meteo parameters in a grid cell
    real           ::  zclwc               ! cloud liquid water density in m3 water/m3 air
    real           ::  xl                  ! volume of water in cloud
    real           ::  xliq                ! zclwc/(cloud cover)
    
    ! Concentrations in a grid cell
    real           ::  ch_so2              ! concentration SO2 in ppb
    real           ::  ch_hno3             ! concentration HNO3 in ppb
    real           ::  ch_nh3              ! concentration NH3 in ppb
    real           ::  ch_so4              ! concentration SO4 in ppb
    real           ::  ch_nh4              ! concentration NH4 in ppb
    real           ::  ch_no3              ! concentration NO3 in ppb
    real           ::  ch_sulf             ! concentration sulfate in ppb ( = ch_so4 xxx is dat inderdaad hetzelfde?)
    real           ::  ch_o3               ! concentration O3 in ppb
    real           ::  ch_h2o2             ! concentration H2O2 in ppb
    ! Dissociation and Henry coefficients in a grid cell
    real           ::  ztr                 ! Temperature-related parameter (1/Kelvin)
                                           ! ztr = (1/T - 1/Tref)
    real           ::  dkh2o               ! dissociation constant water
    real           ::  hkco2               ! dimensionless Henry's constant CO2 
    real           ::  dkco2               ! Dissociation constant CO2
    real           ::  hkso2               ! dimensionless Henry's constant for SO2
    real           ::  dkso2               ! Dissociation constant for SO2
    real           ::  dknh3               ! dissociation constant ammonia
    real           ::  hknh3               ! dimensionless Henry's constant ammonia
    real           ::  hkh2o2              ! dimensionless Henry's constant for hydroperoxide
    real           ::  hko3                ! dimensionless Henry's constant for ozone
    real           ::  dkhso3              ! Dissociation constant for HSO3-
    real           ::  hplus               ! concentration H+ (mol/l)
    !
    real           ::  hlaw                ! Henry's Law constant
    !
    real           ::  gscav_in, gscav_below  ! gaseous scavenging coefficients [1/s]
    real           ::  ascav_in, ascav_below  ! aerosol scavenging coefficients [1/s]
    real           ::  drpvel                 ! drop velocity of rain drops
    real           ::  dtfall
    
    real           ::  delta_cwet(nx,ny,nspec)  ! present drydeposition budget over one layer

    ! meteo data:
    real, pointer        ::    h_m(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::   temp(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::     rh(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::   pres(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::   dens(:,:,:)   ! (lon,lat,lev)
    real, pointer        :: cov_ic(:,:,:)   ! (lon,lat,lev)
    real, pointer        :: cov_bc(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::  iclwc(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::  raini(:,:,:)   ! (lon,lat,lev)
    ! conversion:
    !                        (mm/hr)/(m/s)  = mm/m    s/hr
    real, parameter      ::  m_s__to__mm_hr = 1.0e3 * 3600.0
    
    ! --- begin ---------------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 't', temp, status, check_units='K' )
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_GetPointer( 'rh', rh, status, check_units='%' )
    IF_NOTOK_RETURN(status=1)  
    call LE_Data_GetPointer( 'p', pres, status, check_units='Pa' )
    IF_NOTOK_RETURN(status=1)
    !~ air density:
    call LE_Data_GetPointer( 'dens', dens, status, check_units ='kg/m3' )
    IF_NOTOK_RETURN(status=1)
    !~ in-cloud coverage:
    call LE_Data_GetPointer( 'icc', cov_ic, status, check_units ='1' )
    IF_NOTOK_RETURN(status=1)
    !~ below-cloud-coverage: part of the cell that has no cloud,
    !  but is below a cloud in a higher layer.
    call LE_Data_GetPointer( 'bcc', cov_bc, status, check_units ='1' )
    IF_NOTOK_RETURN(status=1)
    !~ cloud liquid water content:
    call LE_Data_GetPointer( 'iclwc', iclwc, status, check_units ='kg/kg' )
    IF_NOTOK_RETURN(status=1)
    !~ rain intensity:
    call LE_Data_GetPointer( 'raini', raini, status, check_units ='m/s' )
    IF_NOTOK_RETURN(status=1)


    ! initialize

    delta_cwet = 0.0
    ! loop over grid cells
    do ix = 1, nx
    do iy = 1, ny
      
      ! hihest layer with cloud coverage
      ltop = .true.

      ! levels
      do iz = nz, 1, -1
          
        ! check for rain
        if ( raini(ix,iy,iz)*m_s__to__mm_hr >= thr_rain__mm_hr ) then
        
          ! still above cloud ?
          if ( (cov_ic(ix,iy,iz) + cov_bc(ix,iy,iz)) == 0.0 ) cycle
          
          if ( iz > 1 ) then
            dh = h_m(ix,iy,iz) - h_m(ix,iy,iz-1)  ! (m)
          else                              
            dh = h_m(ix,iy,iz)                    ! (m)
          end if
                  
          ! ENVIRON: volrat=regengehalt/rhoh2o=1.0e-7*(rr(k))**100/127   
          volrat=1.0e-7*( abs(raini(ix,iy,iz)*m_s__to__mm_hr)**(100./127.) ) 

          ! Get meteo parameters in current grid cell, to calculate Henry, hplus
          ! zclwc, xl, xliq (defined in module header):
          call get_meteo( temp(ix,iy,iz),cov_bc(ix,iy,iz), dens(ix,iy,iz), &
                          iclwc(ix,iy,iz), cov_ic(ix,iy,iz), &
                          zclwc, xl, xliq )
                          
          ! no cloud coverage (not in-incloud, and not below a cloud that is in a higher layer);
          ! clip values that are slightly outside [0,1] due to rounding errors:
          cov_nc = min( max( 0.0, 1.0 - cov_ic(ix,iy,iz) - cov_bc(ix,iy,iz) ), 1.0 )

          ! Get concentrations in current grid cell (all in ppb);
          ! ch_so2, ch_hno3, ch_nh3, ch_so4, ch_nh4, ch_no3, ch_sulf, ch_o3, ch_h2o2
          ! (defined in module header):
          ! only if gas and sia tracers in model
          if ( i_so2 > 0 .and. i_so4a_f > 0 .and. i_nh3 > 0 .and. &
               i_o3  > 0 .and. i_hno3 > 0 .and. i_nh4a_f > 0 .and. &
               i_h2o2 > 0 .and.i_no3a_f > 0  ) then

            ! convert concentrations (also aerosols!) to ppb
            call get_conc( c(ix,iy,iz,:), dens(ix,iy,iz), &
                           ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2, &
                           ch_so4, ch_sulf, ch_nh4, ch_no3 )


            ! Compute dissociation and Henry coefficients;
            ! dkh2o, hkco2, dkco2, hkso2, dknh3, hknh3, hkh2o2, hko3, dkso2, dkhso3
            ! (defined in module header):
            call dissociation_henry( temp(ix,iy,iz), ztr, &
                                     dkco2, dkso2, dknh3, dkhso3, dkh2o, &
                                     hkco2, hkso2, hknh3, hko3, hkh2o2 )

            ! Get H+ concentration (mol/l):
            !  - input concentrations in ppb
            !    SHOULD 'ch_co2' PARAMETER BE USED HERE?
            !  - minimum value of 1e-10 ?
            !  - sometimes returns [H+] < 0.0 due to numerical errors
            !  - sometimes retunrs extreme large [H+] that lead to pH ~ -4 ...
            hplus = get_hplus(ch_so2, ch_so4, ch_nh3, ch_nh4, ch_hno3, ch_no3, co2, &
                            dkco2, dkso2, dknh3, dkhso3, dkh2o, &
                            hkco2, hkso2, hknh3, &
                            xl)

            ! calculate pH:
            if ( hplus > 0.0 ) then
              ! compute pH ; use minimum value of 1 which is extreme strong acid already ..
              pH(ix,iy,iz) = max( pH_min, -log10(hplus) )
            else
              ! something wrong, [H+] < 0.0 ...
              ! use default value for rain water:
              pH(ix,iy,iz) = pH_rain
            end if
            
          else
          
            ! use default value for rain water:
            pH(ix,iy,iz) = pH_rain

          end if
          
          ! loop over tracers
          do ispec = 1, nspec


            ! check aerosol
            if ( tracer_is_aerosol(ispec) ) then
              
              ! help array with concentration
              ch = c(ix,iy,iz,ispec)
              
              ! aerosol index
              iaerosol = ispec2aerosol(ispec)

              ! diameter enlarged by water agglomeration
              psize = sqrt( Rh_Growth( rh(ix,iy,iz) ) ) * partsize(ix,iy,iz,iaerosol)

              ! calculate scavenging ratios
              call calc_scav_rat(Tmin, raini(ix,iy,iz)*m_s__to__mm_hr, temp(ix,iy,iz), dh, dens(ix,iy,iz), drpvel, status, &
                               prtdia=psize, rhopart=1.e3*rhopart(ispec), &
                               ascav_in=ascav_in, ascav_below=ascav_below )

              
              ! fractional loss term for in and below cloud part of the grid cell
              f_in    = exp(-ascav_in   *dt*60.0)
              f_below = exp(-ascav_below*dt*60.0)
                           
              ! To avoid, too effieciently scavenging due to implicit mixing within the layer
              ! Take rate used as in (Henzing et al., 2006, Size resolved below cloud scavenging of aerosols)
              ! loss terms:
              loss_curr = (f_in**(n  )) * cov_ic(ix,iy,iz) + (f_below**(n  )) * cov_bc(ix,iy,iz) + cov_nc
              loss_prev = (f_in**(n-1)) * cov_ic(ix,iy,iz) + (f_below**(n-1)) * cov_bc(ix,iy,iz) + cov_nc
              ! avoid division by zero ...
              if ( loss_prev > 0.0 ) then
                rate = loss_curr / loss_prev
              else
                rate = 0.0
              end if
                     
              ! testing ...
              if ( rate < 0.0 ) then
                write (gol,'("found negative wet deposition rate:")'); call goErr
                write (gol,'("  cell                   = ",i0,",",i0,",",i0)') ix, iy, iz; call goErr
                write (gol,'("  spec                   = ",i0," (",a,")")') ispec, trim(specname(ispec)); call goErr
                write (gol,*) ' rate                   = ', rate; call goErr
                write (gol,*) ' ch                     = ', ch, ' ', trim(specunit(ispec)); call goErr
                write (gol,*) ' loss_curr, loss_prev   = ', loss_curr, loss_prev; call goErr
                write (gol,*) ' n                      = ', n; call goErr
                write (gol,*) ' f_in, f_below          = ', f_in, f_below; call goErr
                write (gol,*) ' cov_ic, cov_bc, cov_nc = ', cov_ic(ix,iy,iz), cov_bc(ix,iy,iz), cov_nc; call goErr
                TRACEBACK; status=1; return
              end if
              
              ! change:                     
              delta_cwet(ix,iy,ispec) = ch * (1-rate) * dh
              ! update budgets ?
              if ( update_cwet ) then
                ! update wet deposition flux
                cwet(ix,iy,iz,ispec) = cwet(ix,iy,iz,ispec) + delta_cwet(ix,iy,ispec)
              end if             
              
              
              ! update air concentration
              ch = ch * rate
              
              ! save in global array
              c(ix,iy,iz,ispec) = ch

            else
              
              ! skip for not depositing tracers
              if ( henry0(ispec) == 0.0 ) cycle
              
              ! check ...
              if ( pH(ix,iy,iz) < -900.0 ) then
                write (gol,'("no pH defined for this grid cell")'); call goErr
                TRACEBACK; status=1; return
              end if
              
              ! help array for concentration
              ! g tr / mol air  =            *       ppb         * umol tr /mol tr  * kg air / m3 air / (kg air / mol air)
              ch                =    1e-9    * c(ix,iy,iz,ispec) * 1e6              * dens(ix,iy,iz)  / xm_air
              
              ! calculate henry law constant
              hlaw = henry_func( ispec, henry0(ispec), Tfact(ispec), temp(ix,iy,iz), pH(ix,iy,iz) )
              hlaw = hlaw * Rgas * temp(ix,iy,iz)

              ! water concentration in g/m3
              ! g h20 /m3 air   = g/kg * kg h2o/kg air  * kg air / m3 air
              cwat              = 1.e3 * iclwc(ix,iy,iz) * dens(ix,iy,iz)
              ! preciptation water
              pwat = volrat*rhoh2o

              ! rscale defines species behaviour for temperatures below 273K
              if (temp(ix,iy,iz).lt.273. .and. rscale(ispec).gt.0.) then
                 cwat = amax1(0.,iclwc(ix,iy,iz)*(temp(ix,iy,iz) - Tmin)/(273. - Tmin))
              endif

              ! initialize
              delc = 0.

              if (ltop) then 
                 c0(ispec) = 0.                 
              endif 

              convfac = densfac*(273./temp(ix,iy,iz))*(pres(ix,iy,iz)/ 1013.)
              cmin = bdnl(ispec)*convfac
              ch   = amax1(cmin,ch) 
              call calc_scav_rat( Tmin, raini(ix,iy,iz)*m_s__to__mm_hr, temp(ix,iy,iz), dh, dens(ix,iy,iz), drpvel, status, &
                                  cwat=cwat, hlaw=hlaw, ch=ch, difrat=difrat(ispec), rscale=rscale(ispec), &
                                  gscav_in=gscav_in, gscav_below=gscav_below )
              IF_NOTOK_RETURN(status=1)

              totc = ch + c0(ispec)
              totw = cwat + pwat
              ceq = totc/(1. + hlaw*totw/rhoh2o)
              ceq = totc - ceq
              
              ! fractional loss term for in and below cloud part of the grid cell
              f_in    = exp(-gscav_in   *dt*60.0)
              f_below = exp(-gscav_below*dt*60.0)
              
              ! To avoid, too efficient scavenging due to implicit mixing within the layer
              ! Take rate used as in (Henzing et al., 2006, Size resolved below cloud scavenging of aerosols)
              ! loss terms:
              loss_curr = (f_in**(n  )) * cov_ic(ix,iy,iz) + (f_below**(n  )) * cov_bc(ix,iy,iz) + cov_nc
              loss_prev = (f_in**(n-1)) * cov_ic(ix,iy,iz) + (f_below**(n-1)) * cov_bc(ix,iy,iz) + cov_nc
              ! avoid division by zero ...
              if ( loss_prev > 0.0 ) then
                rate = loss_curr / loss_prev
              else
                rate = 0.0
              end if
              
              ! washed out concentration
              delc = (1-rate) * (ceq - c0(ispec) )             
              if (delc.gt.0.) delc = amin1(delc,ch-cmin)

              ! change:
              delta_cwet(ix,iy,ispec)                      =  delc       * dh * 1e9  *       1e-6       / dens(ix,iy,iz)    * xm_air        
              ! update budgets ?
              if ( update_cwet ) then
                ! update deposition budget   ppb m      + umol tr / m3 air * m  *      * mol tr / umol tr / (kg air / m3 air) * kg air / mol air
                cwet(ix,iy,iz,ispec)  = cwet(ix,iy,iz,ispec) +  delta_cwet(ix,iy,ispec)
              end if

              
              ! Update the cell concentration
              ch = ch - delc
              
              ! update in-droplet concentrations, to avoid super saturation
              dtfall = dh/drpvel                             
              delc0 = (ceq - c0(ispec))*(1. - exp(-gscav_below*dtfall)) * cov_bc(ix,iy,iz)         
              if (delc0.gt.0.) delc0 = amin1(delc0,ch-cmin)
              c0(ispec)    = c0(ispec) + delc0

              ! give concentration back to main array
              ! ppb             =            *   umol tr/m3 air  *  mol tr /umol tr  / (kg air / m3 air) * (kg air / mol air)
              c(ix,iy,iz,ispec) =    1e9     *   ch              *       1e-6        / dens(ix,iy,iz)    * xm_air

            end if ! aerosol or gas?
          end do ! tracers   
          
          ! go to lower layer
          ltop = .false.
                                 
        end if ! rain 
        
      end do ! iz
      
    end do ! ix
    end do ! iy
    

    ! ok
    status = 0

  end subroutine LE_WetDepos_Camx_Apply
  
  
  subroutine calc_scav_rat (Tmin, Ri, T, dz, rhoair, drpvel, status, &
                            cwat, hlaw, ch, difrat, rscale, &
                            prtdia, rhopart, &
                            gscav_in, gscav_below, &
                            ascav_in, ascav_below )  
                            
!     calc_scav_ret calculates wet scavenging rates for gases and aerosols.
!     Rates are determined for:
!     1) Uptake of cloud water with dissolved gasses
!     2) Uptake of ambient gasses into precip
!     3) Uptake of cloud water with PM (all PM in cloudy layers is assumed
!     to reside in cloud water)
!     4) Uptake of ambient PM into precip, dependent on particle size and 
!     ice form
!     Super-cooled liquid cloud water is assumed to exist in the temperature
!     range tamin < T < 273K using a ramp function (100% liquid at 273 to 
!     0% liquid at tamin).
!     
  
  use Binas, only : pi
  
  ! --- in/out ----------------
  real, intent(in)    ::  Tmin    ! minimum temperature for liquid cloud water [K]
  real, intent(in)    ::  Ri      ! Rain intensity [mm/hr]
  real, intent(in)    ::  T       ! Temperature [K]
  real, intent(in)    ::  dz      ! Cell depth [m]
  real, intent(in)    ::  rhoair  ! Atmospheric density [kg/m3]

  real, intent(out)   ::  drpvel  ! rain drop fall speed [m/s]
  integer, intent(out)::  status  
  
  ! gas tracer input
  real, intent(in), optional  ::  cwat    ! Cloud Water Content [g/m3]
  real, intent(in), optional  ::  hlaw   ! henry's Law constant
  real, intent(in), optional  ::  ch     ! Gaseous concentration [umol/m3]
  real, intent(in), optional  ::  difrat ! Ratio of H2O to gas diffusivity
  real, intent(in), optional  ::  rscale ! Solubility in liquid clouds below 273.15 K

  ! aerosol tracer input
  real, intent(in), optional  ::  prtdia ! mean aerosol size [m]
  real, intent(in), optional  ::  rhopart ! aerosol density [g/m3]
  
  ! outputs
  real, intent(out), optional ::  gscav_in    ! Gas Scavenging rate in cloud [1/s]
  real, intent(out), optional ::  gscav_below ! Gas Scavenging rate below cloud [1/s]
  real, intent(out), optional ::  ascav_in    ! Aerosol scavening rate in cloud [1/s]
  real, intent(out), optional ::  ascav_below ! Aerosol scavening rate below cloud [1/s]
  
  ! --- const ------------------------
  
  character(len=*), parameter ::  rname = mname//'/Calc_scav_rat'

  real, parameter ::  rhoh2o = 1e-6   ! [g/m3]
  real, parameter ::  difh2o = 2.3e-5 ! [m2/s]
  real, parameter ::  muair  = 1.8e-5 ! [kg/ms]
  real, parameter ::  muh2o  = 1.e-3  ! [kg/ms]
  real, parameter ::  boltz  = 1.38e-23 ! [J/K]
  real, parameter ::  xmfp   = 6.5e-8 ! [m]
  real, parameter ::  cldeff = 0.95   ! neu...vorher 0.9 ??
  
  
  ! --- local ------------------------
  
  real  ::  nuair
  real  ::  drpdia, drpmas
  real  ::  cscav, dscav
  real  ::  eff
  real  ::  cgas, caq
  real  ::  diff
  real  ::  term1, term2, term3
  real  ::  kc
  
  real  ::  reynold
  real  ::  power
  real  ::  scf
  real  ::  difbrwn
  real  ::  schmidt, stoke
  real  ::  top, bot
  real  ::  star
  real  ::  amin1
  real  ::  terma, termb
  real  ::  phi
  
  ! --- begin ------------------------
  
  ! -----Entry point 
  
  if ( present(gscav_in   ) ) gscav_in    = 0.0
  if ( present(gscav_below) ) gscav_below = 0.0
  if ( present(ascav_in   ) ) ascav_in    = 0.0
  if ( present(ascav_below) ) ascav_below = 0.0
  dscav = 0.0
  eff   = 1.0
  
  ! -----Calculate environmental parameters
  nuair = muair/rhoair                  ! air molecular diffusivity [m2/s]
  drpdia = 9.0e-4*(Ri**0.21)            ! rain drop diameter [m]
  drpmas = 1.e6*(pi/6.)*(drpdia)**3     ! rain drop mass [g]
  if (T.gt.260.) then                   ! define no liquid phase temperature
     drpvel = 3.1e3*drpdia              ! rain drop fall speed [m/s]
  else
     drpdia = (29.*1000.*drpmas)**0.56  ! snow diameter [mm]
     drpvel = 0.83*drpdia**0.20         ! snow fall speed [m/s]
     drpdia = 1.e-3*drpdia              ! [m]
  endif
  cscav = 4.2e-7*Ri*cldeff/drpdia       ! cloud scavenging rate [1/s]
  
  
  ! check aerosol of gas scavenging rate  
  if ( present(gscav_in) .and. present(gscav_below) ) then
    
    !     
    ! Gas scavenging
    !

    ! No liquid droplets    
    if (T.le.Tmin) return
    
    !
    cgas = ch
    caq  = 0.
    
    !     
    ! In cloud, partition total gas into equilibrium aqueous and gas phase,
    ! and calculate scavenging rate for dissolved gasses in cloud water
 
    cgas = ch/(1. + hlaw*cwat/rhoh2o)
    caq  = ch - cgas
    dscav = cscav*caq/ch
    
    !     
    ! Below Cloud Calculate scavenging rate for ambient gas dissolving into precip
    !     
      
    if ( T.gt.273. .or. (T.le.273. .and. rscale.eq.0.)) then
       diff        = difh2o/difrat
       term1       = (drpvel*drpdia/nuair)**0.5
       term2       = (nuair/diff)**0.333
       kc          = diff/drpdia*(2. + 0.6*term1*term2)
       gscav_below = 1.67e-6*Ri*kc/(drpdia*drpvel)
    endif
    
    ! Gas scavenging rate  
    !!~ wrong ...
    !gscav_below = gscav_below + dscav
    !gscav_below = gscav_below
    !~ fixed 2021-02:  
    gscav_below = gscav_below
    gscav_in    = dscav

  else if ( present(ascav_in) .and. present(ascav_below) ) then
    
    !     
    ! Aerosol scavenging
    !
    !     
    ! Scavenging rate for aerosol in cloud water is set equal to
    ! cloud water scavenging rate
    ascav_in = cscav
    
    ! Calculate scavenging rate of dry aerosols below cloud as f(size)
    reynold = drpdia*drpvel/(2.*nuair)
    power = amin1(7.6,0.55*prtdia/xmfp)
    scf = 1. + (2.514 + 0.8*exp(-power))*xmfp/prtdia
    difbrwn = boltz*T*scf/(3.*pi*muair*prtdia)
    schmidt = nuair/difbrwn
    stoke = drpvel*prtdia*prtdia*rhopart*scf/(9000.*muair*drpdia)
    top = 1.2 + alog(1. + reynold)/12.
    bot = 1.0 + alog(1. + reynold)
    star = amin1(top/bot,stoke)

    terma = reynold**0.5 * schmidt**0.333
    termb = reynold**0.5 * schmidt**0.5
    term1 = 4./(reynold*schmidt)*(1. + 0.4*terma + 0.16*termb)
    phi = prtdia/drpdia
    term2 = 4.*phi*(muair/muh2o + (1. + 2.*reynold**0.5)*phi)
    term3 = (stoke - star)/(stoke - star + 2./3.)
    term3 = term3**1.5
    eff = term1 + term2 + term3
    if (T.le.273.) eff = amax1(eff,1.e-3)

    ascav_below = cscav*eff/cldeff

  else
    write (gol, '("Aerosol or Gas scavening rate, at least one must be in call statement")') ; call goErr
    TRACEBACK;status=1;return
  end if
  ! ok
  status = 0
  
  end subroutine calc_scav_rat

  ! ***
  
  !-----------------------------------------------------------------------------
  ! Get concentrations in the current grid cell and
  ! convert aerosol concentrations from ug/m3 to ppb

  pure subroutine get_conc( c, dens, &
                       ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2, &
                       ch_so4, ch_sulf, ch_nh4, ch_no3 )
    
    use Binas  , only : xm_air
    use Indices, only : nspec
    use Indices, only : i_so2, i_o3, i_h2o2, i_hno3, i_so4a_f, i_nh3, i_no3a_f, i_nh4a_f
    use Indices, only : specmolm          ! mol mass of elements (kg/mol)

    ! -- in/out ---------------------------------------------------------------

    real, intent(in)  :: c(nspec)     ! concentration (gas in ppb, aerosol in ug/m3)
    real, intent(in)  :: dens         ! air density (kg/m3)
    real, intent(out) :: ch_so2, ch_hno3, ch_nh3, ch_o3, ch_h2o2  ! gas in ppb
    real, intent(out) :: ch_so4, ch_sulf, ch_nh4, ch_no3  ! aerosol in ppb

    ! -- local  ---------------------------------------------------------------

    ! gasesous components:
    ch_so2  = c(i_so2 )  ! ppb
    ch_hno3 = c(i_hno3)  ! ppb
    ch_nh3  = c(i_nh3 )  ! ppb
    ch_o3   = c(i_o3  )  ! ppb
    ch_h2o2 = c(i_h2o2)  ! ppb

    ! aerosols:
    !         (ug aer)/(m3 air) / (kg air)/(m3 air)   (kg air)/(mol air) / (kg aer)/(mol aer)
    ch_so4  =   c(i_so4a_f)     /       dens        *      xm_air        / specmolm(i_so4a_f) ! ug/kg (mol aer)/(mol air) = ppb
    ch_sulf =    ch_so4
    ch_nh4  =   c(i_nh4a_f)     /       dens        *      xm_air        / specmolm(i_nh4a_f)
    ch_no3  =   c(i_no3a_f)     /       dens        *      xm_air        / specmolm(i_no3a_f)

  end subroutine get_conc
  

  ! ***

  !-----------------------------------------------------------------------------
  ! Get meteo parameters in the current grid cell.
  !
  ! Meteo parameters are declared in module header
  ! zclwc               ! cloud liquid water density in m3 water/m3 air
  ! xl                  ! volume of water in cloud
  ! xliq                ! zclwc/(cloud cover)
  !

  pure subroutine get_meteo( temp, cloud, dens, clwc, cov_ic, &
                        zclwc, xl, xliq )
    
    use Binas, only : Rgas
    ! -- in/out ---------------------------------------------------------------

    real, intent(in)  :: temp         ! temperature (K)
    real, intent(in)  :: cloud        ! cloud cover (0-1)
    real, intent(in)  :: dens         ! air density (kg/m3)
    real, intent(in)  :: clwc         ! cloud liquid water content (kg/kg)
    real, intent(in)  :: cov_ic       ! cloud thickness (vertical fraction cell)
    real, intent(out) :: zclwc        ! cloud liquid water density in m3 water/m3 air
    real, intent(out) :: xl           ! volume of water in cloud
    real, intent(out) :: xliq         ! zclwc/(cloud cover)

    ! -- local  ---------------------------------------------------------------

    !real, parameter :: R=8.31451           ! Gas constant xxx use from Binas ??

    ! -- begin  ---------------------------------------------------------------

    ! Compute zclwc = cloud liquid water density in m3 water/m3 air.
    ! The unit in the code for zlwc to be used is m3 water/m3 air;
    ! initially clwc is in kg/kg, so we have to multiply with 10^-3 assuming a 
    ! water density of 1000 kg/m3 and we have to multiply with dens to 
    ! account for the density of air:

    !ES clwc nemen we vast
    !zclwc = 3.e-4*1.e-3*dens
    zclwc = clwc*1.e-3*dens

    ! Compute xl = volume of water in cloud 
    ! xl is used to convert ppb (nmol of species/mol air) concentrations to 
    ! concentrations in water (mol H+/l) 
    ! hplus = c/xl -> xl = c/hplus -> [xl] = (nmol of species/mol air) / (mol H+/l)
    ! ideal gas law: p*V = n*R*T <=> V = n*R*T/p
    ! p = 10^5 hPa
    ! f1 = conversion factor ppb   = 10^-9 mol/mol air
    ! f2 = conversion factor liter = 10^-3 m3
    ! 1/(f1*f2*p) = 10^7
    ! xxx in sulfwetchem: xl = zclwc/cloud*R*temp*10e7 welke is correct??
    ! xxx in getph      : xl = zclwc*R*temp*10e7       welke is correct??
    xl = zclwc*Rgas*1.e+7*temp 
    
    xl = max(xl, xl_threshold)

    ! Compute xliq = zclwcoud
    xliq = zclwc

  end subroutine get_meteo
  
end module LE_WetDepos_Camx
