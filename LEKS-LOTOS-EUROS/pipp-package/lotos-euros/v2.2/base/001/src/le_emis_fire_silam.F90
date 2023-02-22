!#######################################################################
!
! LE_Emis_Fire_SILAM 
!   routines to read and work with daily fire emissions
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if                                                                                
!
#include "le.inc"
!
!#######################################################################


module LE_Emis_Fire_SILAM

  use GO, only : gol, goPr, goErr
  use GO, only : TDate
!  use binas, only : xm_C, xm_O, xm_H

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Emis_Fire_SILAM_Init, LE_Emis_Fire_SILAM_Done
  public  ::  LE_Emis_Fire_SILAM_Setup
  
  ! --- const ------------------------------------

  ! module name:
  character(len=*), parameter ::  mname = 'LE_Emis_Fire_SILAM'


  ! --- var --------------------------------------

  ! profiles:
  real      ::  hour_in_day_index(1:24)
  real      ::  day_in_week_index(1:7)
  real      ::  month_in_year_index(1:12)

  ! height profile:
  integer, parameter  ::  maxlev = 3
  integer             ::  nheights
  real                ::  emis_height(0:maxlev)
  real                ::  emis_height_frac(maxlev)

  ! horizontal mapping:
  character(len=32)   ::  coordinate_of_values
  real                ::  lon_step, lat_step
  real                ::  lon_start, lat_start
  
  ! interval of current emissions:
  type(TDate)                   ::  t1_fire, t2_fire
  
  ! location of emissions files:
  character(len=1024)           ::  emispath_fire
  
  ! allow mising files ?
  logical                       ::  allow_missing


contains


  ! ===============================================================
  

  subroutine LE_Emis_Fire_SILAM_Init( rcF, rckey, status )
  
    use GO, only : AnyDate
    use GO, only : TrcFile, ReadRc
    use Dims, only : nx, ny !, nz, nspec
    use LE_Emis_Fires_Data, only : ffire_release !, ffire_proba
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------
    
    type(TrcFile), intent(in)       ::  rcF
    character(len=*), intent(in)    ::  rckey
    integer, intent(out)            ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Fire_SILAM_Init'
    
    ! --- local ------------------------------
    
    ! --- begin -------------------------------
    
    ! read settings:
    call ReadRc( rcF, trim(rckey)//'.path', emispath_fire, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read settings:
    call ReadRc( rcF, trim(rckey)//'.allow_missing', allow_missing, status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage for emissions:
    allocate( ffire_release(nx,ny) )
    !allocate( ffire_proba  (nx,ny) )
    !! collect:
    !allocate(  emis_ffire(nx,ny,nz,nspec) )
    
    ! no emissions read yet:
    t1_fire = AnyDate()
    t2_fire = AnyDate()

    ! enable data:
    call LE_Data_Enable( 'oro', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Fire_SILAM_Init
  
  
  ! ***
  

  subroutine LE_Emis_Fire_SILAM_Done( status )

    ! --- in/out ------------------------------
    
    integer, intent(out)  ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Fire_SILAM_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------
    
!    write (gol,'("sorry, routine not finished yet ...")'); call goErr
!    TRACEBACK; status=1; return

    ! ok
    status = 0

  end subroutine LE_Emis_Fire_SILAM_Done
  
  
  ! ***
  

  subroutine LE_Emis_Fire_SILAM_Setup( emis_a, t1, t2, status )

    use GO, only : TDate, NewDate, IncrDate
    use Go, only : operator(<), operator(>), operator(+), operator(-)
    use Go, only : Precisely, rTotal, wrtgol
    use Dims, only : nx, ny, nz, nspec

    ! --- in/out ---------------------------

    real, intent(inout)       ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)   ::  t1, t2
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter :: rname = mname//'/LE_Emis_Fire_SILAM_Setup'
    
    ! --- local ----------------------------
    
    ! --- begin ----------------------------

    ! info ...
    write (gol,'("LE:     forest fires ...")'); call goPr

    ! new interval ?
    if ( (t1 > t2_fire) .or. (t2 > t2_fire) ) then
      ! info ...
      write (gol,'("LE:       read emissions for new day ...")'); call goPr
      ! read new data for current day:
      call read_fire_emis( t1, status )
      IF_NOTOK_RETURN(status=1)
      ! info ..
      call wrtgol( 'LE:         valid for : ', t1_fire, ' - ', t2_fire ); call goPr
      ! check ...
      if ( (t1 < t1_fire) .or. (t2 > t2_fire+IncrDate(sec=1)) ) then
        call wrtgol( 'tried to read forest fire emissions for : ', t1, ' - ', t2 ); call goErr
        call wrtgol( 'but data read seems to be valid for     : ', t1_fire, ' - ', t2_fire ); call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! check: intervals of 1 hour only ..
    if ( (.not. Precisely(t1,1.0,'hour')) .or. (.not. Precisely(t2,1.0,'hour')) &
         .or. (rTotal(t2-t1,'hour') > 1.0) ) then
      call wrtgol( 'Forest fires setup only for hourly intervals [0,1] etc, not for: ',t1,' - ',t2 ); call goErr
      TRACEBACK; status=1; return
    end if

    ! add emissions to 'emis_a' for hour 01..24 :
    call make_fire_emis( emis_a, t1%hour+1, status )
    IF_NOTOK_RETURN(status=1)
    
!    write (gol,'("sorry, routine not finished yet ...")'); call goErr
!    TRACEBACK; status=1; return

    ! ok
    status = 0
       
  end subroutine LE_Emis_Fire_SILAM_Setup


  ! ***
  

  subroutine read_fire_emis( t, status )

    use GO, only : TDate
    use GO, only : AnyDate, IsAnyDate, DayNumber, NewDate
    use GO, only : goGetFu
    use GO, only : goSplitLine, goReplace   
    use Grid, only : In_Domain
    use LE_Grid, only : ugg
    use LE_Emis_Fires_Data, only : ffire_release  !, ffire_proba

    ! --- in/out ---------------------------------
    
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/read_fire_emis'
    
    ! --- local ----------------------------------
    
    logical               ::  exist
    integer               ::  fu
    integer               ::  irec
    character(len=1024)   ::  fname
    character(len=1024)   ::  line
    character(len=64)     ::  key
    character(len=1024)   ::  value
    real                  ::  lon, lat, em
    real                  ::  ilonr, ilatr
    logical               ::  inflag
    integer               ::  ix, iy
    real                  ::  h1, h2, frac
    character(len=64)     ::  unit
    integer               ::  ipar
    integer               ::  year, month, day, hour, minu
    real                  ::  rsec, r1, r2, r3
    
    ! --- begin ----------------------------------

    ! clear current fields:
    ffire_release = 0.0
    !ffire_proba   = 0.0

    ! level counter:
    nheights = 0
    
    ! par_str counter:
    ipar = 0
    t1_fire = AnyDate()
    t2_fire = AnyDate()

    ! write filename:
    fname = trim(emispath_fire)
    ! replace some values:
    !  o %{year4}, %{month2}, %{day2}   : year, month and day in 2 or 4 digits
    call goReplace( fname, '%{year4}' , '(i4.4)', t%year , status )
    IF_NOTOK_RETURN(status=1)
    call goReplace( fname, '%{month2}', '(i2.2)', t%month, status )
    IF_NOTOK_RETURN(status=1)
    call goReplace( fname, '%{day2}'  , '(i2.2)', t%day  , status )
    IF_NOTOK_RETURN(status=1)
    !!  o %{iday3}   : daynumber within year, e.g. 001 .. 365
    !call goReplace( fname, '%{iday3}'  , '(i3.3)', DayNumber(t), status )
    !IF_NOTOK_RETURN(status=1)
    
    ! check ...
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      ! warning only ?
      if ( allow_missing ) then
        write (gol,'("WARNING - forest fire emission not found : ",a)') trim(fname); call goPr
        status=0; return
      else
        write (gol,'("forest fire emission not found : ",a)') trim(fname); call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! free file unit:
    call goGetFu( fu, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open:
    open( unit=fu, file=trim(fname), status='old', form='formatted', iostat=status )
    if ( status/=0 ) then
      write (gol,'("opening file : ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! record counter:
    irec = 0
    
    ! loop over records:
    do
    
      ! increase counter:
      irec = irec + 1

      ! read line:
      read (fu,'(a)',iostat=status) line
      ! eof ?
      if (status < 0) exit
      ! error ..
      if (status/=0) then
        write (gol,'("reading record ",i6," from file : ",a)') trim(fname); call goErr
        TRACEBACK; status=1; return
      end if
       
      ! skip empty lines:
      if ( len_trim(line) == 0 ) cycle
      ! skip comment:
      if ( line(1:1) == '#' ) cycle


      ! split in 'key = value'
      call goSplitLine( line, key, '=', value, status )
      IF_NOTOK_RETURN(status=1)
    
      
      ! take appropritate action:
      select case ( trim(key) )
      
        case ( 'AREA_SOURCE_2', 'END_AREA_SOURCE_2' )
          ! skip ...
      
        case ( 'source_name', 'source_sector_name' )
          ! skip ...

        case ( 'release_rate_unit' )
          if ( trim(value) /= 'kg/sec' ) then
            write (gol,'("unsupported ",a," : ",a)') trim(key), trim(value); call goErr
            TRACEBACK; status=1; return
          end if
          
        case ( 'dx'        ) ; read (value,*) lon_step
        case ( 'dy'        ) ; read (value,*) lat_step
        case ( 'lon_start' ) ; read (value,*) lon_start
        case ( 'lat_start' ) ; read (value,*) lat_start
        
        case ( 'nx', 'ny', 'grid_type', 'lon_s_pole', 'lat_s_pole', &
                'lon_pole_stretch', 'lat_pole_stretch', &
                'resol_flag', 'ifReduced', 'earth_flag', 'wind_component', 'reduced_nbr_str' )
          ! skip ...
          
        case ( 'cocktail_composition', 'size_class_split' )
          if ( trim(value) /= 'COMMON' ) then
            write (gol,'("unsupported ",a," : ",a)') trim(key), trim(value); call goErr
            TRACEBACK; status=1; return
          end if
          
        case ( 'emitted_substance' )
          if ( trim(value) /= 'WHOLE_COCKTAIL' ) then
            write (gol,'("unsupported ",a," : ",a)') trim(key), trim(value); call goErr
            TRACEBACK; status=1; return
          end if
        
        case ( 'emitted_size_mode_nbr' )
          ! skip ...
          
        case ( 'hour_in_day_index' )
          read (value,*) hour_in_day_index
          
        case ( 'day_in_week_index' )
          read (value,*) day_in_week_index
          if ( any( day_in_week_index /= 1.0 ) ) then
            write (gol,*) 'unsupported ',trim(key),' : ', trim(key), day_in_week_index; call goErr
            TRACEBACK; status=1; return
          end if
          
        case ( 'month_in_year_index ' )
          read (value,*) month_in_year_index 
          if ( any( month_in_year_index  /= 1.0 ) ) then
            write (gol,*) 'unsupported ',trim(key),' : ', trim(key), month_in_year_index ; call goErr
            TRACEBACK; status=1; return
          end if
          
        case ( 'vertical_unit' )
          if ( trim(value) /= 'm' ) then
            write (gol,'("unsupported ",a," : ",a)') trim(key), trim(value); call goErr
            TRACEBACK; status=1; return
          end if
          
        case ( 'vert_level' )

          ! next level:
          nheights = nheights + 1
          ! check ...
          if ( nheights > maxlev ) then
            write (gol,'("more levels defined in file than maximum : ",i3)') maxlev; call goErr
            TRACEBACK; status=1; return
          end if
          ! extract values:
          read (value,*) unit, h1, h2, frac
          ! check ...
          if ( trim(unit) /= 'HEIGHT_FROM_SURF' ) then
            write (gol,'("unsupported ",a," unit : ",a)') trim(key), trim(unit); call goErr
            TRACEBACK; status=1; return
          end if
          ! store:
          if ( nheights == 1 ) emis_height(0) = 0.0
          emis_height(nheights) = h2
          emis_height_frac(nheights) = frac

  
        case ( 'vertical_distribution' )
          if ( trim(value) /= 'MULTI_LEVEL_FIXED' ) then
            write (gol,'("unsupported ",a," : ",a)') trim(key), trim(value); call goErr
            TRACEBACK; status=1; return
          end if
        
        case ( 'par_str' )
          ! increase counter:
          ipar = ipar + 1
          ! read value:
          read (value,*) year, month, day, hour, minu, rsec, r1, r2, r3, unit
          ! store time:
          if ( ipar == 1 ) then
            t1_fire = NewDate( year=year, month=month, day=day, hour=hour, min=minu, sec=nint(rsec) )
          else if ( ipar == 2 ) then
            t2_fire = NewDate( year=year, month=month, day=day, hour=hour, min=minu, sec=nint(rsec) )
          else
            write (gol,'("not more than 2 `par_str` lines expected")'); call goErr
            TRACEBACK; status=1; return
          end if

          ! check ...
          select case ( trim(unit) )
            case ( 'FIRE_PM_COCKTAIL' )
              ! I guess these are the same ...
            case default
              write (gol,'("unsupported ",a," unit : ",a)') trim(key), trim(unit); call goErr
              TRACEBACK; status=1; return
          end select
          
        case ( 'coordinate_of_values' )
        
          ! store:
          coordinate_of_values = trim(value)
          
        case ( 'val' )
      

          ! fill in correct cell:
          select case ( trim(coordinate_of_values) )

            ! lon lat is provided:
            case ( 'GEOGRAPHICAL' )

              ! extract values:
              read (value,*) lon, lat, em

            case ( 'GRID_INDICES' )
            
              ! extract indices in emission grid:
              read (value,*) ilonr, ilatr, em
              
              ! convert to lon/lat:
              lon = lon_start + lon_step * (ilonr-0.5)
              lat = lat_start + lat_step * (ilatr-0.5)

            case default

              write (gol,'("unsupported coordinate_of_values : ",a)') trim(coordinate_of_values); call goErr
              TRACEBACK; status=1; return

          end select
      
          ! check location:
          call ugg%InDomain(lon,lat,inflag,status)
          IF_NOTOK_RETURN(status=1)

          ! in grid ?
          if ( inflag ) then
            ! find cell location:
            ! ... lon, lat, ix, iy, ...
            ! dummy:
            ix = -999
            iy = -999
            write (gol,'("not implemented yet")'); call goErr
            TRACEBACK; status=1; return
            ! exclude some water pixels near cost, which have extreme high emissions ...
            !if ( inflag .and. (lu_lotos(ix,iy,7) <= 0.1) ) then
              !! testing ...
              !print *, '  c2 add to cell ', ix, iy
            ! add contribution:
            ffire_release(ix,iy) = ffire_release(ix,iy) + em
          end if 
          
        case default

          write (gol,'("unsupported key `",a,"` in forest fire emission file ...")') trim(key); call goErr
          TRACEBACK; status=1; return

      end select

    end do  ! records
    
    ! close:
    close( unit=fu, iostat=status )
    if ( status/=0 ) then
      write (gol,'("closing file : ",a)') trim(fname); call goErr
      TRACEBACK; status=0; return
    end if
    
    ! check ...
    if ( IsAnyDate(t1_fire) .or. IsAnyDate(t2_fire) ) then
      write (gol,'("could not guess time interval for forest fires; need 2 `par_str` lines")'); call goErr
      write (gol,'("  file : ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine read_fire_emis

! ***
  
  subroutine make_fire_emis( emis_a, hour, status )

    use Num, only : IntervalSum
    use indices
    use dims, only : nx, ny, nz, nspec
    use LE_Data, only : LE_Data_GetPointer
    use LE_Emis_Fires_Data
!  
!    use LE_Emis_Fires_Data, only : emis_ffire
!    use LE_Emis_Fires_Data, only : ffire_release
!    use LE_Emis_Fires_Data, only : rEF_CO, rEF_CH4, rEF_NOx, rEF_NO, rEF_NO2, rEF_VOC, rEF_SO2, rEF_NH3, rEF_BC, rEF_PM25
!    use LE_Emis_Fires_Data, only : xfire_OLE, xfire_PAR, xfire_TOL, xfire_XYL, xfire_FORM, xfire_ALD, xfire_ETH, xfire_UNR
!    
!      !use dims, only : volume
      use file_nc

    ! --- in/out ---------------------------------
    
    real, intent(inout)         ::  emis_a(nx,ny,nz,nspec)
    integer, intent(in)         ::  hour  ! 1..24
    integer, intent(out)        ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/make_fire_emis'

    ! --- local ----------------------------------
    
    integer   ::  ix, iy, iz
    integer   ::  ilast
    real      ::  height_frac(nz)
    real      ::  hh(0:nz)
    real      ::  tfrac
    
    ! meteo data:
    real, pointer        ::  oro(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  h_m(:,:,:)   ! (lon,lat,1)

    ! --- begin ----------------------------------
    
    ! point to meteo data:
    call LE_Data_GetPointer( 'oro', oro, status,check_units='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status,check_units='m' )
    IF_NOTOK_RETURN(status=1)
    
    ! init current forest fire emissions to zero:
    
        !call nc_dump( 'ffire_release.nc', ffire_release, 'ffire_release', (/'x','y','z','s'/), status )
        !stop 'break after dump'

    ! time fraction:
    tfrac = hour_in_day_index(hour)  !  ~ 1.0

    ! loop over grid cells:
    do iy = 1, ny
      do ix = 1, nx
 
        ! skip if no fire ...
        if ( ffire_release(ix,iy) <= 0.0 ) cycle

        ! model half level heights above surface:
        hh(0:nz) = (/ 0.0, h_m(ix,iy,:) /)  ! m
        
        ! loop over model layer:
        ilast = 0
        do iz = 1, nz
          ! bottom of layer above emission top ?
          if ( hh(iz-1) > emis_height(nheights) ) then
            ! no emissions here ...
            height_frac(iz) = 0.0
          else
            ! sum fractions of fractions; maximum model height is emission top:
            call IntervalSum( emis_height(0:nheights), emis_height_frac(1:nheights), &
                                hh(iz-1), min(hh(iz),emis_height(nheights)), &
                                height_frac(iz), ilast, status )
            IF_NOTOK_RETURN(status=1)
          end if
        end do
        ! over mountain area's, the model top (above the surface) might be lower than the emission top;
        ! add the remaining emissions to the hightest model layer for this case:
        if ( hh(nz) < emis_height(nheights) ) then
          height_frac(nz) = height_frac(nz) + (1.0 - sum(height_frac(1:nz)))
        end if

!test...        
!print *, 'bbb1 ', ix, iy, height_frac

        ! check the fractions:
        if ( abs( sum(height_frac) - 1.0 ) > 1.0e-2 ) then
          write (gol,'("height fractions should sum to 1.0:")'); call goErr
          write (gol,*) '  emis heights            : ', emis_height(0:nheights); call goErr
          write (gol,*) '  emis emis_height_frac   : ', emis_height_frac(1:nheights); call goErr
          write (gol,*) '  model cell              : ', ix, iy; call goErr
          write (gol,*) '  model heights           : ', oro(ix,iy,1)+hh; call goErr
          write (gol,*) '  model fractions         : ', height_frac; call goErr
          TRACEBACK; status=1; return
        end if

        ! loop over model levels:
        do iz = 1,nz

          ! primary aerosols
          ! unitconversion from kg/sec to ug/min
          !    ug/min                                                     [0-1]          ~1    (kg cocktail)/sec (kg tracer)/(kg cocktail) ug/kg sec/min
   
          if ( i_ppm_f > 0 ) emis_a(ix,iy,iz,i_ppm_f) = emis_a(ix,iy,iz,i_ppm_f) + height_frac(iz) * tfrac * ffire_release(ix,iy) *rEF_PPM_F* 1e9 * 60
          if ( i_ppm_c > 0 ) emis_a(ix,iy,iz,i_ppm_c) = emis_a(ix,iy,iz,i_ppm_c) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_PPM_C*1e9 * 60
!          if ( i_ppm_f > 0 ) emis_a(ix,iy,iz,i_ppm_f) = emis_a(ix,iy,iz,i_ppm_f) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_PM25 * 1e9 * 60
          if ( i_ec_f > 0 ) emis_a(ix,iy,iz,i_ec_f ) = emis_a(ix,iy,iz,i_ec_f ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_EC_F   * 1e9 * 60           
          if ( i_pom_f > 0 ) emis_a(ix,iy,iz,i_pom_f ) = emis_a(ix,iy,iz,i_pom_f ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_POM_F   * 1e9 * 60 
!          ! reactive gasses!
!          ! GASES ; unit conversion from kg/sec to mol/min
          if (i_co  > 0) emis_a(ix,iy,iz,i_co  ) = emis_a(ix,iy,iz,i_co  ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_CO  / specmolm(i_CO ) * 60.0
          if (i_ch4 > 0) emis_a(ix,iy,iz,i_ch4 ) = emis_a(ix,iy,iz,i_ch4 ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_CH4 / specmolm(i_CH4) * 60.0
          if (i_no  > 0) emis_a(ix,iy,iz,i_no  ) = emis_a(ix,iy,iz,i_no  ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_NO  / specmolm(i_NO ) * 60.0
          if (i_no2 > 0) emis_a(ix,iy,iz,i_no2 ) = emis_a(ix,iy,iz,i_no2 ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_NO2 / specmolm(i_NO ) * 60.0
          if (i_so2 > 0) emis_a(ix,iy,iz,i_so2 ) = emis_a(ix,iy,iz,i_so2 ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_SO2 / specmolm(i_SO2) * 60.0
          if (i_nh3 > 0) emis_a(ix,iy,iz,i_nh3 ) = emis_a(ix,iy,iz,i_nh3 ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_NH3 / specmolm(i_NH3) * 60.0
!          ! partioning VOC ; unit conversion from (kg voc)/sec to (mol component)/min
          if (i_eth > 0) emis_a(ix,iy,iz,i_eth ) = emis_a(ix,iy,iz,i_eth ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_VOC * xfire_ETH  * 60.0
          if (i_ole > 0) emis_a(ix,iy,iz,i_ole ) = emis_a(ix,iy,iz,i_ole ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_VOC * xfire_OLE  * 60.0
          if (i_par > 0) emis_a(ix,iy,iz,i_par ) = emis_a(ix,iy,iz,i_par ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_VOC * xfire_PAR  * 60.0
          if (i_tol > 0) emis_a(ix,iy,iz,i_tol ) = emis_a(ix,iy,iz,i_tol ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_VOC * xfire_TOL  * 60.0
          if (i_xyl > 0) emis_a(ix,iy,iz,i_xyl ) = emis_a(ix,iy,iz,i_xyl ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_VOC * xfire_XYL  * 60.0
          if (i_form> 0) emis_a(ix,iy,iz,i_form) = emis_a(ix,iy,iz,i_form) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_VOC * xfire_FORM * 60.0
          if (i_ald > 0) emis_a(ix,iy,iz,i_ald ) = emis_a(ix,iy,iz,i_ald ) + height_frac(iz) * tfrac * ffire_release(ix,iy) * rEF_VOC * xfire_ALD  * 60.0
        end do  ! z


      end do  ! x
    end do ! y
    
!        !! dump ...    
!        !call nc_dump( 'volume.nc', volume, 'volume', (/'x','y','z'/), status )
!        call nc_dump( 'emis_ffire.nc', emis_ffire, 'emis_ffire', (/'x','y','z','s'/), status )
!        stop 'break after dump'
    
    ! ok
    status = 0

  end subroutine make_fire_emis


end module LE_Emis_Fire_SILAM

