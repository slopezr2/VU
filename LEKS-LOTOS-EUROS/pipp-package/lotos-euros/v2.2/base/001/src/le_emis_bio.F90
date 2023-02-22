!###############################################################################
!
! LE_Emis - LOTOS-EUROS emission routines
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

module LE_Emis_Bio

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Emis_Bio_Init, LE_Emis_Bio_Done
  public  ::  LE_Emis_Bio_VOC
  

  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Bio'
  
  !! fraction of the total monoterpene emissions
  !! that are in the form of a-pinene and limone
  !real, parameter :: api_fac=0.21, lim_fac = 0.23

  integer, parameter            ::  nbio_emis = 3
  character(len=10), parameter  ::  bio_emis_names(nbio_emis) = &
                      (/ 'ISOP      ', &
                         'Terp_Light', &
                         'Terp_Store' /)
  integer, parameter            ::  i_iso_emis        = 1
  integer, parameter            ::  i_terp_light_emis = 2
  integer, parameter            ::  i_terp_store_emis = 3

  ! --- types ------------------------------------

  type T_vegetationtype
    character(len=64) ::  name
    character(len=1)  ::  type
    real              ::  biomass
    character(len=64) ::  time_profile
    character(len=64) ::  lat_profile
    character(len=64) ::  source
    real              ::  coeff(nbio_emis)
  end type T_vegetationtype

  ! --- var --------------------------------------

  ! tree class names etc:
  type(T_vegetationtype), allocatable :: vegetation(:)

contains


  ! ==============================================
  
  
  subroutine LE_Emis_Bio_Init( rcF, status )
    
    use GO, only : TrcFile, ReadRc
    use LE_Landuse_Data, only : nveg
    use LE_Data, only : LE_Data_Enable
    
    ! --- in/out ---------------------------------

    type(TrcFile), intent(in) ::  rcF
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Bio_Init'
    
    ! --- local ----------------------------------
    
    character(len=512)              ::  query
    
    
    ! --- begin ----------------------------------
    
    !
    allocate( vegetation(nveg) )
    
    ! initialize tree names etc:
    call ReadRc( rcF, 'le.landuse.vegetation_data', query, status )
    IF_NOTOK_RETURN(status=1)
       
    call read_vegetation_data(query, status)
    IF_NOTOK_RETURN(status=1)

    ! enable data:
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'ssrd', status )
    IF_NOTOK_RETURN(status=1)    
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lat', status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_Bio_Init
  
  
  ! ***
  
  
  subroutine LE_Emis_Bio_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Bio_Done'
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_Bio_Done
  
  
  ! ***


  ! this routine adds biogenic VOC emissions to the emissions array
  ! input: the emission array

  subroutine LE_Emis_Bio_VOC( emis, iday, status )

    use dims, only : runF
    use dims, only : nx, ny, nz, nspec
    use LE_LandUse_Data, only : lu_fracs
    use LE_LandUse_Data, only : nveg, veg_fracs
    use LE_Landuse_Data, only : ilu_arable
    use indices

    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer


    ! --- in/out ---------------------------------
    
    real, intent(inout)       ::  emis(nx,ny,nz,nspec)
    integer, intent(out)      ::  status
    integer                   ::  iday
    
    ! --- const ----------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Bio_VOC'
    
    ! factors for temperature dependency terpene emissions
    real, parameter :: tref = 303.15
    real, parameter :: beta = 0.09
    
    ! --- local ---------------------------------
    
    real    ::  tk
    real    ::  isoprene_emis, monoterp_emis
    real    ::  vegetation_area(nveg)
    real    ::  par
    real    ::  Dmass, T_fac_iso, T_fac_mt, e_mt
    integer ::  ix, iy, iveg
    integer ::  ilu
    real    ::  lat
    real    ::  lai, sai
    real    ::  delta_emis

    ! meteo data:
    real, pointer               ::  temp(:,:,:)   ! (lon,lat,1)
    real, pointer               ::  radd(:,:,:)   ! (lon,lat,1)    
    real, pointer               ::  area(:,:,:)   ! (lon,lat,1)
    real, pointer               ::  lats(:,:,:)   ! (lon,lat,1)

    ! --- begin ----------------------------------

    ! point to meteo data:
    call LE_Data_GetPointer( 't', temp, status, check_units='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'ssrd', radd, status, check_units='J/m2/s' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'area', area, status, check_units ='m2' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lat', lats, status, check_units ='degrees_north' )
    IF_NOTOK_RETURN(status=1)
    

    ! monoterpene and isoprene emissions
    do ix=1,nx
      do iy=1,ny

        ! temperature (K):
        tk = temp(ix,iy,1)

        ! approximation for par =50% of GLOB-RAD and 1W/m2 = 4.6 mueE/m2/s  
        par = 0.5 * 4.6*radd(ix,iy,1)        

        ! calculate the temperature and light dependent scaling factor
        ! for isoprene (not species dependent)
        call get_temp_factor_isoprene(tk, par, T_fac_iso)
        ! idem for monoterpene:
        call get_temp_factor_monoterp(tk, tref, beta, T_fac_mt)
        
        ! leaf-area-index for crops depends on day and latitude:
        ilu = ilu_arable
        lat = lats(ix,iy,1)
        call slai_emberson(iday, lat, ilu, lai, sai)

        ! re-initialize emissions to zero
        isoprene_emis=0.0
        monoterp_emis=0.0

        ! tree coverage fractions:
        vegetation_area(:) = veg_fracs(ix,iy,:)   ! fractional                        

        ! loop over vegetation types:
        do iveg = 1, nveg                                
           
          ! set Dmass (vegetation dry mass in g/m2) ;
          if ( trim(vegetation(iveg)%time_profile) == 'const' .and. &
               trim(vegetation(iveg)%lat_profile ) == 'const' ) then
            
            ! no special treatment of Dmass, take number from table:
            Dmass = vegetation(iveg)%biomass  
            
          else if ( trim(vegetation(iveg)%time_profile) == '0.5*(251+1067.5*dlai)') then

            ! biomass number depends on the lai index
            if ( lai > 0.0 ) then
              ! from old code permanent crops + arable wheat  0.5 * (251 + 1067.5)
              Dmass = 659.25
            else 
              ! from old code permanent crops
              Dmass = 251
            end if

          else if ( trim(vegetation(iveg)%lat_profile) == 'step_1(lat)' ) then
            
            ! mass depends on latitude (method 1):
            if ( lat < 55.0 ) then
              Dmass = 1600.
            else if (lat >= 55.0 .and. lat <= 60.0 ) then
              Dmass = 1400.
            else 
              Dmass = 800.
            endif

          else if ( trim(vegetation(iveg)%lat_profile) == 'step_2(lat)' ) then
            
            ! mass depends on latitude (method 2):
            if (lat > 60.0) then
              Dmass = 500.0
            else
              Dmass = 700.0
            end if

          else 
          
            ! unknown ...
            write (gol,'("unknown time or lat profile")'); call goErr
            write (gol,'("time profile : ",a)') trim(vegetation(iveg)%time_profile) ; call goErr
            write (gol,'("lat_profile : ",a)') trim(vegetation(iveg)%lat_profile)   ; call goErr
            write (gol,'("vegetation_type : ", i4)') iveg ; call goErr
            TRACEBACK; status=1; return            

          end if
          
          ! no vegetation ? then skip:
          if ( Dmass < 0.0 ) cycle
          
          ! add contribution to isoprene emissions:
          if ( vegetation(iveg)%coeff(i_iso_emis) > 0.0 ) then
            isoprene_emis=  isoprene_emis + Dmass*vegetation(iveg)%coeff(i_iso_emis)*vegetation_area(iveg)*T_fac_iso
          endif
          
          ! now calculate the temperature and light dependent scaling factors
          ! for monoterpenes which is species dependent
          if ( vegetation(iveg)%coeff(i_terp_store_emis) < 0.0 ) then 
            ! do exceptions that follow isoprene factor              
            e_mt = vegetation(iveg)%coeff(i_terp_light_emis)
            if(e_mt.lt.0.0) e_mt = 0.0  ! takes care of species that do not emit MT
            monoterp_emis= monoterp_emis + Dmass*e_mt*vegetation_area(iveg)*T_fac_iso
          else
            e_mt = vegetation(iveg)%coeff(i_terp_store_emis)
            monoterp_emis= monoterp_emis + Dmass*e_mt*vegetation_area(iveg)*T_fac_mt
          endif

        enddo              

        ! add the emissions (ug/m2/hr) to the emission array in mol/min!
        if ( i_iso > 0 ) then
          ! mol/min     m2       ug/m2/hr   / (ug/g) / (min/hr) / (   (kg/mol       g/kg )
          delta_emis = area(ix,iy,1) * isoprene_emis /  1e6   /  60.0    / (specmolm(i_iso) * 1e3 )
          emis(ix,iy,1,i_iso ) = emis(ix,iy,1,i_iso ) + delta_emis
        endif 
        if ( i_terp > 0 ) then  
          ! mol/min     m2       ug/m2/hr   / (ug/g) / (min/hr) / (   (kg/mol       g/kg ) 
          delta_emis = area(ix,iy,1) * monoterp_emis/1e6/60./(specmolm(i_terp)*1e3)
          emis(ix,iy,1,i_terp) = emis(ix,iy,1,i_terp) + delta_emis
          !emis(i,j,1,i_api1) = emis(i,j,1,i_api1) + area*(0.5*api_fac)*monoterp_emis/1e6/60./184.0
          !emis(i,j,1,i_api2) = emis(i,j,1,i_api2) + area*(0.5*api_fac)*monoterp_emis/1e6/60./184.0
          !emis(i,j,1,i_lim1) = emis(i,j,1,i_lim1) + area*(0.5*lim_fac)*monoterp_emis/1e6/60./184.0
          !emis(i,j,1,i_lim2) = emis(i,j,1,i_lim2) + area*(0.5*lim_fac)*monoterp_emis/1e6/60./184.0
        endif                  

      enddo  ! ix
    enddo  ! iy

    ! ok
    status = 0

  end subroutine LE_Emis_Bio_VOC


  ! ***
  

  subroutine read_vegetation_data( query, status )

    use GO, only :  goGetFU
    use GO, only :  goVarValue, goSplitString, goMatchValue
    use LE_Landuse_Data, only : nveg
    ! --- in/out ---------------------------------

    character(len=*), intent(in)    ::  query
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   :: rname = mname//'/read_vegetation_data'

    ! maximum number of columns in csv file:
    integer, parameter    ::  maxcol = 10

    ! --- local ----------------------------------

    integer                 ::  iveg
    
    character(len=512)      ::  filename
    character(len=1)        ::  comment
    character(len=1)        ::  sep
    character(len=32)       ::  fileformat
    logical                 ::  exist

    integer                 ::  ifield_nr
    integer                 ::  ifield_name
    integer                 ::  ifield_type
    integer                 ::  ifield_biomass
    integer                 ::  ifield_time_profile
    integer                 ::  ifield_lat_profile
    integer                 ::  ifield_source
    integer, allocatable    ::  icomps(:)   ! (nfield)

    integer                 ::  fu
    integer                 ::  iline
    character(len=1024)     ::  line
    integer                 ::  nheader
    character(len=64)       ::  headers(maxcol)
    character(len=64)       ::  header
    integer                 ::  nfield
    integer                 ::  ifield
    character(len=64)       ::  fields(maxcol)
    character(len=64)       ::  field
    integer                 ::  icomp
    
    integer                 ::  nr
    character(len=64)       ::  name
    character(len=1)        ::  type
    real                    ::  biomass
    character(len=64)       ::  time_profile
    character(len=64)       ::  lat_profile
    character(len=64)       ::  source
    real                    ::  coeff

    ! --- begin ----------------------------------

    ! initialize vegetation type
    do iveg = 1,nveg
      vegetation(iveg)%biomass      = -999.0
      vegetation(iveg)%coeff        = -999.0
      vegetation(iveg)%time_profile = 'const'
      vegetation(iveg)%lat_profile  = 'const'
    enddo
! info ...
    !print *,'scanning tree data'
    ! empty ?
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( len_trim(query) == 0 ) then   ! empty
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else  ! read table from file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! extract filename:
      call goVarValue( trim(query), ';', 'file', '=', filename, status )
      IF_NOTOK_RETURN(status=1)
      ! seperation character:
      sep = ';'
      call goVarValue( trim(query), ';', 'sep', '=', sep, status )
      IF_ERROR_RETURN(status=1)
      ! comment character:
      comment = '#'
      call goVarValue( trim(query), ';', 'comment', '=', comment, status )
      IF_ERROR_RETURN(status=1)
      ! format description:
      fileformat = 'csv'
      call goVarValue( trim(query), ';', 'format', '=', fileformat, status )
      IF_ERROR_RETURN(status=1)
      
      ! file should be present:
      inquire( file=trim(filename), exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! switch per format:
      select case ( trim(fileformat) )
        
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'csv' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
          ! free file unit:
          call goGetFU( fu, status )
          IF_NOTOK_RETURN(status=1)

          ! open file:
          open( fu, file=trim(filename),status = 'old', form='formatted', iostat=status )
          if (status/=0) then
            write (gol,'("opening tree data file : ",a)') trim(filename); call goErr
            TRACEBACK; status=1; return
          end if

          iline = 0

          ! read header line after first comment:
          do
            ! read line:
            read (fu,'(a)',iostat=status) line
            if (status/=0) then
              write (gol,'("reading header line from file : ",a)') trim(filename); call goErr
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
          
          ! set variables from header; first dummy values:
          ifield_nr = -1
          ifield_name = -1
          ifield_type = -1
          ifield_Biomass = -1
          ifield_time_profile = -1
          ifield_lat_profile = -1
          ifield_source = -1
          
          ! storage to map from field index to component index:
          allocate( icomps(nheader) ) ; icomps = -1
          ! loop over fields:
          do ifield = 1, nheader
            ! current:
            header = headers(ifield) 
            ! which column ?
            select case ( trim(header) )
              !~ index nr:
              case ( 'nr', 'Nr', 'NR' )
                ! store field index:
                ifield_nr = ifield
              !~ vegetation name:
              case ( 'name', 'vegetation_name' )
                ! store field index:
                ifield_name = ifield
              !~ vegatation type
              case ( 'type', 'Type' )
                ! store field index:
                ifield_type = ifield
              !~ biomass number
              case ( 'biomass', 'Biomass' )
                ! store field index:
                ifield_biomass = ifield
              !~ time profile
              case ( 'time_prof', 'Time_prof' )
                ! store field index:
                ifield_time_profile = ifield
              !~ lat profile
              case ( 'lat_prof', 'Lat_prof' )
                ! store field index:
                ifield_lat_profile = ifield
              !~ source
              case ( 'source', 'Source' )
                ! store field index:
                ifield_source = ifield
              !~ emitted component coefficient...
              case default
                ! search:
                call goMatchValue( trim(header), bio_emis_names, icomp, status )
                IF_NOTOK_RETURN(status=1)
                ! assign to field index:
                icomps(ifield) = icomp
              !~
            end select  ! header
          end do  ! fields
        
          ! check indices:
          if ( ifield_nr < 0 ) then
            write (gol,'("number index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return
          end if
          
          if ( ifield_name < 0 ) then
            write (gol,'("name index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return            
          end if
          
          if ( ifield_type < 0 ) then
            write (gol,'("type index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return            
          end if
          
          if ( ifield_biomass < 0 ) then
            write (gol,'("biomass index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return            
          end if
          
          if ( ifield_time_profile < 0 ) then
            write (gol,'("time profile index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return            
          end if

          if ( ifield_lat_profile < 0 ) then
            write (gol,'("lat profile index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return            
          end if
          
          if ( ifield_source < 0 ) then
            write (gol,'("source index not defined yet while processing column `",a,"` of line:")') trim(header); call goErr
            write (gol,'("  ",a)') trim(line); call goErr
            TRACEBACK; status=1; return            
          end if
          
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
              write (gol,'("reading line ",i6," from file : ",a)') iline, trim(filename); call goErr
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
              write (gol,'("differs from number of headers ",i2," in file ",a)') nheader, trim(filename); call goErr
              TRACEBACK; status=1; return
            end if
            ! number
            read (fields(ifield_nr),*,iostat=status) nr
            if (status/=0) then
              write (gol,'("reading number index from `",a,"`")') field; call goErr
              TRACEBACK; status=1; return
            end if
            ! name
            read (fields(ifield_name),*,iostat=status) name
            if (status/=0) then
              write (gol,'("reading vegetation name from `",a,"`")') field; call goErr
              TRACEBACK; status=1; return
            end if
            ! type
            read (fields(ifield_type),*,iostat=status) type
            if (status/=0) then
              write (gol,'("reading vegetation type from `",a,"`")') field; call goErr
              TRACEBACK; status=1; return
            end if
            ! biomass
            read (fields(ifield_biomass),*,iostat=status) biomass
            if (status/=0) then
              write (gol,'("reading biomass value from `",a,"`")') field; call goErr
              TRACEBACK; status=1; return
            end if
            ! time_profile
            read (fields(ifield_time_profile),*,iostat=status) time_profile
            if (status/=0) then
              write (gol,'("reading biomass profile from `",a,"`")') field; call goErr
              TRACEBACK; status=1; return
            end if
            ! lat_profile
            read (fields(ifield_lat_profile),*,iostat=status) lat_profile
            if (status/=0) then
              write (gol,'("reading biomass profile from `",a,"`")') field; call goErr
              TRACEBACK; status=1; return
            end if
            ! source
            read (fields(ifield_source),*,iostat=status) source
            if (status/=0) then
              write (gol,'("reading source literature from `",a,"`")') field; call goErr
              TRACEBACK; status=1; return
            end if
            
            ! check type 
            if ( .not. any(type == (/'C','D','H','N'/) ) ) then
              write (gol,'("unsupported vegetation type `",a,"` for vegetation number: ",i4)') type, nr; call goErr
              TRACEBACK; status=1; return
            end if
            
            ! fill in type
            vegetation(nr)%name         = name
            vegetation(nr)%type         = type
            vegetation(nr)%biomass      = biomass
            vegetation(nr)%time_profile = time_profile
            vegetation(nr)%lat_profile  = lat_profile
            vegetation(nr)%source       = source
            
            ! loop over fields:
            do ifield = 1, nfield

              ! current component:
              icomp = icomps(ifield)
              ! no component ? then skip:
              if ( icomp < 0 ) cycle

              ! current field characters:
              field = fields(ifield)            
              ! read factor:
              read (field,*,iostat=status) coeff
              if (status/=0) then
                write (gol,'("reading coefficient from `",a,"`")') field; call goErr
                TRACEBACK; status=1; return
              end if

              vegetation(nr)%coeff(icomp) = coeff
            enddo
            
          end do ! lines
          
          ! clear:
          deallocate( icomps )
     
          ! close
          close( fu, iostat=status )
          if (status/=0) then
            write (gol,'("closing file : ",a)') trim(filename); call goErr
            TRACEBACK; status=1; return
          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          write (gol,'("unsupported format : ",a)') trim(fileformat); call goErr
          TRACEBACK; status=1; return
          
       end select         

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! ok
    status = 0

  end subroutine read_vegetation_data

  
  ! ***


  subroutine get_temp_factor_isoprene(tk, par, iso_factor)
    !-----------------------------------------------------------------------------
    !  compute light dependence for isoprene emissions
    !  formula's (3,4) of David Simpson 1995
    !  Biogenic emissions estimates and uncertainties
    !  J. Geophys. Res. 100, D11
    !      input:
    !      TEMP    Temperature in K
    !      PAR     photosyntetically active ration in micromole photons / (m2 sec )
    !              microE/(m2 sec)
    !-----------------------------------------------------------------------------

    real, parameter :: alpha=0.0027
    real, parameter :: cl1=1.066
    real, parameter :: Rgas = 8.314   ! J/ (K mole)
    real, parameter :: Ts=303.0       ! K standard temperature
    real, parameter :: Tm=314.0       ! K empirical
    real, parameter :: ct1=95000.0    ! J/mole
    real, parameter :: ct2=230000.0   ! J/mole


    real :: tk, par, ct, cl
    real :: arg1, arg2, fac, iso_factor

    fac  = Rgas*tk*ts
    arg1 = ct1*(tk-ts)/fac
    arg2 = ct2*(tk-tm)/fac
    ct = exp(arg1)/(1.0+exp(arg2))
    ct=  max(ct,0.)

    cl = alpha*cl1*par/sqrt(1.0+(alpha*par)**2)

    iso_factor = ct*cl

  end subroutine get_temp_factor_isoprene
   
  ! ***
  
  subroutine get_temp_factor_monoterp(tk, tref, beta, terp_factor)
  
    ! --- in/out ---------------
    real, intent(in)    ::  tk
    real, intent(in)    ::  tref
    real, intent(in)    ::  beta
    real, intent(out)   ::  terp_factor
    
    ! --- begin ----------------
    
    terp_factor = 1.0 * exp(beta*(tk-tref))
    
  end subroutine get_temp_factor_monoterp
  
  ! ***

  !
  !! depricated, now these latitude dependencies are implemented above
  !! following keywords from the vegetation property table
  !
  !subroutine get_Dmass(ix,iy, code,Dmass)
  !
  !use dims, only : runF
  !
  !integer:: ix,iy
  !real :: code, Dmass, lat
  !integer ::  idum
  !
  !! to avoid errors about unused variables ...
  !idum = ix
  !
  !if(code.eq.-100.) then !Sitka and Norway Spruce
  !  Dmass = 1400.
  !  lat = runF%southb + runf%dlat * (iy-1)
  !  if (lat.ge.60.0) Dmass = 800.
  !  if (lat.lt.55.0) Dmass = 1600.
  !endif
  !
  !if(code.eq.-200.) then ! Scotts Pine
  !  Dmass = 700.
  !  lat = runF%southb + runf%dlat * (iy-1)
  !  if (lat.ge.60.0) Dmass = 500.
  !endif
  !
  !end subroutine get_Dmass



  !-------------------------------------------------------------------
  ! rc_lai: compute one-sided leaf area index
  !         based on Emberson xxx reference
  !-------------------------------------------------------------------
  subroutine slai_emberson(day_of_year,lat,lu,lai,sai) 

    implicit none

    ! Input/output variables
    integer, intent (in)      :: day_of_year ! day of year 
    real   , intent (in)      :: lat         ! latitude (degrees)
    integer, intent (in)      :: lu          ! landuse class
    real   , intent(out)      :: lai         ! one-sided leaf area index
    real   , intent(out)      :: sai         ! surface area index  
    integer, parameter :: nlu =9 ! number of land use type 
    
    ! Define type for LAI (leaf area index) parameters:
    type laitype
       integer     :: sgs50     ! start growing season at 50 degrees latitude (days)
       real        :: dsgs      ! shift in start growing season (days/degree latitude)
                                ! dsgs > 0 -> start-growing-season earlier in the south, later in the north
       integer     :: egs50     ! end growing season at 50 degrees latitude (days)
       real        :: degs      ! shift in end growing season (days/degree latitude)
                                ! degs < 0 -> end-growing-season later in the south, earlier in the north
       real        :: laimin    ! leaf area index at start and end of growing season (m2 leaf/m2 ground surface)
                                ! xxx outside growing season LAI = 0
       real        :: laimax    ! maximal leaf area index (m2 leaf/m2 ground surface)
       integer     :: s_lai_len ! length of starting phase of LAI (days)
       integer     :: e_lai_len ! length of end phase of LAI (days)
    end type laitype

    ! lai_par: LAI parameters derived from Emberson xxx
    !          note that for the LAI of arable land and permanent crops the Emberson class "root crop" has been used;
    !          for other parameters, "temp_crop" has been used.
    !
    !          sgs50    dsgs  egs50    degs  laimin  laimax  s_lai_len  e_lai_len
    type(laitype),dimension(nlu), parameter :: lai_par = (/ &                          ! Emberson   DEPAC
       laitype(    0,    0.0,   366,    0.0,    2.0,    3.5,       140,       135), &  ! grass      1. grass
       laitype(  130,    0.0,   250,    0.0,    0.0,    4.2,        35,        65), &  ! root crop  2. arable land
       laitype(  130,    0.0,   250,    0.0,    0.0,    4.2,        35,        65), &  ! root crop  3. permanent crops
       laitype(    0,    0.0,   366,    0.0,    5.0,    5.0,         1,         1), &  ! temp_conif 4. coniferous forest
       laitype(  100,    1.5,   307,   -2.0,    0.0,    4.0,        20,        30), &  ! temp_decid 5. deciduous forest
       laitype( -999, -999.0,  -999, -999.0, -999.0, -999.0,      -999,      -999), &  ! water      6. water
       laitype( -999, -999.0,  -999, -999.0, -999.0, -999.0,      -999,      -999), &  ! urban      7. urban
       laitype(    0,    0.0,   366,    0.0,    2.0,    3.5,       140,       135), &  ! grass      8. other
       laitype( -999, -999.0,  -999, -999.0, -999.0, -999.0,      -999,      -999) /)  ! desert     9. desert

    ! variables from module
    !       LAI
    !        |
    !  laimax|               ---------
    !        |              /         \
    !        |             /           \
    !        |            /             \
    !        |           /               \
    !        |          /                 \
    !  laimin|         /                   \
    !        |        |                     |
    !        |        |                     |   
    !        ---------|------|--------|-----|-----------
    !                sgs    sgs+     egs-   egs
    !                     s_lai_len  e_lai_len
    !
    ! contents of lai_par:
    ! sgs50     ! start growing season at 50 degrees latitude (days)
    ! dsgs      ! shift in start growing season (days/degree latitude)
    !             dsgs > 0 -> start-growing-season earlier in the south, later in the north
    ! egs50     ! end growing season at 50 degrees latitude (days)
    ! degs      ! shift in end growing season (days/degree latitude)
    !             degs < 0 -> end-growing-season later in the south, earlier in the north
    ! laimin    ! leaf area index at start and end of growing season (m2 leaf/m2 ground surface)
    !           ! xxx outside growing season LAI = 0
    ! laimax    ! maximal leaf area index (m2 leaf/m2 ground surface)
    ! s_lai_len ! length of starting phase of LAI (days)
    ! e_lai_len ! length of end phase of LAI (days)

    ! local variables:
    integer :: sgs   ! start growing season at certain latitude (days)
    integer :: egs   ! end growing season at certain latitude (days)
    ! Compute start and end of growing season for current latitude;
    ! dsgs > 0 -> sgs earlier in the south, later in the north
    ! degs < 0 -> egs later in the south, earlier in the north

    type(laitype) :: lai_par1

    lai_par1 = lai_par(lu)

    sgs = int ( 0.5 +  lai_par1%sgs50 + lai_par1%dsgs * (lat-50.0) );
    egs = int ( 0.5 +  lai_par1%egs50 + lai_par1%degs * (lat-50.0) );

    if (missing(lai_par1%laimax)) then
      lai = -999.0
      sai = -999.0
    else
      ! calculation of lai
      if ( day_of_year < sgs .or. day_of_year > egs ) then
        lai = 0.0;
      else
        if (day_of_year <=  sgs + lai_par1%s_lai_len ) then

          ! (lai    - laimin)   (day_of_year - sgs)
          ! ----------------- = -------------------
          ! (laimax - laimin)        s_lai_len

          lai = lai_par1%laimin + (lai_par1%laimax-lai_par1%laimin)*(day_of_year-sgs)/lai_par1%s_lai_len

        elseif (day_of_year >=  egs - lai_par1%e_lai_len) then

          ! (lai    - laimin)   (egs - day_of_year)
          ! ----------------- = -------------------
          ! (laimax - laimin)        e_lai_len

          lai = lai_par1%laimin + (lai_par1%laimax-lai_par1%laimin)*(egs-day_of_year)/lai_par1%e_lai_len
        else
          lai = lai_par1%laimax
        endif
      endif

      ! calculation of sai
      if (lu .eq. 4 .or. lu .eq. 5) then   ! forest
        sai = lai + 1.                      ! EMEP report 1/2003 Simpson et al.
      elseif (lu .eq. 3) then              ! permanent crop
        sai = lai + .5                      ! pers. comm. Roy Wichink Kruit
      elseif (lu .eq. 2) then              ! arable land
        if (day_of_year < sgs .or. day_of_year > egs) then ! EMEP report 1/2003 Simpson et al.
          sai = lai
        elseif (day_of_year <=  sgs + lai_par1%s_lai_len) then
          sai = max(5.0/3.5*lai,lai + 1.5)            ! max statement is used to avoid hiccup in sai value
        elseif (day_of_year >=  egs - lai_par1%e_lai_len) then
          sai = lai + 1.5
        else
          sai = lai + 1.5
        endif
      else                              ! rest 
        sai = lai
      endif
    endif

  end subroutine slai_emberson
  
  
  ! ***
  

  logical function missing(x)

    real, intent(in) :: x

    ! bandwidth for checking (in)equalities of floats
    real, parameter :: EPS = 1.0e-5

    missing = (abs(x + 999.) .le. EPS)

  end function missing

  
end module
