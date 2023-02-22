!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_EDGAR

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile
  
  use EDGAR_Data, only : T_EDGAR_Data

  use LE_Emis_EDGAR_Composition, only : T_Emis_Composition
! use LE_Emis_Scenario   , only : T_Emis_Scenario
  use LE_Emis_EDGAR_HeightDistr, only : T_Emis_HeightDistr
  use LE_Emis_EDGAR_Time_Prof  , only : T_Emis_Time_Prof
  use LE_Emis_EDGAR_TempVar    , only : T_Emis_TempVar

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Emis_EDGAR

  public  ::  LE_Emis_EDGAR_Init, LE_Emis_EDGAR_Done
  public  ::  LE_Emis_EDGAR_Setup


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_EDGAR'
  
  ! maximum number of categories etc:
  integer, parameter   ::  maxcat = 50
  integer, parameter   ::  maxzone = 100

  ! for temporary files that are closed directly after use
  integer, parameter :: u_tmp  = 61

  ! --- local --------------------------------

  ! emission data base
  type T_Emis_EDGAR
    ! label assigned to this emission:
    character(len=32)                ::  label

    ! year for which this data set was made:
    ! (in future, implement this with a time range [t1,t2])
    integer                           ::  year
    integer                           ::  month
    ! last year in current inventory:
    integer                           ::  maxyear
    
    ! settings:
    character(len=512)                ::  rcfile
    character(len=64)                 ::  rckey
    
    ! tracer variables:
    integer                           ::  ntr
    character(len=32), pointer        ::  name_tr(:)
    
    ! emission categories
    integer, pointer                  ::  ncat(:)        ! (ntr)
    character(len=64), pointer        ::  catnames(:,:)  ! (ntr,maxcat)
    integer, pointer                  ::  icat_global(:,:)  ! (ntr,maxcat)
    
    ! emission categories
    integer                           ::  nemiscat        
    character(len=64)                 ::  emiscats(maxcat)  

    ! Geographical zone info
    integer                           ::  nzone
    character(len=4)                  ::  zones(maxzone) 
    integer, pointer                  ::  izone(:,:)  ! (nx,ny) zone index on le grid

    ! skip flags per model tracer:
    logical, allocatable              ::  skip_tracer(:)   ! (nspec)
    integer                           ::  ispec
  
    ! molecule masses:
    real , pointer                   ::  xm(:)  ! kg/mole

    ! base emissions (will be re-read at start of new year)
    ! ~ template for input file:
    character(len=1024)             ::  edg_fname_template
    ! ~ actual data:
    type(T_EDGAR_Data), pointer     ::  edg(:,:)  ! (ntr,maxcat)
    
    ! emisison compositions:
    type(T_Emis_Composition), pointer  ::  emcomp(:)

    ! height distribution:
    logical                           ::  with_height_distribution
    type(T_Emis_HeightDistr)          ::  hdistr

    ! time factors country/snap dependent
    type(T_Emis_Time_Prof)            ::  time_prof
    ! template for input files:
    character(len=1024)               ::  time_prof_template
    character(len=512)                ::  time_prof_query

    ! temperature variations for VOC and CO
    type(T_Emis_TempVar)            ::  temp_var_VOC
    type(T_Emis_TempVar)            ::  temp_var_CO
    
    ! At which time a new emission file should be read in? 
    ! beginning each year, month?...
    character(len=16)               ::  emisset_valid

  end type T_Emis_EDGAR



contains


  subroutine LE_Emis_EDGAR_Init( emd, label, rcF, rckey, t, status )

    use GO        , only : TrcFile
    use GO        , only : TDate
    use LE_Data   , only : LE_Data_Enable

    ! --- in/out ------------------------------

    type(T_Emis_EDGAR), intent(out)     ::  emd
    character(len=*), intent(in)        ::  label
    type(TrcFile), intent(in)           ::  rcF
    character(len=*), intent(in)        ::  rckey
    type(TDate), intent(in)             ::  t         ! start time
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_EDGAR_Init'

    ! --- local ------------------------------
    
    ! --- begin -------------------------------

    ! enable data:
    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'area', status )
    IF_NOTOK_RETURN(status=1)

    ! this is the year for which the currently read data is valid ;
    ! nothing read yet, this will be done during the first call
    ! to the Setup routine:
    emd%year  = -1
    emd%month = -1
    
    ! store label:
    emd%label = trim(label)
    
    ! store info on settings:
    emd%rcfile = trim(rcF%fname)
    emd%rckey = trim(rckey)

    ! read emissions    
    call LE_Emis_EDGAR_Read( emd, t%year, t%month, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok  
    status = 0
  
  end subroutine LE_Emis_EDGAR_Init
  
  
  ! ***
   
    
  subroutine LE_Emis_EDGAR_Done( emd, status )

    ! --- in/out ------------------------------

    type(T_Emis_EDGAR), intent(inout)    ::  emd
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_EDGAR_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------

    ! any data loaded ?
    if ( emd%year > 0 .or. emd%month > 0 ) then
    
      ! cleanup:
      call LE_Emis_EDGAR_Clear( emd, status )
      IF_NOTOK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine LE_Emis_EDGAR_Done
  
  ! ===============================================================


  subroutine LE_Emis_EDGAR_Read( emd, year, month, status )

    use Binas     , only : xm_H, xm_C, xm_N, xm_O, xm_S
    use GO        , only : ReadRc, Init, Done
    use GO        , only : goMatchValue, goSplitString, goReplace, goReadFromLine 
    use GO        , only : goGetFU
    use Indices   , only : nspec_all, specname, nspec
    use EDGAR_Data, only : EDGAR_Data_Init
    use EDGAR_Data, only : EDGAR_Data_Get
    use LE_Emis_EDGAR_Composition, only : LE_Emis_Composition_Init
    use LE_Emis_EDGAR_Heightdistr, only : LE_Emis_Heightdistr_Init
    use LE_Emis_EDGAR_Time_Prof, only : LE_Emis_Time_Prof_Init
    use LE_Emis_EDGAR_TempVar, only : LE_Emis_TempVar_Init
    use LE_Grid, only : ugg
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    type(T_Emis_EDGAR), intent(inout)   ::  emd
    integer, intent(in)                 ::  year
    integer, intent(in)                 ::  month
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_EDGAR_Init'

    ! --- local ------------------------------
    
    type(TrcFile)                  ::  rcF
    character(len=1024)            ::  fname_template
    character(len=256)             ::  tracer_names
    character(len=256)             ::  line
    character(len=1024)            ::  catfile
    character(len=1024)            ::  zonecodefile
    character(len=1024)            ::  emiscodefile
    logical                        ::  exist
    integer                        ::  fu
    integer                        ::  itr
    integer                        ::  icat
    integer                        ::  ispec
    integer                        ::  izone
    integer                        ::  i, j
    character(len=1)               ::  zone_band
    character(len=4)               ::  zone_name
    real                           ::  lon360
    integer                        ::  zone_ilon
    integer                        ::  iemiscat   !to read all emission categories possible
    ! emission composition:
    character(len=512)              ::  query
    ! skip list:
    character(len=32)               ::  emspec
    ! temperature profiles:
    real, allocatable               ::  temp_var(:,:)
    integer                         ::  t0, t1, dt, temp, ncol
    ! PM emission indices
    integer                         ::  iemis_pm25_fossil, iemis_pm25_bio, iemis_pm25, iemis_pm10
    integer                         ::  icat_pm25_fossil, icat_pm25_bio, icat_pm25
    character(len=64)               ::  catname
    
    
    ! --- begin -------------------------------

    ! init file:
    write( gol, '(a)' ) trim(emd%rcfile) ; call GoPr
    call Init( rcF, trim(emd%rcfile), status )
    IF_NOTOK_RETURN(status=1)

    ! max year for which base emissions are available:
    call ReadRc( rcF, trim(emd%rckey)//'.maxyear', emd%maxyear, status )
    IF_NOTOK_RETURN(status=1)

    ! Valid period for the emission set
    call ReadRc( rcF, trim(emd%rckey)//'.timerange_valid', emd%emisset_valid, status )
    IF_NOTOK_RETURN(status=1)

    !--------------------------------------------------------
    write( gol, '(" set names of emitted components ...")' ) ; call GoPr
    !--------------------------------------------------------
    
    ! setup storage for tracer fields:
    allocate( emd%name_tr (nspec_all) )
    ! read name of tracers:
    call ReadRc( rcF, trim(emd%rckey)//'.tracers', tracer_names, status )
    IF_NOTOK_RETURN(status=1)
    ! read values from input line:
    call goSplitString( tracer_names, emd%ntr, emd%name_tr, status )
    IF_NOTOK_RETURN(status=1)

    !--------------------------------------------------------
    write( gol, '(" Read all possible emission categories... " )' ) ; call GoPr
    !--------------------------------------------------------
    ! 
    ! read name of directory with file:
    call ReadRc( rcF, trim(emd%rckey)//'.sourcecodes', emiscodefile, status )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    inquire( file=trim(emiscodefile), exist=exist )
    if ( .not. exist ) then
      write (gol,'("source codes file not found: ",a)') trim(emiscodefile); call goErr
      TRACEBACK; status=1; return
    end if
    ! obtain free file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open:
    open( unit=fu, file=trim(emiscodefile), form='formatted', status='old', iostat=status )
    IF_NOTOK_RETURN(status=1)
    ! init counter:
    iemiscat = 0
    ! loop over lines:
    do
      ! read next line:
      read (fu,'(a)',iostat=status) line
      if (status<0) exit  ! end of file
      IF_NOTOK_RETURN(status=1)
      ! skip empty lines and comment:
      if ( len_trim(line) == 0 ) cycle
      if ( line(1:1) == '#' ) cycle
      ! found; increase counter:
      iemiscat = iemiscat + 1
      ! check ...
      if ( iemiscat > maxcat ) then
        write (gol,'("found more emission categories in file `",a,"` than current maximum of ",i6)') trim(emiscodefile), maxcat; call goErr
        write (gol,'("please increase value of parameter `maxcat` in this module")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! content of line is the category name:
      emd%emiscats(iemiscat) = trim(line)
    end do
    ! store counter:
    emd%nemiscat = iemiscat
    ! close:
    close( unit=fu, iostat=status )
    IF_NOTOK_RETURN(status=1)

    !--------------------------------------------------------
    write( gol, '(" Read category names ...")') ; call GoPr
    !--------------------------------------------------------
    
    ! number of tracers is known now, create storage for category codes:
    allocate( emd%ncat(emd%ntr) )
    allocate( emd%catnames(emd%ntr,maxcat) )
    allocate( emd%icat_global(emd%ntr,maxcat) )
    emd%icat_global = -999
    ! read name of directory with category files:
    call ReadRc( rcF, trim(emd%rckey)//'.sources', fname_template, status )
    IF_NOTOK_RETURN(status=1)
    ! loop over tracers:
    do itr = 1, emd%ntr
      ! name of category file, replace some values:
      catfile = trim(fname_template)
      call goReplace( catfile , '%{tracer}', trim(emd%name_tr(itr)), status )
      IF_NOTOK_RETURN(status=1)   
      ! check ...
      inquire( file=trim(catfile), exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found : ",a)') trim(catfile); call goErr
        TRACEBACK; status=1; return
      end if
      ! obtain free file unit:
      call goGetFU( fu, status )
      IF_NOTOK_RETURN(status=1)
      ! open:
      open( unit=fu, file=trim(catfile), form='formatted', status='old', iostat=status )
      IF_NOTOK_RETURN(status=1)
      ! init counter:
      icat = 0
      ! loop over lines:
      do
        ! read next line:
        read (fu,'(a)',iostat=status) line
        if (status<0) exit  ! end of file
        IF_NOTOK_RETURN(status=1)
        ! skip empty lines and comment:
        if ( len_trim(line) == 0 ) cycle
        if ( line(1:1) == '#' ) cycle
        ! found; increase counter:
        icat = icat + 1
        ! check ...
        if ( icat > maxcat ) then
          write (gol,'("found more categories in file `",a,"` than current maximum of ",i6)') trim(catfile), maxcat; call goErr
          write (gol,'("please increase value of parameter `maxcat` in this module")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! check if already in use ...
        if ( icat > 1 ) then
          ! loop over previous:
          do i = 1, icat-1
            ! compare:
            if ( trim(line) == trim(emd%catnames(itr,i)) ) then
              write (gol,'("category `",a,"` defined twice .. ?")') trim(line); call goErr
              TRACEBACK; status=1; return
            end if
          end do ! i
        end if  ! icat > 1
        ! content of line is the category name:
        emd%catnames(itr,icat) = trim(line)
        ! search index in global list
        call goMatchValue( trim(line), emd%emiscats(1:emd%nemiscat), emd%icat_global(itr,icat), status )
        IF_NOTOK_RETURN(status=1)
      end do
      ! store counter:
      emd%ncat(itr) = icat
      ! close:
      close( unit=fu, iostat=status )
      IF_NOTOK_RETURN(status=1)
    end do

    !-------------------------------------------
    write (gol,' (" Set molecule masses... ")') ; call GoPr
    !-------------------------------------------
    
    allocate( emd%xm(emd%ntr) )
    ! loop over emisisons:
    do itr = 1, emd%ntr
      ! switch:
      select case ( trim(emd%name_tr(itr)) )
        !--> assume full equivalents similar as in TNO emissions,
        !    otherwise explit 'kg(N)' or 'kg(S)' should be used:
        case ( 'SO2' ) ; emd%xm(itr) = xm_S + xm_O * 2  ! kg(SO2)/mole
        case ( 'NOx' ) ; emd%xm(itr) = xm_N + xm_O * 2  ! kg(NO2)/mole
        !<---
        case ( 'NH3' ) ; emd%xm(itr) = xm_N + xm_H * 3  ! kg/mole
        case ( 'CO'  ) ; emd%xm(itr) = xm_C + xm_O      ! kg/mole
        case ( 'CH4' ) ; emd%xm(itr) = xm_C + xm_H * 4  ! kg/mole
        case default   ; emd%xm(itr) = -999.9   ! to be trapped ...
      end select
    end do

    !--------------------------------------------------------
    write(gol, '(" setup base emissions ...")' ) ; call GoPr
    !--------------------------------------------------------
        
    ! storage for emission files for each of the tracer and category pair in use;
    allocate( emd%edg(emd%ntr,maxcat) )
    
    ! read template for name of data file:
    call ReadRc( rcF, trim(emd%rckey)//'.file', emd%edg_fname_template, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over tracers:
    do itr = 1, emd%ntr
      ! loop over categories for this tracer:
      do icat = 1, emd%ncat(itr)
        ! init file:
        call EDGAR_Data_Init( emd%edg(itr,icat), trim(emd%edg_fname_template), &
                                trim(emd%name_tr(itr)), &
                                trim(emd%catnames(itr,icat)), &
                                status )
        IF_NOTOK_RETURN(status=1)

        ! read data
        call EDGAR_Data_Get( emd%edg(itr,icat), min(emd%maxyear,year), month, status )
        IF_NOTOK_RETURN(status=1)
        
      end do  ! cats
    end do  ! tracers
    
    ! Search index for PM-coarse emissions
    iemis_pm10 = -1
    do itr = 1, emd%ntr
      select case (trim(emd%name_tr(itr) ) )
        case ( 'PM10' )
          iemis_pm10 = itr
        case ( 'BC', 'CH4', 'CO', 'NH3', 'NMVOC', 'NOx', 'OC', &
               'PM2.5_bio', 'PM2.5_fossil', 'PM2.5', 'PMcoarse', 'SO2', &
               'CO2_excl_short-cycle_org_C', 'CO2_org_short-cycle_C', 'N2O', &
               'CBIV_ALD2', 'CBIV_ETH', 'CBIV_FORM', 'CBIV_ISOP', 'CBIV_OLE', &
               'CBIV_PAR', 'CBIV_TOL', 'CBIV_XYL' )
          ! not PM10
        case default     
          write (gol,'("do not know whether tracer `",a,"` represents coarse mode aerosol or something else")') trim(emd%name_tr(itr)); call goErr
          TRACEBACK; status=1; return
      end select
    end do
    
    ! found PM10 emission?
    if ( iemis_pm10 > 0 ) then
      ! search for index of fine mode emissions:
      iemis_pm25_fossil = -1
      iemis_pm25_bio    = -1
      iemis_pm25        = -1
      do itr = 1, emd%ntr
        select case ( trim(emd%name_tr(itr) ) )
          case ( 'PM2.5_fossil' )
            iemis_pm25_fossil = itr
          case ( 'PM2.5_bio' )
            iemis_pm25_bio = itr
          case ( 'PM2.5' )
            iemis_pm25 = itr
          case ( 'BC', 'CH4', 'CO', 'NH3', 'NMVOC', 'NOx', 'OC', 'SO2', 'PM10', &
                 'CO2_excl_short-cycle_org_C', 'CO2_org_short-cycle_C', 'N2O' )
            ! not PM2.5
          case default 
            write (gol,'("do not know whether tracer `",a,"` represents fine mode aerosol or something else")') trim(emd%name_tr(itr)); call goErr
          TRACEBACK; status=1; return
        end select
      end do
      
      ! check ..
      if ( ( iemis_pm25_fossil < 0 .or. iemis_pm25_bio < 0 ) .and. iemis_pm25 < 0) then
        write( gol, '("Could not find fine mode pm source in emitted tracers :")' ); call GoErr
        do itr = 1, emd%ntr
          write(gol,'("  ",a)') trim(emd%name_tr(itr)); call goErr
        end do
        write( gol, '("Make clear what to do with PM10 emissions" )' ); call GoErr
        TRACEBACK;status=1;return
      end if
      
      ! found single PM2.5 emitted tracer"
      if ( iemis_pm25 > 0 ) then
        ! info ...
        write( gol, '(" Convert `",a,"` to `pm25_pm10` using `",a,"`")' ) &
                trim(emd%name_tr(iemis_pm10)) , trim(emd%name_tr(iemis_pm25)); call GoPr
        !substract
        do icat = 1, emd%ncat(iemis_pm10)
          catname = emd%catnames(iemis_pm10,icat)

          ! Select category index in PM2.5 emission
          call GoMatchValue( trim(catname), emd%catnames(iemis_pm25,:), icat_pm25, status, quiet = .true. )
          if ( status < 0 ) then
            ! category not found in fine emissions --> all emissions coarse mode??
            write( gol, '("Category name: `",a, "` not found in fine pm emissions, all emission considered as coarse mode ")' ) trim(catname) ; call GoPr
          else if ( status > 0 ) then
            write( gol, '("Error matching category name in PM2.5 emissions: ",a )' ) trim( catname) ; call GoErr
            TRACEBACK;status=1;return
          else 
            ! subtract
            emd%edg(iemis_pm10,icat)%emis = emd%edg(iemis_pm10,icat)%emis - emd%edg(iemis_pm25,icat_pm25)%emis
            ! Error message for PM coarse fraction is negative
            if ( any( emd%edg(iemis_pm10,icat)%emis < 0.0 ) ) then
              write( gol, '("Found larger emission for PM2.5 than for PM10 --> negative PMcoarse emission")' ) ; call GoErr
              write( gol, '("Category: ", a)' ) trim(catname) ; call GoErr
              write (gol, '("IF YOU KNOW WHAT YOU ARE DOING, COMMENT THE FOLLOWING TRACEBACK (#1/3) TO IGNORE THIS ERROR:")'); call goErr
              TRACEBACK;status=1;return
            end if
            
            ! ====================================
            ! ==== Switched of Error handling ====
            ! ==== set negative values to 0   ====
            emd%edg(iemis_pm10,icat)%emis = max( 0.0, emd%edg(iemis_pm10,icat)%emis )
            ! ====================================
          end if
        end do  ! categories
      
      ! found specified fossil and bio pm25 emissions
      else if ( iemis_pm25_fossil > 0 .and. iemis_pm25_bio > 0 )  then   

        ! info ...
        write( gol, '(" Convert `",a,"` to `pm25_pm10` using `",a,"` and `", a,"`")' ) &
                trim(emd%name_tr(iemis_pm10)) , trim(emd%name_tr(iemis_pm25_fossil)), trim(emd%name_tr(iemis_pm25_bio)) ; call GoPr
        !substract
        do icat = 1, emd%ncat(iemis_pm10)
          ! current category name for total pm10:
          catname = emd%catnames(iemis_pm10,icat)

          ! Select category index in PM2.5 emission
          call GoMatchValue( trim(catname), emd%catnames(iemis_pm25_fossil,:), icat_pm25_fossil, status, quiet = .true. )
          if ( status < 0 ) then
            ! category not found in fine emissions --> all emissions coarse mode??
            write( gol, '("Category name: `",a, "` not found in fossil fine pm emissions, emissions considered as coarse or bio ")' ) trim(catname) ; call GoPr
          else if ( status > 0 ) then
            write( gol, '("Error matching category name in PM2.5 emissions: ",a )' ) trim( catname) ; call GoErr
            TRACEBACK;status=1;return
          else 
            ! subtract
            emd%edg(iemis_pm10,icat)%emis = emd%edg(iemis_pm10,icat)%emis - emd%edg(iemis_pm25_fossil,icat_pm25_fossil)%emis
            ! Error message for PM coarse fraction is negative
            if ( any( emd%edg(iemis_pm10,icat)%emis < 0.0 ) ) then
              write( gol, '("Found larger emission for PM2.5 than for PM10 --> negative PMcoarse emission")' ) ; call GoErr
              write( gol, '("Category: ", a)' ) trim(catname) ; call GoErr
              write (gol, '("IF YOU KNOW WHAT YOU ARE DOING, COMMENT THE FOLLOWING TRACEBACK (#2/3) TO IGNORE THIS ERROR:")'); call goErr
              TRACEBACK;status=1;return
            end if
            
            ! ====================================
            ! ==== Switched of Error handling ====
            ! ==== set negative values to 0   ====
            emd%edg(iemis_pm10,icat)%emis = max( 0.0, emd%edg(iemis_pm10,icat)%emis )
            ! ====================================
          end if
          
          ! Select category index in PM2.5 bio emission
          call GoMatchValue( trim(catname), emd%catnames(iemis_pm25_bio,:), icat_pm25_bio, status, quiet = .true. )
          if ( status < 0 ) then
            ! category not found in fine emissions --> all emissions coarse mode??
            write( gol, '("Category name: `",a, "` not found in bio fine pm emissions, emissions considered as coarse or fossil ")' ) trim(catname) ; call GoPr
          else if ( status > 0 ) then
            write( gol, '("Error matching category name in PM2.5 emissions: ",a )' ) trim( catname) ; call GoErr
            TRACEBACK;status=1;return
          else 
            ! subtract
            emd%edg(iemis_pm10,icat)%emis = emd%edg(iemis_pm10,icat)%emis - emd%edg(iemis_pm25_bio,icat_pm25_bio)%emis
            ! Error message for PM coarse fraction is negative
            if ( any( emd%edg(iemis_pm10,icat)%emis < 0.0 ) ) then
              write( gol, '("Found larger emission for PM2.5 than for PM10 --> negative PMcoarse emission")' ) ; call GoErr
              write( gol, '("Category: ", a)' ) trim(catname) ; call GoErr
              write (gol, '("IF YOU KNOW WHAT YOU ARE DOING, COMMENT THE FOLLOWING TRACEBACK (#3/3) TO IGNORE THIS ERROR:")'); call goErr
              TRACEBACK;status=1;return
            end if
            
            ! ====================================
            ! ==== Switched of Error handling ====
            ! ==== set negative values to 0   ====
            emd%edg(iemis_pm10,icat)%emis = max( 0.0, emd%edg(iemis_pm10,icat)%emis )
            ! ====================================
          end if
        end do  ! categories
      end if  ! which pm2.5 emissions found?
      
      ! change name if necessary
      if ( iemis_pm25 > 0 .or. iemis_pm25_bio > 0 .or. iemis_pm25_fossil > 0 ) then
      
        ! PM10 emissions --> only coarse mode emissions 
        emd%name_tr(iemis_pm10) = 'PM25_PM10' 
      end if
    end if  ! pm 10 emission found?
        
    !--------------------------------------------------------
    write( gol, '(" Read geographical zones ...")' ) ; call GoPr
    !--------------------------------------------------------

    ! read name of directory with files:
    call ReadRc( rcF, trim(emd%rckey)//'.zones', zonecodefile, status )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    inquire( file=trim(zonecodefile), exist=exist )
    if ( .not. exist ) then
      write (gol,'("zone codes file not found: ",a)') trim(zonecodefile); call goErr
      TRACEBACK; status=1; return
    end if
    ! obtain free file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open:
    open( unit=fu, file=trim(zonecodefile), form='formatted', status='old', iostat=status )
    IF_NOTOK_RETURN(status=1)
    ! init counter:
    izone = 0
    ! loop over lines:
    do
      ! read next line:
      read (fu,'(a)',iostat=status) line
      if (status<0) exit  ! end of file
      IF_NOTOK_RETURN(status=1)
      ! skip empty lines and comment:
      if ( len_trim(line) == 0 ) cycle
      if ( line(1:1) == '#' ) cycle
      ! found; increase counter:
      izone = izone + 1
      ! check ...
      if ( izone > maxzone ) then
        write (gol,'("found more zones in file `",a,"` than current maximum of ",i6)') trim(zonecodefile), maxzone; call goErr
        write (gol,'("please increase value of parameter `maxzones` in this module")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! content of line is the category name:
      emd%zones(izone) = trim(line)
    end do
    ! store counter:
    emd%nzone = izone
    ! close:
    close( unit=fu, iostat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! fill zone index on le grid:
    ! storage:
    allocate( emd%izone(ugg%nlon,ugg%nlat) )
    ! create zone labels:  [STN](000|015|..|345)
    ! loop over le latitudes:
    do j = 1, ugg%nlat
      ! loop over le longitudes
      do i = 1, ugg%nlon
        ! set south/tropics/north:
        if ( ugg%latitude(i,j) < -30 ) then
          zone_band = 'S'
        else if ( ugg%latitude(i,j) > 30 ) then
          zone_band = 'N'
        else
          zone_band = 'T'
        end if
        ! lon in [0,360] :
        lon360 = modulo( ugg%longitude(i,j), 360.0 )
        ! round to mid of [-7.5,7.5] zones around 0.0, 15.0, ...
        zone_ilon = nint(lon360/15.0)*15
        if ( zone_ilon == 360 ) zone_ilon = 0
        ! fill name:
        write (zone_name,'(a,i3.3)') zone_band, zone_ilon
        ! search in list:
        call goMatchValue( zone_name, emd%zones, emd%izone(i,j), status )
        IF_NOTOK_RETURN(status=1)
      end do
    end do
    
    !--------------------------------------------------------
    write( gol,'(" Setup input of time profiles ...")' ) ; call GoPr
    !--------------------------------------------------------

    call ReadRc( rcF, trim(emd%rckey)//'.time_profiles', emd%time_prof_query, status )
    IF_NOTOK_RETURN(status=1)

    ! store current year/month:    
    emd%year  = year
    emd%month = month
    
    !--------------------------------------------------------
    write (gol,'("    read emission time profiles for year ",i4," ...")') year; call goPr
    !--------------------------------------------------------

    ! copy query:
    line = trim(emd%time_prof_query)
    ! replace some values:
    call goReplace( line, '%{year}', '(i4.4)', emd%year, status )
    IF_NOTOK_RETURN(status=1)

    ! init with just the skip flag, this will allow proper done:
    call LE_Emis_Time_Prof_Init( emd%time_prof, trim(line), &
                                   emd%emiscats(1:emd%nemiscat), &
                                  emd%zones(1:emd%nzone), status )
    IF_NOTOK_RETURN(status=1)


    !--------------------------------------------------------
    write (gol, '(" Read skip list...")') ; call GoPr
    !--------------------------------------------------------

    ! storage per model tracer:
    allocate( emd%skip_tracer(nspec) )
    ! by default skip nothing:
    emd%skip_tracer = .false.

    ! list with model tracer names to be skipped:
    call ReadRc( rcF, trim(emd%rckey)//'.skip.species', query, status )
    IF_NOTOK_RETURN(status=1)
    ! loop over list items:
    do
      ! leave if empty:
      if ( len_trim(query) == 0 ) exit
      ! extract:
      call goReadFromLine( query, emspec, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      ! index:
      call goMatchValue( trim(emspec), specname(1:nspec), ispec, status, quiet=.true. )
      IF_ERROR_RETURN(status=1)
      ! found ?
      if ( ispec > 0 ) then
        ! reset flag:
        emd%skip_tracer(ispec) = .true.
        ! info ...
        write (gol,'("    skip tracer ",a," ...")') trim(specname(ispec)); call goPr
      end if
    end do
    
    !--------------------------------------------------------
    write(gol, '(" Read composition tables")' ) ; call GoPr
    !--------------------------------------------------------
    
    ! storage for composition data
    allocate( emd%emcomp(emd%ntr) )
    ! loop over tracers:
    do itr = 1, emd%ntr
      ! component fractions query:
      call ReadRc( rcF, trim(emd%rckey)//'.composition.'//trim(emd%name_tr(itr)), query, status )
      IF_NOTOK_RETURN(status=1)
      ! special or not ?
      if ( trim(query) == 'special' ) query = 'skip=T'
      ! fill:
      call LE_Emis_Composition_Init( emd%emcomp(itr), query, &
                                    trim(emd%name_tr(itr)), &
                                    emd%catnames(itr,1:emd%ncat(itr)), &
                                    emd%zones(1:emd%nzone), &
                                    specname(1:nspec), status )
      IF_NOTOK_RETURN(status=1)
    end do  ! tracers

    !--------------------------------------------------------
    write(gol, '(" Read height distribution")' ) ; call GoPr
    !--------------------------------------------------------

    ! key to set or read distribution:
    call ReadRc( rcF, trim(emd%rckey)//'.height_distribution', query, status )
    IF_NOTOK_RETURN(status=1)
    ! if not empty than this is enabled ...
    emd%with_height_distribution = len_trim(query) > 0

    ! defined ?
    if ( emd%with_height_distribution ) then
      ! setup:
      call LE_Emis_HeightDistr_Init( emd%hdistr, query, emd%emiscats(1:emd%nemiscat), status )
      IF_NOTOK_RETURN(status=1)
    end if

    !--------------------------------------------------------
    write(gol, '(" Read temperature dependent emission factors")' ) ; call GoPr
    !--------------------------------------------------------

    ! key to set or read distribution:
    call ReadRc( rcF, trim(emd%rckey)//'.temp_var_voc', query, status )
    IF_NOTOK_RETURN(status=1)

    call LE_Emis_TempVar_Init( emd%temp_var_VOC, query, &
                              emd%emiscats(1:emd%nemiscat), status )
    IF_NOTOK_RETURN(status=1)

    call ReadRc( rcF, trim(emd%rckey)//'.temp_var_co', query, status )
    IF_NOTOK_RETURN(status=1)

    call LE_Emis_TempVar_Init( emd%temp_var_CO, query, &
                              emd%emiscats(1:emd%nemiscat), status )
    IF_NOTOK_RETURN(status=1)

    ! close:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)

    
    ! ok
    status = 0

  end subroutine LE_Emis_EDGAR_Read


  ! ***


  subroutine LE_Emis_EDGAR_Clear( emd, status )

    use EDGAR_Data, only : EDGAR_Data_Done
    use LE_Emis_EDGAR_Composition, only : LE_Emis_Composition_Done
    use LE_Emis_EDGAR_Heightdistr, only : LE_Emis_Heightdistr_Done
    use LE_Emis_EDGAR_Time_Prof, only : LE_Emis_Time_Prof_Done
    use LE_Emis_EDGAR_TempVar, only : LE_Emis_TempVar_Done

    ! --- in/out ------------------------------

    type(T_Emis_EDGAR), intent(inout)   ::  emd
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_EDGAR_Done'

    ! --- local -------------------------------
    
    integer   ::  itr
    integer   ::  icat

    ! --- begin -------------------------------
    
    ! Done with emission files for each of the tracer and category pair in use;
    ! loop over tracers:
    do itr = 1, emd%ntr
      ! loop over categories for this tracer:
      do icat = 1, emd%ncat(itr)
        ! close file:
        call EDGAR_Data_Done( emd%edg(itr,icat), status )
        IF_NOTOK_RETURN(status=1)
      end do  ! cats
    end do  ! tracers

    ! loop over tracers:
    do itr = 1, emd%ntr
      ! done with composition:
      call LE_Emis_Composition_Done( emd%emcomp(itr), status )
      IF_NOTOK_RETURN(status=1)
    end do  ! tracers
    deallocate( emd%emcomp )
    
    ! Done with height distribution?
    if ( emd%with_height_distribution ) then
      call LE_Emis_Heightdistr_Done( emd%hdistr, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! Done with time distribution
    call LE_Emis_Time_Prof_Done( emd%time_prof, status )
    IF_NOTOK_RETURN(status=1)

    ! Done with temperature variable emission factors
    call LE_Emis_TempVar_Done( emd%temp_var_VOC, status )
    IF_NOTOK_RETURN(status=1)
    call LE_Emis_TempVar_Done( emd%temp_var_CO, status )
    IF_NOTOK_RETURN(status=1)
    

    ! clear:
    deallocate( emd%name_tr )
    deallocate( emd%ncat )
    deallocate( emd%catnames )
    deallocate( emd%icat_global )
    deallocate( emd%edg )
    deallocate( emd%skip_tracer )
    
    ! ok
    status = 0

  end subroutine LE_Emis_EDGAR_Clear
  

  ! ***


  subroutine LE_Emis_EDGAR_Setup( emd, emis_a, t1, t2, status )

    use Binas     , only : xm_N, xm_S, Avog
    use GO        , only : TDate, Get, operator(+), operator(-), operator(/)
    use GO       , only : GO_Timer_Start, GO_Timer_End
    use GO       , only : calc_DayNumber, days_in_month
    use Num      , only : IntervalSum, Interp_Lin, EXTRAPOL_CONSTANT
    use Grid      , only : AreaOper
    use Dims      , only : nx, ny, nz
    use Indices   , only : nspec, specname, specunit, specmolm
#ifdef with_m7
    ! if it is the M7 mode, we also need to added the number concentraion
    !   nnus    , number concentration nucle. mode soluble 
    !   nais    , number concentration aitken mode soluble 
    !   nacs    , number concentration accum. mode soluble 
    !   ncos    , number concentration coarse mode soluble 
    !   naii    , number concentration aitken mode insoluble 
    !   naci    , number concentration accum. mode insoluble 
    !   ncoi    , number concentration coarse mode insoluble 
    use Indices    , only : ispec_so4ks, ispec_so4as, ispec_ocki, ispec_ocks, ispec_bcki
    use Indices    , only : i_nais, i_nacs, i_naii 
    use indices    , only : tracer_is_m7
    use mo_aero_m7 , only : sigmaln
    use binas      , only: pi
#endif

    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out ---------------------------

    type(T_Emis_EDGAR), intent(inout)     ::  emd
    real, intent(inout)                   ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_EDGAR_Setup'

    ! number of hours in a year of 365 days:
    real, parameter  ::  hour_per_year = 8760.0
    
    ! conversion factors:
    real, parameter  ::  sec_per_min   = 60.0    ! (sec/min)
    real, parameter  ::  min_per_hour  = 60.0    ! (min/hour)
    real, parameter  ::  min_per_year  = min_per_hour * hour_per_year  ! (min/year)
    real, parameter  ::  ug_per_kg     = 1.0e9   ! (ug/kg)
    real, parameter  ::  ug_per_MEg    = 1.0e12  ! (ug/Mg)
    real, parameter  ::  kg_per_MEg    = 1.0e3   ! (kg/Mg)
    real, parameter  ::  mol_per_MEmol = 1.0e6   ! (mol/Mmol)
    
    ! conversion from Mg/year to ug/min :            ug/Mg  / (hour/year)   / (min/hour)
    real, parameter  ::  Mg_per_year_to_ug_per_min = 1.0e12 / hour_per_year /  min_per_hour  ! (ug/mn)/(Mg/year)
    
    ! conversion from Mmol/year to mol/min:          mol/Mmol / (hour/year)   / (min/hour)
    real, parameter   :: Mmol_per_year_to_mol_per_min = 1e6   / hour_per_year / min_per_hour  ! (mol/min)/(Mmol/year)
    
    ! --- local ----------------------------
    
    integer   ::  itr
    integer   ::  icat, icat_glob

    real, allocatable    ::  em(:,:)
    character(len=32)    ::  em_units
    type(TDate)          ::  tmid
    integer              ::  yy, mm, dd, hh
    integer              ::  ix, iy, iz
    integer              ::  icomp
    integer              ::  izone
    integer              ::  ispec
    integer              ::  ilast

    real                 ::  compfrac
    character(len=128)   ::  conversion
    real                 ::  convfac

    real                 ::  temperfrac

    real                 ::  timefrac
    integer              ::  ihour, julian_day
    real, allocatable    ::  hhb(:)
    
    integer              ::  min_per_month

    ! pre-computed height distribution profiles:
    real, allocatable    ::  hd_profiles(:,:,:,:)  ! (nx,ny,nz,ncat)
    real                 ::  emtop
    
    ! help array:
    real, allocatable     ::  delta_emis(:)
    real, allocatable     ::  year_emis(:)

    ! meteo data:
    real, pointer               ::  tsurf(:,:,:)   ! (lon,lat,1)
    real, pointer               ::  h_m(:,:,:)   ! (lon,lat,1)
    real, pointer               ::  area(:,:,:)   ! (lon,lat,1)

#ifdef with_m7
    ! jianbing: these parameters are results discussed with Nick and Bas 
    !densities and radii of emitted species
    real, parameter:: rad_so4ait=0.03e-6 !(m)
    real, parameter:: rad_so4acc=0.075e-6 !(m) 
    ! real, parameter:: rad_ssa=0.079e-6 !(m)
    ! real, parameter:: rad_ssc=0.63e-6 !(m)
    real, parameter:: rad_bcai=0.034e-6 !(m)
    real, parameter:: rad_ocai=0.030e-6  !(m)

    ! real, parameter:: rad_ssac=0.079e-6 !(m)
    ! real, parameter:: rad_ssco=0.63e-6 !(m)
    ! real, parameter:: rad_duai=0.079e-6 !(m)
    ! real, parameter:: rad_duci=0.63e-6 !(m)
    ! real, parameter:: rad_duas=0.079e-6 !(m)
    ! real, parameter:: rad_ducs=0.63e-6 !(m)

    real, parameter:: dens_so4=1841. !kg/m3
    real, parameter:: dens_ss=2165. !kg/m3
    real, parameter:: dens_bc=2000.  !kg/m3
    real, parameter:: dens_oc=2000 !kg/m3
    ! real, parameter:: dens_du=2650 !kg/m3
    real, parameter:: fracsol_oc=0.65

    real, parameter:: convcm=1.e-6
    real, parameter:: kgtoug=1.e9
    real, parameter:: SOA_frac= 1.4 ! ask Bas for this factor 

    ! to pre_calc, instead of calculate them every step  
    real           ::  rad_vol_so4ait 
    real           ::  rad_vol_so4acc 
    real           ::  rad_vol_bcai 
    real           ::  rad_vol_ocai 

    real           ::  vol_so4ait 
    real           ::  vol_so4acc 
    real           ::  vol_bcai 
    real           ::  vol_ocai 
#endif

    ! --- begin ----------------------------
    
    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units ='K')
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m')    
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'area', area, status, check_units ='m2')    
    IF_NOTOK_RETURN(status=1)

#ifdef with_m7
    ! to calc the radius of the volume of different modes 
    rad_vol_so4ait = rad_so4ait * exp(1.5*(sigmaln(2))**2)
    vol_so4ait = 4.0/3.0*pi*rad_vol_so4ait**3  ! m3

    rad_vol_so4acc = rad_so4acc * exp(1.5*(sigmaln(3))**2)
    vol_so4acc = 4.0/3.0*pi*rad_vol_so4acc**3  ! m3
    
    rad_vol_bcai = rad_bcai * exp(1.5*(sigmaln(5))**2) ! should be no difference to ln(2)
    vol_bcai = 4.0/3.0*pi*rad_vol_bcai**3  ! m3

    rad_vol_ocai   = rad_ocai * exp(1.5*(sigmaln(5))**2) ! should be no difference to ln(2)
    vol_ocai = 4.0/3.0*pi*rad_vol_ocai**3  ! m3
    ! radavmass =rad_bcai *exp(1.5*(sigmaln(5))**2)
#endif

    ! extract time parameters:
    tmid = t1 + (t2-t1)/2
    call Get( tmid, year=yy, month=mm, day=dd, hour=hh )
    
    ! time conversion parameter
    min_per_month = days_in_month( tmid ) * 24 * 60

    ! other year then loaded ? OR other month then loaded and emisfile are valid for a month
    if ( yy /= emd%year .or. ( mm /= emd%month .and. trim(emd%emisset_valid) == 'month' )) then
      
      ! Read new emission data
      call LE_Emis_EDGAR_Clear( emd, status )
      IF_NOTOK_RETURN(status=1)
      call LE_Emis_EDGAR_Read( emd, tmid%year, tmid%month, status )
      IF_NOTOK_RETURN(status=1)
            
    end if
  
    ! storage:
    allocate( em(nx,ny) )
    allocate( delta_emis(1:nz) )
    allocate( year_emis(nspec) )

    year_emis=0.0
    
    ! height distribution ?
    if ( emd%with_height_distribution ) then
      ! lookup table:
      allocate( hd_profiles(nx,ny,nz,emd%nemiscat) )
      ! storage for model layer heighs:
      allocate( hhb(0:nz) )
      ! emission top:
      emtop = emd%hdistr%heightb(emd%hdistr%nlay)
      ! init :
      hd_profiles = 0.0
      ! loop over grid cells:
      do iy = 1, ny
        do ix = 1, nx
          ! extra model layer heights:
          hhb(0)    = 0.0  ! m
          hhb(1:nz) = h_m(ix,iy,1:nz)   ! m
          ! loop over categories:
          do icat = 1, emd%nemiscat
            ! loop over model layers:
            ilast = -1
            do iz = 1, nz
              ! fraction of profile in this interval:
              call IntervalSum( emd%hdistr%heightb, emd%hdistr%fraction(icat,:), &
                                 hhb(iz-1), min(hhb(iz),emtop), &
                                 hd_profiles(ix,iy,iz,icat), ilast, status )
              IF_NOTOK_RETURN(status=1)
              ! end ?
              if ( hhb(iz) > emtop ) exit
            end do  ! model layers
          end do  ! categories
        end do ! ix
      end do ! iy
      ! clear:
      deallocate( hhb )
    end if  ! with height distribution
        
    ! set the time factor          
    julian_day = calc_DayNumber('gregorian',yy,mm,dd)
    ihour = (julian_day-1)*24+hh + 1

    ! loop over emitted tracers;
    ! tracers in the EDGAR emission list
    do itr = 1, emd%ntr
      
      ! loop over categories for this tracer:
      do icat = 1, emd%ncat(itr)
      
        ! extract global index:
        icat_glob = emd%icat_global(itr,icat)

        ! get the field for this time interval ;        
        em = emd%edg(itr,icat)%emis
        em_units = emd%edg(itr,icat)%emis_units
                
        ! loop over grid cells:
        do ix = 1, nx
          do iy = 1, ny

            ! extract current zone:
            izone = emd%izone(ix,iy)
            ! time factor:
            timefrac = emd%time_prof%profile(yy,ihour,icat_glob,izone)

            ! set temperature fraction:
            select case ( trim(emd%name_tr(itr)) )
              !~
              case ( 'CO')
                ! temperature fraction is interpolated from table using current temperature;
                ! keep first/last value if outside range:
                ilast = 0
                call Interp_Lin( emd%temp_var_CO%temper(:), emd%temp_var_CO%fraction(icat_glob,:), &
                                     tsurf(ix,iy,1), temperfrac, &
                                     ilast, status, extrapol=EXTRAPOL_CONSTANT )
                IF_NOTOK_RETURN(status=1) 
              !~
              case ( 'NMVOC' )
                ilast = 0
                call Interp_Lin( emd%temp_var_VOC%temper(:), emd%temp_var_VOC%fraction(icat_glob,:), &
                                    tsurf(ix,iy,1), temperfrac, ilast, status, extrapol=EXTRAPOL_CONSTANT )
                IF_NOTOK_RETURN(status=1) 
              !~
! #ifdef with_m7
!               case ( 'NOx', 'SO2', 'CH4', 'NH3', 'PM10', 'SF6', 'PMcoarse', 'PM2.5', 'PM2.5_bio', 'PM2.5_fossil', 'PM25_PM10', 'BC', &
!                     'CBIV_ALD2', 'CBIV_ETH','CBIV_FORM', 'CBIV_ISOP', 'CBIV_OLE', 'CBIV_PAR','CBIV_TOL', 'CBIV_XYL'  )
!                 ! no depencency yet ...
!                 temperfrac = 1.0
!               case ( 'OC' )
!                 temperfrac = 1.0 
!                 SOA = 1.4
! #else 
              case ( 'NOx', 'SO2', 'CH4', 'NH3', 'PM10', 'SF6', 'PMcoarse', 'PM2.5', 'PM2.5_bio', 'PM2.5_fossil', 'PM25_PM10', 'OC', 'BC', &
                    'CBIV_ALD2', 'CBIV_ETH','CBIV_FORM', 'CBIV_ISOP', 'CBIV_OLE', 'CBIV_PAR','CBIV_TOL', 'CBIV_XYL'  )
                ! no depencency yet ...
                temperfrac = 1.0
! #endif               
              !~
              case default
                write (gol,'("no idea if temperature factor should be defined for emitted `",a,"`")') trim(emd%name_tr(itr)); call goErr
                TRACEBACK; status=1; return
            end select
        
            ! loop over components of emitted tracer,
            ! e.g. emitted nox is distributed over no and no2, etc:
            do icomp = 1, emd%emcomp(itr)%ncomp
              ! target tracer:
              ispec = emd%emcomp(itr)%itracer(icomp)
              ! not a tracer in this run ? then skip:
              if ( ispec < 0 ) cycle
              ! explicitly skip ?
              if ( emd%skip_tracer(ispec) ) cycle
              ! extract:
              compfrac = emd%emcomp(itr)%frac(icat,izone,icomp)

              ! check for undefined values ...
              if ( compfrac < 0.0 ) then
                write (gol,'("found undefined component fraction for emission:")'); call goErr
                write (gol,'("  emission   : ",a)') trim(emd%name_tr(itr)); call goErr
                write (gol,'("  category   : ",i6)') emd%catnames(itr,icat); call goErr
                write (gol,'("  country    : ",a)') emd%zones(izone); call goErr
                write (gol,'("  component  : ",a)') specname(ispec); call goErr
                TRACEBACK; status=1; return
              end if

              ! units of emissions and target tracer:
              conversion = trim(em_units)//' ; '//trim(specunit(ispec))
              ! distribute emission (in kg) over components:
              select case ( trim(conversion) )
                !~~~~~~~~~~~~~~~~~~~~~~~~~
                ! convert to mol/min if tracer unit is volume mixing ratio:
                !~~~~~~~~~~~~~~~~~~~~~~~~~
                case ( 'kg m-2 s-1 ; ppb' )
                  ! special: NMVOC composition table is in [mol/(kg NMVOC)]
                  if ( trim(emd%name_tr(itr)) == 'NMVOC' ) then
                    !              m2                s/min 
                    convfac = area(ix,iy,1)  *  sec_per_min   ! (1/min) / (1/m2/s)
                  else
                    ! check ...
                    if ( emd%xm(itr) < 0.0 ) then
                      write (gol,'("no molecule mass defined for emitted tracer `",a,"`")') &
                             trim(emd%name_tr(itr)); call goErr
                      TRACEBACK; status=1; return
                    end if
                    !              m2          / ((kg tr)/mol)        s/min 
                    convfac = area(ix,iy,1)  /  emd%xm(itr)   *  sec_per_min   ! (mol/min) / (kg/m2/s)
                  end if
                !~
                case ( 'kg(N) m-2 s-1 ; ppb' )
                  !              m2         / ((kg N)/mol)         s/min 
                  convfac = area(ix,iy,1) /     xm_N       *   sec_per_min   ! (mol/min) / (kg(N)/m2/s)
                !~
                case ( 'kg(S) m-2 s-1 ; ppb' )
                  !              m2         / ((kg S)/mol)         s/min 
                  convfac = area(ix,iy,1) /     xm_S       *   sec_per_min  ! (mol/min) / (kg(S)/m2/s)
                !~
                case ( 'Mg/m2/year ; ppb' )
                  ! check ...
                  if ( emd%xm(itr) < 0.0 ) then
                    write (gol,'("no molecule mass defined for emitted tracer `",a,"`")') &
                           trim(emd%name_tr(itr)); call goErr
                    TRACEBACK; status=1; return
                  end if
                  !              m2            (kg/Mg)   / ((kg tr)/mol)  /   (min/year)
                  convfac = area(ix,iy,1) * kg_per_MEg /  emd%xm(itr)   /  min_per_year   ! (mol/min) / (Mg/m2/year)
                case ( 'Mg/m2/month ; ppb' )
                  ! check ...
                  if ( emd%xm(itr) < 0.0 ) then
                    write (gol,'("no molecule mass defined for emitted tracer `",a,"`")') &
                           trim(emd%name_tr(itr)); call goErr
                    TRACEBACK; status=1; return
                  end if
                  !              m2            (kg/Mg)   / ((kg tr)/mol)  /   (min/month)
                  convfac = area(ix,iy,1) * kg_per_MEg /  emd%xm(itr)   /  min_per_month   ! (mol/min) / (Mg/m2/month)
                !~
                case ( 'Mmol/m2/year ; ppb' )
                  !              m2            (mol/Mmol)   /   (min/year)
                  convfac = area(ix,iy,1) * mol_per_MEmol /  min_per_year   ! (mol/min) / (Mmol/m2/year)

                case ( 'Mmol/m2/month ; ppb' )
                  !              m2            (mol/Mmol)   /   (min/month )
                  convfac = area(ix,iy,1) * mol_per_MEmol /  min_per_month   ! (mol/min) / (Mmol/m2/month)

                !~~~~~~~~~~~~~~~~~~~~~~~~~
                ! convert to ug/min if targer tracer unit is mass concentration:
                !~~~~~~~~~~~~~~~~~~~~~~~~~
                case ( 'kg m-2 s-1 ; ug/m3' )
                  !              m2         *  (ug/kg)   *     s/min 
                  convfac = area(ix,iy,1) * ug_per_kg  *  sec_per_min  ! (ug/min) / (kg/m2/s)
                !~
                case ( 'kg(S) m-2 s-1 ; ug/m3' )
                  !              m2         *  (ug/kg)    (kg spec)/mol   / ((kg S)/mol) *    s/min 
                  convfac = area(ix,iy,1) * ug_per_kg * specmolm(ispec) /     xm_S     *  sec_per_min   ! (ug/min) / (kg/m2/s)
                !~
                case ( 'Mg/m2/year ; ug/m3' )
                  !              m2         *  (ug/Mg)    /   (min/year)
                  convfac = area(ix,iy,1) * ug_per_MEg  /  min_per_year  ! (ug/min) / (kg/m2/s)
               
                case ( 'Mg/m2/month ; ug/m3' )
                  !              m2         *  (ug/Mg)    /   (min/month)
                  convfac = area(ix,iy,1) * ug_per_MEg  /  min_per_month  ! (ug/min) / (kg/m2/s)
                  
                ! added by jianbing, for the tracer M7 sulfate aerosols
                ! convert to mlc/min 
                case ( 'kg m-2 s-1 ; mlc/cm3' )
                  if ( emd%xm(itr) < 0.0 ) then
                    write (gol,'("no molecule mass defined for emitted tracer so2")'); call goErr
                    TRACEBACK; status=1; return
                  end if
                  ! kg m-2 s-1    m2           s/min       kg/mol      mlc/mol
                  convfac = area(ix,iy,1) *  sec_per_min / emd%xm(itr) * Avog
                  ! write (gol,'("Jianbing to see emdxm is ", f10.7)') emd%xm(itr); call goErr
                  ! write (gol,'("Jianbing to see emdxm is ", f10.7)') Avog; call goErr

                !~~~~~~~~~~~~~~~~~~~~~~~~~
                ! unknown ...
                !~~~~~~~~~~~~~~~~~~~~~~~~~
                case default
                  write (gol,'("unsupported `emis;tracer` units pair : ",a)') trim(conversion); call goErr
                  write (gol,'("  emitted substance : ",a)') trim(emd%name_tr(itr)); call goErr
                  write (gol,'("  tracer component  : ",a)') trim(specname(ispec)); call goErr
                  TRACEBACK; status=1; return
              end select

              ! emission to be added:
              delta_emis = hd_profiles(ix,iy,:,icat_glob) * timefrac * temperfrac * &
                                compfrac * convfac * em(ix,iy)
              year_emis(ispec) = year_emis(ispec) + compfrac * convfac * em(ix,iy)

#ifdef with_m7
              ! if it is the M7 mode, we also need to added the number concentraion
              ! nnus    , number concentration nucle. mode soluble 
              ! nais    , number concentration aitken mode soluble 
              ! nacs    , number concentration accum. mode soluble 
              ! ncos    , number concentration coarse mode soluble 
              ! naii    , number concentration aitken mode insoluble 
              ! naci    , number concentration accum. mode insoluble 
              ! ncoi    , number concentration coarse mode insoluble 

              if ( tracer_is_m7(ispec)  ) then    

                select case ( ispec )       
                ! this is the soluable mode  
                case( ispec_bcki )  ! aitken insoluable              
                  ! mass add 
              emis_a(ix,iy,:,ispec) = emis_a(ix,iy,:,ispec) + delta_emis   ! mol/min, ug/min, mlc/min, #/min
                  ! num add 
                  emis_a(ix,iy,:, i_naii) = emis_a(ix,iy,:, i_naii) + delta_emis/( vol_bcai * dens_bc * kgtoug ) 

                case( ispec_ocki )  ! aitken insoluable  
                  ! mass add 
                  emis_a(ix,iy,:,ispec) =  emis_a(ix,iy,:,ispec) + SOA_frac * delta_emis   ! mol/min, ug/min, mlc/min, #/min
                  ! num add 
                  emis_a(ix,iy,:, i_naii) = emis_a(ix,iy,:, i_naii) + SOA_frac* delta_emis/( vol_ocai * dens_oc * kgtoug ) 

                case( ispec_ocks )  ! aitken soluable 
                  ! mass add
                  emis_a(ix,iy,:,ispec) =  emis_a(ix,iy,:,ispec) + SOA_frac * delta_emis   ! mol/min, ug/min, mlc/min, #/min
                  ! num add 
                  emis_a(ix,iy,:, i_nais) = emis_a(ix,iy,:, i_nais) + SOA_frac* delta_emis/( vol_ocai * dens_oc * kgtoug ) 
                  ! emis_a(ix,iy,:, i_nais) = emis_a(ix,iy,:, i_nais) 

                case( ispec_so4ks ) ! accumulation soluable 
                  ! mass added 
                  emis_a(ix,iy,:,ispec) =  emis_a(ix,iy,:,ispec) +  delta_emis   ! mol/min, ug/min, mlc/min, #/min
                  ! num add 
                  ! emis_a(ix,iy,:, i_nais) = emis_a(ix,iy,:, i_nais) + delta_emis  /( Avog * vol_so4ait * dens_ss * kgtoug) ! this need to be check that, whether unite is ug/m3 
                  emis_a(ix,iy,:, i_nais) = emis_a(ix,iy,:, i_nais) + delta_emis * 96.0e6 /( Avog * vol_so4ait * dens_ss * kgtoug) ! this need to be check that, whether unite is ug/m3 

                case( ispec_so4as ) ! coarse soluable       1/cm3         mlc/min         
                  ! mass add 
                  emis_a(ix,iy,:,ispec) =  emis_a(ix,iy,:,ispec) +  delta_emis   ! mol/min, ug/min, mlc/min, #/min
                  ! num add 
                  ! emis_a(ix,iy,:, i_nacs) = emis_a(ix,iy,:, i_nacs) + delta_emis  /( Avog * vol_so4acc * dens_ss * kgtoug ) ! this need to be check that
                  emis_a(ix,iy,:, i_nacs) = emis_a(ix,iy,:, i_nacs) + delta_emis * 96.0e6/( Avog * vol_so4acc * dens_ss * kgtoug ) ! this need to be check that
                      
                case default 
                  write (gol,'("Jianbing: unsupported ispc ", a )') specname(ispec); call goErr
                    TRACEBACK; status=1; return

                end select

              else 
                ! for gas but not the m7 aerosols 
                emis_a(ix,iy,:,ispec) =  emis_a(ix,iy,:,ispec) + delta_emis   
              end if 
#else 
              ! if not m7, add directly: 
              emis_a(ix,iy,:,ispec) =  emis_a(ix,iy,:,ispec) + delta_emis   ! mol/min, ug/min, mlc/min, #/min

#endif               


            ! delta_emis = u_fact*(smar_Dp(1)+smon_Dp(1) + smar_Dp(2)+smon_Dp(2))*conv 
            ! ! average aersol radius corresponding to mass:
            ! radavmass = rad_ssac * exp(1.5*(sigmaln(3))**2)  ! m
            ! ! aerosol volume:
            ! volavmass = 4.0/3.0*pi*radavmass**3  ! m3
            ! ! add mass:
            ! emis(ix,iy,1,i_ssas) = emis(ix,iy,1,i_ssas) + delta_emis
            ! ! add number:
            ! emis(ix,iy,1,i_nacs) = emis(ix,iy,1,i_nacs) + delta_emis/ (volavmass*dens_ss*kgtoug )

            end do  ! components of emitted tracer

          end do  ! iy
        end do  ! ix          

      end do  ! cats

    end do  ! tracers    
    
    ! clear:
    deallocate( em )
    deallocate( delta_emis )
    deallocate( year_emis )
  
  
    ! ok
    status = 0

  end subroutine LE_Emis_EDGAR_Setup


end module LE_Emis_EDGAR

