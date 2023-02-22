!###############################################################################
!
! LE_Emis_BaseCatIon - routines to read and work with emission file in base format
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

!---------------------------------------------------------------------------------------------------------
! This module contains data and subroutines needed for emissions.
!
! contains routines
! reademis: read base emissions for a limited number of species groups (base_emis), VOC- and time-profiles
! LE_Emis_BaseCatIon_Setup  : compute emissions (emis_a), for each species needed, for a specific time
! alloc_emis: allocate arrays for emissions
! print_emission_output: print emission output table
!
!----------------------------------------------------------------------------------
! FS Nov 17, 2005: extension to different dustributions of emission categories
! A. SNAP10
! B. SNAP_TNO = SNAP10 + subdivision for traffic (71,72,73,74,75)
! C. BTL-DFE: 77 categories
! 
! xxx Data is read from emission distribution files specified in the .rc file 
! (no fixed file names anymore).
!
!
! We have 7 types of files:
! 1. definition of country codes 
!    This file defines the countries; note that all other input must use the same countries 
!    and in the same order as specified here.
! 2. definition of emission categories
!    This file defines the emission categories; note that all other input must use the same categories
!    and in the same order as specified here.
! 3. emission distribution for each country 
! 4. definition of time profiles (distribution of emission in time, for each category)
! 5. definition of VOC-profiles (distribution of VOC-emission over different species, for each category)
! 6. definition of scenario-profiles (used to reduce or increase emissions, for each country, for each category)
! 7. definition of black carbon-profiles (distribution of black carbon, for each country, for each category)
!
! xxx ACTION: check for double definitions in country code and/or emission category ?
!
! For VOC-profiles we have different files for each chemical scheme (currently CBM4, CB99).
! For scenario-profiles, we have different files for each emitted species.
!
! Columns in these files are separated by white space.
! These files can be created and maintained using a spreadsheet (save as TAB separated).
! Each file starts with a header with an arbitrary number of lines, each header line starts with a #
! The strings BEGIN "id_string" and END "id_string" are obligatory.
! The column header (e.g. "code  category_name mon tue wed thu fri sat sun") is also obligatory.
! Category names must be between double quotes if there are spaces in a name.
!
! Example files (using only 1 for all factors):
!
!---------------------------------------------------------------------------------
! # country codes
! # date:  2005-11-17 
! # reference: ...
! # .... 
! # -----------------------
! iso3      country_name                                     
! 'ALB'		 'Albania'
! 'ARM'     'Armenia'
! ...
!
!----------------------------------------------------------------------------------
! # emission distribution BTL-DFE
! # date:  2005-11-17 
! # reference: ...
! # .... 
! # -----------------------
! code      category_name                                     
! 1100 		 "industrie: voedings- en genotmiddelen"
! 1200 		 "industrie: olie raffinaderijen"          
! 1300 		 "industrie: chemische industrie"                
! ...
! 3100 		 "verkeer: wegverkeer                               
! 3110 		 "verkeer: wegverkeer: personenauto's"             
! 3111 		 "verkeer: wegverkeer: personenauto's: benzine"
! ...
!
!----------------------------------------------------------------------------------
! # emission time profiles
! # date:  2005-11-17 
! # reference: ...
! # .... 
! # -----------------------
! BEGIN monthly distribution factors
! code  category_name                            jan feb mar apr may jun jul aug sep oct nov dec
! 1100  "industrie: voedings- en genotmiddelen"  1   1   1   1   1   1   1   1   1   1   1   1
! 1200  "industrie: olie raffinaderijen"         1   1   1   1   1   1   1   1   1   1   1   1
! ....
! END monthly distribution factors
!
! BEGIN daily distribution factors
! code  category_name                            mon tue wed thu fri sat sun
! 1100  "industrie: voedings- en genotmiddelen"  1   1   1   1   1   1   1  
! 1200  "industrie: olie raffinaderijen"         1   1   1   1   1   1   1  
! ....
! END daily distribution factors
!
! BEGIN hourly distribution factors
! code category_name                             1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
! 1100  "industrie: voedings- en genotmiddelen"  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
! 1200  "industrie: olie raffinaderijen"         1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
! ....
! END hourly distribution factors
!
!----------------------------------------------------------------------------------
! Parameters for emission distributions
!
! ndistr                                      : number of distributions
! distr_nam(1:ndistr)                         : name of distribution, e.g. SNAP10, SNAP14, DFE
!
! ncountry                                    : number of countries
! country_code(1:ncountry)                    : country code (e.g. 'ALB', 'BEL')
! country2distr(1:ncountry)                   : index of emission distribution used for each country
! note that for each country there is exactly 1 emission distribution.
!
! ncat(1:ndistr)                              : number of emission categories
! mcat                                        : maximal number of emission categories
! cat_code(1:mcat,1:ndistr)                   : code-number of each category
! cat_nam(1:mcat,1:ndistr)                    : name of each category
!
! imonthdp(1:mcat,12,1:ndistr)                : time dependency factor month, nazoeken definitie XXX
! idaydp(1:mcat,7,1:ndistr)                   : time dependency factor day, nazoeken definitie XXX
! ihourdp(1:mcat,0:24,1:ndistr)               : time dependency factor hour, nazoeken definitie XXX
!                                               timefact=imonthdp*idaydp*ihourdp
!
! nvoc                                        : number of VOC-species that are emitted
! vocprof(1:mcat,1:nvoc,1:ndistr)             : distribution factor voc-split, nazoeken definitie XXX
!
! hstack(1:mcat,1:ndistr)                     : stack height
! frac_bc(1:mcat,1:ncountry,1)                : fraction black carbon
! emisfac(1:mcat,1:ncountry,1:nemis)          : scenario factors 
!------------------------------------------------------------------------------------------------------

module LE_Emis_BaseCatIon

  use GO, only : gol, goPr, goErr

  use dims, only : nx, ny, nz, nspec, runF, outF
  use LE_Logging, only : ident1,ident2
  use LE_Logging, only : u_log, u_err

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  LE_Emis_BaseCatIon_Init, LE_Emis_BaseCatIon_Done
  public  ::  LE_Emis_BaseCatIon_Setup


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_BaseCatIon'

  ! for temporary files that are closed directly after use
  integer, parameter :: u_tmp  = 61


  ! --- local --------------------------------

  character(len=100) :: emispath

  ! emission categories
  integer                         :: mcat           ! maximal number of emission categories 
  ! (max over emission distributions)
  integer,allocatable             :: ncat(:)        ! number of emission categories
  character(len=8), allocatable   :: cat_code(:,:)  ! code-number of each emission category
  character(len=100), allocatable :: cat_nam(:,:)   ! name of each emission category 

  ! define names of emitted species and set pointers
  ! xxx should be made more flexible, also for HM
  integer, parameter  :: nemis  = 4
  character(len=5),dimension(nemis), parameter :: emisnames = &
       (/ 'Ca   ', 'Mg   ', 'K    ', 'Na   ' /)

  integer, private, parameter :: ii_ca =1
  integer, private, parameter :: ii_mg =2
  integer, private, parameter :: ii_k =3
  integer, private, parameter :: ii_na =4
  ! help array for reading emissions
  real :: emis_rd(nemis)

  !HM indices
  integer, parameter :: ii_pb=2
  integer, parameter :: ii_cd=1

  ! data type for the base emissions
  type emisdata
     character(len=3) :: emcountry_code   ! country code, e.g. 'ALB','BEL'
     integer          :: icountry         ! country number (defined by the order of the input file with country codes)
     character(len=8) :: emcat_code       ! emission category code (integer value)
     integer          :: icat             ! index for emission category (defined by the order of the input file with 
     ! category codes and names)
     character(len=1) :: source_type      ! 'a' (area source) or 'p' (point source)
     real             :: lon, lat         ! lon, lat of emission
     real             :: emis(nemis)      ! amount of emission for all emitted species (kg/y)
     integer          :: ix, iy           ! cell indices in LOTOS 
  end type emisdata

  ! arrays with emissions
  type (emisdata), allocatable :: base_emis_bc(:)        ! HM base emissions for country, location, emission category
  real                         :: emis_outside_bc(nemis) ! total of emissions that fall outside the domain

  ! emission distributions
  integer                         :: ndistr              ! number of emission distributions
  character(len=50), allocatable  :: distr_nam(:)        ! name of emission distributions (e.g. SNAP1)

  ! arrays for country code, fraction BC and scenario factors 
  integer                         :: ncountry         ! number of countries
  character(len=3), allocatable   :: country_code(:)  ! country code (e.g. 'ALB', 'BEL')
  integer, allocatable            :: country2distr(:) ! index of emission distribution for each country 
  real, allocatable               :: emisfac(:,:,:)   ! scenario factors for each category, country, emitted species

  ! stack height
  real, allocatable            :: hstack(:,:) ! hstack(1:mcat,ndistr)

  ! time dependency factors
  real, allocatable            :: imonthdp(:,:,:),idaydp(:,:,:),ihourdp(:,:,:)



contains


  !---------------------------------------------------------------------------------------


  subroutine LE_Emis_BaseCatIon_Init( rcF, rckey, status )
  
    use GO                       , only : TrcFile, ReadRc
    
    use dims , only   : runF
    use LE_IO_Tools, only  : io_read_ini, io_read_label, io_read_table1, io_check, io_err1, io_compare_strings,  & 
         io_read_table_country_cat, debopt

    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)       ::  rcF
    character(len=*), intent(in)    ::  rckey
    integer, intent(out)            ::  status
  
    ! --- const ---------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_BaseCatIon_Init'
    
    ! --- local ---------------------------------
    
    character(len=250)      :: fnam
    character(len=500)      :: colheader
    character(len=250)      :: line
    character(len=250)      :: emspecies_all
    character(len=3)        :: country_code_rd
    character(len=50)       :: distr_nam_rd
    integer                 :: iemis
    integer                 :: icountry
    integer                 :: icat
    integer                 :: irec 
    integer                 :: ios
    !integer                 :: indx1
    logical                 :: firstread
    integer                 :: ncol
    character(len=250)      :: long_name
    character(len=400)      :: country_all
    real                    :: fac_unit
    integer                 :: nlines
    integer                 :: lines_bc
    integer                 :: idistr

    ! --- begin ----------------------------------
    
    ! input path:
    call ReadRc( rcF, trim(rckey)//'.aux.dir', emispath, status )
    IF_NOTOK_RETURN(status=1)

    !------------------------------------------------------
    write(*,*) ' Read emission distributions'
    !------------------------------------------------------

    ! Emission distribution file consists of distribution abbreviation (e.g. 'SNAP10') followed by 
    ! an explanatory text
    ! note: the text is not yet used here

    ! open file, read header and count number of lines after header:
    fnam = trim(emispath)//'emis_distributions.txt'
    call io_read_ini(u_tmp,u_log,debopt,fnam,irec,nlines)

    ! allocate array
    ndistr = nlines - 3     ! two extra labels in file + column header
    if (debopt .gt. 0) write(u_log,*) ' number of emission distributions: ',ndistr
    allocate(distr_nam(ndistr))

    ! read emission distribution names
    call io_read_label('BEGIN emission distributions',u_tmp,u_err,fnam,irec)
    read(u_tmp,'(a)') line ! column header
    if (debopt .gt. 0) write(u_log,'(a)') line
    do idistr = 1,ndistr
       read(u_tmp,*,iostat = ios) distr_nam(idistr), long_name
       call io_check(ios,fnam,irec,u_tmp,u_err)
       if (debopt .gt. 0) write(u_log,'(a6,1x,a)') distr_nam(idistr), trim(long_name)
    enddo
    call io_read_label('END emission distributions',u_tmp,u_err,fnam,irec)
    if (debopt .gt. 0) write(u_log,*) ' '

    ! close file
    close (u_tmp)

    !------------------------------------------------------
    write(*,*) ' Read country codes and make indices for each country'
    !------------------------------------------------------

    ! Country file consists of iso3 country_code ('ALB') followed by long country name ('ALBANIA')
    ! note: the long name is not used here

    ! open file, read header and count number of lines after header:
    fnam = trim(emispath)//'country_code.txt'
    call io_read_ini(u_tmp,u_log,debopt,fnam,irec,nlines)

    ! allocate array
    ncountry = nlines - 3     ! two extra labels in file + 1 header line
    if (debopt .gt. 0) write(u_log,*) ' number of countries: ',ncountry
    allocate(country_code(ncountry))

    ! read country codes
    call io_read_label('BEGIN country codes',u_tmp,u_err,fnam,irec)
    read(u_tmp,'(a)') line
    if (debopt .gt. 0) write(u_log,*) line
    do icountry = 1,ncountry
       read(u_tmp,*,iostat = ios) country_code(icountry), long_name
       call io_check(ios,fnam,irec,u_tmp,u_err)
       if (debopt .gt. 0) write(u_log,'(a5,1x,a)') country_code(icountry), trim(long_name)
    enddo
    call io_read_label('END country codes',u_tmp,u_err,fnam,irec)
    if (debopt .gt. 0) write(u_log,*) ' '

    ! concatenate all country codes:
    country_all = ''
    do icountry = 1,ncountry
       country_all = trim(country_all) // ' ' // country_code(icountry)
    enddo

    ! close file
    close (u_tmp)

    !------------------------------------------------------
    write(*,*) ' Read emission distribution for each country'
    !------------------------------------------------------

    ! file consists of iso3 country_code ('ALB') followed by distribution name (e.g. 'SNAP10')

    ! open file, read header:
    fnam = trim(emispath)//'country_distribution.txt'
    call io_read_ini(u_tmp,u_log,debopt,fnam,irec)

    ! allocate array
    allocate(country2distr(ncountry))

    ! read country, distribution:
    call io_read_label('BEGIN country emission distribution',u_tmp,u_err,fnam,irec)
    read(u_tmp,'(a)') line
    if (debopt .gt. 0) write(u_log,*) line
    do icountry = 1,ncountry

       ! read country,distribution name
       read(u_tmp,*,iostat = ios) country_code_rd,distr_nam_rd
       call io_check(ios,fnam,irec,u_tmp,u_err)

       ! check country
       call io_compare_strings('country code',country_code_rd,country_code(icountry),u_err,fnam,irec)

       ! get index of emission distribution:
       call get_nam_in_list(distr_nam_rd,distr_nam,'emission distribution names',country2distr(icountry),u_err,fnam,irec)
       if (debopt .gt. 0) write(u_log,'(a5,1x,a,1x,i5)') &
            country_code(icountry), trim(distr_nam_rd),country2distr(icountry)
    enddo
    call io_read_label('END country emission distribution',u_tmp,u_err,fnam,irec)
    if (debopt .gt. 0) write(u_log,*) ' '

    ! close file
    close (u_tmp)

    !-------------------------------------------
    write(*,*) ' Read definition of emission categories'
    !-------------------------------------------

    ! Emission categories file consists of an integer code followed by the category name 

    ! Initialisation
    allocate(ncat(ndistr))
    mcat = 0

    ! loop over emission distributions
    do idistr = 1,ndistr

       ! open file, read header and count number of lines after header:
       fnam = trim(emispath)//'emis_cat_'//trim(distr_nam(idistr))//'.txt'
       call io_read_ini(u_tmp,u_log,debopt,fnam,irec,nlines)

       ! Get number of emission categories and maximal number:
       ncat(idistr) = nlines - 3   ! 2 extra labels + column header
       if (debopt .gt. 0) write(u_log,*) ' number of emission categories: ',ncat
       mcat = max(mcat,ncat(idistr))

       ! close file
       close (u_tmp)
    enddo

    ! Allocate:
    allocate(cat_code(mcat,ndistr))
    allocate(cat_nam(mcat,ndistr))
    allocate(emisfac(mcat,ncountry,nemis))
    allocate(hstack(mcat,ndistr))


    do idistr = 1,ndistr

       ! open file, read header and count number of lines after header:
       fnam = trim(emispath)//'emis_cat_'//trim(distr_nam(idistr))//'.txt'
       call io_read_ini(u_tmp,u_log,-1,fnam,irec)

       ! read category definitions
       call io_read_label('BEGIN emission categories',u_tmp,u_err,fnam,irec)
       read(u_tmp,'(a)') line ! column header
       if (debopt .gt. 0) write(u_log,'(a)') line
       do icat = 1,ncat(idistr)
          read(u_tmp,*,iostat = ios) cat_code(icat,idistr), cat_nam(icat,idistr)
          call io_check(ios,fnam,irec,u_tmp,u_err)
          if (debopt .gt. 0) write(u_log,'(a,1x,a)') cat_code(icat,idistr), trim(cat_nam(icat,idistr))
       enddo
       call io_read_label('END emission categories',u_tmp,u_err,fnam,irec)
       if (debopt .gt. 0) write(u_log,*) ' '

       ! close file
       close (u_tmp)
    enddo


    ! generate list
    if (debopt > 0) then
       write(u_log,*) ' '
       write(u_log,*) '------------------------------------------------------ '
       write(u_log,*) ' list of country/category'
       write(u_log,*) '------------------------------------------------------ '
       do icountry = 1,ncountry
          idistr = country2distr(icountry)
          do icat = 1,ncat(idistr)
             write(u_log,'(a4,1x,a8,1x,a)') country_code(icountry),cat_code(icat,idistr),trim(cat_nam(icat,idistr))
          enddo
       enddo
    endif

    !--------------------------------------------------
    !write(*,*) ' Read scenario factors'
    !--------------------------------------------------
    ! xxx scenario factors not for other type of runs ?
    ! xxx emisfac also for HM ?
    ! scenario factors for ship emissions ?

    !if (runF%scheme .ne. 'none' .or. &
    !     runF%do_aerosol         .or. &
    !     runF%do_prim            .or. &
    !     runF%do_basecation      .or. &
    !     runF%do_secorg          .or. &
    !     runF%do_sulphur) then

       ! Read emission factors for each emitted species (for each country and emission category);
       ! emisson factors are used to reduce (or increase) the emissions for a specific country
       ! and/or emission category

       ! read file name:
       call ReadRc( rcF, trim(rckey)//'.emisfac.file', fnam, status )
       IF_NOTOK_RETURN(status=1)
       !
       write(*,*) ' Read scenario factors, file = emisfac.txt'
       if (debopt > 0) write(u_log,*) ' reading emission factors from file ',trim(fnam)

       ! Construct list of emitted species names:
       emspecies_all = ''
       do iemis = 1,nemis
          emspecies_all = trim(emspecies_all) // ' ' // trim(emisnames(iemis))
       enddo

       ! read table with emission factors; 
       ! row    -> country/category
       ! column -> emitted species
       firstread = .true.
       call io_read_table_country_cat('emission factors','[-]',ndistr,mcat,ncat,ncountry,nemis,country2distr, &
            'country code category_name '//trim(emspecies_all),country_code,cat_code,cat_nam, &
            firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
            emisfac)

       ! Write info 
       write(u_log,*) ' '
       write(u_log,*) '-------------------------------------------------------------------------- '
       write(u_log,'(a6,1x,a10,40x,2(1x,a8))') 'code','category','minimum','maximum'
       write(u_log,*) '-------------------------------------------------------------------------- '
       do iemis = 1,nemis
          do idistr = 1,ndistr
             write(u_log,*) ' '
             write(u_log,*) '------------------------------------------------------------------------ '
             write(u_log,*) ' emission factors ',trim(emisnames(iemis)),' for ',trim(distr_nam(idistr))
             write(u_log,*) ' '
             write(u_log,*) '------------------------------------------------------------------------ '
             do icat = 1,ncat(idistr)
                write(u_log,'(a,1x,a50,2(1x,f8.2))') cat_code(icat,idistr),cat_nam(icat,idistr), &
                     minval(emisfac(icat,:,iemis)), &
                     maxval(emisfac(icat,:,iemis))
             enddo
          enddo
       enddo

       ! Close file:
       close(u_tmp)
    !endif

    !--------------------------------------------------------
    write(*,*) '  Read stack heights'
    !--------------------------------------------------------

    do idistr = 1,ndistr

       ! read table with stack heights; 
       ! row    -> category
       ! column -> stack height (1 column)
       !
       ! read file name:
       call ReadRc( rcF, trim(rckey)//'.stack_height.file', fnam, status )
       IF_NOTOK_RETURN(status=1)
       !
       write(*,*) '    ',trim(fnam)
       !
       firstread = .true.
       call io_read_table1('stack heights for point sources','[m]',ncat(idistr),1, &
            'code category_name stack_height',cat_code(:,idistr),cat_nam(:,idistr), &
            firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
            hstack(:,idistr), status )
       IF_NOTOK_RETURN(status=1)

       ! xxx hstack traffic = 50 m ?? point sources only

       ! Close file
       close(u_tmp)
    enddo

    !--------------------------------------------------------
    !write(*,*) '  Read time profiles'
    !--------------------------------------------------------

    ! Allocate arrays:
    allocate(imonthdp(mcat,12,ndistr))
    allocate(idaydp(mcat,7,ndistr))
    allocate(ihourdp(mcat,0:24,ndistr))

    do idistr = 1,ndistr

       ! read table with 12 monthly factors
       ! row    -> category
       ! column -> month
       fnam      = trim(emispath)//'time_var_emis_'//trim(distr_nam(idistr))//'.txt'
       write(*,*) '  Read time profiles, file = time_var_emis_'//trim(distr_nam(idistr))//'.txt'
       firstread = .true.
       ncol      = 12
       colheader = 'code category_name jan feb mar apr may jun jul aug sep oct nov dec'
       call io_read_table1('monthly time profiles emissions','[-]',ncat(idistr),ncol, &
            colheader,cat_code(:,idistr),cat_nam(:,idistr), &
            firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
            imonthdp(:,:,idistr), status )
       IF_NOTOK_RETURN(status=1)
       firstread = .false.

       ! read table with 7 daily factors
       ! row    -> category
       ! column -> day in week
       ncol      = 7
       colheader = 'code category_name mon tue wed thu fri sat sun'
       call io_read_table1('daily time profiles emissions','[-]',ncat(idistr),ncol, &
            colheader,cat_code(:,idistr),cat_nam(:,idistr), &
            firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
            idaydp(:,:,idistr), status )
       IF_NOTOK_RETURN(status=1)

       ! read table with 24 hourly factors
       ! row    -> category
       ! column -> hour
       ncol      = 24
       colheader = 'code category_name 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24'
       call io_read_table1('hourly time profiles emissions','[-]',ncat(idistr),ncol, &
            colheader,cat_code(:,idistr),cat_nam(:,idistr), &
            firstread,fnam,irec,u_tmp,u_log,u_err,debopt, &
            ihourdp(:,1:ncol,idistr), status )
       IF_NOTOK_RETURN(status=1)
       ! set hour 0 equal to hour 24
       ihourdp(:,0,idistr)=ihourdp(:,24,idistr)

       if (debopt > 1) then
          write(*,*) 'em/monthdp 7 ',imonthdp(7,:,idistr)
          write(*,*) 'em/idaydp 7 ',idaydp(7,:,idistr)
          write(*,*) 'em/ihourdp 7 ',ihourdp(7,0:24,idistr)
       endif

       ! Close file
       close(u_tmp)
    enddo


    !if (runF%do_basecation) then
       !--------------------------------------------------------
       write(*,*) '  Read emissions for basecations'
       !--------------------------------------------------------

       ! input is in ton/year
       fac_unit = 1000.0

       ! read file name:
       call ReadRc( rcF, trim(rckey)//'.base.file', fnam, status )
       IF_NOTOK_RETURN(status=1)
       !
       ! read emissions for 4 BCs
       call read_base_emis_count(u_tmp,fnam,debopt,lines_bc,status)
       IF_NOTOK_RETURN(status=1)
       allocate(base_emis_bc(lines_bc))
       call read_base_emis(u_tmp,fnam,lines_bc,4,fac_unit,debopt,base_emis_bc,emis_outside_bc,status)   
       IF_NOTOK_RETURN(status=1)
    !endif


  end subroutine LE_Emis_BaseCatIon_Init
  
  
  ! ***
  
  
  subroutine LE_Emis_BaseCatIon_Done( status)
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_BaseCatIon_Done'

    ! --- begin ----------------------------------
    
    ! ok
    status = 0
    
  end subroutine LE_Emis_BaseCatIon_Done
  
  
  ! ***

  !---------------------------------------------------------------------------------------
  subroutine read_base_emis_count(u_file,fnam,debopt,lines,status)

    !
    ! count number of base emissions records
    !

    use LE_IO_Tools, only : io_read_ini
    use LE_Grid, only : ugg
    
    ! Arguments:
    integer,          intent(in)              :: u_file   ! unit number file 
    character(len=*), intent(in)              :: fnam     ! file name
    integer,          intent(in)              :: debopt   ! debug option , if >0 -> screen output
    integer,          intent(out)             :: lines    ! number of emission records
    integer,          intent(out)             :: status    
    
    ! --- const ---------------------------------
    character(len=*), parameter ::  rname = mname//'/read_base_emis_count'


    ! Local 
    integer                 :: irec
    character(len=3)        :: country_code_rd
    character(len=8)        :: cat_code_rd
    character(len=1)        :: type_rd
    real                    :: lat_rd, lon_rd
    logical                 :: inflag
    ! xxx IOSTATUS not yet checked
    !integer                 :: ios  ! IO-status

    ! open file, read header:
    call io_read_ini(u_file,u_log,debopt,fnam,irec)

    ! Count all emission lines that lie within the domain
    write(*,*) ' counting number of emission records ...'
    write(*,*) ' '
    lines = 0
    do
       read (u_file,*,end=99) country_code_rd, cat_code_rd, type_rd, lon_rd, lat_rd

       call ugg%InDomain( lon_rd, lat_rd, inflag, status )
       IF_NOTOK_RETURN(status=1)            

       if (inflag) lines = lines + 1
       ! if (mod(lines,100000) .eq. 0) write(*,'(i8,a5)') lines,' ... '
    enddo
99  continue
    close (u_file)

    write(*,*) ' number of emission records: ',lines

  end subroutine read_base_emis_count

  !---------------------------------------------------------------------------------------
  subroutine read_base_emis(u_file,fnam,lines,ncol,fac_unit,debopt,base_emis_local,emis_outside_local,status)

    !
    ! Read base emissions (country_code,cat_code,type,lon,lat,emis(1:ncol)
    !

    use LE_Grid, only : ugg    
    
    use LE_IO_Tools, only : io_read_ini

    ! Arguments:
    integer,          intent(in)              :: u_file   ! unit number file 
    character(len=*), intent(in)              :: fnam     ! file name
    integer,          intent(in)              :: lines    ! number of emission records
    integer,          intent(in)              :: ncol     ! number of columns with emissions (number of emitted species)
    real   ,          intent(in)              :: fac_unit ! factor to convert input units to kg/year 
    ! (LOTOSEUROS uses kg/year for emissions)
    integer,          intent(in)              :: debopt   ! debug option , if >0 -> screen output
    type (emisdata),  intent(out)             :: base_emis_local(lines)
    real,             intent(out)             :: emis_outside_local(ncol)
    integer,          intent(out)             :: status
    
    ! --- const ---------------------------------
    character(len=*), parameter ::  rname = mname//'/read_base_emis'
    

    ! Local 
    integer                 :: irec
    character(len=3)        :: country_code_rd
    character(len=8)        :: cat_code_rd
    character(len=1)        :: type_rd
    real                    :: lat_rd, lon_rd
    integer                 :: ilin, ix_rd, iy_rd
    logical                 :: inflag
    integer                 :: idistr
    integer                 :: icountry
    integer                 :: icat
    ! xxx IOSTATUS not yet checked
    !integer                 :: ios  ! IO-status

    real :: dum1, dum2, dum3

    ! open file and read header
    call io_read_ini(u_file,u_log,-1,fnam,irec)
    write(*,*) ' loading emission records ...'
    if (debopt>0) write(*,*) ' loading emission records ...'

    ! read until end of file
    ilin = 0
    emis_outside_local = 0.0
    do 

       ! read emissions:
       !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       ! the following generates a segmentation fault. Why ???????
       ! read (u_file,*,end=999) country_code_rd, cat_code_rd, type_rd, lon_rd, lat_rd, emis_rd(1:ncol)
       ! xxx moet anders aantal kolommen doorgeven 
       !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       
       ! Needs attention! if statement does not allow simultaneous runs! XXXX: action
       read (u_file,*,end=999) country_code_rd, cat_code_rd, type_rd, lon_rd, lat_rd, dum1, dum2, dum3, emis_rd(1:4)
       !xxx read extra column with distribution (e.g. SNAP14 or DFE) (not needed but more clear)
       irec = irec + 1

       ! check resolution ...
       if ( irec == 1 ) then
          if ( .not. ugg%type == 'carthesian-regular' ) then
            print *, 'ERROR: check emission file for basecat-ions to fit in rungrid'
            stop
          else  
            if ( (ugg%dlon < 0.50) .or. (ugg%dlat < 0.25) ) then
               ! * base_emis.txt
               if ( (lon_rd == 32.3750) .and. (lat_rd == 53.1875) ) then
                  stop 'ERROR - I think you use base_emis.txt in a zoom run; try to use a zoomed emission file!'
                  ! * rap_uba_20_12_2006.txt
               else if ( (lon_rd == 19.56250) .and. (lat_rd == 42.53125) ) then
                  ! ok, this is a zoomed emission file
               else if ( (lon_rd == 19.55) .and. (lat_rd == 40.73) ) then
                  ! ok, this is the macc emission file
                  ! * other ...
               else
                  print *, 'WARNING - check on emission resolutions during zoom run not implemented for this file ...'
                  print *, 'WARNING -   first point in emission file : ', lon_rd, lat_rd
                  ! stop 'ERROR - please implement the check for this file too !'
               end if
            end if
          end if
       end if

       
       call ugg%InDomain(lon_rd,lat_rd,inflag,status)
       IF_NOTOK_RETURN(status=1)       
            

       if (inflag) then

          ! get the country number icountry (if not found -> stop)
          call get_nam_in_list(country_code_rd,country_code,'country codes',icountry,u_err,fnam,irec)
          idistr = country2distr(icountry)

          ! get the emission category number icat (if not found -> stop)
          call get_nam_in_list(cat_code_rd,cat_code(:,idistr),'emission category codes',icat,u_err,fnam,irec)
          
          ! convert to cell location:
          !call ... lon_rd, lat_rd, ix_rd, iy_rd, ...
          ! dummy:
          ix_rd = -999
          iy_rd = -999
          ! not yet ...
          write (gol,'("not implemented yet")'); call goErr
          TRACEBACK; status=1; return

          ! put data into data/structure base_emis_local  
          ilin = ilin + 1
          if (mod(ilin,100000) .eq. 0) write(*,'(i8,a5)') ilin,' ... '
          base_emis_local(ilin)%emcountry_code = country_code_rd
          base_emis_local(ilin)%icountry     = icountry
          base_emis_local(ilin)%icat         = icat
          base_emis_local(ilin)%source_type  = type_rd
          base_emis_local(ilin)%lon          = lon_rd
          base_emis_local(ilin)%lat          = lat_rd
          base_emis_local(ilin)%ix           = ix_rd
          base_emis_local(ilin)%iy           = iy_rd

          ! multiply emissions by fac_unit and emission factor
          base_emis_local(ilin)%emis(1:ncol) = emis_rd(1:ncol)*fac_unit*emisfac(icat,icountry,1:ncol)
       else

          ! keep total of emission outside domain
          ! note: no emission factors used for emissions outside domain
          emis_outside_local(1:ncol) = emis_outside_local(1:ncol) + emis_rd(1:ncol)*fac_unit
       endif
    enddo
999 continue ! end of file
    close (u_file)

  end subroutine Read_base_emis


  !---------------------------------------------------------------------------------------


  subroutine LE_Emis_BaseCatIon_Setup( emis_a, yy,mm,dd,hh, status )

    use dims, only  : nx, ny, nz,zenith,nspec
    use LE_Time, only : local_time2
    use indices
    
    ! point to meteo data:
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------------

    real, intent(inout)       ::  emis_a(nx,ny,nz,nspec)
    integer, intent(in)       ::  yy, mm, dd, hh
    integer, intent(out)      ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Emis_Dust_Setup'

    ! --- local ----------------------------------

    integer :: yyh, idistr
    character(len=1) :: type

    integer ::  mmh, ddh, hhh, iday, i, k, icat, ix, iy, iz, icountry
    integer ::  itr, ispec
    real    :: fact, hstack1
    !real    :: hemis
    
    ! meteo data:
    real, pointer               ::  tsurf(:,:,:)   ! (lon,lat,1)
    real, pointer               ::    h_m(:,:,:)   ! (lon,lat,lev)

    ! --- begin ----------------------------------
    
    call LE_Data_GetPointer( 'tsurf', tsurf, status, check_units ='K')
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m')
    IF_NOTOK_RETURN(status=1)

    if (.NOT.outF%suppress) print *, ident2,'<constructing emissions>'
    
    ! dummy value to avoid warnings ...
    iz = -999

    !if (runF%do_basecation) then
       do i=1,size(base_emis_bc)
          ix = base_emis_bc(i)%ix 
          iy = base_emis_bc(i)%iy 
          icat     = base_emis_bc(i)%icat
          icountry = base_emis_bc(i)%icountry
          type     = base_emis_bc(i)%source_type
          idistr   = country2distr(icountry)

          ! compute actual time for emission based on the country code
          ! get day-of-the-week number
          call local_time2(base_emis_bc(i)%emcountry_code, yy, mm, dd, hh, yyh, mmh, ddh, hhh, iday, status)
          IF_NOTOK_RETURN(status=1)
          if (type == 'A' .OR. type == 'a') then
             ! area source
             iz      = 1
             hstack1 = 0.0
          else if (type == 'P' .OR. type == 'p') then
             ! point source
             hstack1 = hstack(icat,idistr)
          else
             write(u_err,*) 'source type >>>', trim(type), ' <<<< not supported '
             stop
          endif
          ! determine the layer in which the emission takes place
          if ( hstack1 > 0.0 ) then
             do k = 1, nz
                if ( hstack1 < h_m(ix,iy,k) ) then
                   iz = k
                   exit
                endif
             end do
          end if

          ! compute time factor
          fact=imonthdp(icat,mmh,idistr)*&
               idaydp(icat,iday,idistr)*&
               ihourdp(icat,hhh+1,idistr) ! at local hour 00 (hhh=0), the time profile for the first hour is valid

          do itr = 1, n_basecation
            ispec = ispecs_basecation(itr)
            select case ( ispec )
              case (ispec_Na_f2); emis_a(ix,iy,iz,ispec) = emis_a(ix,iy,iz,ispec) + 0.5* base_emis_bc(i)%emis(ii_na)*fact/8760.0*1.0e9/60.0
              case (ispec_Na_c3); emis_a(ix,iy,iz,ispec) = emis_a(ix,iy,iz,ispec) + 0.5* base_emis_bc(i)%emis(ii_na)*fact/8760.0*1.0e9/60.0
              case (ispec_Ca_f) ; emis_a(ix,iy,iz,ispec) = emis_a(ix,iy,iz,ispec) + 0.5* base_emis_bc(i)%emis(ii_ca)*fact/8760.0*1.0e9/60.0
              case (ispec_Ca_c) ; emis_a(ix,iy,iz,ispec) = emis_a(ix,iy,iz,ispec) + 0.5* base_emis_bc(i)%emis(ii_ca)*fact/8760.0*1.0e9/60.0
              case (ispec_Mg_f) ; emis_a(ix,iy,iz,ispec) = emis_a(ix,iy,iz,ispec) + 0.5* base_emis_bc(i)%emis(ii_mg)*fact/8760.0*1.0e9/60.0
              case (ispec_Mg_c) ; emis_a(ix,iy,iz,ispec) = emis_a(ix,iy,iz,ispec) + 0.5* base_emis_bc(i)%emis(ii_mg)*fact/8760.0*1.0e9/60.0
              case (ispec_K_f ) ; emis_a(ix,iy,iz,ispec) = emis_a(ix,iy,iz,ispec) + 0.5* base_emis_bc(i)%emis(ii_k )*fact/8760.0*1.0e9/60.0
              case (ispec_K_c ) ; emis_a(ix,iy,iz,ispec) = emis_a(ix,iy,iz,ispec) + 0.5* base_emis_bc(i)%emis(ii_k )*fact/8760.0*1.0e9/60.0
              case default
                write (gol,'("unsupported tracer index ",i6)') ispec; call goErr
                TRACEBACK; status=1; return
            end select
          end do

       end do

    !endif

    if (.NOT.outF%suppress) print *, ident2,'  <finished processing emissions>'

    ! ok
    status = 0
    
  end subroutine LE_Emis_BaseCatIon_Setup

  !---------------------------------------------------------------------------------------

  !subroutine alloc_emis
  !
  !implicit none 
  !
  !!only needed for EMEp emissions not here, this is quick and dirty
  !!allocate(height_frac(nz) )
  !
  !print *, 'nothing to be allocated for emissions'
  !
  !return
  !
  !end subroutine alloc_emis

  !---------------------------------------------------------------------------------------
  subroutine get_num_in_list(num,list,listnam,indx,u_err,fnam,irec)

    ! get the index of a given number in a list of numbers; stop if the number is not found.

    use LE_IO_Tools, only : io_err1

    implicit none

    ! arguments
    integer, intent(in)                    :: num     ! input number
    integer, intent(in)                    :: list(:) ! list of numbers
    character(len=*), intent(in)           :: listnam ! name of list (used for error message)
    integer, intent(out)                   :: indx    ! index, such that num = list(indx)
    integer, intent(in)                    :: u_err   ! unit number error file
    character(len=*), optional, intent(in) :: fnam    ! file name (used for error message)
    integer, optional, intent(in)          :: irec    ! record number (used for error message)

    ! local
    logical :: found

    ! Initialisation:
    found = .false.
    indx  = 0

    ! Loop through list:
    do while (.not. found .and. indx .lt. size(list))
       indx  = indx + 1
       found = (num .eq. list(indx))
    enddo

    ! error message, stop:
    if (.not. found) then
       if (present(fnam) .and. present(irec)) call io_err1(fnam,irec,u_err)
       write(u_err,*) ' '
       write(u_err,*) ' cannot find ',num,' in list with ',trim(listnam)
       write(u_err,*) ' '
       stop 
    endif

  end subroutine get_num_in_list

  !---------------------------------------------------------------------------------------
  subroutine get_nam_in_list(nam,list,listnam,indx,u_err,fnam,irec)

    ! get the index of a given name in a list of names; stop if the name is not found.

    use LE_IO_Tools, only : io_err1, io_compress

    implicit none

    ! arguments
    character(len=*), intent(in)           :: nam     ! input name
    character(len=*), intent(in)           :: list(:) ! list of names
    character(len=*), intent(in)           :: listnam ! name of list (used for error message)
    integer, intent(out)                   :: indx    ! index, such that nam = list(indx)
    integer, intent(in)                    :: u_err   ! unit number error file
    character(len=*), optional, intent(in) :: fnam    ! file name (used for error message)
    integer, optional, intent(in)          :: irec    ! record number (used for error message)

    ! local
    logical :: found

    ! Initialisation:
    found = .false.
    indx  = 0

    ! Loop through list:
    do while (.not. found .and. indx .lt. size(list))
       indx  = indx + 1
       !found = (adjustl(nam) .eq. adjustl(list(indx)))
       found = (io_compress(nam) .eq. io_compress(list(indx)))
    enddo

    ! error message, stop:
    if (.not. found) then
       if (present(fnam) .and. present(irec)) call io_err1(fnam,irec,u_err)
       write(u_err,*) ' '
       write(u_err,*) ' cannot find ',trim(nam),' in list with ',trim(listnam)
       write(u_err,*) ' '
       stop 
    endif

  end subroutine get_nam_in_list

  !---------------------------------------------------------------------------------------

end module LE_Emis_BaseCatIon

