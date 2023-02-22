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

module LE_Emis

  use GO, only : gol, goErr, goPr

  use LE_Emis_data   , only : T_Emis_Data
  use LE_Emis_TNO    , only : T_Emis_TNO
  use LE_Emis_EDGAR  , only : T_Emis_EDGAR
  use LE_Emis_MEGAN  , only : T_Emis_MEGAN
#ifdef with_pollen  
  use LE_Emis_Pollen , only : T_Emis_Pollen
#endif
  use LE_Emis_OPS    , only : T_Emis_OPS
  use LE_Emis_Point  , only : T_Emis_Point
  use LE_Emis_GFAS   , only : T_Emis_GFAS

  implicit none


  ! --- in/out --------------------------------

  private

  public  ::  LE_Emis_Init, LE_Emis_Done
  public  ::  LE_Emis_Setup, LE_Emis_Add

  public  ::  emis_a_bio

  public  ::  emis_set, max_emis

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis'

  ! maximum number of supported emissions:
  integer, parameter  ::  max_emis = 15


  ! --- types --------------------------------

  ! info on emisson data:
  type T_Emis_Set
    ! label:
    character(len=64)     ::  name
    ! filetype key:
    character(len=32)     ::  filetype
    ! timer:
    integer               ::  itim
    ! emis data:
    type(T_Emis_Data)     ::  emd
    ! tno data:
    type(T_Emis_TNO)      ::  emt
    ! EDGAR data:
    type(T_Emis_EDGAR)    ::  edgar
    ! MEGAN data:
    type(T_Emis_MEGAN)    ::  megan
#ifdef with_pollen     
    ! pollen data:
    type(T_Emis_Pollen)   ::  pollen
#endif  
    ! ops-mis-data
    type(T_Emis_OPS)      ::  emo
    ! point emis:
    type(T_Emis_Point)    ::  point
    ! GFAS emissions:
    type(T_Emis_GFAS)     ::  gfas
  end type T_Emis_Set



  ! --- local --------------------------------

  ! emission sets and file types:
  integer               ::  n_emis
  type(T_Emis_Set)      ::  emis_set(max_emis)

  ! total bio emis, for output ..
  real, allocatable     ::  emis_a_bio(:,:,:,:)

  ! flags:
  logical               ::  with_pom_to_vbs


contains


  ! ===============================================================


  subroutine LE_Emis_Init( rcF, t, status )

    use GO                       , only : TDate
    use GO                       , only : TrcFile, ReadRc
    use GO                       , only : goSplitString
    use GO                       , only : GO_Timer_Def
    use GO                       , only : TDate
    use dims                     , only : runF
    use dims                     , only : nx,ny,nz
    use Indices                  , only : nspec
    use Indices                  , only : n_seasalt, n_dust, n_basecation, n_hm
    use Indices                  , only : n_vbs_cg
    use LE_Emis_EDGAR            , only : LE_Emis_EDGAR_Init
    use LE_Emis_Bio              , only : LE_Emis_Bio_Init
    use LE_Emis_Bio_NO           , only : LE_Emis_Bio_NO_Init
    use LE_Emis_MEGAN            , only : LE_Emis_MEGAN_Init
#ifdef with_pollen    
    use LE_Emis_Pollen           , only : LE_Emis_Pollen_Init
#endif
    use LE_Emis_SeaSalt          , only : LE_Emis_SeaSalt_Init
    use LE_Emis_Dust_Wind        , only : LE_Emis_Dust_Wind_Init
    use LE_Emis_Dust_Resuspension, only : LE_Emis_Dust_Resuspension_Init
    use LE_Emis_Dust_Agriculture , only : init_agriacti_dust
    use LE_Emis_BaseCatIon       , only : LE_Emis_BaseCatIon_Init
    use LE_Emis_OPS              , only : LE_Emis_OPS_Init
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    type(TRcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Init'

    ! --- local -------------------------------

    character(len=512)          ::  list
    integer                     ::  i_emis
    character(len=64)           ::  rckey
    type(TDate), intent(in)     ::  t       ! start time
    character(len=64)           ::  emis_names(max_emis)

    ! --- begin -------------------------------

    ! info ...
    write (gol,'("initialise emissions")'); call goPr

    ! enable data:
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'p', status )
    IF_NOTOK_RETURN(status=1)

    ! read set names:
    call ReadRc( rcF, 'le.emis.sets', list, status )
    IF_NOTOK_RETURN(status=1)

    ! extract parts:
    call goSplitString( list, n_emis, emis_names, status )
    IF_NOTOK_RETURN(status=1)

    ! loop over emissions:
    do i_emis = 1, n_emis

      ! store name:
      emis_set(i_emis)%name = trim(emis_names(i_emis))

      ! define timer:
      call GO_Timer_Def( emis_set(i_emis)%itim, 'emis_'//trim(emis_set(i_emis)%name), status )
      IF_NOTOK_RETURN(status=1)

      ! first part of rcfiles:
      write (rckey,'("le.emis.",a)') trim(emis_set(i_emis)%name)

      ! read filetype, store for later use:
      call ReadRc( rcF, trim(rckey)//'.filetype', emis_set(i_emis)%filetype, status )
      IF_NOTOK_RETURN(status=1)

      ! setup data:
      select case ( trim(emis_set(i_emis)%filetype) )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'data' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! init data 
          call emis_set(i_emis)%emd%Init( rcf, rckey, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'tno' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! read tables and base emis from .nc file:
          call emis_set(i_emis)%emt%Init( trim(emis_set(i_emis)%name), &
                                            rcF, rckey, t, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'edgar' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! read tables and base emis from .nc file:
          call LE_Emis_EDGAR_Init( emis_set(i_emis)%edgar, &
                                       trim(emis_set(i_emis)%name), &
                                       rcF, rckey, t, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'bio' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! setup biogenic emissions:
          call LE_Emis_Bio_Init( rcF, status  )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'bio-no' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! setup biogenic emissions:
          call LE_Emis_Bio_NO_Init( status  )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'megan' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! setup MEGAN biogenic emissions:
          call LE_Emis_MEGAN_Init( emis_set(i_emis)%megan, &
                                   trim(emis_set(i_emis)%name), &
                                   rcF, rckey, t, status )
          IF_NOTOK_RETURN(status=1)

#ifdef with_pollen
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'silam-pollen' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
          ! info ...
          write (gol,'("initialise pollen emissions")'); call goPr
          
          call LE_Emis_Pollen_Init( emis_set(i_emis)%pollen, &
                                       trim(emis_set(i_emis)%name), &
                                       rcF, rckey, t, status )
          IF_NOTOK_RETURN(status=1)
#endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'point' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
          ! info ...
          write (gol,'("initialise point emissions")'); call goPr
          
          call emis_set(i_emis)%point%Init( trim(emis_set(i_emis)%name), &
                                              rcF, rckey, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'gfas' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
          ! info ...
          write (gol,'("initialise GFAS emissions")'); call goPr
          
          call emis_set(i_emis)%gfas%Init( trim(emis_set(i_emis)%name), &
                                              rcF, rckey, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'sea_salt' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! run includes seasalt ?
          if ( n_seasalt > 0 ) then

            ! Initialise the sea salt emission function
            call LE_Emis_SeaSalt_Init( status )
            IF_NOTOK_RETURN(status=1)

          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-wind' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! dust included ?
          if ( (n_dust > 0) .or. (n_basecation > 0) ) then

            write (gol,'("initialise wind blown dust emissions")'); call goPr

            call LE_Emis_Dust_Wind_Init( rcF, status )
            IF_NOTOK_RETURN(status=1)

          endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-resuspension' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! dust included ?
          if ( (n_dust > 0) .or. (n_basecation > 0) ) then

            write (gol,'("initialise road dust resuspension emissions")'); call goPr

            call LE_Emis_Dust_Resuspension_Init( rcF, rckey, status )
            IF_NOTOK_RETURN(status=1)

          endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-agriculture' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! dust included ?
          if ( (n_dust > 0) .or. (n_basecation > 0) ) then

            write (gol,'("initialise agricultural dust emissions")'); call goPr

            call init_agriacti_dust( status )
            IF_NOTOK_RETURN(status=1)

          endif
  
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'hm_natural' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          ! Nothing to be done
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'basecation' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! base cat ions included ?
          if ( n_basecation > 0 ) then

            ! init base cat ion emissions:
            call LE_Emis_BaseCatIon_Init( rcF, rckey, status  )
            IF_NOTOK_RETURN(status=1)

          endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'ops' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          ! init OPS-emission-data
          call LE_Emis_OPS_Init( emis_set(i_emis)%emo, &
                                  trim(emis_set(i_emis)%name), &
                                  rcF, rckey, t, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! error ...
          write (gol,'("found unsupported emis data type : ",a)') trim(emis_set(i_emis)%filetype); call goErr
          TRACEBACK; status=1; return

      end select


      !--------------------------------------------------------
      ! done
      !--------------------------------------------------------

    end do  ! emission sets

    ! temporary storage for total bio emis:
    allocate( emis_a_bio(nx,ny,nz,nspec), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! vbs settings:
    if ( n_vbs_cg > 0 ) then
      ! re-distribute pom to vbs tracers ?
      call ReadRc( rcF, 'le.vbs.pom_to_vbs', with_pom_to_vbs, status )
      IF_NOTOK_RETURN(status=1)
      ! info ...
      write (gol,'("re-distribute pom emisisons to vbs tracers : ",l1)') with_pom_to_vbs; call goPr
    else
      ! no vbs, so no need to re-distribute:
      with_pom_to_vbs = .false.
    end if

    ! ok
    status = 0

  end subroutine LE_Emis_Init


  ! ***


  subroutine LE_Emis_Done( status )

    use Indices            , only : n_hm, n_basecation, n_seasalt
    use LE_Emis_EDGAR      , only : LE_Emis_EDGAR_Done
    use LE_Emis_Bio        , only : LE_Emis_Bio_Done
    use LE_Emis_Bio_NO     , only : LE_Emis_Bio_NO_Done
    use LE_Emis_MEGAN      , only : LE_Emis_MEGAN_Done
#ifdef with_pollen    
    use LE_Emis_Pollen     , only : LE_Emis_Pollen_Done
#endif
    use LE_Emis_Dust_Wind  , only : LE_Emis_Dust_Wind_Done
    use LE_Emis_SeaSalt    , only : LE_Emis_SeaSalt_Done
    use LE_Emis_BaseCatIon , only : LE_Emis_BaseCatIon_Done
    use LE_Emis_OPS        , only : LE_Emis_OPS_Done

    ! --- in/out ------------------------------

    integer, intent(out)  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Done'

    ! --- local -------------------------------

    integer                ::  i_emis

    ! --- begin -------------------------------

    ! loop over emissions:
    do i_emis = 1, n_emis

      ! setup data:
      select case ( trim(emis_set(i_emis)%filetype) )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'data' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with data 
          call emis_set(i_emis)%emd%Done( status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'tno' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with emission data:
          call emis_set(i_emis)%emt%Done( status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'edgar' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with emission data:
          call LE_Emis_EDGAR_Done( emis_set(i_emis)%edgar, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'bio' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with biogenic emissions:
          call LE_Emis_Bio_Done( status  )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'bio-no' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with biogenic emissions:
          call LE_Emis_Bio_NO_Done( status  )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'megan' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with MEGAN biogenic emissions:
          call LE_Emis_MEGAN_Done(  emis_set(i_emis)%megan, status  )
          IF_NOTOK_RETURN(status=1)

#ifdef with_pollen
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'silam-pollen' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with silam pollen
          call LE_Emis_Pollen_Done( emis_set(i_emis)%pollen, status )
          IF_NOTOK_RETURN(status=1)
#endif
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'point' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with point emissions
          call emis_set(i_emis)%point%Done( status )
          IF_NOTOK_RETURN(status=1)
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'gfas' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with point emissions
          call emis_set(i_emis)%gfas%Done( status )
          IF_NOTOK_RETURN(status=1)
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'sea_salt' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! enabled ?
          if ( n_seasalt > 0 ) then
            call LE_Emis_SeaSalt_Done( status  )
            IF_NOTOK_RETURN(status=1)
          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-wind' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with wind blown dust emissions:
          call LE_Emis_Dust_Wind_Done( status  )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-resuspension' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! nothing to be done ...

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-agriculture' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! nothing to be done ...

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'hm_natural' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          ! Nothing to be done
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'basecation' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! base cat ions included ?
          if ( n_basecation > 0 ) then

            ! done with base cat ion emissions:
            call LE_Emis_BaseCatIon_Done( status  )
            IF_NOTOK_RETURN(status=1)

          endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'ops' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          ! done with OPS-emission data
          call LE_Emis_OPS_Done( emis_set(i_emis)%emo, status )
          IF_NOTOK_RETURN(status=1)
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! error ...
          write (gol,'("found unsupported emis data type : ",a)') trim(emis_set(i_emis)%filetype); call goErr
          TRACEBACK; status=1; return

      end select

    end do  ! emisison sets

    ! clear:
    deallocate( emis_a_bio, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine LE_Emis_Done


  ! ***


  subroutine LE_Emis_Setup( t1, t2, status )

    use GO                       , only : TDate, wrtgol, DayNumber
    use GO                       , only : GO_Timer_Start, GO_Timer_End
    use LE_Logging               , only : ident2
    use dims                     , only : runF, outF
    use dims                     , only : emis_a
    use dims                     , only : nspec
    use Indices                  , only : n_seasalt, n_dust, n_basecation, n_hm
    use Indices                  , only : i_so4a_f, i_so4a_c, i_pom_c, i_pom_f 
#ifdef with_m7    
    use Indices                  , only : n_m7
#endif
    use Indices                  , only : n_pom, ispecs_pom
    use Indices                  , only : n_vbs_pog, ispecs_vbs_pog
    use Indices                  , only : specmolm, specname
    use LE_Emis_EDGAR            , only : LE_Emis_EDGAR_Setup
    use LE_Emis_Bio              , only : LE_Emis_Bio_VOC
    use LE_Emis_Bio_NO           , only : LE_Emis_Bio_NO_Setup
    use LE_Emis_MEGAN            , only : LE_Emis_MEGAN_Setup
#ifdef with_pollen    
    use LE_Emis_Pollen           , only : LE_Emis_Pollen_Setup
#endif
    use LE_Emis_SeaSalt          , only : LE_Emis_SeaSalt_Add
    use LE_Emis_Dust_Wind        , only : LE_Emis_Dust_Wind_Add
    use LE_Emis_Dust_Resuspension, only : roadresp_dust_em
    use LE_Emis_Dust_Agriculture , only : agriacti_dust_em
    use LE_Emis_HM_Natural       , only : get_natural_hm
    use LE_Emis_BaseCatIon       , only : LE_Emis_BaseCatIon_Setup
    use LE_Emis_OPS              , only : LE_Emis_OPS_Setup
#ifdef with_vbs
    use LE_VBS                   , only : vbs_emisfac
#endif
#ifdef with_labeling
    use SA_Labeling              , only : SA_Emis_Reset
#endif

    ! --- in/out ------------------------------

    type(TDate), intent(in)   ::  t1, t2
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Setup'

    ! --- local -------------------------------

    integer                ::  i_emis
    integer                ::  iday

#ifdef with_vbs
    ! vbs stuff:
    integer                ::  ipom, ispec
    integer                ::  ivbs, itr
#endif

    ! --- begin -------------------------------

    ! info ...
    call wrtgol('LE:   setup emissions for ',t1,' - ',t2); call goPr

#ifdef with_m7
    ! check for M7 ...
    if ( n_m7 > 0 ) then
      ! warning ...
      if ( (i_so4a_f > 0) .or. (i_so4a_c > 0) ) then
        write (gol,'("Original LE/M7 code skipped SO4a from SO2 emissions ;")'); call goErr
        write (gol,'("ensure that this is implemented in this version through the rcfile settings")'); call goErr
        write (gol,'("Better is to not have both M7 and SO4a_[fc] ...")'); call goErr
        TRACEBACK; status=1; return
        !write (gol,'("WARNING - continue while testing ...")'); call goPr
      end if
    end if  ! m7 ?
#endif

    ! init collection arrays to zero ; in loop over missions,
    ! contributions are either added to "emis_a" or "emis_a_bio" ;
    ! after the loop, "emis_a_bio" is then added to "emis_a" :
    emis_a     = 0.0
    emis_a_bio = 0.0

#ifdef with_labeling
    ! Reset emissions in Labeling to zero, otherwise emissions are added up
    call SA_Emis_Reset(status)
    IF_NOTOK_RETURN(status=1)
#endif

    ! loop over emissions:
    do i_emis = 1, n_emis

      ! start timing:
      call GO_Timer_Start( emis_set(i_emis)%itim, status )
      IF_NOTOK_RETURN(status=1)

      ! setup data:
      select case ( trim(emis_set(i_emis)%filetype) )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'data' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! done with data 
          call emis_set(i_emis)%emd%Setup( t1, t2, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'tno' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add the actual emissions:
          call emis_set(i_emis)%emt%Setup( emis_a, t1, t2, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'edgar' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add the actual emissions:
          call LE_Emis_EDGAR_Setup( emis_set(i_emis)%edgar, emis_a, t1, t2, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'bio' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! day number in year:
          iday = DayNumber( t1 )
          ! add biogenic emissions to "emis_a_bio":
          call LE_Emis_Bio_VOC( emis_a_bio, iday, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'bio-no' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add biogenic emissions to "emis_a_bio":
          call LE_Emis_Bio_NO_Setup( emis_a_bio, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'megan' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add emissions to "emis_a_bio":
          call LE_Emis_MEGAN_Setup( emis_set(i_emis)%megan, &
                                      emis_a_bio, t1, t2, status )
          IF_NOTOK_RETURN(status=1)

#ifdef with_pollen
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'silam-pollen' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

          ! info ...
          write (gol,'("         <add silam pollen emissions>")'); call goPr
          ! fill the emission array:
          call LE_Emis_Pollen_Setup( emis_set(i_emis)%pollen, emis_a, t1, t2, status )
          IF_NOTOK_RETURN(status=1)
#endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'point' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

          ! info ...
          write (gol,'("         <add point emissions>")'); call goPr
          ! fill the emission array:
          call emis_set(i_emis)%point%Setup( emis_a, t1, t2, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'gfas' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

          ! info ...
          write (gol,'("         add GFAS emissions")'); call goPr
          ! add contribution to the emission array:
          call emis_set(i_emis)%gfas%Setup( emis_a, t1, t2, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'sea_salt' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add sea-salt emissions if necessary:
          if ( n_seasalt > 0 ) then
            call LE_Emis_SeaSalt_Add( emis_a, status )
            IF_NOTOK_RETURN(status=1)
          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-wind' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add dust emissions if necessary:
          if ( (n_dust > 0) .or. (n_basecation > 0) ) then

            call LE_Emis_Dust_Wind_Add( emis_a, t1, status )
            IF_NOTOK_RETURN(status=1)

          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-resuspension' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! local time (necessary in roadresp_dust, agriacti_dust) should be defined in the utils module
          ! by moving subroutine local_time3 in the respective module(s)
          ! ACTION: set correct local time instead of current settings in the modules

          ! add dust emissions if necessary:
          if ( (n_dust > 0) .or. (n_basecation > 0) ) then

            call roadresp_dust_em( emis_a, t1, status )
            IF_NOTOK_RETURN(status=1)

          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'dust-agriculture' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! local time (necessary in roadresp_dust, agriacti_dust) should be defined in the utils module
          ! by moving subroutine local_time3 in the respective module(s)
          ! ACTION: set correct local time instead of current settings in the modules

          ! add dust emissions if necessary:
          if ( (n_dust > 0) .or. (n_basecation > 0) ) then

            call agriacti_dust_em( emis_a, t1, status )
            IF_NOTOK_RETURN(status=1)

          end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'hm_natural' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! heavy metals included ?
          if ( n_hm > 0 ) then

            ! setup dust emissions:
            call get_natural_hm( emis_a, status  )
            IF_NOTOK_RETURN(status=1)

          endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'basecation' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! base cat ions included ?
          if ( n_basecation > 0 ) then

            ! setup cat ion emissions:
            call LE_Emis_BaseCatIon_Setup( emis_a, t1%year, t1%month, t1%day, t1%hour, status  )
            IF_NOTOK_RETURN(status=1)

          endif

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( 'ops' )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! add the actual emissions:
          call LE_Emis_OPS_Setup( emis_set(i_emis)%emo, emis_a, t1, t2, status )
          IF_NOTOK_RETURN(status=1)

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! error ...
          write (gol,'("found unsupported emis data type : ",a)') trim(emis_set(i_emis)%filetype); call goErr
          TRACEBACK; status=1; return

      end select
      
      ! end timing:
      call GO_Timer_End( emis_set(i_emis)%itim, status )
      IF_NOTOK_RETURN(status=1)

    end do   ! emission sets

    ! add "bio" contributions to total:
    emis_a = emis_a + emis_a_bio

#ifdef with_vbs
    !
    ! If vbs is on and pom tracers exist, then primary organic matter (POM) emissions are distributed
    ! amongst the vbs classes according to Shrivastata et al, JGR 2008. Thus, the amount of emissions
    ! into vbs classes is assumed to be 2.5 times the total of POM emissions.
    ! Note that the tracers pom_c and pom_f should actually be named pom_c and pom_f because
    ! it already contains the non-carbon atoms.
    !
    ! redestribute prim.org.carb emissions over vbs tracers ?
    if ( with_pom_to_vbs ) then
      ! info ...
!      write (gol,'("Distribution of primary organic carbon to VBS tracers")'); call goErr
!      write (gol,'("should be implemented using composition tables.")'); call goErr
!      write (gol,'("Also think about whether anthro is different from bio emis.")'); call goErr
!      TRACEBACK; status=1; return
      !write (gol,'("WARNING - continue while testing ...")'); call goPr
      ! loop over prim.org.carbon tracers:
      do ipom = 1, n_pom
        ! source tracer index:
        ispec = ispecs_pom(ipom)
        !
        if ( ispec == i_pom_c ) then
          ! split only fine part of POM
          cycle
        endif
        ! enabled ?
        if ( ispec > 0 ) then
          ! loop over vbs primary organic gasses:
          do ivbs = 1, n_vbs_pog
            ! target tracer index:
            itr = ispecs_vbs_pog(ivbs)
            !      mol/min         mol/min                 ug/min              mass ratio  /  (kg/mol)    kg/ug
            emis_a(:,:,:,itr) = emis_a(:,:,:,itr) + emis_a(:,:,:,ispec) * vbs_emisfac(ivbs)/specmolm(itr)*1.0e-9
          end do  ! vbs tracers
        end if  ! source tracer enabled
        ! reset original emissions to zero, since now distributed over vbs tracers:
        emis_a(:,:,:,ispec) = 0.0
      end do  ! source tracers
    end if   ! pom to vbs
#endif

    ! ok
    status = 0

  end subroutine LE_Emis_Setup


  ! ***


  !
  ! Add emissions emitted during 'dt_min' minutes to concentration array.
  !

  subroutine LE_Emis_Add( c, dt_min, status )

    use dims   , only : nx, ny, nz, nspec
    use dims   , only : runF
    use dims   , only : emis_a
    use binas  , only : Rgas
    use indices, only : specname, specunit
#ifdef with_labeling
    use SA_Labeling, only : SA_Emis_Delta_c
#endif

    ! point to meteo data:
    use LE_Data      , only : LE_Data_GetPointer
    use LE_Meteo_Data, only : volume  ! km3 (zucht)

    ! --- in/out ------------------------------

    real, intent(inout)       ::  c(nx,ny,nz,nspec)
    real, intent(in)          ::  dt_min   ! time step in minutes
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Add'

    ! --- local -------------------------------

    integer   ::  ispec
    integer   ::  status_par

    real, allocatable      ::  ugm3_fac(:,:,:)
    real, allocatable      ::  ncm3_fac(:,:,:)
    real, allocatable      ::  nm3_fac(:,:,:)
    real, allocatable      ::  ppb_fac(:,:,:)
    real, allocatable      ::  unity_fac(:,:,:)

    ! meteo data:
    !real, pointer               ::  volume(:,:,:)   ! (lon,lat,nz)
    real, pointer               ::  temp  (:,:,:)   ! (lon,lat,nz)
    real, pointer               ::  pres  (:,:,:)   ! (lon,lat,nz)    

    ! --- begin ----------------------------------
    
    ! point to meteo data:
    call LE_Data_GetPointer( 't'  , temp  , status, check_units='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'p'  , pres  , status, check_units='Pa' )
    IF_NOTOK_RETURN(status=1)
    !call LE_Data_GetPointer( 'vol', volume, status, check_units='m3' )
    !IF_NOTOK_RETURN(status=1)

    !! info ...
    !write (gol,'("       <add emissions> over ",f6.1," minutes")') dt_min; call goPr

    ! storage:
    allocate( ppb_fac(nx,ny,nz) )
    allocate( ugm3_fac(nx,ny,nz) )
    allocate( ncm3_fac(nx,ny,nz) )
    allocate( nm3_fac(nx,ny,nz) )
    allocate( unity_fac(nx,ny,nz) )
    
    ! conversion factors:
    !            1     K                        Pa                    km3
    ppb_fac = (Rgas * temp(1:nx,1:ny,1:nz)) / (pres(1:nx,1:ny,1:nz)*volume(1:nx,1:ny,1:nz))  ! ppb/mol
    !           km3/m3     /  km3
    ugm3_fac =  1.0e-9     / volume(1:nx,1:ny,1:nz)   ! (ug/m3)/ug
    !           km3/cm3    /  km3
    ncm3_fac =  1.0e-15    / volume(1:nx,1:ny,1:nz)   ! (#/cm3)/#
    !           km3/m3    /  km3
    nm3_fac  =  1.0e-9     / volume(1:nx,1:ny,1:nz)   ! (grns/cm3)/grns
    ! 
    unity_fac = 1.0
    
    ! init parallel status:
    status_par = 0

    !$OMP parallel &
#ifndef __GFORTRAN__
    !$OMP default (none ) &
    !$OMP shared( specname, specunit ) &
#endif
    !$OMP shared( ppb_fac, ugm3_fac, ncm3_fac, nm3_fac, unity_fac ) &
    !$OMP shared( dt_min ) &
    !$OMP shared( c ) &
    !$OMP shared( emis_a ) &
    !$OMP private( ispec ) &
    !$OMP private( status )  &
    !$OMP reduction( + : status_par )
    !$OMP do SCHEDULE( DYNAMIC )
    ! loop over tracers:
    do ispec = 1, nspec

      ! check units ...
      select case ( specunit(ispec) )

        !~ volume mixing ratio's:
        case ( 'ppb' )
          !
          ! unit emis_a for gases   : mol/min
          !                    ppb           min      ppb/mol      mol/min
          c(:,:,:,ispec) = c(:,:,:,ispec) + dt_min * ppb_fac * emis_a(:,:,:,ispec)
#ifdef with_labeling
          ! update label fractions:
          call SA_Emis_Delta_c( ispec, dt_min * ppb_fac, status )
          status_par = status_par + status
#endif

        !~ mass concentrations:
        case ( 'ug/m3' )
          ! unit emis_a for aerosols   : ug/min
          !                    ug/m3           min  (ug/m3)/(ug)      ug/min
          c(:,:,:,ispec) = c(:,:,:,ispec) +  dt_min *  ugm3_fac    * emis_a(:,:,:,ispec)
#ifdef with_labeling
          ! update label fractions:
          call SA_Emis_Delta_c( ispec, dt_min * ugm3_fac, status )
          status_par = status_par + status
#endif

        !~ number concentrations:
        case ( '1/cm3', 'mlc/cm3' )
          ! from PN database, emis_a is in #/min (numbers) or mlc/min (concentrations) :
          !    #/cm3                          min     (#/cm3)/#          #/min
          !  mlc/cm3                          min     (#/cm3)/#          mlc/min
          c(:,:,:,ispec) = c(:,:,:,ispec) + dt_min  *  ncm3_fac  *  emis_a(:,:,:,ispec)

        
#ifdef with_pollen
        ! ~ grains concentrations
        case ( 'grns/m3' )
          ! unit emis_a for pollen grains : grns/min
          !
          ! grns/m3      =    grns/m3     +   min  *  (grns/m3)/grns         grns/min
          c(:,:,:,ispec) = c(:,:,:,ispec) + dt_min *  nm3_fac    * emis_a(:,:,:,ispec)
          
#ifdef with_labeling
          ! update label fractions
          call SA_Emis_Delta_c( ispec, dt_min * unity_fac, status )
          status_par = status_par + status
#endif
                  
#endif
        !~ unknown ...
        case default
          ! no write statements allowed within OpenMP loop ...
          !write (gol,'("unsupported unit `",a,"` for tracer ",i4," (",a,")")'), &
          !                 trim(specunit(ispec)), ispec, trim(specname(ispec)); call goErr
          status_par = status_par + 1
      end select

    end do  ! tracers
    !$OMP end do
    !$OMP end parallel

    ! check status ...
    if ( status_par /= 0 ) then
      write (gol,'("error status returned from parallel loop over tracers : ",i6)') status_par; call goErr
      TRACEBACK; status=1; return
    end if

    ! clear:
    deallocate( ppb_fac )
    deallocate( ugm3_fac )
    deallocate( ncm3_fac )
    deallocate( nm3_fac )
    deallocate( unity_fac )

    ! ok
    status = 0

  end subroutine LE_Emis_Add


end module LE_Emis

