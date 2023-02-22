!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_UDUNITS_NOTOK_RETURN(action) if (status/=UDUNITS_NOERR) then; gol=trim(UDUnits_StrError(status)); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Data_Common

  use GO     , only : gol, goPr, goErr
  use C3PO   , only : T_Levs_Hyb
  use C3PO   , only : T_Grid_Convertors
  use UDUnits, only : UDUNITS_NOERR, UDUnits_StrError

  implicit none


  ! --- in/out -----------------------------------
  
  private
  
  public  ::  LE_Data_Common_Init, LE_Data_Common_Done
  public  ::  LE_Data_CompareUnits
  public  ::  LE_Data_ConversionFactor
  public  ::  LE_Data_ClimatMask
  
  public  ::  levels_type
  public  ::  nlev, nlev_top
  public  ::  mixlayer_surf_top
  public  ::  mixlayer_mix_topmin
  public  ::  mixlayer_dmin, mixlayer_sdofac
  public  ::  mixlayer_top
  public  ::  mixlayer_daloft
  public  ::  hyb, hyb_top
  public  ::  metlevel_select, metlevel_select_top

  public  ::  Grid_Convertors


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'LE_Data_Common'


  ! --- var --------------------------------------
  
  ! current level definition:
  character(len=64)             ::  levels_type

  ! number of levels:
  integer                       ::  nlev
  integer                       ::  nlev_top

  ! settings for "mixlayer" type:
  real                          ::  mixlayer_surf_top
  real                          ::  mixlayer_mix_topmin
  real                          ::  mixlayer_dmin
  real, allocatable             ::  mixlayer_sdofac(:)  ! (1:nlev)
  real, allocatable             ::  mixlayer_top(:)     ! (1:nlev)
  real                          ::  mixlayer_daloft

  ! settings for "hyblevel" type:
  type(T_Levs_Hyb)              ::  hyb        ! (1:nlev)
  type(T_Levs_Hyb)              ::  hyb_top    ! (1:nlev_top)

  ! settings for "metlevel" type;
  ! mapping from original meteo levels to model levels:
  !    (/0,1,3,6,9,12/)   # 5 layers out original 12, coarsening 1,2,3,3,3
  integer, allocatable, target  ::  metlevel_select(:)      ! (0:nlev)
  integer, allocatable, target  ::  metlevel_select_top(:)  ! (0:nlev_ntop)

  ! convertors:
  type(T_Grid_Convertors)       ::  Grid_Convertors



contains


  ! ====================================================================
  ! ===
  ! === module
  ! ===
  ! ====================================================================
  

  subroutine LE_Data_Common_Init( rcF, status )
  
    use GO     , only : TrcFile, ReadRc
    use GO     , only : goSplitString
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)               ::  rcF
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Data_Common_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------
    
    ! * level definitions
    
    ! number of model layers:
    call ReadRc( rcF, 'levels.nz', nlev, status )
    IF_NOTOK_RETURN(status=1)
    
    ! all layers for output of entire atmosphere:
    call ReadRc( rcF, 'levels.nz_top', nlev_top, status )
    IF_NOTOK_RETURN(status=1)
    
    ! type:
    call ReadRc( rcF, 'levels.type', levels_type, status )
    IF_NOTOK_RETURN(status=1)
    ! switch:
    select case ( trim(levels_type) )
      !~
      case ( 'mixlayer' )
        ! read settings:
        call LE_Data_Common_MixLayer_Init( rcF, status )
        IF_NOTOK_RETURN(status=1)
        ! check ...
        if ( nlev_top /= nlev+1 ) then
          write (gol,'("level type `mixlayer` supports one aloft layer only yet")'); call goErr
          TRACEBACK; status=1; return
        end if
      !~
      case ( 'hyblevel' )
        ! fill layer definition from data file:
        call LE_Data_Common_HybLevel_Init( rcF, nlev, hyb, status )
        IF_NOTOK_RETURN(status=1)
        ! layer definition including top layers:
        call LE_Data_Common_HybLevel_Init( rcF, nlev_top, hyb_top, status )
        IF_NOTOK_RETURN(status=1)
      !~
      case ( 'metlevel' )
        ! read settings for both "nlev" and "nlev_top":
        call LE_Data_Common_MetLevel_Init( rcF, status )
        IF_NOTOK_RETURN(status=1)

      case default
        write (gol,'("unsupported level type `",a,"`")') trim(levels_type); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! * regridding
    
    ! init conversions:
    call Grid_Convertors%Init( status )
    IF_NOTOK_RETURN(status=1)
    ! ok
    status = 0
    
  end subroutine LE_Data_Common_Init
  
  
  ! ***
  
  
  subroutine LE_Data_Common_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Common_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! done with conversions:
    call Grid_Convertors%Done( status )
    IF_NOTOK_RETURN(status=1)
    
    ! switch:
    select case ( trim(levels_type) )
      !~
      case ( 'mixlayer' )
        ! clear:
        call LE_Data_Common_MixLayer_Done( status )
        IF_NOTOK_RETURN(status=1)
      !~
      case ( 'hyblevel' )
        ! clear:
        call LE_Data_Common_HybLevel_Done( status )
        IF_NOTOK_RETURN(status=1)

      case ( 'metlevel' )
        ! clear:
        call LE_Data_Common_MetLevel_Done( status )
        IF_NOTOK_RETURN(status=1)
      !
      case default
        write (gol,'("unsupported level type `",a,"`")') trim(levels_type); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0
    
  end subroutine LE_Data_Common_Done
  

  ! ***


  subroutine LE_Data_Common_MixLayer_Init( rcF, status )
  
    use GO     , only : TrcFile, ReadRc
    use GO     , only : goSplitString
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)               ::  rcF
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Data_Common_MixLayer_Init'
    
    ! --- local ------------------------------------------
    
    integer               ::  ntop
    real                  ::  tops(20)
    integer               ::  nfac
    real                  ::  facs(20)
    character(len=1024)   ::  line

    ! --- begin ------------------------------------------
    
    ! thickness of surface layer:
    call ReadRc( rcF, 'mixlayer.surf_top', mixlayer_surf_top, status )
    IF_NOTOK_RETURN(status=1)
    ! minium top of mix layer:
    call ReadRc( rcF, 'mixlayer.mix_topmin', mixlayer_mix_topmin, status )
    IF_NOTOK_RETURN(status=1)
    ! line with values:
    call ReadRc( rcF, 'mixlayer.top', line, status )
    IF_NOTOK_RETURN(status=1)
    ! split:
    call goSplitString( line, ntop, tops, status )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    if ( 3 + ntop /= nlev ) then
      write (gol,'("number of top levels (",i0,") + 3 should be nz (",i0,")")') ntop, nlev; call goErr
      TRACEBACK; status=0; return
    end if
    ! storage for target layer tops:
    allocate( mixlayer_top(nlev), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! fill:
    mixlayer_top(1) = mixlayer_surf_top    ! dummy, not used ...
    mixlayer_top(2) = mixlayer_mix_topmin  ! dummy, not used ...
    mixlayer_top(3) = ( mixlayer_mix_topmin + tops(1) )/2.0  ! dummy, not used ...
    mixlayer_top(4:nlev) = tops(1:ntop)

    ! minimum thickness of residual layer:
    call ReadRc( rcF, 'mixlayer.dmin', mixlayer_dmin, status )
    IF_NOTOK_RETURN(status=1)
    ! line with factors for orography std.dev. added to dmin:
    call ReadRc( rcF, 'mixlayer.sdofac', line, status )
    IF_NOTOK_RETURN(status=1)
    ! split:
    call goSplitString( line, nfac, facs, status )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    if ( nfac /= ntop ) then
      write (gol,'("number of sdo factors (",i0,") should be the same as number of level tops (",i0,")")') nfac, ntop; call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( mixlayer_sdofac(nlev), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! copy:
    mixlayer_sdofac = 0.0
    mixlayer_sdofac(4:nlev) = facs(1:nfac)

    ! thickness of top layer to hold boundary conditions:
    call ReadRc( rcF, 'mixlayer.daloft', mixlayer_daloft, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Data_Common_MixLayer_Init
  
  
  ! ***
  
  
  subroutine LE_Data_Common_MixLayer_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Common_MixLayer_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( mixlayer_top, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( mixlayer_sdofac, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Data_Common_MixLayer_Done
  

  ! ***


  !
  ! Read hybride coefficents.
  ! Example of file content:
  !
  !   #hlev, hyai [Pa], hybi [1]
  !    20,  6759.727  , 0.0000000
  !    19,  8564.624  , 0.0000546672
  !    18, 10558.88   , 0.0005483841
  !    17, 12713.9    , 0.002764719
  !    16, 14922.69   , 0.00903499
  !    15, 16990.62   , 0.02218864
  !    14, 18716.97   , 0.04514613
  !    13, 19919.8    , 0.0807773
  !    12, 20434.16   , 0.1319348
  !    11, 20107.03   , 0.2015201
  !    10, 18798.82   , 0.2919934
  !     9, 16544.59   , 0.3992048
  !     8, 13653.22   , 0.5132797
  !     7, 10471.31   , 0.6265589
  !     6,  7335.165  , 0.7322238
  !     5,  4550.216  , 0.8241854
  !     4,  2356.203  , 0.8977672
  !     3,   895.1935 , 0.9502738
  !     2,   162.0434 , 0.9822381
  !     1,     0.00316, 0.9976301
  !     0,     0.00000, 1.0000000
  !
  ! Each line contains a half-level definition.
  ! If the model is defined on a smaller number of levels,
  ! the top levels are skipped.
  !

  subroutine LE_Data_Common_HybLevel_Init( rcF, n, hb, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : TTextFile, Init, Done, ReadLine
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)               ::  rcF
    integer, intent(in)                     ::  n
    type(T_Levs_Hyb), intent(out)           ::  hb   ! layers (1:n)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Data_Common_HybLevel_Init'
    
    ! --- local ------------------------------------------

    character(len=1024)   ::  filename
    type(TTextFile)       ::  file
    character(len=1024)   ::  line
    integer               ::  hlev
    real, allocatable     ::  hyai(:), hybi(:)
    real, allocatable     ::  hyam(:), hybm(:)
    logical, allocatable  ::  filled(:)

    ! --- begin ------------------------------------------
    
    ! storage:
    allocate( hyai(0:n), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( hybi(0:n), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( hyam(n), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( hybm(n), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( filled(0:n), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! data file:
    call ReadRc( rcF, 'hyblevel.coefficients', filename, status )
    IF_NOTOK_RETURN(status=1)
    
    ! open existing file:
    call Init( file, filename, status, status='old', comment='#' )
    IF_NOTOK_RETURN(status=1)

    ! init flags:
    filled = .false.
    ! loop over lines:
    do
      ! read next line, status<0 if end-of-file:
      call ReadLine( file, line, status  )
      if ( status < 0 ) exit
      IF_NOTOK_RETURN(status=1)
      
      ! half level number:
      call goReadFromLine( line, hlev, status )
      IF_NOTOK_RETURN(status=1)
      ! skip if above current top:
      if ( hlev > n ) cycle
      ! check if it makes sense ...
      if ( hlev < 0 ) then
        write (gol,'("half level in file should be 0,1,.. while found ",i0)') hlev; call goErr
        write (gol,'("input file : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! already filled?
      if ( filled(hlev) ) then
        write (gol,'("half level ",i0," already read")') hlev; call goErr
        write (gol,'("input file : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! read coeff:
      call goReadFromLine( line, hyai(hlev), status )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( line, hybi(hlev), status )
      IF_NOTOK_RETURN(status=1)
      ! set flag:
      filled(hlev) = .true.

    end do ! lines
    
    ! close:
    call Done( file, status )
    IF_NOTOK_RETURN(status=1)
    
    ! not all filled?
    if ( .not. all(filled) ) then
      write (gol,'("not all ",i0," half levels (",i0," incl top) filled:")') n, n+1; call goErr
      write (gol,'("  hlev                   a                   b")'); call goErr
      do hlev = n, 0, -1
        if ( filled(hlev) ) then
          write (gol,'(i6,2f20.1)') hlev, hyai(hlev), hybi(hlev); call goErr
        else
          write (gol,'(i6,2a20)') hlev, '-', '-'; call goErr
        end if
      end do
      write (gol,'("input file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! mid levels:
    hyam = 0.5 * ( hyai(0:n-1) + hyai(1:n) )
    hybm = 0.5 * ( hybi(0:n-1) + hybi(1:n) )

    ! define levels:
    call hb%Init( hyam, hybm, hyai, hybi, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( hyai, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( hybi, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( hyam, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( hybm, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( filled, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Data_Common_HybLevel_Init
  
  
  ! ***
  
  
  subroutine LE_Data_Common_HybLevel_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Common_HybLevel_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    call hyb%Done( status )
    IF_NOTOK_RETURN(status=1)
    call hyb_top%Done( status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Data_Common_HybLevel_Done


  ! ***


  subroutine LE_Data_Common_MetLevel_Init( rcF, status )
  
    use GO     , only : TrcFile, ReadRc
    use GO     , only : goSplitString
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)               ::  rcF
    integer, intent(out)                    ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Data_Common_MetLevel_Init'
    
    ! --- local ------------------------------------------
    
    character(len=1024)   ::  line
    integer               ::  ncomb
    integer               ::  combine(100)
    integer               ::  k

    ! --- begin ------------------------------------------
    
    ! storage:
    allocate( metlevel_select(0:nlev), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( metlevel_select_top(0:nlev_top), stat=status )
    IF_NOTOK_RETURN(status=1)

    ! line with number of meteo layers to combine:
    call ReadRc( rcF, 'metlevel.combine', line, status )
    IF_NOTOK_RETURN(status=1)
    ! defined?
    if ( len_trim(line) > 0 ) then
      ! split:
      call goSplitString( line, ncomb, combine, status )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( ncomb /= nlev ) then
        write (gol,'("number of level combinationss (",i0,") should be the same as number of levels (",i0,")")') ncomb, nlev; call goErr
        TRACEBACK; status=1; return
      end if
      ! fill ranges:
      metlevel_select(0) = 0
      do k = 1, nlev
        metlevel_select(k) = metlevel_select(k-1) + combine(k)
      end do
    else
      ! no combinations, just copy:
      do k = 0, nlev
        metlevel_select(k) = k
      end do
    end if  ! line

    ! line with number of meteo layers to combine:
    call ReadRc( rcF, 'metlevel.combine_top', line, status )
    IF_NOTOK_RETURN(status=1)
    ! defined?
    if ( len_trim(line) > 0 ) then
      ! split:
      call goSplitString( line, ncomb, combine, status )
      IF_NOTOK_RETURN(status=1)
      ! check ...
      if ( ncomb /= nlev_top ) then
        write (gol,'("number of level combinationss (",i0,") should be the same as number of levels incl. top (",i0,")")') ncomb, nlev_top; call goErr
        TRACEBACK; status=1; return
      end if
      ! fill ranges:
      metlevel_select_top(0) = 0
      do k = 1, nlev_top
        metlevel_select_top(k) = metlevel_select_top(k-1) + combine(k)
      end do
    else
      ! no combinations, just copy:
      do k = 0, nlev_top
        metlevel_select_top(k) = k
      end do
    end if  ! line
    
    ! testing ...
    write (gol,'("meteo half level selection:")'); call goPr
    do k = 0, nlev_top
      if ( k <= nlev ) then
        write (gol,'("  model layer ",i3," from ",i3)') k, metlevel_select(k); call goPr
      else
        write (gol,'("  top layer   ",i3," from ",i3)') k, metlevel_select_top(k); call goPr
      end if
    end do

    ! ok
    status = 0
    
  end subroutine LE_Data_Common_MetLevel_Init
  
  
  ! ***
  
  
  subroutine LE_Data_Common_MetLevel_Done( status )
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_Common_MetLevel_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( metlevel_select, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( metlevel_select_top, stat=status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Data_Common_MetLevel_Done

  ! ***
  
  !
  ! Return status:
  !  -1  : units are not equivalent
  !   0  : units are equivalent
  !  >0  : error
  !
  
  subroutine LE_Data_CompareUnits( units1, units2, status, verbose )
  
    use UDUnits, only : UDUnits_Standard

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)            ::  units1, units2
    integer, intent(out)                    ::  status
    logical, intent(in), optional           ::  verbose
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_CompareUnits'

    ! --- local ----------------------------------

    logical               ::  verb
    character(len=64)     ::  units1_standard
    character(len=64)     ::  units2_standard
    
    ! --- begin ----------------------------------
    
    ! shout?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! only further check if values are different ...
    if ( trim(units1) /= trim(units2) ) then

      ! adhoc comparisons to avoid use of udunits;
      ! the following matches are needed to run with open source test data:
      if ( ( (trim(units1) == '(0 - 1)'                    ) .and. (trim(units2) == '1'            ) ) .or. &
           ( (trim(units1) == 'degrees'                    ) .and. (trim(units2) == 'degrees_east' ) ) .or. &
           ( (trim(units1) == 'degrees'                    ) .and. (trim(units2) == 'degrees_north') ) .or. &
           ( (trim(units1) == 'kg kg**-1'                  ) .and. (trim(units2) == 'kg/kg'        ) ) .or. &
           ( (trim(units1) == 'm s**-1'                    ) .and. (trim(units2) == 'm/s'          ) ) .or. &
           ( (trim(units1) == 'm/s'                        ) .and. (trim(units2) == 'm s-1'        ) ) .or. &
           ( (trim(units1) == 'm**3 m**-3'                 ) .and. (trim(units2) == 'm3/m3'        ) ) .or. &
           ( (trim(units1) == 'J m**-2 s**-1'              ) .and. (trim(units2) == 'J/m2/s'       ) ) .or. &
           ( (trim(units1) == 'J/m2/s'                     ) .and. (trim(units2) == 'W/m2'         ) ) .or. &
           ( (trim(units1) == 'm of water equivalent'      ) .and. (trim(units2) == 'm'            ) ) .or. &
           ( (trim(units1) == 'm of water equivalent s**-1') .and. (trim(units2) == 'm/s'          ) )      ) then

        ! equivalent units ...
        
      ! adhoc comparisons to avoid use of udunits;
      ! the following comparisons should provide a warning status
      ! that should be handled by the calling routine:
      else if ( ( (trim(units1) == 'kg/kg') .and. (trim(units2) == 'ppb'  ) ) .or. &
                ( (trim(units1) == 'kg/kg') .and. (trim(units2) == 'ug/m3') )      ) then
      
        ! warning status:
        status = -1; return

      else

        ! convert to standard units:
        call UDUnits_Standard( units1, units1_standard, status )
        if ( status /= UDUNITS_NOERR ) then
          gol=trim(UDUnits_StrError(status)); call goErr
          write (gol,'("could not obtain standard units for `",a,"` for conversion to `",a,"`")') &
                          trim(units1), trim(units2); call goErr
          write (gol,'("no UDUnits library enabled?")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! convert to standard units:
        call UDUnits_Standard( units2, units2_standard, status )
        if ( status /= UDUNITS_NOERR ) then
          gol=trim(UDUnits_StrError(status)); call goErr
          write (gol,'("could not obtain standard units for `",a,"` for conversion from `",a,"`")') &
                          trim(units2), trim(units1); call goErr
          write (gol,'("no UDUnits library enabled?")'); call goErr
          TRACEBACK; status=1; return
        end if

        ! the above calls might return empty units, 
        ! but these are not allowed in the cf standard;
        ! reset them to unity in that case:
        if ( len_trim(units1_standard) == 0 ) units1_standard = '1'
        if ( len_trim(units2_standard) == 0 ) units2_standard = '1'

        !! testing ...
        !write (gol,'("compare units `",a,"` (`",a,"`) with `",a,"` (`",a"`)")') &
        !                trim(units1), trim(units1_standard), &
        !                trim(units2), trim(units2_standard); call goPr

        ! compare:
        if ( trim(units1_standard) /= trim(units2_standard) ) then
          ! info?
          if ( verb ) then
            write (gol,'("units `",a,"` and `",a,"` are not equivalent")') trim(units1), trim(units2); call goPr
            write (gol,'("(standard units `",a,"` and `",a,"`)")') trim(units1_standard), trim(units2_standard); call goPr
          end if
          ! warning status:
          status = -1 ; return
        end if

        ! safety check: m3/m3 is not equivalent to kg/kg ...
        if ( (trim(units1_standard)=='1') .and. (trim(units2_standard)=='1') &
                    .and. (units1(1:1) /= units2(1:1)) ) then
          ! info?
          if ( verb ) then
            write (gol,'("units `",a,"` and `",a,"` might not be equivalent")') trim(units1), trim(units2); call goPr
            write (gol,'("(ambigious standard units `",a,"`)")') trim(units1_standard); call goPr
          end if
          ! warning status:
          status = -1 ; return
        end if

        !! info?
        !if ( verb ) then
        !  write (gol,'("units `",a,"` and `",a,"` are equivalent")') trim(units1), trim(units2); call goPr
        !end if
      
      end if  ! adhoc match

    end if  ! unit descriptions are different
    
    ! ok
    status = 0
    
  end subroutine LE_Data_CompareUnits
  

  ! ***
  
  !
  ! Return status:
  !   0  : units are equivalent, factor computed
  !  >0  : error
  !
  
  subroutine LE_Data_ConversionFactor( units1, units2, factor, status, keys )
  
    use UDUnits, only : UDUnits_Standard
    use UDUnits, only : UDUNITS_NOERR, UDUnits_StrError
    use UDUnits, only : UDUnits_ConversionFactor
    use Binas  , only : xm_O, xm_H, xm_S, xm_N, xm_C
    use Binas  , only : xm_air
    use GO     , only : GoVarValue

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)            ::  units1, units2
    real, intent(out)                       ::  factor
    integer, intent(out)                    ::  status
    character(len=*), intent(in), optional  ::  keys
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_ConversionFactor'

    ! --- local ----------------------------------
    
    logical               ::  convert_mixing_ratio
    character(len=256)    ::  convert
    character(len=32)     ::  tracer
    real                  ::  xm_tracer
    character(len=64)     ::  units1_standard
    character(len=64)     ::  units2_standard
    
    ! --- begin ----------------------------------
    
    ! defaults:
    convert = 'None'
    tracer = 'None'
    ! optional keys provided?
    if ( present(keys) ) then
      ! extract keys:
      call GoVarValue( keys, ';', 'convert', '=', convert, status )
      IF_ERROR_RETURN(status=1)
      call GoVarValue( keys, ';', 'tracer', '=', tracer, status )
      IF_ERROR_RETURN(status=1)
    end if

    ! trap equivalent units that occure in test data;
    ! this to avoid the need for a udunits library:
    if ( (  trim(units1) == trim(units2)                                                           ) .or. &
         ( (trim(units1) == '(0 - 1)'                    ) .and. (trim(units2) == '1'            ) ) .or. &
         ( (trim(units1) == 'degrees'                    ) .and. (trim(units2) == 'degrees_east' ) ) .or. &
         ( (trim(units1) == 'degrees'                    ) .and. (trim(units2) == 'degrees_north') ) .or. &
         ( (trim(units1) == 'kg kg**-1'                  ) .and. (trim(units2) == 'kg/kg'        ) ) .or. &
         ( (trim(units1) == 'm s**-1'                    ) .and. (trim(units2) == 'm/s'          ) ) .or. &
         ( (trim(units1) == 'm**3 m**-3'                 ) .and. (trim(units2) == 'm3/m3'        ) ) .or. &
         ( (trim(units1) == 'J m**-2 s**-1'              ) .and. (trim(units2) == 'J/m2/s'       ) ) .or. &
         ( (trim(units1) == 'm of water equivalent'      ) .and. (trim(units2) == 'm'            ) ) .or. &
         ( (trim(units1) == 'm of water equivalent s**-1') .and. (trim(units2) == 'm/s'          ) )      ) then

      ! no conversion needed:
      factor = 1.0
    
    ! mixing ratio conversion
    else if ( trim(convert) == 'mole/mole to kg/kg' ) then
      
        ! check input units ...
        select case ( trim(units1) )
          case ( 'mole/mole', 'mole mole**-1', 'mole mole-1' )
            ! ok
          case default
            write (gol,'("unsupported input units `",a,"` for conversion `",a,"`")') &
                             trim(units1), trim(convert); call goErr
            TRACEBACK; status=1; return
        end select
      
        ! check output units ...
        select case ( trim(units2) )
          case ( 'kg/kg', 'kg kg**-1', 'kg kg-1' )
            ! ok
          case default
            write (gol,'("unsupported output units `",a,"` for conversion `",a,"`")') &
                             trim(units2), trim(convert); call goErr
            TRACEBACK; status=1; return
        end select
               
        ! mole mass of tracer:
        select case ( trim(tracer) )
          case ( 'FORM'  , 'form'   ) ; xm_tracer = xm_C + xm_H * 2 + xm_O                ! kg/mol
          case ( 'CH4'   , 'ch4'    ) ; xm_tracer = xm_C + xm_H * 4                       ! kg/mol
          case ( 'ETHA'  , 'etha'   ) ; xm_tracer = xm_C * 2 + xm_H * 6                   ! kg/mol
          case ( 'CH3CHO', 'ch3cho' ) ; xm_tracer = xm_C * 2 + xm_H * 4 + xm_O            ! kg/mol
          case ( 'CO'    , 'co'     ) ; xm_tracer = xm_C + xm_O                           ! kg/mol
          case ( 'HNO3'  , 'hno3'   ) ; xm_tracer = xm_H + xm_N + xm_O * 3                ! kg/mol
          case ( 'ISOP'  , 'isop'   ) ; xm_tracer = xm_C * 5 + xm_H * 8                   ! kg/mol
          case ( 'NO'    , 'no'     ) ; xm_tracer = xm_N + xm_O                           ! kg/mol
          case ( 'NO2'   , 'no2'    ) ; xm_tracer = xm_N + xm_O * 2                       ! kg/mol
          case ( 'O3'    , 'o3'     ) ; xm_tracer = xm_O * 3                              ! kg/mol
          case ( 'PAN'   , 'pan'    ) ; xm_tracer = xm_C * 2 + xm_H * 3 + xm_N + xm_O * 5 ! kg/mol
          case ( 'TOL'   , 'tol'    ) ; xm_tracer = xm_C * 7 + xm_H * 8                   ! kg/mol
          case ( 'SO2'   , 'so2'    ) ; xm_tracer = xm_S + xm_O * 2                       ! kg/mol
          case ( 'OH'    , 'oh'     ) ; xm_tracer = xm_O + xm_H                           ! kg/mol
          case default 
            write( gol,'("no mole mass defined for tracer `",a,"`")') trim(tracer); call GoErr
            TRACEBACK;status=1;return
        end select
          
        ! conversion factor:
        !        (kg tr)/(mole tr) / (kg air)/(mole air)
        factor =     xm_tracer     /      xm_air         ! [(kg tr)/(kg air)]/[(mole tr)/(mole air)]
        
        ! testing ..
        write (gol,'("conversion `",a,"` for tracer `",a,"` with factor: ",f12.6)') trim(convert), trim(tracer), factor; call goPr
      
    else

      ! convert to standard units:
      call UDUnits_Standard( units1, units1_standard, status )
      if ( status /= UDUNITS_NOERR ) then
        gol=trim(UDUnits_StrError(status)); call goErr
        write (gol,'("could not obtain standard units for `",a,"` for conversion to `",a,"`")') &
                        trim(units1), trim(units2); call goErr
        write (gol,'("no UDUnits library enabled?")'); call goErr
        TRACEBACK; status=1; return
      end if
      call UDUnits_Standard( units2, units2_standard, status )
      if ( status /= UDUNITS_NOERR ) then
        gol=trim(UDUnits_StrError(status)); call goErr
        write (gol,'("could not obtain standard units for `",a,"` for conversion from `",a,"`")') &
                        trim(units2), trim(units1); call goErr
        write (gol,'("no UDUnits library enabled?")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! the above calls might return empty units, 
      ! but these are not allowed in the cf standard;
      ! reset them to unity in that case:
      if ( len_trim(units1_standard) == 0 ) units1_standard = '1'
      if ( len_trim(units2_standard) == 0 ) units2_standard = '1'

      !! testing ...
      !write (gol,'("compare units `",a,"` (`",a,"`) with `",a,"` (`",a"`)")') &
      !                trim(units1), trim(units1_standard), &
      !                trim(units2), trim(units2_standard); call goPr

      ! safety check: m3/m3 is not equivalent to kg/kg ...
      if ( (trim(units1_standard)=='1') .and. (trim(units2_standard)=='1') &
                  .and. (units1(1:1) /= units2(1:1)) ) then

        write (gol,'("units `",a,"` and `",a,"` might not be equivalent")') trim(units1), trim(units2); call goErr
        write (gol,'("(ambigious standard units `",a,"`)")') trim(units1_standard); call goErr
        TRACEBACK; status=1; return

      else

        ! get conversion factor:
        call UDUnits_ConversionFactor( trim(units1), trim(units2), factor, status )
        if ( status /= UDUNITS_NOERR ) then
          gol=trim(UDUnits_StrError(status)); call goErr
          write (gol,'("could not obtain conversion factor from `",a,"` to `",a,"`")') &
                           trim(units1), trim(units2); call goErr
          TRACEBACK; status=1; return
        end if
        
      end if  ! both x/x ?

    end if  ! unit descriptions are different
    
    ! ok
    status = 0
    
  end subroutine LE_Data_ConversionFactor
  
  
  ! ***
  
  
  !
  ! Convert climatology description 'year day'
  ! into logical array :  (/ .true., .false., .true., .false., .false., .false. /)
  !
  
  subroutine LE_Data_ClimatMask( climat, cmask, status )
  
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)            ::  climat
    logical, intent(out)                    ::  cmask(6)
    integer, intent(out)                    ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Data_ClimatMask'

    ! --- local ----------------------------------

    character(len=1024)   ::  line
    character(len=64)     ::  element

    ! --- begin ----------------------------------

    ! init climat mask:
    cmask = .false.
    ! climat mask?
    if ( len_trim(climat) > 0 ) then
      ! copy:
      line = trim(climat)
      ! loop over elements:
      do
        ! leave if empty:
        if ( len_trim(line) == 0 ) exit
        ! current:
        call goReadFromLine( line, element, status, sep=' ' )
        IF_NOTOK_RETURN(status=1)
        ! switch:
        select case ( trim(element) )
          case ( 'year'  ) ; cmask(1) = .true.   ! no specific years
          case ( 'month' ) ; cmask(2) = .true.   ! no specific months
          case ( 'day'   ) ; cmask(3) = .true.   ! no specific days
          case default
            write (gol,'("unsupported element `",a,"` in climats line: ",a)') &
                             trim(element), trim(climat); call goErr
            TRACEBACK; status=1; return
        end select
      end do  ! climat elements
    end if ! climat defined
    
    ! ok
    status = 0
    
  end subroutine LE_Data_ClimatMask


end module LE_Data_Common
