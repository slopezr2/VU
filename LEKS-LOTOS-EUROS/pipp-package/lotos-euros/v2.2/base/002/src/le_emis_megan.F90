!###############################################################################
!
! NAME
!   LE_Emis_MEGAN - Emission module around MEGAN 2.04
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

module LE_Emis_MEGAN

  use GO, only : gol, goPr, goErr

  use LE_Emis_MEGAN_Input, only : T_LE_Emis_MEGAN_Input

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Emis_MEGAN

  public  ::  LE_Emis_MEGAN_Init, LE_Emis_MEGAN_Done
  public  ::  LE_Emis_MEGAN_Setup


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_MEGAN'


  ! --- local --------------------------------

  ! emission data base
  type T_Emis_MEGAN
      ! label assigned to this emission:
    character(len=32)                         :: label
    ! settings:
    character(len=512)                        ::  rcfile
    character(len=64)                         ::  rckey
    ! emission factor maps:
    real, allocatable                         ::  EF(:,:,:,:)  ! (nx,ny,N_MECH_SPEC,N_MGN_SPC)
    ! LAI input:
    integer                                   ::  LAI_year
    type(T_LE_Emis_MEGAN_Input)               ::  LAI_prev_LMI
    type(T_LE_Emis_MEGAN_Input)               ::  LAI_curr_LMI
    ! mapping from mechanism to LE species:
    integer, allocatable                      ::  mech2le(:)  ! (N_MECH_SPC)
  end type T_Emis_MEGAN



contains


  ! ===============================================================


  subroutine LE_Emis_MEGAN_Init( em, label, rcF, rckey, t, status )

    use GO                 , only : TrcFile, ReadRc
    use GO                 , only : TDate
    use GO                 , only : goMatchValues
    use MEGAN_2_04         , only : N_MGN_SPC, MGN_SPC
    use MEGAN_2_04         , only : N_MECH_SPC, MECH_SPC
    use MEGAN_2_04         , only : N_PFT, PFT_NC_FNAME, PFT_NC_VNAME
    use MEGAN_2_04         , only : MEGAN_2_04_Init
    use MEGAN_2_04         , only : MEGAN_2_04_EF_FROM_PFT
    use MEGAN_2_04         , only : MEGAN_2_04_EF_FROM_MAP
    use LE_Emis_MEGAN_Input, only : LE_Emis_MEGAN_Input_Init, LE_Emis_MEGAN_Input_Done, LE_Emis_MEGAN_Input_Setup
    use Dims               , only : nx, ny
    use Indices

    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    type(T_Emis_MEGAN), intent(out)     ::  em
    character(len=*), intent(in)        ::  label
    type(TrcFile), intent(in)           ::  rcF
    character(len=*), intent(in)        ::  rckey
    type(TDate), intent(in)             ::  t       ! start time
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_MEGAN_Init'

    ! flags to identify how emissions factor maps are set:
    integer, parameter  ::  iEFsrc_PFT = 1
    integer, parameter  ::  iEFsrc_MAP = 2
    
    ! --- local ------------------------------

    character(len=1024)               ::  line
    character(len=1024)               ::  fname
    character(len=32)                 ::  fc_ext
    character(len=512)                ::  vname
    character(len=4)                  ::  comp
    character(len=32)                 ::  component
    character(len=2)                  ::  mt
    integer                           ::  year
    type(T_LE_Emis_MEGAN_Input)       ::  LMI
    integer                           ::  I_PFT
    real, allocatable                 ::  PFT_prcnt(:,:,:)   ! (nx,ny,N_PFT)
    integer                           ::  n
    integer, allocatable              ::  inds(:)
    character(len=4), allocatable     ::  names(:)
    integer                           ::  imap
    integer, allocatable              ::  EFsrc(:)
    integer                           ::  I_MGN
    integer                           ::  i, j
    integer                           ::  I_MECH
    integer                           ::  ispec

    ! --- begin -------------------------------

    call LE_Data_Enable( 'lon', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lat', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'tsurf', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'ssrd', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'tsurf_ema', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'ssrd_ema', status )
    IF_NOTOK_RETURN(status=1)
    
    ! store label:
    em%label = trim(label)
    
    ! store info on settings:
    em%rcfile = trim(rcF%fname)
    em%rckey = trim(rckey)
    
    ! init MEGAN module:
    call MEGAN_2_04_Init( status )
    IF_NOTOK_RETURN(status=1)
    
    ! final storage for emission factors including
    ! conversion to tracers of the chemical mechanism:
    allocate( em%EF(nx,ny,N_MECH_SPC,N_MGN_SPC) )

    ! storage for PFT fractions:
    allocate( PFT_prcnt(nx,ny,N_PFT) )
    ! read filename template for PFT files:
    call ReadRc( rcF, 'le.emis.megan.pft.fname', fname, status )
    IF_NOTOK_RETURN(status=1)
    ! variable names:
    call ReadRc( rcF, 'le.emis.megan.pft.vname', vname, status )
    IF_NOTOK_RETURN(status=1)
    ! PFT fields for fixed year:
    call ReadRc( rcF, 'le.emis.megan.pft.year', year, status )
    IF_NOTOK_RETURN(status=1)
    ! loop over variables:
    do I_PFT = 1, N_PFT
      ! access input files:
      call LE_Emis_MEGAN_Input_Init( LMI, fname, vname, status )
      IF_NOTOK_RETURN(status=1)
      ! get field on LE grid, read and regrid if necessary:
      call LE_Emis_MEGAN_Input_Setup( LMI, status, &
                                       comp=trim(PFT_NC_FNAME(I_PFT)), &
                                       component=trim(PFT_NC_VNAME(I_PFT)), &
                                       year=year, units='%' )
      IF_NOTOK_RETURN(status=1)
      ! store:
      PFT_prcnt(:,:,I_PFT) = LMI%data  ! %
      ! done:
      call LE_Emis_MEGAN_Input_Done( LMI, status )
      IF_NOTOK_RETURN(status=1)
    end do
    
    ! set origin of emission factors:
    allocate( EFsrc(N_MGN_SPC) )
    ! by default from PFT maps:
    EFsrc(:) = iEFsrc_PFT

    ! read list with MEGAN representer species 
    ! for which emisison factor map should be read;
    ! temporary variables:
    allocate( names(N_MGN_SPC) )
    allocate( inds (N_MGN_SPC) )
    ! read line with names:
    call ReadRc( rcF, 'le.emis.megan.ef.maps', line, status )
    IF_NOTOK_RETURN(status=1)
    ! compare with list of representer species names:
    call goMatchValues( line , MGN_SPC, n, names, inds, status )
    IF_NOTOK_RETURN(status=1)
    ! maps needed ?
    if ( n > 0 ) then
      ! reset origin flag:
      do i = 1, n
        EFsrc(inds(i)) = iEFsrc_MAP
     end do
    end if
    ! clear:
    deallocate( names )
    deallocate( inds  )
    
    ! fill maps, loop over MEGAN representer species:
    do I_MGN = 1, N_MGN_SPC
      ! how to fill ?
      select case ( EFsrc(I_MGN) )
        !
        ! from PFT maps:
        !
        case ( iEFsrc_PFT )

          ! compute emission factor maps, convert to tracers of the chemical mechanism:
          call MEGAN_2_04_EF_FROM_PFT( I_MGN, nx, ny, PFT_prcnt, em%EF(:,:,:,I_MGN), status )
          IF_NOTOK_RETURN(status=1)

        !
        ! read:
        !
        case ( iEFsrc_MAP )

          ! read filename template for EF files:
          call ReadRc( rcF, 'le.emis.megan.ef.fname', fname, status )
          IF_NOTOK_RETURN(status=1)
          ! variable names:
          call ReadRc( rcF, 'le.emis.megan.ef.vname', vname, status )
          IF_NOTOK_RETURN(status=1)
          ! EF fields for fixed year:
          call ReadRc( rcF, 'le.emis.megan.ef.year', year, status )
          IF_NOTOK_RETURN(status=1)
          ! set emission specific variable and filename parts and flags:
          select case ( trim(MGN_SPC(I_MGN)) )
            !~ available for v2.0:
            case ( 'ISOP' ) ; component = 'Isoprene'      ; comp = 'iso'  ; mt = ''
            !!~ available for v2.1:
            !case ( 'MYRC' ) ; component = 'myrcene'       ; comp = 'myr'  ; mt = 'MT'
            !case ( 'SABI' ) ; component = 'sabinene'      ; comp = 'sab'  ; mt = 'MT'
            !case ( 'LIMO' ) ; component = 'limonene'      ; comp = 'lim'  ; mt = 'MT'
            !case ( '3CAR' ) ; component = 'carene'        ; comp = 'car'  ; mt = 'MT'
            !case ( 'OCIM' ) ; component = 'ocimene'       ; comp = 'oci'  ; mt = 'MT'
            !case ( 'APIN' ) ; component = 'a_pinene'      ; comp = 'pina' ; mt = 'MT'
            !case ( 'BPIN' ) ; component = 'b_pinene'      ; comp = 'pinb' ; mt = 'MT'
            !case ( 'MBO'  ) ; component = 'methylbutenol' ; comp = 'mbo'  ; mt = ''
            !case ( 'MEOH' ) ; component = 'methanol'      ; comp = 'meoh' ; mt = ''
            !case ( 'NO'   ) ; component = 'nitric_oxide'  ; comp = 'no'   ; mt = ''
            !~ not yet ...
            case default
              write (gol,'("could not set components of EF file for `",a,"` emis")') trim(MGN_SPC(I_MGN)); call goErr
              TRACEBACK; status=1; return
          end select
          ! access input of emission factors:
          call LE_Emis_MEGAN_Input_Init( LMI, fname, vname, status )
          IF_NOTOK_RETURN(status=1)
          ! get field on LE grid, read and regrid if necessary:
          call LE_EMIS_MEGAN_Input_Setup( LMI, status, &
                        year=year, comp=comp, component=component, mt=mt, &
                        units='ug m-2 h-1' )
          IF_NOTOK_RETURN(status=1)
          ! convert to emission factor per components of chemical mechanism:
          call MEGAN_2_04_EF_FROM_MAP( I_MGN, nx, ny, LMI%data, PFT_prcnt, em%EF(:,:,:,I_MGN), status )
          IF_NOTOK_RETURN(status=1)
          ! clear:
          call LE_Emis_MEGAN_Input_Done( LMI, status )
          IF_NOTOK_RETURN(status=1)
          
        case default
          write (gol,'("unsupported EF source number ",i2," for repr.species ",i2)') EFsrc(I_MGN), I_MGN; call goErr
          TRACEBACK; status=1; return
      end select

    end do   ! representer species
    
    ! read filename template for LAI files:
    call ReadRc( rcF, 'le.emis.megan.lai.fname', fname, status )
    IF_NOTOK_RETURN(status=1)
    ! variable names:
    call ReadRc( rcF, 'le.emis.megan.lai.vname', vname, status )
    IF_NOTOK_RETURN(status=1)
    ! LAI fields for fixed year:
    call ReadRc( rcF, 'le.emis.megan.lai.year', em%LAI_year, status )
    IF_NOTOK_RETURN(status=1)
    ! access input for previous month:
    call LE_Emis_MEGAN_Input_Init( em%LAI_prev_LMI, fname, vname, status )
    IF_NOTOK_RETURN(status=1)
    ! access input for current month:
    call LE_Emis_MEGAN_Input_Init( em%LAI_curr_LMI, fname, vname, status )
    IF_NOTOK_RETURN(status=1)
            
    ! mapping from chemical mechanism to LE species:
    allocate( em%mech2le(N_MECH_SPC) )
    ! loop over original species:
    do I_MECH = 1, N_MECH_SPC
      ! try to match:
      select case ( trim(MECH_SPC(I_MECH)) )
        case ( 'ISOP'  ) ; ispec = i_iso     ! Isoprene
        case ( 'NO'    ) ; ispec = i_no      ! Nitric oxide                                  
        case ( 'NO2'   ) ; ispec = i_no2     ! Nitrogen dioxide                              
        case ( 'CO'    ) ; ispec = i_co      ! Carbon monoxide                               
        case ( 'CH3OH' ) ; ispec = -1
        case ( 'ANOL'  ) ; ispec = -1        ! ehtanol ?
        case ( 'ALD2'  ) ; ispec = i_ald     ! Aldehyde             
        case ( 'HCHO'  ) ; ispec = i_form    ! Formaldehyde        
        case ( 'HCOOH' ) ; ispec = -1        ! carbon acid
        case ( 'RCOOH' ) ; ispec = -1        ! carboneous acid
        case ( 'CH4'   ) ; ispec = i_ch4     ! Methane       
        case ( 'C2H6'  ) ; ispec = i_par     ! ethane : assigned to PAR
        case ( 'PAR'   ) ; ispec = i_par     ! Paraffin carbon bond (C-C)    
        case ( 'ETH'   ) ; ispec = i_eth     ! Ethene  
        case ( 'OLET'  ) ; ispec = i_ole     ! Terminal olefin carbon bond (R-C=C)
        case ( 'OLEI'  ) ; ispec = i_ole     ! Terminal olefin carbon bond (R-C=C)
        case ( 'CRES'  ) ; ispec = i_cres    ! Cresol and higher molecular weight phenols    
        case ( 'AONE'  ) ; ispec = -1
        case ( 'TOL'   ) ; ispec = i_tol     ! Toluene and other monoalkyl aromatics    
        case ( 'XYL'   ) ; ispec = i_xyl     ! Xylene and other polyalkyl aromatics          
        case ( 'DMS'   ) ; ispec = -1
        case ( 'NH3'   ) ; ispec = i_nh3     ! Ammonia            
        case ( 'NR'    ) ; ispec = -1
        case ( 'TERP'  ) ; ispec = i_terp    ! terpene
        case default
          write (gol,'("could not map mechanism species `",a,"` to LE species")') trim(MECH_SPC(I_MECH)); call goErr
          TRACEBACK; status=1; return
      end select
      ! store:
      em%mech2le(I_MECH) = ispec
    end do

    ! ok
    status = 0

  end subroutine LE_Emis_MEGAN_Init


  ! ***


  subroutine LE_Emis_MEGAN_Done( em, status )

    use MEGAN_2_04         , only : N_MGN_SPC
    use LE_Emis_MEGAN_Input, only : LE_Emis_MEGAN_Input_Done

    ! --- in/out ------------------------------

    type(T_Emis_MEGAN), intent(inout)     ::  em
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_MEGAN_Done'

    ! --- local -------------------------------

    integer                     ::  imap

    ! --- begin -------------------------------
      
    ! done with input:
    call LE_Emis_MEGAN_Input_Done( em%LAI_prev_LMI, status )
    IF_NOTOK_RETURN(status=1)
    ! done with input:
    call LE_Emis_MEGAN_Input_Done( em%LAI_curr_LMI, status )
    IF_NOTOK_RETURN(status=1)

    ! clear:
    deallocate( em%EF )
    deallocate( em%mech2le  )

    ! ok
    status = 0

  end subroutine LE_Emis_MEGAN_Done


  ! ===============================================================


  subroutine LE_Emis_MEGAN_Setup( em, emis_a, t1, t2, status )

    use GO                 , only : TDate, Get, operator(+), operator(-), operator(/)
    use GO                 , only : DayNumber
    use Grid               , only : AreaOper
    use dims               , only : nx, ny, nz, nspec
    use MEGAN_2_04         , only : N_MGN_SPC, MGN_SPC
    use MEGAN_2_04         , only : N_MECH_SPC
    use MEGAN_2_04         , only : MEGAN_2_04_GET_GAMMA
    use MEGAN_2_04         , only : MEGAN_2_04_GET_RHO
    use LE_Emis_MEGAN_Input, only : LE_Emis_MEGAN_Input_Setup
    use LE_Grid            , only : ugg
    use LE_Data, only : LE_Data_GetPointer

    ! --- in/out ---------------------------

    type(T_Emis_MEGAN), intent(inout)     ::  em
    real, intent(inout)                   ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_MEGAN_Setup'

    ! Parameterization to compute emission activity factor (umol/m2/s)
    ! from SRAD (short wave from sun (W/m2))
    ! assuming 4.766 (umol m-2 s-1) per (W m-2)
    ! assume 1/2 of SRAD is in 400-700nm band
    ! Taken from original MEGAN code.
    real, parameter  ::  ppfd_umol_per_J          = 4.766 ! (umol m-2 s-1) per (W m-2)
    real, parameter  ::  frac_400nm_700nm_in_ssrd = 0.5   ! (1) fraction in 400-700nm band

    ! parameter for unit conversion
    real, parameter   ::  minutes_per_hour = 60.0
    real, parameter   ::  units_per_micro  = 1.0e-6

    ! --- local ----------------------------

    type(TDate)                      ::  tmid
    integer                          ::  year, month, day, hour
    integer                          ::  month_prev
    real, allocatable                ::  PPFD(:,:)           ! (nx,ny)
    real, allocatable                ::  PPFD_dayaver(:,:)   ! (nx,ny)
    real, allocatable                ::  ER_MECH(:,:,:)      ! (nx,ny,N_MECH_SPC)
    integer                          ::  I_MGN
    integer                          ::  I_MECH
    integer                          ::  IDATE   ! YYYYDDD
    integer                          ::  ITIME   ! HHMMSS
 
    real, allocatable                ::  GAMMA(:,:) ! combined factors

    character(len=32)                ::  units

    real                             ::  RHO             ! Production and loss within canopy

    integer                          ::  ispec
    integer                          ::  ix,iy

    ! grid data:
    real, pointer        ::  lons(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  lats(:,:,:)   ! (lon,lat,1)
    ! grid data:
    real, pointer        ::  t2m(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  ssrd(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  t2m_dayaver(:,:,:)   ! (lon,lat,1)
    real, pointer        ::  ssrd_dayaver(:,:,:)   ! (lon,lat,1)
    
   
    ! --- begin ----------------------------

    call LE_Data_GetPointer( 'lon', lons, status, check_units ='degrees_east' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lat', lats, status, check_units ='degrees_north' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'tsurf', t2m, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'ssrd', ssrd, status, check_units ='W/m2' )
    IF_NOTOK_RETURN(status=1)
    
    !
    call LE_Data_GetPointer( 'tsurf_ema', t2m_dayaver, status, check_units ='K' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'ssrd_ema', ssrd_dayaver, status, check_units ='W/m2' )
    IF_NOTOK_RETURN(status=1)
    !

    ! extract time parameters:
    tmid = t1 + (t2-t1)/2
    call Get( tmid, year=year, month=month, day=day, hour=hour )

    ! alternative time values:
    IDATE = t1%year*1000 + DayNumber(t1)     ! YYYYDDD
    ITIME = t1%hour*100*100 + t1%min*100 + t1%sec   ! HHMMSS
    
    allocate( PPFD( ugg%nlon, ugg%nlat) )
    allocate( PPFD_dayaver( ugg%nlon, ugg%nlat) )
    allocate( GAMMA  (ugg%nlon, ugg%nlat) )
    allocate( ER_MECH(ugg%nlon, ugg%nlat,N_MECH_SPC) )
    
    ! Compute PPFD emission activity factors from radiance:
    PPFD         = ssrd(:,:,1)         * ppfd_umol_per_J * frac_400nm_700nm_in_ssrd  ! umol/m2/s
    PPFD_dayaver = ssrd_dayaver(:,:,1) * ppfd_umol_per_J * frac_400nm_700nm_in_ssrd  ! umol/m2/s

    ! Get LAI (m2/m2) for the previous month and the current month
    ! previous month:
    month_prev = t1%month - 1
    if ( month_prev == 0 ) month_prev = 12
    ! get field on LE grid, read and regrid if necessary:
    call LE_Emis_MEGAN_Input_Setup( em%LAI_prev_LMI, status, &
                                     year=em%LAI_year, month=month_prev, &
                                     units='m2 m-2' )
    IF_NOTOK_RETURN(status=1)
    ! similar for current month:
    call LE_Emis_MEGAN_Input_Setup( em%LAI_curr_LMI, status, &
                                     year=em%LAI_year, month=t1%month, &
                                     units='m2 m-2' )
    IF_NOTOK_RETURN(status=1)

    ! init result:
    ER_MECH = 0.0
    ! loop over representer species:
    do I_MGN = 1, N_MGN_SPC

      ! get factors:
      call MEGAN_2_04_GET_GAMMA( I_MGN, IDATE, ITIME, &
                         nx, ny, lons(:,:,1), lats(:,:,1), &
                         em%LAI_prev_LMI%data, em%LAI_curr_LMI%data, &
                         t2m, t2m_dayaver, &
                         PPFD, PPFD_dayaver, &
                         GAMMA, status )
      IF_NOTOK_RETURN(status=1)

      ! get factor for production and loss within plant canopies:
      call MEGAN_2_04_GET_RHO( MGN_SPC(I_MGN), RHO, status )
      IF_NOTOK_RETURN(status=1)
      
      ! loop over species from chemical mechanism:
      do I_MECH = 1, N_MECH_SPC
        ! add contribution to emission rate:
        !       umol/m2/hr          umol/m2/hr          umol/m2/hr              1      1
        ER_MECH(:,:,I_MECH) = ER_MECH(:,:,I_MECH) + em%EF(:,:,I_MECH,I_MGN) * GAMMA * RHO
      end do

    end do   ! representer species

    ! convert from units/m2 to units using area per cell:
    do I_MECH = 1, N_MECH_SPC
      call ugg%AreaOper( ER_MECH(:,:,I_MECH), '*', 'm2', status )
      IF_NOTOK_RETURN(status=1)
    end do

    ! convert from umol/hr to mol/min
    !          umol/hr        1/u       /      (min/hr) :
    ER_MECH = ER_MECH * units_per_micro / minutes_per_hour   !  mol/min
    
    ! add to LE emission array:
    do I_MECH = 1, N_MECH_SPC
      ! target tracer:
      ispec = em%mech2le(I_MECH)
      ! undefined ? then try next:
      if ( ispec < 0 ) cycle
      ! add contribution to surface layer:
      emis_a(:,:,1,ispec) = emis_a(:,:,1,ispec) + ER_MECH(:,:,I_MECH)   ! mol/min
    
    end do
    
    ! clear:
    deallocate( PPFD )
    deallocate( PPFD_dayaver )
    deallocate( GAMMA   )
    deallocate( ER_MECH )

    ! ok
    status = 0

  end subroutine LE_Emis_MEGAN_Setup


end module LE_Emis_MEGAN

