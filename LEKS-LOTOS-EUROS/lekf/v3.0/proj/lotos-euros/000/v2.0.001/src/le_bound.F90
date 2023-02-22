!###############################################################################
!
! LE_Bound - LOTOS-EUROS boundary condition routines
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

module LE_Bound

  use GO, only : gol, goErr, goPr
  
  use LE_Bound_EMEP_Monthly, only : T_LE_Bound_EMEP_Monthly
  use LE_Bound_Data        , only : T_LE_Bound_Data

  implicit none

  ! --- in/out --------------------------------

  private

  public  ::  LE_Bound_Init, LE_Bound_Done
  public  ::  LE_Bound_Setup
  public  ::  LE_Bound_Initial



  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound'

  ! maximum number of bc types:
  integer, parameter    ::  max_bc_type = 10
  
  ! print messages ?
  logical, parameter    ::  verbose = .false.
  !logical, parameter    ::  verbose = .true.
  

  ! --- types --------------------------------
  
  type T_LE_Bound
    ! set name:
    character(len=16)           ::  typ
    ! storage:
    type(T_LE_Bound_EMEP_Monthly)   ::  emep_monthly
    ! storage:
    type(T_LE_Bound_Data)           ::  bdata
  end type


  ! --- local --------------------------------

  ! list with boundary conditions to use:
  integer              ::  n_bc_type
  character(len=16)    ::  le_bound_type(max_bc_type)
  ! storage per type:
  type(T_LE_Bound), allocatable   ::  bnds(:)

  ! boundary conditions from Isaksen model:
  character(len=512)    ::  isak_path
  integer               ::  isak_month_curr

  ! Logan ozone climatology:
  character(len=512)    ::  logan_path

contains


  ! ===============================================================


  subroutine LE_Bound_Init( rcfile, t, status )

    use GO, only : TrcFile, Init, Done, ReadRc
    use GO, only : TDate
    use GO, only : goReadFromLine

    use dims               , only : runF
    use dims               , only : bc_west, bc_east, bc_south, bc_north
    use dims               , only : caloft
    use LE_Bound_Clim_Isak , only : LE_Bound_Clim_Isak_Init, read_bc
    use LE_Bound_Clim_Logan, only : LE_Bound_Clim_Logan_Init
    use LE_Bound_Clim_EMEP , only : LE_Bound_Clim_EMEP_Init
    use LE_Bound_LE        , only : LE_Bound_LE_Init
    use LE_Bound_EMEP_Monthly, only : LE_Bound_EMEP_Monthly_Init
    use LE_Data, only : LE_Data_Enable

    ! --- in/out ------------------------------

    character(len=*), intent(in)    ::  rcfile
    type(TDate), intent(in)         ::  t
    integer, intent(out)            ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_Init'

    ! --- local ------------------------------

    type(TRcFile)         ::  rcF
    character(len=512)    ::  btypes
    character(len=16)     ::  btype
    integer               ::  i_bc_type
    character(len=64)     ::  rckey

    ! --- begin -------------------------------
    
    ! info ...
    write (gol,'(a,": initialize boundary conditions ...")') rname; call goPr

    ! used for regridding:
    call LE_Data_Enable( 'oro', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'h', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 't', status )
    IF_NOTOK_RETURN(status=1)

    ! open rcfile:
    call Init( rcF, rcfile, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read item:
    call ReadRc( rcF, 'le.bound.types', btypes, status )
    IF_NOTOK_RETURN(status=1)

    ! setup list:
    n_bc_type = 0
    do
      ! empty ? then leave
      if ( len_trim(btypes) == 0 ) exit
      ! read part from line
      call goReadFromLine( btypes, btype, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      ! storage ok ?
      if ( n_bc_type == max_bc_type ) then
        write (gol,'("storage too small; increase max_bc_type")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! store:
      n_bc_type = n_bc_type + 1
      le_bound_type(n_bc_type) = trim(btype)
    end do
    
    ! storage:
    allocate( bnds(n_bc_type) )

    ! loop over requested bc types:
    do i_bc_type = 1, n_bc_type
    
      ! store type:
      bnds(i_bc_type)%typ = trim(le_bound_type(i_bc_type))
      
      ! current key:
      write (rckey,'("le.bound.",a)') trim(bnds(i_bc_type)%typ)

      ! which type now ?
      select case ( le_bound_type(i_bc_type) )

        case ( 'clim-isak' )

          ! allocate some boundary condition arrays:
          call LE_Bound_Clim_Isak_Init( status )
          IF_NOTOK_RETURN(status=1)

          ! read path:
          call ReadRc( rcF, 'le.bound.clim-isak.path', isak_path, status )
          IF_NOTOK_RETURN(status=1)

          ! reading boundary conditions
          call read_bc( t%month, isak_path, status )
          IF_NOTOK_RETURN(status=1)
          ! bc read for this month:
          isak_month_curr = t%month

        case ( 'clim-const' )

          ! nothing to init

        case ( 'clim-logan' )

          ! read path:
          call ReadRc( rcF, 'le.bound.clim-logan.path', logan_path, status )
         IF_NOTOK_RETURN(status=1)

          ! read the logan boundary conditions for O3
          call LE_Bound_Clim_Logan_Init( logan_path, status )
          IF_NOTOK_RETURN(status=1)

        case ( 'clim-emep' )

          ! setup the emep boundary conditions:
          call LE_Bound_Clim_EMEP_Init( status )
          IF_NOTOK_RETURN(status=1)

        case ( 'clim-steady' )

          ! nothing to init

       case ( 'emep-monthly' )
          ! init boundary conditions from MOZART:
          call LE_Bound_EMEP_Monthly_Init( bnds(i_bc_type)%emep_monthly, rcfile, rckey, status )
          IF_NOTOK_RETURN(status=1)

       case ( 'data' )
           ! boundary conditions from data variables:
           call bnds(i_bc_type)%bdata%Init( rcF, status )
           IF_NOTOK_RETURN(status=1)

        case ( 'le' )

          ! read boundary conditions from previous (course grid) LE run;
          call LE_Bound_LE_Init( rcfile, status )
          IF_NOTOK_RETURN(status=1)
          
       case default

          ! error ...
          write (gol,'("found unsupported boundary condition type in list : ")'); call goErr
          write (gol,'("  item nr   : ",i6)') i_bc_type; call goErr
          write (gol,'("  type      : ",a)') trim(le_bound_type(i_bc_type)); call goErr
          TRACEBACK; status=1; return

      end select

    end do  ! bc types

    ! close rcfile:
    call Done( rcF, status )
    IF_NOTOK_RETURN(status=1)
    
    ! force all boundary concentrations to zero:
    bc_west  = 0.0
    bc_east  = 0.0
    bc_north = 0.0
    bc_south = 0.0
    caloft   = 0.0
    ! this used to be done every timestep,
    ! but this is now skipped since we need to maintain
    ! an aerosol composition to distribute total-pm10
    ! values over the original components;
    ! this is done by copying the nearby concentrations
    ! into the boundary arrays at the begin of the timestep boundary conditions

    ! ok
    status = 0

  end subroutine LE_Bound_Init


  ! ***


  subroutine LE_Bound_Done( status )

    use LE_Bound_LE        , only : LE_Bound_LE_Done
    use LE_Bound_EMEP_Monthly, only : LE_Bound_EMEP_Monthly_Done
#ifdef with_labeling
    use SA_Labeling        , only : SA_Bound_Done
#endif    

    ! --- in/out ------------------------------

    integer, intent(out)  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_Done'

    ! --- local -------------------------------

    integer               ::  i_bc_type

    ! --- begin -------------------------------

    ! loop over requested bc types:
    do i_bc_type = 1, n_bc_type
      ! which type now ?
      select case ( le_bound_type(i_bc_type) )
        case ( 'clim-isak', 'clim-const', 'clim-logan', 'clim-emep', 'clim-steady' )
          ! nothing to be done
          ! ... or nothing implemented
        case ( 'emep-monthly' )
          ! done with emep-monthly boundary conditions
          call LE_Bound_EMEP_Monthly_Done( bnds(i_bc_type)%emep_monthly, status )
          IF_NOTOK_RETURN(status=1)
       case ( 'data' )
          ! done with boundary conditions from data variables:
          call bnds(i_bc_type)%bdata%Done( status )
          IF_NOTOK_RETURN(status=1)
        case ( 'le' )
          ! done with boundary conditions from previous LOTOS-EUROS run:
          call LE_Bound_LE_Done( status )
          IF_NOTOK_RETURN(status=1)
#ifdef with_labeling
          call SA_Bound_Done( status )
          IF_NOTOK_RETURN(status=1)
#endif                    
        case default
          ! error ...
          write (gol,'("unsupported bc type : ")'); call goErr
          write (gol,'("  nr     : ",i6)') i_bc_type; call goErr
          write (gol,'("  type   : ",a)') trim(le_bound_type(i_bc_type)); call goErr
          TRACEBACK; status=1; return
      end select
    end do  ! bc types
    
    ! clear:
    deallocate( bnds )
    
    ! ok
    status = 0

  end subroutine LE_Bound_Done


  ! ***


  subroutine LE_Bound_Setup( t1, t2, nhour, status )

    use GO, only : TDate, wrtgol, DayNumber
    use GO, only : operator(+), operator(-), operator(/)

    use dims           , only : runF
    use dims           , only : nx, ny, nz, nspec
    use indices        , only : specname
    use dims           , only : bc_west, bc_east, bc_south, bc_north, caloft
    use dims           , only : rk1rk3

    use LE_Bound_Clim_Isak   , only : read_bc, makebc
    use LE_Bound_Clim_Const  , only : LE_Bound_Clim_Const_Setup
    use LE_Bound_Clim_Logan  , only : get_logan
    use LE_Bound_Clim_EMEP   , only : get_emep_bc
    use LE_Bound_Clim_Steady , only : bc_steady
    use LE_Bound_LE          , only : LE_Bound_LE_Get
    use LE_Bound_EMEP_Monthly, only : LE_Bound_EMEP_Monthly_Get
    use LE_Bound_Tools       , only : LE_Bound_Fill
#ifdef with_labeling
    use SA_Labeling          , only : SA_Bound_Setup
#endif    

    ! --- in/out ------------------------------

    type(TDate), intent(in)   ::  t1, t2
    integer, intent(in)       ::  nhour     ! hours since start of forecast or midnight (for t1)
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_Setup'

    ! --- local -------------------------------

    integer               ::  i_bc_type
    integer               ::  iday
    integer               ::  ispec
    logical               ::  bc_filled(nspec)
    logical               ::  filled(nspec)

    ! boundary conditions arrays:
    real, pointer       ::  bc_lons(:)
    real, pointer       ::  bc_lats(:)
    real, pointer       ::  bc_phlev(:,:,:)
    real, pointer       ::  bc_ahlev(:,:,:)
    real, pointer       ::  bc_vmr(:,:,:)
    character(len=32)   ::  bc_unit
    ! converstion:
    real                ::  molemass_tracer, molemass_air
    ! time stamp:
    type(TDate)         ::  tmid

    ! --- begin -------------------------------

    ! info ...
    call wrtgol('LE:   setup boundary conditions for ',t1,' - ',t2); call goPr

    ! nothing filled yet ...
    bc_filled = .false.  ! no boundary conditions filled

    ! loop over requested bc types:
    do i_bc_type = 1, n_bc_type

      ! no fields filled yet:
      nullify( bc_lons  )
      nullify( bc_lats  )
      nullify( bc_phlev )
      nullify( bc_ahlev )
      nullify( bc_vmr   )

      ! which type now ?
      select case ( le_bound_type(i_bc_type) )

        !
        ! Isaksen model
        !

        case ( 'clim-isak' )

          ! chemistry ?
          ! update for new month ?
          if ( t1%month /= isak_month_curr ) then
            call read_bc( t1%month, isak_path, status )
            IF_NOTOK_RETURN(status=1)
            isak_month_curr = t1%month
          end if
          ! fill boundary conditions:
          call makebc( bc_filled, status )
          IF_NOTOK_RETURN(status=1)

        !
        ! contant values
        !

        case ( 'clim-const' )

          ! fill some constant values:
          call LE_Bound_Clim_Const_Setup( bc_filled )


        !
        ! Logan ozone climatology
        !

        case ( 'clim-logan' )

          ! replace O3 bound. cond. by Logan BC for O3
          call get_logan( t1%year, t1%month, bc_filled, status )
          IF_NOTOK_RETURN(status=1)

        !
        ! EMEP boundary conditions
        !

        case ( 'clim-emep' )

          ! day number in year:
          iday = DayNumber( t1 )
          ! setup EMEP boundary conditions:
          call get_emep_bc( t1%year, iday, bc_filled, status )
          IF_NOTOK_RETURN(status=1)

        !
        ! steady state o3/no/no2
        !

        case ( 'clim-steady' )

          ! put BC in steady state
          call bc_steady( rk1rk3, bc_filled )

        !
        ! data variables
        !
        
        case ( 'data' )

          ! setup boundary conditions from data variables:
          call bnds(i_bc_type)%bdata%FillBounds( filled, status )
          IF_NOTOK_RETURN(status=1)
          
          ! combine flags:
          bc_filled = bc_filled .or. filled

        !
        ! LOTOS-EUROS run
        !

        case ( 'le' )

          ! read boundary conditions from previous (course) le run:
          call LE_Bound_LE_Get( t1, bc_filled, status )
          IF_NOTOK_RETURN(status=1)
#ifdef with_labeling
          call SA_Bound_Setup(t1, status)
          IF_NOTOK_RETURN(status=1)
#endif          
        !
        ! ~~ other ...
        !

        case default

          ! error ...
          write (gol,'("unsupported bc type : ")'); call goErr
          write (gol,'("  nr     : ",i6)') i_bc_type; call goErr
          write (gol,'("  type   : ",a)') trim(le_bound_type(i_bc_type)); call goErr
          TRACEBACK; status=1; return

      end select

      ! clear if necessary:
      if ( associated(bc_lons ) ) deallocate( bc_lons  )
      if ( associated(bc_lats ) ) deallocate( bc_lats  )
      if ( associated(bc_phlev) ) deallocate( bc_phlev )
      if ( associated(bc_ahlev) ) deallocate( bc_ahlev )
      if ( associated(bc_vmr  ) ) deallocate( bc_vmr   )

    end do  ! bc types

    ! print info ?
    if ( verbose ) then
      ! some missing ?
      if ( .not. all(bc_filled) ) then
        write (gol,'("WARNING - ")'); call goPr
        write (gol,'("WARNING - boundary conditions not explicitly set for some tracers:")'); call goPr
        do ispec = 1, nspec
          if ( .not. bc_filled(ispec) ) then
            write (gol,'("WARNING -   ",i3," ",a)') ispec, trim(specname(ispec)); call goPr
          end if
        end do
        write (gol,'("WARNING - ")'); call goPr
      end if
    end if

    ! ok
    status = 0

  end subroutine LE_Bound_Setup


  ! ***


  subroutine LE_Bound_Initial( c, t, status )

    use GO, only : TDate, wrtgol, DayNumber

    use dims           , only : runF
    use dims           , only : nx, ny, nz, nspec
    use indices        , only : specname
    use dims           , only : bc_west, bc_east

    use LE_Bound_LE        , only : LE_Bound_LE_Initial
    use LE_Bound_Tools     , only : LE_Bound_Fill_Initial
    use LE_Bound_EMEP_Monthly, only : LE_Bound_EMEP_Monthly_Get

    ! --- in/out ------------------------------

    real, intent(inout)       ::  c(nx,ny,nz,nspec)
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Bound_Initial'

    ! --- local -------------------------------

    integer               ::  i_bc_type
    integer               ::  ispec
    logical               ::  ic_filled(nspec)
    logical               ::  filled(nspec)

    ! boundary conditions:
    real, pointer       ::  bc_lons(:)
    real, pointer       ::  bc_lats(:)
    real, pointer       ::  bc_phlev(:,:,:)
    real, pointer       ::  bc_ahlev(:,:,:)
    real, pointer       ::  bc_vmr(:,:,:)
    character(len=32)   ::  bc_unit
    ! converstion
    real                ::  molemass_tracer, molemass_air

    ! initial concentrations:
    integer             ::  i

    ! --- begin -------------------------------

    ! info ...
    call wrtgol('LE:   setup initial concentrations from boundary conditions for ',t); call goPr

    ! nothing filled yet ...
    ic_filled = .false.  ! no initial concentrations filled

    ! loop over requested bc types:
    do i_bc_type = 1, n_bc_type

      ! no fields filled yet:
      nullify( bc_lons  )
      nullify( bc_lats  )
      nullify( bc_phlev )
      nullify( bc_ahlev )
      nullify( bc_vmr   )

      ! which type now ?
      select case ( trim(le_bound_type(i_bc_type)) )

        !
        ! Isaksen model
        !

        case ( 'clim-isak' )

          ! no initial concentrations

        !
        ! contant values
        !

        case ( 'clim-const' )

          ! no initial concentrations

        !
        ! Logan ozone climatology
        !

        case ( 'clim-logan' )

          ! no initial concentrations

        !
        ! EMEP boundary conditions
        !

        case ( 'clim-emep' )

          ! no initial concentrations

        !
        ! steady state o3/no/no2
        !

        case ( 'clim-steady' )

          ! no initial concentrations

        !
        ! emep monthly
        !
        
        case ( 'emep-monthly' )

          ! loop over tracers:
          
          do ispec = 1, nspec
            ! get field with full grid description:
            call LE_Bound_EMEP_Monthly_Get( bnds(i_bc_type)%emep_monthly, &
                      t, 0, ispec, &
                      bc_lons, bc_lats, bc_phlev, bc_ahlev, &
                      bc_vmr, bc_unit, &
                      molemass_tracer, molemass_air, status, &
                      verbose=verbose )
            select case ( status )
              case ( -1 )
                !! not found ...
                !write (gol,'("    WARNING - EMEP Monthly field not available : ",a)') trim(specname(ispec)); call goPr
                !TRACEBACK; status=1; return
                ! try next field:
                cycle
              case ( 0 )
                ! ok
              case default
                ! error ...
                write (gol,'("reading emep-monthly field")'); call goErr
                TRACEBACK; status=1; return
            end select
            ! convert grids and fill initial concentration:
            call LE_Bound_Fill_Initial( c(:,:,:,ispec), ispec, &
                                        molemass_tracer, molemass_air, &
                                        bc_lons, bc_lats, bc_phlev, &
                                        bc_vmr, bc_unit, &
                                        status )
            IF_NOTOK_RETURN(status=1)

            ! set flag:
            ic_filled(ispec) = .true.

          end do  ! le tracers

        !
        ! data variables
        !
        
        case ( 'data' )

          ! setup initial conditions from data variables,
          ! return array with filled concentrations:
          call bnds(i_bc_type)%bdata%FillConc( c, filled, status )
          IF_NOTOK_RETURN(status=1)
          
          ! combine flags:
          ic_filled = ic_filled .or. filled


        !
        ! LOTOS-EUROS run (nc files)
        !

        case ( 'le' )

          ! read boundary conditions from previous (course) le run:
          call LE_Bound_LE_Initial( t, c, ic_filled, status )
          IF_NOTOK_RETURN(status=1)

        !
        ! ~~ other ...
        !

        case default

          ! error ...
          write (gol,'("unsupported bc type : ")'); call goErr
          write (gol,'("  nr     : ",i6)') i_bc_type; call goErr
          write (gol,'("  type   : ",a)') trim(le_bound_type(i_bc_type)); call goErr
          TRACEBACK; status=1; return

      end select

      ! clear if necessary:
      if ( associated(bc_lons ) ) deallocate( bc_lons  )
      if ( associated(bc_lats ) ) deallocate( bc_lats  )
      if ( associated(bc_phlev) ) deallocate( bc_phlev )
      if ( associated(bc_ahlev) ) deallocate( bc_ahlev )
      if ( associated(bc_vmr  ) ) deallocate( bc_vmr   )

    end do  ! bc types

    ! filled initial concentrations ?
    write (gol,'("fill initial concentrations ...")'); call goPr
    ! loop over specs:
    do ispec = 1, nspec
      ! filled yet ?
      if ( ic_filled(ispec) ) then
        ! filled, for example from previous LE run
        write (gol,'("  ",i3," ",a," already filled ...")') ispec, trim(specname(ispec)); call goPr
      else
        ! fill from interpolation
        write (gol,'("  ",i3," ",a," interpolated between west and east ...")') ispec, trim(specname(ispec)); call goPr
        ! loop over longitudes:
        do i = 1, nx
          c(i,:,:,ispec) = bc_west(:,:,ispec)*(1.0-(i-0.5)/nx) + bc_east(:,:,ispec)*(i-0.5)/nx
        end do
      end if
    end do

    ! ok
    status = 0

  end subroutine LE_Bound_Initial


end module LE_Bound

