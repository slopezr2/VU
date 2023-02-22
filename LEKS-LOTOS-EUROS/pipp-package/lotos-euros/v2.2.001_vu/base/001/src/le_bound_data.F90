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

module LE_Bound_Data

  use GO     , only : gol, goErr, goPr
  use LE_Data, only : T_VarSum

  implicit none
  

  ! --- in/out --------------------------------

  private

  public  ::  T_LE_Bound_Data
  

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Bound_Data'


  ! --- types --------------------------------
  
  type T_LE_Bound_Data
    ! number of tracers:
    integer                       ::  ntracer
    ! mapping to global index:
    integer, allocatable          ::  ispec(:)  ! (ntracer)
    ! info on input variables:
    type(T_VarSum), allocatable   ::  varsum(:) ! (ntracer)
  contains
    procedure :: Init          => Bound_Data_Init
    procedure :: Done          => Bound_Data_Done
    procedure :: FillBounds    => Bound_Data_FillBounds
    procedure :: FillConc      => Bound_Data_FillConc
  end type T_LE_Bound_Data



contains


  ! ====================================================================
  ! ===
  ! === Bound_Data
  ! ===
  ! ====================================================================


  subroutine Bound_Data_Init( self, rcF, status )
  
    use GO     , only : TrcFile, ReadRc
    use Indices, only : specname, n_advected, ispecs_advected
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Bound_Data), intent(out)           ::  self
    type(TrcFile), intent(in)                     ::  rcF
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Bound_Data_Init'
      
    ! --- local ----------------------------------
    
    integer               ::  itr
    integer               ::  ispec
    character(len=512)    ::  line
      
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("setup boundary conditions from data variables ...")'); call goPr
    
    ! maximum storage:
    allocate( self%ispec(n_advected), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%varsum(n_advected), stat=status )
    IF_NOTOK_RETURN(status=1)
    ! no data yet:
    self%ntracer = 0
    
    ! loop over advected tracers:
    do itr = 1, n_advected
      ! spec index:
      ispec = ispecs_advected(itr)
    
      ! read description line:
      call ReadRc( rcF, 'le.bound.data.'//trim(specname(ispec)), line, status )
      IF_NOTOK_RETURN(status=1)
      
      ! undefined? then skip this spec:
      if ( trim(line) == 'None' ) cycle
      
      !! info ...
      !write (gol,'("  fill ",a," from: ",a)') specname(ispec), trim(line); call goPr
      
      ! increase counter:
      self%ntracer = self%ntracer + 1
      ! store model index:
      self%ispec(self%ntracer) = ispec
      ! init weighted sum:
      call self%varsum(self%ntracer)%Init( line, status, verbose=.true. )
      IF_NOTOK_RETURN(status=1)
      ! enable all:
      call self%varsum(self%ntracer)%Enable( status )
      IF_NOTOK_RETURN(status=1)

    end do ! specs

    ! ok
    status = 0
    
  end subroutine Bound_Data_Init


  ! ***


  subroutine Bound_Data_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_LE_Bound_Data), intent(inout)         ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Bound_Data_Done'
      
    ! --- local ----------------------------------
    
    integer       ::  itracer
      
    ! --- begin ----------------------------------

    ! loop over tracers:
    do itracer = 1, self%ntracer
      ! done:
      call self%varsum(itracer)%Done( status )
      IF_NOTOK_RETURN(status=1)
    end do 
    
    ! maximum storage:
    deallocate( self%ispec, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%varsum, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Bound_Data_Done


  ! ***


  ! Copy values from data varaibles to boundary arrays

  subroutine Bound_Data_FillBounds( self, filled, status )
  
    use Binas          , only : xm_air
    use Dims           , only : nx, ny, nz
    use Indices        , only : nspec
    use LE_Bound_Common, only : bc_west, bc_east, bc_south, bc_north
    use LE_Bound_Common, only : caloft
    use LE_Data        , only : LE_Data_GetPointer
    use LE_Data_Common , only : nlev, nlev_top
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Bound_Data), intent(inout)         ::  self
    logical, intent(out)                          ::  filled(:)  ! (nspec)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Bound_Data_FillBounds'
  
    ! --- local ----------------------------------
    
    integer             ::  itracer
    integer             ::  ispec
    integer             ::  n
    integer             ::  i, k
    real                ::  w
    real, pointer       ::  pdata(:,:,:)
    character(len=64)   ::  units
    real                ::  factor
    logical             ::  density_factor
    real, pointer       ::  pdens(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(filled) /= nspec ) then
      write (gol,'("array filled has size (",i0,") while expected (",i0,")")') &
              size(filled), nspec; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no specs filled yet:
    filled = .false.
    
    ! loop over tracers:
    do itracer = 1, self%ntracer
    
      ! target spec:
      ispec = self%ispec(itracer)
            
      ! number of elements:
      call self%varsum(itracer)%Get( status, n=n )
      IF_NOTOK_RETURN(status=1)
      ! any?
      if ( n > 0 ) then
      
        ! init sums:
        bc_west (:  ,:,ispec) = 0.0
        bc_east (:  ,:,ispec) = 0.0
        bc_south(:  ,:,ispec) = 0.0
        bc_north(:  ,:,ispec) = 0.0
        caloft  (:,:,:,ispec) = 0.0

        ! loop over elements:
        do i = 1, n

          ! get pointer to data and weight:
          call self%varsum(itracer)%GetElement( i, status, pdata=pdata, w=w, units=units, &
                                        check_lbo=(/0,0,1/), check_ubo=(/nx+1,ny+1,nlev_top/) )
          IF_NOTOK_RETURN(status=1)

          ! conversion parameters:
          call Bound_Data_SpecConversions( ispec, trim(units),  &
                                             factor, density_factor, status )
          IF_NOTOK_RETURN(status=1)

          ! multiply with density?
          if ( density_factor ) then

            ! pointer to density:
            call LE_Data_GetPointer( 'dens', pdens, status, check_units='kg/m3', &
                                 check_lbo=(/0,0,1/), check_ubo=(/nx+1,ny+1,nz+1/) )
            IF_NOTOK_RETURN(status=1)
            ! copy from data variable data shaped (0:nx+1,0:nx+1,1:nlev_top):
            do k = 1, nlev
              bc_west (:,k,ispec) = bc_west (:,k,ispec) + w * pdata(   0,1:ny,k) * factor * pdens(   0,1:ny,k)
              bc_east (:,k,ispec) = bc_east (:,k,ispec) + w * pdata(nx+1,1:ny,k) * factor * pdens(nx+1,1:ny,k)
              bc_south(:,k,ispec) = bc_south(:,k,ispec) + w * pdata(1:nx,   0,k) * factor * pdens(1:nx,   0,k)
              bc_north(:,k,ispec) = bc_north(:,k,ispec) + w * pdata(1:nx,ny+1,k) * factor * pdens(1:nx,ny+1,k)
            end do
            do k = nlev+1, nlev_top
              caloft(:,:,k,ispec) = caloft(:,:,k,ispec) + w * pdata(1:nx,1:ny,k) * factor * pdens(1:nx,1:ny,k)
            end do

          else

            ! copy from data variable data shaped (0:nx+1,0:nx+1,1:nlev_top):
            do k = 1, nlev
              bc_west (:,k,ispec) = bc_west (:,k,ispec) + w * pdata(   0,1:ny,k) * factor
              bc_east (:,k,ispec) = bc_east (:,k,ispec) + w * pdata(nx+1,1:ny,k) * factor
              bc_south(:,k,ispec) = bc_south(:,k,ispec) + w * pdata(1:nx,   0,k) * factor
              bc_north(:,k,ispec) = bc_north(:,k,ispec) + w * pdata(1:nx,ny+1,k) * factor
            end do
            do k = nlev+1, nlev_top
              caloft(:,:,k,ispec) = caloft(:,:,k,ispec) + w * pdata(1:nx,1:ny,k) * factor
            end do

          end if

        end do ! elements

        ! set flag:
        filled(ispec) = .true.
        
      end if  ! any contribution
    
    end do  ! tracers
    
    ! ok
    status = 0
    
  end subroutine Bound_Data_FillBounds


  ! ***


  subroutine Bound_Data_FillConc( self, c, filled, status )
  
    use Dims   , only : nx, ny, nz
    use Indices, only : nspec
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------------
    
    class(T_LE_Bound_Data), intent(inout)         ::  self
    real, intent(inout)                           ::  c(:,:,:,:)  ! (nx,ny,nz,nspec)
    logical, intent(out)                          ::  filled(:)   ! (nspec)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Bound_Data_FillConc'
  
    ! --- local ----------------------------------
    
    integer             ::  itracer
    integer             ::  ispec
    integer             ::  n
    integer             ::  i
    real                ::  w
    real, pointer       ::  pdata(:,:,:)
    character(len=64)   ::  units
    real                ::  factor
    logical             ::  density_factor
    real, pointer       ::  pdens(:,:,:)
      
    ! --- begin ----------------------------------
    
    ! check ...
    if ( any(shape(c)/=(/nx,ny,nz,nspec/)) ) then
      write (gol,'("array c has shape (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
               shape(c), (/nx,ny,nz,nspec/); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( size(filled) /= nspec ) then
      write (gol,'("array filled has size (",i0,") while expected (",i0,")")') &
              size(filled), nspec; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! no specs filled yet:
    filled = .false.
    
    ! loop over tracers:
    do itracer = 1, self%ntracer
    
      ! target spec:
      ispec = self%ispec(itracer)
            
      ! init sum:
      c(:,:,:,ispec) = 0.0
            
      ! number of elements:
      call self%varsum(itracer)%Get( status, n=n )
      IF_NOTOK_RETURN(status=1)
      ! loop over elements:
      do i = 1, n
      
        ! get pointer to data and weight:
        call self%varsum(itracer)%GetElement( i, status, pdata=pdata, w=w, units=units, &
                                      check_lbo=(/1,1,1/), check_ubo=(/nx,ny,nz/) )
        IF_NOTOK_RETURN(status=1)
        
        ! conversion parameters:
        call Bound_Data_SpecConversions( ispec, trim(units),  &
                                           factor, density_factor, status )
        IF_NOTOK_RETURN(status=1)
        
        ! multiply with density?
        if ( density_factor ) then

          ! pointer to density:
          call LE_Data_GetPointer( 'dens', pdens, status, check_units='kg/m3', &
                               check_lbo=(/1,1,1/), check_ubo=(/nx,ny,nz/) )
          IF_NOTOK_RETURN(status=1)

          ! copy core from data variable:
          c(:,:,:,ispec) = c(:,:,:,ispec) + w * pdata(1:nx,1:ny,1:nz) * factor * pdens(1:nx,1:ny,1:nz)
        
        else

          ! copy core from data variable:
          c(:,:,:,ispec) = c(:,:,:,ispec) + w * pdata(1:nx,1:ny,1:nz) * factor
          
        end if
        
      end do ! elements
      
      ! set flag:
      filled(ispec) = .true.
    
    end do  ! tracers
    
    ! ok
    status = 0
    
  end subroutine Bound_Data_FillConc


  ! ***


  ! factor and flags for conversion from "units" to spec units
  
  subroutine Bound_Data_SpecConversions( ispec, units,  &
                                           factor, density_factor, status )
  
    use Binas  , only : xm_air
    use Indices, only : specunit, specmolm, specname
    use LE_Data, only : LE_Data_CompareUnits
    
    ! --- in/out ---------------------------------
    
    integer, intent(in)                           ::  ispec
    character(len=*), intent(in)                  ::  units
    real, intent(out)                             ::  factor
    logical, intent(out)                          ::  density_factor
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Bound_Data_SpecConversions'
  
    ! --- local ----------------------------------
    
    character(len=256)  ::  conversion
      
    ! --- begin ----------------------------------
    
    ! no conversion by default:
    factor = 1.0
    density_factor = .false.

    ! compare units:
    call LE_Data_CompareUnits( trim(units), trim(specunit(ispec)), status, verbose=.false. )
    IF_ERROR_RETURN(status=1)
    if ( status /= 0 ) then
      ! required conversion:
      conversion = trim(units)//' -> '//trim(specunit(ispec))
      ! switch:
      select case ( trim(conversion) )
        !~ volume mixing ratio:
        case ( 'mole/mole -> ppb', 'mol/mol -> ppb' )
          !         ppb
          factor = 1.0e9
        !~ mass concentration to mass concentration:
        case ( 'kg/m3 -> ug/m3' )
          !        
          factor = 1.0e9
        !~ mass mixing ratio to volume mixing ratio:
        case ( 'kg/kg -> ppb' )
          !  (kg tr)/(kg air) * (kg air)/(mole air) / ((kg tr)/(mole tr)) * 1e9 = ppb
          factor = xm_air / specmolm(ispec) * 1.0e9
        !~ mass mixing ratio to mass concentratin:
        case ( 'kg/kg -> ug/m3' )
          !  (kg tr)/(kg air) * (kg air)/(m3 air) * 1e9 = (ug tr)/(m3 air)
          factor = 1.0e9
          density_factor = .true.
        !~
        case default
          write (gol,'("unsupported conversion `",a,"` needed for spec `",a,"`")') &
                   trim(conversion), trim(specname(ispec)); call goErr
          TRACEBACK; status=1; return
      end select
    end if ! conversion needed

    ! ok
    status = 0
    
  end subroutine Bound_Data_SpecConversions


end module LE_Bound_Data


