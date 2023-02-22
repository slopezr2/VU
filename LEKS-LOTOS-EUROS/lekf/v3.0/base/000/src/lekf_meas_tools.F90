!###############################################################################
!
! Tools for measurement update.
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
!###############################################################################


module LEKF_Meas_Tools

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------

  private

  public    ::  GetSpecApply
  
  public    ::  TGainBoundInfo
  

  ! --- const --------------------------------

  character(len=*), parameter   ::  mname = 'LEKF_Meas_Tools'
    

  ! --- types --------------------------------

  ! gain bounds
  
  type TGainBoundInfo
    !  o list off all cells correlated with a cell (ix,iy)
    integer, allocatable    ::  nc_corr(:,:)
    integer, allocatable    ::  ix_corr(:,:,:)
    integer, allocatable    ::  iy_corr(:,:,:)
    real, allocatable       ::  fc_corr(:,:,:)
  contains
    procedure ::                      gainbound_Init_0d
    procedure ::                      gainbound_Init_2d
    generic   ::  Init            =>  gainbound_Init_0d, &
                                      gainbound_Init_2d
    procedure ::  Done            =>  gainbound_Done
  end type TGainBoundInfo
  

contains


  ! =========================================================================
  ! ===
  ! === tracer selection
  ! ===
  ! =========================================================================

  
  ! Given a keyword, determine analysis weights for species.
  ! This is used to analyze selected species only for certain observations.
  ! Specify one of the keywords: 'gasses','aerosols', 'both'
  ! Return value is logical array of length nspec which is .true. if a specie 
  ! is analyzed and .false. if not.
  
  subroutine GetSpecApply( line, apply, apply_to_sia, status )
  
    use Indices, only : nspec_all, specname
    use Indices, only : specmode, NO_AEROSOL_MODE
    use Indices, only : n_sia, ispecs_sia

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)     ::  line
    logical, intent(out)             ::  apply(:)   ! (nspec_all)
    logical, intent(out)             ::  apply_to_sia
    integer, intent(out)             ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/GetSpecapply'

    ! --- local ----------------------------------
    
    integer                        ::  ispec
    integer                        ::  k
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( size(apply) /= nspec_all ) then
      write (gol,'("output apply should have length ",i6," not ",i6)') nspec_all, size(apply); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! switch:
    select case ( trim(line) )
      case ( 'gasses'   ) ; apply = specmode == NO_AEROSOL_MODE
      case ( 'aerosols' ) ; apply = specmode /= NO_AEROSOL_MODE
      case ( 'both'     ) ; apply = .true.
      case default
        write (gol,'("unsupported spec apply keyword `",a,"`")') trim(line); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! apply to sia? used to enable update of aerh2o:
    apply_to_sia = .false.
    ! loop over sia tracers:
    do k = 1, n_sia
      ! global index:
      ispec = ispecs_sia(k)
      ! enabled?
      if ( apply(ispec) ) then
        apply_to_sia = .true.
        exit
      end if
    end do  ! sia specs

    ! info ...
    write (gol,'("LEKF:       spec weigts:")'); call goPr
    write (gol,'("LEKF:         line  : ",a)') trim(line); call goPr
    write (gol,'("LEKF:         number name (apply)")'); call goPr
    do ispec = 1, nspec_all
      write (gol,'("LEKF:         ",i2," ",a," (",l1,")")') ispec, trim(specname(ispec)), apply(ispec); call goPr
    end do
    write (gol,'("LEKF:         apply to sia : ",l1)') apply_to_sia; call goPr
    
    ! ok
    status = 0
    
  end subroutine GetSpecApply
  

  ! =========================================================================
  ! ===
  ! === gainbound
  ! ===
  ! =========================================================================

  
  subroutine gainbound_Init_0d( self, rho, status )
  
    use dims, only : nx, ny
    
    ! --- in/out ----------------------------
    
    class(TGainBoundInfo), intent(out)    ::  self
    real, intent(in)                      ::  rho
    integer, intent(out)                  ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'gainbound_Init_2d'
    
    ! --- local -----------------------------
    
    real, allocatable   ::  rho_2d(:,:)
    
    ! --- begin -----------------------------
    
    ! storage:
    allocate( rho_2d(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! fill
    rho_2d = rho
    
    ! init:
    call self%Init( rho_2d, status )
    IF_NOTOK_RETURN(status=1)
    
    ! clear:
    deallocate( rho_2d, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine gainbound_Init_0d
  
  ! *  
  
  subroutine gainbound_Init_2d( self, rho, status )
  
    use dims, only : nx, ny
    use dims, only : runF
    
    ! --- in/out ----------------------------
    
    class(TGainBoundInfo), intent(out)    ::  self
    real, intent(in)                      ::  rho(nx,ny)
    integer, intent(out)                  ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'gainbound_Init_2d'
    
    ! --- local -----------------------------
    
    integer                 ::  ix, iy
    integer                 ::  i, j
    integer                 ::  mx_corr
    integer                 ::  icell
    real                    ::  dist
    
    ! --- begin -----------------------------
    
    ! setup counter array:
    allocate( self%nc_corr(nx,ny), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! no cells counted yet:    
    self%nc_corr = 0

    ! loop over all cells:
    do iy = 1, ny
      do ix = 1, nx
        ! loop over all possibly correlated cells:
        do j = 1, ny
          do i = 1, nx
            ! distance in km:
            dist = sqrt( ((i-ix)*runF%dx(iy))**2 + ((j-iy)*runF%dy)**2 )
            ! within range ?
            if ( dist <= 3.5*rho(ix,iy) ) self%nc_corr(ix,iy) = self%nc_corr(ix,iy) + 1
          end do
        end do
        ! end for this cell
      end do
    end do
    
    ! maximum number of correlated cells:
    mx_corr = maxval(self%nc_corr)

    ! indices of correlated cells and weight for all cells:
    allocate( self%ix_corr(nx,ny,mx_corr), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%iy_corr(nx,ny,mx_corr), stat=status )
    IF_NOTOK_RETURN(status=1)
    allocate( self%fc_corr(nx,ny,mx_corr), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over all cells:
    do iy = 1, ny
      do ix = 1, nx
        ! no cells filled yet:
        icell = 0
        ! loop over all possibly correlated cells:
        do j = 1, ny
          do i = 1, nx
            ! distance in km:
            dist = sqrt( ((i-ix)*runF%dx(iy))**2 + ((j-iy)*runF%dy)**2 )
            ! within range ?
            if ( dist <= 3.5*rho(ix,iy) ) then
              ! next index:
              icell = icell + 1
              ! check ...
              if ( icell > mx_corr ) stop 'BUG - icell > mx_corr'
              ! store location and weight:
              self%ix_corr(ix,iy,icell) = i
              self%iy_corr(ix,iy,icell) = j
              ! compute correlation factor:
              if ( rho(ix,iy) == 0.0 ) then
                self%fc_corr(ix,iy,icell) = 1.0
              else
                self%fc_corr(ix,iy,icell) = exp( -0.5*(dist/rho(ix,iy))**2 )
              end if
            end if
          end do
        end do
        ! end for this cell
      end do
    end do
    
    ! ok
    status = 0
  
  end subroutine gainbound_Init_2d


  ! ***
  
  
  subroutine gainbound_Done( self, status )
  
    ! --- in/out ----------------------------
    
    class(TGainBoundInfo), intent(inout)    ::  self
    integer, intent(out)                    ::  status
    
    ! --- const ------------------------------------
    
    character(len=*), parameter  ::  rname = 'gainbound_Done'
    
    ! --- local -----------------------------
    
    ! --- begin -----------------------------
    
    ! clear:
    deallocate( self%nc_corr, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%ix_corr, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%iy_corr, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%fc_corr, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine gainbound_Done
  


end module LEKF_Meas_Tools
