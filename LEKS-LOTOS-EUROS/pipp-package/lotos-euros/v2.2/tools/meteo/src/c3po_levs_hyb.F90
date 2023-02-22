!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#include "c3po.inc"
!
!######################################################################

module C3PO_Levs_Hyb

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Levs_Hyb
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Levs_Hyb'
  

  ! --- types ----------------------------------------
  
  type T_Levs_Hyb
    integer                 ::  nlev                 ! number of levels
    integer                 ::  nv                   ! number of boundary values (2)
    real, allocatable       ::  ap(:)                ! (nlev) hybride coeff
    real, allocatable       ::   b(:)                ! (nlev) hybride coeff
    real, allocatable       ::  ap_bnds(:,:)         ! (nv,nlev) hybride coeff bounds
    real, allocatable       ::   b_bnds(:,:)         ! (nv,nlev) hybride coeff bounds
    real, allocatable       ::  ap_i(:)              ! (nlev) hybride coeff interface
    real, allocatable       ::   b_i(:)              ! (nlev) hybride coeff interface
  contains
    !
    procedure   ::                  Levs_Hyb_Init_ab_i
    procedure   ::                  Levs_Hyb_Init_ab_bnds
    procedure   ::                  Levs_Hyb_Init_coarsen
    generic     ::  Init         => Levs_Hyb_Init_ab_i, &
                                    Levs_Hyb_Init_ab_bnds, &
                                    Levs_Hyb_Init_coarsen
    !
    procedure   ::  Done         => Levs_Hyb_Done
    !
    procedure   ::  Ps_to_Ph     => Levs_Hyb_Ps_to_Ph_1d
    !
    procedure   ::                  Levs_Hyb_Ps_to_dP_r4_2d
    procedure   ::                  Levs_Hyb_Ps_to_dP_r4_3d
    procedure   ::                  Levs_Hyb_Ps_to_dP_r8_2d
    procedure   ::                  Levs_Hyb_Ps_to_dP_r8_3d
    generic     ::  Ps_to_dP     => Levs_Hyb_Ps_to_dP_r4_2d, &
                                    Levs_Hyb_Ps_to_dP_r4_3d, &
                                    Levs_Hyb_Ps_to_dP_r8_2d, &
                                    Levs_Hyb_Ps_to_dP_r8_3d
  end type T_Levs_Hyb

! adhoc ...
#ifdef without_f2003
#define XTYPE type
#else
#define XTYPE class
#endif
    
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** hybride levels
  ! ***
  ! ********************************************************************


  subroutine Levs_Hyb_Init_ab_i( self, ap, b, ap_i, b_i, status  )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(out)  ::  self
    real, intent(in)                ::  ap(:)      ! (  nlev)
    real, intent(in)                ::  b(:)       ! (  nlev)
    real, intent(in)                ::  ap_i(:)    ! (0:nlev)
    real, intent(in)                ::  b_i(:)     ! (0:nlev)
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Init_ab_bnds'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! size:
    self%nlev = size(ap)
    ! fixed:
    self%nv = 2
    
    ! check ...
    if ( size(b) /= self%nlev ) then
      write (gol,'("size of argument b is ",i6," not ",i6)') size(b), self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( size(ap_i) /= self%nlev+1 ) then
      write (gol,'("shape of argument ap_i is ",i6," not ",i6)') size(ap_i), self%nlev+1; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( size(b_i) /=  self%nlev+1 ) then
      write (gol,'("shape of argument b_i is ",i6," not ",i6)') size(b_i), self%nlev+1; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( self%ap(self%nlev) )
    allocate( self%b (self%nlev) )
    allocate( self%ap_bnds(self%nv,self%nlev) )
    allocate( self%b_bnds (self%nv,self%nlev) )
    allocate( self%ap_i(0:self%nlev) )
    allocate( self%b_i (0:self%nlev) )
    
    ! store full level values:
    self%ap   = ap
    self%b    = b
    ! store interface values:
    self%ap_i = ap_i
    self%b_i  = b_i
    ! bounds:
    self%ap_bnds(1,:) = self%ap_i(0:self%nlev-1)
    self%ap_bnds(2,:) = self%ap_i(1:self%nlev  )
    self%b_bnds (1,:) = self%b_i (0:self%nlev-1)
    self%b_bnds (2,:) = self%b_i (1:self%nlev  )

    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Init_ab_i


  ! ***


  subroutine Levs_Hyb_Init_ab_bnds( self, ap, b, ap_bnds, b_bnds, status  )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(out)  ::  self
    real, intent(in)                ::  ap(:)           ! (   nlev)
    real, intent(in)                ::  b(:)            ! (   nlev)
    real, intent(in)                ::  ap_bnds(:,:)    ! (nv,nlev)
    real, intent(in)                ::  b_bnds(:,:)     ! (nv,nlev)
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Init_ab_bnds'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! size:
    self%nlev = size(ap)
    self%nv   = size(ap_bnds,1)
    
    ! check ...
    if ( size(b) /= self%nlev ) then
      write (gol,'("size of argument b is ",i6," not ",i6)') size(b), self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( any( shape(ap_bnds) /=  (/self%nv,self%nlev/) ) ) then
      write (gol,'("shape of argument ap_bnds is ",2i6," not ",2i6)') shape(ap_bnds), self%nv,self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( any( shape(b_bnds) /=  (/self%nv,self%nlev/) ) ) then
      write (gol,'("shape of argument b_bnds is ",2i6," not ",2i6)') shape(b_bnds), self%nv,self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( self%ap(self%nlev) )
    allocate( self%b (self%nlev) )
    allocate( self%ap_bnds(self%nv,self%nlev) )
    allocate( self%b_bnds (self%nv,self%nlev) )
    allocate( self%ap_i(0:self%nlev) )
    allocate( self%b_i (0:self%nlev) )
    
    ! store:
    self%ap = ap
    self%b  = b
    self%ap_bnds = ap_bnds
    self%b_bnds  = b_bnds
    ! fill interface values:
    self%ap_i(0 ) = ap_bnds(1,1)
    self%ap_i(1:) = ap_bnds(2,:)
    self%b_i (0 ) =  b_bnds(1,1)
    self%b_i (1:) =  b_bnds(2,:)

    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Init_ab_bnds


  ! ***


  subroutine Levs_Hyb_Init_coarsen( self, hyb, nn, status, mapping  )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(out)  ::  self
    XTYPE(T_Levs_Hyb), intent(in)   ::  hyb
    integer, intent(in)             ::  nn(:)
    integer, intent(out)            ::  status
    
    integer, intent(out), allocatable, optional   ::  mapping(:,:)  ! (nlev,2)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Init_coarsen'
    
    ! --- local ----------------------------------
    
    integer             ::  nlev, ilev
    integer             ::  k0, k1, k2
    integer             ::  nv
    real, allocatable   ::  ap(:)           ! (   nlev)
    real, allocatable   ::  b(:)            ! (   nlev)
    real, allocatable   ::  ap_bnds(:,:)    ! (nv,nlev)
    real, allocatable   ::  b_bnds(:,:)     ! (nv,nlev)
    
    ! --- begin ----------------------------------
    
    ! expect top-down order ...
    if ( hyb%b(1) > hyb%b(hyb%nlev) ) then
      write (gol,'("expecting top-down input levels, found b range ",f0.0," - ",f0.0)') &
               hyb%b(1), hyb%b(hyb%nlev); call goErr
      TRACEBACK; status=1; return
    end if
    
    !! target number of levels as integer divison (floor):
    !nlev = hyb%nlev / n
    !! offset:
    !k0 = hyb%nlev - nlev * n
    ! target number of levels:
    nlev = size(nn)
    ! offset:
    k0 = hyb%nlev - sum(nn)
    ! copy:
    nv = hyb%nv
    
    ! storage:
    if ( present(mapping) ) then
      if ( allocated(mapping) ) deallocate( mapping )
      allocate( mapping(nlev,2) )
    end if
    
    ! storage:
    allocate( ap(nlev) )
    allocate( b (nlev) )
    allocate( ap_bnds(nv,nlev) )
    allocate( b_bnds (nv,nlev) )
    
    !! show input:
    !print *, 'input level bounds:'
    !do ilev = 1, hyb%nlev
    !  print *, '  ilev ', ilev
    !  print *, '    ap_bnds = ', hyb%ap_bnds(:,ilev)
    !  print *, '    b_bnds  = ', hyb%b_bnds(:,ilev)
    !end do
    !! show output:
    !print *, 'output level bounds:'
    
    ! loop over target layers:
    do ilev = 1, nlev
      ! range of input layers,
      ! take care of reversed order:
      if ( ilev == 1 ) then
        k1 = k0 + 1
      else
        k1 = k0 + sum(nn(1:ilev-1)) + 1
      end if
      k2 = k0 + sum(nn(1:ilev))
      ! store ?
      if ( present(mapping) ) then
        mapping(ilev,1) = k1
        mapping(ilev,2) = k2
      end if
      ! copy bounds:
      ap_bnds(1,ilev) = hyb%ap_bnds(1,k1)
       b_bnds(1,ilev) =  hyb%b_bnds(1,k1)
      ap_bnds(2,ilev) = hyb%ap_bnds(2,k2)
       b_bnds(2,ilev) =  hyb%b_bnds(2,k2)
      ! mid values:
      ap(ilev) = 0.5 * ( ap_bnds(1,ilev) + ap_bnds(2,ilev) )
       b(ilev) = 0.5 * (  b_bnds(1,ilev) +  b_bnds(2,ilev) )
      !! info ...
      !print *, '  ilev ', ilev, k1, k2
      !print *, '    ap_bnds = ', ap_bnds(:,ilev)
      !print *, '    b_bnds  = ', b_bnds(:,ilev)
    end do
    
    !! show input:
    !print *, 'output levels:'
    !do ilev = 1, nlev
    !  print *, ilev, ap(ilev), b(ilev)
    !end do
    
    ! init:
    call self%Init( ap, b, ap_bnds, b_bnds, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( ap )
    deallocate( b  )
    deallocate( ap_bnds )
    deallocate( b_bnds  )

    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Init_coarsen


  ! ***
  
  
  subroutine Levs_Hyb_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(inout)    ::  self
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%ap )
    deallocate( self%b  )
    deallocate( self%ap_bnds )
    deallocate( self%b_bnds  )
    deallocate( self%ap_i )
    deallocate( self%b_i  )
    
    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Done
  
  
  ! ***
  
  
  subroutine Levs_Hyb_Ps_to_Ph_1d( self, ps, ph, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(in)       ::  self
    real, intent(in)                    ::  ps(:)    ! (np) surface pressure [Pa]
    real, intent(out)                   ::  ph(:,:)  ! (np,0:nlev) half.lev.p. [Pa]
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Ps_to_Ph_1d'
    
    ! --- local ----------------------------------
    
    integer      ::  i
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( any( shape(ph) /= (/size(ps),self%nlev+1/) ) ) then
      write (gol,'("unexpected shape of ph : ",2i6)') shape(ph); call goErr
      write (gol,'("  grid size   : ",i6)') size(ps); call goErr
      write (gol,'("  half levels : ",i6," + 1")') self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop:
    do i = 1, size(ps)
      ! apply interface coeffs:
      ph(i,:) = self%ap_i + self%b_i * ps(i)  ! [Pa]
    end do
    
    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Ps_to_Ph_1d
  
  
  ! ***
  
  
  subroutine Levs_Hyb_Ps_to_dP_r4_2d( self, ps, dp, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(in)       ::  self
    real(4), intent(in)                 ::  ps(:  ,:) ! (nx     ,nt) surface pressure [Pa]
    real(4), intent(out)                ::  dp(:,:,:) ! (nx,nlev,nt) delta pressure [Pa]
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Ps_to_dP_r4_2d'
    
    ! --- local ----------------------------------
    
    integer      ::  ilev
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( any( shape(dp) /= (/size(ps,1),self%nlev,size(ps,2)/) ) ) then
      write (gol,'("unexpected shape of dp : ",3i6)') shape(dp); call goErr
      write (gol,'("  sp shape   : ",2i6)') shape(ps); call goErr
      write (gol,'("  levels     : ",i6)') self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop:
    do ilev = 1, self%nlev
      ! apply interface coeffs:
      dp(:,ilev,:) = abs( (self%ap_bnds(2,ilev)-self%ap_bnds(1,ilev)) + &
                          ( self%b_bnds(2,ilev)- self%b_bnds(1,ilev))*ps )
    end do
    
    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Ps_to_dP_r4_2d
  
  
  ! *
  
  
  subroutine Levs_Hyb_Ps_to_dP_r4_3d( self, ps, dp, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(in)       ::  self
    real(4), intent(in)                 ::  ps(:,:  ,:) ! (nx,ny     ,nt) surface pressure [Pa]
    real(4), intent(out)                ::  dp(:,:,:,:) ! (nx,ny,nlev,nt) delta pressure [Pa]
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Ps_to_dP_r4_3d'
    
    ! --- local ----------------------------------
    
    integer      ::  ilev
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( any( shape(dp) /= (/size(ps,1),size(ps,2),self%nlev,size(ps,3)/) ) ) then
      write (gol,'("unexpected shape of dp : ",4i6)') shape(dp); call goErr
      write (gol,'("  sp shape   : ",3i6)') shape(ps); call goErr
      write (gol,'("  levels     : ",i6)') self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop:
    do ilev = 1, self%nlev
      ! apply interface coeffs:
      dp(:,:,ilev,:) = abs( (self%ap_bnds(2,ilev)-self%ap_bnds(1,ilev)) + &
                            ( self%b_bnds(2,ilev)- self%b_bnds(1,ilev))*ps )
    end do
    
    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Ps_to_dP_r4_3d
  
  
  ! ***
  
  
  subroutine Levs_Hyb_Ps_to_dP_r8_2d( self, ps, dp, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(in)       ::  self
    real(8), intent(in)                 ::  ps(:  ,:) ! (nx     ,nt) surface pressure [Pa]
                                                      ! (nx,ny     )
    real(8), intent(out)                ::  dp(:,:,:) ! (nx,nlev,nt) delta pressure [Pa]
                                                      ! (nx,ny,nlev)
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Ps_to_dP_r8_2d'
    
    ! --- local ----------------------------------
    
    integer      ::  ilev
    
    ! --- begin ----------------------------------
    
    ! switch ...
    !~ (nx,nlev,nt)
    if ( all( shape(dp) == (/size(ps,1),self%nlev,size(ps,2)/) ) ) then
    
      ! loop:
      do ilev = 1, self%nlev
        ! apply interface coeffs:
        dp(:,ilev,:) = abs( (self%ap_bnds(2,ilev)-self%ap_bnds(1,ilev)) + &
                            ( self%b_bnds(2,ilev)- self%b_bnds(1,ilev))*ps )
      end do

    !~ (nx,ny,nlev)
    else if ( all( shape(dp) == (/size(ps,1),size(ps,2),self%nlev/) ) ) then
    
      ! loop:
      do ilev = 1, self%nlev
        ! apply interface coeffs:
        dp(:,:,ilev) = abs( (self%ap_bnds(2,ilev)-self%ap_bnds(1,ilev)) + &
                            ( self%b_bnds(2,ilev)- self%b_bnds(1,ilev))*ps )
      end do
    
    else
    
      write (gol,'("unsupported shape of dp : ",3i6)') shape(dp); call goErr
      write (gol,'("  sp shape   : ",2i6)') shape(ps); call goErr
      write (gol,'("  levels     : ",i6)') self%nlev; call goErr
      TRACEBACK; status=1; return

    end if
    
    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Ps_to_dP_r8_2d
  
  
  ! *
  
  
  subroutine Levs_Hyb_Ps_to_dP_r8_3d( self, ps, dp, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(T_Levs_Hyb), intent(in)       ::  self
    real(8), intent(in)                 ::  ps(:,:  ,:) ! (nx,ny     ,nt) surface pressure [Pa]
    real(8), intent(out)                ::  dp(:,:,:,:) ! (nx,ny,nlev,nt) delta pressure [Pa]
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Levs_Hyb_Ps_to_dP_r8_3d'
    
    ! --- local ----------------------------------
    
    integer      ::  ilev
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( any( shape(dp) /= (/size(ps,1),size(ps,2),self%nlev,size(ps,3)/) ) ) then
      write (gol,'("unexpected shape of dp : ",4i6)') shape(dp); call goErr
      write (gol,'("  sp shape   : ",3i6)') shape(ps); call goErr
      write (gol,'("  levels     : ",i6)') self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop:
    do ilev = 1, self%nlev
      ! apply interface coeffs:
      dp(:,:,ilev,:) = abs( (self%ap_bnds(2,ilev)-self%ap_bnds(1,ilev)) + &
                            ( self%b_bnds(2,ilev)- self%b_bnds(1,ilev))*ps )
    end do
    
    ! ok
    status = 0
  
  end subroutine Levs_Hyb_Ps_to_dP_r8_3d
  
  
end module C3PO_Levs_Hyb
