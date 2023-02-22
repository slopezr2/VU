!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_VDiff

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out ----------------------------

  private
  
  public   ::  LE_VDiff_Init, LE_VDiff_Done
  public   ::  LE_VDiff_Apply


  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_VDiff'

  ! local
  character(len=16)    ::  kz_type 
  
contains


  ! ====================================================================
  
  
  subroutine LE_VDiff_Init( rcF, status )

    use LE_Data, only : LE_Data_Enable
    use GO, only : ReadRc, TRcFile

    ! --- in/out ---------------------------------
    
    type(TRcFile), intent(in)        ::  rcF
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ LE_VDiff_Init'

    ! --- local ----------------------------------

    character(len=16)   ::  key
    
    ! --- begin ----------------------------------
    
    ! enable data:
    call LE_Data_Enable( 'h',status )        
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'vol',status )        
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'airm',status )        
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'kz',status )        
    IF_NOTOK_RETURN(status=1)
    
    ! Read Kz type
    call ReadRc( rcF, 'vdiff.kz_type', kz_type, status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_VDiff_Init


  ! ***
  
  
  subroutine LE_VDiff_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ LE_VDiff_Done'
    
    ! --- local ----------------------------------
    
    integer       ::  ibud

    ! --- begin ----------------------------------
    
    ! ok
    status = 0
    
  end subroutine LE_VDiff_Done


  ! ***  


  subroutine LE_VDiff_Apply( c, dt_min, status )

    use LE_Logging , only : ident2
    use Dims       , only : nx, ny, nz
    use LE_Bound   , only : caloft         ! [ppb] or [ug/m3]
    use Indices    , only : nspec
    use LE_Data    , only : LE_Data_GetPointer

    ! --- in/out -------------------------------------

    real, intent(inout)   ::  c(nx,ny,nz,nspec)
    real, intent(in)      ::  dt_min  ! minutes
    integer, intent(out)  ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_VDiff_Apply'
 
    ! --- local --------------------------------------

    real, pointer        ::  h_m(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::  vol(:,:,:)   ! (lon,lat,lev)
    real, pointer        ::  airm(:,:,:)   ! (lon,lat,lev)
    integer              ::  n
    real                 ::  dtvd
    integer              ::  status_par
    integer              ::  i,j
    real, pointer        ::  kz(:,:,:)   ! (lon,lat,lev)          

    ! --- begin ------------------------------------------

    !! info ...
    !write (gol,'(a,"<vertical diffusion>")') ident2; call goPr

    ! pointers to meteo data:
    call LE_Data_GetPointer( 'h', h_m, status, check_units ='m' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'vol', vol, status, check_units ='m3' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'airm', airm, status, check_units ='kg' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'kz', kz, status, check_units = 'm2/s')    
    IF_NOTOK_RETURN(status=1)    

    ! set time step in seconds:
    dtvd = dt_min * 60.0  ! sec

    ! status for parallel computing
    status_par = 0

    !$OMP parallel &
    !$OMP   default( none ) &
    !$OMP   shared( nx, ny, nz ) &
    !$OMP   shared( h_m, Kz ) &
    !$OMP   shared( vol ) &
    !$OMP   shared( airm ) &
    !$OMP   shared( caloft ) &
    !$OMP   shared( c ) &
    !$OMP   shared( dtvd ) &
    !$OMP   private( i, j ) &
    !$OMP   private( status ) &
    !$OMP   reduction( + : status_par )
    !$OMP   do SCHEDULE( DYNAMIC )
    do j = 1, ny
      do i = 1, nx

        ! apply 1D diffusion on all specs:
        call LE_VDiff_Apply_1d_all( i, j, dtvd, h_m(i,j,:), Kz(i,j,1:nz), &
                                     vol(i,j,:), airm(i,j,:), &
                                      caloft(i,j,nz+1,:), c(i,j,:,:), status )
        status_par = status_par + abs(status)

      end do  ! i
    end do ! j
    !$OMP   end do
    !$OMP end parallel

    ! check ...
    if (status_par /= 0 ) then
      write (gol, '("Non-zero return status from OMP loop")'); call goErr
      TRACEBACK; status=1; return
    endif

    ! ok
    status = 0

  end subroutine LE_VDiff_Apply


  ! *

  subroutine LE_VDiff_Apply_1d_all( iloc, jloc, &
                                    dtvd, h_m, Kz, vol, airm, &
                                    caloft, c, status )

    use Binas     , only : xm_air
    use Dims      , only : nz
    use Indices   , only : nspec, specunit, specmolm
    use Indices    , only : i_o3

    ! --- in/out -------------------------------------

    integer, intent(in)   ::  iloc, jloc  ! needed for labeling
    real, intent(in)      ::  dtvd  ! sec
    real, intent(in)      ::  h_m(nz)   ! m
    real, intent(in)      ::  Kz(nz)     ! m2/s
    real, intent(in)      ::  vol(nz)    ! m3
    real, intent(in)      ::  airm(nz)   ! kg
    real, intent(in)      ::  caloft(nspec)  ! [ppb] or [ug/m3]
    real, intent(inout)   ::  c(nz,nspec)    ! [ppb] or [ug/m3]
    integer, intent(out)  ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_VDiff_Apply_1d_all'

    ! --- local --------------------------------------

    integer     ::  l
    real        ::  dh(0:nz)
    real        ::  fac
    real        ::  hh(nz)
    real        ::  diffmat(nz,nz)
    integer     ::  ispec
    real        ::  c_conc(nz), ca_conc
    real        ::  mass_before
    real        ::  mass_after

    ! --- begin ------------------------------------------


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( trim(kz_type) == 'msp' ) then
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! compute the mesh widths between two cell centers
      dh(1) = 0.5 * h_m(2)
      do l = 2, nz-1
        dh(l) = 0.5*( h_m(l+1) - h_m(l-1) )
      end do
      ! copy edges:
      dh(0)  = dh(1)
      dh(nz) = dh(nz-1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else if ( trim(kz_type) == 'normal' ) then
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! set distane between ghost points around interface
      ! from which concentrations are sampled
      ! to 2 * minimum distance of the cell interface to the cell centers:
      !~ first interface:
      dh(1) = 2* min( 0.5*h_m(1) , 0.5*(h_m(2)-h_m(1)) )
      !~ higher interfaces:
      do l = 2, nz-1
        dh(l) = 2 * min( 0.5*(h_m(l) - h_m(l-1)), 0.5*(h_m(l+1)-h_m(l)) )
      end do
      ! copy edges:
      dh(0)  = dh(1)
      dh(nz) = dh(nz-1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
      write( gol, '("Unknown kz type defined: ", a)' ) trim(kz_type) ; call GoErr
      TRACEBACK;status=1;return
    end if
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! compute the depths of the cells
    hh(1) = h_m(1)
    do l = 2, nz
       hh(l) = ( h_m(l) - h_m(l-1) )
    end do

    ! construct diffusion matrix for this grid column
    ! initialize diffusion matrix
    diffmat = 0.0
    do l=1,nz
      ! matrix is built up per flux
      fac = Kz(l)/dh(l)

      ! set weights for cell l
      diffmat(l,l) = diffmat(l,l) - fac/hh(l)
      if (l<nz) diffmat(l,l+1) = diffmat(l,l+1) + fac/hh(l)

      ! set weights for cell l+1
      if (l < nz) then
        diffmat(l+1,l+1)= diffmat(l+1,l+1) - fac/hh(l+1)
        diffmat(l+1,l  )= diffmat(l+1,l  ) + fac/hh(l+1)
      endif
    enddo

    ! make diffmat equal to I-dt*diffmat
    diffmat=-dtvd*diffmat
    do l=1,nz
      diffmat(l,l) = 1.0 + diffmat(l,l)
    enddo

    ! do for all species
    do ispec = 1, nspec
    
      ! check units, convert to (mass) concentration (per volume)
      ! since adjust is based on volume change
      select case ( trim(specunit(ispec)) )
        ! mole (volume) mixing ratio?
        case ( 'ppb' )
          ! conversion to mass concentation in ug/m3:
          !         1.0e-9             ! ((mole c)/(mole air)) / ppb
!          !         1.0e9              ! ug/kg
!          c_conc = c(:,ispec)                        &  ! ppb (mole c)/(mole air)/ppb ug/kg
!                    * specmolm(ispec)/xm_air         &  ! (kg c)/(mole c) (mole air)/(kg air)
!                    * airm / vol                        ! (kg air) / m3
!          ! use top layer properties for aloft:
!          ca_conc = caloft(ispec)                    &  ! (mole c)/(mole air) ug/kg
!                    * specmolm(ispec)/xm_air         &  ! (kg c)/(mole c) (mole air)/(kg air)
!                    * airm(nz) / vol(nz)                ! (kg air) / m3
          c_conc  = c(:,ispec)     
          ca_conc = caloft(ispec)  
          mass_before = sum( 1e-9*c_conc * specmolm(ispec)/xm_air * airm    ) + &
                             1e-9*ca_conc * specmolm(ispec)/xm_air * airm(nz)
        !~ mass concentration
        case ( 'ug/m3', 'grns/m3', 'mlc/cm3', '1/cm3' )
          c_conc  = c(:,ispec)     ! ug/m3
          ca_conc = caloft(ispec)  ! ug/m3
          mass_before = sum( 1e-9*c_conc * vol(:)    ) + &
                             1e-9*ca_conc * vol(nz)
          
        !~
        case default
          write (gol,'("volume concentrations not implemented for units `",a,"`")') &
                            trim(specunit(ispec)); call goErr
          TRACEBACK; status=1; return
      end select
    
      ! apply for this spec:
      call LE_VDiff_Apply_1d( iloc, jloc, ispec, &
                                Kz, hh, dh, dtvd, diffmat, &
                                 ca_conc, c_conc, status )
      IF_NOTOK_RETURN(status=1)
    
      ! convert back using volume and airmass,
      ! no need to reset the aloft:
      select case ( trim(specunit(ispec)) )
        ! mole (volume) mixing ratio?
        case ( 'ppb' )
          ! conversion from mass concentation in ug/m3:
          !         1.0e9              ! ppb / ((mole c)/(mole air))
          !         1.0e-9             ! kg/ug
!          c(:,ispec) = c_conc                        &  ! ug/m3 kg/ug ppb/((mole air)/(mole c))
!                    / specmolm(ispec)*xm_air         &  ! (mole c)/(kg c) (kg air)/(mole air)
!                    / airm * vol                        ! m3 / (kg air)
          mass_after = sum( 1e-9*c_conc * specmolm(ispec)/xm_air * airm    ) + &
                            1e-9*ca_conc * specmolm(ispec)/xm_air * airm(nz)
          
          ! scale with mass difference
          if ( mass_after > 0 ) then
            c(:,ispec) = mass_before/mass_after * c_conc  
          else 
            c(:,ispec) = c_conc
          end if
          
        !~ mass concentration
        case ( 'ug/m3', 'grns/m3', 'mlc/cm3', '1/cm3' )
          mass_after = sum( 1e-9*c_conc * vol(:)    ) + &
                            1e-9*ca_conc * vol(nz)
          if ( mass_after > 0 ) then
            c(:,ispec) = mass_before/mass_after * c_conc    ! ug/m3
          else 
            c(:,ispec) = c_conc
          end if
        !~
        case default
          write (gol,'("volume concentrations not implemented for units `",a,"`")') &
                            trim(specunit(ispec)); call goErr
          TRACEBACK; status=1; return
      end select
                                  
    end do ! spec


    ! ok
    status = 0

  end subroutine LE_VDiff_Apply_1d_all


  ! *


  subroutine LE_VDiff_Apply_1d( iloc, jloc, ispec, &
                                  Kz, hh, dh, dtvd, &
                                  diffmat, caloft, c, status )

    use Dims       , only : nz

    ! --- in/out -------------------------------------

    integer, intent(in)   ::  iloc, jloc, ispec  ! needed for labeling
    real, intent(in)      ::  dtvd  ! sec
    real, intent(in)      ::  diffmat(nz,nz)
    real, intent(in)      ::  Kz(nz)     ! m2/s
    real, intent(in)      ::  hh(nz)     ! m
    real, intent(in)      ::  dh(nz)     ! m
    real, intent(in)      ::  caloft     ! ug/m3
    real, intent(inout)   ::  c(nz)      ! ug/m3
    integer, intent(out)  ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_VDiff_Apply_1d'

    ! --- local --------------------------------------

    real      ::  cnew(nz)
    real      ::  rhs(nz)
    real      ::  diffmat_work(nz,nz)

    ! --- begin ------------------------------------------

    ! initial condition
    cnew(:)   = c(:)
    
    ! set right hand side
    rhs(:) = cnew(:)
    ! apply upper boundary condition
    rhs(nz) = rhs(nz) + Kz(nz)/dh(nz)/hh(nz)*dtvd*caloft

    ! add the deposition velocity to diffmat
    ! diffmat(1,1) = diffmat(1,1) + vd(i,j,k)*dtvd/hh(1)

    ! copy matrix, will be changed:
    diffmat_work = diffmat

    ! solve the tridiagonal system: Diffmat * cnew = rhs
    call TriDSolve( nz, diffmat_work, rhs, cnew, status )
    IF_NOTOK_RETURN(status=1)


    ! copy final solution for current species into "global" conc. array
    c(:) = cnew(:)

    ! ok
    status = 0

  end subroutine LE_VDiff_Apply_1d


  ! ***


  !
  ! Solve tri-diagonal system A*c = b, with A stored as a full matrix
  !

  subroutine TriDSolve( n, A, b, c, status )

    ! --- in/out -------------------------------------

    integer, intent(in)   ::  n
    real, intent(inout)   ::  A(n,n)
    real, intent(in)      ::  b(n)
    real, intent(out)     ::  c(n)
    integer, intent(out)  ::  status

    ! --- const --------------------------------

    character(len=*), parameter   :: rname = mname//'/LE_VDiff_Apply'

    ! --- local --------------------------------------

    integer         ::  l
    real            ::  rhs(n)
    real            ::  pivot

    ! --- begin ------------------------------------------

    ! copy the vector b to the right hand side
    rhs(:) = b(:)

    ! transform into bi-diagonal form
    do l=2,n
      ! check ..
      if (A(l-1,l-1) == 0.0) then
        write (gol,'("singular diffusion matrix at row ",i0)') l; call goErr
        TRACEBACK; status=1; return
      end if
      pivot  = A(l,l-1)/A(l-1,l-1)
      A(l,l) = A(l,l)   - pivot*A(l-1,l)
      if (l /= n) A(l,l+1) = A(l,l+1) - pivot*A(l-1,l+1)
      rhs(l) = rhs(l) - pivot*rhs(l-1)
    end do

    ! perform the backsolve
    c(n) = rhs(n)/A(n,n)
    do l=n-1,1,-1
      c(l) = (rhs(l)-A(l,l+1)*c(l+1))/A(l,l)
    end do
    
    ! ok
    status = 0

  end subroutine TriDSolve


end module LE_VDiff
