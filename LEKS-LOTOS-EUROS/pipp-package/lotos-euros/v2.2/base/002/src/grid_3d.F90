!
! 3d grid transforms
!
!
!  call Fill3D( lli, levi, ps, field, &
!               lliX, leviX, fieldX, &
!               combkey, status )
!
!    o lli,levi   : output field defintions  (in)
!      ps         : surface pressure in hPa  (in)
!      field      : output field             (out)
!
!    o lliX, leviX : input field definitons  (in)
!      fieldX      : input field             (in)
!
!    o combkey : 'mass-aver', 'sum', 'aver'   (in)
!
!    o note that psX is not required ...
!
!  call FillMass      ( m    , lli, levi, sp          , status )
!
!    ! Fill 3D mass array given surface pressure and grid cell aera's.
!
!  call FillMassChange( dm_dt, lli, levi, sp_t1, sp_t2, status )
!
!    ! Fill 3D mass change between two times given different
!    ! surface pressures;
!    ! note that ( a + b sp2 ) - ( a + b sp1 )
!    ! is not the same as  a + b (sp2 - sp1) ...
!
!

module grid_3d

  implicit none
  
  ! --- in/out -----------------------------
  
  private
  
  public  ::  Fill3D
  public  ::  FillMass
  public  ::  FillMassChange
  
  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'grid_3d'  
    
  
  
contains


  ! =========================================
  
    
  subroutine Fill3D( lli, levi, nw, ps, field, &
                     lliX, leviX, fieldX, &
                     combkey, status )
                     
    use grid_type_ll , only : TllGridInfo, FillGrid
    use grid_type_hyb, only : TLevelInfo, FillLevels
  
    ! --- in/out --------------------------------

    type(TllGridInfo), intent(in)        ::  lli
    type(TLevelInfo), intent(in)         ::  levi
    character(len=*), intent(in)         ::  nw
    real, intent(in)                     ::  ps(:,:)        ! Pa
    real, intent(out)                    ::  field(:,:,:)
    type(TllGridInfo), intent(in)        ::  lliX
    type(TLevelInfo), intent(in)         ::  leviX
    real, intent(in)                     ::  fieldX(:,:,:)
    character(len=*), intent(in)         ::  combkey
    integer, intent(out)                 ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  name = mname//'/Fill3D'
 
    ! --- local -------------------------
    
    real, allocatable     ::  field_ll(:,:,:)
    integer               ::  l
  
    ! --- begin -------------------------
    
    ! output horizontal grid, input levels
    allocate( field_ll(lli%nlon,lli%nlat,leviX%nlev) )

    select case ( combkey )
    
      !
      ! mass average
      !
      
      case ( 'mass-aver' )
      
        ! horizontal
        do l = 1, leviX%nlev
          call FillGrid( lli , 'n', field_ll(:,:,l), &
                         lliX, 'n', fieldX(:,:,l), 'area-aver', status )
          if (status<0) then
            write (*,'("ERROR - only part of target grid filled")')
            write (*,'("ERROR in ",a)') name; status=1; return
          end if
          if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if
        end do

        ! vertical
        call FillLevels( levi, nw, ps, field, leviX, field_ll, 'mass-aver', status )
        if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if
      
      !
      ! other (should be supported by FillGrid and FillLevels)
      !
    
      case default
      
        ! horizontal
        do l = 1, leviX%nlev
          call FillGrid( lli , 'n', field_ll(:,:,l), &
                         lliX, 'n', fieldX(:,:,l), combkey, status )
          if (status<0) then
            write (*,'("ERROR - only part of target grid filled")')
            write (*,'("ERROR in ",a)') name; status=1; return
          end if
          if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if
        end do

        ! vertical
        call FillLevels( levi, nw, ps, field, leviX, field_ll, combkey, status )
        if (status/=0) then; write (*,'("ERROR in ",a)') name; status=1; return; end if

    end select
    
    ! done
    deallocate( field_ll )
    
    ! ok
    status = 0
    
  end subroutine Fill3D
  


  ! **************************************************************  

  
  !
  !   p = a + b sp
  ! 


  subroutine FillMass( m, lli, levi, sp, status )

    use Binas        , only : grav
    use grid_type_ll , only : TllGridInfo, AreaOper
    use grid_type_hyb, only : TLevelInfo

    ! --- begin ---------------------------------

    real, intent(out)                    ::  m(:,:,:)   ! kg
    type(TllGridInfo), intent(in)        ::  lli
    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp(:,:)       ! Pa
    integer, intent(out)                 ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/FillMass'

    ! --- local -------------------------

    integer               ::  l

    ! --- begin -------------------------

    ! check shape of target grid:
    if ( (size(m,1) /= lli%nlon ) .or. (size(m,2) /= lli%nlat) .or. &
         (size(m,3) /= levi%nlev) ) then
      write (*,'("ERROR - target array does not match with grid definition:")')
      write (*,'("ERROR -   lli    : ",i3," x ",i3         )') lli%nlon, lli%nlat
      write (*,'("ERROR -   levi   : ",i3                  )') levi%nlev
      write (*,'("ERROR -   ll     : ",i3," x ",i3," x ",i3)') shape(m)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! Pa = kg g / A   ->  kg = A * Pa/g

    ! loop over levels
    do l = 1, levi%nlev
      m(:,:,l) = levi%da(l) + levi%db(l) * sp / grav    ! Pa/g = kg/m2
      call AreaOper( lli, m(:,:,l), '*', 'm2', status )          ! kg
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    end do

    ! ok
    status = 0

  end subroutine FillMass


  ! ***


  ! 
  ! p1 = a + b sp1
  ! p2 = a + b sp2
  ! 
  ! p2 - p1 = b (sp2 - sp1)
  !
  ! m = (p2 - p1)/g * A
  !
  

  subroutine FillMassChange( dm, lli, levi, sp1, sp2, status )

    use Binas        , only : grav
    use grid_type_ll , only : TllGridInfo, AreaOper
    use grid_type_hyb, only : TLevelInfo

    ! --- begin ---------------------------------

    real, intent(out)                    ::  dm(:,:,:)   ! kg
    type(TllGridInfo), intent(in)        ::  lli
    type(TLevelInfo), intent(in)         ::  levi
    real, intent(in)                     ::  sp1(:,:)       ! Pa
    real, intent(in)                     ::  sp2(:,:)       ! Pa
    integer, intent(out)                 ::  status

    ! --- const ---------------------------------

    character(len=*), parameter  ::  rname = mname//'/FillMassChange'

    ! --- local -------------------------

    integer               ::  l

    ! --- begin -------------------------

    ! check shape of target grid:
    if ( (size(dm,1) /= lli%nlon ) .or. (size(dm,2) /= lli%nlat) .or. &
         (size(dm,3) /= levi%nlev) ) then
      write (*,'("ERROR - target array does not match with grid definition:")')
      write (*,'("ERROR -   lli    : ",i3," x ",i3         )') lli%nlon, lli%nlat
      write (*,'("ERROR -   levi   : ",i3                  )') levi%nlev
      write (*,'("ERROR -   ll     : ",i3," x ",i3," x ",i3)') shape(dm)
      write (*,'("ERROR in ",a)') rname; status=1; return
    end if

    ! Pa = kg g / A   ->  kg = A * Pa/g

    ! loop over levels
    do l = 1, levi%nlev
      dm(:,:,l) = abs(levi%db(l)) * ( sp2 - sp1 ) / grav    ! Pa/g = kg/m2
      call AreaOper( lli, dm(:,:,l), '*', 'm2', status )          ! kg
      if (status/=0) then; write (*,'("ERROR in ",a)') rname; status=1; return; end if
    end do


    ! ok
    status = 0

  end subroutine FillMassChange




end module grid_3d
