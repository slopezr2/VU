!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
! C3PO_Grid_Convert - grid conversion tools
!
! DESCRIPTION
!
!   This module provides a data type that can be used to remap fields
!   from one grid definition to the other, and eventually repeat this many times.
!
!   If a certain remapping for a particular combination of grid definitions
!   is requested for the first time, then the conversion data type computes
!   remaping weights and stores them for later re-use. New calls to the convertor
!   will check if weights are available already, and use these to speedup
!   the actual calculations.
!
! USAGE
!
!   use C3PO_Grid_Convert, only : T_Grid_Convertors
!   use Grid_Type_LL     , only : TllGridInfo
!   use C3PO_Grid_Rgg    , only : T_Grid_Rgg
!
!   ! variables:
!   class(T_Grid_Convertors)                 ::  conv
!   ! return status, ok if zero:
!   integer                                  ::  status
!   ! reduced-gaussian-grid definition and field:
!   class(T_Grid_Rgg)                        ::  rgg
!   real, allocatable                        ::  gg(:)
!   ! regular lon/lat grid and field:
!   type(TllGridInfo)                        ::  lli
!   real, allocatable                        ::  ll(:,:)
!
!   ! initialize conversions:
!   call conv%Init( status )
!   if (status/=0) stop
!
!   ! define reduced-guassian-grid :
!   call rgg%Init( ..., status )
!   if (status/=0) stop
!   ! data field defined on this grid:
!   allocate( gg(rgg%npoint) )
!
!   ! define regular lon/lat grid:
!   call Init( lli, ..., status )
!   if (status/=0) stop
!   ! data field defined on this grid:
!   allocate( ll(lli%nlon,lli%nlat) )
!   
!   ! remap field defined on reduced-gaussian-grid to lon/lat grid
!   ! using area-weighted averages:
!   call conv%Rgg_LL_AreaAver( rgg, gg, lli, ll, status )
!   if (status/=0) stop
!
!   ! interpolate from field defined on reduced-gaussian-grid to points
!   ! on lon/lat grid using bi-linear interpolation:
!   call conv%Rgg_LL_Interpol( rgg, gg, lli, ll, status )
!   if (status/=0) stop
!
!   ! clear:
!   call conv%Done( status )
!   if (status/=0) stop
!
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action)  if (status> 0) then; TRACEBACK; action; return; end if
!
#include "c3po.inc"
!
!######################################################################

module C3PO_Grid_Convert

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Grid_Convertors
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Grid_Convert' 
  
  ! number of rea/integer key values:
  integer, parameter ::  nkeyval = 50
  
  ! maximum number of convertors:
  integer,parameter  ::  n_convertor_max = 50


  ! --- type ------------------------------------
  
  type T_Grid_Key
    character(len=32)     ::  ckey
    real                  ::  rkey(nkeyval)
    integer               ::  ikey(nkeyval)
  contains
    procedure   ::  Init    => Grid_Key_Init
    procedure   ::  Done    => Grid_Key_Done
    procedure   ::  Compare => Grid_Key_Compare 
  end type T_Grid_Key
  
  ! *
  
  type T_Grid_Convertor
    character(len=32)     ::  howto
    type(T_Grid_Key)      ::  grid_in
    type(T_Grid_Key)      ::  grid_out
  contains
    procedure   ::  Compare     =>  Grid_Convertor_Compare
  end type T_Grid_Convertor
  
  ! *
  
  type, extends(T_Grid_Convertor) :: T_Grid_Convertor_Cgg_LL
    integer               ::  nlon, nlat
    integer               ::  ni, nj
    integer, allocatable  ::  ncov(:,:)        ! (nlon,nlat)
    integer, allocatable  ::  ii(:,:,:)        ! (nlon,nlat,ncov)
    integer, allocatable  ::  jj(:,:,:)        ! (nlon,nlat,ncov)
    real, allocatable     ::  frac(:,:,:)      ! (nlon,nlat,ncov)
    real, allocatable     ::  A_cg(:,:)        ! (ni,nj) area [m2]
    real, allocatable     ::  A_ll(:,:)        ! (nlon,nlat) area [m2]
  contains
    procedure   ::  Init_AreaRemap    => Grid_Convertor_Cgg_LL_Init_AreaRemap
    procedure   ::  Init_Interpol     => Grid_Convertor_Cgg_LL_Init_Interpol
    procedure   ::  Init_Nearest      => Grid_Convertor_Cgg_LL_Init_Nearest
    procedure   ::  Done              => Grid_Convertor_Cgg_LL_Done
    procedure   ::                       Grid_Convertor_Cgg_LL_AreaAver_2d
    procedure   ::                       Grid_Convertor_Cgg_LL_AreaAver_3d
    generic     ::  AreaAver          => Grid_Convertor_Cgg_LL_AreaAver_2d, &
                                         Grid_Convertor_Cgg_LL_AreaAver_3d
    procedure   ::                       Grid_Convertor_Cgg_LL_Interpol_2d
    procedure   ::                       Grid_Convertor_Cgg_LL_Interpol_3d
    generic     ::  Interpol          => Grid_Convertor_Cgg_LL_Interpol_2d, &
                                         Grid_Convertor_Cgg_LL_Interpol_3d
  end type T_Grid_Convertor_Cgg_LL
  
  ! *
  
  type, extends(T_Grid_Convertor) :: T_Grid_Convertor_Rgg_LL
    integer               ::  nlon, nlat
    integer               ::  np
    integer, allocatable  ::  ncov(:,:)        ! (nlon,nlat)
    integer, allocatable  ::  indx(:,:,:)      ! (nlon,nlat,ncov)
    real, allocatable     ::  frac(:,:,:)      ! (nlon,nlat,ncov)
    real, allocatable     ::  A_gg(:)          ! (npoint)    area [m2]
    real, allocatable     ::  A_ll(:,:)        ! (nlon,nlat) area [m2]
  contains
    procedure   ::  Init_AreaRemap    => Grid_Convertor_Rgg_LL_Init_AreaRemap
    procedure   ::  Init_Interpol     => Grid_Convertor_Rgg_LL_Init_Interpol
    procedure   ::  Done              => Grid_Convertor_Rgg_LL_Done
    procedure   ::                       Grid_Convertor_Rgg_LL_AreaAver_2d
    procedure   ::                       Grid_Convertor_Rgg_LL_AreaAver_3d
    generic     ::  AreaAver          => Grid_Convertor_Rgg_LL_AreaAver_2d, &
                                         Grid_Convertor_Rgg_LL_AreaAver_3d
    procedure   ::                       Grid_Convertor_Rgg_LL_Interpol_2d
    procedure   ::                       Grid_Convertor_Rgg_LL_Interpol_3d
    generic     ::  Interpol          => Grid_Convertor_Rgg_LL_Interpol_2d, &
                                         Grid_Convertor_Rgg_LL_Interpol_3d
  end type T_Grid_Convertor_Rgg_LL
 
  ! *
  
  type, extends(T_Grid_Convertor) :: T_Grid_Convertor_Ugg
    ! source cells:
    integer               ::  ni, nj
    real, allocatable     ::  A_ug_from(:,:)   ! (ni,nj) area [m2]
    ! target cells:
    integer               ::  nlon, nlat
    real, allocatable     ::  A_ug_to(:,:)     ! (nlon,nlat) area [m2]
    ! number of source cells:
    integer, allocatable  ::  ncov(:,:)        ! (nlon,nlat)
    ! source cell indices and weights:
    integer, allocatable  ::  ii(:,:,:)        ! (nlon,nlat,ncov)
    integer, allocatable  ::  jj(:,:,:)        ! (nlon,nlat,ncov)
    real, allocatable     ::  frac(:,:,:)      ! (nlon,nlat,ncov)
    !
    ! ... NEW ......................
    ! range of source cells in 1D arrays:
    integer, allocatable  ::  froms(:,:,:)     ! (nlon,nlat,2)
    ! source cell indices and weight:
    integer               ::  nfrom
    integer               ::  nstep
    integer, pointer      ::  ifrom(:)        ! (total froms)
    integer, pointer      ::  jfrom(:)        ! (total froms)
    real, pointer         ::  wfrom(:)        ! (total froms)
    !
  contains
    procedure   ::  InitMapping            => Grid_Convertor_Ugg_InitMapping
    procedure   ::  AddMapping             => Grid_Convertor_Ugg_AddMapping
    procedure   ::  Done                   => Grid_Convertor_Ugg_Done
    procedure   ::  Init_AreaRemap         => Grid_Convertor_Ugg_Init_AreaRemap
    procedure   ::  Init_Interpol          => Grid_Convertor_Ugg_Init_Interpol
    procedure   ::  Init_Nearest           => Grid_Convertor_Ugg_Init_Nearest
    !procedure   ::  Init_IndexFractions    => Grid_Convertor_Ugg_Init_IndexFractions
    procedure   ::                            Grid_Convertor_Ugg_AreaAver_2d
    procedure   ::                            Grid_Convertor_Ugg_AreaAver_3d
    generic     ::  AreaAver               => Grid_Convertor_Ugg_AreaAver_2d, &
                                              Grid_Convertor_Ugg_AreaAver_3d
    procedure   ::                            Grid_Convertor_Ugg_Interpol_2d
    procedure   ::                            Grid_Convertor_Ugg_Interpol_3d
    generic     ::  Interpol               => Grid_Convertor_Ugg_Interpol_2d, &
                                              Grid_Convertor_Ugg_Interpol_3d
    !procedure   ::  IndexFractions         => Grid_Convertor_Ugg_IndexFractions
  end type T_Grid_Convertor_Ugg
  
  ! *
  
  type T_Grid_Convertors
    ! reduced to regular:
    integer                                      ::  n_rgg_to_ll
    type(T_Grid_Convertor_Rgg_LL), allocatable   ::  rgg_to_ll(:)
    ! cartesian to regular:
    integer                                      ::  n_cgg_to_ll
    type(T_Grid_Convertor_Cgg_LL), allocatable   ::  cgg_to_ll(:)
    ! universal to universal:
    integer                                      ::  n_ugg_to_ugg
    type(T_Grid_Convertor_Ugg), allocatable      ::  ugg_to_ugg(:)
  contains
    !
    procedure   ::  Init            => Grid_Convertors_Init
    procedure   ::  Done            => Grid_Convertors_Done

    procedure   ::  Rgg_LL_Select   => Grid_Convertors_Rgg_LL_Select
    procedure   ::  Cgg_LL_Select   => Grid_Convertors_Cgg_LL_Select
    procedure   ::  Ugg_Select      => Grid_Convertors_Ugg_Select

    ! ***
    
    procedure   ::                     Grid_Convertors_Cgg_LL_AreaAver_2d
    procedure   ::                     Grid_Convertors_Cgg_LL_AreaAver_3d
    procedure   ::                     Grid_Convertors_Cgg_LL_AreaAver_2d_uv
    procedure   ::                     Grid_Convertors_Cgg_LL_AreaAver_3d_uv
    generic     ::  Cgg_LL_AreaAver => Grid_Convertors_Cgg_LL_AreaAver_2d, &
                                       Grid_Convertors_Cgg_LL_AreaAver_3d, &
                                       Grid_Convertors_Cgg_LL_AreaAver_2d_uv, &
                                       Grid_Convertors_Cgg_LL_AreaAver_3d_uv
    !
    procedure   ::                     Grid_Convertors_Cgg_LL_Interpol_2d
    procedure   ::                     Grid_Convertors_Cgg_LL_Interpol_3d
    procedure   ::                     Grid_Convertors_Cgg_LL_Interpol_3d_uv
    generic     ::  Cgg_LL_Interpol => Grid_Convertors_Cgg_LL_Interpol_2d, &
                                       Grid_Convertors_Cgg_LL_Interpol_3d, &
                                       Grid_Convertors_Cgg_LL_Interpol_3d_uv
    !
    procedure   ::  Cgg_LL_Nearest  => Grid_Convertors_Cgg_LL_Nearest_2d
    
    ! ***
    
    procedure   ::                     Grid_Convertors_Ugg_AreaAver_2d
    procedure   ::                     Grid_Convertors_Ugg_AreaAver_3d
    procedure   ::                     Grid_Convertors_Ugg_AreaAver_2d_uv
    procedure   ::                     Grid_Convertors_Ugg_AreaAver_3d_uv
    generic     ::  Ugg_AreaAver    => Grid_Convertors_Ugg_AreaAver_2d, &
                                       Grid_Convertors_Ugg_AreaAver_3d, &
                                       Grid_Convertors_Ugg_AreaAver_2d_uv, &
                                       Grid_Convertors_Ugg_AreaAver_3d_uv
    !
    procedure   ::                     Grid_Convertors_Ugg_Interpol_2d
    procedure   ::                     Grid_Convertors_Ugg_Interpol_2d_uv
    procedure   ::                     Grid_Convertors_Ugg_Interpol_3d
    procedure   ::                     Grid_Convertors_Ugg_Interpol_3d_uv
    generic     ::  Ugg_Interpol    => Grid_Convertors_Ugg_Interpol_2d, &
                                       Grid_Convertors_Ugg_Interpol_2d_uv, &
                                       Grid_Convertors_Ugg_Interpol_3d, &
                                       Grid_Convertors_Ugg_Interpol_3d_uv
    !                                       
    procedure   ::  uGG_Nearest     => Grid_Convertors_Ugg_Nearest_2d    
    !
    procedure   ::  ugg_IndexFractions => Grid_Convertors_Ugg_IndexFractions
    
    ! ***
    
    procedure   ::                     Grid_Convertors_Rgg_LL_AreaAver_2d
    procedure   ::                     Grid_Convertors_Rgg_LL_AreaAver_3d
    generic     ::  Rgg_LL_AreaAver => Grid_Convertors_Rgg_LL_AreaAver_2d, &
                                       Grid_Convertors_Rgg_LL_AreaAver_3d
    !
    procedure   ::                     Grid_Convertors_Rgg_LL_Interpol_2d
    procedure   ::                     Grid_Convertors_Rgg_LL_Interpol_3d
    generic     ::  Rgg_LL_Interpol => Grid_Convertors_Rgg_LL_Interpol_2d, &
                                       Grid_Convertors_Rgg_LL_Interpol_3d
    !
  end type T_Grid_Convertors
  

  ! --- var ------------------------------------
  

contains


  ! ********************************************************************
  ! ***
  ! *** key
  ! ***
  ! ********************************************************************
  
  subroutine Grid_Key_Init( self, ckey, ikey, rkey, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Key), intent(out)        ::  self
    character(len=*), intent(in)          ::  ckey
    integer, intent(in)                   ::  ikey(:)
    real, intent(in)                      ::  rkey(:)
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Key_Init'
    
    ! --- local ----------------------------------
    
    integer   ::  n

    ! --- begin ----------------------------------
    
    ! store:
    self%ckey = trim(ckey)
    
    ! check ...
    n = size(ikey)
    if ( n > size(self%ikey) ) then
      write (gol,'("size of input ikey ",i4," exceeds storage ",i4)') n, size(self%ikey); call goErr
      TRACEBACK; status=1; return
    end if
    ! store:
    self%ikey = 0
    self%ikey(1:n) = ikey
    
    ! check ...
    n = size(rkey)
    if ( n > size(self%rkey) ) then
      write (gol,'("size of input rkey ",i4," exceeds storage ",i4)') n, size(self%rkey); call goErr
      TRACEBACK; status=1; return
    end if
    ! store:
    self%rkey = 0.0
    self%rkey(1:n) = rkey
    
    ! ok
    status = 0

  end subroutine Grid_Key_Init
  
  ! *
  
  subroutine Grid_Key_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Key), intent(inout)        ::  self
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Key_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0

  end subroutine Grid_Key_Done
  
  ! *
  
  ! Compare keys
  
  subroutine Grid_Key_Compare( self, new, match, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Key), intent(in)         ::  self
    class(T_Grid_Key), intent(in)         ::  new
    logical, intent(out)                  ::  match
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Key_Compare'
    
    ! --- local ----------------------------------
    
    integer   ::  n

    ! --- begin ----------------------------------
    
    ! same ?
    match = (trim(self%ckey) == trim(new%ckey)) .and. &
            all(self%ikey == new%ikey) .and. all(self%rkey == new%rkey)
    
    ! ok
    status = 0

  end subroutine Grid_Key_Compare


  ! ********************************************************************
  ! ***
  ! *** convertor
  ! ***
  ! ********************************************************************

  
  subroutine Grid_Convertor_Init( self, howto, grid_in, grid_out, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor), intent(out)    ::  self
    character(len=*), intent(in)            ::  howto
    type(T_Grid_Key), intent(in)            ::  grid_in
    type(T_Grid_Key), intent(in)            ::  grid_out
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Init'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! store:
    self%howto    = trim(howto)
    self%grid_in  = grid_in
    self%grid_out = grid_out

    ! ok
    status = 0

  end subroutine Grid_Convertor_Init
  
  ! *
  
  subroutine Grid_Convertor_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor), intent(inout)    ::  self
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! emtpy:
    self%howto = ''
    ! clear:
    call self%grid_in%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call self%grid_out%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Done
  
  ! *
  
  subroutine Grid_Convertor_Compare( self, howto, grid_in, grid_out, match, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor), intent(in)       ::  self
    character(len=*), intent(in)              ::  howto
    type(T_Grid_Key), intent(in)              ::  grid_in
    type(T_Grid_Key), intent(in)              ::  grid_out
    logical, intent(out)                      ::  match
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Compare'
    
    ! --- local ----------------------------------
    
    logical   ::  match_in, match_out
    
    ! --- begin ----------------------------------
    
    ! compare input grids:
    call self%grid_in%Compare( grid_in, match_in, status )
    IF_NOT_OK_RETURN(status=1)
    ! compare output grids:
    call self%grid_out%Compare( grid_out, match_out, status )
    IF_NOT_OK_RETURN(status=1)
    ! combine:
    match = (trim(self%howto) == trim(howto)) .and. match_in .and. match_out
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Compare


  ! ********************************************************************
  ! ***
  ! *** rgg to ll
  ! ***
  ! ********************************************************************

  
  subroutine Grid_Convertor_Rgg_LL_Init_AreaRemap( self, howto, rgg_key, rgg, lli_key, lli, status )
  
    use num          , only : Interval  
    use grid_tools   , only : pi, ll_area_frac_deg
    use grid_type_ll , only : TllGridInfo, AreaOper
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Rgg_LL), intent(out)     ::  self
    character(len=*), intent(in)                    ::  howto
    type(T_Grid_Key), intent(in)                    ::  rgg_key
    class(T_Grid_Rgg), intent(in)                   ::  rgg
    type(T_Grid_Key), intent(in)                    ::  lli_key
    type(TllGridInfo), intent(in)                   ::  lli
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Rgg_LL_Init_AreaRemap'
    
    ! --- local ----------------------------------
    
    integer              ::  max_ncov_lon, max_ncov_lat, max_ncov
    integer              ::  i, j, ip
    integer              ::  iflag
    integer              ::  gi, gi1, gi2
    integer              ::  gj, gj1, gj2
    integer              ::  nlon
    real, allocatable    ::  gblon(:)  ! (0:ggi%nulon)
    real                 ::  west1, east1, south1, north1
    real                 ::  west2, east2, south2, north2
    real                 ::  frac
    integer              ::  ncov

    ! --- begin ----------------------------------

    ! fill keys:
    call Grid_Convertor_Init( self, howto, rgg_key, lli_key, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! save dimension
    self%nlon = lli%nlon
    self%nlat = lli%nlat
    self%np   = rgg%npoint
    
    ! estimate maximum number of gg cells covering an ll cell:
    max_ncov_lon = ceiling( lli%dlon_deg / minval(abs(rgg%longitude_bnds(2,:)-rgg%longitude_bnds(1,:))) )
    max_ncov_lat = ceiling( lli%dlat_deg / minval(abs(rgg%latitude_bnds (2,:)-rgg%latitude_bnds (1,:))) )
    max_ncov = (max_ncov_lon+1) * (max_ncov_lat+1)
    
    ! allocate arrays
    allocate( self%ncov(lli%nlon,lli%nlat) )
    allocate( self%indx(lli%nlon,lli%nlat,max_ncov) )
    allocate( self%frac(lli%nlon,lli%nlat,max_ncov) )
    ! temporary:
    allocate( gblon(0:rgg%nulon) )
    
    ! zero by default:
    self%ncov = 0
    self%indx = 0
    self%frac = 0.0
    
    gj1 = 0
    gj2 = 0
    
    ! loop over ll cells from north to south (gg direction):
    do j = lli%nlat, 1, -1
      do i = 1, lli%nlon
     
        ! extract boundaries of ll cell:
        west2  = lli%blon_deg(i-1)
        east2  = lli%blon_deg(i)
        if ( east2 < 0.0 ) then
          west2 = west2 + 360.0    ! (0,360)
          east2 = east2 + 360.0    ! (0,360)
        end if
        south2 = lli%blat_deg(j-1)
        north2 = lli%blat_deg(j)

        ! search gg rows including north and south ll lat;
        ! negative lats to get increasing values ...
        !call Interval( -rgg%blat, -north2  , gj1, iflag )
        !if ( iflag /= 0 ) stop 'BUG IN gg2ll_Init : wrong iflag gj1'
        !call Interval( -rgg%blat, -south2, gj2, iflag )
        !if ( iflag /= 0 ) stop 'BUG IN gg2ll_Init : wrong iflag gj2'
        call rgg%GetBandIndex( north2, gj1, status )
        IF_NOT_OK_RETURN(status=1)
        call rgg%GetBandIndex( south2, gj2, status )
        IF_NOT_OK_RETURN(status=1)
        
        gi1 = 0
        gi2 = 0
        ! loop over gg lat rows:
        do gj = gj1, gj2
        
          ! boundary lons
          nlon = rgg%band_nlon(gj)
          do gi = 1, nlon
            ip = rgg%band_i0(gj) + gi
            if ( gi == 1 ) gblon(0) = rgg%longitude_bnds(1,ip)
            gblon(gi) = rgg%longitude_bnds(2,ip)
          end do
          
          ! search cells including west and east bound of ll cell
          if ( west2 < gblon(0) ) then
            call Interval( gblon(0:nlon), west2+360.0, gi1, iflag )
          else if ( west2 > gblon(nlon) ) then
            call Interval( gblon(0:nlon), west2-360.0, gi1, iflag )
          else
            call Interval( gblon(0:nlon), west2      , gi1, iflag )
          end if
          if ( iflag /= 0 ) then
            print *, 'gblon=', gblon(0:nlon)
            print *, 'west2=',west2
            print *, 'iflag=',iflag
            stop 'BUG IN gg2ll_Init : wrong iflag gi1'
          end if
          if ( east2 < gblon(0) ) then
            call Interval( gblon(0:nlon), east2+360.0, gi2, iflag )
          else if ( east2 > gblon(nlon) ) then
            call Interval( gblon(0:nlon), east2-360.0, gi2, iflag )
          else
            call Interval( gblon(0:nlon), east2      , gi2, iflag )
          end if
          if ( iflag /= 0 ) then
            print *, 'gblon=', gblon(0:nlon)
            print *, 'east2=',east2
            stop 'BUG IN gg2ll_Init : wrong iflag gi2'
          end if
        
          ! loop over all gg cells in current row:
          gi = gi1
          do 
          
            ! extract boundaries of gg cell:
            !west1  = (gi-1.5)*rgg%dlon(gj)    ! (0,2pi)
            !east1  = (gi-0.5)*rgg%dlon(gj)    ! (0,2pi)
            !south1 = rgg%blat(gj)
            !north1 = rgg%blat(gj-1)
            ! point index:
            ip = rgg%band_i0(gj) + gi
            ! bounds:
            west1  = minval(rgg%longitude_bnds(:,ip))
            east1  = maxval(rgg%longitude_bnds(:,ip))
            south1 = minval(rgg%latitude_bnds (:,ip))
            north1 = maxval(rgg%latitude_bnds (:,ip))

            ! shift if gg cell is right from [0,360]
            if ( west1 > east2 ) then
              west1 = west1 - 360.0    ! (0,360)
              east1 = east1 - 360.0    ! (0,360)
            end if

            ! determine covarage fraction:
            call ll_area_frac_deg( west1, east1, south1, north1, &
                                   west2, east2, south2, north2, &
                                   frac, status )
            IF_NOT_OK_RETURN(status=1)

            ! fill fraction:
            if ( frac > 0.0 .and. frac <= 1.0 ) then
              ! add contribution:
              ncov = self%ncov(i,j) + 1
              self%ncov(i,j     ) = ncov
              self%indx(i,j,ncov) = rgg%band_i0(gj) + gi
              self%frac(i,j,ncov) = frac
              
            else if ( abs(frac) < 1.0e-4 ) then
              ! almost no coverage ...
              
            else if ( abs(frac-1.0) < 1.0e-4 ) then
              ncov = self%ncov(i,j) + 1
              self%ncov(i,j     ) = ncov
              self%indx(i,j,ncov) = rgg%band_i0(gj) + gi
              self%frac(i,j,ncov) = 1.0
              
            else 
              write (gol,'("BUG - found strange coverage fraction : ",e15.6)') frac; call goErr
              write (gol,'("BUG - source cell : ",4f12.6)') west1, east1, south1, north1; call goErr
              write (gol,'("BUG - target cell : ",4f12.6)') west2, east2, south2, north2; call goErr
              TRACEBACK; status=1; return
            end if

            if ( gi == gi2 ) exit
            gi = gi + 1
            if ( gi == nlon+1 ) gi = 1

          end do  ! rgg cells in row
          
        end do  ! rgg bands
        
        ! check ...
        if ( self%ncov(i,j)  == 0 ) then
          write (gol,'("found target cell not covered by source grid")'); call goErr
          TRACEBACK; status=1; return
        end if

      end do  ! lli i
    end do  ! lli j
    
    ! clear:
    deallocate( gblon )
    
    ! storage for cell area's in rgg:
    allocate( self%A_gg(rgg%npoint) )
    ! copy:
    self%A_gg = rgg%area  ! m2

    ! storage for cell area in lon/lat grid:
    allocate( self%A_ll(lli%nlon,lli%nlat) )
    ! fill:
    call AreaOper( lli, self%A_ll, '=', 'm2', status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertor_Rgg_LL_Init_AreaRemap
  

  ! *

  
  subroutine Grid_Convertor_Rgg_LL_Init_Interpol( self, howto, rgg_key, rgg, lli_key, lli, status )
  
    use grid_type_ll , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Rgg_LL), intent(out)     ::  self
    character(len=*), intent(in)                    ::  howto
    type(T_Grid_Key), intent(in)                    ::  rgg_key
    class(T_Grid_Rgg), intent(in)                   ::  rgg
    type(T_Grid_Key), intent(in)                    ::  lli_key
    type(TllGridInfo), intent(in)                   ::  lli
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Rgg_LL_Init_Interpol'
    
    ! --- local ----------------------------------
    
    integer              ::  max_np
    integer              ::  i, j
    real                 ::  lon, lat

    ! --- begin ----------------------------------

    ! fill keys:
    call Grid_Convertor_Init( self, howto, rgg_key, lli_key, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! save dimension
    self%nlon = lli%nlon
    self%nlat = lli%nlat
    self%np   = rgg%npoint
    
    ! maximum number of surrounding points for bilinear interpolation:
    max_np = 4
    
    ! allocate arrays
    allocate( self%ncov(lli%nlon,lli%nlat) )
    allocate( self%indx(lli%nlon,lli%nlat,max_np) )
    allocate( self%frac(lli%nlon,lli%nlat,max_np) )
    
    ! zero by default:
    self%ncov = 0
    self%indx = 0
    self%frac = 0.0
    
    ! loop over ll cells:
    do j = 1, lli%nlat
    
      ! target latitude:
      lat = lli%lat_deg(j)
    
      ! loop over lon cells    
      do i = 1, lli%nlon

        ! target longitude:
        lon = lli%lon_deg(i)
        
        ! get weights:
        call rgg%GetInterpol( lon, lat, self%ncov(i,j), self%indx(i,j,:), self%frac(i,j,:), status )
        IF_NOT_OK_RETURN(status=1)

      end do  ! lli i

    end do  ! lli j
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Rgg_LL_Init_Interpol
  

  ! *
  

  subroutine Grid_Convertor_Rgg_LL_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Rgg_LL), intent(inout)     ::  self
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Rgg_LL_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    call Grid_Convertor_Done( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! deallocate arrays
    deallocate( self%ncov )
    deallocate( self%indx )
    deallocate( self%frac )
    if ( allocated( self%A_gg) ) deallocate( self%A_gg )
    if ( allocated( self%A_ll) ) deallocate( self%A_ll )
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Rgg_LL_Done
  
  
  ! *


  subroutine Grid_Convertor_Rgg_LL_AreaAver_2d( self, gg, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Rgg_LL), intent(inout)     ::  self
    real, intent(in)                                  ::  gg(:)        ! (npoint)
    real, intent(out)                                 ::  ll(:,:)      ! (nlon,nlat)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Rgg_LL_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer            ::  i, j, k
    integer            ::  ncov, indx
    
    ! --- begin ----------------------------------
    
    ! init output:
    ll = 0.0
    ! loop over target cells:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(i,j)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do k = 1, ncov
            ! source index:
            indx = self%indx(i,j,k)
            ! add contribution:
            ll(i,j) = ll(i,j) + gg(indx) * self%A_gg(indx) * self%frac(i,j,k)
          end do
        end if
        ! area average:
        ll(i,j) = ll(i,j) / self%A_ll(i,j)
      end do ! cells i
    end do ! cells j
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Rgg_LL_AreaAver_2d
  
  
  ! *


  subroutine Grid_Convertor_Rgg_LL_AreaAver_3d( self, gg, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Rgg_LL), intent(inout)     ::  self
    real, intent(in)                                  ::  gg(:,:)    ! (npoint,nlev)
    real, intent(out)                                 ::  ll(:,:,:)  ! (nlon,nlat,nlev)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Rgg_LL_AreaAver_3d'
    
    ! --- local ----------------------------------
    
    integer            ::  i, j, k
    integer            ::  ncov, indx
    
    ! --- begin ----------------------------------
    
    ! init output:
    ll = 0.0
    ! loop over target cells:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(i,j)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do k = 1, ncov
            ! source index:
            indx = self%indx(i,j,k)
            ! add contribution:
            ll(i,j,:) = ll(i,j,:) + gg(indx,:) * self%A_gg(indx) * self%frac(i,j,k)
          end do
        end if
        ! area average:
        ll(i,j,:) = ll(i,j,:) / self%A_ll(i,j)
      end do ! cells i
    end do ! cells j
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Rgg_LL_AreaAver_3d
  
  
  ! *


  subroutine Grid_Convertor_Rgg_LL_Interpol_2d( self, gg, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Rgg_LL), intent(inout)     ::  self
    real, intent(in)                                  ::  gg(:)        ! (npoint)
    real, intent(out)                                 ::  ll(:,:)      ! (nlon,nlat)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Rgg_LL_Interpol_2d'
    
    ! --- local ----------------------------------
    
    integer            ::  i, j, k
    integer            ::  ncov, indx
    
    ! --- begin ----------------------------------
    
    ! init output:
    ll = 0.0
    ! loop over target cells:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(i,j)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do k = 1, ncov
            ! source index:
            indx = self%indx(i,j,k)
            ! add contribution:
            ll(i,j) = ll(i,j) + gg(indx) * self%frac(i,j,k)
          end do
        end if
      end do ! cells i
    end do ! cells j
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Rgg_LL_Interpol_2d
  
  
  ! *


  subroutine Grid_Convertor_Rgg_LL_Interpol_3d( self, gg, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Rgg_LL), intent(inout)     ::  self
    real, intent(in)                                  ::  gg(:,:)        ! (npoint,nlev)
    real, intent(out)                                 ::  ll(:,:,:)      ! (nlon,nlat,nlev)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Rgg_LL_Interpol_3d'
    
    ! --- local ----------------------------------
    
    integer            ::  i, j, k
    integer            ::  ncov, indx
    
    ! --- begin ----------------------------------
    
    ! init output:
    ll = 0.0
    ! loop over target cells:
    do j = 1, self%nlat
      do i = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(i,j)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do k = 1, ncov
            ! source index:
            indx = self%indx(i,j,k)
            ! add contribution:
            ll(i,j,:) = ll(i,j,:) + gg(indx,:) * self%frac(i,j,k)
          end do
        end if
      end do ! cells i
    end do ! cells j
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Rgg_LL_Interpol_3d


  ! ********************************************************************
  ! ***
  ! *** cgg to ll
  ! ***
  ! ********************************************************************

  
  subroutine Grid_Convertor_Cgg_LL_Init_AreaRemap( self, howto, cgg_key, cgg, lli_key, lli, status )
  
    use num          , only : Interval  
    use grid_tools   , only : pi, ll_area_frac_deg
    use grid_type_ll , only : TllGridInfo, AreaOper
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Cgg_LL), intent(out)     ::  self
    character(len=*), intent(in)                    ::  howto
    type(T_Grid_Key), intent(in)                    ::  cgg_key
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    type(T_Grid_Key), intent(in)                    ::  lli_key
    type(TllGridInfo), intent(in)                   ::  lli
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Cgg_LL_Init_AreaRemap'
    
    ! --- local ----------------------------------
    
    integer              ::  max_ncov_lon, max_ncov_lat, max_ncov
    integer              ::  ilon, ilat
    integer              ::  i, i1, i2, idir
    integer              ::  j, j1, j2, jdir
    real                 ::  west1, east1, south1, north1
    real                 ::  west2, east2, south2, north2
    real                 ::  frac
    integer              ::  ncov

    ! --- begin ----------------------------------

    ! fill keys:
    call Grid_Convertor_Init( self, howto, cgg_key, lli_key, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! save dimension
    self%nlon = lli%nlon
    self%nlat = lli%nlat
    self%ni   = cgg%nlon
    self%nj   = cgg%nlat
    
    ! estimate maximum number of gg cells covering an ll cell:
    max_ncov_lon = ceiling( lli%dlon_deg / minval(abs(cgg%longitude_bnds(2,:)-cgg%longitude_bnds(1,:))) )
    max_ncov_lat = ceiling( lli%dlat_deg / minval(abs(cgg%latitude_bnds (2,:)-cgg%latitude_bnds (1,:))) )
    max_ncov = (max_ncov_lon+1) * (max_ncov_lat+1)
    
    ! allocate arrays
    allocate( self%ncov(lli%nlon,lli%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ii  (lli%nlon,lli%nlat,max_ncov), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj  (lli%nlon,lli%nlat,max_ncov), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%frac(lli%nlon,lli%nlat,max_ncov), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! zero by default:
    self%ncov = 0
    self%ii   = 0
    self%jj   = 0
    self%frac = 0.0
    
    ! loop over target cells:
    do ilat = 1, lli%nlat
      do ilon = 1, lli%nlon
     
        ! extract boundaries of ll cell:
        west2  = lli%blon_deg(ilon-1)
        east2  = lli%blon_deg(ilon)
        if ( east2 < 0.0 ) then
          west2 = west2 + 360.0    ! (0,360)
          east2 = east2 + 360.0    ! (0,360)
        end if
        south2 = lli%blat_deg(ilat-1)
        north2 = lli%blat_deg(ilat)

        ! search cgg rows including west and east ll lat:
        call cgg%GetLonIndex( west2, i1, status )
        IF_NOT_OK_RETURN(status=1)
        call cgg%GetLonIndex( east2, i2, status )
        IF_NOT_OK_RETURN(status=1)
        ! direction:
        idir = sign( 1, i2 - i1 )

        ! search cgg rows including north and south ll lat:
        call cgg%GetLatIndex( south2, j1, status )
        IF_NOT_OK_RETURN(status=1)
        call cgg%GetLatIndex( north2, j2, status )
        IF_NOT_OK_RETURN(status=1)
        ! direction:
        jdir = sign( 1, j2 - j1 )
        
        ! loop over source cells:
        do j = j1, j2, jdir
          do i = i1, i2, idir
        
            ! bounds:
            west1  = minval(cgg%longitude_bnds(:,i))
            east1  = maxval(cgg%longitude_bnds(:,i))
            south1 = minval(cgg%latitude_bnds (:,j))
            north1 = maxval(cgg%latitude_bnds (:,j))

            ! shift if gg cell is right from [0,360]
            if ( west1 > east2 ) then
              west1 = west1 - 360.0    ! (0,360)
              east1 = east1 - 360.0    ! (0,360)
            end if

            ! determine covarage fraction:
            call ll_area_frac_deg( west1, east1, south1, north1, &
                                   west2, east2, south2, north2, &
                                   frac, status )
            IF_NOT_OK_RETURN(status=1)

            ! fill fraction:
            if ( frac > 0.0 .and. frac <= 1.0 ) then
              ! add contribution:
              ncov = self%ncov(ilon,ilat) + 1
              self%ncov(ilon,ilat     ) = ncov
              self%ii  (ilon,ilat,ncov) = i
              self%jj  (ilon,ilat,ncov) = j
              self%frac(ilon,ilat,ncov) = frac
              
            else if ( abs(frac) < 1.0e-4 ) then
              ! almost no coverage ...
              
            else if ( abs(frac-1.0) < 1.0e-4 ) then
              ! 100% coverage, but small rounding error
              ncov = self%ncov(ilon,ilat) + 1
              self%ncov(ilon,ilat     ) = ncov
              self%ii  (ilon,ilat,ncov) = i
              self%jj  (ilon,ilat,ncov) = j
              self%frac(ilon,ilat,ncov) = 1.0
              
            else 
              write (gol,'("BUG - found strange coverage fraction : ",e15.6)') frac; call goErr
              write (gol,'("BUG - source cell : ",4f12.6)') west1, east1, south1, north1; call goErr
              write (gol,'("BUG - target cell : ",4f12.6)') west2, east2, south2, north2; call goErr
              TRACEBACK; status=1; return
            end if

          end do  ! cgg i
        end do  ! cgg j
        
        ! check ...
        if ( self%ncov(ilon,ilat)  == 0 ) then
          write (gol,'("found target cell not covered by source grid")'); call goErr
          TRACEBACK; status=1; return
        end if

      end do  ! lli ilon
    end do  ! lli ilat
    
    ! storage for cell area's in cgg:
    allocate( self%A_cg(cgg%nlon,cgg%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    self%A_cg = cgg%area  ! m2

    ! storage for cell area in lon/lat grid:
    allocate( self%A_ll(lli%nlon,lli%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    call AreaOper( lli, self%A_ll, '=', 'm2', status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertor_Cgg_LL_Init_AreaRemap
  

  ! *

  
  subroutine Grid_Convertor_Cgg_LL_Init_Interpol( self, howto, cgg_key, cgg, lli_key, lli, status )
  
    use grid_type_ll , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Cgg_LL), intent(out)     ::  self
    character(len=*), intent(in)                    ::  howto
    type(T_Grid_Key), intent(in)                    ::  cgg_key
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    type(T_Grid_Key), intent(in)                    ::  lli_key
    type(TllGridInfo), intent(in)                   ::  lli
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Cgg_LL_Init_Interpol'
    
    ! --- local ----------------------------------
    
    integer              ::  max_np
    integer              ::  ilon, ilat
    real                 ::  lon, lat

    ! --- begin ----------------------------------

    ! fill keys:
    call Grid_Convertor_Init( self, howto, cgg_key, lli_key, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! save dimension
    self%nlon = lli%nlon
    self%nlat = lli%nlat
    self%ni   = cgg%nlon
    self%nj   = cgg%nlat
    
    ! maximum number of surrounding points for bilinear interpolation:
    max_np = 4
    
    ! allocate arrays
    allocate( self%ncov(lli%nlon,lli%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ii  (lli%nlon,lli%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj  (lli%nlon,lli%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%frac(lli%nlon,lli%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! zero by default:
    self%ncov = 0
    self%ii   = 0
    self%jj   = 0
    self%frac = 0.0
    
    ! loop over target cells:
    do ilat = 1, lli%nlat
      do ilon = 1, lli%nlon
    
        ! target location:
        lon = lli%lon_deg(ilon)
        lat = lli%lat_deg(ilat)
        
        ! get weights:
        call cgg%GetInterpol( lon, lat, self%ncov(ilon,ilat), &
                               self%ii(ilon,ilat,:), self%jj(ilon,ilat,:), &
                               self%frac(ilon,ilat,:), status )
        IF_NOT_OK_RETURN(status=1)

      end do  ! lli i
    end do  ! lli j
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Cgg_LL_Init_Interpol
  

  ! *

  
  subroutine Grid_Convertor_Cgg_LL_Init_Nearest( self, howto, cgg_key, cgg, lli_key, lli, status )
  
    use grid_type_ll , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Cgg_LL), intent(out)     ::  self
    character(len=*), intent(in)                    ::  howto
    type(T_Grid_Key), intent(in)                    ::  cgg_key
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    type(T_Grid_Key), intent(in)                    ::  lli_key
    type(TllGridInfo), intent(in)                   ::  lli
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Cgg_LL_Init_Nearest'
    
    ! --- local ----------------------------------
    
    integer              ::  max_np
    integer              ::  ilon, ilat
    real                 ::  lon, lat
    integer              ::  np
    integer              ::  ii(4), jj(4)
    real                 ::  ff(4)
    integer              ::  ip
    integer              ::  k

    ! --- begin ----------------------------------

    ! fill keys:
    call Grid_Convertor_Init( self, howto, cgg_key, lli_key, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! save dimension
    self%nlon = lli%nlon
    self%nlat = lli%nlat
    self%ni   = cgg%nlon
    self%nj   = cgg%nlat
    
    ! maximum number of surrounding points for bilinear interpolation:
    max_np = 1
    
    ! allocate arrays
    allocate( self%ncov(lli%nlon,lli%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ii  (lli%nlon,lli%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj  (lli%nlon,lli%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%frac(lli%nlon,lli%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! zero by default:
    self%ncov = 0
    self%ii   = 0
    self%jj   = 0
    self%frac = 0.0
    
    ! loop over target cells:
    do ilat = 1, lli%nlat
      do ilon = 1, lli%nlon
    
        ! target location:
        lon = lli%lon_deg(ilon)
        lat = lli%lat_deg(ilat)
        
        ! get weights for bi-linear interpoltion:
        call cgg%GetInterpol( lon, lat, np, ii, jj, ff, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! point with heights weight is most nearby:
        ip = 1
        do k = 2, np
          if ( ff(k) > ff(ip) ) ip = k
        end do
            
        ! single point:
        self%ncov(ilon,ilat) = 1
        ! copy location:
        self%ii  (ilon,ilat,1) = ii(ip)
        self%jj  (ilon,ilat,1) = jj(ip)
        ! full weight:
        self%frac(ilon,ilat,1) = 1.0
        
      end do  ! lli i
    end do  ! lli j
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Cgg_LL_Init_Nearest
  

  ! *
  

  subroutine Grid_Convertor_Cgg_LL_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Cgg_LL), intent(inout)     ::  self
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Cgg_LL_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    call Grid_Convertor_Done( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! deallocate arrays
    deallocate( self%ncov )
    deallocate( self%ii   )
    deallocate( self%jj   )
    deallocate( self%frac )
    if ( allocated( self%A_cg) ) deallocate( self%A_cg )
    if ( allocated( self%A_ll) ) deallocate( self%A_ll )
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Cgg_LL_Done
  
  
  ! ***


  subroutine Grid_Convertor_Cgg_LL_AreaAver_2d( self, cg, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Cgg_LL), intent(inout)     ::  self
    real, intent(in)                                  ::  cg(:,:)      ! (ni,nj)
    real, intent(out)                                 ::  ll(:,:)      ! (nlon,nlat)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Cgg_LL_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer            ::  ilon, ilat
    integer            ::  icov
    integer            ::  i, j
    integer            ::  ncov, indx
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (size(ll,1) /= self%nlon) .or. (size(ll,2) /= self%nlat) ) then
      write (gol,'("shape of target grid (",i0,",",i0,") differs from definition (",i0,",",i0,")")') &
               size(ll,1), size(ll,2), self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init output:
    ll = 0.0
    ! loop over target cells:
    do ilat = 1, self%nlat
      do ilon = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(ilon,ilat)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do icov = 1, ncov
            ! source index:
            i = self%ii(ilon,ilat,icov)
            j = self%jj(ilon,ilat,icov)
            ! add contribution:
            ll(ilon,ilat) = ll(ilon,ilat) + cg(i,j) * self%A_cg(i,j) * self%frac(ilon,ilat,icov)
          end do
        end if
        ! area average:
        ll(ilon,ilat) = ll(ilon,ilat) / self%A_ll(ilon,ilat)
      end do ! cells i
    end do ! cells j
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Cgg_LL_AreaAver_2d


  ! *
  
  
  subroutine Grid_Convertor_Cgg_LL_AreaAver_3d( self, cg, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Cgg_LL), intent(inout)     ::  self
    real, intent(in)                                  ::  cg(:,:,:)      ! (ni,nj,:)
    real, intent(out)                                 ::  ll(:,:,:)      ! (nlon,nlat,:)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Cgg_LL_AreaAver_3d'
    
    ! --- local ----------------------------------
    
    integer            ::  ilon, ilat
    integer            ::  icov
    integer            ::  i, j
    integer            ::  ncov, indx
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (size(ll,1) /= self%nlon) .or. (size(ll,2) /= self%nlat) ) then
      write (gol,'("shape of target grid (",i0,",",i0,") differs from definition (",i0,",",i0,")")') &
               size(ll,1), size(ll,2), self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( size(ll,3) /= size(cg,3) ) then
      write (gol,'("number of levels in target grid (",i0,") differs from source (",i0,")")') &
               size(ll,3), size(cg,3); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init output:
    ll = 0.0
    ! loop over target cells:
    do ilat = 1, self%nlat
      do ilon = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(ilon,ilat)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do icov = 1, ncov
            ! source index:
            i = self%ii(ilon,ilat,icov)
            j = self%jj(ilon,ilat,icov)
            ! add contribution:
            ll(ilon,ilat,:) = ll(ilon,ilat,:) + cg(i,j,:) * self%A_cg(i,j) * self%frac(ilon,ilat,icov)
          end do
        end if
        ! area average:
        ll(ilon,ilat,:) = ll(ilon,ilat,:) / self%A_ll(ilon,ilat)
      end do ! cells i
    end do ! cells j
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Cgg_LL_AreaAver_3d


  ! ***
  

  subroutine Grid_Convertor_Cgg_LL_Interpol_2d( self, cg, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Cgg_LL), intent(inout)     ::  self
    real, intent(in)                                  ::  cg(:,:)      ! (i,j)
    real, intent(out)                                 ::  ll(:,:)      ! (nlon,nlat)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Cgg_LL_Interpol_2d'
    
    ! --- local ----------------------------------
    
    integer            ::  ilon, ilat
    integer            ::  ncov, icov
    integer            ::  i, j
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (size(ll,1) /= self%nlon) .or. (size(ll,2) /= self%nlat) ) then
      write (gol,'("shape of target grid (",i0,",",i0,") differs from definition (",i0,",",i0,")")') &
               size(ll,1), size(ll,2), self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init output:
    ll = 0.0
    ! loop over target cells:
    do ilat = 1, self%nlat
      do ilon = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(ilon,ilat)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do icov = 1, ncov
            ! source indeces:
            i = self%ii(ilon,ilat,icov)
            j = self%jj(ilon,ilat,icov)
            ! add contribution:
            ll(ilon,ilat) = ll(ilon,ilat) + cg(i,j) * self%frac(ilon,ilat,icov)
          end do ! source cells
        end if ! any source cells?
      end do ! lons
    end do ! lats
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Cgg_LL_Interpol_2d


  ! *
  

  subroutine Grid_Convertor_Cgg_LL_Interpol_3d( self, cg, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Cgg_LL), intent(inout)     ::  self
    real, intent(in)                                  ::  cg(:,:,:)      ! (nx,ny,nz)
    real, intent(out)                                 ::  ll(:,:,:)      ! (nlon,nlat,nz)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Cgg_LL_Interpol_3d'
    
    ! --- local ----------------------------------
    
    integer            ::  ilon, ilat
    integer            ::  ncov, icov
    integer            ::  i, j
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (size(ll,1) /= self%nlon) .or. (size(ll,2) /= self%nlat) ) then
      write (gol,'("shape of target grid (",i0,",",i0,") differs from definition (",i0,",",i0,")")') &
               size(ll,1), size(ll,2), self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( size(ll,3) /= size(cg,3) ) then
      write (gol,'("number of levels in target grid (",i0,") differs from source (",i0,")")') &
               size(ll,3), size(cg,3); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init output:
    ll = 0.0
    ! loop over target cells:
    do ilat = 1, self%nlat
      do ilon = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(ilon,ilat)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do icov = 1, ncov
            ! source indeces:
            i = self%ii(ilon,ilat,icov)
            j = self%jj(ilon,ilat,icov)
            ! add contribution:
            ll(ilon,ilat,:) = ll(ilon,ilat,:) + cg(i,j,:) * self%frac(ilon,ilat,icov)
          end do ! source cells
        end if ! any source cells?
      end do ! lons
    end do ! lats
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Cgg_LL_Interpol_3d
  
  
  ! ********************************************************************
  ! ***
  ! *** ugg to ugg
  ! ***
  ! ********************************************************************

  !
  ! Initialize arrays.
  ! Optional nstep is used to extend the "from" arrays if they are too small.
  !

  subroutine Grid_Convertor_Ugg_InitMapping( self, nlon, nlat, status, nstep )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(inout)      ::  self
    integer, intent(in)                             ::  nlon, nlat
    integer, intent(out)                            ::  status

    integer, intent(in), optional                   ::  nstep

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_InitMapping'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store target shape:
    self%nlon = nlon
    self%nlat = nlat
    
    ! no sources yet:
    self%nfrom = 0

    ! step size for extending; 
    ! number of source probably scales with number of target cells:
    self%nstep = nlon*nlat
    ! reset from argument?
    if ( present(nstep) ) self%nstep = nstep
    
    ! init pointer arrays:
    nullify( self%ifrom )
    nullify( self%jfrom )
    nullify( self%wfrom )
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Ugg_InitMapping


  ! *
  
  
  subroutine Grid_Convertor_Ugg_AddMapping( self, ilon, ilat, i, j, w, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(inout)      ::  self
    integer, intent(in)                             ::  ilon, ilat
    integer, intent(in)                             ::  i, j
    real, intent(in)                                ::  w
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_AddMapping'
    
    ! --- local ----------------------------------
    
    integer             ::  n
    integer, pointer    ::  ifrom(:)
    integer, pointer    ::  jfrom(:)
    real, pointer       ::  wfrom(:)

    ! --- begin ----------------------------------
    
    ! new ?
    if ( self%nfrom == 0 ) then
      ! source ranges:
      allocate( self%froms(self%nlon,self%nlat,2), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! initial values:
      self%froms = 0
    end if
    
    ! extra storage needed?
    if ( (self%nfrom == 0) .or. (self%nfrom == size(self%ifrom)) ) then
      ! new size:
      n = self%nfrom + self%nstep
      ! new storage:
      allocate( ifrom(n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( jfrom(n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( wfrom(n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! any current content?
      if ( self%nfrom > 0 ) then
        ! copy:
        ifrom(1:self%nfrom) = self%ifrom
        jfrom(1:self%nfrom) = self%jfrom
        wfrom(1:self%nfrom) = self%wfrom
        ! clear previous storage:
        deallocate( self%ifrom, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( self%jfrom, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( self%wfrom, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if ! copy current storage
      ! reset pointers to new storage:
      self%ifrom => ifrom
      self%jfrom => jfrom
      self%wfrom => wfrom
    end if ! extend storage

    ! increase counter:
    self%nfrom = self%nfrom + 1
    ! set first index?
    if ( self%froms(ilon,ilat,1) == 0 ) then
      ! set start:
      self%froms(ilon,ilat,1) = self%nfrom
    end if
    ! set end index:
    self%froms(ilon,ilat,2) = self%nfrom
    ! store source indices and weight:
    self%ifrom(self%nfrom) = i
    self%jfrom(self%nfrom) = j
    self%wfrom(self%nfrom) = w
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Ugg_AddMapping


  ! *
  
  
  subroutine Grid_Convertor_Ugg_Init_AreaRemap( self, howto, ugg_key_from, ugg_from, &
                                                             ugg_key_to  , ugg_to  , status, debug )
  
    use grid_tools   , only : pi, ll_area_frac_deg
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
    use GO           , only : T_Polygon, T_PlotFile
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(out)        ::  self
    character(len=*), intent(in)                    ::  howto
    type(T_Grid_Key), intent(in)                    ::  ugg_key_from
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    type(T_Grid_Key), intent(in)                    ::  ugg_key_to
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    integer, intent(out)                            ::  status
    
    logical, intent(in), optional ::  debug

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_Init_AreaRemap'
    
    ! --- local ----------------------------------

    !logical             ::  dbg, dbg2
    character(len=256)  ::  gridfromto
    
    integer              ::  max_ncov_lon, max_ncov_lat, max_ncov
    integer              ::  ilon, ilat
    integer              ::  i, i1, i2, idir
    integer              ::  j, j1, j2, jdir
    real                 ::  west1, east1, south1, north1
    real                 ::  west2, east2, south2, north2
    real                 ::  frac
    integer              ::  ncov

    real(8)              ::  bbox(4)
    integer              ::  ibox(4)
    type(T_Polygon)      ::  pg
    real(8)              ::  pg_area
    real                 ::  w
    type(T_PlotFile)     ::  pf

    real(8)              ::  from_area

    ! --- begin ----------------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug

    ! init base class:
    call Grid_Convertor_Init( self, howto, ugg_key_from, ugg_key_to, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init mapping to target array:
    call self%InitMapping( ugg_to%nlon, ugg_to%nlat, status )
    IF_NOT_OK_RETURN(status=1)

    ! copy source dimensions:
    self%ni   = ugg_from%nlon
    self%nj   = ugg_from%nlat
    
    ! storage for cell area's in source grid:
    allocate( self%A_ug_from(ugg_from%nlon,ugg_from%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    self%A_ug_from = ugg_from%area  ! m2

    ! storage for cell area destination grid:
    allocate( self%A_ug_to(ugg_to%nlon,ugg_to%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    self%A_ug_to = ugg_to%area ! m2
    
    ! define in/out gridding
    gridfromto = 'from '//trim(ugg_from%type)//' to '//trim(ugg_to%type)
    ! switch:
    select case ( trim(gridfromto) )

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'from cartesian-regular to cartesian-regular', &
             'from cartesian-regular to cartesian'        , &
             'from cartesian to cartesian-regular'        , &
             'from cartesian to cartesian'                  )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! estimate maximum number of orig cells covering a target cell in each direction:
        max_ncov_lon = ceiling( ugg_to%dlon / minval(abs(ugg_from%longitude_bnds_1d(2,:)-ugg_from%longitude_bnds_1d(1,:))) )
        max_ncov_lat = ceiling( ugg_to%dlat / minval(abs(ugg_from%latitude_bnds_1d (2,:)-ugg_from%latitude_bnds_1d (1,:))) )
        ! combined, allow extra at both sides:
        max_ncov = (max_ncov_lon+2) * (max_ncov_lat+2)

        ! allocate arrays
        allocate( self%ncov(ugg_to%nlon,ugg_to%nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%ii  (ugg_to%nlon,ugg_to%nlat,max_ncov), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%jj  (ugg_to%nlon,ugg_to%nlat,max_ncov), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%frac(ugg_to%nlon,ugg_to%nlat,max_ncov), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! zero by default:
        self%ncov = 0
        self%ii   = 0
        self%jj   = 0
        self%frac = 0.0

        ! loop over target cells:
        do ilat = 1, ugg_to%nlat
          do ilon = 1, ugg_to%nlon

            ! extract boundaries of target cell:
            west2  = ugg_to%longitude_bnds_1d(1,ilon)        
            east2  = ugg_to%longitude_bnds_1d(2,ilon)
            if ( east2 <= 0.0 ) then
              west2 = west2 + 360.0    ! (0,360)
              east2 = east2 + 360.0    ! (0,360)
            end if
            south2 = ugg_to%latitude_bnds_1d(1,ilat)
            north2 = ugg_to%latitude_bnds_1d(2,ilat)

            ! search ugg_from rows including west and east ll lat:
            call ugg_from%GetLonIndex( west2, i1, status )
            IF_NOT_OK_RETURN(status=1)
            call ugg_from%GetLonIndex( east2, i2, status )
            IF_NOT_OK_RETURN(status=1)
            ! direction:
            idir = sign( 1, i2 - i1 )

            ! search ugg_from rows including north and south ll lat:
            call ugg_from%GetLatIndex( south2, j1, status )
            IF_NOT_OK_RETURN(status=1)
            call ugg_from%GetLatIndex( north2, j2, status )
            IF_NOT_OK_RETURN(status=1)
            ! direction:
            jdir = sign( 1, j2 - j1 )

            ! loop over source cells:
            do j = j1, j2, jdir
              do i = i1, i2, idir

                ! bounds:
                west1  = minval(ugg_from%longitude_bnds_1d(:,i))
                east1  = maxval(ugg_from%longitude_bnds_1d(:,i))
                south1 = minval(ugg_from%latitude_bnds_1d (:,j))
                north1 = maxval(ugg_from%latitude_bnds_1d (:,j))

                ! shift if gg cell is right from [0,360]
                if ( west1 > east2 ) then
                  west1 = west1 - 360.0    ! (0,360)
                  east1 = east1 - 360.0    ! (0,360)
                end if

                ! determine covarage fraction:
                call ll_area_frac_deg( west1, east1, south1, north1, &
                                       west2, east2, south2, north2, &
                                       frac, status )
                IF_NOT_OK_RETURN(status=1)

                ! fill fraction:
                if ( abs(frac) < 1.0e-4 ) then
                  ! almost no coverage ...

                else if ( (frac > 0.0) .and. (frac <= 1.0) ) then
                  ! coverage within (0,1]; increase counter:
                  ncov = self%ncov(ilon,ilat) + 1
                  ! check ...
                  if ( ncov > max_ncov ) then
                    write (gol,'("number of covering cells (",i0,") exceeds estimated maximum (",i0,")")') ncov, max_ncov; call goErr
                    write (gol,'("  target resolution        : ",2es12.4)') ugg_to%dlon, ugg_to%dlat; call goErr
                    write (gol,'("  min. source resolution   : ",2es12.4)') &
                                           minval(abs(ugg_from%longitude_bnds_1d(2,:)-ugg_from%longitude_bnds_1d(1,:))), &
                                           minval(abs(ugg_from%latitude_bnds_1d (2,:)-ugg_from%latitude_bnds_1d (1,:))); call goErr
                    write (gol,'("  max. coverage            : ",2i6)') max_ncov_lon, max_ncov_lat; call goErr
                    TRACEBACK; status=1; return
                  end if
                  ! store:
                  self%ncov(ilon,ilat     ) = ncov
                  self%ii  (ilon,ilat,ncov) = i
                  self%jj  (ilon,ilat,ncov) = j
                  self%frac(ilon,ilat,ncov) = frac

                else if ( abs(frac-1.0) < 1.0e-4 ) then
                  ! 100% coverage, but small rounding error; increase counter:
                  ncov = self%ncov(ilon,ilat) + 1
                  ! check ...
                  if ( ncov > max_ncov ) then
                    write (gol,'("number of covering cells (",i0,") exceeds estimated maximum (",i0,")")') ncov, max_ncov; call goErr
                    write (gol,'("  target resolution        : ",2es12.4)') ugg_to%dlon, ugg_to%dlat; call goErr
                    write (gol,'("  min. source resolution   : ",2es12.4)') &
                                           minval(abs(ugg_from%longitude_bnds_1d(2,:)-ugg_from%longitude_bnds_1d(1,:))), &
                                           minval(abs(ugg_from%latitude_bnds_1d (2,:)-ugg_from%latitude_bnds_1d (1,:))); call goErr
                    write (gol,'("  max. coverage            : ",2i6)') max_ncov_lon, max_ncov_lat; call goErr
                    TRACEBACK; status=1; return
                  end if
                  ! store:
                  self%ncov(ilon,ilat     ) = ncov
                  self%ii  (ilon,ilat,ncov) = i
                  self%jj  (ilon,ilat,ncov) = j
                  self%frac(ilon,ilat,ncov) = 1.0

                else 
                  write (gol,'("BUG - found strange coverage fraction : ",e15.6)') frac; call goErr
                  write (gol,'("BUG - source cell : ",4f12.6)') west1, east1, south1, north1; call goErr
                  write (gol,'("BUG - target cell : ",4f12.6)') west2, east2, south2, north2; call goErr
                  TRACEBACK; status=1; return
                end if

              end do  ! ugg_from i
            end do  ! ugg_from j

            ! check ...
            if ( self%ncov(ilon,ilat)  == 0 ) then
              write (gol,'("found target cell not covered by source grid")'); call goErr
              TRACEBACK; status=1; return
            end if

          end do  ! ugg_to ilon
        end do  ! ugg_to ilat
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! generic version ...
        
        !if (dbg) print *, '  ggg1 target ', ugg_to%nlon, ugg_to%nlat
        
        ! loop over target cells:
        do ilat = 1, ugg_to%nlat
          do ilon = 1, ugg_to%nlon
          
            ! get bounding box for target cell: (/west,east,south,north/)
            ! this is used to limit searching for source cells:
            call ugg_to%pg(ilon,ilat)%BoundBox( bbox, status )
            IF_NOT_OK_RETURN(status=1)
            
            ! get range of source cells that might include target cells:
            call ugg_from%GetIndexDomain( real(bbox), ibox, status )
            IF_NOT_OK_RETURN(status=1)
            
            !dbg2 = dbg .and. (ilat==1) .and. (ilon==1)
            !if (dbg2) then
            !  print *, '  ggg1 cell', ilon, ilat
            !  print *, '  ggg1 bbox = ', bbox
            !  print *, '  ggg1 ibox = ', ibox
            !end if
            
            ! init counter:
            ncov = 0
            ! loop over source cells:
            do j = ibox(3), ibox(4)
              do i = ibox(1), ibox(2)

                !! testing ...
                !if (dbg2 ) print *, '  ggg1 source ', i, j
                
                !! intersection, status<0 if no overlap:
                !call ugg_from%pg(i,j)%Intersection( ugg_to%pg(ilon,ilat), pg, status )!, debug=dbg2 )
                !!if (dbg2 ) print *, '  ggg1 intersection status ', status
                !if ( status < 0 ) then
                !  ! no overlap, next source cell:
                !  cycle
                !else if ( status > 0 ) then
                !  write (gol,'("from intersection between source (",i0,",",i0,") and target (",i0,",",i0,")")') &
                !           i, j, ilon, ilat; call goErr
                !  TRACEBACK; status=1; return
                !end if
                
                ! init polygon as copy of input polygon clipped by target:
                call pg%InitClippedPolygon( ugg_from%pg(i,j), ugg_to%pg(ilon,ilat), status )
                IF_NOT_OK_RETURN(status=1)
                ! filled ?
                if ( .not. pg%IsEmpty() ) then
                  ! increase counter:
                  ncov = ncov + 1
                  ! area of intersection:
                  call pg%LonLat_Area( pg_area, status )
                  if ( status /=0 ) then
                    call pf%Init( 'debug-area-remap.py', status )
                    IF_NOT_OK_RETURN(status=1)
                    call ugg_from%pg(i,j)%Plot( 'edges', pf, status, kwargs='color="red"' )
                    IF_NOT_OK_RETURN(status=1)
                    call ugg_to%pg(ilon,ilat)%Plot( 'edges', pf, status, kwargs='color="blue"' )
                    IF_NOT_OK_RETURN(status=1)
                    call pg%Plot( 'edges', pf, status, kwargs='color="green"' )
                    IF_NOT_OK_RETURN(status=1)
                    call pf%Done( status )
                    IF_NOT_OK_RETURN(status=1)
                    TRACEBACK; status=1; return
                  end if
                  IF_NOT_OK_RETURN(status=1)
                  ! skip if zero (might not have been trapped by Intersection):
                  if ( pg_area == 0.0 ) cycle
                  ! weight of contribution is area fraction:
                  w = real( pg_area / ugg_from%area(i,j) )
                  ! round values just above one ...
                  if ( (1.0 < w) .and. (w < 1.0+1.0e-3) ) then
                    w = 1.0
                  end if
                  ! check ...
                  if ( (pg_area <= 0.0) .or. (ugg_from%area(i,j) <=0) .or. &
                        (w <= 0.0) .or. (w > 1.0) ) then
                    !write (gol,'("area of source area")'); call goErr
                    !call ugg_from%pg(i,j)%LonLat_Area( from_area, status, debug=.true. )
                    !IF_NOT_OK_RETURN(status=1)
                    !write (gol,'("area of intersection area")'); call goErr
                    !call pg%LonLat_Area( pg_area, status, debug=.true. )
                    !IF_NOT_OK_RETURN(status=1)
                    write (gol,'("strange mapping weight : ",es12.4)') w; call goErr
                    write (gol,'("  intersection area    : ",es12.4)') pg_area; call goErr
                    write (gol,'("  source area          : ",es12.4)') ugg_from%area(i,j); call goErr
                    write (gol,'("mapping from source (",i0,",",i0,") to target (",i0,",",i0,")")') &
                         i, j, ilon, ilat; call goErr
                    ! ... plot ...
                    call pf%Init( 'debug.py', status )
                    IF_NOT_OK_RETURN(status=1)
                    call ugg_from%pg(i,j)%Plot( 'edges', pf, status, kwargs='color="red"' )
                    IF_NOT_OK_RETURN(status=1)
                    !call ugg_from%pg(i,j)%Plot( 'edge-arrows', pf, status, kwargs='color="red", head_width=0.01' )
                    !IF_NOT_OK_RETURN(status=1)
                    call ugg_to%pg(ilon,ilat)%Plot( 'edges', pf, status, kwargs='color="blue"' )
                    IF_NOT_OK_RETURN(status=1)
                    !call ugg_to%pg(ilon,ilat)%Plot( 'edge-arrows', pf, status, kwargs='color="blue", head_width=0.01' )
                    !IF_NOT_OK_RETURN(status=1)
                    call pg%Plot( 'edges', pf, status, kwargs='color="green"' )
                    IF_NOT_OK_RETURN(status=1)
                    !call pg%Plot( 'edge-arrows', pf, status, kwargs='color="green", head_width=0.01' )
                    !IF_NOT_OK_RETURN(status=1)
                    call pf%Done( status )
                    IF_NOT_OK_RETURN(status=1)
                    ! ..............
                    TRACEBACK; status=1; return
                  end if
                  ! add contribution:
                  call self%AddMapping( ilon, ilat, i, j, w, status )
                  IF_NOT_OK_RETURN(status=1)
                end if
                
                ! clear:
                call pg%Done( status )
                IF_NOT_OK_RETURN(status=1)
                
              end do ! i
            end do ! j
            
            ! check ...
            if ( ncov == 0 ) then
              write (gol,'("no source cells found to cover target")'); call goErr
              write (gol,'("  target cell  : ",2i4)') ilon, ilat; call goErr
              write (gol,'("  source box   : ",4i4)') ibox; call goErr
              call pf%Init( 'debug.py', status )
              IF_NOT_OK_RETURN(status=1)
              do j = ibox(3), ibox(4)
                do i = ibox(1), ibox(2)
                  call ugg_from%pg(i,j)%Plot( 'edges', pf, status, kwargs='color="red"' )
                  IF_NOT_OK_RETURN(status=1)
                  !! info ..
                  !write (gol,'("   source cell: ",2i4)') i,j; call goErr
                  !call ugg_from%pg(i,j)%Intersection( ugg_to%pg(ilon,ilat), pg, status, debug=.true. )
                  !write (gol,'("     intersection status: ",i0)') status; call goErr
                end do
              end do
              call ugg_to%pg(ilon,ilat)%Plot( 'edges', pf, status, kwargs='color="blue"' )
              IF_NOT_OK_RETURN(status=1)
              call pf%Done( status )
              IF_NOT_OK_RETURN(status=1)
              TRACEBACK; status=1; return
            end if

          end do  ! ugg_to ilon
        end do  ! ugg_to ilat

    end select

    ! ok
    status = 0

  end subroutine Grid_Convertor_Ugg_Init_AreaRemap
  

  ! *

  
  subroutine Grid_Convertor_Ugg_Init_Interpol( self, howto, ugg_key_from, ugg_from, ugg_key_to, ugg_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(out)        ::  self
    character(len=*), intent(in)                    ::  howto
    type(T_Grid_Key), intent(in)                    ::  ugg_key_from
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    type(T_Grid_Key), intent(in)                    ::  ugg_key_to
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_Init_Interpol'
    
    ! --- local ----------------------------------
    
    integer              ::  max_np
    integer              ::  ilon, ilat
    real                 ::  lon, lat

    ! --- begin ----------------------------------

    ! fill keys:
    call Grid_Convertor_Init( self, howto, ugg_key_from, ugg_key_to, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init mapping to target array:
    call self%InitMapping( ugg_to%nlon, ugg_to%nlat, status )
    IF_NOT_OK_RETURN(status=1)

    ! save dimension
    self%ni   = ugg_from%nlon
    self%nj   = ugg_from%nlat
    
    ! maximum number of surrounding points for bilinear interpolation:
    max_np = 4
    
    ! allocate arrays
    allocate( self%ncov(ugg_to%nlon,ugg_to%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ii  (ugg_to%nlon,ugg_to%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj  (ugg_to%nlon,ugg_to%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%frac(ugg_to%nlon,ugg_to%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! zero by default:
    self%ncov = 0
    self%ii   = 0
    self%jj   = 0
    self%frac = 0.0
    
    ! loop over target cells:
    do ilat = 1, ugg_to%nlat
      do ilon = 1, ugg_to%nlon
    
        ! target location:
        lon = ugg_to%longitude(ilon,ilat)
        lat = ugg_to%latitude(ilon,ilat)
        
        ! get weights:
        call ugg_from%GetInterpol( lon, lat, self%ncov(ilon,ilat), &
                                     self%ii(ilon,ilat,:), self%jj(ilon,ilat,:), &
                                     self%frac(ilon,ilat,:), status )
        IF_NOT_OK_RETURN(status=1)

      end do  ! lli i
    end do  ! lli j
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Ugg_Init_Interpol


  ! *

  
  subroutine Grid_Convertor_Ugg_Init_Nearest( self, howto, ugg_key_from, ugg_from, ugg_key_to, ugg_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(out)        ::  self
    character(len=*), intent(in)                    ::  howto
    type(T_Grid_Key), intent(in)                    ::  ugg_key_from
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    type(T_Grid_Key), intent(in)                    ::  ugg_key_to
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_Init_Nearest'
    
    ! --- local ----------------------------------
    
    integer              ::  max_np
    integer              ::  ilon, ilat
    real                 ::  lon, lat
    integer              ::  np
    integer              ::  ii(4), jj(4)
    real                 ::  ff(4)
    integer              ::  ip
    integer              ::  k

    ! --- begin ----------------------------------

    ! fill keys:
    call Grid_Convertor_Init( self, howto, ugg_key_from, ugg_key_to, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init mapping to target array:
    call self%InitMapping( ugg_to%nlon, ugg_to%nlat, status )
    IF_NOT_OK_RETURN(status=1)

    ! save dimension
    self%ni   = ugg_from%nlon
    self%nj   = ugg_from%nlat
    
    ! maximum number of surrounding points for bilinear interpolation:
    max_np = 1
    
    ! allocate arrays
    allocate( self%ncov(ugg_to%nlon,ugg_to%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ii  (ugg_to%nlon,ugg_to%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj  (ugg_to%nlon,ugg_to%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%frac(ugg_to%nlon,ugg_to%nlat,max_np), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! zero by default:
    self%ncov = 0
    self%ii   = 0
    self%jj   = 0
    self%frac = 0.0
    
    ! loop over target cells:
    do ilat = 1, ugg_to%nlat
      do ilon = 1, ugg_to%nlon
    
        ! target location:
        lon = ugg_to%longitude(ilon,ilat)
        lat = ugg_to%latitude(ilon,ilat)
        
        ! get weights for bi-linear interpoltion:
        call ugg_from%GetInterpol( lon, lat, np, ii, jj, ff, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! point with heights weight is most nearby:
        ip = 1
        do k = 2, np
          if ( ff(k) > ff(ip) ) ip = k
        end do
            
        ! single point:
        self%ncov(ilon,ilat) = 1
        ! copy location:
        self%ii  (ilon,ilat,1) = ii(ip)
        self%jj  (ilon,ilat,1) = jj(ip)
        ! full weight:
        self%frac(ilon,ilat,1) = 1.0
        
      end do  ! lli i
    end do  ! lli j
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Ugg_Init_Nearest


  ! * 

  
!  subroutine Grid_Convertor_Ugg_Init_IndexFractions( self, howto, ugg_key_from, ugg_from, ugg_key_to, ugg_to, status, &
!                                                             ug_from, ug_to )
!  
!    use C3PO_Grid_Ugg, only : T_Grid_Ugg
!    use grid_tools   , only : pi, ll_area_frac_deg
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_Grid_Convertor_Ugg), intent(out)        ::  self
!    character(len=*), intent(in)                    ::  howto
!    type(T_Grid_Key), intent(in)                    ::  ugg_key_from
!    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
!    type(T_Grid_Key), intent(in)                    ::  ugg_key_to
!    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
!    integer, intent(out)                            ::  status
!
!    integer, intent(in), optional                   ::  ug_from(:,:) ! (nlon_from,nlat_from)
!    real, intent(out), optional                     ::  ug_to(:,:,:) ! (nlon_to,nlat_to,1:max(ug_from))
!    
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_Init_IndexFractions'
!    
!    ! --- local ----------------------------------
!    
!    integer              ::  max_ncov_lon, max_ncov_lat, max_ncov
!    integer              ::  ilon, ilat, k
!    integer              ::  i, i1, i2, idir
!    integer              ::  j, j1, j2, jdir
!    real                 ::  west1, east1, south1, north1
!    real                 ::  west2, east2, south2, north2
!    real                 ::  frac
!    integer              ::  ncov
!    integer              ::  nlev
!
!    ! --- begin ----------------------------------
!
!    ! fill keys:
!    call Grid_Convertor_Init( self, howto, ugg_key_from, ugg_key_to, status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! init mapping to target array:
!    call self%InitMapping( ugg_to%nlon, ugg_to%nlat, status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! save dimension
!    self%ni   = ugg_from%nlon
!    self%nj   = ugg_from%nlat
!
!    ! check ...
!    if ( .not. allocated(ugg_from%longitude_bnds_1d) ) then
!      write (gol,'("bnds_1d not allocated in source grid; not cartesian?")'); call goErr
!      TRACEBACK; status=1; return
!    end if
!
!    ! check ...
!    if ( .not. allocated(ugg_to%longitude_bnds_1d) ) then
!      write (gol,'("bnds_1d not allocated in target grid; not cartesian?")'); call goErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! estimate maximum number of orig cells covering a target cell:
!    max_ncov_lon = ceiling( ugg_to%dlon / minval(abs(ugg_from%longitude_bnds_1d(2,:)-ugg_from%longitude_bnds_1d(1,:))) )
!    max_ncov_lat = ceiling( ugg_to%dlat / minval(abs(ugg_from%latitude_bnds_1d (2,:)-ugg_from%latitude_bnds_1d (1,:))) )
!    max_ncov = (max_ncov_lon+1) * (max_ncov_lat+1)
!    
!    if ( .not. ( present( ug_from ) .and. present( ug_to) ) ) then
!      ! save convertors, for further use allocate arrays
!      allocate( self%ncov(ugg_to%nlon,ugg_to%nlat), stat=status )
!      IF_NOT_OK_RETURN(status=1)
!      allocate( self%ii  (ugg_to%nlon,ugg_to%nlat,max_ncov), stat=status )
!      IF_NOT_OK_RETURN(status=1)
!      allocate( self%jj  (ugg_to%nlon,ugg_to%nlat,max_ncov), stat=status )
!      IF_NOT_OK_RETURN(status=1)
!      allocate( self%frac(ugg_to%nlon,ugg_to%nlat,max_ncov), stat=status )
!      IF_NOT_OK_RETURN(status=1)
!    
!      ! zero by default:
!      self%ncov = 0
!      self%ii   = 0
!      self%jj   = 0
!      self%frac = 0.0
!    else 
!      ncov = 0
!    end if
!    
!    ! loop over target cells:
!    do ilat = 1, ugg_to%nlat
!      do ilon = 1, ugg_to%nlon
!     
!        ! extract boundaries of target cell:
!        west2  = ugg_to%longitude_bnds_1d(1,ilon)        
!        east2  = ugg_to%longitude_bnds_1d(2,ilon)
!        if ( east2 < 0.0 ) then
!          west2 = west2 + 360.0    ! (0,360)
!          east2 = east2 + 360.0    ! (0,360)
!        end if
!        south2 = ugg_to%latitude_bnds_1d(1,ilat)
!        north2 = ugg_to%latitude_bnds_1d(2,ilat)
!
!        ! search ugg_from rows including west and east ll lat:
!        call ugg_from%GetLonIndex( west2, i1, status )
!        IF_NOT_OK_RETURN(status=1)
!        call ugg_from%GetLonIndex( east2, i2, status )
!        IF_NOT_OK_RETURN(status=1)
!        ! direction:
!        idir = sign( 1, i2 - i1 )
!
!        ! search ugg_from rows including north and south ll lat:
!        call ugg_from%GetLatIndex( south2, j1, status )
!        IF_NOT_OK_RETURN(status=1)
!        call ugg_from%GetLatIndex( north2, j2, status )
!        IF_NOT_OK_RETURN(status=1)
!        ! direction:
!        jdir = sign( 1, j2 - j1 )
!        
!        ! loop over source cells:
!        do j = j1, j2, jdir
!          do i = i1, i2, idir
!        
!            ! bounds:
!            west1  = minval(ugg_from%longitude_bnds_1d(:,i))
!            east1  = maxval(ugg_from%longitude_bnds_1d(:,i))
!            south1 = minval(ugg_from%latitude_bnds_1d (:,j))
!            north1 = maxval(ugg_from%latitude_bnds_1d (:,j))
!
!            ! shift if gg cell is right from [0,360]
!            if ( west1 > east2 ) then
!              west1 = west1 - 360.0    ! (0,360)
!              east1 = east1 - 360.0    ! (0,360)
!            end if
!
!            ! Fraction of target cell covered by source
!            call ll_area_frac_deg( west2, east2, south2, north2, &
!                                   west1, east1, south1, north1, &
!                                   frac, status )
!            IF_NOT_OK_RETURN(status=1)
!            
!            if ( .not. ( present( ug_from ) .and. present( ug_to) ) ) then
!              ! fill fractions and indices:
!              if ( frac > 0.0 .and. frac <= 1.0 ) then
!                ! add contribution:
!                ncov = self%ncov(ilon,ilat) + 1
!                self%ncov(ilon,ilat     ) = ncov
!                self%ii  (ilon,ilat,ncov) = i
!                self%jj  (ilon,ilat,ncov) = j
!                self%frac(ilon,ilat,ncov) = frac
!
!              else if ( abs(frac) < 1.0e-4 ) then
!                ! almost no coverage ...
!
!              else if ( abs(frac-1.0) < 1.0e-4 ) then
!                ! 100% coverage, but small rounding error
!                ncov = self%ncov(ilon,ilat) + 1
!                self%ncov(ilon,ilat     ) = ncov
!                self%ii  (ilon,ilat,ncov) = i
!                self%jj  (ilon,ilat,ncov) = j
!                self%frac(ilon,ilat,ncov) = 1.0
!
!              else 
!                write (gol,'("BUG - found strange coverage fraction : ",e15.6)') frac; call goErr
!                write (gol,'("BUG - source cell : ",4f12.6)') west1, east1, south1, north1; call goErr
!                write (gol,'("BUG - target cell : ",4f12.6)') west2, east2, south2, north2; call goErr
!                TRACEBACK; status=1; return
!              end if ! fraction check
!            
!            else 
!              ! calculate now:
!              if ( frac > 0.0 .and. frac <= 1.0 ) then
!                k = ug_from(i,j) 
!                ug_to(ilon,ilat,k) = ug_to(ilon,ilat,k) + ugg_to%area(ilon,ilat) * frac
!                ncov = ncov + 1
!              else if ( abs(frac) < 1.0e-4 ) then
!                ! almost no coverage ...
!
!              else if ( abs(frac-1.0) < 1.0e-4 ) then
!                ! 100% coverage, but small rounding error
!                frac = 1.0
!                k = ug_from(i,j) 
!                ug_to(ilon,ilat,k) = ug_to(ilon,ilat,k) + ugg_to%area(ilon,ilat) * frac
!                ncov = ncov + 1
!              else 
!                write (gol,'("BUG - found strange coverage fraction : ",e15.6)') frac; call goErr
!                write (gol,'("BUG - source cell : ",4f12.6)') west1, east1, south1, north1; call goErr
!                write (gol,'("BUG - target cell : ",4f12.6)') west2, east2, south2, north2; call goErr
!                TRACEBACK; status=1; return
!              end if ! fraction check
!            
!            end if ! present ug_from,ug_to
!            
!          end do  ! ugg_from i
!        end do  ! ugg_from j
!        
!        ! check ...
!        if ( .not. ( present( ug_from ) .and. present( ug_to) ) ) then
!          if ( self%ncov(ilon,ilat)  == 0 ) then
!            write (gol,'("found target cell not covered by source grid")'); call goErr
!            TRACEBACK; status=1; return
!          end if
!        else
!          if ( ncov == 0 ) then
!            write (gol,'("found target cell not covered by source grid")'); call goErr
!            TRACEBACK; status=1; return
!          end if
!        end if
!
!      end do  ! ugg_to ilon
!    end do  ! ugg_to ilat
!    
!    ! area average
!    if ( present( ug_from ) .and. present( ug_to) )  then
!      nlev = size(ug_to,dim=3)
!      do k = 1, nlev
!        ug_to(:,:,k) = ug_to(:,:,k) / ugg_to%area(:,:)
!      end do
!    end if
!    
!    ! storage for cell area in lon/lat grid:
!    allocate( self%A_ug_to(ugg_to%nlon,ugg_to%nlat), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! copy:
!    self%A_ug_to = ugg_to%area ! m2
!    
!    ! ok
!    status = 0
!
!  end subroutine Grid_Convertor_Ugg_Init_IndexFractions

  
  ! *
  
  
  subroutine Grid_Convertor_Ugg_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(inout)     ::  self
    integer, intent(out)                           ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    call Grid_Convertor_Done( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! clear area's?
    if ( allocated(self%A_ug_from) ) then
      deallocate( self%A_ug_from, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%A_ug_to) ) then
      deallocate( self%A_ug_to, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! deallocate arrays
    if ( allocated(self%ncov) ) then
      deallocate( self%ncov, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%ii  , stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%jj  , stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%frac, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! clear generic:
    if ( allocated(self%froms) ) then
      deallocate( self%froms, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%ifrom, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%jfrom, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%wfrom, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0

  end subroutine Grid_Convertor_Ugg_Done
  

  ! ***


  subroutine Grid_Convertor_Ugg_AreaAver_2d( self, ug_from, ug_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(inout)        ::  self
    real, intent(in)                                  ::  ug_from(:,:)      ! (ni,nj)
    real, intent(out)                                 ::  ug_to(:,:)      ! (nlon,nlat)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer            ::  ilon, ilat
    integer            ::  icov
    integer            ::  i, j
    integer            ::  ncov, indx
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (size(ug_to,1) /= self%nlon) .or. (size(ug_to,2) /= self%nlat) ) then
      write (gol,'("shape of target grid (",i0,",",i0,") differs from definition (",i0,",",i0,")")') &
               size(ug_to,1), size(ug_to,2), self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! switch:
    if ( self%nfrom > 0 ) then
    
      ! init output:
      ug_to = 0.0
      ! loop over target cells:
      do ilat = 1, self%nlat
        do ilon = 1, self%nlon
          ! any sources?
          if ( self%froms(ilon,ilat,1) > 0 ) then
            ! loop over source cells:
            do icov = self%froms(ilon,ilat,1), self%froms(ilon,ilat,2)
              ! source index:
              i = self%ifrom(icov)
              j = self%jfrom(icov)
              ! add contribution:
              ug_to(ilon,ilat) = ug_to(ilon,ilat) + ug_from(i,j) * self%A_ug_from(i,j) * self%wfrom(icov)
            end do
          end if
          ! area average:
          ug_to(ilon,ilat) = ug_to(ilon,ilat) / self%A_ug_to(ilon,ilat)
        end do ! cells i
      end do ! cells j
      
    else

      ! init output:
      ug_to = 0.0
      ! loop over target cells:
      do ilat = 1, self%nlat
        do ilon = 1, self%nlon
          ! number of source cells:
          ncov = self%ncov(ilon,ilat)
          ! any source cells ?
          if ( ncov > 0 ) then
            ! loop over source cells:
            do icov = 1, ncov
              ! source index:
              i = self%ii(ilon,ilat,icov)
              j = self%jj(ilon,ilat,icov)
              ! add contribution:
              ug_to(ilon,ilat) = ug_to(ilon,ilat) + ug_from(i,j) * self%A_ug_from(i,j) * self%frac(ilon,ilat,icov)
            end do
          end if
          ! area average:
          ug_to(ilon,ilat) = ug_to(ilon,ilat) / self%A_ug_to(ilon,ilat)
        end do ! cells i
      end do ! cells j

    end if
    
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Ugg_AreaAver_2d


  ! *
  
  
  subroutine Grid_Convertor_Ugg_AreaAver_3d( self, ug_from, ug_to, status, debug )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(inout)        ::  self
    real, intent(in)                                  ::  ug_from(:,:,:)      ! (ni,nj,:)
    real, intent(out)                                 ::  ug_to(:,:,:)      ! (nlon,nlat,:)
    integer, intent(out)                              ::  status
    
    logical, intent(in), optional  ::  debug

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_AreaAver_3d'
    
    ! --- local ----------------------------------
    
    integer            ::  ilon, ilat
    integer            ::  icov
    integer            ::  i, j
    integer            ::  ncov, indx
    
    !logical :: dbg
    
    ! --- begin ----------------------------------
    
    !dbg = .false.
    !if ( present(debug) ) dbg = debug
    
    ! check ...
    if ( (size(ug_to,1) /= self%nlon) .or. (size(ug_to,2) /= self%nlat) ) then
      write (gol,'("shape of target grid (",i0,",",i0,") differs from definition (",i0,",",i0,")")') &
               size(ug_to,1), size(ug_to,2), self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( size(ug_to,3) /= size(ug_from,3) ) then
      write (gol,'("number of levels in target grid (",i0,") differs from source (",i0,")")') &
               size(ug_to,3), size(ug_from,3); call goErr
      TRACEBACK; status=1; return
    end if
    
    !! testing ..
    !if (dbg) then
    !  print *, 'aaa1 target grid ', self%nlon, self%nlat
    !  print *, 'aaa1 nfrom = ', self%nfrom
    !end if
    
    ! switch:
    if ( self%nfrom > 0 ) then
    
      ! init output:
      ug_to = 0.0
      ! loop over target cells:
      do ilat = 1, self%nlat
        do ilon = 1, self%nlon
          !! testing ...
          !if ( ilon == 1 .and. ilat == 1 ) then
          !  print *, 'aaa1 cell ', ilon, ilat
          !  print *, 'aaa1 froms ', self%froms(ilon,ilat,:)
          !  if ( self%froms(ilon,ilat,1) > 0 ) then
          !    do icov = self%froms(ilon,ilat,1), self%froms(ilon,ilat,2)
          !      print *, 'aaa1 from ', icov, self%ifrom(icov), self%jfrom(icov), self%wfrom(icov)
          !    end do
          !  endif
          !  print *, 'aaa1 area ', self%A_ug_to(ilon,ilat)
          !end if
          ! any sources?
          if ( self%froms(ilon,ilat,1) > 0 ) then
            ! loop over source cells:
            do icov = self%froms(ilon,ilat,1), self%froms(ilon,ilat,2)
              ! source index:
              i = self%ifrom(icov)
              j = self%jfrom(icov)
              ! add contribution:
              ug_to(ilon,ilat,:) = ug_to(ilon,ilat,:) + ug_from(i,j,:) * self%A_ug_from(i,j) * self%wfrom(icov)
            end do
          end if
          ! area average:
          ug_to(ilon,ilat,:) = ug_to(ilon,ilat,:) / self%A_ug_to(ilon,ilat)
        end do ! cells i
      end do ! cells j
      
    else

      ! init output:
      ug_to = 0.0
      ! loop over target cells:
      do ilat = 1, self%nlat
        do ilon = 1, self%nlon
          ! number of source cells:
          ncov = self%ncov(ilon,ilat)
          ! any source cells ?
          if ( ncov > 0 ) then
            ! loop over source cells:
            do icov = 1, ncov
              ! source index:
              i = self%ii(ilon,ilat,icov)
              j = self%jj(ilon,ilat,icov)
              ! add contribution:
              ug_to(ilon,ilat,:) = ug_to(ilon,ilat,:) + ug_from(i,j,:) * self%A_ug_from(i,j) * self%frac(ilon,ilat,icov)
            end do
          end if
          ! area average:
          ug_to(ilon,ilat,:) = ug_to(ilon,ilat,:) / self%A_ug_to(ilon,ilat)
        end do ! cells i
      end do ! cells j
      
    end if

    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Ugg_AreaAver_3d


  ! ***
  

  subroutine Grid_Convertor_Ugg_Interpol_2d( self, ug_from, ug_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(inout)        ::  self
    real, intent(in)                                  ::  ug_from(:,:)      ! (i,j)
    real, intent(out)                                 ::  ug_to(:,:)      ! (nlon,nlat)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_Interpol_2d'
    
    ! --- local ----------------------------------
    
    integer            ::  ilon, ilat
    integer            ::  ncov, icov
    integer            ::  i, j
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (size(ug_to,1) /= self%nlon) .or. (size(ug_to,2) /= self%nlat) ) then
      write (gol,'("shape of target grid (",i0,",",i0,") differs from definition (",i0,",",i0,")")') &
               size(ug_to,1), size(ug_to,2), self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! switch:
    if ( self%nfrom > 0 ) then
      write (gol,'("not implemented yet for froms arrays")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init output:
    ug_to = 0.0
    ! loop over target cells:
    do ilat = 1, self%nlat
      do ilon = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(ilon,ilat)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do icov = 1, ncov
            ! source indeces:
            i = self%ii(ilon,ilat,icov)
            j = self%jj(ilon,ilat,icov)
            ! add contribution:
            ug_to(ilon,ilat) = ug_to(ilon,ilat) + ug_from(i,j) * self%frac(ilon,ilat,icov)
          end do ! source cells
        end if ! any source cells?
      end do ! lons
    end do ! lats

    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Ugg_Interpol_2d


  ! *
  

  subroutine Grid_Convertor_Ugg_Interpol_3d( self, ug_from, ug_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertor_Ugg), intent(inout)        ::  self
    real, intent(in)                                  ::  ug_from(:,:,:)      ! (nx,ny,nz)
    real, intent(out)                                 ::  ug_to(:,:,:)      ! (nlon,nlat,nz)
    integer, intent(out)                              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertor_Ugg_Interpol_3d'
    
    ! --- local ----------------------------------
    
    integer            ::  ilon, ilat
    integer            ::  ncov, icov
    integer            ::  i, j
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (size(ug_to,1) /= self%nlon) .or. (size(ug_to,2) /= self%nlat) ) then
      write (gol,'("shape of target grid (",i0,",",i0,") differs from definition (",i0,",",i0,")")') &
               size(ug_to,1), size(ug_to,2), self%nlon, self%nlat; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( size(ug_to,3) /= size(ug_from,3) ) then
      write (gol,'("number of levels in target grid (",i0,") differs from source (",i0,")")') &
               size(ug_to,3), size(ug_from,3); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! switch:
    if ( self%nfrom > 0 ) then
      write (gol,'("not implemented yet for froms arrays")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! init output:
    ug_to = 0.0
    ! loop over target cells:
    do ilat = 1, self%nlat
      do ilon = 1, self%nlon
        ! number of source cells:
        ncov = self%ncov(ilon,ilat)
        ! any source cells ?
        if ( ncov > 0 ) then
          ! loop over source cells:
          do icov = 1, ncov
            ! source indeces:
            i = self%ii(ilon,ilat,icov)
            j = self%jj(ilon,ilat,icov)
            ! add contribution:
            ug_to(ilon,ilat,:) = ug_to(ilon,ilat,:) + ug_from(i,j,:) * self%frac(ilon,ilat,icov)
          end do ! source cells
        end if ! any source cells?
      end do ! lons
    end do ! lats
  
    ! ok
    status = 0
    
  end subroutine Grid_Convertor_Ugg_Interpol_3d

  
  ! ********************************************************************
  ! ***
  ! *** convertors
  ! ***
  ! ********************************************************************

  
  subroutine Grid_Convertors_Init( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(out)     ::  self
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! no conversions defined yet:
    self%n_cgg_to_ll  = 0
    ! allocate maximum storage:
    allocate( self%rgg_to_ll(n_convertor_max), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! no conversions defined yet:
    self%n_rgg_to_ll  = 0
    ! allocate maximum storage:
    allocate( self%cgg_to_ll(n_convertor_max), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! no conversions defined yet:
    self%n_ugg_to_ugg = 0
    allocate( self%ugg_to_ugg(n_convertor_max), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Init
  
  ! *
  
  subroutine Grid_Convertors_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)     ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Done'
    
    ! --- local ----------------------------------
    
    integer   ::  i
    
    ! --- begin ----------------------------------
    
    ! loop over initialized convertors:
    do i = 1, self%n_cgg_to_ll
      ! clear:
      call self%cgg_to_ll(i)%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! clear:
    deallocate( self%cgg_to_ll, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over initialized convertors:
    do i = 1, self%n_rgg_to_ll
      ! clear:
      call self%rgg_to_ll(i)%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! clear:
    deallocate( self%rgg_to_ll, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over initialized convertors:
    do i = 1, self%n_ugg_to_ugg
      ! clear:
      call self%ugg_to_ugg(i)%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! clear:
    deallocate( self%ugg_to_ugg, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Grid_Convertors_Done

  
  ! *
  
  
  subroutine Grid_Convertors_Rgg_LL_Select( self, howto, rgg, lli, iconv, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    character(len=*), intent(in)                    ::  howto
    class(T_Grid_Rgg), intent(in)                   ::  rgg
    type(TllGridInfo), intent(in)                   ::  lli
    integer, intent(out)                            ::  iconv
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Rgg_LL_Select'
    
    ! --- local ----------------------------------
    
    type(T_Grid_Key)    ::  rgg_key, lli_key
    logical             ::  match

    ! --- begin ----------------------------------
    
    ! fill characteristic values:
    call rgg_key%Init( 'rgg', (/ rgg%npoint, rgg%nlat, rgg%nulon /), &
                        (/ rgg%longitude(1), rgg%longitude(rgg%npoint), rgg%latitude(1), rgg%latitude(rgg%npoint) /), &
                        status )
    IF_NOT_OK_RETURN(status=1)
    ! fill characteristic values:
    call lli_key%Init( 'lli', (/ lli%nlon, lli%nlat /), &
                        (/ lli%lon_deg(1), lli%lon_deg(lli%nlon), lli%lat_deg(1), lli%lat_deg(lli%nlat) /), &
                        status )
    IF_NOT_OK_RETURN(status=1)
    ! search:
    match = .false.
    do iconv = 1, self%n_rgg_to_ll
      ! compare:
      call self%rgg_to_ll(iconv)%Compare( howto, rgg_key, lli_key, match, status )
      IF_NOT_OK_RETURN(status=1)
      ! found ?
      if ( match ) exit
    end do
    
    ! not found ? then create new:
    if ( .not. match ) then
      ! info ...
      write (gol,'("define new rgg_to_ll convertor ...")'); call goPr
      ! increase counter:
      self%n_rgg_to_ll = self%n_rgg_to_ll + 1
      ! check ...
      if ( self%n_rgg_to_ll > size(self%rgg_to_ll) ) then
        write (gol,'("maximum number of rgg_to_ll conversions exceeded")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      iconv = self%n_rgg_to_ll
      ! switch:
      select case ( trim(howto) )
        ! area remapping:
        case ( 'area-remap' )
          ! new:
          call self%rgg_to_ll(iconv)%Init_AreaRemap( howto, rgg_key, rgg, lli_key, lli, status )
          IF_NOT_OK_RETURN(status=1)
        ! point interpolation:
        case ( 'interpol' )
          ! new:
          call self%rgg_to_ll(iconv)%Init_Interpol( howto, rgg_key, rgg, lli_key, lli, status )
          IF_NOT_OK_RETURN(status=1)
        ! unknown ...
        case default
          write (gol,'("unsupported grid convertor howto `",a,"`")') trim(howto); call goErr
          TRACEBACK; status=1; return
      end select
    end if        

    ! ok
    status = 0

  end subroutine Grid_Convertors_Rgg_LL_Select

  
  ! *
  
  
  subroutine Grid_Convertors_Cgg_LL_Select( self, howto, cgg, lli, iconv, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    character(len=*), intent(in)                    ::  howto
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    type(TllGridInfo), intent(in)                   ::  lli
    integer, intent(out)                            ::  iconv
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_Select'
    
    ! --- local ----------------------------------
   
    type(T_Grid_Key)    ::  cgg_key, lli_key
    logical             ::  match

    ! --- begin ----------------------------------

    ! fill characteristic values:
    call cgg_key%Init( 'cgg', (/ cgg%nlon, cgg%nlat /), &
                        (/ cgg%longitude(1), cgg%longitude(cgg%nlon), cgg%latitude(1), cgg%latitude(cgg%nlat) /), &
                        status )
    IF_NOT_OK_RETURN(status=1)
    ! fill characteristic values:
    call lli_key%Init( 'lli', (/ lli%nlon, lli%nlat /), &
                        (/ lli%lon_deg(1), lli%lon_deg(lli%nlon), lli%lat_deg(1), lli%lat_deg(lli%nlat) /), &
                        status )
    IF_NOT_OK_RETURN(status=1)
    ! search:
    match = .false.
    do iconv = 1, self%n_cgg_to_ll
      ! compare:
      call self%cgg_to_ll(iconv)%Compare( howto, cgg_key, lli_key, match, status )
      IF_NOT_OK_RETURN(status=1)
      ! found ?
      if ( match ) exit
    end do
    
    ! not found ? then create new:
    if ( .not. match ) then
      ! info ...
      write (gol,'("define new cgg_to_ll convertor ...")'); call goPr
      ! increase counter:
      self%n_cgg_to_ll = self%n_cgg_to_ll + 1
      ! check ...
      if ( self%n_cgg_to_ll > size(self%cgg_to_ll) ) then
        write (gol,'("maximum number of cgg_to_ll conversions exceeded")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      iconv = self%n_cgg_to_ll
      ! switch:
      select case ( trim(howto) )
        ! area remapping:
        case ( 'area-remap' )
          ! new:
          call self%cgg_to_ll(iconv)%Init_AreaRemap( howto, cgg_key, cgg, lli_key, lli, status )
          IF_NOT_OK_RETURN(status=1)
        ! point interpolation:
        case ( 'interpol' )
          ! new:
          call self%cgg_to_ll(iconv)%Init_Interpol( howto, cgg_key, cgg, lli_key, lli, status )
          IF_NOT_OK_RETURN(status=1)
        ! nearest value:
        case ( 'nearest' )
          ! new:
          call self%cgg_to_ll(iconv)%Init_Nearest( howto, cgg_key, cgg, lli_key, lli, status )
          IF_NOT_OK_RETURN(status=1)
        ! unknown ...
        case default
          write (gol,'("unsupported grid convertor howto `",a,"`")') trim(howto); call goErr
          TRACEBACK; status=1; return
      end select
    end if        

    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_Select

  
  ! *
  
  
  subroutine Grid_Convertors_Ugg_Select( self, howto, ugg_from, ugg_to, iconv, status, debug )
  
    use GO           , only : GO_Timer_Start, GO_Timer_End
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    character(len=*), intent(in)                    ::  howto
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    integer, intent(out)                            ::  iconv
    integer, intent(out)                            ::  status
    
    logical, intent(in), optional ::  debug

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_Select'
    
    ! --- local ----------------------------------
   
    type(T_Grid_Key)    ::  ugg_key_from, ugg_key_to
    logical             ::  match

    ! --- begin ----------------------------------
    
    ! fill characteristic values:
    call ugg_key_from%Init( trim(ugg_from%type), (/ ugg_from%nlon, ugg_from%nlat /), &
                        (/ ugg_from%longitude(1,1), ugg_from%longitude(ugg_from%nlon,ugg_from%nlat), &
                           ugg_from%latitude (1,1), ugg_from%latitude (ugg_from%nlon,ugg_from%nlat) /), &
                        status )
    IF_NOT_OK_RETURN(status=1)
    ! fill characteristic values:
    call ugg_key_to%Init( trim(ugg_to%type), (/ ugg_to%nlon, ugg_to%nlat /), &
                        (/ ugg_to%longitude(1,1), ugg_to%longitude(ugg_to%nlon,ugg_to%nlat), &
                           ugg_to%latitude (1,1), ugg_to%latitude (ugg_to%nlon,ugg_to%nlat) /), &
                        status )
    IF_NOT_OK_RETURN(status=1)
    ! search:
    match = .false.
      do iconv = 1, self%n_ugg_to_ugg
        ! compare:
        call self%ugg_to_ugg(iconv)%Compare( howto, ugg_key_from, ugg_key_to, match, status )
        IF_NOT_OK_RETURN(status=1)
        ! found ?
        if ( match ) exit
      end do
        
    ! not found ? then create new:
    if ( .not. match ) then
      ! info ...
      write (gol,'("define new ugg_to_ugg convertor ...")'); call goPr
      ! increase counter:
      self%n_ugg_to_ugg = self%n_ugg_to_ugg + 1
      ! check ...
      if ( self%n_ugg_to_ugg > size(self%ugg_to_ugg) ) then
        write (gol,'("maximum number of ugg_to_ugg conversions exceeded")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      iconv = self%n_ugg_to_ugg
      ! switch:
      select case ( trim(howto) )

        ! area remapping:
        case ( 'area-remap' )
          ! new:          
          call self%ugg_to_ugg(iconv)%Init_AreaRemap( howto, ugg_key_from, ugg_from, &
                                                             ugg_key_to  , ugg_to  , status, debug=debug )
          IF_NOT_OK_RETURN(status=1)

        ! point interpolation:
        case ( 'interpol' )
          ! new:
          call self%ugg_to_ugg(iconv)%Init_Interpol( howto, ugg_key_from, ugg_from, &
                                                            ugg_key_to  , ugg_to  , status )
          IF_NOT_OK_RETURN(status=1)

        ! nearest value:
        case ( 'nearest' )
          ! new:
          call self%ugg_to_ugg(iconv)%Init_Nearest( howto, ugg_key_from, ugg_from, ugg_key_to, ugg_to, status )
          IF_NOT_OK_RETURN(status=1)


        ! unknown ...
        case default
          write (gol,'("unsupported grid convertor howto `",a,"`")') trim(howto); call goErr
          TRACEBACK; status=1; return
      end select
      
    end if    ! new convertor    

    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_Select

  
  ! *
  
  
  subroutine Grid_Convertors_Rgg_LL_AreaAver_2d( self, rgg, gg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Rgg), intent(in)                   ::  rgg
    real, intent(in)                                ::  gg(:)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Rgg_LL_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Rgg_LL_Select( 'area-remap', rgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%rgg_to_ll(iconv)%AreaAver( gg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Rgg_LL_AreaAver_2d

  
  ! *
  
  
  subroutine Grid_Convertors_Rgg_LL_AreaAver_3d( self, rgg, gg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Rgg), intent(in)                   ::  rgg
    real, intent(in)                                ::  gg(:,:)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Rgg_LL_AreaAver_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Rgg_LL_Select( 'area-remap', rgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%rgg_to_ll(iconv)%AreaAver( gg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Rgg_LL_AreaAver_3d

  
  ! *
  
  
  subroutine Grid_Convertors_Rgg_LL_Interpol_2d( self, rgg, gg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Rgg), intent(in)                   ::  rgg
    real, intent(in)                                ::  gg(:)    ! (npoint)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:)  ! (nlon,nlat)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Rgg_LL_Interpol_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Rgg_LL_Select( 'interpol', rgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%rgg_to_ll(iconv)%Interpol( gg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Rgg_LL_Interpol_2d

  
  ! *
  
  
  subroutine Grid_Convertors_Rgg_LL_Interpol_3d( self, rgg, gg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Rgg, only : T_Grid_Rgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Rgg), intent(in)                   ::  rgg
    real, intent(in)                                ::  gg(:,:)    ! (npoint,nlev)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:,:)  ! (nlon,nlat,nlev)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Rgg_LL_Interpol_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Rgg_LL_Select( 'interpol', rgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%rgg_to_ll(iconv)%Interpol( gg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Rgg_LL_Interpol_3d

  
  ! ***
  
  
  subroutine Grid_Convertors_Cgg_LL_AreaAver_2d( self, cgg, cg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    real, intent(in)                                ::  cg(:,:)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Cgg_LL_Select( 'area-remap', cgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%cgg_to_ll(iconv)%AreaAver( cg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_AreaAver_2d

  
  ! *
  
  
  subroutine Grid_Convertors_Cgg_LL_AreaAver_3d( self, cgg, cg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    real, intent(in)                                ::  cg(:,:,:)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Cgg_LL_Select( 'area-remap', cgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%cgg_to_ll(iconv)%AreaAver( cg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_AreaAver_3d

  
  ! ***
  
  
  subroutine Grid_Convertors_Cgg_LL_AreaAver_2d_uv( self, cgg, cg_u, cg_v, &
                                                          lli, ll_u, ll_v, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    real, intent(in)                                ::  cg_u(:,:)
    real, intent(in)                                ::  cg_v(:,:)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll_u(:,:)
    real, intent(out)                               ::  ll_v(:,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv
    integer             ::  shp(2)
    real, allocatable   ::  cg_len(:,:)
    real, allocatable   ::  cg2ll_len(:,:)
    real, allocatable   ::  ll_len(:,:)

    ! --- begin ----------------------------------
    
    ! convertor:
    call Grid_Convertors_Cgg_LL_Select( self, 'area-remap', cgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply to u component:
    call Grid_Convertor_Cgg_LL_AreaAver_2d( self%cgg_to_ll(iconv), cg_u, ll_u, status )
    IF_NOT_OK_RETURN(status=1)
    ! apply to v component:
    call Grid_Convertor_Cgg_LL_AreaAver_2d( self%cgg_to_ll(iconv), cg_v, ll_v, status )
    IF_NOT_OK_RETURN(status=1)

    ! shape of input:
    shp = shape(cg_u)
    ! check ...
    if ( any( shape(cg_v) /= shp ) ) then
      write (gol,'("shape of cg_u (",i0,",",i0,") should match shape of cg_v (",i0,",",i0,")")') &
                      shape(cg_u), shape(cg_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( cg_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! input length:
    cg_len = sqrt( cg_u**2 + cg_v**2 )
    
    ! shape of output:
    shp = shape(ll_u)
    ! check ...
    if ( any( shape(ll_v) /= shp ) ) then
      write (gol,'("shape of ll_u (",i0,",",i0,") should match shape of ll_v (",i0,",",i0,")")') &
                      shape(ll_u), shape(ll_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( ll_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! output length:
    ll_len = sqrt( ll_u**2 + ll_v**2 )

    ! storage for input length mapped to output:
    allocate( cg2ll_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! convert input lenght to output grid:
    call Grid_Convertor_Cgg_LL_AreaAver_2d( self%cgg_to_ll(iconv), cg_len, cg2ll_len, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! scale output:
    !  ll_len_new**2 = (ll_u * cg2ll_len/ll_len)**2 +  (ll_v * cg2ll_len/ll_len)**2
    !                = (ll_u**2 + ll_v**2) * cg2ll_len**2 / ll_len**2
    !                         ll_len**2    * cg2ll_len**2 / ll_len**2
    !                = cg2ll_len**2
    ! only scalable if non-zero:
    where ( ll_len > 0.0 )
      ll_u = ll_u * cg2ll_len / ll_len
      ll_v = ll_v * cg2ll_len / ll_len
    end where
    
    ! clear:
    deallocate( ll_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( cg2ll_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( cg_len, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_AreaAver_2d_uv

  
  ! ***
  
  
  subroutine Grid_Convertors_Cgg_LL_AreaAver_3d_uv( self, cgg, cg_u, cg_v, &
                                                          lli, ll_u, ll_v, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    real, intent(in)                                ::  cg_u(:,:,:)
    real, intent(in)                                ::  cg_v(:,:,:)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll_u(:,:,:)
    real, intent(out)                               ::  ll_v(:,:,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_AreaAver_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  nlev
    integer             ::  k

    ! --- begin ----------------------------------
    
    ! number of levels:
    nlev = size(cg_u,3)
    ! check ...
    if ( any( (/size(cg_v,3),size(ll_u,3),size(ll_v,3)/) /= nlev ) ) then
      write (gol,'("3D fields should have same number of layers")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop:
    do k = 1, nlev
      ! apply to 2D field:
      call self%Cgg_LL_AreaAver( cgg, cg_u(:,:,k), cg_v(:,:,k), &
                                 lli, ll_u(:,:,k), ll_v(:,:,k), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! levels
    
    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_AreaAver_3d_uv
  ! *
  
  
  subroutine Grid_Convertors_Cgg_LL_Interpol_2d( self, cgg, cg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    real, intent(in)                                ::  cg(:,:)  ! (i,j)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:)  ! (nlon,nlat)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_Interpol_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Cgg_LL_Select( 'interpol', cgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%cgg_to_ll(iconv)%Interpol( cg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_Interpol_2d

  
  ! *
  
  
  subroutine Grid_Convertors_Cgg_LL_Interpol_3d( self, cgg, cg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    real, intent(in)                                ::  cg(:,:,:)  ! (i,j,:)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:,:)  ! (nlon,nlat,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_Interpol_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Cgg_LL_Select( 'interpol', cgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%cgg_to_ll(iconv)%Interpol( cg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_Interpol_3d

  
  ! *
  
  
  subroutine Grid_Convertors_Cgg_LL_Interpol_3d_uv( self, cgg, cg_u, cg_v, lli, ll_u, ll_v, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    real, intent(in)                                ::  cg_u(:,:,:)  ! (i,j,:)
    real, intent(in)                                ::  cg_v(:,:,:)  ! (i,j,:)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll_u(:,:,:)  ! (nlon,nlat,:)
    real, intent(out)                               ::  ll_v(:,:,:)  ! (nlon,nlat,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_Interpol_3d_uv'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv
    integer             ::  shp(3)
    real, allocatable   ::  cg_len(:,:,:)
    real, allocatable   ::  cg2ll_len(:,:,:)
    real, allocatable   ::  ll_len(:,:,:)

    ! --- begin ----------------------------------
    
    ! convertor:
    call Grid_Convertors_Cgg_LL_Select( self, 'interpol', cgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply to u-component:
    call Grid_Convertor_Cgg_LL_Interpol_3d( self%cgg_to_ll(iconv), cg_u, ll_u, status )
    IF_NOT_OK_RETURN(status=1)
    ! apply to v-component:
    call Grid_Convertor_Cgg_LL_Interpol_3d( self%cgg_to_ll(iconv), cg_v, ll_v, status )
    IF_NOT_OK_RETURN(status=1)

    ! shape of input:
    shp = shape(cg_u)
    ! check ...
    if ( any( shape(cg_v) /= shp ) ) then
      write (gol,'("shape of cg_u (",i0,",",i0,",",i0,") should match shape of cg_v (",i0,",",i0,",",i0,")")') &
                      shape(cg_u), shape(cg_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( cg_len(shp(1),shp(2),shp(3)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! input length:
    cg_len = sqrt( cg_u**2 + cg_v**2 )
    
    ! shape of output:
    shp = shape(ll_u)
    ! check ...
    if ( any( shape(ll_v) /= shp ) ) then
      write (gol,'("shape of ll_u (",i0,",",i0,",",i0,") should match shape of ll_v (",i0,",",i0,",",i0,")")') &
                      shape(ll_u), shape(ll_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( ll_len(shp(1),shp(2),shp(3)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! output length:
    ll_len = sqrt( ll_u**2 + ll_v**2 )

    ! storage for input length mapped to output:
    allocate( cg2ll_len(shp(1),shp(2),shp(3)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! convert input lenght to output grid:
    call Grid_Convertor_Cgg_LL_Interpol_3d( self%cgg_to_ll(iconv), cg_len, cg2ll_len, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! scale output:
    !  ll_len_new**2 = (ll_u * cg2ll_len/ll_len)**2 +  (ll_v * cg2ll_len/ll_len)**2
    !                = (ll_u**2 + ll_v**2) * cg2ll_len**2 / ll_len**2
    !                         ll_len**2    * cg2ll_len**2 / ll_len**2
    !                = cg2ll_len**2
    ! only scalable if non-zero:
    where ( ll_len > 0.0 )
      ll_u = ll_u * cg2ll_len / ll_len
      ll_v = ll_v * cg2ll_len / ll_len
    end where
    
    ! clear:
    deallocate( ll_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( cg2ll_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( cg_len, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_Interpol_3d_uv

  
  ! *
  
  
  subroutine Grid_Convertors_Cgg_LL_Nearest_2d( self, cgg, cg, lli, ll, status )
  
    use Grid         , only : TllGridInfo
    use C3PO_Grid_Cgg, only : T_Grid_Cgg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Cgg), intent(in)                   ::  cgg
    real, intent(in)                                ::  cg(:,:)  ! (i,j)
    type(TllGridInfo), intent(in)                   ::  lli
    real, intent(out)                               ::  ll(:,:)  ! (nlon,nlat)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Cgg_LL_Nearest_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Cgg_LL_Select( 'nearest', cgg, lli, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply, this is a special form of interpolation (single point):
    call self%cgg_to_ll(iconv)%Interpol( cg, ll, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Cgg_LL_Nearest_2d


  ! ***
  
  
  subroutine Grid_Convertors_Ugg_AreaAver_2d( self, ugg_from, ug_from, ugg_to, ug_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from(:,:)  ! (i,j)
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    real, intent(out)                               ::  ug_to(:,:)  ! (nlon,nlat)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Ugg_Select( 'area-remap', ugg_from, ugg_to, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%ugg_to_ugg(iconv)%AreaAver( ug_from, ug_to, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_AreaAver_2d

  
  ! ***
  
  
  subroutine Grid_Convertors_Ugg_AreaAver_2d_uv( self, ugg_from, ug_from_u, ug_from_v, &
                                                          ugg_to, ug_to_u, ug_to_v, status )
  
    use C3PO_Grid_ugg, only : T_Grid_ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from_u(:,:)
    real, intent(in)                                ::  ug_from_v(:,:)
    type(T_Grid_Ugg), intent(in)                    ::  ugg_to
    real, intent(out)                               ::  ug_to_u(:,:)
    real, intent(out)                               ::  ug_to_v(:,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv
    integer             ::  shp(2)
    real, allocatable   ::  cg_len(:,:)
    real, allocatable   ::  ug_from2to_len(:,:)
    real, allocatable   ::  ug_to_len(:,:)

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Ugg_Select( 'area-remap', ugg_from, ugg_to, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply to u component:
    call self%ugg_to_ugg(iconv)%AreaAver( ug_from_u, ug_to_u, status )
    IF_NOT_OK_RETURN(status=1)
    ! apply to v component:
    call self%ugg_to_ugg(iconv)%AreaAver( ug_from_v, ug_to_v, status )
    IF_NOT_OK_RETURN(status=1)

    ! shape of input:
    shp = shape(ug_from_u)
    ! check ...
    if ( any( shape(ug_from_v) /= shp ) ) then
      write (gol,'("shape of ug_from_u (",i0,",",i0,") should match shape of ug_from_v (",i0,",",i0,")")') &
                      shape(ug_from_u), shape(ug_from_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( cg_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! input length:
    cg_len = sqrt( ug_from_u**2 + ug_from_v**2 )
    
    ! shape of output:
    shp = shape(ug_to_u)
    ! check ...
    if ( any( shape(ug_to_v) /= shp ) ) then
      write (gol,'("shape of ug_to_u (",i0,",",i0,") should match shape of ug_to_v (",i0,",",i0,")")') &
                      shape(ug_to_u), shape(ug_to_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( ug_to_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! output length:
    ug_to_len = sqrt( ug_to_u**2 + ug_to_v**2 )

    ! storage for input length mapped to output:
    allocate( ug_from2to_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! convert input lenght to output grid:
    call self%ugg_to_ugg(iconv)%AreaAver( cg_len, ug_from2to_len, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! scale output:
    !  ug_to_len_new**2 = (ug_to_u * ug_from2to_len/ug_to_len)**2 +  (ug_to_v * ug_from2to_len/ug_to_len)**2
    !                = (ug_to_u**2 + ug_to_v**2) * ug_from2to_len**2 / ug_to_len**2
    !                         ug_to_len**2       * ug_from2to_len**2 / ug_to_len**2
    !                = ug_from2to_len**2
    ! only scalable if non-zero:
    where ( ug_to_len > 0.0 )
      ug_to_u = ug_to_u * ug_from2to_len / ug_to_len
      ug_to_v = ug_to_v * ug_from2to_len / ug_to_len
    end where
    
    ! clear:
    deallocate( ug_to_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( ug_from2to_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( cg_len, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_AreaAver_2d_uv

  
  ! *
  
  
  subroutine Grid_Convertors_Ugg_AreaAver_3d( self, ugg_from, ug_from, ugg_to, ug_to, status, debug )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from(:,:,:)  ! (i,j,:)
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    real, intent(out)                               ::  ug_to(:,:,:)  ! (nlon,nlat,:)
    integer, intent(out)                            ::  status

    logical, intent(in), optional   ::  debug

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_AreaAver_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Ugg_Select( 'area-remap', ugg_from, ugg_to, iconv, status, debug=debug )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%ugg_to_ugg(iconv)%AreaAver( ug_from, ug_to, status, debug=debug )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_AreaAver_3d


  ! ***
  
  
  subroutine Grid_Convertors_Ugg_AreaAver_3d_uv( self, ugg_from, ug_from_u, ug_from_v, &
                                                   ugg_to, ug_to_u, ug_to_v, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from_u(:,:,:)
    real, intent(in)                                ::  ug_from_v(:,:,:)
    type(T_Grid_Ugg), intent(in)                    ::  ugg_to
    real, intent(out)                               ::  ug_to_u(:,:,:)
    real, intent(out)                               ::  ug_to_v(:,:,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_AreaAver_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  nlev
    integer             ::  k

    ! --- begin ----------------------------------
    
    ! number of levels:
    nlev = size(ug_from_u,3)
    ! check ...
    if ( any( (/size(ug_from_v,3),size(ug_to_u,3),size(ug_to_v,3)/) /= nlev ) ) then
      write (gol,'("3D fields should have same number of layers")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! loop:
    do k = 1, nlev
      ! apply to 2D field:
      call self%Ugg_AreaAver( ugg_from, ug_from_u(:,:,k), ug_from_v(:,:,k), &
                              ugg_to, ug_to_u(:,:,k), ug_to_v(:,:,k), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! levels
    
    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_AreaAver_3d_uv

  
  ! *

  subroutine Grid_Convertors_Ugg_Interpol_2d( self, ugg_from, ug_from, ugg_to, ug_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from(:,:)  ! (i,j)
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    real, intent(out)                               ::  ug_to(:,:)  ! (nlon,nlat)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_Interpol_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Ugg_Select( 'interpol', ugg_from, ugg_to, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%ugg_to_ugg(iconv)%Interpol( ug_from, ug_to, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_Interpol_2d
  
  
  ! *
  
  
  subroutine Grid_Convertors_Ugg_Interpol_2d_uv( self, ugg_from, ug_from_u, ug_from_v, &
                                                         ugg_to, ug_to_u, ug_to_v, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from_u(:,:)  ! (i,j,:)
    real, intent(in)                                ::  ug_from_v(:,:)  ! (i,j,:)
    type(T_Grid_Ugg), intent(in)                    ::  ugg_to
    real, intent(out)                               ::  ug_to_u(:,:)  ! (nlon,nlat,:)
    real, intent(out)                               ::  ug_to_v(:,:)  ! (nlon,nlat,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_Interpol_2d_uv'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv
    integer             ::  shp(2)
    real, allocatable   ::  ug_from_len(:,:)
    real, allocatable   ::  ug_from2to_len(:,:)
    real, allocatable   ::  ug_to_len(:,:)

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Ugg_Select( 'interpol', ugg_from, ugg_to, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply to u-component:
    call self%ugg_to_ugg(iconv)%Interpol( ug_from_u, ug_to_u, status )
    IF_NOT_OK_RETURN(status=1)
    ! apply to v-component:
    call self%ugg_to_ugg(iconv)%Interpol( ug_from_v, ug_to_v, status )
    IF_NOT_OK_RETURN(status=1)

    ! shape of input:
    shp = shape(ug_from_u)
    ! check ...
    if ( any( shape(ug_from_v) /= shp ) ) then
      write (gol,'("shape of ug_from_u (",i0,",",i0,",",i0,") should match shape of ug_from_v (",i0,",",i0,",",i0,")")') &
                      shape(ug_from_u), shape(ug_from_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( ug_from_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! input length:
    ug_from_len = sqrt( ug_from_u**2 + ug_from_v**2 )
    
    ! shape of output:
    shp = shape(ug_to_u)
    ! check ...
    if ( any( shape(ug_to_v) /= shp ) ) then
      write (gol,'("shape of ug_to_u (",i0,",",i0,",",i0,") should match shape of ug_to_v (",i0,",",i0,",",i0,")")') &
                      shape(ug_to_u), shape(ug_to_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( ug_to_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! output length:
    ug_to_len = sqrt( ug_to_u**2 + ug_to_v**2 )

    ! storage for input length mapped to output:
    allocate( ug_from2to_len(shp(1),shp(2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! convert input lenght to output grid:
    call self%ugg_to_ugg(iconv)%Interpol( ug_from_len, ug_from2to_len, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! scale output:
    !  ug_to_len_new**2 = (ug_to_u * ug_from2to_len/ug_to_len)**2 +  (ug_to_v * ug_from2to_len/ug_to_len)**2
    !                = (ug_to_u**2 + ug_to_v**2) * ug_from2to_len**2 / ug_to_len**2
    !                         ug_to_len**2    * ug_from2to_len**2 / ug_to_len**2
    !                = ug_from2to_len**2
    ! only scalable if non-zero:
    where ( ug_to_len > 0.0 )
      ug_to_u = ug_to_u * ug_from2to_len / ug_to_len
      ug_to_v = ug_to_v * ug_from2to_len / ug_to_len
    end where
    
    ! clear:
    deallocate( ug_to_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( ug_from2to_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( ug_from_len, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_Interpol_2d_uv
  
  
  ! *
    
  
  subroutine Grid_Convertors_Ugg_Interpol_3d( self, ugg_from, ug_from, ugg_to, ug_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from(:,:,:)  ! (i,j,:)
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    real, intent(out)                               ::  ug_to(:,:,:)  ! (nlon,nlat,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_Interpol_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Ugg_Select( 'interpol', ugg_from, ugg_to, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply:
    call self%ugg_to_ugg(iconv)%Interpol( ug_from, ug_to, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_Interpol_3d

  
  ! *
  
  
  subroutine Grid_Convertors_Ugg_Interpol_3d_uv( self, ugg_from, ug_from_u, ug_from_v, &
                                                         ugg_to, ug_to_u, ug_to_v, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from_u(:,:,:)  ! (i,j,:)
    real, intent(in)                                ::  ug_from_v(:,:,:)  ! (i,j,:)
    type(T_Grid_Ugg), intent(in)                    ::  ugg_to
    real, intent(out)                               ::  ug_to_u(:,:,:)  ! (nlon,nlat,:)
    real, intent(out)                               ::  ug_to_v(:,:,:)  ! (nlon,nlat,:)
    integer, intent(out)                            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_Interpol_3d_uv'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv
    integer             ::  shp(3)
    real, allocatable   ::  ug_from_len(:,:,:)
    real, allocatable   ::  ug_from2to_len(:,:,:)
    real, allocatable   ::  ug_to_len(:,:,:)

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Ugg_Select( 'interpol', ugg_from, ugg_to, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply to u-component:
    call self%ugg_to_ugg(iconv)%Interpol( ug_from_u, ug_to_u, status )
    IF_NOT_OK_RETURN(status=1)
    ! apply to v-component:
    call self%ugg_to_ugg(iconv)%Interpol( ug_from_v, ug_to_v, status )
    IF_NOT_OK_RETURN(status=1)

    ! shape of input:
    shp = shape(ug_from_u)
    ! check ...
    if ( any( shape(ug_from_v) /= shp ) ) then
      write (gol,'("shape of ug_from_u (",i0,",",i0,",",i0,") should match shape of ug_from_v (",i0,",",i0,",",i0,")")') &
                      shape(ug_from_u), shape(ug_from_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( ug_from_len(shp(1),shp(2),shp(3)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! input length:
    ug_from_len = sqrt( ug_from_u**2 + ug_from_v**2 )
    
    ! shape of output:
    shp = shape(ug_to_u)
    ! check ...
    if ( any( shape(ug_to_v) /= shp ) ) then
      write (gol,'("shape of ug_to_u (",i0,",",i0,",",i0,") should match shape of ug_to_v (",i0,",",i0,",",i0,")")') &
                      shape(ug_to_u), shape(ug_to_v); call goErr
      TRACEBACK; status=1; return
    end if
    ! storage:
    allocate( ug_to_len(shp(1),shp(2),shp(3)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! output length:
    ug_to_len = sqrt( ug_to_u**2 + ug_to_v**2 )

    ! storage for input length mapped to output:
    allocate( ug_from2to_len(shp(1),shp(2),shp(3)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! convert input lenght to output grid:
    call self%ugg_to_ugg(iconv)%Interpol( ug_from_len, ug_from2to_len, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! scale output:
    !  ug_to_len_new**2 = (ug_to_u * ug_from2to_len/ug_to_len)**2 +  (ug_to_v * ug_from2to_len/ug_to_len)**2
    !                = (ug_to_u**2 + ug_to_v**2) * ug_from2to_len**2 / ug_to_len**2
    !                         ug_to_len**2    * ug_from2to_len**2 / ug_to_len**2
    !                = ug_from2to_len**2
    ! only scalable if non-zero:
    where ( ug_to_len > 0.0 )
      ug_to_u = ug_to_u * ug_from2to_len / ug_to_len
      ug_to_v = ug_to_v * ug_from2to_len / ug_to_len
    end where
    
    ! clear:
    deallocate( ug_to_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( ug_from2to_len, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( ug_from_len, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_Interpol_3d_uv

  
  ! *
  
  
  subroutine Grid_Convertors_Ugg_Nearest_2d( self, ugg_from, ug_from, ugg_to, ug_to, status )
  
    use C3PO_Grid_Ugg, only : T_Grid_Ugg
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    real, intent(in)                                ::  ug_from(:,:)  ! (i,j)
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    real, intent(out)                               ::  ug_to(:,:)  ! (nlon,nlat)
    integer, intent(out)                            ::  status


    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_Nearest_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  iconv

    ! --- begin ----------------------------------
    
    ! convertor:
    call self%Ugg_Select( 'nearest', ugg_from, ugg_to, iconv, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! apply, this is a special form of interpolation (single point):
    call self%ugg_to_ugg(iconv)%Interpol( ug_from, ug_to, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_Nearest_2d
  
  ! *
  
!  subroutine Grid_Convertors_Ugg_IndexFractions( self, ugg_from, ug_from, ugg_to, ug_to, status, clear )
!  
!    use C3PO_Grid_Ugg, only : T_Grid_Ugg
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_Grid_Convertors), intent(inout)         ::  self
!    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
!    integer, intent(in)                             ::  ug_from(:,:)  ! (i,j)
!    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
!    real, intent(out)                               ::  ug_to(:,:,:)  ! (nlon,nlat,nlev)
!    integer, intent(out)                            ::  status
!    logical, intent(in), optional                   ::  clear
!
!
!    ! --- const --------------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_IndexFractions'
!    
!    ! --- local ----------------------------------
!    
!    integer             ::  iconv
!    logical             ::  clear_convertor
!    
!    ! --- begin ----------------------------------
!    
!    ! calculate field directly and clear the convertor after use
!    clear_convertor = .false.
!    if ( present(clear) ) clear_convertor = clear
!    
!    if (clear_convertor ) then      
!      ! convertor, direct application:
!      call self%Ugg_Select( 'index-fractions', ugg_from, ugg_to, iconv, status, ug_from=ug_from, ug_to=ug_to )
!      IF_NOT_OK_RETURN(status=1)
!    else
!      ! convertor:
!      call self%Ugg_Select( 'index-fractions', ugg_from, ugg_to, iconv, status )
!      IF_NOT_OK_RETURN(status=1)
!      ! apply:
!      call self%ugg_to_ugg(iconv)%IndexFractions( ug_from, ug_to, status )
!      IF_NOT_OK_RETURN(status=1)
!    end if
!    
!    ! ok
!    status = 0
!
!  end subroutine Grid_Convertors_Ugg_IndexFractions
  
  ! *
  
  
  !
  ! Output  ug_to(i,j,k) is the fraction of the output cell (i,j)
  ! that has value k in the input.
  !
  ! Mapping is equivalent to:
  ! - for target value k, split input field in maps with area (in m2) where equal to k
  ! - apply 'area-remap', this gives amount of m2 equal to k
  ! - devide by target area to obtain fraction 0-1
  ! In future, implement using a 'area-remap' convertor.
  ! 
  
  subroutine Grid_Convertors_Ugg_IndexFractions( self, ugg_from, ug_from, ugg_to, ug_to, status, &
                                                          clear )
  
    use C3PO_Grid_Ugg  , only : T_Grid_Ugg
    use C3PO_Grid_Tools, only : ll_area_frac_deg
    use GO             , only : T_Polygon, T_PlotFile
  
    ! --- in/out ---------------------------------
    
    class(T_Grid_Convertors), intent(inout)         ::  self
    class(T_Grid_Ugg), intent(in)                   ::  ugg_from
    integer, intent(in)                             ::  ug_from(:,:)  ! (i,j)
    class(T_Grid_Ugg), intent(in)                   ::  ugg_to
    real, intent(out)                               ::  ug_to(:,:,:)  ! (nlon,nlat,nlev)
    integer, intent(out)                            ::  status

    logical, intent(in), optional                   ::  clear

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Grid_Convertors_Ugg_IndexFractions'
    
    ! --- local ----------------------------------
    
    logical              ::  clear_convertor
    logical              ::  reg2reg
    integer              ::  ilon, ilat
    real(8)              ::  bbox(4)
    integer              ::  ibox(4)
    integer              ::  i, j, k
    type(T_Polygon)      ::  pg
    real(8)              ::  pg_area
    real                 ::  frac
    
    !logical              ::  debug
    !real                 ::  sarea
    !type(T_PlotFile)     ::  pf

    ! --- begin ----------------------------------
    
    ! calculate field directly and clear the convertor after use
    clear_convertor = .true.
    if ( present(clear) ) clear_convertor = clear
    ! check ...
    if ( .not. clear_convertor ) then
      write (gol,'("this version could not store a converter but clear=.false.")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( (minval(ug_from) < 1) .or. (maxval(ug_from) > size(ug_to,3)) ) then
      write (gol,'("range of input values ",i0,",..,",i0," should match output size ",i0,",..,",i0)') &
               minval(ug_from), maxval(ug_from), 1, size(ug_to,3); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( .not. allocated(ugg_from%pg) ) then
      write (gol,'("no polygons in input grid definition; used clear=.true. when creating?")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! regular to regular?
    reg2reg = ( (trim(ugg_from%type) == 'cartesian') .or. (trim(ugg_from%type) == 'cartesian-regular') ) .and. &
              ( (trim(ugg_to%type  ) == 'cartesian') .or. (trim(ugg_to%type  ) == 'cartesian-regular') )
   
    ! loop over target cells:
    do ilat = 1, ugg_to%nlat
      do ilon = 1, ugg_to%nlon
      
        !! testing ...
        !debug = (ilon==1) .and. (ilat==1)
        !if (debug) print *, 'xxx debug target cell ', ilon, ilat
        !sarea = 0.0
        
        !! testing ...
        !if ( debug) then
        !  call pf%Init( 'debug-index-frac.py', status )
        !  IF_NOT_OK_RETURN(status=1)
        !  call ugg_to%pg(ilon,ilat)%Plot( 'edges', pf, status, kwargs='color="blue"' )
        !  IF_NOT_OK_RETURN(status=1)
        !end if
      
        ! init target values:
        ug_to(ilon,ilat,:) = 0.0

        ! get bounding box for target cell: (/west,east,south,north/)
        ! this is used to limit searching for source cells:
        call ugg_to%pg(ilon,ilat)%BoundBox( bbox, status )
        IF_NOT_OK_RETURN(status=1)

        ! get range of source cells that might include target cells:
        call ugg_from%GetIndexDomain( real(bbox), ibox, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! loop over source cells:
        do j = ibox(3), ibox(4)
          do i = ibox(1), ibox(2)

            ! target level is value of source:
            k = ug_from(i,j)
            
            ! two regular grids?
            if ( reg2reg ) then
            
              ! compute fraction of target cell covered by source:
              frac = ll_area_frac_deg( ugg_to%longitude_bnds_1d(1,ilon), ugg_to%longitude_bnds_1d(2,ilon), &
                                       ugg_to%latitude_bnds_1d (1,ilat), ugg_to%latitude_bnds_1d (2,ilat), &
                                       ugg_from%longitude_bnds_1d(1,i) , ugg_from%longitude_bnds_1d(2,i), &
                                       ugg_from%latitude_bnds_1d (1,j) , ugg_from%latitude_bnds_1d (2,j) )
                                          
              ! contribution is fraction of target cell covered by intersection:
              ug_to(ilon,ilat,k) = ug_to(ilon,ilat,k) + frac
            
              !! testing ...
              !print *, 'xxx   source ', i, j, frac, sum(ug_to(ilon,ilat,:))

            else

              !! intersection, status<0 if no overlap:
              !call ugg_from%pg(i,j)%Intersection( ugg_to%pg(ilon,ilat), pg, status )
              !                          !debug= debug .and. all((/i,j/)==(/18,7/)) )
              !                          !debug=all((/i,j,ilon,ilat/)==(/1289,914,296,363/)) )
              !if ( status < 0 ) then
              !  ! no overlap, next source cell:
              !  cycle
              !else if ( status > 0 ) then
              !  write (gol,'("from intersection between source (",i0,",",i0,") and target (",i0,",",i0,")")') &
              !           i, j, ilon, ilat; call goErr
              !  TRACEBACK; status=1; return
              !end if
              
              ! intersection, init pg as copy of input clipped by target:
              call pg%InitClippedPolygon( ugg_from%pg(i,j), ugg_to%pg(ilon,ilat), status )
              IF_NOT_OK_RETURN(status=1)
              ! filled?
              if ( .not. pg%IsEmpty() ) then

                ! area of intersection:
                call pg%LonLat_Area( pg_area, status )  !, debug=debug.and.(j==4) )
                IF_NOT_OK_RETURN(status=1)

                ! testing ...
                !if ( debug ) print *, 'xxx source ', i, j, pg_area
                !if ( status /= 0 ) then
                !  write (gol,'("from intersection between source (",i0,",",i0,") and target (",i0,",",i0,")")') &
                !           i, j, ilon, ilat; call goErr
                !  ! ... plot ...
                !  call pf%Init( 'debug-intersection.py', status )
                !  IF_NOT_OK_RETURN(status=1)
                !  call ugg_from%pg(i,j)%Plot( 'edges', pf, status, kwargs='color="red"' )
                !  IF_NOT_OK_RETURN(status=1)
                !  call ugg_to%pg(ilon,ilat)%Plot( 'edges', pf, status, kwargs='color="blue"' )
                !  IF_NOT_OK_RETURN(status=1)
                !  call pg%Plot( 'edges', pf, status, kwargs='color="green"' )
                !  IF_NOT_OK_RETURN(status=1)
                !  call pf%Done( status )
                !  IF_NOT_OK_RETURN(status=1)
                !  ! ..............
                !  TRACEBACK; status=1; return
                !end if

                ! skip if zero (might not have been trapped by Intersection):
                if ( pg_area == 0.0 ) cycle

                ! contribution is fraction of target cell covered by intersection:
                !  1                   1                     m2   /         m2
                ug_to(ilon,ilat,k) = ug_to(ilon,ilat,k) + pg_area / ugg_to%area(ilon,ilat)

                !! testing ...
                !if (debug) print *, 'xxx   add to ', k, pg_area / ugg_to%area(ilon,ilat), ug_to(ilon,ilat,k), sum(ug_to(ilon,ilat,:))
                !sarea = sarea + pg_area
                !if (debug) print *, 'xxx  source area ', i, j, sarea/ugg_to%area(ilon,ilat)
                !if (debug) print *, '  x  source area ', i, j, pg_area/ugg_to%area(ilon,ilat)

                !! testing ...
                !if ( debug ) then
                !  call ugg_from%pg(i,j)%Plot( 'edges', pf, status, kwargs='color="red"' )
                !  IF_NOT_OK_RETURN(status=1)
                !  !if ( all( (/i,j/) == (/18,4/) ) ) then
                !  !  call pg%Plot( 'edges', pf, status, kwargs='color="green", marker="o", linestyle="-", linewidth="2"' )
                !  !  IF_NOT_OK_RETURN(status=1)
                !  !  call pg%Plot( 'fill', pf, status, kwargs='color="green", alpha=0.5' )
                !  !  IF_NOT_OK_RETURN(status=1)
                !  !end if
                !end if

              end if ! pg filled
              
              ! done:
              call pg%Done( status )
              IF_NOT_OK_RETURN(status=1)
              
            end if ! reg2reg

          end do ! i
        end do ! j
        
        !! testing ..
        !if ( debug) then
        !  call pf%Done( status )
        !  IF_NOT_OK_RETURN(status=1)
        !end if
        
        !! testing ..
        !if (debug) stop

      end do  ! ugg_to ilon
    end do  ! ugg_to ilat
    
    ! ok
    status = 0

  end subroutine Grid_Convertors_Ugg_IndexFractions
  
  ! *


end module C3PO_Grid_Convert

