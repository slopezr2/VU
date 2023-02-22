!#################################################################
!
! NAME
!   GO_LUT - Look-Up-Table
!
! USAGE
!
!   type(T_LUT)   ::  lut
!
!   !
!   ! Init table with 2 axis (max. 7 axis supported).
    !
!   call LUT_Init( lut, 2, status )
!   IF_NOTOK_RETURN(status=1)
!
!   !
!   ! Each axis is defined as a number of equally spaced bins;
!   ! bins are defined by an offset, bin size (scale), 
!   ! and number of bins (n).
!   ! The mid of a bin i is then:
!   !   value_i  =   offset + scale * (i-0.5)       ,  i = 1 .. n
!   !
!   call LUT_Def_Axis( lut, 1, 'name1', 20, 0.0, 0.1, status )
!   IF_NOTOK_RETURN(status=1)
!   call LUT_Def_Axis( lut, 2, 'name2', 40, 10.0, 2.0, status )
!   IF_NOTOK_RETURN(status=1)
!
!   !
!   ! Fill values:
!   !
!   ! loop over first coordinate:
!   do i = 1, 20
!     ! extract bin value:
!     call LUT_Coordinate( lut, 1, i, coor1, status )
!     IF_NOTOK_RETURN(status=1)
!     ! loop over second coordinate:
!     do j = 1, 40
!       ! extract bin value:
!       call LUT_Coordinate( lut, 2, j, coor2, status )
!       IF_NOTOK_RETURN(status=1)
!       ! evaluate expensive function:
!       value = ... coor1, coor2, ...
!       ! store:
!       call LUT_Set_Value( lut, (/i,j/), value, status )
!       IF_NOTOK_RETURN(status=1)
!     end do
!   end do
!
!   ! Get index given single coordinate value:
!   call LUT_Index( lut, 1, 1.12, i, status )
!   IF_NOTOK_RETURN(status=1)
!   call LUT_Index( lut, 1,  3.6, j, status )
!   IF_NOTOK_RETURN(status=1)
!
!   ! Get value given indices:
!   call LUT_Get( lut, (/i,j/), value, status )
!   IF_NOTOK_RETURN(status=1)
!
!   ! clear:
!   call LUT_Done( lut, status )
!   IF_NOTOK_RETURN(status=1)
!
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#################################################################


!module GO_Print
!  implicit none
!  public
!  character(len=1024)    ::  gol
!contains
!  subroutine goPr()
!    write (*,'(a)') trim(gol)
!  end subroutine goPr
!  subroutine goErr()
!    write (*,'("ERROR - ",a)') trim(gol)
!  end subroutine goErr
!end module GO_Print


!#################################################################
    

module GO_LUT

  use GO_Print, only : gol, goPr, goErr
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public    ::  T_LUT
  public    ::  LUT_Init, LUT_Done
  public    ::  LUT_Def_Axis
  public    ::  LUT_Coordinate, LUT_Index
  public    ::  LUT_Set_Value, LUT_Get_Value
  

  ! --- const ------------------------------------
  
  character(len=*), parameter   ::  mname = 'GO_LUT'
  
    
  ! --- types ------------------------------------
  
  ! max 7 dimensions in fortran arrays ...
  integer, parameter      ::  maxdim = 7
  
  ! look-up table:
  type T_LUT
    ! dimension:
    integer                       ::  ndim
    ! axis:
    character(len=64)             ::  dimname(maxdim)
    real                          ::  offset(maxdim)
    real                          ::  scale(maxdim)
    integer                       ::  shp(maxdim)
    logical                       ::  axis_defined(maxdim)
    ! values:
    real, pointer                 ::  values(:,:,:,:,:,:,:)
    logical, pointer              ::  filled(:,:,:,:,:,:,:)
    logical                       ::  defined
    ! meta:
    character(len=256)            ::  name, units, long_name
  end type T_LUT


contains


  ! ===================================================================


  subroutine LUT_Init( lut, ndim, status, &
                         name, long_name, units )
  
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(out)                 ::  lut
    integer, intent(in)                      ::  ndim
    integer, intent(out)                     ::  status
    
    character(len=*), intent(in), optional   ::  name
    character(len=*), intent(in), optional   ::  long_name
    character(len=*), intent(in), optional   ::  units
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Init'
  
    ! --- begin ----------------------------------

    ! check ..
    if ( (ndim < 1) .or. (ndim > maxdim) ) then
      write (gol,'("number of dimensions ",i6," out of range 1 .. ",i6)') ndim, maxdim; call goPr
      TRACEBACK; status=1; return
    end if
    ! store:
    lut%ndim = ndim

    ! single dimensions by default:
    lut%shp = 1
    ! no axis 1..ndim defined yet:
    lut%axis_defined = .true.
    lut%axis_defined(1:ndim) = .false.
   
    ! set overall flag
    lut%defined = .false.
    
    ! meta:
    lut%name       = ''
    lut%long_name  = ''
    lut%units      = ''
    ! fill:
    if ( present(name     ) ) lut%name      = trim(name     )
    if ( present(long_name) ) lut%long_name = trim(long_name)
    if ( present(units    ) ) lut%units     = trim(units    )
    
    ! for safety ...
    nullify( lut%values )
    nullify( lut%filled )

    ! ok
    status = 0
    
  end subroutine LUT_Init
  
  
  ! ***
  

  subroutine LUT_Done( lut, status )
  
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(inout)    ::  lut
    integer, intent(out)          ::  status
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Done'
  
    ! --- begin ----------------------------------

    ! clear:
    if ( associated(lut%values) ) deallocate( lut%values )
    if ( associated(lut%filled) ) deallocate( lut%filled )   
    
    ! ok
    status = 0
    
  end subroutine LUT_Done
  
  
  ! ***
  

  subroutine LUT_Def_Axis( lut, idim, name, n, offset, scale, status )
  
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(inout)      ::  lut
    integer, intent(in)             ::  idim
    character(len=*), intent(in)    ::  name
    integer, intent(in)             ::  n
    real, intent(in)                ::  offset
    real, intent(in)                ::  scale
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Def_Axis'
  
    ! --- begin ----------------------------------
    
    ! check ..
    if ( (idim < 1) .or. (idim > lut%ndim) ) then
      write (gol,'("dimension ",i6," out of range 1 .. ",i6)') idim, lut%ndim; call goPr
      TRACEBACK; status=1; return
    end if

    ! already defined ?
    if ( lut%axis_defined(idim) ) then
      write (gol,'("axis ",i6," already defined")') idim; call goPr
      TRACEBACK; status=1; return
    end if
    
    ! store:
    lut%dimname(idim) = trim(name)
    lut%offset(idim) = offset
    lut%scale (idim) = scale
    lut%shp   (idim) = n
    
    ! reset flag:
    lut%axis_defined(idim) = .true.   
    
    ! ok
    status = 0
    
  end subroutine LUT_Def_Axis
  

  ! ***
  

  subroutine LUT_Coordinate( lut, idim, indx, coor, status )
  
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(in)         ::  lut
    integer, intent(in)             ::  idim
    integer, intent(in)             ::  indx
    real, intent(out)               ::  coor
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Coordinate'
  
    ! --- begin ----------------------------------
    
    ! check ..
    if ( (idim < 1) .or. (idim > lut%ndim) ) then
      write (gol,'("dimension ",i6," out of range 1 .. ",i6)') idim, lut%ndim; call goPr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( (indx < 1) .or. (indx > lut%shp(idim)) ) then
      write (gol,'("axis ",i1," `",a,"` index ",i6," out of range 1 .. ",i6)') &
               idim, trim(lut%dimname(idim)), indx, lut%shp(idim); call goPr
      TRACEBACK; status=1; return
    end if

    ! fill mid of bin:
    coor = lut%offset(idim) + (indx-0.5) * lut%scale(idim)
    
    ! ok
    status = 0
    
  end subroutine LUT_Coordinate
  

  ! ***
  

  subroutine LUT_Index( lut, idim, coor, indx, status )
  
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(in)         ::  lut
    integer, intent(in)             ::  idim
    real, intent(in)                ::  coor
    integer, intent(out)            ::  indx
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Index'
  
    ! --- begin ----------------------------------
    
    ! check ..
    if ( (idim < 1) .or. (idim > lut%ndim) ) then
      write (gol,'("dimension ",i6," out of range 1 .. ",i6)') idim, lut%ndim; call goPr
      TRACEBACK; status=1; return
    end if
    
    ! compute index:
    indx = nint( ( coor - lut%offset(idim) )/lut%scale(idim) + 0.5 )

    ! check ..
    if ( (indx < 1) .or. (indx > lut%shp(idim)) ) then
      write (gol,'("axis ",i1," `",a,"` index ",i6," out of range 1 .. ",i6)') &
               idim, trim(lut%dimname(idim)), indx, lut%shp(idim); call goPr
      write (gol,*) '  coor range : ', lut%offset(idim), lut%offset(idim)+lut%scale(idim)*lut%shp(idim); call goErr
      write (gol,*) '  coordinate : ', coor; call goErr
      write (gol,'("increase the number of bins in the axis definition")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine LUT_Index
  

  ! ***
  
  
  subroutine LUT_Set_Value( lut, indices, value, status )
  
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(inout)      ::  lut
    integer, intent(in)             ::  indices(:)
    real, intent(in)                ::  value
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Set_Value'
    
    ! --- local ----------------------------------
    
    integer     ::  idim
    integer     ::  inds(maxdim)
  
    ! --- begin ----------------------------------
    
    ! check ...
    if ( lut%defined ) then
      write (gol,'("table already defined:")'); call goErr
      write (gol,'("indices in current call : ",7i6)') indices; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( size(indices) /= lut%ndim ) then
      write (gol,'("size of indices ",i6," while ndim ",i6)') size(indices), lut%ndim; call goPr
      TRACEBACK; status=1; return
    end if
    
    ! default indices for maximum number of dimensions:
    inds = 1
    ! loop over dims:
    do idim = 1, lut%ndim
      ! check ...
      if ( .not. lut%axis_defined(idim) ) then
        write (gol,'("axis ",i1,"not defined yet")') idim; call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( (indices(idim) < 1) .or. (indices(idim) > lut%shp(idim)) ) then
        write (gol,'("index ",i6," for axis ",i1," out of range 1 .. ",i6)') &
            indices(idim), idim, lut%shp(idim); call goPr
        TRACEBACK; status=1; return
      end if
      ! store:
      inds(idim) = indices(idim)
    end do
    
    ! value arrays not allocated yet ?
    if ( .not. associated(lut%values) ) then
      ! setup storage:
      allocate( lut%values(lut%shp(1),lut%shp(2),lut%shp(3),lut%shp(4),lut%shp(5),lut%shp(6),lut%shp(7)) )
      ! not filled yet:
      allocate( lut%filled(lut%shp(1),lut%shp(2),lut%shp(3),lut%shp(4),lut%shp(5),lut%shp(6),lut%shp(7)) )
      lut%filled = .false.
    end if
    
    ! check ...
    if ( lut%filled(inds(1),inds(2),inds(3),inds(4),inds(5),inds(6),inds(7)) ) then
      write (gol,'("value already defined for:")'); call goErr
      do idim = 1, lut%ndim
        write (gol,'("  ",i6," ",a)') inds(idim), lut%dimname(idim)
      end do
      TRACEBACK; status=1; return
    end if
    
    ! store:
    lut%values(inds(1),inds(2),inds(3),inds(4),inds(5),inds(6),inds(7)) = value
    
    ! update flags:
    lut%filled(inds(1),inds(2),inds(3),inds(4),inds(5),inds(6),inds(7)) = .true.
    ! overall flag:
    lut%defined = all( lut%filled )
    ! done ?
    if ( lut%defined ) then
      ! clear:
      deallocate( lut%filled )
    end if

    ! ok
    status = 0
    
  end subroutine LUT_Set_Value


  ! ***
  
  
  subroutine LUT_Get_Value( lut, indices, value, status )
  
    ! --- in/out ---------------------------------
    
    type(T_LUT), intent(inout)      ::  lut
    integer, intent(in)             ::  indices(:)
    real, intent(out)               ::  value
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
  
    character(len=*), parameter   ::  rname = mname//'/LUT_Get_Value'
    
    ! --- local ----------------------------------
    
    integer     ::  inds(maxdim)
    integer     ::  idim
  
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. lut%defined ) then
      write (gol,'("table not completely defined yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( size(indices) /= lut%ndim ) then
      write (gol,'("size of indices ",i6," while ndim ",i6)') size(indices), lut%ndim; call goPr
      TRACEBACK; status=1; return
    end if
    
    ! default indices for maximum number of dimensions:
    inds = 1
    ! loop over dims:
    do idim = 1, lut%ndim
      ! check ...
      if ( (inds(idim) < 1) .or. (inds(idim) > lut%shp(idim)) ) then
        write (gol,'("value for coordinate ",i1," `",a,"` out of range:")') idim, trim(lut%dimname(idim)); call goErr
        write (gol,*) '  range : ', lut%offset(idim), lut%offset(idim)+lut%scale(idim)*lut%shp(idim); call goErr
        write (gol,'("increase the number of bins in the axis definition")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! store:
      inds(idim) = indices(idim)
    end do
        
    ! extract:
    value = lut%values(inds(1),inds(2),inds(3),inds(4),inds(5),inds(6),inds(7))

    ! ok
    status = 0
    
  end subroutine LUT_Get_Value


end module GO_LUT


!#######################################################################
!###
!### test
!###
!#######################################################################
!
!#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!
!program test
!
!  use GO_Print, only : gol, goPr, goErr
!  use GO_LUT
!  
!  implicit none
!  
!  character(len=*), parameter   ::  rname = 'test'
!
!  integer         ::  status
!  type(T_LUT)     ::  lut
!  integer         ::  i, j
!  real            ::  x, y, value
!  
!  ! info ...
!  write (gol,'("begin test")'); call goPr
!  
!  ! setup table:
!  call LUT_Init( lut, 2, status )
!  IF_NOTOK_STOP
!  
!  ! define axis:
!  call LUT_Def_Axis( lut, 1, 'xx', 5, 0.0, 10.0, status )
!  IF_NOTOK_STOP
!  call LUT_Def_Axis( lut, 2, 'yy', 4, 0.1, 0.2, status )
!  IF_NOTOK_STOP
!
!  ! fill:
!  print *, ''
!  print *, 'set values ...'
!  do j = 1, 4
!    do i = 1, 5
!      ! coordinates:
!      call LUT_Coordinate( lut, 1, i, x, status )
!      IF_NOTOK_STOP
!      call LUT_Coordinate( lut, 2, j, y, status )
!      IF_NOTOK_STOP
!      ! info ..
!      print *, 'coor ', i, j, x, y, value
!      ! set value:
!      value = x + y
!      ! store:
!      call LUT_Set_Value( lut, (/i,j/), value, status )
!      IF_NOTOK_STOP
!    end do
!  end do
!
!  ! extract:
!  print *, ''
!  print *, 'get values ...'
!  do j = 1, 4
!    y = 0.2 + j*0.05
!    do i = 0, 50
!      x = i*2.0
!      ! get nearest value:
!      call LUT_Get_Value( lut, (/x,y/), value, status )
!      IF_NOTOK_STOP
!      ! info ..
!      print *, 'coor ', x, y, value
!    end do
!  end do
!
!  ! clear:
!  call LUT_Done( lut, status )
!  IF_NOTOK_STOP
!    
!  ! info ...
!  write (gol,'("end test")'); call goPr
!  
!end program test
