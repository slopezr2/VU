!###############################################################################
!
! Radiation lookup table
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################


module LE_Radiation_LUT
 
  use GO, only : gol, goPr, goErr
  use NetCDF , only : NF90_NOERR, NF90_StrError
 
  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  T_Radiation_LUT


  ! --- const ------------------------------------

  character(len=*), parameter ::  mname = 'LE_Radiation_LUT'

  ! lookup tables:  
  integer, parameter :: LUT_NO_FILENAMES = 2
  integer, parameter :: LUT_NO_TABLES = 12

  ! debug ?
  logical, parameter ::  LE_DEBUG = .true.


  ! --- types ----------------------------------

  type T_Radiation_LUT
    ! the LUT's have three indices, in order: Rg, Re, Im. 
    ! CHANGED: first one is 2 * pi * Rg / lambda
    real, allocatable   ::  ext_159(:,:,:)
    real, allocatable   ::  ext_200(:,:,:)
    real, allocatable   ::  a_159(:,:,:)
    real, allocatable   ::  a_200(:,:,:)
    real, allocatable   ::  g_159(:,:,:)
    real, allocatable   ::  g_200(:,:,:)
    real, allocatable   ::  TabInd_Rg(:)
    real, allocatable   ::  TabInd_Re(:)
    real, allocatable   ::  TabInd_Im(:)
    real, allocatable   ::  IndGen_Rg(:)
    real, allocatable   ::  IndGen_Re(:)
    real, allocatable   ::  IndGen_Im(:)
    ! other:
    real, allocatable   ::  opac(:,:) 
    real, allocatable   ::  echam(:,:)
    real, allocatable   ::  segelstein(:,:)
  contains
    procedure   ::  Init      => Radiation_LUT_Init
    procedure   ::  Done      => Radiation_LUT_Done
  end type T_Radiation_LUT
  



contains


  ! ====================================================================
  ! ===
  ! === lookup tables
  ! ===
  ! ====================================================================

  !
  ! Reads in Netcdf look-up tables
  !
  ! output: 
  !   RefracIndexPerSpecie(sw_bands, n_species) 
  !   RefracIndexWater    (sw_bands)
  !
  subroutine Radiation_LUT_Init( self, rcF, status )     

    use GO     , only : TrcFile, ReadRc
    use Indices
    use netcdf

    ! --- in/out ---------------------------------
    
    class(T_Radiation_LUT), intent(out)   ::  self
    type(TrcFile), intent(in)             ::  rcF
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/RefracIndex_LUT_Init'

    ! --- types --------------------------------------

    ! Need the sw bands, I could not find them easily specfied tables
    type le_table 
      character(len=32) :: name
      integer           :: file_no
    end type le_table

    ! --- local ----------------------------------

    ! info per lut:
    type(le_table)      :: le_table_info(LUT_NO_TABLES)
    character(len=256)  :: le_lut_filenames(LUT_NO_FILENAMES)
    integer             :: i, j
    integer             :: ncids(LUT_NO_FILENAMES)
    integer             :: NoDims
    integer             :: ncid, varid
    integer             :: IdDims (NF90_MAX_VAR_DIMS)
    integer             :: LenDims(NF90_MAX_VAR_DIMS)

!    integer                        ::  bandwidth
!    real, allocatable   ::  tmp1(:), tmp2(:)
!    real                            ::  Re, Im

!    integer                         ::  i_aerosol, ispec

    ! --- begin ----------------------------------

    ! info ...
    if ( LE_DEBUG ) then
      write (gol,'(a,": start")') rname; call goPr
    end if

     ! settings:
    call ReadRc( rcF, 'le.radiation.lut.1', le_lut_filenames(1), status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'le.radiation.lut.2', le_lut_filenames(2), status )
    IF_NOTOK_RETURN(status=1)

    ! fill table:
    le_table_info(1) =  le_table('Opac'      , 1)
    le_table_info(2) =  le_table('ECHAM-HAM' , 1)
    le_table_info(3) =  le_table('Segelstein', 1)
    le_table_info(4) =  le_table('ext_159'   , 2)
    le_table_info(5) =  le_table('ext_200'   , 2)
    le_table_info(6) =  le_table('a_159'     , 2)
    le_table_info(7) =  le_table('a_200'     , 2)
    le_table_info(8) =  le_table('g_159'     , 2)
    le_table_info(9) =  le_table('g_200'     , 2)
    le_table_info(10) = le_table('x'         , 2)
    le_table_info(11) = le_table('mr'        , 2)
    le_table_info(12) = le_table('mi'        , 2)

    ! loop over input files:
    do i=1, LUT_NO_FILENAMES
      ! info ... 
      if (LE_DEBUG) then
        write (gol,'(a,": using lookup table: ",a)') rname, trim(le_lut_filenames(i)); call goPr
      end if
      ! open:
      status = NF90_Open( le_lut_filenames(i), NF90_NOWRITE, ncids(i) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! loop over tables:
    do i = 1, LUT_NO_TABLES  

      ! info ... 
      if (LE_DEBUG) then
        write (gol,'(a,": reading lut ",i6)') rname, i; call goPr
      end if

      ! current:
      ncid = ncids( le_table_info(i)%file_no )

      status = nf90_inq_varid(ncid, le_table_info(i)%name, varid)
      IF_NF90_NOT_OK_RETURN(status=1)

      status = nf90_inquire_variable(ncid, varid, ndims=NoDims)
      IF_NF90_NOT_OK_RETURN(status=1)

      status = nf90_inquire_variable(ncid, varid, dimids=IdDims(:NoDims) )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! info:
      if (LE_DEBUG) then
        write (gol,'(a,": reading: ",a," from file ",i0)') rname, &
                trim(le_table_info(i)%name), le_table_info(i)%file_no; call goPr
      end if

      ! loop over dimensions:
      do j=1, NoDims 
        status = nf90_inquire_dimension(ncid, IdDims(j), len=LenDims(j))
        IF_NF90_NOT_OK_RETURN(status=1)
        ! info ..
        if (LE_DEBUG) then
          write (gol,'(a,": dimensions: ",i6)') rname, LenDims(j); call goPr
        end if
      end do
      ! switch:
      select case (i) 
        case (1)
          allocate( self%opac(LenDims(1),LenDims(2)) )
          status=nf90_get_var(ncid, varid, self%opac)
          IF_NF90_NOT_OK_RETURN(status=1)
        case (2) 
          allocate(self%echam(LenDims(1),LenDims(2)))
          status=nf90_get_var(ncid, varid, self%echam)
          IF_NF90_NOT_OK_RETURN(status=1)
        case (3) 
          allocate(self%segelstein(LenDims(1),LenDims(2)))
          status=nf90_get_var(ncid, varid, self%segelstein)
          IF_NF90_NOT_OK_RETURN(status=1)
        case (4) 
          allocate(self%ext_159(LenDims(1),LenDims(2),LenDims(3)))
          status=nf90_get_var(ncid, varid, self%ext_159)
          IF_NF90_NOT_OK_RETURN(status=1)
        case (5) 
          allocate(self%ext_200(LenDims(1),LenDims(2),LenDims(3)))
          status=nf90_get_var(ncid, varid, self%ext_200)
          IF_NF90_NOT_OK_RETURN(status=1)
        case (6) 
          allocate(self%a_159(LenDims(1),LenDims(2),LenDims(3)))
          status=nf90_get_var(ncid, varid, self%a_159)
          IF_NF90_NOT_OK_RETURN(status=1)
        case (7) 
          allocate(self%a_200(LenDims(1),LenDims(2),LenDims(3)))
          status=nf90_get_var(ncid, varid, self%a_200)
          IF_NF90_NOT_OK_RETURN(status=1)
        case (8) 
          allocate(self%g_159(LenDims(1),LenDims(2),LenDims(3)))
          status=nf90_get_var(ncid, varid, self%g_159)
          IF_NF90_NOT_OK_RETURN(status=1)
        case (9) 
          allocate(self%g_200(LenDims(1),LenDims(2),LenDims(3)))
          status=nf90_get_var(ncid, varid, self%g_200)
          print *, 'g_200, nx,ny,nz', LenDims(1),LenDims(2),LenDims(3)
        case (10)
          allocate(self%TabInd_Rg(LenDims(1)))
          status=nf90_get_var(ncid, varid, self%TabInd_Rg)
          IF_NF90_NOT_OK_RETURN(status=1)
          allocate(self%IndGen_Rg(LenDims(1)))
          self%IndGen_Rg = (/(i,i=1,LenDims(1))/)
        case (11)
          allocate(self%TabInd_Re(LenDims(1)))
          status=nf90_get_var(ncid, varid, self%TabInd_Re)
          IF_NF90_NOT_OK_RETURN(status=1)
          allocate(self%IndGen_Re(LenDims(1))); 
          self%IndGen_Re = (/(i,i=1,LenDims(1))/)
        case (12)
          allocate(self%TabInd_Im(LenDims(1)))
          status=nf90_get_var(ncid, varid, self%TabInd_Im)
          IF_NF90_NOT_OK_RETURN(status=1)
          allocate(self%IndGen_Im(LenDims(1)))
          self%IndGen_Im = (/(i,i=1,LenDims(1))/)
      end select

    end do  ! lookup tables

    ! loop over input files:
    do i = 1, LUT_NO_FILENAMES
      ! close:
      status = NF90_Close( ncids(i) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    !print *, 'refracindexperspecie', RefracIndexPerSpecie
    ! info ...
    if ( LE_DEBUG ) then
      write (gol,'(a,": end")') rname; call goPr
    end if

    ! ok
    status = 0

  end subroutine Radiation_LUT_Init
  
  
  ! ***


  subroutine Radiation_LUT_Done( self, status )     

    ! --- in/out ---------------------------------
    
    class(T_Radiation_LUT), intent(inout)   ::  self
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/Radiation_LUT_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! clear lut:
    deallocate( self%ext_159   )
    deallocate( self%ext_200   )
    deallocate( self%a_159     )
    deallocate( self%a_200     )
    deallocate( self%g_159     )
    deallocate( self%g_200     )
    deallocate( self%TabInd_Rg )
    deallocate( self%TabInd_Re )
    deallocate( self%TabInd_Im )
    deallocate( self%IndGen_Rg )
    deallocate( self%IndGen_Re )
    deallocate( self%IndGen_Im )
    
    ! clear:
    deallocate( self%opac, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%echam, stat=status )
    IF_NOTOK_RETURN(status=1)
    deallocate( self%segelstein, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Radiation_LUT_Done


end module LE_Radiation_LUT
