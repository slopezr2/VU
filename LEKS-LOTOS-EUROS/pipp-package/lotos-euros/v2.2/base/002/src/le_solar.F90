!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Solar

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  zangle_init
  public  ::  solardec
  public  ::  zangle


  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LE_Solar'


  ! --- var -------------------------------
  
  real  ::  soldec
  

contains

  subroutine zangle_init( status )

    use LE_Data, only : LE_Data_Enable
  
    ! --- in/out ---
    integer, intent(out)                  ::  status
    
    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/zangle_init'
      
    ! --- begin ---

    call LE_Data_Enable( 'lon', status )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_Enable( 'lat', status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine zangle_init
  
  ! this routine computes for each grid cell the cosine
  ! of the solar zenith angle

  subroutine zangle( ih, minutes, status )

    use Dims, only : nx, ny, zenith, coszen
    use Dims, only : zenith, coszen
    use LE_Data, only : LE_Data_GetPointer
    
    ! --- in/out ---------------------------------

    integer, intent(in)   ::  ih
    real, intent(in)      ::  minutes
    integer, intent(out)  ::  status

    
    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/zangle'
    ! --- local ----------------------------------

    real      ::  tod, lambda, phi
    integer   ::  i, j
    
    real, pointer    ::  lons(:,:,:)   ! (lon,lat,1)                                
    real, pointer    ::  lats(:,:,:)   ! (lon,lat,1)                                
    ! --- begin ---------------------------------

    call LE_Data_GetPointer( 'lon', lons, status, check_units ='degrees_east' )
    IF_NOTOK_RETURN(status=1)
    call LE_Data_GetPointer( 'lat', lats, status, check_units ='degrees_north' )
    IF_NOTOK_RETURN(status=1)

    write (gol,'("computation of zenith angle ...")'); call goPr
    
    tod  = float(ih) + minutes/60.0

    ! loop over all grid cells in hor. plane
    do i = 1, nx
      do j = 1, ny
        lambda = lons(i,j,1)
        phi = lats(i,j,1)
        ! compute solar angle for this point
        call zenangle(lambda,phi,tod,zenith(i,j) )
      end do
    end do

    ! set cosine of solar angle
    coszen = max(0.0, cos(zenith) )
    
    ! ok
    status = 0

  end subroutine
  
  
  ! ***


  subroutine zenangle( lambda, phi, tod, zen )

    use Binas, only : deg2rad
    
    ! --- in/out ---------------------------------

    real, intent(in)    :: lambda   ! longitude [degrees]
    real, intent(in)    :: phi      ! latitude [degrees]
    real, intent(in)    :: tod      ! time-of-day [hours] ; 2.5 is 02:30
    real, intent(out)   :: zen      ! zenit angle [radians]
    
    ! --- const ---------------------------------
    
    real, parameter :: conv=.26179939

    ! --- local ----------------------------------

    real :: ss,cc,t0
    
    ! --- begin ---------------------------------

    ss = sin(soldec) * sin(phi*deg2rad)
    cc = cos(soldec) * cos(phi*deg2rad)
    t0= (tod + lambda*24./360. - 12.67)*conv
    !coszen=max(0.0, ss+cc*cos(t0) )
    zen = acos( ss+cc*cos(t0) )

  end subroutine zenangle
  
  
  ! ***
  

  subroutine solardec( im, id )

    ! --- in/out ---------------------------------

    integer, intent(in)   :: im, id
    
    ! --- const ----------------------------------
    
    integer, parameter :: ndays(12) =  &
           (/0,31,59,90,120,151,181,212,243,273,304,334/)

    ! --- in/out ---------------------------------
    
    real                :: dd
    
    ! --- begin ----------------------------------
    
    write (gol,'("computation of solar declination ...")'); call goPr

    ! calculate delta (needed for calculation of solar angle)
    dd = 2.0 * 3.1415927 * (ndays(im)+id) / 365.0
    soldec = 0.006918 - 0.399912*cos(    dd) + 0.070257*sin(    dd) &
                      - 0.006758*cos(2.0*dd) + 0.000907*sin(2.0*dd) &
                      - 0.002697*cos(3.0*dd) + 0.00148 *sin(3.0*dd)

  end subroutine solardec


end module LE_Solar

