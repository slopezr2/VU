!###############################################################################
!
! Soil composition of base-cat-ions.
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

module LE_LandUse_BCatIon

  use GO, only : gol, goPr, goErr
 
  implicit none


  ! --- in/out -----------------------------
  
  private
  
  public  ::  nbcat
  public  ::  ibcat_Ca, ibcat_Mg, ibcat_K, ibcat_Na
  public  ::  bcat_frac

  public  ::  LandUse_BCatIon_Init, LandUse_BCatIon_Done


  ! --- const -----------------------------

  character(len=*), parameter   ::  mname = 'LandUse_BCatIon'

  ! number of base cat ions in database:
  integer, parameter :: nbcat = 4
  ! indices:
  integer, parameter :: ibcat_Ca = 1
  integer, parameter :: ibcat_Mg = 2
  integer, parameter :: ibcat_K  = 3
  integer, parameter :: ibcat_Na = 4


  ! --- var --------------------------------------

  ! basecation fraction in topsoil:
  real, allocatable   ::  bcat_frac(:,:,:)  ! nx,ny,nbcat



contains


  ! ========================================================================
  

  subroutine LandUse_BCatIon_Init( rcF, status )
  
    use GO  , only : TrcFile
    use Dims, only : nx, ny
    use Indices, only : n_basecation
  
    ! --- in/out ---------------------------------
    
    type(TrcFile), intent(in)   ::  rcF
    integer, intent(out)        ::  status

    ! --- const --------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_LandUse_BCatIon_Init'
    
    ! --- local ------------------------------------------
    
    ! --- begin ------------------------------------------
    
   ! ! enabled ? only Na is not enough ...Adaptation for additional Na classes, should be done with better check!!!
    if ( n_basecation > 5 ) then

      ! storage:
      allocate( bcat_frac(nx,ny,nbcat) )

      ! fill:
      call get_soil_composition( status )
      IF_NOTOK_RETURN(status=1)
      
    end if  ! enabled

    ! ok
    status = 0
    
  end subroutine LandUse_BCatIon_Init
  
  
  ! ***
  
  
  subroutine LandUse_BCatIon_Done( status )

    use Indices, only : n_basecation
  
    ! --- in/out ---------------------------------
    
    integer, intent(out)      ::  status
    
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_LandUse_BCatIon_Done'
    
    ! --- begin ----------------------------------
    
    ! enabled ? only Na is not enough ...
    if ( n_basecation > 5 ) then

      ! clear:
      deallocate( bcat_frac )
      
    end if  ! enabled

    ! ok
    status = 0
    
  end subroutine LandUse_BCatIon_Done
  
  
  ! ***


  subroutine get_soil_composition( status )

    use dims, only : nx, ny

    ! --- in/out ---------------------------------

    integer, intent(out)    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   :: rname = mname//'/get_tree_data'

    ! the full dimensions are defined by the size of
    ! the meteorological input which is given on a
    ! resolution of 0.5x0.25 degree lon/lat
    integer, parameter :: nx_met=140, ny_met=140

    ! for temporary files that are closed directly after use
    integer, parameter :: u_tmp  = 61
   
    ! --- local ----------------------------------

    integer     ::  k, zrec
    real        ::  sc(nx_met, ny_met, nbcat)
    integer     ::  reclen_fac_bin_io

    ! --- begin ----------------------------------

    print *, "do something here on interpol and name"
    stop

    ! Initialise the record length of a 4-byte real value
    ! in case of binary i/o (grads files).
    inquire (iolength=reclen_fac_bin_io) real(1.0,kind=4)

    open (u_tmp,file='ca_content.dat',form='unformatted',recl=nx_met*ny_met*reclen_fac_bin_io, &
        access='direct',status='old')
    zrec = 1 ! skip first useless record
    do k=1,nbcat
      zrec = zrec+1
      print *, zrec
      read (u_tmp, rec=zrec) sc(:,:,k)
      ! no data is -999, thus put them to zero
      where(sc(:,:,k).lt.0.0) sc(:,:,k) = 0.0
      !MSP call interpol( sc(:,:,k), bcat_frac(:,:,k), 'constant') ! selects the area used in this simulation and interpolates, if necessary
    enddo
    ! convert bcatfrac from % to fraction:
    bcat_frac = bcat_frac / 100.

    close(u_tmp)

    ! ok
    status = 0

  end subroutine get_soil_composition


end module LE_LandUse_BCatIon
