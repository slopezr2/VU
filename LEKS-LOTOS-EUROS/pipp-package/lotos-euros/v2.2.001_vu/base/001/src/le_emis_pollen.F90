!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_Pollen

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile
  
  use Pollen_Data, only : T_Pollen_Data

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Emis_Pollen

  public  ::  LE_Emis_Pollen_Init, LE_Emis_Pollen_Done
  public  ::  LE_Emis_Pollen_Setup


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Pollen'

  ! maximum number of emitted tracers:
  integer, parameter            ::  max_emis = 10

  ! --- local --------------------------------

  ! emission data base
  type T_Emis_Pollen
    ! label assigned to this emission:
    character(len=32)                 ::  label
    
    ! number of emitted tracers
    integer                           ::  ntr
    
    ! input file:
    type(T_Pollen_Data)               ::  emp(max_emis)
    
    ! use Pollen landuse ?
    logical                           ::  with_pol_landuse
    ! scaling method:
    character(len=16)                 ::  scaling
  
  end type T_Emis_Pollen

contains


  ! ===============================================================


  subroutine LE_Emis_Pollen_Init( ems, label, rcF, rckey, t, status )

    use GO         , only : TrcFile, ReadRc
    use GO         , only : goMatchValue, goSplitString
    use GO         , only : TDate
    use Pollen_Data, only : Pollen_Data_Init
    use indices    , only : specname

    ! --- in/out ------------------------------

    type(T_Emis_Pollen), intent(out)    ::  ems
    character(len=*), intent(in)        ::  label
    type(TrcFile), intent(in)           ::  rcF
    character(len=*), intent(in)        ::  rckey
    type(TDate), intent(in)             ::  t
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Pollen_Init'

    ! --- local ------------------------------
    
    character(len=1024)         ::  fname_template
    character(len=512)          ::  list
    character(len=32)           ::  tracer_names(max_emis)
    integer                     ::  itr

    ! --- begin -------------------------------

    ! store label:
    ems%label = trim(label)
             
    ! read name of tracer:  IN FUTURE MULTIPLE !
    call ReadRc( rcF, trim(rckey)//'.tracers', list, status )
    IF_NOTOK_RETURN(status=1)    
    
    ! target tracer index
    call goSplitString( list, ems%ntr, tracer_names,  status )
    IF_NOTOK_RETURN(status=1)
        
    do itr = 1, ems%ntr
      ! set tracer name in pollen data type
      ems%emp(itr)%tracer = trim(tracer_names(itr) )
      
      ! match tracer name to global tracer index
      call goMatchValue( trim(tracer_names(itr) ), specname, ems%emp(itr)%itr_glob, status )
      IF_NOTOK_RETURN(status=1)
      
      ! read in tracer dependent emission information
      call Pollen_Data_Init(ems%emp(itr), rcF, trim(rckey), t, status )
      IF_NOTOK_RETURN(status=1)
    end do
        
    ! ok
    status = 0

  end subroutine LE_Emis_Pollen_Init


  ! ***


  subroutine LE_Emis_Pollen_Done( ems, status )
    
    use Pollen_Data, only : Pollen_Data_Done
    
    ! --- in/out ------------------------------

    type(T_Emis_Pollen), intent(inout)    ::  ems
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Pollen_Done'

    ! --- local -------------------------------
    
    integer                       ::  itr
    
    ! --- begin -------------------------------
    
    do itr = 1, ems%ntr 
      ! close file:
      call Pollen_Data_Done( ems%emp(itr), status )
      IF_NOTOK_RETURN(status=1)
    end do
	
    ! ok
    status = 0

  end subroutine LE_Emis_Pollen_Done


  ! ***


  subroutine LE_Emis_Pollen_Setup( ems, emis_a, t1, t2, status )

    use GO          , only : TDate
    use Dims        , only : nx, ny, nz, nspec
    use Pollen_Data , only : Pollen_Data_Get    


 
    ! --- in/out ---------------------------

    type(T_Emis_Pollen), intent(inout)    ::  ems
    real, intent(inout)                   ::  emis_a(nx,ny,nz,nspec)
    type(TDate), intent(in)               ::  t1, t2
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/LE_Emis_Pollen_Setup'
    
    ! --- local ----------------------------
    
    integer                             ::  itr
    integer                             ::  ix,iy
    real                                ::  delta_emis_a
     
    ! --- begin ----------------------------
     
    ! is pollen in the model anyway ?
    do itr = 1, ems%ntr
      
      call Pollen_Data_Get( ems%emp(itr), t1, t2, status )
      IF_NOTOK_RETURN(status=1)
      
      emis_a(:,:,1,ems%emp(itr)%itr_glob) = ems%emp(itr)%emis
    end do 
    
    ! ok
    status = 0

  end subroutine LE_Emis_Pollen_Setup
  
  
  ! ***
 
end module LE_Emis_Pollen


 
