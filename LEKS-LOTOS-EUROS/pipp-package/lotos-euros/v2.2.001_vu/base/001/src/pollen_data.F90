!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
!###############################################################################

module Pollen_Data

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile
  use GO, only : TDate
  
  use Pollen_Data_Birch, only : T_Pollen_Data_Birch
  use Pollen_Data_Olive, only : T_Pollen_Data_Olive
  use Pollen_Data_Grass, only : T_Pollen_Data_Grass
 
  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Pollen_Data

  public  ::  Pollen_Data_Init, Pollen_Data_Done
  public  ::  Pollen_Data_Get

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'Pollen_Data'
  

  ! --- local --------------------------------

  ! emission data base
  type T_Pollen_Data

    character(len=32)                 ::  tracer
    integer                           ::  itr_glob
    
    type(T_Pollen_Data_Birch)         ::  Birch_pollen
    type(T_Pollen_Data_Olive)         ::  Olive_pollen
    type(T_Pollen_Data_Grass)         ::  Grass_pollen

    real, allocatable                 ::  emis(:,:)
    
    
  end type T_Pollen_Data



contains


  ! ===============================================================


  subroutine Pollen_Data_Init( emp, rcF, rckey, t, status )
    
    use dims, only : nx, ny

    use Pollen_Data_Birch, only : Pollen_Data_Birch_Init
    use Pollen_Data_Olive, only : Pollen_Data_Olive_Init
    use Pollen_Data_Grass, only : Pollen_Data_Grass_Init

    ! --- in/out ------------------------------

    type(T_Pollen_Data), intent(inout)  ::  emp
    type(TrcFile), intent(in)           ::  rcF
    character(len=*), intent(in)        ::  rckey
    type(TDate), intent(in)             ::  t
    integer, intent(out)                ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Init'

    ! --- local ------------------------------

    ! --- begin -------------------------------
    
    allocate ( emp%emis(nx,ny) )
    
    if ( emp%tracer == 'pol_b' ) then
      emp%Birch_pollen%tracer = 'pol_b' 
      call Pollen_Data_Birch_Init( emp%Birch_pollen, rcF, rckey, t, status  )
      IF_NOTOK_RETURN(status=1)

    else if ( emp%tracer == 'pol_o' ) then
      emp%Olive_Pollen%tracer = 'pol_o'
      call Pollen_Data_Olive_Init( emp%Olive_pollen, rcF, rckey, t, status  )
      IF_NOTOK_RETURN(status=1)

    else if ( emp%tracer == 'pol_g' ) then
      emp%Grass_Pollen%tracer = 'pol_g' 
      call Pollen_Data_Grass_Init( emp%Grass_pollen, rcF, rckey, t, status  )
      IF_NOTOK_RETURN(status=1)

    else 
      write( gol, '("Unknown pollen tracer: ", a)' ) emp%tracer ; call goErr
      TRACEBACK;status=1;return
    end if
    
    ! ok
    status = 0

  end subroutine Pollen_Data_Init

  
  ! ***


  subroutine Pollen_Data_Done( emp, status )

    use Pollen_Data_Birch, only : Pollen_Data_Birch_Done
    use Pollen_Data_Olive, only : Pollen_Data_Olive_Done
    use Pollen_Data_Grass, only : Pollen_Data_Grass_Done

    ! --- in/out ------------------------------

    type(T_Pollen_Data), intent(inout)    ::  emp
    integer, intent(out)                  ::  status

    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Done'

    ! --- local -------------------------------

    ! --- begin -------------------------------
    
    deallocate( emp%emis )
    
    ! clear ?
    if ( emp%tracer == 'pol_b' ) then
      call Pollen_Data_Birch_Done( emp%Birch_pollen, status  )
      IF_NOTOK_RETURN(status=1)
    
    else if ( emp%tracer == 'pol_o' ) then
      call Pollen_Data_Olive_Done( emp%Olive_pollen, status  )
      IF_NOTOK_RETURN(status=1)

    else if ( emp%tracer == 'pol_g' ) then
      call Pollen_Data_Grass_Done( emp%Grass_pollen, status  )
      IF_NOTOK_RETURN(status=1)

    else 
      write( gol, '("Unknown pollen tracer: ", a)' ) emp%tracer ; call goErr
      TRACEBACK;status=1;return
    end if
    
    !  ok
    status = 0

  end subroutine Pollen_Data_Done


  ! ***

  
  subroutine Pollen_Data_Get( emp, t1, t2, status )
    
    use GO, only : TDate
    use Pollen_Data_Birch, only : Pollen_Data_Birch_Get
    use Pollen_Data_Olive, only : Pollen_Data_Olive_Get
    use Pollen_Data_Grass, only : Pollen_Data_Grass_Get

    
    ! --- in/out ---------------------------

    type(T_Pollen_Data), intent(inout)    ::  emp
    type(TDate), intent(in)               ::  t1, t2   
    integer, intent(out)                  ::  status
    
    ! --- const -------------------------------

    character(len=*), parameter ::  rname = mname//'/Pollen_Data_Get'
    
    ! --- local ----------------------------
        
    ! --- begin ----------------------------
                                      
    ! 
    if ( emp%tracer == 'pol_b' ) then
      call Pollen_Data_Birch_Get( emp%Birch_pollen, t1, t2, status  )
      IF_NOTOK_RETURN(status=1)
      
      emp%emis = emp%Birch_pollen%emis

    else if ( emp%tracer == 'pol_o' ) then
      call Pollen_Data_Olive_Get( emp%Olive_pollen, t1, t2, status  )
      IF_NOTOK_RETURN(status=1)

      emp%emis = emp%Olive_pollen%emis
      
    else if ( emp%tracer == 'pol_g' ) then
      call Pollen_Data_Grass_Get( emp%Grass_pollen, t1, t2, status  )
      IF_NOTOK_RETURN(status=1)

      emp%emis = emp%Grass_pollen%emis

    else 
      write( gol, '("Unknown pollen tracer: ", a)' ) emp%tracer ; call goErr
      TRACEBACK;status=1;return
    end if
    
     ! ok
    status = 0

  end subroutine Pollen_Data_Get     
  
end module Pollen_Data
