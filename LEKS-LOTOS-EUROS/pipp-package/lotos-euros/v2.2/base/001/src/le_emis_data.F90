!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Emis_Data

  use GO, only : gol, goPr, goErr

  use LE_Data_Variable, only : LEN_VARNAME

  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  T_Emis_Data


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Data'
  

  ! --- types --------------------------------

  ! emission data base
  type T_Emis_Data
    ! list with data variables:
    integer                             ::  nvar
    character(len=LEN_VARNAME)          ::  varnames(10)   ! just a maximum ...
    ! target spec indices:
    integer, allocatable                ::  ispecs(:)   ! (nvar)
  contains
    procedure :: Init          => Emis_Data_Init
    procedure :: Done          => Emis_Data_Done
    procedure :: Setup         => Emis_Data_Setup
  end type T_Emis_Data



contains


  ! ====================================================================
  ! ===
  ! === Emis Data
  ! ===
  ! ====================================================================


  subroutine Emis_Data_Init( self, rcF, rckey, status )
  
    use GO     , only : TrcFile, ReadRc
    use GO     , only : goSplitString, goMatchValue
    use LE_Data, only : LE_Data_Enable
    use Indices, only : specname
    
    ! --- in/out ---------------------------------
    
    class(T_Emis_Data), intent(out)           ::  self
    type(TrcFile), intent(in)                 ::  rcF
    character(len=*), intent(in)              ::  rckey
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Emis_Data_Init'
      
    ! --- local ----------------------------------
    
    character(len=1024)           ::  line
    integer                       ::  ivar
    character(len=32)             ::  spec
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("setup emissions from data variables ...")'); call goPr
    
    ! read line with supported variables:
    call ReadRc( rcF, trim(rckey)//'.vars', line, status )
    IF_NOTOK_RETURN(status=1)
    ! split:
    call goSplitString( line, self%nvar, self%varnames, status )
    IF_NOTOK_RETURN(status=1)
    
    ! storage:
    allocate( self%ispecs(self%nvar), stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over variables:
    do ivar = 1, self%nvar

      ! enable data:
      call LE_Data_Enable( trim(self%varnames(ivar)), status )
      IF_NOTOK_RETURN(status=1)

      ! target spec:
      call ReadRc( rcF, trim(rckey)//'.'//trim(self%varnames(ivar))//'.spec', spec, status )
      IF_NOTOK_RETURN(status=1)
      ! global index:
      call goMatchValue( spec , specname, self%ispecs(ivar), status )
      IF_NOTOK_RETURN(status=1)

    end do
    
    ! ok
    status = 0
    
  end subroutine Emis_Data_Init


  ! ***


  subroutine Emis_Data_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Emis_Data), intent(inout)         ::  self
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Emis_Data_Done'
      
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%ispecs, stat=status )
    IF_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Emis_Data_Done


  ! ***


  subroutine Emis_Data_Setup( self, t1, t2, status )
  
    use GO     , only : TDate
    use LE_Data, only : LE_Data_GetPointer
    use Dims   , only : nx, ny
    use dims   , only : emis_a
    use Indices, only : specname, specunit, specmolm
  
    ! --- in/out ---------------------------------
    
    class(T_Emis_Data), intent(inout)         ::  self
    type(TDate), intent(in)                   ::  t1, t2
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Emis_Data_Setup'
    
    ! conversions:
    real, parameter   ::  sec_per_min = 60.0
      
    ! --- local ----------------------------------
    
    integer             ::  ivar
    character(len=64)   ::  units
    integer             ::  ispec
    real, pointer       ::  area(:,:,:)
    real, pointer       ::  emis(:,:,:)
    character(len=256)  ::  conversion
    
    ! --- begin ----------------------------------

    ! grid cell area:
    call LE_Data_GetPointer( 'area', area, status, check_units='m2' )
    IF_NOTOK_RETURN(status=1)
    
    ! loop over variables:
    do ivar = 1, self%nvar
    
      ! info ...
      write (gol,'(" add `",a,"` emissions ...")') trim(self%varnames(ivar)); call goPr

      ! pointer to emission array:
      call LE_Data_GetPointer( trim(self%varnames(ivar)), emis, status, &
                                 units=units )
      IF_NOTOK_RETURN(status=1)
      
      ! index of target spec:
      ispec = self%ispecs(ivar)
      
      ! required conversion:
      select case ( specunit(ispec) )
        !~ emission for tracers expressed as volume mixing ratio's:
        case ( 'ppb' )
          ! unit emis_a for gases   : mol/min
          conversion = trim(units)//' -> mol/min'
        !~ unknown ...
        case default
          write (gol,'("unsupported unit `",a,"` for tracer ",i4," (",a,")")') &
                           trim(specunit(ispec)), ispec, trim(specname(ispec)); call goErr
          TRACEBACK; status=1; return
      end select

      ! fill emission array:
      select case ( trim(conversion) )

        !~ emission of trace gasses:
        case ( 'kg/m2/s -> mol/min' )
          !       mol/min                     mol/min               kg/m2/s   *         m2        /     (kg/mol)    *  (s/min)
          emis_a(1:nx,1:ny,1,ispec) = emis_a(1:nx,1:ny,1,ispec) + emis(:,:,1) * area(1:nx,1:ny,1) / specmolm(ispec) * sec_per_min

        !~ unknown ...
        case default
          write (gol,'("unsupported conversion `",a,"` for tracer ",i4," (",a,") emitted from `",a,"`")') &
                           trim(conversion), ispec, trim(specname(ispec)), &
                           trim(self%varnames(ivar)); call goErr
          TRACEBACK; status=1; return
      end select

    end do
    
    ! ok
    status = 0
    
  end subroutine Emis_Data_Setup


end module LE_Emis_Data

