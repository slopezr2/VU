!###############################################################################
!
! NAME
!
!   LE_Output_Common  -  LOTOS-EUROS output tools and global variables
!
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Output_Common

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

#ifdef with_netcdf  
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  implicit none


  ! --- in/out -----------------------------
  
  public
  
  public  ::  T_LE_Output_Common


  ! --- const ------------------------------
    
  character(len=*), parameter   ::  mname = 'LE_Output_Common'
  
    
  ! --- types ------------------------------
  
  type T_LE_Output_Common
    ! output directory:
    character(len=256)          ::  outdir
    ! description:
    character(len=64)           ::  author
    character(len=64)           ::  institution
    character(len=64)           ::  version
    ! experiment:
    character(len=32)           ::  model
    character(len=64)           ::  expid
    ! CF convention:
    character(len=32)           ::  CF_convention
    ! t0
    type(Tdate)                 ::  t0
    ! replace existing files?
    logical                     ::  replace
    ! deflate level:
    integer                     ::  deflate_level
  contains
    procedure :: Init                       =>  LE_Output_Common_Init
    procedure :: Done                       =>  LE_Output_Common_Done
    procedure :: PutOut_GlobalAttributes    =>  LE_Output_Common_PutOut_GlobalAttributes
  end type T_LE_Output_Common

  
contains


  ! ====================================================
  
  
  subroutine LE_Output_Common_Init( self, rcF, rckey, status )

    use GO     , only : TrcFile
    use GO     , only : TDate
    use GO     , only : NewDate
    use GO     , only : GoTranslate
    
    use LE_Config, only : outputdir
  
    ! --- in/out --------------------------------
    
    class(T_LE_Output_Common), intent(out)  ::  self
    type(TrcFile), intent(in)               ::  rcF
    character(len=*), intent(in)            ::  rckey
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Common_Init'
    
    ! --- local ---------------------------------
    
    character(len=32)     ::  stime
    integer               ::  minu, seco
    integer               ::  yy_s, mm_s, dd_s, hh_s

    ! --- begin ---------------------------------
    
    ! copy output directory:
    self%outdir = trim(outputdir)
    
    ! description:
    call rcF%Get( trim(rckey)//'.author'     , self%author     , status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( trim(rckey)//'.institution', self%institution, status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( trim(rckey)//'.version'    , self%version    , status )
    IF_NOTOK_RETURN(status=1)

    ! read keys for filename etc:
    call rcF%Get( trim(rckey)//'.model', self%model, status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( trim(rckey)//'.expid', self%expid, status )
    IF_NOTOK_RETURN(status=1)

    ! CF convention: "CF-1.5"
    call rcF%Get( 'cf.convention', self%CF_convention, status )
    IF_NOTOK_RETURN(status=1)

    ! reference time for out put file (seconds from January first of the start year of simulation)   
    call rcF%Get( 'timerange.start' , stime, status )
    IF_NOTOK_RETURN(status=1)
    call goTranslate( stime, '/-:', ' ', status )
    IF_NOTOK_RETURN(status=1)
    read (stime,*,iostat=status) yy_s, mm_s, dd_s, hh_s, minu, seco
    if ( status /= 0 ) then
      write (gol,'("could not read start time from : ",a)') trim(stime); call goErr
      TRACEBACK; status=1; return
    end if
    self%t0 = NewDate( time6=(/yy_s,01,01,00,00,00/) )

    ! replace existing files?
    call rcF%Get( trim(rckey)//'.replace', self%replace, status )
    IF_NOTOK_RETURN(status=1)
    
    ! deflate level:
    call rcF%Get( trim(rckey)//'.deflate_level', self%deflate_level, status )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LE_Output_Common_Init
  
  
  ! ***
  

  subroutine LE_Output_Common_Done( self, status )
  
    ! --- in/out --------------------------------
    
    class(T_LE_Output_Common), intent(inout)  ::  self
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Common_Done'
    
    ! --- begin ---------------------------------
    
    ! dummy to avoid warnings about unused variables ...
    self%outdir = ''
    
    ! ok
    status = 0
    
  end subroutine LE_Output_Common_Done


  ! ====================================================
  
  ! for backwards compatibility ...  
  
  subroutine leoc_common_Init( leoc, rcF, rckey, status )

    use GO     , only : TrcFile
    use GO     , only : TDate
    use GO     , only : NewDate
    use GO     , only : GoTranslate
    
    use LE_Config, only : outputdir
  
    ! --- in/out --------------------------------
    
    type(T_LE_Output_Common), intent(out)   ::  leoc
    type(TrcFile), intent(in)               ::  rcF
    character(len=*), intent(in)            ::  rckey
    integer, intent(out)                    ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/leoc_common_Init'
    
    ! --- local ---------------------------------
    
    character(len=32)     ::  stime
    integer               ::  minu, seco
    integer               ::  yy_s, mm_s, dd_s, hh_s

    ! --- begin ---------------------------------
    
    ! copy output directory:
    leoc%outdir = trim(outputdir)
    
    ! description:
    call rcF%Get( trim(rckey)//'.author'     , leoc%author     , status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( trim(rckey)//'.institution', leoc%institution, status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( trim(rckey)//'.version'    , leoc%version    , status )
    IF_NOTOK_RETURN(status=1)

    ! read keys for filename etc:
    call rcF%Get( trim(rckey)//'.model', leoc%model, status )
    IF_NOTOK_RETURN(status=1)
    call rcF%Get( trim(rckey)//'.expid', leoc%expid, status )
    IF_NOTOK_RETURN(status=1)

    ! CF convention: "CF-1.5"
    call rcF%Get( 'cf.convention', leoc%CF_convention, status )
    IF_NOTOK_RETURN(status=1)

    ! reference time for out put file (seconds from January first of the start year of simulation)   
    call rcF%Get( 'timerange.start' , stime, status )
    IF_NOTOK_RETURN(status=1)
    call goTranslate( stime, '/-:', ' ', status )
    IF_NOTOK_RETURN(status=1)
    read (stime,*,iostat=status) yy_s, mm_s, dd_s, hh_s, minu, seco
    if ( status /= 0 ) then
      write (gol,'("could not read start time from : ",a)') trim(stime); call goErr
      TRACEBACK; status=1; return
    end if
    leoc%t0 = NewDate( time6=(/yy_s,01,01,00,00,00/) )
    
    ! ok
    status = 0
    
  end subroutine leoc_common_Init
  
  
  ! ***
  

  subroutine leoc_common_Done( leoc, status )
  
    ! --- in/out --------------------------------
    
    type(T_LE_Output_Common), intent(inout)   ::  leoc
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/leoc_common_Done'
    
    ! --- begin ---------------------------------
    
    ! dummy to avoid warnings about unused variables ...
    leoc%outdir = ''
    
    ! ok
    status = 0
    
  end subroutine leoc_common_Done
  
  
  ! ***
  
  
  subroutine LE_Output_Common_PutOut_GlobalAttributes( leoc, ncid, status )

#ifdef with_netcdf  
    use NetCDF , only : NF90_Put_Att, NF90_GLOBAL
#endif
  
    ! --- in/out ------------------------------
    
    class(T_LE_Output_Common), intent(in)     ::  leoc
    integer, intent(in)                       ::  ncid
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/LE_Output_Common_PutOut_GlobalAttributes'
    
    ! --- begin ---------------------------------
    
#ifndef with_netcdf
    ! check ...
    stop 'not compiled with netcdf support'
#endif
    
#ifdef with_netcdf  
    ! data origin:    
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'author', trim(leoc%author) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'institution', trim(leoc%institution) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'version', trim(leoc%version) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! model and experiment info:
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'model', trim(leoc%model) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'expid', trim(leoc%expid) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! CF convention
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'Conventions', trim(leoc%CF_convention) )
    IF_NF90_NOTOK_RETURN(status=1)
#else
    ! dummy to avoid warnings about unused variables ...
    print *, 'LE_Output_Common_PutOut_GlobalAttributes: ', ncid, trim(leoc%model)
#endif

    ! ok
    status = 0
    
  end subroutine LE_Output_Common_PutOut_GlobalAttributes

  
end module le_output_common
