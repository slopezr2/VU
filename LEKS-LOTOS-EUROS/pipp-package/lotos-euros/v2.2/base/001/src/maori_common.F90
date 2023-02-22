!#######################################################################
!
! NAME
!
!   MAORI - Model And Output Routine Interface
!
! DESCRIPTION
!
!   Mode independent output routines.
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "mdf.inc"
!
!#######################################################################

module MAORI_Common

  use GO, only : gol, goPr, goErr
  use GO, only : TDate
  
  use MAORI_Param, only : MAORI_LEN_FILE, MAORI_LEN_NAME, MAORI_LEN_LONGNAME

#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif

  implicit none


  ! --- in/out ---------------------------

  private

  public  ::  MAORI_Common_Init, MAORI_Common_Done
#ifdef with_netcdf
  public  ::  MAORI_Common_NF90_GlobalAttributes
#endif

  public  ::  com_path
  public  ::  com_model_name
  public  ::  com_model_version
  public  ::  com_author_name
  public  ::  com_author_inst
  public  ::  com_exp_id
  public  ::  com_data_version
  public  ::  com_t0


  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'MAORI_Common'


  ! --- types ------------------------------
  
  character(len=MAORI_LEN_FILE)       ::  com_path
  character(len=MAORI_LEN_NAME)       ::  com_model_name
  character(len=MAORI_LEN_NAME)       ::  com_model_version
  character(len=MAORI_LEN_LONGNAME)   ::  com_author_name
  character(len=MAORI_LEN_LONGNAME)   ::  com_author_inst
  character(len=MAORI_LEN_LONGNAME)   ::  com_exp_id
  character(len=MAORI_LEN_LONGNAME)   ::  com_data_version
  type(TDate)                         ::  com_t0



contains


  ! ==================================================================


  subroutine MAORI_Common_Init( rcF, status )

    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine
    use GO, only : NewDate, wrtgol
    
    ! --- in/out ---------------------------------

    type(TRcFile), intent(in)         ::  rcF
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Common_Init'

    ! --- local -----------------------------------
    
    character(len=20)   ::  sdatetime
    integer             ::  year, month, day, hour, minu, sec

    ! --- begin -----------------------------------

    ! read output path:
    call ReadRc( rcF, 'maori.path', com_path, status )
    IF_NOTOK_RETURN(status=1)

    ! read author and institute:
    call ReadRc( rcF, 'maori.author.name', com_author_name, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'maori.author.inst', com_author_inst, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read model name etc:
    call ReadRc( rcF, 'maori.model.name', com_model_name, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'maori.model.version', com_model_version, status )
    IF_NOTOK_RETURN(status=1)
    
    ! experiment:
    call ReadRc( rcF, 'maori.experiment.id', com_exp_id, status )
    IF_NOTOK_RETURN(status=1)
    
    ! data:
    call ReadRc( rcF, 'maori.data.version', com_data_version, status )
    IF_NOTOK_RETURN(status=1)

    ! times computed as days since t0
    ! ~ read date as 'yyyy-dd-mm hh:mm:ss'
    call ReadRc( rcF, 'maori.time.t0', sdatetime, status )
    IF_NOTOK_RETURN(status=1)
    ! ~ extract:
    call goReadFromLine( sdatetime, year, status, sep='-', default=0 )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( sdatetime, month, status, sep='-', default=0 )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( sdatetime, day, status, sep=' ', default=0 )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( sdatetime, hour, status, sep=':', default=0 )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( sdatetime, minu, status, sep=':', default=0 )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( sdatetime, sec, status, default=0 )
    IF_NOTOK_RETURN(status=1)
    ! ~ fill:
    com_t0 = NewDate( year=year, month=month, day=day, hour=hour, min=minu, sec=sec )
    
    ! ok
    status = 0

  end subroutine MAORI_Common_Init


  ! ***


  subroutine MAORI_Common_Done( status )

    ! --- in/out ---------------------------------

    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/MAORI_Common_Done'

    ! --- begin ----------------------------------
    
    ! ok
    status = 0

  end subroutine MAORI_Common_Done


  ! ***
  
  
#ifdef with_netcdf  

  subroutine MAORI_Common_NF90_GlobalAttributes( ncid, status )

    use NetCDF , only : NF90_Put_Att, NF90_GLOBAL
  
    ! --- in/out ------------------------------
    
    integer, intent(in)                       ::  ncid
    integer, intent(out)                      ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MAORI_Common_NF90_GlobalAttributes'
    
    ! --- begin ---------------------------------
    
    ! data origin:    
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'author_name', trim(com_author_name) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'author_institution', trim(com_author_inst) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! model info:
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'model_name'   , trim(com_model_name) )
    IF_NF90_NOTOK_RETURN(status=1)
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'model_version', trim(com_model_version) )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! experiment:
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'experiment_id', trim(com_exp_id) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! data version:
    status = NF90_Put_Att( ncid, NF90_GLOBAL, 'data_version', trim(com_data_version) )
    IF_NF90_NOTOK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine MAORI_Common_NF90_GlobalAttributes

#endif

  
end module MAORI_Common
