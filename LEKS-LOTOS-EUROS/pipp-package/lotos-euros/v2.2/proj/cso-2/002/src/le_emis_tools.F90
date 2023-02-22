!###############################################################################
!
! LE_Emis_Tools
!
! Conversion routines etc.
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

module LE_Emis_Tools

  use GO, only : gol, goErr, goPr
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  emis_a_bio
  
  public  ::  Get_ShortSNAP
  public  ::  ShortSNAP_to_Code
  public  ::  MDF_Get_StrArr


  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_Tools'
  
  
  ! --- var ----------------------------------

  ! total bio emis, for output ..
  real, allocatable     ::  emis_a_bio(:,:,:,:)
  
  
contains



  ! ===============================================================
  ! ===
  ! === snap conversions
  ! ===
  ! ===============================================================
  
  !
  ! Conversions:
  !   '1'      ->   1
  !   '07.01'  ->  71
  
  subroutine Get_ShortSNAP( field, shortsnap, status )
  
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  field
    integer, intent(out)              ::  shortsnap
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/Get_ShortSNAP'
    
    ! --- local ----------------------------------
    
    character(len=256)      ::  line
    integer                 ::  snap0, snap1

    ! --- begin ----------------------------------
    
    ! main and sub ? then dot included:
    if ( index(field,'.') > 0 ) then
      ! copy:
      line = trim(field)
      ! extract main number:
      call goReadFromLine( line, snap0, status, sep='.' )
      IF_NOTOK_RETURN(status=1)
      ! extract sub number:
      call goReadFromLine( line, snap1, status )
      IF_NOTOK_RETURN(status=1)
      ! convert:
      if ( snap1 == 0 ) then
        shortsnap = snap0
      else
        shortsnap = snap0 * 10 + snap1
      end if
    !~ single number
    else
      ! read number directly:
      read (field,'(i2)',iostat=status) shortsnap
      if (status/=0) then
        write (gol,'("reading short snap code from `",a,"`")') trim(field); call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! ok
    status = 0
    
  end subroutine Get_ShortSNAP
  
  !
  ! Conversions:
  !    1 -> '01:00'
  !   71 -> '07.02'
  !
  
  subroutine ShortSNAP_to_Code( shortsnap, code, status )
  
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    integer, intent(in)               ::  shortsnap
    character(len=*), intent(out)     ::  code
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter ::  rname = mname//'/ShortSNAP_to_Code'
    
    ! --- local ----------------------------------
    
    integer                 ::  snap0, snap1

    ! --- begin ----------------------------------
    
    ! get snap levels:
    if ( (shortsnap >= 70) .and. (shortsnap <= 79) ) then
      snap0 = 7
      snap1 = shortsnap - 70
    else
      snap0 = shortsnap
      snap1 = 0
    end if
    
    ! fill code:
    write (code,'(i2.2,".",i2.2)') snap0, snap1
    
    ! ok
    status = 0
    
  end subroutine ShortSNAP_to_Code
  
  
  ! ***


  ! Read 2d character variable in a variable with oversized length

  subroutine MDF_Get_StrArr( hid, varname, nchar, nval, values, status )
  
    use MDF, only : MDF_Inq_VarID, MDF_Get_Var

    ! --- in/out -----------------------

    integer, intent(in)             ::  hid
    character(len=*), intent(in)    ::  varname
    integer, intent(in)             ::  nchar, nval
    character(len=*), intent(out)   ::  values(nval)
    integer, intent(out)            ::  status

    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/MDF_Get_StrArr'
    
    ! --- local -------------------------

    integer                 ::  varid
    character(len=nchar)    ::  tmp_values(nval)
    integer                 ::  ival
    integer                 ::  l

    ! --- begin -------------------------

    ! get variable id given name:      
    call MDF_Inq_VarID( hid, varname, varid, status )
    IF_NOTOK_RETURN(status=1)
    
    ! read data:
    call MDF_Get_Var( hid, varid, tmp_values, status )
    IF_NOTOK_RETURN(status=1)
    
    ! copy:
    do ival = 1, nval
      ! seach <nul> characters:
      l = index( tmp_values(ival), char(0) )
      ! found a <nul> character ?
      if ( l < 1 ) then
        ! no <nul>, just copy:
        values(ival) = trim(tmp_values(ival))
      else if ( l == 1 ) then
        ! first character is <nul>, thus empty:
        values(ival) = ''
      else
        ! found a <nul> somewhere in the string; copy part before:
        values(ival) = trim(tmp_values(ival)(1:l-1))
      end if
    end do

    ! ok
    status = 0

  end subroutine MDF_Get_StrArr
    

end module LE_Emis_Tools
