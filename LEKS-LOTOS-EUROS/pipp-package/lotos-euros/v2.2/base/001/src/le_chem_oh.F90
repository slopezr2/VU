!###############################################################################
!
! LE_Chem_OH - OH fields for use in sulphur-only or methane-only chemistries
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

module LE_Chem_OH

  use GO, only : gol, goErr, goPr
  
#ifdef with_netcdf
  use NetCDF, only : NF90_StrError, NF90_NOERR
#endif  
  
  implicit none
  
  ! --- in/out --------------------------------
  
  private
  
  public  ::  LE_Chem_OH_Init, LE_Chem_OH_Done
  public  ::  LE_Chem_OH_Setup
  

  ! --- const --------------------------------

  ! module name:
  character(len=*), parameter ::  mname = 'LE_Chem_OH'


  ! --- local ---------------------------------
  
  ! OH archive description:    
  character(len=512)        ::  arch_path
  character(len=512)        ::  arch_format
  character(len=512)        ::  arch_type
  character(len=32)         ::  arch_model
  character(len=32)         ::  arch_runid


contains


  ! ===============================================================
  

  subroutine LE_Chem_OH_Init( rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goVarValue

    ! --- in/out ------------------------------
    
    type(TrcFile), intent(in)       ::  rcF
    integer, intent(out)            ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Chem_OH_Init'
    
    ! --- local ------------------------------
    
    character(len=512)    ::  key
    
    ! --- begin -------------------------------

    ! OH archive description:
    call ReadRc( rcF, 'le.chem.OH.path', arch_path, status )
    IF_NOTOK_RETURN(status=1)
    call ReadRc( rcF, 'le.chem.OH.key', key, status )
    IF_NOTOK_RETURN(status=1)

    ! extract parts from key:
    !
    !   format=netcdf;model=LE;runid=base;type=conc
    !
    arch_format = 'netcdf'
      call goVarValue( key, ';', 'format' , '=', arch_format, status )
      IF_ERROR_RETURN(status=1)
    arch_model = 'LE'
      call goVarValue( key, ';', 'model', '=', arch_model, status )
      IF_ERROR_RETURN(status=1)
    arch_runid = 'base'
      call goVarValue( key, ';', 'runid', '=', arch_runid, status )
      IF_ERROR_RETURN(status=1)
    arch_type = 'conc'
      call goVarValue( key, ';', 'type' , '=', arch_type , status )
      IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine LE_Chem_OH_Init
  
  
  ! ***
  

  subroutine LE_Chem_OH_Done( status )
  
    ! --- in/out ------------------------------
    
    integer, intent(out)  ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Chem_OH_Done'
    
    ! --- local -------------------------------
    
    ! --- begin -------------------------------

    ! ok
    status = 0

  end subroutine LE_Chem_OH_Done
  
  
  ! ***
  
  
  subroutine LE_Chem_OH_Setup( t, status )
  
    use GO, only : TDate, NewDate, wrtgol, operator(/=)

    use dims, only : nx, ny, nz, nspec
    use dims, only : OH
    
#ifdef with_netcdf
    use NetCDF, only : NF90_NoWrite
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
#endif

    ! --- in/out ------------------------------
    
    type(TDate), intent(in)   ::  t
    integer, intent(out)      ::  status
  
    ! --- const -------------------------------
    
    character(len=*), parameter ::  rname = mname//'/LE_Chem_OH_Setup'
    
    ! --- local -------------------------------
    
    character(len=512)    ::  fname
    type(TDate)           ::  tday
    logical               ::  exist
    integer               ::  ncid
    integer               ::  varid
    integer               ::  irec
    integer               ::  time6(6)
    type(TDate)           ::  trec
    character(len=32)     ::  unit
    real                  ::  unitconv
    
    ! --- begin -------------------------------

    ! info ...
    call wrtgol('LE:   setup OH fields for ',t); call goPr
    
    ! read OH field:
    select case ( arch_format )
    
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( 'netcdf' )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! daily output files, valid for hours [0,24)
        tday = t
      
        ! set file name:  /archive/LE_base_conc_20050102.nc
        write (fname,'(a,"/",3(a,"_"),i4.4,2i2.2,".nc")') &
            trim(arch_path), trim(arch_model), trim(arch_runid), trim(arch_type), &
            tday%year, tday%month, tday%day
        
        ! check ...
        inquire( file=trim(fname), exist=exist )
        if ( .not. exist ) then
          write (gol,'("OH input file not found:")'); call goErr
          write (gol,'("  ",a)') trim(fname); call goErr
          TRACEBACK; status=1; return
        end if

#ifdef with_netcdf
        ! open netcdf file:
        status = NF90_Open( trim(fname), NF90_NoWrite, ncid )
        IF_NF90_NOTOK_RETURN(status=1)
        
        ! record number: 1=00:00, 2=01:00 ...
        irec = tday%hour+1
         
        ! date variable:
        status =  NF90_Inq_VarID( ncid, 'date', varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! read date:
        status = NF90_Get_Var( ncid, varid, time6, start=(/1,irec/), count=(/6,1/) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! check ...
        trec = NewDate( time6=time6 )
        if ( trec /= t ) then
          write (gol,'("unexpected time in record:")'); call goErr
          call wrtgol( '  model time     : ', t ); call goErr
          call wrtgol( '  record time    : ', trec ); call goErr
          write (gol,'("  record number  : ",i4)') irec; call goErr
          TRACEBACK; status=1; return
        end if
        
        ! info ...
        write (gol,'("WARNING - copy OH surface field to all levels ..")')

        ! OH variable:
        status = NF90_Inq_VarID( ncid, 'oh', varid )
        IF_NF90_NOTOK_RETURN(status=1)
        ! read 3D OH field:
        status = NF90_Get_Var( ncid, varid, OH, &
                        start=(/1,1,1,irec/), count=(/nx,ny,1,1/) )
        IF_NF90_NOTOK_RETURN(status=1)
        ! get unit:
        status = NF90_Get_Att( ncid, varid, 'units', unit )
        IF_NF90_NOTOK_RETURN(status=1)
        ! conversion factor; guess OH is in ppb ...
        select case ( unit )
          case ( 'mole mole-1', 'mole mole**-1' )
            unitconv = 1.0e9   ! mole/mole -> ppb
          case default
            write (gol,'("unsupported OH unit : ",a)') trim(unit); call goErr
            TRACEBACK; status=1; return
        end select
        ! convert:
        OH = OH * unitconv
        
        ! close:
        status = NF90_Close( ncid )
        IF_NF90_NOTOK_RETURN(status=1)
        
#else
        stop 'not compiled with netcdf support'
#endif

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        write (gol,'("unsupported OH file format : ",a)') trim(arch_format); call goErr
        TRACEBACK; status=1; return
        
    end select
    
    ! ok
    status = 0

  end subroutine LE_Chem_OH_Setup
  
  
end module LE_Chem_OH

