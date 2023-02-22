!###############################################################################
!
! LE_Input
!
! Read concentrations from conc-3d output files.
! This could be used to run the model as a postprocessing tool
! to create extra output.
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN (action) if (status> 0) then; TRACEBACK; action; return; end if
!
!#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################
  
module LE_Input

  use GO, only : gol, goPr, goErr
  use Indices, only : nspec_all
  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_LE_Input
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'LE_Input'

  ! --- types --------------------------------
  
  type T_LE_Input
    ! filename temmplates for 'conc' files,
    ! time values to be inserted:
    character(len=1024)                 ::  filename_template
    !
  contains
    procedure :: Init            => LE_Input_Init
    procedure :: Done            => LE_Input_Done
    procedure :: PutIn           => LE_Input_PutIn
  end type T_LE_Input


  ! --- var ----------------------------------
  
  logical       ::  with_input_conc
  logical       ::  put_zero
  integer      ::  n_non_zero
  character(len=1024)            ::  list
  character(len=16)            ::  list_non_zero(nspec_all)
  
contains


  ! ====================================================================
  ! ===
  ! === Data
  ! ===
  ! ====================================================================


  subroutine LE_Input_Init( self, rcF, rcbase, status )
  
    use GO     , only : TrcFile
    use GO              , only : goSplitString
    ! --- in/out ---------------------------------
    
    class(T_LE_Input), intent(out)              ::  self
    type(TrcFile), intent(in)                   ::  rcF
    character(len=*), intent(in)                ::  rcbase
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LE_Input_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! flag:
    call rcF%get( rcbase//'.conc.apply', with_input_conc, status )
    IF_NOT_OK_RETURN(status=1)
    
        ! enabled?
    if ( with_input_conc ) then
      ! flag put variables in zero:
       call rcF%get( rcbase//'.conc.put_zero', put_zero, status )
       IF_NOT_OK_RETURN(status=1)
       if (put_zero) then
        ! read non_zero_variables:
        call rcF%Get( rcbase//'.conc.non_zero_variables', list, status )
        IF_NOT_OK_RETURN(status=1)
        call goSplitString( list, n_non_zero, list_non_zero, status )
        IF_NOT_OK_RETURN(status=1)
       end if

      ! filename template:
      call rcF%get( rcbase//'.conc.conc-3d.filenames', self%filename_template, status )
      IF_NOT_OK_RETURN(status=1)
      ! info ..
      write (gol,'(a,": read concentrations from input files:")') rname; call goPr
      write (gol,'(a,":   ",a)') rname, trim(self%filename_template); call goPr
      
    else
      ! info ..
      write (gol,'(a,": no concentrations read input files ...")') rname; call goPr
    end if
    
    ! ok
    status = 0
    
  end subroutine LE_Input_Init


  ! ***


  subroutine LE_Input_Done( self, status )

    ! --- in/out ---------------------------------
    
    class(T_LE_Input), intent(inout)              ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Input_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! ok
    status = 0
    
  end subroutine LE_Input_Done


  ! ***


  subroutine LE_Input_PutIn( self, t, c, status )
  
    use GO              , only : goc
    use GO              , only : TDate, wrtgol
    use GO              , only : goFormat
    use GO              , only : goMatchValue
    use GO              , only : goSplitString
    use C3PO            , only : T_File_Nc
    use Dims            , only : nx, ny, nz
    use LE_Grid         , only : dom
    use Indices         , only : nspec, specname, specunit

    ! --- in/out ---------------------------------
    
    class(T_LE_Input), intent(in)           ::  self
    type(TDate), intent(in)                 ::  t
    real, intent(inout)                     ::  c(nx,ny,nz,nspec)
    integer, intent(out)                    ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/LE_Input_PutIn'
    
    ! --- local ----------------------------------
    
    character(len=1024)       ::  filename
    type(T_File_Nc)           ::  f
    integer                   ::  nlon, nlat, nlev, nt
    integer                   ::  nvar, ivar
    character(len=64)         ::  varname
    integer                   ::  nwarning
    integer                   ::  ispec,isnonzero
    real, allocatable         ::  values(:,:,:)   ! (nlon,nlat,nlev)
    integer                   ::  off(2)
    integer                   ::  irec
    character(len=64)         ::  units
    character(len=128)        ::  conversion
    
    
    
    ! --- begin ----------------------------------
    
    ! enabled?
    if ( with_input_conc ) then
       
      ! info ...
      call wrtgol( rname//': read concentrations for ', t ); call goPr

      ! source file, evaluate time templates:
      call goFormat( t, self%filename_template, filename, status )
      IF_NOT_OK_RETURN(status=1)

      ! open:
      call f%Open( filename, status )
      IF_NOT_OK_RETURN(status=1)

      ! dims:
      call f%Inquire_Dimension( 'longitude', status, length=nlon )
      IF_NOT_OK_RETURN(status=1)
      call f%Inquire_Dimension( 'latitude', status, length=nlat )
      IF_NOT_OK_RETURN(status=1)
      call f%Inquire_Dimension( 'level', status, length=nlev )
      IF_NOT_OK_RETURN(status=1)
      call f%Inquire_Dimension( 'time', status, length=nt )
      IF_NOT_OK_RETURN(status=1)

      ! check, size should be the global shape ..
      if ( any( (/nlon,nlat,nlev/) /= (/dom%glb_shp(1),dom%glb_shp(2),nz/) ) ) then
        write (gol,'("input concentrations of shape (",i0,2(",",i0),") do not match with global shape (",i0,2(",",i0),")")') &
                       nlon,nlat,nlev, dom%glb_shp(1),dom%glb_shp(2),nz; call goErr
        TRACEBACK; status=1; return
      end if
      
      ! storage:
      allocate( values(nx,ny,nz), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! offset in global grid:
      call dom%Get( status, off=off )
      IF_NOT_OK_RETURN(status=1)
      
      ! index of time record:
      call f%Inq_TimeRecord( 'time', t, irec, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! init counter:
      nwarning = 0
      
      ! number of variables:
      call f%Inquire( status, nvar=nvar )
      IF_NOT_OK_RETURN(status=1)
      
      ! loop over variables:
      do ivar = 1, nvar

        ! name:
        call f%Inq_Var( ivar, status, name=varname )
        IF_NOT_OK_RETURN(status=1)
        ! info ..
        write (gol,'(a,":   variable `",a,"`")') rname, trim(varname); call goPr

        ! switch:
        select case ( trim(varname) )

          ! coordinates:
          case ( 'longitude', 'longitude_bnds', 'longitude_crnr', &
                 'latitude', 'latitude_bnds', 'latitude_crnr', &
                 'level', 'time', &
                 'altitude', 'T', 'Q', 'hp', 'halt' )
            ! skip

          ! tracer or other
          case default
            ! try if this is a tracer:
            call goMatchValue( trim(varname), specname, ispec, status, quiet=.true. )
            !~ no match?
            if ( status < 0 ) then
              ! no tracer ..
              write (gol,'(a,":     WARNING - unsupported variable name `",a,"` in input file:")') &
                                 rname, trim(varname); call goPr
              nwarning = nwarning + 1
            !~ match:
            else if ( status == 0 ) then
              ! could be accumulated ...
              if ( ispec <= nspec ) then
                ! info ...
                write (gol,'(a,":     read ...")') rname; call goPr
                ! read slab:
                call f%Get_Var( 'var_name='//trim(varname), values, units, status, &
                                  start=(/off(1)+1,off(2)+1,1,irec/), count=(/nx,ny,nz,1/) )
                IF_NOT_OK_RETURN(status=1)
                ! check:
                if ( trim(units) /= trim(specunit(ispec)) ) then
                  ! converstion:
                  conversion = trim(units)//' -> '//trim(specunit(ispec))
                  ! convert:
                  select case ( trim(conversion) )
                    !~ volume mixing ratios, mass concentration:
                    case ( 'mole mole-1 -> ppb', 'kg m-3 -> ug/m3' )
                      ! apply factor:
                      values = values * 1e9
                    !~
                    case default
                      write (gol,'("unsupported conversion `",a,"`")') trim(conversion); call goErr
                      TRACEBACK; status=1; return
                  end select ! conversions
                end if ! units differ
                ! store:
                if (put_zero) then
                    call goMatchValue( trim(varname), list_non_zero, isnonzero, status, quiet=.true. )
                    if (status<0) then 
                        c(1:nx,1:ny,1:nz,ispec) = 0.0
                          write (gol,'(a,": Put variable zero `",a,"`")') &
                                 rname, trim(varname); call goPr
                    else if ( status == 0 ) then
                        c(1:nx,1:ny,1:nz,ispec) = values
                    end if       
                else
                     c(1:nx,1:ny,1:nz,ispec) = values
                end if    
                
              else
                ! info ...
                write (gol,'(a,":     skip; accumulated tracer ...")') rname; call goPr
              end if   ! ispec <= nspec

            !~ error ...
            else
              TRACEBACK; status=1; return
            end if

        end select ! varname

      end do  ! ivar

      ! close:
      call f%Close( status )
      IF_NOT_OK_RETURN(status=1)
      
      ! warnings?
      if ( nwarning > 0 ) then
        write (gol,'("unsupported variables in file: ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! clear:
      deallocate( values, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if ! with input

    ! ok
    status = 0
    
  end subroutine LE_Input_PutIn


end module LE_Input


