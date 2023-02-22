!######################################################################
!
! C3PO - CF Convention Compliance Python Objects
!
! Extract region from NetCDF file with grid data.
! For most grids this could be done with using CDO tools:
!   cdo merge -sellonlatbox west,east,south,north input.nc output.nc
! However, CDO cannot handle reduced Gaussian grid.
!
! This program is based on the Python routine:
!    c3po.file.ExtractDomain
! Since this is rather slow, it was re-implemented in Fortran.
!
!
!######################################################################
!
#define TRACEBACK write (*,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; stop 1; end if
!
!######################################################################

program C3PO_Select

  use GO, only : gol, goPr, goErr
  
  use C3PO_Select_Main, only : ParseArguments
  use C3PO_Select_Main, only : ExtractDomain

  implicit none
  
  
  ! --- const --------------------------------------------------
  
  character(len=*), parameter  ::  rname = 'C3PO_Select'
  
  
  ! --- local ---------------------------------------------------
  
  integer                         ::  status
  
  real                            ::  domain(4) ! (w,e,s,n) [deg]
  character(len=1024)             ::  input_filename
  character(len=1024)             ::  output_filename
  
  
  ! --- begin ---------------------------------------------------
  
  ! get arguments:
  call ParseArguments( domain, input_filename, output_filename, status )
  ! usage displayed ?
  if ( status == -1 ) stop 0
  ! check:
  IF_NOT_OK_RETURN(status=1)
      
  ! info ..
  write (gol,'("")'); call goPr
  write (gol,'("C3PO - Select domain from gridded file.")'); call goPr
  write (gol,'("")'); call goPr
  write (gol,'("domain      : ",4f8.2)') domain; call goPr
  write (gol,'("input  file : ",a)') trim(input_filename); call goPr
  write (gol,'("output file : ",a)') trim(output_filename); call goPr
  write (gol,'("")'); call goPr

  ! extract:
  call ExtractDomain( domain, input_filename, output_filename, status )
  IF_NOT_OK_RETURN(status=1)
  
  ! info ...
  write (gol,'("")'); call goPr
  write (gol,'("End")'); call goPr
  write (gol,'("")'); call goPr

end program C3PO_Select

