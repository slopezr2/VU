!###############################################################################
!
! NAME
!   LE_Emis_OPS_PM_Composition - PM composition file for OPS emissions, also used for country/SNAP information
!
! Distribution of PM emissions over different size distributions
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

module LE_Emis_OPS_PM_Composition

  use GO, only : gol, goPr, goErr
  use GO, only : TrcFile


  implicit none


  ! --- in/out -------------------------------

  private

  public  ::  LE_Emis_OPS_PM_Composition_Init
  public  ::  LE_Emis_OPS_PM_Composition_Done

  ! --- const --------------------------------

  character(len=*), parameter ::  mname = 'LE_Emis_OPS_PM_Composition'
  
contains  
  
  ! ========================
  
  subroutine LE_Emis_OPS_PM_Composition_Init( query, max_country,max_cat,country_ISO3,country_codes, ncountry, cat_codes, ncat, pm_comp, status )
  
    use GO, only : GoMatchValue, goVarValue, goGetFu
    
    ! --- in/out ---
    
    character(len=*), intent(in)              ::  query
    integer, intent(in)                       ::  max_country, max_cat
    character(len=3), intent(inout)           ::  country_ISO3(max_country)
    integer, intent(inout)                    ::  country_codes(max_country)
    integer, intent(out)                      ::  ncountry
    
    character(len=*), intent(inout)           ::  cat_codes(max_cat)
    integer, intent(out)                      ::  ncat
    real, intent(inout)                       ::  pm_comp(max_country,max_cat,6)
    integer, intent(out)                      ::  status
    
    ! --- const ---
    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_PM_Composition_Init'
    
    ! --- local ---
    logical             ::  exist
    character(len=512)  ::  filename
    integer             ::  fu
    character(len=256)  ::  line
    integer             ::  iline
    character(len=1)    ::  comment
    
    real                ::  class1, class2, class3
    real                ::  class4, class5, class6
    character(len=3)    ::  ISO3
    integer             ::  country_code
    character(len=8)    ::  cat_code
    integer             ::  cat_nr
    integer             ::  icountry, icat
    integer             ::  ii_country, ii_cat
    logical             ::  found
    
    ! --- begin ---
    
    ! extract filename, by default empty:
    filename = ''
    call goVarValue( trim(query), ';', 'file', '=', filename, status )
    IF_ERROR_RETURN(status=1)
    ! comment character:
    comment = '#'
    call goVarValue( trim(query), ';', 'comment', '=', comment, status )
    IF_ERROR_RETURN(status=1)

    ! file should be present:
    inquire( file=trim(filename), exist=exist  )
    if ( .not. exist ) then
      write (gol, '("file not found : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! new file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)

    ! open file:
    open( fu, file=trim(filename), status='old', form='formatted', iostat=status )
    if (status/=0) then
      write (gol,'("opening file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! line counter
    iline = 0
    
    ! read headerline
    do
      ! read line:
      read (fu,'(a)',iostat=status) line
      if (status/=0) then
        write (gol,'("reading header line from file : ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! empty ? then skip:
      if ( len_trim(line) == 0 ) cycle
      ! comment ? then skip:
      if ( line(1:1) == comment ) cycle
      ! found non-comment line, leave loop:
      exit
    end do
    
    ! init
    ncountry = 0
    ncat = 0
    ! loop over records 
    do 
      ! increase record counter: 
      iline = iline + 1
      ! try to read line:
      read (fu,'(i5,i1,6f7.1,2x,a3,a,i2)',iostat=status) country_code,cat_nr , class1, class2, class3, class4, class5, class6, ISO3, cat_code
      
      ! eof ?
      if (status<0) exit
      ! error ?
      if (status>0) then
        write (gol,'("reading line ",i6," from file : ",a)') iline, trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! fill in country_codes
      found = .false.
      do icountry = 1, max_country
        if ( country_codes(icountry) == country_code ) then
          ii_country = icountry
          found = .true.
          exit
        end if
      end do
      
      ! new country
      if ( .not. found ) then
        ncountry = ncountry + 1
        ii_country = ncountry
        country_ISO3(ncountry) = trim(ISO3)
        country_codes(ncountry) = country_code
        
      end if
      
      found = .false.
      ! fill in category codes
      if ( cat_code == '0' ) then
        cat_code = '10'
      end if

      do icat = 1, max_cat
        if ( cat_code == cat_codes(icat) ) then
          ii_cat = icat
          found = .true.
          exit
        end if
      end do
      
      ! new category
      if ( .not. found ) then
        ncat = ncat + 1
        ii_cat = ncat
        cat_codes(ncat) = cat_code
      end if
      
      ! fill in fractions per class
      pm_comp(ii_country,ii_cat,1) = 1e-2*class1
      pm_comp(ii_country,ii_cat,2) = 1e-2*class2
      pm_comp(ii_country,ii_cat,3) = 1e-2*class3
      pm_comp(ii_country,ii_cat,4) = 1e-2*class4
      pm_comp(ii_country,ii_cat,5) = 1e-2*class5
      pm_comp(ii_country,ii_cat,6) = 1e-2*class6
      
    end do    
      
    ! ok
    status = 0
    
  end subroutine LE_Emis_OPS_PM_Composition_Init
  
  subroutine LE_Emis_OPS_PM_Composition_Done( status )

    ! --- in/out ---
    integer, intent(out)                      ::  status
    
    ! --- const ---
    character(len=*), parameter ::  rname = mname//'/LE_Emis_OPS_PM_Composition_Done'
    
    ! --- local ---
    
    ! --- begin ---
    
    ! ok
    status = 0
  
  end subroutine LE_Emis_OPS_PM_Composition_Done
  
  
end module LE_Emis_OPS_PM_Composition
