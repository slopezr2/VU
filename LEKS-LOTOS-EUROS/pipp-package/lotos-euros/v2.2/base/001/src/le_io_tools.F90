!###############################################################################
!
! This module contains subroutines that are used for Input/Output.
!
! io_read_string      : read string from file, preceded by an identification
! io_read_ini         : open file and read header; option to count the number of lines after the header
! io_read_table1      : read table from file
! io_read_label       : read label from file and stop if label is not correct
! io_read_comments    : read block of comment lines (lines starting with #)
! io_read_count_lines : read and count lines from current file pointer on
! io_compare_strings  : compare two strings
! io_compress         : function that deletes leading spaces or tabs and 
!                       replaces multiple spaces or tabs by one space
! io_err1             : write start of error message (with file name and record number)
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

module LE_IO_Tools

  use GO, only : gol, goPr, goErr

  implicit none


  ! --- in/out -----------------------------

  public
  

  ! --- const ------------------------------

  character(len=*), parameter   ::  mname = 'LE_IO_Tools'
  
  ! debug option 
  ! debopt = 0 -> no output
  ! debopt = 1 -> output on read processes etc. to log file
  integer, parameter  :: debopt = 0


contains


!---------------------------------------------------------------------------------------
subroutine io_read_string(u_in,u_err,id,string,fnam,irec)

! read string from file.
! input record is 
! identification: string
! it is an error if the identification is not present in the line currently read
!

integer, intent(in)           :: u_in  ! unit number input file
integer, intent(in)           :: u_err ! unit number error file
character(len=*), intent(in)  :: id    ! identification string for file
character(len=*), intent(out) :: string! string read from file
character(len=*), intent(in)  :: fnam  ! file name
integer                       :: irec  ! record number

! local
character(len=250) :: line
character(len=250) :: id_read
integer            :: indx

! read line:
read(u_in,'(a)') line
irec = irec + 1

! Get position of ":" and write error message if not present:
indx = index(line,':')
if (indx .eq. 0) then
   call io_err1(fnam,irec,u_err)
   write(u_err,*) ' expecting to read record of type "identification: string" '
   write(u_err,*) ' cannot find ":" in line read'
   write(u_err,*) ' line read: ',trim(line)
   write(u_err,*) ' '
   stop
endif

! Split into identification and file name:
id_read = adjustl(line(1:indx-1))
string  = line(indx+1:)

! Check identification:
if (id .ne. id_read) then
   call io_err1(fnam,irec,u_err)
   write(u_err,*) ' identification not correct '
   write(u_err,*) ' identification expected: ',trim(id)
   write(u_err,*) ' identification read    : ',trim(id_read)
   write(u_err,*) ' '
   stop
endif

end subroutine io_read_string

!---------------------------------------------------------------------------------------
subroutine io_read_ini(u_file,u_out,debopt,fnam,irec,nlines)

!
! initialise reading of a formatted file
! 1. open file
! 2. read header (lines starting with #)
! 3. write header to screen
! 4. if optional argument nlines is present: count the number of lines after the header
! 
! The file pointer is left at the last line of the header 
! (next read-statement is first non-header line).
!
! it is a fatal error if the file contains only header lines or is empty.

! arguments:
integer, intent(in)            :: u_file ! unit number file
integer, intent(in)            :: u_out  ! unit number output file
integer, intent(in)            :: debopt ! debug option (if > 0: echo header to file)
character(len=*),intent(in)    :: fnam   ! file name
integer, intent(out)           :: irec   ! record number of last header line 
!                                         (where file pointer is left)
!integer, intent(out), optional :: nlines ! number of lines after the header
! Salford compiler produces strange run-time error if intent(out) is given.
integer, optional :: nlines ! number of lines after the header

! const:
character(len=*), parameter ::  rname = mname//'/io_read_ini'

! Local
!logical                     :: eof    ! end of file
!logical                     :: eoh    ! end of header
!character(len=256)          :: line   ! line read from file 
!integer                     :: status
logical                     ::  exist

!-------------
! Write info to screen:
if (debopt > 0) then
   write(u_out,*) ' '
   write(u_out,'(a)') ' --------------------------------------------------------------------------------------------'
   write(u_out,'(a)') ' opening file ',trim(fnam)
   write(u_out,'(a)') ' --------------------------------------------------------------------------------------------'
   write(u_out,*) ' '
endif

! check ...
inquire( file=trim(fnam), exist=exist )
if ( .not. exist ) then
  write (gol,'("file not found : ",a)') trim(fnam); call goErr
  TRACEBACK; stop
end if

! Open file
open(unit = u_file, file = fnam, status = 'OLD')
rewind(u_file)
irec = 0

! Read header:
call io_read_comments(u_file,u_out,debopt,irec)

! Count number of lines:
if (present(nlines)) then
   call io_read_count_lines(u_file,nlines)

   ! Reset file pointer to end of header
   rewind(u_file)
   irec = 0
   call io_read_comments(u_file,u_out,-1,irec)
endif

end subroutine io_read_ini


!---------------------------------------------------------------------------------------

!************************************************************************
! io_read_table1
!              reads a table with column header line and row names.
!              The input file should contain the following lines:
!
!                 BEGIN id_string
!                 UNIT unit
!                 # header line 1
!                 # header line 2
!                 # ...
!                 h_code         h_rownam    h_1         h_2              h_ncol
!                 rowcode(1)     rownam(1)   val(1,1)    val(1,2)    ...  val(1,ncol)
!                 rowcode(2)     rownam(2)   val(2,1)    val(2,2)    ...  val(2,ncol)
!                 ...........................................
!                 ...........................................
!                 rowcode(nrow)  rownam(2)   val(nrow,1) val(nrow,2) ...  val(nrow,ncol)
!                 END id_string
! 
! - Each row name must be in single quotes.
! - The column header line may not contain quotes.
!
! - It is a fatal error if id_string read from file does not match with the input argument
!   (both at BEGIN and END of table).
! - It is a fatal error if the unit read is not the same as the input argument.
! - It is a fatal error if the column header line (line with h_'s) read from file does not 
!   match with the input argument (the number of spaces between each column header is irrelevant).
! - It is also a fatal error if the code and row names read from file do not match with the input arguments 
!   and are not in the same order. 
!
! Example:
!                 BEGIN daily distribution factors
!                 UNIT [-]
!                 # emission distribution over days in the week
!                 # date: 2005-11-23
!                 # reference: ...
!                 code   category  mon tue wed thu fri sat  sun
!                 1100   industry  1.1 1.1 1.1 1.1 1.1 0.75 0.75
!                 1100   traffic   1.2 1.1 1.1 1.0 1.0 0.8  0.8
!                 ...........................................
!                 ...........................................
!                 END daily distribution factors
!
!************************************************************************
     subroutine io_read_table1(id_string,unit,nrow,ncol, &
                               colheader,rowcode,rownam, &
                               firstread,fnam,irec,u_in,u_out,u_err,debopt, &
                               mat, status )

      use LE_Emis_Tools, only : Get_ShortSNAP, ShortSNAP_to_Code

!CONSTANTS
    
      character(len=*), parameter ::  rname = mname//'/io_read_table1'

!VARIABLES
! 
! INPUT
!   id_string      : identification string for matrix
!   unit           : unit to be used for values (read from line after BEGIN id_string)
!                    unit is surrounded by brackets []; dimensionless -> [-]
!   nrow           : number of rows (excluding header row)
!   ncol           : number of data-columns to read (excluding
!                    columns with code and name)
!   colheader      : column header line
!   rowcode        : row codes
!   rownam         : row names
!   firstread      : first matrix is read (open file)
!   fnam           : file name
!   u_in           : logical unit number input file
!   u_out          : logical unit number output file
!   u_err          : logical unit number error file
!   debopt         : debug option (screen output); if > 0: echo input to screen
      integer, intent(in)          :: ncol
      integer, intent(in)          :: nrow
      character(len=*), intent(in) :: fnam
      character(len=*), intent(in) :: id_string
      logical, intent(in)          :: firstread
      integer, intent(in)          :: u_in
      integer, intent(in)          :: u_out
      integer, intent(in)          :: u_err
      integer, intent(in)          :: debopt
      character(len=*), intent(in) :: rowcode(nrow)
      character(len=*), intent(in) :: rownam(nrow)
      character(len=*), intent(in) :: unit
      character(len=*), intent(in) :: colheader

! INPUT/OUTPUT
!   irec           : record number
      integer, intent(inout)     :: irec

! OUTPUT
!   mat            : matrix with data for component
      real, intent(out)             :: mat(nrow,ncol)
      integer, intent(out)          :: status

! LOCAL
!    line          : line read from file
!
      character(len=500) :: line
      character(len=50)  :: fmt
      character(len=500) :: unit_rd, rownam_rd
      integer            :: shortsnap
      character(len=8)   :: code_rd
      integer            :: irow
      integer            :: ios   ! IO status
      real, allocatable,dimension(:) :: rdumx

!***********************************************************************
!
      allocate(rdumx(ncol))
!
!---- If first matrix is read, open file and read header:
      IF (firstread) THEN
         CALL io_read_ini(u_in,u_out,debopt,fnam,irec)
      ENDIF
!      print*,'     ',trim(fnam)
!
!---- Read "BEGIN id_string":
      call io_read_label('BEGIN '//id_string,u_in,u_err,fnam,irec)
      if (debopt .gt. 0) write(u_out,'(a)') 'BEGIN '//trim(id_string)
!
!---- Read line with UNIT:
      read(u_in,'(a)',IOSTAT=ios) line
      call io_check(ios,fnam,irec,u_in,u_err)
      line = adjustl(line)
      if (line(1:4) .ne. 'UNIT') then
         call io_err1(fnam,irec,u_err)
         write(u_err,*) ' keyword UNIT not found in table '
         write(u_err,'(a,a)') ' line read: ', trim(line)
         stop
      else

!------- Read unit and check:
         read(line(5:),'(a)') unit_rd
         if (io_compress(unit_rd) .ne. io_compress(unit)) then
            call io_err1(fnam,irec,u_err)
            write(u_err,*) ' incorrect unit found in table '
            write(u_err,*) ' unit read    : ', trim(unit_rd)
            write(u_err,*) ' unit expected: ', trim(unit)
            stop
         else
           if (debopt .gt. 0)  write(u_out,'(a)') 'UNIT ' // trim(unit)
         endif
      endif
!
!---- Read comment lines:
      call io_read_comments(u_in,u_out,debopt,irec)
!
!---- Read column header line:
      call io_read_label(colheader,u_in,u_err,fnam,irec)
      if (debopt .gt. 0)  write(u_out,'(a)') colheader
!
!---- Construct format for writing:
      if (debopt .gt. 0) then
         fmt = '(i6,1x,a50,   (e12.5))'
         write(fmt(12:14),'(i3)') ncol
      endif

      ! loop over rows:
      do irow = 1,nrow

         ! read row with code, name, values from file:
         read (u_in,*,IOSTAT=ios) code_rd, rownam_rd, rdumx(1:ncol)
         call io_check(ios,fnam,irec,u_in,u_err)
         
         ! snap or gnfr?
         if ( index('0123456789',code_rd(1:1)) > 0 ) then
           ! convert from '1' or '01.00' to 1:
           call Get_ShortSNAP( code_rd, shortsnap, status )
           IF_NOTOK_RETURN(status=1)
           ! convert from 1 to '01.00'
           call ShortSNAP_to_Code( shortsnap, code_rd, status )
           IF_NOTOK_RETURN(status=1)
         end if

         ! check row code:
         if ( trim(code_rd) /= trim(rowcode(irow)) ) then
            call io_err1(fnam,irec,u_err)
            write (gol,*) ' inconsistency in code '; call goErr
            write (gol,*) ' code read from file: ', code_rd; call goErr
            write (gol,*) ' code expected      : ', rowcode(irow); call goErr
            TRACEBACK; status=1; return
         endif

         ! check row name:
         if (io_compress(rownam_rd) .ne. io_compress(rownam(irow))) then
            write(*,*) 'io_read_table1---',trim(rownam_rd),'---'
            write(*,*) 'io_read_table1---',trim(rownam(irow)),'---'
            write (gol,*) ' inconsistency in row name '; call goErr
            write (gol,*) ' row name read from file: ', trim(rownam_rd); call goErr
            write (gol,*) ' row name expected      : ', trim(rownam(irow)); call goErr
            TRACEBACK; status=1; return
         endif

         ! copy data:
         mat(irow,:) = rdumx(:)

      enddo
!
!---- Read "END id_string":
      call io_read_label('END '//id_string,u_in,u_err,fnam,irec)
      if (debopt .gt. 0) write(u_out,'(a)') 'END '//trim(id_string)

      deallocate(rdumx)
      
      ! ok
      status = 0

      end subroutine io_read_table1


!---------------------------------------------------------------------------------------

!************************************************************************
! io_read_table_country_cat
!              reads a table with column header line and rows consisting of country/category and real values.
!              The input file should contain the following lines:
!
!                 BEGIN id_string
!                 UNIT unit
!                 # header line 1
!                 # header line 2
!                 # ...
!                 h_country        h_cat_code       h_cat_nam          h_1           h_2                h_ncol
!                 country_code(1)  category_code(1) category_name(1)   mat(1,1,1)    mat(1,1,2)    ...  mat(1,1,ncol)
!                 country_code(1)  category_code(2) category_name(2)   mat(2,1,1)    mat(2,1,2)    ...  mat(2,1,ncol)
!                 ...
!                 country_code(2)  category_code(1) category_name(1)   mat(1,2,1)    mat(1,2,2)    ...  mat(1,2,ncol)
!                 country_code(2)  category_code(2) category_name(2)   mat(2,2,1)    mat(2,2,2)    ...  mat(2,2,ncol)
!                 ...
!                 END id_string
! 
! Output are floating point values in array mat(1:cat,1:ncountry,1:ncol).
!
! - Each row name must be in single quotes.
! - The column header line may not contain quotes.
! - Emission categories may differ for each country.
!
! - It is a fatal error if id_string read from file does not match with the input argument
!   (both at BEGIN and END of table).
! - It is a fatal error if the unit read is not the same as the input argument.
! - It is a fatal error if the column header line (line with h_'s) read from file does not 
!   match with the input argument (the number of spaces between each column header is irrelevant).
! - It is also a fatal error if the country names and emission category code and names read from file 
!   do not match with the input arguments and are not in the same order. 
!
! Example:
!                 BEGIN daily distribution factors
!                 UNIT [-]
!                 # emission distribution over days in the week
!                 # date: 2005-11-23
!                 # reference: ...
!                 country code   category                         mon tue wed thu fri sat  sun
!                 'ALB'   1      industry                         1.1 1.1 1.1 1.1 1.1 0.75 0.75
!                 'ALB'   2      xxx                              1.1 1.1 1.1 1.1 1.1 0.75 0.75
!                 ...........................................
!                 'NLD'   1100   'industrie, voedingsmiddelen'    1.1 1.1 1.1 1.1 1.1 0.75 0.75
!  
!                 ...........................................
!                 END daily distribution factors
!
!************************************************************************
     subroutine io_read_table_country_cat(id_string,unit,ndistr,mcat,ncat,ncountry,ncol,country2distr, &
                               colheader,country_code,cat_code,cat_nam, &
                               firstread,fnam,irec,u_in,u_out,u_err,debopt, &
                               mat)

!CONSTANTS

!VARIABLES
! 
! INPUT
!   id_string      : identification string for matrix
!   unit           : unit to be used for values (read from line after BEGIN id_string)
!                    unit is surrounded by brackets []; dimensionless -> [-]
!   ndistr         : number of emission categories
!   ncountry       : number of countries
!   mcat           : maximal number of emission categories 
!   ncat           : number of emission categories (for each emission distribution)
!   ncol           : number of data-columns to read (excluding columns with country, cat_code and cat_name)
!   country2distr  : index of emission distribution used for each country
!   colheader      : column header line
!   country_code   : country codes
!   cat_code       : emission category codes
!   cat_nam        : emission category names
!   firstread      : first matrix is read (open file)
!   fnam           : file name
!   u_in           : logical unit number input file
!   u_out          : logical unit number output file
!   u_err          : logical unit number error file
!   debopt         : debug option (screen output); if > 0: echo input to screen
      integer, intent(in)          :: ndistr
      integer, intent(in)          :: ncountry
      integer, intent(in)          :: mcat
      integer, intent(in)          :: ncat(ndistr)
      integer, intent(in)          :: ncol
      integer, intent(in)          :: country2distr(ncountry)
      character(len=*), intent(in) :: fnam
      character(len=*), intent(in) :: id_string
      logical, intent(in)          :: firstread
      integer, intent(in)          :: u_in
      integer, intent(in)          :: u_out
      integer, intent(in)          :: u_err
      integer, intent(in)          :: debopt
      character(len=*), intent(in) :: country_code(ncountry)
      character(len=*), intent(in) :: cat_code(mcat,ndistr)
      character(len=*), intent(in) :: cat_nam(mcat,ndistr)
      character(len=*), intent(in) :: unit
      character(len=*), intent(in) :: colheader

! INPUT/OUTPUT
!   irec           : record number
      integer, intent(inout)     :: irec

! OUTPUT
!   mat            : matrix with data for component
      real, intent(out)             :: mat(mcat,ncountry,ncol)

! LOCAL
!    line          : line read from file
!
      character(len=250) :: line
      character(len=50)  :: fmt
      character(len=250) :: unit_rd, cat_nam_rd
      character(len=3)   :: country_code_rd
      character(len=8)   :: cat_code_rd
      integer            :: idistr, icountry, icat
      integer            :: ios   ! IO status
      real, allocatable,dimension(:)    :: rdumx
!      real               :: rdum
      integer,dimension(8) :: idum8

!***********************************************************************
!
      allocate(rdumx(ncol))
!
!---- If first matrix is read, open file and read header:
      IF (firstread) THEN
         CALL io_read_ini(u_in,u_out,debopt,fnam,irec)
      ENDIF
!
!---- Read "BEGIN id_string":
      call io_read_label('BEGIN '//id_string,u_in,u_err,fnam,irec)
      if (debopt .gt. 0) write(u_out,'(a)') 'BEGIN '//trim(id_string)
      ! print '(a,a)','     ',trim(fnam)
      ! print '(a,i3)','1b1, ncol = ',ncol
!
!---- Read line with UNIT:
      read(u_in,'(a)',IOSTAT=ios) line
      call io_check(ios,fnam,irec,u_in,u_err)
      line = adjustl(line)
      if (line(1:4) .ne. 'UNIT') then
         call io_err1(fnam,irec,u_err)
         write(u_err,*) ' keyword UNIT not found in table '
         write(u_err,'(a,a)') ' line read: ', trim(line)
         stop
      else

!------- Read unit and check:
         read(line(5:),'(a)') unit_rd
         call io_compare_strings('unit',unit_rd,unit,u_err,fnam,irec)
         if (debopt .gt. 0)  write(u_out,'(a)') 'UNIT ' // trim(unit)
      endif
!
!---- Read comment lines:
      call io_read_comments(u_in,u_out,debopt,irec)
!
!---- Read column header line:
      call io_read_label(colheader,u_in,u_err,fnam,irec)
      if (debopt .gt. 0)  write(u_out,'(a)') colheader
      
!---- Construct format for writing:
      if (debopt .gt. 0) then
         fmt = '(a5,1x,a8,1x,a50,   (e12.5))'
         write(fmt(18:20),'(i3)') ncol
! write(*,*) ' io_mod2/500: ',fmt
      endif

      ! Loop over countries:
      do icountry = 1,ncountry

         ! Get corresponding emission distribution:
         idistr = country2distr(icountry)

         ! Loop over countries and emission categories:
         do icat = 1,ncat(idistr)
            ! Read row with code, name, values from file:
            if ( ncol == 8 ) then
               read(u_in,*,IOSTAT=ios) country_code_rd,cat_code_rd,cat_nam_rd,idum8(1:ncol)
               mat(icat,icountry,1:ncol) = real(idum8(1:ncol))
               ! print '(a,a3,i3,x,a,8i4)','1b1a: ',country_code_rd,cat_code_rd,trim(cat_nam_rd),idum8(1:ncol)
            else
               read(u_in,*,IOSTAT=ios) country_code_rd,cat_code_rd,cat_nam_rd,rdumx(1:ncol)
               mat(icat,icountry,1:ncol) = rdumx(1:ncol)
               ! print '(a,a3,i3,x,a,f12.4)','1b1a: ',country_code_rd,cat_code_rd,trim(cat_nam_rd),rdumx(1:ncol)
            end if
            call io_check(ios,fnam,irec,u_in,u_err)

            ! Check country code, emission category code and name:
            call io_compare_strings('country code',country_code_rd,country_code(icountry),u_err,fnam,irec)
            call io_compare_strings('emission category code',cat_code_rd,cat_code(icat,idistr),u_err,fnam,irec)
            call io_compare_strings('emission category name',cat_nam_rd,cat_nam(icat,idistr),u_err,fnam,irec)
            if (debopt .gt. 0) write(u_out,fmt) country_code_rd,cat_code_rd,cat_nam_rd,mat(icat,icountry,:)

         enddo
      enddo
!
!---- Read "END id_string":
      call io_read_label('END '//id_string,u_in,u_err,fnam,irec)
      if (debopt .gt. 0) write(u_out,'(a)') 'END '//trim(id_string)
      ! print*,'END '//trim(id_string)

      deallocate(rdumx)

     end subroutine io_read_table_country_cat

      
      ! ***
      
      
     subroutine io_read_table_country_cat_data( id_string, unit, ncat, ncountry, ncol, &
                               colheader, country_code, cat_code, cat_nam, &
                               firstread,fnam, irec, &
                               u_in, u_out, u_err, debopt, mat)

      !CONSTANTS

      !VARIABLES
      ! 
      ! INPUT
      !   id_string      : identification string for matrix
      !   unit           : unit to be used for values (read from line after BEGIN id_string)
      !                    unit is surrounded by brackets []; dimensionless -> [-]
      !   ncountry       : number of countries
      !   ncat           : number of emission categories
      !   ncol           : number of data-columns to read (excluding columns with country, cat_code and cat_name)
      !   colheader      : column header line
      !   country_code   : country codes
      !   cat_code       : emission category codes
      !   cat_nam        : emission category names
      !   firstread      : first matrix is read (open file)
      !   fnam           : file name
      !   u_in           : logical unit number input file
      !   u_out          : logical unit number output file
      !   u_err          : logical unit number error file
      !   debopt         : debug option (screen output); if > 0: echo input to screen
      integer, intent(in)          :: ncountry
      integer, intent(in)          :: ncat
      integer, intent(in)          :: ncol
      character(len=*), intent(in) :: fnam
      character(len=*), intent(in) :: id_string
      logical, intent(in)          :: firstread
      integer, intent(in)          :: u_in
      integer, intent(in)          :: u_out
      integer, intent(in)          :: u_err
      integer, intent(in)          :: debopt
      character(len=*), intent(in) :: country_code(ncountry)
      character(len=*), intent(in) :: cat_code(ncat)
      character(len=*), intent(in) :: cat_nam(ncat)
      character(len=*), intent(in) :: unit
      character(len=*), intent(in) :: colheader

      ! INPUT/OUTPUT
      !   irec           : record number
      integer, intent(inout)     :: irec

      ! OUTPUT
      !   mat            : matrix with data for component
      real, intent(out)             :: mat(ncat,ncountry,ncol)

      ! LOCAL
      !    line          : line read from file
      !
      character(len=250) :: line
      character(len=50)  :: fmt
      character(len=250) :: unit_rd, cat_nam_rd
      character(len=3)   :: country_code_rd
      character(len=8)   :: cat_code_rd
      integer            :: icountry, icat
      integer            :: ios   ! IO status

      !***********************************************************************
      !
      !---- If first matrix is read, open file and read header:
      IF (firstread) THEN
         CALL io_read_ini(u_in,u_out,debopt,fnam,irec)
      ENDIF

      !---- Read "BEGIN id_string":
      call io_read_label('BEGIN '//id_string,u_in,u_err,fnam,irec)
      if (debopt .gt. 0) write(u_out,'(a)') 'BEGIN '//trim(id_string)

      !---- Read line with UNIT:
      read(u_in,'(a)',IOSTAT=ios) line
      call io_check(ios,fnam,irec,u_in,u_err)
      line = adjustl(line)
      if (line(1:4) .ne. 'UNIT') then
         call io_err1(fnam,irec,u_err)
         write(u_err,*) ' keyword UNIT not found in table '
         write(u_err,'(a,a)') ' line read: ', trim(line)
         stop
      else

         !------- Read unit and check:
         read(line(5:),'(a)') unit_rd
         call io_compare_strings('unit',unit_rd,unit,u_err,fnam,irec)
         if (debopt .gt. 0)  write(u_out,'(a)') 'UNIT ' // trim(unit)
      endif

      !---- Read comment lines:
      call io_read_comments(u_in,u_out,debopt,irec)

      !---- Read column header line:
      call io_read_label(colheader,u_in,u_err,fnam,irec)
      if (debopt .gt. 0)  write(u_out,'(a)') colheader
      
      !---- Construct format for writing:
      if (debopt .gt. 0) then
         fmt = '(a5,1x,a8,1x,a50,   (e12.5))'
         write(fmt(18:20),'(i3)') ncol
      endif

      ! Loop over countries:
      do icountry = 1,ncountry

         ! Loop over countries and emission categories:
         do icat = 1,ncat
            ! Read row with code, name, values from file:
            read (u_in,*,IOSTAT=ios) country_code_rd, cat_code_rd, cat_nam_rd, mat(icat,icountry,:)
            call io_check(ios,fnam,irec,u_in,u_err)

            ! Check country code, emission category code and name:
            call io_compare_strings('country code',country_code_rd,country_code(icountry),u_err,fnam,irec)
            call io_compare_strings('emission category code',cat_code_rd,cat_code(icat),u_err,fnam,irec)
            call io_compare_strings('emission category name',cat_nam_rd,cat_nam(icat),u_err,fnam,irec)
            if (debopt .gt. 0) write(u_out,fmt) country_code_rd,cat_code_rd,cat_nam_rd,mat(icat,icountry,:)

         enddo
      enddo

      !---- Read "END id_string":
      call io_read_label('END '//id_string,u_in,u_err,fnam,irec)
      if (debopt .gt. 0) write(u_out,'(a)') 'END '//trim(id_string)

    end subroutine io_read_table_country_cat_data


!---------------------------------------------------------------------------------------


subroutine io_read_label(label,u_in,u_err,fnam,irec)

!
! read label from file and stop if label is not the same as input argument
!
character(len=*), intent(in)     :: label   ! label to be read from file
integer, intent(in)              :: u_in    ! logical unit number input file
integer, intent(in)              :: u_err   ! logical unit number error file
character(len=*), intent(in)     :: fnam    ! file name
integer, intent(inout), optional :: irec    ! record number

! Local
character(len=500) :: line            ! line read from file
integer            :: ios             ! IO-status

! Read line:
read(u_in,'(a)',IOSTAT=ios) line
call io_check(ios,fnam,irec,u_in,u_err)

! Check 
!write(*,*) 'io_read_label---',trim(line),'---',len(line),len_trim(line)
!write(*,*) 'io_read_label---',trim(label),'---',len(label),len_trim(label)
if (io_compress(line) .ne. io_compress(label)) then
   call io_err1(fnam,irec,u_err)
   write(u_err,'(a,a)') ' label read     : ', trim(line)
   write(u_err,'(a,a)') ' label expected : ', trim(label)
   stop
endif

end subroutine io_read_label

!---------------------------------------------------------------------------------------
subroutine io_read_comments(u_file,u_out,debopt,irec)

!
! read block of comment lines (lines starting with #)
! the file pointer is left at the first line after this block that is not a comment line
!

! arguments:
integer, intent(in)           :: u_file ! unit number file
integer, intent(in)           :: u_out  ! unit number output file
integer, intent(in)           :: debopt ! debug option (if > 0: echo header to file)
integer, intent(inout)        :: irec   ! record number

! Local
logical                     :: eof    ! end of file
logical                     :: eoh    ! end of header
character(len=250)          :: line   ! line read from file 
character(len=250)          :: line2  ! help line
integer                     :: ios    ! IO-status

!-------------
! Read/write header:
eof = .false.
eoh = .false.
do while (.not. eof .and. .not. eoh)

   ! Read line and check for end of file:
   read(u_file,'(a)',iostat = ios) line
   eof = (ios < 0)
   if (.not. eof) then

      irec = irec + 1

      ! Check for end of header:
      ! (non-empty line with another character as # as first character)
      line2 = adjustl(line)
      eoh = (len_trim(line) > 0 .and. line2(1:1) .ne. '#')

      ! Echo line:
      if (.not. eoh .and. debopt > 0) write(u_out,'(a)') trim(line)
   endif
enddo

! backspace one record to set the file pointer to the last line of the header:
backspace(u_file)
irec = irec - 1

end subroutine io_read_comments

!---------------------------------------------------------------------------------------
subroutine io_read_count_lines(u_file,nlines)

!
! read and count lines until end of file (from the current file record)
! Note: empty lines and comment lines (starting with a #) are not counted 
!
! Arguments:
integer, intent(in)  :: u_file ! unit number file
integer, intent(out) :: nlines ! number of lines read

! Local
logical :: eof ! end of file
integer :: ios ! IO-status
character(len=250) :: line

nlines = 0

! Read  and count lines:
eof = .false.
do while (.not. eof) 
   read(u_file,'(a)',iostat = ios) line
   eof = (ios < 0)
   if (.not. eof) then
      line = adjustl(line)
      if (len_trim(line) > 0 .and. line(1:1) .ne. '#') nlines = nlines + 1
   endif
enddo

end subroutine io_read_count_lines
!---------------------------------------------------------------------------------------
subroutine io_compare_strings(str_nam,str_rd,str_exp,u_err,fnam,irec)

!
! compare two strings, one read from file (fnam,irec), the other expected by the program
! Fatal error if strings are not equal.
!

! Argument:
character(len=*), intent(in) :: str_nam ! name of strings to compare (for error message)
character(len=*), intent(in) :: str_rd  ! string read from file
character(len=*), intent(in) :: str_exp ! string expected by program
integer, intent(in)          :: u_err   ! unit number error file
character(len=*), intent(in) :: fnam    ! file name (used for error message)
integer, intent(in)          :: irec    ! record number (used for error message)

! Local

if (io_compress(str_rd) .ne. io_compress(str_exp)) then
   call io_err1(fnam,irec,u_err)
   write(u_err,*) ' fatal error while reading ',trim(str_nam)
   write(u_err,*) trim(str_nam),' read from file: ',trim(str_rd)
   write(u_err,*) trim(str_nam),' expected      : ',trim(str_exp)
   stop
endif

end subroutine io_compare_strings

!---------------------------------------------------------------------------------------
subroutine io_compare_ints(int_nam,int_rd,int_exp,u_err,fnam,irec)

!
! compare two integers, one read from file (fnam,irec), the other expected by the program
! Fatal error if integers are not equal.
!

! Argument:
character(len=*), intent(in) :: int_nam ! name of integers to compare (for error message)
integer, intent(in)          :: int_rd  ! integer read from file
integer, intent(in)          :: int_exp ! integer expected by program
integer, intent(in)          :: u_err   ! unit number error file
character(len=*), intent(in) :: fnam    ! file name (used for error message)
integer, intent(in)          :: irec    ! record number (used for error message)

! Local

if (int_rd .ne. int_exp) then
   call io_err1(fnam,irec,u_err)
   write(u_err,*) ' fatal error while reading ',trim(int_nam)
   write(u_err,*) trim(int_nam),' read from file: ',int_rd
   write(u_err,*) trim(int_nam),' expected      : ',int_exp
   stop
endif

end subroutine io_compare_ints

!---------------------------------------------------------------------------------------
function io_compress(line)
!
! compress line (delete leading spaces and replace multiple spaces by one space)
!

! Argument:
character(len=*), intent(in) :: line
character(len = len(line))   :: io_compress

! Local
character(len=len(line))   :: line_out         ! output line (with spaces compressed)
character(len=1)           :: c                ! current character from line
integer, parameter         :: iachar_tab = 9   ! ASCII number of TAB
                                               ! xxx is this system dependent ?
logical                    :: prev_space       ! previous character was white space 
integer                    :: ii

!
! current character is white space and previous character was white space --> skip
! current character is white space and previous character was alfa-numeric --> output one space
! current character is alfa-numeric --> output one character
!

! initialisation:
line_out    = ''
ii          = 0
prev_space  = .false.

! Loop 
do while (ii .lt. len_trim(line))

   ! Check whether next character is white space
   ii = ii + 1
   c  = line(ii:ii)
   if (c .eq. ' ' .or. c .eq. char(iachar_tab)) then

      ! add space to line_out, only if this is first space:
      if (.not. prev_space) line_out  = trim(line_out) // ' '

      ! set logical for next cycle of loop
      prev_space = .true.
   else

      ! add character to line_out:
      line_out    = trim(line_out) // c

      ! set logical for next cycle of loop
      prev_space = .false.
   endif
   ! write(*,*) '478: ','---',trim(line_out)
enddo

! output:
io_compress = adjustl(line_out)

end function io_compress

!---------------------------------------------------------------------------------------
subroutine io_err1(fnam,irec,u_err)
!
! issue start of error message to error output
!

! arguments:
character(len=*), intent(in)     :: fnam    ! file name
integer,          intent(in)     :: irec    ! record number
integer,          intent(in)     :: u_err   ! unit number error file

write(u_err,*) ' '
write(u_err,*) ' !----------------------------------------------------------'
write(u_err,*) ' ERROR while reading file ', trim(fnam)
write(u_err,*) ' record number: ',irec
write(u_err,*) ' '

end subroutine io_err1

!---------------------------------------------------------------------------------------
subroutine io_check(ios,fnam,irec,u_in,u_err,eof)

!
! check IO read statement
!
! - increase record number with 1
! - check IO status (error or end of file)
! - issue error message if necessary
! - error while reading --> stop
! - end of file --> if eof is present    : return eof
!                   if eof is not present: stop
!

! arguments:
integer,           intent(in)     :: ios     ! IO status
character(len=*),  intent(in)     :: fnam    ! file name
integer,           intent(inout)  :: irec    ! record number
integer,           intent(in)     :: u_in    ! unit number input file
integer,           intent(in)     :: u_err   ! unit number error file
!logical,intent(out), optional     :: eof     ! end of file has been reached
! Salford compiler produces strange run-time error if intent(out) is given.
logical, optional                 :: eof     ! end of file has been reached

! Local:
character(len=250) :: line

!write(*,*) 'io_mod 780: ', ios,fnam,irec,u_in,u_err

! increase record number
irec = irec + 1

! Check for end of file:
if (present(eof)) then
   eof = (ios < 0)
else

   ! error if end of file:
   if (ios .lt. 0) then
      call io_err1(fnam,irec,u_err)
      write(u_err,*) ' end of file reached '
      write(u_err,*) ' IO-status :',ios
      stop
   endif
endif

! Check for error:
if (ios .gt. 0) then
   call io_err1(fnam,irec,u_err)
   write(u_err,*) ' error occurred when reading record'
   write(u_err,*) ' IO-status :',ios
   write(u_err,*) ' '

   ! echo erroneous record and stop:
   backspace(u_in)
   read(u_in,'(a)') line
   write(u_err,*) ' '
   write(u_err,*) '  erroneous record:'
   write(u_err,*) ' '
   write(u_err,'(a)') trim(line)
   write(u_err,*) ' '
   write(u_err,*) ' '
   stop
endif

end subroutine io_check

subroutine rdlab(labexp,iin,iout)

      integer iin,iout
      character*(*) labexp
      character*79  labin,string

!   Subroutine RDLAB ReaDs a LABel and stops if this label is not
!   as expected.
!   Subroutine RDLAB
!   Author      : Ferd Sauter,Arthur Beusen
!   Date        : 1988-01-27
!   Last update : 2006-02-23 by MSP
!   Call        : CALL RDLAB(LABEXP,IIN,IOUT)
! *
! *   Arguments :
! *  (i)  IIN    ( I*4  ): unit number input file
! *  (i)  IOUT   ( I*4  ): unit number output file
! *  (i)  LABEXP ( L*4  ): expected label.
!=*
! *   Local variables :
! *      LABIN: input label.
! *
! *   No subroutines called.
!-*******************************************************************
! *
 1000 continue
      read(iin,'(a)',end=2001,err=2002) string
      if (string .eq. ' ') goto 1000
      call rddspa(string,labin)

      if (labin .ne. labexp) then
         call rderr1(iin,iout)
         write(*   ,2100) labexp,labin
         write(iout,2100) labexp,labin
         stop
      endif
      goto 2003 !everything is OK

! *-- Error handling:
2001  call rderr(1,iin,iout)
      goto 2003
2002  call rderr(2,iin,iout)
      goto 2003
2003  continue
! *-- Error message:
 2100 format(' No label  : ',a,/,&
             ' Label read: ',a,/,&
             /,/,1x,77('*'),/)

      return
end subroutine rdlab

SUBROUTINE RDDSPA(TXTIN,TXTOUT)
      ! removes spaces at beginning of string
      INTEGER       I,IBEG
      CHARACTER*(*) TXTIN,TXTOUT
      IBEG   = 1
      TXTOUT = TXTIN
      DO 100 I = 1,LEN(TXTIN)
         IF (TXTIN(I:I) .NE. ' ') THEN
            IBEG = I
            GOTO 110
         ENDIF
  100 CONTINUE
  110 CONTINUE
      TXTOUT = TXTIN(IBEG:)
      RETURN
END SUBROUTINE RDDSPA

subroutine rderr(kerr,iin,iout)

!%*   Subroutine RDERR writes ERRor messages which occur when ReaDing
!%*   a file and stops.

! *-- Declarations:
      integer       kerr,iin,iout
      character*79  line

! *   SUBROUTINE RDERR
! *   Author      : Ferd Sauter
! *   Date        : 1989-08-03
! *   Last update : 1990-01-03
! *   CALL        : CALL RDERR(KERR,IIN,IOUT)
! *
! *
! *   Arguments :
! * (i) KERR  = 1: error occurred when reading record
!             = 2: end of file reached
! * (i) IIN   : unit number of the input file (may not be connected to
! *             standard input, which is normally unit 5)
! * (i) IOUT  : unit number of the output file
! *             if IOUT = 0, output is only written to screen
! *             Otherwise output is to a file and to the screen
!-*******************************************************************
! *-- Write heading of error message (including file name):
      call rderr1(iin,iout)
      if (kerr .eq. 1) then
! *----- End of file:
         write(*   ,2100)
         if (iout .ne. 0) write(iout,2100)
      elseif (kerr .eq. 2) then
! *----- Get erroneous record:
         backspace(iin)
         read(iin,'(a)') line
         write(*   ,2300) line
         if (iout .ne. 0) write(iout,2300) line
      endif
! *-- Write finishing line and stop:
      write(*   ,2400)
      if (iout .ne. 0) write(iout,2400)
      stop

! *-- Error messages:
 2100 format(' End of file reached.',/,/)
 2300 format(' Erroneous record: ',/,1x,a,/,/,&
        ' (Note: it is possible that the error is due to missing data',&
        /,' on the previous record)',/)
 2400 format(1x,77('*'),/)

end subroutine rderr

subroutine rderr1(iin,iout)

! *-- Declarations:
      integer       iin,iout
      character*79  fname

!%*   Subroutine RDERR1 writes the first lines of an ERRor message
!%*   (including the file name).
! *   Author      : Ferd Sauter
! *   Date        : 1989-08-03
! *   Last update : 1990-01-03
! *   CALL        : CALL RDERR1(IIN,IOUT)

! *   Arguments :
! * (i) IIN   : unit number of the input file (may not be connected to
! *             standard input, which is normally unit 5)
! * (i) IOUT  : unit number of the output file
! *             if IOUT = 0, output is only written to screen
! *             Otherwise output is to a file and to the screen
!-*******************************************************************

! *-- Get filename:
      inquire(unit = iin, err = 1000, name = fname)
      write(*   ,2100)
! FS1 disabled, gives problems when iin is standard input
!     write(*   ,*   ) fname
!     if (iout .ne. 0) write(iout,2100) fname
! FS2
      goto 2000

 1000 continue
      write(*   ,2200)
      if (iout .ne. 0) write(iout,2200)

 2000 continue

! *-- Error message:
 2100 format(/,1x,5('*'),' Error ',65('*'),/,/,&
       ' RDERR1 ERROR when reading file ')
 2200 format(/,1x,5('*'),' Error ',65('*'),/,/,&
       ' RDERR1 ERROR when reading file with unknown file name.',/,/)

      return
end subroutine rderr1


!---------------------------------------------------------------------------------------
end module LE_IO_Tools
