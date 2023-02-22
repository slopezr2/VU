!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "lekf.inc"
!
!###############################################################################

module LEKF_meas_hour

use LEKF_meas

implicit none

integer   ::  out_unit
integer   ::  out_record


contains

subroutine initmeas(yy,mm,dd,hh)

! opens the measurement file (if present)
! end sets some parameters

use dims, only : nx, ny, runF, locF
use utils, only : indomain
use indices, only : i_o3
use LEKF_state, only :  m_s, maxmeas
use constants, only : cum_days
use units, only : u_tmp
use sysdep

!input
integer :: yy, mm, dd, hh

!local
integer :: i, icnt, h_help, dummy

! open the file that specifies the measurements
!open(u_tmp,file='/user4/lotos/measurements/meas/grnd/emep/ok.txt')
open(u_tmp,file='/linuxMA5b/lotos/input/measurements/ozon/2003/o3_emep03.txt')

! read the number of measurements
read (u_tmp,*) nmeas

print *, '>>>KF<<< number of stations in the file ', nmeas
print *, '>>>KF<<< the following station are within the domain:'

! counter for the actual number of measurements
icnt = 0

! read all locations
do i=1,nmeas
  icnt = icnt + 1
  if (icnt > maxmeas) stop 'increase the maximum number of measurements'
  read(u_tmp,*) meas(icnt)%code, &
            meas(icnt)%lon, meas(icnt)%lat, dummy, meas(icnt)%irec
!  print *, meas(icnt)%code, meas(icnt)%lon, meas(icnt)%lat, meas(icnt)%irec
           
  meas(i)%name='unknown'

  ! check if meas is in domain
  call indomain(meas(icnt)%lon, meas(icnt)%lat, runF%dlon, runF%dlat, runF%southb, runF%westb, &
                meas(icnt)%ix,  meas(icnt)%iy, meas(icnt)%indomain)

  if (meas(icnt)%indomain) then
     print *, '>>>KF<<<',meas(icnt)%code, meas(icnt)%name, &
                         meas(icnt)%lat, meas(icnt)%lon
     ! the index of the measurement of the first hour is the icnt-th element
     ! from the start in the measurement part of the state vector
     meas(icnt)%index = m_s-1 + icnt
     meas(icnt)%rho   = 100.0
     select case (meas(icnt)%code(1:2))
       case ('NL','nl')
         meas(icnt)%href  = 3.6
       case ('GB','gb')
         meas(icnt)%href  = 3.0
       case default
         meas(icnt)%href  = 2.0
     end select
  else
     ! reset meas number
     icnt = icnt - 1
  endif
enddo

! set number of measurements in domain
nmeas = icnt
print *, '>>>KF<<<    number of meas within the domain ', nmeas

close(u_tmp)

! open the meas file
!open(meas_unit,file='/user4/lotos/measurements/meas/grnd/emep/emep97.dat', &
!     form='unformatted', &
!     recl=emepmeas*reclen_fac_bin_io ,access='direct',status='old')
open(meas_unit,file='/linuxMA5b/lotos/input/measurements/ozon/2003/o3_emep03.dat', &
     form='unformatted', recl=emepmeas*reclen_fac_bin_io ,access='direct',status='old') 

out_unit = u_meas_min+2
open(out_unit,file=trim(locF%outpath)//'o3_stat.dat', &
     form='unformatted', recl=1*reclen_fac_bin_io ,access='direct',status='new') 

! determine the record where to start
! assumption: binary file starts with the first hour of the year in GMT
! the record is set to 
meas_record = (cum_days(mm) + dd - 1)*24 + hh
out_record = 0
end subroutine

subroutine measupdate(nmodes,mm,dd,hh)

! implementation of measurement update
! routine expects a measurement list +
! a list specifying the corresponding element nr.
! of the statevector.

use LEKF_state
use LEKF_reduce
use units

! input: the actual number of modes
integer :: nmodes

!local copy of the meas + its sd
real :: c, sd, z, z1(nmeas), z2(nmeas), temp_xb(emepmeas),temp_x(emepmeas), temp_meas(emepmeas)
integer :: i, ind, ix, iy, ihist
integer :: yy, mm, dd, hh
real    :: rho
integer :: ind2


if (nmeas > 0) then
print *, '<measurement updating>'

do i=1,nmeas
  if (.NOT.meas(i)%indomain) stop 'measurement not in domain ???!!'
  if (meas(i)%indomain .AND. meas(i)%c >= 0) then
    ind = meas(i)%index
    z1(i) = x(ind)-meas(i)%c
  endif
enddo

do i=1,nmeas

! ! reduce rank if necessary
! if (nmodes == maxmodesnoise) call reduce_rank(nmodes)

  if (meas(i)%c > 0.0 .AND. meas(i)%indomain) then
      ind = meas(i)%index
      c   = meas(i)%c
      sd  = meas(i)%sd
      ix  = meas(i)%ix
      iy  = meas(i)%iy
      rho = meas(i)%rho
      if (ind > nstate) then
        print *, 'index too large! nstate = ',nstate,' index = ',meas(i)%index
        stop
      else if (nmodes > maxmodesnoise) then
        print *, 'nmodes too large! maxmodes = ',maxmodesnoise,'nmodes = ', nmodes
        stop
      endif
      call meas_update(c,sd,rho,ix,iy,ind,nmodes,4)
   endif

enddo

! compute the final residues and write to file
do i=1,nmeas
  if (meas(i)%c > 0) then
     ind = meas(i)%index
     z2(i) = x(ind)-meas(i)%c
     write(u_meas_min+3,fmt='(3i3,a7,6f10.3)') mm,dd,hh, meas(i)%code, meas(i)%c, xb(ind), x(ind), xb(ind)-meas(i)%c, z1(i), z2(i)
  else
     ind = meas(i)%index
     c=-999.
     write(u_meas_min+3,fmt='(3i3,a7,6f10.3)') mm,dd,hh, meas(i)%code, c ,  xb(ind), x(ind), c,c,c !-999.0 waar noodazakelijk
  endif
enddo
!write binary file like the original emep data file .... than same control file to be used
   temp_meas = -999.
do i=1,nmeas
  if (meas(i)%c > 0) then
     ind = meas(i)%index ! in state vector
     ind2  = meas(i)%irec  ! in original meas file
   temp_xb(ind2) = xb(ind)
   temp_x(ind2) = x(ind)
   temp_meas(ind2) = meas(i)%c
  endif
enddo
do ind2=1,emepmeas
   out_record = out_record + 1
   write (out_unit,rec=out_record) temp_meas(ind2)
   out_record = out_record + 1
   write (out_unit,rec=out_record) temp_xb(ind2)
   out_record = out_record + 1
   write (out_unit,rec=out_record) temp_x(ind2)
enddo

! flush for convenience of viewing intermediate results
call flush(u_meas_min)

print *, '<finished measurement updating>'

! one more time: reduce rank if necessary
!if (nmodes > maxmodes) call reduce_rank(nmodes)

else

print *, '<no measurement updating: no measurements available>'

endif

return
end subroutine

subroutine getmeas

! local
real :: c(emepmeas)
integer  ::  i, ind

meas_record = meas_record + 1
read (meas_unit,rec=meas_record) (c(i), i=1,emepmeas)

! set meas
do i=1,nmeas
   ind       = meas(i)%irec
   meas(i)%c = c(ind)
   ! compute the s.d.
   ! (10% of the concentration with a minimum of 1 and a maximum of 5)
   meas(i)%sd = min(5.0, max(1.0,0.1*meas(i)%c))
!  print *, 'read meas #',i,' value = ', c(ind)
   !convert to ug/m3!
   meas(i)%sd = meas(i)%sd / 2.
   meas(i)%c = meas(i)%c / 2.
enddo

end subroutine

end module
