!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "le.inc"
!
!###############################################################################

module LE_Country_Data

  !
  ! country name table 
  !
  character(len=3)  :: country_name(100)
  !
  data country_name(  1)/'AL '/
  data country_name(  2)/'AT '/
  data country_name(  3)/'BE '/
  data country_name(  4)/'BG '/
  data country_name(  5)/'FCS'/
  data country_name(  6)/'DK '/
  data country_name(  7)/'FI '/
  data country_name(  8)/'FR '/
  data country_name(  9)/'FGD'/
  data country_name( 10)/'FFR'/
  data country_name( 11)/'GR '/
  data country_name( 12)/'HU '/
  data country_name( 13)/'IS '/
  data country_name( 14)/'IE '/
  data country_name( 15)/'IT '/
  data country_name( 16)/'LU '/
  data country_name( 17)/'NL '/
  data country_name( 18)/'NO '/
  data country_name( 19)/'PL '/
  data country_name( 20)/'PT '/
  data country_name( 21)/'RO '/
  data country_name( 22)/'ES '/
  data country_name( 23)/'SE '/
  data country_name( 24)/'CH '/
  data country_name( 25)/'TR '/
  data country_name( 26)/'FSU'/
  data country_name( 27)/'GB '/
  data country_name( 28)/'VUL'/
  data country_name( 29)/'REM'/
  data country_name( 30)/'BAS'/
  data country_name( 31)/'NOS'/
  data country_name( 32)/'ATL'/
  data country_name( 33)/'MED'/
  data country_name( 34)/'BLS'/
  data country_name( 35)/'NAT'/
  data country_name( 36)/'RUO'/
  data country_name( 37)/'RUP'/
  data country_name( 38)/'RUA'/
  data country_name( 39)/'BY '/
  data country_name( 40)/'UA '/
  data country_name( 41)/'MD '/
  data country_name( 42)/'RUR'/
  data country_name( 43)/'EE '/
  data country_name( 44)/'LV '/
  data country_name( 45)/'LT '/
  data country_name( 46)/'CZ '/
  data country_name( 47)/'SK '/
  data country_name( 48)/'SI '/
  data country_name( 49)/'HR '/
  data country_name( 50)/'BA '/
  data country_name( 51)/'YU '/
  data country_name( 52)/'MK '/
  data country_name( 53)/'KZ '/
  data country_name( 54)/'GE '/
  data country_name( 55)/'CY '/
  data country_name( 56)/'AM '/
  data country_name( 57)/'MT '/
  data country_name( 58)/'ASI'/
  data country_name( 59)/'LI '/
  data country_name( 60)/'DE '/
  data country_name( 61)/'RU '/
  data country_name( 62)/'MC '/
  data country_name( 63)/'NOA'/
  data country_name( 64)/'EU '/
  data country_name( 65)/'US '/
  data country_name( 66)/'CA '/
  data country_name( 67)/'BIC'/
  data country_name( 68)/'KG '/
  data country_name( 69)/'AZ '/
  data country_name( 70)/'???'/
  data country_name( 71)/'???'/
  data country_name( 72)/'???'/
  data country_name( 73)/'???'/
  data country_name( 74)/'???'/
  data country_name( 75)/'???'/
  data country_name( 76)/'???'/
  data country_name( 77)/'???'/
  data country_name( 78)/'???'/
  data country_name( 79)/'???'/
  data country_name( 80)/'???'/
  data country_name( 81)/'???'/
  data country_name( 82)/'???'/
  data country_name( 83)/'???'/
  data country_name( 84)/'???'/
  data country_name( 85)/'???'/
  data country_name( 86)/'???'/
  data country_name( 87)/'???'/
  data country_name( 88)/'???'/
  data country_name( 89)/'???'/
  data country_name( 90)/'???'/
  data country_name( 91)/'???'/
  data country_name( 92)/'???'/
  data country_name( 93)/'???'/
  data country_name( 94)/'???'/
  data country_name( 95)/'???'/
  data country_name( 96)/'???'/
  data country_name( 97)/'???'/
  data country_name( 98)/'???'/
  data country_name( 99)/'???'/
  data country_name(100)/'???'/


end module LE_Country_Data
