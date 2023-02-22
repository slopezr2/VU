C-----------------------------------------------------------------------
C
C   SUBROUTINE: SOLARANGLE
C
C   Description: To calculate the solar zenith angle.  This will give
C                SIN(BETA), not the BETA.
C
C   Call: None
C
C   Require: None
C
C   Input:
C            1) Day of year
C            2) Latitude
C            3) Hour
C
C   Output: CALCBETA (Solar zenith angle)
C
C   Created by Tan 11/15/06  (bases on xxxx's program)
C
C-----------------------------------------------------------------------
      SUBROUTINE SOLARANGLE( DAY, SHOUR, LAT, SINbeta)

      IMPLICIT NONE

C...  Define parameters
      INTEGER :: DAY                  ! DOY or julian day
      REAL    :: SHOUR                ! Solar hour
      REAL    :: LAT                  ! Latitude
      REAL    :: BETA                 ! solar elevation angle
      REAL    :: sindelta, cosdelta, A, B
      REAL    :: SINbeta

      REAL    :: PI, D2RAD, RAD2D

C...  Constants
      PI           = 3.14159
      D2RAD        = PI/180.0
      RAD2D        = 180.0/PI

C...  Calculation
      sindelta = -SIN(0.40907) * COS( 6.28*(DAY+10)/365 )
      cosdelta = (1-sindelta**2.)**0.5

      A = SIN( LAT*D2RAD ) * sindelta
      B = COS( LAT*D2RAD ) * cosdelta

      SINbeta = A + B * COS(2*PI*(SHOUR-12)/24)  ! This will be transfered
                                                 ! to gamma_p function

      BETA = ASIN(sinbeta)*RAD2D    ! This is not used.

      RETURN
      END SUBROUTINE SOLARANGLE
C-----------------------------------------------------------------------

