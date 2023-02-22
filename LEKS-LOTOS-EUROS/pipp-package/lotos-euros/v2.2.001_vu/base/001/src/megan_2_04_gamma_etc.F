C=======================================================================
C     MODULE GAMMA_ETC
C
C     This module contain the function to calculate
C     GAMMA_P, GAMMA_T, GAMMA_L, GAMMA_A for BVOCs.
C
C     CONTAINS: 1)GAMMA_LAI
C               2)GAMMA_P
C               3)GAMMA_TISOP
C               4)GAMMA_TNISP
C               5)GAMMA_A
C               6)GAMMA_S
C
C     Note:
C
C     Requirement:
C
C     CALL: SOLARANGLE
C
C     Created by Tan 11/21/06 for MEGAN v2.0
C
C     History:
C     08/01/07 Guenther A. - Move to MEGANv2.02 with modification to
C                            correct calculation of GAMMA_P
C
C=======================================================================

      MODULE GAMMA_ETC

      IMPLICIT NONE

C...  Program I/O parameters

C...  External parameters

      CONTAINS
C***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Scientific algorithm
!
!             Emission = [EF][GAMMA][RHO]
!           where [EF]    = emission factor (ug/m2h)
!                 [GAMMA] = emission activity factor (non-dimension)
!                 [RHO]   = production and loss within plant canopies
!                           (non-dimensino)
!                 Assumption: [RHO] = 1 (11/27/06) (See PDT_LOT_CP.inc)
!
!             GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
!           where [GAMMA_CE]  = canopy correction factor
!                 [GAMMA_age] = leaf age correction factor
!                 [GAMMA_SM]  = soil moisture correction factor
!                 Assumption: [GAMMA_SM]  = 1 (11/27/06)
!
!             GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
!           where [GAMMA_LAI] = leaf area index factor
!                 [GAMMA_P]   = PPFD emission activity factor
!                 [GAMMA_T]   = temperature response factor
!
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!        Derivation:
!             Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P]
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
!             Emission = [EF][GAMMA_ect]{ (1-LDF) + [LDF][GAMMA_P] }
!           where LDF = light dependent function (non-dimension)
!
!     For ISOPRENE
!                 Assumption: LDF = 1 for isoprene            (11/27/06)
!
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!
!     For NON-ISOPRENE
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
!                        { (1-LDF) + [LDF][GAMMA_P] }
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C=======================================================================
C...  Begin module
C=======================================================================


C-----------------------------------------------------------------------
C.....1) Calculate GAM_L (GAMMA_LAI)
C-----------------------------------------------------------------------
C                            0.49[LAI]
C             GAMMA_LAI = ----------------    (non-dimension)
C                         (1+0.2LAI^2)^0.5
C
C     SUBROUTINE GAMMA_LAI returns the GAMMA_LAI values
C-----------------------------------------------------------------------
      SUBROUTINE GAMMA_LAI( NCOLS, NROWS,
     &                      LAI, GAM_L )

      IMPLICIT NONE

      INTEGER NCOLS, NROWS
      REAL LAI(NCOLS,NROWS), GAM_L(NCOLS,NROWS)

      GAM_L = (0.49*LAI) / ( (1+0.2*(LAI**2))**0.5 )

      RETURN
      END SUBROUTINE GAMMA_LAI
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C.....2) Calculate GAM_P (GAMMA_P)
C-----------------------------------------------------------------------
C             GAMMA_P = 0.0         a<=0, a>=180, sin(a) <= 0.0
C
C             GAMMA_P = sin(a)[ 2.46*(1+0.0005(Pdaily-400))*PHI - 0.9*PHI^2 ]
C                                   0<a<180, sin(a) > 0.0
C           where PHI    = above canopy PPFD transmission (non-dimension)
C                 Pdaily = daily average above canopy PPFD (umol/m2s)
C                 a      = solar angle (degree)
C
C                 Note: AAA = 2.46*BBB*PHI - 0.9*PHI^2
C                       BBB = (1+0.0005(Pdaily-400))
C                       GAMMA_P = sin(a)*AAA
C
C                       Pac
C             PHI = -----------
C                   sin(a)*Ptoa
C           where Pac  = above canopy PPFD (umol/m2s)
C                 Ptoa = PPFD at the top of atmosphere (umol/m2s)
C
C             Pac =  SRAD * 4.766 mmmol/m2-s * 0.5
C
C             Ptoa = 3000 + 99*cos[2*3.14-( DOY-10)/365 )]
C           where DOY = day of year
C
C     SUBROUTINE GAMMA_P returns the GAMMA_P values
C-----------------------------------------------------------------------
      SUBROUTINE GAMMA_P( JDATE, JTIME, LAT, LONG,
     &                    NCOLS, NROWS,
     &                    PPFD, D_PPFD, GAM_P )

      IMPLICIT NONE

      CHARACTER*16  :: FUNCNAME = 'GAMMA_P'

      INTEGER NCOLS, NROWS, JDATE, JTIME

      INTEGER DAY(NCOLS,NROWS)    ! DAY is DOY (JDATE)
      REAL    HOUR(NCOLS,NROWS)   ! HOUR is solar hour
      REAL    LAT(NCOLS,NROWS)    ! Latitude
      REAL    LONG(NCOLS,NROWS)   ! Longitude
      REAL    PPFD(NCOLS,NROWS)   ! Photosynthetic Photon Flus Density
      REAL    D_PPFD(NCOLS,NROWS) ! Daily PPFD
      REAL    GAM_P(NCOLS,NROWS)  ! GAMMA_P

C...  Local parameters
      CHARACTER*256 MESG          ! message buffer

      REAL    AAA, BBB

      REAL    BETA                ! Solar zenith angle
      REAL    SINbeta             ! sin(beta)
      REAL    SZANGLE(NCOLS,NROWS)! Solar zenith angle array
      
      REAL    Ptoa(NCOLS,NROWS), Pac(NCOLS,NROWS), PHI(NCOLS,NROWS)

      INTEGER I, J

      REAL    LONG_ZONE           ! Longitude in [-180,180]

C...  Constants
      REAL        PI
      PARAMETER ( PI = 3.14159 )
      REAL        D2RAD
      PARAMETER ( D2RAD = PI/180.0 )
      REAL        RAD2D
      PARAMETER ( RAD2D = 180.0/PI )

C...  Convert date and time format to local time
      ! DAY is Julian day
      DAY(:,:)  =  MOD(JDATE,1000)
      ! Convert from XXXXXX format to XX.XX (solar hour)
      ! HOUR = 0 -> 23.xx
      DO I = 1, NCOLS
         DO J = 1, NROWS

            !! Solar hour : ASSUMES THAT LONG in [-180,180]
            !HOUR(I,J) = JTIME/10000 + LONG(I,J)/15

            ! convert longitude to [-180,180] if necessary:
            LONG_ZONE = LONG(I,J)
            if ( LONG_ZONE > 180.0 ) LONG_ZONE = LONG_ZONE - 360.0
            ! Solar hour :
            HOUR(I,J) = JTIME/10000 + LONG_ZONE/15.0
            
            IF ( HOUR(I,J) .LT. 0.0 ) THEN
               HOUR(I,J) = HOUR(I,J) + 24.0
               DAY(I,J) = DAY(I,J) - 1
            ELSEIF ( HOUR(I,J) .GE. 24.0 ) THEN
               HOUR(I,J) = HOUR(I,J) - 24.0
               DAY(I,J) = DAY(I,J) + 1
            ENDIF
               
         ENDDO
      ENDDO

C...  Begin estimating gamma_p
      ! getting solar radiation
      Pac = PPFD

C...  Initialize parameters
      SZANGLE = 0.
      Ptoa = 0.
      PHI = 0.

      DO I = 1, NCOLS
         DO J = 1, NROWS
            ! Get solar elevation angle
            CALL SOLARANGLE( DAY(I,J), HOUR(I,J), LAT(I,J), SINbeta )
            BETA = ASIN(SINbeta)*RAD2D            ! Degree
            SZANGLE(I,J) =  BETA                  ! Degree
            !AJS: problem with underflow for too low solar angles;
            !AJS: changed threshold to skip almost zero angles too:
            !IF (SINbeta .LE. 0.0) THEN
            IF (SINbeta .LE. 0.01) THEN
               GAM_P(I,J) = 0.0
            !ELSEIF (SINbeta .GT. 0.0) THEN
            ELSE
               Ptoa(I,J) = 3000.0 + 99.0 * COS( 2*3.14*(DAY(I,J)-10)/365 )

               PHI(I,J) = Pac(I,J)/(SINbeta * Ptoa(I,J))

               BBB = 1 + 0.0005*( D_PPFD(I,J)-400 )
               AAA = ( 2.46 * BBB * PHI(I,J) ) - ( 0.9 * PHI(I,J)**2 )

               GAM_P(I,J) = SINbeta * AAA
            !ELSE
            !   MESG = 'Error: Solar angle is invalid'
            !   CALL M3EXIT(FUNCNAME,JDATE,JTIME,MESG,2)

            ENDIF

!            IF (BETA .LE. 0.0 .OR. BETA .GE. 180.0) THEN
!               GAM_P(I,J) = 0.0
!
!            ELSEIF (BETA .GT. 0.0 .AND. BETA .LT. 180.0) THEN
!               Ptoa(I,J) = 3000.0 + 99.0 * COS( 2*3.14*(DAY(I,J)-10)/365 )
!
!               PHI(I,J) = Pac(I,J)/(SINbeta * Ptoa(I,J))
!
!               BBB = 1 + 0.0005*( D_PPFD(I,J)-400 )
!               AAA = 2.46 * PHI(I,J)**3 * 0.9 * BBB
!
!               GAM_P(I,J) = SINbeta * AAA
!
!            ELSE
!               MESG = 'Error: Solar angle is invalid'
!               CALL M3EXIT(FUNCNAME,JDATE,JTIME,MESG,2)
!
!            ENDIF

            ! Screening the unforced errors
            ! IF solar elevation angle is less than 1 THEN
            ! gamma_p can not be greater than 0.1.
            IF (BETA .LT. 1.0 .AND. GAM_P(I,J) .GT. 0.1) THEN
               GAM_P(I,J) = 0.0
            ENDIF

         ENDDO    ! End loop for NROWS
      ENDDO       ! End loop for NCOLS

      RETURN
      END SUBROUTINE GAMMA_P
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C.....3) Calculate GAM_T (GAMMA_T) for isoprene
C-----------------------------------------------------------------------
C                          Eopt*CT2*exp(CT1*x)
C             GAMMA_T =  ------------------------
C                        [CT2-CT1*(1-exp(CT2*x))]
C           where x      = [ (1/Topt)-(1/Thr) ] / 0.00831
C                 Eopt   = 1.75*exp(0.08(Tdaily-297)
C                 CT1    = 80
C                 CT2    = 200
C                 Thr    = hourly average air temperature (K)
C                 Tdaily = daily average air temperature (K)
C                 Topt   = 313 + 0.6(Tdaily-297)
C
C                 Note: AAA = Eopt*CT2*exp(CT1*x)
C                       BBB = [CT2-CT1*(1-exp(CT2*x))]
C                       GAMMA_T = AAA/BBB
C
C     SUBROUTINE GAMMA_TISOP returns the GAMMA_T value for isoprene
C-----------------------------------------------------------------------
      SUBROUTINE GAMMA_TISOP( NCOLS, NROWS, TEMP, D_TEMP, GAM_T )

      IMPLICIT NONE

      INTEGER NCOLS, NROWS

      REAL    TEMP(NCOLS,NROWS)              ! hourly surface temperature
      REAL    D_TEMP(NCOLS,NROWS)            ! daily surface temperature
      REAL    GAM_T(NCOLS,NROWS)             ! GAMMA_T

C...  Local parameters
      REAL    Eopt(NCOLS,NROWS), Topt(NCOLS,NROWS), X(NCOLS,NROWS)
      REAL    AAA(NCOLS,NROWS), BBB(NCOLS,NROWS)
      REAL, PARAMETER :: CT1 = 80.0
      REAL, PARAMETER :: CT2 = 200.0

      Eopt = 1.75 * exp(0.08*(D_TEMP-297.0))
      Topt = 313.0 + ( 0.6*(D_TEMP-297.0) )
      X = ( (1/Topt)-(1/TEMP) ) / 0.00831

      AAA = Eopt*CT2*exp(CT1*X)
      BBB = (  CT2-CT1*( 1-exp(CT2*X) )  )

      GAM_T = AAA/BBB

      RETURN
      END SUBROUTINE GAMMA_TISOP
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C.....4) Calculate GAM_T (GAMMA_T) for non-isoprene
C-----------------------------------------------------------------------
C
C             GAMMA_T =  exp[BETA*(T-Ts)]
C           where BETA   = temperature dependent parameter
C                 Ts     = standard temperature (normally 303K, 30C)
C
C     SUBROUITINE GAMMA_TNISP returns the GAMMA_T value for non-isoprene
C-----------------------------------------------------------------------
      SUBROUTINE GAMMA_TNISP( NCOLS, NROWS, SPCNAM, TEMP,
     &                        GAM_T                               )

      use GO, only : goMatchValue
      
      IMPLICIT NONE

      INCLUDE 'megan_2_04_TEMPD_PRM.inc'
      !INCLUDE 'megan_2_04_TEMPD_PRM_AVE.inc'
      !INCLUDE 'megan_2_04_TEMPD_PRM_PSD.inc'
      !INCLUDE 'megan_2_04_TEMPD_PRM_MSD.inc'

      !INTEGER       INDEX1
      !EXTERNAL      INDEX1

      integer        ::  status
      CHARACTER*16  SPCNAM
      INTEGER       NCOLS, NROWS
      INTEGER       SPCNUM                             ! Species number
      REAL           TEMP(NCOLS,NROWS), GAM_T(NCOLS,NROWS)
      REAL, PARAMETER :: Ts = 303.0

      !! use search routine from IOAPI library:
      !SPCNUM = INDEX1(SPCNAM,N_TDF_SPC,TDF_SPC)
      ! use search routine from GO library:
      call goMatchValue( SPCNAM, TDF_SPC, SPCNUM, status )
      if (status/=0) stop 'ERROR GAMMA_TNISP (__FILE__, line __LINE__)'
      ! fill:
      GAM_T = exp( TDF_PRM(SPCNUM)*(TEMP-Ts) )

      RETURN
      END SUBROUTINE GAMMA_TNISP
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C.....5) Calculate GAM_A (GAMMA_age)
C-----------------------------------------------------------------------
C
C             GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
C           where Fnew = new foliage fraction
C                 Fgro = growing foliage fraction
C                 Fmat = mature foliage fraction
C                 Fold = old foliage fraction
C                 Anew = relative emission activity for new foliage
C                 Agro = relative emission activity for growing foliage
C                 Amat = relative emission activity for mature foliage
C                 Aold = relative emission activity for old foliage
C
C
C             For foliage fraction
C             Case 1) LAIc = LAIp
C             Fnew = 0.0  , Fgro = 0.1  , Fmat = 0.8  , Fold = 0.1
C
C             Case 2) LAIp > LAIc
C             Fnew = 0.0  , Fgro = 0.0
C             Fmat = 1-Fold
C             Fold = (LAIp-LAIc)/LAIp
C
C             Case 3) LAIp < LAIc
C             Fnew = 1-(LAIp/LAIc)                       t <= ti
C                  = (ti/t) * ( 1-(LAIp/LAIc) )          t >  ti
C
C             Fmat = LAIp/LAIc                           t <= tm
C                  = (LAIp/LAIc) +
C                      ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )  t >  tm
C
C             Fgro = 1 - Fnew - Fmat
C             Fold = 0.0
C
C           where
C             ti = 5 + (0.7*(300-Tt))                   Tt <= 303
C                = 2.9                                  Tt >  303
C             tm = 2.3*ti
C
C             t  = length of the time step (days)
C             ti = number of days between budbreak and the induction of
C                  emission
C             tm = number of days between budbreak and the initiation of
C                  peak emissions rates
C             Tt = average temperature (K) near top of the canopy during
C                  current time period (daily ave temp for this case)
C
C
C             For relative emission activity
C             Case 1) Constant
C             Anew = 1.0  , Agro = 1.0  , Amat = 1.0  , Aold = 1.0
C
C             Case 2) Monoterpenes
C             Anew = 2.0  , Agro = 1.8  , Amat = 0.95 , Aold = 1.0
C
C             Case 3) Sesquiterpenes
C             Anew = 0.4  , Agro = 0.6  , Amat = 1.075, Aold = 1.0
C
C             Case 4) Methanol
C             Anew = 3.0  , Agro = 2.6  , Amat = 0.85 , Aold = 1.0
C
C             Case 5) Isoprene
C             Anew = 0.05 , Agro = 0.6  , Amat = 1.125, Aold = 1.0
C
C     SUBROUTINE GAMMA_A returns GAMMA_A
C-----------------------------------------------------------------------
      SUBROUTINE GAMMA_A( JDATE, JTIME, NCOLS, NROWS, SPC_NAME,
     &                    LAIARp, LAIARc, TSTLEN,
     &                    D_TEMP, GAM_A )

      IMPLICIT NONE

      INCLUDE 'megan_2_04_REL_EM_ACT.inc'

      CHARACTER*16  :: FUNCNAME = 'GAMMA_A'

C...  External parameters
      INTEGER     NCOLS, NROWS, JDATE, JTIME
      CHARACTER*8 SPC_NAME
      REAL        D_TEMP(NCOLS,NROWS)
      REAL        LAIARp(NCOLS,NROWS), LAIARc(NCOLS,NROWS)
      INTEGER     TSTLEN
      REAL        GAM_A(NCOLS,NROWS)

C...  Local parameters
      REAL  Fnew, Fgro, Fmat, Fold

      INTEGER        AINDX          ! relative emission acitivity index
      CHARACTER*256  MESG           ! message buffer

      REAL     LAIp                 ! LAI at previous time step
      REAL     LAIc                 ! LAI at current time step
      INTEGER  t                    ! time step
      REAL     ti                   ! number of days between budbreak
                                    ! and the induction of emission
      REAL     tm                   ! number of days between budbreak
                                    ! and the initiation of peak
                                    ! emissions rates
      REAL     Tt                   ! average temperature (K)
                                    ! daily ave temp

      INTEGER  I, J                 ! loop index

C...  Choose relative emission activity
      SELECT CASE ( TRIM(SPC_NAME) )
      CASE ('ACTO','ACTA','FORM','CH4','NO','CO')
         AINDX = 1

      CASE ('MYRC','SABI','LIMO','3CAR','OCIM','BPIN','APIN','OMTP')
         AINDX = 2

      CASE ('FARN','BCAR','OSQT')
         AINDX = 3

      CASE ('MEOH')
         AINDX = 4

      CASE ( 'ISOP','MBO' )
         AINDX = 5

      CASE DEFAULT
!         MESG = 'Error: Chemical species, invalid variable: '
!     &          //TRIM(SPC_NAME)
!         CALL M3EXIT(FUNCNAME,JDATE,JTIME,MESG,2)
         write (*,'("Error: Chemical species, invalid variable: ",a)') TRIM(SPC_NAME)
         write (*,'("in ",a," (",a,", line",i5,")")') FUNCNAME, __FILE__, __LINE__
         call Exit(2)

      ENDSELECT

      t = TSTLEN
      DO I = 1, NCOLS
         DO J = 1, NROWS
            LAIc = LAIARc(I,J)
            LAIp = LAIARp(I,J)
            Tt   = D_TEMP(I,J)
C...  Calculate foliage fraction
            IF (LAIp .EQ. LAIc) THEN
               Fnew = 0.0
               Fgro = 0.1
               Fmat = 0.8
               Fold = 0.1

            ELSEIF (LAIp .GT. LAIc) THEN
               Fnew = 0.0
               Fgro = 0.0
               Fold = ( LAIp-LAIc ) / LAIp
               Fmat = 1-Fold

            ELSEIF (LAIp .LT. LAIc) THEN
C              Calculate ti and tm
               IF (Tt .LE. 303.0) THEN
                  ti = 5.0 + 0.7*(300-Tt)
               ELSEIF (Tt .GT. 303.0) THEN
                  ti = 2.9
               ENDIF
               tm = 2.3*ti

C              Calculate Fnew and Fmat, then Fgro and Fold
C              Fnew
               IF (t .LE. ti) THEN
                  Fnew = 1.0 - (LAIp/LAIc)
               ELSEIF (t .GT. ti) THEN
                  Fnew = (ti/t) * ( 1-(LAIp/LAIc) )
               ENDIF

C              Fmat
               IF (t .LE. tm) THEN
                  Fmat = LAIp/LAIc
               ELSEIF (t .GT. tm) THEN
                  Fmat = (LAIp/LAIc) + ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )
               ENDIF

               Fgro = 1.0 - Fnew - Fmat
               Fold = 0.0
         
            ENDIF

C...  Calculate GAMMA_A
      GAM_A(I,J) = Fnew*Anew(AINDX) + Fgro*Agro(AINDX) +
     &             Fmat*Amat(AINDX) + Fold*Aold(AINDX)

         ENDDO    ! End loop for NROWS
      ENDDO       ! End loop for NCOLS

      RETURN
      END SUBROUTINE GAMMA_A


C-----------------------------------------------------------------------
C.....6) Calculate GAM_SMT (GAMMA_SM)
C-----------------------------------------------------------------------
C
C             GAMMA_SM =     1.0   (non-dimension)
C
C
C     SUBROUTINE GAMMA_S returns the GAMMA_SM values
C-----------------------------------------------------------------------
      SUBROUTINE GAMMA_S( NCOLS, NROWS, GAM_S )

      IMPLICIT NONE

      INTEGER NCOLS, NROWS
      REAL    GAM_S(NCOLS,NROWS)

      GAM_S = 1.0

      RETURN
      END SUBROUTINE GAMMA_S


C=======================================================================
      END MODULE GAMMA_ETC
