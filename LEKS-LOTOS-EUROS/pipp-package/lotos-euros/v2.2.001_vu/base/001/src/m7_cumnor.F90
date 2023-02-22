SUBROUTINE m7_cumnor ( arg, RESULT, ccum )
  !
  !*******************************************************************************
  !
  !! CUMNOR computes the cumulative normal distribution.
  !
  !
  !     the integral from -infinity to x of
  !          (1/sqrt(2*pi)) exp(-u*u/2) du
  !
  !  Author:
  !  -------
  !  Original source:
  !
  !    W. J. Cody    Mathematics and Computer Science Division
  !                  Argonne National Laboratory
  !                  Argonne, IL 60439
  !
  !    DCDFLIB is attributed to Barry Brown, James Lovato, and Kathy Russell
  !            bwb@odin.mda.uth.tmc.edu.
  !
  !    Adopted to ECHAM/M7:
  !
  !    Philip Stier  (MPI-MET)                    2001
  !
  !
  !  Reference:
  !  ----------
  !
  !    W D Cody, 
  !    "ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special 
  !    Function Routines and Test Drivers"
  !    ACM Transactions on Mathematical Software,
  !    Volume 19, 1993, pages 22-32.
  !
  !  Parameters:
  !
  !     ARG --> Upper limit of integration.
  !                                        X is REAL
  !
  !     RESULT <-- Cumulative normal distribution.
  !                                        RESULT is REAL
  !
  !     CCUM <-- Complement of Cumulative normal distribution.
  !                                        CCUM is REAL
  !
  !
  ! Original Comments:
  !
  !
  ! This function evaluates the normal distribution function:
  !
  !                              / x
  !                     1       |       -t*t/2
  !          P(x) = ----------- |      e       dt
  !                 sqrt(2 pi)  |
  !                             /-oo
  !
  !   The main computation evaluates near-minimax approximations
  !   derived from those in "Rational Chebyshev approximations for
  !   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
  !   This transportable program uses rational functions that
  !   theoretically approximate the normal distribution function to
  !   at least 18 significant decimal digits.  The accuracy achieved
  !   depends on the arithmetic system, the compiler, the intrinsic
  !   functions, and proper selection of the machine-dependent
  !   constants.
  !
  !  Explanation of machine-dependent constants.
  !
  !   MIN   = smallest machine representable number.
  !
  !   EPS   = argument below which anorm(x) may be represented by
  !           0.5  and above which  x*x  will not underflow.
  !           A conservative value is the largest machine number X
  !           such that   1.0 + X = 1.0   to machine precision.
  !
  !  Error returns
  !
  !  The program returns  ANORM = 0     for  ARG .LE. XLOW.
  !
  !  Author: 
  !
  !    W. J. Cody
  !    Mathematics and Computer Science Division
  !    Argonne National Laboratory
  !    Argonne, IL 60439
  !
  !  Latest modification: March 15, 1992
  !
  REAL, PARAMETER, DIMENSION ( 5 ) :: a = (/ &
       2.2352520354606839287d00, &
       1.6102823106855587881d02, &
       1.0676894854603709582d03, &
       1.8154981253343561249d04, &
       6.5682337918207449113d-2 /)
  
  real arg
  REAL, PARAMETER, DIMENSION ( 4 ) :: b = (/ &
       4.7202581904688241870d01, &
       9.7609855173777669322d02, &
       1.0260932208618978205d04, &
       4.5507789335026729956d04 /)
  REAL, PARAMETER, DIMENSION ( 9 ) :: c = (/ &
       3.9894151208813466764d-1, &
       8.8831497943883759412d00, &
       9.3506656132177855979d01, &
       5.9727027639480026226d02, &
       2.4945375852903726711d03, &
       6.8481904505362823326d03, &
       1.1602651437647350124d04, &
       9.8427148383839780218d03, &
       1.0765576773720192317d-8 /)
  REAL ccum
  REAL, PARAMETER, DIMENSION ( 8 ) :: d = (/ &
       2.2266688044328115691d01, &
       2.3538790178262499861d02, &
       1.5193775994075548050d03, &
       6.4855582982667607550d03, &
       1.8615571640885098091d04, &
       3.4900952721145977266d04, &
       3.8912003286093271411d04, &
       1.9685429676859990727d04 /)
  REAL del
!@@@ REAL dpmpar
  REAL eps
  INTEGER i
  REAL min
  REAL, PARAMETER, DIMENSION ( 6 ) :: p = (/ &
       2.1589853405795699d-1, &
       1.274011611602473639d-1, &
       2.2235277870649807d-2, &
       1.421619193227893466d-3, &
       2.9112874951168792d-5, &
       2.307344176494017303d-2 /)
  REAL, PARAMETER, DIMENSION ( 5 ) :: q = (/ &
       1.28426009614491121d00, &
       4.68238212480865118d-1, &
       6.59881378689285515d-2, &
       3.78239633202758244d-3, &
       7.29751555083966205d-5 /)
  REAL RESULT
  REAL, PARAMETER :: root32 = 5.656854248d0
  REAL, PARAMETER :: sixten = 16.0
  REAL temp
  REAL, PARAMETER :: sqrpi = 3.9894228040143267794d-1
  REAL, PARAMETER :: thrsh = 0.66291d0
  REAL x
  REAL xden
  REAL xnum
  REAL y
  REAL xsq
  !
  !  Machine dependent constants
  !
  eps = EPSILON ( 1.0d0 ) * 0.5d0
  !
  !@@@ Simplified calculation of the smallest machine representable number
  !    (Higher accuracy than needed!)
  !
  !@@@ min = dpmpar(2)

  min = epsilon ( 1.0D0)
  

  x = arg
  y = ABS ( x )

  IF ( y <= thrsh ) THEN
     !
     !  Evaluate  anorm  for  |X| <= 0.66291
     !
     IF ( y > eps ) THEN
        xsq = x * x
     ELSE
        xsq = 0.0
     END IF

     xnum = a(5) * xsq
     xden = xsq
     DO i = 1, 3
        xnum = ( xnum + a(i) ) * xsq
        xden = ( xden + b(i) ) * xsq
     END DO
     RESULT = x * ( xnum + a(4) ) / ( xden + b(4) )
     temp = RESULT
     RESULT = 0.5 + temp
     ccum = 0.5 - temp
     !
     !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
     !
  ELSE IF ( y <= root32 ) THEN

     xnum = c(9) * y
     xden = y
!CDIR UNROLL=7
     DO i = 1, 7
        xnum = ( xnum + c(i) ) * y
        xden = ( xden + d(i) ) * y
     END DO
     RESULT = ( xnum + c(8) ) / ( xden + d(8) )
     xsq = AINT ( y * sixten ) / sixten
     del = ( y - xsq ) * ( y + xsq )
     RESULT = EXP(-xsq*xsq*0.5) * EXP(-del*0.5) * RESULT
     ccum = 1.0 - RESULT

     IF ( x > 0.0 ) THEN
        temp = RESULT
        RESULT = ccum
        ccum = temp
     END IF
     !
     !  Evaluate  anorm  for |X| > sqrt(32).
     !
  ELSE

     RESULT = 0.0
     xsq = 1.0 / ( x * x )
     xnum = p(6) * xsq
     xden = xsq
     DO i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
     END DO

     RESULT = xsq * ( xnum + p(5) ) / ( xden + q(5) )
     RESULT = ( sqrpi - RESULT ) / y
     xsq = AINT ( x * sixten ) / sixten
     del = ( x - xsq ) * ( x + xsq )
     RESULT = EXP ( - xsq * xsq * 0.5 ) * EXP ( - del * 0.5 ) * RESULT
     ccum = 1.0 - RESULT

     IF ( x > 0.0 ) THEN
        temp = RESULT
        RESULT = ccum
        ccum = temp
     END IF

  END IF

  IF ( RESULT < min ) THEN
     RESULT = 0.0d0
  END IF

  IF ( ccum < min ) THEN
     ccum = 0.0d0
  END IF

END SUBROUTINE m7_cumnor
