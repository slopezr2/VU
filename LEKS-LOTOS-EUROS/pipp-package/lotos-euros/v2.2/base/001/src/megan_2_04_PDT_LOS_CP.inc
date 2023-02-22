C=======================================================================
C  PDT_LOS_CP.EXT
C  This include file contains "produciton and loss within canopy"
C  factors.
C
C
C  MEGAN v2.02
C  INPUT version 200
C
C  History:
C  Who          When       What
C  ---------------------------------------------------------------------
C  Tan          12/02/06 - Creates this file
C  Tan          08/14/07 - Move from MEGAN v2.0 to MEGAN v2.02 with no update.
C=======================================================================

      INTEGER        N_RHO_SPC
      PARAMETER     (N_RHO_SPC  = 20)
      CHARACTER*16   RHO_SPC(N_RHO_SPC)
      REAL           RHO_FCT(N_RHO_SPC)
      INTEGER        RHO_MAP(N_RHO_SPC)

      DATA     RHO_SPC(  1)      , RHO_FCT(  1), RHO_MAP(  1)
     &       / 'ISOP            ',   1.0       , 1           /

      DATA     RHO_SPC(  2)      , RHO_FCT(  2), RHO_MAP(  2)
     &       / 'MYRC            ',   1.0       , 2           /

      DATA     RHO_SPC(  3)      , RHO_FCT(  3), RHO_MAP(  3)
     &       / 'SABI            ',   1.0       , 3           /

      DATA     RHO_SPC(  4)      , RHO_FCT(  4), RHO_MAP(  4)
     &       / 'LIMO            ',   1.0       , 4           /

      DATA     RHO_SPC(  5)      , RHO_FCT(  5), RHO_MAP(  5)
     &       / '3CAR            ',   1.0       , 5           /

      DATA     RHO_SPC(  6)      , RHO_FCT(  6), RHO_MAP(  6)
     &       / 'OCIM            ',   1.0       , 6           /

      DATA     RHO_SPC(  7)      , RHO_FCT(  7), RHO_MAP(  7)
     &       / 'BPIN            ',   1.0       , 7           /

      DATA     RHO_SPC(  8)      , RHO_FCT(  8), RHO_MAP(  8)
     &       / 'APIN            ',   1.0       , 8           /

      DATA     RHO_SPC(  9)      , RHO_FCT(  9), RHO_MAP(  9)
     &       / 'OMTP            ',   1.0       , 9           /

      DATA     RHO_SPC( 10)      , RHO_FCT( 10), RHO_MAP( 10)
     &       / 'FARN            ',   1.0       , 10          /

      DATA     RHO_SPC( 11)      , RHO_FCT( 11), RHO_MAP( 11)
     &       / 'BCAR            ',   1.0       , 11          /

      DATA     RHO_SPC( 12)      , RHO_FCT( 12), RHO_MAP( 12)
     &       / 'OSQT            ',   1.0       , 12          /

      DATA     RHO_SPC( 13)      , RHO_FCT( 13), RHO_MAP( 13)
     &       / 'MBO             ',   1.0       , 13          /

      DATA     RHO_SPC( 14)      , RHO_FCT( 14), RHO_MAP( 14)
     &       / 'MEOH            ',   1.0       , 14          /

      DATA     RHO_SPC( 15)      , RHO_FCT( 15), RHO_MAP( 15)
     &       / 'ACTO            ',   1.0       , 15          /

      DATA     RHO_SPC( 16)      , RHO_FCT( 16), RHO_MAP( 16)
     &       / 'CH4             ',   1.0       , 16          /

      DATA     RHO_SPC( 17)      , RHO_FCT( 17), RHO_MAP( 17)
     &       / 'NO              ',   1.0       , 17          /

      DATA     RHO_SPC( 18)      , RHO_FCT( 18), RHO_MAP( 18)
     &       / 'ACTA            ',   1.0       , 18          /

      DATA     RHO_SPC( 19)      , RHO_FCT( 19), RHO_MAP( 19)
     &       / 'FORM            ',   1.0       , 19          /

      DATA     RHO_SPC( 20)      , RHO_FCT( 20), RHO_MAP( 20)
     &       / 'CO              ',   1.0       , 20          /
