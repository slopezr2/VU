C=======================================================================
C  LD_FCT.EXT
C  This include file contains "light dependent" factors.
C
C
C
C  MEGAN v2.02
C  INPUT version 210
C
C  History:
C  Who          When       What
C  ---------------------------------------------------------------------
C  Tan          12/02/06 - Creates this file
C  Guenther A.  08/11/07 - Move from MEGAN v2.0 to MEGAN v2.02 with updates.
C                          See the update document.
C=======================================================================

      INTEGER        N_LDF_SPC
      PARAMETER     (N_LDF_SPC  = 20)
      CHARACTER*16   LDF_SPC(N_LDF_SPC)
      REAL           LDF_FCT(N_LDF_SPC)
      INTEGER        LDF_MAP(N_LDF_SPC)

      DATA     LDF_SPC(  1)      , LDF_FCT(  1), LDF_MAP(  1)
     &       / 'ISOP            ', 0.9999      , 1            /

      DATA     LDF_SPC(  2)      , LDF_FCT(  2), LDF_MAP(  2)
     &       / 'MYRC            ', 0.05        , 2            /

      DATA     LDF_SPC(  3)      , LDF_FCT(  3), LDF_MAP(  3)
     &       / 'SABI            ', 0.1         , 3            /

      DATA     LDF_SPC(  4)      , LDF_FCT(  4), LDF_MAP(  4)
     &       / 'LIMO            ', 0.05        , 4            /

      DATA     LDF_SPC(  5)      , LDF_FCT(  5), LDF_MAP(  5)
     &       / '3CAR            ', 0.05        , 5            /

      DATA     LDF_SPC(  6)      , LDF_FCT(  6), LDF_MAP(  6)
     &       / 'OCIM            ', 0.8         , 6            /

      DATA     LDF_SPC(  7)      , LDF_FCT(  7), LDF_MAP(  7)
     &       / 'BPIN            ', 0.1         , 7            /

      DATA     LDF_SPC(  8)      , LDF_FCT(  8), LDF_MAP(  8)
     &       / 'APIN            ', 0.1         , 8            /

      DATA     LDF_SPC(  9)      , LDF_FCT(  9), LDF_MAP(  9)
     &       / 'OMTP            ', 0.1         , 9            /

      DATA     LDF_SPC( 10)      , LDF_FCT( 10), LDF_MAP( 10)
     &       / 'FARN            ', 0.5         , 10           /

      DATA     LDF_SPC( 11)      , LDF_FCT( 11), LDF_MAP( 11)
     &       / 'BCAR            ', 0.5         , 11           /

      DATA     LDF_SPC( 12)      , LDF_FCT( 12), LDF_MAP( 12)
     &       / 'OSQT            ', 0.5         , 12           /

      DATA     LDF_SPC( 13)      , LDF_FCT( 13), LDF_MAP( 13)
     &       / 'MBO             ', 0.9999      , 13           /

      DATA     LDF_SPC( 14)      , LDF_FCT( 14), LDF_MAP( 14)
     &       / 'MEOH            ', 0.75        , 14           /

      DATA     LDF_SPC( 15)      , LDF_FCT( 15), LDF_MAP( 15)
     &       / 'ACTO            ', 0.25        , 15           /

      DATA     LDF_SPC( 16)      , LDF_FCT( 16), LDF_MAP( 16)
     &       / 'CH4             ', 0.75        , 16           /

      DATA     LDF_SPC( 17)      , LDF_FCT( 17), LDF_MAP( 17)
     &       / 'NO              ', 0.0         , 17           /

      DATA     LDF_SPC( 18)      , LDF_FCT( 18), LDF_MAP( 18)
     &       / 'ACTA            ', 0.5         , 18           /

      DATA     LDF_SPC( 19)      , LDF_FCT( 19), LDF_MAP( 19)
     &       / 'FORM            ', 0.5         , 19           /

      DATA     LDF_SPC( 20)      , LDF_FCT( 20), LDF_MAP( 20)
     &       / 'CO              ', 0.5         , 20           /
