C=======================================================================
C  TEMPD_PRM.EXT
C  This include file contains "temperature dependent" parameter for
C  light-independent emissions.
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

      INTEGER        N_TDF_SPC
      PARAMETER     (N_TDF_SPC  = 20)
      CHARACTER*16   TDF_SPC(N_TDF_SPC)
      REAL           TDF_PRM(N_TDF_SPc)
      INTEGER        TDF_MAP(N_TDF_SPc)

      DATA     TDF_SPC(  1)      , TDF_PRM(  1), TDF_MAP(  1) 
     &       / 'ISOP            ', 0.09        , 1           /

      DATA     TDF_SPC(  2)      , TDF_PRM(  2), TDF_MAP(  2) 
     &       / 'MYRC            ', 0.1         , 2           /

      DATA     TDF_SPC(  3)      , TDF_PRM(  3), TDF_MAP(  3) 
     &       / 'SABI            ', 0.1         , 3           /

      DATA     TDF_SPC(  4)      , TDF_PRM(  4), TDF_MAP(  4) 
     &       / 'LIMO            ', 0.1         , 4           /

      DATA     TDF_SPC(  5)      , TDF_PRM(  5), TDF_MAP(  5) 
     &       / '3CAR            ', 0.1         , 5           /

      DATA     TDF_SPC(  6)      , TDF_PRM(  6), TDF_MAP(  6) 
     &       / 'OCIM            ', 0.1         , 6           /

      DATA     TDF_SPC(  7)      , TDF_PRM(  7), TDF_MAP(  7) 
     &       / 'BPIN            ', 0.1         , 7           /

      DATA     TDF_SPC(  8)      , TDF_PRM(  8), TDF_MAP(  8) 
     &       / 'APIN            ', 0.1         , 8           /

      DATA     TDF_SPC(  9)      , TDF_PRM(  9), TDF_MAP(  9) 
     &       / 'OMTP            ', 0.1         , 9           /

      DATA     TDF_SPC( 10)      , TDF_PRM( 10), TDF_MAP( 10) 
     &       / 'FARN            ', 0.17        , 10          /

      DATA     TDF_SPC( 11)      , TDF_PRM( 11), TDF_MAP( 11) 
     &       / 'BCAR            ', 0.17        , 11          /

      DATA     TDF_SPC( 12)      , TDF_PRM( 12), TDF_MAP( 12) 
     &       / 'OSQT            ', 0.17        , 12          /

      DATA     TDF_SPC( 13)      , TDF_PRM( 13), TDF_MAP( 13) 
     &       / 'MBO             ', 0.09        , 13          /

      DATA     TDF_SPC( 14)      , TDF_PRM( 14), TDF_MAP( 14) 
     &       / 'MEOH            ', 0.08        , 14          /

      DATA     TDF_SPC( 15)      , TDF_PRM( 15), TDF_MAP( 15) 
     &       / 'ACTO            ', 0.11        , 15          /

      DATA     TDF_SPC( 16)      , TDF_PRM( 16), TDF_MAP( 16) 
     &       / 'CH4             ', 0.05        , 16          /

      DATA     TDF_SPC( 17)      , TDF_PRM( 17), TDF_MAP( 17) 
     &       / 'NO              ', 0.11        , 17          /

      DATA     TDF_SPC( 18)      , TDF_PRM( 18), TDF_MAP( 18) 
     &       / 'ACTA            ', 0.13        , 18          /

      DATA     TDF_SPC( 19)      , TDF_PRM( 19), TDF_MAP( 19) 
     &       / 'FORM            ', 0.09        , 19          /

      DATA     TDF_SPC( 20)      , TDF_PRM( 20), TDF_MAP( 20) 
     &       / 'CO              ', 0.09        , 20          /
