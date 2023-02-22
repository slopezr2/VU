C=======================================================================
C  REL_EM_ACT.EXT
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

      INTEGER        N_CAT
      PARAMETER     (N_CAT  = 5)
      REAL           Anew(N_CAT)
      REAL           Agro(N_CAT)
      REAL           Amat(N_CAT)
      REAL           Aold(N_CAT)

      DATA    Anew(  1),  Agro(  1),  Amat(  1),  Aold(  1)
     &     /  1.0      ,  1.0      ,  1.0      ,  1.0       /

      DATA    Anew(  2),  Agro(  2),  Amat(  2),  Aold(  2)
     &     /  2.0      ,  1.8      ,  0.95     ,  1.0       /

      DATA    Anew(  3),  Agro(  3),  Amat(  3),  Aold(  3)
     &     /  0.4      ,  0.6      ,  1.075    ,  1.0       /

      DATA    Anew(  4),  Agro(  4),  Amat(  4),  Aold(  4)
     &     /  3.0      ,  2.6      ,  0.85     ,  1.0       /

      DATA    Anew(  5),  Agro(  5),  Amat(  5),  Aold(  5)
     &     /  0.05     ,  0.6      ,  1.125    ,  1.0       /
