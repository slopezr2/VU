C=======================================================================
C  SPC_CBM4.EXT
C  This include file contains CBM4 species and their MW.
C
C
C  Mechanism Name: CBMZ
C  MEGAN v2.02
C  INPUT version 200
C
C  History:
C  Who          When       What
C  ---------------------------------------------------------------------
C  Tan          12/02/06 - Creates this file
C  Tan          08/14/07 - Move from MEGAN v2.0 to MEGAN v2.02 with no update.
C  Arjo Segers  2012/08  - Created based on CBMZ version.
C                          Added 'TERP' component for terpenes.
C=======================================================================

      CHARACTER*16   MECH_SPC( N_CBM4_SPC )  ! Mechanism species name
      REAL           MECH_MWT( N_CBM4_SPC )  ! Mechanism species molecular weight


      DATA  MECH_SPC(  1), MECH_MWT(  1) / 'ISOP            ', 68.12   /
      DATA  MECH_SPC(  2), MECH_MWT(  2) / 'NO              ', 30.01   /
      DATA  MECH_SPC(  3), MECH_MWT(  3) / 'NO2             ', 44.01   /
      DATA  MECH_SPC(  4), MECH_MWT(  4) / 'CO              ', 28.01   /
      DATA  MECH_SPC(  5), MECH_MWT(  5) / 'CH3OH           ', 32.04   /
      DATA  MECH_SPC(  6), MECH_MWT(  6) / 'ANOL            ', 46.07   /
      DATA  MECH_SPC(  7), MECH_MWT(  7) / 'ALD2            ', 44.05   /
      DATA  MECH_SPC(  8), MECH_MWT(  8) / 'HCHO            ', 30.03   /
      DATA  MECH_SPC(  9), MECH_MWT(  9) / 'HCOOH           ', 46.03   /
      DATA  MECH_SPC( 10), MECH_MWT( 10) / 'RCOOH           ', 60.05   /
      DATA  MECH_SPC( 11), MECH_MWT( 11) / 'CH4             ', 16.04   /
      DATA  MECH_SPC( 12), MECH_MWT( 12) / 'C2H6            ', 30.07   /
      DATA  MECH_SPC( 13), MECH_MWT( 13) / 'PAR             ', 13.00   /
      DATA  MECH_SPC( 14), MECH_MWT( 14) / 'ETH             ', 28.05   /
      DATA  MECH_SPC( 15), MECH_MWT( 15) / 'OLET            ', 27.00   /
      DATA  MECH_SPC( 16), MECH_MWT( 16) / 'OLEI            ', 26.00   /
      DATA  MECH_SPC( 17), MECH_MWT( 17) / 'CRES            ', 108.14  /
      DATA  MECH_SPC( 18), MECH_MWT( 18) / 'AONE            ', 58.08   /
      DATA  MECH_SPC( 19), MECH_MWT( 19) / 'TOL             ', 92.14   /
      DATA  MECH_SPC( 20), MECH_MWT( 20) / 'XYL             ', 106.17  /
      DATA  MECH_SPC( 21), MECH_MWT( 21) / 'DMS             ', 62.14   /
      DATA  MECH_SPC( 22), MECH_MWT( 22) / 'NH3             ', 17.03   /
      DATA  MECH_SPC( 23), MECH_MWT( 23) / 'NR              ', 1.00    /
      ! terpenes (C10H6) :
      DATA  MECH_SPC( 24), MECH_MWT( 24) / 'TERP            ', 136.00  /
