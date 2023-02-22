C=======================================================================
C  MAP_CV2CBM4.EXT
C  This include file contains conversion table for 138 speciated species
C  to CBM4 species
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
C  Arjo Segers  2012/08  - Based on mapping for CBZ. 
C                          Replaced mapping from '*terp*' to two ISOP molecules
C                          by mapping to single TERP molecule
C                          (monoterpene consists of two isoprene units)
C=======================================================================

      INTEGER        N_SCON_SPC
      PARAMETER     (N_SCON_SPC = 173)        ! Number of map species

      CHARACTER*16   SPMH_NAM( N_SCON_SPC )   ! speciated species name
      INTEGER        SPMH_MAP( N_SCON_SPC )   ! speciated species name
                                              ! mapped to SPC_SPCAT.EXT
      CHARACTER*16   MECH_NAM( N_SCON_SPC )   ! mechanism species
      INTEGER        MECH_MAP( N_SCON_SPC )   ! mechanism species mapped
                                              ! to SPC_CBMZ.EXT
      REAL           CONV_FAC( N_SCON_SPC )   ! conversion factor


      DATA  SPMH_NAM(  1)     , SPMH_MAP(  1)    / 'isoprene        ', 1              /
      DATA  MECH_NAM(  1)     , MECH_MAP(  1)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  1)    / 1.0            /
      DATA  SPMH_NAM(  2)     , SPMH_MAP(  2)    / 'myrcene         ', 2              /
      DATA  MECH_NAM(  2)     , MECH_MAP(  2)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  2)    / 2.0            /
      DATA  SPMH_NAM(  3)     , SPMH_MAP(  3)    / 'sabinene        ', 3              /
      DATA  MECH_NAM(  3)     , MECH_MAP(  3)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  3)    / 2.0            /
      DATA  SPMH_NAM(  4)     , SPMH_MAP(  4)    / 'limonene        ', 4              /
      DATA  MECH_NAM(  4)     , MECH_MAP(  4)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  4)    / 2.0            /
      DATA  SPMH_NAM(  5)     , SPMH_MAP(  5)    / 'carene_3        ', 5              /
      DATA  MECH_NAM(  5)     , MECH_MAP(  5)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  5)    / 2.0            /
      DATA  SPMH_NAM(  6)     , SPMH_MAP(  6)    / 'ocimene_t_b     ', 6              /
      DATA  MECH_NAM(  6)     , MECH_MAP(  6)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  6)    / 2.0            /
      DATA  SPMH_NAM(  7)     , SPMH_MAP(  7)    / 'pinene_b        ', 7              /
      DATA  MECH_NAM(  7)     , MECH_MAP(  7)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  7)    / 2.0            /
      DATA  SPMH_NAM(  8)     , SPMH_MAP(  8)    / 'pinene_a        ', 8              /
      DATA  MECH_NAM(  8)     , MECH_MAP(  8)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  8)    / 2.0            /
      DATA  SPMH_NAM(  9)     , SPMH_MAP(  9)    / '2met_styrene    ', 9              /
      DATA  MECH_NAM(  9)     , MECH_MAP(  9)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(  9)    / 2.0            /
      DATA  SPMH_NAM( 10)     , SPMH_MAP( 10)    / 'cymene_p        ', 10             /
      DATA  MECH_NAM( 10)     , MECH_MAP( 10)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 10)    / 2.0            /
      DATA  SPMH_NAM( 11)     , SPMH_MAP( 11)    / 'cymene_o        ', 11             /
      DATA  MECH_NAM( 11)     , MECH_MAP( 11)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 11)    / 2.0            /
      DATA  SPMH_NAM( 12)     , SPMH_MAP( 12)    / 'phellandrene_a  ', 12             /
      DATA  MECH_NAM( 12)     , MECH_MAP( 12)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 12)    / 2.0            /
      DATA  SPMH_NAM( 13)     , SPMH_MAP( 13)    / 'thujene_a       ', 13             /
      DATA  MECH_NAM( 13)     , MECH_MAP( 13)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 13)    / 2.0            /
      DATA  SPMH_NAM( 14)     , SPMH_MAP( 14)    / 'terpinene_a     ', 14             /
      DATA  MECH_NAM( 14)     , MECH_MAP( 14)    / 'TERP            ', 24             /
      DATA  CONV_FAC( 14)    / 1.0            /
      DATA  SPMH_NAM( 15)     , SPMH_MAP( 15)    / 'terpinene_g     ', 15             /
      DATA  MECH_NAM( 15)     , MECH_MAP( 15)    / 'TERP            ', 24             /
      DATA  CONV_FAC( 15)    / 1.0            /
      DATA  SPMH_NAM( 16)     , SPMH_MAP( 16)    / 'terpinolene     ', 16             /
      DATA  MECH_NAM( 16)     , MECH_MAP( 16)    / 'TERP            ', 24             /
      DATA  CONV_FAC( 16)    / 1.0            /
      DATA  SPMH_NAM( 17)     , SPMH_MAP( 17)    / 'phellandrene_b  ', 17             /
      DATA  MECH_NAM( 17)     , MECH_MAP( 17)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 17)    / 2.0            /
      DATA  SPMH_NAM( 18)     , SPMH_MAP( 18)    / 'camphene        ', 18             /
      DATA  MECH_NAM( 18)     , MECH_MAP( 18)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 18)    / 2.0            /
      DATA  SPMH_NAM( 19)     , SPMH_MAP( 19)    / 'bornene         ', 19             /
      DATA  MECH_NAM( 19)     , MECH_MAP( 19)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 19)    / 2.0            /
      DATA  SPMH_NAM( 20)     , SPMH_MAP( 20)    / 'fenchene_a      ', 20             /
      DATA  MECH_NAM( 20)     , MECH_MAP( 20)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 20)    / 2.0            /
      DATA  SPMH_NAM( 21)     , SPMH_MAP( 21)    / 'ocimene_al      ', 21             /
      DATA  MECH_NAM( 21)     , MECH_MAP( 21)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 21)    / 2.0            /
      DATA  SPMH_NAM( 22)     , SPMH_MAP( 22)    / 'ocimene_c_b     ', 22             /
      DATA  MECH_NAM( 22)     , MECH_MAP( 22)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 22)    / 2.0            /
      DATA  SPMH_NAM( 23)     , SPMH_MAP( 23)    / 'tricyclene      ', 23             /
      DATA  MECH_NAM( 23)     , MECH_MAP( 23)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 23)    / 2.0            /
      DATA  SPMH_NAM( 24)     , SPMH_MAP( 24)    / 'estragole       ', 24             /
      DATA  MECH_NAM( 24)     , MECH_MAP( 24)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 24)    / 2.0            /
      DATA  SPMH_NAM( 25)     , SPMH_MAP( 25)    / 'camphor         ', 25             /
      DATA  MECH_NAM( 25)     , MECH_MAP( 25)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 25)    / 2.0            /
      DATA  SPMH_NAM( 26)     , SPMH_MAP( 26)    / 'fenchone        ', 26             /
      DATA  MECH_NAM( 26)     , MECH_MAP( 26)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 26)    / 2.0            /
      DATA  SPMH_NAM( 27)     , SPMH_MAP( 27)    / 'piperitone      ', 27             /
      DATA  MECH_NAM( 27)     , MECH_MAP( 27)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 27)    / 2.0            /
      DATA  SPMH_NAM( 28)     , SPMH_MAP( 28)    / 'thujone_a       ', 28             /
      DATA  MECH_NAM( 28)     , MECH_MAP( 28)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 28)    / 2.0            /
      DATA  SPMH_NAM( 29)     , SPMH_MAP( 29)    / 'thujone_b       ', 29             /
      DATA  MECH_NAM( 29)     , MECH_MAP( 29)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 29)    / 2.0            /
      DATA  SPMH_NAM( 30)     , SPMH_MAP( 30)    / 'cineole_1_8     ', 30             /
      DATA  MECH_NAM( 30)     , MECH_MAP( 30)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 30)    / 2.0            /
      DATA  SPMH_NAM( 31)     , SPMH_MAP( 31)    / 'borneol         ', 31             /
      DATA  MECH_NAM( 31)     , MECH_MAP( 31)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 31)    / 2.0            /
      DATA  SPMH_NAM( 32)     , SPMH_MAP( 32)    / 'linalool        ', 32             /
      DATA  MECH_NAM( 32)     , MECH_MAP( 32)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 32)    / 2.0            /
      DATA  SPMH_NAM( 33)     , SPMH_MAP( 33)    / 'terpineol_4     ', 33             /
      DATA  MECH_NAM( 33)     , MECH_MAP( 33)    / 'TERP            ', 24             /
      DATA  CONV_FAC( 33)    / 1.0            /
      DATA  SPMH_NAM( 34)     , SPMH_MAP( 34)    / 'terpineol_a     ', 34             /
      DATA  MECH_NAM( 34)     , MECH_MAP( 34)    / 'TERP            ', 24             /
      DATA  CONV_FAC( 34)    / 1.0            /
      DATA  SPMH_NAM( 35)     , SPMH_MAP( 35)    / 'linalool_OXD_c  ', 35             /
      DATA  MECH_NAM( 35)     , MECH_MAP( 35)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 35)    / 2.0            /
      DATA  SPMH_NAM( 36)     , SPMH_MAP( 36)    / 'linalool_OXD_t  ', 36             /
      DATA  MECH_NAM( 36)     , MECH_MAP( 36)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 36)    / 2.0            /
      DATA  SPMH_NAM( 37)     , SPMH_MAP( 37)    / 'ionone_b        ', 37             /
      DATA  MECH_NAM( 37)     , MECH_MAP( 37)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 37)    / 3.0            /
      DATA  SPMH_NAM( 38)     , SPMH_MAP( 38)    / 'bornyl_ACT      ', 38             /
      DATA  MECH_NAM( 38)     , MECH_MAP( 38)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 38)    / 2.0            /
      DATA  SPMH_NAM( 39)     , SPMH_MAP( 39)    / 'farnescene_a    ', 39             /
      DATA  MECH_NAM( 39)     , MECH_MAP( 39)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 39)    / 3.0            /
      DATA  SPMH_NAM( 40)     , SPMH_MAP( 40)    / 'caryophyllene_b ', 40             /
      DATA  MECH_NAM( 40)     , MECH_MAP( 40)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 40)    / 3.0            /
      DATA  SPMH_NAM( 41)     , SPMH_MAP( 41)    / 'acoradiene      ', 41             /
      DATA  MECH_NAM( 41)     , MECH_MAP( 41)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 41)    / 3.0            /
      DATA  SPMH_NAM( 42)     , SPMH_MAP( 42)    / 'aromadendrene   ', 42             /
      DATA  MECH_NAM( 42)     , MECH_MAP( 42)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 42)    / 3.0            /
      DATA  SPMH_NAM( 43)     , SPMH_MAP( 43)    / 'bergamotene_a   ', 43             /
      DATA  MECH_NAM( 43)     , MECH_MAP( 43)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 43)    / 3.0            /
      DATA  SPMH_NAM( 44)     , SPMH_MAP( 44)    / 'bergamotene_b   ', 44             /
      DATA  MECH_NAM( 44)     , MECH_MAP( 44)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 44)    / 3.0            /
      DATA  SPMH_NAM( 45)     , SPMH_MAP( 45)    / 'bisabolene_a    ', 45             /
      DATA  MECH_NAM( 45)     , MECH_MAP( 45)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 45)    / 3.0            /
      DATA  SPMH_NAM( 46)     , SPMH_MAP( 46)    / 'bisabolene_b    ', 46             /
      DATA  MECH_NAM( 46)     , MECH_MAP( 46)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 46)    / 3.0            /
      DATA  SPMH_NAM( 47)     , SPMH_MAP( 47)    / 'bourbonene_b    ', 47             /
      DATA  MECH_NAM( 47)     , MECH_MAP( 47)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 47)    / 3.0            /
      DATA  SPMH_NAM( 48)     , SPMH_MAP( 48)    / 'cadinene_d      ', 48             /
      DATA  MECH_NAM( 48)     , MECH_MAP( 48)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 48)    / 3.0            /
      DATA  SPMH_NAM( 49)     , SPMH_MAP( 49)    / 'cadinene_g      ', 49             /
      DATA  MECH_NAM( 49)     , MECH_MAP( 49)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 49)    / 3.0            /
      DATA  SPMH_NAM( 50)     , SPMH_MAP( 50)    / 'cedrene_a       ', 50             /
      DATA  MECH_NAM( 50)     , MECH_MAP( 50)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 50)    / 3.0            /
      DATA  SPMH_NAM( 51)     , SPMH_MAP( 51)    / 'copaene_a       ', 51             /
      DATA  MECH_NAM( 51)     , MECH_MAP( 51)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 51)    / 3.0            /
      DATA  SPMH_NAM( 52)     , SPMH_MAP( 52)    / 'cubebene_a      ', 52             /
      DATA  MECH_NAM( 52)     , MECH_MAP( 52)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 52)    / 3.0            /
      DATA  SPMH_NAM( 53)     , SPMH_MAP( 53)    / 'cubebene_b      ', 53             /
      DATA  MECH_NAM( 53)     , MECH_MAP( 53)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 53)    / 3.0            /
      DATA  SPMH_NAM( 54)     , SPMH_MAP( 54)    / 'elemene_b       ', 54             /
      DATA  MECH_NAM( 54)     , MECH_MAP( 54)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 54)    / 3.0            /
      DATA  SPMH_NAM( 55)     , SPMH_MAP( 55)    / 'farnescene_b    ', 55             /
      DATA  MECH_NAM( 55)     , MECH_MAP( 55)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 55)    / 3.0            /
      DATA  SPMH_NAM( 56)     , SPMH_MAP( 56)    / 'germacrene_B    ', 56             /
      DATA  MECH_NAM( 56)     , MECH_MAP( 56)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 56)    / 3.0            /
      DATA  SPMH_NAM( 57)     , SPMH_MAP( 57)    / 'germacrene_D    ', 57             /
      DATA  MECH_NAM( 57)     , MECH_MAP( 57)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 57)    / 3.0            /
      DATA  SPMH_NAM( 58)     , SPMH_MAP( 58)    / 'gurjunene_b     ', 58             /
      DATA  MECH_NAM( 58)     , MECH_MAP( 58)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 58)    / 3.0            /
      DATA  SPMH_NAM( 59)     , SPMH_MAP( 59)    / 'humulene_a      ', 59             /
      DATA  MECH_NAM( 59)     , MECH_MAP( 59)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 59)    / 3.0            /
      DATA  SPMH_NAM( 60)     , SPMH_MAP( 60)    / 'humulene_g      ', 60             /
      DATA  MECH_NAM( 60)     , MECH_MAP( 60)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 60)    / 3.0            /
      DATA  SPMH_NAM( 61)     , SPMH_MAP( 61)    / 'isolongifolene  ', 61             /
      DATA  MECH_NAM( 61)     , MECH_MAP( 61)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 61)    / 3.0            /
      DATA  SPMH_NAM( 62)     , SPMH_MAP( 62)    / 'longifolene     ', 62             /
      DATA  MECH_NAM( 62)     , MECH_MAP( 62)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 62)    / 3.0            /
      DATA  SPMH_NAM( 63)     , SPMH_MAP( 63)    / 'longipinene     ', 63             /
      DATA  MECH_NAM( 63)     , MECH_MAP( 63)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 63)    / 3.0            /
      DATA  SPMH_NAM( 64)     , SPMH_MAP( 64)    / 'muurolene_a     ', 64             /
      DATA  MECH_NAM( 64)     , MECH_MAP( 64)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 64)    / 3.0            /
      DATA  SPMH_NAM( 65)     , SPMH_MAP( 65)    / 'muurolene_g     ', 65             /
      DATA  MECH_NAM( 65)     , MECH_MAP( 65)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 65)    / 3.0            /
      DATA  SPMH_NAM( 66)     , SPMH_MAP( 66)    / 'selinene_b      ', 66             /
      DATA  MECH_NAM( 66)     , MECH_MAP( 66)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 66)    / 3.0            /
      DATA  SPMH_NAM( 67)     , SPMH_MAP( 67)    / 'selinene_d      ', 67             /
      DATA  MECH_NAM( 67)     , MECH_MAP( 67)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 67)    / 3.0            /
      DATA  SPMH_NAM( 68)     , SPMH_MAP( 68)    / 'nerolidol_c     ', 68             /
      DATA  MECH_NAM( 68)     , MECH_MAP( 68)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 68)    / 3.0            /
      DATA  SPMH_NAM( 69)     , SPMH_MAP( 69)    / 'nerolidol_t     ', 69             /
      DATA  MECH_NAM( 69)     , MECH_MAP( 69)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 69)    / 3.0            /
      DATA  SPMH_NAM( 70)     , SPMH_MAP( 70)    / 'cedrol          ', 70             /
      DATA  MECH_NAM( 70)     , MECH_MAP( 70)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 70)    / 3.0            /
      DATA  SPMH_NAM( 71)     , SPMH_MAP( 71)    / 'MBO_2m3e2ol     ', 71             /
      DATA  MECH_NAM( 71)     , MECH_MAP( 71)    / 'OLET            ', 15             /
      DATA  CONV_FAC( 71)    / 1.0            /
      DATA  SPMH_NAM( 72)     , SPMH_MAP( 72)    / 'MBO_2m3e2ol     ', 71             /
      DATA  MECH_NAM( 72)     , MECH_MAP( 72)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 72)    / 3.0            /
      DATA  SPMH_NAM( 73)     , SPMH_MAP( 73)    / 'methanol        ', 72             /
      DATA  MECH_NAM( 73)     , MECH_MAP( 73)    / 'CH3OH           ', 5              /
      DATA  CONV_FAC( 73)    / 1.0            /
      DATA  SPMH_NAM( 74)     , SPMH_MAP( 74)    / 'acetone         ', 73             /
      DATA  MECH_NAM( 74)     , MECH_MAP( 74)    / 'AONE            ', 18             /
      DATA  CONV_FAC( 74)    / 1.0            /
      DATA  SPMH_NAM( 75)     , SPMH_MAP( 75)    / 'methane         ', 74             /
      DATA  MECH_NAM( 75)     , MECH_MAP( 75)    / 'CH4             ', 11             /
      DATA  CONV_FAC( 75)    / 1.0            /
      DATA  SPMH_NAM( 76)     , SPMH_MAP( 76)    / 'ammonia         ', 75             /
      DATA  MECH_NAM( 76)     , MECH_MAP( 76)    / 'NH3             ', 22             /
      DATA  CONV_FAC( 76)    / 1.0            /
      DATA  SPMH_NAM( 77)     , SPMH_MAP( 77)    / 'nitrous_OXD     ', 76             /
      DATA  MECH_NAM( 77)     , MECH_MAP( 77)    / 'NR              ', 23             /
      DATA  CONV_FAC( 77)    / 1.0            /
      DATA  SPMH_NAM( 78)     , SPMH_MAP( 78)    / 'nitric_OXD      ', 77             /
      DATA  MECH_NAM( 78)     , MECH_MAP( 78)    / 'NO              ', 2              /
      DATA  CONV_FAC( 78)    / 1.0            /
      DATA  SPMH_NAM( 79)     , SPMH_MAP( 79)    / 'acetaldehyde    ', 78             /
      DATA  MECH_NAM( 79)     , MECH_MAP( 79)    / 'ALD2            ', 7              /
      DATA  CONV_FAC( 79)    / 1.0            /
      DATA  SPMH_NAM( 80)     , SPMH_MAP( 80)    / 'ethanol         ', 79             /
      DATA  MECH_NAM( 80)     , MECH_MAP( 80)    / 'ANOL            ', 6              /
      DATA  CONV_FAC( 80)    / 1.0            /
      DATA  SPMH_NAM( 81)     , SPMH_MAP( 81)    / 'formic_acid     ', 80             /
      DATA  MECH_NAM( 81)     , MECH_MAP( 81)    / 'HCOOH           ', 9              /
      DATA  CONV_FAC( 81)    / 1.0            /
      DATA  SPMH_NAM( 82)     , SPMH_MAP( 82)    / 'formaldehyde    ', 81             /
      DATA  MECH_NAM( 82)     , MECH_MAP( 82)    / 'HCHO            ', 8              /
      DATA  CONV_FAC( 82)    / 1.0            /
      DATA  SPMH_NAM( 83)     , SPMH_MAP( 83)    / 'acetic_acid     ', 82             /
      DATA  MECH_NAM( 83)     , MECH_MAP( 83)    / 'RCOOH           ', 10             /
      DATA  CONV_FAC( 83)    / 1.0            /
      DATA  SPMH_NAM( 84)     , SPMH_MAP( 84)    / 'MBO_3m2e1ol     ', 83             /
      DATA  MECH_NAM( 84)     , MECH_MAP( 84)    / 'ALD2            ', 7              /
      DATA  CONV_FAC( 84)    / 1.0            /
      DATA  SPMH_NAM( 85)     , SPMH_MAP( 85)    / 'MBO_3m2e1ol     ', 83             /
      DATA  MECH_NAM( 85)     , MECH_MAP( 85)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 85)    / 3.0            /
      DATA  SPMH_NAM( 86)     , SPMH_MAP( 86)    / 'MBO_3m3e1ol     ', 84             /
      DATA  MECH_NAM( 86)     , MECH_MAP( 86)    / 'HCHO            ', 8              /
      DATA  CONV_FAC( 86)    / 1.0            /
      DATA  SPMH_NAM( 87)     , SPMH_MAP( 87)    / 'MBO_3m3e1ol     ', 84             /
      DATA  MECH_NAM( 87)     , MECH_MAP( 87)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 87)    / 4.0            /
      DATA  SPMH_NAM( 88)     , SPMH_MAP( 88)    / 'benzaldehyde    ', 85             /
      DATA  MECH_NAM( 88)     , MECH_MAP( 88)    / 'TOL             ', 19             /
      DATA  CONV_FAC( 88)    / 1.0            /
      DATA  SPMH_NAM( 89)     , SPMH_MAP( 89)    / 'butanone_2      ', 86             /
      DATA  MECH_NAM( 89)     , MECH_MAP( 89)    / 'AONE            ', 18             /
      DATA  CONV_FAC( 89)    / 1.0            /
      DATA  SPMH_NAM( 90)     , SPMH_MAP( 90)    / 'butanone_2      ', 86             /
      DATA  MECH_NAM( 90)     , MECH_MAP( 90)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 90)    / 1.0            /
      DATA  SPMH_NAM( 91)     , SPMH_MAP( 91)    / 'decanal         ', 87             /
      DATA  MECH_NAM( 91)     , MECH_MAP( 91)    / 'ALD2            ', 7              /
      DATA  CONV_FAC( 91)    / 1.0            /
      DATA  SPMH_NAM( 92)     , SPMH_MAP( 92)    / 'decanal         ', 87             /
      DATA  MECH_NAM( 92)     , MECH_MAP( 92)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 92)    / 8.0            /
      DATA  SPMH_NAM( 93)     , SPMH_MAP( 93)    / 'dodecene_1      ', 88             /
      DATA  MECH_NAM( 93)     , MECH_MAP( 93)    / 'OLET            ', 15             /
      DATA  CONV_FAC( 93)    / 1.0            /
      DATA  SPMH_NAM( 94)     , SPMH_MAP( 94)    / 'dodecene_1      ', 88             /
      DATA  MECH_NAM( 94)     , MECH_MAP( 94)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 94)    / 10.0           /
      DATA  SPMH_NAM( 95)     , SPMH_MAP( 95)    / 'geranyl_acetone ', 89             /
      DATA  MECH_NAM( 95)     , MECH_MAP( 95)    / 'ISOP            ', 1              /
      DATA  CONV_FAC( 95)    / 3.0            /
      DATA  SPMH_NAM( 96)     , SPMH_MAP( 96)    / 'heptanal        ', 90             /
      DATA  MECH_NAM( 96)     , MECH_MAP( 96)    / 'ALD2            ', 7              /
      DATA  CONV_FAC( 96)    / 1.0            /
      DATA  SPMH_NAM( 97)     , SPMH_MAP( 97)    / 'heptanal        ', 90             /
      DATA  MECH_NAM( 97)     , MECH_MAP( 97)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 97)    / 5.0            /
      DATA  SPMH_NAM( 98)     , SPMH_MAP( 98)    / 'heptane         ', 91             /
      DATA  MECH_NAM( 98)     , MECH_MAP( 98)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 98)    / 7.0            /
      DATA  SPMH_NAM( 99)     , SPMH_MAP( 99)    / 'hexane          ', 92             /
      DATA  MECH_NAM( 99)     , MECH_MAP( 99)    / 'PAR             ', 13             /
      DATA  CONV_FAC( 99)    / 6.0            /
      DATA  SPMH_NAM(100)     , SPMH_MAP(100)    / 'met_benzoate    ', 93             /
      DATA  MECH_NAM(100)     , MECH_MAP(100)    / 'TOL             ', 19             /
      DATA  CONV_FAC(100)    / 1.0            /
      DATA  SPMH_NAM(101)     , SPMH_MAP(101)    / 'met_heptenone   ', 94             /
      DATA  MECH_NAM(101)     , MECH_MAP(101)    / 'AONE            ', 18             /
      DATA  CONV_FAC(101)    / 1.0            /
      DATA  SPMH_NAM(102)     , SPMH_MAP(102)    / 'met_heptenone   ', 94             /
      DATA  MECH_NAM(102)     , MECH_MAP(102)    / 'PAR             ', 13             /
      DATA  CONV_FAC(102)    / 3.0            /
      DATA  SPMH_NAM(103)     , SPMH_MAP(103)    / 'met_heptenone   ', 94             /
      DATA  MECH_NAM(103)     , MECH_MAP(103)    / 'OLET            ', 15             /
      DATA  CONV_FAC(103)    / 1.0            /
      DATA  SPMH_NAM(104)     , SPMH_MAP(104)    / 'neryl_acetone   ', 95             /
      DATA  MECH_NAM(104)     , MECH_MAP(104)    / 'AONE            ', 18             /
      DATA  CONV_FAC(104)    / 1.0            /
      DATA  SPMH_NAM(105)     , SPMH_MAP(105)    / 'neryl_acetone   ', 95             /
      DATA  MECH_NAM(105)     , MECH_MAP(105)    / 'PAR             ', 13             /
      DATA  CONV_FAC(105)    / 8.0            /
      DATA  SPMH_NAM(106)     , SPMH_MAP(106)    / 'neryl_acetone   ', 95             /
      DATA  MECH_NAM(106)     , MECH_MAP(106)    / 'OLEI            ', 16             /
      DATA  CONV_FAC(106)    / 2.0            /
      DATA  SPMH_NAM(107)     , SPMH_MAP(107)    / 'nonanal         ', 96             /
      DATA  MECH_NAM(107)     , MECH_MAP(107)    / 'ALD2            ', 7              /
      DATA  CONV_FAC(107)    / 1.0            /
      DATA  SPMH_NAM(108)     , SPMH_MAP(108)    / 'nonanal         ', 96             /
      DATA  MECH_NAM(108)     , MECH_MAP(108)    / 'PAR             ', 13             /
      DATA  CONV_FAC(108)    / 7.0            /
      DATA  SPMH_NAM(109)     , SPMH_MAP(109)    / 'nonenal         ', 97             /
      DATA  MECH_NAM(109)     , MECH_MAP(109)    / 'ALD2            ', 7              /
      DATA  CONV_FAC(109)    / 1.0            /
      DATA  SPMH_NAM(110)     , SPMH_MAP(110)    / 'nonenal         ', 97             /
      DATA  MECH_NAM(110)     , MECH_MAP(110)    / 'PAR             ', 13             /
      DATA  CONV_FAC(110)    / 6.0            /
      DATA  SPMH_NAM(111)     , SPMH_MAP(111)    / 'nonenal         ', 97             /
      DATA  MECH_NAM(111)     , MECH_MAP(111)    / 'OLEI            ', 16             /
      DATA  CONV_FAC(111)    / 1.0            /
      DATA  SPMH_NAM(112)     , SPMH_MAP(112)    / 'octanal         ', 98             /
      DATA  MECH_NAM(112)     , MECH_MAP(112)    / 'ALD2            ', 7              /
      DATA  CONV_FAC(112)    / 1.0            /
      DATA  SPMH_NAM(113)     , SPMH_MAP(113)    / 'octanal         ', 98             /
      DATA  MECH_NAM(113)     , MECH_MAP(113)    / 'PAR             ', 13             /
      DATA  CONV_FAC(113)    / 6.0            /
      DATA  SPMH_NAM(114)     , SPMH_MAP(114)    / 'octanol         ', 99             /
      DATA  MECH_NAM(114)     , MECH_MAP(114)    / 'PAR             ', 13             /
      DATA  CONV_FAC(114)    / 8.0            /
      DATA  SPMH_NAM(115)     , SPMH_MAP(115)    / 'octenol_1e3ol   ', 100            /
      DATA  MECH_NAM(115)     , MECH_MAP(115)    / 'PAR             ', 13             /
      DATA  CONV_FAC(115)    / 6.0            /
      DATA  SPMH_NAM(116)     , SPMH_MAP(116)    / 'octenol_1e3ol   ', 100            /
      DATA  MECH_NAM(116)     , MECH_MAP(116)    / 'OLET            ', 15             /
      DATA  CONV_FAC(116)    / 1.0            /
      DATA  SPMH_NAM(117)     , SPMH_MAP(117)    / 'oxopentanal     ', 101            /
      DATA  MECH_NAM(117)     , MECH_MAP(117)    / 'PAR             ', 13             /
      DATA  CONV_FAC(117)    / 3.0            /
      DATA  SPMH_NAM(118)     , SPMH_MAP(118)    / 'oxopentanal     ', 101            /
      DATA  MECH_NAM(118)     , MECH_MAP(118)    / 'ALD2            ', 7              /
      DATA  CONV_FAC(118)    / 1.0            /
      DATA  SPMH_NAM(119)     , SPMH_MAP(119)    / 'pentane         ', 102            /
      DATA  MECH_NAM(119)     , MECH_MAP(119)    / 'PAR             ', 13             /
      DATA  CONV_FAC(119)    / 5.0            /
      DATA  SPMH_NAM(120)     , SPMH_MAP(120)    / 'phenyl_CCO      ', 103            /
      DATA  MECH_NAM(120)     , MECH_MAP(120)    / 'ALD2            ', 7              /
      DATA  CONV_FAC(120)    / 1.0            /
      DATA  SPMH_NAM(121)     , SPMH_MAP(121)    / 'phenyl_CCO      ', 103            /
      DATA  MECH_NAM(121)     , MECH_MAP(121)    / 'TOL             ', 19             /
      DATA  CONV_FAC(121)    / 1.0            /
      DATA  SPMH_NAM(122)     , SPMH_MAP(122)    / 'pyruvic_acid    ', 104            /
      DATA  MECH_NAM(122)     , MECH_MAP(122)    / 'HCOOH           ', 9              /
      DATA  CONV_FAC(122)    / 1.0            /
      DATA  SPMH_NAM(123)     , SPMH_MAP(123)    / 'pyruvic_acid    ', 104            /
      DATA  MECH_NAM(123)     , MECH_MAP(123)    / 'AONE            ', 18             /
      DATA  CONV_FAC(123)    / 1.0            /
      DATA  SPMH_NAM(124)     , SPMH_MAP(124)    / 'terpinyl_ACT_a  ', 105            /
      DATA  MECH_NAM(124)     , MECH_MAP(124)    / 'TERP            ', 24             /
      DATA  CONV_FAC(124)    / 1.0            /
      DATA  SPMH_NAM(125)     , SPMH_MAP(125)    / 'tetradecene_1   ', 106            /
      DATA  MECH_NAM(125)     , MECH_MAP(125)    / 'PAR             ', 13             /
      DATA  CONV_FAC(125)    / 12.0           /
      DATA  SPMH_NAM(126)     , SPMH_MAP(126)    / 'tetradecene_1   ', 106            /
      DATA  MECH_NAM(126)     , MECH_MAP(126)    / 'OLET            ', 15             /
      DATA  CONV_FAC(126)    / 1.0            /
      DATA  SPMH_NAM(127)     , SPMH_MAP(127)    / 'toluene         ', 107            /
      DATA  MECH_NAM(127)     , MECH_MAP(127)    / 'TOL             ', 19             /
      DATA  CONV_FAC(127)    / 1.0            /
      DATA  SPMH_NAM(128)     , SPMH_MAP(128)    / 'carbon_monoxide ', 108            /
      DATA  MECH_NAM(128)     , MECH_MAP(128)    / 'CO              ', 4              /
      DATA  CONV_FAC(128)    / 1.0            /
      DATA  SPMH_NAM(129)     , SPMH_MAP(129)    / 'butene          ', 109            /
      DATA  MECH_NAM(129)     , MECH_MAP(129)    / 'OLET            ', 15             /
      DATA  CONV_FAC(129)    / 1.0            /
      DATA  SPMH_NAM(130)     , SPMH_MAP(130)    / 'butene          ', 109            /
      DATA  MECH_NAM(130)     , MECH_MAP(130)    / 'PAR             ', 13             /
      DATA  CONV_FAC(130)    / 2.0            /
      DATA  SPMH_NAM(131)     , SPMH_MAP(131)    / 'ethane          ', 110            /
      DATA  MECH_NAM(131)     , MECH_MAP(131)    / 'C2H6            ', 12             /
      DATA  CONV_FAC(131)    / 1.0            /
      DATA  SPMH_NAM(132)     , SPMH_MAP(132)    / 'ethene          ', 111            /
      DATA  MECH_NAM(132)     , MECH_MAP(132)    / 'ETH             ', 14             /
      DATA  CONV_FAC(132)    / 1.0            /
      DATA  SPMH_NAM(133)     , SPMH_MAP(133)    / 'hydrogen_cyanide', 112            /
      DATA  MECH_NAM(133)     , MECH_MAP(133)    / 'NR              ', 23             /
      DATA  CONV_FAC(133)    / 1.0            /
      DATA  SPMH_NAM(134)     , SPMH_MAP(134)    / 'propane         ', 113            /
      DATA  MECH_NAM(134)     , MECH_MAP(134)    / 'PAR             ', 13             /
      DATA  CONV_FAC(134)    / 3.0            /
      DATA  SPMH_NAM(135)     , SPMH_MAP(135)    / 'propene         ', 114            /
      DATA  MECH_NAM(135)     , MECH_MAP(135)    / 'OLET            ', 15             /
      DATA  CONV_FAC(135)    / 1.0            /
      DATA  SPMH_NAM(136)     , SPMH_MAP(136)    / 'propene         ', 114            /
      DATA  MECH_NAM(136)     , MECH_MAP(136)    / 'PAR             ', 13             /
      DATA  CONV_FAC(136)    / 1.0            /
      DATA  SPMH_NAM(137)     , SPMH_MAP(137)    / 'carbon_2s       ', 115            /
      DATA  MECH_NAM(137)     , MECH_MAP(137)    / 'NR              ', 23             /
      DATA  CONV_FAC(137)    / 1.0            /
      DATA  SPMH_NAM(138)     , SPMH_MAP(138)    / 'carbonyl_s      ', 116            /
      DATA  MECH_NAM(138)     , MECH_MAP(138)    / 'NR              ', 23             /
      DATA  CONV_FAC(138)    / 1.0            /
      DATA  SPMH_NAM(139)     , SPMH_MAP(139)    / 'diallyl_2s      ', 117            /
      DATA  MECH_NAM(139)     , MECH_MAP(139)    / 'DMS             ', 21             /
      DATA  CONV_FAC(139)    / 1.0            /
      DATA  SPMH_NAM(140)     , SPMH_MAP(140)    / 'diallyl_2s      ', 117            /
      DATA  MECH_NAM(140)     , MECH_MAP(140)    / 'PAR             ', 13             /
      DATA  CONV_FAC(140)    / 2.0            /
      DATA  SPMH_NAM(141)     , SPMH_MAP(141)    / 'diallyl_2s      ', 117            /
      DATA  MECH_NAM(141)     , MECH_MAP(141)    / 'OLET            ', 15             /
      DATA  CONV_FAC(141)    / 2.0            /
      DATA  SPMH_NAM(142)     , SPMH_MAP(142)    / '2met_2s         ', 118            /
      DATA  MECH_NAM(142)     , MECH_MAP(142)    / 'DMS             ', 21             /
      DATA  CONV_FAC(142)    / 1.0            /
      DATA  SPMH_NAM(143)     , SPMH_MAP(143)    / '2met_s          ', 119            /
      DATA  MECH_NAM(143)     , MECH_MAP(143)    / 'DMS             ', 21             /
      DATA  CONV_FAC(143)    / 1.0            /
      DATA  SPMH_NAM(144)     , SPMH_MAP(144)    / 'met_chloride    ', 120            /
      DATA  MECH_NAM(144)     , MECH_MAP(144)    / 'NR              ', 23             /
      DATA  CONV_FAC(144)    / 1.0            /
      DATA  SPMH_NAM(145)     , SPMH_MAP(145)    / 'met_bromide     ', 121            /
      DATA  MECH_NAM(145)     , MECH_MAP(145)    / 'NR              ', 23             /
      DATA  CONV_FAC(145)    / 1.0            /
      DATA  SPMH_NAM(146)     , SPMH_MAP(146)    / 'met_iodide      ', 122            /
      DATA  MECH_NAM(146)     , MECH_MAP(146)    / 'NR              ', 23             /
      DATA  CONV_FAC(146)    / 1.0            /
      DATA  SPMH_NAM(147)     , SPMH_MAP(147)    / 'hydrogen_s      ', 123            /
      DATA  MECH_NAM(147)     , MECH_MAP(147)    / 'NR              ', 23             /
      DATA  CONV_FAC(147)    / 1.0            /
      DATA  SPMH_NAM(148)     , SPMH_MAP(148)    / 'met_mercaptan   ', 124            /
      DATA  MECH_NAM(148)     , MECH_MAP(148)    / 'PAR             ', 13             /
      DATA  CONV_FAC(148)    / 1.0            /
      DATA  SPMH_NAM(149)     , SPMH_MAP(149)    / 'met_propenyl_2s ', 125            /
      DATA  MECH_NAM(149)     , MECH_MAP(149)    / 'DMS             ', 21             /
      DATA  CONV_FAC(149)    / 1.0            /
      DATA  SPMH_NAM(150)     , SPMH_MAP(150)    / 'met_propenyl_2s ', 125            /
      DATA  MECH_NAM(150)     , MECH_MAP(150)    / 'OLEI            ', 16             /
      DATA  CONV_FAC(150)    / 1.0            /
      DATA  SPMH_NAM(151)     , SPMH_MAP(151)    / 'PPPP_2s         ', 126            /
      DATA  MECH_NAM(151)     , MECH_MAP(151)    / 'DMS             ', 21             /
      DATA  CONV_FAC(151)    / 1.0            /
      DATA  SPMH_NAM(152)     , SPMH_MAP(152)    / 'PPPP_2s         ', 126            /
      DATA  MECH_NAM(152)     , MECH_MAP(152)    / 'PAR             ', 13             /
      DATA  CONV_FAC(152)    / 2.0            /
      DATA  SPMH_NAM(153)     , SPMH_MAP(153)    / 'PPPP_2s         ', 126            /
      DATA  MECH_NAM(153)     , MECH_MAP(153)    / 'OLEI            ', 16             /
      DATA  CONV_FAC(153)    / 1.0            /
      DATA  SPMH_NAM(154)     , SPMH_MAP(154)    / '2met_nonatriene ', 127            /
      DATA  MECH_NAM(154)     , MECH_MAP(154)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(154)    / 2.0            /
      DATA  SPMH_NAM(155)     , SPMH_MAP(155)    / 'met_salicylate  ', 128            /
      DATA  MECH_NAM(155)     , MECH_MAP(155)    / 'TOL             ', 19             /
      DATA  CONV_FAC(155)    / 1.0            /
      DATA  SPMH_NAM(156)     , SPMH_MAP(156)    / 'indole          ', 129            /
      DATA  MECH_NAM(156)     , MECH_MAP(156)    / 'TOL             ', 19             /
      DATA  CONV_FAC(156)    / 1.0            /
      DATA  SPMH_NAM(157)     , SPMH_MAP(157)    / 'jasmone         ', 130            /
      DATA  MECH_NAM(157)     , MECH_MAP(157)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(157)    / 2.0            /
      DATA  SPMH_NAM(158)     , SPMH_MAP(158)    / 'met_jasmonate   ', 131            /
      DATA  MECH_NAM(158)     , MECH_MAP(158)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(158)    / 3.0            /
      DATA  SPMH_NAM(159)     , SPMH_MAP(159)    / '3met_3DCTT      ', 132            /
      DATA  MECH_NAM(159)     , MECH_MAP(159)    / 'ISOP            ', 1              /
      DATA  CONV_FAC(159)    / 3.0            /
      DATA  SPMH_NAM(160)     , SPMH_MAP(160)    / 'hexanal         ', 133            /
      DATA  MECH_NAM(160)     , MECH_MAP(160)    / 'ALD2            ', 7              /
      DATA  CONV_FAC(160)    / 1.0            /
      DATA  SPMH_NAM(161)     , SPMH_MAP(161)    / 'hexanal         ', 133            /
      DATA  MECH_NAM(161)     , MECH_MAP(161)    / 'PAR             ', 13             /
      DATA  CONV_FAC(161)    / 4.0            /
      DATA  SPMH_NAM(162)     , SPMH_MAP(162)    / 'hexanol_1       ', 134            /
      DATA  MECH_NAM(162)     , MECH_MAP(162)    / 'PAR             ', 13             /
      DATA  CONV_FAC(162)    / 6.0            /
      DATA  SPMH_NAM(163)     , SPMH_MAP(163)    / 'hexenal_c3      ', 135            /
      DATA  MECH_NAM(163)     , MECH_MAP(163)    / 'ALD2            ', 7              /
      DATA  CONV_FAC(163)    / 1.0            /
      DATA  SPMH_NAM(164)     , SPMH_MAP(164)    / 'hexenal_c3      ', 135            /
      DATA  MECH_NAM(164)     , MECH_MAP(164)    / 'PAR             ', 13             /
      DATA  CONV_FAC(164)    / 3.0            /
      DATA  SPMH_NAM(165)     , SPMH_MAP(165)    / 'hexenal_c3      ', 135            /
      DATA  MECH_NAM(165)     , MECH_MAP(165)    / 'OLEI            ', 16             /
      DATA  CONV_FAC(165)    / 1.0            /
      DATA  SPMH_NAM(166)     , SPMH_MAP(166)    / 'hexenal_t2      ', 136            /
      DATA  MECH_NAM(166)     , MECH_MAP(166)    / 'ALD2            ', 7              /
      DATA  CONV_FAC(166)    / 1.0            /
      DATA  SPMH_NAM(167)     , SPMH_MAP(167)    / 'hexenal_t2      ', 136            /
      DATA  MECH_NAM(167)     , MECH_MAP(167)    / 'PAR             ', 13             /
      DATA  CONV_FAC(167)    / 6.0            /
      DATA  SPMH_NAM(168)     , SPMH_MAP(168)    / 'hexenal_t2      ', 136            /
      DATA  MECH_NAM(168)     , MECH_MAP(168)    / 'OLEI            ', 16             /
      DATA  CONV_FAC(168)    / 1.0            /
      DATA  SPMH_NAM(169)     , SPMH_MAP(169)    / 'hexenol_c3      ', 137            /
      DATA  MECH_NAM(169)     , MECH_MAP(169)    / 'PAR             ', 13             /
      DATA  CONV_FAC(169)    / 5.0            /
      DATA  SPMH_NAM(170)     , SPMH_MAP(170)    / 'hexenol_c3      ', 137            /
      DATA  MECH_NAM(170)     , MECH_MAP(170)    / 'OLEI            ', 16             /
      DATA  CONV_FAC(170)    / 1.0            /
      DATA  SPMH_NAM(171)     , SPMH_MAP(171)    / 'hexenyl_ACT_c3  ', 138            /
      DATA  MECH_NAM(171)     , MECH_MAP(171)    / 'AONE            ', 18             /
      DATA  CONV_FAC(171)    / 1.0            /
      DATA  SPMH_NAM(172)     , SPMH_MAP(172)    / 'hexenyl_ACT_c3  ', 138            /
      DATA  MECH_NAM(172)     , MECH_MAP(172)    / 'PAR             ', 13             /
      DATA  CONV_FAC(172)    / 3.0            /
      DATA  SPMH_NAM(173)     , SPMH_MAP(173)    / 'hexenyl_ACT_c3  ', 138            /
      DATA  MECH_NAM(173)     , MECH_MAP(173)    / 'OLEI            ', 16             /
      DATA  CONV_FAC(173)    / 1.0            /

