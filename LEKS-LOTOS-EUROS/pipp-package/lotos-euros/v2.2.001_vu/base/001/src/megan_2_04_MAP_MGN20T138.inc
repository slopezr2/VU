C=======================================================================
C  MAP_MGN20T138.EXT
C  This include file contains conversion table for MEGAN species to
C  134 species
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

      INTEGER        N_SMAP_SPC
      PARAMETER     (N_SMAP_SPC = 138)        ! Number of map species

      CHARACTER*16   SPCA_NAM( N_SMAP_SPC )   ! speciated species name
      INTEGER        SPCA_MAP( N_SMAP_SPC )   ! speciated species name
                                              ! mapped to SPCAT_SPC.EXT
      CHARACTER*16   MG20_NAM( N_SMAP_SPC )   ! MEGAN species
      INTEGER        MG20_MAP( N_SMAP_SPC )   ! MEGAN species mapped to
                                              ! MGN_SPC.EXT

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C _a  = alpha     , _b  = beta      , _c  = cis       , _al = allo      ,
C _g  = gamma     , _d  = delta     , _t  = trans     , _m  = methyl    ,
C _p  = para      , _o  = ortho     , _e  = ene       , _ol = ol        ,
C met = methyl    , 2met= dimethyl  , MBO = methylbutenol               ,
C 2s  = disulfide , s   = sulfide   , OXD = oxide     , ACT = acetate   ,
C PPPP= propenylpropyl              , DCTT= decatetraene                ,
C CCO = acetaldehyde
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      DATA  SPCA_NAM(  1)     , SPCA_MAP(  1), MG20_NAM(  1)     , MG20_MAP(  1)
     &    / 'isoprene        ', 1            , 'ISOP            ', 1             /
      DATA  SPCA_NAM(  2)     , SPCA_MAP(  2), MG20_NAM(  2)     , MG20_MAP(  2)
     &    / 'myrcene         ', 2            , 'MYRC            ', 2             /
      DATA  SPCA_NAM(  3)     , SPCA_MAP(  3), MG20_NAM(  3)     , MG20_MAP(  3)
     &    / 'sabinene        ', 3            , 'SABI            ', 3             /
      DATA  SPCA_NAM(  4)     , SPCA_MAP(  4), MG20_NAM(  4)     , MG20_MAP(  4)
     &    / 'limonene        ', 4            , 'LIMO            ', 4             /
      DATA  SPCA_NAM(  5)     , SPCA_MAP(  5), MG20_NAM(  5)     , MG20_MAP(  5)
     &    / 'carene_3        ', 5            , '3CAR            ', 5             /
      DATA  SPCA_NAM(  6)     , SPCA_MAP(  6), MG20_NAM(  6)     , MG20_MAP(  6)
     &    / 'ocimene_t_b     ', 6            , 'OCIM            ', 6             /
      DATA  SPCA_NAM(  7)     , SPCA_MAP(  7), MG20_NAM(  7)     , MG20_MAP(  7)
     &    / 'pinene_b        ', 7            , 'BPIN            ', 7             /
      DATA  SPCA_NAM(  8)     , SPCA_MAP(  8), MG20_NAM(  8)     , MG20_MAP(  8)
     &    / 'pinene_a        ', 8            , 'APIN            ', 8             /
      DATA  SPCA_NAM(  9)     , SPCA_MAP(  9), MG20_NAM(  9)     , MG20_MAP(  9)
     &    / '2met_styrene    ', 9            , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 10)     , SPCA_MAP( 10), MG20_NAM( 10)     , MG20_MAP( 10)
     &    / 'cymene_p        ', 10           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 11)     , SPCA_MAP( 11), MG20_NAM( 11)     , MG20_MAP( 11)
     &    / 'cymene_o        ', 11           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 12)     , SPCA_MAP( 12), MG20_NAM( 12)     , MG20_MAP( 12)
     &    / 'phellandrene_a  ', 12           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 13)     , SPCA_MAP( 13), MG20_NAM( 13)     , MG20_MAP( 13)
     &    / 'thujene_a       ', 13           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 14)     , SPCA_MAP( 14), MG20_NAM( 14)     , MG20_MAP( 14)
     &    / 'terpinene_a     ', 14           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 15)     , SPCA_MAP( 15), MG20_NAM( 15)     , MG20_MAP( 15)
     &    / 'terpinene_g     ', 15           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 16)     , SPCA_MAP( 16), MG20_NAM( 16)     , MG20_MAP( 16)
     &    / 'terpinolene     ', 16           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 17)     , SPCA_MAP( 17), MG20_NAM( 17)     , MG20_MAP( 17)
     &    / 'phellandrene_b  ', 17           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 18)     , SPCA_MAP( 18), MG20_NAM( 18)     , MG20_MAP( 18)
     &    / 'camphene        ', 18           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 19)     , SPCA_MAP( 19), MG20_NAM( 19)     , MG20_MAP( 19)
     &    / 'bornene         ', 19           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 20)     , SPCA_MAP( 20), MG20_NAM( 20)     , MG20_MAP( 20)
     &    / 'fenchene_a      ', 20           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 21)     , SPCA_MAP( 21), MG20_NAM( 21)     , MG20_MAP( 21)
     &    / 'ocimene_al      ', 21           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 22)     , SPCA_MAP( 22), MG20_NAM( 22)     , MG20_MAP( 22)
     &    / 'ocimene_c_b     ', 22           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 23)     , SPCA_MAP( 23), MG20_NAM( 23)     , MG20_MAP( 23)
     &    / 'tricyclene      ', 23           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 24)     , SPCA_MAP( 24), MG20_NAM( 24)     , MG20_MAP( 24)
     &    / 'estragole       ', 24           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 25)     , SPCA_MAP( 25), MG20_NAM( 25)     , MG20_MAP( 25)
     &    / 'camphor         ', 25           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 26)     , SPCA_MAP( 26), MG20_NAM( 26)     , MG20_MAP( 26)
     &    / 'fenchone        ', 26           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 27)     , SPCA_MAP( 27), MG20_NAM( 27)     , MG20_MAP( 27)
     &    / 'piperitone      ', 27           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 28)     , SPCA_MAP( 28), MG20_NAM( 28)     , MG20_MAP( 28)
     &    / 'thujone_a       ', 28           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 29)     , SPCA_MAP( 29), MG20_NAM( 29)     , MG20_MAP( 29)
     &    / 'thujone_b       ', 29           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 30)     , SPCA_MAP( 30), MG20_NAM( 30)     , MG20_MAP( 30)
     &    / 'cineole_1_8     ', 30           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 31)     , SPCA_MAP( 31), MG20_NAM( 31)     , MG20_MAP( 31)
     &    / 'borneol         ', 31           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 32)     , SPCA_MAP( 32), MG20_NAM( 32)     , MG20_MAP( 32)
     &    / 'linalool        ', 32           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 33)     , SPCA_MAP( 33), MG20_NAM( 33)     , MG20_MAP( 33)
     &    / 'terpineol_4     ', 33           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 34)     , SPCA_MAP( 34), MG20_NAM( 34)     , MG20_MAP( 34)
     &    / 'terpineol_a     ', 34           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 35)     , SPCA_MAP( 35), MG20_NAM( 35)     , MG20_MAP( 35)
     &    / 'linalool_OXD_c  ', 35           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 36)     , SPCA_MAP( 36), MG20_NAM( 36)     , MG20_MAP( 36)
     &    / 'linalool_OXD_t  ', 36           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 37)     , SPCA_MAP( 37), MG20_NAM( 37)     , MG20_MAP( 37)
     &    / 'ionone_b        ', 37           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 38)     , SPCA_MAP( 38), MG20_NAM( 38)     , MG20_MAP( 38)
     &    / 'bornyl_ACT      ', 38           , 'OMTP            ', 9             /
      DATA  SPCA_NAM( 39)     , SPCA_MAP( 39), MG20_NAM( 39)     , MG20_MAP( 39)
     &    / 'farnescene_a    ', 39           , 'FARN            ', 10            /
      DATA  SPCA_NAM( 40)     , SPCA_MAP( 40), MG20_NAM( 40)     , MG20_MAP( 40)
     &    / 'caryophyllene_b ', 40           , 'BCAR            ', 11            /
      DATA  SPCA_NAM( 41)     , SPCA_MAP( 41), MG20_NAM( 41)     , MG20_MAP( 41)
     &    / 'acoradiene      ', 41           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 42)     , SPCA_MAP( 42), MG20_NAM( 42)     , MG20_MAP( 42)
     &    / 'aromadendrene   ', 42           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 43)     , SPCA_MAP( 43), MG20_NAM( 43)     , MG20_MAP( 43)
     &    / 'bergamotene_a   ', 43           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 44)     , SPCA_MAP( 44), MG20_NAM( 44)     , MG20_MAP( 44)
     &    / 'bergamotene_b   ', 44           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 45)     , SPCA_MAP( 45), MG20_NAM( 45)     , MG20_MAP( 45)
     &    / 'bisabolene_a    ', 45           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 46)     , SPCA_MAP( 46), MG20_NAM( 46)     , MG20_MAP( 46)
     &    / 'bisabolene_b    ', 46           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 47)     , SPCA_MAP( 47), MG20_NAM( 47)     , MG20_MAP( 47)
     &    / 'bourbonene_b    ', 47           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 48)     , SPCA_MAP( 48), MG20_NAM( 48)     , MG20_MAP( 48)
     &    / 'cadinene_d      ', 48           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 49)     , SPCA_MAP( 49), MG20_NAM( 49)     , MG20_MAP( 49)
     &    / 'cadinene_g      ', 49           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 50)     , SPCA_MAP( 50), MG20_NAM( 50)     , MG20_MAP( 50)
     &    / 'cedrene_a       ', 50           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 51)     , SPCA_MAP( 51), MG20_NAM( 51)     , MG20_MAP( 51)
     &    / 'copaene_a       ', 51           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 52)     , SPCA_MAP( 52), MG20_NAM( 52)     , MG20_MAP( 52)
     &    / 'cubebene_a      ', 52           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 53)     , SPCA_MAP( 53), MG20_NAM( 53)     , MG20_MAP( 53)
     &    / 'cubebene_b      ', 53           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 54)     , SPCA_MAP( 54), MG20_NAM( 54)     , MG20_MAP( 54)
     &    / 'elemene_b       ', 54           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 55)     , SPCA_MAP( 55), MG20_NAM( 55)     , MG20_MAP( 55)
     &    / 'farnescene_b    ', 55           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 56)     , SPCA_MAP( 56), MG20_NAM( 56)     , MG20_MAP( 56)
     &    / 'germacrene_B    ', 56           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 57)     , SPCA_MAP( 57), MG20_NAM( 57)     , MG20_MAP( 57)
     &    / 'germacrene_D    ', 57           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 58)     , SPCA_MAP( 58), MG20_NAM( 58)     , MG20_MAP( 58)
     &    / 'gurjunene_b     ', 58           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 59)     , SPCA_MAP( 59), MG20_NAM( 59)     , MG20_MAP( 59)
     &    / 'humulene_a      ', 59           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 60)     , SPCA_MAP( 60), MG20_NAM( 60)     , MG20_MAP( 60)
     &    / 'humulene_g      ', 60           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 61)     , SPCA_MAP( 61), MG20_NAM( 61)     , MG20_MAP( 61)
     &    / 'isolongifolene  ', 61           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 62)     , SPCA_MAP( 62), MG20_NAM( 62)     , MG20_MAP( 62)
     &    / 'longifolene     ', 62           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 63)     , SPCA_MAP( 63), MG20_NAM( 63)     , MG20_MAP( 63)
     &    / 'longipinene     ', 63           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 64)     , SPCA_MAP( 64), MG20_NAM( 64)     , MG20_MAP( 64)
     &    / 'muurolene_a     ', 64           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 65)     , SPCA_MAP( 65), MG20_NAM( 65)     , MG20_MAP( 65)
     &    / 'muurolene_g     ', 65           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 66)     , SPCA_MAP( 66), MG20_NAM( 66)     , MG20_MAP( 66)
     &    / 'selinene_b      ', 66           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 67)     , SPCA_MAP( 67), MG20_NAM( 67)     , MG20_MAP( 67)
     &    / 'selinene_d      ', 67           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 68)     , SPCA_MAP( 68), MG20_NAM( 68)     , MG20_MAP( 68)
     &    / 'nerolidol_c     ', 68           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 69)     , SPCA_MAP( 69), MG20_NAM( 69)     , MG20_MAP( 69)
     &    / 'nerolidol_t     ', 69           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 70)     , SPCA_MAP( 70), MG20_NAM( 70)     , MG20_MAP( 70)
     &    / 'cedrol          ', 70           , 'OSQT            ', 12            /
      DATA  SPCA_NAM( 71)     , SPCA_MAP( 71), MG20_NAM( 71)     , MG20_MAP( 71)
     &    / 'MBO_2m3e2ol     ', 71           , 'MBO             ', 13            /
      DATA  SPCA_NAM( 72)     , SPCA_MAP( 72), MG20_NAM( 72)     , MG20_MAP( 72)
     &    / 'methanol        ', 72           , 'MEOH            ', 14            /
      DATA  SPCA_NAM( 73)     , SPCA_MAP( 73), MG20_NAM( 73)     , MG20_MAP( 73)
     &    / 'acetone         ', 73           , 'ACTO            ', 15            /
      DATA  SPCA_NAM( 74)     , SPCA_MAP( 74), MG20_NAM( 74)     , MG20_MAP( 74)
     &    / 'methane         ', 74           , 'CH4             ', 16            /
      DATA  SPCA_NAM( 75)     , SPCA_MAP( 75), MG20_NAM( 75)     , MG20_MAP( 75)
     &    / 'ammonia         ', 75           , 'NO              ', 17            /
      DATA  SPCA_NAM( 76)     , SPCA_MAP( 76), MG20_NAM( 76)     , MG20_MAP( 76)
     &    / 'nitrous_OXD     ', 76           , 'NO              ', 17            /
      DATA  SPCA_NAM( 77)     , SPCA_MAP( 77), MG20_NAM( 77)     , MG20_MAP( 77)
     &    / 'nitric_OXD      ', 77           , 'NO              ', 17            /
      DATA  SPCA_NAM( 78)     , SPCA_MAP( 78), MG20_NAM( 78)     , MG20_MAP( 78)
     &    / 'acetaldehyde    ', 78           , 'ACTA            ', 18            /
      DATA  SPCA_NAM( 79)     , SPCA_MAP( 79), MG20_NAM( 79)     , MG20_MAP( 79)
     &    / 'ethanol         ', 79           , 'ACTA            ', 18            /
      DATA  SPCA_NAM( 80)     , SPCA_MAP( 80), MG20_NAM( 80)     , MG20_MAP( 80)
     &    / 'formic_acid     ', 80           , 'FORM            ', 19            /
      DATA  SPCA_NAM( 81)     , SPCA_MAP( 81), MG20_NAM( 81)     , MG20_MAP( 81)
     &    / 'formaldehyde    ', 81           , 'FORM            ', 19            /
      DATA  SPCA_NAM( 82)     , SPCA_MAP( 82), MG20_NAM( 82)     , MG20_MAP( 82)
     &    / 'acetic_acid     ', 82           , 'FORM            ', 19            /
      DATA  SPCA_NAM( 83)     , SPCA_MAP( 83), MG20_NAM( 83)     , MG20_MAP( 83)
     &    / 'MBO_3m2e1ol     ', 83           , 'CO              ', 20            /
      DATA  SPCA_NAM( 84)     , SPCA_MAP( 84), MG20_NAM( 84)     , MG20_MAP( 84)
     &    / 'MBO_3m3e1ol     ', 84           , 'CO              ', 20            /
      DATA  SPCA_NAM( 85)     , SPCA_MAP( 85), MG20_NAM( 85)     , MG20_MAP( 85)
     &    / 'benzaldehyde    ', 85           , 'CO              ', 20            /
      DATA  SPCA_NAM( 86)     , SPCA_MAP( 86), MG20_NAM( 86)     , MG20_MAP( 86)
     &    / 'butanone_2      ', 86           , 'CO              ', 20            /
      DATA  SPCA_NAM( 87)     , SPCA_MAP( 87), MG20_NAM( 87)     , MG20_MAP( 87)
     &    / 'decanal         ', 87           , 'CO              ', 20            /
      DATA  SPCA_NAM( 88)     , SPCA_MAP( 88), MG20_NAM( 88)     , MG20_MAP( 88)
     &    / 'dodecene_1      ', 88           , 'CO              ', 20            /
      DATA  SPCA_NAM( 89)     , SPCA_MAP( 89), MG20_NAM( 89)     , MG20_MAP( 89)
     &    / 'geranyl_acetone ', 89           , 'CO              ', 20            /
      DATA  SPCA_NAM( 90)     , SPCA_MAP( 90), MG20_NAM( 90)     , MG20_MAP( 90)
     &    / 'heptanal        ', 90           , 'CO              ', 20            /
      DATA  SPCA_NAM( 91)     , SPCA_MAP( 91), MG20_NAM( 91)     , MG20_MAP( 91)
     &    / 'heptane         ', 91           , 'CO              ', 20            /
      DATA  SPCA_NAM( 92)     , SPCA_MAP( 92), MG20_NAM( 92)     , MG20_MAP( 92)
     &    / 'hexane          ', 92           , 'CO              ', 20            /
      DATA  SPCA_NAM( 93)     , SPCA_MAP( 93), MG20_NAM( 93)     , MG20_MAP( 93)
     &    / 'met_benzoate    ', 93           , 'CO              ', 20            /
      DATA  SPCA_NAM( 94)     , SPCA_MAP( 94), MG20_NAM( 94)     , MG20_MAP( 94)
     &    / 'met_heptenone   ', 94           , 'CO              ', 20            /
      DATA  SPCA_NAM( 95)     , SPCA_MAP( 95), MG20_NAM( 95)     , MG20_MAP( 95)
     &    / 'neryl_acetone   ', 95           , 'CO              ', 20            /
      DATA  SPCA_NAM( 96)     , SPCA_MAP( 96), MG20_NAM( 96)     , MG20_MAP( 96)
     &    / 'nonanal         ', 96           , 'CO              ', 20            /
      DATA  SPCA_NAM( 97)     , SPCA_MAP( 97), MG20_NAM( 97)     , MG20_MAP( 97)
     &    / 'nonenal         ', 97           , 'CO              ', 20            /
      DATA  SPCA_NAM( 98)     , SPCA_MAP( 98), MG20_NAM( 98)     , MG20_MAP( 98)
     &    / 'octanal         ', 98           , 'CO              ', 20            /
      DATA  SPCA_NAM( 99)     , SPCA_MAP( 99), MG20_NAM( 99)     , MG20_MAP( 99)
     &    / 'octanol         ', 99           , 'CO              ', 20            /
      DATA  SPCA_NAM(100)     , SPCA_MAP(100), MG20_NAM(100)     , MG20_MAP(100)
     &    / 'octenol_1e3ol   ', 100          , 'CO              ', 20            /
      DATA  SPCA_NAM(101)     , SPCA_MAP(101), MG20_NAM(101)     , MG20_MAP(101)
     &    / 'oxopentanal     ', 101          , 'CO              ', 20            /
      DATA  SPCA_NAM(102)     , SPCA_MAP(102), MG20_NAM(102)     , MG20_MAP(102)
     &    / 'pentane         ', 102          , 'CO              ', 20            /
      DATA  SPCA_NAM(103)     , SPCA_MAP(103), MG20_NAM(103)     , MG20_MAP(103)
     &    / 'phenyl_CCO      ', 103          , 'CO              ', 20            /
      DATA  SPCA_NAM(104)     , SPCA_MAP(104), MG20_NAM(104)     , MG20_MAP(104)
     &    / 'pyruvic_acid    ', 104          , 'CO              ', 20            /
      DATA  SPCA_NAM(105)     , SPCA_MAP(105), MG20_NAM(105)     , MG20_MAP(105)
     &    / 'terpinyl_ACT_a  ', 105          , 'CO              ', 20            /
      DATA  SPCA_NAM(106)     , SPCA_MAP(106), MG20_NAM(106)     , MG20_MAP(106)
     &    / 'tetradecene_1   ', 106          , 'CO              ', 20            /
      DATA  SPCA_NAM(107)     , SPCA_MAP(107), MG20_NAM(107)     , MG20_MAP(107)
     &    / 'toluene         ', 107          , 'CO              ', 20            /
      DATA  SPCA_NAM(108)     , SPCA_MAP(108), MG20_NAM(108)     , MG20_MAP(108)
     &    / 'carbon_monoxide ', 108          , 'CO              ', 20            /
      DATA  SPCA_NAM(109)     , SPCA_MAP(109), MG20_NAM(109)     , MG20_MAP(109)
     &    / 'butene          ', 109          , 'CO              ', 20            /
      DATA  SPCA_NAM(110)     , SPCA_MAP(110), MG20_NAM(110)     , MG20_MAP(110)
     &    / 'ethane          ', 110          , 'CO              ', 20            /
      DATA  SPCA_NAM(111)     , SPCA_MAP(111), MG20_NAM(111)     , MG20_MAP(111)
     &    / 'ethene          ', 111          , 'CO              ', 20            /
      DATA  SPCA_NAM(112)     , SPCA_MAP(112), MG20_NAM(112)     , MG20_MAP(112)
     &    / 'hydrogen_cyanide', 112          , 'CO              ', 20            /
      DATA  SPCA_NAM(113)     , SPCA_MAP(113), MG20_NAM(113)     , MG20_MAP(113)
     &    / 'propane         ', 113          , 'CO              ', 20            /
      DATA  SPCA_NAM(114)     , SPCA_MAP(114), MG20_NAM(114)     , MG20_MAP(114)
     &    / 'propene         ', 114          , 'CO              ', 20            /
      DATA  SPCA_NAM(115)     , SPCA_MAP(115), MG20_NAM(115)     , MG20_MAP(115)
     &    / 'carbon_2s       ', 115          , 'CO              ', 20            /
      DATA  SPCA_NAM(116)     , SPCA_MAP(116), MG20_NAM(116)     , MG20_MAP(116)
     &    / 'carbonyl_s      ', 116          , 'CO              ', 20            /
      DATA  SPCA_NAM(117)     , SPCA_MAP(117), MG20_NAM(117)     , MG20_MAP(117)
     &    / 'diallyl_2s      ', 117          , 'CO              ', 20            /
      DATA  SPCA_NAM(118)     , SPCA_MAP(118), MG20_NAM(118)     , MG20_MAP(118)
     &    / '2met_2s         ', 118          , 'CO              ', 20            /
      DATA  SPCA_NAM(119)     , SPCA_MAP(119), MG20_NAM(119)     , MG20_MAP(119)
     &    / '2met_s          ', 119          , 'CO              ', 20            /
      DATA  SPCA_NAM(120)     , SPCA_MAP(120), MG20_NAM(120)     , MG20_MAP(120)
     &    / 'met_chloride    ', 120          , 'CO              ', 20            /
      DATA  SPCA_NAM(121)     , SPCA_MAP(121), MG20_NAM(121)     , MG20_MAP(121)
     &    / 'met_bromide     ', 121          , 'CO              ', 20            /
      DATA  SPCA_NAM(122)     , SPCA_MAP(122), MG20_NAM(122)     , MG20_MAP(122)
     &    / 'met_iodide      ', 122          , 'CO              ', 20            /
      DATA  SPCA_NAM(123)     , SPCA_MAP(123), MG20_NAM(123)     , MG20_MAP(123)
     &    / 'hydrogen_s      ', 123          , 'CO              ', 20            /
      DATA  SPCA_NAM(124)     , SPCA_MAP(124), MG20_NAM(124)     , MG20_MAP(124)
     &    / 'met_mercaptan   ', 124          , 'CO              ', 20            /
      DATA  SPCA_NAM(125)     , SPCA_MAP(125), MG20_NAM(125)     , MG20_MAP(125)
     &    / 'met_propenyl_2s ', 125          , 'CO              ', 20            /
      DATA  SPCA_NAM(126)     , SPCA_MAP(126), MG20_NAM(126)     , MG20_MAP(126)
     &    / 'PPPP_2s         ', 126          , 'CO              ', 20            /
      DATA  SPCA_NAM(127)     , SPCA_MAP(127), MG20_NAM(127)     , MG20_MAP(127)
     &    / '2met_nonatriene ', 127          , 'CO              ', 20            /
      DATA  SPCA_NAM(128)     , SPCA_MAP(128), MG20_NAM(128)     , MG20_MAP(128)
     &    / 'met_salicylate  ', 128          , 'CO              ', 20            /
      DATA  SPCA_NAM(129)     , SPCA_MAP(129), MG20_NAM(129)     , MG20_MAP(129)
     &    / 'indole          ', 129          , 'CO              ', 20            /
      DATA  SPCA_NAM(130)     , SPCA_MAP(130), MG20_NAM(130)     , MG20_MAP(130)
     &    / 'jasmone         ', 130          , 'CO              ', 20            /
      DATA  SPCA_NAM(131)     , SPCA_MAP(131), MG20_NAM(131)     , MG20_MAP(131)
     &    / 'met_jasmonate   ', 131          , 'CO              ', 20            /
      DATA  SPCA_NAM(132)     , SPCA_MAP(132), MG20_NAM(132)     , MG20_MAP(132)
     &    / '3met_3DCTT      ', 132          , 'CO              ', 20            /
      DATA  SPCA_NAM(133)     , SPCA_MAP(133), MG20_NAM(133)     , MG20_MAP(133)
     &    / 'hexanal         ', 133          , 'CO              ', 20            /
      DATA  SPCA_NAM(134)     , SPCA_MAP(134), MG20_NAM(134)     , MG20_MAP(134)
     &    / 'hexanol_1       ', 134          , 'CO              ', 20            /
      DATA  SPCA_NAM(135)     , SPCA_MAP(135), MG20_NAM(135)     , MG20_MAP(135)
     &    / 'hexenal_c3      ', 135          , 'CO              ', 20            /
      DATA  SPCA_NAM(136)     , SPCA_MAP(136), MG20_NAM(136)     , MG20_MAP(136)
     &    / 'hexenal_t2      ', 136          , 'CO              ', 20            /
      DATA  SPCA_NAM(137)     , SPCA_MAP(137), MG20_NAM(137)     , MG20_MAP(137)
     &    / 'hexenol_c3      ', 137          , 'CO              ', 20            /
      DATA  SPCA_NAM(138)     , SPCA_MAP(138), MG20_NAM(138)     , MG20_MAP(138)
     &    / 'hexenyl_ACT_c3  ', 138          , 'CO              ', 20            /

