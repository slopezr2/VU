C=======================================================================
C  SPC_SPCAT.EXT
C  This include file contains MEGAN speciated species and their MW.
C
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

      INTEGER        N_SPCA_SPC
      PARAMETER     (N_SPCA_SPC = 138)        ! Number of speciated species

      CHARACTER*16   SPCA_SPC( N_SPCA_SPC )   ! speciated species name
      REAL           SPCA_MWT( N_SPCA_SPC )   ! Mechanism species molecular weight

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C _a  = alpha     , _b  = beta      , _c  = cis       , _al = allo      ,
C _g  = gamma     , _d  = delta     , _t  = trans     , _m  = methyl    ,
C _p  = para      , _o  = ortho     , _e  = ene       , _ol = ol        ,
C met = methyl    , 2met= dimethyl  , MBO = methylbutenol               ,
C 2s  = disulfide , s   = sulfide   , OXD = oxide     , ACT = acetate   ,
C PPPP= propenylpropyl              , DCTT= decatetraene                ,
C CCO = acetaldehyde
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Isoprene
      DATA     SPCA_SPC(  1), SPCA_MWT(  1) / 'isoprene        ', 68.12  /

C MTP
      DATA     SPCA_SPC(  2), SPCA_MWT(  2) / 'myrcene         ', 136.23 /
      DATA     SPCA_SPC(  3), SPCA_MWT(  3) / 'sabinene        ', 136.23 /
      DATA     SPCA_SPC(  4), SPCA_MWT(  4) / 'limonene        ', 136.23 /
      DATA     SPCA_SPC(  5), SPCA_MWT(  5) / 'carene_3        ', 136.23 /
      DATA     SPCA_SPC(  6), SPCA_MWT(  6) / 'ocimene_t_b     ', 136.23 /
      DATA     SPCA_SPC(  7), SPCA_MWT(  7) / 'pinene_b        ', 136.23 /
      DATA     SPCA_SPC(  8), SPCA_MWT(  8) / 'pinene_a        ', 136.23 /
C Other MT
      DATA     SPCA_SPC(  9), SPCA_MWT(  9) / '2met_styrene    ', 132.20 /
      DATA     SPCA_SPC( 10), SPCA_MWT( 10) / 'cymene_p        ', 134.22 /
      DATA     SPCA_SPC( 11), SPCA_MWT( 11) / 'cymene_o        ', 134.22 /
      DATA     SPCA_SPC( 12), SPCA_MWT( 12) / 'phellandrene_a  ', 136.23 /
      DATA     SPCA_SPC( 13), SPCA_MWT( 13) / 'thujene_a       ', 136.23 /
      DATA     SPCA_SPC( 14), SPCA_MWT( 14) / 'terpinene_a     ', 136.23 /
      DATA     SPCA_SPC( 15), SPCA_MWT( 15) / 'terpinene_g     ', 136.23 /
      DATA     SPCA_SPC( 16), SPCA_MWT( 16) / 'terpinolene     ', 136.23 /
      DATA     SPCA_SPC( 17), SPCA_MWT( 17) / 'phellandrene_b  ', 136.23 /
      DATA     SPCA_SPC( 18), SPCA_MWT( 18) / 'camphene        ', 136.23 /
      DATA     SPCA_SPC( 19), SPCA_MWT( 19) / 'bornene         ', 136.23 /
      DATA     SPCA_SPC( 20), SPCA_MWT( 20) / 'fenchene_a      ', 136.23 /
      DATA     SPCA_SPC( 21), SPCA_MWT( 21) / 'ocimene_al      ', 136.23 /
      DATA     SPCA_SPC( 22), SPCA_MWT( 22) / 'ocimene_c_b     ', 136.23 /
      DATA     SPCA_SPC( 23), SPCA_MWT( 23) / 'tricyclene      ', 136.23 /
      DATA     SPCA_SPC( 24), SPCA_MWT( 24) / 'estragole       ', 148.20 /
      DATA     SPCA_SPC( 25), SPCA_MWT( 25) / 'camphor         ', 152.23 /
      DATA     SPCA_SPC( 26), SPCA_MWT( 26) / 'fenchone        ', 152.23 /
      DATA     SPCA_SPC( 27), SPCA_MWT( 27) / 'piperitone      ', 152.23 /
      DATA     SPCA_SPC( 28), SPCA_MWT( 28) / 'thujone_a       ', 152.23 /
      DATA     SPCA_SPC( 29), SPCA_MWT( 29) / 'thujone_b       ', 152.23 /
      DATA     SPCA_SPC( 30), SPCA_MWT( 30) / 'cineole_1_8     ', 154.25 /
      DATA     SPCA_SPC( 31), SPCA_MWT( 31) / 'borneol         ', 154.25 /
      DATA     SPCA_SPC( 32), SPCA_MWT( 32) / 'linalool        ', 154.25 /
      DATA     SPCA_SPC( 33), SPCA_MWT( 33) / 'terpineol_4     ', 154.25 /
      DATA     SPCA_SPC( 34), SPCA_MWT( 34) / 'terpineol_a     ', 154.25 /
      DATA     SPCA_SPC( 35), SPCA_MWT( 35) / 'linalool_OXD_c  ', 170.25 /
      DATA     SPCA_SPC( 36), SPCA_MWT( 36) / 'linalool_OXD_t  ', 170.25 /
      DATA     SPCA_SPC( 37), SPCA_MWT( 37) / 'ionone_b        ', 192.30 /
      DATA     SPCA_SPC( 38), SPCA_MWT( 38) / 'bornyl_ACT      ', 196.29 /

C SQT
      DATA     SPCA_SPC( 39), SPCA_MWT( 39) / 'farnescene_a    ', 204.35 /
      DATA     SPCA_SPC( 40), SPCA_MWT( 40) / 'caryophyllene_b ', 204.35 /
C Other SQT
      DATA     SPCA_SPC( 41), SPCA_MWT( 41) / 'acoradiene      ', 204.35 /
      DATA     SPCA_SPC( 42), SPCA_MWT( 42) / 'aromadendrene   ', 204.35 /
      DATA     SPCA_SPC( 43), SPCA_MWT( 43) / 'bergamotene_a   ', 204.35 /
      DATA     SPCA_SPC( 44), SPCA_MWT( 44) / 'bergamotene_b   ', 204.35 /
      DATA     SPCA_SPC( 45), SPCA_MWT( 45) / 'bisabolene_a    ', 204.35 /
      DATA     SPCA_SPC( 46), SPCA_MWT( 46) / 'bisabolene_b    ', 204.35 /
      DATA     SPCA_SPC( 47), SPCA_MWT( 47) / 'bourbonene_b    ', 204.35 /
      DATA     SPCA_SPC( 48), SPCA_MWT( 48) / 'cadinene_d      ', 204.35 /
      DATA     SPCA_SPC( 49), SPCA_MWT( 49) / 'cadinene_g      ', 204.35 /
      DATA     SPCA_SPC( 50), SPCA_MWT( 50) / 'cedrene_a       ', 204.35 /
      DATA     SPCA_SPC( 51), SPCA_MWT( 51) / 'copaene_a       ', 204.35 /
      DATA     SPCA_SPC( 52), SPCA_MWT( 52) / 'cubebene_a      ', 204.35 /
      DATA     SPCA_SPC( 53), SPCA_MWT( 53) / 'cubebene_b      ', 204.35 /
      DATA     SPCA_SPC( 54), SPCA_MWT( 54) / 'elemene_b       ', 204.35 /
      DATA     SPCA_SPC( 55), SPCA_MWT( 55) / 'farnescene_b    ', 204.35 /
      DATA     SPCA_SPC( 56), SPCA_MWT( 56) / 'germacrene_B    ', 204.35 /
      DATA     SPCA_SPC( 57), SPCA_MWT( 57) / 'germacrene_D    ', 204.35 /
      DATA     SPCA_SPC( 58), SPCA_MWT( 58) / 'gurjunene_b     ', 204.35 /
      DATA     SPCA_SPC( 59), SPCA_MWT( 59) / 'humulene_a      ', 204.35 /
      DATA     SPCA_SPC( 60), SPCA_MWT( 60) / 'humulene_g      ', 204.35 /
      DATA     SPCA_SPC( 61), SPCA_MWT( 61) / 'isolongifolene  ', 204.35 /
      DATA     SPCA_SPC( 62), SPCA_MWT( 62) / 'longifolene     ', 204.35 /
      DATA     SPCA_SPC( 63), SPCA_MWT( 63) / 'longipinene     ', 204.35 /
      DATA     SPCA_SPC( 64), SPCA_MWT( 64) / 'muurolene_a     ', 204.35 /
      DATA     SPCA_SPC( 65), SPCA_MWT( 65) / 'muurolene_g     ', 204.35 /
      DATA     SPCA_SPC( 66), SPCA_MWT( 66) / 'selinene_b      ', 204.35 /
      DATA     SPCA_SPC( 67), SPCA_MWT( 67) / 'selinene_d      ', 204.35 /
      DATA     SPCA_SPC( 68), SPCA_MWT( 68) / 'nerolidol_c     ', 222.37 /
      DATA     SPCA_SPC( 69), SPCA_MWT( 69) / 'nerolidol_t     ', 222.37 /
      DATA     SPCA_SPC( 70), SPCA_MWT( 70) / 'cedrol          ', 222.37 /

C VOC
      DATA     SPCA_SPC( 71), SPCA_MWT( 71) / 'MBO_2m3e2ol     ', 86.13  /
      DATA     SPCA_SPC( 72), SPCA_MWT( 72) / 'methanol        ', 32.04  /
      DATA     SPCA_SPC( 73), SPCA_MWT( 73) / 'acetone         ', 58.08  /
      DATA     SPCA_SPC( 74), SPCA_MWT( 74) / 'methane         ', 16.04  /
C Ammonia, NO2, and NO
      DATA     SPCA_SPC( 75), SPCA_MWT( 75) / 'ammonia         ', 17.03  /
      DATA     SPCA_SPC( 76), SPCA_MWT( 76) / 'nitrous_OXD     ', 44.01  /
      DATA     SPCA_SPC( 77), SPCA_MWT( 77) / 'nitric_OXD      ', 30.01  /
C Acetaldehyde + ethanol
      DATA     SPCA_SPC( 78), SPCA_MWT( 78) / 'acetaldehyde    ', 44.05  /
      DATA     SPCA_SPC( 79), SPCA_MWT( 79) / 'ethanol         ', 46.07  /
C Formic acid + formaldehyde + acetic acid
      DATA     SPCA_SPC( 80), SPCA_MWT( 80) / 'formic_acid     ', 46.03  /
      DATA     SPCA_SPC( 81), SPCA_MWT( 81) / 'formaldehyde    ', 30.03  /
      DATA     SPCA_SPC( 82), SPCA_MWT( 82) / 'acetic_acid     ', 60.05  /
C Other VC
      DATA     SPCA_SPC( 83), SPCA_MWT( 83) / 'MBO_3m2e1ol     ', 86.13  /
      DATA     SPCA_SPC( 84), SPCA_MWT( 84) / 'MBO_3m3e1ol     ', 86.13  /
      DATA     SPCA_SPC( 85), SPCA_MWT( 85) / 'benzaldehyde    ', 106.12 /
      DATA     SPCA_SPC( 86), SPCA_MWT( 86) / 'butanone_2      ', 72.11  /
      DATA     SPCA_SPC( 87), SPCA_MWT( 87) / 'decanal         ', 156.27 /
      DATA     SPCA_SPC( 88), SPCA_MWT( 88) / 'dodecene_1      ', 168.32 /
      DATA     SPCA_SPC( 89), SPCA_MWT( 89) / 'geranyl_acetone ', 194.31 /
      DATA     SPCA_SPC( 90), SPCA_MWT( 90) / 'heptanal        ', 114.19 /
      DATA     SPCA_SPC( 91), SPCA_MWT( 91) / 'heptane         ', 100.20 /
      DATA     SPCA_SPC( 92), SPCA_MWT( 92) / 'hexane          ', 86.18  /
      DATA     SPCA_SPC( 93), SPCA_MWT( 93) / 'met_benzoate    ', 136.15 /
      DATA     SPCA_SPC( 94), SPCA_MWT( 94) / 'met_heptenone   ', 126.20 /
      DATA     SPCA_SPC( 95), SPCA_MWT( 95) / 'neryl_acetone   ', 194.31 /
      DATA     SPCA_SPC( 96), SPCA_MWT( 96) / 'nonanal         ', 142.24 /
      DATA     SPCA_SPC( 97), SPCA_MWT( 97) / 'nonenal         ', 140.22 /
      DATA     SPCA_SPC( 98), SPCA_MWT( 98) / 'octanal         ', 128.21 /
      DATA     SPCA_SPC( 99), SPCA_MWT( 99) / 'octanol         ', 130.23 /
      DATA     SPCA_SPC(100), SPCA_MWT(100) / 'octenol_1e3ol   ', 128.21 /
      DATA     SPCA_SPC(101), SPCA_MWT(101) / 'oxopentanal     ', 100.12 /
      DATA     SPCA_SPC(102), SPCA_MWT(102) / 'pentane         ', 72.15  /
      DATA     SPCA_SPC(103), SPCA_MWT(103) / 'phenyl_CCO      ', 120.15 /
      DATA     SPCA_SPC(104), SPCA_MWT(104) / 'pyruvic_acid    ', 88.06  /
      DATA     SPCA_SPC(105), SPCA_MWT(105) / 'terpinyl_ACT_a  ', 196.29 /
      DATA     SPCA_SPC(106), SPCA_MWT(106) / 'tetradecene_1   ', 196.37 /
      DATA     SPCA_SPC(107), SPCA_MWT(107) / 'toluene         ', 92.14  /
      DATA     SPCA_SPC(108), SPCA_MWT(108) / 'carbon_monoxide ', 28.01  /
      DATA     SPCA_SPC(109), SPCA_MWT(109) / 'butene          ', 56.11  /
      DATA     SPCA_SPC(110), SPCA_MWT(110) / 'ethane          ', 30.07  /
      DATA     SPCA_SPC(111), SPCA_MWT(111) / 'ethene          ', 28.05  /
      DATA     SPCA_SPC(112), SPCA_MWT(112) / 'hydrogen_cyanide', 27.03  /
      DATA     SPCA_SPC(113), SPCA_MWT(113) / 'propane         ', 44.10  /
      DATA     SPCA_SPC(114), SPCA_MWT(114) / 'propene         ', 42.08  /
      DATA     SPCA_SPC(115), SPCA_MWT(115) / 'carbon_2s       ', 76.14  /
      DATA     SPCA_SPC(116), SPCA_MWT(116) / 'carbonyl_s      ', 60.08  /
      DATA     SPCA_SPC(117), SPCA_MWT(117) / 'diallyl_2s      ', 146.28 /
      DATA     SPCA_SPC(118), SPCA_MWT(118) / '2met_2s         ', 94.20  /
      DATA     SPCA_SPC(119), SPCA_MWT(119) / '2met_s          ', 62.14  /
      DATA     SPCA_SPC(120), SPCA_MWT(120) / 'met_chloride    ', 50.49  /
      DATA     SPCA_SPC(121), SPCA_MWT(121) / 'met_bromide     ', 94.94  /
      DATA     SPCA_SPC(122), SPCA_MWT(122) / 'met_iodide      ', 141.94 /
      DATA     SPCA_SPC(123), SPCA_MWT(123) / 'hydrogen_s      ', 34.08  /
      DATA     SPCA_SPC(124), SPCA_MWT(124) / 'met_mercaptan   ', 48.11  /
      DATA     SPCA_SPC(125), SPCA_MWT(125) / 'met_propenyl_2s ', 120.24 /
      DATA     SPCA_SPC(126), SPCA_MWT(126) / 'PPPP_2s         ', 148.29 /
      DATA     SPCA_SPC(127), SPCA_MWT(127) / '2met_nonatriene ', 150.26 /
      DATA     SPCA_SPC(128), SPCA_MWT(128) / 'met_salicylate  ', 152.15 /
      DATA     SPCA_SPC(129), SPCA_MWT(129) / 'indole          ', 117.15 /
      DATA     SPCA_SPC(130), SPCA_MWT(130) / 'jasmone         ', 164.24 /
      DATA     SPCA_SPC(131), SPCA_MWT(131) / 'met_jasmonate   ', 224.30 /
      DATA     SPCA_SPC(132), SPCA_MWT(132) / '3met_3DCTT      ', 218.38 /
      DATA     SPCA_SPC(133), SPCA_MWT(133) / 'hexanal         ', 100.16 /
      DATA     SPCA_SPC(134), SPCA_MWT(134) / 'hexanol_1       ', 102.17 /
      DATA     SPCA_SPC(135), SPCA_MWT(135) / 'hexenal_c3      ', 98.14  /
      DATA     SPCA_SPC(136), SPCA_MWT(136) / 'hexenal_t2      ', 98.14  /
      DATA     SPCA_SPC(137), SPCA_MWT(137) / 'hexenol_c3      ', 100.16 /
      DATA     SPCA_SPC(138), SPCA_MWT(138) / 'hexenyl_ACT_c3  ', 142.20 /
