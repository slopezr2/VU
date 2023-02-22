C=======================================================================
C  EFFS_MGN20T138.EXT
C  This include file contains EF fractions for speciation from 20 MEGAN
C  categories to 138 species.  The values in this file have to be in the
C  same order as in MAP_MGN20T138.EXT
C
C  MEGAN v2.02
C  INPUT version 210
C
C  History:
C  Who          When       What
C  ---------------------------------------------------------------------
C  Tan          12/02/06 - Creates this file
C  Guenther A.  08/11/07 - Move from MEGAN v2.0 to MEGAN v2.02 with update on
C                          Nitrogen gas.
C=======================================================================
                                                                                                     
      INTEGER        N_EFFS_SPC
      PARAMETER     (N_EFFS_SPC = 138)        ! Number of chemical species
                                                                                                     
      REAL           EFFS_BT( N_EFFS_SPC )    ! speciation factor for
                                              ! broadleaf
      REAL           EFFS_NT( N_EFFS_SPC )    ! speciation factor for
                                              ! needleleaf
      REAL           EFFS_SB( N_EFFS_SPC )    ! speciation factor for
                                              ! shrub
      REAL           EFFS_HB( N_EFFS_SPC )    ! speciation factor for
                                              ! herb/grass/crop


      DATA    EFFS_BT(  1), EFFS_NT(  1), EFFS_SB(  1), EFFS_HB(  1)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT(  2), EFFS_NT(  2), EFFS_SB(  2), EFFS_HB(  2)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT(  3), EFFS_NT(  3), EFFS_SB(  3), EFFS_HB(  3)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT(  4), EFFS_NT(  4), EFFS_SB(  4), EFFS_HB(  4)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT(  5), EFFS_NT(  5), EFFS_SB(  5), EFFS_HB(  5)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT(  6), EFFS_NT(  6), EFFS_SB(  6), EFFS_HB(  6)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT(  7), EFFS_NT(  7), EFFS_SB(  7), EFFS_HB(  7)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT(  8), EFFS_NT(  8), EFFS_SB(  8), EFFS_HB(  8)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT(  9), EFFS_NT(  9), EFFS_SB(  9), EFFS_HB(  9)
     &      / 0.011468    , 0.005543    , 0.009242    , 0.010395     /
      DATA    EFFS_BT( 10), EFFS_NT( 10), EFFS_SB( 10), EFFS_HB( 10)
     &      / 0.057339    , 0.055432    , 0.046211    , 0.041580     /
      DATA    EFFS_BT( 11), EFFS_NT( 11), EFFS_SB( 11), EFFS_HB( 11)
     &      / 0.034404    , 0.016630    , 0.027726    , 0.031185     /
      DATA    EFFS_BT( 12), EFFS_NT( 12), EFFS_SB( 12), EFFS_HB( 12)
     &      / 0.045872    , 0.055432    , 0.046211    , 0.041580     /
      DATA    EFFS_BT( 13), EFFS_NT( 13), EFFS_SB( 13), EFFS_HB( 13)
     &      / 0.011468    , 0.033259    , 0.036969    , 0.041580     /
      DATA    EFFS_BT( 14), EFFS_NT( 14), EFFS_SB( 14), EFFS_HB( 14)
     &      / 0.057339    , 0.055432    , 0.046211    , 0.041580     /
      DATA    EFFS_BT( 15), EFFS_NT( 15), EFFS_SB( 15), EFFS_HB( 15)
     &      / 0.057339    , 0.055432    , 0.046211    , 0.041580     /
      DATA    EFFS_BT( 16), EFFS_NT( 16), EFFS_SB( 16), EFFS_HB( 16)
     &      / 0.057339    , 0.066519    , 0.055453    , 0.062370     /
      DATA    EFFS_BT( 17), EFFS_NT( 17), EFFS_SB( 17), EFFS_HB( 17)
     &      / 0.057339    , 0.166297    , 0.092421    , 0.103950     /
      DATA    EFFS_BT( 18), EFFS_NT( 18), EFFS_SB( 18), EFFS_HB( 18)
     &      / 0.172018    , 0.249446    , 0.184843    , 0.145530     /
      DATA    EFFS_BT( 19), EFFS_NT( 19), EFFS_SB( 19), EFFS_HB( 19)
     &      / 0.010321    , 0.004989    , 0.008318    , 0.008316     /
      DATA    EFFS_BT( 20), EFFS_NT( 20), EFFS_SB( 20), EFFS_HB( 20)
     &      / 0.003440    , 0.001663    , 0.002773    , 0.004158     /
      DATA    EFFS_BT( 21), EFFS_NT( 21), EFFS_SB( 21), EFFS_HB( 21)
     &      / 0.011468    , 0.005543    , 0.009242    , 0.010395     /
      DATA    EFFS_BT( 22), EFFS_NT( 22), EFFS_SB( 22), EFFS_HB( 22)
     &      / 0.045872    , 0.022173    , 0.036969    , 0.041580     /
      DATA    EFFS_BT( 23), EFFS_NT( 23), EFFS_SB( 23), EFFS_HB( 23)
     &      / 0.011468    , 0.005543    , 0.009242    , 0.010395     /
      DATA    EFFS_BT( 24), EFFS_NT( 24), EFFS_SB( 24), EFFS_HB( 24)
     &      / 0.003440    , 0.001663    , 0.002773    , 0.004158     /
      DATA    EFFS_BT( 25), EFFS_NT( 25), EFFS_SB( 25), EFFS_HB( 25)
     &      / 0.034404    , 0.033259    , 0.046211    , 0.041580     /
      DATA    EFFS_BT( 26), EFFS_NT( 26), EFFS_SB( 26), EFFS_HB( 26)
     &      / 0.011468    , 0.005543    , 0.009242    , 0.010395     /
      DATA    EFFS_BT( 27), EFFS_NT( 27), EFFS_SB( 27), EFFS_HB( 27)
     &      / 0.003440    , 0.001663    , 0.002773    , 0.004158     /
      DATA    EFFS_BT( 28), EFFS_NT( 28), EFFS_SB( 28), EFFS_HB( 28)
     &      / 0.011468    , 0.027716    , 0.046211    , 0.041580     /
      DATA    EFFS_BT( 29), EFFS_NT( 29), EFFS_SB( 29), EFFS_HB( 29)
     &      / 0.002294    , 0.005543    , 0.009242    , 0.010395     /
      DATA    EFFS_BT( 30), EFFS_NT( 30), EFFS_SB( 30), EFFS_HB( 30)
     &      / 0.057339    , 0.011086    , 0.036969    , 0.041580     /
      DATA    EFFS_BT( 31), EFFS_NT( 31), EFFS_SB( 31), EFFS_HB( 31)
     &      / 0.008028    , 0.003880    , 0.006470    , 0.006237     /
      DATA    EFFS_BT( 32), EFFS_NT( 32), EFFS_SB( 32), EFFS_HB( 32)
     &      / 0.137615    , 0.066519    , 0.110906    , 0.124740     /
      DATA    EFFS_BT( 33), EFFS_NT( 33), EFFS_SB( 33), EFFS_HB( 33)
     &      / 0.006881    , 0.003326    , 0.005545    , 0.006237     /
      DATA    EFFS_BT( 34), EFFS_NT( 34), EFFS_SB( 34), EFFS_HB( 34)
     &      / 0.034404    , 0.016630    , 0.027726    , 0.031185     /
      DATA    EFFS_BT( 35), EFFS_NT( 35), EFFS_SB( 35), EFFS_HB( 35)
     &      / 0.006881    , 0.003326    , 0.005545    , 0.006237     /
      DATA    EFFS_BT( 36), EFFS_NT( 36), EFFS_SB( 36), EFFS_HB( 36)
     &      / 0.034404    , 0.016630    , 0.027726    , 0.031185     /
      DATA    EFFS_BT( 37), EFFS_NT( 37), EFFS_SB( 37), EFFS_HB( 37)
     &      / 0.002294    , 0.001109    , 0.001848    , 0.002079     /
      DATA    EFFS_BT( 38), EFFS_NT( 38), EFFS_SB( 38), EFFS_HB( 38)
     &      / 0.001147    , 0.002772    , 0.002773    , 0.002079     /
      DATA    EFFS_BT( 39), EFFS_NT( 39), EFFS_SB( 39), EFFS_HB( 39)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT( 40), EFFS_NT( 40), EFFS_SB( 40), EFFS_HB( 40)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT( 41), EFFS_NT( 41), EFFS_SB( 41), EFFS_HB( 41)
     &      / 0.018570    , 0.015949    , 0.019157    , 0.021858     /
      DATA    EFFS_BT( 42), EFFS_NT( 42), EFFS_SB( 42), EFFS_HB( 42)
     &      / 0.007428    , 0.006380    , 0.007663    , 0.010929     /
      DATA    EFFS_BT( 43), EFFS_NT( 43), EFFS_SB( 43), EFFS_HB( 43)
     &      / 0.083565    , 0.143541    , 0.095785    , 0.098361     /
      DATA    EFFS_BT( 44), EFFS_NT( 44), EFFS_SB( 44), EFFS_HB( 44)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 45), EFFS_NT( 45), EFFS_SB( 45), EFFS_HB( 45)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 46), EFFS_NT( 46), EFFS_SB( 46), EFFS_HB( 46)
     &      / 0.055710    , 0.119617    , 0.067050    , 0.076503     /
      DATA    EFFS_BT( 47), EFFS_NT( 47), EFFS_SB( 47), EFFS_HB( 47)
     &      / 0.027855    , 0.023923    , 0.028736    , 0.032787     /
      DATA    EFFS_BT( 48), EFFS_NT( 48), EFFS_SB( 48), EFFS_HB( 48)
     &      / 0.013928    , 0.011962    , 0.014368    , 0.016393     /
      DATA    EFFS_BT( 49), EFFS_NT( 49), EFFS_SB( 49), EFFS_HB( 49)
     &      / 0.009285    , 0.007974    , 0.009579    , 0.010929     /
      DATA    EFFS_BT( 50), EFFS_NT( 50), EFFS_SB( 50), EFFS_HB( 50)
     &      / 0.005571    , 0.004785    , 0.005747    , 0.005464     /
      DATA    EFFS_BT( 51), EFFS_NT( 51), EFFS_SB( 51), EFFS_HB( 51)
     &      / 0.009285    , 0.007974    , 0.009579    , 0.010929     /
      DATA    EFFS_BT( 52), EFFS_NT( 52), EFFS_SB( 52), EFFS_HB( 52)
     &      / 0.013928    , 0.011962    , 0.014368    , 0.016393     /
      DATA    EFFS_BT( 53), EFFS_NT( 53), EFFS_SB( 53), EFFS_HB( 53)
     &      / 0.009285    , 0.007974    , 0.009579    , 0.010929     /
      DATA    EFFS_BT( 54), EFFS_NT( 54), EFFS_SB( 54), EFFS_HB( 54)
     &      / 0.018570    , 0.015949    , 0.019157    , 0.021858     /
      DATA    EFFS_BT( 55), EFFS_NT( 55), EFFS_SB( 55), EFFS_HB( 55)
     &      / 0.278552    , 0.239234    , 0.287356    , 0.218579     /
      DATA    EFFS_BT( 56), EFFS_NT( 56), EFFS_SB( 56), EFFS_HB( 56)
     &      / 0.009285    , 0.007974    , 0.009579    , 0.010929     /
      DATA    EFFS_BT( 57), EFFS_NT( 57), EFFS_SB( 57), EFFS_HB( 57)
     &      / 0.027855    , 0.023923    , 0.028736    , 0.032787     /
      DATA    EFFS_BT( 58), EFFS_NT( 58), EFFS_SB( 58), EFFS_HB( 58)
     &      / 0.004643    , 0.003987    , 0.004789    , 0.005464     /
      DATA    EFFS_BT( 59), EFFS_NT( 59), EFFS_SB( 59), EFFS_HB( 59)
     &      / 0.139276    , 0.199362    , 0.172414    , 0.163934     /
      DATA    EFFS_BT( 60), EFFS_NT( 60), EFFS_SB( 60), EFFS_HB( 60)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 61), EFFS_NT( 61), EFFS_SB( 61), EFFS_HB( 61)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 62), EFFS_NT( 62), EFFS_SB( 62), EFFS_HB( 62)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 63), EFFS_NT( 63), EFFS_SB( 63), EFFS_HB( 63)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 64), EFFS_NT( 64), EFFS_SB( 64), EFFS_HB( 64)
     &      / 0.013928    , 0.011962    , 0.014368    , 0.016393     /
      DATA    EFFS_BT( 65), EFFS_NT( 65), EFFS_SB( 65), EFFS_HB( 65)
     &      / 0.046425    , 0.039872    , 0.047893    , 0.054645     /
      DATA    EFFS_BT( 66), EFFS_NT( 66), EFFS_SB( 66), EFFS_HB( 66)
     &      / 0.185701    , 0.079745    , 0.114943    , 0.109290     /
      DATA    EFFS_BT( 67), EFFS_NT( 67), EFFS_SB( 67), EFFS_HB( 67)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 68), EFFS_NT( 68), EFFS_SB( 68), EFFS_HB( 68)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 69), EFFS_NT( 69), EFFS_SB( 69), EFFS_HB( 69)
     &      / 0.004643    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 70), EFFS_NT( 70), EFFS_SB( 70), EFFS_HB( 70)
     &      / 0.001857    , 0.001595    , 0.001916    , 0.005464     /
      DATA    EFFS_BT( 71), EFFS_NT( 71), EFFS_SB( 71), EFFS_HB( 71)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT( 72), EFFS_NT( 72), EFFS_SB( 72), EFFS_HB( 72)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT( 73), EFFS_NT( 73), EFFS_SB( 73), EFFS_HB( 73)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT( 74), EFFS_NT( 74), EFFS_SB( 74), EFFS_HB( 74)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT( 75), EFFS_NT( 75), EFFS_SB( 75), EFFS_HB( 75)
     &      / 0.000000    , 0.000000    , 0.000000    , 0.000000     /
      DATA    EFFS_BT( 76), EFFS_NT( 76), EFFS_SB( 76), EFFS_HB( 76)
     &      / 0.000000    , 0.000000    , 0.000000    , 0.000000     /
      DATA    EFFS_BT( 77), EFFS_NT( 77), EFFS_SB( 77), EFFS_HB( 77)
     &      / 1.000000    , 1.000000    , 1.000000    , 1.000000     /
      DATA    EFFS_BT( 78), EFFS_NT( 78), EFFS_SB( 78), EFFS_HB( 78)
     &      / 0.500000    , 0.500000    , 0.500000    , 0.500000     /
      DATA    EFFS_BT( 79), EFFS_NT( 79), EFFS_SB( 79), EFFS_HB( 79)
     &      / 0.500000    , 0.500000    , 0.500000    , 0.500000     /
      DATA    EFFS_BT( 80), EFFS_NT( 80), EFFS_SB( 80), EFFS_HB( 80)
     &      / 0.285714    , 0.285714    , 0.285714    , 0.285714     /
      DATA    EFFS_BT( 81), EFFS_NT( 81), EFFS_SB( 81), EFFS_HB( 81)
     &      / 0.428571    , 0.428571    , 0.428571    , 0.428571     /
      DATA    EFFS_BT( 82), EFFS_NT( 82), EFFS_SB( 82), EFFS_HB( 82)
     &      / 0.285714    , 0.285714    , 0.285714    , 0.285714     /
      DATA    EFFS_BT( 83), EFFS_NT( 83), EFFS_SB( 83), EFFS_HB( 83)
     &      / 0.000516    , 0.000516    , 0.000516    , 0.000516     /
      DATA    EFFS_BT( 84), EFFS_NT( 84), EFFS_SB( 84), EFFS_HB( 84)
     &      / 0.000516    , 0.000516    , 0.000516    , 0.000516     /
      DATA    EFFS_BT( 85), EFFS_NT( 85), EFFS_SB( 85), EFFS_HB( 85)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT( 86), EFFS_NT( 86), EFFS_SB( 86), EFFS_HB( 86)
     &      / 0.001032    , 0.001032    , 0.001032    , 0.001032     /
      DATA    EFFS_BT( 87), EFFS_NT( 87), EFFS_SB( 87), EFFS_HB( 87)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT( 88), EFFS_NT( 88), EFFS_SB( 88), EFFS_HB( 88)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT( 89), EFFS_NT( 89), EFFS_SB( 89), EFFS_HB( 89)
     &      / 0.003095    , 0.003095    , 0.003095    , 0.003095     /
      DATA    EFFS_BT( 90), EFFS_NT( 90), EFFS_SB( 90), EFFS_HB( 90)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT( 91), EFFS_NT( 91), EFFS_SB( 91), EFFS_HB( 91)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT( 92), EFFS_NT( 92), EFFS_SB( 92), EFFS_HB( 92)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT( 93), EFFS_NT( 93), EFFS_SB( 93), EFFS_HB( 93)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT( 94), EFFS_NT( 94), EFFS_SB( 94), EFFS_HB( 94)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT( 95), EFFS_NT( 95), EFFS_SB( 95), EFFS_HB( 95)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT( 96), EFFS_NT( 96), EFFS_SB( 96), EFFS_HB( 96)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT( 97), EFFS_NT( 97), EFFS_SB( 97), EFFS_HB( 97)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT( 98), EFFS_NT( 98), EFFS_SB( 98), EFFS_HB( 98)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT( 99), EFFS_NT( 99), EFFS_SB( 99), EFFS_HB( 99)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT(100), EFFS_NT(100), EFFS_SB(100), EFFS_HB(100)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT(101), EFFS_NT(101), EFFS_SB(101), EFFS_HB(101)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT(102), EFFS_NT(102), EFFS_SB(102), EFFS_HB(102)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(103), EFFS_NT(103), EFFS_SB(103), EFFS_HB(103)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(104), EFFS_NT(104), EFFS_SB(104), EFFS_HB(104)
     &      / 0.002064    , 0.002064    , 0.002064    , 0.002064     /
      DATA    EFFS_BT(105), EFFS_NT(105), EFFS_SB(105), EFFS_HB(105)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(106), EFFS_NT(106), EFFS_SB(106), EFFS_HB(106)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(107), EFFS_NT(107), EFFS_SB(107), EFFS_HB(107)
     &      / 0.002064    , 0.002064    , 0.002064    , 0.002064     /
      DATA    EFFS_BT(108), EFFS_NT(108), EFFS_SB(108), EFFS_HB(108)
     &      / 0.619067    , 0.619067    , 0.619067    , 0.619067     /
      DATA    EFFS_BT(109), EFFS_NT(109), EFFS_SB(109), EFFS_HB(109)
     &      / 0.036112    , 0.036112    , 0.036112    , 0.036112     /
      DATA    EFFS_BT(110), EFFS_NT(110), EFFS_SB(110), EFFS_HB(110)
     &      / 0.002064    , 0.002064    , 0.002064    , 0.002064     /
      DATA    EFFS_BT(111), EFFS_NT(111), EFFS_SB(111), EFFS_HB(111)
     &      / 0.134131    , 0.134131    , 0.134131    , 0.134131     /
      DATA    EFFS_BT(112), EFFS_NT(112), EFFS_SB(112), EFFS_HB(112)
     &      / 0.004127    , 0.004127    , 0.004127    , 0.004127     /
      DATA    EFFS_BT(113), EFFS_NT(113), EFFS_SB(113), EFFS_HB(113)
     &      / 0.001032    , 0.001032    , 0.001032    , 0.001032     /
      DATA    EFFS_BT(114), EFFS_NT(114), EFFS_SB(114), EFFS_HB(114)
     &      / 0.082542    , 0.082542    , 0.082542    , 0.082542     /
      DATA    EFFS_BT(115), EFFS_NT(115), EFFS_SB(115), EFFS_HB(115)
     &      / 0.000310    , 0.000310    , 0.000310    , 0.000310     /
      DATA    EFFS_BT(116), EFFS_NT(116), EFFS_SB(116), EFFS_HB(116)
     &      / 0.000619    , 0.000619    , 0.000619    , 0.000619     /
      DATA    EFFS_BT(117), EFFS_NT(117), EFFS_SB(117), EFFS_HB(117)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(118), EFFS_NT(118), EFFS_SB(118), EFFS_HB(118)
     &      / 0.000310    , 0.000310    , 0.000310    , 0.000310     /
      DATA    EFFS_BT(119), EFFS_NT(119), EFFS_SB(119), EFFS_HB(119)
     &      / 0.001238    , 0.001238    , 0.001238    , 0.001238     /
      DATA    EFFS_BT(120), EFFS_NT(120), EFFS_SB(120), EFFS_HB(120)
     &      / 0.001032    , 0.001032    , 0.001032    , 0.001032     /
      DATA    EFFS_BT(121), EFFS_NT(121), EFFS_SB(121), EFFS_HB(121)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT(122), EFFS_NT(122), EFFS_SB(122), EFFS_HB(122)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(123), EFFS_NT(123), EFFS_SB(123), EFFS_HB(123)
     &      / 0.000516    , 0.000516    , 0.000516    , 0.000516     /
      DATA    EFFS_BT(124), EFFS_NT(124), EFFS_SB(124), EFFS_HB(124)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(125), EFFS_NT(125), EFFS_SB(125), EFFS_HB(125)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(126), EFFS_NT(126), EFFS_SB(126), EFFS_HB(126)
     &      / 0.000103    , 0.000103    , 0.000103    , 0.000103     /
      DATA    EFFS_BT(127), EFFS_NT(127), EFFS_SB(127), EFFS_HB(127)
     &      / 0.020636    , 0.020636    , 0.020636    , 0.020636     /
      DATA    EFFS_BT(128), EFFS_NT(128), EFFS_SB(128), EFFS_HB(128)
     &      / 0.002064    , 0.002064    , 0.002064    , 0.002064     /
      DATA    EFFS_BT(129), EFFS_NT(129), EFFS_SB(129), EFFS_HB(129)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT(130), EFFS_NT(130), EFFS_SB(130), EFFS_HB(130)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT(131), EFFS_NT(131), EFFS_SB(131), EFFS_HB(131)
     &      / 0.002064    , 0.002064    , 0.002064    , 0.002064     /
      DATA    EFFS_BT(132), EFFS_NT(132), EFFS_SB(132), EFFS_HB(132)
     &      / 0.000206    , 0.000206    , 0.000206    , 0.000206     /
      DATA    EFFS_BT(133), EFFS_NT(133), EFFS_SB(133), EFFS_HB(133)
     &      / 0.003095    , 0.003095    , 0.003095    , 0.003095     /
      DATA    EFFS_BT(134), EFFS_NT(134), EFFS_SB(134), EFFS_HB(134)
     &      / 0.003095    , 0.003095    , 0.003095    , 0.003095     /
      DATA    EFFS_BT(135), EFFS_NT(135), EFFS_SB(135), EFFS_HB(135)
     &      / 0.025794    , 0.025794    , 0.025794    , 0.025794     /
      DATA    EFFS_BT(136), EFFS_NT(136), EFFS_SB(136), EFFS_HB(136)
     &      / 0.015477    , 0.015477    , 0.015477    , 0.015477     /
      DATA    EFFS_BT(137), EFFS_NT(137), EFFS_SB(137), EFFS_HB(137)
     &      / 0.015477    , 0.015477    , 0.015477    , 0.015477     /
      DATA    EFFS_BT(138), EFFS_NT(138), EFFS_SB(138), EFFS_HB(138)
     &      / 0.015477    , 0.015477    , 0.015477    , 0.015477     /
