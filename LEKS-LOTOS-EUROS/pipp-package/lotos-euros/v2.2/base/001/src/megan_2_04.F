!###############################################################################
!
! MEGAN 2.04 module
!
! This module provides access to entities of the MEGAN code.
! Constants are imported from the include files, and
! subroutines are formed from pieced of the original main programs.
!
! FILE NAMES
!
!   The source files of the original code were renamed by:
!    - prepending 'megan_2_04_' to the names ;
!    - changing extension '.EXT' to '.inc' which is more
!      common for Fortran include files.
!
!   The following table shows the filenames under this module
!   and the original names within the downloaded source tree:
!
!     MEGAN 2.04 module                MEGAN_V2_0/Code/MEGANv2.04/src/
!     ----------------------------     -----------------------------------
!
!     megan_2_04.F                  <- EMPROC/megan.F
!     megan_2_04_gamma_etc.F         ~ EMPROC/gamma_etc.F
!     megan_2_04_solarangle.F        = EMPROC/solarangle.F
!                                      EMPROC/checkmem.f
!                                      EMPROC/ckgdioapi2d.f
!     megan_2_04_TEMPD_PRM.inc       = EMPROC/INCLDIR/TEMPD_PRM.EXT
!     megan_2_04_SPC_MGN.inc         = EMPROC/INCLDIR/SPC_MGN.EXT
!     megan_2_04_REL_EM_ACT.inc      = EMPROC/INCLDIR/REL_EM_ACT.EXT
!                                      EMPROC/INCLDIR/CONST.EXT
!     megan_2_04_LD_FCT.inc          = EMPROC/INCLDIR/LD_FCT.EXT
!     megan_2_04_PDT_LOS_CP.inc      = EMPROC/INCLDIR/PDT_LOS_CP.EXT
!
!     megan_2_04.F                  <- MECHCONV/mg2mech.F
!                                      MECHCONV/checkmem.f
!                                      MECHCONV/parsline.f
!
!                                      MECHCONV/INCDIR/SPC_MGN.EXT
!                                      MECHCONV/INCDIR/EF_MGN20.EXT
!
!     megan_2_04_SPC_SPCAT.inc       = MECHCONV/INCDIR/SPC_SPCAT.EXT
!     megan_2_04_EFFS_MGN20T138.inc  = MECHCONV/INCDIR/MAP_MGN20T138.EXT
!     megan_2_04_MAP_MGN20T138.inc   = MECHCONV/INCDIR/EFFS_MGN20T138.EXT
!
!     megan_2_04_ALLMECH.inc         = MECHCONV/INCDIR/ALLMECH.EXT
!
!     megan_2_04_SPC_CBMZ.inc        = MECHCONV/INCDIR/SPC_CBMZ.EXT
!     megan_2_04_MAP_CV2CBMZ.inc     = MECHCONV/INCDIR/MAP_CV2CBMZ.EXT
!
!                                      MECHCONV/INCDIR/SPC_SAPRC99.EXT
!                                      MECHCONV/INCDIR/MAP_CV2SAPRC99.EXT
!
!                                      MECHCONV/INCDIR/SPC_SAPRCII.EXT
!                                      MECHCONV/INCDIR/MAP_CV2SAPRCII.EXT
!
!                                      MECHCONV/INCDIR/SPC_RACM.EXT
!                                      MECHCONV/INCDIR/MAP_CV2RACM.EXT
!
!                                      MECHCONV/INCDIR/SPC_RADM2.EXT
!                                      MECHCONV/INCDIR/MAP_CV2RADM2.EXT
!
!     megan_2_04_index1.f            = ioapi-3.1/ioapi/index1.f
!
!   Legend:
!      =  same content
!      ~  minor changes, e.g. not using the IOAPI library
!     <-  new file formed from parts of original codes
!
! HISTORY
!   2012-08, Arjo Segers
!     Created module to access entities from the original MEGAN 2.04 
!     include files. Added subroutines to convert to the tracers
!     of the chemical meachanism (here CBZ).
!
!###############################################################################
!
!  (from the original commen in "megan.F")
!
!   This program computes biogenic emission using input emission 
!   capacity maps and MCIP output variables.
!   The emission capacity map (ECMAP) are gridded in netCDF-IOAPI format
!   with all the daily average PPFD and daily average temperature.
!
!   Note: The projection and input grids of the two files must be 
!   identical.
!
!   Requirement:
!      Requires libnetcdf.a and libioapi.a to compile
!
!      setenv ECMAP <input netCDF-ioapi emission capacity map>
!      setenv EROUT <megan emission output>
!      setenv TEMPVAR <temperature variable>
!      setenv SRADVAR <solar radiation variable>
!      setenv SDATE <start date>
!      setenv STIME <start time>
!      setenv RLENG <run length>
!
!      If SDATE,STIME and RLENG are not set, use values from MET input
!
!   CALL:
!      CHECKMEM
!      CKGDIOAPI2D
!      MODULE GAMMA_ETC
!         GAMMA_LAI
!         GAMMA_P
!         GAMMA_TISOP
!         GAMMA_TNISP
!         GAMMA_A
!         GAMMA_S
!
!   History:
!     Created by Jack Chen 11/04
!     Modified by Tan 11/21/06 for MEGAN v2.0
!     08/14/07 Tan    - Move to MEGANv2.02 with no update
!
!   Scientific algorithm:
!
!             Emission = [EF][GAMMA][RHO]
!           where [EF]    = emission factor (ug/m2h)
!                 [GAMMA] = emission activity factor (non-dimension)
!                 [RHO]   = production and loss within plant canopies
!                           (non-dimensino)
!                 Assumption: [RHO] = 1  (11/27/06)   (See PDT_LOT_CP.inc)
!
!             GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
!           where [GAMMA_CE]  = canopy correction factor
!                 [GAMMA_age] = leaf age correction factor
!                 [GAMMA_SM]  = soil moisture correction factor
!                 Assumption: [GAMMA_SM]  = 1  (11/27/06)
!             GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
!           where [GAMMA_LAI] = leaf area index factor
!                 [GAMMA_P]   = PPFD emission activity factor
!                 [GAMMA_T]   = temperature response factor
!
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age]
!        Derivation:
!             Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P]
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
!           where LDF = light dependent function (non-dimension)
!                               (See LD_FCT.inc)
!
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age]*
!                        { (1-LDF) + [LDF][GAMMA_P] }
!
!!###############################################################################

      module MEGAN_2_04

        implicit none
        
        
        ! --- in/out --------------------------------
        
        private
        
        public    ::  N_PFT
        public    ::  PFT_NC_FNAME, PFT_NC_VNAME
        
        public    ::  N_MGN_SPC         ! number of MEGAN representer species
        public    ::  MGN_SPC           ! names of MEGAN representer species
        
        public    ::  N_MECH_SPC  ! number of output
        public    ::  MECH_SPC    ! (N_MECH_SPC) Mechanism species name (16 char)

        public    ::  MEGAN_2_04_Init, MEGAN_2_04_Done
        public    ::  MEGAN_2_04_EF_FROM_PFT
        public    ::  MEGAN_2_04_EF_FROM_MAP

        public    ::  MEGAN_2_04_GET_GAMMA, MEGAN_2_04_GET_RHO
        

        ! --- includes ------------------------------
        
        INCLUDE 'megan_2_04_SPC_MGN.inc'        ! MEGAN species (20 species)
        INCLUDE 'megan_2_04_EF_MGN20.inc'       ! EF by PFT!

        INCLUDE 'megan_2_04_PDT_LOS_CP.inc'      !  Production and loss within canopy (RHO)
        INCLUDE 'megan_2_04_LD_FCT.inc'          !  Light dependent factor (LDF)

        INCLUDE 'megan_2_04_SPC_SPCAT.inc'      ! Speciated species (138 species)
        INCLUDE 'megan_2_04_MAP_MGN20T138.inc'  ! Speciation map (20 to 138)
        INCLUDE 'megan_2_04_EFFS_MGN20T138.inc' ! EF fraction for Speciation

        ! ** translation to chemical mechanisms
        
        INCLUDE 'megan_2_04_ALLMECH.inc'       !  Number of mechanism species

        !INCLUDE 'megan_2_04_SPC_SAPRCII.EXT'   !  SAPRCII species and their MW
        !INCLUDE 'megan_2_04_MAP_CV2SAPRCII.EXT'!  Conversion map and factors for SAPRCII
        !integer, parameter    ::  N_MECH_SPC = N_SAPRCII_SPC

        !INCLUDE 'megan_2_04_SPC_SAPRC99.EXT'   !  SAPRC99 species and their MW
        !INCLUDE 'megan_2_04_MAP_CV2SAPRC99.EXT'!  Conversion map and factors for SAPRC99
        !integer, parameter    ::  N_MECH_SPC = N_SAPRC99_SPC

        !INCLUDE 'megan_2_04_SPC_RADM2.EXT'     ! RADM2 species and their MW
        !INCLUDE 'megan_2_04_MAP_CV2RADM2.EXT'  !  Conversion map and factors for RADM2
        !integer, parameter    ::  N_MECH_SPC = N_RADM2_SPC

        !INCLUDE 'megan_2_04_SPC_RACM.EXT'      !  RACM species and their MW
        !INCLUDE 'megan_2_04_MAP_CV2RACM.EXT'   !  Conversion map and factors for RACM
        !integer, parameter    ::  N_MECH_SPC = N_RACM_SPC

        !INCLUDE 'megan_2_04_SPC_CBMZ.inc'      !  CBMZ species and their MW
        !INCLUDE 'megan_2_04_MAP_CV2CBMZ.inc'   !  Conversion map and factors for CBMZ
        !integer, parameter    ::  N_MECH_SPC = N_CBMZ_SPC

        INCLUDE 'megan_2_04_SPC_CBM4.inc'      !  CBM4 species and their MW
        INCLUDE 'megan_2_04_MAP_CV2CBM4.inc'   !  Conversion map and factors for CBM4
        integer, parameter    ::  N_MECH_SPC = N_CBM4_SPC
        

        ! --- const ---------------------------------

        ! number of PFT's:
        integer, parameter   ::  N_PFT = 4

        ! indices:
        integer, parameter   :: I_PFT_BT = 1  ! Broadleaf Tree cover
        integer, parameter   :: I_PFT_NT = 2  ! Needleleaf Tree cover
        integer, parameter   :: I_PFT_SB = 3  ! ShruB cover
        integer, parameter   :: I_PFT_HB = 4  ! HerBaceous cover

        ! parameters for .nc files (PFT 2.1) :
        !  ~ in filenames:
        character(len=*), parameter  ::  PFT_NC_FNAME(1:N_PFT) =
     &         (/ 'btr', 
     &            'ntr', 
     &            'shr', 
     &            'hrb' /)
        !  ~ in variable names:
        character(len=*), parameter  ::  PFT_NC_VNAME(1:N_PFT) =
     &         (/ 'Broadleaf_tree ', 
     &            'Needleleaf_tree', 
     &            'Shrub          ', 
     &            'Herbaceous     ' /)

        ! Length of the time step used for aging:
        integer, parameter   :: TSTLEN = 30   ! days


        ! --- var --------------------------------------
        
        ! new storage for per-pft emission factors and speciation,
        ! to allow looping over pft's:
        real, allocatable         ::  EF_PFT  (:,:)  ! (N_MGN_SPC ,N_PFT)  [ug/m2/hr]
        real, allocatable         ::  EFFS_PFT(:,:)  ! (N_SPCA_SPC,N_PFT)  [(g SPCA)/(g MGN)]
        

      contains
      

        ! =================================================================
      

        subroutine MEGAN_2_04_Init( status )

          ! ---in/out --------------------------------

          integer, intent(out)      ::  status
          
          ! --- const ---------------------------------
          
          ! --- local ---------------------------------

          ! --- begin ---------------------------------
          
          ! storage for per-pft emisison factors in new array:
          allocate( EF_PFT(N_MGN_SPC,N_PFT) )
          ! copy:
          EF_PFT(:,I_PFT_BT) = EF_BT
          EF_PFT(:,I_PFT_NT) = EF_NT
          EF_PFT(:,I_PFT_SB) = EF_SB
          EF_PFT(:,I_PFT_HB) = EF_HB
          
          ! storage for per-pft release of intermediate species:
          allocate( EFFS_PFT(N_SPCA_SPC,N_PFT) )
          ! copy:
          EFFS_PFT(:,I_PFT_BT) = EFFS_BT
          EFFS_PFT(:,I_PFT_NT) = EFFS_NT
          EFFS_PFT(:,I_PFT_SB) = EFFS_SB
          EFFS_PFT(:,I_PFT_HB) = EFFS_HB

          ! ok
          status = 0
          
        end subroutine MEGAN_2_04_Init


        ! ***
      

        subroutine MEGAN_2_04_Done( status )
        
          ! ---in/out --------------------------------

          integer, intent(out)      ::  status
          
          ! --- const ---------------------------------
          
          ! --- local ---------------------------------

          ! --- begin ---------------------------------
          
          ! clear:
          deallocate( EF_PFT )
          deallocate( EFFS_PFT )

          ! ok
          status = 0
          
        end subroutine MEGAN_2_04_Done


        ! ***
        
        
        ! 
        ! Compute Emisison Factors from Plant-Functional-Type fractions
        ! for a given MEGAN representar species.
        ! Return values for species of the chemical mechanism.
        !
        
        
        subroutine MEGAN_2_04_EF_FROM_PFT( I_MGN, NCOLS, NROWS, PFT_prcnt, EF_MECH, status )
        
          ! --- in/out ---------------------------
          
          ! index of MGN representer species:
          integer, intent(in)       ::  I_MGN
          ! grid dimensions:
          integer                   ::  NCOLS, NROWS
          ! PFT fraction in % per cell
          real, intent(in)          ::  PFT_prcnt(NCOLS,NROWS,N_PFT) 
          ! emisison factor (umol/m2/hr) for this MEGAN representar
          ! species derived from Plan-Function-Type maps:
          real, intent(out)         ::  EF_MECH(NCOLS,NROWS,N_MECH_SPC)
          ! return status:
          integer, intent(out)      ::  status
          
          ! --- local -----------------------------
          
          integer            ::  i, j
          integer            ::  S
          integer            ::  I_PFT
          integer            ::  I_SPCA
          integer            ::  I_MECH
          real               ::  EF_SPCA(NCOLS,NROWS,N_SPCA_SPC)
          
          ! --- begin -----------------------------
          
          ! init to zero:
          EF_SPCA = 0.0
          
          ! loop over distributions from megan (20) to intermediate (138) species:
          do S = 1, N_SMAP_SPC
          
            ! skip if origin is not current MEGAN species ...
            if ( MG20_MAP(S) /= I_MGN ) cycle
            
            ! target intermediate species:
            I_SPCA = SPCA_MAP(S)
            
            ! loop over Plant-Functional-Types:
            do I_PFT = 1, N_PFT
              ! add contribution of this PFT for each cell:
              EF_SPCA(:,:,I_SPCA) = EF_SPCA(:,:,I_SPCA) +
     &                                EF_PFT(I_MGN ,I_PFT) *     ! (ug mgn)/m2/hr
     &                              EFFS_PFT(I_SPCA,I_PFT) *     ! (g spca)/(g mgn)
     &                               PFT_prcnt(:,:,I_PFT)*0.01   ! 0-1
            end do
            
          end do  ! distribution from MGN20 to SPCA138
          
          ! loop over grid cells:
          do j = 1, NROWS
            do i = 1, NCOLS
              ! convert to chemical mechanism:
              !                    ug/m2/hr        umol/m2/hr
              call SPCA_TO_MECH( EF_SPCA(i,j,:), EF_MECH(i,j,:) )
            end do
          end do
          
          ! ok
          status = 0
          
        end subroutine MEGAN_2_04_EF_FROM_PFT


        ! ***
        
        
        ! 
        ! Compute Emisison Factors for species in the chemical meachism
        ! given a map of grid-cell-average emission factors 
        ! and the (underlying?) maps of Plant-Functional-Type fractions.
        ! Return values for species of the chemical mechanism.
        !
        
        
        subroutine MEGAN_2_04_EF_FROM_MAP( I_MGN, NCOLS, NROWS, EF_MGN, PFT_prcnt, EF_MECH, status )
        
          ! --- in/out ---------------------------
          
          ! index of MGN representer species:
          integer, intent(in)       ::  I_MGN
          ! grid dimensions:
          integer                   ::  NCOLS, NROWS
          ! grid cell average emission factors (ug/m2/hr):
          real, intent(in)          ::  EF_MGN(NCOLS,NROWS)
          ! PFT fraction in % per cell
          real, intent(in)          ::  PFT_prcnt(NCOLS,NROWS,N_PFT) 
          ! emisison factor (umol/m2/hr) for this MEGAN representar
          ! species derived from Plan-Function-Type maps:
          real, intent(out)         ::  EF_MECH(NCOLS,NROWS,N_MECH_SPC)
          ! return status:
          integer, intent(out)      ::  status
          
          ! --- local -----------------------------
          
          integer            ::  i, j
          integer            ::  S
          integer            ::  I_PFT
          integer            ::  I_SPCA
          integer            ::  I_MECH
          real               ::  PFT_total(NCOLS,NROWS)
          real               ::  EF_SPCA(NCOLS,NROWS,N_SPCA_SPC)
          
          ! --- begin -----------------------------
          
          ! add all fractions:
          PFT_total = sum(PFT_prcnt,3)
          
          ! distribution over 138 intermediate species, init to zero:
          EF_SPCA = 0.0
          
          ! loop over distributions from megan (20) to intermediate (138) species:
          do S = 1, N_SMAP_SPC
          
            ! skip if origin is not current MEGAN species ...
            if ( MG20_MAP(S) /= I_MGN ) cycle
            
            ! target intermediate species:
            I_SPCA = SPCA_MAP(S)
            
            ! loop over Plant-Functional-Types:
            do I_PFT = 1, N_PFT
              ! add relative contribution of this PFT for each cell,
              ! devide later on by sum over all PFT fractions:
              EF_SPCA(:,:,I_SPCA) = EF_SPCA(:,:,I_SPCA) +   ! (ug/m2/hr)*%
     &                                  EF_MGN(:,:) *          ! (ug/m2/hr)
     &                              EFFS_PFT(I_SPCA,I_PFT) *   ! (g spca)/(g mgn)
     &                               PFT_prcnt(:,:,I_PFT)      ! %
            end do
            ! normalize, result is already zero if all PFT fractions are zero:
            where ( PFT_total > 0.0 )
              !  (ug/m2/hr)       =    (ug/m2/hr)*%     /     %
              EF_SPCA(:,:,I_SPCA) = EF_SPCA(:,:,I_SPCA) / PFT_total
            end where
            
          end do  ! distribution from MGN20 to SPCA138
          
          ! loop over grid cells:
          do j = 1, NROWS
            do i = 1, NCOLS
              ! convert to chemical mechanism:
              !                    ug/m2/hr        umol/m2/hr
              call SPCA_TO_MECH( EF_SPCA(i,j,:), EF_MECH(i,j,:) )
            end do
          end do
          
          ! ok
          status = 0
          
        end subroutine MEGAN_2_04_EF_FROM_MAP


        ! ***
        
        
        !
        ! Convert from g of intermediate species 
        ! to mol of chemical mechanism species.
        !
        
        subroutine SPCA_TO_MECH( SPCA_MASS, MECH_MOLE )
        
          ! --- in/out ---------------------------
          
          real, intent(in)          ::  SPCA_MASS(N_SPCA_SPC)   ! species mass (g)
          real, intent(out)         ::  MECH_MOLE(N_MECH_SPC)   ! mechanism (mol)
          
          ! --- local -----------------------------
          
          integer            ::  S
          integer            ::  I_SPCA
          integer            ::  I_MECH
          
          ! --- begin -----------------------------

          ! init result to zero:
          MECH_MOLE = 0.0
          ! lumping to mechanism species
          do S = 1, N_SCON_SPC
            ! source and target indices:
            I_SPCA = SPMH_MAP(S)         ! Mapping value for SPCAT
            I_MECH = MECH_MAP(S)         ! Mapping value for mechanism
            ! add contribution:
            !    (umol mech)/m2/hr
            MECH_MOLE(I_MECH) = MECH_MOLE(I_MECH) + 
     &                SPCA_MASS(I_SPCA) / SPCA_MWT(I_SPCA) * CONV_FAC(S)
                      !   ug/m2/hr    /    (g/mol)    (mol mech)/(mol spca)
          end do ! end species loop
          
        end subroutine SPCA_TO_MECH
        
        
        !
        ! ***
        !
        
        
        subroutine MEGAN_2_04_GET_GAMMA( I_MGN, IDATE, ITIME, 
     &                        NCOLS, NROWS, LONG, LAT, 
     &                        LAIp, LAIc, TEMP, D_TEMP, PPFD, D_PPFD,
     &                        GAMMA, status )
        
          ! Module containing gamma functions
          use GAMMA_ETC, only : GAMMA_TISOP, GAMMA_TNISP
          use GAMMA_ETC, only : GAMMA_LAI, GAMMA_P, GAMMA_A, GAMMA_S

          ! --- in/out ---------------------------
          
          ! index of MGN representer species:
          integer, intent(in)       ::  I_MGN
          ! date and time:
          integer, intent(in)       ::  IDATE   ! YYYYDDD
          integer, intent(in)       ::  ITIME   ! HHMMSS
          ! grid dimensions:
          integer                   ::  NCOLS, NROWS
          ! grid cell locations:
          real, intent(in)          ::  LONG(NCOLS,NROWS)  ! degree east
          real, intent(in)          ::  LAT (NCOLS,NROWS)  ! degree north
          ! leaf-area-indices
          real, intent(in)          ::  LAIp(NCOLS,NROWS)  ! m2/m2 ; previous month
          real, intent(in)          ::  LAIc(NCOLS,NROWS)  ! m2/m2 ; current month
          ! temperature
          real, intent(in)          ::    TEMP(NCOLS,NROWS)  ! K, current hour
          real, intent(in)          ::  D_TEMP(NCOLS,NROWS)  ! K, daily average
          ! radiance:
          real, intent(in)          ::    PPFD(NCOLS,NROWS)  ! W/m2, current hour
          real, intent(in)          ::  D_PPFD(NCOLS,NROWS)  ! W/m2, daily average
          ! combined factors:
          real, intent(out)         ::  GAMMA(NCOLS,NROWS)  ! 1
          ! return status:
          integer, intent(out)      ::  status
          
          ! --- const -----------------------------
          
          character(len=*), parameter  ::  rname = 'MEGAN_2_04_GET_GAMMA'
          
          ! --- local -----------------------------
          
          real           ::  GAM_PHO(NCOLS,NROWS) ! light correction factor
          real           ::  GAM_TMP(NCOLS,NROWS) ! temperature correction factor
          real           ::  GAM_LHT(NCOLS,NROWS) ! LAI correction factor
          real           ::  GAM_AGE(NCOLS,NROWS) ! leaf age correction factor
          real           ::  GAM_SMT(NCOLS,NROWS) ! Soil moilture correction factor
          real           ::  LDF                  ! Light dependent factor

          ! --- begin -----------------------------
          
          ! Select algorithms for different chemical species
          ! Due to the differences in the calculation for gamma,
          ! this logical condition will check for species and choose
          ! (call) the right algorithms to calculate the gamma and
          ! return back the main.
          SELECT CASE ( trim(MGN_SPC(I_MGN)) )
            CASE ( 'ISOP' )
              CALL GAMMA_TISOP( NCOLS, NROWS, TEMP, D_TEMP, GAM_TMP )
            CASE (  'MBO','MYRC','SABI','LIMO','3CAR','OCIM','BPIN',
     &              'APIN','FARN','BCAR','MEOH','ACTO','ACTA','FORM',
     &              'CH4',  'NO','OMTP','OSQT',  'CO' )
              CALL GAMMA_TNISP( NCOLS, NROWS, MGN_SPC(I_MGN), TEMP, GAM_TMP )
            CASE DEFAULT
              write (*,'("ERROR - unsupported representer species for Gamma_temper : ",a)') trim(MGN_SPC(I_MGN))
              write (*,'("ERROR in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
              status=1; return
          END SELECT
          
          ! leaf-area-index dependent:
          CALL GAMMA_LAI( NCOLS, NROWS, LAIc, GAM_LHT )

          ! photolysis factor:
          CALL GAMMA_P( IDATE, ITIME, LAT, LONG, NCOLS, NROWS,
     &                    PPFD, D_PPFD, GAM_PHO )
          ! sometimes negatives for low solar angles, remove:
          GAM_PHO = max( 0.0, GAM_PHO )
          
          ! aging:
          CALL GAMMA_A( IDATE, ITIME, NCOLS, NROWS, MGN_SPC(I_MGN),
     &                    LAIp, LAIc, TSTLEN, D_TEMP, GAM_AGE )

          ! soil-moisture:
          CALL GAMMA_S( NCOLS, NROWS, GAM_SMT )

          ! get Light-Depended-Factor given name of MEGAN species:
          call MEGAN_2_04_GET_LDF( MGN_SPC(I_MGN), LDF, status )
          if (status/=0) then
            write (*,'("ERROR in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
            status=1; return
          end if

          ! combined factor:
          GAMMA = GAM_TMP * GAM_AGE * GAM_LHT * GAM_SMT * ( (1.0-LDF) + GAM_PHO*LDF )

          ! ok
          status = 0
          
        end subroutine MEGAN_2_04_GET_GAMMA

        !
        ! ***
        !
        
        ! Return Light-Depended-Factor given name of MEGAN species.
        
        subroutine MEGAN_2_04_GET_LDF( NAME, LDF, status )
        
          ! --- in/out ---------------------------
          
          character(len=*), intent(in)   ::  NAME  ! MGN name
          real, intent(out)              ::  LDF   ! Light-Depended-Factor
          integer, intent(out)           ::  status
          
          ! --- const -----------------------------
          
          character(len=*), parameter  ::  rname = 'MEGAN_2_04_GET_LDF'
          
          ! --- external --------------------------
          
          integer, external     ::  INDEX1
          
          ! --- local -----------------------------
          
          integer        ::  NMAP
          
          ! --- begin -----------------------------
          
          ! search for index in list of MEGAN names:
          NMAP = INDEX1( NAME, N_LDF_SPC, LDF_SPC )
          if ( NMAP <= 0 ) then
            write (*,'("ERROR - could not find entry for `",a,"` in LDF_SPC list")') trim(name)
            write (*,'("ERROR in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
            status=1; return
          end if

          ! extract corresponding entry:
          LDF = LDF_FCT(NMAP)
          
          ! ok
          status = 0
          
        end subroutine MEGAN_2_04_GET_LDF
        
        
        !
        ! ***
        !
        

        ! Return factor for production and loss within plant canopies.
        
        subroutine MEGAN_2_04_GET_RHO( NAME, RHO, status )
        
          ! --- in/out ---------------------------
          
          character(len=*), intent(in)   ::  NAME  ! MGN name
          real, intent(out)              ::  RHO   ! Light-Depended-Factor
          integer, intent(out)           ::  status
          
          ! --- const -----------------------------
          
          character(len=*), parameter  ::  rname = 'MEGAN_2_04_GET_RHO'
          
          ! --- external --------------------------
          
          integer, external     ::  INDEX1
          
          ! --- local -----------------------------
          
          integer        ::  NMAP
          
          ! --- begin -----------------------------
          
          ! search for index in list of MEGAN names:
          NMAP = INDEX1( NAME, N_RHO_SPC, RHO_SPC )
          if ( NMAP <= 0 ) then
            write (*,'("ERROR - could not find entry for `",a,"` in RHO_SPC list")') trim(name)
            write (*,'("ERROR in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__
            status=1; return
          end if

          ! extract corresponding entry:
          RHO = RHO_FCT(NMAP)
          
          ! ok
          status = 0
          
        end subroutine MEGAN_2_04_GET_RHO


      end module MEGAN_2_04

