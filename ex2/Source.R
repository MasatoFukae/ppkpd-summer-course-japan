###############################################################################
# PROGRAM NAME        : Source.R
# DESCRIPTION         : To load necessary R packages and functions
# COMPOUND NAME       : Not Applicable
# PROGRAM VERSION     : 01
# CREATION DATE       : 28-July-2025
# PROGRAMMER          : pkpd_seminor
# EXTERNAL FILES USED : None
###############################################################################

library(tidyverse)
library(readxl)
library(mrgsolve)
library(ggpubr)
library(pracma)  # for correlation statistics in plots
library(drc)
library(patchwork)  # for combining plots
library(gridExtra)  # for arranging multiple plots
library(GGally)
library(MatchIt)

code <- '
  $PROB
  TCE_PKPD
  
  $PLUGIN nm-vars autodec tad
  
  $GLOBAL
  mrg::tadose tad_cmt_1(1);
  mrg::tadose tad_cmt_2(2);
  capture IL6_AUC_TAU = 0;
  
  $PARA
  // PK
  TVCL = 0.02 // L/h
  TVVC = 3 // L
  TVQ  = 0.05 // L/h
  TVVP = 3 // L
  TVKa = 0.188/24 // Absorption rate constant (1/h)
  TVF1 = 0.4 // Bioavailablity for SC
  
  BSA     =   1.88 // Median BSA (m2) - use individual value for simulation
  BSA_REF =   1.88 // Reference BSA (m2)
  BSA_EFF =   0.75 // BSA effect on clearance
  BSA_V_EFF  = 1.00  // BSA effect on volume of distribution
  
  // PD (IL6)
  TVEMAX = 150 // pg/mL/hour
  TVEC50 = 3 // ug/mL
  TVHILL = 2.5 // unitless
  TVIMAX = 0.75 // unitless
  TVIC50 = 8000 // pg*hour/mL
  TVKDEG = 0.3 // 1/hour
  TVKPRIM = 4 // unitless
  TVIL6BASE = 0.695 // pg/mL
  TVMTT = 0.5 // hour
  
  // PD (Tumor Growth)
  TUMOR_REF = 100 // Reference tumor burden (%)
  TUMOR_EFF = -1  // Tumor burden effect on clearance
  BASE_TUMOR = 100 //  Reference tumor burden
  KG_TUMOR = 0.0014   // Tumor growth rate constant (1/h)
  KMAX_TUMOR  = 0.0005 // Maximum kill rate (1/h)
  TVKILL_TUMOR_INTACT = 0.0001 // kill rate by intact (1/h)
  EC50_TUMOR_KILL = 3 // Concentration for half-maximal kill rate (μg/L)
  HILL_TUMOR_KILL = 2 // Hill coefficient for tumor kill
  
  Dose = 3  // Dose - - use individual value for simulation
  Cycle = 1 // Number of Cycle
  
  $CMT @annotated
  ABS  : Absorption compartment
  CENT : Central compartment
  PERI : Peripheral compartment
  AUC : Area Under the Curve
  IL6_AUC : Area Under the Curve for IL6
  IL6_release : IL6 release from the source
  IL6_trans1 : IL6 transition compartment 1
  IL6_trans2 : IL6 transition compartment 2
  IL6_trans3 : IL6 transition compartment 3
  IL6_trans4 : IL6 transition compartment 4
  IL6_trans5 : IL6 transition compartment 5
  IL6_circ : IL6 in circulation
  TUMOR : Tumor growth
  
  $MAIN
  // PK
  CL = TVCL * pow(BSA/BSA_REF, BSA_EFF)   * pow(BASE_TUMOR/TUMOR_REF, TUMOR_EFF) * EXP(ETA(1));
  VC = TVVC * pow(BSA/BSA_REF, BSA_V_EFF) * EXP(ETA(2));
  Q  = TVQ  * pow(BSA/BSA_REF, BSA_EFF)   * EXP(ETA(3));
  VP = TVVP * pow(BSA/BSA_REF, BSA_V_EFF) * EXP(ETA(4));
  Ka = TVKa * EXP(ETA(14));
  F1 = TVF1 * EXP(ETA(15));
  
  KE = CL / VC;
  K12 = Q / VC;
  K21 = Q / VP;
  CENT_0 = 0;
  PERI_0 = 0;
  
  // PD (IL6)
  EMAX = TVEMAX * EXP(ETA(5));
  EC50 = TVEC50 * EXP(ETA(6));
  HILL = TVHILL * EXP(ETA(7));
  IMAX = TVIMAX * EXP(ETA(8));
  IC50 = TVIC50 * EXP(ETA(9));
  KDEG = TVKDEG * EXP(ETA(10));
  KPRIM = TVKPRIM * EXP(ETA(11));
  IL6BASE = TVIL6BASE * EXP(ETA(12));
  MTT = TVMTT * EXP(ETA(13));
  IL6_release_0 = IL6BASE;
  IL6_trans1_0 = IL6BASE;
  IL6_trans2_0 = IL6BASE;
  IL6_trans3_0 = IL6BASE;
  IL6_trans4_0 = IL6BASE;
  IL6_trans5_0 = IL6BASE;
  IL6_circ_0 = IL6BASE;
  IL6_AUC_0 = 0;
  
  KTR = 5 / MTT;
  
  capture tad1 = tad_cmt_1.tad(self);
  capture tad2 = tad_cmt_2.tad(self);
  
  if(NEWIND <= 1){
    NDOSE = 0;
  }
  if(NEWIND == 2 && EVID == 1) {
    NDOSE = NDOSE + 1;
  }
  
  // PD (Tumor growth)
  TUMOR_0 = BASE_TUMOR;
  KILL_TUMOR_INTACT = TVKILL_TUMOR_INTACT * EXP(ETA(16));
  
  $OMEGA
  0.1 // CL OMEGA, 1
  0.1 // VC OMEGA, 2
  0.1 // Q OMEGA, 3
  0.1 // VP OMEGA, 4
  0.2 // EMAX OMEGA, 5
  0.1 // EC50 OMEGA, 6
  0.1 // HILL OMEGA, 7
  0 // IMAX OMEGA, 8
  0.1 // IC50 OMEGA, 9
  0.1 // KDEG OMEGA, 10
  0.1 // KPRIM OMEGA, 11
  0.1 // IL6BASE OMEGA, 12
  0.1 // MTT OMEGA, 13
  0.277  // BSV on Ka, 14
  0.3    // BSV on F1, 15
  0.3    // BSV on tumor kill, 16
  
  $SIGMA
  0.03
  0.03
  0.09
  
  $ODE
  // PK
  dxdt_ABS  = - Ka * ABS; 
  dxdt_CENT = - KE * CENT - K12 * CENT + K21 * PERI + Ka * ABS;
  dxdt_PERI = K12  * CENT - K21 * PERI;
  dxdt_AUC = CENT/VC;
  
  // PD (IL6)
  CONC = CENT / VC;
  dxdt_IL6_AUC = IL6_circ;
  CONC_C = CONC/pow(Cycle, 1.5);
  
  RL = EMAX * pow(CONC_C, HILL) / (pow(EC50, HILL) + pow(CONC_C, HILL));
  IH = IMAX * IL6_AUC / (IC50/pow(KPRIM, NDOSE - 1) + IL6_AUC);
  dxdt_IL6_release = KDEG * IL6BASE + RL * (1 - IH) - KTR * IL6_release;
  dxdt_IL6_trans1 = KTR * IL6_release - KTR * IL6_trans1;
  dxdt_IL6_trans2 = KTR * IL6_trans1 - KTR * IL6_trans2;
  dxdt_IL6_trans3 = KTR * IL6_trans2 - KTR * IL6_trans3;
  dxdt_IL6_trans4 = KTR * IL6_trans3 - KTR * IL6_trans4;
  dxdt_IL6_trans5 = KTR * IL6_trans4 - KTR * IL6_trans5;
  dxdt_IL6_circ = KTR * IL6_trans5 - KDEG * IL6_circ;
  
  // PD (Tumor growth)
  KILL_RATE = KILL_TUMOR_INTACT + (KMAX_TUMOR * pow(CONC, HILL_TUMOR_KILL))/(pow(EC50_TUMOR_KILL, HILL_TUMOR_KILL) + pow(CONC, HILL_TUMOR_KILL));
    dxdt_TUMOR = KG_TUMOR * TUMOR * (1 - TUMOR/100) - KILL_RATE * TUMOR;
  
  $TABLE
  capture CP = (CENT/VC);
  capture CP_OBS = CP * (1 + EPS(1));
  capture CYK_IPRED = IL6_circ;
  capture CYK_OBS = IL6_circ * (1 + EPS(2));
  capture CAVE = AUC/TIME;
  capture TUMOR_IPRED = TUMOR;
  capture TUMOR_OBS = TUMOR * (1 + EPS(3));
  
  $CAPTURE
  NEWIND CL VC Q VP EMAX EC50 HILL IMAX IC50 KDEG KPRIM IL6BASE MTT BSA BSA, BASE_TUMOR, Dose, Cycle
  
'

###############################################################################
#
# End of program
#
###############################################################################