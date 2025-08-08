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
TVQ = 0.05 // L/h
TVVP = 3 // L

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

$CMT @annotated
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

$MAIN
// PK
CL = TVCL * EXP(ETA(1));
VC = TVVC * EXP(ETA(2));
Q = TVQ * EXP(ETA(3));
VP = TVVP * EXP(ETA(4));
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

$OMEGA
0.1 // CL OMEGA
0.1 // VC OMEGA
0.1 // Q OMEGA
0.1 // VP OMEGA
0.2 // EMAX OMEGA
0.1 // EC50 OMEGA
0.1 // HILL OMEGA
0 // IMAX OMEGA
0.1 // IC50 OMEGA
0.1 // KDEG OMEGA
0.1 // KPRIM OMEGA
0.1 // IL6BASE OMEGA
0 // MTT OMEGA

$SIGMA
0.01
0.01

$ODE
// PK
dxdt_CENT = -KE * CENT - K12 * CENT + K21 * PERI;
dxdt_PERI = K12 * CENT - K21 * PERI;
dxdt_AUC = CENT/VC;

// PD (IL6)
CONC = CENT / VC;
dxdt_IL6_AUC = IL6_circ;

RL = EMAX * pow(CONC, HILL) / (pow(EC50, HILL) + pow(CONC, HILL));
IH = IMAX * IL6_AUC / (IC50/pow(KPRIM, NDOSE - 1) + IL6_AUC);
dxdt_IL6_release = KDEG * IL6BASE + RL * (1 - IH) - KTR * IL6_release;
dxdt_IL6_trans1 = KTR * IL6_release - KTR * IL6_trans1;
dxdt_IL6_trans2 = KTR * IL6_trans1 - KTR * IL6_trans2;
dxdt_IL6_trans3 = KTR * IL6_trans2 - KTR * IL6_trans3;
dxdt_IL6_trans4 = KTR * IL6_trans3 - KTR * IL6_trans4;
dxdt_IL6_trans5 = KTR * IL6_trans4 - KTR * IL6_trans5;
dxdt_IL6_circ = KTR * IL6_trans5 - KDEG * IL6_circ;


$TABLE
capture IPRED_CP = (CENT/VC);
capture DV_CP = IPRED_CP * (1 + EPS(1));
capture IPRED_IL6 = IL6_circ;
capture DV_IL6 = IPRED_IL6 * (1 + EPS(2));
capture CAVE = AUC/TIME;

$CAPTURE
NEWIND CL VC Q VP EMAX EC50 HILL IMAX IC50 KDEG KPRIM IL6BASE MTT