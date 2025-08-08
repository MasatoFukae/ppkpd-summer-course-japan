$PROB PKPD (IPP)
$INPUT ID DOSE CYCLE TIME DAY DVID DV MDV EVID AMT RATE DUR ADDL II IL6BSL BSLFL 
G1 G2 G3 G2_3 ICL IV1 IQ IV2

$DATA pkpd02.csv IGNORE=@ IGNORE(DVID.EQ.1) IGNORE(BSLFL.EQ.1)

$SUBROUTINES ADVAN13 TOL=6
$MODEL
COMP=(CENT)
COMP=(PERI)
COMP=(PKAUCCOMP)
COMP=(IL6AUCCOMP)
COMP=(IL6RELEASE)
COMP=(IL6TRANS1)
COMP=(IL6TRANS2)
COMP=(IL6TRANS3)
COMP=(IL6TRANS4)
COMP=(IL6TRANS5)
COMP=(IL6CIRC)

$PK
; Number of dose---------------------
IF(NEWIND.LE.1) NDOSE = 1                        ; NDOSE=0 (B1 method), NDOSE=1 (B2 method) 
IF(NEWIND.EQ.2.AND.EVID.EQ.1) NDOSE = NDOSE + 1
NEWIND2 = NEWIND

; PK--------------------------------------
CL = ICL ; (L/h) 
V1 = IV1 ; (L)
Q  = IQ ; (L/h) 
V2 = IV2 ; (L)

S1 = V1
S2 = V2
 
K12 = Q/V1
K21 = Q/V2
K = CL/V1

; PD-------------------------------------
TVEMAX = THETA(1)
EMAX = TVEMAX * EXP(ETA(1)) ; (pg/mL/h)

TVEC50 = THETA(2)
EC50 = TVEC50 * EXP(ETA(2)) ; (ug/mL)

TVHILL = THETA(3)
HILL = TVHILL * EXP(ETA(3))

TVIMAX = THETA(4)
IMAX = TVIMAX * EXP(ETA(4))

TVIC50 = THETA(5)
IC50 = TVIC50 * EXP(ETA(5)) ; (pg*h/mL)

TVKDEG = THETA(6)
KDEG = TVKDEG * EXP(ETA(6)) ; (/h)

TVKPRIM = THETA(7)
KPRIM = TVKPRIM * EXP(ETA(7))

TVMTT = THETA(8)
MTT = TVMTT * EXP(ETA(8)) ; (h)

PROPERR = THETA(9)
ADDERR = THETA(10)

; Assume baseline deviates from observed value by the same variability of residual error                                        
REWT = SQRT(ADDERR**2 + (PROPERR**2)*(IL6BSL**2))       
IL6BASE = IL6BSL + REWT*ETA(9)

A_0(4) = 0
A_0(5) = IL6BASE
A_0(6) = IL6BASE
A_0(7) = IL6BASE
A_0(8) = IL6BASE
A_0(9) = IL6BASE
A_0(10) = IL6BASE
A_0(11) = IL6BASE

KTR = 5/MTT

$DES
; PK--------------------------------------------
CP = A(1)/S1 ; Drug concentration
PKAUC = A(3) ; cumulative drug AUC

DADT(1) = - K*A(1) - K12*A(1) + K21*A(2)
DADT(2) = K12*A(1) - K21*A(2)
DADT(3) = CP

; PD--------------------------------------------
IL6 = A(11) ; IL6 concentration
IL6AUC = A(4) ; cumulative IL6 AUC

RL =     ; Stimulation effect on IL6 release 
IH =      ; Inhibition effect (negative feedback) on IL6 release 

DADT(4) =    ; Cumulative IL6 AUC
DADT(5) =    ; IL6 release
DADT(6) =    ; Transit compartment 1
DADT(7) =    ; Transit compartment 2
DADT(8) =    ; Transit compartment 3
DADT(9) =    ; Transit compartment 4
DADT(10) =  ; Transit compartment 5
DADT(11) =   ; Plasma IL6

$ERROR
CONC = A(1)/S1 ; drug concentration
IL6CONC = A(11) ; IL6 concentration

IPRED = A(11)
IRES  = DV-IPRED
W = SQRT(IPRED*IPRED*PROPERR**2 + ADDERR**2)
IWRES = IRES/(W+1E-33)
Y = IPRED + W*EPS(1)   

$THETA
(0, 150)      ; TVEMAX
(0, 3)        ; TVEC50
(0, 2.5)      ; TVHILL
(0, 0.75, 1)     ; TVIMAX
(0, 8000)     ; TVIC50
(0, 0.3)      ; TVKDEG
(0, 4)        ; TVKPRIM
(0, 0.5)      ; TVMTT

(0, 0.1)     ; PROP_ERR_IL6
(0 FIXED)      ; ADD_ERR_IL6

$OMEGA
0 FIXED ; IIV_EMAX
0.1 ; IIV_EC50
0.1 ; IIV_HILL
0 FIXED ; IIV_IMAX
0 FIXED ; IIV_IC50
0 FIXED ; IIV_KDEG
0 FIXED ; IIV_KPRIM
0 FIXED ; IIV_MTT
1 FIXED ; Dummy (IL6BASE)


$SIGMA
1 FIXED              ; DUMMY

$EST METHOD=0 POSTHOC MAXEVAL=9999 PRINT=1 ETASTYPE=1

$COV PRINT=E MATRIX=S

$TABLE
ID DOSE CYCLE TIME DAY DVID MDV IPRED IWRES CWRES NPDE
NOPRINT ONEHEADER ESAMPLE=1000
FILE=sdtab001 FORMAT=,1PE15.8

$TABLE
ID DOSE TIME CL V1 Q V2 EMAX EC50 HILL IMAX IC50 KDEG KPRIM IL6BASE MTT ETAS(1:LAST)
NEWIND2 NDOSE PKAUC IL6AUC RL IH CP IL6
NOPRINT NOAPPEND ONEHEADER
FILE=patab001 FORMAT=,1PE15.8

;---------------------------------------------
; THE END OF PROGRAM
;---------------------------------------------
