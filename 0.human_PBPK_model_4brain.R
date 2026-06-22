human.4B.code <- '
$PLUGIN Rcpp mrgx

$PARAM @annotated
//-----Physiological parameters, L was set as 1/5000 of Q, 96 of Baxter T cell PBPK
Q_Lung    :  312000   :  ml/h Dr. Shah Platform PBPK
Q_Spleen  :  4290   :  ml/h Dr. Shah Platform PBPK
Q_GI      :  43274   :  ml/h Dr. Shah Platform PBPK
Q_Liver  :  70824   :  ml/h Dr. Shah Platform PBPK
Q_Kidney :  54600   :  ml/h Dr. Shah Platform PBPK
Q_Brain  :  35568   :  ml/h Dr. Shah Platform PBPK
Q_Tumor  :  33.84   :  ml/h 1996 paper Baxter
Q_Other  :  103409.76   :  ml/h (373-all)



// from the sympcy paper https://www.sciencedirect.com/science/article/pii/S1347436716300088#bib20

Q_bulk  : 5.25 :ml/min 0.3*1/4
Q_SCSF_sink : 7.98 :0.3*0.38 ml/min https://scholar.google.com/scholar_lookup?title=Spinal%20CSF%20absorption%20in%20healthy%20individuals&publication_year=2004&author=M.%20Edsbagge&author=M.%20Tisell&author=L.%20Jacobsson&author=C.%20Wikkelso
Q_CCSF_sink : 13.02 :0.3*0.62 ml/min, https://scholar.google.com/scholar_lookup?title=Spinal%20CSF%20absorption%20in%20healthy%20individuals&publication_year=2004&author=M.%20Edsbagge&author=M.%20Tisell&author=L.%20Jacobsson&author=C.%20Wikkelso
Q_SC : 7.182: ml/min, 0.3*3/4
Q_CS : 15.162: ml/min, sympcy, https://scholar.google.com/scholar_lookup?title=The%20distribution%20of%20medication%20along%20the%20spinal%20canal%20after%20chronic%20intrathecal%20administration&publication_year=1993&author=J.S.%20Kroin&author=A.%20Ali&author=M.%20York&author=R.D.%20Penn


L_Lung  :  624   :  ml/h
L_Spleen  :  8.58   :  ml/h
L_GI  :  86.55   :  ml/h
L_Liver  :  141.65   :  ml/h
L_Kidney  :  109.2   :  ml/h
L_Brain  :  71.14   :  ml/h
L_Tumor  :  1.8   :  ml/h
L_Other  :  206.82   :  ml/h
L_LN  :  1249.73   :  ml/h sum of all the lymphatic flow rates

V_Blood  :  5300   :  ml
V_Lung_v  :  55.48   :  ml
V_Lung_e  :  499.32   :  ml
V_Spleen_v  :  9.79   :  ml
V_Spleen_e  :  180.01   :  ml
V_GI_v  :  36.22   :  ml
V_GI_e  :  1212.08   :  ml
V_Liver_v  :  206.37   :  ml
V_Liver_e  :  1669.73   :  ml
V_Kidney_v  :  115.63   :  ml
V_Kidney_e  :  205.57   :  ml
V_Brain_v  :  58.4   :  ml
V_Brain_e  :  1401.6   :  ml
V_Other_v  :  5095.84   :  ml (others has to be all organs sum which are not included)
V_Other_e  :  61954.66   :  ml (others has to be all organs sum which are not included)
V_LN  :  274   :  ml, from platform paper
//V_CCSF : 160 : ml https://www.mdpi.com/1999-4923/10/1/14
// from the sympcy paper https://www.sciencedirect.com/science/article/pii/S1347436716300088#bib20
V_CCSF : 130 : ml, from https://content.iospress.com/articles/restorative-neurology-and-neuroscience/rnn00229
V_SCSF : 30 : ml, https://content.iospress.com/articles/restorative-neurology-and-neuroscience/rnn00229



// TABLE 2 PARAMETERS AOOXISATED WITH PBPK model 
// These J values and Kdep value are estimated from the PBPK model (model 9, RUN4)
J_Lung  :  0.0000187 :
J_Spleen  :  0.0000508 :
J_Liver  :  0.462 :
J_Kidney  :  0.000000201 :
J_Tumor  :  14368:
// Other J values, from 2019 Shah paper
J_Brain  :  0.164   :  1/h
J_GI  :  0.00000000939   :  1/h
J_Other  :  0.00000000000317   :  1/h

//brain J value
J_CSFBB : 0.00104 :1/h
J_BCSFB : 0.00000261 :1/h


// Other parameters
Density_CAR  :  3496.8 : overal density of CARs on CAR-T cells (numbers/CAR-T cell), Ag^CAR, EGFR # got from experiment 
TargetCell   :  1e8 : Cells/mL, tumor cells/mL of tumor tissue, tumor cell concentration inside tumor extravacular space




// parameters------
// parameters associated with PBPK-PD model 
Kdep_Liver  :  506.4: The 1 st oder of depletion of CAR-T cells from the liver
kon_pop: 0.000000000000027: //EGFR
koff_pop : 0.677:
Kg_Exp_Tumor_pop: 0.00394 : The 1st order exponential growth rate of tumors (1/h), Kg_Exp^(Tumor)
Kkill_max_pop: 0.0316 :The 1st order maximum rate of tuomr growth inhibition (1/h), Kkill_max^(CAR-T), BCMA
IC50_pop: 341.3 : The number of ‘CAR-Target Complexes per tumor cell’ required to achieve 50% of the maximum rate of tumor growth inhibition (number/cell), KC^(CAR-T)_50
Kexp_CART_max_pop: 0.105 : The maximum 1st order rate constant for in vivo CAR-T cell expansion (1/h), T^(Act)_Max
EC50_pop: 0.000000756 : The number of ‘CAR-Target Complexes per tumor cell’ required to achieve 50% of the maximum rate of CAR-T cell expansion (number/cell), EC^Act_50
tau_kill_pop :16.1: The transit time parameter associated with signal transduction of killing signal (h), τ
gamma_act_pop:0.555 : Curve fitting parameter used to determine the exposure-response relationship
gamma_kill_pop: 0.0119 :
psi_pop : 20: //
Kg_linear_pop : 4.291:
TV_max: 142200: The maximum achievable tumor burden value (mm3)



$MAIN
//---------------------------DIFFERENTIAL EQUATION-----------------------------

double Density_TAA = 48093; // changed from experiment we measured 48092.7
double kon = kon_pop;
double koff = koff_pop;
double Kg_Exp_Tumor = Kg_Exp_Tumor_pop;
double Kkill_max =  Kkill_max_pop;
double IC50 =  IC50_pop;
double Kexp_CART_max = Kexp_CART_max_pop;
double EC50 =  EC50_pop;
double tau_kill =  tau_kill_pop; //The transit time parameter associated with signal transduction of killing signal (h), τ
double gamma_act = gamma_act_pop; // from monolix fitting after adding 2d experiment measured density
double gamma_kill = gamma_kill_pop;
double psi = psi_pop;
double Kg_linear= Kg_linear_pop;

$CMT C_Blood C_V_Lung C_E_Lung C_V_Spleen C_E_Spleen C_V_Kidney C_E_Kidney 
    C_V_Brain C_E_Brain C_CCSF C_SCSF C_V_Other C_E_Other C_V_Tumor A_E_Tumor 
     C_V_GI C_E_GI  C_V_Liver C_E_Liver C_LN  K1 K2 K3 K4 C_E_Tumor TE V_Tumor  

$ODE

//—————--------------------- model-------------
//-------------disposition of CAR-Ts in Blood 
// 11.3% of total volume is vascular volume, calculated from 1996 PBPK paper
double V_Tumor_v  =  0.113* V_Tumor;  
double V_Tumor_e  =  0.887* V_Tumor;


dxdt_C_Blood = (-(Q_Lung+L_Lung)*C_Blood + 
                  (Q_GI-L_GI+Q_Spleen-L_Spleen +Q_Liver-L_Liver)*C_V_Liver + 
                  (Q_Kidney-L_Kidney)*C_V_Kidney + 
                  (Q_Brain-L_Brain)*C_V_Brain + 
                  L_LN*C_LN + 
                  (Q_Other-L_Other)*C_V_Other + 
                  (Q_Tumor-L_Tumor)*C_V_Tumor)
                  /V_Blood;
// dxdt_N_Blood = C_Blood * V_Blood;                  

//--------------Disposition of CAR-Ts in Lungs 
// Vascular space
dxdt_C_V_Lung = ((Q_Lung + L_Lung)*C_Blood - J_Lung*C_V_Lung*V_Lung_v -
                Q_Lung*C_V_Lung)/V_Lung_v;
// Extravascular Space
dxdt_C_E_Lung = (J_Lung*C_V_Lung*V_Lung_v - L_Lung*C_E_Lung)/V_Lung_e;

//-------------------Disposition of CAR-Ts in Spleen 
// Vascular space
dxdt_C_V_Spleen = (Q_Spleen*C_V_Lung - J_Spleen*C_V_Spleen*V_Spleen_v - 
                  (Q_Spleen-L_Spleen)*C_V_Spleen)/V_Spleen_v;
// Extravascular Space
dxdt_C_E_Spleen = (J_Spleen*C_V_Spleen*V_Spleen_v - L_Spleen*C_E_Spleen)/V_Spleen_e;

//---------------------Disposition of CAR-Ts in Kidney 
// Vascular space
dxdt_C_V_Kidney = (Q_Kidney*C_V_Lung - J_Kidney*C_V_Kidney*V_Kidney_v - 
                  (Q_Kidney-L_Kidney)*C_V_Kidney)/V_Kidney_v;
// Extravascular Space
dxdt_C_E_Kidney = (J_Kidney*C_V_Kidney*V_Kidney_v - L_Kidney*C_E_Kidney)/V_Kidney_e;

//-----------------Disposition of CAR-Ts in Brain 
// Vascular space
//dxdt_C_V_Brain = (Q_Brain*C_V_Lung - J_Brain*C_V_Brain*V_Brain_v - 
//                  (Q_Brain-L_Brain)*C_V_Brain)/V_Brain_v;

//dxdt_C_V_Brain = ((Q_Brain + Q_Tumor)*C_V_Lung - J_Brain*C_V_Brain*V_Brain_v - 
//                  (Q_Brain-L_Brain)*C_V_Brain)/V_Brain_v;

dxdt_C_V_Brain = (
                (Q_Brain + Q_Tumor)*C_V_Lung 
                - J_Brain*C_V_Brain*V_Brain_v 
                - J_BCSFB*C_V_Brain * V_Brain_v
                - (Q_Brain-L_Brain)*C_V_Brain
                + Q_CCSF_sink * C_CCSF
                + Q_SCSF_sink * C_SCSF
                )/V_Brain_v;                  
                  
// Extravascular Space
//dxdt_C_E_Brain = (J_Brain*C_V_Brain*V_Brain_v - L_Brain*C_E_Brain)/V_Brain_e;

//dxdt_C_E_Brain = (J_Brain*C_V_Brain*V_Brain_v - L_Brain*C_E_Brain - Q_Tumor * C_E_Brain)/V_Brain_e;


dxdt_C_E_Brain = (J_Brain * C_V_Brain * V_Brain_v 
                  + J_CSFBB * C_CCSF * V_CCSF
                  - L_Brain * C_E_Brain
                  - Q_Tumor * C_E_Brain
                  - Q_bulk * C_E_Brain
                  )/V_Brain_e;

// cranial CSF space
dxdt_C_CCSF = (Q_SC * C_SCSF
             + Q_bulk * C_E_Brain
             + J_BCSFB * C_V_Brain * V_Brain_v
             - Q_CS * C_CCSF
             - Q_CCSF_sink * C_CCSF
             - J_CSFBB * C_CCSF * V_CCSF)/V_CCSF;


// spinal CSF space
dxdt_C_SCSF = (Q_CS * C_CCSF 
             - Q_SC * C_SCSF
             - Q_SCSF_sink * C_SCSF)/V_SCSF;

//------------------Disposition of CAR-Ts in Others 
// Vascular space
dxdt_C_V_Other = (Q_Other*C_V_Lung - J_Other*C_V_Other*V_Other_v - 
                  (Q_Other-L_Other)*C_V_Other)/V_Other_v;
// Extravascular Space
dxdt_C_E_Other = (J_Other*C_V_Other*V_Other_v - L_Other*C_E_Other)/V_Other_e;

//---------------------Disposition of CAR-Ts in GI 
// Vascular space
dxdt_C_V_GI = (Q_GI*C_V_Lung - J_GI*C_V_GI*V_GI_v - (Q_GI-L_GI)*C_V_GI)/V_GI_v;
// Extravascular Space
dxdt_C_E_GI = (J_GI*C_V_GI*V_GI_v - L_GI*C_E_GI)/V_GI_e;

//----------------Disposition of CAR-Ts in Liver 
// Vascular space
dxdt_C_V_Liver = (Q_Liver*C_V_Lung  - J_Liver*C_V_Liver*V_Liver_v 
                  - (Q_GI-L_GI + Q_Spleen - L_Spleen + Q_Liver - L_Liver)*C_V_Liver
                  + (Q_GI-L_GI)*C_V_GI
                  + (Q_Spleen-L_Spleen)*C_V_Spleen)/V_Liver_v;
// Extravascular Space
dxdt_C_E_Liver = (J_Liver*C_V_Liver*V_Liver_v - L_Liver*C_E_Liver)/V_Liver_e - 
                  Kdep_Liver * C_E_Liver;

// ----------------------------Disposition of CAR-Ts in LN 
dxdt_C_LN = (L_Liver*C_E_Liver + L_GI*C_E_GI + L_Tumor*(A_E_Tumor/V_Tumor_e - TE/Density_CAR) 
            + L_Other*C_E_Other + L_Brain*C_E_Brain + L_Kidney*C_E_Kidney + 
            L_Spleen*C_E_Spleen + L_Lung*C_E_Lung - L_LN*C_LN)/V_LN;

// ----------------------------Disposition of CAR-Ts in Tumors

// Vascular space in tumor

dxdt_C_V_Tumor = (Q_Tumor*C_E_Brain - J_Tumor*C_V_Tumor*V_Tumor_v - 
                  (Q_Tumor-L_Tumor)*C_V_Tumor)/V_Tumor_v; //unit is #/ml
//dxdt_C_V_Tumor = (Q_Tumor*C_V_Brain - J_Tumor*C_V_Tumor*V_Tumor_v - 
//                  (Q_Tumor-L_Tumor)*C_V_Tumor);



// TE equals concentration of complexes inside tumor Extravascular
// CplxPT eqaulra CAR-Target complexes per tumor cell

double CplxPT = TE/TargetCell; // = TE * V_Tumor_e/ (TargetCell * V_Tumor_e)

double kg = Kg_Exp_Tumor*(1-V_Tumor/TV_max) / pow((1+pow((Kg_Exp_Tumor/Kg_linear*V_Tumor),psi)), 1/psi);

// Kact_CART: expansion of activated CAR-T cells
double Kact_CART = Kexp_CART_max * pow(CplxPT,gamma_act) / (pow(CplxPT,gamma_act) + pow(EC50,gamma_act)); //equation 26 or 4?

//tumor killing functions 

double Kkill = Kkill_max * pow(CplxPT,gamma_kill) / (pow(CplxPT,gamma_kill) + pow(IC50,gamma_kill)); // equation 25

dxdt_C_E_Tumor = (J_Tumor*C_V_Tumor*V_Tumor_v 
                  - L_Tumor*(C_E_Tumor- TE / Density_CAR) 
                  + Kact_CART*C_E_Tumor* V_Tumor_e) / V_Tumor_e;          

dxdt_K1 = (1/tau_kill)*(Kkill - K1);

dxdt_K2 = (1/tau_kill)*(K1 - K2);

dxdt_K3 = (1/tau_kill)*(K2 - K3);

dxdt_K4 = (1/tau_kill)*(K3 - K4);



// equation 1 from paper, interaction between CAR-T cells and tumor cells
dxdt_TE = kon*((A_E_Tumor /V_Tumor_e)*Density_CAR - TE)*(TargetCell*Density_TAA - TE) - koff*TE ; //unit is #/ml

//amounts (# cells/mL) of CAR-T cells in the extravascular space and vascular of the tumor 
dxdt_A_E_Tumor = J_Tumor*C_V_Tumor*V_Tumor_v 
                 - L_Tumor*(A_E_Tumor/V_Tumor_e - TE/Density_CAR) 
                 + Kact_CART * A_E_Tumor ; //unit is #, equation 27
      

dxdt_V_Tumor =  kg* V_Tumor - K4 * V_Tumor; // tumor size, equation 32


$TABLE
capture OutputCART =  C_Blood;//convert #/mL
capture OutputCART_brain_vascular = C_V_Brain; // #/ml
capture OutputCART_tumor =  (A_E_Tumor + C_V_Tumor *V_Tumor_v)/V_Tumor; //#/mL
capture OutputCART_tumor_ex =  A_E_Tumor /V_Tumor_e;
capture OutputCART_SCSF =  C_SCSF;//convert #/mL

capture OutputCART_tumor_vascular =  C_V_Tumor; //#/mL
capture OutputVolume = V_Tumor * 1000;//convert ml to mm3

capture CplxPCART = TE * V_Tumor_e / A_E_Tumor; // equation 2 
'

