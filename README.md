# Ocean2020
Ocean Speciation UROP 2020 

---- Numbering System

   Major Cations    
   Na = 0;       K = 1;     Mg = 2;       Ca = 3;      Sr = 4

   ************************************************************************************

   Minor Cations
   H = 5;      Li = 6;      Rb = 7;       Cs = 8;      NH4 = 9;    Ba = 10;   Mn = 11;    
Fe(II)= 12;     Co = 13;     Ni = 14;   Cu(II)= 15;      Zn = 16;   UO2 = 17;   Be = 18;    
   Cd = 19;  Pb(II)= 20;   Cu(I)= 21;      La = 22;      Ce = 23;    Pr = 24;   Nd = 25;      
   Pm = 26;     Sm = 27;     Eu = 28;      Gd = 29;      Tb = 30;    Dy = 31;   Ho = 32;      
   Er = 33;     Tm = 34;     Yb = 35;      Lu = 36;       Y = 37;    Al = 38;   Ga = 39;      
   In = 40; Fe(III)= 41;

   ************************************************************************************

   Major Anions   
   Cl = 0;     SO4 = 1;     CO3 = 2;     HCO3 = 3;    Br = 4;     F = 5   B(OH)4 = 6

   ************************************************************************************

   Minor Anions
   OH = 7;   HSO4 = 8;     HS = 9;      I = 10;    ClO3  = 11;  ClO4 = 12;    BrO3 = 13;     
   CNS = 14;   NO2 = 15;   NO3 = 16;  H2PO4 = 17;     HPO4 = 18;   PO4 = 19;  H2AsO4 = 20;  
   HAsO4 = 21;  AsO4 = 22;  HSO3 = 23;    SO3 = 24;  Acetate = 25;   
   
   ************************************************************************************

---- Calccomp.ipynb

This function calculates the ionic strength and composition of background electrolyte i.e standard seawater [Millero et al., 2008]. On a molar concentration basis, seven chemical elements account for 99.9% of the dissolved species in seawater. SHIAMI defines the background seawater electrolyte containing 12 elements whose concentrations exceed 1 mol/kg: 

INPUT - S, Tc, P, PAR1, PAR2, PAR1TYPE, PAR2TYPE, SIL, PO4, pHSCALE

S  = Salinity (ppm)
Tc      (degr. C) : scalar or vector of size n x 1
Pressure    (atm) : scalar or vector of size n x 1
PAR1  (some unit) : scalar or vector of size n x 1
PAR2  (some unit) : scalar or vector of size n x 1
PAR1TYPE       () : scalar or vector of size n x 1(*)
PAR2TYPE       () : scalar or vector of size n x 1(*)
SIL    (umol/kgSW): scalar or vector of size n x 1
PO4    (umol/kgSW): scalar or vector of size n x 1
pHSCALE           : scalar or vector of size n x 1(*,*)

(*) Each element must be an integer, indicating that PAR1 (or PAR2) is of type: 
1 = Total Alkalinity
2 = DIC (Total CO2)
3 = pH
4 = pCO2
5 = fCO2

(*,*) Each element must be an integer, indicating that the pH-input (PAR1 or PAR2, if any) is at:
1 = Total scale
2 = Seawater scale
3 = Free scale
4 = NBS scale

OUTPUT - I, mT_cat, mT_an, OH

I              () : Ionic strength (I = 19.924.*S./(1000-1.005*S), S = salinity)
mT_cat   (mol/kg) : Total concentration  of major cations in seawater (Na, K, Mg, Ca, Sr)
mT_an    (mol\kg) : Total concentration  of major anions in seawater  (Cl, SO4, CO3, HCO3, Br, F, B(OH)4)
OH       (mol\kg) : Total concentration  of OH anion in seawater (OH = 10**(-14+pH))

---- CreatePitzerV3.ipynb

The T dependence coefficients that are needed to generate the Pitzer binary interaction at a given T. This was primarily done to reduce the number of lines of code in the "BinaryV3" function and improve readability.

---- BinaryV3.ipynb

This function calculates  binary interaction parameters beta_0, beta_1, beta_2, C_phi at a given temperature and produces arrays to store them. For HSO4, two additional parameters are required C0_HSO4, C1_HSO4. It makes use of coefficients provided in PitzerV3.mat. Alternatively, these coefficients can be found in "CreatePitzerV3.ipynb" file.

INPUT - T

T       (degr. C) : scalar or vector of size n x 1

OUTPUT - beta_0, beta_1, beta_2, C_phi, C0_HSO4, C1_HSO4

These account from short-range interactions which are divided into interactions between dissimilar ions of like charges (Na+ – Mg2+) and opposite 
charges (-Cl – SO42-)

---- TernaryV3.ipynb

This function calculates mixing parameter, theta_Mc and theta_Xa (which account for interactions of like-charged ions) and triplet interaction parameters, phi_AAC and phi_CCA, at a given temperature. The function ouputs these as arrays. These interaction parameters are needed to calculate free activity coefficients using Pitzer equations.

INPUT - T

T       (degr. C) : scalar or vector of size n x 1

OUTPUT - theta_Mc, theta_Xa, phi_CCA, phi_AAC

theta_Mc          : Interactions of like-charged ions between cations
theta_Xa          : Interactions of like-charged ions between anions
phi_CCA           : Interactions between cation-cation-anion 
phi_AAC           : Interactions between anion-anion-cation

---- GammaFV3.ipynb

This function calculates  free ion activity coefficients (free gamma), for the full set of ions  in Millero Model (Millero & Pierrot 1998) and some ion pairs whose pitzer parameters is available, at a given temperature(s), ionic strength(s) and water composition(s).

Free gammas are calculated using a truncated form of the general Pitzer equation (Pitzer,1993) as the higher-order electrostatic terms are relevant only in brines with much higher ionic strength than ordinary seawater. ---> Hain et al 2015//SUPPLIMENTARY INFORMATION. 

INPUT - Tc, I, mT_cat, mT_an

Tc      (degr. C) : Scalar or vector of size n x 1
I              () : Ionic strength (I = 19.924.*S./(1000-1.005*S), S = salinity)
mT_cat   (mol/kg) : Total concentration  of major cations in seawater (Na, K, Mg, Ca, Sr)
mT_an    (mol\kg) : Total concentration  of major anions in seawater  (Cl, SO4, CO3, HCO3, Br, F, B(OH)4)

OUTPUT - gF_cat, gF_an, gIP1, gdivIP, gtrivIP

gF_cat   : an array containing free gammas for the following cations (Free cation activity coefficients): 

   Major Cations                  Minor Cations 
   Na, K, Mg, Ca, Sr              H, Li, Rb, Cs, NH4, Ba, Mn, Fe(II), Co, Ni, 
                                  Cu(II), Zn, UO2, Be, Cd, Pb(II), Cu(I)
         
   REE
   La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Y
                           
   Trivalent
   Al, Ga, In, Fe(III)

gF_an   : an array containing free gammas for the following anions (Free anion activity coefficients):

   Major Anions                    Minor Anions
   Cl, SO4, CO3, HCO3              OH, HSO4, HS,  I, ClO3, ClO4, BrO3, CNS, NO2, NO3, H2PO4, HPO4,
   Br, F, B(OH)4                   PO4, H2AsO4, HAsO4, AsO4, HSO3, SO3, Acetate

   For charged ion pairs whose pitzer parameters are available 

gIP1    : an array containing free gammas for the following ion pairs (Ion pair activity coefficients): 
      
   1)MgOH    2)MgF    3)CaF    4)CuCl    5)CuCl2-   6)CuCl32-                                
           
gdivIP  : an array containing free gammas for the following divalent ion pairs:

   Fe(CO3)   Fe(CO3)2   CuCl2      CuCO3   CuCO32   CuHCO3 
   CuOH      CuOH2      UO2(CO3)3  BeOH3   BeOH4    PbCl+    
   PbCl2     PbCl3-     PbCl42-

   gtrivIP : an array containing the following trivalent ion pairs:
           *all trivalent ion pairs required* 

---- GammaNV3.ipynb

This function calculates  free ion coefficients for the following neutral solutes: NH3,  B(OH)3, H3PO4, H2S, SO2, HF, CO2

INPUT - Tc, I, mT_cat, mT_an

Tc      (degr. C) : scalar or vector of size n x 1
I              () : Ionic strength (I = 19.924.*S./(1000-1.005*S), S = salinity)
mT_cat   (mol/kg) : Total concentration  of major cations in seawater (Na, K, Mg, Ca, Sr)
mT_an    (mol\kg) : Total concentration  of major anions in seawater  (Cl, SO4, CO3, HCO3, Br, F, B(OH)4)

OUTPUT - gN

gN                : Free ion coefficients for the neutral solutes: NH3, B(OH)3, H3PO4, H2S, SO2, HF, CO2

---- GammaIP.ipynb

This function calculates and stores activity coefficients for all ion pairs that are needed to convert thermodynamic equilibrium constants to stoichiometric constants. 

INPUT - Tc, I, mT_cat, mT_an
Tc      (degr. C) : scalar or vector of size n x 1
I              () : Ionic strength (I = 19.924.*S./(1000-1.005*S), S = salinity)
mT_cat   (mol/kg) : Total concentration  of major cations in seawater (Na, K, Mg, Ca, Sr)
mT_an    (mol\kg) : Total concentration  of major anions in seawater  (Cl, SO4, CO3, HCO3, Br, F, B(OH)4)

OUTPUT - gmjIP, mnCatIP, gdivIP, gtrivIP

gmjIP          () : Neutral ion pairs made from monovalent cations and anions   (NaCl, NaHCO3, NaBr, NaF, 
                                                                                 NaB(OH)4, NaOH, KCl, KBr, KOH)

-              () : -1 ion pairs made from monovalent cation and divalent anion (NaSO4, NaCO3, KSO4)

-              () : +1 ion pairs made from divalent cation and monovalent anion (MgCl, MgHCO3, MgF, MgB(OH)4,              
                                                                                 MgOH, CaCl, CaHCO3, CaF, CaB(OH)4,  
                                                                                 CaOH, SrCl, SrHCO3, SrF, SrB(OH)4, SrOH)
mnCatIP        () : (    HCl,     HSO4,     HCO3,     HHCO3,     HBr,     HF)
-                   (   LiCl,    LiSO4,    LiCO3,    LiHCO3,    LiBr,    LiF)
-                   (   RbCl,    RbSO4,    RbCO3,    RbHCO3,    RbBr,    RbF)
-                   (   CsCl,    CsSO4,    CsCO3,    CsHCO3,    CsBr,    CsF)
-                   (  NH4Cl,   NH4SO4,   NH4CO3,   NH4HCO3,   NH4Br,   NH4F)
-                   (   BaCl,    BaSO4,    BaCO3,    BaHCO3,    BaBr,    BaF)
-                   (Cu(I)Cl, Cu(I)SO4, Cu(I)CO3, Cu(I)HCO3, Cu(I)Br, Cu(I)F)

-                   (    HCl2,     H(SO4)2,     H(CO3)2,     H(HCO3)2,     HBr2,     HF2)
-                   (   LiCl2,    Li(SO4)2,    Li(CO3)2,    Li(HCO3)2,    LiBr2,    LiF2)
-                   (   RbCl2,    Rb(SO4)2,    Rb(CO3)2,    Rb(HCO3)2,    RbBr2,    RbF2)
-                   (   CsCl2,    Cs(SO4)2,    Cs(CO3)2,    Cs(HCO3)2,    CsBr2,    CsF2)
-                   (  NH4Cl2,   NH4(SO4)2,   NH4(CO3)2,   NH4(HCO3)2,   NH4Br2,   NH4F2)
-                   (   BaCl2,    Ba(SO4)2,    Ba(CO3)2,    Ba(HCO3)2,    BaBr2,    BaF2)
-                   (Cu(I)Cl2, Cu(I)(SO4)2, Cu(I)(CO3)2, Cu(I)(HCO3)2, Cu(I)Br2, Cu(I)F2)

-                   (    HCl3,     H(SO4)3,     H(CO3)3,     H(HCO3)3,     HBr3,     HF3)
-                   (   LiCl3,    Li(SO4)3,    Li(CO3)3,    Li(HCO3)3,    LiBr3,    LiF3)
-                   (   RbCl3,    Rb(SO4)3,    Rb(CO3)3,    Rb(HCO3)3,    RbBr3,    RbF3)
-                   (   CsCl3,    Cs(SO4)3,    Cs(CO3)3,    Cs(HCO3)3,    CsBr3,    CsF3)
-                   (  NH4Cl3,   NH4(SO4)3,   NH4(CO3)3,   NH4(HCO3)3,   NH4Br3,   NH4F3)
-                   (   BaCl3,    Ba(SO4)3,    Ba(CO3)3,    Ba(HCO3)3,    BaBr3,    BaF3)
-                   (Cu(I)Cl3, Cu(I)(SO4)3, Cu(I)(CO3)3, Cu(I)(HCO3)3, Cu(I)Br3, Cu(I)F3)

gdivIP         () :

gtrivIP        () :

---- Thermodatabase1.ipynb

List of all the thermodynamic equilibrium constants used in the model as functions of T where possible. We later use activity coefficients determined using our pitzer and empirical values  model to convert to conditional constants.

The log Ks provided are thermodynamic association constants for infinite dilution. 
The log K is assumed to apply at all ionic strengths because they are in terms of activities. 
All ionic strength corrections go into the activity coefficients of the aqueous species.
 
Effectively, PHREEQC automatically uses: 
        log K(conditional) = log K(thermodynamic) * gamma(reactants)/gamma(products) 

If you use conditional constants for seawater, the mass action equations  are written with molalities rather than activities and all ionic strength corrections go into the log K. 

The enthalpy of reaction is used in the Vant Hoff equation to determine the temperature dependence of the equilibrium constant. Internally, all enthalpy calculations are performed in the units of kJ/mol. Default units are kJ/mol. 

For some species, analytical expression for the temperature dependence of log K is defined. 
If defined, the analytical expression takes precedence over log_k and the Vant Hoff equation to determine the temperature dependence of the equilibrium constant.  

---- Kthermo.ipynb

Calculates...
           Major cations - anions
           Minor cations - Major anions
           Minor anions -  Major cations
           DIVALENTS - Major anions

INPUT - Tc

Tc      (degr. C) : scalar or vector of size n x 1

OUTPUT - mjKthermo, mnCat_Kthermo, divKthermo, trivKthermo, mnAn_Kthermo

mjKthermo      () : Association constants for Major Cat - Major An

mnCat_Kthermo  () : Association constants for Minor Cat - Major An

divKthermo     () : Dissociation constants for divalents

trivKthermo    () : Association constants for trivalents

---- Kcond.ipynb

This function converts thermodynamic equilibrium constants (K) into conditional/stoichiometric equilibrium constants (K*)

INPUT - Tc, I, mT_cat, mT_an

Tc      (degr. C) : scalar or vector of size n x 1
I              () : Ionic strength (I = 19.924.*S./(1000-1.005*S), S = salinity)
mT_cat   (mol/kg) : Total concentration  of major cations in seawater (Na, K, Mg, Ca, Sr)
mT_an    (mol\kg) : Total concentration  of major anions in seawater  (Cl, SO4, CO3, HCO3, Br, F, B(OH)4)

OUTPUT - mjKcond, mnCat_Kcond, divKcond, trivKcond

mjKcond        () : Association constants for Major Cat - Major An

mnCat_Kcond    () : Association constants for Minor Cat - Major An

divKcond       () : Dissociation constants for divalents

trivKcond      () : Association constants for trivalents

---- Guess_AlphaV3.ipynb

Allows iteration until convergance by setting all alpha_cat and alpha_an to 1

---- AlphaV3.ipynb

This function calculates the fraction of free ions in a solution of fixed composition and ionic strength (given) needed for correcting free activity coefficients for ion complexation. 
                   
alpha_cat, alpha_an, iter = alphasV3(mT_cat, mT_an)
alpha_cat  .......... [M]f/[M]T
alpha_an   .......... [X]f/[X]T

Assumptions:
1) 1:1 complexes only. Higher order complexes ignored
2) Trivalents form strong hydroxyl complexes
3) Br/Cl- ion pairs weak therefore ignored
                  AlphaCl = 1.0  or [Cl]f = [Cl]T
                  AlphaBr = 1.0  or [Br]f = [Br]T
5) Assume K forms only KSO4- Complex
6) The speciation of minor components of seawater can bemade without making any iterations, since minor components do not affect the speciation of the major componenents 

INPUT - Tc, I, mT_cat, mT_an, OH

Tc      (degr. C) : scalar or vector of size n x 1
I              () : Ionic strength (I = 19.924.*S./(1000-1.005*S), S = salinity)
mT_cat   (mol/kg) : Total concentration  of major cations in seawater (Na, K, Mg, Ca, Sr)
mT_an    (mol\kg) : Total concentration  of major anions in seawater  (Cl, SO4, CO3, HCO3, Br, F, B(OH)4)
OH       (mol\kg) : Total concentration  of OH anion in seawater

OUTPUT - alpha_cat, alpha_an, mF_an, mF_cat, itera

alpha_cat      () : 

alpha_an       () : 

mF_an          () : 

mF_cat         () : 

itera          () : 

---- GammaTV3.ipynb

---- Speciation.ipynb

INPUT - Tc, I, mT_cat, mT_an, OH

Tc      (degr. C) : scalar or vector of size n x 1
I              () : Ionic strength (I = 19.924.*S./(1000-1.005*S), S = salinity)
mT_cat   (mol/kg) : Total concentration  of major cations in seawater (Na, K, Mg, Ca, Sr)
mT_an    (mol\kg) : Total concentration  of major anions in seawater  (Cl, SO4, CO3, HCO3, Br, F, B(OH)4)
OH       (mol\kg) : Total concentration  of OH anion in seawater

OUTPUT - mjcat_spec, mjcat_FREE, mjan_spec, mjan_FREE, mncat_spec, mncat_FREE, div_spec, div_FREE, triv_spec, triv_FREE

mjcat_spec     () : Major cations speciation (NaCl, NaSO4, NaCO3, NaHCO3, NaBr, NaF, NaB(OH)4, NaOH)
                                             ( KCl,  KSO4,  KCO3,  KHCO3,  KBr,  KF,  KB(OH)4,  KOH)
                                             (MgCl, MgSO4, MgCO3, MgHCO3, MgBr, MgF, MgB(OH)4, MgOH)
                                             (CaCl, CaSO4, CaCO3, CaHCO3, CaBr, CaF, CaB(OH)4, CaOH)
                                             (SrCl, SrSO4, SrCO3, SrHCO3, SrBr, SrF, SrB(OH)4, SrOH)

mjcat_FREE     () : Free major cations       (Na, K, Mg, Ca, Sr)

mjcat_spec and mjcat_FREE make up Table 3

mjan_spec      () : Major anion speciation   (    NaCl,     KCl,     MgCl,     CaCl,     SrCl) 
                                             (   NaSO4,    KSO4,    MgSO4,    CaSO4,    SrSO4)
                                             (   NaCO3,    KCO3,    MgCO3,    CaCO3,    SrCO3)
                                             (  NaHCO3,   KHCO3,   MgHCO3,   CaHCO3,   SrHCO3)
                                             (    NaBr,     KBr,     MgBr,     CaBr,     SrBr)
                                             (     NaF,      KF,      MgF,      CaF,      SrF)
                                             (NaB(OH)4, KB(OH)4, MgB(OH)4, CaB(OH)4, SrB(OH)4)
                                             (    NaOH,     KOH,     MgOH,     CaOH,     SrOH)

mjan_FREE      () : Free major anions        (Cl, SO4, CO3, HCO3, Br, F, B(OH)4, OH)

mjan_spec and mjan_FREE make up Table 4

mncat_spec     () : Minor cations speciation (  LiCl,     LiSO4,     LiCO3,     LiHCO3,   LiBr,   LiF,   LiOH)
                                             (  RbCl,     RbSO4,     RbCO3,     RbHCO3,   RbBr,   RbF,   RbOH)
                                             (  CsCl,     CsSO4,     CsCO3,     CsHCO3,   CsBr,   CsF,   CsOH)
                                             ( NH4Cl,    NH4SO4,    NH4CO3,    NH4HCO3,  NH4Br,  NH4F,  NH4OH)
                                             (  BaCl,     BaSO4,     BaCO3,     BaHCO3,   BaBr,   BaF,   BaOH)
                                             (  MnCl,     MnSO4,     MnCO3,     MnHCO3,   MnBr,   MnF,   MnOH)
                                             (  LaCl,     LaSO4,     LaCO3,     LaHCO3,   LaBr,   LaF,   LaOH)

-                                            ( LiCl2,  Li(SO4)2,  Li(CO3)2,  Li(HCO3)2,  LiBr2,  LiF2,  LiOH2)
-                                            ( RbCl2,  Rb(SO4)2,  Rb(CO3)2,  Rb(HCO3)2,  RbBr2,  RbF2,  RbOH2)
                                             ( CsCl2,  Cs(SO4)2,  Cs(CO3)2,  Cs(HCO3)2,  CsBr2,  CsF2,  CsOH2)
                                             (NH4Cl2, NH4(SO4)2, NH4(CO3)2, NH4(HCO3)2, NH4Br2, NH4F2, NH4OH2)
                                             ( BaCl2,  Ba(SO4)2,  Ba(CO3)2,  Ba(HCO3)2,  BaBr2,  BaF2,  BaOH2)
                                             ( MnCl2,  Mn(SO4)2,  Mn(CO3)2,  Mn(HCO3)2,  MnBr2,  MnF2,  MnOH2)
                                             ( LaCl2,  La(SO4)2,  La(CO3)2,  La(HCO3)2,  LaBr2,  LaF2,  LaOH2)

                                             ( LiCl3,  Li(SO4)3,  Li(CO3)3,  Li(HCO3)3,  LiBr3,  LiF3,  LiOH3)
                                             ( RbCl3,  Rb(SO4)3,  Rb(CO3)3,  Rb(HCO3)3,  RbBr3,  RbF3,  RbOH3)
                                             ( CsCl3,  Cs(SO4)3,  Cs(CO3)3,  Cs(HCO3)3,  CsBr3,  CsF3,  CsOH3)
                                             (NH4Cl3, NH4(SO4)3, NH4(CO3)3, NH4(HCO3)3, NH4Br3, NH4F3, NH4OH3)
                                             ( BaCl3,  Ba(SO4)3,  Ba(CO3)3,  Ba(HCO3)3,  BaBr3,  BaF3,  BaOH3)
                                             ( MnCl3,  Mn(SO4)3,  Mn(CO3)3,  Mn(HCO3)3,  MnBr3,  MnF3,  MnOH3)
                                             ( LaCl3,  La(SO4)3,  La(CO3)3,  La(HCO3)3,  LaBr3,  LaF3,  LaOH3)

mncat_FREE     () : Free minor cations       (Li, Rb, Cs, NH4, Ba, Mn, La)

div_spec       () : Divalent speciation      ( Fe(II)Cl,    Fe(II)SO4,    Fe(II)CO3,    Fe(II)HCO3,  Fe(II)F,  Fe(II)OH)
                                             (     CoCl,        CoSO4,        CoCO3,        CoHCO3,      CoF,      CoOH)
                                             (     NiCl,        NiSO4,        NiCO3,        NiHCO3,      NiF,      NiOH)
                                             ( Cu(II)Cl,    Cu(II)SO4,    Cu(II)CO3,    Cu(II)HCO3,  Cu(II)F,  Cu(II)OH)
                                             (     ZnCl,        ZnSO4,        ZnCO3,        ZnHCO3,      ZnF,      ZnOH)
                                             (    UO2Cl,       UO2SO4,       UO2CO3,       UO2HCO3,     UO2F,     UO2OH)
                                             (     BeCl,        BeSO4,        BeCO3,        BeHCO3,      BeF,      BeOH)
                                             (     CdCl,        CdSO4,        CdCO3,        CdHCO3,      CdF,      CdOH)
                                             ( Pb(II)Cl,    Pb(II)SO4,    Pb(II)CO3,    Pb(II)HCO3,  Pb(II)F,  Pb(II)OH)
                                             (  Cu(I)Cl,     Cu(I)SO4,     Cu(I)CO3,     Cu(I)HCO3,   Cu(I)F,   Cu(I)OH)

-                                            (Fe(II)Cl2, Fe(II)(SO4)2, Fe(II)(CO3)2, Fe(II)(HCO3)2, Fe(II)F2, Fe(II)OH2)
                                             (    CoCl2,     Co(SO4)2,     Co(CO3)2,     Co(HCO3)2,     CoF2,     CoOH2)
                                             (    NiCl2,     Ni(SO4)2,     Ni(CO3)2,     Ni(HCO3)2,     NiF2,     NiOH2)
                                             (Cu(II)Cl2, Cu(II)(SO4)2, Cu(II)(CO3)2, Cu(II)(HCO3)2, Cu(II)F2, Cu(II)OH2)
                                             (    ZnCl2,     Zn(SO4)2,     Zn(CO3)2,     Zn(HCO3)2,     ZnF2,     ZnOH2)
                                             (   UO2Cl2,    UO2(SO4)2,    UO2(CO3)2,    UO2(HCO3)2,    UO2F2,    UO2OH2)
                                             (    BeCl2,     Be(SO4)2,     Be(CO3)2,     Be(HCO3)2,     BeF2,     BeOH2)
                                             (    CdCl2,     Cd(SO4)2,     Cd(CO3)2,     Cd(HCO3)2,     CdF2,     CdOH2)
                                             (Pb(II)Cl2, Pb(II)(SO4)2, Pb(II)(CO3)2, Pb(II)(HCO3)2, Pb(II)F2, Pb(II)OH2)
                                             ( Cu(I)Cl2,  Cu(I)(SO4)2,  Cu(I)(CO3)2,  Cu(I)(HCO3)2,  Cu(I)F2,  Cu(I)OH2)

                                             (Fe(II)Cl3, Fe(II)(SO4)3, Fe(II)(CO3)3, Fe(II)(HCO3)3, Fe(II)F3, Fe(II)OH3)
                                             (    CoCl3,     Co(SO4)3,     Co(CO3)3,     Co(HCO3)3,     CoF3,     CoOH3)
                                             (    NiCl3,     Ni(SO4)3,     Ni(CO3)3,     Ni(HCO3)3,     NiF3,     NiOH3)
                                             (Cu(II)Cl3, Cu(II)(SO4)3, Cu(II)(CO3)3, Cu(II)(HCO3)3, Cu(II)F3, Cu(II)OH3)
                                             (    ZnCl3,     Zn(SO4)3,     Zn(CO3)3,     Zn(HCO3)3,     ZnF3,     ZnOH3)
                                             (   UO2Cl3,    UO2(SO4)3,    UO2(CO3)3,    UO2(HCO3)3,    UO2F3,    UO2OH3)
                                             (    BeCl3,     Be(SO4)3,     Be(CO3)3,     Be(HCO3)3,     BeF3,     BeOH3)
                                             (    CdCl3,     Cd(SO4)3,     Cd(CO3)3,     Cd(HCO3)3,     CdF3,     CdOH3)
                                             (Pb(II)Cl3, Pb(II)(SO4)3, Pb(II)(CO3)3, Pb(II)(HCO3)3, Pb(II)F3, Pb(II)OH3)
                                             ( Cu(I)Cl3,  Cu(I)(SO4)3,  Cu(I)(CO3)3,  Cu(I)(HCO3)3,  Cu(I)F3,  Cu(I)OH3)

                                             (Fe(II)Cl4, Fe(II)(SO4)4, Fe(II)(CO3)4, Fe(II)(HCO3)4, Fe(II)F4, Fe(II)OH4)
                                             (    CoCl4,     Co(SO4)4,     Co(CO3)4,     Co(HCO3)4,     CoF4,     CoOH4)
                                             (    NiCl4,     Ni(SO4)4,     Ni(CO3)4,     Ni(HCO3)4,     NiF4,     NiOH4)
                                             (Cu(II)Cl4, Cu(II)(SO4)4, Cu(II)(CO3)4, Cu(II)(HCO3)4, Cu(II)F4, Cu(II)OH4)
                                             (    ZnCl4,     Zn(SO4)4,     Zn(CO3)4,     Zn(HCO3)4,     ZnF4,     ZnOH4)
                                             (   UO2Cl4,    UO2(SO4)4,    UO2(CO3)4,    UO2(HCO3)4,    UO2F4,    UO2OH4)
                                             (    BeCl4,     Be(SO4)4,     Be(CO3)4,     Be(HCO3)4,     BeF4,     BeOH4)
                                             (    CdCl4,     Cd(SO4)4,     Cd(CO3)4,     Cd(HCO3)4,     CdF4,     CdOH4)
                                             (Pb(II)Cl4, Pb(II)(SO4)4, Pb(II)(CO3)4, Pb(II)(HCO3)4, Pb(II)F4, Pb(II)OH4)
                                             ( Cu(I)Cl4,  Cu(I)(SO4)4,  Cu(I)(CO3)4,  Cu(I)(HCO3)4,  Cu(I)F4,  Cu(I)OH4)

div_FREE       () : Free divalent mn anions  (Fe(II), Co, Ni, Cu(II), Zn, UO2, Be, Cd, Pb(II), Cu(I))

triv_spec      () : Divalent speciation      (      LaCl,         LaSO4,         LaCO3,         LaHCO3,       LaF,       LaOH]
                                             (      CeCl,         CeSO4,         CeCO3,         CeHCO3,       CeF,       CeOH]
                                             (      PrCl,         PrSO4,         PrCO3,         PrHCO3,       PrF,       PrOH]
                                             (      NdCl,         NdSO4,         NdCO3,         NdHCO3,       NdF,       NdOH]
                                             (      PmCl,         PmSO4,         PmCO3,         PmHCO3,       PmF,       PmOH]
                                             (      SmCl,         SmSO4,         SmCO3,         SmHCO3,       SmF,       SmOH]
                                             (      EuCl,         EuSO4,         EuCO3,         EuHCO3,       EuF,       EuOH]
                                             (      GdCl,         GdSO4,         GdCO3,         GdHCO3,       GdF,       GdOH]
                                             (      TbCl,         TbSO4,         TbCO3,         TbHCO3,       TbF,       TbOH]
                                             (      DyCl,         DySO4,         DyCO3,         DyHCO3,       DyF,       DyOH]
                                             (      HoCl,         HoSO4,         HoCO3,         HoHCO3,       HoF,       HoOH]
                                             (      ErCl,         ErSO4,         ErCO3,         ErHCO3,       ErF,       ErOH]
                                             (      TmCl,         TmSO4,         TmCO3,         TmHCO3,       TmF,       TmOH]
                                             (      YbCl,         YbSO4,         YbCO3,         YbHCO3,       YbF,       YbOH]
                                             (      LuCl,         LuSO4,         LuCO3,         LuHCO3,       LuF,       LuOH]
                                             (       YCl,          YSO4,          YCO3,          YHCO3,        YF,        YOH]
                                             (      AlCl,         AlSO4,         AlCO3,         AlHCO3,       AlF,       AlOH]
                                             (      GaCl,         GaSO4,         GaCO3,         GaHCO3,       GaF,       GaOH]
                                             (      InCl,         InSO4,         InCO3,         InHCO3,       InF,       InOH]
                                             ( Fe(III)Cl,    Fe(III)SO4,    Fe(III)CO3,    Fe(III)HCO3,  Fe(III)F,  Fe(III)OH]

-                                            (     LaCl2,      La(SO4)2,      La(CO3)2,      La(HCO3)2,      LaF2,      LaOH2]
                                             (     CeCl2,      Ce(SO4)2,      Ce(CO3)2,      Ce(HCO3)2,      CeF2,      CeOH2]
                                             (     PrCl2,      Pr(SO4)2,      Pr(CO3)2,      Pr(HCO3)2,      PrF2,      PrOH2]
                                             (     NdCl2,      Nd(SO4)2,      Nd(CO3)2,      Nd(HCO3)2,      NdF2,      NdOH2]
                                             (     PmCl2,      Pm(SO4)2,      Pm(CO3)2,      Pm(HCO3)2,      PmF2,      PmOH2]
                                             (     SmCl2,      Sm(SO4)2,      Sm(CO3)2,      Sm(HCO3)2,      SmF2,      SmOH2]
                                             (     EuCl2,      Eu(SO4)2,      Eu(CO3)2,      Eu(HCO3)2,      EuF2,      EuOH2]
                                             (     GdCl2,      Gd(SO4)2,      Gd(CO3)2,      Gd(HCO3)2,      GdF2,      GdOH2]
                                             (     TbCl2,      Tb(SO4)2,      Tb(CO3)2,      Tb(HCO3)2,      TbF2,      TbOH2]
                                             (     DyCl2,      Dy(SO4)2,      Dy(CO3)2,      Dy(HCO3)2,      DyF2,      DyOH2]
                                             (     HoCl2,      Ho(SO4)2,      Ho(CO3)2,      Ho(HCO3)2,      HoF2,      HoOH2]
                                             (     ErCl2,      Er(SO4)2,      Er(CO3)2,      Er(HCO3)2,      ErF2,      ErOH2]
                                             (     TmCl2,      Tm(SO4)2,      Tm(CO3)2,      Tm(HCO3)2,      TmF2,      TmOH2]
                                             (     YbCl2,      Yb(SO4)2,      Yb(CO3)2,      Yb(HCO3)2,      YbF2,      YbOH2]
                                             (     LuCl2,      Lu(SO4)2,      Lu(CO3)2,      Lu(HCO3)2,      LuF2,      LuOH2]
                                             (      YCl2,       Y(SO4)2,       Y(CO3)2,       Y(HCO3)2,       YF2,       YOH2]
                                             (     AlCl2,      Al(SO4)2,      Al(CO3)2,      Al(HCO3)2,      AlF2,      AlOH2]
                                             (     GaCl2,      Ga(SO4)2,      Ga(CO3)2,      Ga(HCO3)2,      GaF2,      GaOH2]
                                             (     InCl2,      In(SO4)2,      In(CO3)2,      In(HCO3)2,      InF2,      InOH2]
                                             (Fe(III)Cl2, Fe(III)(SO4)2, Fe(III)(CO3)2, Fe(III)(HCO3)2, Fe(III)F2, Fe(III)OH2]

                                             (     LaCl3,      La(SO4)3,      La(CO3)3,      La(HCO3)3,      LaF3,      LaOH3]
                                             (     CeCl3,      Ce(SO4)3,      Ce(CO3)3,      Ce(HCO3)3,      CeF3,      CeOH3]
                                             (     PrCl3,      Pr(SO4)3,      Pr(CO3)3,      Pr(HCO3)3,      PrF3,      PrOH3]
                                             (     NdCl3,      Nd(SO4)3,      Nd(CO3)3,      Nd(HCO3)3,      NdF3,      NdOH3]
                                             (     PmCl3,      Pm(SO4)3,      Pm(CO3)3,      Pm(HCO3)3,      PmF3,      PmOH3]
                                             (     SmCl3,      Sm(SO4)3,      Sm(CO3)3,      Sm(HCO3)3,      SmF3,      SmOH3]
                                             (     EuCl3,      Eu(SO4)3,      Eu(CO3)3,      Eu(HCO3)3,      EuF3,      EuOH3]
                                             (     GdCl3,      Gd(SO4)3,      Gd(CO3)3,      Gd(HCO3)3,      GdF3,      GdOH3]
                                             (     TbCl3,      Tb(SO4)3,      Tb(CO3)3,      Tb(HCO3)3,      TbF3,      TbOH3]
                                             (     DyCl3,      Dy(SO4)3,      Dy(CO3)3,      Dy(HCO3)3,      DyF3,      DyOH3]
                                             (     HoCl3,      Ho(SO4)3,      Ho(CO3)3,      Ho(HCO3)3,      HoF3,      HoOH3]
                                             (     ErCl3,      Er(SO4)3,      Er(CO3)3,      Er(HCO3)3,      ErF3,      ErOH3]
                                             (     TmCl3,      Tm(SO4)3,      Tm(CO3)3,      Tm(HCO3)3,      TmF3,      TmOH3]
                                             (     YbCl3,      Yb(SO4)3,      Yb(CO3)3,      Yb(HCO3)3,      YbF3,      YbOH3]
                                             (     LuCl3,      Lu(SO4)3,      Lu(CO3)3,      Lu(HCO3)3,      LuF3,      LuOH3]
                                             (      YCl3,       Y(SO4)3,       Y(CO3)3,       Y(HCO3)3,       YF3,       YOH3]
                                             (     AlCl3,      Al(SO4)3,      Al(CO3)3,      Al(HCO3)3,      AlF3,      AlOH3]
                                             (     GaCl3,      Ga(SO4)3,      Ga(CO3)3,      Ga(HCO3)3,      GaF3,      GaOH3]
                                             (     InCl3,      In(SO4)3,      In(CO3)3,      In(HCO3)3,      InF3,      InOH3]
                                             (Fe(III)Cl3, Fe(III)(SO4)3, Fe(III)(CO3)3, Fe(III)(HCO3)3, Fe(III)F3, Fe(III)OH3]

                                             (     LaCl4,      La(SO4)4,      La(CO3)4,      La(HCO3)4,      LaF4,      LaOH4]
                                             (     CeCl4,      Ce(SO4)4,      Ce(CO3)4,      Ce(HCO3)4,      CeF4,      CeOH4]
                                             (     PrCl4,      Pr(SO4)4,      Pr(CO3)4,      Pr(HCO3)4,      PrF4,      PrOH4]
                                             (     NdCl4,      Nd(SO4)4,      Nd(CO3)4,      Nd(HCO3)4,      NdF4,      NdOH4]
                                             (     PmCl4,      Pm(SO4)4,      Pm(CO3)4,      Pm(HCO3)4,      PmF4,      PmOH4]
                                             (     SmCl4,      Sm(SO4)4,      Sm(CO3)4,      Sm(HCO3)4,      SmF4,      SmOH4]
                                             (     EuCl4,      Eu(SO4)4,      Eu(CO3)4,      Eu(HCO3)4,      EuF4,      EuOH4]
                                             (     GdCl4,      Gd(SO4)4,      Gd(CO3)4,      Gd(HCO3)4,      GdF4,      GdOH4]
                                             (     TbCl4,      Tb(SO4)4,      Tb(CO3)4,      Tb(HCO3)4,      TbF4,      TbOH4]
                                             (     DyCl4,      Dy(SO4)4,      Dy(CO3)4,      Dy(HCO3)4,      DyF4,      DyOH4]
                                             (     HoCl4,      Ho(SO4)4,      Ho(CO3)4,      Ho(HCO3)4,      HoF4,      HoOH4]
                                             (     ErCl4,      Er(SO4)4,      Er(CO3)4,      Er(HCO3)4,      ErF4,      ErOH4]
                                             (     TmCl4,      Tm(SO4)4,      Tm(CO3)4,      Tm(HCO3)4,      TmF4,      TmOH4]
                                             (     YbCl4,      Yb(SO4)4,      Yb(CO3)4,      Yb(HCO3)4,      YbF4,      YbOH4]
                                             (     LuCl4,      Lu(SO4)4,      Lu(CO3)4,      Lu(HCO3)4,      LuF4,      LuOH4]
                                             (      YCl4,       Y(SO4)4,       Y(CO3)4,       Y(HCO3)4,       YF4,       YOH4]
                                             (     AlCl4,      Al(SO4)4,      Al(CO3)4,      Al(HCO3)4,      AlF4,      AlOH4]
                                             (     GaCl4,      Ga(SO4)4,      Ga(CO3)4,      Ga(HCO3)4,      GaF4,      GaOH4]
                                             (     InCl4,      In(SO4)4,      In(CO3)4,      In(HCO3)4,      InF4,      InOH4]
                                             (Fe(III)Cl4, Fe(III)(SO4)4, Fe(III)(CO3)4, Fe(III)(HCO3)3, Fe(III)F4, Fe(III)OH4]

triv_FREE      () : Free trivalent mn anions (La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Y, Al, Ga, In, Fe(III))

---- PyCO2SYS.ipynb

Parameter info - generate all combinations of marine carbonate system parameters

defpars     = (Total Alkalinity, Total CO2, Input pH, Pressure CO2 input, Fugacity of CO2 Input, Input CO3, Input HCO3, Input CO2)
defpartypes = (               1,         2,        3,                  4,                     5,         6,          7,         8)

---- Literature_Graphs.ipynb

Creation of graphs using data from Cantrell and Byrne 1987, Stanley and Byrne 1990 and:
Earth Syst. Sci. Data, 8, 325–340, 2016 www.earth-syst-sci-data.net/8/325/2016/ doi:10.5194/essd-8-325-2016
