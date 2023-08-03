string ExtrSysPathOmega = "_Train110357";
//Old ones: "_Train109354";
string ExtrSysPathXi = "_Train109827";
// old ones: "_Train100720", "_Sel3July"
string Extryear = "LHC22o_pass4_MinBias_Train108123";
// old ones: "LHC22o_pass4_Train99659", "LHC22o_pass4_Train89684", "LHC22o_pass3_Train75538", "LHC22r_pass3_Train67853"
TString ExtrSPathInEffOmega = "effChiaraOmega_inelgt0_lhc23e1b_3aug.root";
TString ExtrSPathInEffXi = "effChiaraXi_inelgt0_lhc23e1b_3aug.root";
// old ones: "eff_LHC22o_pass4_Sel23June.root", "eff6June.root"
TString PathInSyst = "SystematicErrors/SpectraSystematics-Omega-13TeV-FT0M-0-100_21July.root";

Int_t ExtrParticle = 8;
Bool_t ExtrisBkgParab = 1;
Bool_t ExtrUseTwoGauss = 1;
Int_t ExtrevFlag = 1;   // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
Int_t ExtrMultType = 1; // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A