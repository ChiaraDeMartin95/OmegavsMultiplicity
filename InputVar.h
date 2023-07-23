TString ExtrSysPath = "_Train100720";  
//old ones: "_Sel3July"
TString Extryear = "LHC22o_pass4_Train99659";  
//old ones: "LHC22o_pass4_Train89684", "LHC22o_pass3_Train75538", "LHC22r_pass3_Train67853"
TString ExtrSPathInEff = "eff_LHC22o_pass4_Train100720.root";  
//old ones: "eff_LHC22o_pass4_Sel23June.root", "eff6June.root"
TString PathInSyst = "SystematicErrors/SpectraSystematics-Omega-13TeV-FT0M-0-100_21July.root";

Bool_t ExtrisBkgParab = 1;
Bool_t ExtrUseTwoGauss = 1;
Int_t ExtrevFlag = 1;  // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
Int_t ExtrMultType = 1;  // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A