const Int_t numMult = 5; // number of multiplicity sub-intervals
const Int_t numPart = 7; // number of particles
const Int_t numFiles = 8;
const Int_t numTopoVar = 10;
const Int_t numChoice = 5; //mean, sigma, purity, yield, significance

Int_t MultiplicityPerc[numMult + 1] = {0, 10, 20, 40, 60, 100};
Int_t Color[] = {634, 628, 797,815,418, 429, 867,601,1};

const Float_t massParticle[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
TString Spart[numPart] = {"K0s", "Lambda", "AntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPos"};
TString SpartType[numPart] = {"K0s", "Lambda", "Lambda", "Xi", "Xi", "Omega", "Omega"};
TString NamePart[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};

TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};

TString SIsBkgParab[2] = {"_BkgRetta", "_BkgParab"};

TString TypeHisto[numChoice] = {"Mean", "Sigma", "Purity", "Yield", "Significance"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "Yield/#sigma_{Yield}"};

TString TopoVar[numTopoVar] = {"casccospa", "dcacascdau", "dcabachtopv", "dcapostopv", "dcanegtopv", "lambdamasswin", "rejcomp", "nsigmatpcKa", "cascradius"};
TString TopoVarSigned[numTopoVar] = {"casccospa > ", "dcacascdau < ", "dcabachtopv > ", "dcapostopv > ", "dcanegtopv > ", "|mV0 - mLambda| < ", "|mCasc-mXi| > ", "|n#sigma_{Ka}^{TPC}| < ", "rXi > "};

TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";