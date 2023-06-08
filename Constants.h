const Int_t numMult = 5; //10 // number of multiplicity sub-intervals
const Int_t numPart = 7;  // number of particles
const Int_t numFiles = 8;
const Int_t numTopoVar = 13;
const Int_t numChoice = 7;    // mean, sigma, purity, yield, significance, efficiency, yieldCorr
const Int_t numParticles = 3; // part, apart, sum

//Int_t MultiplicityPerc[numMult + 1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 60, 100};
Int_t MultiplicityPerc[numMult + 1] = {0, 10, 20, 30, 40, 100};
Int_t ColorMult[] = {634, 628, 807, 797, 815, 418, 429, 867, 856, 601, 1};
Float_t SizeMult[] = {2, 2, 2.8, 2.5, 2.8, 2, 2.8, 2.5, 2, 2, 2};
Float_t SizeMultRatio[] = {1, 1, 1.8, 1.5, 1.8, 1, 1.8, 1.5, 1, 1, 1};
Int_t MarkerMult[] = {20, 21, 33, 34, 29, 24, 27, 28, 25, 25, 25};
Float_t ScaleFactor[] = {4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 1};
TString sScaleFactor[] = {" (x2^{12})", " (x2^{11})", " (x2^{10})", " (x2^{9})", " (x2^{8})", " (x2^{7})",
                                     " (x2^{6})", " (x2^{5})", " (x2^{4})", " (x2^{3})", ""};

Int_t Color[] = {634, 628, 797, 815, 418, 429, 867, 601, 1};
Int_t ColorPart[] = {634, kBlue + 2};

TString SMultType[] = {"", "FT0M", "FV0A"};

const Float_t massParticle[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
TString Spart[numPart] = {"K0s", "Lambda", "AntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPos"};
TString SpartType[numPart] = {"K0s", "Lambda", "Lambda", "Xi", "Xi", "Omega", "Omega"};
TString NamePart[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};

TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};

TString SIsBkgParab[2] = {"_BkgRetta", "_BkgParab"};

TString TypeHisto[numChoice] = {"Mean", "Sigma", "Purity", "Yield", "Significance", "Efficiency", "YieldCorr"};
TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "Yield/#sigma_{Yield}", "Efficiency x acceptance", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}"};

TString TopoVar[numTopoVar] = {"casccospa", "dcacascdau", "dcabachtopv", "dcapostopv", "dcanegtopv", "lambdamasswin", "rejcomp", "nsigmatpcKa", "cascradius", "v0radius", "dcav0dau", "v0cospa", "casclifetime"};
TString TopoVarSigned[numTopoVar] = {"casccospa > ", "dcacascdau < ", "dcabachtopv > ", "dcapostopv > ", "dcanegtopv > ", "|mV0 - mLambda| < ", "|mCasc-mXi| > ", "|n#sigma_{Ka}^{TPC}| < ", "rXi > ", "rV0 > ", "dcav0dau < ", "v0cospa > ", "c#tau < "};

TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";
TString TitleYYield = "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}";