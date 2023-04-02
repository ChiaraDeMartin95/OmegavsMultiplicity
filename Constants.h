const Int_t numMult = 5; // number of multiplicity sub-intervals
const Int_t numPart = 7; // number of particles

Int_t MultiplicityPerc[numMult + 1] = {0, 10, 20, 40, 60, 100};

const Float_t massParticle[numPart] = {0.497611, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
TString Spart[numPart] = {"K0s", "Lambda", "AntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus"};
TString SpartType[numPart] = {"K0s", "Lambda", "Lambda", "Xi", "Xi", "Omega", "Omega"};
TString NamePart[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
TString IsOneOrTwoGauss[2] = {"_OneGaussFit", ""};

TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";