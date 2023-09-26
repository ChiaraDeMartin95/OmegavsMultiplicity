#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TROOT.h>

#include "THStack.h"
#include "TCanvas.h"
#include <TLatex.h>

#include <iostream>
#include <fstream>

TFile* fileData;
TFile* fileOut;

TH3F* hPtCascadePlusTrueAssoc;
TH3F* hPtCascadeMinusTrueAssoc;
TH3F* hPtCascadeSumAssoc;

TH3F* hPtCascadePlusTrue;
TH3F* hPtCascadeMinusTrue;
TH3F* hPtCascadeSumTrue;

TH3F* hPtCascadePlusTrueRec;
TH3F* hPtCascadeMinusTrueRec;
TH3F* hPtCascadeRecSum;

const Double_t massOmega = 1.67245;
Double_t gBinwidth = 1;

// Binning used in MB analysis
const Int_t nPtBinsMB = 7;
Double_t xBinsMB[nPtBinsMB+1] = {1.00, 1.40, 1.80, 2.30, 2.80, 3.30, 3.80, 4.80};

const Int_t nBinsChiara = 19;
Double_t xbinsChiara[nBinsChiara+1] = {0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 8.00};

const Int_t nBinsUpasana = 3;
Double_t xbinsUpasana[nBinsUpasana+1] = {0.6, 1.4, 2.0, 3.0};


