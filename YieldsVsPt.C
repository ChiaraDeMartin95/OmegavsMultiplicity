#include <Riostream.h>
#include <string>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TCutG.h>
#include "TFitResult.h"
#include "TLegend.h"
#include "/Users/mbp-cdm-01/Desktop/AssegnoRicerca/Run3Analyses/OmegavsMult/Constants.h"
#include "/Users/mbp-cdm-01/Desktop/AssegnoRicerca/Run3Analyses/OmegavsMult/InputVar.h"

void StyleCanvas(TCanvas *canvas, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(LMargin);
  canvas->SetRightMargin(RMargin);
  canvas->SetTopMargin(TMargin);
  canvas->SetBottomMargin(BMargin);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  // gStyle->SetPalette(55, 0);
}

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX,
                TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)
    histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);
}

const Float_t UpperLimitLSBOmega = 1.655; // upper limit of fit of left sidebands for omega
const Float_t LowerLimitRSBOmega = 1.689; // lower limit of fit of right sidebands for omega
const Float_t UpperLimitLSBXi = 1.302;    // upper limit of fit of left sidebands for Xi
const Float_t LowerLimitRSBXi = 1.34;     // lower limit of fit of right sidebands for Xi

Bool_t reject = 1;
Double_t fparab(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[3] == 0)
  {
    LimInf = 0.474;
    LimSup = 0.520;
  }
  else if (par[3] == 3 || par[3] == 4 || par[3] == 5)
  {
    LimInf = UpperLimitLSBXi; // 1.31
    LimSup = LowerLimitRSBXi; // 1.335
  }
  else if (par[3] == 6 || par[3] == 7 || par[3] == 8)
  {
    LimInf = UpperLimitLSBOmega; // 1.66
    LimSup = LowerLimitRSBOmega; // 1.692
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

Double_t fretta(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[2] == 0)
  {
    LimInf = 0.47;
    LimSup = 0.530;
  }
  else if (par[2] == 3 || par[2] == 4 || par[2] == 5)
  {
    LimInf = UpperLimitLSBXi; // 1.31
    LimSup = LowerLimitRSBXi; // 1.335
  }
  else if (par[2] == 6 || par[2] == 7 || par[2] == 8)
  {
    LimInf = UpperLimitLSBOmega; // 1.65
    LimSup = LowerLimitRSBOmega; // 1.69
  }
  if (reject && x[0] > LimInf && x[0] < LimSup)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0];
}

TString titlePt = "p_{T} (GeV/c)";
TString titleYield = "1/N_{ev} dN/dp_{T}";

TString TitleInvMass[numPart] = {"(#pi^{+}, #pi^{-})", "(p, #pi^{-})", "(#bar{p}, #pi^{-})", "(#Lambda, #pi^{-})", "(#bar{Lambda}, #pi^{+})", "(#Lambda, #pi)", "(#Lambda, K^{-})", "(#bar{#Lambda}, K^{+})", "(#Lambda, K)"};
TString SInvMass = "invariant mass (GeV/#it{c}^{2})";
TString namehisto[numPart] = {"h3dMassK0Short", "", "", "hCascMinusInvMassvsPt", "hCascPlusInvMassvsPt", "hCascMinusInvMassvsPt", "hCascMinusInvMassvsPt", "hCascPlusInvMassvsPt", "hCascMinusInvMassvsPt"};

// fit ranges
Float_t min_range_signal[numPart] = {0.46, 1.105, 1.105, 1.3, 1.3, 1.3, 1.65, 1.65, 1.65}; // gauss fit range
Float_t max_range_signal[numPart] = {0.535, 1.125, 1.125, 1.335, 1.335, 1.335, 1.69, 1.69, 1.69};
Float_t liminf[numPart] = {0.45, 1.1153, 1.1153, 1.29, 1.29, 1.29, 1.63, 1.63, 1.63}; // bkg and total fit range
Float_t limsup[numPart] = {0.545, 1.1168, 1.1168, 1.352, 1.352, 1.352, 1.71, 1.71, 1.71};

// visualisation ranges
Float_t LowMassRange[numPart] = {0.48, 1.09, 1.09, 1.31, 1.31, 1.31, 1.655, 1.655, 1.655}; // range to compute approximate yield (signal + bkg)
Float_t UpMassRange[numPart] = {0.51, 1.14, 1.14, 1.33, 1.33, 1.33, 1.685, 1.685, 1.685};
Float_t gaussDisplayRangeLow[numPart] = {0.42, 1.09, 1.09, 1.29, 1.29, 1.29, 1.63, 1.63, 1.63}; // display range of gauss functions (from total fit)
Float_t gaussDisplayRangeUp[numPart] = {0.57, 1.14, 1.14, 1.35, 1.35, 1.35, 1.71, 1.71, 1.71};
Float_t bkgDisplayRangeLow[numPart] = {0.42, 1.09, 1.09, 1.28, 1.28, 1.28, 1.626, 1.626, 1.626}; // display range of bkg function (from total fit)
Float_t bkgDisplayRangeUp[numPart] = {0.57, 1.14, 1.14, 1.36, 1.36, 1.36, 1.72, 1.72, 1.72};
Float_t histoMassRangeLow[numPart] = {0, 0, 0, 1.28, 1.28, 1.28, 1.626, 1.626, 1.626}; // display range of mass histograms
Float_t histoMassRangeUp[numPart] = {0, 0, 0, 1.36, 1.36, 1.36, 1.72, 1.72, 1.72};

void YieldsVsPt(
    Int_t part = ExtrParticle,
    Bool_t isMB = 1,
    Int_t mul = 0,
    Bool_t isSysStudy = 1,
    string SysPath = "",
    TString StringPathInEvt = "",
    Bool_t isYAxisMassZoomed = 0,
    Int_t MassRebin = 2,
    Int_t MultType = ExtrMultType, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
    TString year = Extryear,
    TString SPathIn = "OutputFilesCascPPTask/LHC22o_pass4_MinBias/AnalysisResults",
    TString SPathInEvt = "AnalysisResultsCascQATask/LHC22o_pass4_MinBias/AnalysisResultsEvts",
    TString OutputDir = "Yields",
    Bool_t UseTwoGauss = ExtrUseTwoGauss,
    Bool_t isBkgParab = ExtrisBkgParab,
    Bool_t isMeanFixedPDG = 0,
    Float_t sigmacentral = 4.2,
    Int_t evFlag = ExtrevFlag, // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1,
    Int_t SysSigExtr = 0       // 0: default, 1: pol1 bkg, 2: tighter range + pol2, 3: tighter range + pol1
)
{

  if (part == 3 || part == 4 || part == 5)
    SysPath = ExtrSysPathXi;
  else if (part == 6 || part == 7 || part == 8)
    SysPath = ExtrSysPathOmega;

  Float_t UpperLimitLSB = 0;
  Float_t LowerLimitRSB = 0;
  if (part == 3 || part == 4 || part == 5)
  {
    UpperLimitLSB = UpperLimitLSBXi;
    LowerLimitRSB = LowerLimitRSBXi;
  }
  else if (part == 6 || part == 7 || part == 8)
  {
    UpperLimitLSB = UpperLimitLSBOmega;
    LowerLimitRSB = LowerLimitRSBOmega;
  }

  if (SysSigExtr == 1)
    isBkgParab = 0;

  if (mul > numMult)
  {
    cout << "Multiplciity out of range" << endl;
    return;
  }
  if (MultType == 0 && (part == 5 || part == 8))
  { // pos + neg cascades
    cout << "No backward compatibility for this case" << endl;
    return;
  }
  if (part < 3)
  {
    cout << "this value of part is not specified, choose 3 - 4 (Xi) or 5 - 6  (Omega) " << endl;
    return;
  }

  SPathIn += "_" + year + "_" + SpartType[part];
  if (isSysStudy)
    SPathIn += SysPath;
  SPathIn += ".root";

  SPathInEvt += "_" + year;
  if (SysPath.find("_Train109827") != string::npos)
    SPathInEvt += "_AllRuns";
  SPathInEvt += StringPathInEvt;
  // SPathInEvt += "_run528531";
  SPathInEvt += ".root";
  TFile *filein = new TFile(SPathIn, "");
  if (!filein)
  {
    cout << "FileIn not available" << endl;
    return;
  }

  TFile *fileinEvt = new TFile(SPathInEvt, "");
  if (!fileinEvt)
  {
    cout << "FileInEvt not available" << endl;
    return;
  }

  TDirectoryFile *dir;
  TDirectoryFile *dirEvt;
  TDirectoryFile *dirCasc;
  TH3F *h3;
  TH3F *h3Plus;
  TH2F *h2;
  TH2F *h2Bis;
  TH1F *hEvents;
  TH3F *hEventsFT0M3D;
  TH2F *hEventsFT0M2D;
  TH2F *hEventsFV0A2D;
  TH1F *hEventsFT0M;
  TH1F *hEventsFV0A;
  Int_t MultLowBin = 0;
  Int_t MultUpBin = 0;

  // if (SysPath == "_Train100720")
  if (SysPath.find("_Train100720") != string::npos)
  {
    if (part == 3 || part == 4 || part == 5)
      dir = (TFile *)filein->Get("lf-cascpostprocessing_id4873");
    else if (part == 6 || part == 7 || part == 8)
      dir = (TFile *)filein->Get("lf-cascpostprocessing_id4477");
  }
  else if (SysPath.find("_Train109354") != string::npos)
  {
    dir = (TFile *)filein->Get("lf-cascpostprocessing_id4477");
  }
  else
    dir = (TFile *)filein->Get("lf-cascpostprocessing");
  if (!dir)
  {
    cout << "dir not available" << endl;
    return;
  }

  TH1F *histoCandidateSelections = (TH1F *)dir->Get("CascadeSelectionSummary");
  if (!histoCandidateSelections)
  {
    cout << "histoCandidateSelections not found" << endl;
    return;
  }

  if (MultType == 0)
  {
    h2 = (TH2F *)dir->Get(namehisto[part]);
    if (!h2)
    {
      cout << "h2 cascade not avilable " << endl;
      return;
    }
  }
  else
  {
    if (MultType == 1)
    {
      h3 = (TH3F *)dir->Get(namehisto[part] + "_FT0M");
      h3Plus = (TH3F *)dir->Get("hCascPlusInvMassvsPt_FT0M");
    }
    else if (MultType == 2)
    {
      h3 = (TH3F *)dir->Get(namehisto[part] + "_FV0A");
      h3Plus = (TH3F *)dir->Get("hCascPlusInvMassvsPt_FV0A");
    }
    if (!h3)
    {
      cout << "h3 cascade not avilable " << endl;
      return;
    }
    if (!h3Plus)
    {
      cout << "h3 positive cascade not avilable " << endl;
      return;
    }
    if (part == 5 || part == 8)
      h3->Add(h3Plus);
    if (isMB)
    {
      MultLowBin = h3->GetXaxis()->FindBin(100 + 0.001);
      MultUpBin = h3->GetXaxis()->FindBin(100 - 0.001);
    }
    else
    {
      MultLowBin = h3->GetXaxis()->FindBin(MultiplicityPerc[mul] + 0.001);
      MultUpBin = h3->GetXaxis()->FindBin(MultiplicityPerc[mul + 1] - 0.001);
    }
    h3->GetXaxis()->SetRange(MultLowBin, MultUpBin);
    h2 = (TH2F *)h3->Project3D("zyoe");
  }

  Double_t NEvents = 0;
  dirEvt = (TDirectoryFile *)fileinEvt->Get("lf-cascqaanalysis");
  hEvents = (TH1F *)dirEvt->Get("hNEvents");
  if (!hEvents)
  {
    cout << "hEvents not available " << endl;
    return;
  }

  if (SysPath.find("_Train100720") != string::npos)
  {
    hEventsFT0M2D = (TH2F *)dirEvt->Get("hCentFT0M");
    if (!hEventsFT0M2D)
    {
      cout << "hCentFT0M not available " << endl;
      return;
    }
  }
  else if (SysPath.find("_run528531") != string::npos || SysPath.find("_Train109") != string::npos || SysPath.find("_Train110") != string::npos)
  {
    hEventsFT0M3D = (TH3F *)dirEvt->Get("hFT0Mglobal");
    if (!hEventsFT0M3D)
    {
      cout << "hFT0Mglobal not available " << endl;
      return;
    }
    hEventsFT0M2D = (TH2F *)hEventsFT0M3D->Project3D("zxoe");
    if (!hEventsFT0M2D)
    {
      cout << "hCentFT0M not available " << endl;
      return;
    }
    hEventsFT0M2D->SetName("hCentFT0M");
  }

  if (SysPath.find("_Train100720") != string::npos)
  {
    hEventsFV0A2D = (TH2F *)dirEvt->Get("hCentFV0A");
    if (!hEventsFV0A2D)
    {
      cout << "hCentFV0A not available " << endl;
      return;
    }
  }
  else if (SysPath.find("_run528531") != string::npos || SysPath.find("_Train109") != string::npos || SysPath.find("_Train110") != string::npos)
  {
    hEventsFV0A = (TH1F *)dirEvt->Get("hCentFV0A");
    if (!hEventsFV0A)
    {
      cout << "hCentFV0A not available " << endl;
      return;
    }
  }

  if (evFlag == 2)
  { // INEL > 1
    hEventsFT0M = (TH1F *)hEventsFT0M2D->ProjectionX("hEventsFT0M", hEventsFT0M2D->GetYaxis()->FindBin(2.), hEventsFT0M2D->GetYaxis()->FindBin(2.));
    if (SysPath.find("_Train100720") != string::npos)
      hEventsFV0A = (TH1F *)hEventsFV0A2D->ProjectionX("hEventsFV0A", hEventsFV0A2D->GetYaxis()->FindBin(2.), hEventsFV0A2D->GetYaxis()->FindBin(2.));
  }
  else if (evFlag == 1)
  { // INEL > 0
    hEventsFT0M = (TH1F *)hEventsFT0M2D->ProjectionX("hEventsFT0M", hEventsFT0M2D->GetYaxis()->FindBin(1.), hEventsFT0M2D->GetYaxis()->FindBin(2.));
    if (SysPath.find("_Train100720") != string::npos)
      hEventsFV0A = (TH1F *)hEventsFV0A2D->ProjectionX("hEventsFV0A", hEventsFV0A2D->GetYaxis()->FindBin(1.), hEventsFV0A2D->GetYaxis()->FindBin(2.));
  }
  else
  { // INEL
    hEventsFT0M = (TH1F *)hEventsFT0M2D->ProjectionX("hEventsFT0M", hEventsFT0M2D->GetYaxis()->FindBin(0.), hEventsFT0M2D->GetYaxis()->FindBin(2.));
    if (SysPath.find("_Train100720") != string::npos)
      hEventsFV0A = (TH1F *)hEventsFV0A2D->ProjectionX("hEventsFV0A", hEventsFV0A2D->GetYaxis()->FindBin(0.), hEventsFV0A2D->GetYaxis()->FindBin(2.));
  }

  if (isMB == 1)
  {
    // NEvents = hEvents->GetBinContent(1);
    for (Int_t b = 1; b <= hEventsFT0M->GetNbinsX(); b++)
    {
      if (hEventsFT0M->GetBinCenter(b) > 0 && hEventsFT0M->GetBinCenter(b) < 100)
      {
        NEvents += hEventsFT0M->GetBinContent(b);
      }
    }
  }
  else
  {
    if (MultType == 1)
    {
      for (Int_t b = 1; b <= hEventsFT0M->GetNbinsX(); b++)
      {
        if (hEventsFT0M->GetBinCenter(b) > MultiplicityPerc[mul] && hEventsFT0M->GetBinCenter(b) < MultiplicityPerc[mul + 1])
        {
          NEvents += hEventsFT0M->GetBinContent(b);
        }
      }
    }
    else if (MultType == 2)
    {
      for (Int_t b = 1; b <= hEventsFV0A->GetNbinsX(); b++)
      {
        if (hEventsFV0A->GetBinCenter(b) > MultiplicityPerc[mul] && hEventsFV0A->GetBinCenter(b) < MultiplicityPerc[mul + 1])
        {
          NEvents += hEventsFV0A->GetBinContent(b);
        }
      }
    }
  }
  cout << "NEvents" << NEvents << endl;

  const Int_t numPt = 19; // 17 for omega, 19 for Xi
  Int_t numPtMax = 19;
  if (part == 6 || part == 7 || part == 8)
    numPtMax = 17;
  // default for 22m_pass4_Train79153
  // Xi
  /*
   Float_t binpt[numPt + 1] = {0.4, 0.6, 0.8, 1.0,
                              1.2, 1.4, 1.6, 1.8, 2.0,
                              2.2, 2.4, 2.6, 2.8, 3.0,
                              3.5, 4.0, 4.5, 5.0, 6.0, 8.0};
                              */
  // Omega
  /*
  Float_t binpt[numPt + 1] = {0.8, 1.0,
                              1.2, 1.4, 1.6, 1.8, 2.0,
                              2.2, 2.4, 2.6, 2.8, 3.0,
                              3.5, 4.0, 4.5, 5.0, 6.0, 8.0};
*/
  Float_t binpt[numPt + 1] = {0.4, 0.6, 0.8, 1.0,
                              1.2, 1.4, 1.6, 1.8, 2.0,
                              2.2, 2.4, 2.6, 2.8, 3.0,
                              3.5, 4.0, 4.5, 5.0, 6.0, 8.0};
  TString SPt[numPt] = {""};
  TH1F *hInvMass[numPt];

  TCanvas *canvas[4];
  for (Int_t c = 0; c < 4; c++)
  {
    canvas[c] = new TCanvas(Form("canvas_%i", c), Form("canvas%i", c), 1800, 1400);
    canvas[c]->Divide(3, 3);
    StyleCanvas(canvas[c], 0.15, 0.05, 0.05, 0.15);
  }
  TCanvas *canvasMass = new TCanvas("canvasMass", "canvasMass", 800, 1800);
  canvasMass->Divide(2, 3);
  StyleCanvas(canvasMass, 0.15, 0.05, 0.05, 0.15);

  TH1F *histoCountsPerEvent = new TH1F("histoCountsPerEvent", "histoCountsPerEvent", numPt, binpt);
  TH1F *histoMean = new TH1F("histoMean", "histoMean", numPt, binpt);
  TH1F *histoSigma = new TH1F("histoSigma", "histoSigma", numPt, binpt);
  TH1F *histoSigmaNarrow = new TH1F("histoSigmaNarrow", "histoSigmaNarrow", numPt, binpt);
  TH1F *histoPurity = new TH1F("histoPurity", "histoPurity", numPt, binpt);
  TH1F *histoYield = new TH1F("histoYield", "histoYield", numPt, binpt);
  TH1F *histoSignificance = new TH1F("histoSignificance", "histoSignificance", numPt, binpt);

  Float_t counts = 0;
  Float_t errcount = 0;
  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      if (binpt[pt] < 0.8)
        continue;
    }
    SPt[pt] = Form("%.1f < p_{T} < %.1f", binpt[pt], binpt[pt + 1]);
    cout << "Analysed pt interval: " << binpt[pt] << "-" << binpt[pt + 1] << endl;
    cout << binpt[pt] << endl;

    hInvMass[pt] = (TH1F *)h2->ProjectionY(Form("hInvMass_pt%i", pt), h2->GetXaxis()->FindBin(binpt[pt] + 0.001), h2->GetXaxis()->FindBin(binpt[pt + 1] - 0.001));
    hInvMass[pt]->Rebin(MassRebin);
    StyleHisto(hInvMass[pt], 0, 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), 1, 20, TitleInvMass[part] + " " + SInvMass, "Counts", SPt[pt] + " GeV/#it{c}", 1, histoMassRangeLow[part], histoMassRangeUp[part], 1.4, 1.6, 0.7);
    if (isYAxisMassZoomed)
    {
      if (part >= 6 && part <= 8)
        hInvMass[pt]->GetYaxis()->SetRangeUser(0, 2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(1.65)));
      else if (part >= 3 && part <= 5)
        hInvMass[pt]->GetYaxis()->SetRangeUser(0, 2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->FindBin(1.29)));
    }
    if (pt < 9)
      canvas[0]->cd(pt + 1);
    else if (pt < 18)
      canvas[1]->cd(pt + 1 - 9);
    else if (pt < 27)
      canvas[2]->cd(pt + 1 - 18);
    else
      canvas[3]->cd(pt + 1 - 27);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.1);
    gPad->SetBottomMargin(0.2);
    // hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);
    hInvMass[pt]->Draw("e same");

    counts = 0;
    errcount = 0;
    for (Int_t bmass = hInvMass[pt]->GetXaxis()->FindBin(LowMassRange[part]); bmass <= hInvMass[pt]->GetXaxis()->FindBin(UpMassRange[part]); bmass++)
    {
      counts += hInvMass[pt]->GetBinContent(bmass);
    }
    errcount = sqrt(counts);
    histoCountsPerEvent->SetBinContent(pt + 1, counts / NEvents / histoCountsPerEvent->GetBinWidth(pt + 1));
    histoCountsPerEvent->SetBinError(pt + 1, errcount / NEvents / histoCountsPerEvent->GetBinWidth(pt + 1));
  }

  // fits

  TF1 **functionsFirst = new TF1 *[numPt];
  TF1 **functionsSecond = new TF1 *[numPt];
  TF1 **functions1 = new TF1 *[numPt];
  TF1 **functions2 = new TF1 *[numPt];
  TF1 **bkg1 = new TF1 *[numPt];
  TF1 **bkg2 = new TF1 *[numPt];
  TF1 **bkgretta = new TF1 *[numPt]; // initial bkg fit
  TF1 **bkgparab = new TF1 *[numPt]; // initial bkg fit
  TF1 **total = new TF1 *[numPt];
  TF1 **totalbis = new TF1 *[numPt];

  Double_t parTwoGaussParab[numPt + 1][9];
  Double_t parTwoGaussRetta[numPt + 1][8];
  Double_t parOneGaussParab[numPt + 1][6];
  Double_t parOneGaussRetta[numPt + 1][5];

  TFitResultPtr fFitResultPtr0[numPt];
  TFitResultPtr fFitResultPtr1[numPt];

  Float_t mean[numPt] = {0};
  Float_t errmean[numPt] = {0};
  Float_t sigma[numPt] = {0};
  Float_t errsigma[numPt] = {0};
  Float_t sigmaNarrow[numPt] = {0};
  Float_t errsigmaNarrow[numPt] = {0};
  Float_t LowLimit[numPt] = {0};
  Float_t UpLimit[numPt] = {0};
  Float_t b[numPt] = {0};
  Float_t errb[numPt] = {0};
  Float_t SSB[numPt] = {0};
  Float_t errSSB[numPt] = {0};
  Float_t entries_range[numPt] = {0};
  Float_t Yield[numPt] = {0};
  Float_t ErrYield[numPt] = {0};
  Float_t TotYield = 0;
  Float_t TotSigBkg = 0;

  TLine *lineP3Sigma[numPt];
  TLine *lineM3Sigma[numPt];

  Float_t bTest[numPt] = {0};
  Float_t errbTest[numPt] = {0};
  Float_t SignalTest[numPt] = {0};
  Float_t errSignalTest[numPt] = {0};
  Float_t YieldTest[numPt] = {0};
  Float_t ErrYieldTest[numPt] = {0};
  TH1F *hYieldTest[numPt];
  TH1F *hYieldRelErrorTest[numPt];
  TH1F *hYieldRelErrorTestRelative[numPt];
  Float_t LowLimitTest[numPt] = {0};
  Float_t UpLimitTest[numPt] = {0};
  Float_t LowBin0[numPt] = {0};
  Float_t UpBin0[numPt] = {0};
  Float_t LowBin[numPt] = {0};
  Float_t UpBin[numPt] = {0};
  Int_t numMassInt = 10;
  TF1 *LineAt1 = new TF1("LineAt1", "[0]+[1]*x", 0, numMassInt);
  LineAt1->SetParameter(0, 1);
  TF1 *LineAt995 = new TF1("LineAt995", "[0]+[1]*x", 0, numMassInt);
  LineAt995->SetParameter(0, 0.995);

  TLine *lineBkgLimitA[numPt];
  TLine *lineBkgLimitB[numPt];
  TLine *lineBkgLimitC[numPt];
  TLine *lineBkgLimitD[numPt];

  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      if (binpt[pt] < 0.8)
        continue;
    }
    if (part == 6 || part == 7 || part == 8)
    { // Omega
      liminf[part] = 1.63;
      limsup[part] = 1.71;
      if (SysSigExtr == 2)
      {
        liminf[part] = 1.635;
        limsup[part] = 1.705;
      }
      if (SysSigExtr == 3)
      {
        liminf[part] = 1.626;
        limsup[part] = 1.715;
      }
    }
    if (part == 3 || part == 4 || part == 5)
    { // Xi
      liminf[part] = 1.29;
      limsup[part] = 1.355;
    }

    if (part == 6 || part == 7 || part == 8)
    {
      if (binpt[pt] < 0.6)
        continue;
    }

    if (pt < 9)
      canvas[0]->cd(pt + 1);
    else if (pt < 18)
      canvas[1]->cd(pt + 1 - 9);
    else if (pt < 27)
      canvas[2]->cd(pt + 1 - 18);
    else
      canvas[3]->cd(pt + 1 - 27);

    functionsFirst[pt] = new TF1(Form("1f_%i", pt), "gaus", min_range_signal[part], max_range_signal[part]);
    functionsFirst[pt]->SetLineColor(881);
    functionsFirst[pt]->SetParameter(1, massParticle[part]);
    functionsFirst[pt]->SetParName(0, "norm");
    functionsFirst[pt]->SetParName(1, "mean");
    functionsFirst[pt]->SetParName(2, "sigma");
    functionsFirst[pt]->SetParLimits(1, min_range_signal[part], max_range_signal[part]);
    functionsFirst[pt]->SetParLimits(2, 0.001, 0.1);
    functionsFirst[pt]->SetParLimits(0, 0, 1.1 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));

    functionsSecond[pt] = new TF1(Form("2f_%i", pt), "gaus", min_range_signal[part], max_range_signal[part]);
    functionsSecond[pt]->SetLineColor(867);
    functionsSecond[pt]->SetParameter(1, massParticle[part]);
    functionsSecond[pt]->SetParName(0, "norm");
    functionsSecond[pt]->SetParName(1, "mean");
    functionsSecond[pt]->SetParName(2, "sigma");
    functionsSecond[pt]->SetParLimits(1, min_range_signal[part], max_range_signal[part]);
    functionsSecond[pt]->SetParLimits(2, 0.001, 0.15);
    functionsSecond[pt]->SetParLimits(0, 0, 1.1 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));

    functions1[pt] = new TF1(Form("1f_%i_final", pt), "gaus", gaussDisplayRangeLow[part], gaussDisplayRangeUp[part]);
    functions1[pt]->SetLineColor(kRed); // 867
    functions1[pt]->SetParName(0, "norm");
    functions1[pt]->SetParName(1, "mean");
    functions1[pt]->SetParName(2, "sigma");

    functions2[pt] = new TF1(Form("2f_%i_final", pt), "gaus", gaussDisplayRangeLow[part], gaussDisplayRangeUp[part]);
    functions2[pt]->SetLineColor(kMagenta); // 891
    functions2[pt]->SetParName(0, "norm");
    functions2[pt]->SetParName(1, "mean");
    functions2[pt]->SetParName(2, "sigma");

    bkg1[pt] = new TF1(Form("bkg1%i", pt), "pol1", bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
    bkg1[pt]->SetLineColor(418);
    bkg1[pt]->SetLineStyle(2);

    bkg2[pt] = new TF1(Form("bkg2%i", pt), "pol2", bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
    bkg2[pt]->SetLineColor(1);
    bkg2[pt]->SetLineStyle(2);

    bkgretta[pt] = new TF1(Form("retta%i", pt), fretta, liminf[part], limsup[part], 3);
    bkgretta[pt]->SetLineColor(kGreen + 3);
    bkgretta[pt]->FixParameter(2, part);

    bkgparab[pt] = new TF1(Form("parab%i", pt), fparab, liminf[part], limsup[part], 4);
    bkgparab[pt]->SetLineColor(kAzure + 7);
    bkgparab[pt]->FixParameter(3, part);

    Bool_t UseTwoGaussUpdated = 1;
    if (UseTwoGauss)
    {
      cout << "\n\e[35mFit with two gauss \e[39m"
           << " Pt: " << binpt[pt] << "-" << binpt[pt + 1] << endl;

      if (isBkgParab)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+pol2(6)", liminf[part], limsup[part]);
      else
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+gaus(3)+pol1(6)", liminf[part], limsup[part]);
      total[pt]->SetLineColor(597);
      total[pt]->SetParName(0, "norm");
      total[pt]->SetParName(1, "mean");
      total[pt]->SetParName(2, "sigma");
      total[pt]->SetParName(3, "norm2");
      total[pt]->SetParName(4, "mean2");
      total[pt]->SetParName(5, "sigma2");

      cout << "\n\n fit gauss1 " << endl;
      hInvMass[pt]->Fit(functionsFirst[pt], "RB");
      cout << "\n\n fit gauss2 " << endl;
      hInvMass[pt]->Fit(functionsSecond[pt], "RB");

      bkg1[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkg2[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkgparab[pt]->SetRange(liminf[part], limsup[part]);
      bkgretta[pt]->SetRange(liminf[part], limsup[part]);
      total[pt]->SetRange(liminf[part], limsup[part]);

      cout << "\n\n fit bkg " << endl;
      if (isBkgParab)
        hInvMass[pt]->Fit(bkgparab[pt], "RB0");
      else
        hInvMass[pt]->Fit(bkgretta[pt], "RB0");

      functionsFirst[pt]->GetParameters(&parTwoGaussParab[pt][0]);
      functionsFirst[pt]->GetParameters(&parTwoGaussRetta[pt][0]);
      functionsSecond[pt]->GetParameters(&parTwoGaussParab[pt][3]);
      functionsSecond[pt]->GetParameters(&parTwoGaussRetta[pt][3]);
      if (isBkgParab)
      {
        bkgparab[pt]->GetParameters(&parTwoGaussParab[pt][6]);
        total[pt]->SetParameters(parTwoGaussParab[pt]);
      }
      else
      {
        bkgretta[pt]->GetParameters(&parTwoGaussRetta[pt][6]);
        total[pt]->SetParameters(parTwoGaussRetta[pt]);
      }

      cout << "\n\n fit total " << endl;
      if (Spart[part] == "XiNeg" || Spart[part] == "XiPos" || Spart[part] == "Xi")
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.31, 1.335);
        total[pt]->SetParLimits(2, 0.0012, 0.010);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.31, 1.335);
        total[pt]->SetParLimits(5, 0.001, 0.01);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, massParticle[part]);
          total[pt]->FixParameter(4, massParticle[part]);
        }
      }
      else if (Spart[part] == "OmegaNeg" || Spart[part] == "OmegaPos" || Spart[part] == "Omega")
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.66, 1.68);
        total[pt]->SetParLimits(2, 0.002, 0.01);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.66, 1.68);
        total[pt]->SetParLimits(5, 0.001, 0.01);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, massParticle[part]);
          total[pt]->FixParameter(4, massParticle[part]);
        }
      }
      else if (Spart[part] == "K0s")
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 0.485, 0.505);
        total[pt]->SetParLimits(2, 0.001, 0.01);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(4, 0.485, 0.505);
        total[pt]->SetParLimits(5, 0.001, 0.015);
        // total[pt]->SetParLimits(6, -2000, 2000);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, massParticle[part]);
          total[pt]->FixParameter(4, massParticle[part]);
        }
        cout << "max value " << hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()) << endl;
      }

      fFitResultPtr0[pt] = hInvMass[pt]->Fit(total[pt], "SRB+"); // per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0
      // la gaussiana più larga deve esserte quella più bassa
      if (total[pt]->GetParameter(2) > total[pt]->GetParameter(5))
      {
        if (total[pt]->GetParameter(0) > total[pt]->GetParameter(3))
          UseTwoGaussUpdated = kFALSE;
      }
      else
      {
        if (total[pt]->GetParameter(0) < total[pt]->GetParameter(3))
          UseTwoGaussUpdated = kFALSE;
      }

      cout << "UseTwoGauss = " << UseTwoGaussUpdated << endl;

      totalbis[pt] = (TF1 *)total[pt]->Clone();
      fFitResultPtr1[pt] = fFitResultPtr0[pt];

      functions1[pt]->FixParameter(0, total[pt]->GetParameter(0));
      functions1[pt]->FixParameter(1, total[pt]->GetParameter(1));
      functions1[pt]->FixParameter(2, total[pt]->GetParameter(2));
      functions2[pt]->FixParameter(0, total[pt]->GetParameter(3));
      functions2[pt]->FixParameter(1, total[pt]->GetParameter(4));
      functions2[pt]->FixParameter(2, total[pt]->GetParameter(5));

      totalbis[pt]->FixParameter(0, 0);
      totalbis[pt]->FixParameter(1, 0);
      totalbis[pt]->FixParameter(2, 0);
      totalbis[pt]->FixParameter(3, 0);
      totalbis[pt]->FixParameter(4, 0);
      totalbis[pt]->FixParameter(5, 0);

      if (isBkgParab)
      {
        bkg2[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg2[pt]->FixParameter(1, total[pt]->GetParameter(7));
        bkg2[pt]->FixParameter(2, total[pt]->GetParameter(8));
      }
      else
      {
        bkg1[pt]->FixParameter(0, total[pt]->GetParameter(6));
        bkg1[pt]->FixParameter(1, total[pt]->GetParameter(7));
      }

      if (UseTwoGaussUpdated)
      {

        if (pt < 9)
          canvas[0]->cd(pt + 1);
        else if (pt < 18)
          canvas[1]->cd(pt + 1 - 9);
        else if (pt < 27)
          canvas[2]->cd(pt + 1 - 18);
        else
          canvas[3]->cd(pt + 1 - 27);
        hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);
        hInvMass[pt]->Draw("same e");
        functions1[pt]->Draw("same");
        functions2[pt]->Draw("same");
        if (isBkgParab)
          bkg2[pt]->Draw("same");
        else
          bkg1[pt]->Draw("same");

        TMatrixDSym cov = fFitResultPtr0[pt]->GetCovarianceMatrix();
        Double_t cov_mean = cov[1][4];
        Double_t cov_sigma = cov[2][5];
        mean[pt] = (functions1[pt]->GetParameter(1) + functions2[pt]->GetParameter(1)) / 2;
        errmean[pt] = (total[pt]->GetParError(1) + total[pt]->GetParError(4)) / 2;
        sigma[pt] = (functions1[pt]->GetParameter(2) + functions2[pt]->GetParameter(2)) / 2;
        sigmaNarrow[pt] = functions1[pt]->GetParameter(2);
        errsigmaNarrow[pt] = functions1[pt]->GetParError(2);
        if (functions1[pt]->GetParameter(2) > functions2[pt]->GetParameter(2))
        {
          sigmaNarrow[pt] = functions2[pt]->GetParameter(2);
          errsigmaNarrow[pt] = functions2[pt]->GetParError(2);
        }
        errsigma[pt] = sqrt(pow(total[pt]->GetParError(2), 2) + pow(total[pt]->GetParError(5), 2) + 2 * cov_sigma) / 2;
      }
    }
    if (!UseTwoGaussUpdated || !UseTwoGauss)
    {
      cout << "\n\e[36mFit with one gauss only: \e[39m"
           << " Pt: " << binpt[pt] << "-" << binpt[pt + 1] << endl;

      if (isBkgParab)
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+pol2(3)", liminf[part], limsup[part]);
      else
        total[pt] = new TF1(Form("total%i", pt), "gaus(0)+pol1(3)", liminf[part], limsup[part]);
      total[pt]->SetLineColor(7);
      total[pt]->SetParName(0, "norm");
      total[pt]->SetParName(1, "mean");
      total[pt]->SetParName(2, "sigma");

      cout << "\n\n fit gauss " << endl;
      hInvMass[pt]->Fit(functionsFirst[pt], "RB");

      bkg1[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkg2[pt]->SetRange(bkgDisplayRangeLow[part], bkgDisplayRangeUp[part]);
      bkgparab[pt]->SetRange(liminf[part], limsup[part]);
      bkgretta[pt]->SetRange(liminf[part], limsup[part]);
      total[pt]->SetRange(liminf[part], limsup[part]);

      cout << "\n\n fit bkg " << endl;
      if (isBkgParab)
        hInvMass[pt]->Fit(bkgparab[pt], "RB0");
      else
        hInvMass[pt]->Fit(bkgretta[pt], "RB0");

      functionsFirst[pt]->GetParameters(&parOneGaussParab[pt][0]);
      functionsFirst[pt]->GetParameters(&parOneGaussRetta[pt][0]);
      if (isBkgParab)
      {
        bkgparab[pt]->GetParameters(&parOneGaussParab[pt][3]);
        total[pt]->SetParameters(parOneGaussParab[pt]);
      }
      else
      {
        bkgretta[pt]->GetParameters(&parOneGaussRetta[pt][3]);
        total[pt]->SetParameters(parOneGaussRetta[pt]);
      }

      cout << "\n\n fit total " << endl;
      if (Spart[part] == "XiNeg" || Spart[part] == "XiPos" || Spart[part] == "Xi")
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.31, 1.335);
        total[pt]->SetParLimits(2, 0.0012, 0.010);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, massParticle[part]);
        }
      }
      else if (Spart[part] == "OmegaNeg" || Spart[part] == "OmegaPos" || Spart[part] == "Omega")
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 1.66, 1.68);
        total[pt]->SetParLimits(2, 0.001, 0.02);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, massParticle[part]);
        }
      }
      else if (Spart[part] == "K0s")
      {
        total[pt]->SetParLimits(0, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()));
        total[pt]->SetParLimits(1, 0.485, 0.505);
        total[pt]->SetParLimits(2, 0.001, 0.01);
        if (isMeanFixedPDG)
        {
          total[pt]->FixParameter(1, massParticle[part]);
        }
      }
      cout << "max value " << hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()) << endl;
      fFitResultPtr0[pt] = hInvMass[pt]->Fit(total[pt], "SRB+"); // per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0

      totalbis[pt] = (TF1 *)total[pt]->Clone();
      fFitResultPtr1[pt] = fFitResultPtr0[pt];

      functions1[pt]->FixParameter(0, total[pt]->GetParameter(0));
      functions1[pt]->FixParameter(1, total[pt]->GetParameter(1));
      functions1[pt]->FixParameter(2, total[pt]->GetParameter(2));

      totalbis[pt]->FixParameter(0, 0);
      totalbis[pt]->FixParameter(1, 0);
      totalbis[pt]->FixParameter(2, 0);

      if (isBkgParab)
      {
        bkg2[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg2[pt]->FixParameter(1, total[pt]->GetParameter(4));
        bkg2[pt]->FixParameter(2, total[pt]->GetParameter(5));
      }
      else
      {
        bkg1[pt]->FixParameter(0, total[pt]->GetParameter(3));
        bkg1[pt]->FixParameter(1, total[pt]->GetParameter(4));
      }
      if (pt < 9)
        canvas[0]->cd(pt + 1);
      else if (pt < 18)
        canvas[1]->cd(pt + 1 - 9);
      else if (pt < 27)
        canvas[2]->cd(pt + 1 - 18);
      else
        canvas[3]->cd(pt + 1 - 27);
      hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);
      hInvMass[pt]->Draw("same e");

      if (isBkgParab)
        bkg2[pt]->Draw("same");
      else
        bkg1[pt]->Draw("same");
      functions1[pt]->Draw("same");

      mean[pt] = total[pt]->GetParameter(1);
      errmean[pt] = total[pt]->GetParError(1);
      sigma[pt] = total[pt]->GetParameter(2);
      sigmaNarrow[pt] = sigma[pt];
      errsigma[pt] = total[pt]->GetParError(2);
      errsigmaNarrow[pt] = errsigma[pt];
    }

    cout << "\nMean: " << mean[pt] << " +/- " << errmean[pt] << endl;
    cout << "Sigma: " << sigma[pt] << " +/- " << errsigma[pt] << endl;

    TLine *linebkgFitLL = new TLine(liminf[part], 0, liminf[part], hInvMass[pt]->GetMaximum()); // low limit of left SB
    TLine *linebkgFitRR = new TLine(limsup[part], 0, limsup[part], hInvMass[pt]->GetMaximum()); // upper limit of right SB
    linebkgFitLL->SetLineColor(kBlue);
    linebkgFitRR->SetLineColor(kBlue);
    TLine *linebkgFitLR = new TLine(UpperLimitLSB, 0, UpperLimitLSB, hInvMass[pt]->GetMaximum()); // upper limit of left SB
    TLine *linebkgFitRL = new TLine(LowerLimitRSB, 0, LowerLimitRSB, hInvMass[pt]->GetMaximum()); // lower limit of right SB
    linebkgFitLR->SetLineColor(kBlue);
    linebkgFitRL->SetLineColor(kBlue);
    lineBkgLimitA[pt] = new TLine(liminf[part], 0, liminf[part], hInvMass[pt]->GetMaximum());
    lineBkgLimitB[pt] = new TLine(UpperLimitLSB, 0, UpperLimitLSB, hInvMass[pt]->GetMaximum());
    lineBkgLimitC[pt] = new TLine(limsup[part], 0, limsup[part], hInvMass[pt]->GetMaximum());
    lineBkgLimitD[pt] = new TLine(LowerLimitRSB, 0, LowerLimitRSB, hInvMass[pt]->GetMaximum());
    lineBkgLimitA[pt]->SetLineColor(kViolet + 1);
    lineBkgLimitB[pt]->SetLineColor(kViolet + 1);
    lineBkgLimitC[pt]->SetLineColor(kViolet + 1);
    lineBkgLimitD[pt]->SetLineColor(kViolet + 1);
    lineBkgLimitA[pt]->Draw("same");
    lineBkgLimitB[pt]->Draw("same");
    lineBkgLimitC[pt]->Draw("same");
    lineBkgLimitD[pt]->Draw("same");

    // linebkgFitLL->Draw("same");
    // linebkgFitRR->Draw("same");
    // linebkgFitLR->Draw("same");
    // linebkgFitRL->Draw("same");

    if (isBkgParab)
      bkgparab[pt]->Draw("same");
    else
      bkgretta[pt]->Draw("same");

    if (pt < 9)
      canvas[0]->cd(pt + 1);
    else if (pt < 18)
      canvas[1]->cd(pt + 1 - 9);
    else if (pt < 27)
      canvas[2]->cd(pt + 1 - 18);
    else
      canvas[3]->cd(pt + 1 - 27);

    LowLimit[pt] = hInvMass[pt]->GetXaxis()->GetBinLowEdge(hInvMass[pt]->GetXaxis()->FindBin(mean[pt] - sigmacentral * sigma[pt]));
    UpLimit[pt] = hInvMass[pt]->GetXaxis()->GetBinUpEdge(hInvMass[pt]->GetXaxis()->FindBin(mean[pt] + sigmacentral * sigma[pt]));

    lineP3Sigma[pt] = new TLine(UpLimit[pt], 0, UpLimit[pt], hInvMass[pt]->GetMaximum());
    lineM3Sigma[pt] = new TLine(LowLimit[pt], 0, LowLimit[pt], hInvMass[pt]->GetMaximum());
    lineP3Sigma[pt]->Draw("same");
    lineM3Sigma[pt]->Draw("same");

    b[pt] = 0;
    errb[pt] = 0;
    if (isBkgParab)
    {
      b[pt] = bkg2[pt]->Integral(LowLimit[pt], UpLimit[pt]);
      errb[pt] = totalbis[pt]->IntegralError(LowLimit[pt], UpLimit[pt], fFitResultPtr1[pt]->GetParams(),
                                             (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
    }
    else
    {
      b[pt] = bkg1[pt]->Integral(LowLimit[pt], UpLimit[pt]);
      errb[pt] = totalbis[pt]->IntegralError(LowLimit[pt], UpLimit[pt], fFitResultPtr1[pt]->GetParams(),
                                             (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
    }
    b[pt] = b[pt] / hInvMass[pt]->GetBinWidth(1);
    errb[pt] = errb[pt] / hInvMass[pt]->GetBinWidth(1);

    entries_range[pt] = 0;
    for (Int_t l = hInvMass[pt]->GetXaxis()->FindBin(mean[pt] - sigmacentral * sigma[pt]); l <= hInvMass[pt]->GetXaxis()->FindBin(mean[pt] + sigmacentral * sigma[pt]); l++)
    { // I inlcude bins where the limits lie
      entries_range[pt] += hInvMass[pt]->GetBinContent(l);
    }

    Yield[pt] = entries_range[pt] - b[pt];
    ErrYield[pt] = sqrt(entries_range[pt] + pow(errb[pt], 2));
    TotYield += Yield[pt];
    TotSigBkg += entries_range[pt];

    SSB[pt] = (entries_range[pt] - b[pt]) / entries_range[pt];
    errSSB[pt] = SSB[pt] * sqrt(1. / entries_range[pt] + pow(errb[pt] / b[pt], 2));

    // Study of integration range ******************
    hYieldTest[pt] = new TH1F(Form("hYieldTest%d", pt), "", numMassInt, 0, numMassInt);
    hYieldRelErrorTest[pt] = new TH1F(Form("hYieldRelErrorTest%d", pt), "", numMassInt, 0, numMassInt);
    hYieldRelErrorTestRelative[pt] = new TH1F(Form("hYieldRelErrorTestRelative%d", pt), "", numMassInt, 0, numMassInt);

    LowBin0[pt] = hInvMass[pt]->GetXaxis()->FindBin(mean[pt] - 3 * sigma[pt]);
    UpBin0[pt] = hInvMass[pt]->GetXaxis()->FindBin(mean[pt] + 3 * sigma[pt]);

    Float_t NsigmaTest = 0;
    for (Int_t i = 0; i < numMassInt; i++)
    {
      NsigmaTest = 3 + 0.2 * i;
      SignalTest[pt] = 0;
      // LowBin[pt] = LowBin0[pt] - i;
      LowBin[pt] = hInvMass[pt]->GetXaxis()->FindBin(mean[pt] - NsigmaTest * sigma[pt]);
      // UpBin[pt] = UpBin0[pt] + i;
      UpBin[pt] = hInvMass[pt]->GetXaxis()->FindBin(mean[pt] + NsigmaTest * sigma[pt]);
      LowLimitTest[pt] = hInvMass[pt]->GetXaxis()->GetBinLowEdge(LowBin[pt]);
      UpLimitTest[pt] = hInvMass[pt]->GetXaxis()->GetBinUpEdge(UpBin[pt]);
      for (Int_t l = LowBin[pt]; l <= UpBin[pt]; l++)
      {
        SignalTest[pt] += hInvMass[pt]->GetBinContent(l);
      }
      if (isBkgParab)
      {
        bTest[pt] = bkg2[pt]->Integral(LowLimitTest[pt], UpLimitTest[pt]);
        errbTest[pt] = totalbis[pt]->IntegralError(LowLimitTest[pt], UpLimitTest[pt], fFitResultPtr1[pt]->GetParams(),
                                                   (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
      }
      else
      {
        bTest[pt] = bkg1[pt]->Integral(LowLimitTest[pt], UpLimitTest[pt]);
        errbTest[pt] = totalbis[pt]->IntegralError(LowLimitTest[pt], UpLimitTest[pt], fFitResultPtr1[pt]->GetParams(),
                                                   (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
      }
      bTest[pt] = bTest[pt] / hInvMass[pt]->GetBinWidth(1);
      errbTest[pt] = errbTest[pt] / hInvMass[pt]->GetBinWidth(1);
      YieldTest[pt] = SignalTest[pt] - bTest[pt];
      ErrYieldTest[pt] = sqrt(SignalTest[pt] + pow(errbTest[pt], 2));
      hYieldTest[pt]->SetBinContent(i + 1, YieldTest[pt]);
      hYieldRelErrorTest[pt]->SetBinContent(i + 1, ErrYieldTest[pt] / YieldTest[pt]);
      hYieldRelErrorTestRelative[pt]->SetBinContent(i + 1, ErrYieldTest[pt] / YieldTest[pt] / hYieldRelErrorTest[pt]->GetBinContent(1));
      // hYieldTest[pt]->SetBinError(i + 1, ErrYieldTest[pt]);
      hYieldTest[pt]->SetBinError(i + 1, 0);
      // hYieldTest[pt]->GetXaxis()->SetBinLabel(i + 1, Form("%.2f", (UpLimitTest[pt] - LowLimitTest[pt]) / sigma[pt] / 2));
      hYieldTest[pt]->GetXaxis()->SetBinLabel(i + 1, Form("%.1f", NsigmaTest));
      // hYieldRelErrorTest[pt]->GetXaxis()->SetBinLabel(i + 1, Form("%.2f", (UpLimitTest[pt] - LowLimitTest[pt]) / sigma[pt] / 2));
      hYieldRelErrorTest[pt]->GetXaxis()->SetBinLabel(i + 1, Form("%.1f", NsigmaTest));
      // hYieldRelErrorTestRelative[pt]->GetXaxis()->SetBinLabel(i + 1, Form("%.2f", (UpLimitTest[pt] - LowLimitTest[pt]) / sigma[pt] / 2));
      hYieldRelErrorTestRelative[pt]->GetXaxis()->SetBinLabel(i + 1, Form("%.1f", NsigmaTest));
    }
    hYieldTest[pt]->Scale(1. / hYieldTest[pt]->GetBinContent(numMassInt));

    for (Int_t i = 0; i < numMassInt; i++)
    {
      // hYieldTest[pt]->SetBinContent(i + 1, TMath::Abs(1 - hYieldTest[pt]->GetBinContent(i + 1)));
      //  hYieldTest[pt]->SetBinError(i + 1, );
    }

    //*********************************************

    histoYield->SetBinContent(pt + 1, Yield[pt] / NEvents / histoYield->GetBinWidth(pt + 1));
    histoYield->SetBinError(pt + 1, ErrYield[pt] / NEvents / histoCountsPerEvent->GetBinWidth(pt + 1));

    histoMean->SetBinContent(pt + 1, mean[pt]);
    histoMean->SetBinError(pt + 1, errmean[pt]);

    histoSigma->SetBinContent(pt + 1, sigma[pt]);
    histoSigma->SetBinError(pt + 1, errsigma[pt]);

    histoSigmaNarrow->SetBinContent(pt + 1, sigmaNarrow[pt]);
    histoSigmaNarrow->SetBinError(pt + 1, errsigmaNarrow[pt]);

    histoPurity->SetBinContent(pt + 1, SSB[pt]);
    histoPurity->SetBinError(pt + 1, errSSB[pt]);

    histoSignificance->SetBinContent(pt + 1, Yield[pt] / ErrYield[pt]);
    histoSignificance->SetBinError(pt + 1, 0);
  }

  TotYield = TotYield / NEvents;
  TotSigBkg = TotSigBkg / NEvents;
  histoYield->SetLineColor(kRed);

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 1000, 800);
  if (part == 6 || part == 7 || part == 8)
    histoCountsPerEvent->GetYaxis()->SetRangeUser(0, 0.00007);
  histoCountsPerEvent->Draw("same");
  histoYield->DrawClone("same");

  TLegend *legendPt = new TLegend(0.5, 0.6, 0.9, 0.88);
  legendPt->SetTextAlign(21);
  legendPt->SetFillColor(0);

  TLegend *legendYield = new TLegend(0.7, 0.7, 0.9, 0.9);
  legendYield->AddEntry(histoCountsPerEvent, "w/o bkg subtraction", "pl");
  legendYield->AddEntry(histoYield, "w/ bkg subtraction", "pl");
  legendYield->Draw("");

  TCanvas *canvasSummary = new TCanvas("canvasSummary", "canvasSummary", 1700, 1000);
  canvasSummary->Divide(3, 2);
  TCanvas *canvasSummaryBis = new TCanvas("canvasSummaryBis", "canvasSummaryBis", 1950, 700);
  canvasSummaryBis->Divide(3, 1);

  canvasSummary->cd(1);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoMean, gaussDisplayRangeLow[part], gaussDisplayRangeUp[part], 1, 1, titlePt, "#mu (GeV/c^{2})", "histoMean", 0, 0, 0, 1.4, 1.4, 1.2);
  histoMean->Draw("");
  canvasSummary->cd(2);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoSigma, 0, 0.010, 1, 1, titlePt, "#sigma (GeV/c^{2})", "histoSigma", 0, 0, 0, 1.4, 1.4, 1.2);
  StyleHisto(histoSigmaNarrow, 0, 0.010, kAzure + 1, 1, titlePt, "#sigma (GeV/c^{2})", "histoSigma", 0, 0, 0, 1.4, 1.4, 1.2);
  if (part >= 6)
    histoSigma->GetYaxis()->SetRangeUser(0, 0.02);
  if (part >= 6)
    histoSigmaNarrow->GetYaxis()->SetRangeUser(0, 0.02);
  histoSigma->Draw("");
  // histoSigmaNarrow->Draw("");
  canvasSummary->cd(3);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoPurity, 0, 1, 1, 1, titlePt, "S / (S+B)", "histoPurity", 0, 0, 0, 1.4, 1.4, 1.2);
  histoPurity->Draw("");
  canvasSummary->cd(4);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoYield, 0, 1.2 * histoYield->GetBinContent(histoYield->GetMaximumBin()), 1, 1, titlePt, titleYield, "histoYield", 0, 0, 0, 1.4, 1.4, 1.2);
  histoYield->Draw("same");

  /*
    canvasSummaryBis->cd(1);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.14);
    StyleHisto(histoSignificance, 0, 1.2 * histoSignificance->GetBinContent(histoSignificance->GetMaximumBin()), 1, 1, titlePt, "Yield/#sigma_{Yield}", "histoSignificance", 0, 0, 0, 1.4, 1.4, 1.2);
    histoSignificance->Draw("same");
  */
  canvasSummaryBis->cd(1);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.04);
  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      if (binpt[pt] < 0.8)
        continue;
    }
    if (pt % 2 == 0)
      continue;
    StyleHisto(hYieldTest[pt], 0.95, 1.05, ColorPt[pt], 1, "N_{sigma}", "Normalised raw yield", "", 0, 0, 0, 1.4, 1.8, 1.2);
    legendPt->AddEntry(hYieldTest[pt], Form("%.1f < p_{T} < %.1f GeV/#it{c}", binpt[pt], binpt[pt + 1]), "pl");
    hYieldTest[pt]->Draw("same");
  }
  legendPt->Draw("");
  LineAt1->SetLineColor(1);
  LineAt1->SetLineStyle(10);
  LineAt1->Draw("same");
  LineAt995->SetLineColor(1);
  LineAt995->SetLineStyle(10);
  LineAt995->Draw("same");

  canvasSummaryBis->cd(2);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.04);
  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      if (binpt[pt] < 0.8)
        continue;
    }
    if (pt % 2 == 0)
      continue;
    StyleHisto(hYieldRelErrorTest[pt], 0., 0.03, ColorPt[pt], 1, "N_{sigma}", "Rel error (%)", "", 0, 0, 0, 1.4, 1.8, 1.2);
    hYieldRelErrorTest[pt]->Draw("same");
  }
  legendPt->Draw("");

  canvasSummaryBis->cd(3);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.04);
  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      if (binpt[pt] < 0.8)
        continue;
    }
    if (pt % 2 == 0)
      continue;
    StyleHisto(hYieldRelErrorTestRelative[pt], 1., 1.6, ColorPt[pt], 1, "N_{sigma}", "Normalised rel error (%)", "", 0, 0, 0, 1.4, 1.8, 1.2);
    hYieldRelErrorTestRelative[pt]->Draw("same");
  }
  legendPt->Draw("");

  Int_t index = 0;
  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      if (binpt[pt] < 0.8)
        continue;
    }

    if (pt == 0)
      index = 1;
    else if (pt == 3)
      index = 2;
    else if (pt == 6)
      index = 3;
    else if (pt == 9)
      index = 4;
    else if (pt == 12)
      index = 5;
    else if (pt == 15)
      index = 6;
    else
      continue;

    canvasMass->cd(index);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.18);

    hInvMass[pt]->GetXaxis()->SetRangeUser(histoMassRangeLow[part], histoMassRangeUp[part]);
    hInvMass[pt]->Draw("");
    functions1[pt]->Draw("same");
    functions2[pt]->Draw("same");
    bkg2[pt]->SetLineColor(1);
    bkg2[pt]->SetLineStyle(2);
    bkg1[pt]->SetLineColor(1);
    bkg1[pt]->SetLineStyle(2);
    if (isBkgParab)
      bkg2[pt]->Draw("same");
    else
      bkg1[pt]->Draw("same");
    lineP3Sigma[pt]->Draw("same");
    lineM3Sigma[pt]->Draw("same");
  }

  TString Soutputfile;
  Soutputfile = OutputDir + "/Yields_" + Spart[part] + "_" + year;
  Soutputfile += IsOneOrTwoGauss[UseTwoGauss];
  Soutputfile += SIsBkgParab[isBkgParab];
  if (isMB)
    Soutputfile += "_Mult0-100";
  else
    Soutputfile += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  if (isSysStudy)
    Soutputfile += SysPath;
  // Soutputfile += "_Test";
  Soutputfile += "_" + EventType[evFlag];
  Soutputfile += SSysSigExtr[SysSigExtr];

  // save canvases
  canvas[0]->SaveAs(Soutputfile + ".pdf(");
  canvas[1]->SaveAs(Soutputfile + ".pdf");
  canvas[2]->Close();
  canvas[3]->Close();
  // canvas[2]->SaveAs(Soutputfile + ".pdf");
  // canvas[3]->SaveAs(Soutputfile + ".pdf");
  canvasYield->SaveAs(Soutputfile + ".pdf");
  canvasSummary->SaveAs(Soutputfile + ".pdf");
  canvasSummaryBis->SaveAs(Soutputfile + ".pdf)");
  canvasMass->SaveAs(Soutputfile + "_MassPlot.pdf");
  canvasSummary->SaveAs(Soutputfile + "_Nsigma.pdf");

  TFile *outputfile = new TFile(Soutputfile + ".root", "RECREATE");
  outputfile->WriteTObject(histoCandidateSelections);
  for (Int_t i = 0; i < 4; i++)
  {
    outputfile->WriteTObject(canvas[i]);
  }
  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      if (binpt[pt] < 0.8)
        continue;
    }
    hInvMass[pt]->Write();
  }
  outputfile->WriteTObject(histoCountsPerEvent);
  outputfile->WriteTObject(histoYield);
  outputfile->WriteTObject(histoMean);
  outputfile->WriteTObject(histoSigma);
  outputfile->WriteTObject(histoSigmaNarrow);
  outputfile->WriteTObject(histoPurity);
  outputfile->WriteTObject(histoSignificance);
  outputfile->Close();
  cout << "\nA partire dal file:\n"
       << SPathIn << "\n"
       << SPathInEvt << endl;
  cout << "\nHo creato il file: " << Soutputfile << ".root" << endl;

  cout << "Total raw yield (signal only) " << TotYield << endl;
  cout << "Total raw yield (signal+bkg within 3sigmas) " << TotSigBkg << endl;
  cout << "Total number of analysed events " << NEvents << endl;

  TH1F *hEventsFT0MFinal = new TH1F("hEventsFT0MFinal", "hEventsFT0MFinal", 100, 0, 100);
  for (Int_t b = 1; b <= hEventsFT0MFinal->GetNbinsX(); b++)
  {
    counts = 0;
    for (Int_t c = 1; c <= hEventsFT0M->GetNbinsX(); c++)
    {
      if (hEventsFT0M->GetBinCenter(c) > hEventsFT0MFinal->GetXaxis()->GetBinLowEdge(b) && hEventsFT0M->GetBinCenter(c) < hEventsFT0MFinal->GetXaxis()->GetBinUpEdge(b))
      {
        counts += hEventsFT0M->GetBinContent(c);
      }
    }
    hEventsFT0MFinal->SetBinContent(b, counts);
  }

  TH1F *hEventsFV0AFinal = new TH1F("hEventsFV0AFinal", "hEventsFV0AFinal", 100, 0, 100);
  for (Int_t b = 1; b <= hEventsFV0AFinal->GetNbinsX(); b++)
  {
    counts = 0;
    for (Int_t c = 1; c <= hEventsFV0A->GetNbinsX(); c++)
    {
      if (hEventsFV0A->GetBinCenter(c) > hEventsFV0AFinal->GetXaxis()->GetBinLowEdge(b) && hEventsFV0A->GetBinCenter(c) < hEventsFV0AFinal->GetXaxis()->GetBinUpEdge(b))
      {
        counts += hEventsFV0A->GetBinContent(c);
      }
    }
    hEventsFV0AFinal->SetBinContent(b, counts);
  }

  TCanvas *canvasFV0A = new TCanvas("canvasFV0A", "canvasFV0A", 1000, 800);
  StyleCanvas(canvasFV0A, 0.14, 0.05, 0.11, 0.15);
  StyleHisto(hEventsFV0AFinal, 0, 1.2 * hEventsFV0AFinal->GetBinContent(hEventsFV0AFinal->GetMaximumBin()), 1, 1, "FV0A Multiplicity percentile", "Counts", "hEventsFV0AFinal", 0, 0, 0, 1.4, 1.4, 1.2);
  hEventsFV0AFinal->Draw("same");
  canvasFV0A->SaveAs(Soutputfile + "_FV0A.pdf");
  canvasFV0A->SaveAs(Soutputfile + "_FV0A.png");
  // canvasFV0A->Close();

  TCanvas *canvasFT0M = new TCanvas("canvasFT0M", "canvasFT0M", 1000, 800);
  StyleCanvas(canvasFT0M, 0.14, 0.05, 0.11, 0.15);
  StyleHisto(hEventsFT0MFinal, 0, 1.2 * hEventsFT0MFinal->GetBinContent(hEventsFT0MFinal->GetMaximumBin()), 1, 1, "FT0M Multiplicity percentile", "Counts", "hEventsFT0MFinal", 0, 0, 0, 1.4, 1.4, 1.2);
  hEventsFT0MFinal->Draw("same");
  canvasFT0M->SaveAs(Soutputfile + "_FT0M.pdf");
  canvasFT0M->SaveAs(Soutputfile + "_FT0M.png");
  // canvasFT0M->Close();
}
