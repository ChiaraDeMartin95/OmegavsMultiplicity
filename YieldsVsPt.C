#include <Riostream.h>
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

Bool_t reject;
Double_t fparab(Double_t *x, Double_t *par)
{
  Float_t LimInf = 0;
  Float_t LimSup = 0;
  if (par[3] == 0)
  {
    LimInf = 0.474;
    LimSup = 0.520;
  }
  else if (par[3] == 3 || par[3] == 4)
  {
    LimInf = 1.310;
    LimSup = 1.335;
  }
  else if (par[2] == 5 || par[2] == 6)
  {
    LimInf = 1.66;
    LimSup = 1.682;
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
  else if (par[2] == 3 || par[2] == 3)
  {
    LimInf = 1.310;
    LimSup = 1.335;
  }
  else if (par[2] == 5 || par[2] == 6)
  {
    LimInf = 1.66;  // 1.65
    LimSup = 1.682; // 1.69
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

TString TitleInvMass[numPart] = {"(#pi^{+}, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(p, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(#bar{p}, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(#Lambda, #pi^{-}) invariant mass (GeV/#it{c}^{2})"};
TString namehisto[numPart] = {"h3dMassK0Short", "", "", "hCascMinusInvMassvsPt", "hCascPlusInvMassvsPt", "hCascMinusInvMassvsPt", "hCascPlusInvMassvsPt"};
Float_t LowLimitMass[numPart] = {0.42, 1.09, 1.09, 1.29, 1.29, 1.61, 1.61}; // 0.44
Float_t UpLimitMass[numPart] = {0.57, 1.14, 1.14, 1.35, 1.35, 1.73, 1.73};  // 0.55
Float_t LowMassRange[numPart] = {0.48, 1.09, 1.09, 1.31, 1.31, 1.655, 1.655};
Float_t UpMassRange[numPart] = {0.51, 1.14, 1.14, 1.33, 1.33, 1.685, 1.685};

Float_t min_range_signal[numPart] = {0.46, 1.105, 1.105, 1.31, 1.31, 1.65, 1.65}; // estremi region fit segnale (gaussiane)
Float_t max_range_signal[numPart] = {0.535, 1.125, 1.125, 1.334, 1.334, 1.69, 1.69};
Float_t min_histo[numPart] = {0.42, 1.09, 1.09, 1.29, 1.29, 1.62, 1.62}; // estremi del range visualizzazione bkg
Float_t max_histo[numPart] = {0.57, 1.14, 1.14, 1.35, 1.35, 1.72, 1.72};
Float_t liminf[numPart] = {0.45, 1.1153, 1.1153, 1.29, 1.29, 1.63, 1.63}; // estremi regione fit del bkg e total
Float_t limsup[numPart] = {0.545, 1.1168, 1.1168, 1.35, 1.35, 1.71, 1.71};

// histogram limits
// Omega 22o_pass3: 1.62-1.72
// Omega 22m_pass4: 1.64 - 1.69

void YieldsVsPt(Bool_t isSysStudy = 1,
                TString SysPath = "",
                Int_t part = 5,
                Int_t MassRebin = 2,
                Int_t MultType = 0, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
                Bool_t isMB = 1,
                Int_t mul = 0,
                TString year = "LHC22m_pass4_Train79153" /*"LHC22o_pass3_Train75538" /*"LHC22r_pass3_Train67853"*/,
                TString SPathIn = "OutputFilesCascPPTask/LHC22m_pass4/AnalysisResults",
                TString SPathInEvt = "AnalysisResultsCascQATask/LHC22m_pass4/AnalysisResultsEvts",
                TString OutputDir = "Yields",
                Bool_t UseTwoGauss = 0,
                Bool_t isBkgParab = 0,
                Bool_t isMeanFixedPDG = 0,
                Float_t sigmacentral = 3)
{

  if (mul > numMult)
  {
    cout << "Multiplciity out of range" << endl;
    return;
  }
  if (part < 3)
  {
    cout << "this value of part is not specified, choose 3 - 4 (Xi) or 5 - 6  (Omega) " << endl;
    return;
  }
  if (part == 5 || part == 6)
  { // Omega
    if (year == "LHC22m_pass4_Train79153")
    {
      LowLimitMass[part] = 1.63;
      UpLimitMass[part] = 1.71;
    }
  }

  SPathIn += "_" + year + "_" + SpartType[part];
  if (isSysStudy)
    SPathIn += SysPath;
  SPathIn += ".root";

  SPathInEvt += "_" + year + ".root";
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
  TH2F *h2;
  TH2F *h2Bis;
  TH1F *hEvents;
  TH1F *hEventsFT0M;
  TH1F *hEventsFV0A;
  Int_t MultLowBin = 0;
  Int_t MultUpBin = 0;

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
      h3 = (TH3F *)dir->Get(namehisto[part] + "_FT0M");
    else if (MultType == 2)
      h3 = (TH3F *)dir->Get(namehisto[part] + "_FV0A");
    if (!h3)
    {
      cout << "h3 cascade not avilable " << endl;
      return;
    }
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
  hEventsFT0M = (TH1F *)dirEvt->Get("hCentFT0M");
  if (!hEventsFT0M)
  {
    cout << "hCentFT0M not available " << endl;
    return;
  }
  hEventsFV0A = (TH1F *)dirEvt->Get("hCentFV0A");
  if (!hEventsFV0A)
  {
    cout << "hCentFV0A not available " << endl;
    return;
  }
  if (isMB == 1)
    NEvents = hEvents->GetBinContent(1);
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

  const Int_t numPt = 19; // default: 19; topo: 11;

  // Xi
  //  Float_t binpt[numPt + 1] = {0.4, 0.6, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3, 4}; // 10
  //  Float_t binpt[numPt + 1] = {0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0};
  // Float_t binpt[numPt + 1] = {0.6, 1.0, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 4.0};

  // Float_t binpt[numPt + 1] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
  // 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
  // 2.2, 2.6, 3.0, 3.5, 4.0, 4.5, 5.0};

  // Omega 22o_pass3
  // Float_t binpt[numPt + 1] = {0.4, 0.6, 0.8, 1.0,
  //                          1.2, 1.4, 1.6, 1.8, 2.0,
  //                          2.5, 3.5, 5.0};

  // Omega 22m_pass4
  // Float_t binpt[numPt + 1] = {0.4, 0.6, 0.8, 1.0, 1.1,
  //                          1.2, 1.3, 1.4, 1.5, 1.6,
  //                        1.7, 1.8, 1.9, 2.0, 2.1,
  //                      2.2, 2.3, 2.4, 2.5, 2.7,
  //                    2.9, 3.1, 3.3, 3.5, 4.0,
  //                  4.5, 5.0, 6.0, 8.0};
  // 28

  // 19
  // default
  Float_t binpt[numPt + 1] = {0.4, 0.6, 0.8, 1.0,
                              1.2, 1.4, 1.6, 1.8, 2.0,
                              2.2, 2.4, 2.6, 2.8, 3.0,
                              3.5, 4.0, 4.5, 5.0, 6.0, 8.0};

  // for topo studies
  // Float_t binpt[numPt + 1] = {0.4, 0.8,
  //                          1.2, 1.6, 2.0,
  //                      2.4, 2.8, 3.2, 4.0, 5.0, 6.0, 8.0};

  TString SPt[numPt] = {""};
  TH1F *hInvMass[numPt];

  TCanvas *canvas[4];
  for (Int_t c = 0; c < 4; c++)
  {
    canvas[c] = new TCanvas(Form("canvas_%i", c), Form("canvas%i", c), 1800, 1400);
    canvas[c]->Divide(3, 3);
    StyleCanvas(canvas[c], 0.15, 0.05, 0.05, 0.15);
  }

  TH1F *histoCountsPerEvent = new TH1F("histoCountsPerEvent", "histoCountsPerEvent", numPt, binpt);
  TH1F *histoMean = new TH1F("histoMean", "histoMean", numPt, binpt);
  TH1F *histoSigma = new TH1F("histoSigma", "histoSigma", numPt, binpt);
  TH1F *histoPurity = new TH1F("histoPurity", "histoPurity", numPt, binpt);
  TH1F *histoYield = new TH1F("histoYield", "histoYield", numPt, binpt);
  TH1F *histoSignificance = new TH1F("histoSignificance", "histoSignificance", numPt, binpt);

  Float_t counts = 0;
  Float_t errcount = 0;
  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 5 || part == 6)
    {
      if (binpt[pt] < 0.6)
        continue;
      // if (binpt[pt] > 4)
      // continue;
    }
    SPt[pt] = Form("%.1f < p_{T} < %.1f", binpt[pt], binpt[pt + 1]);
    cout << "Analysed pt interval: " << binpt[pt] << "-" << binpt[pt + 1] << endl;
    cout << binpt[pt] << endl;

    hInvMass[pt] = (TH1F *)h2->ProjectionY(Form("hInvMass_pt%i", pt), h2->GetXaxis()->FindBin(binpt[pt] + 0.001), h2->GetXaxis()->FindBin(binpt[pt + 1] - 0.001));
    hInvMass[pt]->Rebin(MassRebin);
    StyleHisto(hInvMass[pt], 0, 1.2 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), 1, 20, TitleInvMass[part], "Counts", SPt[pt], 1, LowLimitMass[part], UpLimitMass[part], 1.4, 1.4, 0.7);
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
  Float_t b[numPt] = {0};
  Float_t errb[numPt] = {0};
  Float_t SSB[numPt] = {0};
  Float_t errSSB[numPt] = {0};
  Float_t entries_range[numPt] = {0};
  Float_t Yield[numPt] = {0};
  Float_t ErrYield[numPt] = {0};
  Float_t TotYield = 0;
  Float_t TotSigBkg = 0;

  for (Int_t pt = 0; pt < numPt; pt++)
  {
    if (part == 5 || part == 6)
    { // Omega
      if (year == "LHC22m_pass4_Train79153")
      {
        min_histo[part] = 1.645; // estremi del range visualizzazione bkg
        max_histo[part] = 1.695;
        liminf[part] = 1.645; // estremi regione fit del bkg e total
        limsup[part] = 1.695;
        if (binpt[pt] > 2)
        {
          min_histo[part] = 1.64;
          liminf[part] = 1.64;
          max_histo[part] = 1.7;
          limsup[part] = 1.7;
        }
        UpLimitMass[part] = 1.71;
        if (binpt[pt] < 1.2)
        {
          min_histo[part] = 1.635;
          liminf[part] = 1.635;
          max_histo[part] = 1.7;
          limsup[part] = 1.7;
        }
        // retta or pol2
        if (binpt[pt] < 1. || binpt[pt] > 4.)
          isBkgParab = 0;
        else
          isBkgParab = 1;
        if (part == 5 || part == 6)
        {
          // if (binpt[pt] == 1.4 || binpt[pt] == 2. || binpt[pt] == 2.2) isBkgParab = 0;
          if (binpt[pt] >= 1.99 && binpt[pt] <= 2.1)
            isBkgParab = 0;
          if (binpt[pt] >= 1.59 && binpt[pt] <= 1.61)
            isBkgParab = 0;
        }
        if (part == 6)
        {
          if (binpt[pt] >= 0.9 && binpt[pt] <= 2.1)
            isBkgParab = 0;
        }
      }
    }
    if (part == 5 || part == 6)
    {
      if (binpt[pt] < 0.6)
        continue;
    }
    // special settings
    if (part == 0)
    {
      if (pt == 1)
        limsup[part] = 0.54;
      if (binpt[pt] > 1.5)
      {
        liminf[part] = 0.42;
        limsup[part] = 0.57;
      }
    }
    // end of special settings

    // UseTwoGauss = 1;
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

    functions1[pt] = new TF1(Form("1f_%i_final", pt), "gaus", LowLimitMass[part], UpLimitMass[part]);
    functions1[pt]->SetLineColor(kRed); // 867
    functions1[pt]->SetParName(0, "norm");
    functions1[pt]->SetParName(1, "mean");
    functions1[pt]->SetParName(2, "sigma");

    functions2[pt] = new TF1(Form("2f_%i_final", pt), "gaus", LowLimitMass[part], UpLimitMass[part]);
    functions2[pt]->SetLineColor(kMagenta); // 891
    functions2[pt]->SetParName(0, "norm");
    functions2[pt]->SetParName(1, "mean");
    functions2[pt]->SetParName(2, "sigma");

    bkg1[pt] = new TF1(Form("bkg1%i", pt), "pol1", min_histo[part], max_histo[part]);
    bkg1[pt]->SetLineColor(418);

    bkg2[pt] = new TF1(Form("bkg2%i", pt), "pol2", min_histo[part], max_histo[part]);
    bkg2[pt]->SetLineColor(kBlue);

    bkgretta[pt] = new TF1(Form("retta%i", pt), fretta, liminf[part], limsup[part], 3);
    bkgretta[pt]->SetLineColor(418);
    bkgretta[pt]->FixParameter(2, part);

    bkgparab[pt] = new TF1(Form("parab%i", pt), fparab, liminf[part], limsup[part], 4);
    bkgparab[pt]->SetLineColor(kBlue);
    bkgparab[pt]->FixParameter(3, part);

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

      bkg1[pt]->SetRange(min_histo[part], max_histo[part]);
      bkg2[pt]->SetRange(min_histo[part], max_histo[part]);
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
        total[pt]->SetParLimits(1, 1.318, 1.326);
        total[pt]->SetParLimits(2, 0.0012, 0.010);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.318, 1.326);
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
        total[pt]->SetParLimits(2, 0.002, 0.02);
        total[pt]->SetParLimits(3, 0.08 * hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin()), hInvMass[pt]->GetBinContent(hInvMass[pt]->GetMaximumBin())); // maximum was wothout 0.3
        total[pt]->SetParLimits(4, 1.66, 1.68);
        total[pt]->SetParLimits(5, 0.001, 0.02);
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
          UseTwoGauss = kFALSE;
      }
      else
      {
        if (total[pt]->GetParameter(0) < total[pt]->GetParameter(3))
          UseTwoGauss = kFALSE;
      }

      cout << "UseTwoGauss = " << UseTwoGauss << endl;

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

      if (part == 3)
        UseTwoGauss = 0;

      if (UseTwoGauss)
      {
        if (pt < 9)
          canvas[0]->cd(pt + 1);
        else if (pt < 18)
          canvas[1]->cd(pt + 1 - 9);
        else if (pt < 27)
          canvas[2]->cd(pt + 1 - 18);
        else
          canvas[3]->cd(pt + 1 - 27);
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
        errsigma[pt] = sqrt(pow(total[pt]->GetParError(2), 2) + pow(total[pt]->GetParError(5), 2) + 2 * cov_sigma) / 2;
      }
    }
    if (!UseTwoGauss)
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

      bkg1[pt]->SetRange(min_histo[part], max_histo[part]);
      bkg2[pt]->SetRange(min_histo[part], max_histo[part]);
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
        total[pt]->SetParLimits(1, 1.318, 1.326);
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
      if (isBkgParab)
        bkg2[pt]->Draw("same");
      else
        bkg1[pt]->Draw("same");
      functions1[pt]->Draw("same");

      mean[pt] = total[pt]->GetParameter(1);
      errmean[pt] = total[pt]->GetParError(1);
      sigma[pt] = total[pt]->GetParameter(2);
      errsigma[pt] = total[pt]->GetParError(2);
    }

    // cout << mean[pt] - sigmacentral * sigma[pt] << "-" << mean[pt] + sigmacentral * sigma[pt] << endl;
    //  Compute yield
    b[pt] = 0;
    errb[pt] = 0;
    if (isBkgParab)
    {
      b[pt] = bkg2[pt]->Integral(mean[pt] - sigmacentral * sigma[pt], mean[pt] + sigmacentral * sigma[pt]);
      errb[pt] = totalbis[pt]->IntegralError(mean[pt] - sigmacentral * sigma[pt], mean[pt] + sigmacentral * sigma[pt], fFitResultPtr1[pt]->GetParams(),
                                             (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
    }
    else
    {
      b[pt] = bkg1[pt]->Integral(mean[pt] - sigmacentral * sigma[pt], mean[pt] + sigmacentral * sigma[pt]);
      errb[pt] = totalbis[pt]->IntegralError(mean[pt] - sigmacentral * sigma[pt], mean[pt] + sigmacentral * sigma[pt], fFitResultPtr1[pt]->GetParams(),
                                             (fFitResultPtr1[pt]->GetCovarianceMatrix()).GetMatrixArray());
    }
    b[pt] = b[pt] / hInvMass[pt]->GetBinWidth(1);
    errb[pt] = errb[pt] / hInvMass[pt]->GetBinWidth(1);

    entries_range[pt] = 0;
    for (Int_t l = hInvMass[pt]->GetXaxis()->FindBin(mean[pt] - sigmacentral * sigma[pt]); l <= hInvMass[pt]->GetXaxis()->FindBin(mean[pt] + sigmacentral * sigma[pt]); l++)
    {
      entries_range[pt] += hInvMass[pt]->GetBinContent(l);
    }

    Yield[pt] = entries_range[pt] - b[pt];
    ErrYield[pt] = sqrt(entries_range[pt] + pow(errb[pt], 2));
    TotYield += Yield[pt];
    TotSigBkg += entries_range[pt];

    SSB[pt] = (entries_range[pt] - b[pt]) / entries_range[pt];
    errSSB[pt] = SSB[pt] * sqrt(1. / entries_range[pt] + pow(errb[pt] / b[pt], 2));

    histoYield->SetBinContent(pt + 1, Yield[pt] / NEvents / histoYield->GetBinWidth(pt + 1));
    histoYield->SetBinError(pt + 1, ErrYield[pt] / NEvents / histoCountsPerEvent->GetBinWidth(pt + 1));

    histoMean->SetBinContent(pt + 1, mean[pt]);
    histoMean->SetBinError(pt + 1, errmean[pt]);

    histoSigma->SetBinContent(pt + 1, sigma[pt]);
    histoSigma->SetBinError(pt + 1, errsigma[pt]);

    histoPurity->SetBinContent(pt + 1, SSB[pt]);
    histoPurity->SetBinError(pt + 1, errSSB[pt]);

    histoSignificance->SetBinContent(pt + 1, Yield[pt] / ErrYield[pt]);
    histoSignificance->SetBinError(pt + 1, 0);
  }

  TotYield = TotYield / NEvents;
  TotSigBkg = TotSigBkg / NEvents;
  histoYield->SetLineColor(kRed);

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 1000, 800);
  if (part == 5 || part == 6)
    histoCountsPerEvent->GetYaxis()->SetRangeUser(0, 0.00007);
  histoCountsPerEvent->Draw("same");
  histoYield->DrawClone("same");
  TLegend *legendYield = new TLegend(0.7, 0.7, 0.9, 0.9);
  legendYield->AddEntry(histoCountsPerEvent, "w/o bkg subtraction", "pl");
  legendYield->AddEntry(histoYield, "w/ bkg subtraction", "pl");
  legendYield->Draw("");

  TCanvas *canvasSummary = new TCanvas("canvasSummary", "canvasSummary", 1000, 800);
  canvasSummary->Divide(2, 2);
  TCanvas *canvasSummaryBis = new TCanvas("canvasSummaryBis", "canvasSummaryBis", 1000, 800);
  canvasSummaryBis->Divide(2, 2);

  canvasSummary->cd(1);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoMean, LowLimitMass[part], UpLimitMass[part], 1, 1, titlePt, "#mu (GeV/c^{2})", "histoMean", 0, 0, 0, 1.4, 1.4, 1.2);
  histoMean->Draw("");
  canvasSummary->cd(2);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoSigma, 0, 0.010, 1, 1, titlePt, "#sigma (GeV/c^{2})", "histoSigma", 0, 0, 0, 1.4, 1.4, 1.2);
  if (part >= 5)
    histoSigma->GetYaxis()->SetRangeUser(0, 0.02);
  histoSigma->Draw("");
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
  canvasSummaryBis->cd(1);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.14);
  StyleHisto(histoSignificance, 0, 1.2 * histoSignificance->GetBinContent(histoSignificance->GetMaximumBin()), 1, 1, titlePt, "Yield/#sigma_{Yield}", "histoSignificance", 0, 0, 0, 1.4, 1.4, 1.2);
  histoSignificance->Draw("same");

  TString Soutputfile;
  Soutputfile = OutputDir + "/Yields_" + Spart[part] + "_" + year;
  Soutputfile += IsOneOrTwoGauss[UseTwoGauss];
  // Soutputfile += SIsBkgParab[isBkgParab];
  if (isMB)
    Soutputfile += "_Mult0-100";
  else
    Soutputfile += Form("_Mult%i-%i", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  // Soutputfile += "_DefSel";
  if (isSysStudy)
    Soutputfile += SysPath;
  Soutputfile += "_FewPtBins";

  // save canvases
  canvas[0]->SaveAs(Soutputfile + ".pdf(");
  canvas[1]->SaveAs(Soutputfile + ".pdf");
  canvas[2]->SaveAs(Soutputfile + ".pdf");
  canvas[3]->SaveAs(Soutputfile + ".pdf");
  canvasYield->SaveAs(Soutputfile + ".pdf");
  canvasSummary->SaveAs(Soutputfile + ".pdf");
  canvasSummaryBis->SaveAs(Soutputfile + ".pdf)");

  TFile *outputfile = new TFile(Soutputfile + ".root", "RECREATE");
  outputfile->WriteTObject(histoCandidateSelections);
  for (Int_t i = 0; i < 4; i++)
  {
    outputfile->WriteTObject(canvas[i]);
  }
  outputfile->WriteTObject(histoCountsPerEvent);
  outputfile->WriteTObject(histoYield);
  outputfile->WriteTObject(histoMean);
  outputfile->WriteTObject(histoSigma);
  outputfile->WriteTObject(histoPurity);
  outputfile->WriteTObject(histoSignificance);
  outputfile->Close();
  cout << "\nA partire dal file:\n"
       << SPathIn << "\n"
       << SPathInEvt << endl;
  cout << "\nHo creato il file: " << Soutputfile << ".root" << endl;

  cout << "Total raw yield (signal only) " << TotYield << endl;
  cout << "Total raw yield (signal+bkg within 3sigmas) " << TotSigBkg << endl;
}
