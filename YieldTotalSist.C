#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TSpline.h>
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include "Constants.h"
#include "ErrRatioCorr.C"
#include "/Users/mbp-cdm-01/Desktop/AssegnoRicerca/Run3Analyses/OmegavsMult/InputVar.h"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(1.5);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->SetTitle(title);
}

void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); // 1.2
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void SetFont(TH1F *histo)
{
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelFont(43);
}
void SetTickLength(TH1F *histo, Float_t TickLengthX, Float_t TickLengthY)
{
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}

void SetHistoTextSize(TH1F *histo, Float_t XSize, Float_t XLabelSize, Float_t XOffset, Float_t XLabelOffset, Float_t YSize, Float_t YLabelSize, Float_t YOffset, Float_t YLabelOffset)
{
  histo->GetXaxis()->SetTitleSize(XSize);
  histo->GetXaxis()->SetLabelSize(XLabelSize);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetXaxis()->SetLabelOffset(XLabelOffset);
  histo->GetYaxis()->SetTitleSize(YSize);
  histo->GetYaxis()->SetLabelSize(YLabelSize);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetYaxis()->SetLabelOffset(YLabelOffset);
}

void StyleCanvas(TCanvas *canvas, Float_t TopMargin, Float_t BottomMargin, Float_t LeftMargin, Float_t RightMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(TopMargin);
  gPad->SetLeftMargin(LeftMargin);
  gPad->SetBottomMargin(BottomMargin);
  gPad->SetRightMargin(RightMargin);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}

void StylePad(TPad *pad, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}

// take spectra in input
// produces ratio of spectra wrt 0-100% multiplciity class

void YieldTotalSist(Bool_t is900GeV = 0,
                    Int_t part = 8,
                    TString PathInSyst = ExtrPathInSyst,
                    TString SysPath = "",
                    TString OutputDir = "SystematicErrors/",
                    TString year = Extryear,
                    Bool_t isSysStudy = 1,
                    Int_t MultType = ExtrMultType, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
                    Bool_t UseTwoGauss = ExtrUseTwoGauss,
                    Bool_t isBkgParab = ExtrisBkgParab,
                    Int_t evFlag = ExtrevFlag // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
)
{

  if (is900GeV)
  {
    year = Extryear900GeV;
    PathInSyst = ExtrPathInSyst900GeV;
  }
  if (part == 3 || part == 4 || part == 5)
    SysPath = ExtrSysPathXi;
  else if (part == 6 || part == 7 || part == 8)
  {
    if (is900GeV)
      SysPath = ExtrSysPathOmega900GeV;
    else
      SysPath = ExtrSysPathOmega;
  }

  gStyle->SetOptStat(0);

  Int_t numMultEff = numMult;
  if (is900GeV)
    numMultEff = numMult900GeV;

  // multiplicity related variables
  TString Smolt[numMultEff + 1];
  TString SmoltBis[numMultEff + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TString PathInSistTopo;
  TString PathInSistMB;
  TString PathInSistSigExt;
  TString PathInSistPileUp;
  TString PathInSistOOBPileUp;
  TString PathInSistMCEff;
  TString PathInSistTotalInput;
  TFile *fileIn[numMultEff + 1];
  TFile *fileInSist[numMultEff + 1];
  TFile *fileInSistTopo[numMultEff + 1];
  TFile *fileInSistMB[numMultEff + 1];
  TFile *fileInSistSigExt[numMultEff + 1];
  TFile *fileInSistPileUp[numMultEff + 1];
  TFile *fileInSistOOBPileUp[numMultEff + 1];
  TFile *fileInSistMCEff[numMultEff + 1];
  TFile *fileInSistTotalInput[numMultEff + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "TotalSysError_" + year + "_";
  stringout += Spart[part];
  if (isSysStudy)
    stringout += SysPath;
  stringout += "_" + EventType[evFlag];
  stringoutpdf = stringout;
  stringout += ".root";

  // canvases
  TCanvas *canvasTotalSistRelErr = new TCanvas("canvasTotalSistRelErr", "canvasTotalSistRelErr", 700, 900);
  StyleCanvas(canvasTotalSistRelErr, 0.05, 0.15, 0.2, 0.02);

  TH1F *fHistSpectrumStat[numMultEff + 1];
  TH1F *fHistSpectrumSist[numMultEff + 1];
  TH1F *fHistSpectrumSistPtCorr[numMultEff + 1];

  TH1F *fHistSpectrumRelErrSist_Topo[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist_MB[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist_SigExt[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist_PileUp[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist_OOBPileUp[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist_MCEff[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist_SignalLoss[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist_AnchRun[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist_TotalInput[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSist[numMultEff + 1];
  TH1F *fHistSpectrumRelErrSistPtCorr[numMultEff + 1];

  TH1F *fHistSpectrumStatScaled[numMultEff + 1];
  TH1F *fHistSpectrumSistScaled[numMultEff + 1];
  TH1F *fHistSpectrumStatScaledB[numMultEff + 1];
  TH1F *fHistSpectrumSistScaledB[numMultEff + 1];
  TH1F *fHistSpectrumSistScaledForLegend[numMultEff + 1];
  TH1F *fHistSpectrumStatMultRatio[numMultEff + 1];
  TH1F *fHistSpectrumSistMultRatio[numMultEff + 1];

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendAllMult = new TLegend(0.22, 0.5, 0.73, 0.7);
  legendAllMult->SetHeader(SMultType[MultType] + " Multiplicity Percentile");
  legendAllMult->SetNColumns(3);
  legendAllMult->SetFillStyle(0);
  TLegendEntry *lheaderAllMult = (TLegendEntry *)legendAllMult->GetListOfPrimitives()->First();
  lheaderAllMult->SetTextSize(0.04);

  TLegend *LegendTitle = new TLegend(0.18, 0.75, 0.4, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(13);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  LegendTitle->AddEntry("", NamePart[part] + ", |y| < 0.5", "");

  TLegend *legendStatBoxK0s = new TLegend(0.25, 0.33, 0.47, 0.42);
  legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[numMultEff], "stat. error", "pe");
  legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[numMultEff], "syst. error", "ef");

  TLine *lineat1Mult = new TLine(0, 1, 8, 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  // get spectra in multiplicity classes
  for (Int_t m = numMultEff; m >= 0; m--)
  {
    PathIn = "Yields/YieldEffCorr" + year + "_";
    PathIn += Spart[part];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[isBkgParab];
    if (m == numMultEff)
    {
      Smolt[m] += "_Mult0-100";
      SmoltBis[m] += "0#minus100";
    }
    else
    {
      if (is900GeV)
      {
        Smolt[m] += Form("_Mult%.1f-%.1f", MultiplicityPerc900GeV[m], MultiplicityPerc900GeV[m + 1]);
        SmoltBis[m] += Form("%.1f#minus%.1f", MultiplicityPerc900GeV[m], MultiplicityPerc900GeV[m + 1]);
      }
      else
      {
        Smolt[m] += Form("_Mult%.1f-%.1f", MultiplicityPerc[m], MultiplicityPerc[m + 1]);
        SmoltBis[m] += Form("%.1f#minus%.1f", MultiplicityPerc[m], MultiplicityPerc[m + 1]);
      }
    }
    PathIn += Smolt[m];
    if (isSysStudy)
      PathIn += SysPath;
    PathIn += "_" + EventType[evFlag];
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;

    fileIn[m] = TFile::Open(PathIn);
    fHistSpectrumStat[m] = (TH1F *)fileIn[m]->Get("histoYieldCorrELCorr");
    fHistSpectrumStat[m]->SetName("histoSpectrumStat" + Smolt[m]);
    if (!fHistSpectrumStat[m])
    {
      cout << " no hist spectrum stat" << endl;
      return;
    }

    // 0. Syst ERROR TOTAL from input file
    PathInSistTotalInput = PathInSyst;
    fileInSistTotalInput[m] = TFile::Open(PathInSistTotalInput);
    fHistSpectrumRelErrSist_TotalInput[m] = (TH1F *)fileInSistTotalInput[m]->Get("hSyst900GeV");
    fHistSpectrumRelErrSist_TotalInput[m]->SetName("histoRelSistTotalInput" + Smolt[m]);

    if (!fHistSpectrumRelErrSist_TotalInput[m])
    {
      cout << " no hist spectrum rel err sist total input" << endl;
      return;
    }

    if (!is900GeV)
    {
      // 1. Syst ERROR TOPOLOGICAL AND KINEMATIC SELECTIONS
      PathInSistTopo = PathInSyst;
      cout << "PathInSistTopo: " << PathInSistTopo << endl;
      fileInSistTopo[m] = TFile::Open(PathInSistTopo);
      fHistSpectrumRelErrSist_Topo[m] = (TH1F *)fileInSistTopo[m]->Get("hSystMultiTrial");
      fHistSpectrumRelErrSist_Topo[m]->SetName("histoRelSistTopo" + Smolt[m]);

      if (!fHistSpectrumRelErrSist_Topo[m])
      {
        cout << " no hist spectrum rel err sist" << endl;
        return;
      }

      // 2. SystErrorMB
      PathInSistMB = PathInSyst;
      cout << "PathInSistMB: " << PathInSistMB << endl;
      fileInSistMB[m] = TFile::Open(PathInSistMB);
      fHistSpectrumRelErrSist_MB[m] = (TH1F *)fileInSistMB[m]->Get("hSystMatBudg");
      fHistSpectrumRelErrSist_MB[m]->SetName("histoRelSistMB" + Smolt[m]);

      if (!fHistSpectrumRelErrSist_MB[m])
      {
        cout << " no hist spectrum rel err sist" << endl;
        return;
      }

      // 3. SysError Signal Extraction
      PathInSistSigExt = PathInSyst;
      cout << "PathInSistSigExt: " << PathInSistSigExt << endl;
      fileInSistSigExt[m] = TFile::Open(PathInSistSigExt);
      fHistSpectrumRelErrSist_SigExt[m] = (TH1F *)fileInSistSigExt[m]->Get("hBkgSyst");
      fHistSpectrumRelErrSist_SigExt[m]->SetName("histoRelSistSigExt" + Smolt[m]);

      if (!fHistSpectrumRelErrSist_SigExt[m])
      {
        cout << " no hist spectrum rel err sist" << endl;
        return;
      }

      // 4. Non-mult dependent efficiency
      PathInSistMCEff = PathInSyst;
      cout << "PathInSistMCEff: " << PathInSistMCEff << endl;
      fileInSistMCEff[m] = TFile::Open(PathInSistMCEff);
      fHistSpectrumRelErrSist_MCEff[m] = (TH1F *)fileInSistMCEff[m]->Get("hSystEffMultDep");
      fHistSpectrumRelErrSist_MCEff[m]->SetName("histoRelSistMCEff" + Smolt[m]);

      if (!fHistSpectrumRelErrSist_MCEff[m])
      {
        cout << " no hist spectrum rel err sist" << endl;
        return;
      }

      // 5. SystError IB PileUp
      // PathInSistPileUp = "SystematicErrors/SysErrorPileUp" + year + "_";
      PathInSistPileUp = PathInSyst;
      cout << "PathInSistPileUp: " << PathInSistPileUp << endl;
      fileInSistPileUp[m] = TFile::Open(PathInSistPileUp);
      fHistSpectrumRelErrSist_PileUp[m] = (TH1F *)fileInSistPileUp[m]->Get("hSystOOB");
      fHistSpectrumRelErrSist_PileUp[m]->SetName("histoRelSistPileUp" + Smolt[m]);

      if (!fHistSpectrumRelErrSist_PileUp[m])
      {
        cout << " no hist spectrum rel err sist" << endl;
        return;
      }
      for (Int_t i = 1; i <= fHistSpectrumRelErrSist_PileUp[m]->GetNbinsX(); i++)
      {
        fHistSpectrumRelErrSist_PileUp[m]->SetBinContent(i, 0);
      }

      // 6. SystError OB PileUp
      // PathInSistOOBPileUp = "SystematicErrors/SysErrorOOBPileUp" + year + "_";
      PathInSistOOBPileUp = PathInSyst;
      cout << "PathInSistOOBPileUp: " << PathInSistOOBPileUp << endl;
      fileInSistOOBPileUp[m] = TFile::Open(PathInSistOOBPileUp);
      fHistSpectrumRelErrSist_OOBPileUp[m] = (TH1F *)fileInSistOOBPileUp[m]->Get("hSystOOB");
      fHistSpectrumRelErrSist_OOBPileUp[m]->SetName("histoRelSistOOBPileUp" + Smolt[m]);

      if (!fHistSpectrumRelErrSist_OOBPileUp[m])
      {
        cout << " no hist spectrum rel err sist" << endl;
        return;
      }
      for (Int_t i = 1; i <= fHistSpectrumRelErrSist_OOBPileUp[m]->GetNbinsX(); i++)
      {
        // fHistSpectrumRelErrSist_OOBPileUp[m]->SetBinContent(i, 0);
      }
    }

    // 7. SystError mult dependence of signal loss
    Float_t RelErrSignalLoss = 0.025;
    if (is900GeV)
      RelErrSignalLoss = 0.025;
    fHistSpectrumRelErrSist_SignalLoss[m] = (TH1F *)fHistSpectrumStat[m]->Clone("hSigLossSyst");
    for (Int_t i = 1; i <= fHistSpectrumRelErrSist_SignalLoss[m]->GetNbinsX(); i++)
    {
      fHistSpectrumRelErrSist_SignalLoss[m]->SetBinContent(i, RelErrSignalLoss);
      fHistSpectrumRelErrSist_SignalLoss[m]->SetBinError(i, 0);
    }

    // 7. SystError anchoring run
    Float_t RelErrAnch = 0.04;
    if (is900GeV)
      RelErrAnch = 0.1;
    fHistSpectrumRelErrSist_AnchRun[m] = (TH1F *)fHistSpectrumStat[m]->Clone("hSigLossSyst");
    for (Int_t i = 1; i <= fHistSpectrumRelErrSist_AnchRun[m]->GetNbinsX(); i++)
    {
      fHistSpectrumRelErrSist_AnchRun[m]->SetBinContent(i, RelErrAnch);
      fHistSpectrumRelErrSist_AnchRun[m]->SetBinError(i, 0);
    }

    Int_t Bin_Topo = 0;
    Int_t Bin_MB = 0;
    Int_t Bin_SigExt = 0;
    Int_t Bin_PileUp = 0;
    Int_t Bin_OOBPileUp = 0;
    Int_t Bin_SignalLoss = 0;
    Int_t Bin_MCEff = 0;
    Int_t Bin_AnchRun = 0;
    Int_t Bin_TotalInput = 0;
    Int_t MinBinPtEff = MinBinPt[part];
    if (is900GeV)
      MinBinPtEff = MinBinPt900GeV[part];
    //*********** SUM IN QUADRATURE OF ALL SOURCES OF UNCERTAINTY *************//
    fHistSpectrumRelErrSist[m] = (TH1F *)fHistSpectrumStat[m]->Clone("histoSpectrumRelErrSist" + Smolt[m]);
    fHistSpectrumRelErrSistPtCorr[m] = (TH1F *)fHistSpectrumStat[m]->Clone("histoSpectrumRelErrSistPtCorr" + Smolt[m]);
    for (Int_t i = 1; i <= fHistSpectrumRelErrSist[m]->GetNbinsX(); i++)
    {
      if ((fHistSpectrumRelErrSist[m]->GetXaxis()->GetBinLowEdge(i)) < MinBinPtEff)
      {
        continue;
      }
      cout << "bin center: " << fHistSpectrumRelErrSist[m]->GetBinCenter(i) << endl;
      if (!is900GeV)
      {
        Bin_Topo = fHistSpectrumRelErrSist_Topo[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));
        Bin_MB = fHistSpectrumRelErrSist_MB[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));
        Bin_SigExt = fHistSpectrumRelErrSist_SigExt[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));
        Bin_PileUp = fHistSpectrumRelErrSist_PileUp[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));
        Bin_OOBPileUp = fHistSpectrumRelErrSist_OOBPileUp[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));
        Bin_MCEff = fHistSpectrumRelErrSist_MCEff[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));
      }
      Bin_SignalLoss = fHistSpectrumRelErrSist_SignalLoss[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));
      Bin_AnchRun = fHistSpectrumRelErrSist_AnchRun[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));
      Bin_TotalInput = fHistSpectrumRelErrSist_TotalInput[m]->FindBin(fHistSpectrumRelErrSist[m]->GetBinCenter(i));

      if (is900GeV)
      {
        fHistSpectrumRelErrSist[m]->SetBinContent(i, TMath::Sqrt(
                                                         TMath::Power(fHistSpectrumRelErrSist_SignalLoss[m]->GetBinContent(Bin_SignalLoss), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_TotalInput[m]->GetBinContent(Bin_MCEff), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_AnchRun[m]->GetBinContent(Bin_AnchRun), 2)));
        fHistSpectrumRelErrSistPtCorr[m]->SetBinContent(i, TMath::Sqrt(
                                                               TMath::Power(fHistSpectrumRelErrSist_SignalLoss[m]->GetBinContent(Bin_SignalLoss), 2) +
                                                               TMath::Power(fHistSpectrumRelErrSist_TotalInput[m]->GetBinContent(Bin_MCEff), 2) +
                                                               TMath::Power(fHistSpectrumRelErrSist_AnchRun[m]->GetBinContent(Bin_AnchRun), 2)));
      }
      else
      {
        fHistSpectrumRelErrSist[m]->SetBinContent(i, TMath::Sqrt(
                                                         TMath::Power(fHistSpectrumRelErrSist_Topo[m]->GetBinContent(Bin_Topo), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_MB[m]->GetBinContent(Bin_MB), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_SigExt[m]->GetBinContent(Bin_SigExt), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_PileUp[m]->GetBinContent(Bin_PileUp), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_OOBPileUp[m]->GetBinContent(Bin_OOBPileUp), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_SignalLoss[m]->GetBinContent(Bin_SignalLoss), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_MCEff[m]->GetBinContent(Bin_MCEff), 2) +
                                                         TMath::Power(fHistSpectrumRelErrSist_AnchRun[m]->GetBinContent(Bin_AnchRun), 2)));
        fHistSpectrumRelErrSistPtCorr[m]->SetBinContent(i, TMath::Sqrt(
                                                               TMath::Power(fHistSpectrumRelErrSist_MB[m]->GetBinContent(Bin_MB), 2) +
                                                               TMath::Power(fHistSpectrumRelErrSist_SigExt[m]->GetBinContent(Bin_SigExt), 2) +
                                                               TMath::Power(fHistSpectrumRelErrSist_PileUp[m]->GetBinContent(Bin_PileUp), 2) +
                                                               TMath::Power(fHistSpectrumRelErrSist_OOBPileUp[m]->GetBinContent(Bin_OOBPileUp), 2) +
                                                               TMath::Power(fHistSpectrumRelErrSist_SignalLoss[m]->GetBinContent(Bin_SignalLoss), 2) +
                                                               TMath::Power(fHistSpectrumRelErrSist_MCEff[m]->GetBinContent(Bin_MCEff), 2) +
                                                               TMath::Power(fHistSpectrumRelErrSist_AnchRun[m]->GetBinContent(Bin_AnchRun), 2)));
      }
      cout << "stat. " << fHistSpectrumStat[m]->GetBinError(i) << endl;
      cout << "syst. " << fHistSpectrumRelErrSist[m]->GetBinContent(i) << endl;
      cout << "syst. pt corr " << fHistSpectrumRelErrSistPtCorr[m]->GetBinContent(i) << endl;
    }

    fHistSpectrumSist[m] = (TH1F *)fHistSpectrumStat[m]->Clone("histoSpectrumSist" + Smolt[m]);
    fHistSpectrumSistPtCorr[m] = (TH1F *)fHistSpectrumStat[m]->Clone("histoSpectrumSistPtCorr" + Smolt[m]);
    for (Int_t i = 1; i <= fHistSpectrumSist[m]->GetNbinsX(); i++)
    {
      fHistSpectrumSist[m]->SetBinError(i, fHistSpectrumStat[m]->GetBinContent(i) * fHistSpectrumRelErrSist[m]->GetBinContent(fHistSpectrumRelErrSist[m]->FindBin(fHistSpectrumSist[m]->GetBinCenter(i))));
      fHistSpectrumSistPtCorr[m]->SetBinError(i, fHistSpectrumStat[m]->GetBinContent(i) * fHistSpectrumRelErrSistPtCorr[m]->GetBinContent(fHistSpectrumRelErrSistPtCorr[m]->FindBin(fHistSpectrumSist[m]->GetBinCenter(i))));
    }
  } // end loop on mult

  // draw total syst. rel errors in multiplicity classes
  Float_t LimSupSpectra = 1;
  Float_t LimInfSpectra = 0.2 * 1e-3;
  Float_t xTitle = 30;
  Float_t xOffset = 2;
  Float_t yTitle = 30;
  Float_t yOffset = 2;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.01;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.042;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 8);
  // for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)  hDummy->SetBinContent(i, 1e-12);
  canvasTotalSistRelErr->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, LimInfSpectra, LimSupSpectra, 1, 1, TitleXPt, "Rel. syst. uncertainty", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, 8);
  hDummy->Draw("same");

  for (Int_t m = numMultEff; m >= 0; m--)
  {
    if (m == numMultEff)
    {
      ColorMult[m] = ColorMB;
      MarkerMult[m] = MarkerMB;
      SizeMult[m] = SizeMB;
    }
    fHistSpectrumRelErrSist[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumRelErrSist[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumRelErrSist[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumRelErrSist[m]->SetMarkerSize(SizeMult[m]);
    fHistSpectrumRelErrSist[m]->SetLineStyle(10);
    fHistSpectrumRelErrSist[m]->Draw("same h");
    fHistSpectrumRelErrSistPtCorr[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumRelErrSistPtCorr[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumRelErrSistPtCorr[m]->SetMarkerStyle(22);
    fHistSpectrumRelErrSistPtCorr[m]->SetMarkerSize(SizeMult[m]);
    fHistSpectrumRelErrSistPtCorr[m]->Draw("same h");
    for (Int_t i = 1; i <= fHistSpectrumRelErrSist[m]->GetNbinsX(); i++)
    {
      cout << "bin center: " << fHistSpectrumRelErrSist[m]->GetBinCenter(i) << endl;
      cout << "syst. " << fHistSpectrumRelErrSist[m]->GetBinContent(i) << endl;
      cout << "syst. pt corr " << fHistSpectrumRelErrSistPtCorr[m]->GetBinContent(i) << endl;
    }
    legendAllMult->AddEntry(fHistSpectrumRelErrSist[m], SmoltBis[m] + "%", "pef");
  } // end loop on mult
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  TFile *fileout = new TFile(stringout, "RECREATE");

  for (Int_t m = numMultEff; m >= 0; m--)
  {
    fHistSpectrumStat[m]->Write();
    fHistSpectrumRelErrSist[m]->Write();
    fHistSpectrumSist[m]->Write();
    fHistSpectrumRelErrSistPtCorr[m]->Write();
    fHistSpectrumSistPtCorr[m]->Write();
  }

  fileout->Close();
  canvasTotalSistRelErr->SaveAs(stringoutpdf + ".pdf");
  canvasTotalSistRelErr->SaveAs(stringoutpdf + ".png");
  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
  cout << "WARNING!! file of syst error is always the 0-100% one" << endl;
}
