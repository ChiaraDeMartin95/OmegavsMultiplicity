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

void YieldTotalSist(Int_t part = 6,
                    Int_t ChosenMult = numMult,
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

  if (part == 3 || part == 4 || part == 5)
    SysPath = ExtrSysPathXi;
  else if (part == 6 || part == 7 || part == 8)
    SysPath = ExtrSysPathOmega;

  gStyle->SetOptStat(0);
  if (ChosenMult > numMult)
  {
    cout << "Chosen Mult outside of available range" << endl;
    return;
  }

  // multiplicity related variables
  TString Smolt[numMult + 1];
  TString SmoltBis[numMult + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TString PathInSistTopo;
  TString PathInSistMB;
  TString PathInSistSigExt;
  TString PathInSistPileUp;
  TString PathInSistOOBPileUp;
  TString PathInSistMCEff;
  TFile *fileIn[numMult + 1];
  TFile *fileInSist[numMult + 1];
  TFile *fileInSistTopo[numMult + 1];
  TFile *fileInSistMB[numMult + 1];
  TFile *fileInSistSigExt[numMult + 1];
  TFile *fileInSistPileUp[numMult + 1];
  TFile *fileInSistOOBPileUp[numMult + 1];
  TFile *fileInSistMCEff[numMult + 1];

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

  TH1F *fHistSpectrumStat[numMult + 1];
  TH1F *fHistSpectrumSist[numMult + 1];
  TH1F *fHistSpectrumSistPtCorr[numMult + 1];

  TH1F *fHistSpectrumRelErrSist_Topo[numMult + 1];
  TH1F *fHistSpectrumRelErrSist_MB[numMult + 1];
  TH1F *fHistSpectrumRelErrSist_SigExt[numMult + 1];
  TH1F *fHistSpectrumRelErrSist_PileUp[numMult + 1];
  TH1F *fHistSpectrumRelErrSist_OOBPileUp[numMult + 1];
  TH1F *fHistSpectrumRelErrSist_MCEff[numMult + 1];
  TH1F *fHistSpectrumRelErrSist_SignalLoss[numMult + 1];
  TH1F *fHistSpectrumRelErrSist[numMult + 1];
  TH1F *fHistSpectrumRelErrSistPtCorr[numMult + 1];

  TH1F *fHistSpectrumStatScaled[numMult + 1];
  TH1F *fHistSpectrumSistScaled[numMult + 1];
  TH1F *fHistSpectrumStatScaledB[numMult + 1];
  TH1F *fHistSpectrumSistScaledB[numMult + 1];
  TH1F *fHistSpectrumSistScaledForLegend[numMult + 1];
  TH1F *fHistSpectrumStatMultRatio[numMult + 1];
  TH1F *fHistSpectrumSistMultRatio[numMult + 1];

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
  legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[numMult], "stat. error", "pe");
  legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[numMult], "syst. error", "ef");

  TLine *lineat1Mult = new TLine(0, 1, 8, 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  // get spectra in multiplicity classes
  for (Int_t m = numMult; m >= 0; m--)
  {
    PathIn = "Yields/YieldEffCorr" + year + "_";
    PathIn += Spart[part];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[isBkgParab];
    if (m == numMult)
    {
      Smolt[m] += "_Mult0-100";
      SmoltBis[m] += "0#minus100";
    }
    else
    {
      Smolt[m] += Form("_Mult%.1f-%.1f", MultiplicityPerc[m], MultiplicityPerc[m + 1]);
      SmoltBis[m] += Form("%.1f#minus%.1f", MultiplicityPerc[m], MultiplicityPerc[m + 1]);
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
      //fHistSpectrumRelErrSist_OOBPileUp[m]->SetBinContent(i, 0);
    }

    // 7. SystError mult dependence of signal loss
    fHistSpectrumRelErrSist_SignalLoss[m] = (TH1F *)fHistSpectrumRelErrSist_OOBPileUp[m]->Clone("hSigLossSyst");
    for (Int_t i = 1; i <= fHistSpectrumRelErrSist_SignalLoss[m]->GetNbinsX(); i++)
    {
      fHistSpectrumRelErrSist_SignalLoss[m]->SetBinContent(i, 0.01);
    }

    //*********** SUM IN QUADRATURE OF ALL SOURCES OF UNCERTAINTY *************//
    fHistSpectrumRelErrSist[m] = (TH1F *)fHistSpectrumRelErrSist_Topo[m]->Clone("histoSpectrumRelErrSist" + Smolt[m]);
    fHistSpectrumRelErrSistPtCorr[m] = (TH1F *)fHistSpectrumRelErrSist_Topo[m]->Clone("histoSpectrumRelErrSistPtCorr" + Smolt[m]);
    for (Int_t i = 1; i <= fHistSpectrumRelErrSist[m]->GetNbinsX(); i++)
    {
      fHistSpectrumRelErrSist[m]->SetBinContent(i, TMath::Sqrt(
                                                       TMath::Power(fHistSpectrumRelErrSist_Topo[m]->GetBinContent(i), 2) +
                                                       TMath::Power(fHistSpectrumRelErrSist_MB[m]->GetBinContent(i), 2) +
                                                       TMath::Power(fHistSpectrumRelErrSist_SigExt[m]->GetBinContent(i), 2) +
                                                       TMath::Power(fHistSpectrumRelErrSist_PileUp[m]->GetBinContent(i), 2) +
                                                       TMath::Power(fHistSpectrumRelErrSist_OOBPileUp[m]->GetBinContent(i), 2) +
                                                       TMath::Power(fHistSpectrumRelErrSist_SignalLoss[m]->GetBinContent(i), 2) +
                                                       TMath::Power(fHistSpectrumRelErrSist_MCEff[m]->GetBinContent(i), 2)));
      fHistSpectrumRelErrSistPtCorr[m]->SetBinContent(i, TMath::Sqrt(
                                                             TMath::Power(fHistSpectrumRelErrSist_MB[m]->GetBinContent(i), 2) +
                                                             TMath::Power(fHistSpectrumRelErrSist_SigExt[m]->GetBinContent(i), 2) +
                                                             TMath::Power(fHistSpectrumRelErrSist_PileUp[m]->GetBinContent(i), 2) +
                                                             TMath::Power(fHistSpectrumRelErrSist_OOBPileUp[m]->GetBinContent(i), 2) +
                                                             TMath::Power(fHistSpectrumRelErrSist_SignalLoss[m]->GetBinContent(i), 2) +
                                                             TMath::Power(fHistSpectrumRelErrSist_MCEff[m]->GetBinContent(i), 2)));
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
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasTotalSistRelErr->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, LimInfSpectra, LimSupSpectra, 1, 1, TitleXPt, "Rel. syst. uncertainty", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, 8);
  hDummy->Draw("same");

  for (Int_t m = numMult; m >= 0; m--)
  {
    if (m == numMult)
    {
      ColorMult[m] = ColorMB;
      MarkerMult[m] = MarkerMB;
      SizeMult[m] = SizeMB;
    }
    fHistSpectrumRelErrSist[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumRelErrSist[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumRelErrSist[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumRelErrSist[m]->SetMarkerSize(SizeMult[m]);
    fHistSpectrumRelErrSist[m]->Draw("same");
    fHistSpectrumRelErrSistPtCorr[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumRelErrSistPtCorr[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumRelErrSistPtCorr[m]->SetMarkerStyle(22);
    fHistSpectrumRelErrSistPtCorr[m]->SetMarkerSize(SizeMult[m]);
    fHistSpectrumRelErrSistPtCorr[m]->Draw("same");

    legendAllMult->AddEntry(fHistSpectrumRelErrSist[m], SmoltBis[m] + "%", "pef");
  } // end loop on mult
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  TFile *fileout = new TFile(stringout, "RECREATE");

  for (Int_t m = numMult; m >= 0; m--)
  {
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
