#include "Riostream.h"
#include "string"
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
// fits them to get pt-integrated yields

void CompareYieldFitF(
    Int_t part = 8,
    TString SysPath = "",
    TString OutputDir = "CompareFitFunctions/",
    TString year = Extryear,
    Bool_t isBkgParab = ExtrisBkgParab,
    Bool_t isSysStudy = 1,
    Int_t MultType = ExtrMultType, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
    Bool_t UseTwoGauss = ExtrUseTwoGauss,
    Int_t evFlag = ExtrevFlag // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
)
{

  if (part == 3 || part == 4 || part == 5)
    SysPath = ExtrSysPathXi;
  else if (part == 6 || part == 7 || part == 8)
    SysPath = ExtrSysPathOmega;

  // multiplicity related variables
  TString Smolt[numMult + 1];
  TString SmoltBis[numMult + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn = "";
  TFile *fileIn;

  gStyle->SetOptStat(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 900, 700);
  StyleCanvas(canvasYield, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasChi2 = new TCanvas("canvasChi2", "canvasChi2", 900, 700);
  StyleCanvas(canvasChi2, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasTemp = new TCanvas("canvasTemp", "canvasTemp", 900, 700);
  StyleCanvas(canvasTemp, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasFracExtrYield = new TCanvas("canvasFracExtrYield", "canvasFracExtrYield", 900, 700);
  StyleCanvas(canvasFracExtrYield, 0.05, 0.15, 0.2, 0.02);

  TLegend *legendAllMult = new TLegend(0.22, 0.03, 0.73, 0.28);
  legendAllMult->SetHeader(SMultType[MultType] + " Multiplicity Percentile");
  legendAllMult->SetNColumns(3);
  legendAllMult->SetFillStyle(0);
  TLegendEntry *lheaderAllMult = (TLegendEntry *)legendAllMult->GetListOfPrimitives()->First();
  lheaderAllMult->SetTextSize(0.04);

  TLegend *LegendTitle = new TLegend(0.54, 0.75, 0.95, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  LegendTitle->AddEntry("", NamePart[part] + ", |y| < 0.5", "");

  TLegend *LegendPub = new TLegend(0.58, 0.65, 0.95, 0.73);
  LegendPub->SetFillStyle(0);
  LegendPub->SetTextAlign(32);
  LegendPub->SetTextSize(0.025);

  TLine *lineat1Mult = new TLine(0, 1, numMult + 1, 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  TLine *lineat0Mult = new TLine(0, 0, numMult + 1, 0);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  TLegend *legendfit = new TLegend(0.25, 0.25, 0.4, 0.45);
  legendfit->SetFillStyle(0);
  legendfit->SetTextSize(0.04);

  TLegend *legendfitSummary = new TLegend(0.23, 0.75, 0.42, 0.93);
  legendfitSummary->SetFillStyle(0);
  legendfitSummary->SetTextSize(0.04);
  legendfitSummary->SetTextAlign(13);

  // fit spectra
  Float_t LimInfYield = 0;
  Float_t LimSupYield = 0.01;
  Float_t YoffsetYield = 2;

  TH1F *hYield[numfittipo];
  TH1F *hChi2[numfittipo];
  TH1F *hFracExtrYield[numfittipo];
  TH1F *hTemp[numfittipo];
  TH1F *hYieldPubStat;
  TH1F *hYieldPubSist;

  for (Int_t f = 0; f < numfittipo; f++)
  {
    if (f > 4)
      continue;
    PathIn = "PtIntegratedYields/";
    PathIn += "PtIntegratedYields_" + year;
    PathIn += Spart[part];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[isBkgParab];
    if (isSysStudy)
      PathIn += SysPath;
    PathIn += "_" + nameFitFile[f];
    PathIn += "_" + EventType[evFlag];
    PathIn += ".root";
    cout << "PathIn: " << PathIn << endl;
    fileIn = TFile::Open(PathIn, "READ");

    hYield[f] = (TH1F *)fileIn->Get("hYield");
    hYield[f]->SetName(Form("hYield_%d", f));
    hChi2[f] = (TH1F *)fileIn->Get("hChi2");
    hChi2[f]->SetName(Form("hChi2_%d", f));
    hTemp[f] = (TH1F *)fileIn->Get("hTemp");
    hTemp[f]->SetName(Form("hTemp_%d", f));
    hFracExtrYield[f] = (TH1F *)fileIn->Get("hFracExtrYield");
    hFracExtrYield[f]->SetName(Form("hFracExtrYield_%d", f));
    if (f == 0)
    {
      hYieldPubStat = (TH1F *)fileIn->Get("hYieldPubStat");
      hYieldPubSist = (TH1F *)fileIn->Get("hYieldPubSist");
    }

    StyleHistoYield(hYield[f], LimInfYield, LimSupYield, ColorFit[f], 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
    StyleHistoYield(hYieldPubStat, LimInfYield, LimSupYield, kAzure + 7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
    StyleHistoYield(hYieldPubSist, LimInfYield, LimSupYield, kAzure + 7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);

    if (f == 0)
      LegendPub->AddEntry(hYieldPubSist, "Eur.Phys.J.C 80 (2020) 167, 2020", "pl");
    legendfitSummary->AddEntry(hYield[f], nameFit[f] + " fit", "pl");

    canvasYield->cd();
    hYield[f]->Draw("same e");
    hYieldPubStat->Draw("same e0x0");
    hYieldPubSist->SetFillStyle(0);
    hYieldPubSist->Draw("same e2");

    StyleHistoYield(hChi2[f], 0, 20, ColorFit[f], 22, SMultType[MultType] + " Multiplicity Percentile", "Chi2/NDF", "", 2, 1.15, YoffsetYield);
    canvasChi2->cd();
    hChi2[f]->Draw("same e");

    StyleHistoYield(hTemp[f], 0, 1, ColorFit[f], 22, SMultType[MultType] + " Multiplicity Percentile", "T parameter", "", 2, 1.15, YoffsetYield);
    canvasTemp->cd();
    hTemp[f]->Draw("same e");

    StyleHistoYield(hFracExtrYield[f], 0, 1, ColorFit[f], 22, SMultType[MultType] + " Multiplicity Percentile", "FracExtrYield", "", 2, 1.15, YoffsetYield);
    canvasFracExtrYield->cd();
    hFracExtrYield[f]->Draw("same e");
  }

  canvasYield->cd();
  LegendTitle->Draw("");
  LegendPub->Draw("");
  legendfitSummary->Draw("same");

  canvasChi2->cd();
  LegendTitle->Draw("");
  legendfitSummary->Draw("same");

  canvasTemp->cd();
  LegendTitle->Draw("");
  legendfitSummary->Draw("same");

  canvasFracExtrYield->cd();
  LegendTitle->Draw("");
  legendfitSummary->Draw("same");

  // compute systematic uncertainty associated to choice of fit function
  TCanvas *canvasRatioToTsallis = new TCanvas("canvasRatioToTsallis", "canvasRatioToTsallis", 800, 600);
  StyleCanvas(canvasRatioToTsallis, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasNSigmaToTsallis = new TCanvas("canvasNSigmaToTsallis", "canvasNSigmaToTsallis", 800, 600);
  StyleCanvas(canvasNSigmaToTsallis, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasRelUncFitChoice = new TCanvas("canvasRelUncFitChoice", "canvasRelUncFitChoice", 800, 600);
  StyleCanvas(canvasRelUncFitChoice, 0.05, 0.15, 0.2, 0.02);

  TH1F *hRatioToTsallis[numfittipo];
  TH1F *hNSigmaToTsallis[numfittipo];
  TH1F *hRatioToTsallisMax;
  TH1F *hRelUncFitChoice;
  for (Int_t f = 0; f < numfittipo; f++)
  {
    if (f > 3)
      continue;
    hRatioToTsallis[f] = (TH1F *)hYield[f]->Clone(Form("hRatioToTsallis_%d", f));
    hNSigmaToTsallis[f] = (TH1F *)hYield[f]->Clone(Form("hNSigmaToTsallis_%d", f));
    hRatioToTsallis[f]->Divide(hYield[3]);
    ErrRatioCorr(hYield[f], hYield[3], hRatioToTsallis[f], 1); // I assume full correlation between the two yields
    StyleHistoYield(hRatioToTsallis[f], 0.9, 1.1, ColorFit[f], 22, SMultType[MultType] + " Multiplicity Percentile", "Yield / Yield_{Tsallis}", "", 2, 1.15, YoffsetYield);
    if (f == 0)
      hRatioToTsallisMax = (TH1F *)hRatioToTsallis[f]->Clone("hRatioToTsallisMax");
    for (Int_t i = 0; i < hRatioToTsallis[f]->GetNbinsX(); i++)
    {
      if (TMath::Abs(hRatioToTsallis[f]->GetBinContent(i + 1) - 1) > TMath::Abs(hRatioToTsallisMax->GetBinContent(i + 1) - 1))
        hRatioToTsallisMax->SetBinContent(i + 1, hRatioToTsallis[f]->GetBinContent(i + 1));
    }

    canvasRatioToTsallis->cd();
    if (f != 3)
      hRatioToTsallis[f]->Draw("same e");
    for (Int_t i = 0; i < hRatioToTsallis[f]->GetNbinsX(); i++)
    {
      hNSigmaToTsallis[f]->SetBinContent(i + 1, hRatioToTsallis[f]->GetBinContent(i + 1) / hRatioToTsallis[f]->GetBinError(i + 1));
      hNSigmaToTsallis[f]->SetBinError(i + 1, 0);
    }
    StyleHistoYield(hNSigmaToTsallis[f], -10, 10, ColorFit[f], 22, SMultType[MultType] + " Multiplicity Percentile", "N_{#sigma}", "", 2, 1.15, YoffsetYield);
    canvasNSigmaToTsallis->cd();
    if (f != 3)
      hNSigmaToTsallis[f]->Draw("same e");
  }

  hRelUncFitChoice = (TH1F *)hRatioToTsallisMax->Clone("hRelUncFitChoice");
  for (Int_t b = 0; b < hRelUncFitChoice->GetNbinsX(); b++)
  {
    hRelUncFitChoice->SetBinContent(b + 1, 1. / 2 * TMath::Abs(1 - hRatioToTsallisMax->GetBinContent(b + 1)));
    hRelUncFitChoice->SetBinError(b + 1, 0);
  }

  canvasRatioToTsallis->cd();
  StyleHistoYield(hRatioToTsallisMax, 0.9, 1.1, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "Yield / Yield_{Tsallis}", "", 2, 1.15, YoffsetYield);
  LegendTitle->Draw("");
  legendfitSummary->Draw("same");
  lineat1Mult->Draw("same");
  // hRatioToTsallisMax->Draw("same");

  canvasNSigmaToTsallis->cd();
  LegendTitle->Draw("");
  legendfitSummary->Draw("same");
  lineat0Mult->Draw("same");

  canvasRelUncFitChoice->cd();
  StyleHistoYield(hRelUncFitChoice, 0, 0.1, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "Relative Uncertainty", "", 2, 1.15, YoffsetYield);
  hRelUncFitChoice->Draw("same");
  LegendTitle->Draw("");

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "CompareFitF_" + year;
  stringout += Spart[part];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[isBkgParab];
  if (isSysStudy)
    stringout += SysPath;
  stringout += "_" + EventType[evFlag];
  stringoutpdf = stringout;
  stringout += ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");
  hRelUncFitChoice->Write();
  fileout->Close();

  canvasYield->SaveAs(stringoutpdf + "_YieldsvsPerc.pdf");
  canvasChi2->SaveAs(stringoutpdf + "_Chi2vsPerc.pdf");
  canvasTemp->SaveAs(stringoutpdf + "_TempvsPerc.pdf");
  canvasFracExtrYield->SaveAs(stringoutpdf + "_FracExtrYieldvsPerc.pdf");
  canvasYield->SaveAs(stringoutpdf + "_YieldsvsPerc.png");
  canvasChi2->SaveAs(stringoutpdf + "_Chi2vsPerc.png");
  canvasTemp->SaveAs(stringoutpdf + "_TempvsPerc.png");
  canvasFracExtrYield->SaveAs(stringoutpdf + "_FracExtrYieldvsPerc.png");
  canvasRatioToTsallis->SaveAs(stringoutpdf + "_RatioToTsallis.pdf");
  canvasNSigmaToTsallis->SaveAs(stringoutpdf + "_NSigmaToTsallis.pdf");
  canvasRelUncFitChoice->SaveAs(stringoutpdf + "_RelUncFitChoice.pdf");
  canvasRatioToTsallis->SaveAs(stringoutpdf + "_RatioToTsallis.png");
  canvasNSigmaToTsallis->SaveAs(stringoutpdf + "_NSigmaToTsallis.png");
  canvasRelUncFitChoice->SaveAs(stringoutpdf + "_RelUncFitChoice.png");

  cout << "\nStarting from the files (for the different fit types): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
