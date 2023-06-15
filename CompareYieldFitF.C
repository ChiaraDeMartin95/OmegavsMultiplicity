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
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include </data/dataalice/cdemart/AliPhysicsChiara/AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C>

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

void CompareFitF(
    Int_t part = 5,
    TString SysPath = "_Sel6June",
    TString OutputDir = "CompareFitFunctions/",
    TString year = "LHC22m_pass4_Train79153",
    Bool_t isSysStudy = 1,
    Int_t MultType = 1, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
    Bool_t UseTwoGauss = 0)
{

  // multiplicity related variables
  TString Smolt[numMult + 1];
  TString SmoltBis[numMult + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn = "";
  TFile *fileIn;

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "CompareFitF_" + year;
  stringout += Spart[part];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  if (isSysStudy)
    stringout += SysPath;
  stringout += "_" + nameFitFile[typefit];
  stringoutpdf = stringout;
  stringout += ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");

  gStyle->SetOptStat(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

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

  TLegend *LegendPub = new TLegend(0.54, 0.65, 0.95, 0.73);
  LegendPub->SetFillStyle(0);
  LegendPub->SetTextAlign(33);
  LegendPub->SetTextSize(0.025);

  TLine *lineat1Mult = new TLine(0, 1, 8, 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  TLegend *legendfit = new TLegend(0.25, 0.25, 0.4, 0.45);
  legendfit->SetFillStyle(0);
  legendfit->SetTextSize(0.04);

  TLegend *legendfitSummary = new TLegend(0.73, 0.65, 0.92, 0.75);
  legendfitSummary->SetFillStyle(0);
  legendfitSummary->SetTextSize(0.04);
  legendfitSummary->SetTextAlign(33);
  legendfitSummary->AddEntry("", nameFit[typefit] + " fit", "");

  // fit spectra
  Float_t LimInfYield = 0;
  Float_t LimSupYield = 0.003;
  Float_t YoffsetYield = 2;

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 900, 700);
  StyleCanvas(canvasYield, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasChi2 = new TCanvas("canvasChi2", "canvasChi2", 900, 700);
  StyleCanvas(canvasChi2, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasTemp = new TCanvas("canvasTemp", "canvasTemp", 900, 700);
  StyleCanvas(canvasTemp, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasFracExtrYield = new TCanvas("canvasFracExtrYield", "canvasFracExtrYield", 900, 700);
  StyleCanvas(canvasFracExtrYield, 0.05, 0.15, 0.2, 0.02);

  for (Int_t f = 0; m <= numfittipo 0; m--)
  {
    PathIn = "PtIntegratedYields/";
    PathIn += "PtIntegratedYields_" + year;
    PathIn += Spart[part];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    if (isSysStudy)
      PathIn += SysPath;
    PathIn += "_" + nameFitFile[typefit];
    PathIn += ".root";
    fileIn = TFile::Open(PathIn, "READ");

    hYield[f] = (TH1F *)fileIn->Get("hYield");
    hYield[f]->SetName(Form("hYield_%d", f));
    hChi2[f][f] = (TH1F *)fileIn->Get("hChi2[f]");
    hChi2[f][f]->SetName(Form("hChi2[f]_%d", f));
    hTemp[f][f] = (TH1F *)fileIn->Get("hTemp[f]");
    hTemp[f][f]->SetName(Form("hTemp[f]_%d", f));
    hFracExtrYield[f][f] = (TH1F *)fileIn->Get("hFracExtrYield[f]");
    hFracExtrYield[f][f]->SetName(Form("hFracExtrYield[f]_%d", f));
    if (f == 0)
    {
      hYieldPubStat = (TH1F *)fileIn->Get("hYieldPubStat");
      hYieldPubSist = (TH1F *)fileIn->Get("hYieldPubSist");
    }

    StyleHistoYield(hYield[f], LimInfYield, LimSupYield, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
    StyleHistoYield(hYieldPubStat, LimInfYield, LimSupYield, kAzure + 7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
    StyleHistoYield(hYieldPubSist, LimInfYield, LimSupYield, kAzure + 7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);

    if (f == 0)
      LegendPub->AddEntry(hYieldPubSist, "Eur.Phys.J.C 80 (2020) 167, 2020", "pl");
    canvasYield->cd();
    hYield[f]->Draw("e");
    hYieldPubStat->Draw("same e0x0");
    hYieldPubSist->SetFillStyle(0);
    hYieldPubSist->Draw("same e2");

    StyleHistoYield(hChi2[f], 0, 20, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "Chi2/NDF", "", 2, 1.15, YoffsetYield);
    canvasChi2->cd();
    hChi2[f]->Draw("e");

    StyleHistoYield(hTemp[f], 0, 1, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "T parameter", "", 2, 1.15, YoffsetYield);
    canvasTemp->cd();
    hTemp[f]->Draw("e");

    StyleHistoYield(hFracExtrYield[f], 0, 1, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "FracExtrYield", "", 2, 1.15, YoffsetYield);
    canvasFracExtrYield->cd();
    hFracExtrYield[f]->Draw("e");
  }

  canvasYield->cd();
  LegendTitle->Draw("");
  LegendPub->Draw("");

  canvasChi2->cd();
  LegendTitle->Draw("");
  legendfitSummary->Draw("");

  canvasChi2->cd();
  LegendTitle->Draw("");
  legendfitSummary->Draw("");

  canvasChi2->cd();
  LegendTitle->Draw("");
  legendfitSummary->Draw("");

  canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
  canvasYield->SaveAs(stringoutpdf + "_YieldsvsPerc.pdf");
  canvasChi2->SaveAs(stringoutpdf + "_Chi2vsPerc.pdf");
  canvasTemp->SaveAs(stringoutpdf + "_TempvsPerc.pdf");
  canvasFracExtrYield->SaveAs(stringoutpdf + "_FracExtrYieldvsPerc.pdf");
  canvasPtSpectra->SaveAs(stringoutpdf + ".png");
  canvasYield->SaveAs(stringoutpdf + "_YieldsvsPerc.png");
  canvasChi2->SaveAs(stringoutpdf + "_Chi2vsPerc.png");
  canvasTemp->SaveAs(stringoutpdf + "_TempvsPerc.png");
  canvasFracExtrYield->SaveAs(stringoutpdf + "_FracExtrYieldvsPerc.png");
  fileout->Close();

  cout << "\nStarting from the files (for the different fit types): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}