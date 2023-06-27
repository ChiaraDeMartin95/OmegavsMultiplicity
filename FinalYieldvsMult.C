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
void StyleTGraphErrors(TGraphAsymmErrors *tgraph, Int_t color, Int_t style, Float_t mSize, Int_t linestyle)
{
  tgraph->SetLineColor(color);
  tgraph->SetLineWidth(3);
  tgraph->SetMarkerColor(color);
  tgraph->SetMarkerStyle(style);
  tgraph->SetMarkerSize(mSize);
  tgraph->SetLineStyle(linestyle);
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

void FinalYieldvsMult(
    Int_t part = 8,
    TString SysPath = "_Sel23June" /*"_Sel6June"*/,
    Bool_t isBkgParab = 1,
    TString OutputDir = "PtIntegratedYields/",
    TString year = "LHC22o_pass4_Train89684" /*"LHC22m_pass4_Train79153"*/,
    Bool_t isSysStudy = 1,
    Int_t MultType = 1, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
    Bool_t UseTwoGauss = 1)
{

  // multiplicity related variables
  TString Smolt[numMult + 1];
  TString SmoltBis[numMult + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn = "";
  TString PathInSistFitChoice = "";
  TFile *fileIn;
  TFile *fileInSistFitChoice;

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "FinalYieldvsMult_" + year;
  stringout += Spart[part];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[isBkgParab];
  if (isSysStudy)
    stringout += SysPath;
  stringoutpdf = stringout;
  stringout += ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");

  gStyle->SetOptStat(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 900, 700);
  StyleCanvas(canvasYield, 0.05, 0.15, 0.2, 0.02);

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

  TLegend *legendfitSummary = new TLegend(0.23, 0.75, 0.42, 0.93);
  legendfitSummary->SetFillStyle(0);
  legendfitSummary->SetTextSize(0.04);
  legendfitSummary->SetTextAlign(13);

  // fit spectra
  Float_t LimInfYield = 0;
  Float_t LimSupYield = 0.01;
  Float_t YoffsetYield = 2;

  TH1F *hYield;
  TH1F *hYieldSist;
  TH1F *hYieldSistTotal;
  TH1F *hYieldSistFitChoice;
  TH1F *hYieldSistRelFitChoice;
  TH1F *hChi2;
  TH1F *hFracExtrYield;
  TH1F *hTemp;
  TH1F *hYieldPubStat;
  TH1F *hYieldPubSist;

  PathIn = "PtIntegratedYields/";
  PathIn += "PtIntegratedYields_" + year;
  PathIn += Spart[part];
  PathIn += IsOneOrTwoGauss[UseTwoGauss];
  PathIn += SIsBkgParab[isBkgParab];
  if (isSysStudy)
    PathIn += SysPath;
  PathIn += "_Levi";
  PathIn += ".root";
  cout << "PathIn: " << PathIn << endl;
  fileIn = TFile::Open(PathIn, "READ");

  PathInSistFitChoice = "CompareFitFunctions/CompareFitF_" + year;
  PathInSistFitChoice += Spart[part];
  PathInSistFitChoice += IsOneOrTwoGauss[UseTwoGauss];
  PathInSistFitChoice += SIsBkgParab[isBkgParab];
  if (isSysStudy)
    PathInSistFitChoice += SysPath;
  PathInSistFitChoice += ".root";
  cout << "PathInSistFitChoice: " << PathInSistFitChoice << endl;
  fileInSistFitChoice = TFile::Open(PathInSistFitChoice, "READ");

  hYield = (TH1F *)fileIn->Get("hYield");
  hYieldSist = (TH1F *)fileIn->Get("hYieldSist");
  hYieldSistRelFitChoice = (TH1F *)fileInSistFitChoice->Get("hRelUncFitChoice");
  hYieldSistFitChoice = (TH1F *)hYieldSist->Clone("hYieldSistFitChoice");
  for (Int_t i = 1; i <= hYieldSistFitChoice->GetNbinsX(); i++)
  {
    hYieldSistFitChoice->SetBinError(i, hYieldSistRelFitChoice->GetBinContent(i) * hYieldSist->GetBinContent(i));
  }
  hYieldSistTotal = (TH1F *)hYieldSist->Clone("hYieldSistTotal");
  for (Int_t i = 1; i <= hYieldSistTotal->GetNbinsX(); i++)
  {
    hYieldSistTotal->SetBinError(i, TMath::Sqrt(hYieldSist->GetBinError(i) * hYieldSist->GetBinError(i) + hYieldSistFitChoice->GetBinError(i) * hYieldSistFitChoice->GetBinError(i)));
  }
  hYieldPubStat = (TH1F *)fileIn->Get("hYieldPubStat");
  hYieldPubSist = (TH1F *)fileIn->Get("hYieldPubSist");

  StyleHistoYield(hYield, LimInfYield, LimSupYield, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldSistTotal, LimInfYield, LimSupYield, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);

  StyleHistoYield(hYieldPubStat, LimInfYield, LimSupYield, kAzure + 7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldPubSist, LimInfYield, LimSupYield, kAzure + 7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);

  LegendPub->AddEntry(hYieldPubSist, "Eur.Phys.J.C 80 (2020) 167, 2020", "pl");

  canvasYield->cd();
  hYieldSistTotal->SetFillStyle(0);
  hYieldSistTotal->Draw("same e2");
  hYield->Draw("same e");
  hYieldPubStat->Draw("same e0x0");
  hYieldPubSist->SetFillStyle(0);
  hYieldPubSist->Draw("same e2");
  LegendTitle->Draw("");
  LegendPub->Draw("");
  legendfitSummary->Draw("same");

  canvasYield->SaveAs(stringoutpdf + "_YieldsvsPerc.pdf");
  canvasYield->SaveAs(stringoutpdf + "_YieldsvsPerc.png");

  //*********** Draw yield vs dNdeta ************//
  Float_t YieldsErrorsStat[numMult + 1] = {0};
  Float_t YieldsErrorsSist[numMult + 1] = {0};
  Float_t Yields[numMult + 1] = {0};
  for (Int_t i = 0; i < numMult; i++)
  {
    YieldsErrorsStat[i] = hYield->GetBinError(i + 1);
    YieldsErrorsSist[i] = hYieldSistTotal->GetBinError(i + 1);
    Yields[i] = hYield->GetBinContent(i + 1);
    cout << "dNdeta " << dNdEtaRun2[i] << endl;
    cout << "Yields[" << i << "] = " << Yields[i] << endl;
  }

  Float_t xTitle = 40;
  Float_t xOffset = 2;

  Float_t yTitle = 30;
  Float_t yOffset = 0.5; // 1.6

  Float_t xLabel = 27;
  Float_t yLabel = 27;
  Float_t xLabelOffset = 0.03;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.03;
  Float_t tickXL = 0.06;
  Float_t tickYL = 0.045;

  Int_t ColorDiff = 1;
  Int_t MarkerType = 20;
  Float_t MarkerSize = 1.5;
  Int_t LineStyle = 1;
  Float_t UpperValueX = 30;
  Float_t Up = 0.025;
  Float_t Low = 1e-4;

  TCanvas *canvasYieldvsdNdeta = new TCanvas("canvasYieldvsdNdeta", "canvasYieldvsdNdeta", 1500, 1500);
  StyleCanvas(canvasYieldvsdNdeta, 0.02, 0.15, 0.1, 0.02);
  TH1F *histoYieldDummy = new TH1F("histoYieldDummy", "histoYieldDummy", 100, 0, UpperValueX);
  StyleHistoYield(histoYieldDummy, Low, Up, 1, 1, TitleXMult, TitleYYieldPtInt, "", 1, 1.15, yOffset);
  SetFont(histoYieldDummy);
  SetHistoTextSize(histoYieldDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  histoYieldDummy->GetYaxis()->SetDecimals(kTRUE);
  histoYieldDummy->GetYaxis()->SetTitleSize(40);
  histoYieldDummy->GetYaxis()->SetTitleOffset(1.75);
  SetTickLength(histoYieldDummy, tickX, tickY);
  histoYieldDummy->Draw("");

  TGraphAsymmErrors *ghistoYield;
  TGraphAsymmErrors *ghistoYieldSist;
  ghistoYield = new TGraphAsymmErrors(numMult, dNdEtaRun2, Yields, 0, 0, YieldsErrorsStat, YieldsErrorsStat);
  ghistoYieldSist = new TGraphAsymmErrors(numMult, dNdEtaRun2, Yields, dNdEtaRun2ErrorL, dNdEtaRun2ErrorR, YieldsErrorsSist, YieldsErrorsSist);
  StyleTGraphErrors(ghistoYield, ColorDiff, MarkerType, MarkerSize, LineStyle);
  StyleTGraphErrors(ghistoYieldSist, ColorDiff, MarkerType, MarkerSize, LineStyle);
  ghistoYield->SetFillColor(ColorDiff);
  ghistoYieldSist->SetFillStyle(0);
  ghistoYieldSist->SetFillColor(ColorDiff);
  ghistoYieldSist->Draw("same p2");
  //  ghistoYieldSistMultUnCorr->SetFillStyle(3001);
  //  ghistoYieldSistMultUnCorr->SetFillColor(ColorDiff);
  //  ghistoYieldSistMultUnCorr->DrawClone("same p2");
  ghistoYield->Draw("same e");

  // get published yields vs mult
  TFile *yieldPub = new TFile("PublishedYield13TeV/HEPData-ins1748157-v1-Table_8cO.root", "");
  TDirectoryFile *dirYieldPub = (TDirectoryFile *)yieldPub->Get("Table 8cO");
  TH1F *hYieldvsMultPub = (TH1F *)dirYieldPub->Get("Hist1D_y1");
  TH1F *hYieldvsMultPubStatErr = (TH1F *)dirYieldPub->Get("Hist1D_y1_e1");
  TH1F *hYieldvsMultPubSistErr = (TH1F *)dirYieldPub->Get("Hist1D_y1_e2");
  TH1F *hYieldvsMultPubStat = (TH1F *)hYieldvsMultPub->Clone("hYieldPubStat");
  TH1F *hYieldvsMultPubSist = (TH1F *)hYieldvsMultPub->Clone("hYieldPubSist");
  for (Int_t b = 1; b <= hYieldvsMultPubStat->GetNbinsX(); b++)
  {
    hYieldvsMultPubStat->SetBinError(b, hYieldvsMultPubStatErr->GetBinContent(b));
    hYieldvsMultPubSist->SetBinError(b, hYieldvsMultPubSistErr->GetBinContent(b));
  }
  hYieldvsMultPubStat->SetMarkerStyle(20);
  hYieldvsMultPubStat->SetMarkerColor(kAzure + 7);
  hYieldvsMultPubStat->SetLineColor(kAzure + 7);
  hYieldvsMultPubSist->SetMarkerStyle(20);
  hYieldvsMultPubSist->SetMarkerColor(kAzure + 7);
  hYieldvsMultPubSist->SetLineColor(kAzure + 7);
  canvasYieldvsdNdeta->cd();
  hYieldvsMultPubStat->Draw("same e0x0");
  hYieldvsMultPubSist->SetFillStyle(0);
  hYieldvsMultPubSist->Draw("same e2");
  LegendTitle->Draw("");
  LegendPub->Draw("");
  canvasYieldvsdNdeta->SaveAs("YieldvsdNdeta.pdf");

  fileout->Close();

  cout << "\nStarting from the files (for the different fit types): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
