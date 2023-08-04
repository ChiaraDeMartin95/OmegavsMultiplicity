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
    Int_t dNdEtaFlag = 1,
    Int_t part = 8,
    TString SysPath = "",
    Bool_t isBkgParab = ExtrisBkgParab,
    TString OutputDir = "PtIntegratedYields/",
    TString year = Extryear,
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
  stringout += "_" + EventType[evFlag];
  stringoutpdf = stringout;
  stringout += ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");

  gStyle->SetOptStat(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 900, 700);
  StyleCanvas(canvasYield, 0.05, 0.15, 0.2, 0.02);

  TCanvas *canvasRelErrorYield = new TCanvas("canvasRelErrorYield", "canvasRelErrorYield", 900, 700);
  StyleCanvas(canvasRelErrorYield, 0.05, 0.15, 0.2, 0.02);

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

  TLegend *legendErrorSummary = new TLegend(0.23, 0.75, 0.42, 0.93);
  legendErrorSummary->SetFillStyle(0);
  legendErrorSummary->SetTextSize(0.04);
  legendErrorSummary->SetTextAlign(13);

  // fit spectra
  Float_t LimInfYield = 0;
  Float_t LimSupYield = 0.01;
  Float_t YoffsetYield = 2;

  TH1F *hYield;
  TH1F *hYieldSist;
  TH1F *hYieldSistTotal;
  TH1F *hYieldSistFitChoice;
  TH1F *hYieldSistTopoSel;
  TH1F *hYieldSistRelFitChoice;
  TH1F *hYieldSistRelTopoSel;
  TH1F *hYieldSistRelTotal;
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
  PathIn += "_" + EventType[evFlag];
  PathIn += ".root";
  cout << "PathIn: " << PathIn << endl;
  fileIn = TFile::Open(PathIn, "READ");

  PathInSistFitChoice = "CompareFitFunctions/CompareFitF_" + year;
  PathInSistFitChoice += Spart[part];
  PathInSistFitChoice += IsOneOrTwoGauss[UseTwoGauss];
  PathInSistFitChoice += SIsBkgParab[isBkgParab];
  if (isSysStudy)
    PathInSistFitChoice += SysPath;
  PathInSistFitChoice += "_" + EventType[evFlag];
  PathInSistFitChoice += ".root";
  cout << "PathInSistFitChoice: " << PathInSistFitChoice << endl;
  fileInSistFitChoice = TFile::Open(PathInSistFitChoice, "READ");

  hYield = (TH1F *)fileIn->Get("hYield");
  hYield->SetName("hYieldStrange");
  hYieldSist = (TH1F *)fileIn->Get("hYieldSist");
  hYieldSistTopoSel = (TH1F *)hYieldSist->Clone("hYieldSistTopoSel");
  hYieldSistRelFitChoice = (TH1F *)fileInSistFitChoice->Get("hRelUncFitChoice");
  hYieldSistFitChoice = (TH1F *)hYieldSist->Clone("hYieldSistFitChoice");
  Float_t RelErrorSystYieldTopoSel = RelSystYieldTopoSel;
  Float_t ErrorSystYieldTopoSel[numMult + 1] = {0};
  for (Int_t i = 1; i <= hYieldSistTopoSel->GetNbinsX(); i++)
  {
    ErrorSystYieldTopoSel[i] = RelErrorSystYieldTopoSel * hYieldSistTopoSel->GetBinContent(i);
    hYieldSistTopoSel->SetBinError(i, ErrorSystYieldTopoSel[i]);
  }
  for (Int_t i = 1; i <= hYieldSistFitChoice->GetNbinsX(); i++)
  {
    hYieldSistFitChoice->SetBinError(i, hYieldSistRelFitChoice->GetBinContent(i) * hYieldSist->GetBinContent(i));
  }
  hYieldSistTotal = (TH1F *)hYieldSist->Clone("hYieldSistTotal");
  for (Int_t i = 1; i <= hYieldSistTotal->GetNbinsX(); i++)
  {
    // hYieldSistTotal->SetBinError(i, TMath::Sqrt(hYieldSist->GetBinError(i) * hYieldSist->GetBinError(i) + hYieldSistFitChoice->GetBinError(i) * hYieldSistFitChoice->GetBinError(i)));
    hYieldSistTotal->SetBinError(i, TMath::Sqrt(
                                        hYieldSist->GetBinError(i) * hYieldSist->GetBinError(i) +
                                        hYieldSistTopoSel->GetBinError(i) * hYieldSistTopoSel->GetBinError(i) +
                                        hYieldSistFitChoice->GetBinError(i) * hYieldSistFitChoice->GetBinError(i)));
  }

  hYieldSistRelTotal = (TH1F *)hYieldSistTotal->Clone("hYieldSistRelTotal");
  hYieldSistRelTopoSel = (TH1F *)hYieldSistTopoSel->Clone("hYieldSistRelTopoSel");
  for (Int_t i = 1; i <= hYieldSistRelTotal->GetNbinsX(); i++)
  {
    hYieldSistRelTotal->SetBinContent(i, hYieldSistTotal->GetBinError(i) / hYieldSistTotal->GetBinContent(i));
    hYieldSistRelTotal->SetBinError(i, 0);
    hYieldSistRelTopoSel->SetBinContent(i, hYieldSistTopoSel->GetBinError(i) / hYieldSistTopoSel->GetBinContent(i));
    hYieldSistRelTopoSel->SetBinError(i, 0);
  }

  hYieldPubStat = (TH1F *)fileIn->Get("hYieldPubStat");
  hYieldPubSist = (TH1F *)fileIn->Get("hYieldPubSist");

  StyleHistoYield(hYield, LimInfYield, LimSupYield, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldSistTotal, LimInfYield, LimSupYield, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);

  StyleHistoYield(hYieldPubStat, LimInfYield, LimSupYield, kAzure + 7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldPubSist, LimInfYield, LimSupYield, kAzure + 7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);

  StyleHistoYield(hYieldSistRelTotal, 0, 0.3, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "Rel. error", "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldSistRelTopoSel, 0, 0.3, kGreen + 2, 22, SMultType[MultType] + " Multiplicity Percentile", "Rel. error", "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldSistRelFitChoice, 0, 0.3, kViolet, 22, SMultType[MultType] + " Multiplicity Percentile", "Rel. error", "", 2, 1.15, YoffsetYield);

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

  canvasYield->SaveAs(stringoutpdf + "_YieldsvsPerc.pdf");
  canvasYield->SaveAs(stringoutpdf + "_YieldsvsPerc.png");

  canvasRelErrorYield->cd();
  hYieldSistRelTotal->Draw("same");
  hYieldSistRelTopoSel->Draw("same");
  hYieldSistRelFitChoice->Draw("same");
  LegendTitle->Draw("");
  legendErrorSummary->AddEntry(hYieldSistRelTotal, "Total", "pl");
  legendErrorSummary->AddEntry(hYieldSistRelTopoSel, "Topological selection", "pl");
  legendErrorSummary->AddEntry(hYieldSistRelFitChoice, "Fit choice", "pl");
  legendErrorSummary->Draw("same");
  canvasRelErrorYield->SaveAs(stringoutpdf + "_RelErrorYieldsvsPerc.pdf");
  canvasRelErrorYield->SaveAs(stringoutpdf + "_RelErrorYieldsvsPerc.png");

  //*********** Draw yield vs dNdeta ************//
  Float_t dNdEta[numMult + 1] = {0};
  Float_t dNdEtaErrorL[numMult + 1] = {0};
  Float_t dNdEtaErrorR[numMult + 1] = {0};
  Float_t dNdEtaMB = 0;
  Float_t dNdEtaMBErrorL = 0;
  Float_t dNdEtaMBErrorR = 0;
  Float_t dNdEtaToMB[numMult + 1] = {0};
  Float_t dNdEtaErrorLToMB[numMult + 1] = {0};
  Float_t dNdEtaErrorRToMB[numMult + 1] = {0};

  Float_t YieldsErrorsStat[numMult + 1] = {0};
  Float_t YieldsErrorsSist[numMult + 1] = {0};
  Float_t Yields[numMult + 1] = {0};
  Float_t YieldMB = 0;
  Float_t YieldMBErrorSist = 0;
  Float_t YieldMBErrorStat = 0;
  Float_t YieldsToMB[numMult + 1] = {0};
  Float_t YieldsErrorsSistToMB[numMult + 1] = {0};
  Float_t YieldsErrorsStatToMB[numMult + 1] = {0};

  for (Int_t i = 0; i < numMult; i++)
  {
    YieldsErrorsStat[i] = hYield->GetBinError(i + 1);
    YieldsErrorsSist[i] = hYieldSistTotal->GetBinError(i + 1);
    Yields[i] = hYield->GetBinContent(i + 1);
    if (dNdEtaFlag == 0)
    { // Run2
      dNdEta[i] = dNdEtaRun2[i];
      dNdEtaErrorL[i] = dNdEtaRun2ErrorL[i];
      dNdEtaErrorR[i] = dNdEtaRun2ErrorR[i];
    }
    else if (dNdEtaFlag == 1)
    { // Run3 Nicolò estimate
      dNdEta[i] = dNdEtaRun3[i];
      dNdEtaErrorL[i] = dNdEtaRun3ErrorL[i];
      dNdEtaErrorR[i] = dNdEtaRun3ErrorR[i];
    }
    cout << "dNdeta " << dNdEta[i] << endl;
    cout << "Yields[" << i << "] = " << Yields[i] << endl;
  }

  YieldMBErrorStat = hYield->GetBinError(numMult + 1);
  YieldMBErrorSist = hYieldSistTotal->GetBinError(numMult + 1);
  YieldMB = hYield->GetBinContent(numMult + 1);
  cout << "YieldMB = " << YieldMB << endl;
  cout << "YieldMBErrorStat = " << YieldMBErrorStat << endl;

  if (dNdEtaFlag == 0)
  { // Run2
    dNdEtaMB = dNdEtaRun2MB;
    dNdEtaMBErrorL = dNdEtaRun2MBErrorL;
    dNdEtaMBErrorR = dNdEtaRun2MBErrorR;
  }
  else if (dNdEtaFlag == 1)
  { // Run3 Nicolò estimate
    dNdEtaMB = dNdEtaRun3MB;
    dNdEtaMBErrorL = dNdEtaRun3MBErrorL;
    dNdEtaMBErrorR = dNdEtaRun3MBErrorR;
  }
  cout << "\n\e[35mdNdeta MB " << dNdEtaMB << "\n\e[39m" << endl;
  for (Int_t i = 0; i < numMult; i++)
  {
    dNdEtaToMB[i] = dNdEta[i] / dNdEtaMB;
    dNdEtaErrorLToMB[i] = dNdEtaToMB[i] * sqrt(pow(dNdEtaErrorL[i] / dNdEta[i], 2) + pow(dNdEtaMBErrorL / dNdEtaMB, 2));
    dNdEtaErrorRToMB[i] = dNdEtaToMB[i] * sqrt(pow(dNdEtaErrorR[i] / dNdEta[i], 2) + pow(dNdEtaMBErrorR / dNdEtaMB, 2));
    cout << "dNdeta  " << dNdEta[i] << endl;
    cout << "dNdeta / dNdetaMB " << dNdEtaToMB[i] << endl;
    YieldsToMB[i] = Yields[i] / YieldMB;
    YieldsErrorsStatToMB[i] = YieldsToMB[i] * sqrt(pow(YieldsErrorsStat[i] / Yields[i], 2) + pow(YieldMBErrorStat / YieldMB, 2));
    YieldsErrorsSistToMB[i] = YieldsToMB[i] * sqrt(pow(YieldsErrorsSist[i] / Yields[i], 2) + pow(YieldMBErrorSist / YieldMB, 2));
    cout << "Yields / YieldMB " << YieldsToMB[i] << endl;
  }

  Float_t YieldPubMB = 0;
  Float_t YieldPubMBErrStat = 0;
  Float_t YieldPubMBErrSist = 0;

  if (part == 3)
  { // Xi
    YieldPubMB = YieldXiMB13TeV;
    YieldPubMBErrStat = YieldXiMB13TeVErrStat;
    YieldPubMBErrSist = YieldXiMB13TeVErrSist;
  }
  if (part == 8)
  { // Omega
    YieldPubMB = YieldOmegaMB13TeV;
    YieldPubMBErrStat = YieldOmegaMB13TeVErrStat;
    YieldPubMBErrSist = YieldOmegaMB13TeVErrSist;
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
  Float_t UpperValueX = 35;
  Float_t Up = 0.025;
  Float_t Low = 1e-4;

  // fit points and get ratio
  TF1 *fitPubYieldvsMult = new TF1("fitPubYieldvsMult", "pol1", 0, 35);
  fitPubYieldvsMult->SetLineColor(kAzure + 7);
  fitPubYieldvsMult->SetLineStyle(2);
  TF1 *fitRun3YieldvsMult = new TF1("fitRun3YieldvsMult", "pol1", 0, 35);
  fitRun3YieldvsMult->SetLineColor(1);
  fitRun3YieldvsMult->SetLineStyle(2);

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
  histoYieldDummy->DrawClone("");

  TGraphAsymmErrors *ghistoYield;
  TGraphAsymmErrors *ghistoYieldSist;
  ghistoYield = new TGraphAsymmErrors(numMult, dNdEta, Yields, 0, 0, YieldsErrorsStat, YieldsErrorsStat);
  ghistoYieldSist = new TGraphAsymmErrors(numMult, dNdEta, Yields, dNdEtaErrorL, dNdEtaErrorR, YieldsErrorsSist, YieldsErrorsSist);
  StyleTGraphErrors(ghistoYield, ColorDiff, MarkerType, MarkerSize, LineStyle);
  StyleTGraphErrors(ghistoYieldSist, ColorDiff, MarkerType, MarkerSize, LineStyle);
  ghistoYield->SetFillColor(ColorDiff);
  ghistoYieldSist->SetFillStyle(0);
  ghistoYieldSist->SetFillColor(ColorDiff);
  ghistoYield->Fit(fitRun3YieldvsMult, "R0");
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
  hYieldvsMultPubStat->Fit("fitPubYieldvsMult", "R0");
  hYieldvsMultPubStat->Draw("same e0x0");
  hYieldvsMultPubSist->SetFillStyle(0);
  hYieldvsMultPubSist->Draw("same e2");
  fitRun3YieldvsMult->Draw("same");
  fitPubYieldvsMult->Draw("same");
  LegendTitle->Draw("");
  LegendPub->Draw("");
  canvasYieldvsdNdeta->SaveAs(stringoutpdf + "_OmegaYieldvsdNdeta.pdf");
  canvasYieldvsdNdeta->SaveAs(stringoutpdf + "_OmegaYieldvsdNdeta.png");

  // Yield/Yield MB vs dNdeta/dNdeta MB
  TCanvas *canvasYieldvsdNdetaToMB = new TCanvas("canvasYieldvsdNdetaToMB", "canvasYieldvsdNdetaToMB", 1500, 1500);
  StyleCanvas(canvasYieldvsdNdetaToMB, 0.02, 0.15, 0.1, 0.02);
  histoYieldDummy->GetXaxis()->SetTitle(TitleXMultToMB);
  histoYieldDummy->GetYaxis()->SetTitle(TitleYYieldPtIntToMB);
  histoYieldDummy->GetXaxis()->SetRangeUser(0, 3.5);
  histoYieldDummy->GetYaxis()->SetRangeUser(0, 5);
  histoYieldDummy->DrawClone("");

  Float_t YieldsPubToMB[10] = {0};
  Float_t YieldsPubErrorsStatToMB[10] = {0};
  Float_t YieldsPubErrorsSistToMB[10] = {0};
  Float_t dNdEtaPubToMB[10] = {0};
  Float_t dNdEtaErrorPubToMB[10] = {0};
  for (Int_t i = 0; i < hYieldvsMultPubStat->GetNbinsX(); i++)
  {
    if (hYieldvsMultPubStat->GetBinContent(i + 1) == 0)
      continue;
    dNdEtaPubToMB[i] = hYieldvsMultPubStat->GetBinCenter(i + 1) / dNdEtaRun2MB;
    dNdEtaErrorPubToMB[i] = hYieldvsMultPubStat->GetBinWidth(i + 1) / dNdEtaRun2MB / 2;
    YieldsPubToMB[i] = hYieldvsMultPubStat->GetBinContent(i + 1) / YieldPubMB;
    YieldsPubErrorsStatToMB[i] = YieldsPubToMB[i] * sqrt(pow(hYieldvsMultPubStat->GetBinError(i + 1) / hYieldvsMultPubStat->GetBinContent(i + 1), 2) + pow(YieldPubMBErrStat / YieldPubMB, 2));
    YieldsPubErrorsSistToMB[i] = YieldsPubToMB[i] * sqrt(pow(hYieldvsMultPubSist->GetBinError(i + 1) / hYieldvsMultPubStat->GetBinContent(i + 1), 2) + pow(YieldPubMBErrSist / YieldPubMB, 2));
    cout << hYieldvsMultPubStat->GetBinCenter(i + 1) << " dNdetaMB Run 2 " << dNdEtaRun2MB << endl;
  }

  TGraphAsymmErrors *ghistoYieldPubToMB;
  TGraphAsymmErrors *ghistoYieldPubSistToMB;
  ghistoYieldPubToMB = new TGraphAsymmErrors(numMult, dNdEtaPubToMB, YieldsPubToMB, 0, 0, YieldsPubErrorsStatToMB, YieldsPubErrorsStatToMB);
  ghistoYieldPubSistToMB = new TGraphAsymmErrors(numMult, dNdEtaPubToMB, YieldsPubToMB, dNdEtaErrorPubToMB, dNdEtaErrorPubToMB, YieldsPubErrorsSistToMB, YieldsPubErrorsSistToMB);
  StyleTGraphErrors(ghistoYieldPubToMB, kAzure + 7, MarkerType, MarkerSize, LineStyle);
  StyleTGraphErrors(ghistoYieldPubSistToMB, kAzure + 7, MarkerType, MarkerSize, LineStyle);
  ghistoYieldPubToMB->SetFillColor(kAzure + 7);
  ghistoYieldPubSistToMB->SetFillStyle(0);
  ghistoYieldPubSistToMB->SetFillColor(kAzure + 7);
  ghistoYieldPubSistToMB->Draw("same p2");
  ghistoYieldPubToMB->Draw("same e");

  TGraphAsymmErrors *ghistoYieldToMB;
  TGraphAsymmErrors *ghistoYieldSistToMB;
  ghistoYieldToMB = new TGraphAsymmErrors(numMult, dNdEtaToMB, YieldsToMB, 0, 0, YieldsErrorsStatToMB, YieldsErrorsStatToMB);
  ghistoYieldSistToMB = new TGraphAsymmErrors(numMult, dNdEtaToMB, YieldsToMB, dNdEtaErrorLToMB, dNdEtaErrorRToMB, YieldsErrorsSistToMB, YieldsErrorsSistToMB);
  StyleTGraphErrors(ghistoYieldToMB, ColorDiff, MarkerType, MarkerSize, LineStyle);
  StyleTGraphErrors(ghistoYieldSistToMB, ColorDiff, MarkerType, MarkerSize, LineStyle);
  ghistoYieldToMB->SetFillColor(ColorDiff);
  ghistoYieldSistToMB->SetFillStyle(0);
  ghistoYieldSistToMB->SetFillColor(ColorDiff);
  ghistoYieldSistToMB->Draw("same p2");
  ghistoYieldToMB->Draw("same e");

  canvasYieldvsdNdetaToMB->SaveAs(stringoutpdf + "_OmegaYieldvsdNdetaToMB.pdf");
  canvasYieldvsdNdetaToMB->SaveAs(stringoutpdf + "_OmegaYieldvsdNdetaToMB.png");

  // Ratio between fit functions
  TCanvas *canvasRatio = new TCanvas("canvasRatio", "canvasRatio", 1500, 1500);
  TH1F *hRatio = (TH1F *)histoYieldDummy->Clone("hRatio");
  for (Int_t b = 1; b <= hRatio->GetNbinsX(); b++)
  {
    hRatio->SetBinContent(b, fitRun3YieldvsMult->Eval(hRatio->GetBinCenter(b)) / fitPubYieldvsMult->Eval(hRatio->GetBinCenter(b)));
    hRatio->SetBinError(b, 0);
  }
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerColor(1);
  hRatio->SetLineColor(1);
  hRatio->SetMarkerSize(1.5);
  hRatio->GetXaxis()->SetRangeUser(5, UpperValueX);
  hRatio->GetYaxis()->SetRangeUser(0., 1.);
  hRatio->GetYaxis()->SetTitle("Ratio");
  hRatio->Draw();

  // Ratio to pions
  // Get Run3 pion yield vs mult
  TString SDirName = "Run3Pions/Run3PiKPr";
  TFile *filePionsPos;
  TH1F *hTempPionsPos;
  TFile *filePionsNeg;
  TH1F *hTempPionsNeg;
  TString sfilePionPos;
  TString sfilePionNeg;
  Float_t PionYield = 0;
  Float_t PionYieldError = 0;
  TH1F *hYieldPions = (TH1F *)fileIn->Get("hYield");
  hYieldPions->SetName("hYieldPions");
  for (Int_t i = 0; i < numMult + 1; i++)
  {
    PionYield = 0;
    PionYieldError = 0;
    sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", MultiplicityPerc[i], MultiplicityPerc[i + 1]);
    sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", MultiplicityPerc[i], MultiplicityPerc[i + 1]);
    if (MultiplicityPerc[i] == 20 || MultiplicityPerc[i] == 25)
    {
      sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", 20.0, 30.0);
      sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", 20.0, 30.0);
    }
    if (MultiplicityPerc[i] == 30 || MultiplicityPerc[i] == 35)
    {
      sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", 30.0, 40.0);
      sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", 30.0, 40.0);
    }
    if (MultiplicityPerc[i] == 50)
    {
      sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", 50.0, 60.0);
      sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", 50.0, 60.0);
    }
    if (MultiplicityPerc[i] == 70)
    {
      sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", 70.0, 80.0);
      sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", 70.0, 80.0);
    }
    if (i == numMult)
    {
      // sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", 0.0, 100.0);
      // sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", 0.0, 100.0);
      sfilePionPos = "Yields_Pos_Pi_Levi.root";
      sfilePionNeg = "Yields_Neg_Pi_Levi.root";
    }

    filePionsPos = new TFile(SDirName + "/" + sfilePionPos, "");
    filePionsNeg = new TFile(SDirName + "/" + sfilePionNeg, "");
    hTempPionsPos = (TH1F *)filePionsPos->Get("IntegratedWithLevi");
    hTempPionsPos->SetName("hTempPionsPos");
    hTempPionsNeg = (TH1F *)filePionsNeg->Get("IntegratedWithLevi");
    hTempPionsNeg->SetName("hTempPionsNeg");
    PionYield = hTempPionsPos->GetBinContent(1) + hTempPionsNeg->GetBinContent(1);
    PionYieldError = pow(hTempPionsPos->GetBinContent(2), 2) + pow(hTempPionsNeg->GetBinContent(2), 2);
    if (MultiplicityPerc[i] == 50)
    {
      sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", 60.0, 70.0);
      sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", 60.0, 70.0);
      filePionsPos = new TFile(SDirName + "/" + sfilePionPos, "");
      filePionsNeg = new TFile(SDirName + "/" + sfilePionNeg, "");
      hTempPionsPos = (TH1F *)filePionsPos->Get("IntegratedWithLevi");
      hTempPionsPos->SetName("hTempPionsPos");
      hTempPionsNeg = (TH1F *)filePionsNeg->Get("IntegratedWithLevi");
      hTempPionsNeg->SetName("hTempPionsNeg");
      PionYield += hTempPionsPos->GetBinContent(1) + hTempPionsNeg->GetBinContent(1);
      PionYieldError += pow(hTempPionsPos->GetBinContent(2), 2) + pow(hTempPionsNeg->GetBinContent(2), 2);
      PionYield = PionYield / 2;
    }
    if (MultiplicityPerc[i] == 70)
    {
      sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", 80.0, 90.0);
      sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", 80.0, 90.0);
      filePionsPos = new TFile(SDirName + "/" + sfilePionPos, "");
      filePionsNeg = new TFile(SDirName + "/" + sfilePionNeg, "");
      hTempPionsPos = (TH1F *)filePionsPos->Get("IntegratedWithLevi");
      hTempPionsPos->SetName("hTempPionsPos");
      hTempPionsNeg = (TH1F *)filePionsNeg->Get("IntegratedWithLevi");
      hTempPionsNeg->SetName("hTempPionsNeg");
      PionYield += hTempPionsPos->GetBinContent(1) + hTempPionsNeg->GetBinContent(1);
      PionYieldError += pow(hTempPionsPos->GetBinContent(2), 2) + pow(hTempPionsNeg->GetBinContent(2), 2);
      sfilePionPos = Form("Yields_Pos_Pi_%.2f_%.2f_Levi.root", 90.0, 100.0);
      sfilePionNeg = Form("Yields_Neg_Pi_%.2f_%.2f_Levi.root", 90.0, 100.0);
      filePionsPos = new TFile(SDirName + "/" + sfilePionPos, "");
      filePionsNeg = new TFile(SDirName + "/" + sfilePionNeg, "");
      hTempPionsPos = (TH1F *)filePionsPos->Get("IntegratedWithLevi");
      hTempPionsPos->SetName("hTempPionsPos");
      hTempPionsNeg = (TH1F *)filePionsNeg->Get("IntegratedWithLevi");
      hTempPionsNeg->SetName("hTempPionsNeg");
      PionYield += hTempPionsPos->GetBinContent(1) + hTempPionsNeg->GetBinContent(1);
      PionYieldError += pow(hTempPionsPos->GetBinContent(2), 2) + pow(hTempPionsNeg->GetBinContent(2), 2);
      PionYield = PionYield / 3;
    }
    PionYieldError = sqrt(PionYieldError);
    if (MultiplicityPerc[i] == 50)
    {
      PionYieldError = PionYieldError / 2;
    }
    if (MultiplicityPerc[i] == 70)
    {
      PionYieldError = PionYieldError / 3;
    }
    hYieldPions->SetBinContent(i + 1, PionYield);
    hYieldPions->SetBinError(i + 1, PionYieldError);
  }

  TCanvas *canvasYieldPions = new TCanvas("canvasYieldPions", "canvasYieldPions", 900, 700);
  StyleCanvas(canvasYieldPions, 0.05, 0.15, 0.2, 0.02);

  StyleHistoYield(hYieldPions, 0, 10, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  // StyleHistoYield(hYieldSistTotalPions, LimInfYield, LimSupYield, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  hYieldPions->Draw("same");
  canvasYieldPions->SaveAs(stringoutpdf + "_PionYield.pdf");
  canvasYieldPions->SaveAs(stringoutpdf + "_PionYield.png");

  // Omega over pion ratio
  TCanvas *canvasRatioOmToPi = new TCanvas("canvasRatioOmToPi", "canvasRatioOmToPi", 900, 700);
  StyleCanvas(canvasRatioOmToPi, 0.05, 0.15, 0.2, 0.02);
  TH1F *hYieldRatio = (TH1F *)fileIn->Get("hYield");
  hYieldRatio->SetName("hYieldRatio");
  hYieldRatio->Sumw2();
  hYieldRatio->Divide(hYieldPions);
  StyleHistoYield(hYieldRatio, 0, 0.001, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "#Omega / #pion", "", 2, 1.15, 1.2);
  hYieldRatio->Draw("same");

  canvasRatioOmToPi->SaveAs(stringoutpdf + "_RatioOmToPi.pdf");
  canvasRatioOmToPi->SaveAs(stringoutpdf + "_RatioOmToPi.png");

  fileout->Close();

  cout << "\nStarting from the files (for the different fit types): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
