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

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRan\
ge,
                Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize)
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
// produces ratio of two plots (efficiencies)

Float_t YLow = {0};
Float_t YUp = {0};

Float_t YLowRatio = {0.99};
Float_t YUpRatio = {1.01};

void EfficiencyRatio(TString year0 = "GapTriggeredAnchoredTo22o",
                     TString year1 = "GapTriggeredNotAnchored",
                     TString Sfilein0 = "Efficiency/effChiaraOmega_inelgt0_gt13tev_27aug.root",
                     TString Sfilein1 = "Efficiency/effChiaraOmega_inelgt0_gapTriggered_14aug.root",
                     Int_t Choice = 0, // 0:effxacc, 1:signal loss
                     Int_t ChosenPart = 8,
                     TString OutputDir = "Efficiency/")
{
  
  Int_t ChosenType = -1;
  TString TypeHisto = "EfficiencyvsPt";
  TString Spart[numPart] = {"K0S", "Lam", "ALam", "XiMin", "XiPlu", "Xi", "OmMin", "OmPlu", "Om"};
  TString NamePart[numPart] = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
  TString TitleY[numChoice] = {"Mean (GeV/#it{c}^{2})", "Sigma (GeV/#it{c}^{2})", "S/(S+B)", "1/#it{N}_{evt} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}", "Efficiency"};
  TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";

  Float_t YLow[numPart] = {0};
  Float_t YUp[numPart] = {0};
  Float_t YLowRatio = 0;
  Float_t YUpRatio = 2;
  if (Choice==1){
    YLowRatio = 0.9;
    YUpRatio = 1;
  }

  Int_t color0 = kRed + 2;
  Int_t color1 = kBlue + 2;

  TH1F *histo0[numPart];
  TH1F *histo1[numPart];
  TH1F *histoRatio[numPart];
  TCanvas *canvas[numPart];
  TPad *pad1[numPart];
  TPad *pad2[numPart];

  TString Sfileout = "";

  Bool_t isGen = 0;

  TF1 *pol0fit = new TF1("pol0fit", "pol0", 0, 10);
  pol0fit->SetLineColor(kBlack);
  pol0fit->SetLineWidth(1);
  pol0fit->SetLineStyle(2);

  TFile *filein0 = new TFile(Sfilein0, "");
  if (!filein0)
  {
    cout << "No input file n.0" << endl;
    return;
  }
  TFile *filein1 = new TFile(Sfilein1, "");
  if (!filein1)
  {
    cout << "No input file n.1" << endl;
    return;
  }

  // Sfileout = OutputDir + "Compare" + TypeHisto + "_" + Spart[ChosenType] + "_" + year0 + "vs" + year1;
  Sfileout = OutputDir + "Compare" + TypeHisto + "_";
  Sfileout += Spart[ChosenType];
  Sfileout += year0 + "vs" + year1;
  cout << "Output file: " << Sfileout << endl;

  for (Int_t part = 0; part < numPart; part++)
  {
    if (part != ChosenPart)
      continue;
    cout << "\n\e[35mParticle:\e[39m " << Spart[part] << endl;

    // TString inputName = TypeHisto + "_" + Spart[part];
    //  histo0[part] = (TH1F *)filein0->Get(inputName);
    TString inputName = "hEffCascSum";
    if (Choice==1) inputName = "hEffCascSumSignal";
    TString namedir = "effAcc";
    if (Choice==1) namedir = "effSignal";
    TDirectoryFile *dir = (TDirectoryFile *)filein0->Get(namedir);
    if (!dir)
      return;
    histo0[part] = (TH1F *)dir->Get(inputName);
    if (!histo0[part])
    {
      cout << "No histo name: " << inputName << " in file0" << endl;
      return;
    }
    histo0[part]->SetName(inputName + "_file0");
    // histo1[part] = (TH1F *)filein1->Get(inputName);
    dir = (TDirectoryFile *)filein1->Get(namedir);
    histo1[part] = (TH1F *)dir->Get(inputName);
    if (!histo1[part])
    {
      cout << "No histo name: " << inputName << " in file1" << endl;
      return;
    }
    histo1[part]->SetName(inputName + "_file1");
    //histo0[part]->Smooth();
    //histo1[part]->Smooth();

    // Ratios
    histoRatio[part] = (TH1F *)histo0[part]->Clone(inputName + "_Ratio");
    if (histo0[part]->GetNbinsX() != histo1[part]->GetNbinsX())
    {
      cout << "The number of bins of the two histograms are different " << endl;
      return;
    }
    histoRatio[part]->Divide(histo1[part]);

    for (Int_t b = 1; b <= histoRatio[part]->GetNbinsX(); b++)
    {
      // cout << "Num: " << histo0[part]->GetBinContent(b) << endl;
      // cout << "Denom " << histo1[part]->GetBinContent(b) << endl;
      // cout << "Ratio " << histoRatio[part]->GetBinContent(b) << endl;
    }

    histoRatio[part]->Fit("pol0fit", "R+");
    canvas[part] = new TCanvas("canvas" + Spart[part], "canvas" + Spart[part], 1000, 800);
    StyleCanvas(canvas[part], 0.15, 0.05, 0.05, 0.15);
    pad1[part] = new TPad("pad1" + Spart[part], "pad1" + Spart[part], 0, 0.36, 1, 1);
    pad2[part] = new TPad("pad2" + Spart[part], "pad2" + Spart[part], 0, 0.01, 1, 0.35);
    StylePad(pad1[part], 0.15, 0.05, 0.05, 0.01);
    StylePad(pad2[part], 0.15, 0.05, 0.03, 0.2);

    TLegend *legend;
    if (Spart[part] == "XiPlu" && Choice == 0)
      legend = new TLegend(0.5, 0.25, 0.8, 0.45);
    // else if (Spart[part] == "K0S" && Choice == 1)
    // legend = new TLegend(0.5, 0.5, 0.8, 0.7);
    else if (part <= 2 && Choice == 2)
      legend = new TLegend(0.5, 0.25, 0.8, 0.45);
    else
      legend = new TLegend(0.5, 0.75, 0.8, 0.9);
    legend->AddEntry("", NamePart[part], "");

    YUp[part] = 0.2;
    if (Choice==1) {
      YLow[part] = 0.8;
      YUp[part] = 1;
    }
    TString TitleY = "Efficiency x Acc.";
    if (Choice==1) TitleY = "Signal loss";
    StyleHisto(histo0[part], YLow[part], YUp[part], color0, 33, "", TitleY, "", 0, 0, 0, 1.5, 1.5, 2);
    StyleHisto(histo1[part], YLow[part], YUp[part], color1, 33, "", TitleY, "", 0, 0, 0, 1.5, 1.5, 2);
    StyleHisto(histoRatio[part], YLowRatio, YUpRatio, color0, 33, TitleXPt, "Ratio to non-anch.", "", 0, 0, 0, 1.5, 1.5, 2);
    histoRatio[part]->GetXaxis()->SetLabelSize(0.08);
    histoRatio[part]->GetXaxis()->SetTitleSize(0.08);
    histoRatio[part]->GetXaxis()->SetTitleOffset(1.2);
    histoRatio[part]->GetYaxis()->SetLabelSize(0.08);
    histoRatio[part]->GetYaxis()->SetTitleSize(0.08);
    histoRatio[part]->GetYaxis()->SetTitleOffset(0.8);

    canvas[part]->cd();
    pad1[part]->Draw();
    pad1[part]->cd();
    histo0[part]->Draw("same");
    histo1[part]->Draw("same");
    legend->AddEntry(histo0[part], year0, "pl");
    legend->AddEntry(histo1[part], year1, "pl");
    legend->Draw("");

    canvas[part]->cd();
    pad2[part]->Draw();
    pad2[part]->cd();
    histoRatio[part]->Draw("same");

    if (part == 0)
      canvas[part]->SaveAs(Sfileout + ".pdf(");
    else if (part == numPart - 1)
      canvas[part]->SaveAs(Sfileout + ".pdf)");
    else
      canvas[part]->SaveAs(Sfileout + ".pdf");
  }

  cout << "\nI started from the files: " << endl;
  cout << Sfilein0 << "\n"
       << Sfilein1 << endl;

  cout << "\nI created the file: " << endl;
  cout << Sfileout << endl;
}
