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

Double_t SetEfficiencyError(Int_t k, Int_t n)
{
  return sqrt(((Double_t)k + 1) * ((Double_t)k + 2) / (n + 2) / (n + 3) - pow((Double_t)(k + 1), 2) / pow(n + 2, 2));
}

void EstimateEventsMC(Float_t Factor = 100, Int_t NumBins = 1, Int_t part = 5)
{

  Float_t NumEvents = 0;
  gStyle->SetOptStat(0);
  //TString PathIn = "Efficiency/effChiaraOmega_inelgt0_gt13tev_27aug.root";
  TString PathIn = "Efficiency/effChiaraXi_inelgt0_lhc23f4b2_26sep.root";
  if (PathIn == "Efficiency/effChiaraOmega_inelgt0_gt13tev_27aug.root") NumEvents = 1.7; //x 10^6
  else if (PathIn == "Efficiency/effChiaraXi_inelgt0_lhc23f4b2_26sep.root") NumEvents = 20; //x 10^6
  TFile *fileMC = new TFile(PathIn, "");
  if (!fileMC)
  {
    cout << "No fileMC" << endl;
    return;
  }
  TString dirName = "effAcc";
  TDirectoryFile *dir = (TDirectoryFile *)fileMC->Get(dirName);
  TString inputNameEff = "hEffCasc";
  if (part == 3 || part == 6)
    inputNameEff += "Minus";
  else if (part == 4 || part == 7)
    inputNameEff += "Plus";
  else if (part == 5 || part == 8)
    inputNameEff += "Sum";
  TH1F *histoEff = (TH1F *)dir->Get(inputNameEff);

  if (!histoEff)
  {
    cout << "No histo name: " << inputNameEff << endl;
    return;
  }

  Int_t Denom[histoEff->GetNbinsX()];
  Int_t Num[histoEff->GetNbinsX()];
  Double_t ExError[histoEff->GetNbinsX()];
  Int_t TargetDenom[histoEff->GetNbinsX()];
  Int_t TargetNum[histoEff->GetNbinsX()];
  Double_t TargetError[histoEff->GetNbinsX()];
  TH1F *heffStatErr = (TH1F *)histoEff->Clone("heffStatErr");
  TH1F *heffStatErrEx = (TH1F *)histoEff->Clone("heffStatErrEx");
  TH1F *heffStatErrTarget = (TH1F *)histoEff->Clone("heffStatErrTarget");
  TH1F *heffStatErrEstimate = (TH1F *)histoEff->Clone("heffStatErrEstimate");
  for (Int_t iBin = 1; iBin <= heffStatErr->GetNbinsX(); iBin++)
  {
    if (histoEff->GetBinContent(iBin) != 0)
    {
      heffStatErr->SetBinContent(iBin, histoEff->GetBinError(iBin) / histoEff->GetBinContent(iBin));
      heffStatErr->SetBinError(iBin, 0);

      Denom[iBin - 1] = 30000; //number of generated particles in the pT bin with less particles
      Num[iBin - 1] = histoEff->GetBinContent(iBin) * Denom[iBin - 1];
      ExError[iBin - 1] = SetEfficiencyError(Num[iBin - 1], Denom[iBin - 1]);
      heffStatErrEx->SetBinContent(iBin, ExError[iBin - 1] / histoEff->GetBinContent(iBin));
      heffStatErrEx->SetBinError(iBin, 0);

      TargetDenom[iBin - 1] = Denom[iBin - 1] * Factor;
      TargetNum[iBin - 1] = histoEff->GetBinContent(iBin) * TargetDenom[iBin - 1];
      TargetError[iBin - 1] = SetEfficiencyError((float)TargetNum[iBin - 1] / NumBins, (float)TargetDenom[iBin - 1] / NumBins);
      heffStatErrTarget->SetBinContent(iBin, TargetError[iBin - 1] / histoEff->GetBinContent(iBin));
      heffStatErrTarget->SetBinError(iBin, 0);

      heffStatErrEstimate->SetBinContent(iBin, heffStatErr->GetBinContent(iBin) * 1. / sqrt(Factor / NumBins));
      heffStatErrEstimate->SetBinError(iBin, 0);
    }
    else
    {
      heffStatErr->SetBinContent(iBin, 0);
      heffStatErr->SetBinError(iBin, 0);
      heffStatErrEx->SetBinContent(iBin, 0);
      heffStatErrEx->SetBinError(iBin, 0);
      heffStatErrTarget->SetBinContent(iBin, 0);
      heffStatErrTarget->SetBinError(iBin, 0);
      heffStatErrEstimate->SetBinContent(iBin, 0);
      heffStatErrEstimate->SetBinError(iBin, 0);
    }
  }

  TCanvas *cEff = new TCanvas("cEff", "cEff", 800, 600);
  StyleCanvas(cEff, 0.1, 0.15, 0.15, 0.05);
  StyleHisto(histoEff, 0.0, 1.0, kBlack, 20, "p_{T} (GeV/c)", "Efficiency", "");
  StyleHisto(heffStatErr, 0.0, 0.5, kBlack, 20, "p_{T} (GeV/c)", "Efficiency relative error", "");
  StyleHisto(heffStatErrEx, 0.0, 0.5, kRed, 20, "p_{T} (GeV/c)", "Efficiency relative error", "");
  StyleHisto(heffStatErrTarget, 0.0, 0.5, kBlue, 20, "p_{T} (GeV/c)", "Efficiency relative error", "");
  StyleHisto(heffStatErrEstimate, 0.0, 0.5, kGreen, 20, "p_{T} (GeV/c)", "Efficiency relative error", "");
  heffStatErr->Draw("same");
  //heffStatErrEx->Draw("same");
  //heffStatErrTarget->Draw("same");
  heffStatErrEstimate->Draw("same");

  TLegend *legEff = new TLegend(0.2, 0.65, 0.35, 0.85);
  legEff->SetBorderSize(0);
  legEff->SetFillStyle(0);
  legEff->SetTextFont(43);
  legEff->SetTextSize(20);
  if (PathIn == "Efficiency/effChiaraOmega_inelgt0_gt13tev_27aug.root") {
    legEff->AddEntry(heffStatErr, Form("%.1f x 10^{6} events", NumEvents), "p");
    legEff->AddEntry(heffStatErrEstimate, Form("%.1f x 10^{6} events, rel. error in each eta bin (%i bins)", Factor * NumEvents, NumBins), "p");
  }
  else if (PathIn == "Efficiency/effChiaraXi_inelgt0_lhc23f4b2_26sep.root") {
    legEff->AddEntry(heffStatErr, Form("%.1f x 10^{6} events", NumEvents), "p");
    legEff->AddEntry(heffStatErrEstimate, Form("%.1f x 10^{6} events, rel. error in each eta bin (%i bins)", Factor * NumEvents, NumBins), "p");
  }
  //legEff->AddEntry(heffStatErrEx, "3.0 x 10^{4} generated Xis per p_{T} interval", "p");
  //legEff->AddEntry(heffStatErrTarget, Form("x %.1f generated Xis, rel. error in each eta bin (%i bins)", Factor, NumBins), "p");
  
  legEff->Draw("same");
}
