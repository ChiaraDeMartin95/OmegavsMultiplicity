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
#include "DATA/Constants.h"
#include "DATA/ErrRatioCorr.C"
#include "DATA/InputVar.h"

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

TSpline3 *sp3;
Double_t spline(Double_t *x, Double_t *p)
{
  Double_t xx = x[0];
  return sp3->Eval(xx);
}

TString titlePt = "p_{T} (GeV/c)";
TString titleYield = "1/N_{ev} dN/dp_{T}";

void YieldEffCorr(
    TString sysname = "Default",
    Int_t part = 8,
    Bool_t isMB = 1,
    Int_t mul = 0,
    TString SPathInEff = "Eff-LHC22oapass4-",
    string SysPath = "",
    TString OutputDir = "CorrSpectra",
    Bool_t isBkgParab = ExtrisBkgParab,
    Bool_t isSysStudy = 1,
    Int_t MultType = ExtrMultType, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0M
    Bool_t UseTwoGauss = ExtrUseTwoGauss,
    Int_t evFlag = ExtrevFlag // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
) {

  TString Sfileout = OutputDir + "/YieldEffCorr" + Spart[part];
  //Sfileout += IsOneOrTwoGauss[UseTwoGauss];
  //Sfileout += SIsBkgParab[isBkgParab];
  if (isMB)
    Sfileout += "_Mult0100_";
  else
    Sfileout += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  //if (isSysStudy)
    //Sfileout += SysPath;
  //Sfileout += "_" + EventType[evFlag];
  // Sfileout += "_Test";
  Sfileout += sysname;

  TString SPathIn = "DATA/Yields/Yields_" + Spart[part];
  //SPathIn += "_" + year;
  //SPathIn += IsOneOrTwoGauss[UseTwoGauss];
  //SPathIn += SIsBkgParab[isBkgParab];
  if (isMB)
    SPathIn += "_Mult0100_";
  else
    SPathIn += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  //if (isSysStudy)
    //SPathIn += SysPath;
  //SPathIn += "_" + EventType[evFlag];
  // SPathIn += "_FewPtBins";
  SPathIn += sysname;
  SPathIn += ".root";

  SPathInEff+= sysname;
  SPathInEff += ".root";

  cout << "SPathIn: " << SPathIn << endl;
  TFile *filein = new TFile(SPathIn, "");
  if (!filein)
  {
    cout << "No yield file" << endl;
    return;
  }

  TString SPathInEffFinal = "MC/Eff/" + SPathInEff;
  cout << "SPathInEff: " << SPathInEffFinal << endl;
  TFile *fileinEff = new TFile(SPathInEffFinal, "");
  if (!fileinEff)
  {
    cout << "No efficiency file" << endl;
    return;
  }

  TH1F *histoYield;
  TH1F *histoEff;
  TH1F *histoYieldCorr;

  TString inputName = "histoYield";
  histoYield = (TH1F *)filein->Get(inputName);
  histoYield->Sumw2();
  if (!histoYield)
  {
    cout << "No histo name: " << inputName << endl;
    return;
  }

  TString dirName = "effAcc";
  TDirectoryFile *dir = (TDirectoryFile *)fileinEff->Get(dirName);
  TString inputNameEff = "hEffCascSum";
  histoEff = (TH1F *)dir->Get(inputNameEff);
  histoEff->Sumw2();
  if (!histoEff)
  {
    cout << "No histo name: " << inputNameEff << endl;
    return;
  }

  histoYieldCorr = (TH1F *)histoYield->Clone("histoYieldCorr");

  Float_t RelErr = 0;
  for (Int_t i = 1; i <= histoYield->GetNbinsX(); i++)
  {
    if (histoYield->GetXaxis()->GetBinLowEdge(i) < MinBinPt[part])
    {
      histoYieldCorr->SetBinContent(i, 0);
      histoYieldCorr->SetBinError(i, 0);
      continue;
    }
    cout << "\npt: " << histoYield->GetBinCenter(i) << endl;
    cout << histoYield->GetBinContent(i) << " +- " << histoYield->GetBinError(i) << endl;
    cout << histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))) << " +- " << histoEff->GetBinError(histoEff->FindBin(histoYield->GetBinCenter(i))) << endl;
    histoYieldCorr->SetBinContent(i, histoYield->GetBinContent(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
    RelErr = sqrt(pow(histoYield->GetBinError(i) / histoYield->GetBinContent(i), 2) + pow(histoEff->GetBinError(histoEff->FindBin(histoYield->GetBinCenter(i))) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))), 2));
    // histoYieldCorr->SetBinError(i, histoYield->GetBinError(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
    histoYieldCorr->SetBinError(i, RelErr * histoYieldCorr->GetBinContent(i));
    cout << histoYieldCorr->GetBinContent(i) << " +- " << histoYieldCorr->GetBinError(i) << endl;
  }

  TFile *fileout = new TFile(Sfileout + ".root", "RECREATE");
  histoYieldCorr->Write();
  fileout->Close();

  cout << "\nA partire dai file: \n"
       << "\e[35mRaw yield: \e[0m" << SPathIn << "\n"
       << "\e[35mEfficiency: \e[0m" << SPathInEffFinal << "\n"
       << endl;
}