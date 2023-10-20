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
#include "/Users/mbp-cdm-01/Desktop/AssegnoRicerca/Run3Analyses/OmegavsMult/Constants.h"
#include "ErrRatioCorr.C"
#include "/Users/mbp-cdm-01/Desktop/AssegnoRicerca/Run3Analyses/OmegavsMult/InputVar.h"

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

void YieldEffCorr(Int_t part = 8,
                  Bool_t isMB = 1,
                  Int_t mul = 0,
                  TString SPathInEff = "",
                  string SysPath = "",
                  TString OutputDir = "Yields",
                  TString year = Extryear,
                  Int_t BkgType = ExtrBkgType, // 0: pol1, 1:pol2, 2:pol3
                  Bool_t isSysStudy = 1,
                  Int_t MultType = ExtrMultType, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0M
                  Bool_t UseTwoGauss = ExtrUseTwoGauss,
                  Int_t evFlag = ExtrevFlag // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
)
{

  if (part == 3 || part == 4 || part == 5)
    SysPath = ExtrSysPathXi;
  else if (part == 6 || part == 7 || part == 8)
    SysPath = ExtrSysPathOmega;

  if (part == 3 || part == 4 || part == 5)
    SPathInEff = ExtrSPathInEffXi;
  else if (part == 6 || part == 7 || part == 8)
    SPathInEff = ExtrSPathInEffOmega;

  TString SfileSignalLossMB = "";
  TString SfileSignalLoss = "";
  if (part == 3 || part == 4 || part == 5)
  {
    SfileSignalLossMB = ExtrSfileSignalLossMBXi;
    SfileSignalLoss = ExtrSfileSignalLossXi;
  }
  else if (part == 6 || part == 7 || part == 8)
  {
    SfileSignalLossMB = ExtrSfileSignalLossMBOmega;
    SfileSignalLoss = ExtrSfileSignalLossOmega;
  }

  if (mul > numMult)
  {
    cout << "Multiplciity out of range" << endl;
    return;
  }
  if (MultType == 0 && (part == 5 || part == 8))
  { // pos + neg cascades
    cout << "No backward compatibility for this case" << endl;
    return;
  }
  if (part < 3)
  {
    return;
  }

  TString Sfileout = OutputDir + "/YieldEffCorr" + year + "_" + Spart[part];
  Sfileout += IsOneOrTwoGauss[UseTwoGauss];
  Sfileout += SIsBkgParab[BkgType];
  if (isMB)
    Sfileout += "_Mult0-100";
  else
    Sfileout += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  if (isSysStudy)
    Sfileout += SysPath;
  Sfileout += "_" + EventType[evFlag];
  // Sfileout += "_Test";

  TString SPathIn;
  SPathIn = "Yields/Yields_" + Spart[part];
  SPathIn += "_" + year;
  SPathIn += IsOneOrTwoGauss[UseTwoGauss];
  SPathIn += SIsBkgParab[BkgType];
  if (isMB)
    SPathIn += "_Mult0-100";
  else
    SPathIn += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  if (isSysStudy)
    SPathIn += SysPath;
  SPathIn += "_" + EventType[evFlag];
  // SPathIn += "_FewPtBins";
  SPathIn += ".root";

  cout << "SPathIn: " << SPathIn << endl;
  TFile *filein = new TFile(SPathIn, "");
  if (!filein)
  {
    cout << "No yield file" << endl;
    return;
  }

  TString SPathInEffFinal = "Efficiency/" + SPathInEff;
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

  TString dirName = "";
  if (SysPath == "_Sel23June")
    dirName = "effCascade";
  else if (SysPath == "_Train100720" || SysPath.find("_Train109") != string::npos || SysPath.find("_Train110") != string::npos)
    dirName = "effAcc";
  else
    dirName = "effOmega";
  TDirectoryFile *dir = (TDirectoryFile *)fileinEff->Get(dirName);
  TString inputNameEff = "hEff";
  if (SysPath == "_Sel23June" || SysPath == "_Train100720" || SysPath.find("_Train109") != string::npos || SysPath.find("_Train110") != string::npos)
    inputNameEff += "Casc";
  else
  {
    if (part >= 3 && part <= 5)
      inputNameEff += "Xi";
    else if (part >= 6 && part <= 8)
      inputNameEff += "Omega";
  }
  if (part == 3 || part == 6)
    inputNameEff += "Minus";
  else if (part == 4 || part == 7)
    inputNameEff += "Plus";
  else if (part == 5 || part == 8)
    inputNameEff += "Sum";
  histoEff = (TH1F *)dir->Get(inputNameEff);
  histoEff->Sumw2();
  histoEff->Smooth();
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
    cout << histoYield->GetBinContent(i) << " +- " << histoYield->GetBinError(i) << " (rel.error: " << histoYield->GetBinError(i) / histoYield->GetBinContent(i) << ")" << endl;
    cout << histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))) << " +- " << histoEff->GetBinError(histoEff->FindBin(histoYield->GetBinCenter(i))) << endl;
    histoYieldCorr->SetBinContent(i, histoYield->GetBinContent(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
    RelErr = sqrt(pow(histoYield->GetBinError(i) / histoYield->GetBinContent(i), 2) + pow(histoEff->GetBinError(histoEff->FindBin(histoYield->GetBinCenter(i))) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))), 2));
    // histoYieldCorr->SetBinError(i, histoYield->GetBinError(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
    histoYieldCorr->SetBinError(i, RelErr * histoYieldCorr->GetBinContent(i));
    cout << histoYieldCorr->GetBinContent(i) << " +- " << histoYieldCorr->GetBinError(i) << endl;
  }

  /*
  cout << "Check: " << endl;
  cout << histoYield->GetNbinsX() << endl;
  cout << histoEff->GetNbinsX() << endl;
  for (Int_t i = 1; i <= histoEff->GetNbinsX(); i++)
  {
    cout << histoEff->GetBinContent(i) << endl;
  }
  for (Int_t i = 1; i <= histoYield->GetNbinsX(); i++)
  {
    cout << histoYield->GetBinCenter(i) << endl;
    cout << histoEff->GetBinCenter(i) << endl;
    histoYieldCorr->SetBinContent(i, histoYield->GetBinContent(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
  }
  */

  // EVENT NORMALISATION CORRECTION
  TFile *fileinNormMB = new TFile(SfileinNormMB, "");
  TFile *fileinNorm = new TFile(SfileinNorm, "");
  TH1F *histoYieldCorrELCorr = (TH1F *)histoYieldCorr->Clone("histoYieldCorrELCorr");
  if (!fileinNorm)
  {
    cout << "No normalisation file" << endl;
    return;
  }
  if (!fileinNormMB)
  {
    cout << "No normalisation file" << endl;
    return;
  }
  TCanvas *cNormMB = (TCanvas *)fileinNormMB->Get("canvasEventCorr");
  if (!cNormMB)
  {
    cout << "No normalisation canvas" << endl;
    return;
  }
  TH1F *histoNormMB = (TH1F *)cNormMB->GetPrimitive("eventCorr");
  if (!histoNormMB)
  {
    cout << "No normalisation histo MB" << endl;
    return;
  }
  histoNormMB->SetName("histoNormMB");
  TCanvas *cNorm = (TCanvas *)fileinNorm->Get("canvasEventCorr");
  if (!cNorm)
  {
    cout << "No normalisation canvas" << endl;
    return;
  }
  TH1F *histoNorm = (TH1F *)cNorm->GetPrimitive("eventCorr");
  if (!histoNorm)
  {
    cout << "No normalisation histo" << endl;
    return;
  }
  histoNorm->SetName("histoNormMB");

  Float_t RelErrSpectrum = 0;
  Float_t EvLoss = histoNorm->GetBinContent(histoNorm->FindBin(MultiplicityPerc[mul] + 0.0001));
  Float_t EvLossMB = histoNormMB->GetBinContent(1);
  Float_t RelErrEL = histoNorm->GetBinError(histoNorm->FindBin(MultiplicityPerc[mul] + 0.0001)) / histoNorm->GetBinContent(histoNorm->FindBin(MultiplicityPerc[mul] + 0.0001));
  Float_t RelErrELMB = histoNormMB->GetBinError(1) / histoNormMB->GetBinContent(1);
  for (Int_t i = 1; i <= histoYieldCorrELCorr->GetNbinsX(); i++)
  {
    if ((histoYieldCorrELCorr->GetXaxis()->GetBinLowEdge(i)) < MinBinPt[part])
    {
      continue;
    }
    RelErrSpectrum = histoYieldCorrELCorr->GetBinError(i) / histoYieldCorrELCorr->GetBinContent(i);
    if (isMB == 1)
    {
      histoYieldCorrELCorr->SetBinContent(i, histoYieldCorrELCorr->GetBinContent(i) * EvLossMB);
      histoYieldCorrELCorr->SetBinError(i, sqrt(pow(RelErrSpectrum, 2) + pow(RelErrELMB, 2)) * histoYieldCorrELCorr->GetBinContent(i));
      cout << "\nNorm factor: " << EvLossMB << endl;
      cout << "spectrum bef. corr: " << histoYieldCorr->GetBinContent(i) << " +- " << histoYieldCorr->GetBinError(i) << " rel error: " << histoYieldCorr->GetBinError(i) / histoYieldCorr->GetBinContent(i) << endl;
      cout << "spectrum after corr: " << histoYieldCorrELCorr->GetBinContent(i) << " +- " << histoYieldCorrELCorr->GetBinError(i) << " rel error: " << histoYieldCorrELCorr->GetBinError(i) / histoYieldCorrELCorr->GetBinContent(i) << endl;
      cout << "Ratio corr / uncorr: " << histoYieldCorrELCorr->GetBinContent(i) / histoYieldCorr->GetBinContent(i) << endl;
    }
    else
    {
      histoYieldCorrELCorr->SetBinContent(i, histoYieldCorrELCorr->GetBinContent(i) * EvLoss);
      histoYieldCorrELCorr->SetBinError(i, sqrt(pow(RelErrSpectrum, 2) + pow(RelErrEL, 2)) * histoYieldCorrELCorr->GetBinContent(i));
      cout << "\nNorm factor: " << EvLoss << endl;
      cout << "Ratio corr / uncorr: " << histoYieldCorrELCorr->GetBinContent(i) / histoYieldCorr->GetBinContent(i) << endl;
    }
  }

  // SIGNAL LOSS CORRECTION
  TFile *fileinSignalLossMB = new TFile(SfileSignalLossMB, "");
  TFile *fileinSignalLoss = new TFile(SfileSignalLoss, "");
  if (!fileinSignalLoss)
  {
    cout << "No signal loss file" << endl;
    return;
  }
  if (!fileinSignalLossMB)
  {
    cout << "No signal loss file MB" << endl;
    return;
  }
  TDirectory *dirSLMB = (TDirectory *)fileinSignalLossMB->Get("effSignal");
  TH1F *histoSignalLossMB = (TH1F *)dirSLMB->Get("hEffCascSumSignal");
  if (!histoSignalLossMB)
  {
    cout << "No signal loss histo MB" << endl;
    return;
  }
  histoSignalLossMB->SetName("histoSignalLossMB");

  TDirectory *dirSL = (TDirectory *)fileinSignalLoss->Get("SignalLossFinalHists");
  TString histoName = "";
  if (MultiplicityPerc[mul] < 30)
    histoName = "hEffCascSumSignal_0";
  else if (MultiplicityPerc[mul] < 70)
    histoName = "hEffCascSumSignal_1";
  else
    histoName = "hEffCascSumSignal_2";
  TH1F *histoSignalLoss = (TH1F *)dirSL->Get(histoName);
  if (!histoSignalLoss)
  {
    cout << "No signal loss histo" << endl;
    return;
  }

  TH1F *histoSignalLossFinal;
  if (isMB == 1)
    histoSignalLossFinal = (TH1F *)histoSignalLossMB->Clone("histoSignalLoss");
  else
    histoSignalLossFinal = (TH1F *)histoSignalLoss->Clone("histoSignalLoss");

  TF1 *fitSLpol0 = new TF1("fitSLpol0", "pol0", 0, 10);
  histoSignalLossFinal->Fit(fitSLpol0, "R+");
  cout << "Chi2/NDF " << fitSLpol0->GetChisquare() / fitSLpol0->GetNDF() << endl;

  TCanvas *cSignalLoss = new TCanvas("cSignalLoss", "cSignalLoss", 800, 600);
  cSignalLoss->cd();
  histoSignalLossFinal->Draw();

  RelErrSpectrum = 0;
  Float_t SigLoss = 0;
  Float_t RelErrSigLoss = 0;
  Float_t SigLossPol0 = 0;
  Float_t RelErrSigLossPol0 = 0;
  Float_t SigLossFinal = 0;
  Float_t RelErrSigLossFinal = 0;
  for (Int_t i = 1; i <= histoYieldCorrELCorr->GetNbinsX(); i++)
  {
    if (histoYieldCorrELCorr->GetXaxis()->GetBinLowEdge(i) < MinBinPt[part])
      continue;
    RelErrSpectrum = histoYieldCorrELCorr->GetBinError(i) / histoYieldCorrELCorr->GetBinContent(i);
    SigLoss = histoSignalLossFinal->GetBinContent(histoSignalLossFinal->FindBin(histoYieldCorrELCorr->GetXaxis()->GetBinCenter(i)));
    SigLossPol0 = fitSLpol0->GetParameter(0);
    RelErrSigLoss = histoSignalLossFinal->GetBinError(histoSignalLossFinal->FindBin(histoYieldCorrELCorr->GetXaxis()->GetBinCenter(i))) / SigLoss;
    RelErrSigLossPol0 = fitSLpol0->GetParError(0) / fitSLpol0->GetParameter(0);
    if (isMB == 1 || MultiplicityPerc[mul] >= 70)
    {
      SigLossFinal = SigLoss;
      RelErrSigLossFinal = RelErrSigLoss;
    }
    else
    {
      SigLossFinal = SigLossPol0;
      RelErrSigLossFinal = RelErrSigLossPol0;
    }
    cout << "\npt interval: " << histoYieldCorrELCorr->GetXaxis()->GetBinLowEdge(i) << "-" << histoYieldCorrELCorr->GetXaxis()->GetBinUpEdge(i) << endl;
    cout << "Signal loss factor: " << SigLoss << ", rel error: " << RelErrSigLoss << endl;
    cout << "Signal loss from fit: " << fitSLpol0->Eval(histoYieldCorrELCorr->GetXaxis()->GetBinCenter(i)) << ", rel error: " << fitSLpol0->GetParError(0) / fitSLpol0->GetParameter(0) << endl;
    cout << "Chosen value: " << SigLossFinal << ", rel error: " << RelErrSigLossFinal << endl;
    cout << "Spectrum bef. SL corr: " << histoYieldCorrELCorr->GetBinContent(i) << " +- " << histoYieldCorrELCorr->GetBinError(i) << " rel error: " << histoYieldCorrELCorr->GetBinError(i) / histoYieldCorrELCorr->GetBinContent(i) << endl;
    histoYieldCorrELCorr->SetBinContent(i, histoYieldCorrELCorr->GetBinContent(i) / SigLossFinal);
    histoYieldCorrELCorr->SetBinError(i, sqrt(pow(RelErrSpectrum, 2) + pow(RelErrSigLoss, 2)) * histoYieldCorrELCorr->GetBinContent(i));
    cout << "Spectrum after SL corr: " << histoYieldCorrELCorr->GetBinContent(i) << " +- " << histoYieldCorrELCorr->GetBinError(i) << " rel error: " << histoYieldCorrELCorr->GetBinError(i) / histoYieldCorrELCorr->GetBinContent(i) << endl;
  }

  // CORRECTION FOR ANCHORING - EFFICIENCY
  // at the moment not applied
  TFile *fileAnchFactor = new TFile(SfileinAnchoring, "");
  if (!fileAnchFactor)
  {
    cout << "No anchoring file" << endl;
    return;
  }
  TH1F *hAnchFactor = (TH1F *)fileAnchFactor->Get("hFinalAnch");
  TCanvas *canvasAnchFactor = new TCanvas("canvasAnchFactor", "canvasAnchFactor", 800, 600);
  hAnchFactor->Draw();

  Float_t AnchFactor = 0;
  for (Int_t i = 1; i <= histoYieldCorrELCorr->GetNbinsX(); i++)
  {
    if (histoYieldCorrELCorr->GetXaxis()->GetBinLowEdge(i) < MinBinPt[part])
      continue;
    AnchFactor = hAnchFactor->GetBinContent(hAnchFactor->FindBin(histoYieldCorrELCorr->GetBinCenter(i)));
    cout << "Spectrum bef. anch. factor corr: " << histoYieldCorrELCorr->GetBinContent(i) << " +- " << histoYieldCorrELCorr->GetBinError(i) << " rel error: " << histoYieldCorrELCorr->GetBinError(i) / histoYieldCorrELCorr->GetBinContent(i) << endl;
    histoYieldCorrELCorr->SetBinContent(i, histoYieldCorrELCorr->GetBinContent(i) * AnchFactor);
    histoYieldCorrELCorr->SetBinError(i, histoYieldCorrELCorr->GetBinError(i) * AnchFactor);
    cout << "Spectrum after anch. factor corr: " << histoYieldCorrELCorr->GetBinContent(i) << " +- " << histoYieldCorrELCorr->GetBinError(i) << " rel error: " << histoYieldCorrELCorr->GetBinError(i) / histoYieldCorrELCorr->GetBinContent(i) << endl;
  }

  Int_t partC = 0;
  if (part == 3 || part == 6)
    partC = 0;
  else if (part == 4 || part == 7)
    partC = 1;
  StyleHisto(histoYieldCorr, 0, 1.3 * histoYieldCorr->GetBinContent(histoYieldCorr->GetMaximumBin()), ColorPart[partC], 33, TitleXPt, titleYield, "", 0, 0, 0, 1.5, 1.5, 2);
  if (part == 6 || part == 7)
    histoYieldCorr->GetYaxis()->SetRangeUser(0, 8 * 1e-4);
  else if (part == 8)
    histoYieldCorr->GetYaxis()->SetRangeUser(0, 16 * 1e-4);
  // histoYieldCorr->GetXaxis()->SetLabelSize(0.08);
  // histoYieldCorr->GetXaxis()->SetTitleSize(0.08);
  // histoYieldCorr->GetXaxis()->SetTitleOffset(1.2);
  // histoYieldCorr->GetYaxis()->SetLabelSize(0.08);
  // histoYieldCorr->GetYaxis()->SetTitleSize(0.08);
  // histoYieldCorr->GetYaxis()->SetTitleOffset(0.8);

  TCanvas *canvas = new TCanvas("canvas" + Spart[part], "canvas" + Spart[part], 1000, 800);
  StyleCanvas(canvas, 0.15, 0.05, 0.05, 0.15);
  canvas->cd();
  histoYieldCorrELCorr->Draw("");
  canvas->SaveAs(Sfileout + ".pdf");
  canvas->SaveAs(Sfileout + ".png");

  // COMPARE TO PUBLISHED YIELDS
  TString SPathInPub = "PublishedYield13TeV/Results-Omega";
  TFile *fileinPub;
  TString histoPubName;
  TH1F *histoPub;
  TH1F *histoPubSyst;
  TH1F *histoPubErr;
  TH1F *histoPubErrSyst;
  TH1F *histoPubFinal;
  TH1F *histoPubFinalSyst;
  TDirectory *dirPub;
  TString SdirPub = "";
  Float_t LimSupComp[numPart] = {0};
   if (part >= 3 && part <= 5)
    LimSupComp[part] = 6;
  else if (part >= 6 && part <= 8)
    LimSupComp[part] = 3;

  if (part >= 3 && part <= 5)
  {
    SPathInPub = "PublishedYield13TeV/HEPData-ins1748157-v1-Table_3.root";
    SdirPub = "Table 3";
    histoPubName = "Hist1D_y11";
  }
  else if (part >= 6 && part <= 8)
  {
    SPathInPub = "PublishedYield13TeV/HEPData-ins1748157-v1-Table_4.root";
    SdirPub = "Table 4";
    histoPubName = "Hist1D_y6";
    // histoPubName = "Hist1D_y4"; // dndeta similar to Run 3 in 30-40
  }
  cout << SPathInPub << endl;
  fileinPub = new TFile(SPathInPub, "");
  if (!fileinPub)
  {
    cout << "No pub yield file " << endl;
    return;
  }
  dirPub = (TDirectory *)fileinPub->Get(SdirPub);
  histoPub = (TH1F *)dirPub->Get(histoPubName);
  histoPub->SetName("histoPub");
  histoPubSyst = (TH1F *)dirPub->Get(histoPubName);
  histoPub->SetName("histoPubSyst");
  histoPubErr = (TH1F *)dirPub->Get(histoPubName + "_e1");
  histoPubErrSyst = (TH1F *)dirPub->Get(histoPubName + "_e2");
  for (Int_t j = 1; j <= histoPub->GetNbinsX(); j++)
  {
    histoPub->SetBinError(j, histoPubErr->GetBinContent(j));
    histoPubSyst->SetBinError(j, histoPubErrSyst->GetBinContent(j));
  }
  histoPubFinal = (TH1F *)histoPub->Clone("histoPubFinal");
  histoPubFinalSyst = (TH1F *)histoPubSyst->Clone("histoPubFinalSyst");
  if (!(part == 5 || part == 8))
  {
    histoPubFinal->Scale(1. / 2);
    histoPubFinalSyst->Scale(1. / 2);
  }

  for (Int_t j = 1; j <= histoPubFinal->GetNbinsX(); j++)
  {
    cout << "PUB error: " << histoPubFinal->GetBinError(j) / histoPubFinal->GetBinContent(j) << endl;
  }

  StyleHisto(histoPubFinal, 0, 1.3 * histoPubFinal->GetBinContent(histoPubFinal->GetMaximumBin()), kAzure + 7, 33, TitleXPt, titleYield, "", 0, 0, 0, 1.5, 1.5, 2);
  StyleHisto(histoPubFinalSyst, 0, 1.3 * histoPubFinal->GetBinContent(histoPubFinal->GetMaximumBin()), kAzure + 7, 33, TitleXPt, titleYield, "", 0, 0, 0, 1.5, 1.5, 2);
  TSpline3 *splinePub;
  TF1 *fsplinePub;
  splinePub = new TSpline3(histoPubFinal, "Spline");
  sp3 = (TSpline3 *)splinePub->Clone("SplineClone");
  fsplinePub = new TF1("fSpline", spline, 0, 5.5);

  TCanvas *canvasComp = new TCanvas("canvasComp" + Spart[part], "canvasComp" + Spart[part], 1000, 800);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  pad1->Draw();
  pad1->cd();
  StyleHisto(histoYieldCorrELCorr, 0, 1.3 * histoYieldCorrELCorr->GetBinContent(histoYieldCorrELCorr->GetMaximumBin()), ColorPart[partC], 33, TitleXPt, titleYield, "", 0, 0, 0, 1.5, 1.5, 2);
  TH1F *histoYieldCorrELCorrClone = (TH1F *)histoYieldCorrELCorr->Clone("histoYieldCorrELCorrClone");
  histoYieldCorrELCorrClone->GetXaxis()->SetRangeUser(MinBinPt[part], LimSupComp[part]);
  histoYieldCorrELCorrClone->Draw("");
  histoPubFinal->Draw("same");
  histoPubFinalSyst->SetFillStyle(0);
  histoPubFinalSyst->Draw("same e2");
  fsplinePub->SetLineColor(kAzure + 7);
  fsplinePub->SetLineWidth(2);
  fsplinePub->Draw("same");

  canvasComp->cd();
  padL1->Draw();
  padL1->cd();
  TH1F *histoRatioToPub = (TH1F *)histoYieldCorrELCorr->Clone("histoRatioToPub");
  TH1F *histoRatioToPubSyst = (TH1F *)histoYieldCorrELCorr->Clone("histoRatioToPubSyst");
  Float_t ALow = 0;
  Float_t AUp = 0;
  Float_t RelYield = 0;
  Float_t RelYieldPub = 0;
  Float_t RelYieldPubSyst = 0;

  for (Int_t b = 1; b <= histoYieldCorrELCorr->GetNbinsX(); b++)
  {
    ALow = histoYieldCorrELCorr->GetBinLowEdge(b);
    AUp = histoYieldCorrELCorr->GetBinLowEdge(b + 1);
    RelYield = histoYieldCorrELCorr->GetBinError(b) / histoYieldCorrELCorr->GetBinContent(b);
    RelYieldPub = histoPubFinal->GetBinError(histoPubFinal->FindBin(AUp - 0.001)) / histoPubFinal->GetBinContent(histoPubFinal->FindBin(AUp - 0.001));
    RelYieldPubSyst = histoPubFinalSyst->GetBinError(histoPubFinalSyst->FindBin(AUp - 0.001)) / histoPubFinalSyst->GetBinContent(histoPubFinalSyst->FindBin(AUp - 0.001));
    cout << "rel yield pub sist " << RelYieldPubSyst << endl;
    histoRatioToPub->SetBinContent(b, histoYieldCorrELCorr->GetBinContent(b) * histoYieldCorrELCorr->GetBinWidth(b) / fsplinePub->Integral(ALow, AUp));
    histoRatioToPub->SetBinError(b, histoRatioToPub->GetBinContent(b) * sqrt(pow(RelYield, 2) + pow(RelYieldPub, 2)));
    histoRatioToPubSyst->SetBinContent(b, histoYieldCorrELCorr->GetBinContent(b) * histoYieldCorrELCorr->GetBinWidth(b) / fsplinePub->Integral(ALow, AUp));
    histoRatioToPubSyst->SetBinError(b, histoRatioToPubSyst->GetBinContent(b) * sqrt(pow(0.06, 2) + pow(RelYieldPubSyst, 2)));
  }
  StyleHisto(histoRatioToPub, 0, 2, ColorPart[partC], 33, TitleXPt, "Ratio to published", "", 0, 0, 0, 1.5, 1.5, 2);
  StyleHisto(histoRatioToPubSyst, 0, 2, ColorPart[partC], 33, TitleXPt, "Ratio to published", "", 0, 0, 0, 1.5, 1.5, 2);
  histoRatioToPub->GetYaxis()->SetTitleSize(0.1);
  histoRatioToPub->GetYaxis()->SetTitleOffset(0.5);
  histoRatioToPub->GetYaxis()->SetLabelSize(0.07);
  histoRatioToPub->GetXaxis()->SetTitleSize(0.1);
  histoRatioToPub->GetXaxis()->SetTitleOffset(0.9);
  histoRatioToPub->GetXaxis()->SetLabelSize(0.07);
  histoRatioToPub->GetYaxis()->SetRangeUser(0.4, 1.6);
  histoRatioToPub->GetXaxis()->SetRangeUser(MinBinPt[part],  LimSupComp[part]);
  histoRatioToPub->Draw("same ex0");
  histoRatioToPubSyst->SetFillStyle(0);
  histoRatioToPubSyst->Draw("same e2");
  TF1 *lineAt1 = new TF1("lineAt1", "1", MinBinPt[part],  LimSupComp[part]);
  lineAt1->SetLineColor(1);
  lineAt1->SetLineWidth(1);
  lineAt1->SetLineStyle(7);
  lineAt1->Draw("same");

  if (isMB == 1)
  {
    canvasComp->SaveAs(Sfileout + "_CompPub.pdf");
    canvasComp->SaveAs(Sfileout + "_CompPub.png");
  }

  TFile *fileout = new TFile(Sfileout + ".root", "RECREATE");
  histoYieldCorr->Write();
  histoYieldCorrELCorr->Write();
  fileout->Close();

  cout << "\nA partire dai file: \n"
       << "\e[35mRaw yield: \e[0m" << SPathIn << "\n"
       << "\e[35mEfficiency: \e[0m" << SPathInEffFinal << "\n"
       << "\e[35mEvent correction: \e[0m" << SfileinNorm << " (for mult classes)\n"
       << "\e[35mEvent correction: \e[0m" << SfileinNormMB << " (for MB)\n"
       << "\eho creato il file: " << Sfileout + ".root\n"
       << endl;
      
}