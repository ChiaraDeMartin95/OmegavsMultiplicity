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

TString titlePt = "p_{T} (GeV/c)";
TString titleYield = "1/N_{ev} dN/dp_{T}";

void YieldEffCorr(Int_t part = 8,
                  Bool_t isMB = 1,
                  Int_t mul = 0,
                  TString SPathInEff = "",
                  string SysPath = "",
                  TString OutputDir = "Yields",
                  TString year = Extryear,
                  Bool_t isBkgParab = ExtrisBkgParab,
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
  Sfileout += SIsBkgParab[isBkgParab];
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
  SPathIn += SIsBkgParab[isBkgParab];
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
  if (!histoEff)
  {
    cout << "No histo name: " << inputNameEff << endl;
    return;
  }

  histoYieldCorr = (TH1F *)histoYield->Clone("histoYieldCorr");

  Float_t RelErr = 0;
  for (Int_t i = 1; i <= histoYield->GetNbinsX(); i++)
  {
    if (part >= 6 && part <= 8 && histoYield->GetXaxis()->GetBinLowEdge(i) < MinBinPt[part])
      continue;
    cout << "\npt: " << histoYield->GetBinCenter(i) << endl;
    cout << histoYield->GetBinContent(i) << " +- " << histoYield->GetBinError(i) << endl;
    cout << histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))) << " +- " << histoEff->GetBinError(histoEff->FindBin(histoYield->GetBinCenter(i))) << endl;
    cout << histoYieldCorr->GetBinContent(i) << " +- " << histoYieldCorr->GetBinError(i) << endl;
    histoYieldCorr->SetBinContent(i, histoYield->GetBinContent(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
    RelErr = sqrt(pow(histoYield->GetBinError(i) / histoYield->GetBinContent(i), 2) + pow(histoEff->GetBinError(histoEff->FindBin(histoYield->GetBinCenter(i))) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))), 2));
    // histoYieldCorr->SetBinError(i, histoYield->GetBinError(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
    histoYieldCorr->SetBinError(i, RelErr * histoYieldCorr->GetBinContent(i));
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

  Float_t Err = 0;
  Float_t ErrEL = histoNorm->GetBinError(histoNorm->FindBin(MultiplicityPerc[mul] + 0.0001)) / histoNorm->GetBinContent(histoNorm->FindBin(MultiplicityPerc[mul] + 0.0001));
  Float_t ErrELMB = histoNormMB->GetBinError(1) / histoNormMB->GetBinContent(1);
  for (Int_t i = 1; i <= histoYieldCorrELCorr->GetNbinsX(); i++)
  {
    if (histoYieldCorrELCorr->GetXaxis()->GetBinLowEdge(i) < MinBinPt[part])
      continue;
    Err = histoYieldCorrELCorr->GetBinError(i) / histoYieldCorrELCorr->GetBinContent(i);
    if (isMB == 1)
    {
      histoYieldCorrELCorr->SetBinContent(i, histoYieldCorrELCorr->GetBinContent(i) * histoNormMB->GetBinContent(1));
      histoYieldCorrELCorr->SetBinError(i, sqrt(pow(Err, 2) + pow(ErrELMB, 2)) * histoYieldCorrELCorr->GetBinContent(i));
      cout << "Norm factor: " << histoNormMB->GetBinContent(1) << endl;
      cout << histoYieldCorr->GetBinContent(i) << " +- " << histoYieldCorr->GetBinError(i) << endl;
      cout << histoYieldCorrELCorr->GetBinContent(i) << " +- " << histoYieldCorrELCorr->GetBinError(i) << endl;
      cout << histoYieldCorrELCorr->GetBinContent(i) / histoYieldCorr->GetBinContent(i) << endl;
    }
    else
    {
      histoYieldCorrELCorr->SetBinContent(i, histoYieldCorrELCorr->GetBinContent(i) * histoNorm->GetBinContent(histoNorm->FindBin(MultiplicityPerc[mul] + 0.0001)));
      histoYieldCorrELCorr->SetBinError(i, sqrt(pow(Err, 2) + pow(ErrEL, 2)) * histoYieldCorrELCorr->GetBinContent(i));
      cout << "Norm factor: " << histoNorm->GetBinContent(mul + 1) << endl;
      cout << histoYieldCorrELCorr->GetBinContent(i) / histoYieldCorr->GetBinContent(i) << endl;
    }
  }

  // SIGNAL LOSS CORRECTION

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
  histoYieldCorr->Draw("");
  canvas->SaveAs(Sfileout + ".pdf");
  canvas->SaveAs(Sfileout + ".png");

  TString SPathInPub = "PublishedYield13TeV/Results-Omega";
  TFile *fileinPub;
  TString histoPubName;
  TH1F *histoPub;
  TH1F *histoPubErr;
  TH1F *histoPubFinal;
  TDirectory *dirPub;
  TString SdirPub = "";

  /*
  for (Int_t i = 0; i <= 1; i++)
  { // loop over particle and antiparticle and sum them if needed

    if (part == 6 && i == 1)
      continue;
    if (part == 7 && i == 0)
      continue;
    SPathInPub = "PublishedYield13TeV/Results-Omega";
    if (i == 0)
      SPathInPub += "Minus";
    else if (i == 1)
      SPathInPub += "Plus";
    SPathInPub += "-V0M-000to100_WithV0refitAndImprovedDCA.root";
    cout << SPathInPub << endl;
    fileinPub = new TFile(SPathInPub, "");
    if (!fileinPub)
    {
      cout << "No pub yield file " << endl;
      return;
    }
    histoPubName = "fHistPtOmega";
    if (i == 0)
      histoPubName += "Minus";
    else if (i == 1)
      histoPubName += "Plus";
    histoPub = (TH1F *)fileinPub->Get(histoPubName);
    if (part == 6 || part == 7)
      histoPubFinal = (TH1F *)histoPub->Clone("histoPubFinal");
    else
    {
      if (i == 0)
        histoPubFinal = (TH1F *)histoPub->Clone("histoPubFinal");
      else
        histoPubFinal->Add(histoPub);
    }
  }
*/
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
  histoPubErr = (TH1F *)dirPub->Get(histoPubName + "_e1");
  for (Int_t j = 1; j <= histoPub->GetNbinsX(); j++)
  {
    histoPub->SetBinError(j, histoPubErr->GetBinContent(j));
  }
  histoPubFinal = (TH1F *)histoPub->Clone("histoPubFinal");
  if (!(part == 5 || part == 8))
    histoPubFinal->Scale(1. / 2);

  StyleHisto(histoPubFinal, 0, 1.3 * histoPubFinal->GetBinContent(histoPubFinal->GetMaximumBin()), kAzure + 7, 33, TitleXPt, titleYield, "", 0, 0, 0, 1.5, 1.5, 2);

  TCanvas *canvasComp = new TCanvas("canvasComp" + Spart[part], "canvasComp" + Spart[part], 1000, 800);
  StyleCanvas(canvasComp, 0.15, 0.05, 0.05, 0.15);
  canvasComp->cd();
  histoYieldCorr->GetYaxis()->SetRangeUser(0, 1.3 * histoPubFinal->GetBinContent(histoPubFinal->GetMaximumBin()));
  histoYieldCorr->GetXaxis()->SetRangeUser(MinBinPt[part], MaxBinPt[part]);
  histoYieldCorr->Draw("same");
  histoPubFinal->Draw("same");
  if (isMB == 1)
  {
    canvasComp->SaveAs(Sfileout + "_CompPub.pdf");
    canvasComp->SaveAs(Sfileout + "_CompPub.png");
  }

  TFile *fileout = new TFile(Sfileout + ".root", "RECREATE");
  histoYieldCorr->Write();
  histoYieldCorrELCorr->Write();
  fileout->Close();

  cout << "\nA partire dal file: \n"
       << SPathIn << " per lo yield raw\n e dal file: " << SPathInEffFinal << " per l'efficienza, \nho creato il file: " << Sfileout << ".pdf and .png and .root" << endl;
}