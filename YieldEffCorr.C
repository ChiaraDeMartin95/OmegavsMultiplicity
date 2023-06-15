#include <Riostream.h>
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

void YieldEffCorr(Int_t part = 6,
                  TString SPathInEff = "eff6June.root",
                  TString SysPath = "_Sel6June",
                  TString OutputDir = "Yields",
                  TString year = "LHC22m_pass4_Train79153",
                  Bool_t isSysStudy = 1,
                  Int_t MultType = 1, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0M
                  Bool_t isMB = 1,
                  Int_t mul = 0,
                  Bool_t UseTwoGauss = 0)
{

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
  if (isMB)
    Sfileout += "_Mult0-100";
  else
    Sfileout += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  if (isSysStudy)
    Sfileout += SysPath;
  //Sfileout += "_Test";

  TString SPathIn;
  SPathIn = "Yields/Yields_" + Spart[part];
  SPathIn += "_" + year;
  SPathIn += IsOneOrTwoGauss[UseTwoGauss];
  if (isMB)
    SPathIn += "_Mult0-100";
  else
    SPathIn += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  if (isSysStudy)
    SPathIn += SysPath;
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

  TDirectoryFile *dir = (TDirectoryFile *)fileinEff->Get("effOmega");
  TString inputNameEff = "hEff";
  if (part >= 3 && part <= 5)
    inputNameEff += "Xi";
  else if (part >= 6 && part <= 8)
    inputNameEff += "Omega";
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
  // histoYieldCorr->Divide(histoEff);
  for (Int_t i = 1; i <= histoEff->GetNbinsX(); i++)
  {
    cout << "\npt: " << histoYield->GetBinCenter(i) << endl;
    cout << histoYield->GetBinContent(i) << " +- " << histoYield->GetBinError(i) << endl;
    cout << histoEff->GetBinContent(i) << " +- " << histoEff->GetBinError(i) << endl;
    cout << histoYieldCorr->GetBinContent(i) << " +- " << histoYieldCorr->GetBinError(i) << endl;
    histoYieldCorr->SetBinContent(i, histoYield->GetBinContent(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
    histoYieldCorr->SetBinError(i, histoYield->GetBinError(i) / histoEff->GetBinContent(histoEff->FindBin(histoYield->GetBinCenter(i))));
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
  TH1F *histoPubFinal;
  for (Int_t i = 0; i <= 1; i++)
  { // loop over particle and antiparticle and sum them if needed
    {
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
  }

  StyleHisto(histoPub, 0, 1.3 * histoPubFinal->GetBinContent(histoPubFinal->GetMaximumBin()), kAzure + 7, 33, TitleXPt, titleYield, "", 0, 0, 0, 1.5, 1.5, 2);

  TCanvas *canvasComp = new TCanvas("canvasComp" + Spart[part], "canvasComp" + Spart[part], 1000, 800);
  StyleCanvas(canvasComp, 0.15, 0.05, 0.05, 0.15);
  canvasComp->cd();
  histoYieldCorr->Draw("same");
  histoPubFinal->Draw("same");
  canvasComp->SaveAs(Sfileout + "_CompPub.pdf");
  canvasComp->SaveAs(Sfileout + "_CompPub.png");

  TFile *fileout = new TFile(Sfileout + ".root", "RECREATE");
  histoYieldCorr->Write();

  cout << "\nA partire dal file: \n"
       << SPathIn << " per lo yield raw\n e dal file: " << SPathInEffFinal << " per l'efficienza, \nho creato il file: " << Sfileout << ".pdf and .png and .root" << endl;
}