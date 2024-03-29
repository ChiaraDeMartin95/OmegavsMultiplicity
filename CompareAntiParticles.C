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

Float_t YLowMean[numPart] = {0.485, 1.110, 1.110, 1.316, 1.316, 1.316, 1.664, 1.664, 1.664};
Float_t YUpMean[numPart] = {0.51, 1.130, 1.130, 1.327, 1.327, 1.327, 1.68, 1.68, 1.68};
Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002};
Float_t YUpSigma[numPart] = {0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};
Float_t YLowPurity[numPart] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0};

Float_t YLowRatio[numChoice] = {0.98, 0.8, 0.8, 0.6, 0.8, 0.8, 0.8};
Float_t YUpRatio[numChoice] = {1.02, 1.2, 1.2, 1.4, 1.2, 1.2, 1.2};

void CompareAntiParticles(Int_t part = 6,
                          Int_t Choice = 0,
                          string SysPath = "",
                          TString OutputDir = "Yields",
                          TString year = Extryear,
                          Bool_t isBkgParab = ExtrisBkgParab,
                          Bool_t isSysStudy = 1,
                          Int_t MultType = ExtrMultType, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0M
                          Bool_t isMB = 1,
                          Int_t mul = 0,
                          Bool_t UseTwoGauss = ExtrUseTwoGauss,
                          Int_t evFlag = ExtrevFlag // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
)
{

  if (mul > numMult)
  {
    cout << "Multiplciity out of range" << endl;
    return;
  }
  if (part < 3 || part == 4 || part == 5 || part == 7 || part == 8)
  {
    cout << "this value of part is not specified, choose 3 (Xi) or 6 (Omega) " << endl;
    return;
  }

  if (part == 3)
    SysPath = ExtrSysPathXi;
  else if (part == 6)
    SysPath = ExtrSysPathOmega;

  TString SPathIn;
  TString SPathInFinal[numParticles];

  TH1F *histo[numParticles];
  TH1F *histoRatio[numParticles];
  TH1F *histoParticle;

  TCanvas *canvas;
  TPad *pad1;
  TPad *pad2;
  canvas = new TCanvas("canvas" + Spart[part], "canvas" + Spart[part], 1000, 800);
  StyleCanvas(canvas, 0.15, 0.05, 0.05, 0.15);
  pad1 = new TPad("pad1" + Spart[part], "pad1" + Spart[part], 0, 0.36, 1, 1);
  pad2 = new TPad("pad2" + Spart[part], "pad2" + Spart[part], 0, 0.01, 1, 0.35);
  StylePad(pad1, 0.15, 0.05, 0.06, 0.01);
  StylePad(pad2, 0.15, 0.05, 0.03, 0.2);

  TF1 *lineMass = new TF1("pol0", "pol0", 0, 8);
  lineMass->SetParameter(0, massParticle[part]);
  lineMass->SetLineColor(kBlack);
  lineMass->SetLineStyle(7);

  TF1 *lineAt1 = new TF1("pol0", "pol0", 0, 8);
  lineAt1->SetParameter(0, 1);
  lineAt1->SetLineColor(kBlack);
  lineAt1->SetLineStyle(1);
  lineAt1->SetLineWidth(1);

  TF1 *lineAt09 = new TF1("pol0", "pol0", 0, 8);
  lineAt09->SetParameter(0, 0.9);
  lineAt09->SetLineColor(kBlack);
  lineAt09->SetLineStyle(7);
  lineAt09->SetLineWidth(1);

  TString Sfileout = "";
  // cout << "Do you want to compare Mean (=0), Sigma (=1), Purity (=2) or Yield per event (=3) or significance (=4)?" << endl;
  // cin >> Choice;
  cout << Choice << " " << TypeHisto[Choice] << endl;
  if (Choice > (numChoice - 1))
  {
    cout << "Option not implemented" << endl;
    return;
  }
  Sfileout = OutputDir + "/CompareAntiParticle" + year + "_" + TypeHisto[Choice] + "_" + Spart[part];
  if (isSysStudy)
    Sfileout += SysPath;
  Sfileout += "_" + EventType[evFlag];

  TLegend *legend;
  if (Choice == 0 || Choice == 1)
    legend = new TLegend(0.2, 0.65, 0.5, 0.9);
  else if (Choice == 2)
    legend = new TLegend(0.2, 0.15, 0.5, 0.4);
  else if (Choice == 3 || Choice == 6)
    legend = new TLegend(0.6, 0.65, 0.9, 0.9);
  else
    legend = new TLegend(0.6, 0.65, 0.9, 0.9);
  legend->SetTextSize(0.035);
  // legend->AddEntry("", NamePart[part], "");

  TFile *filein[numParticles];
  Int_t numParticlesEff = 0;
  Int_t numParticlesEffBis = 0;
  for (Int_t ifile = 0; ifile < numParticles; ifile++)
  {
    if (ifile > 1)
      continue;
    numParticlesEff++;
    if (ifile == 0)
      SPathIn = "Yields/Yields_" + Spart[part];
    else if (ifile == 1)
      SPathIn = "Yields/Yields_" + Spart[part + 1];
    SPathIn += "_" + year;
    SPathIn += IsOneOrTwoGauss[UseTwoGauss];
    SPathIn += SIsBkgParab[isBkgParab];
    if (isMB)
      SPathIn += "_Mult0-100";
    else
      SPathIn += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
    if (Choice == 5)
    {
      if (part == 3)
        SPathIn = "Efficiency/" + ExtrSPathInEffXi;
      else if (part == 6)
        SPathIn = "Efficiency/" + ExtrSPathInEffOmega;
    }
    else if (Choice == 6)
    {
      if (ifile == 0)
        SPathIn = "Yields/YieldEffCorr" + year + "_" + Spart[part];
      else if (ifile == 1)
        SPathIn = "Yields/YieldEffCorr" + year + "_" + Spart[part + 1];
      SPathIn += IsOneOrTwoGauss[UseTwoGauss];
      SPathIn += SIsBkgParab[isBkgParab];
      if (isMB)
        SPathIn += "_Mult0-100";
      else
        SPathIn += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
    }
    if (isSysStudy && Choice != 5)
      SPathIn += SysPath;
    if (Choice != 5)
      SPathIn += "_" + EventType[evFlag];
    SPathInFinal[ifile] = SPathIn;
    if (Choice != 5)
      SPathInFinal[ifile] += ".root";

    cout << "Getting file..." << SPathInFinal[ifile] << endl;
    filein[ifile] = new TFile(SPathInFinal[ifile], "");
    // cout << filein[ifile] << endl;
    if (!filein[ifile])
    {
      cout << "No input file n." << ifile << endl;
      break;
    }
  }

  cout << "\n\e[35mParticle:\e[39m " << Spart[part] << endl;
  cout << "numParticlesEff " << numParticlesEff << endl;
  for (Int_t ifile = 0; ifile < numParticlesEff; ifile++)
  {
    numParticlesEffBis++;
    if (Choice == 0)
    {
      YLow[part] = YLowMean[part];
      YUp[part] = YUpMean[part];
    }
    else if (Choice == 1)
    {
      YLow[part] = YLowSigma[part];
      YUp[part] = YUpSigma[part];
    }
    else if (Choice == 2)
    {
      YLow[part] = YLowPurity[part];
      YUp[part] = 1;
    }
    else if (Choice == 4)
    {
      YUp[part] = 100;
    }
    else if (Choice == 5)
    {
      YLow[part] = 0;
      YUp[part] = 0.2;
    }

    TString inputName;
    TString dirName = "";
    if (SysPath == "_Sel23June")
      dirName = "effCascade";
    else if (SysPath == "_Train100720" || SysPath.find("_Train109") != string::npos || SysPath.find("_Train110") != string::npos)
      dirName = "effAcc";
    else
      dirName = "effOmega";
    if (Choice == 5)
    {
      TDirectoryFile *dir = (TDirectoryFile *)filein[ifile]->Get(dirName);
      if (SysPath == "_Sel23June" || SysPath == "_Train100720" || SysPath.find("_Train109") != string::npos || SysPath.find("_Train110") != string::npos)
        inputName += "hEffCasc";
      else
        inputName = "hEffOmega";
      if (ifile == 0)
        inputName += "Minus";
      else if (ifile == 1)
        inputName += "Plus";
      histo[ifile] = (TH1F *)dir->Get(inputName);
    }
    else
    {
      inputName = "histo" + TypeHisto[Choice];
      histo[ifile] = (TH1F *)filein[ifile]->Get(inputName);
    }
    histo[ifile]->Sumw2();
    if (!histo[ifile])
    {
      cout << "No histo name: " << inputName << " in file " << ifile << endl;
      return;
    }
    if (ifile == 0)
      histoParticle = (TH1F *)histo[ifile]->Clone("histoParticle");

    // Ratios
    histoRatio[ifile] = (TH1F *)histo[ifile]->Clone(Form("hRatio_%i", ifile));
    histoRatio[ifile]->Divide(histoParticle);
    // ErrRatioCorr(histo[ifile], histoParticle, histoRatio[ifile], 1);

    if (Choice == 3 || Choice == 6)
    {
      YLow[part] = 0;
      if (ifile == 0)
        YUp[part] = 1.4 * histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin());
      else
      {
        if (histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin()) > histo[ifile - 1]->GetBinContent(histo[ifile - 1]->GetMaximumBin()))
        {
          YUp[part] = 1.4 * histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin());
        }
      }
    }
    if (Choice != 3 && Choice != 6)
    {
      YLow[part] += 10e-5;
      if (Choice != 2)
        YUp[part] -= 10e-5;
    }

    StyleHisto(histo[ifile], YLow[part], YUp[part], ColorPart[ifile], 33, "", TitleY[Choice], "", 0, 0, 0, 1.5, 1.5, 2);
    StyleHisto(histoRatio[ifile], YLowRatio[Choice], YUpRatio[Choice], ColorPart[ifile], 33, TitleXPt, "Antiparticle/particle", "", 0, 0, 0, 1.5, 1.5, 2);
    histoRatio[ifile]->GetXaxis()->SetLabelSize(0.08);
    histoRatio[ifile]->GetXaxis()->SetTitleSize(0.08);
    histoRatio[ifile]->GetXaxis()->SetTitleOffset(1.2);
    histoRatio[ifile]->GetYaxis()->SetLabelSize(0.08);
    histoRatio[ifile]->GetYaxis()->SetTitleSize(0.08);
    histoRatio[ifile]->GetYaxis()->SetTitleOffset(0.8);
    histoRatio[ifile]->GetXaxis()->SetRangeUser(MinBinPt[part], MaxBinPt[part]);

    canvas->cd();
    pad1->Draw();
    pad1->cd();
    histo[ifile]->GetXaxis()->SetRangeUser(MinBinPt[part], MaxBinPt[part]);
    histo[ifile]->Draw("same");
    if (Choice == 0)
      lineMass->DrawClone("same");
    if (ifile == 0)
      legend->AddEntry(histo[ifile], NamePart[part], "pl");
    else
      legend->AddEntry(histo[ifile], NamePart[part + 1], "pl");
    if (ifile == numParticlesEffBis - 1)
      legend->Draw("");

    canvas->cd();
    pad2->Draw();
    pad2->cd();
    if (ifile != 0)
    {
      histoRatio[ifile]->Draw("same");
      lineAt1->Draw("same");
    }

    cout << "\nA partire dal file: " << SPathInFinal[ifile] << endl;
  }
  canvas->SaveAs(Sfileout + ".pdf");
  canvas->SaveAs(Sfileout + ".png");
  cout << "\nHo creato il file: " << Sfileout << ".pdf and .png" << endl;
}
