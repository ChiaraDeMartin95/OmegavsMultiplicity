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

TString TitleInvMass[numPart] = {"(#pi^{+}, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(p, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(#bar{p}, #pi^{-}) invariant mass (GeV/#it{c}^{2})", "(#Lambda, #pi^{-}) invariant mass (GeV/#it{c}^{2})"};
TString namehisto[numPart] = {"h3dMassK0Short", "", "", "hCascMinusInvMassvsPt", "hCascPlusInvMassvsPt", "hCascMinusInvMassvsPt", "hCascPlusInvMassvsPt"};
Float_t YLowMean[numPart] = {0.485, 1.110, 1.110, 1.316, 1.316, 1.316, 1.664, 1.664, 1.664};
Float_t YUpMean[numPart] = {0.51, 1.130, 1.130, 1.327, 1.327, 1.327, 1.68, 1.68, 1.68};
Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002};
Float_t YUpSigma[numPart] = {0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015};
Float_t YLowPurity[numPart] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0};

TString CutLabelSummary[25] = {"MassWin", "y", "EtaDau", "dcapostopv", "dcanegtopv", "dcabachtopv",
                               "CascCosPA", "V0CosPA", "DCACascDau", "DCAV0Dau", "rCasc", "rV0", "DCAV0ToPV",
                               "LambdaMass", "TPCPr", "TPCPi", "TOFPr", "TOFPi", "TPCBach",
                               "TOFBach", "proplifetime", "rejcomp", "ptthrtof"};
// cascospa, dca casc dau, dcabachtopv, dcapostopv, dcanegtopv, lambdamass, rejcomp, nsigmatpcKa, cascradius, v0radius, dcav0dau, v0cospa, casclifetime, bachbaryon, dcabachbar, dcav0topv
Int_t Bin[] = {7, 9, 6, 4, 5, 14, 22, 19, 11, 12, 10, 8, 21, 24, 25, 13};

// Float_t YLowRatio[numChoice] = {0.98, 0, 0.9, 0.6, 0.6};
// Float_t YUpRatio[numChoice] = {1.02, 1.2, 1.5, 1.2, 1.2};
Float_t YLowRatio[numChoice] = {0.98, 0.8, 0.9, 0.8, 0.9};
Float_t YUpRatio[numChoice] = {1.02, 1.2, 1.6, 1.1, 1.2};

void CompareYields(Int_t itopovar = 2, // cospa, dcacascdau, dcabachtopv, dcapostopv, dcanegtopv, lambdamass, rejcomp, nsigmatpcKa, cascradius, v0radius, dcav0dau, v0cospa, casclifetime, cosbachbaryon, dcabachbar, dcav0topv
                   Int_t Choice = 0,
                   TString SysPath = "",
                   Bool_t isBkgParab = 1,
                   TString OutputDir = "CompareTopo",
                   TString year = "LHC22o_pass4_MinBias_Train108123" /*"LHC22o_pass4_Train89684" /*"LHC22m_pass4_Train79153"*/,
                   Int_t part = 5,
                   Bool_t isSysStudy = 1,
                   Int_t MultType = 1, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0M
                   Bool_t isMB = 1,
                   Int_t mul = 0,
                   Bool_t UseTwoGauss = 1,
                   Int_t evFlag = ExtrevFlag // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
)
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
    cout << "this value of part is not specified, choose 3 - 4 (Xi) or 5 - 6  (Omega) " << endl;
  }

  TString SPathIn;
  TString SPathInFinal[numFiles];
  SPathIn = "Yields/Yields_" + Spart[part] + "_" + year;
  SPathIn += IsOneOrTwoGauss[UseTwoGauss];
  SPathIn += SIsBkgParab[isBkgParab];
  if (isMB)
    SPathIn += "_Mult0-100";
  else
    SPathIn += Form("_Mult%.1f-%.1f", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
  SPathIn += "_run528531";
  SPathIn += "_" + TopoVar[itopovar];

  TH1F *histo[numFiles];
  TH1F *histoRatio[numFiles];
  TH1F *histoLoosest;
  Float_t CutValue[numFiles];

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
  Sfileout = OutputDir + "/Compare" + TypeHisto[Choice] + "_" + Spart[part] + "_" + TopoVar[itopovar];

  TLegend *legend;
  if (Choice == 0 || Choice == 1)
    legend = new TLegend(0.2, 0.65, 0.5, 0.9);
  else if (Choice == 2)
    legend = new TLegend(0.2, 0.15, 0.5, 0.4);
  else if (Choice == 3)
    legend = new TLegend(0.6, 0.65, 0.9, 0.9);
  else if (Choice == 4)
    legend = new TLegend(0.6, 0.65, 0.9, 0.9);
  legend->SetTextSize(0.035);
  legend->AddEntry("", NamePart[part], "");

  TFile *filein[numFiles];
  // Int_t numFilesEff = -1;
  Int_t numFilesEff = 0;
  Int_t numFilesEffBis = 0;
  for (Int_t ifile = 0; ifile < numFiles; ifile++)
  {
    numFilesEff++;
    SPathInFinal[ifile] = SPathIn + Form("%i", ifile);
    SPathInFinal[ifile] += "_" + EventType[evFlag];
    SPathInFinal[ifile] += ".root";
    filein[ifile] = new TFile(SPathInFinal[ifile], "");
    // cout << "Getting file..." << SPathInFinal[ifile] << endl;
    // cout << filein[ifile] << endl;
    if (!filein[ifile])
    {
      cout << "No input file n." << ifile << endl;
      break;
    }
  }

  cout << "\n\e[35mParticle:\e[39m " << Spart[part] << endl;
  cout << "numFilesEff " << numFilesEff << endl;
  for (Int_t ifile = 0; ifile < numFilesEff; ifile++)
  {
    if (itopovar == 1 && ifile > 7)
      continue; // dcacascdau
    if (itopovar == 2 && ifile > 5)
      continue; // dcabachtopv
    if (itopovar == 3 && ifile > 7)
      continue; // dcapostopv
    if (itopovar == 4 && ifile > 7)
      continue; // dcanegtopv
    if (itopovar == 5 && ifile > 7)
      continue; // lambdamasswin
    if (itopovar == 6 && ifile > 7)
      continue; // rejcomp
    if (itopovar == 7 && ifile > 4)
      continue; // nsigmatpcKa
    if (itopovar == 8 && ifile > 4)
      continue; // cascradius
    if (itopovar == 9 && ifile > 7)
      continue; // v0radius
    if (itopovar == 10 && ifile > 7)
      continue; // dcav0dau
    if (itopovar == 11 && ifile > 6)
      continue; // v0cospa
    if (itopovar == 12 && ifile > 6)
      continue; // cascproplifetime
    if (itopovar == 13 && ifile > 7)
      continue; // cosbachbaryon
    if (itopovar == 14 && ifile > 7)
      continue; // dcabachbaryon
    if (itopovar == 15 && ifile > 5)
      continue; // dcav0topv
    numFilesEffBis++;
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

    TH1F *histoCuts = (TH1F *)filein[ifile]->Get("CascadeSelectionSummary");
    if (!histoCuts)
    {
      cout << "histoCuts not present" << endl;
      return;
    }
    CutValue[ifile] = histoCuts->GetBinContent(Bin[itopovar]);

    TString inputName = "histo" + TypeHisto[Choice];
    histo[ifile] = (TH1F *)filein[ifile]->Get(inputName);
    if (!histo[ifile])
    {
      cout << "No histo name: " << inputName << " in file " << ifile << endl;
      return;
    }
    if (ifile == 0)
      histoLoosest = (TH1F *)histo[ifile]->Clone("LoosestSelection");

    // Ratios
    histoRatio[ifile] = (TH1F *)histo[ifile]->Clone(Form("hRatio_%i", ifile));
    histoRatio[ifile]->Divide(histoLoosest);
    ErrRatioCorr(histo[ifile], histoLoosest, histoRatio[ifile], 0);

    for (Int_t b = 1; b <= histoRatio[ifile]->GetNbinsX(); b++)
    {
      // cout << "Num: " << histo->GetBinContent(b) << endl;
      // cout << "Denom " << histoLoosest->GetBinContent(b) << endl;
      // cout << "Ratio " << histoRatio->GetBinContent(b) << endl;
    }

    if (Choice == 3)
    {
      YLow[part] = 0;
      if (ifile == 0)
        YUp[part] = 1.2 * histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin());
      else
      {
        if (histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin()) > histo[ifile - 1]->GetBinContent(histo[ifile - 1]->GetMaximumBin()))
        {
          YUp[part] = 1.2 * histo[ifile]->GetBinContent(histo[ifile]->GetMaximumBin());
        }
      }
    }
    if (Choice != 3)
    {
      YLow[part] += 10e-5;
      if (Choice != 2)
        YUp[part] -= 10e-5;
    }

    StyleHisto(histo[ifile], YLow[part], YUp[part], Color[ifile], 33, "", TitleY[Choice], "", 0, 0, 0, 1.5, 1.5, 2);
    StyleHisto(histoRatio[ifile], YLowRatio[Choice], YUpRatio[Choice], Color[ifile], 33, TitleXPt, "Ratio to loosest selection", "", 0, 0, 0, 1.5, 1.5, 2);
    histoRatio[ifile]->GetXaxis()->SetLabelSize(0.08);
    histoRatio[ifile]->GetXaxis()->SetTitleSize(0.08);
    histoRatio[ifile]->GetXaxis()->SetTitleOffset(1.2);
    histoRatio[ifile]->GetYaxis()->SetLabelSize(0.08);
    histoRatio[ifile]->GetYaxis()->SetTitleSize(0.08);
    histoRatio[ifile]->GetYaxis()->SetTitleOffset(0.8);

    canvas->cd();
    pad1->Draw();
    pad1->cd();
    histo[ifile]->Draw("same");
    if (Choice == 0)
      lineMass->DrawClone("same");
    legend->AddEntry(histo[ifile], TopoVarSigned[itopovar] + Form("%.4f", CutValue[ifile]), "pl");
    if (ifile == numFilesEffBis - 1)
      legend->Draw("");

    canvas->cd();
    pad2->Draw();
    pad2->cd();
    if (ifile != 0)
    {
      histoRatio[ifile]->Draw("same");
      lineAt1->Draw("same");
      lineAt09->Draw("same");
    }

    cout << "\nA partire dal file: " << SPathInFinal[ifile] << endl;
  }
  canvas->SaveAs(Sfileout + ".pdf");
  canvas->SaveAs(Sfileout + ".png");
  cout << "\nHo creato il file: " << Sfileout << ".pdf and .png" << endl;
}
