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
#include "InputVar.h"

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

// take spectra in input
// produces ratio of spectra wrt 0-100% multiplciity class

Float_t YLowMean[numPart] = {0.485, 1.110, 1.110, 1.316, 1.316, 1.316, 1.66, 1.66, 1.66};
Float_t YUpMean[numPart] = {0.51, 1.130, 1.130, 1.327, 1.327, 1.327, 1.68, 1.68, 1.68};
Float_t YLowSigma[numPart] = {0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0, 0.0, 0.0};
Float_t YUpSigma[numPart] = {0.025, 0.015, 0.015, 0.015, 0.015, 0.015, 0.008, 0.008, 0.008};
Float_t YLowPurity[numPart] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

Float_t YLow[numPart] = {0};
Float_t YUp[numPart] = {0};

Float_t YLowRatio[numChoice] = {0.99, 0.4, 0.8, 0.9, 0.8, 0.8, 0.8};
Float_t YUpRatio[numChoice] = {1.01, 1.6, 1.2, 1.2, 1.2, 1.2, 1.2};

const Int_t numRuns = 13;
// multiplicity related variables
// LHC22o_pass4
// TString Srun[numRuns+1] = {"_Run528534", "_Run528531", "_Run528463", "_Run528461", "_Run528448", "_Run528381", "_Run528379", "_Run528232", ""};
// TString SrunBis[numRuns+1] = {"528534", "528531", "528463", "528461", "528448", "528381", "528379", "528232", "All runs"};

// LHC22o_pass4_MinBias
TString SrunBis[numRuns + 1] = {"528531", "528461", "528292", "527899", "527895", "527871", "527850", "527240", "527109", "527057", "527041", "526964", "526641", "All runs"};
TString Srun[numRuns + 1] = {"_Run528531", "_Run528461", "_Run528292", "_Run527899", "_Run527895", "_Run527871", "_Run527850", "_Run527240", "_Run527109", "_Run527057", "_Run527041", "_Run526964", "_Run526641", ""};

void CorrectionFactorAnch(Int_t part = 8,
                          Int_t ChosenRun = 13,
                          Int_t Choice = 0,
                          Bool_t isBkgParab = ExtrisBkgParab,
                          TString SysPath = "",
                          TString OutputDir = "RunByRunComparison/",
                          TString year = Extryear,
                          Bool_t isSysStudy = 1,
                          Int_t MultType = ExtrMultType, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
                          Bool_t UseTwoGauss = ExtrUseTwoGauss,
                          Int_t evFlag = ExtrevFlag // 0: INEL - 1; 1: INEL > 0; 2: INEL > 1
)
{

  TH1F *fHistSpectrumAnch;
  if (part == 3 || part == 4 || part == 5)
    SysPath = ExtrSysPathXi;
  else if (part == 6 || part == 7 || part == 8)
    SysPath = ExtrSysPathOmega;

  gStyle->SetOptStat(0);
  if (ChosenRun > numRuns)
  {
    cout << "Chosen Mult outside of available range" << endl;
    return;
  }
  cout << Choice << " " << TypeHisto[Choice] << endl;
  if (Choice > (numChoice - 1))
  {
    cout << "Option not implemented" << endl;
    return;
  }
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
  else if (Choice == 3)
  {
    if (part == 6 || part == 7 || part == 8)
    {
      YUp[part] = 1e-4;
      YLow[part] = 1e-9;
    }
    else if (part == 3 || part == 4 || part == 5)
    {
      YUp[part] = 1e-3;
      YLow[part] = 1e-7;
    }
  }

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TFile *fileIn[numRuns + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "AnchoringFactor_" + year;
  stringout += "_" + TypeHisto[Choice];
  stringout += "_" + Spart[part];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[isBkgParab];
  if (isSysStudy)
    stringout += SysPath;
  stringoutpdf += "_" + EventType[evFlag];
  stringoutpdf = stringout;
  stringout += ".root";

  // canvases
  TCanvas *canvasPtSpectra = new TCanvas("canvasPtSpectra", "canvasPtSpectra", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TH1F *fHistSpectrum[numRuns + 1];
  TH1F *fHistSpectrumScaled[numRuns + 1];
  TString sScaleFactorFinal[numRuns + 1];
  TH1F *fHistSpectrumMultRatio[numRuns + 1];

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendAllMult = new TLegend(0.22, 0.03, 0.73, 0.28);
  legendAllMult->SetHeader("Runs");
  legendAllMult->SetNColumns(3);
  legendAllMult->SetFillStyle(0);
  TLegendEntry *lheaderAllMult = (TLegendEntry *)legendAllMult->GetListOfPrimitives()->First();
  lheaderAllMult->SetTextSize(0.04);

  TLegend *LegendTitle = new TLegend(0.54, 0.75, 0.95, 0.92);
  LegendTitle->SetFillStyle(0);
  LegendTitle->SetTextAlign(33);
  LegendTitle->SetTextSize(0.04);
  LegendTitle->AddEntry("", "#bf{ALICE Work In Progress}", "");
  LegendTitle->AddEntry("", "pp, #sqrt{#it{s}} = 13.6 TeV", "");
  LegendTitle->AddEntry("", NamePart[part] + ", |y| < 0.5", "");

  TLine *lineat1Mult = new TLine(0, 1, 8, 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  // get spectra in multiplicity classes
  for (Int_t m = numRuns; m >= 0; m--)
  {
    if (SrunBis[m] != "526641" && SrunBis[m] != "527041" && SrunBis[m] != "All runs")
      continue;
    PathIn = "Yields/Yields_";
    PathIn += Spart[part];
    PathIn += "_" + year;
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[isBkgParab];
    PathIn += "_Mult0-100";
    if (isSysStudy)
      PathIn += SysPath;
    // PathIn += "_Run" + Srun[m];
    PathIn += Srun[m];
    PathIn += "_" + EventType[evFlag];
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;

    fileIn[m] = TFile::Open(PathIn);
    fHistSpectrum[m] = (TH1F *)fileIn[m]->Get("histo" + TypeHisto[Choice]);
    fHistSpectrum[m]->SetName("histoSpectrum_" + SrunBis[m]);
    if (!fHistSpectrum[m])
    {
      cout << " no hist " << endl;
      return;
    }
    if (SrunBis[m] == "526641")
    {
      fHistSpectrumAnch = (TH1F *)fHistSpectrum[m]->Clone("fHistSpectrumAnc_");
      fHistSpectrumAnch->Scale(0.75);
    }
    else if (SrunBis[m] == "527041")
    {
      fHistSpectrum[m]->Scale(0.25);
      fHistSpectrumAnch->Add(fHistSpectrum[m]);
      fHistSpectrum[m]->Scale(1. / 0.25);
    }
  } // end loop on mult

  // draw spectra in multiplicity classes
  Float_t xTitle = 15;
  Float_t xOffset = 4;
  Float_t yTitle = 30;
  Float_t yOffset = 2;

  Float_t xLabel = 30;
  Float_t yLabel = 30;
  Float_t xLabelOffset = 0.05;
  Float_t yLabelOffset = 0.01;

  Float_t tickX = 0.03;
  Float_t tickY = 0.042;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 8);
  for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
    hDummy->SetBinContent(i, 1e-12);
  canvasPtSpectra->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, YLow[part], YUp[part], 1, 1, TitleXPt, TitleY[Choice], "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, 8);
  pad1->Draw();
  pad1->cd();
  if (Choice == 3)
    gPad->SetLogy();
  hDummy->Draw("same");

  for (Int_t m = numRuns; m >= 0; m--)
  {
    if (SrunBis[m] != "526641" && SrunBis[m] != "527041" && SrunBis[m] != "All runs")
      continue;
    fHistSpectrumScaled[m] = (TH1F *)fHistSpectrum[m]->Clone("fHistSpectrumScaled_" + SrunBis[m]);
    for (Int_t b = 1; b <= fHistSpectrum[m]->GetNbinsX(); b++)
    {
      cout << "bin " << b << " " << fHistSpectrum[m]->GetBinContent(b) << "+-" << fHistSpectrum[m]->GetBinError(b) << endl;
      cout << "bin " << b << " " << fHistSpectrumScaled[m]->GetBinContent(b) << "+-" << fHistSpectrumScaled[m]->GetBinError(b) << endl;
    }
    SetFont(fHistSpectrumScaled[m]);
    fHistSpectrumScaled[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumScaled[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumScaled[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumScaled[m]->SetMarkerSize(0.6 * SizeMult[m]);
    fHistSpectrumScaled[m]->GetYaxis()->SetRangeUser(YLow[part], YUp[part]);
    fHistSpectrumScaled[m]->Draw("same e0x0");
    sScaleFactorFinal[m] = "";
    legendAllMult->AddEntry(fHistSpectrumScaled[m], SrunBis[m] + sScaleFactorFinal[m] + " ", "pef");
  } // end loop on mult
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  // Compute and draw spectra ratios
  Float_t LimSupMultRatio = 5.1;
  Float_t LimInfMultRatio = 1e-2;
  Float_t YoffsetSpectraRatio = 1.1;
  Float_t xTitleR = 35;
  Float_t xOffsetR = 1;
  Float_t yTitleR = 30;
  Float_t yOffsetR = 2;

  Float_t xLabelR = 25;
  Float_t yLabelR = 25;
  Float_t xLabelOffsetR = 0.02;
  Float_t yLabelOffsetR = 0.04;

  TString TitleYSpectraRatio = "Ratio to " + SrunBis[ChosenRun];
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, 8);
  for (Int_t i = 1; i <= hDummyRatio->GetNbinsX(); i++)
    hDummyRatio->SetBinContent(i, 1e-12);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, YLowRatio[Choice], YUpRatio[Choice], 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickX, tickY);
  hDummyRatio->GetXaxis()->SetRangeUser(0, 8);
  canvasPtSpectra->cd();
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("same");
  lineat1Mult->Draw("same");

  TF1 *pol0Anch = new TF1("pol0Anch", "pol0", 0.8, 8);
  TH1F *hFinalAnch;
  for (Int_t m = numRuns; m >= 0; m--)
  {
    if (SrunBis[m] != "All runs")
      continue;
    fHistSpectrumMultRatio[m] = (TH1F *)fHistSpectrumAnch->Clone("fHistSpectrumMultRatio");
    fHistSpectrumMultRatio[m]->Divide(fHistSpectrum[ChosenRun]);
    ErrRatioCorr(fHistSpectrumAnch, fHistSpectrum[ChosenRun], fHistSpectrumMultRatio[m], 0);
    for (Int_t b = 1; b <= fHistSpectrum[m]->GetNbinsX(); b++)
    {
      // cout << "bin " << b << " " << fHistSpectrum[m]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrum[ChosenRun]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumMultRatio[m]->GetBinContent(b) << endl;
    }
    fHistSpectrumMultRatio[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumMultRatio[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumMultRatio[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumMultRatio[m]->SetMarkerSize(SizeMultRatio[m]);
    // fHistSpectrumMultRatio[m]->Smooth();
    fHistSpectrumMultRatio[m]->Draw("same ep");
    fHistSpectrumMultRatio[m]->Fit(pol0Anch, "R+");
    hFinalAnch = (TH1F *)fHistSpectrumMultRatio[m]->Clone("hFinalAnch");
    for (Int_t b = 1; b <= hFinalAnch->GetNbinsX(); b++)
    {
      if (hFinalAnch->GetBinCenter(b) > 2) hFinalAnch->SetBinContent(b, pol0Anch->GetParameter(0));
      hFinalAnch->SetBinError(b, 0);
    }

  } // end loop on mult

  TFile *fileout = new TFile(stringout, "RECREATE");
  fileout->WriteTObject(hFinalAnch);
  fileout->Close();
  canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
  canvasPtSpectra->SaveAs(stringoutpdf + ".png");
  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
