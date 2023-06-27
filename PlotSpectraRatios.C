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

void PlotSpectraRatios(Int_t part = 8,
                       Int_t ChosenMult = numMult,
                       TString SysPath = "_Sel23June" /*"_Sel6June"*/,
                       TString OutputDir = "PtSpectraMultClasses/",
                       Bool_t isBkgParab = 1,
                       TString year = "LHC22o_pass4_Train89684" /*"LHC22m_pass4_Train79153"*/,
                       Bool_t isSysStudy = 1,
                       Int_t MultType = 1, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
                       Bool_t UseTwoGauss = 1)
{

  gStyle->SetOptStat(0);
  if (ChosenMult > numMult)
  {
    cout << "Chosen Mult outside of available range" << endl;
    return;
  }

  // multiplicity related variables
  TString Smolt[numMult + 1];
  TString SmoltBis[numMult + 1];
  TString sScaleFactorFinal[numMult + 1];
  Float_t ScaleFactorFinal[numMult + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TString PathInSist;
  PathInSist = "SystematicErrors/TotalSysError_" + year + "_";
  PathInSist += Spart[part];
  // PathInSist += Smolt[numMult];
  //if (isSysStudy) PathInSist += SysPath;
  PathInSist += ".root";

  TFile *fileIn[numMult + 1];
  TFile *fileInSist[numMult + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "PlotSpectraRatios_" + year;
  stringout += Spart[part];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[isBkgParab];
  if (isSysStudy)
    stringout += SysPath;
  stringoutpdf = stringout;
  stringout += ".root";
  TFile *fileout = new TFile(stringout, "RECREATE");

  // canvases
  TCanvas *canvasPtSpectra = new TCanvas("canvasPtSpectra", "canvasPtSpectra", 700, 900);
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
  TPad *padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

  StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
  StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

  TH1F *fHistSpectrumStat[numMult + 1];
  TH1F *fHistSpectrumSist[numMult + 1];
  TH1F *fHistSpectrumRelErrSist[numMult + 1];
  TH1F *fHistSpectrumStatScaled[numMult + 1];
  TH1F *fHistSpectrumSistScaled[numMult + 1];
  TH1F *fHistSpectrumStatScaledB[numMult + 1];
  TH1F *fHistSpectrumSistScaledB[numMult + 1];
  TH1F *fHistSpectrumSistScaledForLegend[numMult + 1];
  TH1F *fHistSpectrumStatMultRatio[numMult + 1];
  TH1F *fHistSpectrumSistMultRatio[numMult + 1];

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *legendAllMult = new TLegend(0.22, 0.03, 0.73, 0.28);
  legendAllMult->SetHeader(SMultType[MultType] + " Multiplicity Percentile");
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

  TLegend *legendStatBoxK0s = new TLegend(0.25, 0.33, 0.47, 0.42);
  legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[numMult], "stat. error", "pe");
  legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[numMult], "syst. error", "ef");

  TLine *lineat1Mult = new TLine(0, 1, 8, 1);
  lineat1Mult->SetLineColor(1);
  lineat1Mult->SetLineStyle(2);

  // get spectra in multiplicity classes
  for (Int_t m = numMult; m >= 0; m--)
  {
    PathIn = "Yields/YieldEffCorr" + year + "_";
    PathIn += Spart[part];
    PathIn += IsOneOrTwoGauss[UseTwoGauss];
    PathIn += SIsBkgParab[isBkgParab];
    if (m == numMult)
    {
      Smolt[m] += "_Mult0-100";
      SmoltBis[m] += "0#minus100";
    }
    else
    {
      Smolt[m] += Form("_Mult%.1f-%.1f", MultiplicityPerc[m], MultiplicityPerc[m + 1]);
      SmoltBis[m] += Form("%.1f#minus%.1f", MultiplicityPerc[m], MultiplicityPerc[m + 1]);
    }
    PathIn += Smolt[m];
    if (isSysStudy)
      PathIn += SysPath;
    PathIn += ".root";
    cout << "Path in : " << PathIn << endl;

    fileIn[m] = TFile::Open(PathIn);
    fHistSpectrumStat[m] = (TH1F *)fileIn[m]->Get("histoYieldCorr");
    fHistSpectrumStat[m]->SetName("histoSpectrumStat_" + Smolt[m]);
    if (!fHistSpectrumStat[m])
    {
      cout << " no hist spectrum stat" << endl;
      return;
    }

    cout << "PathInSist: " << PathInSist << endl;
    fileInSist[m] = TFile::Open(PathInSist);
    // fHistSpectrumRelErrSist[m] = (TH1F *)fileInSist[m]->Get("hTotalSyst");
    fHistSpectrumRelErrSist[m] = (TH1F *)fileInSist[m]->Get("histoSpectrumRelErrSist" + Smolt[m]);
    fHistSpectrumRelErrSist[m]->SetName("histoRelSist_" + Smolt[m]);

    if (!fHistSpectrumRelErrSist[m])
    {
      cout << " no hist spectrum rel err sist" << endl;
      return;
    }
    fHistSpectrumSist[m] = (TH1F *)fHistSpectrumStat[m]->Clone("histoSpectrumSist_" + Smolt[m]);
    for (Int_t i = 1; i <= fHistSpectrumSist[m]->GetNbinsX(); i++)
    {
      fHistSpectrumSist[m]->SetBinError(i, fHistSpectrumStat[m]->GetBinContent(i) * fHistSpectrumRelErrSist[m]->GetBinContent(i));
    }
  } // end loop on mult

  // draw spectra in multiplicity classes
  // Float_t LimSupSpectra = 9.99;
  Float_t LimSupSpectra = 9999.99;
  Float_t LimInfSpectra = 0.2 * 1e-7;
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
  StyleHistoYield(hDummy, LimInfSpectra, LimSupSpectra, 1, 1, TitleXPt, TitleYYield, "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  hDummy->GetXaxis()->SetRangeUser(0, 8);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();
  hDummy->Draw("same");

  for (Int_t m = numMult; m >= 0; m--)
  {
    ScaleFactorFinal[m] = ScaleFactor[m];
    if (m == numMult)
    {
      ColorMult[m] = ColorMB;
      MarkerMult[m] = MarkerMB;
      SizeMult[m] = SizeMB;
      ScaleFactorFinal[m] = ScaleFactorMB;
    }
    fHistSpectrumSistScaled[m] = (TH1F *)fHistSpectrumSist[m]->Clone("fHistSpectrumSistScaled_" + Smolt[m]);
    fHistSpectrumStatScaled[m] = (TH1F *)fHistSpectrumStat[m]->Clone("fHistSpectrumStatScaled_" + Smolt[m]);
    fHistSpectrumStatScaled[m]->Scale(ScaleFactorFinal[m]);
    fHistSpectrumSistScaled[m]->Scale(ScaleFactorFinal[m]);
    for (Int_t b = 1; b <= fHistSpectrumStat[m]->GetNbinsX(); b++)
    {
      // cout << "bin " << b << " " << fHistSpectrumStat[m]->GetBinContent(b) << "+-" << fHistSpectrumStat[m]->GetBinError(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumSist[m]->GetBinContent(b) << "+-" << fHistSpectrumSist[m]->GetBinError(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumStatScaled[m]->GetBinContent(b) << "+-" << fHistSpectrumStatScaled[m]->GetBinError(b) << endl;
    }
    fHistSpectrumStatScaled[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumStatScaled[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumStatScaled[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumStatScaled[m]->SetMarkerSize(SizeMult[m]);
    fHistSpectrumSistScaled[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumSistScaled[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumSistScaled[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumSistScaled[m]->SetMarkerSize(SizeMult[m]);
    fHistSpectrumStatScaled[m]->Draw("same e0x0");
    fHistSpectrumSistScaled[m]->SetFillStyle(0);
    fHistSpectrumSistScaled[m]->Draw("same e2");
    fHistSpectrumSistScaledForLegend[m] = (TH1F *)fHistSpectrumSistScaled[m]->Clone("fHistSpectrumSistScaledForLegend_" + Smolt[m]);
    sScaleFactorFinal[m] = Form(" (x2^{%i})", int(log2(ScaleFactorFinal[m])));
    legendAllMult->AddEntry(fHistSpectrumSistScaledForLegend[m], SmoltBis[m] + "%" + sScaleFactorFinal[m] + " ", "pef");
  } // end loop on mult
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  // Compute and draw spectra ratios
  // Float_t LimSupMultRatio = 5.1;
  Float_t LimSupMultRatio = 9;
  Float_t LimInfMultRatio = 1.1 * 1e-1;
  Float_t YoffsetSpectraRatio = 1.1;
  Float_t xTitleR = 35;
  Float_t xOffsetR = 1;
  Float_t yTitleR = 30;
  Float_t yOffsetR = 2;

  Float_t xLabelR = 25;
  Float_t yLabelR = 25;
  Float_t xLabelOffsetR = 0.02;
  Float_t yLabelOffsetR = 0.04;

  TString TitleYSpectraRatio = "Ratio to " + SmoltBis[ChosenMult] + "%";
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, 8);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, LimInfMultRatio, LimSupMultRatio, 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickX, tickY);
  hDummyRatio->GetXaxis()->SetRangeUser(0, 8);
  canvasPtSpectra->cd();
  padL1->Draw();
  padL1->cd();
  gPad->SetLogy();
  hDummyRatio->Draw("same");

  for (Int_t m = numMult; m >= 0; m--)
  {
    fHistSpectrumStatMultRatio[m] = (TH1F *)fHistSpectrumStat[m]->Clone("fHistSpectrumStatMultRatio_" + Smolt[m]);
    fHistSpectrumSistMultRatio[m] = (TH1F *)fHistSpectrumSist[m]->Clone("fHistSpectrumSistMultRatio_" + Smolt[m]);
    fHistSpectrumStatMultRatio[m]->Divide(fHistSpectrumStat[ChosenMult]);
    fHistSpectrumSistMultRatio[m]->Divide(fHistSpectrumSist[ChosenMult]);
    ErrRatioCorr(fHistSpectrumStat[m], fHistSpectrumStat[ChosenMult], fHistSpectrumStatMultRatio[m], 0);
    ErrRatioCorr(fHistSpectrumSist[m], fHistSpectrumSist[ChosenMult], fHistSpectrumSistMultRatio[m], 0);
    fHistSpectrumStatMultRatio[m]->GetYaxis()->SetRangeUser(LimInfMultRatio, LimSupMultRatio);
    fHistSpectrumSistMultRatio[m]->GetYaxis()->SetRangeUser(LimInfMultRatio, LimSupMultRatio);
    for (Int_t b = 1; b <= fHistSpectrumStat[m]->GetNbinsX(); b++)
    {
      // cout << "bin " << b << " " << fHistSpectrumStat[m]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumStat[ChosenMult]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumStatMultRatio[m]->GetBinContent(b) << endl;
      // cout << "bin " << b << " " << fHistSpectrumSistMultRatio[m]->GetBinError(b) << endl;
    }
    fHistSpectrumStatMultRatio[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumStatMultRatio[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumStatMultRatio[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumStatMultRatio[m]->SetMarkerSize(SizeMultRatio[m]);
    fHistSpectrumSistMultRatio[m]->SetMarkerColor(ColorMult[m]);
    fHistSpectrumSistMultRatio[m]->SetLineColor(ColorMult[m]);
    fHistSpectrumSistMultRatio[m]->SetMarkerStyle(MarkerMult[m]);
    fHistSpectrumSistMultRatio[m]->SetMarkerSize(SizeMultRatio[m]);

    if (m != ChosenMult)
    {
      fHistSpectrumStatMultRatio[m]->Draw("same e0x0");
      fHistSpectrumSistMultRatio[m]->SetFillStyle(0);
      fHistSpectrumSistMultRatio[m]->Draw("same e2");
    }

  } // end loop on mult
  fileout->Close();
  canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
  canvasPtSpectra->SaveAs(stringoutpdf + ".png");
  cout << "\nStarting from the files (for the different mult): " << PathIn << " and PathInSist: " << PathInSist << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
  cout << "WARNING!! file of syst error is always the 0-100% one" << endl;
}
