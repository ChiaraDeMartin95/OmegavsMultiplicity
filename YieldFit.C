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
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include </data/dataalice/cdemart/AliPhysicsChiara/AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C>

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
// fits them to get pt-integrated yields

void YieldFit(Int_t typefit = 3, //mT scaling, Boltzmann, Fermi-Direc, Levi
	       Int_t part = 8,
               TString SysPath = "_Sel23June"/*"_Sel6June"*/,
	      Bool_t isBkgParab =1,
               TString OutputDir = "PtIntegratedYields/",
	      TString year = "LHC22o_pass4_Train89684"/*"LHC22m_pass4_Train79153"*/,
               Bool_t isSysStudy = 1,
               Int_t MultType = 1, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0A
               Bool_t UseTwoGauss = 1) 
{

  // multiplicity related variables
  TString Smolt[numMult + 1];
  TString SmoltBis[numMult + 1];

  TString SErrorSpectrum[3] = {"stat.", "syst. uncorr.", "syst. corr."};

  // filein
  TString PathIn;
  TFile *fileIn[numMult + 1];

  //fileinSist
  TString PathInSist;
  PathInSist = "SystematicErrors/TotalSysError_" + year + "_";
  PathInSist += Spart[part];
  // PathInSist += Smolt[numMult];
  //if (isSysStudy) PathInSist += SysPath;
  PathInSist += ".root";
  TFile *fileInSist[numMult + 1];

  // fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = OutputDir + "PtIntegratedYields_" + year;
  stringout += Spart[part];
  stringout += IsOneOrTwoGauss[UseTwoGauss];
  stringout += SIsBkgParab[isBkgParab];

  if (isSysStudy)
    stringout += SysPath;
  stringout += "_"+nameFitFile[typefit];
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
  TH1F *fHistSpectrumStatScaled[numMult + 1];
  TH1F *fHistSpectrumSistScaled[numMult + 1];
  TH1F *fHistSpectrumSistScaledForLegend[numMult + 1];
  TH1F *fHistSpectrumRatioFit[numMult + 1];
  TH1F *fHistSpectrumSistMultRatio[numMult + 1];

  gStyle->SetOptStat(0);

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

  TLegend *LegendPub = new TLegend(0.54, 0.65, 0.95, 0.73);
  LegendPub->SetFillStyle(0);
  LegendPub->SetTextAlign(23);
  LegendPub->SetTextSize(0.025);

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
    fHistSpectrumStat[m]->SetName("histoSpectrumSist_" + Smolt[m]);
    if (!fHistSpectrumStat[m])
    {
      cout << " no hist spectrum stat" << endl;
      return;
    }
    fileInSist[m] = TFile::Open(PathInSist);
    fHistSpectrumSist[m] = (TH1F *)fileInSist[m]->Get("histoSpectrumSist" + Smolt[m]);

    if (!fHistSpectrumSist[m])
    {
      cout << " no hist spectrum sist" << endl;
      return;
    }
  } // end loop on mult

  AliPWGFunc pwgfunc;
  TLegend *legendfit = new TLegend(0.25, 0.25, 0.4, 0.45);
  legendfit->SetFillStyle(0);
  legendfit->SetTextSize(0.04);
  TString namepwgfunc[numMult + 1];

  TLegend *legendfitSummary = new TLegend(0.73, 0.65, 0.92, 0.75);
  legendfitSummary->SetFillStyle(0);
  legendfitSummary->SetTextSize(0.04);
  legendfitSummary->SetTextAlign(33);
  legendfitSummary->AddEntry("", nameFit[typefit] + " fit", "");

  // fit spectra
  Int_t ColorFit[numfittipo + 1] = {860, 881, 868, 628, 419};
  TFitResultPtr fFitResultPtr0[numMult + 1];
  TF1 *fit_pwgfunc[numMult + 1];
  TF1 *fit_pwgfunc_Scaled[numMult + 1];
  TF1 *fit_pwgfunc_ScaledBis[numMult + 1];
  TF1 *fit_pwgfuncBis[numMult + 1];
  TH1 *hhout[numMult + 1];
  Float_t   Yield[numMult + 1] = {0};
  Float_t   YieldExtr[numMult + 1] = {0};
  Float_t   YieldErrStat[numMult + 1] = {0};
  Float_t   YieldErrSistHi[numMult + 1] =  {0};
  Float_t   YieldErrSistLow[numMult + 1] = {0};
  Float_t   Mean[numMult + 1] =  {0};
  Float_t   MeanErrStat[numMult + 1] = {0};
  Float_t   MeanErrSistHi[numMult + 1]= {0};
  Float_t   MeanErrSistLow[numMult + 1]= {0};
  Float_t   Chi2NDF[numMult + 1] = {0};
  Float_t   Temp[numMult + 1] = {0};
  Float_t   TempError[numMult + 1] = {0};

  Double_t LowRange[numMult+1]= {0, 0, 0, 0, 0, 0};
  Double_t UpRange[numMult+1]= {4,4,4,4,4,4};
  if (typefit==3){
    for (Int_t b=0; b<= numMult; b++){
      UpRange[b] = 8;
    }
  }
  for (Int_t b=0; b<= numMult; b++){
    UpRange[b] = 3;
  }

  TString Titlehhout[9] = {"kYield",
                           "kYieldStat",
                           "kYieldSysHi",
                           "kYieldSysLo",
                           "kMean",
                           "kMeanStat",
                           "kMeanSysHi",
                           "kMeanSysLo",
                           "kExtra"};

  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  Int_t factor = 1;

  for (Int_t m = numMult; m >= 0; m--)
  {
    namepwgfunc[m] = Form("fitpwgfunc_m%i_fit%i", m, typefit);
    cout << "\n***********Fitting pt spectra with: " << nameFit[typefit] << endl;
    if (typefit == 0)
      fit_pwgfunc[m] = pwgfunc.GetMTExp(massParticle[part], 0.1, 0.04 * factor, namepwgfunc[m]); // mass, T, norm, name
    if (typefit == 1)
    {
      fit_pwgfunc[m] = pwgfunc.GetBoltzmann(massParticle[part], 0.1, 0.04 * factor, namepwgfunc[m]);
    }
    if (typefit == 2)
      fit_pwgfunc[m] = pwgfunc.GetFermiDirac(massParticle[part], 0.1, 0.04 * factor, namepwgfunc[m]);
    if (typefit == 3)
    {
      fit_pwgfunc[m] = pwgfunc.GetLevi(massParticle[part], 0.1, 0.03, 0.04 * factor, namepwgfunc[m]);                          //norm, n, T, mass (but the function must be called with these parameters in inverse order)
      fit_pwgfunc[m]->SetParLimits(0, 0, fHistSpectrumStat[m]->GetBinContent(fHistSpectrumStat[m]->GetMaximumBin()) * 0.5 * 10); //norm
      fit_pwgfunc[m]->SetParLimits(1, 2, 30);                                                                                    // n
      fit_pwgfunc[m]->SetParLimits(2, 0.1, 10);                                                                                  // T
      fit_pwgfunc[m]->SetParameter(2, 0.7);
    }

    fit_pwgfunc[m]->SetLineColor(ColorMult[m]);
    fit_pwgfunc[m]->SetLineStyle(7);
    fit_pwgfunc[m]->SetLineWidth(2);
    fit_pwgfunc[m]->SetRange(LowRange[m], UpRange[m]);
    fit_pwgfuncBis[m]= (TF1*) fit_pwgfunc[m]->Clone(namepwgfunc[m]+ "_Bis");
    fFitResultPtr0[m] = fHistSpectrumStat[m]->Fit(fit_pwgfuncBis[m],"SR0I");

    cout << "Calling YieldMean macro" << endl;
    //    hhout[m] = YieldMean(fHistSpectrumStat[m], fHistSpectrumSist[m], fit_pwgfunc[m], 0, 20, 0.01, 0.1, "0qI", "log.root", LowRange[m], UpRange[m]);
    hhout[m] = YieldMean(fHistSpectrumStat[m], fHistSpectrumStat[m], fit_pwgfunc[m], 0, 20, 0.01, 0.1, "0qI", "log.root", LowRange[m], UpRange[m]);
    cout << "End of call " << endl;

    hhout[m]->SetLineColor(ColorFit[typefit]);
    hhout[m]->SetName("hhout_" + nameFit[typefit] + Form("_m%i", m));
    hhout[m]->GetYaxis()->SetRangeUser(0, 2);
    for (Int_t b = 1; b <= hhout[m]->GetNbinsX(); b++)
    {
      hhout[m]->GetXaxis()->SetBinLabel(b, Titlehhout[b - 1]);
    }

    cout << "m " << m << " typefit " << typefit << endl;
    Yield[m] = hhout[m]->GetBinContent(1);           // yield (spectra + extrapolated one)
    YieldExtr[m] = hhout[m]->GetBinContent(9);       // extrapolated yield
    YieldErrStat[m] = hhout[m]->GetBinContent(2);    // stat error
    YieldErrSistHi[m] = hhout[m]->GetBinContent(3);  // syst error
    YieldErrSistLow[m] = hhout[m]->GetBinContent(4); // syst error
    Mean[m] = hhout[m]->GetBinContent(5);            // mean pt
    MeanErrStat[m] = hhout[m]->GetBinContent(6);     // stat error
    MeanErrSistHi[m] = hhout[m]->GetBinContent(7);   // syst error
    MeanErrSistLow[m] = hhout[m]->GetBinContent(8);  // syst error
    cout << "************************************" << endl;
    cout << "Multiplicity class: " << SmoltBis[m] << endl;
    cout << "Yield: " << Yield[m] << " +- " << YieldErrStat[m] << " (stat.) +" << YieldErrSistHi[m] << " - " << YieldErrSistLow[m] << " (syst.) " << endl;
    cout << "Mean: " << Mean[m] << " +- " << MeanErrStat[m] << endl;

    Chi2NDF[m] = fit_pwgfunc[m]->GetChisquare()/fit_pwgfunc[m]->GetNDF();
    if (typefit == 3){
      Temp[m] = fit_pwgfunc[m]->GetParameter(2);   
      TempError[m] = fit_pwgfunc[m]->GetParError(2);
    } else {
      Temp[m] = fit_pwgfunc[m]->GetParameter(1);   
      TempError[m] = fit_pwgfunc[m]->GetParError(1);
    }

  } // end loop on multiplicity classes

  // draw spectra in multiplicity classes
  TString sScaleFactorFinal[numMult+1] = {};
  Float_t ScaleFactorFinal[numMult + 1];

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

  Float_t tickX = 0.035;
  Float_t tickY = 0.035;

  TH1F *hDummy = new TH1F("hDummy", "hDummy", 1000, 0, 8);
  canvasPtSpectra->cd();
  SetFont(hDummy);
  StyleHistoYield(hDummy, LimInfSpectra, LimSupSpectra, 1, 1, TitleXPt, TitleYYield, "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummy, tickX, tickY);
  pad1->Draw();
  pad1->cd();
  gPad->SetLogy();
  hDummy->Draw("");

  for (Int_t m = numMult; m >= 0; m--) {
    ScaleFactorFinal[m] = ScaleFactor[m];
    if (m == numMult)      {
      ColorMult[m] = ColorMB;
      MarkerMult[m] = MarkerMB;
      SizeMult[m] = SizeMB;
      ScaleFactorFinal[m] = ScaleFactorMB;
    }
    fHistSpectrumSistScaled[m] = (TH1F *)fHistSpectrumSist[m]->Clone("fHistSpectrumSistScaled_" + Smolt[m]);
    fHistSpectrumStatScaled[m] = (TH1F *)fHistSpectrumStat[m]->Clone("fHistSpectrumStatScaled_" + Smolt[m]);
    fHistSpectrumStatScaled[m]->Scale(ScaleFactorFinal[m]);
    fHistSpectrumSistScaled[m]->Scale(ScaleFactorFinal[m]);
    fit_pwgfunc_Scaled[m] = (TF1*)fit_pwgfunc[m]->Clone(namepwgfunc[m] + "_Scaled");
    fit_pwgfunc_ScaledBis[m] = (TF1*)fit_pwgfuncBis[m]->Clone(namepwgfunc[m] + "_ScaledBis");
    fit_pwgfunc_Scaled[m]->SetLineColor(ColorMult[m]);
    fit_pwgfunc_Scaled[m]->SetParameter(0, fit_pwgfunc[m]->GetParameter(0) * ScaleFactorFinal[m]);
    fit_pwgfunc_ScaledBis[m]->SetParameter(0, fit_pwgfuncBis[m]->GetParameter(0) * ScaleFactorFinal[m]);
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
    fit_pwgfunc_Scaled[m]->Draw("same");
    //fit_pwgfunc_ScaledBis[m]->Draw("same");
    if (m == numMult)
      legendfit->AddEntry(fit_pwgfunc_Scaled[m], nameFit[typefit], "l");
    fHistSpectrumSistScaledForLegend[m] = (TH1F *)fHistSpectrumSistScaled[m]->Clone("fHistSpectrumSistScaledForLegend_" + Smolt[m]);
    sScaleFactorFinal[m] = Form(" (x2^{%i})", int(log2(ScaleFactor[m])));
    legendAllMult->AddEntry(fHistSpectrumSistScaledForLegend[m], SmoltBis[m] + "%" + sScaleFactorFinal[m] + " ", "pef");
  } // end loop on mult
  legendfit->Draw("");
  LegendTitle->Draw("");
  legendAllMult->Draw("");

  // Compute and draw spectra ratios
  Float_t LimSupMultRatio = 1.9;
  Float_t LimInfMultRatio = 1e-2;
  Float_t YoffsetSpectraRatio = 1.1;
  Float_t xTitleR = 25;
  Float_t xOffsetR = 4;
  Float_t yTitleR = 30;
  Float_t yOffsetR = 2;

  Float_t xLabelR = 25;
  Float_t yLabelR = 25;
  Float_t xLabelOffsetR = 0.025;
  Float_t yLabelOffsetR = 0.04;

  Float_t tickXR = 0.045;
  Float_t tickYR = 0.045;

  TF1 * rettaUno= new TF1("rettaUno", "pol0",0,8);
  rettaUno->FixParameter(0,1);

  TString TitleYSpectraRatio = "Data/fit";
  canvasPtSpectra->cd();
  TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 1000, 0, 8);
  SetFont(hDummyRatio);
  StyleHistoYield(hDummyRatio, LimInfMultRatio, LimSupMultRatio, 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
  SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
  SetTickLength(hDummyRatio, tickXR, tickYR);
  padL1->Draw();
  padL1->cd();
  hDummyRatio->Draw("");

  for (Int_t m = numMult; m >= 0; m--)
  {
    fHistSpectrumRatioFit[m] = (TH1F *)fHistSpectrumStat[m]->Clone(Form("HistRatioFit_m%i_typefit%i", m, typefit));
    for (Int_t b=1; b<= fHistSpectrumStat[m]->GetNbinsX(); b++){
      fHistSpectrumRatioFit[m]->SetBinContent(b, fHistSpectrumStat[m]->GetBinContent(b)*fHistSpectrumStat[m]->GetBinWidth(b)/fit_pwgfunc[m]->Integral(fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b), fHistSpectrumStat[m]->GetXaxis()->GetBinUpEdge(b)));
      fHistSpectrumRatioFit[m]->SetBinError(b, fHistSpectrumStat[m]->GetBinError(b)*fHistSpectrumStat[m]->GetBinWidth(b)/fit_pwgfunc[m]->Integral(fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b), fHistSpectrumStat[m]->GetXaxis()->GetBinUpEdge(b)));
  }
    //    fHistSpectrumRatioFit[m]->Divide(fit_pwgfunc[m]);
    rettaUno->SetLineColor(kBlack);
    fHistSpectrumRatioFit[m]->Draw("samee");
    rettaUno->Draw("same");

    SetFont(fHistSpectrumRatioFit[m]);
    StyleHistoYield(fHistSpectrumRatioFit[m], LimInfMultRatio, LimSupMultRatio, ColorMult[m], MarkerMult[m], TitleXPt, TitleYSpectraRatio, "", SizeMultRatio[m], 1.15, YoffsetSpectraRatio);
    SetHistoTextSize(fHistSpectrumRatioFit[m], xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
    SetTickLength(fHistSpectrumRatioFit[m], tickXR, tickYR);

    fHistSpectrumRatioFit[m]->Draw("same e0x0");

  } // end loop on mult


  Float_t LimInfYield = 0;
  Float_t LimSupYield = 0.01;
  Float_t YoffsetYield = 2;

  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 900, 700);
  StyleCanvas(canvasYield, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasChi2 = new TCanvas("canvasChi2", "canvasChi2", 900, 700);
  StyleCanvas(canvasChi2, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasTemp = new TCanvas("canvasTemp", "canvasTemp", 900, 700);
  StyleCanvas(canvasTemp, 0.05, 0.15, 0.2, 0.02);
  TCanvas *canvasFracExtrYield = new TCanvas("canvasFracExtrYield", "canvasFracExtrYield", 900, 700);
  StyleCanvas(canvasFracExtrYield, 0.05, 0.15, 0.2, 0.02);

  TH1F *hYield = new TH1F("hYield", "hYield", numMult+1, 0, numMult+1);
  TH1F *hYieldSist = new TH1F("hYieldSist", "hYieldSist", numMult+1, 0, numMult+1);
  TH1F *hChi2 = new TH1F("hChi2", "hChi2", numMult+1, 0, numMult+1);
  TH1F *hFracExtrYield = new TH1F("hFracExtrYield", "hFracExtrYield", numMult+1, 0, numMult+1);
  TH1F *hTemp = new TH1F("hTemp", "hTemp", numMult+1, 0, numMult+1);
  for (Int_t m = numMult; m >= 0; m--){
    hYield->GetXaxis()->SetBinLabel(m+1, SmoltBis[m] + "%");
    hYield->SetBinContent(m+1, Yield[m]);
    hYield->SetBinError(m+1, YieldErrStat[m]);
    hYieldSist->GetXaxis()->SetBinLabel(m+1, SmoltBis[m] + "%");
    hYieldSist->SetBinContent(m+1, Yield[m]);
    hYieldSist->SetBinError(m+1, 1./2 * (YieldErrSistHi[m] + YieldErrSistLow[m]));
    hChi2->GetXaxis()->SetBinLabel(m+1, SmoltBis[m] + "%");
    hChi2->SetBinContent(m+1, Chi2NDF[m]);
    hChi2->SetBinError(m+1, 0);
    hTemp->GetXaxis()->SetBinLabel(m+1, SmoltBis[m] + "%");
    hTemp->SetBinContent(m+1, Temp[m]);
    hTemp->SetBinError(m+1, TempError[m]);
    hFracExtrYield->GetXaxis()->SetBinLabel(m+1, SmoltBis[m] + "%");
    hFracExtrYield->SetBinContent(m+1, YieldExtr[m]/Yield[m]);
    hFracExtrYield->SetBinError(m+1, 0);
    cout << "\n\n************************************" << endl;
    cout << "Multiplicity class: " << SmoltBis[m] << endl;
    cout << "Yield: " << Yield[m] << " +- " << YieldErrStat[m] << endl;
    cout << "Mean: " << Mean[m] << " +- " << MeanErrStat[m] << endl;

  } // end loop on mult
  
  //taken from HEP data; INEL > 0, Table 11a
  TH1F *hYieldPubStat = new TH1F("hYieldPubStat", "hYieldPubStat", numMult+1, 0, numMult+1);
  TH1F *hYieldPubSist = new TH1F("hYieldPubSist", "hYieldPubSist", numMult+1, 0, numMult+1);
  if (part==6 || part==7) {
    hYieldPubStat->SetBinContent(numMult+1, 0.0024070792/2); 
    hYieldPubSist->SetBinContent(numMult+1, 0.0024070792/2);
  }
  else if (part==8) {
    hYieldPubStat->SetBinContent(numMult+1, 0.0024070792); 
    hYieldPubSist->SetBinContent(numMult+1, 0.0024070792);
  }
  if (part==6 || part==7){
    hYieldPubStat->SetBinError(numMult+1, 0.0000630865/2);
    hYieldPubSist->SetBinError(numMult+1, 0.0002274319/2);
  }
  else if (part==8){
    hYieldPubStat->SetBinError(numMult+1, 0.0000630865);
    hYieldPubSist->SetBinError(numMult+1, 0.0002274319);
  }

  for (Int_t m = numMult; m >= 0; m--){
    if (m==numMult) continue;
    hYieldPubStat->SetBinContent(m+1, -999);
    hYieldPubStat->SetBinError(m+1, -999);
    hYieldPubSist->SetBinContent(m+1, -999);
    hYieldPubSist->SetBinError(m+1, -999);
  }
  
  StyleHistoYield(hYield, LimInfYield, LimSupYield, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldSist, LimInfYield, LimSupYield, 1, 22, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldPubStat, LimInfYield, LimSupYield, kAzure+7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  StyleHistoYield(hYieldPubSist, LimInfYield, LimSupYield, kAzure+7, 33, SMultType[MultType] + " Multiplicity Percentile", TitleYYieldPtInt, "", 2, 1.15, YoffsetYield);
  LegendPub->AddEntry(hYieldPubSist, "Eur.Phys.J.C 80 (2020) 167, 2020", "pl");
  canvasYield->cd();
  hYield->Draw("e");
  hYieldSist->SetFillStyle(0);
  hYieldSist->Draw("same e2");
  hYieldPubStat->Draw("same e0x0");
  hYieldPubSist->SetFillStyle(0);
  hYieldPubSist->Draw("same e2");
  LegendTitle->Draw("");
  LegendPub->Draw("");

  StyleHistoYield(hChi2, 0, 20, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "Chi2/NDF", "", 2, 1.15, YoffsetYield);
  canvasChi2->cd();
  hChi2->Draw("e");
  LegendTitle->Draw("");
  legendfitSummary->Draw("");

  StyleHistoYield(hTemp, 0, 1, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "T parameter", "", 2, 1.15, YoffsetYield);
  canvasTemp->cd();
  hTemp->Draw("e");
  LegendTitle->Draw("");
  legendfitSummary->Draw("");

  StyleHistoYield(hFracExtrYield, 0, 1, 1, 22, SMultType[MultType] + " Multiplicity Percentile", "FracExtrYield", "", 2, 1.15, YoffsetYield);
  canvasFracExtrYield->cd();
  hFracExtrYield->Draw("e");
  LegendTitle->Draw("");
  legendfitSummary->Draw("");

  canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
  canvasYield->SaveAs(stringoutpdf+ "_YieldsvsPerc.pdf");
  canvasChi2->SaveAs(stringoutpdf+ "_Chi2vsPerc.pdf");
  canvasTemp->SaveAs(stringoutpdf+ "_TempvsPerc.pdf");
  canvasFracExtrYield->SaveAs(stringoutpdf+ "_FracExtrYieldvsPerc.pdf");
  canvasPtSpectra->SaveAs(stringoutpdf + ".png");
  canvasYield->SaveAs(stringoutpdf+ "_YieldsvsPerc.png");
  canvasChi2->SaveAs(stringoutpdf+ "_Chi2vsPerc.png");
  canvasTemp->SaveAs(stringoutpdf+ "_TempvsPerc.png");
  canvasFracExtrYield->SaveAs(stringoutpdf+ "_FracExtrYieldvsPerc.png");

  fileout->cd();
  hYield->Write();
  hYieldSist->Write();
  hYieldPubStat->Write();
  hYieldPubSist->Write();
  hChi2->Write();
  hTemp->Write();
  hFracExtrYield->Write();
  fileout->Close();

  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;

  cout << "\nI have created the file:\n " << stringout << endl;
}
