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

TSpline3 *sp3;
Double_t spline(Double_t *x, Double_t *p)
{
    Double_t xx = x[0];
    return sp3->Eval(xx);
}

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

Int_t color = kRed + 2;
Int_t color0 = kRed + 2;
Int_t color1 = kBlue + 2;

void PseudoEfficiency(Int_t ChosenPart = 5,
                      Int_t MultType = 0, // 0: no mult for backward compatibility, 1: FT0M, 2: FV0M
                      Bool_t isMB = 1,
                      Int_t mul = 0,
                      TString year = "LHC22r_pass3_Train67853",
                      TString OutputDir = "Yields",
                      Bool_t UseTwoGauss = 0,
                      Bool_t isBkgParab = 0,
                      Bool_t isMeanFixedPDG = 0,
                      TString SPublishedYieldForPseudoEff = "PublishedYield13TeV/HEPData-ins1748157-v1-Table")
{

    Bool_t isboth = 0;
    if (ChosenPart != 0)
    {
        cout << "Do you want to analyse both particles and antiparticles? (type 1 for both)" << endl;
        cin >> isboth;
    }

    TH1F *histoPub[numPart];
    TH1F *histoPubError[numPart];
    TH1F *histoPubWithError[numPart];
    TH1F *histoNum[numPart];
    TH1F *histoRatioToPub[numPart];
    TDirectoryFile *dir;
    TString FileName[numPart] = {"1", "2", "2", "3", "3", "4", "4"};
    TString HistoNumber[numPart] = {"11", "11", "11", "11", "11", "6", "6"};
    TSpline3 *splinePub[numPart];
    TF1 *fsplinePub[numPart];
    TCanvas *canvasYields[numPart];
    TCanvas *canvasRatioToPub = new TCanvas("canvasRatioToPub", "canvasRatioToPub", 1000, 800);
    StyleCanvas(canvasRatioToPub, 0.15, 0.05, 0.05, 0.15);

    TLegend *legendRatioToPub;
    legendRatioToPub = new TLegend(0.25, 0.76, 0.35, 0.91);
    //legendRatioToPub->SetMargin(0);
    legendRatioToPub->SetTextSize(0.045);

    Float_t YLowRatioToPub = 0;
    Float_t YUpRatioToPub[numPart] = {0.2, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05};
    Float_t LowerPub = 0;
    Float_t UpperPub = 0;

    for (Int_t part = 0; part < numPart; part++)
    {
        if (!isboth && part != ChosenPart)
            continue;
        else if (isboth)
        {
            if (ChosenPart == 1 && (part != 1 && part != 2))
                continue;
            if (ChosenPart == 3 && (part != 3 && part != 4))
                continue;
            if (ChosenPart == 5 && (part != 5 && part != 6))
                continue;
        }

        TString SfileIn;
        SfileIn = OutputDir + "/Yields_" + Spart[part] + "_" + year;
        SfileIn += IsOneOrTwoGauss[UseTwoGauss];
        if (isMB)
            SfileIn += "_Mult0-100";
        else
            SfileIn += Form("_Mult%i-%i", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);
        SfileIn += ".root";

        TFile *fileIn = new TFile(SfileIn, "");
        if (!fileIn)
            return;
        cout << "File Run 3 data: " << SfileIn << endl;
        histoNum[part] = (TH1F *)fileIn->Get("histoYield");
        if (!histoNum[part])
        {
            cout << "histoNum not found " << endl;
            return;
        }

        TString SfilePub = SPublishedYieldForPseudoEff + "_" + FileName[part] + ".root";
        TFile *filePub = new TFile(SfilePub, "");
        if (!filePub)
        {
            cout << "File " << SfilePub << " not available " << endl;
            return;
        }
        cout << "Input file with published yields: " << SfilePub << endl;
        dir = (TDirectoryFile *)filePub->Get("Table " + FileName[part]);
        if (!dir)
        {
            cout << "Input dir not available " << endl;
            return;
        }

        histoPub[part] = (TH1F *)dir->Get("Hist1D_y" + HistoNumber[part]);
        if (!histoPub[part])
        {
            cout << "Published histo not found" << endl;
            return;
        }
        histoPub[part]->SetName("histoYieldPub" + Spart[part]);

        histoPubError[part] = (TH1F *)dir->Get("Hist1D_y" + HistoNumber[part] + "_e1");
        if (!histoPubError[part])
        {
            cout << "Published histo of STAT. ERRORS not found" << endl;
            return;
        }
        histoPubError[part]->SetName("histoYieldPubError" + Spart[part]);

        if (part != 0)
        {
            histoPub[part]->Scale(1. / 2); // particle and antiparticle yields are summed
            histoPubError[part]->Scale(1. / 2);
        }

        histoPubWithError[part] = (TH1F *)histoPub[part]->Clone("histoYieldPubWithError" + Spart[part]);
        for (Int_t b = 1; b <= histoPub[part]->GetNbinsX(); b++)
        {
            if (b == 1)
                LowerPub = histoPubWithError[part]->GetXaxis()->GetBinCenter(b);
            if (b == histoPub[part]->GetNbinsX())
                UpperPub = histoPubWithError[part]->GetXaxis()->GetBinCenter(b);
            histoPubWithError[part]->SetBinError(b, histoPubError[part]->GetBinContent(b));
        }
        splinePub[part] = new TSpline3(histoPubWithError[part], "Spline" + Spart[part]);
        sp3 = (TSpline3 *)splinePub[part]->Clone("SplineClone" + Spart[part]);
        fsplinePub[part] = new TF1("fSpline" + Spart[part], spline, 0, 10);
        for (Int_t b = 1; b <= histoPub[part]->GetNbinsX(); b++)
        {
            // cout << fsplinePub[part]->Eval(histoPub[part]->GetBinCenter(b)) << " vs "
            //<< histoPub[part]->GetBinContent(b) << endl;
        }

        canvasYields[part] = new TCanvas("canvasYields" + Spart[part], "canvasYields" + Spart[part], 1000, 800);
        StyleCanvas(canvasYields[part], 0.15, 0.05, 0.05, 0.15);
        histoPubWithError[part]->SetLineColor(kRed);
        histoPubWithError[part]->SetMarkerColor(kRed);
        histoPubWithError[part]->Draw("same");
        splinePub[part]->Draw("same");
        histoNum[part]->Draw("same");

        // RATIO
        histoRatioToPub[part] = (TH1F *)histoNum[part]->Clone("RatioToPub" + Spart[part]);
        for (Int_t b = 1; b <= histoNum[part]->GetNbinsX(); b++)
        {
            if (histoNum[part]->GetXaxis()->GetBinLowEdge(b) < LowerPub)
            {
                histoRatioToPub[part]->SetBinContent(b, 0);
                histoRatioToPub[part]->SetBinError(b, 0);
                continue;
            }
            if (histoNum[part]->GetXaxis()->GetBinUpEdge(b) > UpperPub)
            {
                histoRatioToPub[part]->SetBinContent(b, 0);
                histoRatioToPub[part]->SetBinError(b, 0);
                continue;
            }
            Float_t ALow = histoNum[part]->GetXaxis()->GetBinLowEdge(b);
            Float_t AUp = histoNum[part]->GetXaxis()->GetBinUpEdge(b);

            if (fsplinePub[part]->Integral(ALow, AUp) <= 10e-8)
            {
                cout << "integral of spline smaller or equal to 10e-8" << endl;
                continue;
            }

            Float_t Numerator = histoNum[part]->GetBinContent(b) * histoNum[part]->GetBinWidth(b);
            Float_t SplineIntegralError = 0;
            // if (part != 0)
            //   histoPubError[part]->Scale(1.2);
            for (Int_t b = 1; b <= histoPubError[part]->GetNbinsX(); b++)
            {
                if (histoPubError[part]->GetBinCenter(b) < ALow)
                    continue;
                if (histoPubError[part]->GetBinCenter(b) > AUp)
                    continue;
                SplineIntegralError += pow(histoPubError[part]->GetBinContent(b), 2);
            }
            if (SplineIntegralError == 0)
            {
                cout << "Spline error mode 2" << endl;
                for (Int_t b = 1; b <= histoPubError[part]->GetNbinsX(); b++)
                {
                    if (histoPubError[part]->GetXaxis()->GetBinLowEdge(b) > AUp)
                        continue;
                    if (histoPubError[part]->GetXaxis()->GetBinUpEdge(b) < ALow)
                        continue;
                    SplineIntegralError += pow(histoPubError[part]->GetBinContent(b), 2);
                }
            }

            SplineIntegralError = sqrt(SplineIntegralError);

            histoRatioToPub[part]->SetBinContent(b, Numerator / fsplinePub[part]->Integral(ALow, AUp));
            histoRatioToPub[part]->SetBinError(b, sqrt(
                                                      pow(histoNum[part]->GetBinError(b) / histoNum[part]->GetBinContent(b), 2) +
                                                      pow(SplineIntegralError / fsplinePub[part]->Integral(ALow, AUp), 2)) *
                                                      histoRatioToPub[part]->GetBinContent(b));

            cout << "\nhistoNum bin center: " << histoNum[part]->GetBinCenter(b) << " Alow: " << ALow << " AUp " << AUp << endl;
            cout << "histoNum bin content: " << histoNum[part]->GetBinContent(b) << " +- " << histoNum[part]->GetBinError(b) << endl;
            cout << "spline integral: " << fsplinePub[part]->Integral(ALow, AUp) << " +- " << SplineIntegralError << endl;
            cout << "histoRatio bin content: " << histoRatioToPub[part]->GetBinContent(b) << " +- " << histoRatioToPub[part]->GetBinError(b) << endl;
        }

        canvasRatioToPub->cd();
        if (isboth)
        {
            if (part == ChosenPart)
                color = color0;
            else
                color = color1;
        }

        StyleHisto(histoRatioToPub[part], YLowRatioToPub, YUpRatioToPub[part], color, 33, TitleXPt, "Ratio to published yield at 13 TeV", "", 0, 0, 0, 1.5, 1.5, 2);
        legendRatioToPub->AddEntry(histoRatioToPub[part], NamePart[part], "pl");
        if (part != ChosenPart)
            legendRatioToPub->AddEntry("", year, "");

        histoRatioToPub[part]->Draw("same");
        legendRatioToPub->Draw("");
    }

    TString Soutputfile = OutputDir + "/PseudoEfficiency";
    if (isboth)
        Soutputfile += SpartType[ChosenPart];
    else
        Soutputfile += Spart[ChosenPart];
    Soutputfile += "_" + year;
    Soutputfile += IsOneOrTwoGauss[UseTwoGauss];
    if (isMB)
        Soutputfile += "_Mult0-100";
    else
        Soutputfile += Form("_Mult%i-%i", MultiplicityPerc[mul], MultiplicityPerc[mul + 1]);

    canvasRatioToPub->SaveAs(Soutputfile + ".pdf");
}
