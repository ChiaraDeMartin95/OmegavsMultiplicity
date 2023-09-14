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

void CalculatedNdeta(TString PathIn = "hCentProfileFT0M.root"/*"profileFT0M.root"*/)
{

  cout << "\nStarting from the files (for the different mult): " << PathIn << endl;
  TFile *file = new TFile(PathIn.Data());
  //TCanvas *c = (TCanvas *)file->Get("profileFT0M");
  TProfile *hprof = (TProfile *)file->Get("hCentProfileFT0M_pub_corr");
  if (!hprof)
  {
    cout << "Error: histogram not found!" << endl;
    return;
  }

  Float_t dndeta[numMult];
  Float_t dndetaMB = 0;
  Float_t dummy[numMult];
  Int_t counter[numMult];
  for (int m = 0; m < numMult; m++)
  {
    if (m != numMult)
      cout << "\nBin " << MultiplicityPerc[m] << " - " << MultiplicityPerc[m + 1] << endl;
    else
      cout << "\nBin " << 0 << " - " << 100 << endl;
    dummy[m] = 1;
    dndeta[m] = 0;
    counter[m] = 0;
    // cout << hprof->GetBinContent(1) << endl;
    for (Int_t i = hprof->FindBin(MultiplicityPerc[m] + 0.00001); i <= hprof->FindBin(MultiplicityPerc[m + 1] - 0.00001); i++)
    {
      // cout << "Bin " << i << " : " << hprof->GetBinWidth(i) << " " << hprof->GetBinCenter(i) << endl;
      if (hprof->GetBinContent(i) == 0)
        continue;
      dndeta[m] += hprof->GetBinContent(i);
      counter[m]++;
    }
    dndeta[m] = dndeta[m] / counter[m];
    cout << "dNdeta = " << dndeta[m] << endl;
  }


  for (int m = 0; m < numMult; m++)
  {
    dndetaMB += dndeta[m] * (MultiplicityPerc[m + 1] - MultiplicityPerc[m]);
  }
  dndetaMB /= 100;
  cout << "dNdetaMB = " << dndetaMB << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  TH1F * hDummy = new TH1F("hDummy", "hDummy", 100, 0, 20);
  hDummy->GetXaxis()->SetTitle("dN/d#eta");
  hDummy->GetYaxis()->SetTitle("Arbitrary value");
  hDummy->GetYaxis()->SetRangeUser(0.9, 1.1);
  TGraph *g = new TGraph(numMult, dndeta, dummy);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.5);
  hDummy->Draw("");
  g->Draw("p");
}
