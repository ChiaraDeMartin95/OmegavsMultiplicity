
double ErrorInRatio(Double_t A, Double_t Aerr, Double_t B, Double_t Berr);
void DivideAndComputeRogerBarlow(TH1F *h1, TH1F *h2);
double PassRogerBarlowCriterion(int nsigmas, Double_t dev, Double_t RBsigma);
TH1F *makeSystPlotsV0s(int num = 1);

void MultiTrial()
{

    const int trials = 500;
    TH1F *h[trials];
    for (int i = 0; i < trials; i++){
      h[i] = makeSystPlotsV0s(i+1);
    }

    const int bins = h[0]->GetNbinsX();
    TH1F* hPtDev[bins];
    TF1* fgaus[bins];

    for (int i = 0; i < bins; i++){
      //some pt bins have different binning for the deviation histo
      if (i==2) hPtDev[i] = new TH1F(Form("hPtDev%i", i), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(i+1), h[0]->GetBinLowEdge(i+2)), 15, -0.5, +0.5);
      else if (i==3) hPtDev[i] = new TH1F(Form("hPtDev%i", i), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(i+1), h[0]->GetBinLowEdge(i+2)), 15, -0.5, +0.5);
      else if (i==4) hPtDev[i] = new TH1F(Form("hPtDev%i", i), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(i+1), h[0]->GetBinLowEdge(i+2)), 20, -0.5, +0.5);
      else if (i==(bins-1)) hPtDev[i] = new TH1F(Form("hPtDev%i", i), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(i+1), h[0]->GetBinLowEdge(i+2)), 25, -0.5, +0.5);
      else hPtDev[i] = new TH1F(Form("hPtDev%i", i), Form("p_{T} bin [%.1f-%.1f] GeV/c;Y_{sys}/Y_{def} - 1;Counts", h[0]->GetBinLowEdge(i+1), h[0]->GetBinLowEdge(i+2)), 50, -0.5, +0.5);
      hPtDev[i]->SetStats(0);
      fgaus[i] = new TF1(Form("fgaus%i", i), "gaus", -0.5, +0.5);
    }

    for (int i = 1; i < trials; i++){
      for (int j = 0; j < bins; j++){
        hPtDev[j]->Fill(h[i]->GetBinContent(j+1));
      }
    }

    for (int i = 0; i < bins; i++){
      fgaus[i]->SetParameter(0, hPtDev[i]->GetMaximum());
      hPtDev[i]->Fit(fgaus[i], "R");
    }

    TH1F *hSystMultiTrial = (TH1F *)h[0]->Clone("hSystMultiTrial");
    hSystMultiTrial->Reset();
    for (int i = 0; i < bins; i++) {
      hSystMultiTrial->SetBinContent(i + 1, fgaus[i]->GetParameter(2));
      hSystMultiTrial->SetBinError(i + 1, 0);
    }

    TFile* Write = new TFile("SystMultiTrial.root", "RECREATE");
    for (int i = 0; i < bins; i++){
      hPtDev[i]->Write();
      fgaus[i]->Write();
    }
    hSystMultiTrial->Write();

    TLatex *ltx = new TLatex();
    ltx->SetTextSize(0.05);
    ltx->SetTextColor(kRed);

    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1200);
    c2->Divide(3, 4);
    for (int i = 2; i < bins; i++){
      c2->cd(i-2+1);
      c2->cd(i - 2 + 1)->SetBottomMargin(0.15);
      hPtDev[i]->GetYaxis()->SetRangeUser(0., hPtDev[i]->GetMaximum() * 1.2);
      hPtDev[i]->GetXaxis()->SetTitleSize(0.06);
      hPtDev[i]->GetXaxis()->SetTitleOffset(1.);
      hPtDev[i]->Draw("EP");
      hPtDev[i]->SetMarkerStyle(kFullCircle);
      hPtDev[i]->Draw("EP SAME");
      fgaus[i]->Draw("same");

      ltx->DrawLatexNDC(0.6, 0.8, Form("#mu = %.3f", fgaus[i]->GetParameter(1)));
      ltx->DrawLatexNDC(0.6, 0.7, Form("#sigma = %.3f", fgaus[i]->GetParameter(2)));
      ltx->DrawLatexNDC(0.6, 0.6, Form("#chi^{2}/ndf = %.1f", fgaus[i]->GetChisquare() / fgaus[i]->GetNDF()));
    }

    c2->SaveAs("images/MultiTrial.pdf");

}

//-------------------------------------------------------------------------------
//------------------------- DEFINE FUNCTIONS ------------------------------------
//-------------------------------------------------------------------------------

//---------------------------------------------------------------
TH1F* makeSystPlotsV0s(int num = 1)
{

    TFile *fdef = TFile::Open("CorrSpectra/YieldEffCorrOmega_Mult0100_Default.root");
    TFile *fvaried = TFile::Open(Form("CorrSpectra/YieldEffCorrOmega_Mult0100_%i.root", num));

    TH1F *hVariedCut = (TH1F *)fvaried->Get("histoYieldCorr");
    hVariedCut->SetName(Form("hVarCutSet%i", num));
    TH1F *hDefault = (TH1F *)fdef->Get("histoYieldCorr");
    hDefault->SetName(Form("hDefault"));

    TH1F *hDev = (TH1F *)hDefault->Clone("hDev");
    hDev->Reset();

    DivideAndComputeRogerBarlow(hVariedCut, hDefault);

    for (int i = 1; i <= hDev->GetNbinsX(); i++)
    {
      double dev = hVariedCut->GetBinContent(i) - 1;
      double err = hVariedCut->GetBinError(i);

      hDev->SetBinContent(i, PassRogerBarlowCriterion(1, dev, err));
  }

  hDev->GetYaxis()->SetRangeUser(0., 0.5);

  //Return Max Dev Histo
  return hDev;

}

//---------------------------------------------------------------------------------------------------
double ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ) {
    //Error in a Ratio
    if(B!=0) {
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( TMath::Abs(errorfromtop - errorfrombottom) );
    }
    return 1.;
}

//----------------------------------------------------------------------------------------------------
void DivideAndComputeRogerBarlow( TH1F* h1, TH1F *h2 ){
  //Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX();
  Double_t lh2NBins = h2->GetNbinsX();

  if( lh1NBins != lh2NBins ){
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return;
  }

  Double_t lSigmaDelta[100];
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    //Computation of roger barlow sigma_{delta}
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    //Computation of relationship to h2 for plotting in ratio plot
    if ( h2->GetBinContent(i) > 1e-12 ){
      lSigmaDelta[i] /= h2->GetBinContent(i);
    }else{
      lSigmaDelta[i] = 0;
    }
  }
  //Regular Division
  h1 -> Divide( h2 );
  //Replace Erorrs
  for( Int_t i=1; i<h1->GetNbinsX()+1; i++){
    h1->SetBinError(i, lSigmaDelta[i]);
  }
}

//----------------------------------------------------------------------------------------------------
double PassRogerBarlowCriterion(int nsigmas, Double_t dev, Double_t RBsigma){

  if (TMath::Abs(dev)>(nsigmas*RBsigma)) {return dev;}
  else {return 0.;}
}