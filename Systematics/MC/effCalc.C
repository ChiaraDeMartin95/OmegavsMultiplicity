#include "effCalc.h"

/*
.L effCalc.C

// Xi

effCalc("/Users/rnepeiv/workLund/PhD_work/run3omega/postprocessingo2physics/results/chiara/xi/lhc23f4b2/AnalysisResults.root", "/Users/rnepeiv/workLund/PhD_work/run3omega/efficiencyCalculation/results/effChiaraXi_inelgt0_lhc23f4b2_24aug.root", 0, 1);

// Omega

effCalc("/Users/rnepeiv/workLund/PhD_work/run3omega/postprocessingo2physics/results/chiara/omega/lhc23f4b2/AnalysisResults.root", "/Users/rnepeiv/workLund/PhD_work/run3omega/efficiencyCalculation/results/effChiaraOmega_inelgt0_lhc23f4b2_24aug.root", 1, 1);

effCalc("/Users/rnepeiv/workLund/PhD_work/run3omega/postprocessingo2physics/results/chiara/omega/gapTriggered27aug/AnalysisResults.root", "/Users/rnepeiv/workLund/PhD_work/run3omega/efficiencyCalculation/results/effChiaraOmega_inelgt0_gt13tev_27aug.root", 1, 1);
*/

TFile* FindFileFresh(const Char_t* fileName){
  // Find file
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName);
  if(file) {
    file->Close();
    delete file;
  }

  file = TFile::Open(fileName, "READ");

  if(!file)
    cout << "File : " << fileName << " was not found" << endl;

  return file;
}

void PrintBinning(TH1* hist)
{
  const Int_t nBinsX = hist->GetXaxis()->GetNbins();

  cout << "const Int_t nBins = " << nBinsX << ";" << endl
       << "Double_t bins[nBins+1] = " << endl
       << "{";
  for(Int_t bin = 1; bin <= nBinsX; bin++) {
    
    printf("%.2f, ", hist->GetXaxis()->GetBinLowEdge(bin));
    if (bin%10 == 0) 
      cout << endl;
  }
  
  printf("%.2f", hist->GetXaxis()->GetBinUpEdge(nBinsX));
  cout << "}" << endl;  
}

void SetErrorDivide(TH1D* histFinal, TH1D* hist1, TH1D* hist2)
{
  const Int_t nBins = hist1->GetXaxis()->GetNbins();
  for(Int_t bin = 1; bin <= nBins; bin++) {
    Int_t k = hist1->GetBinContent(bin);
    Int_t n = hist2->GetBinContent(bin);
    Double_t errorInBin = sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
    histFinal->SetBinError(bin, errorInBin);
  }
}

void effCalc(const Char_t* fileNameData, const Char_t* outputFileName, Bool_t isOmega = 1, Bool_t is13tevBinning = 1){
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat("me");
  // Get Data
  fileData = FindFileFresh(fileNameData);

  // Create output file
  fileOut = new TFile(outputFileName, "RECREATE");

  if(isOmega){
    hPtCascadePlusTrueAssoc = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtOmegaPlusTrueAssocWithSelColl");
    hPtCascadeMinusTrueAssoc = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtOmegaMinusTrueAssocWithSelColl");
    hPtCascadePlusTrue = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtOmegaPlusTrue");
    hPtCascadeMinusTrue = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtOmegaMinusTrue");
  } else {
    hPtCascadePlusTrueAssoc = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtXiPlusTrueAssocWithSelColl");
    hPtCascadeMinusTrueAssoc = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtXiMinusTrueAssocWithSelColl");
    hPtCascadePlusTrue = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtXiPlusTrue");
    hPtCascadeMinusTrue = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtXiMinusTrue");
  }

  hPtCascadePlusTrueRec = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtCascPlusTrueRec");
  hPtCascadeMinusTrueRec = (TH3F*)fileData->Get("lf-cascpostprocessing/hPtCascMinusTrueRec");

  hPtCascadeSumAssoc = (TH3F*)hPtCascadePlusTrueAssoc->Clone("SumAssoc");
  hPtCascadeSumAssoc->Add(hPtCascadeMinusTrueAssoc);
  hPtCascadeSumAssoc->SetName("hPtSumTrueAssoc");

  hPtCascadeRecSum = (TH3F*)hPtCascadePlusTrueRec->Clone("SumRec");
  hPtCascadeRecSum->Add(hPtCascadeMinusTrueRec);
  hPtCascadeRecSum->SetName("hPtRecSum");

  hPtCascadeSumTrue = (TH3F*)hPtCascadePlusTrue->Clone("SumTrue");
  hPtCascadeSumTrue->Add(hPtCascadeMinusTrue);
  hPtCascadeSumTrue->SetName("hPtSumTrue");

  TDirectory* dir = fileOut->mkdir("supportHistos");
  dir->cd();

  Double_t leftCentrFT0 = 0;
  Double_t rightCentrFT0 = 300;

  Double_t leftRap = -0.5;
  Double_t rightRap = 0.5;

  // True lvl Sum Not Assoc
  TH3F* hPtCascadeCentrSumTrue = (TH3F*)hPtCascadeSumTrue->Clone("SumTrueClone");
  hPtCascadeCentrSumTrue->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrSumTrue->Write();

  TH2F* hProfilePtCascadeZSumTrue = static_cast<TH2F*>(hPtCascadeCentrSumTrue->Project3D("xy"));
  hProfilePtCascadeZSumTrue->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZSumTrue->Write();

  TH1D* hProfilePtCascadeYZSumTrue = (hProfilePtCascadeZSumTrue->ProjectionY());
  hProfilePtCascadeYZSumTrue->Sumw2();
  hProfilePtCascadeYZSumTrue->Write();

  // True lvl Plus Not Assoc
  TH3F* hPtCascadeCentrPlusTrue = (TH3F*)hPtCascadePlusTrue->Clone("PlusTrueClone");
  hPtCascadeCentrPlusTrue->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrPlusTrue->Write();

  TH2F* hProfilePtCascadeZPlusTrue = static_cast<TH2F*>(hPtCascadeCentrPlusTrue->Project3D("xy"));
  hProfilePtCascadeZPlusTrue->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZPlusTrue->Write();

  TH1D* hProfilePtCascadeYZPlusTrue = (hProfilePtCascadeZPlusTrue->ProjectionY());
  hProfilePtCascadeYZPlusTrue->Sumw2();
  hProfilePtCascadeYZPlusTrue->Write();

  // True lvl Minus Not Assoc
  TH3F* hPtCascadeCentrMinusTrue = (TH3F*)hPtCascadeMinusTrue->Clone("MinusTrueClone");
  hPtCascadeCentrMinusTrue->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrMinusTrue->Write();

  TH2F* hProfilePtCascadeZMinusTrue = static_cast<TH2F*>(hPtCascadeCentrMinusTrue->Project3D("xy"));
  hProfilePtCascadeZMinusTrue->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZMinusTrue->Write();

  TH1D* hProfilePtCascadeYZMinusTrue = (hProfilePtCascadeZMinusTrue->ProjectionY());
  hProfilePtCascadeYZMinusTrue->Sumw2();
  hProfilePtCascadeYZMinusTrue->Write();

  // True lvl Sum Assoc
  TH3F* hPtCascadeCentrSum = (TH3F*)hPtCascadeSumAssoc->Clone("SumAssocClone");
  hPtCascadeCentrSum->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrSum->Write();

  TH2F* hProfilePtCascadeZSum = static_cast<TH2F*>(hPtCascadeCentrSum->Project3D("xy"));
  hProfilePtCascadeZSum->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZSum->Write();

  TH1D* hProfilePtCascadeYZSum = (hProfilePtCascadeZSum->ProjectionY());
  hProfilePtCascadeYZSum->Sumw2();
  hProfilePtCascadeYZSum->Write();

  // True lvl Plus Assoc
  TH3F* hPtCascadeCentrPlus = (TH3F*)hPtCascadePlusTrueAssoc->Clone("PlusTrueAssocClone");
  hPtCascadeCentrPlus->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrPlus->Write();

  TH2F* hProfilePtCascadeZPlus = static_cast<TH2F*>(hPtCascadeCentrPlus->Project3D("xy"));
  hProfilePtCascadeZPlus->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZPlus->Write();

  TH1D* hProfilePtCascadeYZPlus = (hProfilePtCascadeZPlus->ProjectionY());
  hProfilePtCascadeYZPlus->Sumw2();
  hProfilePtCascadeYZPlus->Write();

  // True lvl Minus Assoc
  TH3F* hPtCascadeCentrMinus = (TH3F*)hPtCascadeMinusTrueAssoc->Clone("MinusTrueAssocClone");
  hPtCascadeCentrMinus->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrMinus->Write();

  TH2F* hProfilePtCascadeZMinus = static_cast<TH2F*>(hPtCascadeCentrMinus->Project3D("xy"));
  hProfilePtCascadeZMinus->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZMinus->Write();

  TH1D* hProfilePtCascadeYZMinus = (hProfilePtCascadeZMinus->ProjectionY());
  hProfilePtCascadeYZMinus->Sumw2();
  hProfilePtCascadeYZMinus->Write();

  // Rec lvl Sum
  TH3F* hPtCascadeCentrRecSum = (TH3F*)hPtCascadeRecSum->Clone("SumRecClone");
  hPtCascadeCentrRecSum->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrRecSum->Write();

  TH2F* hProfilePtCascadeZRecSum = static_cast<TH2F*>(hPtCascadeCentrRecSum->Project3D("xy"));
  hProfilePtCascadeZRecSum->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZRecSum->Write();

  TH1D* hProfilePtCascadeYZRecSum = (hProfilePtCascadeZRecSum->ProjectionY());
  hProfilePtCascadeYZRecSum->Sumw2();
  hProfilePtCascadeYZRecSum->Write();

  // Rec lvl Plus
  TH3F* hPtCascadeCentrRecPlus = (TH3F*)hPtCascadePlusTrueRec->Clone("PlusTrueRecClone");
  hPtCascadeCentrRecPlus->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrRecPlus->Write();

  TH2F* hProfilePtCascadeZRecPlus = static_cast<TH2F*>(hPtCascadeCentrRecPlus->Project3D("xy"));
  hProfilePtCascadeZRecPlus->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZRecPlus->Write();

  TH1D* hProfilePtCascadeYZRecPlus = (hProfilePtCascadeZRecPlus->ProjectionY());
  hProfilePtCascadeYZRecPlus->Sumw2();
  hProfilePtCascadeYZRecPlus->Write();

  // Rec lvl Minus
  TH3F* hPtCascadeCentrRecMinus = (TH3F*)hPtCascadeMinusTrueRec->Clone("MinusTrueRecClone");
  hPtCascadeCentrRecMinus->GetZaxis()->SetRangeUser(leftCentrFT0, rightCentrFT0); // MB eff
  hPtCascadeCentrRecMinus->Write();

  TH2F* hProfilePtCascadeZRecMinus = static_cast<TH2F*>(hPtCascadeCentrRecMinus->Project3D("xy"));
  hProfilePtCascadeZRecMinus->GetXaxis()->SetRangeUser(leftRap, rightRap); // Set rapidity interval
  hProfilePtCascadeZRecMinus->Write();

  TH1D* hProfilePtCascadeYZRecMinus = (hProfilePtCascadeZRecMinus->ProjectionY());
  hProfilePtCascadeYZRecMinus->Sumw2();
  hProfilePtCascadeYZRecMinus->Write();

  if(is13tevBinning){
    hProfilePtCascadeYZRecPlus = (TH1D*)hProfilePtCascadeYZRecPlus->Rebin(nBinsChiara, "hProfilePtCascadeYZRecPlus", xbinsChiara);
    hProfilePtCascadeYZPlus = (TH1D*)hProfilePtCascadeYZPlus->Rebin(nBinsChiara, "hProfilePtCascadeYZPlus", xbinsChiara);
    hProfilePtCascadeYZPlusTrue = (TH1D*)hProfilePtCascadeYZPlusTrue->Rebin(nBinsChiara, "hProfilePtCascadeYZPlusTrue", xbinsChiara);

    hProfilePtCascadeYZRecMinus = (TH1D*)hProfilePtCascadeYZRecMinus->Rebin(nBinsChiara, "hProfilePtCascadeYZRecMinus", xbinsChiara);
    hProfilePtCascadeYZMinus = (TH1D*)hProfilePtCascadeYZMinus->Rebin(nBinsChiara, "hProfilePtCascadeYZMinus", xbinsChiara);
    hProfilePtCascadeYZMinusTrue = (TH1D*)hProfilePtCascadeYZMinusTrue->Rebin(nBinsChiara, "hProfilePtCascadeYZMinusTrue", xbinsChiara);

    hProfilePtCascadeYZRecSum = (TH1D*)hProfilePtCascadeYZRecSum->Rebin(nBinsChiara, "hProfilePtCascadeYZRecSum", xbinsChiara);
    hProfilePtCascadeYZSum = (TH1D*)hProfilePtCascadeYZSum->Rebin(nBinsChiara, "hProfilePtCascadeYZSum", xbinsChiara);
    hProfilePtCascadeYZSumTrue = (TH1D*)hProfilePtCascadeYZSumTrue->Rebin(nBinsChiara, "hProfilePtCascadeYZSumTrue", xbinsChiara);
  } else {
    hProfilePtCascadeYZRecPlus = (TH1D*)hProfilePtCascadeYZRecPlus->Rebin(nBinsUpasana, "hProfilePtCascadeYZRecPlus", xbinsUpasana);
    hProfilePtCascadeYZPlus = (TH1D*)hProfilePtCascadeYZPlus->Rebin(nBinsUpasana, "hProfilePtCascadeYZPlus", xbinsUpasana);
    hProfilePtCascadeYZPlusTrue = (TH1D*)hProfilePtCascadeYZPlusTrue->Rebin(nBinsUpasana, "hProfilePtCascadeYZPlusTrue", xbinsUpasana);

    hProfilePtCascadeYZRecMinus = (TH1D*)hProfilePtCascadeYZRecMinus->Rebin(nBinsUpasana, "hProfilePtCascadeYZRecMinus", xbinsUpasana);
    hProfilePtCascadeYZMinus = (TH1D*)hProfilePtCascadeYZMinus->Rebin(nBinsUpasana, "hProfilePtCascadeYZMinus", xbinsUpasana);
    hProfilePtCascadeYZMinusTrue = (TH1D*)hProfilePtCascadeYZMinusTrue->Rebin(nBinsUpasana, "hProfilePtCascadeYZMinusTrue", xbinsUpasana);

    hProfilePtCascadeYZRecSum = (TH1D*)hProfilePtCascadeYZRecSum->Rebin(nBinsUpasana, "hProfilePtCascadeYZRecSum", xbinsUpasana);
    hProfilePtCascadeYZSum = (TH1D*)hProfilePtCascadeYZSum->Rebin(nBinsUpasana, "hProfilePtCascadeYZSum", xbinsUpasana);
    hProfilePtCascadeYZSumTrue = (TH1D*)hProfilePtCascadeYZSumTrue->Rebin(nBinsUpasana, "hProfilePtCascadeYZSumTrue", xbinsUpasana);
  }

  TDirectory* dirEffAcc = fileOut->mkdir("effAcc");
  dirEffAcc->cd();

  // Efficiency x Acceptance
  // Eff Plus
  TH1D* hEffCascadePlus = (TH1D*)hProfilePtCascadeYZRecPlus->Clone("EffPlus");
  hEffCascadePlus->SetTitle("EfficiencyForCascadePlotPlus");
  hEffCascadePlus->SetName("hEffCascPlus");
  hEffCascadePlus->GetYaxis()->SetTitle("Efficiency");
  hEffCascadePlus->Divide(hProfilePtCascadeYZPlus);
  SetErrorDivide(hEffCascadePlus, hProfilePtCascadeYZRecPlus, hProfilePtCascadeYZPlus);
  hEffCascadePlus->SetStats(0);
  hEffCascadePlus->Write();

  // Eff Minus
  TH1D* hEffCascadeMinus = (TH1D*)hProfilePtCascadeYZRecMinus->Clone("EffMinus");
  hEffCascadeMinus->SetTitle("EfficiencyForCascadePlotMinus");
  hEffCascadeMinus->SetName("hEffCascMinus");
  hEffCascadeMinus->GetYaxis()->SetTitle("Efficiency");
  hEffCascadeMinus->Divide(hProfilePtCascadeYZMinus);
  SetErrorDivide(hEffCascadeMinus, hProfilePtCascadeYZRecMinus, hProfilePtCascadeYZMinus);
  hEffCascadeMinus->SetStats(0);
  hEffCascadeMinus->Write();

  // Eff Sum
  TH1D* hEffCascadeSum = (TH1D*)hProfilePtCascadeYZRecSum->Clone("EffSum");
  hEffCascadeSum->SetTitle("EfficiencyForCascadePlotSum");
  hEffCascadeSum->SetName("hEffCascSum");
  hEffCascadeSum->GetYaxis()->SetTitle("Efficiency");
  hEffCascadeSum->Divide(hProfilePtCascadeYZSum);
  SetErrorDivide(hEffCascadeSum, hProfilePtCascadeYZRecSum, hProfilePtCascadeYZSum);
  hEffCascadeSum->SetStats(0);
  hEffCascadeSum->Write();

  TDirectory* dirEffSignal = fileOut->mkdir("effSignal");
  dirEffSignal->cd();

  // Signal loss
  // Eff Plus
  TH1D* hEffCascadePlusSignal = (TH1D*)hProfilePtCascadeYZPlus->Clone("SignalLossPlus");
  hEffCascadePlusSignal->SetTitle("EfficiencyForCascadePlotPlus");
  hEffCascadePlusSignal->SetName("hEffCascPlusSignal");
  hEffCascadePlusSignal->GetYaxis()->SetTitle("Signal Loss");
  hEffCascadePlusSignal->Divide(hProfilePtCascadeYZPlusTrue);
  SetErrorDivide(hEffCascadePlusSignal, hProfilePtCascadeYZPlus, hProfilePtCascadeYZPlusTrue);
  hEffCascadePlusSignal->SetStats(0);
  hEffCascadePlusSignal->Write();

  // Eff Minus
  TH1D* hEffCascadeMinusSignal = (TH1D*)hProfilePtCascadeYZMinus->Clone("SignalLossMinus");
  hEffCascadeMinusSignal->SetTitle("EfficiencyForCascadePlotMinus");
  hEffCascadeMinusSignal->SetName("hEffCascMinusSignal");
  hEffCascadeMinusSignal->GetYaxis()->SetTitle("Signal Loss");
  hEffCascadeMinusSignal->Divide(hProfilePtCascadeYZMinusTrue);
  SetErrorDivide(hEffCascadeMinusSignal, hProfilePtCascadeYZMinus, hProfilePtCascadeYZMinusTrue);
  hEffCascadeMinusSignal->SetStats(0);
  hEffCascadeMinusSignal->Write();

  // Eff Sum
  TH1D* hEffCascadeSumSignal = (TH1D*)hProfilePtCascadeYZSum->Clone("SignalLossSum");
  hEffCascadeSumSignal->SetTitle("EfficiencyForCascadePlotSum");
  hEffCascadeSumSignal->SetName("hEffCascSumSignal");
  hEffCascadeSumSignal->GetYaxis()->SetTitle("Signal Loss");
  hEffCascadeSumSignal->Divide(hProfilePtCascadeYZSumTrue);
  SetErrorDivide(hEffCascadeSumSignal, hProfilePtCascadeYZSum, hProfilePtCascadeYZSumTrue);
  hEffCascadeSumSignal->SetStats(0);
  hEffCascadeSumSignal->Write();

  //PrintBinning(hProfilePtOmegaSumYZRec);

  fileOut->Close();
  gROOT->SetBatch(kFALSE);
}


