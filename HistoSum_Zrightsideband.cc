#include "ostream"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TMath.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "Fit/FitResult.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooChebychev.h"
#include "RooPlot.h"
#include "RooDataHist.h"

using namespace RooFit;
void HistoSum_Zrightsideband(){
 gROOT->Reset();

//  gStyle->SetOptTitle(kFALSE);
//  gStyle->SetOptStat(000000000);

 TFile *f1 = TFile::Open("outputs/hist_data.root");

 gDirectory->ls();
 TH1F* d_z = (TH1F*)f1->Get("k11");
 
 RooRealVar di_mu("di_mu", "di_mu",  60,120);

// Assigning Histogram dataset to di_mu object
//---------------------------------------------

RooDataHist data("data", "dataset",RooArgSet(di_mu) ,d_z);

// Background Model and parameters (Chebychev polynomial PDF)
//-----------------------------------------------------------

   RooRealVar a0("a0","a0",0.5,-10.,10.) ;
   RooRealVar a1("a1","a1",0.6,-10.,10.) ;
   RooChebychev bkg("bkg","Background", di_mu ,RooArgSet(a0,a1)) ;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   S a m p l e   d a t a ,   F i t   m o d e l  <|>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



TCanvas *c6 = new TCanvas("c6", "c6");

   // Fit in the left and right sideband regions
  // -------------------------------------------

c6->cd();
//   di_mu.setRange("left", 60., 80.);
   di_mu.setRange("right", 103., 120.);
   
   RooFitResult* r2 = bkg.fitTo(data,Range("right"),Save()) ;
   r2->Print();
   
// Plotting the left and right sidebands
// -------------------------------------
   
   RooPlot * frame2 = di_mu.frame(Title("Fit in right sideband"));
   data.plotOn(frame2);
   bkg.plotOn(frame2, VisualizeError(*r2));
   bkg.plotOn(frame2);
   bkg.paramOn(frame2, Layout(0.6,0.9,0.9));
   data.statOn(frame2, Layout(0.6,0.9,0.7));
   frame2->Draw();
   
gPad->SetGrid();
c6->Draw(); 
 
c6->SaveAs("outputs/RightSideband.png");
}  
