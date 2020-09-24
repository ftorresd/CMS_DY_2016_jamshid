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
void HistoSum_Zall(){
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

// S I G N A L    P A R T :
// =======================
// =======================

// Defining Model and its parameters for Crystal Ball Function PDF
// ***************************************************************

RooRealVar mean("mean", "mean of the PDF", 91.019);
RooRealVar cbsigma("cbsigma", "width of the CB PDF", 1.9181);
//RooRealVar cbsig("cbsig", "cbsignal", 1000., 0., 1000000.);
RooRealVar n("n", "n", 0.77582);
RooRealVar alpha("alpha", "alpha", 1.7615);

// Crystal Ball  PDF
// *****************

RooCBShape cb("cb", "cb", di_mu, mean, cbsigma, alpha, n);

// Defining Model and Parameters for Gaussian PDF
// **********************************************

RooRealVar gaussigma("gaussigma", "Width of the Gaussian PDF", 7.00);

// Gaussian PDF
// ************

RooGaussian gauss("gauss", " Gaussiam PDF", di_mu, mean, gaussigma);

//....................................
// Composite Signal Model (CB+Gauss) :
// ...................................

RooRealVar f("f", "signal fraction", 0.7337);
RooAddPdf comp("comp", "cb+gauss", RooArgList(cb, gauss), f);


//================================
//================================
// B A C K G R O U N D   P A R T :
//================================
//================================

// Background Model and parameters (Chebychev polynomial PDF)
//-----------------------------------------------------------

   RooRealVar a0("a0","a0", -1.1519) ;
   RooRealVar a1("a1","a1", 0.24363) ;

// Chebychev PDF
// ************

   RooChebychev bkg("bkg","Background", di_mu ,RooArgSet(a0,a1)) ;


// Associate nsig/nbkg as expected number of events with sig and bkg in the range SignalRange
// -----------------------------------------------------------------------------------------

RooRealVar nsig("nsig","number of signal events in signalRange",0.89, 0., 1.) ;


// *************************************************
// Designing the composite Model (Background+Signal)
// *************************************************

RooAddPdf model("model", "Signal+Background", RooArgList(comp,bkg), nsig);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//   S a m p l e   d a t a ,   F i t   m o d e l  <|>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



TCanvas *c6 = new TCanvas("c6", "c6");


c6->cd();

  

RooFitResult* r1 = model.fitTo(data,Save());
r1->Print();

// Full Plot :-
//________________

RooPlot *frame = di_mu.frame(Title("Signal Yield"));
data.plotOn(frame);
model.plotOn(frame, LineColor(kMagenta));
model.paramOn(frame, Layout(0.6, 0.9, 0.7));
// data.statOn(frame);


// Background Fit :-
// ________________

  

model.plotOn(frame, Components("bkg"),LineColor(2));

model.plotOn(frame, Components("comp"), LineColor(4));


frame->Draw();


c6->Draw(); 
 
c6->SaveAs("outputs/Zfull1_mass.png");
}  
