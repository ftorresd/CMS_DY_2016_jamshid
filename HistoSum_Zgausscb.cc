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
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooChebychev.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooCBShape.h"

using namespace RooFit;

void HistoSum_Zgausscb(){
 gROOT->Reset();

// gStyle->SetOptTitle(kFALSE);
 // gStyle->SetOptStat(000000000);

// Accessing the ROOT file containing the binned data (Histogram)
// *************************************************************

  TFile *f0 = TFile::Open("outputs/new_hist.root");
  gDirectory->ls();
  TH1F* h_z = (TH1F*)f0->Get("t11");

// Redefining the observable with which the histogram "h11"  is filled i.e "dimu":
// *******************************************************************************

  RooRealVar dimu("dimu", "dimu", 60., 120.);

// Assigning histogram dataset to dimu object
// ******************************************

  RooDataHist data("data", "dataset", RooArgSet(dimu), h_z);

// Defining Model and its parameters for Crystal Ball Function PDF
// ***************************************************************

RooRealVar mean("mean", "mean of the PDF", 91., 60., 120.);
RooRealVar cbsigma("cbsigma", "width of the CB PDF", 1.0, 0., 2.6);
RooRealVar cbsig("cbsig", "cbsignal", 1000., 0., 1000000.);
RooRealVar n("n", "n", 2.5, 0., 5.);
RooRealVar alpha("alpha", "alpha", 2.5, 0., 5.);

// Crystal Ball  PDF
//*******************

RooCBShape cb("cb", "cb", dimu, mean, cbsigma, alpha, n);

// Defining Model and Parameters for Gaussian PDF
// **********************************************

// RooRealVar gaussmean("gaussmean", "Mean of the Gaussian PDF", 91., 60., 120.);
RooRealVar gaussigma("gaussigma", "Width of the Gaussian PDF", 2.5, 0., 7.);

// Gaussian PDF
// ************

RooGaussian gauss("gauss", " Gaussiam PDF", dimu, mean, gaussigma);

// Designing the composite Model
// *****************************

RooRealVar f("f", "signal fraction", 0.8, 0., 1.);
RooAddPdf comp("comp", "cb+gauss", RooArgList(cb, gauss), f);


//------------------------------------------------------
// **  S A M P L E   D A T A   ,   F I T    M O D E L  **
//-------------------------------------------------------


 TCanvas *c = new TCanvas("c", "c");
 c->cd();
 c->SetTicks(1, 1);
 RooFitResult *r = comp.fitTo(data, Save());
 r -> Print();

// Plotting the Model
//******************

RooPlot *frame = dimu.frame(Title("CB+Gauss Fit of #mu^{+}#mu^{-} Invariant mass"));
data.plotOn(frame);
comp.plotOn(frame);
comp.paramOn(frame, Layout(0.6, .9, .9));
frame->Draw();
c->Draw();
c->SaveAs("outputs/Zgausscb_mass.png");
}  
