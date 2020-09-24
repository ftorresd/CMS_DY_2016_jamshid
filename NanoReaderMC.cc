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

#include "include/progressbar.hpp"

void NanoReaderMC()
{
    bool _debug = false; // should print stuff?

    gROOT->Reset();

    TFile *file = new TFile("RunIISummer16NanoAODv6_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_NANOAODSIM_PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1.root");

    TTreeReader reader("Events", file);
    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_mass(reader, "Muon_mass");
    TTreeReaderArray<UInt_t> nMuon(reader, "nMuon");
    TTreeReaderArray<Int_t> Muon_charge(reader, "Muon_charge");
    TTreeReaderArray<Bool_t> Muon_softId(reader, "Muon_softId");
    TTreeReaderValue<Bool_t> HLT_IsoMu24(reader, "HLT_IsoMu24");
    TTreeReaderValue<Bool_t> HLT_IsoMu27(reader, "HLT_IsoMu27");

    // define the readers of the LHEPart info
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    TTreeReaderValue<UInt_t> nLHEPart(reader, "nLHEPart");
    TTreeReaderArray<Float_t> LHEPart_pt(reader, "LHEPart_pt");
    TTreeReaderArray<Float_t> LHEPart_eta(reader, "LHEPart_eta");
    TTreeReaderArray<Float_t> LHEPart_phi(reader, "LHEPart_phi");
    TTreeReaderArray<Int_t> LHEPart_pdgId(reader, "LHEPart_pdgId");

    // define the total number of events counter
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int nEvents = reader.GetEntries(0);
    Double_t ptCut = 5;
    Int_t nMuonCut = 2;
    Int_t count = 0;

    // Defining Histograms
    // ^^^^^^^^^^^^^^^^^^^

    TH1F *t1 = new TH1F("t1", "Muon0_pt", 50, 0, 150);
    TH1F *t2 = new TH1F("t2", "Muon0_eta", 50, -10, 10);
    TH1F *t3 = new TH1F("t3", "Muon0_phi", 51, -6, 6);
    TH1F *t4 = new TH1F("t4", "Muon0_charge", 50, -10, 10);
    TH1F *t5 = new TH1F("t5", "Muon1_pt", 50, 0, 60);
    TH1F *t6 = new TH1F("t6", "Muon1_eta", 50, -10, 10);
    TH1F *t7 = new TH1F("t7", "Muon1_phi", 51, -6, 6);
    TH1F *t8 = new TH1F("t8", "Muon1_charge", 50, -10, 10);
    TH1F *t9 = new TH1F("t9", "Muon0_mass", 20, 0, 1);
    TH1F *t10 = new TH1F("t10", "Muon1_mass", 20, 0, 1);
    TH1F *t11 = new TH1F("t11", "di-#mu Invariant Mass", 50, 60, 120);

    int maxEvents = -1;
    if (maxEvents == -1)
    {
        maxEvents = nEvents;
    }
    std::cout << "--> Max Events to be processed: " << maxEvents << "\n\n" << std::endl;

    int eventsReadCounter = 0;
    int goodLHEMassEvents = 0;
    std::cout << "Starting event loop..." << std::endl;
    progressbar bar(maxEvents);
    while (reader.Next() && eventsReadCounter < maxEvents)
    {
        bar.update();
        eventsReadCounter++;

        if (_debug)
            std::cout << "==================== Start GEN level filtering ==========================" << std::endl;

        // At least two particles
        // ^^^^^^^^^^^^^^^^^^^^^^

        if (*nLHEPart < 2)
            continue;

        // Loop over LHE particles and select the first two muons
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        TLorentzVector LHEMuon_plus;
        TLorentzVector LHEMuon_minus;
        Int_t nLHEMuons = 0;
        for (UInt_t i = 0; i < *nLHEPart; i++)
        {

            // Check for positive muons
            // ^^^^^^^^^^^^^^^^^^^^^^^^

            if (LHEPart_pdgId[i] == 13)
            {
                LHEMuon_plus.SetPtEtaPhiM(LHEPart_pt[i], LHEPart_eta[i], LHEPart_phi[i], 105.6583755 / 1000.);
                nLHEMuons++;
            }

            // Check for negative muons
            // ^^^^^^^^^^^^^^^^^^^^^^^^

            if (LHEPart_pdgId[i] == -13)
            {
                LHEMuon_minus.SetPtEtaPhiM(LHEPart_pt[i], LHEPart_eta[i], LHEPart_phi[i], 105.6583755 / 1000.);
                nLHEMuons++;
            }
        }

        // At least two muons
        // ^^^^^^^^^^^^^^^^^^

        if (nLHEMuons < 2)
            continue;

        // Build dimuon pair
        // ^^^^^^^^^^^^^^^^^

        float LHEdimuon_mass = (LHEMuon_plus + LHEMuon_minus).Mag();

        // If the events has gone so far and the mass is within the range (60 to 120), the event is analyzed
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        if (LHEdimuon_mass < 60. || LHEdimuon_mass > 120.)
        {
            continue;
        }

        if (_debug)
            std::cout << "==================== End GEN level filtering ==========================" << std::endl;

        goodLHEMassEvents++;

        int n = nMuon[0];

        Float_t ptMax = Muon_pt[0];
        if (ptMax < ptCut || n < nMuonCut)
        {
            continue;
        }
        // count++;
        if (_debug)
            std::cout << "==============================================" << std::endl;
        if (_debug)
            std::cout << "Event Number: " << reader.GetCurrentEntry() << std::endl;
        if (_debug)
            std::cout << "Number of Muons: " << n << std::endl;
        TLorentzVector pMuon[n];
        // for (Int_t i = 0; i < n; i++)
        for (Int_t i = 0; i < 2; i++) // you only need to loop over the first two muons
        {
            // pMuon[i].SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
            pMuon[i].SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 105.6583755 / 1000.); // just to be sure, lets give that the PDG muon mass. Most probably maskes no difference.

            /*std::cout << Muon_charge[i]     << "       ";
            std::cout << Muon_pt[i]         << "       ";
            std::cout << pMuon[i].E()       << "       ";
            std::cout << pMuon[i].Px()      << "       ";
            std::cout << pMuon[i].Py()      << "       ";
            std::cout << pMuon[i].Pz()      << std::endl; */

            //     A P P L Y I N G    T H E    C U T S :
            //     *************************************
        }
        Double_t dimu = (pMuon[0] + pMuon[1]).Mag();
        Double_t cut0 = 60;
        Double_t cut1 = 120;
        Double_t zchar = Muon_charge[0] + Muon_charge[1];
        Double_t cut2 = 0;
        Double_t high_pt0 = Muon_pt[0];
        Double_t high_pt1 = Muon_pt[1];
        Double_t cut4 = 26;
        Double_t cut5 = 3;
        Bool_t softId0 = Muon_softId[0];
        Bool_t softId1 = Muon_softId[1];
        if (dimu < cut0 || dimu > cut1 || zchar != cut2 || high_pt0 < cut4 || high_pt1 < cut5 || *HLT_IsoMu24 == false || softId0 == false || softId1 == false)
        {
            continue;
        }
        count++;

        if (_debug)
            std::cout << " Muon[0] pT = " << Muon_pt[0] << " Muon[0] eta =  " << Muon_eta[0] << " Muon[0] phi = " << Muon_phi[0] << " Muon[0] charge = " << Muon_charge[0] << std::endl;
        if (_debug)
            std::cout << " Muon[1] pT = " << Muon_pt[1] << " Muon[1]  eta = " << Muon_eta[1] << " Muon[1] phi = " << Muon_phi[1] << " Muon[1] charge = " << Muon_charge[1] << std::endl;
        if (_debug)
            std::cout << " Largest pT of Muon = " << high_pt0 << std::endl;
        if (_debug)
            std::cout << " Di_Muon Invariant Mass from Z decay = " << dimu << std::endl;
        if (_debug)
            std::cout << " Electric charge of Z = " << zchar << std::endl;
        if (_debug)
            std::cout << " Pass HLT_IsoMu24 ? = " << *HLT_IsoMu24 << std::endl;
        if (_debug)
            std::cout << " Muon[0]  ID = " << softId0 << std::endl;
        if (_debug)
            std::cout << " Muon[1]  ID = " << softId1 << std::endl;

        //  std::cout<< " Largest pT of Muon = "<< high_pt <<std::endl;

        // Filling Histograms
        // ^^^^^^^^^^^^^^^^^^

        t1->Fill(Muon_pt[0]);
        t2->Fill(Muon_eta[0]);
        t3->Fill(Muon_phi[0]);
        t4->Fill(Muon_charge[0]);
        t5->Fill(Muon_pt[1]);
        t6->Fill(Muon_eta[1]);
        t7->Fill(Muon_phi[1]);
        t8->Fill(Muon_charge[1]);
        t9->Fill(pMuon[0].Mag());
        t10->Fill(pMuon[1].Mag());
        t11->Fill(dimu);
    }

    // float eff = float(count) / float(nEvents);
    float eff = float(count) / float(goodLHEMassEvents);
    std::cout << "\n" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Number of events on the dataset: " << nEvents << std::endl;
    std::cout << "Total processed events: " << eventsReadCounter << std::endl;
    std::cout << "Total processed events that pass the LHE mass filter: " << goodLHEMassEvents << std::endl;
    // std::cout << "Number of events that passed the cut Pt > " << ptCut
    //           << " GeV and number of muons > " << nMuonCut << ": " << count << std::endl;
    std::cout << "Number of events that passed full selection: " << count << std::endl;
    std::cout << "Selection efficiency: " << eff << std::endl;
    std::cout << "==============================================" << std::endl;


    TFile *hist_sv = new TFile("outputs/new_hist.root", "RECREATE");

    // Defining Canvases
    // ^^^^^^^^^^^^^^^^^

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();
    t1->Draw();
    t1->Write();
    c1->SaveAs("outputs/new_mu_pt0mc.root");
    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->cd();
    t2->Draw();
    t2->Write();
    c2->SaveAs("outputs/new_mu_eta0mc.root");

    TCanvas *c3 = new TCanvas("c3", "c3");
    c3->cd();
    t3->Draw();
    t3->Write();
    c3->SaveAs("outputs/new_mu_phi0mc.root");

    TCanvas *c4 = new TCanvas("c4", "c4");
    c4->cd();
    t4->Draw();
    t4->Write();

    c4->SaveAs("outputs/new_mu_charge0mc.root");

    TCanvas *c5 = new TCanvas("c5", "c5");
    c5->cd();
    t5->Draw();
    t5->Write();
    c5->SaveAs("outputs/new_mu_pt1mc.root");

    TCanvas *c6 = new TCanvas("c6", "c6");
    c6->cd();
    t6->Draw();
    t6->Write();
    c6->SaveAs("outputs/new_mu_eta1mc.root");

    TCanvas *c7 = new TCanvas("c7", "c7");
    c7->cd();
    t7->Draw();
    t7->Write();
    c7->SaveAs("outputs/new_mu_phi1mc.root");

    TCanvas *c8 = new TCanvas("c8", "c8");
    c8->cd();
    t8->Draw();
    t8->Write();
    c8->SaveAs("outputs/new_mu_charge1mc.root");

    TCanvas *c9 = new TCanvas("c9", "c9");
    c9->cd();
    t9->Draw();
    t9->Write();
    c9->SaveAs("outputs/new_mu_0_massmc.root");

    TCanvas *c10 = new TCanvas("c10", "c10");
    c10->cd();
    t10->Draw();
    t10->Write();
    c10->SaveAs("outputs/new_mu_1_massmc.root");

    TCanvas *c11 = new TCanvas("c11", "di-muon-mass");
    c11->cd();
    t11->Draw();
    t11->Write();
    c11->SaveAs("outputs/new_di_mu-massmc.root");

    hist_sv->Write();
    hist_sv->Close();
}
