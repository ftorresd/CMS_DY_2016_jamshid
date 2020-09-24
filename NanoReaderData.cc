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
#include <string>


#include "include/json.hpp"
#include "include/progressbar.hpp"


using nlohmann::json;

void NanoReaderData()
{
    bool _debug = false; // should print stuff?

    gROOT->Reset();

    // TFile *file = new TFile("SingleMuon_Run2016G-Nano25Oct2019-v1_NANOAOD.root");
    TFile *file = new TFile("SingleMuon_Run2016H-Nano25Oct2019-v1_NANOAOD.root");
    TTreeReader reader("Events", file);
    TTreeReaderArray<Float_t> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<Float_t> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<Float_t> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<Float_t> Muon_mass(reader, "Muon_mass");
    TTreeReaderArray<UInt_t> nMuon(reader, "nMuon");
    TTreeReaderArray<Int_t> Muon_charge(reader, "Muon_charge");
    TTreeReaderArray<Bool_t> Muon_softId(reader, "Muon_softId");
    TTreeReaderValue<Bool_t> HLT_IsoMu24(reader, "HLT_IsoMu24");
    TTreeReaderValue<UInt_t> run(reader, "run");
    TTreeReaderValue<UInt_t> luminosityBlock(reader, "luminosityBlock");

    int nEvents = reader.GetEntries(0);
    Double_t ptCut = 5;
    Int_t nMuonCut = 2;
    Int_t count = 0;

    // Defining Histograms
    TH1F *k1 = new TH1F("k1", "Muon0_pt", 50, 0, 150);
    TH1F *k2 = new TH1F("k2", "Muon0_eta", 50, -10, 10);
    TH1F *k3 = new TH1F("k3", "Muon0_phi", 51, -6, 6);
    TH1F *k4 = new TH1F("k4", "Muon0_charge", 50, -10, 10);
    TH1F *k5 = new TH1F("k5", "Muon1_pt", 50, 0, 60);
    TH1F *k6 = new TH1F("k6", "Muon1_eta", 50, -10, 10);
    TH1F *k7 = new TH1F("k7", "Muon1_phi", 51, -6, 6);
    TH1F *k8 = new TH1F("k8", "Muon1_charge", 50, -10, 10);
    TH1F *k9 = new TH1F("k9", "Muon0_mass", 20, 0, 1);
    TH1F *k10 = new TH1F("k10", "Muon1_mass", 20, 0, 1);
    TH1F *k11 = new TH1F("k11", "di-#mu Invariant Mass", 50, 60, 120);

    // std::ifstream ifs("jamshid.json");
    std::ifstream ifs("jamshid_Run2016H.json");
    json good_runs_lumis = json::parse(ifs);

    int maxEvents = -1;
    if (maxEvents == -1) {
        maxEvents = nEvents;
    }
    std::cout << "--> Max Events to be processed: " << maxEvents << "\n\n" << std::endl;

    int eventsReadCounter = 0;
    std::cout << "Starting event loop..." << std::endl;
    progressbar bar(maxEvents);
    while (reader.Next() && eventsReadCounter < maxEvents)
    {
        bar.update();
        eventsReadCounter++;

        bool is_good_run_lumi = false;
        // bool is_good_run_lumi = true;
        string test_run = std::to_string(*run);
        int test_lumi = *luminosityBlock;
        if (good_runs_lumis.find(test_run) != good_runs_lumis.end())
        {

            for (auto const &interval : good_runs_lumis[test_run])
            {

                int low = interval.front();
                int high = interval.back();
                if (test_lumi >= low && test_lumi <= high)
                {
                    is_good_run_lumi = true;
                }
            }
        }
        if (is_good_run_lumi == false)
            continue;

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
            /* std::cout << Muon_charge[i]     << "       ";
            std::cout << Muon_pt[i]         << "       ";
            std::cout << pMuon[i].E()       << "       ";
            std::cout << pMuon[i].Px()      << "       ";
            std::cout << pMuon[i].Py()      << "       ";
            std::cout << pMuon[i].Pz()      << std::endl; */
        }
        // Accessing to pMuon[0] and pMuon[1] since they are the highest pt Muons in the event
        if (_debug)
            std::cout << " Muon 0 PT = " << Muon_pt[0] << " Muon 0 eta =  " << Muon_eta[0] << " Muon 0 Phi = " << Muon_phi[0] << " Muon 0 Charge = " << Muon_charge[0] << std::endl;
        if (_debug)
            std::cout << " Muon 1 PT = " << Muon_pt[1] << " Muon 1 eta = " << Muon_eta[1] << " Muon 1 Phi = " << Muon_phi[1] << " Muon 1 Charge = " << Muon_charge[1] << std::endl;
        // Invariant Masses of Muon[0] & Muon[1]
        /*
      std::cout<<" Muon 0 Px = "<< pMuon[0].Px()<<" Muon 0 Py = "<<pMuon[0].Py()<<" Muon 0 Pz = "<<pMuon[0].Pz()<<" Muon 0 E = "          <<pMuon[0].E()<<std::endl;
      std::cout<<" Muon 0 invariant Mass = "<<pMuon[0].Mag()<<std::endl;
      std::cout<<" Muon 1 Px = "<< pMuon[1].Px()<<" Muon 1 Py = "<<pMuon[1].Py()<<" Muon 1 Pz = "<<pMuon[1].Pz()<<" Muon 1 E = "
      <<pMuon[1].E()<<std::endl;
      std::cout<<" Muon 1 invariant Mass = "<<pMuon[1].Mag()<<std::endl;
*/
        // Invariant mass of Di-Muon
        //  std::cout<<" Di-Muon Invariant Mass = "<<(pMuon[0]+pMuon[1]).Mag()<<std::endl;

        Double_t di_mu = (pMuon[0] + pMuon[1]).Mag();
        Double_t cut2 = 60;
        Double_t cut3 = 120;
        Double_t zchar = Muon_charge[0] + Muon_charge[1];
        Double_t cut4 = 0;
        Double_t cut5 = 26;
        Double_t cut6 = 3;
        Double_t high_pt0 = Muon_pt[0];
        Double_t high_pt1 = Muon_pt[1];
        Double_t eta0 = Muon_eta[0];
        Double_t eta1 = Muon_eta[1];
        Double_t cut7 = 2.4;
        Double_t cut8 = -2.4;
        Bool_t softId0 = Muon_softId[0];
        Bool_t softId1 = Muon_softId[1];


        if (di_mu < cut2 || di_mu > cut3 || zchar != cut4 || high_pt0 < cut5 || eta0 > cut7 || eta1 > cut7 || eta0 < cut8 || eta1 < cut8 || high_pt1 < cut6 || *HLT_IsoMu24 == false || softId0 == false || softId1 == false)
        {
            continue;
        }
        count++;

        if (_debug)
        {
            std::cout << " Di_Muon Invariant mass from Z decay = " << di_mu << std::endl;
            std::cout << " Electric charge of Z = " << zchar << std::endl;
            std::cout << " Largest pT of Muon = " << high_pt0 << std::endl;
            std::cout << " Pass HLT_IsoMu24 ? = " << *HLT_IsoMu24 << std::endl;
            std::cout << " Muon[0] eta = " << eta0 << std::endl;
            std::cout << " Muon[1] eta = " << eta1 << std::endl;
        }

        // Filling Histograms
        k1->Fill(Muon_pt[0]);
        k2->Fill(Muon_eta[0]);
        k3->Fill(Muon_phi[0]);
        k4->Fill(Muon_charge[0]);
        k5->Fill(Muon_pt[1]);
        k6->Fill(Muon_eta[1]);
        k7->Fill(Muon_phi[1]);
        k8->Fill(Muon_charge[1]);
        k9->Fill(pMuon[0].Mag());
        k10->Fill(pMuon[1].Mag());
        //  k11->Fill((pMuon[0]+pMuon[1]).Mag());
        k11->Fill(di_mu);
    }

    // float eff = float(count) / float(nEvents);
    float eff = float(count) / float(eventsReadCounter);
    std::cout << "\n" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Number of events on the dataset: " << nEvents << std::endl;
    std::cout << "Total processed events: " << eventsReadCounter << std::endl;
    // std::cout << "Number of events that passed the cut Pt > " << ptCut
    //           << " GeV and number of muons > " << nMuonCut << ": " << count << std::endl;
    std::cout << "Number of events that passed full selection: " << count << std::endl;
    std::cout << "Selection efficiency (useless): " << eff << std::endl;
    std::cout << "==============================================" << std::endl;

    TFile *hist_data = new TFile("outputs/hist_data.root", "RECREATE");
    // Defining Canvases
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();
    k1->Draw();
    k1->Write();

    c1->SaveAs("outputs/mu_pt0.root");

    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->cd();
    k2->Draw();
    k2->Write();

    c2->SaveAs("outputs/mu_eta0.root");

    TCanvas *c3 = new TCanvas("c3", "c3");
    c3->cd();
    k3->Draw();
    k3->Write();

    c3->SaveAs("outputs/mu_phi0.root");

    TCanvas *c4 = new TCanvas("c4", "c4");
    c4->cd();
    k4->Draw();
    k4->Write();

    c4->SaveAs("outputs/mu_charge0.root");

    TCanvas *c5 = new TCanvas("c5", "c5");
    c5->cd();
    k5->Draw();
    k5->Write();

    c5->SaveAs("outputs/mu_pt1.root");

    TCanvas *c6 = new TCanvas("c6", "c6");
    c6->cd();
    k6->Draw();
    k6->Write();

    c6->SaveAs("outputs/mu_eta1.root");

    TCanvas *c7 = new TCanvas("c7", "c7");
    c7->cd();
    k7->Draw();
    k7->Write();

    c7->SaveAs("outputs/mu_phi1.root");

    TCanvas *c8 = new TCanvas("c8", "c8");
    c8->cd();
    k8->Draw();
    k8->Write();

    c8->SaveAs("outputs/mu_charge1.root");

    TCanvas *c9 = new TCanvas("c9", "c9");
    c9->cd();
    k9->Draw();
    k9->Write();

    c9->SaveAs("outputs/mu_0_mass.root");

    TCanvas *c10 = new TCanvas("c10", "c10");
    c10->cd();
    k10->Draw();
    k10->Write();

    c10->SaveAs("outputs/mu_1_mass.root");

    TCanvas *c11 = new TCanvas("c11", "di-muon-mass");
    c11->cd();
    k11->Draw();
    k11->Write();

    c11->SaveAs("outputs/di_mu-mass.root");

    hist_data->Write();
    hist_data->Close();
}
