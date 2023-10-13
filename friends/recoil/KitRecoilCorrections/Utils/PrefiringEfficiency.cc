#ifndef Prefiring_Efficiency
#define Prefiring_Efficiency

#include <TBenchmark.h> // class to track macro running statistics
#include <TFile.h> // file handle class
#include <TH1D.h> // histogram class
#include <TROOT.h> // access to gROOT, entry point to ROOT system
#include <TStyle.h> // class to handle ROOT plotting styles
#include <TSystem.h> // interface to OS
#include <TTree.h> // class to access ntuples
#include <fstream>
#include <fstream> // functions for file I/O
#include <iomanip> // functions to format standard I/O
#include <iostream> // standard I/O
#include <sstream>
#include <sstream> // class for parsing strings
#include <stdexcept>
#include <string> // C++ string class
#include <vector> // STL vector class

#include <TF1.h>
#include <TH2F.h>

// The bacon particle interfaces
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"

// need the loose muon id
#include "../Utils/LeptonIDCuts.hh"

//
// the muon prefiring calculation is take from
//  https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/plugins/L1PrefiringWeightProducer.cc
//  and the twiki page https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
//

class PrefiringEfficiency {
public:
    PrefiringEfficiency() {}
    ~PrefiringEfficiency() {}
    PrefiringEfficiency(TString fname, TString era, TString fmuonname)
    {
        this->fname = fname;
        this->era = era;
        f = new TFile(fname);
        assert(f);
        loadJetMap();
        loadPhotonMap();
        fmuon = new TFile(fmuonname);
        assert(fmuon);
        loadMuonMap();
    }

    void loadJetMap()
    {
        cout << "Loading Prefiring Jet map for " << era.Data() << endl;
        h_prefmap_jet_ = (TH2F*)f->Get("L1prefiring_jetpt_" + era);
        return;
    }
    void loadPhotonMap()
    {
        cout << "Loading Prefiring Photon map for " << era.Data() << endl;
        h_prefmap_photon_ = (TH2F*)f->Get("L1prefiring_photonpt_" + era);
        return;
    }
    void loadMuonMap()
    {
        TString dataeraMuon = "20172018";
        TString paramName = "L1prefiring_muonparam_0.0To0.2_" + dataeraMuon;
        parametrization0p0To0p2_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_0.2To0.3_" + dataeraMuon;
        parametrization0p2To0p3_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_0.3To0.55_" + dataeraMuon;
        parametrization0p3To0p55_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_0.55To0.83_" + dataeraMuon;
        parametrization0p55To0p83_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_0.83To1.24_" + dataeraMuon;
        parametrization0p83To1p24_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_1.24To1.4_" + dataeraMuon;
        parametrization1p24To1p4_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_1.4To1.6_" + dataeraMuon;
        parametrization1p4To1p6_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_1.6To1.8_" + dataeraMuon;
        parametrization1p6To1p8_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_1.8To2.1_" + dataeraMuon;
        parametrization1p8To2p1_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_2.1To2.25_" + dataeraMuon;
        parametrization2p1To2p25_ = (TF1*)fmuon->Get(paramName);
        paramName = "L1prefiring_muonparam_2.25To2.4_" + dataeraMuon;
        parametrization2p25To2p4_ = (TF1*)fmuon->Get(paramName);

        if (parametrization0p0To0p2_ == nullptr || parametrization0p2To0p3_ == nullptr || parametrization0p3To0p55_ == nullptr || parametrization0p55To0p83_ == nullptr || parametrization0p83To1p24_ == nullptr || parametrization1p24To1p4_ == nullptr || parametrization1p4To1p6_ == nullptr || parametrization1p6To1p8_ == nullptr || parametrization1p8To2p1_ == nullptr || parametrization2p1To2p25_ == nullptr || parametrization2p25To2p4_ == nullptr) {
            assert(0);
        }
        return;
    }

    void setObjects(TClonesArray* photons, TClonesArray* jets, TClonesArray* muons)
    {
        scArr = photons;
        jetArr = jets;
        muonArr = muons;
        return;
    }

    void computeJetsOnly(float& main, float& up, float& down)
    {
        // loop through jets
        float uncTot = 0;
        for (Int_t ij = 0; ij < jetArr->GetEntriesFast(); ij++) {
            const baconhep::TJet* jet = (baconhep::TJet*)((*jetArr)[ij]);
            if (jet->pt < 20.0)
                continue;
            if (!etaCut(jet->eta))
                continue;
            if (!((1 - jet->chEmFrac - jet->neuEmFrac - jet->chHadFrac - jet->neuHadFrac) < 0.50))
                continue;

            double prate_central = 0., prate_up = 0., prate_down = 0.;
            getPrefiringRateEcal(jet->eta, jet->pt, h_prefmap_jet_, prate_central, prate_up, prate_down);
            main *= 1 - prate_central;
            up *= 1 - prate_up;
            down *= 1 - prate_down;
        }
        return;
    }

    void computePhotonsOnly(float& main, float& up, float& down)
    {
        float uncTot = 0;
        for (Int_t ip = 0; ip < scArr->GetEntriesFast(); ip++) {
            const baconhep::TPhoton* photon = (baconhep::TPhoton*)((*scArr)[ip]);
            if (photon->pt < 20.0)
                continue;
            if (!etaCut(photon->eta))
                continue;
            double prate_central = 0., prate_up = 0., prate_down = 0.;
            getPrefiringRateEcal(photon->eta, photon->pt, h_prefmap_photon_, prate_central, prate_up, prate_down);
            main *= 1 - prate_central;
            up *= 1 - prate_up;
            down *= 1 - prate_down;
        }
        return;
    }

    void computeEcalsOnly(float& main, float& up, float& down)
    {
        computePhotonsOnly(main, up, down);

        for (Int_t ij = 0; ij < jetArr->GetEntriesFast(); ij++) {
            const baconhep::TJet* jet = (baconhep::TJet*)((*jetArr)[ij]);
            double pt_jet = jet->pt;
            double eta_jet = jet->eta;
            double phi_jet = jet->phi;
            if (pt_jet < 20.)
                continue;
            if (fabs(eta_jet) < 2.)
                continue;
            if (fabs(eta_jet) > 3.)
                continue;
            if (!((1 - jet->chEmFrac - jet->neuEmFrac - jet->chHadFrac - jet->neuHadFrac) < 0.50))
                continue;

            //Loop over photons to remove overlap
            double nonprefiringprobfromoverlappingphotons = 1.;
            double nonprefiringprobfromoverlappingphotons_up = 1.;
            double nonprefiringprobfromoverlappingphotons_down = 1.;
            bool foundOverlappingPhotons = false;
            for (Int_t ip = 0; ip < scArr->GetEntriesFast(); ip++) {
                const baconhep::TPhoton* photon = (baconhep::TPhoton*)((*scArr)[ip]);
                double pt_gam = photon->pt;
                double eta_gam = photon->eta;
                double phi_gam = photon->phi;
                if (pt_gam < 20.)
                    continue;
                if (fabs(eta_gam) < 2.)
                    continue;
                if (fabs(eta_gam) > 3.)
                    continue;
                double dR = toolbox::deltaR(eta_jet, phi_jet, eta_gam, phi_gam);
                if (dR > 0.4)
                    continue;
                double prefiringprob_gam = 0., prefiringprob_gam_up = 0., prefiringprob_gam_down = 0.;
                getPrefiringRateEcal(eta_gam, pt_gam, h_prefmap_photon_, prefiringprob_gam, prefiringprob_gam_up, prefiringprob_gam_down);
                nonprefiringprobfromoverlappingphotons *= (1. - prefiringprob_gam);
                nonprefiringprobfromoverlappingphotons_up *= (1. - prefiringprob_gam_up);
                nonprefiringprobfromoverlappingphotons_down *= (1. - prefiringprob_gam_down);
                foundOverlappingPhotons = true;
            }
            double prefiringprobfromoverlappingjet = 0.0;
            double prefiringprobfromoverlappingjet_up = 0.0;
            double prefiringprobfromoverlappingjet_down = 0.0;
            getPrefiringRateEcal(eta_jet, pt_jet, h_prefmap_jet_, prefiringprobfromoverlappingjet, prefiringprobfromoverlappingjet_up, prefiringprobfromoverlappingjet_down);
            if (!foundOverlappingPhotons) {
                main *= 1 - prefiringprobfromoverlappingjet;
                up *= 1 - prefiringprobfromoverlappingjet_up;
                down *= 1 - prefiringprobfromoverlappingjet_down;
            }
            //If overlapping photons have a non prefiring rate larger than the jet, then replace these weights by the jet one
            else if (nonprefiringprobfromoverlappingphotons > 1 - prefiringprobfromoverlappingjet) {
                if (nonprefiringprobfromoverlappingphotons > 0.) {
                    main *= (1 - prefiringprobfromoverlappingjet) / nonprefiringprobfromoverlappingphotons;
                    up *= (1 - prefiringprobfromoverlappingjet_up) / nonprefiringprobfromoverlappingphotons_up;
                    down *= (1 - prefiringprobfromoverlappingjet_down) / nonprefiringprobfromoverlappingphotons_down;
                } else {
                    main = 0.;
                    up = 0.;
                    down = 0.;
                }
            }
            //Last case: if overlapping photons have a non prefiring rate smaller than the jet, don't consider the jet in the event weight, and do nothing.
        }
        return;
    }

    void computeMuonsOnly(float& main, float& up, float& down, float& up_stat, float& down_stat, float& up_syst, float& down_syst)
    {
        for (Int_t im = 0; im < muonArr->GetEntriesFast(); im++) {
            const baconhep::TMuon* muon = (baconhep::TMuon*)((*muonArr)[im]);
            // Remove crappy tracker muons which would not have prefired the L1 trigger
            // if (pt < 5 || !muon.isLooseMuon())
            if (muon->pt < 5 || passMuonLooseID(muon))
                continue;
            double prate_central = 0., prate_up = 0., prate_down = 0., prate_statup = 0., prate_statdown = 0., prate_systup = 0., prate_systdown = 0.;
            getPrefiringRateMuon(muon->eta, muon->phi, muon->pt, prate_central, prate_up, prate_down, prate_statup, prate_statdown, prate_systup, prate_systdown);
            main *= 1 - prate_central;
            up *= 1 - prate_up;
            down *= 1 - prate_down;
            up_stat *= 1 - prate_statup;
            down_stat *= 1 - prate_statdown;
            up_syst *= 1 - prate_systup;
            down_syst *= 1 - prate_systdown;
        }
        return;
    }

    void getPrefiringRateEcal(double eta,
        double pt,
        TH2F* h_prefmap,
        double& prate_central,
        double& prate_up,
        double& prate_down)
    {
        //Check pt is not above map overflow
        int nbinsy = h_prefmap->GetNbinsY();
        double maxy = h_prefmap->GetYaxis()->GetBinLowEdge(nbinsy + 1);
        if (pt >= maxy)
            pt = maxy - 0.01;
        int thebin = h_prefmap->FindBin(eta, pt);

        double prefrate = h_prefmap->GetBinContent(thebin);

        double statuncty = h_prefmap->GetBinError(thebin);
        double systuncty = 0.2 * prefrate;

        prate_central = prefrate;
        prate_up = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
        prate_down = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));

        CheckRange(prate_central);
        CheckRange(prate_up);
        CheckRange(prate_down);

        return;
    }

    void getPrefiringRateMuon(double eta,
        double phi,
        double pt,
        double& prate_central,
        double& prate_up,
        double& prate_down,
        double& prate_statup,
        double& prate_statdown,
        double& prate_systup,
        double& prate_systdown)
    {
        double prefrate = 0.;
        double statuncty = 0.;
        if (std::abs(eta) < 0.2) {
            prefrate = parametrization0p0To0p2_->Eval(pt);
            statuncty = parametrization0p0To0p2_->GetParError(2);
        } else if (std::abs(eta) < 0.3) {
            prefrate = parametrization0p2To0p3_->Eval(pt);
            statuncty = parametrization0p2To0p3_->GetParError(2);
        } else if (std::abs(eta) < 0.55) {
            prefrate = parametrization0p3To0p55_->Eval(pt);
            statuncty = parametrization0p3To0p55_->GetParError(2);
        } else if (std::abs(eta) < 0.83) {
            prefrate = parametrization0p55To0p83_->Eval(pt);
            statuncty = parametrization0p55To0p83_->GetParError(2);
        } else if (std::abs(eta) < 1.24) {
            prefrate = parametrization0p83To1p24_->Eval(pt);
            statuncty = parametrization0p83To1p24_->GetParError(2);
        } else if (std::abs(eta) < 1.4) {
            prefrate = parametrization1p24To1p4_->Eval(pt);
            statuncty = parametrization1p24To1p4_->GetParError(2);
        } else if (std::abs(eta) < 1.6) {
            prefrate = parametrization1p4To1p6_->Eval(pt);
            statuncty = parametrization1p4To1p6_->GetParError(2);
        } else if (std::abs(eta) < 1.8) {
            prefrate = parametrization1p6To1p8_->Eval(pt);
            statuncty = parametrization1p6To1p8_->GetParError(2);
        } else if (std::abs(eta) < 2.1) {
            prefrate = parametrization1p8To2p1_->Eval(pt);
            statuncty = parametrization1p8To2p1_->GetParError(2);
        } else if (std::abs(eta) < 2.25) {
            prefrate = parametrization2p1To2p25_->Eval(pt);
            statuncty = parametrization2p1To2p25_->GetParError(2);
        } else if (std::abs(eta) < 2.4) {
            prefrate = parametrization2p25To2p4_->Eval(pt);
            statuncty = parametrization2p25To2p4_->GetParError(2);
        } else {
            //std::cout << "Muon outside of |eta| <= 2.4. Prefiring weight set to 0. Muon pt " << pt << " eta " << eta << std::endl;
        }
        double systuncty = 0.2 * prefrate;

        prate_central = prefrate;
        prate_up = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
        prate_down = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
        prate_systup = std::min(1., prefrate + systuncty);
        prate_systdown = std::max(0., prefrate - systuncty);
        prate_statup = std::min(1., prefrate + statuncty);
        prate_statdown = std::max(0., prefrate - statuncty);

        CheckRange(prate_central);
        CheckRange(prate_up);
        CheckRange(prate_down);
        CheckRange(prate_systup);
        CheckRange(prate_systdown);
        CheckRange(prate_statup);
        CheckRange(prate_statdown);
    }

private:
    TString fname;
    TString era;
    TH2F* h_prefmap_photon_;
    TH2F* h_prefmap_jet_;

    TClonesArray* scArr;
    TClonesArray* jetArr;
    TClonesArray* muonArr;
    TFile* f;
    TFile* fmuon;

    TF1* parametrization0p0To0p2_;
    TF1* parametrization0p2To0p3_;
    TF1* parametrization0p3To0p55_;
    TF1* parametrization0p55To0p83_;
    TF1* parametrization0p83To1p24_;
    TF1* parametrization1p24To1p4_;
    TF1* parametrization1p4To1p6_;
    TF1* parametrization1p6To1p8_;
    TF1* parametrization1p8To2p1_;
    TF1* parametrization2p1To2p25_;
    TF1* parametrization2p25To2p4_;
    TF1* parametrizationHotSpot_;

    bool etaCut(double eta)
    {
        if (fabs(eta) < 2 || fabs(eta) > 3)
            return false;
        return true;
    }

    bool oob(double eff)
    {
        if (eff < 0 || eff > 1)
            return true;
        return false;
    }

    void CheckRange(double& prate)
    {
        if (prate > 1.) {
            std::cout << "Found a prefiring probability > 1. Setting to 1." << std::endl;
            prate = 1.;
        }
    }
};

#endif
