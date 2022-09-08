//================================================================================================
//
// Perform fits to recoil against Z->mm, Z->ee, W->emet, and W->mmet events
//
//  * Outputs a ROOT file of functions parametrizing the distribution of recoil components
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TLorentzVector.h" // 4-vector class
#include <TF1.h> // 1D function
#include <TFile.h> // file handle class
#include <TFitResult.h> // class to handle fit results
#include <TGraphErrors.h> // graph class
#include <TLatex.h>
#include <TTree.h> // class to access ntuples
#include <fstream> // standard I/O
#include <iostream> // standard I/O
#include <sstream>

#include "Utils/CPlot.hh" // helper class for plots
#include "Utils/MitStyleRemix.hh" // style settings for drawing


#include "Math/Minimizer.h"
#include "Math/MinimizerOptions.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooRealIntegral.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "RooProdPdf.h"
#include "RooPolynomial.h"

#include "RooCrystalBall.h"

#endif

using namespace RooFit;
using namespace std;

bool doLikelihoodScan = false;
bool doLog = false; // true for data; false for MC
bool do_keys = 0;
bool do_5TeV = 0;
int etaBinCategory = 0;

//=== FUNCTION DECLARATIONS ======================================================================================
//--------------------------------------------------------------------------------------------------

// perform fit of recoil component
void performFit(const vector<TH1D*> hv, const vector<TH1D*> hbkgv, const Double_t* ptbins, const Int_t nbins,
    const Int_t model, const Bool_t sigOnly,
    const vector<RooDataSet> lDataSet, const vector<RooRealVar> lVar,
    TCanvas* c, const char* plabel, const char* xlabel,
    Double_t* mean1Arr, Double_t* mean1ErrArr,
    Double_t* mean2Arr, Double_t* mean2ErrArr,
    Double_t* mean3Arr, Double_t* mean3ErrArr,
    Double_t* sigma0Arr, Double_t* sigma0ErrArr,
    Double_t* sigma1Arr, Double_t* sigma1ErrArr,
    Double_t* sigma2Arr, Double_t* sigma2ErrArr,
    Double_t* sigma3Arr, Double_t* sigma3ErrArr,
    Double_t* frac2Arr, Double_t* frac2ErrArr,
    Double_t* frac3Arr, Double_t* frac3ErrArr,
    Double_t* chi2Arr, Double_t* chi2ErrArr,
    Double_t* alphaLArr, Double_t* alphaLErrArr,
    Double_t* alphaRArr, Double_t* alphaRErrArr,
    Double_t* nLArr,     Double_t* nLErrArr,
    Double_t* nRArr,     Double_t* nRErrArr,
    RooWorkspace* workspace,
    const char* outputname,
    int etaBinCategory, bool do_keys, bool doCrystalBall, Double_t luminosity);

//=== MAIN MACRO =================================================================================================

void fitRecoil_master(
    Int_t pfu1model = 2, // u1 model (0=> CrystalBall + Gaussian, 1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
    Int_t pfu2model = 2, // u2 model (0=> CrystalBall + Gaussian, 1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
    Bool_t sigOnly = 1, // 1=> signal event only, 0=> add background
    Bool_t doElectron = 0, // 0=> Muon channels, 1=> Electron channels
    Double_t inputType = 0, // Which dataset to run on 0=> Data, 1=> Z->ll MC, 2=> W+ jets MC, 3=> W- jets MC
    std::string inputDirectory = "/ceph/moh/CROWN_samples/EarlyRun3_V12",
    std::string metVar = "met_uncorrected",
    std::string metPhiVar = "metphi_uncorrected",
    TString outputDir = "./", // output directory
    Double_t lumi = 1)
{
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);

    //--------------------------------------------------------------------------------------------------------------
    // Settings and files
    //==============================================================================================================

    TString outputname("");
    outputname += TString(metVar);
    if (inputType==0) { outputname += TString("_data"); } else if (inputType==1) {outputname += TString("_ZllMC");}
    else if (inputType==2) { outputname += TString("_WposMC"); } else if (inputType==3) { outputname += TString("_WnegMC"); };
    if (pfu1model==0) {outputname += TString("_CB");} else if (pfu1model==1) {outputname += TString("_single");}
    else if (pfu1model==2) {outputname += TString("_double");} else if (pfu1model==3) {outputname += TString("_triple");};
    if (doElectron) { outputname += TString("_elec"); } else { outputname += TString("_muon"); } ;
    if (sigOnly) { outputname += TString("_sigOnly"); } else { outputname += TString("_sigAndBck"); };
    CPlot::sOutDir = outputDir + TString("/") + outputname + TString("/plots");
    outputDir =  outputDir + TString("/") + outputname;
    gSystem->mkdir(outputDir, kTRUE); 

    // preserving the fine binning at low pT but the higher-pT bins (>75 GeV have been adjusted to be slightly wider)
    //Double_t ptbins[] = { 0, 5.0, 25, 100, 1000 };
    //Double_t ptbins[] = {45, 47};
    Double_t ptbins[] = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000};
    //Double_t ptbins[] = {0,3.0,6.0,12.5,20,27.5,37.5,50,60,75,90,125,200,1000}; 
    Int_t nbins = sizeof(ptbins) / sizeof(Double_t) - 1;

    TString formulaPFu1mean("pol2");
    TString formulaPFu2mean("pol2");
    TString formulaPFu1meanScale("pol2");
    TString formulaPFu2meanScale("pol2");

    // Instantiate booleans (modified later)
    bool sigIsData = false;
    bool doCrystalBall = false;

    vector<TChain*> fnamev;
    vector<Bool_t> isBkgv;

    TChain *bkgv_chain = new TChain();
    TChain *sig_chain = new TChain();
    TChain *sig_friend_chain = new TChain();
    TChain *sig_friend_chain2 = new TChain();
    TChain *bkgv_friend_chain = new TChain();
    TChain *bkgv_friend_chain2 = new TChain();

    TString infilename = inputDirectory + "/ntuples/2022/";
    TString infilename_friend = infilename;
    TString infilename_friend2 = infilename;
    infilename_friend = infilename_friend.ReplaceAll("ntuples", "friends/crosssection");
    infilename_friend2 = infilename_friend2.ReplaceAll("ntuples", "friends/lep_corr").ReplaceAll("moh", "jdriesch");
    TString bkgvdirname;
    TString bkgvdirname2;
    if (doElectron)
        bkgvdirname = inputDirectory + "/ntuples/2022/TTTo2L2Nu_CP5_13p6TeV_powheg-pythia8-Run3Winter22MiniAOD-122X/ee/";
    if (!doElectron)
        bkgvdirname = inputDirectory + "/ntuples/2022/TTTo2L2Nu_CP5_13p6TeV_powheg-pythia8-Run3Winter22MiniAOD-122X/mm/";
    TString bkgvdirname_friend = bkgvdirname;
    TString bkgvdirname_friend2 = bkgvdirname;
    bkgvdirname_friend = bkgvdirname_friend.ReplaceAll("ntuples", "friends/crosssection");
    bkgvdirname_friend2 = bkgvdirname_friend2.ReplaceAll("ntuples", "friends/lep_corr").ReplaceAll("moh", "jdriesch");

    //Data samples
    if(inputType == 0) {
        sigIsData = true;
        if(!doElectron) {
            sig_chain->Add(infilename + "SingleMuon_Run2022C-PromptReco-v1/mm/*.root/ntuple");
            sig_friend_chain->Add(infilename_friend + "SingleMuon_Run2022C-PromptReco-v1/mm/*.root/ntuple");
            sig_friend_chain2->Add(infilename_friend2 + "SingleMuon_Run2022C-PromptReco-v1/mm/*.root/ntuple");

            sig_chain->Add(infilename + "Muon_Run2022C-PromptReco-v1/mm/*.root/ntuple");
            sig_friend_chain->Add(infilename_friend + "Muon_Run2022C-PromptReco-v1/mm/*.root/ntuple");
            sig_friend_chain2->Add(infilename_friend2 + "Muon_Run2022C-PromptReco-v1/mm/*.root/ntuple");
        }
        if(doElectron) {
            sig_chain->Add(infilename + "EGamma_Run2022C-PromptReco-v1/ee/*.root/ntuple");
            sig_friend_chain->Add(infilename_friend + "EGamma_Run2022C-PromptReco-v1/ee/*.root/ntuple");
            sig_friend_chain2->Add(infilename_friend2 + "EGamma_Run2022C-PromptReco-v1/ee/*.root/ntuple");
        }
    }
    //Z->ll SIMULATION SAMPLE
    if (inputType == 1) {
        if(doElectron) {
            sig_chain->Add(infilename + "DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/ee/*.root/ntuple");
            sig_friend_chain->Add(infilename_friend + "DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/ee/*.root/ntuple");
            sig_friend_chain2->Add(infilename_friend2 + "DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/ee/*.root/ntuple");
        } else {
            sig_chain->Add(infilename + "DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/mm/*.root/ntuple");
            sig_friend_chain->Add(infilename_friend + "DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/mm/*.root/ntuple");
            sig_friend_chain2->Add(infilename_friend2 + "DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/mm/*.root/ntuple");
        }    
    }
    // W SIMULATION SAMPLE
    if (inputType > 1) {
        if (!doElectron) {
            sig_chain->Add(infilename + "WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/mmet/*.root/ntuple");
            sig_friend_chain->Add(infilename_friend + "WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/mmet/*.root/ntuple");
            sig_friend_chain2->Add(infilename_friend2 + "WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/mmet/*.root/ntuple");
        }
        if (doElectron) {
            sig_chain->Add(infilename + "WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/emet/*.root/ntuple");
            sig_friend_chain->Add(infilename_friend + "WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/emet/*.root/ntuple");
            sig_friend_chain2->Add(infilename_friend2 + "WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X/emet/*.root/ntuple");
        }
    }

    bkgv_chain->Add(bkgvdirname + "*.root/ntuple");
    bkgv_friend_chain->Add(bkgvdirname_friend + "*.root/ntuple");
    bkgv_friend_chain2->Add(bkgvdirname_friend2 + "*.root/ntuple");

    bkgv_chain->AddFriend(bkgv_friend_chain);
    bkgv_chain->AddFriend(bkgv_friend_chain2);
    sig_chain->AddFriend(sig_friend_chain);
    sig_chain->AddFriend(sig_friend_chain2);

    //pushing the chains to the event loop
    fnamev.push_back(sig_chain); isBkgv.push_back(kFALSE);
    if (!sigOnly) {
        fnamev.push_back(bkgv_chain); isBkgv.push_back(kTRUE);
    }

    const Double_t MASS_LOW = 60;
    const Double_t MASS_HIGH = 120;
    const Double_t PT_CUT = 25;
    const Double_t ETA_CUT = 2.4;
    Double_t lep_MASS = 0.1057; //Muon Mass
    if (doElectron) {
        lep_MASS = 0.000511; //Electron Mass
    }

    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================

    char hname[100];
    vector<TH1D*> hPTv;
    vector<TH1D*> hPFu1v, hPFu1Bkgv;
    vector<TH1D*> hPFu2v, hPFu2Bkgv;
    vector<RooDataSet> lDataSetU1;
    vector<RooDataSet> lDataSetU2;

    vector<RooRealVar> vu1Var;
    vector<RooRealVar> vu2Var;
    vector<RooRealVar> vptVar;

    RooWorkspace pdfsU1("pdfsU1");
    RooWorkspace pdfsU2("pdfsU2");

    for (Int_t ibin = 0; ibin < nbins; ibin++) {

        // Puppi-PF
        int range = 100;
        if (ptbins[ibin] > 80)
            range = 125;
        if (ptbins[ibin] > 150)
            range = 150;

        sprintf(hname, "hPFu1_%i", ibin);
        hPFu1v.push_back(new TH1D(hname, "", 100, -range - ptbins[ibin], range - ptbins[ibin]));
        hPFu1v[ibin]->Sumw2();
        sprintf(hname, "hPFu1Bkg_%i", ibin);
        hPFu1Bkgv.push_back(new TH1D(hname, "", 100, -range - ptbins[ibin], range - ptbins[ibin]));
        hPFu1Bkgv[ibin]->Sumw2();

        sprintf(hname, "hPFu2_%i", ibin);
        hPFu2v.push_back(new TH1D(hname, "", 100, -range, range));
        hPFu2v[ibin]->Sumw2();
        sprintf(hname, "hPFu2Bkg_%i", ibin);
        hPFu2Bkgv.push_back(new TH1D(hname, "", 100, -range, range));
        hPFu2Bkgv[ibin]->Sumw2();

        sprintf(hname, "hPT_%i", ibin);
        hPTv.push_back(new TH1D(hname, "", 100, -range + ptbins[ibin], range + ptbins[ibin]));
        //hPTv[ibin]->Sumw2();

        std::stringstream name;
        name << "u_" << ibin;

        RooRealVar u1Var(name.str().c_str(), name.str().c_str(), 0, -range - ptbins[ibin], range - ptbins[ibin]);
        RooRealVar u2Var(name.str().c_str(), name.str().c_str(), 0, -range, range);
        u1Var.setBins(100);
        u2Var.setBins(100);

        vu1Var.push_back(u1Var);
        vu2Var.push_back(u2Var);

        //RooRealVar ptVar(name.str().c_str(), name.str().c_str(), 0, -range - ptbins[ibin], range - ptbins[ibin]);
        //ptVar.setBins(100);
        //ptVar.push_back(ptVar)

        sprintf(hname, "hDataSetU1_%i", ibin);
        RooDataSet dataSetU1(hname, hname, RooArgSet(u1Var));
        lDataSetU1.push_back(dataSetU1);
        sprintf(hname, "hDataSetU2_%i", ibin);
        RooDataSet dataSetU2(hname, hname, RooArgSet(u2Var));
        lDataSetU2.push_back(dataSetU2);
    }

    TFile* infile = 0;
    TTree* intree = 0;

    UInt_t category;
    UInt_t runNum;
    Float_t genVPt, genVPhi, genVy, genVMass;
    Float_t lep1_pt_uncorr, lep1_eta, lep1_phi, lep1_mass;
    Float_t lep2_pt_uncorr, lep2_eta, lep2_phi, lep2_mass;
    Double_t lep1_pt, lep2_pt;
    Float_t scale1fb;
    Float_t met, metPhi;
    Int_t q1, q2;
    Float_t scale1fbSumW, genWeight;
    Bool_t trg_match_1, trg_match_2;
    TLorentzVector *dilep = 0, *lep1 = 0, *lep2 = 0, *lep1_raw = 0, *lep2_raw = 0, *genlep1 = 0, *genlep2 = 0;

    for (UInt_t ifile = 0; ifile < fnamev.size(); ifile++) {

        intree = fnamev.at(ifile);
        bool isData = sigIsData && !isBkgv[ifile];

        intree->SetBranchAddress("run", &runNum);
        intree->SetBranchAddress("genbosonpt", &genVPt); // GEN boson pT (signal MC)
        intree->SetBranchAddress("genbosonphi", &genVPhi); // GEN boson phi (signal MC)
        intree->SetBranchAddress("genbosonrapidity", &genVy); // GEN boson rapidity (signal MC)
        intree->SetBranchAddress("genbosonmass", &genVMass); // GEN boson mass (signal MC)
        if (!isData){
            intree->SetBranchAddress("scale1fb_sumw", &scale1fbSumW);
            intree->SetBranchAddress("genweight", &genWeight);
        }
        intree->SetBranchAddress(TString(metVar), &met); // Uncorrected MET
        intree->SetBranchAddress(TString(metPhiVar), &metPhi); // phi(MET)
        intree->SetBranchAddress("q_1", &q1); // charge of tag lepton
        intree->SetBranchAddress("pt_1", &lep1_pt_uncorr); // uncorrected pt of tag lepton
        intree->SetBranchAddress("pt_1_corr", &lep1_pt); // corrected pt of tag lepton
        intree->SetBranchAddress("eta_1", &lep1_eta); // eta of tag lepton
        intree->SetBranchAddress("phi_1", &lep1_phi); // phi of tag lepton
        intree->SetBranchAddress("mass_1", &lep1_mass); // mass of tag lepton
        if (doElectron)
            intree->SetBranchAddress("trg_single_ele27_1", &trg_match_1);
        else
            intree->SetBranchAddress("trg_single_mu24_1", &trg_match_1);

        if (inputType<2) { // Z specific variables
            intree->SetBranchAddress("q_2", &q2); // charge of probe lepton
            intree->SetBranchAddress("pt_2", &lep2_pt_uncorr); // uncorrected pt of probe lepton
            intree->SetBranchAddress("pt_2_corr", &lep2_pt); // corrected pt of probe lepton
            intree->SetBranchAddress("eta_2", &lep2_eta); // eta of probe lepton
            intree->SetBranchAddress("phi_2", &lep2_phi); // phi of probe lepton
            intree->SetBranchAddress("mass_2", &lep2_mass); // mass of probe lepton
            if (doElectron)
                intree->SetBranchAddress("trg_single_ele27_2", &trg_match_2);
            else
                intree->SetBranchAddress("trg_single_mu24_2", &trg_match_2);
        }

        //
        // Loop over events
        //
        int iterator = 1;
        for (Int_t ientry = 0; ientry < intree->GetEntries(); ientry += iterator) {
            intree->GetEntry(ientry);

            if (isData) { //proper scaling
                scale1fb=1;
            } else {
                scale1fb = scale1fbSumW*genWeight;
            } 
            
            /////////
            /// Creating leptons 
            //
            TLorentzVector mu1, mu2;
            TLorentzVector dl;
            double mll, etall, ptll;

            // cout << "368" << endl;
            // cout << lep_MASS <<endl;
            // cout << lep1_mass << endl;
            // cout << lep2_mass << endl;
            // cout << lep1_pt << endl;
            // cout << lep1_eta << endl;
            // cout << lep1_phi << endl;
            // cout << lep2_pt << endl;
            // cout << lep2_eta << endl;
            // cout << lep2_phi << endl;

            if (inputType < 2 ) { //Z events
                mu1.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, lep1_mass);
                mu2.SetPtEtaPhiM(lep2_pt, lep2_eta, lep2_phi, lep2_mass);

                dl = mu1 + mu2;
                mll = dl.M();
                etall = dl.Rapidity();
                ptll = dl.Pt();
            } else { //W Events
                mu1.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, lep1_mass);
            }


            /////////
            /// Applying Cuts
            //
            if (inputType < 2 ){ //Z events
                if (!isBkgv[ifile]) {
                    if (mll < MASS_LOW || mll > MASS_HIGH)
                        continue;
                }
                if (mu1.Pt() < PT_CUT || mu2.Pt() < PT_CUT)
                    continue;
                if (fabs(mu1.Eta()) > ETA_CUT || fabs(mu2.Eta()) > ETA_CUT)
                    continue;
                if (q1*q2 > 0)
                    continue;
                if (!(trg_match_1 || trg_match_2))
                    continue;
            } else { //W events
                if (inputType == 2 && q1 < 0)
                    continue;
                if (inputType == 3 && q2 > 0)
                    continue;
                if (mu1.Pt() < PT_CUT)
                    continue;
                if (fabs(lep1_eta) > ETA_CUT)
                    continue;
                if (!trg_match_1)
                    continue;
            }

            //// TODO: implement this to do eta bin categorizing
            if (!isBkgv[ifile]) {
                if (!isData) {
                    if (etaBinCategory == 1 && fabs(genVy) > 0.5)
                        continue;
                    if (etaBinCategory == 2 && (fabs(genVy) <= 0.5 || fabs(genVy) >= 1))
                        continue;
                    if (etaBinCategory == 3 && fabs(genVy) < 1)
                        continue;
                } else {
                    if (etaBinCategory == 1 && fabs(etall) > 0.5)
                        continue;
                    if (etaBinCategory == 2 && (fabs(etall) <= 0.5 || fabs(etall) >= 1))
                        continue;
                    if (etaBinCategory == 3 && fabs(etall) < 1)
                        continue;
                }
            }

            /////////
            /// Boson pt binning
            //
            Int_t ipt = -1;
            if (inputType < 2) {
                for (Int_t ibin = 0; ibin < nbins; ibin++) {
                    if (ptll > ptbins[ibin] && ptll <= ptbins[ibin + 1])
                        ipt = ibin;
                }
            } else {
                for (Int_t ibin = 0; ibin < nbins; ibin++) {
                    if (genVPt > ptbins[ibin] && genVPt <= ptbins[ibin + 1])
                        ipt = ibin;
                }
            }
            if (ipt < 0)
                continue;

            /////////
            /// RECOIL Filling
            //
            double pU1 = 0.;
            double pU2 = 0.;
            double pUX = 0.;
            double pUY = 0.;
            double pU  = 0.;
            double pSin = 0.;
            double pCos = 0.;
            
            TVector2 vLepRaw1, vLepRaw2;
            TVector2 vLepCor1, vLepCor2;
            if (inputType < 2) { //Z Events
                vLepRaw1.Set(lep1_pt_uncorr * cos(lep1_phi), lep1_pt_uncorr * sin(lep1_phi));
                vLepRaw2.Set(lep2_pt_uncorr * cos(lep2_phi), lep2_pt_uncorr * sin(lep2_phi));

                vLepCor1.Set((lep1_pt) * cos(lep1_phi), (lep1_pt) * sin(lep1_phi));
                vLepCor2.Set((lep2_pt) * cos(lep2_phi), (lep2_pt) * sin(lep2_phi));


                TVector2 vMetCorr((met)*cos(metPhi), (met)*sin(metPhi));
                Double_t corrMetWithLepton = (vMetCorr + vLepRaw1 + vLepRaw2 - vLepCor1 - vLepCor2).Mod();
                Double_t corrMetWithLeptonPhi = (vMetCorr + vLepRaw1 + vLepRaw2 - vLepCor1 - vLepCor2).Phi();
                pUX = corrMetWithLepton * cos(corrMetWithLeptonPhi) + dl.Pt() * cos(dl.Phi());
                pUY = corrMetWithLepton * sin(corrMetWithLeptonPhi) + dl.Pt() * sin(dl.Phi());
                pU = sqrt(pUX * pUX + pUY * pUY);
                pCos = -(pUX * cos(dl.Phi()) + pUY * sin(dl.Phi())) / pU;
                pSin = (pUX * sin(dl.Phi()) - pUY * cos(dl.Phi())) / pU;

            } else {
                vLepRaw1.Set(lep1_pt_uncorr * cos(lep1_phi), lep1_pt_uncorr * sin(lep1_phi));
                vLepCor1.Set(lep1_pt * cos(lep1_phi), lep1_pt * sin(lep1_phi));

                TVector2 vMetCorr((met)*cos(metPhi), (met)*sin(metPhi));
                Double_t corrMetWithLepton = (vMetCorr + vLepRaw1 - vLepCor1).Mod();
                Double_t corrMetWithLeptonPhi = (vMetCorr + vLepRaw1 - vLepCor1).Phi();
                pUX = corrMetWithLepton * cos(corrMetWithLeptonPhi) + mu1.Pt() * cos(lep1_phi);
                pUY = corrMetWithLepton * sin(corrMetWithLeptonPhi) + mu1.Pt() * sin(lep1_phi);
                pU = sqrt(pUX * pUX + pUY * pUY);
                pCos = -(pUX * cos(genVPhi) + pUY * sin(genVPhi)) / pU;
                pSin = (pUX * sin(genVPhi) - pUY * cos(genVPhi)) / pU;
            }

            pU1 = pU * pCos; // U1 in data
            pU2 = pU * pSin; // U2 in data

            // cout << "boson pt: " << dl.Pt() << endl;
            // cout << "pu1: " << pU1 << endl;
            // cout << "pu2: " << pU2 << endl; 


            if (isBkgv[ifile]) {
                hPFu1Bkgv[ipt]->Fill(pU1, scale1fb * lumi);
                hPFu2Bkgv[ipt]->Fill(pU2, scale1fb * lumi);
                
            } else {
                if (isData) {
                    hPFu1v[ipt]->Fill(pU1, 1);
                    hPFu2v[ipt]->Fill(pU2, 1);
                }
                if(!isData) {
                    hPFu1v[ipt]->Fill(pU1, scale1fb * lumi);
                    hPFu2v[ipt]->Fill(pU2, scale1fb * lumi);
                }

                if (inputType<2) {
                    hPTv[ipt]->Fill(ptll);
                } else {
                    hPTv[ipt]->Fill(genVPt);
                }

                // this is the dataset for the RooKey
                // clean the under/overflow
                int range = 100;
                if (ptbins[ipt] > 80)
                    range = 125;
                if (ptbins[ipt] > 150)
                    range = 150;

                if (pU1 < (-range - ptbins[ipt]))
                    continue;
                if (pU1 > (range - ptbins[ipt]))
                    continue;
                if (pU2 < (-range))
                    continue;
                if (pU2 > (range))
                    continue;

                vu1Var[ipt].setVal(pU1);
                vu2Var[ipt].setVal(pU2);

                lDataSetU1[ipt].add(RooArgSet(vu1Var[ipt])); // need to add the weights
                lDataSetU2[ipt].add(RooArgSet(vu2Var[ipt]));
            }
        }
        delete infile;
        infile = 0, intree = 0;
    }

    TCanvas* c = MakeCanvas("c", "c", 800, 800);

    Double_t xval[nbins], xerr[nbins];
    for (Int_t ibin = 0; ibin < nbins; ibin++) {
        xval[ibin] = 0.5 * (ptbins[ibin + 1] + ptbins[ibin]);
        xerr[ibin] = 0.5 * (ptbins[ibin + 1] - ptbins[ibin]);
    }

    Double_t dataMean1[nbins], dataMeanError1[nbins]; 
    Double_t dataMean2[nbins], dataMeanError2[nbins];
    Double_t dataRMS1[nbins], dataRMS2[nbins];
    Double_t ptMean[nbins];

    for (Int_t ibin = 0; ibin < nbins; ibin++) {
        dataMean1[ibin] = -1 * hPFu1v[ibin]->GetMean() / hPTv[ibin]->GetMean();
        dataMean2[ibin] = -1 * hPFu2v[ibin]->GetMean() / hPTv[ibin]->GetMean();

        ptMean[ibin] = hPTv[ibin]->GetMean();

        dataRMS1[ibin] =  hPFu1v[ibin]->GetRMS();
        dataRMS2[ibin] =  hPFu2v[ibin]->GetRMS();
    }

    TGraph *grPFu1mean=0;
    TGraph *grPFu2mean=0;
    TGraph *grPFu1meanerr=0;
    TGraph *grPFu2meanerr=0;

    grPFu1mean = new TGraph(nbins,ptMean,dataMean1);
    grPFu1mean ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu1mean ->GetYaxis()->SetRangeUser(0., 2.5);
    grPFu1mean ->SetName("grPFu1mean");
    CPlot plotPFu1mean("pfu1mean","","Dimuon p_{T} [GeV/c]","<u_{#parallel}>/<p_{T}>");
    plotPFu1mean.AddTextBox("CMS        Hadronic Recoil Response   (13.0 TeV)", 0.1, 0.92, 0.95, 0.97, 0, kBlack, -1);
    plotPFu1mean.AddGraph(grPFu1mean,"",kBlack,kOpenCircle);
    //latexLabel.DrawLatex(0.20, 0.8, label);
    plotPFu1mean.Draw(c,kTRUE,"png");

    grPFu2mean = new TGraph(nbins,ptMean,dataMean2);
    grPFu2mean ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu2mean ->GetYaxis()->SetRangeUser(-.1, .1);
    grPFu2mean ->SetName("grPFu2mean");
    CPlot plotPFu2mean("pfu2mean","","Dimuon p_{T} [GeV/c]","<u_{#perp}  >/<p_{T}>");
    plotPFu2mean.AddTextBox( "CMS           Hadronic Recoil Response   (13.0 TeV)", 0.1, 0.92, 0.95, 0.97, 0, kBlack, -1);
    plotPFu2mean.AddGraph(grPFu2mean,"",kBlack,kOpenCircle);
    //latexLabel.DrawLatex(0.20, 0.8, label);
    plotPFu2mean.Draw(c,kTRUE,"png");

    grPFu1meanerr = new TGraph(nbins,ptMean,dataRMS1);
    grPFu1meanerr ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu1meanerr ->GetYaxis()->SetRangeUser(0, 30);
    grPFu1meanerr ->SetName("grPFu1meanerr");
    CPlot plotPFu1meanerr("pfu1meanerr","","Dimuon p_{T} [GeV/c]","RMS(u_{#parallel})");
    plotPFu1meanerr.AddTextBox("CMS     Hadronic Recoil Resolution  (13.0 TeV)", 0.1, 0.92, 0.95, 0.97, 0, kBlack, -1);
    plotPFu1meanerr.AddGraph(grPFu1meanerr,"",kBlack,kOpenCircle);
    //latexLabel.DrawLatex(0.20, 0.8, label);
    plotPFu1meanerr.Draw(c,kTRUE,"png");

    grPFu2meanerr = new TGraph(nbins,ptMean,dataRMS2);
    grPFu2meanerr ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu2meanerr ->GetYaxis()->SetRangeUser(0, 30);
    grPFu2meanerr ->SetName("grPFu2meanerr");
    CPlot plotPFu2meanerr("pfu2meanerr","","Dimuon p_{T} [GeV/c]","RMS(u_{#perp}  )");
    plotPFu2meanerr.AddTextBox("CMS        Hadronic Recoil Resolution  (13.0 TeV)", 0.1, 0.92, 0.95, 0.97, 0, kBlack, -1);
    plotPFu2meanerr.AddGraph(grPFu2meanerr,"",kBlack,kOpenCircle);
    //latexLabel.DrawLatex(0.20, 0.8, label);
    plotPFu2meanerr.Draw(c,kTRUE,"png");

    cout << "Mean Plots Created" << endl;

    //
    // Arrays and graphs to store fit results
    //
    Double_t pfu1Mean[nbins], pfu1MeanErr[nbins];
    Double_t pfu1Mean2[nbins], pfu1Mean2Err[nbins];
    Double_t pfu1Mean3[nbins], pfu1Mean3Err[nbins];
    Double_t pfu1MeanScale[nbins], pfu1MeanErrScale[nbins];
    Double_t pfu1Mean2Scale[nbins], pfu1Mean2ErrScale[nbins];
    Double_t pfu1Mean3Scale[nbins], pfu1Mean3ErrScale[nbins];
    Double_t pfu1Sigma0[nbins], pfu1Sigma0Err[nbins];
    Double_t pfu1Sigma1[nbins], pfu1Sigma1Err[nbins];
    Double_t pfu1Sigma2[nbins], pfu1Sigma2Err[nbins];
    Double_t pfu1Sigma3[nbins], pfu1Sigma3Err[nbins];
    Double_t pfu1Frac2[nbins], pfu1Frac2Err[nbins];
    Double_t pfu1Frac3[nbins], pfu1Frac3Err[nbins];
    Double_t pfu1chi2[nbins], pfu1chi2Err[nbins];
    Double_t pfu1alphaLArr[nbins], pfu1alphaLErrArr[nbins];
    Double_t pfu1alphaRArr[nbins], pfu1alphaRErrArr[nbins];
    Double_t pfu1nLArr[nbins],     pfu1nLErrArr[nbins];
    Double_t pfu1nRArr[nbins],     pfu1nRErrArr[nbins];

    TGraphErrors *grPFu1chi2=0;
    ;

    Double_t pfu2Mean[nbins], pfu2MeanErr[nbins];
    Double_t pfu2Mean2[nbins], pfu2Mean2Err[nbins];
    Double_t pfu2Mean3[nbins], pfu2Mean3Err[nbins];
    Double_t pfu2Sigma0[nbins], pfu2Sigma0Err[nbins];
    Double_t pfu2Sigma1[nbins], pfu2Sigma1Err[nbins];
    Double_t pfu2Sigma2[nbins], pfu2Sigma2Err[nbins];
    Double_t pfu2Sigma3[nbins], pfu2Sigma3Err[nbins];
    Double_t pfu2Frac2[nbins], pfu2Frac2Err[nbins];
    Double_t pfu2Frac3[nbins], pfu2Frac3Err[nbins];
    Double_t pfu2chi2[nbins], pfu2chi2Err[nbins];
    Double_t pfu2alphaLArr[nbins], pfu2alphaLErrArr[nbins];
    Double_t pfu2alphaRArr[nbins], pfu2alphaRErrArr[nbins];
    Double_t pfu2nLArr[nbins],     pfu2nLErrArr[nbins];
    Double_t pfu2nRArr[nbins],     pfu2nRErrArr[nbins];

    TGraphErrors *grPFu2chi2=0;

    char outpdfname[50];
    sprintf(outpdfname, "%s/%s.root", outputDir.Data(), "pdfsU1");

    doCrystalBall = false;
    if (pfu1model==0) {
        doCrystalBall = true;
    }

    // Fitting PF-MET u1
    performFit(hPFu1v, hPFu1Bkgv, ptbins, nbins, pfu1model, sigOnly,
        lDataSetU1, vu1Var,
        c, "pfu1", "u_{#parallel} [GeV]",
        pfu1Mean, pfu1MeanErr,
        pfu1Mean2, pfu1Mean2Err,
        pfu1Mean3, pfu1Mean3Err,
        pfu1Sigma0, pfu1Sigma0Err,
        pfu1Sigma1, pfu1Sigma1Err,
        pfu1Sigma2, pfu1Sigma2Err,
        pfu1Sigma3, pfu1Sigma3Err,
        pfu1Frac2, pfu1Frac2Err,
        pfu1Frac3, pfu1Frac3Err,
        pfu1chi2, pfu1chi2Err,
        pfu1alphaLArr, pfu1alphaLErrArr,
        pfu1alphaRArr, pfu1alphaRErrArr,
        pfu1nLArr,     pfu1nLErrArr,
        pfu1nRArr,     pfu1nRErrArr,
        &pdfsU1,
        outpdfname,
        etaBinCategory, do_keys, doCrystalBall, lumi);

    pdfsU1.writeToFile(outpdfname, kFALSE);

    grPFu1chi2 = new TGraphErrors(nbins,xval,pfu1chi2,xerr,pfu1chi2Err);
    grPFu1chi2 ->GetYaxis()->SetRangeUser(0., 10.);
    grPFu1chi2 ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu1chi2 ->SetName("grPFu1chi2");
    CPlot plotPFu1chi2("pfu1chi2","","p_{T} [GeV/c]","#chi^{2}(u_{#parallel})");
    plotPFu1chi2.AddGraph(grPFu1chi2,"",kBlack,kOpenCircle);
    //latexLabel.DrawLatex(0.20, 0.8, label);
    plotPFu1chi2.Draw(c,kTRUE,"png");

    sprintf(outpdfname, "%s/%s.root", outputDir.Data(), "pdfsU2");

    doCrystalBall = false;
    if (pfu2model==0) {
        doCrystalBall = true;
    }

    // Fitting PF-MET u2
    performFit(hPFu2v, hPFu2Bkgv, ptbins, nbins, pfu2model, sigOnly,
        lDataSetU2, vu2Var,
        c, "pfu2", "u_{#perp  } [GeV/c]",
        pfu2Mean, pfu2MeanErr,
        pfu2Mean2, pfu2Mean2Err,
        pfu2Mean3, pfu2Mean3Err,
        pfu2Sigma0, pfu2Sigma0Err,
        pfu2Sigma1, pfu2Sigma1Err,
        pfu2Sigma2, pfu2Sigma2Err,
        pfu2Sigma3, pfu2Sigma3Err,
        pfu2Frac2, pfu2Frac2Err,
        pfu2Frac3, pfu2Frac3Err,
        pfu2chi2, pfu2chi2Err,
        pfu2alphaLArr, pfu2alphaLErrArr,
        pfu2alphaRArr, pfu2alphaRErrArr,
        pfu2nLArr,     pfu2nLErrArr,
        pfu2nRArr,     pfu2nRErrArr,
        &pdfsU2,
        outpdfname,
        etaBinCategory, do_keys, doCrystalBall, lumi);

    pdfsU2.writeToFile(outpdfname, kFALSE);

    grPFu2chi2 = new TGraphErrors(nbins,xval,pfu2chi2,xerr,pfu2chi2Err);
    grPFu2chi2 ->GetYaxis()->SetRangeUser(0., 10.);
    grPFu2chi2 ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu2chi2 ->SetName("grPFu2chi2");
    CPlot plotPFu2chi2("pfu2chi2","","p_{T} [GeV/c]","#chi^{2}(u_{#perp} ) [GeV]");
    plotPFu2chi2.AddGraph(grPFu2chi2,"",kBlack,kOpenCircle);
    //latexLabel.DrawLatex(0.20, 0.8, label);
    plotPFu2chi2.Draw(c,kTRUE,"png");

    cout << endl;
    cout << "  <> Output saved in " << outputDir << "/" << endl;
    cout << endl;

    return;
}

//--------------------------------------------------------------------------------------------------
void performFit(const vector<TH1D*> hv, const vector<TH1D*> hbkgv, const Double_t* ptbins, const Int_t nbins,
    const Int_t model, const Bool_t sigOnly,
    vector<RooDataSet> lDataSet, vector<RooRealVar> lVar,
    TCanvas* c, const char* plabel, const char* xlabel,
    Double_t* mean1Arr, Double_t* mean1ErrArr,
    Double_t* mean2Arr, Double_t* mean2ErrArr,
    Double_t* mean3Arr, Double_t* mean3ErrArr,
    Double_t* sigma0Arr, Double_t* sigma0ErrArr,
    Double_t* sigma1Arr, Double_t* sigma1ErrArr,
    Double_t* sigma2Arr, Double_t* sigma2ErrArr,
    Double_t* sigma3Arr, Double_t* sigma3ErrArr,
    Double_t* frac2Arr, Double_t* frac2ErrArr,
    Double_t* frac3Arr, Double_t* frac3ErrArr,
    Double_t* chi2Arr, Double_t* chi2ErrArr,

    Double_t* alphaLArr, Double_t* alphaLErrArr,
    Double_t* alphaRArr, Double_t* alphaRErrArr,
    Double_t* nLArr,     Double_t* nLErrArr,
    Double_t* nRArr,     Double_t* nRErrArr,

    RooWorkspace* wksp,
    const char* outputname,
    int etaBinCategory, bool do_keys, bool doCrystalBall, 
    Double_t luminosity)
{
    char lumi[50];
    char pname[50];
    char ylabel[50];
    char xlabel_new[50];
    char binlabel[50];
    char binYlabel[50];
    char nsigtext[50];
    char nbkgtext[50];

    char mean1text[50];
    char mean2text[50];
    char mean3text[50];
    char sig0text[50];
    char sig1text[50];
    char sig2text[50];
    char sig3text[50];

    char alphaLtext[50];
    char alphaRtext[50];
    char nLtext[50];
    char nRtext[50];

    char frac2text[50];

    double frac2_ini = 0;
    double frac3_ini = 0;

    double sigma1_ini = 0;
    double sigma2_ini = 0;
    double sigma3_ini = 0;

    for (Int_t ibin = 0; ibin < nbins; ibin++) {

        TH1D* hvLOG = (TH1D*)hv[ibin]->Clone();
        for (Int_t i = 0; i < hv[ibin]->GetNbinsX() + 1; i++) {
            if (hv[ibin]->GetBinContent(i) == 0) {
                hvLOG->SetBinContent(i, 0);
                hvLOG->SetBinError(i, 0);
            } else {
                hvLOG->SetBinContent(i, TMath::Log10(hv[ibin]->GetBinContent(i)));
                hvLOG->SetBinError(i, 0);
            }
        }

        TH1D* hbkgvLOG = (TH1D*)hbkgv[ibin]->Clone();
        for (Int_t i = 0; i < hbkgv[ibin]->GetNbinsX() + 1; i++) {
            if (hbkgv[ibin]->GetBinContent(i) == 0) {
                hbkgvLOG->SetBinContent(i, 0);
                hbkgvLOG->SetBinError(i, 0);
            } else {
                hbkgvLOG->SetBinContent(i, TMath::Log10(hbkgv[ibin]->GetBinContent(i)));
                hbkgvLOG->SetBinError(i, 0);
            }
        }

        std::stringstream name;
        // unfortunately have to give each variable individual names for each bin
        name << "u_" << ibin;
        RooRealVar u(name.str().c_str(), name.str().c_str(), hv[ibin]->GetXaxis()->GetXmin(), hv[ibin]->GetXaxis()->GetXmax());
        name.str("");
        u.setBins(100);

        //    HERE THE LOG histo
        RooDataHist dataHist("dataHist", "dataHist", RooArgSet(u), hv[ibin]);
        RooDataHist dataHistLog("dataHistLog", "dataHistLog", RooArgSet(u), hvLOG);

        //
        // Set up background histogram templates
        //
        RooDataHist bkgHist("bkgHist", "bkgHist", RooArgSet(u), hbkgv[ibin]);
        RooDataHist bkgHistLog("bkgHistLog", "bkgHistLog", RooArgSet(u), hbkgvLOG);
        //    RooHistPdf bkg("bkg","bkg",u,bkgHist,0);
        name.str("");
        name << "bkg_" << ibin << std::endl;
        RooHistPdf bkg(name.str().c_str(), name.str().c_str(), u, bkgHist, 0);
        name.str("");
        name.str("");
        name << "bkgLog_" << ibin << std::endl;
        RooHistPdf bkgLog(name.str().c_str(), name.str().c_str(), u, bkgHistLog, 0);
        name.str("");

        //
        // Set up fit parameters
        //
        name.str("");
        name << "mean1_" << ibin;
        RooRealVar mean1(name.str().c_str(), name.str().c_str(),
            hv[ibin]->GetMean(),
            hv[ibin]->GetXaxis()->GetXmin() + 50,
            hv[ibin]->GetXaxis()->GetXmax() - 50);
        

        name.str("");
        name << "mean2_" << ibin;
        RooRealVar mean2(name.str().c_str(), name.str().c_str(),
            hv[ibin]->GetMean(),
            hv[ibin]->GetXaxis()->GetXmin() + 50,
            hv[ibin]->GetXaxis()->GetXmax() - 50);
        name.str("");
        name << "mean3_" << ibin;
        RooRealVar mean3(name.str().c_str(), name.str().c_str(),
            hv[ibin]->GetMean() * 0.85,
            hv[ibin]->GetXaxis()->GetXmin() + 50,
            hv[ibin]->GetXaxis()->GetXmax() - 50);

        name.str("");
        name << "sigma1_" << ibin;
        RooRealVar sigma1(name.str().c_str(), name.str().c_str(), 0.3 * (hv[ibin]->GetRMS()), 0.1 * hv[ibin]->GetRMS(), 2.5 * (hv[ibin]->GetRMS()));
        name.str("");
        name << "sigma2_" << ibin;
        RooRealVar sigma2(name.str().c_str(), name.str().c_str(), 1.0 * (hv[ibin]->GetRMS()), 0.1 * hv[ibin]->GetRMS(), 4.5 * (hv[ibin]->GetRMS()));
        name.str("");
        name << "sigma3_" << ibin;
        RooRealVar sigma3(name.str().c_str(), name.str().c_str(), 2.0 * (hv[ibin]->GetRMS()), 0.1 * hv[ibin]->GetRMS(), 9 * hv[ibin]->GetRMS());
        // fraction
        name.str("");
        name << "frac2_" << ibin;
        RooRealVar frac2(name.str().c_str(), name.str().c_str(), 0.5, 0.15, 0.85);
        name.str("");
        name << "frac3_" << ibin;
        RooRealVar frac3(name.str().c_str(), name.str().c_str(), 0.05, 0., 0.15);

        if (model == 2) {
            frac2.setVal(0.5);
            frac2.setRange(0., 1.);
        }

        if (string(plabel) == string("pfu2")) {
            mean1.setVal(0);
            mean1.setRange(-5., 5.);
            mean1.setConstant(kTRUE);
            mean2.setVal(0);
            mean2.setRange(-5., 5.);
            mean2.setConstant(kTRUE);
            mean3.setVal(0);
            mean3.setRange(-5., 5.);
            mean3.setConstant(kTRUE);
        }

        name.str("");
        name << "gauss1_" << ibin;
        RooGaussian gauss1(name.str().c_str(), name.str().c_str(), u, mean1, sigma1);
        name.str("");
        name << "gauss2_" << ibin;
        RooGaussian gauss2(name.str().c_str(), name.str().c_str(), u, mean2, sigma2);
        name.str("");
        name.str("");
        name << "gauss3_" << ibin;
        RooGaussian gauss3(name.str().c_str(), name.str().c_str(), u, mean3, sigma3);
        name.str("");

        RooGaussian constGauss1("constGauss1", "constGauss1", mean1, RooConst(hv[ibin]->GetMean()), RooConst(0.15 * hv[ibin]->GetRMS()));
        RooGaussian constGauss2("constGauss2", "constGauss2", mean2, RooConst(hv[ibin]->GetMean()), RooConst(0.15 * hv[ibin]->GetRMS()));
        RooGaussian constGauss3("constGauss3", "constGauss3", mean3, RooConst(0.85 * hv[ibin]->GetMean()), RooConst(0.15 * hv[ibin]->GetRMS()));
        
        //cout << " Mean of const: " << hv[ibin]->GetMean() << ", Std Dev of const: " << 0.001 * hv[ibin]->GetRMS() << endl;
        RooGaussian constGauss_forCB("constGauss_forCB", "constGauss_forCB", mean1,  RooConst(hv[ibin]->GetMean()), RooConst(0.15 * hv[ibin]->GetRMS()));

        name.str("");
        name << "sigmaL_" << ibin;
        RooRealVar sigmaL(name.str().c_str(), name.str().c_str(),  (hv[ibin]->GetRMS()), 0.1 * hv[ibin]->GetRMS(), 2.5 * (hv[ibin]->GetRMS()));
        name.str("");
        name << "sigmaR_" << ibin;
        RooRealVar sigmaR(name.str().c_str(), name.str().c_str(),  (hv[ibin]->GetRMS()), 0.1 * hv[ibin]->GetRMS(), 2.5 * (hv[ibin]->GetRMS()));
        name.str("");
        name << "alphaL_" << ibin;
        RooRealVar alphaL(name.str().c_str(), name.str().c_str(), 5, 0, 100);
        name.str("");
        name << "alphaR_" << ibin;
        RooRealVar alphaR(name.str().c_str(), name.str().c_str(), 5, 0, 100);
        name.str("");
        name << "nL_" << ibin;
        RooRealVar nL(name.str().c_str(), name.str().c_str(), 10, 1, 1000);
        name.str("");
        name << "nR_" << ibin;
        RooRealVar nR(name.str().c_str(), name.str().c_str(), 10, 1, 1000);

        name.str("");
        name << "size_" << ibin;
        RooRealVar size(name.str().c_str(), name.str().c_str(), 861947, 861946, 861948);

        name.str("");
        name << "doubleCB_" << ibin;
        RooCrystalBall doubleCB(name.str().c_str(), name.str().c_str(), u, mean1, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
        name.str("");


        if (string(plabel) == string("pfu2")) {
            mean1.setVal(0);
            mean1.setRange(-5., 5.);
            mean1.setConstant(kTRUE);
        }

        //
        // Define formula for overall width (sigma0)
        //
        char formula[100];
        RooArgList params;

        if (!doCrystalBall) {
            if(model==1) {
                sprintf(formula,"sigma1_%d", ibin);
            } else if(model==2) {
                sprintf(formula, "(1.-frac2_%d)*sigma1_%d + frac2_%d*sigma2_%d", ibin, ibin, ibin, ibin);
                params.add(frac2);
                params.add(sigma1);
                params.add(sigma2);
            } else if(model==3) {
                sprintf(formula, "(1.-frac2_%d-frac3_%d)*sigma1_%d + frac2_%d*sigma2_%d + frac3_%d*sigma3_%d", ibin, ibin, ibin, ibin, ibin, ibin, ibin);
                params.add(frac2);
                params.add(frac3);
                params.add(sigma1);
                params.add(sigma2);
                params.add(sigma3);
            }

        }
        if (doCrystalBall) {
            sprintf(formula, "sigmaL_%d + sigmaR_%d", ibin, ibin);
            params.add(sigmaL);
            params.add(sigmaR);
        }

        RooFormulaVar sigma0("sigma0", formula, params);

        //
        // Construct fit model
        //
        RooArgList shapes;
        RooArgList fracs;
        if (!doCrystalBall) {
            if (model >= 3)
                shapes.add(gauss3);
            if (model >= 2)
                shapes.add(gauss2);
            shapes.add(gauss1);

            if (model >= 3)
                fracs.add(frac3);
            if (model >= 2)
                fracs.add(frac2);
        } 
        if (doCrystalBall) {
            //shapes.add(doubleCB);

            shapes.add(gauss1);  //consider swapping the order here
            shapes.add(doubleCB);
            fracs.add(frac2);
        }

        name.str("");
        name << "sig_" << ibin;
        RooAddPdf sig(name.str().c_str(), name.str().c_str(), shapes, fracs);
        name.str("");

        // name.str("");
        // name << "sigCB_" << ibin;
        // RooCrystalBall sigCB(name.str().c_str(), name.str().c_str(), u, mean1, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
        // name.str("");

        RooArgList parts;
        if(doCrystalBall) {
            parts.add(sig);
        }
        if(!doCrystalBall) {
            parts.add(sig);
        }

        if (!sigOnly && !doLog)
            parts.add(bkg);
        if (!sigOnly && doLog)
            parts.add(bkgLog);

        RooArgList yields;

        name.str("");
        name << "nsig_" << ibin;
        RooRealVar nsig(name.str().c_str(), name.str().c_str(), 0.98 * (hv[ibin]->Integral()), 0., 1.1 * hv[ibin]->Integral()); // just to be sure that doesn't it the boundary
        name.str("");
        name << "nbkg_" << ibin;
        RooRealVar nbkg(name.str().c_str(), name.str().c_str(), (hbkgv[ibin]->Integral()), 0.1 * (hbkgv[ibin]->Integral()), 0.25 * (hv[ibin]->Integral()));
        RooRealVar* lAbkgFrac = new RooRealVar("AbkgFrac", "AbkgFrac", (hv[ibin]->Integral() / (hv[ibin]->Integral() + hbkgv[ibin]->Integral())), 0.9, .98);

        if (sigOnly) {
            yields.add(nsig);
            nbkg.setVal(0);
        } else {
            RooFormulaVar* sigbkgFrac = new RooFormulaVar("bkgfrac", "@0", RooArgSet(*lAbkgFrac));
            yields.add(*sigbkgFrac);
        }

        name.str("");
        name << "modelpdf_" << ibin << std::endl;
        RooAddPdf modelpdf(name.str().c_str(), name.str().c_str(), parts, yields);
        name.str("");

        RooFitResult* fitResultLOG = 0;
        if (doLog) {
            if(!doCrystalBall) {
                fitResultLOG = modelpdf.fitTo(dataHistLog,
                    NumCPU(4),
                    Minimizer("Minuit2", "minimize"),
                    ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //				    ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
            }
            if(doCrystalBall) {
                fitResultLOG = modelpdf.fitTo(dataHistLog,
                    NumCPU(4),
                    Minimizer("Minuit2", "minimize"),
                    //ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //				    ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
            }
        }
        
        RooFitResult* fitResult = 0;

        if(!doCrystalBall) {
            fitResult = modelpdf.fitTo(dataHist,
                NumCPU(4),
                Minimizer("Minuit2", "minimize"),
                ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                RooFit::Minos(),
                RooFit::Strategy(2),
                RooFit::Save());
            int nTries = 0;

            do {
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    //				 Minimizer("Minuit2","minimize"),
                    Minimizer("Minuit2", "scan"),
                    ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());

                // nTries++;

                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    //				 Minimizer("Minuit2","minimize"),
                    Minimizer("Minuit2", "migrad"),
                    ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    RooFit::Hesse(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    //				 Minimizer("Minuit2","minimize"),
                    Minimizer("Minuit2", "improve"),
                    ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());

                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "minimize"),
                    ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                nTries++;
            } while ((fitResult->status() > 0 || fitResult->covQual() < 3) && nTries < 10);
        }
        if (doCrystalBall) {
            fitResult = modelpdf.fitTo(dataHist,
                NumCPU(4),
                Minimizer("Minuit2", "minimize"),
                //ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                RooFit::Minos(),
                RooFit::Strategy(2),
                RooFit::Save());
            int nTries = 0;

            do {
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    //				 Minimizer("Minuit2","minimize"),
                    Minimizer("Minuit2", "scan"),
                    //ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());

                // nTries++;

                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    //				 Minimizer("Minuit2","minimize"),
                    Minimizer("Minuit2", "migrad"),
                    //ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                    RooFit::Hesse(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    //				 Minimizer("Minuit2","minimize"),
                    Minimizer("Minuit2", "improve"),
                    //ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //				 ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());

                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "minimize"),
                    //ExternalConstraints(RooArgSet(constGauss1, constGauss2, constGauss3)), //ExternalConstraints(constGauss2), ExternalConstraints(constGauss3),
                    //			       ExternalConstraints(constGauss2),ExternalConstraints(constGauss3),
                    ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                nTries++;
            } while ((fitResult->status() > 0 || fitResult->covQual() < 3) && nTries < 10);
        }

        char rname[100];
        if (string(plabel) == string("pfu1"))
            sprintf(rname, "fitResultU1_%i", ibin);
        if (string(plabel) == string("pfu2"))
            sprintf(rname, "fitResultU2_%i", ibin);

        fitResult->SetName(rname);

        TFile* lFile = TFile::Open(outputname, "UPDATE");
        fitResult->Write();
        lFile->Write();
        lFile->Close();

        c->SetFillColor(kWhite);
        if (fitResult->status() > 0)
            c->SetFillColor(kYellow);

        wksp->import(u);
        wksp->import(modelpdf);
        if(!doCrystalBall)
            wksp->import(sig);
        if(doCrystalBall)
            wksp->import(sig);  // wksp->import(sigCB);
        if (!doLog)
            wksp->import(bkg);
        if (doLog)
            wksp->import(bkgLog);

        mean1Arr[ibin] = mean1.getVal();
        mean1ErrArr[ibin] = mean1.getError();
        if(!doCrystalBall) {
            frac2_ini = frac2.getVal();
            frac3_ini = frac3.getVal();
            sigma1_ini = sigma1.getVal();
            sigma2_ini = sigma2.getVal();
            sigma3_ini = sigma3.getVal();

            if (model >= 2) {
                mean2Arr[ibin] = mean2.getVal();
                mean2ErrArr[ibin] = mean2.getError();
                sigma0Arr[ibin] = sigma0.getVal();
                sigma0ErrArr[ibin] = sigma0.getPropagatedError(*fitResult);
                sigma2Arr[ibin] = sigma2.getVal();
                sigma2ErrArr[ibin] = sigma2.getError();
                frac2Arr[ibin] = frac2.getVal();
                frac2ErrArr[ibin] = frac2.getError();
            }
            if (model >= 3) {
                mean3Arr[ibin] = mean3.getVal();
                mean3ErrArr[ibin] = mean3.getError();
                sigma3Arr[ibin] = sigma3.getVal();
                sigma3ErrArr[ibin] = sigma3.getError();
                frac3Arr[ibin] = frac3.getVal();
                frac3ErrArr[ibin] = frac3.getError();
            }
        }
        if (doCrystalBall) {
            sigma0Arr[ibin] = sigma0.getVal();
            sigma0ErrArr[ibin] = sigma0.getPropagatedError(*fitResult);
            sigma1Arr[ibin] = sigmaL.getVal();
            sigma1ErrArr[ibin] = sigmaL.getError();
            sigma2Arr[ibin] = sigmaR.getVal();
            sigma2ErrArr[ibin] = sigmaR.getError();

            alphaLArr[ibin] = alphaL.getVal();
            alphaLErrArr[ibin] = alphaL.getError();
            alphaRArr[ibin] = alphaR.getVal();
            alphaRErrArr[ibin] = alphaR.getError();

            nLArr[ibin] = nL.getVal();
            nLErrArr[ibin] = nL.getError();
            nRArr[ibin] = nR.getVal();
            nRErrArr[ibin] = nR.getError();

            //Addition of Gaussian
            sigma3Arr[ibin] = sigma1.getVal();
            sigma3ErrArr[ibin] = sigma1.getError();
            frac2Arr[ibin] = frac2.getVal();
            frac2ErrArr[ibin] = frac2.getError();
        }

        std::cout << "Plot Fit results " << std::endl;
        //
        // Plot fit results
        //
        RooPlot* frame = u.frame(Bins(100));
        dataHist.plotOn(frame, MarkerStyle(kFullCircle), MarkerSize(0.8), DrawOption("ZP"));
        modelpdf.plotOn(frame);

        if (!doLog)
            name.str("");
        name << "bkg_" << ibin;
        if (doLog)
            name.str("");
        name << "bkgLog_" << ibin;
        if (!sigOnly)
            modelpdf.plotOn(frame, Components(bkg), FillColor(kRed), DrawOption("F"));

        if(!doCrystalBall) {
            name.str("");
            name << "gauss1_" << ibin;
            if (model >= 2)
                sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kRed));
            name.str("");
            name << "gauss2_" << ibin;
            if (model >= 2)
                sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kMagenta));
            name.str("");
            name << "gauss3_" << ibin;
            if (model >= 3)
                sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kGreen + 2));
        }
        if(doCrystalBall) {
            //name.str("");
            //name << "doubleCB_" << ibin;
            //sigCB.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kRed));

            name.str("");
            name << "doubleCB_" << ibin;
            sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kRed));

            name.str("");
            name << "gauss1_" << ibin;
            sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kGreen));
        }

        
        if (!doLog)
            name.str("");
        name << "bkg_" << ibin;
        if (doLog)
            name.str("");
        name << "bkgLog_" << ibin;

        // draw the curve
        if (!sigOnly)
            modelpdf.plotOn(frame, FillColor(kGray), VisualizeError(*fitResult, 1), RooFit::Components(modelpdf)); // 1 sigma band
        if (!sigOnly)
            modelpdf.plotOn(frame, RooFit::LineColor(kGray + 2));
        if (!doCrystalBall) {
            sig.plotOn(frame, FillColor(7), VisualizeError(*fitResult, 1), RooFit::Components(sig)); // 1 sigma band
            sig.plotOn(frame, RooFit::LineColor(kBlue));
        } 
        if (doCrystalBall) {
            //sigCB.plotOn(frame, FillColor(7), VisualizeError(*fitResult, 1), RooFit::Components(sig)); // 1 sigma band
            //sigCB.plotOn(frame, RooFit::LineColor(kBlue));

            sig.plotOn(frame, FillColor(7), VisualizeError(*fitResult, 1), RooFit::Components(sig)); // 1 sigma band
            sig.plotOn(frame, RooFit::LineColor(kBlue));
        }  

        // redraw the data
        dataHist.plotOn(frame, MarkerStyle(kFullCircle), MarkerSize(0.8), DrawOption("ZP"));

        if (do_keys) {

            //      lDataSet[ibin].Print();
            name.str("");
            name << "key_" << ibin;
            //      RooKeysPdf * pdf_keys = new RooKeysPdf(name.str().c_str(),name.str().c_str(), u, dataHist, RooKeysPdf::NoMirror, 2);
            RooKeysPdf pdf_keys(name.str().c_str(), name.str().c_str(), lVar[ibin], lDataSet[ibin], RooKeysPdf::NoMirror, 2);

            RooPlot* xframe = lVar[ibin].frame(Title(Form("%s Zp_{T}=%.1f - %.1f GeV/c ", plabel, ptbins[ibin], ptbins[ibin + 1])));

            lDataSet[ibin].plotOn(xframe);
            TCanvas* c = new TCanvas("validatePDF", "validatePDF", 800, 800);
            c->cd();
            pdf_keys.plotOn(xframe, LineColor(kBlue));
            xframe->Draw();

            c->SaveAs(Form("%s_%d_dataset.png", plabel, ibin));

            name.str("");

            pdf_keys.plotOn(frame, LineColor(kRed));

            //      wksp->import(lDataSet[ibin]);
            wksp->import(pdf_keys);
            //      wksp->import(lVar[ibin],RooFit::RecycleConflictNodes(),RooFit::Silence());
            //      wksp->import(pdf_keys,RooFit::RecycleConflictNodes(),RooFit::Silence());

            wksp->Print();
        }

        int sizeParam = 0;
        if (string(plabel) == string("pfu1"))
            sizeParam = 8; // 3 means + 3 sigma + 2 frac
        if (string(plabel) == string("pfu2"))
            sizeParam = 5; // 0 means + 3 sigma + 2 frac

        TString nameRooHist = Form("h_%s", dataHist.GetName());
        TString nameRooCurve = Form("sig_%d_Norm[u_%d]", ibin, ibin);

        RooHist* hist = frame->getHist(nameRooHist.Data());
        RooCurve* fitCurve = frame->getCurve(nameRooCurve.Data());

        RooHist* hist_pull = hist->makePullHist(*fitCurve);
        hist_pull->SetTitle("");
        hist_pull->GetYaxis()->SetTitle("pull");
        hist_pull->GetYaxis()->SetRangeUser(-5., 5.);
        hist_pull->SetMarkerColor(kAzure);
        hist_pull->SetLineColor(kAzure);
        hist_pull->SetFillColor(kAzure);
        //    hist_pull->GetYaxis()->SetTitleFont(42);
        //    hist_pull->GetXaxis()->SetTitleFont(42);
        hist_pull->GetYaxis()->SetTitleSize(0.055);
        hist_pull->GetYaxis()->SetTitleOffset(1.600);
        hist_pull->GetYaxis()->SetLabelOffset(0.014);
        hist_pull->GetYaxis()->SetLabelSize(0.050);
        hist_pull->GetYaxis()->SetLabelFont(42);
        hist_pull->GetXaxis()->SetTitleSize(0.055);
        hist_pull->GetXaxis()->SetTitleOffset(1.300);
        hist_pull->GetXaxis()->SetLabelOffset(0.014);
        hist_pull->GetXaxis()->SetLabelSize(0.050);
        hist_pull->GetXaxis()->SetLabelFont(42);

        chi2Arr[ibin] = frame->chiSquare(nameRooCurve.Data(), nameRooHist.Data(), sizeParam);
        chi2ErrArr[ibin] = 0;
        if (chi2Arr[ibin] > 10) {
            chi2Arr[ibin] = 0;
            chi2ErrArr[ibin] = 200;
        } // just a larger number so that is easy to notice on the plot
        //    cout << " chi2Arr[ibin]=" << chi2Arr[ibin] << " chi2ErrArr[ibin]=" << chi2ErrArr[ibin] << endl;

        // if(do_5TeV) sprintf(lumi,"CMS                               27.4 pb^{-1} (5 TeV)");
        // if(!do_5TeV) sprintf(lumi,"CMS                               2.2 fb^{-1} (13 TeV)");
        sprintf(lumi,"CMS Preliminary                    %.1f fb^{-1} (13 TeV)", luminosity);
        //sprintf(lumi,"CMS Preliminary                     92.8 fb^{-1} (13 TeV)");
        sprintf(pname, "%sfit_%i", plabel, ibin);
        sprintf(ylabel, "Events / %.1f GeV", hv[ibin]->GetBinWidth(1));
        if (string(plabel) == string("pfu1")) {
            sprintf(xlabel_new, "u_{#parallel}");
        } else {
            sprintf(xlabel_new, "u_{#perp}");
        }
        sprintf(binlabel, "p_{T}(Z) = %.1f - %.1f GeV/c", ptbins[ibin], ptbins[ibin + 1]);

        if (etaBinCategory == 1)
            sprintf(binYlabel, "|y| < 0.5");
        if (etaBinCategory == 2)
            sprintf(binYlabel, "0.5 < |y| < 1");
        if (etaBinCategory == 3)
            sprintf(binYlabel, "|y| > 1");

        if (sigOnly) {
            sprintf(nsigtext, "N_{evts} = %i", (Int_t)hv[ibin]->Integral());
        } else {
            sprintf(nsigtext, "N_{sig}/(N_{bkg}+N_{sig}) = %.3f #pm %.3f", lAbkgFrac->getVal(), lAbkgFrac->getError());
            //      sprintf(nsigtext,"N_{sig} = %.1f #pm %.1f",nsig.getVal(),nsig.getError());
            //      sprintf(nbkgtext,"N_{bkg} = %.1f #pm %.1f",nbkg.getVal(),nbkg.getError());
        }
        if (doCrystalBall) {
            sprintf(mean1text, "#mu_{1} = %.1f #pm %.1f", mean1Arr[ibin], mean1ErrArr[ibin]);
            sprintf(sig1text, "#sigma_{L} = %.1f #pm %.1f", sigma1Arr[ibin], sigma1ErrArr[ibin]);
            sprintf(sig2text, "#sigma_{R} = %.1f #pm %.1f", sigma2Arr[ibin], sigma2ErrArr[ibin]);
            sprintf(sig0text, "#sigma_{0} = %.1f #pm %.1f", sigma0Arr[ibin], sigma0ErrArr[ibin]);
            sprintf(alphaLtext, "#alpha_{L} = %.1f #pm %.1f", alphaLArr[ibin], alphaLErrArr[ibin]);
            sprintf(alphaRtext, "#alpha_{R} = %.1f #pm %.1f", alphaRArr[ibin], alphaRErrArr[ibin]);
            sprintf(nLtext, "n_{L} = %.1f #pm %.1f", nLArr[ibin], nLErrArr[ibin]);
            sprintf(nRtext, "n_{R} = %.1f #pm %.1f", nRArr[ibin], nRErrArr[ibin]);

            sprintf(sig3text, "sigma_{G} = %.1f #pm %.1f", sigma3Arr[ibin], sigma3ErrArr[ibin]);
            sprintf(frac2text, "fraction = %.1f #pm %.1f", frac2Arr[ibin], frac2ErrArr[ibin]);
        }
        if (!doCrystalBall) {
            sprintf(mean1text, "#mu_{1} = %.1f #pm %.1f", mean1Arr[ibin], mean1ErrArr[ibin]);
            sprintf(sig1text, "#sigma = %.1f #pm %.1f", sigma1Arr[ibin], sigma1ErrArr[ibin]);
            if (model >= 2) {
                sprintf(mean2text, "#mu_{2} = %.1f #pm %.1f", mean2Arr[ibin], mean2ErrArr[ibin]);
                //      sprintf(mean2text,"#mu_{2} = #mu_{1} ");
                sprintf(sig0text, "#sigma = %.1f #pm %.1f", sigma0Arr[ibin], sigma0ErrArr[ibin]);
                //sprintf(sig1text, "#sigma_{1} = %.1f #pm %.1f", sigma1Arr[ibin], sigma1ErrArr[ibin]);
                sprintf(sig2text, "#sigma_{2} = %.1f #pm %.1f", sigma2Arr[ibin], sigma2ErrArr[ibin]);
            }
            if (model >= 3) {
                sprintf(mean3text, "#mu_{3} = %.1f #pm %.1f", mean3Arr[ibin], mean3ErrArr[ibin]);
                sprintf(sig3text, "#sigma_{3} = %.1f #pm %.1f", sigma3Arr[ibin], sigma3ErrArr[ibin]);
            }
        }

        ///////////
        /////////// Draw Linear
        ///////////
        ///////////

        TCanvas* cLin = MakeCanvas("cLin", "cLin", 800, 800);
        cLin->Divide(1, 2, 0, 0);
        cLin->cd(1)->SetPad(0, 0.3, 1.0, 1.0);
        cLin->cd(1)->SetTopMargin(0.1);
        cLin->cd(1)->SetBottomMargin(0.01);
        cLin->cd(1)->SetLeftMargin(0.15);
        cLin->cd(1)->SetRightMargin(0.07);
        cLin->cd(1)->SetTickx(1);
        cLin->cd(1)->SetTicky(1);
        cLin->cd(2)->SetPad(0, 0, 1.0, 0.3);
        cLin->cd(2)->SetTopMargin(0.05);
        cLin->cd(2)->SetBottomMargin(0.45);
        cLin->cd(2)->SetLeftMargin(0.15);
        cLin->cd(2)->SetRightMargin(0.07);
        cLin->cd(2)->SetTickx(1);
        cLin->cd(2)->SetTicky(1);

        CPlot plot(pname, frame, "", xlabel_new, ylabel);
        //    pad1->cd();
        plot.AddTextBox(lumi, 0.1, 0.92, 0.95, 0.97, 0, kBlack, -1);
        plot.AddTextBox(binlabel, 0.21, 0.80, 0.51, 0.85, 0, kBlack, -1);
        plot.AddTextBox(xlabel_new, 0.75, 0.1, 1.00, 0, 0, kBlack, -1);
        if (etaBinCategory != 0)
            plot.AddTextBox(binYlabel, 0.21, 0.85, 0.51, 0.9, 0, kBlack, -1);
        if (sigOnly)
            plot.AddTextBox(nsigtext, 0.21, 0.78, 0.51, 0.73, 0, kBlack, -1);
        //    else        plot.AddTextBox(0.21,0.78,0.51,0.68,0,kBlack,-1,2,nsigtext,nbkgtext);
        else
        plot.AddTextBox(nsigtext, 0.21, 0.78, 0.51, 0.73, 0, kBlack, -1); // this print the fraction now
        if(!doCrystalBall) {
            if (model == 1)
                plot.AddTextBox(0.70, 0.90, 0.95, 0.80, 0, kBlack, -1, 2, mean1text, sig1text);
            else if (model == 2)
                plot.AddTextBox(0.70, 0.90, 0.95, 0.70, 0, kBlack, -1, 5, mean1text, mean2text, sig0text, sig1text, sig2text);
            //    else if(model==3) plot.AddTextBox(0.70,0.90,0.95,0.65,0,kBlack,-1,7,mean1text,mean2text,mean3text,sig0text,sig1text,sig2text,sig3text);
            else if (model == 3)
                plot.AddTextBox(0.70, 0.90, 0.95, 0.65, 0, kBlack, -1, 7, mean1text, mean2text, mean3text, sig0text, sig1text, sig2text, sig3text);
        }
        if(doCrystalBall)
            plot.AddTextBox(0.70, 0.90, 0.95, 0.70, 0, kBlack, -1, 10, mean1text, sig0text, sig1text, sig2text, alphaLtext, alphaRtext, nLtext, nRtext, sig1text, frac2text);
            //plot.AddTextBox(0.70, 0.90, 0.95, 0.70, 0, kBlack, -1, 8, mean1text, sig0text, sig1text, sig2text, alphaLtext, alphaRtext, nLtext, nRtext);
        plot.Draw(cLin, kFALSE, "png", 1);
        plot.Draw(cLin, kFALSE, "pdf", 1);

        cLin->cd(2);
        hist_pull->Draw("A3 L ");
        TLine* lineZero = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0, hv[ibin]->GetXaxis()->GetXmax(), 0);
        lineZero->SetLineColor(kBlack);
        lineZero->Draw("same");
        TLine* lineZero1SigmaM = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0, hv[ibin]->GetXaxis()->GetXmax(), 0);
        lineZero1SigmaM->SetLineColor(11);
        lineZero1SigmaM->Draw("same");
        TLine* lineZero1SigmaP = new TLine(hv[ibin]->GetXaxis()->GetXmin(), 0, hv[ibin]->GetXaxis()->GetXmax(), 0);
        lineZero1SigmaP->SetLineColor(11);
        lineZero1SigmaP->Draw("same");

        plot.Draw(cLin, kTRUE, "png", 1);
        plot.Draw(cLin, kTRUE, "pdf", 1);

        TCanvas* c1 = MakeCanvas("c1", "c1", 800, 800);
        c1->Divide(1, 2, 0, 0);
        c1->cd(1)->SetPad(0, 0.3, 1.0, 1.0);
        c1->cd(1)->SetTopMargin(0.1);
        c1->cd(1)->SetBottomMargin(0.01);
        c1->cd(1)->SetLeftMargin(0.15);
        c1->cd(1)->SetRightMargin(0.07);
        c1->cd(1)->SetTickx(1);
        c1->cd(1)->SetTicky(1);
        c1->cd(2)->SetPad(0, 0, 1.0, 0.3);
        c1->cd(2)->SetTopMargin(0.05);
        c1->cd(2)->SetBottomMargin(0.45);
        c1->cd(2)->SetLeftMargin(0.15);
        c1->cd(2)->SetRightMargin(0.07);
        c1->cd(2)->SetTickx(1);
        c1->cd(2)->SetTicky(1);

        sprintf(pname, "%sfitlog_%i", plabel, ibin);
        plot.SetYRange(0.1, 10 * hv[ibin]->GetMaximum());
        //plot.SetYRange(0.001, 10 * hv[ibin]->GetMaximum());
        plot.SetName(pname);
        plot.SetLogy();
        plot.Draw(c1, kFALSE, "png", 1);
        plot.Draw(c1, kFALSE, "pdf", 1);

        c1->cd(2);
        hist_pull->SetTitle("");
        hist_pull->GetYaxis()->SetTitle("pull");
        hist_pull->GetYaxis()->SetRangeUser(-5., 5.);
        hist_pull->SetMarkerColor(kAzure);
        hist_pull->SetLineColor(kAzure);
        hist_pull->SetFillColor(kAzure);    
        hist_pull->Draw("A3 L ");
        lineZero->SetLineColor(kBlack);
        lineZero->Draw("same");
        lineZero1SigmaM->SetLineColor(11);
        lineZero1SigmaM->Draw("same");
        lineZero1SigmaP->SetLineColor(11);
        lineZero1SigmaP->Draw("same");
        plot.Draw(c1, kTRUE, "png", 1);
        plot.Draw(c1, kTRUE, "pdf", 1);

        // reset color canvas
        c->SetFillColor(kWhite);
    }

    return;
}