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

//#include "RooCrystalBall.h"

#endif

using namespace RooFit;
using namespace std;

int etaBinCategory = 0;
bool do_keys = 0;

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

void fitRecoil(
    Int_t pfu1model = 2, // u1 model (0=> CrystalBall + Gaussian, 1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
    Int_t pfu2model = 2, // u2 model (0=> CrystalBall + Gaussian, 1 => single Gaussian, 2 => double Gaussian, 3 => triple Gaussian)
    Bool_t sigOnly = 1, // 1=> signal event only, 0=> add background
    Bool_t doElectron = 0, // 0=> Muon channels, 1=> Electron channels
    Double_t inputType = 0, // Which dataset to run on 0=> Data, 1=> Z->ll MC, 2=> W+ jets MC, 3=> W- jets MC
    std::string inputDirectory = "/ceph/jdriesch/CROWN_samples/Run3V07",
    std::string metVar = "pfmet_uncorrected",
    std::string metPhiVar = "pfmetphi_uncorrected",
    TString outputDir = "./", // output directory
    Double_t lumi = 1,
    Int_t rapbin = -1,
    TString hist_file = "/work/jdriesch/earlyrun3/Z_early_Run3/output/earlyRun3_2022_Zpt.root")
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
    TString zrapbin = "";
    if (rapbin > -1) {
        outputname += TString::Format("_zrap%d", rapbin);
        zrapbin = TString::Format("_zrap%d", rapbin);
    }
    CPlot::sOutDir = outputDir + TString("/") + outputname + TString("/plots");
    outputDir =  outputDir + TString("/") + outputname;
    gSystem->mkdir(outputDir, kTRUE);

    // preserving the fine binning at low pT but the higher-pT bins (>75 GeV have been adjusted to be slightly wider)
    Double_t ptbins[] = {0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,65,70,75,80,90,100,125,150,1000};
    Int_t nbins = sizeof(ptbins) / sizeof(Double_t) - 1;

    // Instantiate booleans (modified later)
    bool sigIsData = (inputType == 0);
    bool doCrystalBall = false;

    //--------------------------------------------------------------------------------------------------------------
    // Main analysis code
    //==============================================================================================================

    char hname_out[100];
    vector<TH1D*> hPTv;
    vector<TH1D*> hPFu1v, hPFu1Bkgv;
    vector<TH1D*> hPFu2v, hPFu2Bkgv;
    vector<RooDataSet> lDataSetU1;
    vector<RooDataSet> lDataSetU2;

    vector<RooRealVar> vu1Var;
    vector<RooRealVar> vu2Var;

    RooWorkspace pdfsU1("pdfsU1");
    RooWorkspace pdfsU2("pdfsU2");

    TFile *fin = TFile::Open(hist_file);
    cout << ">>> Input file: " << hist_file << endl;

    for (Int_t ibin = 0; ibin < nbins; ibin++) {

        TString uP1_var, uP2_var, bosonpt_var, hnamein_data, hnamein_DY, hnamein_W, hnamein_TT, hnamein_ST, hnamein_VV, hnamein_DYtau, hnamein_Wtau, hnamein_pt;
        if (inputType < 2) { //Z events
            uP1_var = "uP1_uncorrected";
            if (metVar == "pfmet_uncorrected") {
                uP1_var = "pf"+uP1_var;
            }
            bosonpt_var = "bosonpt";

            hnamein_data = TString::Format("data#mm_corr-355862-357482-data#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            cout << "string: " << hnamein_data << endl;
            hnamein_DY = TString::Format("DY#mm_corr-355862-357482-DY#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_W = TString::Format("W#mm_corr-355862-357482-W#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_TT = TString::Format("TT#mm_corr-355862-357482-TT#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_ST = TString::Format("ST#mm_corr-355862-357482-ST#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_VV = TString::Format("VV#mm_corr-355862-357482-VV#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_DYtau = TString::Format("DYtau#mm_corr-355862-357482-DYtau#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_Wtau = TString::Format("Wtau#mm_corr-355862-357482-Wtau#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            if (inputType == 0) {
                hnamein_pt = TString::Format("data#mm_corr-355862-357482-data#Nominal#%s_zpt%d%s", bosonpt_var.Data(), ibin, zrapbin.Data());
            } else if (inputType == 1) {
                hnamein_pt = TString::Format("DY#mm_corr-355862-357482-DY#Nominal#%s_zpt%d%s", bosonpt_var.Data(), ibin, zrapbin.Data());
            }
        } else {
            if (inputType == 2) {
                uP1_var = "uP1_uncorrected_pos";
                bosonpt_var = "bosonpt_pos";
            } else if (inputType == 3) {
                uP1_var = "uP1_uncorrected_neg";
                bosonpt_var = "bosonpt_neg";
            }
            if (metVar == "pfmet_uncorrected") {
                uP1_var = "pf"+uP1_var;
            }
            hnamein_data = TString::Format("data#mmet_corr-355862-357482-data#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_DY = TString::Format("DY#mmet_corr-355862-357482-DY#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_W = TString::Format("W#mmet_corr-355862-357482-W#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_TT = TString::Format("TT#mmet_corr-355862-357482-TT#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_ST = TString::Format("ST#mmet_corr-355862-357482-ST#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_VV = TString::Format("VV#mmet_corr-355862-357482-VV#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_DYtau = TString::Format("DYtau#mmet_corr-355862-357482-DYtau#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_Wtau = TString::Format("Wtau#mmet_corr-355862-357482-Wtau#Nominal#%s_zpt%d%s", uP1_var.Data(), ibin, zrapbin.Data());
            hnamein_pt = TString::Format("W#mmet_corr-355862-357482-W#Nominal#%s_zpt%d%s", bosonpt_var.Data(), ibin, zrapbin.Data());
        }
        TString hname = TString::Format("huP1_%d", ibin);
        TH1D *huP1_data = (TH1D*)fin->Get(hnamein_data)->Clone(hname);
        TH1D *huP1_DY = (TH1D*)fin->Get(hnamein_DY)->Clone(hname+"_DY");
        TH1D *huP1_W = (TH1D*)fin->Get(hnamein_W)->Clone(hname+"_W");
        TH1D *huP1_TT = (TH1D*)fin->Get(hnamein_TT)->Clone(hname+"_TT");
        TH1D *huP1_ST = (TH1D*)fin->Get(hnamein_ST)->Clone(hname+"_ST");
        TH1D *huP1_VV = (TH1D*)fin->Get(hnamein_VV)->Clone(hname+"_VV");
        TH1D *huP1_DYtau = (TH1D*)fin->Get(hnamein_DYtau)->Clone(hname+"_DYtau");
        TH1D *huP1_Wtau = (TH1D*)fin->Get(hnamein_Wtau)->Clone(hname+"_Wtau");
        TH1D *huP1_bkg = (TH1D*)huP1_W->Clone(hname+"_bkg");
        huP1_bkg->Add(huP1_TT);
        huP1_bkg->Add(huP1_ST);
        huP1_bkg->Add(huP1_VV);
        huP1_bkg->Add(huP1_DYtau);
        huP1_bkg->Add(huP1_Wtau);

        if (inputType < 2) { //Z events
            uP2_var = "uP2_uncorrected";
            if (metVar == "pfmet_uncorrected") {
                uP2_var = "pf"+uP2_var;
            }

            hnamein_data = TString::Format("data#mm_corr-355862-357482-data#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_DY = TString::Format("DY#mm_corr-355862-357482-DY#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_W = TString::Format("W#mm_corr-355862-357482-W#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_TT = TString::Format("TT#mm_corr-355862-357482-TT#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_ST = TString::Format("ST#mm_corr-355862-357482-ST#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_VV = TString::Format("VV#mm_corr-355862-357482-VV#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_DYtau = TString::Format("DYtau#mm_corr-355862-357482-DYtau#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_Wtau = TString::Format("Wtau#mm_corr-355862-357482-Wtau#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
        } else {
            if (inputType == 2) {
                uP2_var = "uP2_uncorrected_pos";
            } else if (inputType == 3) {
                uP2_var = "uP2_uncorrected_neg";
            }
            if (metVar == "pfmet_uncorrected") {
                uP2_var = "pf"+uP2_var;
            }
            hnamein_data = TString::Format("data#mmet_corr-355862-357482-data#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_DY = TString::Format("DY#mmet_corr-355862-357482-DY#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_W = TString::Format("W#mmet_corr-355862-357482-W#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_TT = TString::Format("TT#mmet_corr-355862-357482-TT#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_ST = TString::Format("ST#mmet_corr-355862-357482-ST#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_VV = TString::Format("VV#mmet_corr-355862-357482-VV#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_DYtau = TString::Format("DYtau#mmet_corr-355862-357482-DYtau#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
            hnamein_Wtau = TString::Format("Wtau#mmet_corr-355862-357482-Wtau#Nominal#%s_zpt%d%s", uP2_var.Data(), ibin, zrapbin.Data());
        }
        hname = TString::Format("huP2_%d", ibin);
        TH1D *huP2_data = (TH1D*)fin->Get(hnamein_data)->Clone(hname);
        TH1D *huP2_DY = (TH1D*)fin->Get(hnamein_DY)->Clone(hname+"_DY");
        TH1D *huP2_W = (TH1D*)fin->Get(hnamein_W)->Clone(hname+"_W");
        TH1D *huP2_TT = (TH1D*)fin->Get(hnamein_TT)->Clone(hname+"_TT");
        TH1D *huP2_ST = (TH1D*)fin->Get(hnamein_ST)->Clone(hname+"_ST");
        TH1D *huP2_VV = (TH1D*)fin->Get(hnamein_VV)->Clone(hname+"_VV");
        TH1D *huP2_DYtau = (TH1D*)fin->Get(hnamein_DYtau)->Clone(hname+"_DYtau");
        TH1D *huP2_Wtau = (TH1D*)fin->Get(hnamein_Wtau)->Clone(hname+"_Wtau");
        TH1D *huP2_bkg = (TH1D*)huP2_W->Clone(hname+"_bkg");
        huP2_bkg->Add(huP2_TT);
        huP2_bkg->Add(huP2_ST);
        huP2_bkg->Add(huP2_VV);
        huP2_bkg->Add(huP2_DYtau);
        huP2_bkg->Add(huP2_Wtau);

        hname = TString::Format("hPt_%d", ibin);
        TH1D *hPT = (TH1D*)fin->Get(hnamein_pt)->Clone(hname);

        if (inputType == 0) {
            hPFu1v.push_back(huP1_data);
            hPFu2v.push_back(huP2_data);
        } else if (inputType == 1) {
            hPFu1v.push_back(huP1_DY);
            hPFu2v.push_back(huP2_DY);
        } else {
            hPFu1v.push_back(huP1_W);
            hPFu2v.push_back(huP2_W);
        }

        hPFu1Bkgv.push_back(huP1_bkg);
        hPFu2Bkgv.push_back(huP2_bkg);

        hPTv.push_back(hPT);

        std::stringstream name;
        name << "u_" << ibin;

        int range = 250;
        if (ptbins[ibin] > 140)
            range = 500;

        RooRealVar u1Var(name.str().c_str(), name.str().c_str(), 0, -range - ptbins[ibin], range - ptbins[ibin]);
        RooRealVar u2Var(name.str().c_str(), name.str().c_str(), 0, -range, range);
        u1Var.setBins(100);
        u2Var.setBins(100);

        vu1Var.push_back(u1Var);
        vu2Var.push_back(u2Var);

        sprintf(hname_out, "hDataSetU1_%i", ibin);
        RooDataSet dataSetU1(hname_out, hname_out, RooArgSet(u1Var));
        lDataSetU1.push_back(dataSetU1);
        sprintf(hname_out, "hDataSetU2_%i", ibin);
        RooDataSet dataSetU2(hname_out, hname_out, RooArgSet(u2Var));
        lDataSetU2.push_back(dataSetU2);
    }

    cout << ">>> Input histograms loaded" << endl;

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
    plotPFu1mean.Draw(c,kTRUE,"png");

    grPFu2mean = new TGraph(nbins,ptMean,dataMean2);
    grPFu2mean ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu2mean ->GetYaxis()->SetRangeUser(-.1, .1);
    grPFu2mean ->SetName("grPFu2mean");
    CPlot plotPFu2mean("pfu2mean","","Dimuon p_{T} [GeV/c]","<u_{#perp}  >/<p_{T}>");
    plotPFu2mean.AddTextBox( "CMS           Hadronic Recoil Response   (13.0 TeV)", 0.1, 0.92, 0.95, 0.97, 0, kBlack, -1);
    plotPFu2mean.AddGraph(grPFu2mean,"",kBlack,kOpenCircle);
    plotPFu2mean.Draw(c,kTRUE,"png");

    grPFu1meanerr = new TGraph(nbins,ptMean,dataRMS1);
    grPFu1meanerr ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu1meanerr ->GetYaxis()->SetRangeUser(0, 30);
    grPFu1meanerr ->SetName("grPFu1meanerr");
    CPlot plotPFu1meanerr("pfu1meanerr","","Dimuon p_{T} [GeV/c]","RMS(u_{#parallel})");
    plotPFu1meanerr.AddTextBox("CMS     Hadronic Recoil Resolution  (13.0 TeV)", 0.1, 0.92, 0.95, 0.97, 0, kBlack, -1);
    plotPFu1meanerr.AddGraph(grPFu1meanerr,"",kBlack,kOpenCircle);
    plotPFu1meanerr.Draw(c,kTRUE,"png");

    grPFu2meanerr = new TGraph(nbins,ptMean,dataRMS2);
    grPFu2meanerr ->GetXaxis()->SetRangeUser(0., 150.);
    grPFu2meanerr ->GetYaxis()->SetRangeUser(0, 30);
    grPFu2meanerr ->SetName("grPFu2meanerr");
    CPlot plotPFu2meanerr("pfu2meanerr","","Dimuon p_{T} [GeV/c]","RMS(u_{#perp}  )");
    plotPFu2meanerr.AddTextBox("CMS        Hadronic Recoil Resolution  (13.0 TeV)", 0.1, 0.92, 0.95, 0.97, 0, kBlack, -1);
    plotPFu2meanerr.AddGraph(grPFu2meanerr,"",kBlack,kOpenCircle);
    plotPFu2meanerr.Draw(c,kTRUE,"png");

    cout << ">>> Mean Plots Created" << endl;

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

    char outpdfname[100];
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
        if (!sigOnly) {
            hv[ibin]->Add(hbkgv[ibin], -1.0);
        }

        TH1D* hvLOG = (TH1D*)hv[ibin]->Clone();
        for (Int_t i = 0; i < hv[ibin]->GetNbinsX() + 1; i++) {
            if (hv[ibin]->GetBinContent(i) <= 0) {
                hv[ibin]->SetBinContent(i, 0);
                hv[ibin]->SetBinError(i, 0);
                hvLOG->SetBinContent(i, 0);
                hvLOG->SetBinError(i, 0);
            } else {
                hvLOG->SetBinContent(i, TMath::Log10(hv[ibin]->GetBinContent(i)));
                hvLOG->SetBinError(i, 0);
            }
        }

        std::stringstream name;
        // unfortunately have to give each variable individual names for each bin
        name << "u_" << ibin;
        RooRealVar u(name.str().c_str(), name.str().c_str(), hv[ibin]->GetXaxis()->GetXmin(), hv[ibin]->GetXaxis()->GetXmax());
        name.str("");
        u.setBins(100);

        // data hist
        RooDataHist dataHist("dataHist", "dataHist", RooArgSet(u), hv[ibin]);

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
        RooRealVar sigma1(name.str().c_str(), name.str().c_str(), 0.7 * (hv[ibin]->GetRMS()), 0.5 * hv[ibin]->GetRMS(), 2.5 * (hv[ibin]->GetRMS()));
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
        //RooCrystalBall doubleCB(name.str().c_str(), name.str().c_str(), u, mean1, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
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
            shapes.add(gauss1);  //consider swapping the order here
            //shapes.add(doubleCB);
            fracs.add(frac2);
        }

        name.str("");
        name << "sig_" << ibin;
        RooAddPdf sig(name.str().c_str(), name.str().c_str(), shapes, fracs);
        name.str("");

        RooArgList parts;
        if(doCrystalBall) {
            parts.add(sig);
        }
        if(!doCrystalBall) {
            parts.add(sig);
        }

        RooArgList yields;

        name.str("");
        name << "nsig_" << ibin;
        RooRealVar nsig(name.str().c_str(), name.str().c_str(), 0.98 * (hv[ibin]->Integral()), 0., 1.1 * hv[ibin]->Integral()); // just to be sure that doesn't it the boundary
        name.str("");
        name << "nbkg_" << ibin;
        RooRealVar nbkg(name.str().c_str(), name.str().c_str(), (hbkgv[ibin]->Integral()), 0.1 * (hbkgv[ibin]->Integral()), 0.25 * (hv[ibin]->Integral()));
        double bkg_frac = (hv[ibin]->Integral() / (hv[ibin]->Integral() + hbkgv[ibin]->Integral()));
        RooRealVar* lAbkgFrac = new RooRealVar("AbkgFrac", "AbkgFrac", bkg_frac, 0.95*bkg_frac, 1.05*bkg_frac);

        yields.add(nsig);
        nbkg.setVal(0);

        name.str("");
        name << "modelpdf_" << ibin << std::endl;
        RooAddPdf modelpdf(name.str().c_str(), name.str().c_str(), parts, yields);
        name.str("");

        RooFitResult* fitResult = 0;

        if(!doCrystalBall) {
            auto externalConstraintsArgSet = RooArgSet(constGauss1, constGauss2, constGauss3);

            if (model == 2) {
                externalConstraintsArgSet = RooArgSet(constGauss1, constGauss2);
            } else if (model == 1) {
                externalConstraintsArgSet = RooArgSet(constGauss1);
            }


            fitResult = modelpdf.fitTo(dataHist,
                NumCPU(4),
                Minimizer("Minuit2", "minimize"),
                ExternalConstraints(externalConstraintsArgSet),
                RooFit::Minos(),
                RooFit::Strategy(2),
                RooFit::Save());
            int nTries = 0;

            do {
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "scan"),
                    ExternalConstraints(externalConstraintsArgSet),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "migrad"),
                    ExternalConstraints(externalConstraintsArgSet),
                    RooFit::Hesse(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "improve"),
                    ExternalConstraints(externalConstraintsArgSet),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "minimize"),
                    ExternalConstraints(externalConstraintsArgSet),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                nTries++;
            } while ((fitResult->status() > 0 || fitResult->covQual() < 3) && nTries < 5);
        }
        if (doCrystalBall) {
            fitResult = modelpdf.fitTo(dataHist,
                NumCPU(4),
                Minimizer("Minuit2", "minimize"),
                ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                RooFit::Minos(),
                RooFit::Strategy(2),
                RooFit::Save());
            int nTries = 0;

            do {
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "scan"),
                    ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "migrad"),
                    ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                    RooFit::Hesse(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "improve"),
                    ExternalConstraints(RooArgSet(constGauss_forCB, constGauss1)),
                    RooFit::Minos(),
                    RooFit::Strategy(2),
                    RooFit::Save());
                fitResult = modelpdf.fitTo(dataHist,
                    NumCPU(4),
                    Minimizer("Minuit2", "minimize"),
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
            wksp->import(sig);

        mean1Arr[ibin] = mean1.getVal();
        mean1ErrArr[ibin] = mean1.getError();
        sigma1Arr[ibin] = sigma1.getVal();
        sigma1ErrArr[ibin] = sigma1.getError();
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
            name.str("");
            name << "doubleCB_" << ibin;
            sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kRed));

            name.str("");
            name << "gauss1_" << ibin;
            sig.plotOn(frame, Components(name.str().c_str()), LineStyle(kDashed), LineColor(kGreen));
        }

        // draw the curve
        if (!doCrystalBall) {
            sig.plotOn(frame, FillColor(7), VisualizeError(*fitResult, 1), RooFit::Components(sig)); // 1 sigma band
            sig.plotOn(frame, RooFit::LineColor(kBlue));
        } 
        if (doCrystalBall) {
            sig.plotOn(frame, FillColor(7), VisualizeError(*fitResult, 1), RooFit::Components(sig)); // 1 sigma band
            sig.plotOn(frame, RooFit::LineColor(kBlue));
        }  

        // redraw the data
        dataHist.plotOn(frame, MarkerStyle(kFullCircle), MarkerSize(0.8), DrawOption("ZP"));

        if (do_keys) {
            name.str("");
            name << "key_" << ibin;
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
            wksp->import(pdf_keys);
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
        hist_pull->GetXaxis()->SetRangeUser(frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
        hist_pull->GetYaxis()->SetTitle("pull");
        hist_pull->GetYaxis()->SetRangeUser(-5., 5.);
        hist_pull->SetMarkerColor(kAzure);
        hist_pull->SetLineColor(kAzure);
        hist_pull->SetFillColor(kAzure);
        hist_pull->GetYaxis()->SetTitleFont(42);
        hist_pull->GetXaxis()->SetTitleFont(42);
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

        sprintf(lumi,"CMS Preliminary                    %.1f fb^{-1} (13 TeV)", luminosity);
        sprintf(pname, "%sfit_%i", plabel, ibin);
        sprintf(ylabel, "Events / %.1f GeV", hv[ibin]->GetBinWidth(1));
        if (string(plabel) == string("pfu1")) {
            sprintf(xlabel_new, "u_{#parallel}");
        } else {
            sprintf(xlabel_new, "u_{#perp}");
        }
        sprintf(binlabel, "p_{T}(V) = %.1f - %.1f GeV/c", ptbins[ibin], ptbins[ibin + 1]);

        if (etaBinCategory == 1)
            sprintf(binYlabel, "|y| < 0.5");
        if (etaBinCategory == 2)
            sprintf(binYlabel, "0.5 < |y| < 1");
        if (etaBinCategory == 3)
            sprintf(binYlabel, "|y| > 1");

        sprintf(nsigtext, "N_{evts} = %.0f", (double)hv[ibin]->Integral());

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
            sprintf(sig1text, "#sigma_{1} = %.1f #pm %.1f", sigma1Arr[ibin], sigma1ErrArr[ibin]);
            if (model >= 2) {
                sprintf(mean2text, "#mu_{2} = %.1f #pm %.1f", mean2Arr[ibin], mean2ErrArr[ibin]);
                sprintf(sig0text, "#sigma = %.1f #pm %.1f", sigma0Arr[ibin], sigma0ErrArr[ibin]);
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
        plot.AddTextBox(nsigtext, 0.21, 0.78, 0.51, 0.73, 0, kBlack, -1);

        if(!doCrystalBall) {
            if (model == 1)
                plot.AddTextBox(0.70, 0.90, 0.95, 0.80, 0, kBlack, -1, 2, mean1text, sig1text);
            else if (model == 2)
                plot.AddTextBox(0.70, 0.90, 0.95, 0.70, 0, kBlack, -1, 4, mean1text, mean2text, sig1text, sig2text);
            else if (model == 3)
                plot.AddTextBox(0.70, 0.90, 0.95, 0.65, 0, kBlack, -1, 6, mean1text, mean2text, mean3text, sig1text, sig2text, sig3text);  // sig0text
        }
        if(doCrystalBall)
            plot.AddTextBox(0.70, 0.90, 0.95, 0.70, 0, kBlack, -1, 10, mean1text, sig0text, sig1text, sig2text, alphaLtext, alphaRtext, nLtext, nRtext, sig1text, frac2text);
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
