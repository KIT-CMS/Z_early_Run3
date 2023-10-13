#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TTree.h"
#include <sstream>
#include <string>
#include <vector>

#include <math.h>
#include <stdio.h>
//
// ** apply MET corrections **
//
//

using namespace std;

class METXYCorrector {

public:
    METXYCorrector(string iNameZDat, int iSeed = 0xDEADBEEF);
    METXYCorrector(string iNameZDat1, string iPrefix, int iSeed = 0xDEADBEEF);
    void loadXYCorrection(string iNameFile);

    void CorrectMETXY(double& pfmet, double& pfmetphi, int nPV, bool isData);
    void CalcU1U2(double iMet, double iMPhi, double iLepPt, double iLepPhi, double iGenPt, double iGenPhi, double& pU1, double& pU2);

private:
    TF1* f1_metX_vs_nPV_data;
    TF1* f1_metY_vs_nPV_data;
    TF1* f1_metX_vs_nPV_mc;
    TF1* f1_metY_vs_nPV_mc;
};

//-----------------------------------------------------------------------------------------------------------------------------------------
// constructors, but some of the parts are pointless lol

METXYCorrector::METXYCorrector(string iNameZDat, std::string iPrefix, int iSeed)
{
}

METXYCorrector::METXYCorrector(string iNameZ, int iSeed)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------
void METXYCorrector::loadXYCorrection(std::string iFName)
{
    TFile* lFile = new TFile(iFName.c_str());
    f1_metX_vs_nPV_data = (TF1*)lFile->Get("f1_metX_vs_nPV_data");
    f1_metY_vs_nPV_data = (TF1*)lFile->Get("f1_metY_vs_nPV_data");
    f1_metX_vs_nPV_mc = (TF1*)lFile->Get("f1_metX_vs_nPV_mc");
    f1_metY_vs_nPV_mc = (TF1*)lFile->Get("f1_metY_vs_nPV_mc");
    lFile->Delete();
}

void METXYCorrector::CorrectMETXY(double& pfmet, double& pfmetphi, int nPV, bool isData)
{
    double pfmetx = pfmet * TMath::Cos(pfmetphi);
    double pfmety = pfmet * TMath::Sin(pfmetphi);

    nPV = (nPV > 10) ? 10 : nPV;

    //std::cout << "pfmetx " << pfmetx << " pfmety " << pfmety << " nPV " << nPV;

    if (isData) {
        pfmetx -= f1_metX_vs_nPV_data->Eval(nPV);
        pfmety -= f1_metY_vs_nPV_data->Eval(nPV);
        //std::cout << " corrx " << f1_metX_vs_nPV_data->Eval(nPV) << " corry " << f1_metY_vs_nPV_data->Eval(nPV);
    } else {
        pfmetx -= f1_metX_vs_nPV_mc->Eval(nPV);
        pfmety -= f1_metY_vs_nPV_mc->Eval(nPV);
    }

    //std::cout << " after corr pfmetx " << pfmetx << " pfmety " << pfmety << std::endl;

    pfmet = TMath::Sqrt(pfmetx*pfmetx + pfmety*pfmety);
    pfmetphi = TMath::ATan2(pfmety, pfmetx);
}

void METXYCorrector::CalcU1U2(double iMet, double iMPhi, double iLepPt, double iLepPhi, double iGenPt, double iGenPhi, double& pU1, double& pU2)
{
    // calculate the u1 and u2 based on vphi

    double pUX = iMet * cos(iMPhi) + iLepPt * cos(iLepPhi);
    double pUY = iMet * sin(iMPhi) + iLepPt * sin(iLepPhi);
    double pU = sqrt(pUX * pUX + pUY * pUY);
    double pCos = -(pUX * cos(iGenPhi) + pUY * sin(iGenPhi)) / pU;
    double pSin = (pUX * sin(iGenPhi) - pUY * cos(iGenPhi)) / pU;
    pU1 = pU * pCos; // uparallel
    pU2 = pU * pSin; // uperp
}

