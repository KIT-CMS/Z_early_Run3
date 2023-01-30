#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from TagAndProbeFitter import TagAndProbeFitter
import ROOT
ROOT.gROOT.SetBatch()
# kInfo = 1001, kWarning = 2001, ...
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")
ROOT.gROOT.LoadMacro('RooCMSShape.cc+')

def hist_fitter(
    inFName,
    run_tag,
    lumi_online,
    doTemplate = True
):
    tnpNomFitSig = [
        # "meanP[-0.0, -5.0, 5.0]",
        # "sigmaLP[0.9, 0.05, 5.0]",
        # "sigmaRP[0.9, 0.05, 5.0]",
        
        # "meanF[-0.0, -5.0, 5.0]",
        # "sigmaF[0.9, 0.05, 5.0]",

        # "RooCrystalBall::sigResPass()",

        "meanP[-0.0, -5.0, 5.0]", "sigmaP[0.9, 0.05, 5.0]",
        "meanF[-0.0, -5.0, 5.0]", "sigmaF[0.9, 0.05, 5.0]",
        "Gaussian::sigResPass(x, meanP, sigmaP)",
        "Gaussian::sigResFail(x, meanF, sigmaF)",
    ]

    tnpNomFitBkg = [
        "acmsP[60., 50., 190.]", "betaP[0.05, 0.01, 0.08]",
        "gammaP[0.1, -2, 2]", "peakP[91.2]",
        "acmsF[60., 50., 190.]", "betaF[0.05, 0.01, 0.08]",
        "gammaF[0.1, -2, 2]", "peakF[91.2]",
        "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
    ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpNomFitSig)
    tnpWorkspace.extend(tnpNomFitBkg)

    run_tag_name = "-{}-data".format(run_tag) if run_tag != None else ""
    infile = ROOT.TFile(inFName, "read")
    hP = infile.Get('data#mm-HLT2{}#Nominal#m_toFit'.format(run_tag_name))
    hF = infile.Get('data#mm-HLT1{}#Nominal#m_toFit'.format(run_tag_name))
    hP.Scale(lumi_online)
    hF.Scale(lumi_online)
    hP.Print()
    hF.Print()

    fitter = TagAndProbeFitter("test")
    fitter.set_histograms(hP, hF)
    fitter.set_fit_range(60, 120)

    histZLineShapeP = infile.Get('DY#mm-HLT2-DY#Nominal#m_toFit')
    histZLineShapeF = infile.Get('DY#mm-HLT1-DY#Nominal#m_toFit')
    histZLineShapeP.Scale(lumi_online)
    histZLineShapeF.Scale(lumi_online)
    fitter.set_gen_shapes(histZLineShapeP, histZLineShapeF)
    histZLineShapeP.Print()
    histZLineShapeF.Print()

    histZMC = infile.Get('DY#mm-DY#Nominal#m_toFit')
    histZMC.Scale(lumi_online)
    nZ_MC = histZMC.Integral(0, -1)
    histZMC.Print()

    infile.Close()

    # set workspace
    fitter.set_workspace(tnpWorkspace, doTemplate)

    # fit
    outFName = "./fits/fitResults_all.root"
    if run_tag != None:
        outFName = "./fits/fitResults_{}.root".format(run_tag.replace("-","_"))
    os.makedirs(os.path.dirname(outFName), exist_ok=True)

    lumi_Z, lumi_Z_err = fitter.fit(outFName, lumi_online, nZ_MC)

    return lumi_Z, lumi_Z_err

if __name__ == "__main__":
    output = {}
    inFName = "../output/earlyRun3_2022_mm_runPlot0_PUPPI.root"
    lumi_Z, lumi_Z_err = hist_fitter(inFName, None, 1.)
    output["all"] = [lumi_Z, lumi_Z_err]

    # inFName = "../output/earlyRun3_2022_mm_runPlot1_PUPPI.root"
    # runs = {
    #     "355872-356075": 104.472137129,
    #     "356076-356323": 166.562837288,
    #     "356371-356378": 49.549307419,
    #     "356381": 131.203129037,
    #     "356383-356433": 110.396733126,
    #     "356434-356446": 94.701628537
    # }

    inFName = "../output/earlyRun3_2022_mm_runPlot1_each.root"
    run_list = [355872,355892,355912,355913,355921,355933,355942,355988,355989,356043,356071,356074,356075,356076,356077,356135,356309,356316,356321,356322,356323,356371,356375,356378,356381,356383,356385,356386,356426,356428,356433,356434,356435,356446]
    lumi_list = [21.029407135,1.976640736,9.007156928,5.725368440,22.873241336,0.724271840,1.371972162,2.780268121,1.265671750,5.567402786,15.248489527,2.999361674,13.902884694,17.133845292,49.787188261,2.612524012,10.957679862,9.292921771,7.684875420,2.060615508,67.033187162,1.569606629,6.515209807,41.464490983,131.203129037,2.867059944,2.467573112,8.706233768,2.563759568,45.822680251,47.969426483,1.724343482,0.441028997,92.536256058]

    # for run, run_lumi in runs.items():
    for irun, run in enumerate(run_list):
        run_lumi = lumi_list[irun]
        lumi_Z, lumi_Z_err = hist_fitter(inFName, str(run), run_lumi)
        output[run] = [lumi_Z, lumi_Z_err]

    import json
    j = json.dumps(output, indent=4)
    print(j)
    with open('lumi_Z.json', 'w') as fp:
        fp.write(j)
        fp.close()

    sys.exit(0)
