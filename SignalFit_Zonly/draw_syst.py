"""
script to make systmatic tempalate comparisons
"""
import ROOT
import os,sys,math
from collections import OrderedDict
import json
import re
import numpy as np

from modules.Binnings import mass_bins_w

from CMSPLOTS.myFunction import DrawHistos, DrawConfig

ROOT.gROOT.SetBatch(True)
# ROOT.gStyle.SetPalette(77)
palette = ROOT.TColor.GetPalette()

def MakeSystPlot(version: str, ifilename: str, hist_name: str, syst_group: str, systs: list, suffix: str = "", x_label: str = ""):
    """
    compare the systmatic tempalate
    """

    out_dir = f"systs_{version}"
    if not os.path.exists("plots/"+out_dir):
        os.makedirs("plots/"+out_dir)

    # Add nominal
    if "" not in systs:
        systs = systs + [""]

    # Set label
    signal = None
    if "_dy_" in hist_name:
        signal = "Z"
    elif "_w_" in hist_name:
        if "_pos" in hist_name:
            signal = "W^{+}"
        elif "_neg" in hist_name:
            signal = "W^{-}"
    elif "_qcd_" in hist_name:
        if "_pos" in hist_name:
            signal = "QCD #mu^{+}"
        elif "_neg" in hist_name:
            signal = "QCD #mu^{-}"

    print("")
    print("#"*50)
    print("signal:    ", signal)
    print("hist_name: ", hist_name)
    print("syst_group:", syst_group)
    print("systs:     ", systs)
    print("suffix:    ", suffix)
    print("ifile:     ", ifilename)
    print("#"*50)
    

    ifile = ROOT.TFile(ifilename)
    outputname = f"{out_dir}/systs_{version}_{hist_name.replace('h_', '')}_{syst_group}{suffix}"

    hists = []
    drawoptions = []
    legendoptions = []
    ratiooptions = []
    labels = []
    bin_width = None
    nsysts = len(systs)-1
    for idx, syst in enumerate(systs):
        if syst == "":
            hist = ifile.Get(f"{hist_name}")
            hist.SetBinContent(0, hist.GetBinContent(1))
            hist.SetBinError(0, 0)
            hist.SetBinContent(hist.GetNbinsX()+1, hist.GetBinContent(hist.GetNbinsX()))
            hist.SetBinError(hist.GetNbinsX()+1, 0)
            hist.Scale(1.0, "width")
            if "_qcd_" in hist_name:
                # hist.SetFillColorAlpha(ROOT.kBlack, 0.5)
                hist.SetMarkerSize(0)
                hist.SetMarkerColorAlpha(ROOT.kBlack, 0.5)
                hist.SetLineColorAlpha(ROOT.kBlack, 0.5)
                hist.SetLineStyle(1)
                hist.SetLineWidth(1)
            else:
                hist.SetFillColorAlpha(ROOT.kBlack, 0.5)
                hist.SetMarkerSize(0)
                hist.SetMarkerColorAlpha(ROOT.kBlack, 0.5)
                hist.SetLineColorAlpha(0, 0)
                hist.SetLineStyle(1)
                hist.SetLineWidth(0)
            if not bin_width:
                bin_width = hist.GetBinWidth(1)
                assert bin_width == hist.GetBinWidth(2)

            hists.append(hist)
            if "_qcd_" in hist_name:
                drawoptions.append("HIST")
                legendoptions.append("L")
            else:
                drawoptions.append("E2")
                legendoptions.append("F")
            ratiooptions.append("HIST")
            labels.append(f"{signal}: nominal (stat. unc.)")
        else:
            for updn in ["Up", "Down"]:
                icolor = int(idx*(255./nsysts))
                istyle = 1 if updn == "Up" else 2
                hist = ifile.Get(f"{hist_name}_{syst}{updn}")
                hist.Scale(1.0, "width")
                hist.SetLineColor(palette.At(icolor))
                hist.SetLineStyle(istyle)
                hist.SetLineWidth(1)

                hists.append(hist)
                drawoptions.append("HIST")
                legendoptions.append("l")
                ratiooptions.append("HIST")
                labels.append(syst+updn)


    if "m_toFit" in hist_name:
        xlabel = "m_{ll} [GeV]"
    elif "mt" in hist_name:
        xlabel = "m_{T} [GeV]"
    elif "met" in hist_name:
        xlabel = "MET [GeV]"
    if x_label != "":
        xlabel = x_label

    ymaxs = {
        "W^{+}": 3.5e6*1.3,
        "W^{-}": 3.5e6*1.3,
        "Z": 1.0e6*1.3,
        "QCD #mu^{+}": 3.5e5*1.3,
        "QCD #mu^{-}": 3.5e5*1.3,
    }

    yrmin=0.95
    yrmax=1.05
    if "_qcd_" in hist_name:
        yrmin=0.85
        yrmax=1.15
        if "QCDSyst" in syst_group:
            yrmin=0.5
            yrmax=1.5
    if "PDF" in syst_group:
        yrmin=0.995
        yrmax=1.005

    drawconfigs = DrawConfig(
        xmin = hists[0].GetXaxis().GetBinLowEdge(1),
        xmax = hists[0].GetXaxis().GetBinLowEdge(hists[0].GetNbinsX()+1),
        xlabel = xlabel,
        ymin = 0,
        ymax = ymaxs[signal] / bin_width,
        ylabel = f"Events / GeV",
        outputname = outputname,
        noCMS=False,
        dology=False,
        addOverflow=False,
        addUnderflow=False,
        showratio=True,
        ratiobase=len(hists)-1,
        yrmin=yrmin,
        yrmax=yrmax,
        yrlabel = "#scale[0.6]{Syst/Nominal}",
        legendPos=[0.45, 0.72, 0.88, 0.88],
    )
    DrawHistos(
        hists,
        labels,
        drawconfigs.xmin,
        drawconfigs.xmax,
        drawconfigs.xlabel,
        drawconfigs.ymin,
        drawconfigs.ymax,
        drawconfigs.ylabel,
        drawconfigs.outputname,
        dology=drawconfigs.dology,
        dologx=drawconfigs.dologx,
        showratio=drawconfigs.showratio,
        yrmax = drawconfigs.yrmax,
        yrmin = drawconfigs.yrmin,
        yrlabel = drawconfigs.yrlabel,
        donormalize=drawconfigs.donormalize,
        ratiobase=drawconfigs.ratiobase,
        legendNCols = 2,
        legendTextSize = 0.02,
        legendPos = drawconfigs.legendPos,
        redrawihist = drawconfigs.redrawihist,
        extraText = drawconfigs.extraText,
        noCMS = drawconfigs.noCMS,
        addOverflow = drawconfigs.addOverflow,
        addUnderflow = drawconfigs.addUnderflow,
        nMaxDigits = drawconfigs.nMaxDigits,
        hratiopanel=None,
        drawoptions=drawoptions,
        legendoptions=legendoptions,
        ratiooptions=ratiooptions,
        showpull=False,
        hpulls=None,
        W_ref = 600,
        is5TeV = False
    )

    return 1

if __name__  == "__main__":
    version = "v12"
    syst_groups = {
        f"root/{version}/output_shapes_pfmt_corr_wbin0_systAll.root": {
            "h_dy_m_toFit": {
                "SF": [
                    "SFTrk",
                    "SFSta",
                    "SFID",
                    "SFIso",
                    "SFTrg",
                    "SFPrefire",
                ],
                "SFTrk": [
                    "SFTrk"
                ],
                "PDF": [
                    f"LHEPdfWeight{ipdf}" for ipdf in range(1, 101)
                ],
                "QCDScale": [
                    "LHEPdfWeightAlphaS",
                    "LHEScaleWeightMUF",
                    "LHEScaleWeightMUR",
                    "LHEScaleWeightMUFMUR",
                ],
            },

            "h_w_pfmt_corr_pos": {
                "SF": [
                    "SFTrk",
                    "SFSta",
                    "SFID",
                    "SFIso",
                    "SFTrg",
                    "SFPrefire",
                ],
                "SFTrk": [
                    "SFTrk"
                ],
                "PDF": [
                    f"LHEPdfWeight{ipdf}" for ipdf in range(1, 101)
                ],
                "QCDScale": [
                    "LHEPdfWeightAlphaS",
                    "LHEScaleWeightMUF",
                    "LHEScaleWeightMUR",
                    "LHEScaleWeightMUFMUR",
                ],
                "RecoilSyst": [
                    "RecoilDoubleGauss",
                    "RecoilSigOnlyFit",
                    "RecoilZrap0",
                    "RecoilZrap1",
                    "RecoilZrap2",
                ],
                "RecoilZrap": [
                    "RecoilZrap0",
                    "RecoilZrap1",
                    "RecoilZrap2",
                ],
                "RecoilStat": [
                    f"RecoilStat{istat}" for istat in range(6)  # range(15)
                ],
            },

            "h_w_pfmt_corr_neg": {
                "SF": [
                    "SFTrk",
                    "SFSta",
                    "SFID",
                    "SFIso",
                    "SFTrg",
                    "SFPrefire",
                ],
                "SFTrk": [
                    "SFTrk"
                ],
                "PDF": [
                    f"LHEPdfWeight{ipdf}" for ipdf in range(1, 101)
                ],
                "QCDScale": [
                    "LHEPdfWeightAlphaS",
                    "LHEScaleWeightMUF",
                    "LHEScaleWeightMUR",
                    "LHEScaleWeightMUFMUR",
                ],
                "RecoilSyst": [
                    "RecoilDoubleGauss",
                    "RecoilSigOnlyFit",
                    "RecoilZrap0",
                    "RecoilZrap1",
                    "RecoilZrap2",
                ],
                "RecoilZrap": [
                    "RecoilZrap0",
                    "RecoilZrap1",
                    "RecoilZrap2",
                ],
                "RecoilStat": [
                    f"RecoilStat{istat}" for istat in range(6)  # range(15)
                ],
            },

            "h_qcd_pfmt_corr_pos": {
                "QCDSyst": [
                    "mcScale",
                    "Pol1shape",
                ],
                "QCDStat": [
                    f"bin{i}shape" for i in range(1, len(mass_bins_w))
                ],
            },
            
            "h_qcd_pfmt_corr_neg": {
                "QCDSyst": [
                    "mcScale",
                    "Pol1shape",
                ],
                "QCDStat": [
                    f"bin{i}shape" for i in range(1, len(mass_bins_w))
                ],
            },
        },
    }

    for ifile, hist_infos in syst_groups.items():
        for hist_name, hist_info in hist_infos.items():
            for syst_group, systs in hist_info.items():
                MakeSystPlot(version, ifile, hist_name, syst_group, systs)
