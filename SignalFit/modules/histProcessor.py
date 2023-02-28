import ROOT
import re
import numpy as np
import os
from CMSPLOTS.myFunction import AddOverflowsTH1, RebinHisto

ROOT.gROOT.SetBatch(True)

def DoRebin(hist, mass_bins):
    """
    rebin the histogram to be used for the mass fit
    """
    hist = RebinHisto(hist, mass_bins, hist.GetName())

    # clear underflow
    hist.SetBinContent(0, 0.)
    hist.SetBinError(0, 0.)

    return hist


def SaveHistToFile(hist, ofile):
    """
    post processing for histograms
    """
    # hist.SetDirectory(ofile)
    ofile.cd()
    hist.Write()


def mergeEWK(hist_name, hists):
    hout = None
    isFirst = True
    for h in hists:
        if isFirst:
            hout = h.Clone(hist_name)
            isFirst = False
        else:
            hout.Add(h)

    return hout


def renameHist(hist, base_name, syst = None):
    outname = base_name
    if syst is not None and syst != "Nominal":
        outname += "_"+syst
        outname = outname.replace("Dn", "Down")
    out = hist.Clone(outname)
    return out


def ProcessHists(
    ifile, foutput,
    processes = ["EWK", "TT", "DY", "data"],
    ewk = ["Wtau", "DYtau", "VV", "ST"],
    fit_variable = "m_toFit",
    bins = None,
    systs = []):
    """
    process all the hists in ifile and save them into ofile
    """

    finput = ROOT.TFile.Open(ifile)
    keys = finput.GetListOfKeys()
    hnames = []
    for key in keys:
        hname = str(key.GetName())
        if "data" not in hname:
            continue
        if "mm_corr" not in hname and "mmet_corr" not in hname:
            continue
        if "pfmtcut" in hname:
            continue
        if not hname.endswith(fit_variable):
            continue
        # this is illegal...
        if 'LepCorr' in hname:
            continue
        hnames.append(hname.replace("data", "PROCREPLACE"))

    hists_all = {}
    for hname in hnames:
        for syst in ["Nominal"]+systs:
            hname_syst = hname
            if ("met" in fit_variable or "mt" in fit_variable) and syst.startswith("Recoil"):
                # this is too hacky...
                if "RecoilDoubleGauss" in syst:
                    hname_syst = hname_syst\
                        .replace("#pfmet_corr", "#pfmet_corr_double")\
                        .replace("#pfmt_corr", "#pfmt_corr_double")\
                        .replace("#met_corr", "#met_corr_double")\
                        .replace("#mt_corr", "#mt_corr_double")\
                        .replace("#ptOverPfmt_corr", "#ptOverPfmt_corr_double")
                elif "RecoilSigOnlyFit" in syst:
                    hname_syst = hname_syst\
                        .replace("#pfmet_corr", "#pfmet_corr_sigOnly")\
                        .replace("#pfmt_corr", "#pfmt_corr_sigOnly")\
                        .replace("#met_corr", "#met_corr_sigOnly")\
                        .replace("#mt_corr", "#mt_corr_sigOnly")\
                        .replace("#ptOverPfmt_corr", "#ptOverPfmt_corr_sigOnly")
                elif "RecoilZrap0" in syst:
                    hname_syst = hname_syst\
                        .replace("#pfmet_corr", "#pfmet_corr_zrap0")\
                        .replace("#pfmt_corr", "#pfmt_corr_zrap0")\
                        .replace("#met_corr", "#met_corr_zrap0")\
                        .replace("#mt_corr", "#mt_corr_zrap0")\
                        .replace("#ptOverPfmt_corr", "#ptOverPfmt_corr_zrap0")
                elif "RecoilZrap1" in syst:
                    hname_syst = hname_syst\
                        .replace("#pfmet_corr", "#pfmet_corr_zrap1")\
                        .replace("#pfmt_corr", "#pfmt_corr_zrap1")\
                        .replace("#met_corr", "#met_corr_zrap1")\
                        .replace("#mt_corr", "#mt_corr_zrap1")\
                        .replace("#ptOverPfmt_corr", "#ptOverPfmt_corr_zrap1")
                elif "RecoilZrap2" in syst:
                    hname_syst = hname_syst\
                        .replace("#pfmet_corr", "#pfmet_corr_zrap2")\
                        .replace("#pfmt_corr", "#pfmt_corr_zrap2")\
                        .replace("#met_corr", "#met_corr_zrap2")\
                        .replace("#mt_corr", "#mt_corr_zrap2")\
                        .replace("#ptOverPfmt_corr", "#ptOverPfmt_corr_zrap2")
                else:
                    assert syst.startswith("RecoilStat")
                    syst_replace = syst.replace("RecoilStat", "stat")
                    hname_syst = hname_syst\
                        .replace("#pfmet_corr", f"#pfmet_corr_{syst_replace}")\
                        .replace("#pfmt_corr", f"#pfmt_corr_{syst_replace}")\
                        .replace("#met_corr", f"#met_corr_{syst_replace}")\
                        .replace("#mt_corr", f"#mt_corr_{syst_replace}")\
                        .replace("#ptOverPfmt_corr", f"#ptOverPfmt_corr_{syst_replace}")
            else:
                hname_syst = hname_syst.replace("#Nominal#", f"#{syst}#")

            for p in processes:
                if p == "data" and syst != "Nominal":
                    continue
                if p == "EWK" and ("LHEPdfWeight" in syst or "LHEScaleWeight" in syst):
                    continue

                name_base = f"h_{p.lower()}_{fit_variable}"
                h_renamed = None
                if p == "EWK":
                    hists = []
                    for pewk in ewk:
                        h_tmp = finput.Get(hname_syst.replace("PROCREPLACE", pewk))
                        hists.append(h_tmp)
                    h_ewk = mergeEWK("h_ewk", hists)
                    h_renamed = renameHist(h_ewk, name_base, syst)
                else:
                    h_tmp = finput.Get(hname_syst.replace("PROCREPLACE", p))
                    h_renamed = renameHist(h_tmp, name_base, syst)

                if bins is not None:
                    h_renamed = DoRebin(h_renamed, bins)

                if syst in ["RecoilDoubleGauss", "RecoilSigOnlyFit", "RecoilZrap0", "RecoilZrap1", "RecoilZrap2"]:
                    h_renamed_up = h_renamed.Clone(h_renamed.GetName()+"Up")
                    h_renamed_dn = hists_all[h_renamed.GetName().replace(f"_{syst}", "")].Clone(h_renamed.GetName()+"Down")
                    h_renamed_dn.Scale(2.)
                    h_renamed_dn.Add(h_renamed_up, -1.)
                    SaveHistToFile(h_renamed_up, foutput)
                    SaveHistToFile(h_renamed_dn, foutput)
                elif syst in [f"LHEPdfWeight{i}" for i in range(1, 101)]:
                    h_renamed_no = hists_all[h_renamed.GetName().replace(f"_{syst}", "")].Clone(h_renamed.GetName()+"Nominal")
                    n_renamed_no = h_renamed_no.Integral(0, -1)
                    h_renamed_up = h_renamed.Clone(h_renamed.GetName()+"Up")
                    n_renamed_up = h_renamed_up.Integral(0, -1)
                    h_renamed_dn = hists_all[h_renamed.GetName().replace(f"_{syst}", "")].Clone(h_renamed.GetName()+"Down")
                    h_renamed_dn.Scale(2.)
                    h_renamed_dn.Add(h_renamed_up, -1.)
                    n_renamed_dn = h_renamed_dn.Integral(0, -1)
                    #h_renamed_up.Scale(n_renamed_no/n_renamed_up)
                    #h_renamed_dn.Scale(n_renamed_no/n_renamed_dn)
                    SaveHistToFile(h_renamed_up, foutput)
                    SaveHistToFile(h_renamed_dn, foutput)
                    #assert h_renamed_up.Integral(0, -1) - h_renamed_no.Integral(0, -1) < 1.e-6, f"{h_renamed_up.Integral(0, -1)} {h_renamed_no.Integral(0, -1)}"
                    #assert h_renamed_dn.Integral(0, -1) - h_renamed_no.Integral(0, -1) < 1.e-6, f"{h_renamed_dn.Integral(0, -1)} {h_renamed_no.Integral(0, -1)}"
                elif syst in ["LHEPdfWeightAlphaSUp", "LHEPdfWeightAlphaSDown", "LHEScaleWeightMUFUp", "LHEScaleWeightMUFDown", "LHEScaleWeightMURUp", "LHEScaleWeightMURDown", "LHEScaleWeightMUFMURUp", "LHEScaleWeightMUFMURDown"]:
                    h_renamed_no = hists_all[h_renamed.GetName().replace(f"_{syst}", "")].Clone(h_renamed.GetName()+"Nominal")
                    n_renamed_no = h_renamed_no.Integral(0, -1)
                    h_renamed_sy = h_renamed.Clone(h_renamed.GetName())
                    n_renamed_sy = h_renamed_sy.Integral(0, -1)
                    #h_renamed_sy.Scale(n_renamed_no/n_renamed_sy)
                    SaveHistToFile(h_renamed_sy, foutput)
                    #assert h_renamed_sy.Integral(0, -1) - h_renamed_no.Integral(0, -1) < 1.e-6, f"{h_renamed_sy.Integral(0, -1)} {h_renamed_no.Integral(0, -1)}"
                else:
                    hists_all[h_renamed.GetName()] = h_renamed
                    SaveHistToFile(h_renamed, foutput)
    return 1


def ProcessHistsQCD(ifile, foutput, fit_variable = None, bins = None, systs = []):
    finput = ROOT.TFile.Open(ifile)
    keys = finput.GetListOfKeys()
    hnames = []
    for key in keys:
        hname = str(key.GetName())
        if "h_QCD_Extrapolated" not in hname:
            continue
        if not hname.endswith(fit_variable):
            continue
        hnames.append(hname)

    hists_all = {}
    for hname in hnames:
        for syst in ["Nominal"]+systs:
            hname_syst = hname
            if syst != "Nominal":
                hname_syst += "_"+syst

            name_base = f"h_qcd_{fit_variable}"
            h_tmp = finput.Get(hname_syst)
            h_renamed = renameHist(h_tmp, name_base, syst)
            if bins is not None:
                h_renamed = DoRebin(h_renamed, bins)

            if "mcScale" in syst:
                h_renamed_no = hists_all[h_renamed.GetName().replace(f"_{syst}", "")].Clone(h_renamed.GetName()+"Nominal")
                n_renamed_no = h_renamed_no.Integral(0, -1)
                h_renamed_sy = h_renamed.Clone(h_renamed.GetName())
                n_renamed_sy = h_renamed_sy.Integral(0, -1)
                h_renamed_sy.Scale(n_renamed_no/n_renamed_sy)
                SaveHistToFile(h_renamed_sy, foutput)
                assert h_renamed_sy.Integral(0, -1) - h_renamed_no.Integral(0, -1) < 1.e-6, f"{h_renamed_sy.Integral(0, -1)} {h_renamed_no.Integral(0, -1)}"
            else:
                hists_all[h_renamed.GetName()] = h_renamed
                SaveHistToFile(h_renamed, foutput)

    return 1


def ProcessHistsAll(ifile, ifile_qcd, ofile, mass_bins_w, mass_bins_z, fit_variable_w = "pfmet_corr"):
    outdir = ofile.rpartition('/')[0]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    foutput = ROOT.TFile(ofile, "RECREATE")

    # Z histograms
    ProcessHists(
        ifile, foutput,
        processes = ["EWK", "TT", "DY", "data"],
        ewk = ["Wtau", "DYtau", "VV", "ST"],
        fit_variable = "m_toFit",
        bins = mass_bins_z,
        systs = [
            "SFTrkUp",
            "SFTrkDn",
            "SFStaUp",
            "SFStaDn",
            "SFIDUp",
            "SFIDDn",
            "SFIsoUp",
            "SFIsoDn",
            "SFTrgUp",
            "SFTrgDn",
        ] + [
            f"LHEPdfWeight{ipdf}" for ipdf in range(1, 101)
        ] + [
            "LHEPdfWeightAlphaSUp",
            "LHEPdfWeightAlphaSDown",
            "LHEScaleWeightMUFUp",
            "LHEScaleWeightMUFDown",
            "LHEScaleWeightMURUp",
            "LHEScaleWeightMURDown",
            "LHEScaleWeightMUFMURUp",
            "LHEScaleWeightMUFMURDown",
        ] + [
            "LepCorrUp",
            "LepCorrDn",
        ] + [
            "puWeightUp",
            "puWeightDn",
        ]

    )

    # W histograms
    ProcessHists(
        ifile, foutput,
        processes = ["EWK", "TT", "DY", "W", "data"],
        ewk = ["Wtau", "DYtau", "VV", "ST"],
        fit_variable = fit_variable_w+"_pos",
        bins = mass_bins_w,
        systs = [
            "SFTrkUp",
            "SFTrkDn",
            "SFStaUp",
            "SFStaDn",
            "SFIDUp",
            "SFIDDn",
            "SFIsoUp",
            "SFIsoDn",
            "SFTrgUp",
            "SFTrgDn",
            "RecoilDoubleGauss",
            "RecoilSigOnlyFit",
            "RecoilZrap0",
            "RecoilZrap1",
            "RecoilZrap2",
        ] + [
            f"RecoilStat{istat}{updn}" for istat in range(15) for updn in ["Up", "Down"]
        ] + [
            f"LHEPdfWeight{ipdf}" for ipdf in range(1, 101)
        ] + [
            "LHEPdfWeightAlphaSUp",
            "LHEPdfWeightAlphaSDown",
            "LHEScaleWeightMUFUp",
            "LHEScaleWeightMUFDown",
            "LHEScaleWeightMURUp",
            "LHEScaleWeightMURDown",
            "LHEScaleWeightMUFMURUp",
            "LHEScaleWeightMUFMURDown",
        ] + [
            "puWeightUp",
            "puWeightDn",
        ]
    )

    ProcessHists(
        ifile, foutput,
        processes = ["EWK", "TT", "DY", "W", "data"],
        ewk = ["Wtau", "DYtau", "VV", "ST"],
        fit_variable = fit_variable_w+"_neg",
        bins = mass_bins_w,
        systs = [
            "SFTrkUp",
            "SFTrkDn",
            "SFStaUp",
            "SFStaDn",
            "SFIDUp",
            "SFIDDn",
            "SFIsoUp",
            "SFIsoDn",
            "SFTrgUp",
            "SFTrgDn",
            "RecoilDoubleGauss",
            "RecoilSigOnlyFit",
            "RecoilZrap0",
            "RecoilZrap1",
            "RecoilZrap2",
        ] + [
            f"RecoilStat{istat}{updn}" for istat in range(15) for updn in ["Up", "Down"]
        ] + [
            f"LHEPdfWeight{ipdf}" for ipdf in range(1, 101)
        ] + [
            "LHEPdfWeightAlphaSUp",
            "LHEPdfWeightAlphaSDown",
            "LHEScaleWeightMUFUp",
            "LHEScaleWeightMUFDown",
            "LHEScaleWeightMURUp",
            "LHEScaleWeightMURDown",
            "LHEScaleWeightMUFMURUp",
            "LHEScaleWeightMUFMURDown",
        ] + [
            "puWeightUp",
            "puWeightDn",
        ]
    )


    ProcessHistsQCD(
        ifile_qcd, foutput,
        fit_variable = fit_variable_w+"_pos",
        bins = mass_bins_w,
        systs = [
            "mcScaleUp",
            "mcScaleDown",
            "Pol1shapeUp",
            "Pol1shapeDown",
        ] + [
            f"bin{i}shapeUp" for i in range(1, len(mass_bins_w))
        ] + [
            f"bin{i}shapeDown" for i in range(1, len(mass_bins_w))
        ]
    )

    ProcessHistsQCD(
        ifile_qcd, foutput,
        fit_variable = fit_variable_w+"_neg",
        bins = mass_bins_w,
        systs = [
            "mcScaleUp",
            "mcScaleDown",
            "Pol1shapeUp",
            "Pol1shapeDown",
        ] + [
            f"bin{i}shapeUp" for i in range(1, len(mass_bins_w))
        ] + [
            f"bin{i}shapeDown" for i in range(1, len(mass_bins_w))
        ]
    )

    foutput.Close()

    return 1
