"""
script to make postfit comparisons
"""
import sys
from modules.postFitScripts import MakePostPlot, result2json, GetPOIValue, ComparePOIs, DumpGroupImpacts
from modules.CombineHarvester.plotImpacts import plotImpacts
from modules.Binnings import mass_bins_w, mass_bins_w_ptOverMt, mass_bins_z, mass_bins_w_scan
from modules.Systematics import syst_groups
from modules.Utils import FormatTable
import ROOT
import numpy as np
from collections import OrderedDict
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

doInclusive = True

# boolean flag to config if the pull
# distribution should be included in the plots
showPULL = True
if doInclusive:
    sqrtS = "13p6TeV"
    do5TeV = False
    lumi_unc = 0.03

    version = sys.argv[1]
    doAsimov = ("Asimov" in version)

    w_fit_vars = [
        "pfmt_corr",
        # "pfmet_corr",
        # "ptOverPfmt_corr",
    ]

    nscans = len(mass_bins_w_scan)


    for fit_variable_w in w_fit_vars:
        fit_var_short = fit_variable_w.replace("_corr", "")

        starts_mT = []
        vals_lep_pos = []
        vals_lep_neg = []
        vals_leplep  = []
        vals_qcd_pos = []
        vals_qcd_neg = []
        errs_lep_pos = []
        errs_lep_neg = []
        errs_leplep  = []
        errs_qcd_pos = []
        errs_qcd_neg = []

        for idx in range(nscans):

            for syst_name, systs in syst_groups.items():

                workdir = f"cards/{version}/{fit_variable_w}/scan_wbin{idx}_syst{syst_name}/"
                filename = workdir + f"card_mu.root"

                print(filename)

                mass_bins = mass_bins_w_scan[idx]
                if mass_bins[0] > 45:
                    continue

                if w_fit_vars == "ptOverPfmt_corr":
                    mass_bins = mass_bins_w_ptOverMt

                _                = MakePostPlot(filename, "muplus",  "prefit",  mass_bins,   f"{version}_{fit_var_short}_wbin{idx}_syst{syst_name}", showPULL,                                  is5TeV=do5TeV)
                _                = MakePostPlot(filename, "muminus", "prefit",  mass_bins,   f"{version}_{fit_var_short}_wbin{idx}_syst{syst_name}", showPULL, startbin = len(mass_bins),       is5TeV=do5TeV)
                _                = MakePostPlot(filename, "mumu",    "prefit",  mass_bins_z, f"{version}_{fit_var_short}_wbin{idx}_syst{syst_name}", showPULL, startbin = len(mass_bins)*2 - 1, is5TeV=do5TeV)

                if syst_name != "All":
                    continue

                nevts = OrderedDict()
                nevts['muplus']  = MakePostPlot(filename, "muplus",  "postfit", mass_bins,   f"{version}_{fit_var_short}_wbin{idx}_syst{syst_name}", showPULL,                                  is5TeV=do5TeV)
                nevts['muminus'] = MakePostPlot(filename, "muminus", "postfit", mass_bins,   f"{version}_{fit_var_short}_wbin{idx}_syst{syst_name}", showPULL, startbin = len(mass_bins),       is5TeV=do5TeV)
                nevts['mumu']    = MakePostPlot(filename, "mumu",    "postfit", mass_bins_z, f"{version}_{fit_var_short}_wbin{idx}_syst{syst_name}", showPULL, startbin = len(mass_bins)*2 - 1, is5TeV=do5TeV)

                result2json(filename, "lepplus_sig_mu",  f"{workdir}/impacts_lepplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json")
                result2json(filename, "lepminus_sig_mu", f"{workdir}/impacts_lepminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json")
                result2json(filename, "leplep_sig_mu",   f"{workdir}/impacts_leplep_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json")
                result2json(filename, "qcd_muplus_mu",   f"{workdir}/impacts_qcd_muplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json")
                result2json(filename, "qcd_muminus_mu",  f"{workdir}/impacts_qcd_muminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json")
                result2json(filename, "Wplus_sig_sumxsec",  f"{workdir}/impacts_Wplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",  "nuisance_impact_sumpois")
                result2json(filename, "Wminus_sig_sumxsec", f"{workdir}/impacts_Wminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json", "nuisance_impact_sumpois")
                result2json(filename, "Winc_sig_sumxsec",   f"{workdir}/impacts_Winc_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",   "nuisance_impact_sumpois")
                result2json(filename, "Zinc_sig_sumxsec",   f"{workdir}/impacts_Zinc_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",   "nuisance_impact_sumpois")
                result2json(filename, "WplusZRatio_ratiometaratio",  f"{workdir}/impacts_WplusZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",  "nuisance_impact_ratiometapois")
                result2json(filename, "WminusZRatio_ratiometaratio", f"{workdir}/impacts_WminusZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json", "nuisance_impact_ratiometapois")
                result2json(filename, "WZRatio_ratiometaratio",      f"{workdir}/impacts_WZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",      "nuisance_impact_ratiometapois")
                result2json(filename, "WchgRatio_ratiometaratio",    f"{workdir}/impacts_WchgRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",    "nuisance_impact_ratiometapois")
                #result2json(filename, "WchgAsym_chargemetaasym",  f"{workdir}/impacts_WchgAsym_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",  "nuisance_impact_chargemetapois")

                plotImpacts(f"{workdir}/impacts_lepplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",      f"impacts_lepplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_lepminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",     f"impacts_lepminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_leplep_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",       f"impacts_leplep_mll_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_qcd_muplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",   f"impacts_qcd_muplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_qcd_muminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",  f"impacts_qcd_muminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_Wplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",        f"impacts_Wplus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_Wminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",       f"impacts_Wminus_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_Winc_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",         f"impacts_Winc_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_Zinc_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",         f"impacts_Zinc_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_WplusZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",  f"impacts_WplusZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_WminusZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json", f"impacts_WminusZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_WZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",      f"impacts_WZRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                plotImpacts(f"{workdir}/impacts_WchgRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",    f"impacts_WchgRatio_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")
                #plotImpacts(f"{workdir}/impacts_WchgAsym_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}.json",  f"impacts_WchgAsym_{version}_{fit_var_short}_wbin{idx}_syst{syst_name}")

                impacts = OrderedDict()
                _                   = DumpGroupImpacts(filename, "qcd_muplus_mu")
                _                   = DumpGroupImpacts(filename, "qcd_muminus_mu")
                impacts['lepplus']  = DumpGroupImpacts(filename, "lepplus_sig_mu")
                impacts['lepminus'] = DumpGroupImpacts(filename, "lepminus_sig_mu")
                impacts['leplep']   = DumpGroupImpacts(filename, "leplep_sig_mu")
                impacts['Wp']       = DumpGroupImpacts(filename, "Wplus_sig_sumxsec",           "nuisance_group_impact_sumpois")
                impacts['Wm']       = DumpGroupImpacts(filename, "Wminus_sig_sumxsec",          "nuisance_group_impact_sumpois")
                impacts['Winc']     = DumpGroupImpacts(filename, "Winc_sig_sumxsec",            "nuisance_group_impact_sumpois")
                impacts['Zinc']     = DumpGroupImpacts(filename, "Zinc_sig_sumxsec",            "nuisance_group_impact_sumpois")
                impacts['WpOverZ']  = DumpGroupImpacts(filename, "WplusZRatio_ratiometaratio",  "nuisance_group_impact_ratiometapois")
                impacts['WmOverZ']  = DumpGroupImpacts(filename, "WminusZRatio_ratiometaratio", "nuisance_group_impact_ratiometapois")
                impacts['WOverZ']   = DumpGroupImpacts(filename, "WZRatio_ratiometaratio",      "nuisance_group_impact_ratiometapois")
                impacts['WpOverWm'] = DumpGroupImpacts(filename, "WchgRatio_ratiometaratio",    "nuisance_group_impact_ratiometapois")
                impacts['WAsym']    = DumpGroupImpacts(filename, "WchgAsym_chargemetaasym",  "nuisance_group_impact_chargemetapois")

                impacts_subset = OrderedDict()
                impacts_subset['lepplus']  = DumpGroupImpacts(filename, "lepplus_sig_mu")
                impacts_subset['lepminus'] = DumpGroupImpacts(filename, "lepminus_sig_mu")
                impacts_subset['leplep']   = DumpGroupImpacts(filename, "leplep_sig_mu")
                impacts_subset['WpOverZ']  = DumpGroupImpacts(filename, "WplusZRatio_ratiometaratio",  "nuisance_group_impact_ratiometapois")
                impacts_subset['WmOverZ']  = DumpGroupImpacts(filename, "WminusZRatio_ratiometaratio", "nuisance_group_impact_ratiometapois")
                impacts_subset['WpOverWm'] = DumpGroupImpacts(filename, "WchgRatio_ratiometaratio",    "nuisance_group_impact_ratiometapois")

                ## print out the nevts information
                print("\n\n")
                print(FormatTable(nevts, caption=f"Event yield {fit_var_short} wbin{idx}_syst{syst_name}", label = f"tab:yield", simple=True))
                print("\n")
                print(FormatTable(nevts, caption=f"Event yield {fit_var_short} wbin{idx}_syst{syst_name}", label = f"tab:yield"))
                print("\n\n")
                print(FormatTable(impacts, caption=f"Systematic uncertainties in percentage", label = f"tab:impacts", precision=2, simple=True))
                print("\n")
                print(FormatTable(impacts, caption=f"Systematic uncertainties in percentage", label = f"tab:impacts", precision=2))
                print("\n\n")
                print(FormatTable(impacts_subset, caption=f"Systematic uncertainties in percentage", label = f"tab:impacts", precision=2, tdrStyle=True))
                starts_mT.append( mass_bins[0])
                val_lep_pos, err_lep_pos = GetPOIValue(filename, "lepplus_sig_mu")
                val_lep_neg, err_lep_neg = GetPOIValue(filename, "lepminus_sig_mu")
                val_leplep,  err_leplep  = GetPOIValue(filename, "leplep_sig_mu")
                val_qcd_pos, err_qcd_pos = GetPOIValue(filename, "qcd_muplus_mu")
                val_qcd_neg, err_qcd_neg = GetPOIValue(filename, "qcd_muminus_mu")
                vals_lep_pos.append( val_lep_pos )
                errs_lep_pos.append( err_lep_pos )
                vals_lep_neg.append( val_lep_neg )
                errs_lep_neg.append( err_lep_neg )
                vals_leplep .append( val_leplep  )
                errs_leplep .append( err_leplep  )
                vals_qcd_pos.append( val_qcd_pos )
                errs_qcd_pos.append( err_qcd_pos )
                vals_qcd_neg.append( val_qcd_neg )
                errs_qcd_neg.append( err_qcd_neg )

        if nscans > 1:
            # compare POIs
            vals_mT = np.array(starts_mT)
            vals_lep_pos = np.array(vals_lep_pos)
            errs_lep_pos = np.array(errs_lep_pos)
            vals_lep_neg = np.array(vals_lep_neg)
            errs_lep_neg = np.array(errs_lep_neg)
            vals_leplep  = np.array(vals_leplep)
            errs_leplep  = np.array(errs_leplep)
            vals_qcd_pos = np.array(vals_qcd_pos)
            errs_qcd_pos = np.array(errs_qcd_pos)
            vals_qcd_neg = np.array(vals_qcd_neg)
            errs_qcd_neg = np.array(errs_qcd_neg)

            # scale back to 1
            vals_lep_pos = vals_lep_pos / vals_lep_pos[0]
            vals_lep_neg = vals_lep_neg / vals_lep_neg[0]
            vals_leplep  = vals_leplep  / vals_leplep [0]
            vals_qcd_pos = vals_qcd_pos / vals_qcd_pos[0]
            vals_qcd_neg = vals_qcd_neg / vals_qcd_neg[0]

            labels = [
                "W^{+} #rightarrow l^{+}#nu",
                "W^{-} #rightarrow l^{-}#bar{#nu}",
                "Z #rightarrow l^{+}l^{-}",
                "QCD l^{+}",
                "QCD l^{-}"
            ]
            markers = [20, 22, 23, 34, 33]
            colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kCyan+1]
            ComparePOIs(
                vals_mT,
                [vals_lep_pos, vals_lep_neg, vals_leplep, vals_qcd_pos, vals_qcd_neg],
                [errs_lep_pos, errs_lep_neg, errs_leplep, errs_qcd_pos, errs_qcd_neg],
                labels, colors, markers,
                f"poi_vs_wbin_{version}_{fit_var_short}",
                is5TeV = False
            )



    # HERE
    if not doAsimov and False:
        workdir = f"cards/{version}/pfmet_pfmt/"
        filename = workdir + f"card_mu.root"

        mass_bins = mass_bins_w

        nevts = OrderedDict()
        _                      = MakePostPlot(filename, "muplus_pfmet",  "prefit",  mass_bins,   f"{version}_pfmet_pfmt", showPULL, x_label = "MET [GeV]",                                     is5TeV=do5TeV)
        _                      = MakePostPlot(filename, "muminus_pfmet", "prefit",  mass_bins,   f"{version}_pfmet_pfmt", showPULL, x_label = "MET [GeV]",    startbin = len(mass_bins),       is5TeV=do5TeV)
        _                      = MakePostPlot(filename, "muplus_pfmt",   "prefit",  mass_bins,   f"{version}_pfmet_pfmt", showPULL, x_label = "m_{T} [GeV]",  startbin = len(mass_bins)*2 - 1, is5TeV=do5TeV)
        _                      = MakePostPlot(filename, "muminus_pfmt",  "prefit",  mass_bins,   f"{version}_pfmet_pfmt", showPULL, x_label = "m_{T} [GeV]",  startbin = len(mass_bins)*3 - 2, is5TeV=do5TeV)
        _                      = MakePostPlot(filename, "mumu",          "prefit",  mass_bins_z, f"{version}_pfmet_pfmt", showPULL, x_label = "m_{ll} [GeV]", startbin = len(mass_bins)*4 - 3, is5TeV=do5TeV)
        nevts['muplus_pfmet']  = MakePostPlot(filename, "muplus_pfmet",  "postfit", mass_bins,   f"{version}_pfmet_pfmt", showPULL, x_label = "MET [GeV]",                                     is5TeV=do5TeV)
        nevts['muminus_pfmet'] = MakePostPlot(filename, "muminus_pfmet", "postfit", mass_bins,   f"{version}_pfmet_pfmt", showPULL, x_label = "MET [GeV]",    startbin = len(mass_bins),       is5TeV=do5TeV)
        nevts['muplus_pfmt']   = MakePostPlot(filename, "muplus_pfmt",   "postfit", mass_bins,   f"{version}_pfmet_pfmt", showPULL, x_label = "m_{T} [GeV]",  startbin = len(mass_bins)*2 - 1, is5TeV=do5TeV)
        nevts['muminus_pfmt']  = MakePostPlot(filename, "muminus_pfmt",  "postfit", mass_bins,   f"{version}_pfmet_pfmt", showPULL, x_label = "m_{T} [GeV]",  startbin = len(mass_bins)*3 - 2, is5TeV=do5TeV)
        nevts['mumu']          = MakePostPlot(filename, "mumu",          "postfit", mass_bins_z, f"{version}_pfmet_pfmt", showPULL, x_label = "m_{ll} [GeV]", startbin = len(mass_bins)*4 - 3, is5TeV=do5TeV)

        result2json(filename, "lepplus_sig_mu",  f"{workdir}/impacts_lepplus_{version}_pfmet_pfmt.json")
        result2json(filename, "lepminus_sig_mu", f"{workdir}/impacts_lepminus_{version}_pfmet_pfmt.json")
        result2json(filename, "leplep_sig_mu",   f"{workdir}/impacts_leplep_{version}_pfmet_pfmt.json")
        result2json(filename, "qcd_muplus_mu",   f"{workdir}/impacts_qcd_muplus_{version}_pfmet_pfmt.json")
        result2json(filename, "qcd_muminus_mu",  f"{workdir}/impacts_qcd_muminus_{version}_pfmet_pfmt.json")
        result2json(filename, "Wplus_sig_sumxsec",  f"{workdir}/impacts_Wplus_{version}_pfmet_pfmt.json",  "nuisance_impact_sumpois")
        result2json(filename, "Wminus_sig_sumxsec", f"{workdir}/impacts_Wminus_{version}_pfmet_pfmt.json", "nuisance_impact_sumpois")
        result2json(filename, "Winc_sig_sumxsec",   f"{workdir}/impacts_Winc_{version}_pfmet_pfmt.json",   "nuisance_impact_sumpois")
        result2json(filename, "Zinc_sig_sumxsec",   f"{workdir}/impacts_Zinc_{version}_pfmet_pfmt.json",   "nuisance_impact_sumpois")
        result2json(filename, "WplusZRatio_ratiometaratio",  f"{workdir}/impacts_WplusZRatio_{version}_pfmet_pfmt.json",  "nuisance_impact_ratiometapois")
        result2json(filename, "WminusZRatio_ratiometaratio", f"{workdir}/impacts_WminusZRatio_{version}_pfmet_pfmt.json", "nuisance_impact_ratiometapois")
        result2json(filename, "WZRatio_ratiometaratio",      f"{workdir}/impacts_WZRatio_{version}_pfmet_pfmt.json",      "nuisance_impact_ratiometapois")
        result2json(filename, "WchgRatio_ratiometaratio",    f"{workdir}/impacts_WchgRatio_{version}_pfmet_pfmt.json",    "nuisance_impact_ratiometapois")
        #result2json(filename, "WchgAsym_chargemetaasym",  f"{workdir}/impacts_WchgAsym_{version}_pfmet_pfmt.json",  "nuisance_impact_chargemetapois")

        plotImpacts(f"{workdir}/impacts_lepplus_{version}_pfmet_pfmt.json",      f"impacts_lepplus_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_lepminus_{version}_pfmet_pfmt.json",     f"impacts_lepminus_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_leplep_{version}_pfmet_pfmt.json",       f"impacts_leplep_mll_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_qcd_muplus_{version}_pfmet_pfmt.json",   f"impacts_qcd_muplus_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_qcd_muminus_{version}_pfmet_pfmt.json",  f"impacts_qcd_muminus_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_Wplus_{version}_pfmet_pfmt.json",        f"impacts_Wplus_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_Wminus_{version}_pfmet_pfmt.json",       f"impacts_Wminus_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_Winc_{version}_pfmet_pfmt.json",         f"impacts_Winc_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_Zinc_{version}_pfmet_pfmt.json",         f"impacts_Zinc_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_WplusZRatio_{version}_pfmet_pfmt.json",  f"impacts_WplusZRatio_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_WminusZRatio_{version}_pfmet_pfmt.json", f"impacts_WminusZRatio_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_WZRatio_{version}_pfmet_pfmt.json",      f"impacts_WZRatio_{version}_pfmet_pfmt")
        plotImpacts(f"{workdir}/impacts_WchgRatio_{version}_pfmet_pfmt.json",    f"impacts_WchgRatio_{version}_pfmet_pfmt")
        #plotImpacts(f"{workdir}/impacts_WchgAsym_{version}_pfmet_pfmt.json",  f"impacts_WchgAsym_{version}_pfmet_pfmt")

        impacts = OrderedDict()
        _                   = DumpGroupImpacts(filename, "qcd_muplus_mu")
        _                   = DumpGroupImpacts(filename, "qcd_muminus_mu")
        impacts['lepplus']  = DumpGroupImpacts(filename, "lepplus_sig_mu")
        impacts['lepminus'] = DumpGroupImpacts(filename, "lepminus_sig_mu")
        impacts['leplep']   = DumpGroupImpacts(filename, "leplep_sig_mu")
        impacts['Wp']       = DumpGroupImpacts(filename, "Wplus_sig_sumxsec",           "nuisance_group_impact_sumpois")
        impacts['Wm']       = DumpGroupImpacts(filename, "Wminus_sig_sumxsec",          "nuisance_group_impact_sumpois")
        impacts['Winc']     = DumpGroupImpacts(filename, "Winc_sig_sumxsec",            "nuisance_group_impact_sumpois")
        impacts['Zinc']     = DumpGroupImpacts(filename, "Zinc_sig_sumxsec",            "nuisance_group_impact_sumpois")
        impacts['WpOverZ']  = DumpGroupImpacts(filename, "WplusZRatio_ratiometaratio",  "nuisance_group_impact_ratiometapois")
        impacts['WmOverZ']  = DumpGroupImpacts(filename, "WminusZRatio_ratiometaratio", "nuisance_group_impact_ratiometapois")
        impacts['WOverZ']   = DumpGroupImpacts(filename, "WZRatio_ratiometaratio",      "nuisance_group_impact_ratiometapois")
        impacts['WpOverWm'] = DumpGroupImpacts(filename, "WchgRatio_ratiometaratio",    "nuisance_group_impact_ratiometapois")
        #impacts['WAsym']    = DumpGroupImpacts(filename, "WchgAsym_chargemetaasym",  "nuisance_group_impact_chargemetapois")

        impacts_subset = OrderedDict()
        impacts_subset['lepplus']  = DumpGroupImpacts(filename, "lepplus_sig_mu")
        impacts_subset['lepminus'] = DumpGroupImpacts(filename, "lepminus_sig_mu")
        impacts_subset['leplep']   = DumpGroupImpacts(filename, "leplep_sig_mu")
        impacts_subset['WpOverZ']  = DumpGroupImpacts(filename, "WplusZRatio_ratiometaratio",  "nuisance_group_impact_ratiometapois")
        impacts_subset['WmOverZ']  = DumpGroupImpacts(filename, "WminusZRatio_ratiometaratio", "nuisance_group_impact_ratiometapois")
        impacts_subset['WpOverWm'] = DumpGroupImpacts(filename, "WchgRatio_ratiometaratio",    "nuisance_group_impact_ratiometapois")

        ## print out the nevts information
        print("\n\n")
        print(FormatTable(nevts, caption=f"Event yield pfmet_pfmt", label = f"tab:yield", simple=True))
        print("\n")
        print(FormatTable(nevts, caption=f"Event yield pfmet_pfmt", label = f"tab:yield"))
        print("\n\n")
        print(FormatTable(impacts, caption=f"Systematic uncertainties in percentage", label = f"tab:impacts", precision=2, simple=True))
        print("\n")
        print(FormatTable(impacts, caption=f"Systematic uncertainties in percentage", label = f"tab:impacts", precision=2))
        print("\n\n")
        print(FormatTable(impacts_subset, caption=f"Systematic uncertainties in percentage", label = f"tab:impacts", precision=2, tdrStyle=True))
