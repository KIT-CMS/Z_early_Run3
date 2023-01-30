import sys
import ROOT
import re
import numpy as np
# from CMSPLOTS.myFunction import AddOverflowsTH1, RebinHisto
# from modules.qcdExtrapolater import ExtrapolateQCD
from modules.cardMaker import MakeWJetsCards, MakeZJetsCards, GenerateRunCommand, MakeXSecCard
from modules.histProcessor import ProcessHistsAll
from modules.Binnings import mass_bins_w, mass_bins_w_ptOverMt, mass_bins_z, mass_bins_w_scan
from modules.Systematics import syst_groups
from collections import OrderedDict

ROOT.gROOT.SetBatch(True)


def RunPreparations(f_input, fqcd_input, f_output, lepname = "mu", is5TeV = False, mass_bins_w = mass_bins_w, mass_bins_z = mass_bins_z, outdir_card = "cards", applyLFU = False, fit_variable_w = "pfmet_corr", systs=None):

    ProcessHistsAll(f_input, fqcd_input, f_output, mass_bins_w, mass_bins_z, fit_variable_w)

    card_plus, card_minus, card_z, card_xsec_plus, card_xsec_minus, card_xsec_z = None, None, None, None, None, None

    # generate card based on the signal and qcd templates
    card_plus  = MakeWJetsCards(f_output, f_output, lepname+"plus",  rebinned = False, nMTBins = len(mass_bins_w)-1, is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU, fit_variable_w=fit_variable_w, systs=systs)
    card_minus = MakeWJetsCards(f_output, f_output, lepname+"minus", rebinned = False, nMTBins = len(mass_bins_w)-1, is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU, fit_variable_w=fit_variable_w, systs=systs)
    card_z     = MakeZJetsCards(f_output,           lepname+lepname, rebinned = False,                               is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU,                                systs=systs)

    # generate xsec hists and cards
    card_xsec_plus  = MakeXSecCard(lepname+"plus",  is5TeV, outdir = outdir_card, applyLFU=applyLFU)
    card_xsec_minus = MakeXSecCard(lepname+"minus", is5TeV, outdir = outdir_card, applyLFU=applyLFU)
    card_xsec_z     = MakeXSecCard(lepname+lepname, is5TeV, outdir = outdir_card, applyLFU=applyLFU)

    return card_plus, card_minus, card_z, card_xsec_plus, card_xsec_minus, card_xsec_z

if __name__  == "__main__":
    applyLFU = True

    card_muplus = None
    card_muminus = None
    card_zmumu = None

    version = sys.argv[1]
    doAsimov = ("Asimov" in version)

    w_fit_vars = [
        "pfmt_corr",
        "pfmet_corr",
        # "ptOverPfmt_corr",
    ]

    f_input = "root/earlyRun3_2022_combined.root"
    fqcd_input = "root/QCD/qcdshape_extrapolated_combined.root"

    for fit_variable_w in w_fit_vars:
        print(f">>> fit_variable_w: {fit_variable_w}")
        for idx, _mass_bins_w in mass_bins_w_scan.items():
            _mass_bins_w_toUse = _mass_bins_w
            if fit_variable_w == "ptOverPfmt_corr":
                _mass_bins_w_toUse = mass_bins_w_ptOverMt

            for syst_name, systs in syst_groups.items():
                f_output = f"root/test_{version}/output_shapes_{fit_variable_w}_wbin{idx}_syst{syst_name}.root"

                card_muplus, card_muminus, card_zmumu, card_xsec_muplus, card_xsec_muminus, card_xsec_mumu = RunPreparations(
                    f_input,
                    fqcd_input,
                    f_output,
                    "mu",
                    outdir_card = f"cards/test_{version}/{fit_variable_w}/scan_wbin{idx}_syst{syst_name}",
                    mass_bins_w = _mass_bins_w_toUse,
                    mass_bins_z = mass_bins_z,
                    applyLFU=applyLFU,
                    is5TeV=False,
                    fit_variable_w=fit_variable_w,
                    systs=systs
                )

                # mu channel fit
                GenerateRunCommand(
                    "card_mu",
                    [card_muplus, card_muminus, card_zmumu],
                    ["muplus", "muminus", "mumu"],
                    [card_xsec_muplus, card_xsec_muminus, card_xsec_mumu],
                    applyLFU=applyLFU,
                    doAsimov=doAsimov,
                )