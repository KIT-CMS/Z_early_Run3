import numpy as np
from collections import OrderedDict
from modules.Binnings import mass_bins_w

"""
Systematic uncertainties
"""

syst_groups = {
    "All": None,
    "SFTrk": ["SFTrk"],
    #"QCD_ScaledMC": ["mcScale"],

    "SF": [
        "SFTrk",
        "SFSta",
        "SFID",
        "SFIso",
        "SFTrg",
        "SFPrefire",
    ],
    "PDF": [
        f"LHEPdfWeight{ipdf}" for ipdf in range(1, 101)
    ],
    #"QCDScale": [
    #    "LHEPdfWeightAlphaS",
    #    "LHEScaleWeightMUF",
    #    "LHEScaleWeightMUR",
    #    "LHEScaleWeightMUFMUR",
    #],
    "PartonShower": [
        "PSWeightISR",
        "PSWeightFSR",
        # "PSWeightISRFSR",
    ],
    #"RecoilSyst": [
    #    "RecoilDoubleGauss",
    #    "RecoilSigOnlyFit",
    #    "RecoilZrap0",
    #    "RecoilZrap1",
    #    "RecoilZrap2",
    #],
    #"RecoilZrap": [
    #    "RecoilZrap0",
    #    "RecoilZrap1",
    #    "RecoilZrap2",
    #],
    #"RecoilStat": [
    #    f"RecoilStat{istat}" for istat in range(6)  # range(15)
    #],

    #"QCDSyst": [
    #    "mcScale",
    #    # "Pol1shape",
    #],
    #"QCDStat": [
    #    [
    #        f"QCDp_bin{i}shape" for i in range(1, len(mass_bins_w))
    #    ] +
    #    [
    #        f"QCDm_bin{i}shape" for i in range(1, len(mass_bins_w))
    #    ]
    #],
    "Momentum": [
        "LepCorr",
    ],
    "Pileup": [
        "puWeight",
    ],

}
