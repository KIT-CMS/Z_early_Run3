import os, sys, glob

base_path_dummy = "/ceph/moh/CROWN_samples/NTUPLE_VERSION/ntuples/2018/"
"""ntuple_dirs = [
    "EarlyRun3_V09",
    "EarlyRun3_V09_CutMuonIso0p15To0p20",
    "EarlyRun3_V09_CutMuonIso0p20To0p25",
    "EarlyRun3_V09_CutMuonIso0p25To0p30",
    "EarlyRun3_V09_CutMuonIso0p30To0p35",
    "EarlyRun3_V09_CutMuonIso0p35To0p40",
    "EarlyRun3_V09_CutMuonIso0p40To0p45",
    "EarlyRun3_V09_CutMuonIso0p45To0p50",
    "EarlyRun3_V09_CutMuonIso0p50To0p55",
    "EarlyRun3_V09_CutMuonIso0p55To0p60",
    "EarlyRun3_V09_CutMuonIso0p60To1p00",
]"""

# ntuple_dirs = [
#     "EarlyRun3_V10",
#     "EarlyRun3_V10_CutMuonIso0p15To0p20",
#     "EarlyRun3_V10_CutMuonIso0p20To0p25",
#     "EarlyRun3_V10_CutMuonIso0p25To0p30",
#     "EarlyRun3_V10_CutMuonIso0p30To0p35",
#     "EarlyRun3_V10_CutMuonIso0p35To0p40",
#     "EarlyRun3_V10_CutMuonIso0p40To0p45",
#     "EarlyRun3_V10_CutMuonIso0p45To0p50",
#     "EarlyRun3_V10_CutMuonIso0p50To0p55",
#     "EarlyRun3_V10_CutMuonIso0p55To0p60",
#     "EarlyRun3_V10_CutMuonIso0p60To1p00",
# ]

ntuple_dirs = ["EarlyRun3_V12"]

channels = [
    "mm",
    "mmet",
    "ee",
    "emet",
]

samples = [
    "EGamma_Run2018*-UL2018",
    "SingleMuon_Run2018*-UL2018",

    #"SingleMuon_Run2018D-UL2018",

    "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",

    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X",
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X",

    "WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
    "WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
    "WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
]

for ntuple_dir in ntuple_dirs:
    base_path = base_path_dummy.replace("NTUPLE_VERSION", ntuple_dir)
    outdir = "./inputs_"+ntuple_dir+"/"
    for sample in samples:
        sample_dir = sample.replace("*", "")
        for ch in channels:
            outdir_full = outdir+sample_dir+"/"+ch+"/"
            if not os.path.exists(outdir_full):
                os.makedirs(outdir_full)

            files = glob.glob(base_path+sample+"/"+ch+"/*.root")

            fout = open(outdir_full+"input_files.txt", "w")
            for file in files:
                fout.write(file+"\n")
            fout.close()