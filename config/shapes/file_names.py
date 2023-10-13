# 2018 test files
common_files_2018 = {
    "DY": [
        "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
    ],
    "TT": [
        "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X",
        "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X",
    ],
    "W": [
        "WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
        "WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
        "WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9-106X",
    ],
}

common_files_2022 = {
    "DY": [
        "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv11-126X",
        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv11-126X",
    ],
    "W": [
        "WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv11-126X",
    ],
    "TT": [
        "TT_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv10-124X",
    ],
    "ST": [
        "TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8_Run3Summer22NanoAODv11-126X",
        "TbarWplus_DR_AtLeastOneLepton_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv11-126X",
        "TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8_Run3Summer22NanoAODv11-126X",
        "TWminus_DR_AtLeastOneLepton_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer22NanoAODv11-126X",
    ],
    "VV": [
        "WW_TuneCP5_13p6TeV_pythia8_Run3Summer22NanoAODv11-126X",
        "WZ_TuneCP5_13p6TeV_pythia8_Run3Summer22NanoAODv11-126X",
        "ZZ_TuneCP5_13p6TeV_pythia8_Run3Summer22NanoAODv11-126X",
        "EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8_RunIIAutumn18NanoAODv5-Nano1June2019",
    ],
    "DYtau": [
        "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv11-126X",
        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv11-126X",
    ],
    "Wtau": [
        "WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer22NanoAODv11-126X",
    ],
}

files = {
    "2022": {
        "mm": dict(
            {
                "data": [
                    "SingleMuon_Run2022C-PromptNanoAODv10",
                    "Muon_Run2022C-PromptNanoAODv10",
                    # "Muon_Run2022D-PromptReco-v1",
                    # "Muon_Run2022D-PromptReco-v2",
                ],
            },
            **common_files_2022
        ),
        "ee": dict(
            {
                "data": [
                    "EGamma_Run2022C-PromptReco-v1",
                    "EGamma_Run2022D-PromptReco-v1",
                    "EGamma_Run2022D-PromptReco-v2",
                ],
            },
            **common_files_2022
        ),
        "mmet": dict(
            {
                "data": [
                    # "SingleMuon_Run2022C-PromptReco-v1",
                    "SingleMuon_Run2022C-PromptNanoAODv10",
                    "Muon_Run2022C-PromptNanoAODv10",
                    # "SingleMuon_Run2022C-PromptReco-v1",
                    # "Muon_Run2022D-PromptReco-v1",
                    # "Muon_Run2022D-PromptReco-v2",
                ],
            },
            **common_files_2022
        ),
        "emet": dict(
            {
                "data": [
                    "EGamma_Run2022C-PromptReco-v1",
                    "EGamma_Run2022D-PromptReco-v1",
                    "EGamma_Run2022D-PromptReco-v2",
                ],
            },
            **common_files_2022
        ),
    },

    "2018": {
        "mm": dict(
            {
                "data": [
                    # "SingleMuon_Run2018A-UL2018",
                    # "SingleMuon_Run2018B-UL2018",
                    # "SingleMuon_Run2018C-UL2018",
                    "SingleMuon_Run2018D-UL2018",
                ],
            },
            **common_files_2018
        ),
        "ee": dict(
            {
                "data": [
                    # "EGamma_Run2018A-UL2018",
                    # "EGamma_Run2018B-UL2018",
                    # "EGamma_Run2018C-UL2018",
                    "EGamma_Run2018D-UL2018",
                ],
            },
            **common_files_2018
        ),
        "mmet": dict(
            {
                "data": [
                    # "SingleMuon_Run2018A-UL2018",
                    # "SingleMuon_Run2018B-UL2018",
                    # "SingleMuon_Run2018C-UL2018",
                    "SingleMuon_Run2018D-UL2018",
                ],
            },
            **common_files_2018
        ),
        "emet": dict(
            {
                "data": [
                    # "EGamma_Run2018A-UL2018",
                    # "EGamma_Run2018B-UL2018",
                    # "EGamma_Run2018C-UL2018",
                    "EGamma_Run2018D-UL2018",
                ],
            },
            **common_files_2018
        ),
    },
}

