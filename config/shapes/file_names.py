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
        "DYJetsToLL_M-10to50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X",
        "DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X",
    ],
    "W": [
        "WtoLNu_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X",
    ],
    "TT": [
        "TTTo2J1L1Nu_CP5_13p6TeV_powheg-pythia8-Run3Winter22MiniAOD-122X",
        "TTTo2L2Nu_CP5_13p6TeV_powheg-pythia8-Run3Winter22MiniAOD-122X",
    ],
    "ST": [
        "TbarBQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8-Run3Winter22MiniAOD-122X",
        "TBbarQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8-Run3Winter22MiniAOD-122X",
        "TbarWplus_DR_AtLeastOneLepton_CP5_13p6TeV_powheg-pythia8-Run3Winter22MiniAOD-122X",
        "TWminus_DR_AtLeastOneLepton_CP5_13p6TeV_powheg-pythia8-Run3Winter22MiniAOD-122X",
    ],
    "VV": [
        "WW_TuneCP5_13p6TeV-pythia8-Run3Winter22MiniAOD-122X",
        "WZ_TuneCP5_13p6TeV-pythia8-Run3Winter22MiniAOD-122X",
        "ZZ_TuneCP5_13p6TeV-pythia8-Run3Winter22MiniAOD-122X",
    ],
    "DYtau": [
        "DYJetsToLL_M-10to50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X",
        "DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X",
    ],
    "Wtau": [
        "WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8-Run3Winter22MiniAOD-122X",
    ],
}

files = {
    "2022": {
        "mm": dict(
            {
                "data": [
                    "SingleMuon_Run2022C-PromptNanoAODv10_v1-v1",
                    # "SingleMuon_Run2022C-PromptReco-v1",
                    # "Muon_Run2022C-PromptReco-v1",
                    "Muon_Run2022C-PromptNanoAODv10_v1-v1",
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
                    "SingleMuon_Run2022C-PromptNanoAODv10_v1-v1",
                    "Muon_Run2022C-PromptNanoAODv10_v1-v1",
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

