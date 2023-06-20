import os, sys
import json
from collections import OrderedDict
import uproot
import numpy as np
import pandas as pd
import tqdm
import gc

def read_ntuples(path, treename = "ntuple"):
    df = uproot.concatenate(f'{path}:{treename}',
                            filter_name = ["genweight", "*pt*", "*eta*", "*phi*", "*mass*", "*pdgId*", "*dR*", "is_*"],
                            file_handler = uproot.MultithreadedFileSource,
                            num_workers = 24,
                            library='pd')
    df.loc[:, "w"] = np.sign(df.genweight)
    df.loc[:, "w2"] = df.w * df.w
    return df

def preselection(df, signal, channel):
    assert signal in ["Z", "Wpos", "Wneg"]
    assert channel in ["mu", "e"]
    if signal == "Z":
        assert len(df.index) == len(df.loc[df.is_dy_mm, :].index) + len(df.loc[df.is_dy_ee, :].index)
    else:
        assert len(df.index) == len(df.loc[df.is_w_m, :].index) + len(df.loc[df.is_w_e, :].index)

    sel = False
    if signal == "Z":
        if channel == "mu":
            sel = (df.is_dy_mm & (df.genlep_pdgId_1 == 13))
        elif channel == "e":
            sel = (df.is_dy_ee & (df.genlep_pdgId_1 == 11))
    elif signal == "Wpos":
        if channel == "mu":
            sel = (df.is_w_m & (df.genlep_pdgId_2 == -13))
        elif channel == "e":
            sel = (df.is_w_e & (df.genlep_pdgId_2 == -11))
    elif signal == "Wneg":
        if channel == "mu":
            sel = (df.is_w_m & (df.genlep_pdgId_1 == 13))
        elif channel == "e":
            sel = (df.is_w_e & (df.genlep_pdgId_1 == 11))

    out = df.loc[sel, :]
    assert len(out.index) > 0
    if signal == "Wpos":
        assert out.genlep_pt_2.min() > 0.
    elif signal == "Wneg":
        assert out.genlep_pt_1.min() > 0.
    else:
        assert out.genlep_pt_1.min() > 0.
        assert out.genlep_pt_2.min() > 0.

    return out

def count(df, signal, prefix = "genlep", lepton_idx = "1", pt_cut = 25., abseta_cut = 2.4):
    num_pw, num_nw, den_pw, den_nw = None, None, None, None
    num_sel = None
    if signal == "Z":
        # always take the true Z mass before FSR
        mass_cut = (
            (df["genlepPreFSR_dilepton_mass"] > 60.) &
            (df["genlepPreFSR_dilepton_mass"] < 120.)
        )
        num_sel = (
            (df[f"{prefix}_pt_1"] > pt_cut) &
            (abs(df[f"{prefix}_eta_1"]) < abseta_cut) &
            (df[f"{prefix}_pt_2"] > pt_cut) &
            (abs(df[f"{prefix}_eta_2"]) < abseta_cut)
        )
        num_pw = len(df.loc[((df.w > 0) & mass_cut & num_sel), :].index)
        num_nw = len(df.loc[((df.w < 0) & mass_cut & num_sel), :].index)
        den_pw = len(df.loc[((df.w > 0) & mass_cut), :].index)
        den_nw = len(df.loc[((df.w < 0) & mass_cut), :].index)
    else:
        num_sel = (
            (df[f"{prefix}_pt_{lepton_idx}"] > pt_cut) &
            (abs(df[f"{prefix}_eta_{lepton_idx}"]) < abseta_cut)
        )
        num_pw = len(df.loc[((df.w > 0) & num_sel), :].index)
        num_nw = len(df.loc[((df.w < 0) & num_sel), :].index)
        den_pw = len(df.loc[(df.w > 0), :].index)
        den_nw = len(df.loc[(df.w < 0), :].index)

    return num_pw, num_nw, den_pw, den_nw

def calc_acc(num_pw, num_nw, den_pw, den_nw):
    # Error propagation: https://arxiv.org/pdf/0908.0130.pdf

    assert den_pw != 0
    assert den_nw != 0
    assert den_pw > den_nw

    a_pw = (num_pw) / den_pw
    v_pw = a_pw*(1.-a_pw)/den_pw
    a_nw = (num_nw) / den_nw
    v_nw = a_nw*(1.-a_nw)/den_nw

    acc = (num_pw - num_nw) / (den_pw - den_nw)
    err = np.sqrt(num_pw**2 * v_pw + num_nw**2 * v_nw) / (den_pw - den_nw)

    return acc, err


if __name__ == '__main__':

    paths = {
        "Z": "/ceph/moh/CROWN_samples/Run3V03_gen/ntuples/2022/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/DYtoLL_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root",
        "Wpos": "/ceph/moh/CROWN_samples/Run3V03_gen/ntuples/2022/WtoLNu_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/WtoLNu_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root",
        "Wneg": "/ceph/moh/CROWN_samples/Run3V03_gen/ntuples/2022/WtoLNu_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X/mm/WtoLNu_NoTau_CP5_13p6TeV_amcatnloFXFX-pythia8-Run3Winter22MiniAOD-122X_*.root",
    }

    channels = {
        "mu": {
            "pt_cut": 25.,
            "abseta_cut": 2.4,
        },
        # "e": {
        #     "pt_cut": 29.,
        #     "abseta_cut": 2.5,
        # }
    }

    types = [
        "genlep",
        "genlepPreFSR",
        "genDressed",
    ]

    infos = OrderedDict({})
    for signal, path in paths.items():
        print(f"\n>>> Start {signal}")
        infos[signal] = OrderedDict({})
        df_tmp = read_ntuples(paths[signal])
        sys.stdout.flush()
        for channel, cuts in channels.items():
            print(f">>>\t channel: {channel}")
            sys.stdout.flush()
            infos[signal][channel] = OrderedDict({})
            infos[signal][channel]["signal"] = signal
            infos[signal][channel]["channel"] = channel
            infos[signal][channel]["df"] = preselection(df_tmp, signal, channel)
            print(">>>\t --> preselection done!")
            sys.stdout.flush()
            lepton_idx = "2" if signal == "Wpos" else "1"
            for t in types:
                print(f">>>\t\t type: {t}")
                num_pw, num_nw, den_pw, den_nw = count(infos[signal][channel]["df"], signal, t, lepton_idx, cuts["pt_cut"], cuts["abseta_cut"])
                acc, err = calc_acc(num_pw, num_nw, den_pw, den_nw)
                infos[signal][channel]["pt_cut"] = cuts["pt_cut"]
                infos[signal][channel]["abseta_cut"] = cuts["abseta_cut"]
                infos[signal][channel][t+"_num_pw"] = num_pw
                infos[signal][channel][t+"_num_nw"] = num_nw
                infos[signal][channel][t+"_den_pw"] = den_pw
                infos[signal][channel][t+"_den_nw"] = den_nw
                infos[signal][channel][t+"_acc"] = acc
                infos[signal][channel][t+"_err"] = err
                print(">>>\t\t --> acc calc done!")
                sys.stdout.flush()

    for signal, path in paths.items():
        for channel, cuts in channels.items():
            del infos[signal][channel]["df"]

    with open('acceptance_new.json', 'w') as f:
        f.write(json.dumps(infos, indent=4))

    print(">>> --> all done!")
    sys.stdout.flush()


